using Downloads, HTTP, SHA, Printf
using DataFrames, Dates, Statistics, Plots, LinearAlgebra, Distributions
using LocalProjections

const GDRIVE_FILE_ID  = "1MF7Am-3utva0NDH0-FE_cvAQOiPN0372"
const EXPECTED_SHA256 = "6dd2f041036209a7abfee1fe29430c97e685478d073178ab25f5b57951b5cb34"
const LOCAL_PATH      = "cps_00010.dat"

# Download large file from Google Drive
function download_gdrive_large(file_id::String, dest::String;
                                expected_hash::Union{String,Nothing}=nothing)
    if isfile(dest)
        first_bytes = String(read(dest, min(20, filesize(dest))))
        if startswith(first_bytes, "<!") || startswith(first_bytes, "<h")
            println("Existing file is HTML (bad download) -- deleting and retrying.")
            rm(dest)
        else
            println("File already exists: $dest -- skipping download.")
        end
    end

    if !isfile(dest)
        base_url = "https://drive.google.com/uc?export=download&id=$(file_id)"
        resp     = HTTP.get(base_url)
        body     = String(resp.body)
        uuid_m   = match(r"name=\"uuid\" value=\"([^\"]+)\"", body)
        real_url = if !isnothing(uuid_m)
            "https://drive.usercontent.google.com/download?id=$(file_id)&export=download&confirm=t&uuid=$(uuid_m.captures[1])"
        else
            "https://drive.usercontent.google.com/download?id=$(file_id)&export=download&confirm=t"
        end
        Downloads.download(real_url, dest)
    end

    if !isnothing(expected_hash)
        actual = bytes2hex(open(sha256, dest))
        actual != expected_hash && error("Hash mismatch!\n  expected: $expected_hash\n  actual:   $actual")
    end
end

# Occupation classification (OCC1990 -> 9 categories)
const OCC_LABELS = Dict(
    1 => "Managerial",
    2 => "Professional_specialty",
    3 => "High_tech",
    4 => "Sales",
    5 => "Administrative_support",
    6 => "Service",
    7 => "Farming_forestry_construction",
    8 => "Precision_production_repair",
    9 => "Machine_operators_transport",
)

function classify_occ1990(occ::Int)::Union{Int,Missing}
    occ in 3:37                                           && return 1
    occ in 43:200                                         && return 2
    occ in 203:235                                        && return 3
    occ in 243:283                                        && return 4
    occ in 303:389                                        && return 5
    occ in 405:469                                        && return 6
    (occ in 473:498 || occ in 558:599 || occ in 614:617) && return 7
    (occ in 503:549 || occ in 628:699)                   && return 8
    (occ in 703:799 || occ in 803:889)                   && return 9
    return missing
end

# Parse fixed-width CPS data
function parse_cps(path::String)::DataFrame
    rows = Vector{NamedTuple}()

    open(path, "r") do io
        for line in eachline(io)
            length(line) < 46 && continue

            year_s    = strip(line[1:4])
            month_s   = strip(line[5:6])
            wtfinl_s  = strip(line[7:20])
            empstat_s = strip(line[21:22])
            occ_s     = strip(line[23:25])
            uhrs_s    = strip(line[26:28])
            earnwt_s  = strip(line[29:38])
            earnwk_s  = strip(line[39:46])

            any(isempty, (year_s, month_s, wtfinl_s, empstat_s, earnwt_s)) && continue

            empstat = tryparse(Int, empstat_s); isnothing(empstat) && continue
            empstat in (10, 12, 20, 21, 22)   || continue

            yr  = parse(Int, year_s)
            mth = parse(Int, month_s)
            yr < 1994 || yr > 2023 || (yr == 2023 && mth > 3) && continue

            occ    = tryparse(Int,     occ_s);    isnothing(occ)    && continue
            wtfinl = tryparse(Float64, wtfinl_s); isnothing(wtfinl) && continue
            earnwt = tryparse(Float64, earnwt_s); isnothing(earnwt) && continue

            wtfinl /= 10_000.0
            earnwt /= 10_000.0
            wtfinl > 0 || continue

            occ_cat = classify_occ1990(occ)
            ismissing(occ_cat) && continue

            is_employed = empstat in (10, 12)
            uhrs        = missing
            earnweek    = missing

            if is_employed && earnwt > 0
                isempty(uhrs_s) || isempty(earnwk_s) && continue
                u = tryparse(Int,     uhrs_s)
                e = tryparse(Float64, earnwk_s)
                if !isnothing(u) && !isnothing(e)
                    e /= 100.0
                    if e > 0 && u > 0 && u <= 168
                        uhrs     = Float64(u)
                        earnweek = e
                    end
                end
            end

            push!(rows, (
                year        = yr,
                month       = mth,
                occ_cat     = occ_cat,
                wtfinl      = wtfinl,
                earnwt      = earnwt,
                is_employed = is_employed,
                uhrs        = uhrs,
                earnweek    = earnweek,
            ))
        end
    end

    return DataFrame(rows)
end

# Aggregate to quarterly panel
function aggregate_quarterly(df::DataFrame)::DataFrame
    df[!, :quarter] = df.year .+ div.((df.month .- 1), 3) ./ 4

    agg = combine(
        groupby(df, [:quarter, :occ_cat]),
        [:earnweek, :uhrs, :wtfinl, :earnwt, :is_employed] => (
            (ew, uh, wf, wt, emp) -> begin
                W_person   = sum(wf)
                emp_wf     = sum(Float64.(emp) .* wf)
                unemp_rate = W_person > 0 ? 1.0 - emp_wf / W_person : missing

                emp_mask = emp .== true
                has_earn = emp_mask .& .!ismissing.(ew) .& .!ismissing.(uh)
                W_earn   = sum(wt[has_earn])

                if W_earn > 0
                    mean_income = log(sum(skipmissing(ew[has_earn]) .* wt[has_earn]) / W_earn)
                    mean_hours  = log(sum(skipmissing(Float64.(uh[has_earn])) .* wt[has_earn]) / W_earn)
                    mean_wage   = mean_income - mean_hours
                else
                    mean_income = missing
                    mean_hours  = missing
                    mean_wage   = missing
                end

                return (mean_income = mean_income,
                        mean_hours  = mean_hours,
                        mean_wage   = mean_wage,
                        unemp_rate  = unemp_rate)
            end
        ) => AsTable
    )

    sort!(agg, [:occ_cat, :quarter])
    return agg
end

# Moving average smoother
function moving_avg(x::AbstractVector, n::Int)
    out = Vector{Union{Missing,Float64}}(missing, length(x))
    for i in n:length(x)
        window = x[i-n+1:i]
        if !any(ismissing, window)
            out[i] = mean(skipmissing(window))
        end
    end
    return out
end

# Bias correction delta
function lp_biascorr_delta(irs::AbstractVector, w::AbstractMatrix)::Vector{Float64}
    H     = length(irs)
    T, k  = size(w)
    delta = zeros(H)

    for j in 1:k
        wj      = w[:, j]
        y_ar    = wj[2:end]
        x_ar    = [ones(T-1) wj[1:end-1]]
        rho_hat = (x_ar' * x_ar) \ (x_ar' * y_ar)
        rho     = clamp(rho_hat[2], -0.99, 0.99)

        bias0 = (1 + rho) / (T * (1 - rho))
        for h in 1:H
            delta[h] += rho^h * bias0
        end
    end

    return delta
end

# Run local projection IRF
function run_lp_irf(df_merged::DataFrame, outcome::Symbol;
                    nlag::Int             = 2,
                    nhorz::Int            = 8,
                    controls              = (:ln_gdp_diff, :pi_p, :Interest),
                    alpha::Float64        = 0.32,
                    bias_corr::Bool       = true,
                    bootstrap::Bool       = false,
                    boot_num::Int         = 500,
                    boot_blocklength::Int = 4)

    if bias_corr && !bootstrap
        @warn "bias_corr=true but bootstrap=false: delta-method CI is indicative only."
    end

    r = lp(df_merged, outcome;
           xnames  = (:markup_shock,),
           wnames  = controls,
           nlag    = nlag,
           nhorz   = nhorz,
           minhorz = 1)
    f        = irf(r, outcome, :markup_shock)
    coef_raw = collect(f.B)
    ses      = collect(f.SE)
    horizons = collect(f.minhorz : f.minhorz + length(f.B) - 1)

    # Bias correction per outcome
    delta = zeros(length(coef_raw))
    if bias_corr
        ctrl_cols = [controls...]
        w_mat = Matrix{Float64}(undef, nrow(df_merged), length(ctrl_cols))
        for (j, col) in enumerate(ctrl_cols)
            w_mat[:, j] = df_merged[!, col]
        end
        delta = lp_biascorr_delta(coef_raw, w_mat)
    end
    coef = coef_raw .- delta

    # Delta-method CI
    cv   = quantile(Normal(), 1 - alpha / 2)
    lo68 = coef .- cv .* ses
    hi68 = coef .+ cv .* ses

    # Bootstrap
    cis_boot = nothing
    ses_boot = nothing
    if bootstrap
        T        = nrow(df_merged)
        boot_mat = zeros(boot_num, length(coef_raw))
        Threads.@threads for b in 1:boot_num
            if boot_blocklength == 1
                idx_b = rand(1:T, T)
            else
                nblk   = ceil(Int, T / boot_blocklength)
                starts = rand(1:(T - boot_blocklength + 1), nblk)
                idx_b  = vcat([s:s+boot_blocklength-1 for s in starts]...)[1:T]
            end
            df_b = df_merged[idx_b, :]
            try
                r_b = lp(df_b, outcome;
                          xnames  = (:markup_shock,),
                          wnames  = controls,
                          nlag    = nlag,
                          nhorz   = nhorz,
                          minhorz = 1)
                f_b     = irf(r_b, outcome, :markup_shock)
                coef_b  = collect(f_b.B)
                if bias_corr
                    w_b = Matrix{Float64}(undef, nrow(df_b), length(controls))
                    for (j, col) in enumerate(controls)
                        w_b[:, j] = df_b[!, col]
                    end
                    delta_b = lp_biascorr_delta(coef_b, w_b)
                    coef_b  -= delta_b
                end
                boot_mat[b, :] = coef_b
            catch
                boot_mat[b, :] = coef
            end
        end
        ses_boot     = vec(std(boot_mat, dims=1))
        pseudo_truth = vec(mean(boot_mat, dims=1))
        cis_boot     = boot_ci(pseudo_truth, coef, ses_boot,
                               boot_mat, repeat(ses_boot', boot_num, 1), alpha)
    end

    return (coef     = coef,
            coef_raw = coef_raw,
            delta    = delta,
            ses      = ses,
            ses_boot = ses_boot,
            lo68     = lo68,
            hi68     = hi68,
            horizons = horizons,
            cis_boot = cis_boot)
end

# Plot IRF helper
function plot_irf(res, title_str::String, ylabel_str::String,
                  alpha::Float64, bootstrap::Bool)
    if bootstrap && !isnothing(res.cis_boot)
        lo_plot  = res.cis_boot[1, :, 3]
        hi_plot  = res.cis_boot[2, :, 3]
        ci_label = "$(round(Int,(1-alpha)*100))% bootstrap CI (Hall-t)"
    else
        lo_plot  = res.lo68
        hi_plot  = res.hi68
        ci_label = "$(round(Int,(1-alpha)*100))% delta-method CI"
    end

    pl = plot(res.horizons, lo_plot;
              fillrange = hi_plot,
              fillalpha = 0.25,
              fillcolor = :steelblue,
              linealpha = 0,
              label     = ci_label,
              title     = title_str,
              xlabel    = "Horizon (quarters)",
              ylabel    = ylabel_str,
              legend    = :topright)
    plot!(pl, res.horizons, res.coef;
          lw=2, color=:steelblue, label="point estimate (bc)")
    plot!(pl, res.horizons, res.coef_raw;
          lw=1, ls=:dash, color=:gray, label="point estimate (raw)", alpha=0.6)
    hline!([0]; color=:black, ls=:dot, label="")
    return pl
end

# Main analysis
function run_labor_lp_analysis(markup_shocks::AbstractVector,
                                shock_dates::AbstractVector{Date},
                                macro_df::DataFrame;
                                output_suffix::String  = "v1",
                                nhorz::Int             = 8,
                                nlag::Int              = 2,
                                ma_window::Int         = 4,
                                controls               = (:ln_gdp_diff, :pi_p, :Interest),
                                alpha::Float64         = 0.32,
                                bias_corr::Bool        = true,
                                bootstrap::Bool        = false,
                                boot_num::Int          = 500,
                                boot_blocklength::Int  = 4)

    download_gdrive_large(GDRIVE_FILE_ID, LOCAL_PATH; expected_hash=EXPECTED_SHA256)

    df_raw = parse_cps(LOCAL_PATH)
    df_agg = aggregate_quarterly(df_raw)

    shock_df = DataFrame(
        quarter      = year.(shock_dates) .+ (quarterofyear.(shock_dates) .- 1) ./ 4,
        markup_shock = markup_shocks
    )
    unique!(shock_df, :quarter)

    macro_q = select(macro_df, :observation_date, controls...)
    macro_q[!, :quarter] = year.(macro_q.observation_date) .+ (quarterofyear.(macro_q.observation_date) .- 1) ./ 4
    select!(macro_q, Not(:observation_date))
    unique!(macro_q, :quarter)
    dropmissing!(macro_q)

    save_dir = joinpath(pwd(), "results", output_suffix, "labor_lp")
    mkpath(save_dir)

    outcomes = [:mean_income, :mean_hours, :mean_wage, :unemp_rate]
    outcome_labels = Dict(
        :mean_income => "log Income (EARNWEEK)",
        :mean_hours  => "log Hours (UHRSWORKT)",
        :mean_wage   => "log Wage (EARNWEEK/UHRSWORKT)",
        :unemp_rate  => "Unemployment Rate",
    )

    all_results = Dict()

    for occ_id in 1:9
        occ_label = OCC_LABELS[occ_id]
        df_occ = filter(r -> r.occ_cat == occ_id, df_agg)
        select!(df_occ, :quarter, :occ_cat, outcomes...)
        unique!(df_occ, :quarter)
        sort!(df_occ, :quarter)

        for col in outcomes
            df_occ[!, Symbol(string(col)*"_ma")] = moving_avg(df_occ[!, col], ma_window)
        end
        dropmissing!(df_occ)

        df_merged = innerjoin(df_occ, shock_df; on=:quarter)
        df_merged = innerjoin(df_merged, macro_q;  on=:quarter)
        dropmissing!(df_merged)
        sort!(df_merged, :quarter)

        if nrow(df_merged) < nlag + nhorz + 5
            @warn "Occupation $occ_id: too few observations ($(nrow(df_merged))), skipping."
            continue
        end

        occ_plots = []

        for outcome in outcomes
            ma_col = Symbol(string(outcome)*"_ma")
            res = run_lp_irf(df_merged, ma_col;
                             nlag             = nlag,
                             nhorz            = nhorz,
                             controls         = controls,
                             alpha            = alpha,
                             bias_corr        = bias_corr,
                             bootstrap        = bootstrap,
                             boot_num         = boot_num,
                             boot_blocklength = boot_blocklength)
            all_results[(occ_id, outcome)] = res

            pl = plot_irf(res,
                          "Occ $occ_id: $(outcome_labels[outcome])",
                          "Response",
                          alpha, bootstrap)
            push!(occ_plots, pl)
        end

        # Save plots
        occ_dir = joinpath(save_dir, "occ$(occ_id)_$(occ_label)")
        mkpath(occ_dir)
        ncols = 2
        for (i, pl) in enumerate(occ_plots)
            savefig(pl, joinpath(occ_dir, "irf_$(outcomes[i]).png"))
        end
    end

    return all_results, df_agg, shock_df
end