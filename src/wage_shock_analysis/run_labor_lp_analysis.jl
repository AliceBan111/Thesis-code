"""
run_labor_lp_analysis.jl

LP analysis: markup shock → income (EARNWEEK), hours (UHRSWORKT),
             wage (EARNWEEK / UHRSWORKT), unemployment rate
By 9 OCC1990 occupation categories.
Macro controls: ln_gdp_diff, pi_p, Interest

Data layout (cps_00010.dat):
    YEAR       1-4
    MONTH      5-6
    WTFINL     7-20   (implied decimals: 4)  ← person weight for unemployment
    EMPSTAT    21-22
    OCC1990    23-25
    UHRSWORKT  26-28
    EARNWT     29-38  (implied decimals: 4)  ← earnings weight for wages/hours
    EARNWEEK   39-46  (implied decimals: 2)
"""

using Downloads, HTTP, SHA, Printf
using DataFrames, Dates, Statistics, Plots
using LocalProjections

# ─────────────────────────────────────────────────────────────────────────────
# GOOGLE DRIVE CREDENTIALS
# ─────────────────────────────────────────────────────────────────────────────
const GDRIVE_FILE_ID  = "1MF7Am-3utva0NDH0-FE_cvAQOiPN0372"
const EXPECTED_SHA256 = "6dd2f041036209a7abfee1fe29430c97e685478d073178ab25f5b57951b5cb34"
const LOCAL_PATH      = "cps_00010.dat"

# ─────────────────────────────────────────────────────────────────────────────
# 1. GOOGLE DRIVE DOWNLOAD
# ─────────────────────────────────────────────────────────────────────────────
function download_gdrive_large(file_id::String, dest::String;
                                expected_hash::Union{String,Nothing}=nothing)
    if isfile(dest)
        first_bytes = String(read(dest, min(20, filesize(dest))))
        if startswith(first_bytes, "<!") || startswith(first_bytes, "<h")
            println("Existing file is HTML (bad download) — deleting and retrying.")
            rm(dest)
        else
            println("File already exists: $dest — skipping download.")
        end
    end

    if !isfile(dest)
        base_url = "https://drive.google.com/uc?export=download&id=$(file_id)"
        println("Fetching confirmation page...")
        resp = HTTP.get(base_url)
        body = String(resp.body)
        uuid_m = match(r"name=\"uuid\" value=\"([^\"]+)\"", body)
        if !isnothing(uuid_m)
            uuid     = uuid_m.captures[1]
            real_url = "https://drive.usercontent.google.com/download" *
                       "?id=$(file_id)&export=download&confirm=t&uuid=$(uuid)"
        else
            real_url = "https://drive.usercontent.google.com/download" *
                       "?id=$(file_id)&export=download&confirm=t"
        end
        println("Downloading (please wait)...")
        Downloads.download(real_url, dest)
        println("Download complete.")
    end

    if !isnothing(expected_hash)
        println("Verifying hash...")
        actual = bytes2hex(open(sha256, dest))
        if actual != expected_hash
            error("Hash mismatch!\n  expected: $expected_hash\n  actual:   $actual")
        end
        println("✓ Hash verified.")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 2. OCCUPATION CLASSIFICATION  (OCC1990 → 9 categories)
# ─────────────────────────────────────────────────────────────────────────────
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

# ─────────────────────────────────────────────────────────────────────────────
# 3. PARSE FIXED-WIDTH MICRODATA
# ─────────────────────────────────────────────────────────────────────────────
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

            any(isempty, (year_s, month_s, wtfinl_s,
                          empstat_s, occ_s, earnwt_s)) && continue

            empstat = tryparse(Int, empstat_s); isnothing(empstat) && continue
            empstat in (10, 12, 20, 21, 22)   || continue
            yr  = parse(Int, year_s)
            mth = parse(Int, month_s)
            yr < 1994                          && continue
            yr > 2023                          && continue
            (yr == 2023 && mth > 3)            && continue

            occ    = tryparse(Int,     occ_s);   isnothing(occ)    && continue
            wtfinl = tryparse(Float64, wtfinl_s); isnothing(wtfinl) && continue
            earnwt = tryparse(Float64, earnwt_s); isnothing(earnwt) && continue

            wtfinl /= 10_000.0
            earnwt /= 10_000.0

            wtfinl > 0 || continue

            occ_cat = classify_occ1990(occ)
            ismissing(occ_cat) && continue

            is_employed = empstat in (10, 12)

            uhrs     = missing
            earnweek = missing

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
                year        = parse(Int, year_s),
                month       = parse(Int, month_s),
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

# ─────────────────────────────────────────────────────────────────────────────
# 4. AGGREGATE TO QUARTERLY PANEL
# ─────────────────────────────────────────────────────────────────────────────
function aggregate_quarterly(df::DataFrame)::DataFrame
    df[!, :quarter] = df.year .+ ((df.month .- 1) .÷ 3) ./ 4

    agg = combine(
        groupby(df, [:quarter, :occ_cat]),
        [:earnweek, :uhrs, :wtfinl, :earnwt, :is_employed] => (
            (ew, uh, wf, wt, emp) -> begin
                W_person   = sum(wf)
                emp_wf     = sum(Float64.(emp) .* wf)
                unemp_rate = W_person > 0 ? 1.0 - emp_wf / W_person : missing

                emp_mask = emp .== true
                has_earn = emp_mask .& .!ismissing.(ew) .& .!ismissing.(uh)
                W_earn = sum(wt[has_earn])

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

# ─────────────────────────────────────────────────────────────────────────────
# 6. LP HELPER
# ─────────────────────────────────────────────────────────────────────────────
function run_lp_irf(df_merged::DataFrame, outcome::Symbol;
                    nlag::Int  = 2,
                    nhorz::Int = 8,
                    controls   = (:ln_gdp_diff, :pi_p, :Interest))
    r = lp(df_merged, outcome;
           xnames  = (:markup_shock,),
           wnames  = controls,
           nlag    = nlag,
           nhorz   = nhorz,
           minhorz = 1)
    f = irf(r, outcome, :markup_shock)

    coef     = collect(f.B)
    se       = collect(f.SE)
    lo68     = coef .- 1.0 .* se
    hi68     = coef .+ 1.0 .* se
    horizons = collect(f.minhorz : f.minhorz + length(f.B) - 1)
    return coef, lo68, hi68, horizons
end

# ─────────────────────────────────────────────────────────────────────────────
# 7. SANITY CHECK — print unemployment rate stats after aggregation
# ─────────────────────────────────────────────────────────────────────────────
function check_unemp_rates(df_agg::DataFrame)
    println("\n── Unemployment rate sanity check ──")
    println(@sprintf("%-5s  %-30s  %8s  %8s  %8s  %8s",
                     "Occ", "Label", "Mean", "Std", "Min", "Max"))
    for occ_id in 1:9
        df_occ = filter(r -> r.occ_cat == occ_id, df_agg)
        dropmissing!(df_occ, :unemp_rate)
        @printf("%-5d  %-30s  %8.4f  %8.4f  %8.4f  %8.4f\n",
                occ_id, OCC_LABELS[occ_id],
                mean(df_occ.unemp_rate), std(df_occ.unemp_rate),
                minimum(df_occ.unemp_rate), maximum(df_occ.unemp_rate))
    end
    println()
end

# ─────────────────────────────────────────────────────────────────────────────
# 8. MAIN ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────
function run_labor_lp_analysis(markup_shocks::AbstractVector,
                                shock_dates::AbstractVector{Date},
                                macro_df::DataFrame;
                                output_suffix::String = "v1",
                                nhorz::Int            = 8,
                                nlag::Int             = 2,
                                controls              = (:ln_gdp_diff, :pi_p,
                                                         :Interest))

    # ── Download & parse CPS ─────────────────────────────────────────────────
    download_gdrive_large(GDRIVE_FILE_ID, LOCAL_PATH; expected_hash=EXPECTED_SHA256)

    println("Parsing microdata...")
    df_raw = parse_cps(LOCAL_PATH)
    println("  $(nrow(df_raw)) observations loaded (employed + unemployed).")

    df_agg = aggregate_quarterly(df_raw)
    check_unemp_rates(df_agg)

    # ── Shock series → quarter index ─────────────────────────────────────────
    shock_df = DataFrame(
        quarter      = year.(shock_dates) .+ (quarterofyear.(shock_dates) .- 1) ./ 4,
        markup_shock = markup_shocks
    )
    unique!(shock_df, :quarter)

    # ── Macro controls → quarter index ───────────────────────────────────────
    macro_q = select(macro_df, :observation_date, controls...)
    macro_q[!, :quarter] = year.(macro_q.observation_date) .+ 
                           (quarterofyear.(macro_q.observation_date) .- 1) ./ 4
    select!(macro_q, Not(:observation_date))
    unique!(macro_q, :quarter)
    dropmissing!(macro_q)

    # ── Output directory ─────────────────────────────────────────────────────
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
        println("\n── Occupation $occ_id: $occ_label ──")

        df_occ = filter(r -> r.occ_cat == occ_id, df_agg)
        unique!(df_occ, :quarter)
        sort!(df_occ, :quarter)

        # ── Three-way merge: CPS + markup shock + macro controls
        df_merged = innerjoin(df_occ, shock_df; on=:quarter)
        df_merged = innerjoin(df_merged, macro_q; on=:quarter)
        dropmissing!(df_merged)
        sort!(df_merged, :quarter)

        if nrow(df_merged) < nlag + nhorz + 5
            @warn "Occupation $occ_id: too few obs ($(nrow(df_merged))), skipping."
            continue
        end
        println("  $(nrow(df_merged)) quarters after merge.")

        # ── Run LP for each outcome ───────────────────────────────────────────
        occ_plots = []
        for outcome in outcomes
            coef, lo, hi, horz = run_lp_irf(df_merged, outcome;
                                             nlag=nlag, nhorz=nhorz,
                                             controls=controls)
            all_results[(occ_id, outcome)] = (coef=coef, lo=lo, hi=hi, horizons=horz)

            pl = plot(horz, coef;
                      ribbon    = (coef .- lo, hi .- coef),
                      fillalpha = 0.25, lw=2,
                      label     = outcome_labels[outcome],
                      title     = "Occ $occ_id: $(outcome_labels[outcome])",
                      xlabel    = "Horizon (quarters)",
                      ylabel    = "Response",
                      legend    = :topright)
            hline!([0]; color=:black, ls=:dot, label="")
            push!(occ_plots, pl)
        end

        # ── Consistency check: IRF_wage ≈ IRF_income - IRF_hours ─────────────
        wage_keys = (:mean_income, :mean_hours, :mean_wage)
        if all(k -> haskey(all_results, (occ_id, k)), wage_keys)
            ic           = all_results[(occ_id, :mean_income)].coef
            hc           = all_results[(occ_id, :mean_hours)].coef
            wc           = all_results[(occ_id, :mean_wage)].coef
            implied_wage = ic .- hc

            pl_check = plot(all_results[(occ_id, :mean_wage)].horizons, wc;
                            lw=2, label="IRF_wage (direct)", color=:blue,
                            title="Occ $occ_id: consistency check",
                            xlabel="Horizon (quarters)", ylabel="Response")
            plot!(pl_check, all_results[(occ_id, :mean_income)].horizons, implied_wage;
                  lw=2, ls=:dash, label="IRF_income − IRF_hours", color=:red)
            hline!([0]; color=:black, ls=:dot, label="")
            push!(occ_plots, pl_check)

            max_dev = maximum(abs.(wc .- implied_wage))
            println("  Consistency check max deviation: $(round(max_dev, digits=6))")
        end

        # ── Save ──────────────────────────────────────────────────────────────
        occ_dir = joinpath(save_dir, "occ$(occ_id)_$(occ_label)")
        mkpath(occ_dir)
        fig = plot(occ_plots...; layout=(2, 3), size=(1400, 900),
                   plot_title="Markup Shock → Labor Outcomes, Occ $occ_id")
        savefig(fig, joinpath(occ_dir, "irf_panel.png"))
        println("  Saved: $(joinpath(occ_dir, "irf_panel.png"))")
        display(fig)
    end

    println("\n✓ All done. Results in: $save_dir")
    return all_results, df_agg, shock_df, macro_q
end