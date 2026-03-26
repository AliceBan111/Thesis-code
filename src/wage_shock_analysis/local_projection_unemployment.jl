function run_unemployment_lp_analysis(markup_shocks, shock_dates; output_suffix="method2", nhorz=8, nlag=2)

    # 1. Date Alignment
    shock_df = DataFrame(date = shock_dates, markup_shock = markup_shocks)
    shock_df[!, :quarter] = Dates.year.(shock_df.date) .+ (Dates.quarterofyear.(shock_df.date) .- 1) ./ 4

    # 2. CPS Data Procurement
    url = "https://www.dropbox.com/scl/fi/iuhysf8d9u3jwyszmydoi/cps_00008.dat?rlkey=s778eox3hhtsiagm8islncwk3&st=a7b18oue&dl=0"
    local_path = "cps_00008.dat"
    expected_hash = "1e78220d8b6f0bed74ee51ee6233378f94a809952530338f21b207faf1f41acf"

    if !isfile(local_path)
        Downloads.download(url, local_path)
    end

    actual_hash = bytes2hex(open(sha256, local_path))
    if actual_hash != expected_hash
        error("Integrity error: CPS data hash mismatch!")
    end

    # 3. Processing Microdata
    # Include both employed (10,12) and unemployed (20,21,22) for labor force
    df_working = DataFrame(YEAR=Int[], MONTH=Int[], EARNWT=Float64[], EMPSTAT=Int[], OCC=Int[])

    open(local_path, "r") do io
        for line in eachline(io)
            length(line) < 29 && continue

            empstat_str = strip(line[7:8])
            isempty(empstat_str) && continue
            empstat = tryparse(Int, empstat_str)
            isnothing(empstat) && continue

            # Include employed (10,12) and unemployed (20,21,22)
            if empstat in [10, 12, 20, 21, 22]
                occ_str    = strip(line[9:11])
                earnwt_str = strip(line[12:21])

                (isempty(occ_str) || isempty(earnwt_str)) && continue

                occ    = tryparse(Int,     occ_str);    isnothing(occ)    && continue
                earnwt = tryparse(Float64, earnwt_str); isnothing(earnwt) && continue

                earnwt /= 10000.0

                earnwt > 0 || continue

                push!(df_working, (
                    parse(Int, strip(line[1:4])),
                    parse(Int, strip(line[5:6])),
                    earnwt, empstat, occ
                ))
            end
        end
    end

    # 4. Occupation Classification
    # Note: unemployed persons may have OCC=0, so we classify based on last known occupation
    # Only persons with a valid OCC code can be assigned to a collar group
    function classify_collar(occ::Int)
        if occ in vcat(473:498, 503:549, 558:599, 614:617, 628:699, 703:799, 803:889)
            return 0   # Blue-collar
        elseif occ in vcat(3:37, 43:200, 203:235, 243:283, 303:389, 405:469)
            return 1   # White-collar
        else
            return missing
        end
    end

    df_working[!, :collar] = [classify_collar(occ) for occ in df_working.OCC]
    df_working = dropmissing(df_working, :collar)
    df_working[!, :quarter] = df_working.YEAR .+ ((df_working.MONTH .- 1) .÷ 3) ./ 4
    df_working[!, :employed] = [empstat in [10, 12] ? 1.0 : 0.0 for empstat in df_working.EMPSTAT]

    # 5. Aggregate to quarterly unemployment rate by collar
    # u_h,t = (labor_force - employed) / labor_force  (weighted)
    df_q = combine(
        groupby(df_working, [:quarter, :collar]),
        [:employed, :EARNWT] => (
            (emp, wt) -> begin
                labor_force = sum(wt)
                employed    = sum(emp .* wt)
                unemp_rate  = (labor_force - employed) / labor_force
                return unemp_rate
            end
        ) => :unemp_rate
    )

    df_white_u = sort(filter(row -> row.collar == 1, df_q), :quarter)
    df_blue_u  = sort(filter(row -> row.collar == 0, df_q), :quarter)

    # 6. Merge with shock data
    unique!(df_white_u, :quarter)
    unique!(df_blue_u,  :quarter)
    unique!(shock_df,   :quarter)

    df_merged = innerjoin(
        innerjoin(
            df_white_u[:, [:quarter, :unemp_rate]],
            df_blue_u[:,  [:quarter, :unemp_rate]],
            on=:quarter, makeunique=true
        ),
        shock_df[:, [:quarter, :markup_shock]],
        on=:quarter
    )
    rename!(df_merged, :unemp_rate => :white_unemp, :unemp_rate_1 => :blue_unemp)
    df_merged = dropmissing(df_merged)
    sort!(df_merged, :quarter)

    # 7. Moving average smoothing
    function moving_avg(x, n)
        result = Vector{Union{Missing, Float64}}(missing, length(x))
        for i in n:length(x)
            result[i] = mean(x[i-n+1:i])
        end
        return result
    end

    df_merged[!, :white_unemp_ma] = moving_avg(df_merged.white_unemp, 4)
    df_merged[!, :blue_unemp_ma]  = moving_avg(df_merged.blue_unemp,  4)
    df_merged_clean = dropmissing(df_merged, [:white_unemp_ma, :blue_unemp_ma])

    # 8. Local Projections
    r_white_u = lp(df_merged_clean, :white_unemp_ma,
                   xnames  = (:markup_shock,),
                   wnames  = (:markup_shock,),
                   nlag    = nlag,
                   nhorz   = nhorz,
                   minhorz = 1)

    r_blue_u = lp(df_merged_clean, :blue_unemp_ma,
                  xnames  = (:markup_shock,),
                  wnames  = (:markup_shock,),
                  nlag    = nlag,
                  nhorz   = nhorz,
                  minhorz = 1)

    irf_white_u = irf(r_white_u, :white_unemp_ma, :markup_shock)
    irf_blue_u  = irf(r_blue_u,  :blue_unemp_ma,  :markup_shock)

    # 9. Extract IRF
    function extract_irf(f)
        coef     = collect(f.B)
        se       = collect(f.SE)
        lo       = coef .- 1.645 .* se
        hi       = coef .+ 1.645 .* se
        horizons = collect(f.minhorz : (f.minhorz + length(f.B) - 1))
        return coef, lo, hi, horizons
    end

    white_u_coef, white_u_lo, white_u_hi, horizons_wu = extract_irf(irf_white_u)
    blue_u_coef,  blue_u_lo,  blue_u_hi,  horizons_bu = extract_irf(irf_blue_u)

    # 10. Save & Plot
    save_path = joinpath(pwd(), "results", output_suffix, "unemployment_lp")
    mkpath(save_path)

    plt_white_u = plot(horizons_wu, white_u_coef,
                       ribbon     = (white_u_coef .- white_u_lo, white_u_hi .- white_u_coef),
                       fillalpha  = 0.25, color = :blue, lw = 2,
                       label      = "White-collar",
                       title      = "LP-IRF: Markup Shock → White-collar Unemployment ($output_suffix)",
                       xlabel     = "Horizon (Quarters)",
                       ylabel     = "Response of Unemployment Rate",
                       size       = (900, 500))
    hline!([0], color = :black, ls = :dot, label = "")
    savefig(plt_white_u, joinpath(save_path, "lp_irf_white_unemp.png"))

    plt_blue_u = plot(horizons_bu, blue_u_coef,
                      ribbon     = (blue_u_coef .- blue_u_lo, blue_u_hi .- blue_u_coef),
                      fillalpha  = 0.25, color = :green, lw = 2,
                      label      = "Blue-collar",
                      title      = "LP-IRF: Markup Shock → Blue-collar Unemployment ($output_suffix)",
                      xlabel     = "Horizon (Quarters)",
                      ylabel     = "Response of Unemployment Rate",
                      size       = (900, 500))
    hline!([0], color = :black, ls = :dot, label = "")
    savefig(plt_blue_u, joinpath(save_path, "lp_irf_blue_unemp.png"))

    plt_compare_u = plot(horizons_wu, white_u_coef,
                         ribbon    = (white_u_coef .- white_u_lo, white_u_hi .- white_u_coef),
                         fillalpha = 0.2, color = :blue, lw = 2, label = "White-collar",
                         title     = "LP-IRF Comparison: Markup Shock → Unemployment ($output_suffix)",
                         xlabel    = "Horizon (Quarters)",
                         ylabel    = "Response of Unemployment Rate",
                         size      = (1000, 600))
    plot!(plt_compare_u, horizons_bu, blue_u_coef,
          ribbon    = (blue_u_coef .- blue_u_lo, blue_u_hi .- blue_u_coef),
          fillalpha = 0.2, color = :green, lw = 2, label = "Blue-collar")
    hline!([0], color = :black, ls = :dot, label = "")
    savefig(plt_compare_u, joinpath(save_path, "lp_irf_unemp_comparison.png"))

    println("Unemployment LP plots saved to: $save_path")
    display(plt_white_u); display(plt_blue_u); display(plt_compare_u)

    return df_merged_clean, r_white_u, r_blue_u, irf_white_u, irf_blue_u
end