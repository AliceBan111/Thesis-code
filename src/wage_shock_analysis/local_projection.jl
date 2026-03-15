function run_wage_lp_analysis(markup_shocks, shock_dates; output_suffix="method2", nhorz=8, nlag=4)

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
    df_working = DataFrame(YEAR=Int[], MONTH=Int[], EARNWT=Float64[], EMPSTAT=Int[], OCC=Int[], EARNWEEK=Float64[])

    open(local_path, "r") do io
        for line in eachline(io)
            length(line) < 29 && continue 
            
            empstat_str = strip(line[7:8])
            isempty(empstat_str) && continue
            empstat = tryparse(Int, empstat_str)
            isnothing(empstat) && continue
            
            if empstat == 10 || empstat == 12
                occ_str    = strip(line[9:11])
                earn_str   = strip(line[22:29])
                earnwt_str = strip(line[12:21])

                (isempty(occ_str) || isempty(earn_str) || isempty(earnwt_str)) && continue
                
                occ    = tryparse(Int,     occ_str);    isnothing(occ)    && continue
                earn   = tryparse(Float64, earn_str);   isnothing(earn)   && continue
                earnwt = tryparse(Float64, earnwt_str); isnothing(earnwt) && continue

                earn   /= 100.0
                earnwt /= 10000.0

                if occ != 0 && earn != 9999.99 && earnwt > 0
                    push!(df_working, (
                        parse(Int, strip(line[1:4])),
                        parse(Int, strip(line[5:6])),
                        earnwt, empstat, occ, earn
                    ))
                end
            end
        end
    end

    # 4. Occupation Classification
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

    # 5. Aggregate Wages & Log Growth
    df_q = combine(
        groupby(df_working, [:quarter, :collar]),
        [:EARNWEEK, :EARNWT] => ((w, wt) -> sum(w .* wt) / sum(wt)) => :avg_wage
    )

    df_white = sort(filter(row -> row.collar == 1, df_q), :quarter)
    df_blue  = sort(filter(row -> row.collar == 0, df_q), :quarter)

    df_white[!, :dlog_wage] = [missing; diff(log.(df_white.avg_wage))]
    df_blue[!, :dlog_wage]  = [missing; diff(log.(df_blue.avg_wage))]

    # 6. Merge into one panel
    unique!(df_white, :quarter)
    unique!(df_blue,  :quarter)
    unique!(shock_df, :quarter)

    df_merged = innerjoin(
        innerjoin(
            df_white[:, [:quarter, :dlog_wage]],
            df_blue[:,  [:quarter, :dlog_wage]],
            on=:quarter, makeunique=true
        ),
        shock_df[:, [:quarter, :markup_shock]],
        on=:quarter
    )
    rename!(df_merged, :dlog_wage => :white_wage, :dlog_wage_1 => :blue_wage)
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

    df_merged[!, :white_wage_ma] = moving_avg(df_merged.white_wage, 4)
    df_merged[!, :blue_wage_ma]  = moving_avg(df_merged.blue_wage,  4)
    df_merged_clean = dropmissing(df_merged, [:white_wage_ma, :blue_wage_ma])

    # 8. Local Projections
    r_white = lp(df_merged_clean, :white_wage_ma,
             xnames  = (:markup_shock,),
             wnames  = (:markup_shock,),
             nlag    = 2,
             nhorz   = nhorz,
             minhorz = 1)

    r_blue = lp(df_merged_clean, :blue_wage_ma,
                xnames  = (:markup_shock,),
                wnames  = (:markup_shock,),
                nlag    = 2,
                nhorz   = nhorz,
                minhorz = 1)

    irf_white = irf(r_white, :white_wage_ma, :markup_shock)
    irf_blue  = irf(r_blue,  :blue_wage_ma,  :markup_shock)

    # 9. Extract IRF
    function extract_irf(f)
        coef     = collect(f.B)
        se       = collect(f.SE)
        lo       = coef .- 1.645 .* se
        hi       = coef .+ 1.645 .* se
        horizons = collect(f.minhorz : (f.minhorz + length(f.B) - 1))
        return coef, lo, hi, horizons
    end

    white_coef, white_lo, white_hi, horizons_w = extract_irf(irf_white)
    blue_coef,  blue_lo,  blue_hi,  horizons_b = extract_irf(irf_blue)

    # 10. Save & Plot
    save_path = joinpath(pwd(), "results", output_suffix, "wage_heterogeneity_lp")
    mkpath(save_path)

    plt_white = plot(horizons_w, white_coef,
                     ribbon = (white_coef .- white_lo, white_hi .- white_coef),
                     fillalpha = 0.25, color = :blue, lw = 2,
                     label = "White-collar",
                     title = "LP-IRF: Markup Shock → White-collar Wage ($output_suffix)",
                     xlabel = "Horizon (Quarters)", ylabel = "Response of Δlog Wage",
                     size = (900, 500))
    hline!([0], color = :black, ls = :dot, label = "")
    savefig(plt_white, joinpath(save_path, "lp_irf_white_collar.png"))

    plt_blue = plot(horizons_b, blue_coef,
                    ribbon = (blue_coef .- blue_lo, blue_hi .- blue_coef),
                    fillalpha = 0.25, color = :green, lw = 2,
                    label = "Blue-collar",
                    title = "LP-IRF: Markup Shock → Blue-collar Wage ($output_suffix)",
                    xlabel = "Horizon (Quarters)", ylabel = "Response of Δlog Wage",
                    size = (900, 500))
    hline!([0], color = :black, ls = :dot, label = "")
    savefig(plt_blue, joinpath(save_path, "lp_irf_blue_collar.png"))

    plt_compare = plot(horizons_w, white_coef,
                       ribbon = (white_coef .- white_lo, white_hi .- white_coef),
                       fillalpha = 0.2, color = :blue, lw = 2, label = "White-collar",
                       title = "LP-IRF Comparison: Markup Shock → Wages ($output_suffix)",
                       xlabel = "Horizon (Quarters)", ylabel = "Response of Δlog Wage",
                       size = (1000, 600))
    plot!(plt_compare, horizons_b, blue_coef,
          ribbon = (blue_coef .- blue_lo, blue_hi .- blue_coef),
          fillalpha = 0.2, color = :green, lw = 2, label = "Blue-collar")
    hline!([0], color = :black, ls = :dot, label = "")
    savefig(plt_compare, joinpath(save_path, "lp_irf_comparison.png"))

    println("LP plots saved to: $save_path")
    display(plt_white)
    display(plt_blue)
    display(plt_compare)

    return df_merged_clean, r_white, r_blue, irf_white, irf_blue
end