"""
    run_wage_heterogeneity_analysis(markup_shocks, shock_dates; output_suffix="method2")

Performs wage heterogeneity analysis by:
1. Downloading and processing CPS microdata.
2. Classifying occupations into White-collar and Blue-collar.
3. Calculating correlation and dynamic responses between markup shocks and wage growth.
4. Generating time-series, scatter, and lag-response plots.
"""
function run_wage_heterogeneity_analysis(markup_shocks, shock_dates; output_suffix="method2")
    # 1. Date Alignment
    shock_df = DataFrame(date = shock_dates, markup_shock = markup_shocks)
    shock_df[!, :quarter] = Dates.year.(shock_df.date) .+ (Dates.quarterofyear.(shock_df.date) .- 1) ./ 4

    # 2. CPS Data Procurement
    url = "https://www.dropbox.com/scl/fi/6m7ccu2c2n58f6zcrjd7k/cps_00006.dat?rlkey=n937jeomgs4be3k9i70oe2zde&st=v1dw2ea7&dl=1"
    local_path = "cps_00006.dat"
    expected_hash = "34fcc6371a47152ae5fc48f83421cf28336acaaa5c2b4185914eef8fef131539"

    if !isfile(local_path)
        Downloads.download(url, local_path)
    end

    # Integrity Check
    actual_hash = bytes2hex(open(sha256, local_path))
    if actual_hash != expected_hash
        error("Integrity error: CPS data hash mismatch!")
    end

    # 3. Processing Microdata
    df_working = DataFrame(YEAR=Int[], MONTH=Int[], EARNWT=Float64[], EMPSTAT=Int[], OCC=Int[], EARNWEEK=Float64[])

    open(local_path, "r") do io
        for line in eachline(io)
            empstat = parse(Int, strip(line[7:8]))
            if empstat == 10 || empstat == 12
                occ     = parse(Int, strip(line[9:11]))
                earn    = parse(Float64, strip(line[22:29])) / 100.0
                earnwt  = parse(Float64, strip(line[12:21])) / 10000.0
                
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

    df_working.EARNWEEK .= df_working.EARNWEEK ./ 100
    df_working.EARNWT .= df_working.EARNWT ./ 10000

    # Occupation Classification Logic
    function classify_collar(occ::Int)
        if occ in vcat(473:498, 503:549, 558:599, 614:617, 628:699, 703:799, 803:889)
            return 0 # Blue-collar
        elseif occ in vcat(3:37, 43:200, 203:235, 243:283, 303:389, 405:469)
            return 1 # White-collar
        else
            return missing
        end
    end

    df_working[!, :collar] = [classify_collar(occ) for occ in df_working.OCC]
    df_working = dropmissing(df_working, :collar)
    df_working[!, :quarter] = df_working.YEAR .+ ((df_working.MONTH .- 1) .÷ 3) ./ 4

    # 4. Aggregate Wages & Calculate Growth
    df_q = combine(groupby(df_working, [:quarter, :collar]), 
           [:EARNWEEK, :EARNWT] => ((w, wt) -> sum(w .* wt) / sum(wt)) => :avg_wage)

    # Split and calculate log growth
    df_white = filter(row -> row.collar == 1, df_q)
    df_blue  = filter(row -> row.collar == 0, df_q)
    
    df_white[!, :dlog_wage] = [missing; diff(log.(df_white.avg_wage))]
    df_blue[!, :dlog_wage]  = [missing; diff(log.(df_blue.avg_wage))]

    # 5. Merge with Shocks
    unique!(df_white, :quarter)
    unique!(df_blue, :quarter)
    unique!(shock_df, :quarter)

    df_plot = innerjoin(
        innerjoin(df_white[:, [:quarter, :dlog_wage]], df_blue[:, [:quarter, :dlog_wage]], on=:quarter, makeunique=true),
        shock_df[:, [:quarter, :markup_shock]], on=:quarter
    )
    rename!(df_plot, :dlog_wage => :white_wage, :dlog_wage_1 => :blue_wage)
    df_plot = dropmissing(df_plot)

    # 6. Save path
    save_path = joinpath(pwd(), "results", output_suffix, "wage_heterogeneity")
    mkpath(save_path)
    
    # 7. Visualization & Correlation
    # Plot 1: Time Series
    plt1 = plot(df_plot.quarter, [df_plot.white_wage df_plot.blue_wage], 
                label=["White-collar" "Blue-collar"], color=[:blue :green], lw=2, size=(1000,600),
                title="Wage Growth vs Markup Shocks ($output_suffix)")
    plt1_twin = twinx(plt1)
    plot!(plt1_twin, df_plot.quarter, df_plot.markup_shock, label="Markup Shock", color=:red, ls=:dash, alpha=0.6)
    
    savefig(plt1, joinpath(save_path, "wage_vs_shocks_timeseries.png"))

    # Plot 2: Dynamic Correlation (Pseudo-IRF)
    max_lag = 4
    white_irf = [cor(df_plot.markup_shock[1:end-l], df_plot.white_wage[l+1:end]) for l in 0:max_lag]
    blue_irf  = [cor(df_plot.markup_shock[1:end-l], df_plot.blue_wage[l+1:end]) for l in 0:max_lag]
    
    plt2 = plot(0:max_lag, [white_irf blue_irf], marker=[:circle :diamond], 
                label=["White-collar" "Blue-collar"], title="Wage Correlation Response ($output_suffix)",
                xlabel="Lags (Quarters)", ylabel="Correlation")
    hline!([0], color=:black, ls=:dot, label="")
    savefig(plt2, joinpath(save_path, "wage_correlation_irf.png"))

    println("Plots successfully saved to: $save_path")
    
    display(plt1)
    display(plt2)

    return df_plot, (white_irf, blue_irf)
end