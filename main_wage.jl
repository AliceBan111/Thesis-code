# Include utility functions
include("src/ToolBox/load_var.jl")
#include("src/data/build_macro_data.jl")
include("src/data/build_macro_data_with_TFP.jl")
#include("src/data/build_markup_data_method2.jl")
include("src/data/build_markup_data_method3.jl")
#include("src/data/build_markup_data_method2_no_hp.jl")
#include("src/data/build_markup_data_method2_no_hp_latest.jl")
#include("src/models/var_sr.jl")
include("src/models/var_sr_v2.jl")
include("src/models/var.jl")
include("src/wage_shock_analysis/extract_shocks.jl")
#include("src/wage_shock_analysis/wage_heterogeneity.jl")
#include("src/wage_shock_analysis/local_projection_wage.jl")
#include("src/wage_shock_analysis/local_projection_unemployment.jl")
include("src/wage_shock_analysis/run_labor_lp_analysis.jl")
#include("src/wage_shock_analysis/run_labor_lp_analysis_v2.jl")
include("src/wage_shock_analysis/wald_test_occ.jl")
include("src/wage_shock_analysis/plot_occ_irf_comparison.jl")
include("src/wage_shock_analysis/dynamic_multiplier.jl")

using CSV, DataFrames, Statistics, Plots, Dates, Downloads, SHA, LocalProjections

cd(@__DIR__)

# 1. Load macro data
# start_date_method2 = Date("1987-01-01")
# end_date_method2   = Date("2012-12-31")

start_date_method2_latest = Date("1997-01-01")
end_date_method2_latest   = Date("2023-12-31")

start_date_method3 = Date("1976-01-01")
end_date_method3   = Date("2024-12-31")

# macro_df_method2 = build_macro_data(start_date_method2, end_date_method2)
# macro_df_method2_no_hp = build_macro_data(start_date_method2, end_date_method2)
macro_df_method2_latest = build_macro_data_with_TFP(start_date_method2_latest, end_date_method2_latest)
macro_df_method3 = build_macro_data_with_TFP(start_date_method3, end_date_method3)


# 2. Build markup data
# df_method2 = build_markup_method2(macro_df_method2)
# df_method2_no_hp = build_markup_method2_no_hp(macro_df_method2_no_hp)
df_method2_no_hp_latest = build_markup_method2_no_hp_latest(macro_df_method2_latest, "data/KLEMS_latest.xlsx")
df_method3 = build_markup_method3(macro_df_method3)

# 2.5 Plot VAR inputs
plot_VAR_inputs(df_method2_no_hp_latest; tag="method2_no_hp_latest", use_growth  = true)
plot_VAR_inputs(df_method2_no_hp_latest; tag="method2_no_hp_latest", use_growth  = false)
plot_VAR_inputs(df_method3; tag="method3", use_growth  = false)


# 3. Estimate VAR + Sign Restrictions
# res_method2 = estimate_VAR_SR(df_method2; method_name="method2")
# res_method2_no_hp_level = estimate_VAR_SR(df_method2_no_hp; method_name="method2_no_hp", use_growth = false)


res_method2_no_hp_latest_growth = estimate_VAR_SR(df_method2_no_hp_latest; method_name="method2_no_hp_latest", use_growth  = true)
res_method2_no_hp_latest_level = estimate_VAR_SR(df_method2_no_hp_latest; method_name="method2_no_hp_latest", use_growth  = false)
res_method3 = estimate_VAR_SR(df_method3; method_name="method3", use_growth  = false)

# resids_growth  = res_method2_no_hp_latest_growth[:resids]
# ARCH_LM_test(resids_growth; lags=4)

# resids_level  = res_method2_no_hp_latest_level[:resids]
# ARCH_LM_test(resids_level; lags=4)

# varnames = ["GDP growth", "Inflation", "Unemployment", "Markup", "Interest rate", "TFP"]
# plot_resid_diagnostics(resids_growth, varnames)


# 4. extract shocks
# markup_shocks_m2_growth = extract_structural_shocks(res_method2_no_hp_latest_growth, shock_index=1)
# n_shocks_growth = length(markup_shocks_m2)

# markup_shocks_m2_level = extract_structural_shocks(res_method2_no_hp_latest_level, shock_index=1)
# n_shocks_level = length(markup_shocks_m2_level)

markup_shocks_m3_level = extract_structural_shocks(res_method3, shock_index=1)
n_shocks_level = length(markup_shocks_m3_level)

# Get the corresponding dates for the shocks
# shock_dates = last(df_method2_no_hp_latest.observation_date, n_shocks)

# println("Shocks length: ", n_shocks)
# println("Dates length: ", length(shock_dates))

# shock_dates_level = last(df_method2_no_hp_latest.observation_date, n_shocks_level)

# println("Shocks length: ", n_shocks_level)
# println("Dates length: ", length(shock_dates_level))

shock_dates_level = last(df_method3.observation_date, n_shocks_level)

println("Shocks length: ", n_shocks_level)
println("Dates length: ", length(shock_dates_level))

plot(
    shock_dates_level, 
    markup_shocks_m3_level,
    seriestype = :line,
    lw = 2,
    color = :blue,
    xlabel = "Date",
    ylabel = "Markup Shock",
    title = "Time Series of Markup Shocks (Method 3 Level)",
    legend = false,
    grid = true
)

# 可选：添加水平线 y=0，方便观察正负冲击
hline!([0]; color=:black, ls=:dash)

# 分组
shock_before1994 = [s for (s,d) in zip(markup_shocks_m3_level, shock_dates_level) if year(d) < 1994]
shock_after1994  = [s for (s,d) in zip(markup_shocks_m3_level, shock_dates_level) if year(d) >= 1994]

# 统计函数
function shock_stats(shocks::Vector{Float64})
    n   = length(shocks)
    mean_s = mean(shocks)
    std_s  = std(shocks)
    min_s  = minimum(shocks)
    max_s  = maximum(shocks)
    median_s = median(shocks)
    return n, mean_s, std_s, min_s, max_s, median_s
end

# 计算
n1, mean1, std1, min1, max1, median1 = shock_stats(shock_before1994)
n2, mean2, std2, min2, max2, median2 = shock_stats(shock_after1994)

# 打印结果
W = 80
println("Markup Shock Statistics")
println("="^W)
@printf("%-15s  %6s  %10s  %10s  %10s  %10s  %10s\n",
        "Period", "N", "Mean", "Std", "Min", "Max", "Median")
println("-"^W)
@printf("%-15s  %6d  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n",
        "<1994", n1, mean1, std1, min1, max1, median1)
@printf("%-15s  %6d  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n",
        ">=1994", n2, mean2, std2, min2, max2, median2)
println("="^W)

# 5. Run heterogeneity analysis
# df_results, irf_data = run_wage_heterogeneity_analysis(
#     markup_shocks_m2,
#     shock_dates;
#     output_suffix="method2"
# )

# df_results, irf_data = run_wage_heterogeneity_analysis(
#     markup_shocks_m2_level,
#     shock_dates_level;
#     output_suffix="method2_no_hp_level"
# )

# 6. Local Projection
# df_lp_growth, r_white_g, r_blue_g, irf_white_g, irf_blue_g = run_wage_lp_analysis(
#     markup_shocks_m2,
#     shock_dates;
#     output_suffix = "method2_latest_growth",
#     nhorz = 8,
#     nlag  = 4
# )

# df_lp_level, r_white_l, r_blue_l, irf_white_l, irf_blue_l = run_wage_lp_analysis(
#     markup_shocks_m2_level,
#     shock_dates_level;
#     output_suffix = "method2_latest_level",
#     nhorz = 8,
#     nlag  = 2
# )

# df_unemp_level, r_white_u, r_blue_u, irf_white_u, irf_blue_u = run_unemployment_lp_analysis(
#     markup_shocks_m2_level,
#     shock_dates_level;
#     output_suffix = "method2_latest_level",
#     nhorz = 8,
#     nlag  = 2
# )

results, df_agg, shock_df, macro_q = run_labor_lp_analysis(
    markup_shocks_m3_level,
    shock_dates_level,
    macro_df_method3;
    output_suffix = "method3",
    nhorz         = 8,
    nlag          = 2
)

results, df_agg, shock_df = run_labor_lp_analysis(
    markup_shocks_m3_level, shock_dates_level, macro_df_method3;
    output_suffix    = "method3_boot",
    nhorz            = 8,
    nlag             = 2,
    bias_corr        = true,
    bootstrap        = true,
    boot_num         = 500,
    boot_blocklength = 4)


wald_results = run_wald_tests(df_agg, shock_df, macro_q;
                               nlag     = 2,   
                               nhorz    = 8,
                               controls = (:ln_gdp_diff, :pi_p, :Interest))


# 7. Dynamic Multiplier
# White-collar
println("=== White-collar ===")
phi_white, cumw_white, cumu_white = compute_dynamic_multiplier(r_white_l, r_white_u; nhorz=8)

# Blue-collar
println("\n=== Blue-collar ===")
phi_blue, cumw_blue, cumu_blue = compute_dynamic_multiplier(r_blue_l, r_blue_u; nhorz=8)

println("\nMultiplier ratio Φ_w(white) / Φ_w(blue) at each horizon:")
for h in 1:8
    println("H=$h  ratio = $(round(phi_white[h] / phi_blue[h], sigdigits=4))")
end


println("=== Unemployment IRF (white, nlag=4) ===")
for h in 1:8
    b = r_white_u.B[1, 1, h]
    println("h=$h  β=$b  cumsum=$(sum([r_white_u.B[1,1,k] for k in 1:h]))")
end

println("=== Unemployment IRF (white, nlag=2) ===")
for h in 1:8
    b = r_white_u.B[1, 1, h]
    println("h=$h  β=$b  cumsum=$(sum([r_white_u.B[1,1,k] for k in 1:h]))")
end

println("\n=== Wage IRF (white, nlag=2) ===")
for h in 1:8
    b = r_white_l.B[1, 1, h]
    println("h=$h  β=$b  cumsum=$(sum([r_white_l.B[1,1,k] for k in 1:h]))")
end

println("\n=== Wage IRF (white, nlag=4) ===")
for h in 1:8
    b = r_white_l.B[1, 1, h]
    println("h=$h  β=$b  cumsum=$(sum([r_white_l.B[1,1,k] for k in 1:h]))")
end


using Plots
plot(shock_dates_level, markup_shocks_m3_level,
     title="Markup shock series",
     xlabel="Date", ylabel="Shock", label="")
hline!([0]; color=:black, ls=:dot, label="")

println("Mean:   ", mean(markup_shocks_m3_level))
println("Std:    ", std(markup_shocks_m3_level))
println("Ratio:  ", abs(mean(markup_shocks_m3_level)) / std(markup_shocks_m3_level))

for occ_id in 1:9
    df_occ = filter(r -> r.occ_cat == occ_id, df_agg)
    df_occ = innerjoin(df_occ, shock_df; on=:quarter)
    df_occ = innerjoin(df_occ, macro_q;  on=:quarter)
    dropmissing!(df_occ)
    println("Occ $occ_id ($(OCC_LABELS[occ_id])): $(nrow(df_occ)) quarters")
end

println(describe(df_agg[:, [:mean_income, :mean_hours, :mean_wage]]))

# 1. 确认时间范围
for occ_id in 1:9
    df_occ = filter(r -> r.occ_cat == occ_id, df_agg)
    dropmissing!(df_occ, :mean_hours)
    println("Occ $occ_id: $(nrow(df_occ)) quarters, $(minimum(df_occ.quarter)) - $(maximum(df_occ.quarter))")
end

# 2. 确认 shock 时间范围
println("Shock dates: $(minimum(shock_dates_level)) - $(maximum(shock_dates_level))")
println("Shock quarters: $(minimum(shock_df.quarter)) - $(maximum(shock_df.quarter))")

# 3. merge之后实际有多少季度
df_test = filter(r -> r.occ_cat == 3, df_agg)
dropmissing!(df_test, :mean_hours)
df_test[!, :mean_hours_ma] = moving_avg(df_test.mean_hours, 4)
dropmissing!(df_test)
df_test = innerjoin(df_test, shock_df; on=:quarter)
df_test = innerjoin(df_test, macro_q;  on=:quarter)
dropmissing!(df_test)
println("After merge: $(nrow(df_test)) quarters, $(minimum(df_test.quarter)) - $(maximum(df_test.quarter))")

# 4. 看 shock 序列和 hours 的相关性
df_merged_check = innerjoin(
    filter(r -> r.occ_cat == 3, df_agg),
    shock_df; on=:quarter)
dropmissing!(df_merged_check, [:mean_hours, :markup_shock])
println("Correlation (shock, hours) Occ 3: ", 
        cor(df_merged_check.mean_hours, df_merged_check.markup_shock))


# shock 序列的基本统计
println("Mean:   ", mean(markup_shocks_m3_level))
println("Std:    ", std(markup_shocks_m3_level))

# shock 在 CPS 时间范围内的子集
shock_sub = filter(r -> r.quarter >= 1994.0 && r.quarter <= 2023.0, shock_df)
println("Shock obs in CPS window: ", nrow(shock_sub))

# CPS hours 在该窗口内和 shock 的相关性（直接，不做MA）
df_corr = innerjoin(filter(r -> r.occ_cat == 3, df_agg), shock_sub; on=:quarter)
dropmissing!(df_corr, [:mean_hours, :markup_shock])
println("Correlation in window: ", cor(df_corr.mean_hours, df_corr.markup_shock))
println("Shock std in window:   ", std(shock_sub.markup_shock))

df_ht = filter(r -> r.occ_cat == 3, df_agg)
dropmissing!(df_ht, :mean_hours)
sort!(df_ht, :quarter)
df_ht[!, :hours_raw] = exp.(df_ht.mean_hours)
println(df_ht[1:20, [:quarter, :hours_raw]])

tfp_shock_df = DataFrame(
    quarter      = year.(macro_df_method3.observation_date) .+ 
                   (quarterofyear.(macro_df_method3.observation_date) .- 1) ./ 4,
    markup_shock = macro_df_method3.tfp_util  # 直接用tfp_util作为shock
)

df_corr = innerjoin(filter(r -> r.occ_cat == 3, df_agg), tfp_shock_df; on=:quarter)
dropmissing!(df_corr, [:mean_hours, :markup_shock])
println("Correlation (TFP shock, hours) Occ 3: ",
        cor(df_corr.mean_hours, df_corr.markup_shock))

# 看 aggregate hours 的时序波动
df_all = combine(groupby(df_agg, :quarter),
                 :mean_hours => mean => :agg_hours)
sort!(df_all, :quarter)

# 和 FRED 的实际工时对比
# 正常的 aggregate hours 在2008-2009应该有明显下降
println(df_all[80:100, :])  # 大概2008-2012年附近

df_corr2 = innerjoin(filter(r -> r.occ_cat == 3, df_agg), tfp_shock_df; on=:quarter)
dropmissing!(df_corr2, [:mean_wage, :mean_income, :markup_shock])
println("Correlation (TFP, wage)   Occ 3: ", cor(df_corr2.mean_wage,   df_corr2.markup_shock))
println("Correlation (TFP, income) Occ 3: ", cor(df_corr2.mean_income, df_corr2.markup_shock))

# 同时看 income 的时序波动，确认它有周期信号
df_inc = combine(groupby(df_agg, :quarter), :mean_income => mean => :agg_income)
sort!(df_inc, :quarter)
println("\nagg_income std: ", std(dropmissing(df_inc).agg_income))
println(dropmissing(df_inc)[70:90, :])  # 2008附近

# 找2008-2010附近的行
df_inc_full = combine(groupby(df_agg, :quarter), :mean_income => mean => :agg_income)
sort!(df_inc_full, :quarter)
dropmissing!(df_inc_full)

# 找2007-2011
idx = findall(q -> 2007.0 <= q <= 2011.0, df_inc_full.quarter)
println(df_inc_full[idx, :])

# 同时看全样本的 min/max
println("\nFull sample:")
println("min: ", minimum(df_inc_full.agg_income), " at ", 
        df_inc_full.quarter[argmin(df_inc_full.agg_income)])
println("max: ", maximum(df_inc_full.agg_income), " at ",
        df_inc_full.quarter[argmax(df_inc_full.agg_income)])