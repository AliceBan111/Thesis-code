# Include utility functions
include("src/ToolBox/load_var.jl")
include("src/data/build_macro_data.jl")
include("src/data/build_macro_data_with_TFP.jl")
include("src/data/build_markup_data_method2.jl")
include("src/data/build_markup_data_method2_no_hp.jl")
include("src/data/build_markup_data_method2_no_hp_latest.jl")
include("src/models/var_sr.jl")
include("src/models/var.jl")
include("src/wage_shock_analysis/extract_shocks.jl")
include("src/wage_shock_analysis/wage_heterogeneity.jl")
include("src/wage_shock_analysis/local_projection.jl")

using CSV, DataFrames, Statistics, Plots, Dates, Downloads, SHA, LocalProjections

cd(@__DIR__)

# 1. Load macro data
# start_date_method2 = Date("1987-01-01")
# end_date_method2   = Date("2012-12-31")

start_date_method2_latest = Date("1997-01-01")
end_date_method2_latest   = Date("2023-12-31")

# macro_df_method2 = build_macro_data(start_date_method2, end_date_method2)
# macro_df_method2_no_hp = build_macro_data(start_date_method2, end_date_method2)
macro_df_method2_latest = build_macro_data_with_TFP(start_date_method2_latest, end_date_method2_latest)

# 2. Build markup data
# df_method2 = build_markup_method2(macro_df_method2)
# df_method2_no_hp = build_markup_method2_no_hp(macro_df_method2_no_hp)
df_method2_no_hp_latest = build_markup_method2_no_hp_latest(macro_df_method2_latest, "data/KLEMS_latest.xlsx")

# 2.5 Plot VAR inputs
plot_VAR_inputs(df_method2_no_hp_latest; tag="method2_no_hp_latest", use_growth  = true)
plot_VAR_inputs(df_method2_no_hp_latest; tag="method2_no_hp_latest", use_growth  = false)

# 3. Estimate VAR + Sign Restrictions
# res_method2 = estimate_VAR_SR(df_method2; method_name="method2")
# res_method2_no_hp_level = estimate_VAR_SR(df_method2_no_hp; method_name="method2_no_hp", use_growth = false)


res_method2_no_hp_latest_growth = estimate_VAR_SR(df_method2_no_hp_latest; method_name="method2_no_hp_latest", use_growth  = true)
res_method2_no_hp_latest_level = estimate_VAR_SR(df_method2_no_hp_latest; method_name="method2_no_hp_latest", use_growth  = false)
# 4. extract shocks
markup_shocks_m2_growth = extract_structural_shocks(res_method2_no_hp_latest_growth, shock_index=1)
n_shocks_growth = length(markup_shocks_m2)

markup_shocks_m2_level = extract_structural_shocks(res_method2_no_hp_latest_level, shock_index=1)
n_shocks_level = length(markup_shocks_m2_level)

# Get the corresponding dates for the shocks
shock_dates = last(df_method2_no_hp_latest.observation_date, n_shocks)

println("Shocks length: ", n_shocks)
println("Dates length: ", length(shock_dates))

shock_dates_level = last(df_method2_no_hp_latest.observation_date, n_shocks_level)

println("Shocks length: ", n_shocks_level)
println("Dates length: ", length(shock_dates_level))

# 5. Run heterogeneity analysis
df_results, irf_data = run_wage_heterogeneity_analysis(
    markup_shocks_m2,
    shock_dates;
    output_suffix="method2"
)

df_results, irf_data = run_wage_heterogeneity_analysis(
    markup_shocks_m2_level,
    shock_dates_level;
    output_suffix="method2_no_hp_level"
)

# 6. Local Projection
df_lp_growth, r_white_g, r_blue_g, irf_white_g, irf_blue_g = run_wage_lp_analysis(
    markup_shocks_m2,
    shock_dates;
    output_suffix = "method2_latest_growth",
    nhorz = 8,
    nlag  = 4
)

df_lp_level, r_white_l, r_blue_l, irf_white_l, irf_blue_l = run_wage_lp_analysis(
    markup_shocks_m2_level,
    shock_dates_level;
    output_suffix = "method2_latest_level",
    nhorz = 8,
    nlag  = 4
)



# 1. 基本统计
println("Mean:   ", mean(markup_shocks_m2_level))
println("Std:    ", std(markup_shocks_m2_level))
println("Min:    ", minimum(markup_shocks_m2_level))
println("Max:    ", maximum(markup_shocks_m2_level))

# 2. 画图：时序 + 直方图
p1 = plot(shock_dates_level, markup_shocks_m2_level,
    title = "Markup Shock (level)", xlabel = "Date", ylabel = "Shock",
    legend = false, linecolor = :blue)
hline!([0], linestyle = :dash, color = :black)

p2 = histogram(markup_shocks_m2_level,
    title = "Shock Distribution", xlabel = "Shock value",
    bins = 20, legend = false)

plot(p1, p2, layout = (2,1), size = (800, 600))


# 3. 检查是否近似均值为零、无自相关（结构冲击应满足）
using HypothesisTests
println(OneSampleTTest(markup_shocks_m2_level))   # 均值是否为零


# 先看一下结果对象的结构
dump(r_white_l)
# 或者
fieldnames(typeof(r_white_l))

println("xnames: ", r_white_l.xnames)
println("ynames: ", r_white_l.ynames)
println("B size: ", size(r_white_l.B))
println("B = ", r_white_l.B)

# B: (6, 1, 8) -> B[x_idx, y_idx, horz]
# V: 需要先确认维度
println("V size: ", size(r_white_l.V))


shock_idx = 1  # markup_shock

betas = [r_white_l.B[shock_idx, 1, h] for h in 1:8]
ses   = [sqrt(r_white_l.V[shock_idx, shock_idx, h]) for h in 1:8]

println("\nWhite-collar LP results:")
println("h  |  β            |  SE           |  t-stat  |  MDE(95%)")
for h in 1:8
    t   = betas[h] / ses[h]
    mde = 1.96 * ses[h]
    println("$h  |  $(round(betas[h], sigdigits=4))  |  $(round(ses[h], sigdigits=4))  |  $(round(t, sigdigits=3))  |  $(round(mde, sigdigits=4))")
end

# 联合 Wald 检验（使用完整协方差矩阵）
using LinearAlgebra, Distributions
b_vec = betas
# 跨 horizon 的联合协方差：取各 horizon 对角块拼成块对角
# 保守做法：只用各自对角元素（忽略跨horizon相关）
V_diag = Diagonal([r_white_l.V[shock_idx, shock_idx, h] for h in 1:8])
W = b_vec' * inv(V_diag) * b_vec
df = 8
p_wald = 1 - cdf(Chisq(df), W)
println("\nJoint Wald statistic: ", round(W, sigdigits=4))
println("p-value (χ²($df)):    ", round(p_wald, sigdigits=4))


println("ynames: ", r_blue_l.ynames)
println("B size: ", size(r_blue_l.B))

shock_idx = 1

betas_blue = [r_blue_l.B[shock_idx, 1, h] for h in 1:8]
ses_blue   = [sqrt(r_blue_l.V[shock_idx, shock_idx, h]) for h in 1:8]

println("\nBlue-collar LP results:")
println("h  |  β            |  SE           |  t-stat  |  MDE(95%)")
for h in 1:8
    t   = betas_blue[h] / ses_blue[h]
    mde = 1.96 * ses_blue[h]
    println("$h  |  $(round(betas_blue[h], sigdigits=4))  |  $(round(ses_blue[h], sigdigits=4))  |  $(round(t, sigdigits=3))  |  $(round(mde, sigdigits=4))")
end

using LinearAlgebra, Distributions
V_diag_blue = Diagonal([r_blue_l.V[shock_idx, shock_idx, h] for h in 1:8])
W_blue = betas_blue' * inv(V_diag_blue) * betas_blue
df = 8
p_wald_blue = 1 - cdf(Chisq(df), W_blue)
println("\nJoint Wald statistic: ", round(W_blue, sigdigits=4))
println("p-value (χ²($df)):    ", round(p_wald_blue, sigdigits=4))