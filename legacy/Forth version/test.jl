cd(@__DIR__)
include("sign restriction_forth.jl")
using DataFrames, Plots

# 提取历史冲击序列（Historical Structural Shocks）

# 从SR结果中获取识别矩阵
# SRout[:Ball] 包含所有被接受的旋转矩阵
# 我们使用中位数对应的旋转矩阵

# 获取中位数对应的旋转索引
n_accepted = size(SRout[:Ball], 3)
median_idx = div(n_accepted, 2)

# 提取中位数旋转矩阵
B_median = SRout[:Ball][:, :, median_idx]

# 计算历史结构冲击
# 结构冲击 = B^(-1) * 残差
residuals = VAR[:residuals]  # T x n 矩阵
structural_shocks = residuals * inv(B_median)'  # T x n 矩阵

# 提取价格加成冲击（第1列，因为这是你识别的第一个冲击）
markup_shocks = structural_shocks[:, 1]

# 创建包含日期和冲击的DataFrame
shock_df = DataFrame(
    date = final_df.observation_date[1:size(markup_shocks, 1)],
    markup_shock = markup_shocks
)

# 找出最大的正向和负向冲击
sorted_shocks = sort(shock_df, :markup_shock, rev=true)

println("\n" * "="^60)
println("历史价格加成冲击分析")
println("="^60)

println("\n前10个最大正向价格加成冲击:")
println(first(sorted_shocks, 10))

println("\n前10个最大负向价格加成冲击:")
println(last(sorted_shocks, 10))

# 绘制历史冲击时间序列
shock_plot = plot(shock_df.date, shock_df.markup_shock,
    title="Historical Price Markup Shocks",
    xlabel="Date",
    ylabel="Shock Size",
    label="Markup Shock",
    linewidth=1.5,
    legend=:topright)
hline!([0], linestyle=:dash, color=:black, label="")
display(shock_plot)

# 统计信息
println("\n冲击统计:")
println("均值: ", round(mean(markup_shocks), digits=4))
println("标准差: ", round(std(markup_shocks), digits=4))
println("最大值: ", round(maximum(markup_shocks), digits=4), 
        " (日期: ", shock_df.date[argmax(markup_shocks)], ")")
println("最小值: ", round(minimum(markup_shocks), digits=4),
        " (日期: ", shock_df.date[argmin(markup_shocks)], ")")

# 识别显著冲击（超过1个标准差）
threshold = std(markup_shocks)
significant_shocks = shock_df[abs.(shock_df.markup_shock) .> threshold, :]
sort!(significant_shocks, :markup_shock, rev=true)

println("\n显著冲击（|冲击| > 1个标准差）:")
println(significant_shocks)

# 保存结果
CSV.write("output/historical_markup_shocks.csv", shock_df)
println("\n历史冲击序列已保存至 output/historical_markup_shocks.csv")