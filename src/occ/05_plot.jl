using CSV
using DataFrames
using CairoMakie

# 设置路径（请确保与之前脚本生成的路径一致）
const BASE_DIR = normpath(joinpath(@__DIR__, "..", "..", "result", "occ"))
const OUTPUT_DIR = joinpath(BASE_DIR, "analysis")
const EXTREMES_CSV_PATH = joinpath(OUTPUT_DIR, "occupation_irf_extremes_turning_points.csv")

"""
    visualize_extremes_and_turning_points(csv_path, output_dir)

读取极值与拐点数据，为每个 Outcome 生成一张“双面板联动”可视化图表。
"""
function visualize_extremes_and_turning_points(csv_path::String, output_dir::String)
    if !isfile(csv_path)
        error("找不到数据文件: $csv_path \n请先运行之前提取极值的脚本。")
    end

    df = CSV.read(csv_path, DataFrame)
    mkpath(output_dir)

    # 遍历每个宏观变量 (Outcome) 生成单独的图
    for sub_df in groupby(df, :outcome)
        outcome_name = first(sub_df.outcome)
        
        # 按照最大冲击的实际数值（从小到大，即负面冲击最深的在最下面）进行排序
        sort!(sub_df, :max_abs_shock)
        
        occupations = sub_df.occupation
        n = length(occupations)
        y_positions = 1:n

        # 初始化画布
        fig = Figure(size = (1200, 600), fontsize = 14)

        # ==========================================
        # 面板 1：最大冲击幅度 (Bar Chart)
        # ==========================================
        ax1 = Axis(fig[1, 1], 
            title = "Maximum Shock Magnitude ($outcome_name)",
            xlabel = "Shock Amplitude (β)",
            yticks = (y_positions, occupations),
            yreversed = false # 负面冲击最严重的在底部，视觉上更有重量感
        )

        # 根据数值正负赋予不同颜色（红色代表受损，蓝色代表正向响应）
        bar_colors = [val < 0 ? :indianred : :steelblue for val in sub_df.max_abs_shock]
        
        barplot!(ax1, y_positions, sub_df.max_abs_shock, 
                 direction = :x, color = bar_colors)

        # 添加一条 0 轴的辅助线
        vlines!(ax1, 0, color = :black, linewidth = 1.5, linestyle = :dash)

        # ==========================================
        # 面板 2：极值与拐点的时间轴 (Timeline)
        # ==========================================
        ax2 = Axis(fig[1, 2], 
            title = "Timeline: Extremes & Turning Points",
            xlabel = "Horizon (Months)",
            yticks = (y_positions, occupations),
            yticklabelsvisible = false, # 隐藏 y 轴标签，因为左侧已经有了
            limits = (-1, 37, nothing, nothing),
            xticks = 0:6:36 # 每半年显示一个刻度
        )

        # 为每个职业画一条浅色的水平辅助线，增强对齐感
        hlines!(ax2, y_positions, color = (:gray, 0.2), linestyle = :dash)

        # 逐行绘制极值点和拐点
        for (i, row) in enumerate(eachrow(sub_df))
            # 1. 绘制拐点 (Turning Points) - 橙色圆点
            if row.turning_point_horizons != "None"
                # 将 "6, 12, 24" 这种字符串解析为整数数组
                tps_strs = split(row.turning_point_horizons, ",")
                tps = parse.(Int, strip.(tps_strs))
                
                scatter!(ax2, tps, fill(i, length(tps)), 
                         color = :orange, markersize = 12, marker = :circle)
            end

            # 2. 绘制极值发生月份 (Max Shock Horizon) - 红色五角星
            scatter!(ax2, [row.max_shock_horizon], [i], 
                     color = :red, markersize = 20, marker = :star5)
        end

        # ==========================================
        # 添加全局图例
        # ==========================================
        legend_elements = [
            MarkerElement(color = :red, marker = :star5, markersize = 16),
            MarkerElement(color = :orange, marker = :circle, markersize = 12)
        ]
        Legend(fig[1, 3], legend_elements, ["Max Shock Month", "Turning Point / Reversal"], "Legend")

        # 联动左右两个面板的 Y 轴，确保它们绝对对齐
        linkyaxes!(ax1, ax2)

        # 调整布局间距
        colgap!(fig.layout, 15)

        # 保存图片
        save_path_png = joinpath(output_dir, "viz_extremes_timeline_$(outcome_name).png")
        save_path_pdf = joinpath(output_dir, "viz_extremes_timeline_$(outcome_name).pdf")
        save(save_path_png, fig)
        save(save_path_pdf, fig)
        
        println("Generated Dual-Panel plot for [$outcome_name] -> $save_path_png")
    end
end

# 执行绘图
visualize_extremes_and_turning_points(EXTREMES_CSV_PATH, OUTPUT_DIR)