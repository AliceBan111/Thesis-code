using DataFrames
using CSV
using XLSX
using StatFiles
using Statistics
using Plots
using GLM


"""
    calculate_cumulative_irf_single(input_file, target_outcome; data_dir)

Reads the merged IRF trajectories CSV, filters for `target_outcome`, sums betas
over horizons 1-36 for each occupational group, and returns a 2-column DataFrame:
  :Occupational_Group  |  Symbol(target_outcome)
"""

function calculate_cumulative_irf_single(input_file::String, target_outcome::String;
                                         data_dir::String = "../../result/occ/analysis",
                                         standardize::Bool = true)

    filepath = joinpath(data_dir, input_file)
    isfile(filepath) || error("File not found: $filepath")

    df = CSV.read(filepath, DataFrame; missingstring=["", "NA", "N/A"])

    # Filter for target outcome
    filter!(row -> row.outcome == target_outcome, df)
    nrow(df) > 0 || error("Outcome '$target_outcome' not found in $input_file")

    # Identify group columns
    exclude = Set([:outcome, :horizon])
    grp_cols = [c for c in propertynames(df) if !(c in exclude)]
    isempty(grp_cols) && error("No group columns found in $input_file")

    # Wide -> Long
    df_long = stack(df, grp_cols,
                    variable_name = :Occupational_Group,
                    value_name = :Beta)

    # Keep horizons 1-36
    filter!(row -> !ismissing(row.horizon) &&
                   row.horizon >= 1 &&
                   row.horizon <= 36, df_long)

    # Sum cumulative IRF
    df_cum = combine(groupby(df_long, :Occupational_Group),
                     :Beta => (b -> sum(skipmissing(b))) => Symbol(target_outcome))

    # Standardization (z-score)
    if standardize
        col = Symbol(target_outcome)
        vals = df_cum[!, col]

        μ = mean(skipmissing(vals))
        σ = std(skipmissing(vals))

        if σ == 0
            df_cum[!, col] .= 0.0
        else
            df_cum[!, col] = (vals .- μ) ./ σ
        end
    end

    sort!(df_cum, :Occupational_Group)
    return df_cum
end

"""
计算职业组特征（Abstract, Routine, Manual）的加权平均值
作为绘图的 Y 轴数据
"""
function prepare_y_axis_data(
    dta_path::String = "../../data/occ1990dd_task_alm.dta",
    xlsx_path::String = "../../result/mapping/mapping_done.xlsx"
)
    # ---------------------------------------------------------
    # 1. 读取 DTA 文件 (任务特征数据)
    # ---------------------------------------------------------
    isfile(dta_path) || error("未找到 DTA 文件: $dta_path")
    df_tasks = DataFrame(load(dta_path))
    
    # 确保列名匹配并只保留需要的列
    task_cols = [:occ1990dd, :task_abstract, :task_routine, :task_manual]
    select!(df_tasks, task_cols)
    
    # ---------------------------------------------------------
    # 2. 读取 Excel 文件 (映射和权重数据)
    # ---------------------------------------------------------
    isfile(xlsx_path) || error("未找到 Excel 文件: $xlsx_path")
    
    # 读取第一个 sheet，并转换为 DataFrame
    # 假设你的第一行是表头
    df_mapping = DataFrame(XLSX.readtable(xlsx_path, "Sheet1"))
    
    select!(df_mapping, 2 => :occ1990dd, 6 => :group, 8 => :weights)
    
    # 确保格式正确（比如从 Excel 读进来的可能是 Any 类型，转为具体类型）
    df_mapping.group = convert.(Int, df_mapping.group)
    df_mapping.weights = convert.(Float64, df_mapping.weights)
    
    # ---------------------------------------------------------
    # 3. 合并数据表 (Inner Join)
    # ---------------------------------------------------------
    df_merged = innerjoin(df_tasks, df_mapping, on = :occ1990dd)
    
    # 剔除存在缺失值的行，确保加权平均计算不出错
    dropmissing!(df_merged)
    
    # ---------------------------------------------------------
    # 4. 按 group 计算加权平均 (Weighted Mean)
    # ---------------------------------------------------------
    # 加权平均公式: sum(x * weight) / sum(weight)
    calc_wmean(x, w) = sum(w) == 0 ? 0.0 : sum(x .* w) / sum(w)
    
    df_y = combine(groupby(df_merged, :group),
        [:task_abstract, :weights] => calc_wmean => :task_abstract,
        [:task_routine, :weights]  => calc_wmean => :task_routine,
        [:task_manual, :weights]   => calc_wmean => :task_manual
    )
    
    # ---------------------------------------------------------
    # 5. 将 Group (1-9) 映射成与 X 轴数据相同的名字
    # ---------------------------------------------------------
    # 这一步非常关键，它让你之后可以直接按 :Occupational_Group 把 X轴和 Y轴合并！
    group_map = Dict(
        1 => "group1_Managerial",
        2 => "group2_Professional_specialty",
        3 => "group3_High_tech",
        4 => "group4_Sales",
        5 => "group5_Administrative_support",
        6 => "group6_Service",
        7 => "group7_Farming_forestry_construction",
        8 => "group8_Precision_production_repair",
        9 => "group9_Machine_operators_transport"
    )
    
    # 根据 group 数字生成对应的组名，作为新的列
    df_y.Occupational_Group = [group_map[g] for g in df_y.group]
    
    # 去掉不再需要的数字 group 列，调整列顺序
    select!(df_y, :Occupational_Group, :task_abstract, :task_routine, :task_manual)
    
    return df_y
end

function plot_correlation_with_ci(df_merged::DataFrame, x_col::Symbol, y_col::Symbol;
                                  out_file::String = "correlation_plot.png",
                                  x_label::String = "Cumulative IRF Beta",
                                  y_label::String = "Occupation Characteristic")
    
    # 1. 复制数据，避免修改原表
    df_plot = copy(df_merged)
    
    # 2. 在画图这一步，实现长名字到短名字的映射
    plot_name_map = Dict(
        "group1_Managerial" => "Managerial",
        "group2_Professional_specialty" => "Prof. specialty",
        "group3_High_tech" => "High Tech",
        "group4_Sales" => "Sales",
        "group5_Administrative_support" => "Admin",
        "group6_Service" => "Service",
        "group7_Farming_forestry_construction" => "Constr., Extract., Farm",
        "group8_Precision_production_repair" => "Prod, Repair",
        "group9_Machine_operators_transport" => "Machineists, Transp."
    )
    
    # 生成用于图表显示的 Plot_Label 列
    # 使用 get() 函数提供默认值防报错：如果找不到匹配，就用原名
    df_plot.Plot_Label = [get(plot_name_map, name, name) for name in df_plot.Occupational_Group]

    # 3. 剔除含有缺失值的行
    dropmissing!(df_plot, [x_col, y_col])
    
    # 4. 使用 GLM 拟合 OLS 线性模型 (Y ~ X)
    fm = term(y_col) ~ term(x_col)
    model = lm(fm, df_plot)
    
    # 生成平滑的 X 轴数据（左右稍微多延伸一点），用于画拟合线和阴影
    x_min, x_max = extrema(df_plot[!, x_col])
    y_min, y_max = extrema(df_plot[!, y_col])
    x_span = x_max - x_min == 0 ? 1.0 : (x_max - x_min)
    y_span = y_max - y_min == 0 ? 1.0 : (y_max - y_min)

    padding = x_span * 0.1 
    x_range = range(x_min - padding, x_max + padding * 2.5, length=100)
    
    df_pred = DataFrame()
    df_pred[!, x_col] = x_range
    
    # 预测出拟合值 (prediction) 和 90% 置信区间的上下限 (lower, upper)
    pred = predict(model, df_pred, interval=:confidence, level=0.9)
    
    # 5. 开始绘图 (Plots.jl)
    # 设置全局学术风格
    default(fontfamily="sans-serif", framestyle=:box, grid=false, tickfontsize=9, guidefontsize=11)
    
    # 第 1 层：画阴影和拟合线 (铺在底层)
    p = plot(x_range, pred.prediction,
             ribbon = (pred.prediction .- pred.lower, pred.upper .- pred.prediction),
             fillalpha = 0.2, fillcolor = :gray,
             linecolor = :blue, linewidth = 2,
             legend = false,
             xlims = (x_min - padding, x_max + padding * 2.5), # 扩展 X 轴画幅
             left_margin = 8Plots.mm,   # 解决 Y 轴标题被切
             bottom_margin = 6Plots.mm, # 解决 X 轴标题被切
             right_margin = 5Plots.mm)
             
    # 第 2 层：画散点图 (红色空心圆，盖在拟合线上)
    scatter!(p, df_plot[!, x_col], df_plot[!, y_col],
             markershape = :circle,
             markercolor = :white,       # 内部空心
             markerstrokecolor = :red,   # 边缘红色
             markerstrokewidth = 1.5,
             markersize = 6,
             label = "")
             
    # 第 3 层：加上文本标签 (短名字)
    # 计算一个小小的 Y 轴偏移量，防止文字和圆圈重合
    sort!(df_plot, x_col)
    
    for i in 1:nrow(df_plot)
        row = df_plot[i, :]
        
        # 默认：放在圆圈右上方
        x_shift = x_span * 0.015
        y_shift = y_span * 0.02
        v_align = :bottom
        
        # 防重叠逻辑：如果和前一个点靠得太近（水平方向和垂直方向都很近）
        if i > 1
            prev_row = df_plot[i-1, :]
            if abs(row[x_col] - prev_row[x_col]) < (x_span * 0.15) && 
               abs(row[y_col] - prev_row[y_col]) < (y_span * 0.15)
                # 就把当前文字放到圆圈右下方，强行错开
                y_shift = -y_span * 0.02
                v_align = :top
            end
        end
        
        # 对齐方式设为 :left，配合 x_shift，让文字整体位于点的右侧
        annotate!(p, row[x_col] + x_shift, row[y_col] + y_shift, 
                  text(row[:Plot_Label], 8, :left, v_align, :black))
    end
    
    # 设置标签和标题
    xlabel!(p, x_label)
    ylabel!(p, y_label)
    
    # 6. 保存图表
    savefig(p, out_file)
    println("Plots succesfully saved: ", out_file)
    
    return p
end



# ============================================================
#  MAIN – LOOP OVER ALL OUTCOMES × ALL O*NET TABLES
# ============================================================

function main()
    outcome_list = [
        "employment", "hourly_rate", "income_share",
        "inequality", "median", "unemployment", "income", "hours"
    ]

    output_dir = "../../result/occ/correlation/task"
    mkpath(output_dir)

    # 1. 提取全局的 Y 轴数据 (所有 Task 的特征和权重)
    # 因为不同 outcome 共享相同的职业特征数据，所以只需提取一次，提高运行效率
    println("正在准备 Y 轴特征数据...")
    df_y = prepare_y_axis_data()

    # 定义我们要遍历的三个 task 及其在图表上的 Y 轴显示标签
    task_list = [
        (:task_abstract, "Abstract Task Weight"),
        (:task_routine,  "Routine Task Weight"),
        (:task_manual,   "Manual Task Weight")
    ]

    # 2. 循环遍历每一个 outcome
    for outcome in outcome_list
        println("\n--- 开始处理 Outcome: $outcome ---")
        
        # 2.1 提取当前 outcome 的 X 轴数据 (累积 IRF Beta)
        # 注意：这里假设你的 input_file 是 "merged_occ_irf_trajectories.csv"
        df_x = calculate_cumulative_irf_single("merged_occ_irf_trajectories.csv", outcome)
        
        # 2.2 将 X 轴与 Y 轴数据按照 Occupational_Group 合并
        df_merged = innerjoin(df_x, df_y, on = :Occupational_Group)
        
        # 创建一个空数组，用于存放当前 outcome 生成的 3 个 plot 对象
        plots_array = []
        
        # 2.3 遍历三个任务特征，分别做回归和画图
        for (task_col, task_label) in task_list
            
            # 单独的图片保存路径 (由于你给的画图函数内置了保存，我们让它顺便保存单张图)
            single_plot_path = joinpath(output_dir, "$(outcome)_$(task_col).png")
            
            # 调用你的绘图函数，生成当前 task 的散点图
            p = plot_correlation_with_ci(
                df_merged, 
                Symbol(outcome), 
                task_col;
                out_file = single_plot_path,
                x_label = "Cumulative IRF Beta ($outcome)",
                y_label = task_label
            )
            
            # 把生成的 plot 对象塞进数组里，留着后面拼装
            push!(plots_array, p)
        end
        
        # 2.4 将三张图拼合在一起，并保存为一个 PDF 文件
        pdf_path = joinpath(output_dir, "$(outcome)_all_tasks.pdf")
        
        # 使用 layout 将 3 个 plot 纵向排列 (3 行 1 列)
        # size 参数设置了整个 PDF 的画布大小，宽 700 高 1200，保证纵向空间充足，标签不重叠
        combined_plot = plot(plots_array..., 
                             layout = (3, 1), 
                             size = (700, 1200),
                             plot_title = "Correlation: $outcome vs Tasks")
        
        # 导出包含这 3 个子图的最终 PDF
        savefig(combined_plot, pdf_path)
    end

    println("\n✓ Pipeline complete.")
end

cd(@__DIR__)
main()