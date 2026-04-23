using DataFrames
using CSV
using XLSX
using StatFiles
using Statistics
using Plots
using GLM

"""
    calculate_cumulative_irf_single(input_file, target_outcome; data_dir)
读取 IRF 轨迹数据，按 outcome 过滤，计算 horizon 1-36 的累积 beta，并进行标准化 (Z-score)。
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
    prepare_y_axis_data()
计算职业组特征（Abstract, Routine, Manual）的加权平均值，作为绘图的 Y 轴数据。
"""
function prepare_y_axis_data(
    dta_path::String = "../../data/occ1990dd_task_alm.dta",
    xlsx_path::String = "../../result/mapping/mapping_done.xlsx"
)
    # 1. 读取 DTA 文件 (任务特征数据)
    isfile(dta_path) || error("未找到 DTA 文件: $dta_path")
    df_tasks = DataFrame(load(dta_path))
    task_cols = [:occ1990dd, :task_abstract, :task_routine, :task_manual]
    select!(df_tasks, task_cols)
    
    # 2. 读取 Excel 文件 (映射和权重数据)
    isfile(xlsx_path) || error("未找到 Excel 文件: $xlsx_path")
    df_mapping = DataFrame(XLSX.readtable(xlsx_path, "Sheet1"))
    
    # 使用编号提取：B列(2)->occ1990dd, F列(6)->group, H列(8)->weights
    select!(df_mapping, 2 => :occ1990dd, 6 => :group, 8 => :weights)
    df_mapping.group = convert.(Int, df_mapping.group)
    df_mapping.weights = convert.(Float64, df_mapping.weights)
    
    # 3. 合并数据表 (Inner Join)
    df_merged = innerjoin(df_tasks, df_mapping, on = :occ1990dd)
    dropmissing!(df_merged)
    
    # 4. 按 group 计算加权平均 (Weighted Mean)
    calc_wmean(x, w) = sum(w) == 0 ? 0.0 : sum(x .* w) / sum(w)
    df_y = combine(groupby(df_merged, :group),
        [:task_abstract, :weights] => calc_wmean => :task_abstract,
        [:task_routine, :weights]  => calc_wmean => :task_routine,
        [:task_manual, :weights]   => calc_wmean => :task_manual
    )
    
    # 5. 将 Group (1-9) 映射成与 X 轴数据相同的名字
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
    df_y.Occupational_Group = [group_map[g] for g in df_y.group]
    select!(df_y, :Occupational_Group, :task_abstract, :task_routine, :task_manual)
    
    return df_y
end

"""
    plot_correlation_with_ci()
绘制带有 90% 置信区间阴影的回归散点图（防文字重叠、去图例、扩大边距）
"""
function plot_correlation_with_ci(df_merged::DataFrame, x_col::Symbol, y_col::Symbol;
                                  x_label::String = "Cumulative IRF Beta",
                                  y_label::String = "Occupation Characteristic")
    
    df_plot = copy(df_merged)
    
    # 长名字到短名字的映射
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
    df_plot.Plot_Label = [get(plot_name_map, name, name) for name in df_plot.Occupational_Group]

    dropmissing!(df_plot, [x_col, y_col])
    
    # 使用 GLM 拟合 OLS 线性模型 (Y ~ X)
    fm = term(y_col) ~ term(x_col)
    model = lm(fm, df_plot)
    
    x_min, x_max = extrema(df_plot[!, x_col])
    y_min, y_max = extrema(df_plot[!, y_col])
    
    x_span = x_max - x_min == 0 ? 1.0 : (x_max - x_min)
    y_span = y_max - y_min == 0 ? 1.0 : (y_max - y_min)
    
    padding = x_span * 0.1 
    x_range = range(x_min - padding, x_max + padding * 2.5, length=100)
    
    df_pred = DataFrame()
    df_pred[!, x_col] = x_range
    
    # 预测出拟合值和 90% 置信区间
    pred = predict(model, df_pred, interval=:confidence, level=0.9)
    
    # 设置全局学术风格
    default(fontfamily="sans-serif", framestyle=:box, grid=false, 
            tickfontsize=9, guidefontsize=11)
    
    # 第 1 层：画阴影和拟合线 (隐藏图例，留出边距防止被切)
    p = plot(x_range, pred.prediction,
             ribbon = (pred.prediction .- pred.lower, pred.upper .- pred.prediction),
             fillalpha = 0.2, fillcolor = :gray,
             linecolor = :blue, linewidth = 2,
             legend = false, 
             xlims = (x_min - padding, x_max + padding * 2.5),
             left_margin = 8Plots.mm,
             bottom_margin = 6Plots.mm,
             right_margin = 5Plots.mm)
             
    # 第 2 层：画散点图 (隐藏 label)
    scatter!(p, df_plot[!, x_col], df_plot[!, y_col],
             markershape = :circle,
             markercolor = :white,
             markerstrokecolor = :red,
             markerstrokewidth = 1.5,
             markersize = 6,
             label = "")
             
    # 第 3 层：加上防重叠的文本标签
    sort!(df_plot, x_col)
    
    for i in 1:nrow(df_plot)
        row = df_plot[i, :]
        x_shift = x_span * 0.015
        y_shift = y_span * 0.02
        v_align = :bottom
        
        # 简单防重叠逻辑
        if i > 1
            prev_row = df_plot[i-1, :]
            if abs(row[x_col] - prev_row[x_col]) < (x_span * 0.15) && 
               abs(row[y_col] - prev_row[y_col]) < (y_span * 0.15)
                y_shift = -y_span * 0.02
                v_align = :top
            end
        end
        
        annotate!(p, row[x_col] + x_shift, row[y_col] + y_shift, 
                  text(row[:Plot_Label], 8, :left, v_align, :black))
    end
    
    xlabel!(p, x_label)
    ylabel!(p, y_label)
    
    return p
end

"""
    main()
循环执行所有 outcome，并将同一 outcome 的三个 task 合并成一个 PDF
"""
function main()
    outcome_list = [
        "employment", "hourly_rate", "income_share",
        "inequality", "median", "unemployment", "income", "hours"
    ]

    output_dir = "../../result/occ/correlation/task"
    mkpath(output_dir)

    println("正在准备 Y 轴特征数据...")
    df_y = prepare_y_axis_data()

    task_list = [
        (:task_abstract, "Abstract Task Weight"),
        (:task_routine,  "Routine Task Weight"),
        (:task_manual,   "Manual Task Weight")
    ]

    for outcome in outcome_list
        println("\n--- 开始处理 Outcome: $outcome ---")
        
        # 1. 提取 X 轴
        df_x = calculate_cumulative_irf_single("merged_occ_irf_trajectories.csv", outcome)
        
        # 2. 合并数据
        df_merged = innerjoin(df_x, df_y, on = :Occupational_Group)
        
        plots_array = []
        
        # 3. 循环画 3 张子图
        for (task_col, task_label) in task_list
            p = plot_correlation_with_ci(
                df_merged, 
                Symbol(outcome), 
                task_col;
                x_label = "Cumulative IRF Beta ($outcome)",
                y_label = task_label
            )
            push!(plots_array, p)
        end
        
        # 4. 拼合保存 PDF
        pdf_path = joinpath(output_dir, "$(outcome)_all_tasks.pdf")
        png_path = joinpath(output_dir, "$(outcome)_all_tasks.png")
        
        combined_plot = plot(plots_array..., 
                             layout = (3, 1), 
                             size = (700, 1200),
                             plot_title = "Correlation: $outcome vs Tasks")
        
        savefig(combined_plot, pdf_path)
        savefig(combined_plot, png_path)
        println(">>> 成功生成合并版 PDF: $(outcome)_all_tasks.pdf <<<")
    end
    
    println("\n全部 8 个 Outcome 处理完毕！")
end

# 运行主程序
cd(@__DIR__)
main()