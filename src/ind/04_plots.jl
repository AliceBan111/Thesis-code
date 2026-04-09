# =============================================================================
# 05_plots_unified.jl
# IRF plots for all industry groups (统一动态版 / 无基准组方案)
# 特性：动态Y轴刻度对齐、动态标题映射、13个行业专属颜色映射、自动网格布局
# =============================================================================

using CSV, DataFrames, DataFramesMeta
using CairoMakie
using Printf

# =============================================================================
# 0. 全局配置与动态映射
# =============================================================================
# 默认的行业分组标签 (13个组)
const IND_LABELS = Dict(
    1  => "Agriculture_forestry_fishing",  2  => "Mining",
    3  => "Construction",                  4  => "Manufacturing_nondurable",
    5  => "Manufacturing_durable",         6  => "Transportation_utilities",
    7  => "Wholesale_trade",               8  => "Retail_trade",
    9  => "Finance_insurance_realestate",  10 => "Business_repair_services",
    11 => "Personal_entertainment_services", 12 => "Professional_related_services",
    13 => "Public_administration",
)

# 13 个行业专属颜色
const GROUP_COLORS = [
    :steelblue, :tomato, :seagreen, :darkorange,
    :mediumpurple, :saddlebrown, :hotpink, :teal, :goldenrod,
    :navy, :crimson, :limegreen, :violet
]

# 动态获取当前变量在图表上显示的文本名称
function get_plot_labels(variant::Symbol)
    labels = Dict(
        :hourly_rate      => "Response of Log Hourly Wage",
        :unemployment     => "Response of Unemployment Rate (%)",
        :hours            => "Response of Hours Worked",
        :income           => "Response of Log Income",
        :employment       => "Response of Log Employment",
        :inequality       => "Response of Inequality",
        :median           => "Response of Median Wage",
        :income_share_var => "Response of Income Share"
    )
    return get(labels, variant, "Response of $(string(variant))")
end

# 动态计算全局 Y 轴刻度 (以 95% CI 的上下限为准，加 10% 的 padding)
function get_global_ylims(irfs::Dict{Int, DataFrame})
    g_min, g_max = Inf, -Inf
    
    for g in keys(IND_LABELS)
        !haskey(irfs, g) && continue
        df = irfs[g]
        # 排除可能存在的 NaN
        valid_min = minimum(filter(!isnan, df[!, :ci_lo95]))
        valid_max = maximum(filter(!isnan, df[!, :ci_hi95]))
        
        g_min = min(g_min, valid_min)
        g_max = max(g_max, valid_max)
    end
    
    padding = (g_max - g_min) * 0.1
    # 如果全平 (比如全是0)，给个默认范围
    if padding == 0 
        padding = 0.1
    end
    
    return (g_min - padding, g_max + padding)
end

# =============================================================================
# 1. 绘制绝对 IRF (beta_abs)
# =============================================================================
function plot_absolute_irf(irfs::Dict{Int, DataFrame}, pw_bh::DataFrame, variant::Symbol, output_dir::String)
    
    y_label = get_plot_labels(variant)
    ylims_abs = get_global_ylims(irfs)
    
    # 因为有13个行业，我们使用 4列 x 4行 的网格，调大画布尺寸
    fig = Figure(resolution = (1600, 1200))
    
    group_ids = sort(collect(keys(IND_LABELS)))
    for (idx, g) in enumerate(group_ids)
        !haskey(irfs, g) && continue
        df  = irfs[g]
        
        # 4列布局
        row = ceil(Int, idx / 4)
        col = mod1(idx, 4)
        
        # 获取当前行业的专属颜色
        base_color = GROUP_COLORS[g]
        
        ax = Axis(fig[row, col],
                  title  = get(IND_LABELS, g, "Group $g"),
                  xlabel = "Horizon (months)",
                  ylabel = col == 1 ? y_label : "") 
        
        ylims!(ax, ylims_abs[1], ylims_abs[2])
        hlines!(ax, [0.0], color = :black, linewidth = 1)
        
        # 使用专属颜色绘制置信区间 (带透明度)
        band!(ax, df.horizon, df.ci_lo90, df.ci_hi90, color = (base_color, 0.3))
        band!(ax, df.horizon, df.ci_lo68, df.ci_hi68, color = (base_color, 0.5))
        
        # 使用专属颜色绘制主线
        lines!(ax, df.horizon, df.beta_abs, color = base_color, linewidth = 2.5)
        
        # 显著性散点 (基于 pointwise BH 修正) 
        # 注: 散点保持红色或黑色五角星效果最好，能突出显著点
        sig_sub = filter(r -> r.ind_group == g && r.bh_reject, pw_bh)
        if nrow(sig_sub) > 0
            sig_points = innerjoin(sig_sub, df, on=:horizon, makeunique=true)
            scatter!(ax, sig_points.horizon, sig_points.beta_abs, 
                     color = :red, marker = :star5, markersize = 12)
        end
    end
    
    output_path = joinpath(output_dir, "irf_absolute.pdf")
    save(output_path, fig)
    println("Saved absolute IRFs to: $output_path")
    return fig
end

# =============================================================================
# 2. 绘制持久性摘要图 (Persistence Classification)
# =============================================================================
function plot_persistence_summary(sig_table::DataFrame, variant::Symbol, output_dir::String)
    
    # 颜色映射 (对应 03_significance_tests 里的分类标签)
    type_colors = Dict(
        "Persistent"         => :darkred,
        "Temporary"          => :steelblue,
        "Delayed_persistent" => :darkorange,
        "Not_significant"    => :gray
    )
    
    labels = [get(IND_LABELS, r.ind_group, "Group $(r.ind_group)") for r in eachrow(sig_table)]
    
    # 高度稍作增加以容纳13个标签
    fig = Figure(resolution = (1000, 700))
    ax  = Axis(fig[1, 1],
               yticks = (1:length(labels), labels),
               xlabel = "p-value",
               title  = "Persistence Classification: $(get_plot_labels(variant))")
               
    for (i, row) in enumerate(eachrow(sig_table))
        c = get(type_colors, row.persistence_type, :gray)
        scatter!(ax, [row.p_short], [i]; color = c, marker = :circle, markersize = 12, label = i==1 ? "Short-run" : "")
        scatter!(ax, [row.p_long],  [i]; color = c, marker = :diamond, markersize = 12, label = i==1 ? "Long-run" : "")
        lines!(ax, [row.p_short, row.p_long], [i, i]; color = (c, 0.5))
    end
    
    vlines!(ax, [0.05]; color = :black, linestyle = :dash, linewidth = 1)
    
    Legend(fig[1, 2],
           [MarkerElement(color = v, marker = :circle, markersize = 12) for v in values(type_colors)],
           collect(keys(type_colors)),
           "Persistence type", framevisible = false)
           
    output_path = joinpath(output_dir, "persistence_summary.pdf")
    save(output_path, fig)
    println("Saved persistence summary to: $output_path")
    return fig
end

# =============================================================================
# 3. MAIN RUNNER (统一调用入口)
# =============================================================================
function run_plots(irfs::Dict{Int, DataFrame}, pw_bh::DataFrame, sig_table::DataFrame, variant::Symbol)
    println("\n=== Generating Plots for $variant ===")
    
    # 动态获取对应的输出文件夹 (该函数应在外部或02文件中已定义)
    output_dir = get_output_dir_ind(variant)
    
    # 由于是无基准组设计，只有绝对效应绘图
    plot_absolute_irf(irfs, pw_bh, variant, output_dir)
    plot_persistence_summary(sig_table, variant, output_dir)
    
    println("All plots generated successfully!\n")
end