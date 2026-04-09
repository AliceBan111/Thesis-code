# =============================================================================
# 05_plots_occ.jl
# IRF plots for all occupational groups (统一动态版)
# 特性：动态Y轴刻度对齐、动态标题映射、完全解耦、职业专属颜色映射
# =============================================================================

using CSV, DataFrames, DataFramesMeta
using CairoMakie
using Printf

# =============================================================================
# 0. 全局配置与动态映射
# =============================================================================
const OCC_LABELS = Dict(
    1 => "Managerial",
    2 => "Professional_specialty",
    3 => "High_tech",
    4 => "Sales",
    5 => "Administrative_support",
    6 => "Service",
    7 => "Farming_forestry_construction",
    8 => "Precision_production_repair",
    9 => "Machine_operators_transport",
)

# 9 个职业专属颜色
const GROUP_COLORS = [
    :steelblue, :tomato, :seagreen, :darkorange,
    :mediumpurple, :saddlebrown, :hotpink, :teal, :goldenrod
]

# 动态获取当前变量在图表上显示的文本名称
function get_plot_labels(variant::Symbol)
    labels = Dict(
        :hourly_rate  => "Response of Log Hourly Wage",
        :unemployment => "Response of Unemployment Rate (%)",
        :hours        => "Response of Hours Worked",
        :income       => "Response of Log Income",
        :emp          => "Response of Log Employment",
        :inequality   => "Response of Inequality",
        :median       => "Response of Median Wage"
    )
    return get(labels, variant, "Response of $(string(variant))")
end

# 动态计算全局 Y 轴刻度 (以 95% CI 的上下限为准，加 10% 的 padding)
function get_global_ylims(irfs::Dict{Int, DataFrame}; is_relative::Bool = false)
    g_min, g_max = Inf, -Inf
    
    col_min = is_relative ? :ci_lo95_theta : :ci_lo95
    col_max = is_relative ? :ci_hi95_theta : :ci_hi95
    
    groups_to_check = is_relative ? (2:9) : (1:9)
    
    for g in groups_to_check
        !haskey(irfs, g) && continue
        df = irfs[g]
        # 排除可能存在的 NaN
        valid_min = minimum(filter(!isnan, df[!, col_min]))
        valid_max = maximum(filter(!isnan, df[!, col_max]))
        
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
# 1. 绘制绝对 IRF (βh + θh_g)
# =============================================================================
function plot_absolute_irf(irfs::Dict{Int, DataFrame}, pw_bh::DataFrame, variant::Symbol, output_dir::String)
    
    y_label = get_plot_labels(variant)
    ylims_abs = get_global_ylims(irfs, is_relative=false)
    
    fig = Figure(resolution = (1400, 900))
    
    for (idx, g) in enumerate(1:9)
        !haskey(irfs, g) && continue
        df  = irfs[g]
        row = ceil(Int, idx / 3)
        col = mod1(idx, 3)
        
        # 获取当前职业的专属颜色
        base_color = GROUP_COLORS[g]
        
        ax = Axis(fig[row, col],
                  title  = get(OCC_LABELS, g, "Group $g"),
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
        sig_sub = filter(r -> r.occ_group == g && r.bh_reject, pw_bh)
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
# 2. 绘制相对 IRF (θh_g) - 仅限 Group 2-9
# =============================================================================
function plot_relative_irf(irfs::Dict{Int, DataFrame}, variant::Symbol, output_dir::String)
    
    y_label = "Relative to Managers: " * get_plot_labels(variant)
    ylims_rel = get_global_ylims(irfs, is_relative=true)
    
    fig = Figure(resolution = (1400, 900))
    
    for (idx, g) in enumerate(2:9)
        !haskey(irfs, g) && continue
        df  = irfs[g]
        
        # 获取当前职业的专属颜色
        base_color = GROUP_COLORS[g]
        
        plot_idx = idx - 1 
        row = ceil(Int, plot_idx / 3)
        col = mod1(plot_idx, 3)
        
        ax = Axis(fig[row, col],
                  title  = get(OCC_LABELS, g, "Group $g"),
                  xlabel = "Horizon (months)",
                  ylabel = col == 1 ? y_label : "")
        
        ylims!(ax, ylims_rel[1], ylims_rel[2])
        hlines!(ax, [0.0], color = :black, linewidth = 1)
        
        # 使用专属颜色绘制置信区间 (带透明度)
        band!(ax, df.horizon, df.ci_lo90_theta, df.ci_hi90_theta, color = (base_color, 0.3))
        band!(ax, df.horizon, df.ci_lo68_theta, df.ci_hi68_theta, color = (base_color, 0.5))
        
        # 使用专属颜色绘制主线
        lines!(ax, df.horizon, df.theta, color = base_color, linewidth = 2.5)
    end
    
    output_path = joinpath(output_dir, "irf_relative.pdf")
    save(output_path, fig)
    println("Saved relative IRFs to: $output_path")
    return fig
end

# =============================================================================
# 3. 绘制持久性摘要图 (Persistence Classification)
# =============================================================================
function plot_persistence_summary(sig_table::DataFrame, variant::Symbol, output_dir::String)
    
    type_colors = Dict(
        "Persistent" => :darkred,
        "Temporary"  => :steelblue,
        "Delayed"    => :darkorange,
        "None"       => :gray
    )
    
    labels = [get(OCC_LABELS, r.occ_group, "Group $(r.occ_group)") for r in eachrow(sig_table)]
    
    fig = Figure(resolution = (900, 600))
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
# 4. MAIN RUNNER (统一调用入口)
# =============================================================================
function run_plots(irfs::Dict{Int, DataFrame}, pw_bh::DataFrame, sig_table::DataFrame, variant::Symbol)
    println("\n=== Generating Plots for $variant ===")
    
    # 使用你写好的全局路径函数！
    output_dir = get_output_dir(variant)
    
    plot_absolute_irf(irfs, pw_bh, variant, output_dir)
    plot_relative_irf(irfs, variant, output_dir)
    plot_persistence_summary(sig_table, variant, output_dir)
    
    println("All plots generated successfully!\n")
end