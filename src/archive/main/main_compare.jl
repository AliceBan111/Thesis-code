# src/wage_shock_analysis/06_decomposition.jl

# ── Step 1: 运行三套LP，保存结果 ──────────────────────────────

# earnweek
include("src/data/data_prep_occ_income.jl")
include("src/wage_shock_analysis/02_lp_estimation_occ_income.jl")
cd(@__DIR__)
panel_earn = main()
# describe(DataFrame(shock = panel_earn.shock))
# std(panel_earn.shock)
results_earn, boot_earn, coef_earn, irfs_earn = run_estimation(panel_earn)

# hourly rate
include("src/data/data_prep_occ_hourlyrate.jl")
include("src/wage_shock_analysis/02_lp_estimation_occ_hourlyrate.jl")
cd(@__DIR__)
panel_hourly = main()
results_hourly, boot_hourly, coef_hourly, irfs_hourly = run_estimation(panel_hourly)

# hours
include("src/data/data_prep_occ_hours.jl")
include("src/wage_shock_analysis/02_lp_estimation_occ_hours.jl")
cd(@__DIR__)
panel_hours = main()
results_hours, boot_hours, coef_hours, irfs_hours = run_estimation(panel_hours)

# ── Step 2: 提取 beta，合并 ────────────────────────────────────
# 假设你的 results_df 格式是 occ_group, horizon, beta_abs
function extract_beta(irfs::Dict{Int,DataFrame}, g::Int)
    return select(irfs[g], :horizon, :beta_abs)
end

df_decomp = innerjoin(
    rename(extract_beta(irfs_earn,   1), :beta_abs => :beta_earn),
    rename(extract_beta(irfs_hourly, 1), :beta_abs => :beta_wage),
    rename(extract_beta(irfs_hours,  1), :beta_abs => :beta_hours),
    on = :horizon
)

# ── Step 3: 计算分解 ───────────────────────────────────────────
@transform!(df_decomp,
    :wage_share  = :beta_wage  ./ :beta_earn,
    :hours_share = :beta_hours ./ :beta_earn,
    :residual    = :beta_earn .- :beta_wage .- :beta_hours
)

# ── Step 4: 看 h=12 截面 ───────────────────────────────────────
# 对所有9个职业组合并
occ_labels = Dict(
    1 => "Managerial",
    2 => "Professional_specialty",
    3 => "High_tech",
    4 => "Sales",
    5 => "Administrative_support",
    6 => "Service",
    7 => "Farming_forestry_construction",
    8 => "Precision_production_repair",
    9 => "Machine_operators_transport"
)

df_all = DataFrame()

for g in 1:9
    df_g = innerjoin(
        rename(select(irfs_earn[g],   :horizon, :beta_abs), :beta_abs => :beta_earn),
        rename(select(irfs_hourly[g], :horizon, :beta_abs), :beta_abs => :beta_wage),
        rename(select(irfs_hours[g],  :horizon, :beta_abs), :beta_abs => :beta_hours),
        on = :horizon
    )
    @transform!(df_g,
        :occ_group       = occ_labels[g],
        :wage_share      = :beta_wage  ./ :beta_earn,
        :hours_share     = :beta_hours ./ :beta_earn,
        :decomp_residual = :beta_earn .- :beta_wage .- :beta_hours
    )
    append!(df_all, df_g)
end

# 看 h=12 截面
df_h12 = @chain df_all begin
    @subset(:horizon .== 12)
    @select(:occ_group, :beta_earn, :beta_wage, :beta_hours,
            :wage_share, :hours_share, :decomp_residual)
    @orderby(:beta_earn)
end

println(df_h12)

# ── Step 6: 加置信带版本（用boot_earn的CI）─────────────────────
# 这个更严谨，建议用这个替代Step 5
using CairoMakie
const PLOT_DIR = joinpath(@__DIR__, "result", "occ", "analysis")
mkpath(PLOT_DIR)

focus_occs = [
    3 => "High_tech",
    4 => "Sales",
    8 => "Precision_production_repair"
]

fig = Figure(size = (700, 900))

for (i,(g, label)) in enumerate(focus_occs)
    df_earn_g   = irfs_earn[g]
    df_hourly_g = irfs_hourly[g]
    df_hours_g  = irfs_hours[g]

    sort!(df_earn_g,   :horizon)
    sort!(df_hourly_g, :horizon)
    sort!(df_hours_g,  :horizon)

    hs = df_earn_g.horizon

    # fig = Figure(size = (700, 420))
    ax  = Axis(fig[i, 1],
        title  = label,
        xlabel = "Horizon (months)",
        ylabel = "Cumulative response (log points)",
        xticks = 0:6:36,
    )
    ylims!(ax, -0.02, 0.01)

    # earnweek: 黑色 + 68% 带
    band!(ax, hs, df_earn_g.ci_lo68, df_earn_g.ci_hi68,
          color = (:black, 0.12))
    lines!(ax, hs, df_earn_g.beta_abs,
           color = :black, linewidth = 2, label = "Earnweek")

    # hourly wage: 蓝色 + 68% 带
    band!(ax, hs, df_hourly_g.ci_lo68, df_hourly_g.ci_hi68,
          color = (:blue, 0.12))
    lines!(ax, hs, df_hourly_g.beta_abs,
           color = :blue, linewidth = 2, label = "Hourly wage")

    # hours: 红色 + 68% 带
    band!(ax, hs, df_hours_g.ci_lo68, df_hours_g.ci_hi68,
          color = (:red, 0.12))
    lines!(ax, hs, df_hours_g.beta_abs,
           color = :red, linewidth = 2, label = "Hours")

    hlines!(ax, [0.0], color = :gray, linewidth = 1, linestyle = :dash)

    if i == 1
        axislegend(ax, position = :lt, orientation = :horizontal)
    end
end
save(joinpath(PLOT_DIR, "decomp_ci_all_focus_occs.pdf"), fig)


# ── Step 7: h=12 截面 bar chart（所有occupation对比）────────────

df_h12_sorted = sort(df_h12, :beta_earn)
occs  = df_h12_sorted.occ_group
xs    = 1:nrow(df_h12_sorted)
width = 0.25

fig = Figure(size = (900, 480))
ax  = Axis(fig[1, 1],
    title    = "Earnweek decomposition at h=12",
    ylabel   = "Cumulative response (log points)",
    xticks   = (xs, occs),
    xticklabelrotation = π/4,
    xticklabelsize = 11,
)

barplot!(ax, xs .- width, df_h12_sorted.beta_earn,
         width = width, color = :black,  label = "Earnweek")
barplot!(ax, xs,          df_h12_sorted.beta_wage,
         width = width, color = :blue,   label = "Hourly wage")
barplot!(ax, xs .+ width, df_h12_sorted.beta_hours,
         width = width, color = :red,    label = "Hours")

hlines!(ax, [0.0], color = :gray, linewidth = 1, linestyle = :dash)
axislegend(ax, position = :rt)

save(joinpath(PLOT_DIR, "decomp_h12_all_occs.pdf"), fig)
println("h=12 bar chart saved.")


# ── Step 8: 每个occupation一张图，wage_share + hours_share两条线 ──

fig_all = Figure(size = (1400, 900))

for (idx, g) in enumerate(1:9)
    label = occ_labels[g]
    row   = ceil(Int, idx / 3)
    col   = mod1(idx, 3)

    df_g = @subset(df_all, :occ_group .== label)
    sort!(df_g, :horizon)
    hs = df_g.horizon

    wage_share  = clamp.(df_g.beta_wage  ./ df_g.beta_earn, -3.0, 3.0)
    hours_share = clamp.(df_g.beta_hours ./ df_g.beta_earn, -3.0, 3.0)

    ax = Axis(fig_all[row, col],
        title  = label,
        xlabel = row == 3 ? "Horizon (months)" : "",
        ylabel = col == 1 ? "Share of earnweek response" : "",
        xticks = 0:6:36,
        yticks = -3:1.5:3,
        limits = (nothing, (-3.0, 3.0)),
    )

    lines!(ax, hs, wage_share,  color = :blue, linewidth = 2, label = "Wage share")
    lines!(ax, hs, hours_share, color = :red,  linewidth = 2, label = "Hours share")
    hlines!(ax, [0.0], color = :gray, linewidth = 1, linestyle = :dash)

    # 只在第一张图加legend
    idx == 1 && axislegend(ax, position = :rb)
end

save(joinpath(PLOT_DIR, "share_decomp_all_occs.pdf"), fig_all)
println("Share decomposition plot saved.")