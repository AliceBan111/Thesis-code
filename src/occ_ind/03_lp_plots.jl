# =============================================================================
# 04_plots_occ_ind.jl
# IRF visualisation for occ-within-industry LP results
#
# Produces one PDF per (variant × merged industry):
#   result/occ_ind/{variant}/{ind_label}/irf_occupations.pdf
#
# Layout: 3×3 grid, one panel per occupation group.
# Each panel shows:
#   - 90% CI band (light)
#   - 68% CI band (darker)
#   - Point estimate line
#   - Zero reference line
#
# Entry point:
#   run_visualization(all_results, variant)
#   where all_results = Dict{Int, DataFrame} returned by run_lp()
# =============================================================================

using CairoMakie
using CSV, DataFrames
using Printf

# =============================================================================
# LABELS & COLOURS  (must match 02_lp_estimation_occ_ind.jl)
# =============================================================================
const PLOT_OCC_LABELS = Dict(
    1 => "Managerial",
    2 => "Professional_specialty",
    3 => "High_tech",
    4 => "Sales",
    5 => "Administrative_support",
    6 => "Service",
    7 => "Farming_forestry_constr",
    8 => "Precision_prod_repair",
    9 => "Machine_ops_transport",
)

# const PLOT_IND_LABELS = Dict(
#     1 => "Energy_intensive",
#     2 => "Manufacturing_Construction",
#     3 => "Trade",
#     4 => "Services",
# )

# const PLOT_IND_LABELS = Dict(
#     1  => "Agriculture_forestry_fishing",
#     2  => "Mining",
#     3  => "Construction",
#     4  => "Manufacturing_nondurable",
#     5  => "Manufacturing_durable",
#     6  => "Transportation_utilities",
#     7  => "Wholesale_trade",
#     8  => "Retail_trade",
#     9  => "Finance_insurance_realestate",
#     10 => "Business_repair_services",
#     11 => "Personal_entertainment_services",
#     12 => "Professional_related_services",
# )

const PLOT_IND_LABELS = Dict(
    1 => "Energy_intensive",
    2 => "Manufacturing_Construction",
    3 => "Trade",
    4 => "Services",
)

const PLOT_OCC_COLORS = [
    :steelblue,
    :tomato,
    :seagreen,
    :darkorange,
    :mediumpurple,
    :saddlebrown,
    :hotpink,
    :teal,
    :goldenrod,
]

const VARIANT_YLABELS = Dict(
    :hourly_rate      => "Response of Log Hourly Wage",
    :hours            => "Response of Log Hours Worked",
    :income           => "Response of Log Weekly Income",
    :income_share_var => "Response of Income Share",
    :inequality       => "Response of Log Income Ratio (P75/P25)",
    :median           => "Response of Median Log Income",
    :unemployment     => "Response of Unemployment Rate",
    :employment       => "Response of Log Employment",
)

const H_MAX_PLOT = 36

finite_values(v) = Float64[Float64(x) for x in skipmissing(v) if isfinite(Float64(x))]

# =============================================================================
# OUTPUT DIRECTORY  (mirrors data prep)
# =============================================================================
function get_plot_output_dir(variant::Symbol, ind_group::Int)::String
    ind_label = get(PLOT_IND_LABELS, ind_group, "industry_$(ind_group)")
    dir = joinpath(@__DIR__, "../..", "result", "occ_ind",
                   string(variant), ind_label)
    mkpath(dir)
    return dir
end

# =============================================================================
# GLOBAL Y-LIMITS across all occ groups (for consistent axis scaling)
# =============================================================================
function get_global_ylims(irfs::Dict{Int,DataFrame})
    g_min =  Inf
    g_max = -Inf

    for df in values(irfs)
        nrow(df) == 0 && continue
        for col in (:ci_lo90, :ci_hi90, :beta)
            vals = finite_values(df[!, col])
            isempty(vals) && continue
            g_min = min(g_min, minimum(vals))
            g_max = max(g_max, maximum(vals))
        end
    end

    (!isfinite(g_min) || !isfinite(g_max)) && return (-0.01, 0.01)

    span    = g_max - g_min
    padding = span == 0 ? 0.1 : span * 0.1
    return (g_min - padding, g_max + padding)
end

# =============================================================================
# SINGLE IRF PANEL
# =============================================================================
function plot_single_irf!(ax::Axis, df::DataFrame, occ_g::Int, ylims_all)
    base_color = PLOT_OCC_COLORS[occ_g]

    ylims!(ax, ylims_all[1], ylims_all[2])
    hlines!(ax, [0.0]; color=:black, linewidth=1, linestyle=:dash)

    line_df = dropmissing(df, [:horizon, :beta])
    if nrow(line_df) > 0
        keep_line = isfinite.(Float64.(line_df.horizon)) .& isfinite.(Float64.(line_df.beta))
        line_df = line_df[keep_line, :]
    end

    if nrow(line_df) == 0
        y_mid = (ylims_all[1] + ylims_all[2]) / 2
        text!(ax, [H_MAX_PLOT / 2], [y_mid];
              text=["No data"], align=(:center, :center), color=:gray50)
        return
    end

    band_df = dropmissing(df, [:horizon, :beta, :ci_lo90, :ci_hi90, :ci_lo68, :ci_hi68])
    if nrow(band_df) > 0
        keep_band =
            isfinite.(Float64.(band_df.horizon)) .&
            isfinite.(Float64.(band_df.beta)) .&
            isfinite.(Float64.(band_df.ci_lo90)) .&
            isfinite.(Float64.(band_df.ci_hi90)) .&
            isfinite.(Float64.(band_df.ci_lo68)) .&
            isfinite.(Float64.(band_df.ci_hi68))
        band_df = band_df[keep_band, :]
    end

    if nrow(band_df) > 0
        hs_band = Float64.(band_df.horizon)
        band!(ax, hs_band, Float64.(band_df.ci_lo90), Float64.(band_df.ci_hi90);
              color=(base_color, 0.20))
        band!(ax, hs_band, Float64.(band_df.ci_lo68), Float64.(band_df.ci_hi68);
              color=(base_color, 0.40))
    end

    hs    = Float64.(line_df.horizon)
    betas = Float64.(line_df.beta)
    lines!(ax, hs, betas; color=base_color, linewidth=2.5)
end

# =============================================================================
# 3×3 IRF GRID for one (variant × industry)
# =============================================================================
function plot_irf_grid(irfs::Dict{Int,DataFrame},
                       variant::Symbol,
                       ind_group::Int,
                       output_dir::String)

    fig = Figure(size=(1350, 1050), fontsize=12)

    ind_label  = get(PLOT_IND_LABELS, ind_group, "industry_$(ind_group)")
    ylabel_str = get(VARIANT_YLABELS, variant, "Response of $(string(variant))")
    title_text = "$(ind_label) — $(ylabel_str)"

    Label(fig[0, 1:3], title_text; fontsize=18, font=:bold, tellwidth=false)

    ylims_all  = get_global_ylims(irfs)
    occ_groups = sort(collect(keys(irfs)))

    for (idx, g) in enumerate(occ_groups)
        row_pos = ceil(Int, idx / 3)
        col_pos = mod1(idx, 3)
        label   = get(PLOT_OCC_LABELS, g, "group_$(g)")

        ax = Axis(fig[row_pos, col_pos];
            title   = "($(g)) $(label)",
            xlabel  = "Horizon (months)",
            ylabel  = col_pos == 1 ? ylabel_str : "",
            xticks  = 0:12:H_MAX_PLOT,
        )

        plot_single_irf!(ax, irfs[g], g, ylims_all)
    end

    # Fill empty grid cells if fewer than 9 occ groups
    n_filled = length(occ_groups)
    for idx in (n_filled+1):9
        row_pos = ceil(Int, idx / 3)
        col_pos = mod1(idx, 3)
        Box(fig[row_pos, col_pos]; color=:white, strokecolor=:white)
    end

    output_path = joinpath(output_dir, "irf_occupations.pdf")
    save(output_path, fig)
    println("  Saved: $output_path")
    return output_path
end

# =============================================================================
# EXTRACT IRFs from coefficients DataFrame
# (mirrors extract_irf in 02_lp_estimation_occ_ind.jl but works from CSV too)
# =============================================================================
function extract_irfs_from_results(results_df::DataFrame)::Dict{Int,DataFrame}
    irfs = Dict{Int,DataFrame}()
    coef_names = unique(results_df.coef_name)
    for cname in coef_names
        m = match(r"^shock_g(\d+)$", cname)
        isnothing(m) && continue
        g   = parse(Int, m.captures[1])
        sub = filter(r -> r.coef_name == cname, results_df)
        irfs[g] = sub[:, [:horizon, :beta, :ci_lo90, :ci_hi90, :ci_lo68, :ci_hi68]]
    end
    return irfs
end

# =============================================================================
# MAIN ENTRY POINT
#
# all_results :: Dict{Int, DataFrame}
#   key   = merged ind_group (1-4)
#   value = full lp_coefficients DataFrame for that industry
# =============================================================================
function run_visualization(all_results::Dict{Int,DataFrame}, variant::Symbol)

    println("="^60)
    println("Plotting IRFs | variant: $variant")
    println("="^60)

    for ind_group in sort(collect(keys(all_results)))
        results_df = all_results[ind_group]
        ind_label  = get(PLOT_IND_LABELS, ind_group, "industry_$(ind_group)")

        @info "Plotting $ind_label..."

        irfs       = extract_irfs_from_results(results_df)
        output_dir = get_plot_output_dir(variant, ind_group)

        plot_irf_grid(irfs, variant, ind_group, output_dir)
    end

    println("All plots saved.")
end
