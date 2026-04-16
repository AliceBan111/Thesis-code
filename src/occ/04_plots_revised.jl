# =============================================================================
# 04_visualization.jl
# IRF visualization — extremes, turning points, and dual-panel plots
#
# Design alignment with 02_lp_estimation.jl:
#   [1] No baseline group. All 9 groups plotted equally.
#   [2] IRF DataFrames have columns: horizon, beta, ci_lo90, ci_hi90,
#       ci_lo68, ci_hi68.  The point estimate column is "beta" (not "theta").
#   [3] Coefficient names follow the pattern "shock_g{g}".
#   [4] OCC_LABELS covers groups 1-9 (all included).
# =============================================================================

using CSV
using DataFrames
using CairoMakie

# =============================================================================
# PATHS
# =============================================================================
const BASE_DIR    = normpath(joinpath(@__DIR__, "..", "..", "result", "occ"))
const OUTPUT_DIR  = joinpath(BASE_DIR, "analysis")
const EXTREMES_CSV_PATH = joinpath(OUTPUT_DIR, "occupation_irf_extremes_turning_points.csv")

# =============================================================================
# EXTRACT EXTREMES AND TURNING POINTS FROM IRF DATA
# =============================================================================
# For each (group, outcome) combination, computes:
#   - max_abs_shock   : the beta value at the horizon of largest |beta|
#   - max_shock_horizon: the horizon at which |beta| is maximized
#   - turning_point_horizons: horizons where the IRF changes sign
#
# Inputs:
#   irf_dict : Dict{Int, DataFrame}  mapping group id -> IRF DataFrame
#   outcome  : string label for the outcome variable (for labeling)
#
# Returns a DataFrame with one row per group.

function extract_extremes_and_turning_points(irf_dict::Dict{Int, DataFrame},
                                              outcome::String)::DataFrame

    rows = []

    for g in 1:9
        !haskey(irf_dict, g) && continue
        df  = sort(irf_dict[g], :horizon)
        occ = get(OCC_LABELS, g, "group_$g")

        betas = Float64.(df.beta)

        # Peak absolute response
        peak_idx      = argmax(abs.(betas))
        max_abs_shock = betas[peak_idx]
        max_horizon   = df.horizon[peak_idx]

        # Turning points: horizons where the IRF changes sign
        # (identified as consecutive pairs with opposite-sign betas)
        turning_horizons = Int[]
        for i in 2:length(betas)
            if sign(betas[i]) != sign(betas[i-1]) && betas[i] != 0.0 && betas[i-1] != 0.0
                push!(turning_horizons, df.horizon[i])
            end
        end

        tp_str = isempty(turning_horizons) ?
                     "None" :
                     join(string.(turning_horizons), ", ")

        push!(rows, (
            outcome               = outcome,
            occ_group             = g,
            occupation            = occ,
            max_abs_shock         = max_abs_shock,
            max_shock_horizon     = max_horizon,
            turning_point_horizons = tp_str,
        ))
    end

    return DataFrame(rows)
end

# =============================================================================
# DUAL-PANEL VISUALIZATION: EXTREMES AND TURNING POINTS
# =============================================================================
# For each outcome, produces a figure with two side-by-side panels:
#
#   Left panel  : horizontal bar chart of peak beta values
#                 (red = negative response, blue = positive response)
#   Right panel : timeline showing peak horizon (star) and turning points (dot)
#
# Groups are sorted by peak beta magnitude (ascending) so the most affected
# occupation appears at the bottom — a common convention in labor economics.

function visualize_extremes_and_turning_points(csv_path::String,
                                                output_dir::String)

    if !isfile(csv_path)
        error("Data file not found: $csv_path\nPlease run the extremes extraction step first.")
    end

    df = CSV.read(csv_path, DataFrame)
    mkpath(output_dir)

    for sub_df in groupby(df, :outcome)

        outcome_name = first(sub_df.outcome)

        # Sort by peak beta ascending (largest negative at bottom)
        sort!(sub_df, :max_abs_shock)

        occupations = sub_df.occupation
        n           = nrow(sub_df)
        y_positions = 1:n

        fig = Figure(size = (1300, max(400, 60 * n + 100)), fontsize = 13)

        # ── Left panel: bar chart of peak betas ──────────────────────────────
        ax1 = Axis(fig[1, 1],
            title    = "Peak IRF Response — $outcome_name",
            xlabel   = "Peak Beta (shock_g)",
            yticks   = (collect(y_positions), collect(occupations)),
            yreversed = false,
        )

        bar_colors = [v < 0 ? :indianred : :steelblue
                      for v in sub_df.max_abs_shock]

        barplot!(ax1, collect(y_positions), sub_df.max_abs_shock,
                 direction = :x, color = bar_colors)

        vlines!(ax1, 0.0, color = :black, linewidth = 1.5, linestyle = :dash)

        # ── Right panel: timeline of peak horizon and turning points ─────────
        ax2 = Axis(fig[1, 2],
            title             = "Response Timeline (Horizon in Months)",
            xlabel            = "Horizon",
            yticks            = (collect(y_positions), collect(occupations)),
            yticklabelsvisible = false,
            limits            = (-1, H_MAX + 1, nothing, nothing),
            xticks            = 0:6:H_MAX,
        )

        # Horizontal reference lines per occupation
        hlines!(ax2, collect(y_positions), color = (:gray, 0.25), linestyle = :dash)

        for (i, row) in enumerate(eachrow(sub_df))

            # Turning points — orange circles
            if row.turning_point_horizons != "None"
                tp_strs = split(row.turning_point_horizons, ",")
                tps     = parse.(Int, strip.(tp_strs))
                scatter!(ax2, tps, fill(i, length(tps)),
                         color = :orange, markersize = 12, marker = :circle)
            end

            # Peak horizon — red star
            scatter!(ax2, [row.max_shock_horizon], [i],
                     color = :crimson, markersize = 20, marker = :star5)
        end

        # ── Legend ────────────────────────────────────────────────────────────
        legend_elements = [
            MarkerElement(color = :crimson, marker = :star5,  markersize = 16),
            MarkerElement(color = :orange,  marker = :circle, markersize = 12),
        ]
        Legend(fig[1, 3], legend_elements,
               ["Peak Response Horizon", "Sign Reversal (Turning Point)"],
               "Legend")

        linkyaxes!(ax1, ax2)
        colgap!(fig.layout, 15)

        # ── Save ──────────────────────────────────────────────────────────────
        save_png = joinpath(output_dir, "viz_extremes_timeline_$(outcome_name).png")
        save_pdf = joinpath(output_dir, "viz_extremes_timeline_$(outcome_name).pdf")
        save(save_png, fig)
        save(save_pdf, fig)

        println("Saved dual-panel plot for [$outcome_name] -> $save_png")
    end
end

# =============================================================================
# IRF PANEL PLOT — all 9 groups on a single figure
# =============================================================================
# Plots IRFs for all groups in a grid layout (3 columns).
# Shaded bands show 90% and 68% bootstrap confidence intervals.
# Optionally overlays pointwise significance stars from 03_significance_tests.jl.
#
# Inputs:
#   irfs     : Dict{Int, DataFrame} from run_lp / extract_irf
#   variant  : Symbol (:hourly_rate, :income, etc.) for file naming
#   pw_table : (optional) pointwise significance table from build_pointwise_table

function plot_irf_grid(irfs::Dict{Int, DataFrame},
                        variant::Symbol;
                        pw_table::Union{DataFrame, Nothing} = nothing,
                        output_dir::String = get_output_dir(variant))

    n_groups = length(irfs)
    n_cols   = 3
    n_rows   = ceil(Int, n_groups / n_cols)

    fig = Figure(size = (380 * n_cols, 280 * n_rows + 60), fontsize = 12)

    Label(fig[0, :],
          "Impulse Response Functions — $(string(variant)) (all 9 occupation groups)",
          fontsize = 16, font = :bold)

    sorted_groups = sort(collect(keys(irfs)))

    for (idx, g) in enumerate(sorted_groups)
        row_pos = ceil(Int, idx / n_cols)
        col_pos = mod1(idx, n_cols)

        df    = sort(irfs[g], :horizon)
        label = get(OCC_LABELS, g, "group_$g")

        ax = Axis(fig[row_pos, col_pos],
            title   = "($g) $label",
            xlabel  = "Horizon (months)",
            ylabel  = "Beta",
            xticks  = 0:12:H_MAX,
        )

        horizons = Float64.(df.horizon)
        betas    = Float64.(df.beta)
        ci_lo90  = Float64.(df.ci_lo90)
        ci_hi90  = Float64.(df.ci_hi90)
        ci_lo68  = Float64.(df.ci_lo68)
        ci_hi68  = Float64.(df.ci_hi68)

        # 90% CI band (lighter)
        band!(ax, horizons, ci_lo90, ci_hi90,
              color = (:steelblue, 0.20))

        # 68% CI band (darker)
        band!(ax, horizons, ci_lo68, ci_hi68,
              color = (:steelblue, 0.40))

        # Point estimates
        lines!(ax, horizons, betas,
               color = :steelblue, linewidth = 2)

        # Zero reference line
        hlines!(ax, 0.0, color = :black, linewidth = 1.0, linestyle = :dash)

        # Optional: significance stars at BH-rejected horizons
        if !isnothing(pw_table)
            sub_pw = subset(pw_table,
                            :occ_group  => g_col -> g_col .== g,
                            :bh_reject  => r     -> r .== true)
            if nrow(sub_pw) > 0
                sig_h = Float64.(sub_pw.horizon)
                sig_b = Float64.(sub_pw.beta)
                scatter!(ax, sig_h, sig_b .+ 0.005,
                         color = :red, markersize = 8, marker = :star5)
            end
        end
    end

    save_png = joinpath(output_dir, "irf_grid_$(string(variant)).png")
    save_pdf = joinpath(output_dir, "irf_grid_$(string(variant)).pdf")
    save(save_png, fig)
    save(save_pdf, fig)

    println("Saved IRF grid plot -> $save_png")
    return fig
end

# =============================================================================
# MAIN RUNNER
# =============================================================================
# Typical usage from main.jl:
#
#   results, irfs = run_lp(panel, :hourly_rate)
#   sig_table, pw_table = run_significance_tests(irfs, boot_store, coef_names, :hourly_rate)
#
#   # Extract and save extremes CSV (call once across all variants)
#   all_extremes = vcat([extract_extremes_and_turning_points(irfs, string(v))
#                        for (v, irfs) in all_irfs]...)
#   CSV.write(EXTREMES_CSV_PATH, all_extremes)
#
#   # Render plots
#   visualize_extremes_and_turning_points(EXTREMES_CSV_PATH, OUTPUT_DIR)
#   plot_irf_grid(irfs, :hourly_rate; pw_table = pw_table)

function run_visualization(irfs::Dict{Int, DataFrame},
                            variant::Symbol;
                            pw_table::Union{DataFrame, Nothing} = nothing)

    println("\n=== Visualization for variant: $variant ===")

    output_dir = get_output_dir(variant)

    # 1. Compute extremes and turning points for this variant
    extremes_df = extract_extremes_and_turning_points(irfs, string(variant))

    extremes_path = joinpath(output_dir, "irf_extremes_turning_points.csv")
    CSV.write(extremes_path, extremes_df)
    println("Extremes table saved -> $extremes_path")

    # 2. Dual-panel timeline visualization (reads from the just-saved CSV)
    visualize_extremes_and_turning_points(extremes_path, output_dir)

    # 3. IRF grid across all 9 groups
    plot_irf_grid(irfs, variant;
                  pw_table   = pw_table,
                  output_dir = output_dir)

    println("Visualization complete for $variant.")
end
