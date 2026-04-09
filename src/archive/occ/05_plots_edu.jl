# =============================================================================
# 05_plots_edu.jl
# IRF plots for education groups
#   Fig 1: Absolute IRF for each group (1×3 panel)
#   Fig 2: Pairwise difference IRFs (1×3 panel, all pairs)
# Significance markers:
#   ★  = group's own β_h significantly ≠ 0 (from pw, p < 0.05)
#   △  = pairwise difference significantly ≠ 0 (from pd, p < 0.05)
# =============================================================================

include("03_significance_tests_edu.jl")

using CairoMakie

const EDU_COLORS = Dict(
    1 => :steelblue,
    2 => :tomato,
    3 => :seagreen,
)

const PAIRS = [(1,2), (1,3), (2,3)]

# =============================================================================
# 1. ABSOLUTE IRF — one subplot per group
# =============================================================================
function plot_irf_edu(all_irfs::Dict{Int, DataFrame},
                      pw::DataFrame;
                      sig_level::Float64   = 0.05,
                      output_path::String  = joinpath(OUTPUT_DIR, "irf_edu.pdf"))

    fig = Figure(resolution = (1200, 420))

    for (col, g) in enumerate(1:3)
        !haskey(all_irfs, g) && continue
        df    = all_irfs[g]
        label = get(EDU_LABELS, g, "Group $g")
        color = EDU_COLORS[g]

        ax = Axis(fig[1, col];
                  title     = label,
                  xlabel    = "Horizon (months)",
                  ylabel    = col == 1 ? "Log real income change" : "",
                  titlesize = 13,
                  xgridvisible = false,
                  ygridvisible = false)

        # ── CI bands ────────────────────────────────────────────────────────
        band!(ax, df.horizon, df.ci_lo95, df.ci_hi95;
              color = (color, 0.12))
        band!(ax, df.horizon, df.ci_lo90, df.ci_hi90;
              color = (color, 0.25))
        band!(ax, df.horizon, df.ci_lo68, df.ci_hi68;
              color = (color, 0.40))

        # ── Point estimate ───────────────────────────────────────────────────
        lines!(ax, df.horizon, df.beta_abs;
               color = color, linewidth = 2.0)

        # ── Zero line ────────────────────────────────────────────────────────
        hlines!(ax, [0.0]; color = :black, linewidth = 0.8, linestyle = :dash)

        # ── ★ own significance ───────────────────────────────────────────────
        pw_g = filter(r -> r.edu_group == g && r.p_value < sig_level, pw)
        if nrow(pw_g) > 0
            sig_df = innerjoin(select(pw_g, :horizon),
                               select(df, :horizon, :beta_abs),
                               on = :horizon)
            scatter!(ax, sig_df.horizon, sig_df.beta_abs;
                     color = color, marker = :star5, markersize = 10,
                     label = "own sig.")
        end
    end

    # ── Legend note ──────────────────────────────────────────────────────────
    Label(fig[2, 1:3],
          "Shaded bands: 68% / 90% / 95% block-bootstrap CIs  (B=$N_BOOT, block=$BLOCK_SIZE months).  " *
          "★ = own β_h significantly ≠ 0 at $(Int(round(sig_level*100)))% level.",
          fontsize = 9, tellwidth = false)

    rowsize!(fig.layout, 2, Relative(0.08))
    save(output_path, fig)
    println("Saved: $output_path")
    return fig
end

# =============================================================================
# 2. PAIRWISE DIFFERENCE IRF — one subplot per pair
# =============================================================================
function plot_diff_irf_edu(pd::DataFrame;
                            sig_level::Float64  = 0.05,
                            output_path::String = joinpath(OUTPUT_DIR, "irf_diff_edu.pdf"))

    fig = Figure(resolution = (1200, 420))

    for (col, (g1, g2)) in enumerate(PAIRS)
        pd_pair = filter(r -> r.g1 == g1 && r.g2 == g2, pd)
        nrow(pd_pair) == 0 && continue

        l1    = get(EDU_LABELS, g1, "Group $g1")
        l2    = get(EDU_LABELS, g2, "Group $g2")
        title = "$l1 − $l2"
        color = EDU_COLORS[g1]

        ax = Axis(fig[1, col];
                  title     = title,
                  xlabel    = "Horizon (months)",
                  ylabel    = col == 1 ? "Difference in log real income change" : "",
                  titlesize = 12,
                  xgridvisible = false,
                  ygridvisible = false)

        # ── CI bands ────────────────────────────────────────────────────────
        band!(ax, pd_pair.horizon, pd_pair.ci_lo95, pd_pair.ci_hi95;
              color = (color, 0.12))
        band!(ax, pd_pair.horizon, pd_pair.ci_lo90, pd_pair.ci_hi90;
              color = (color, 0.25))

        # ── Difference point estimate ────────────────────────────────────────
        lines!(ax, pd_pair.horizon, pd_pair.diff;
               color = color, linewidth = 2.0)

        # ── Zero line ────────────────────────────────────────────────────────
        hlines!(ax, [0.0]; color = :black, linewidth = 0.8, linestyle = :dash)

        # ── △ pairwise significance ──────────────────────────────────────────
        sig_pair = filter(r -> r.p_value < sig_level, pd_pair)
        if nrow(sig_pair) > 0
            scatter!(ax, sig_pair.horizon, sig_pair.diff;
                     color = color, marker = :utriangle, markersize = 10)
        end
    end

    Label(fig[2, 1:3],
          "Difference = group row − group col.  " *
          "Shaded bands: 90% / 95% block-bootstrap CIs.  " *
          "△ = difference significantly ≠ 0 at $(Int(round(sig_level*100)))% level.",
          fontsize = 9, tellwidth = false)

    rowsize!(fig.layout, 2, Relative(0.08))
    save(output_path, fig)
    println("Saved: $output_path")
    return fig
end

# =============================================================================
# 3. MAIN
# =============================================================================
function run_plots_edu(all_irfs::Dict{Int, DataFrame},
                       pw::DataFrame,
                       pd::DataFrame)

    println("\n=== Generating plots ===")
    plot_irf_edu(all_irfs, pw)
    plot_diff_irf_edu(pd)
    println("All plots saved to: $OUTPUT_DIR")
end