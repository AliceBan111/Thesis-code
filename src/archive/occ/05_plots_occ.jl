# =============================================================================
# 05_plots.jl
# IRF plots for all occupational groups
# =============================================================================

using CSV, DataFrames, DataFramesMeta
using CairoMakie
using Printf

include("03_significance_tests_occ.jl")

# Color palette: 9 distinguishable colors
const GROUP_COLORS = [
    :steelblue, :tomato, :seagreen, :darkorange,
    :mediumpurple, :saddlebrown, :hotpink, :teal, :goldenrod
]

# =============================================================================
# 1. PLOT ABSOLUTE IRF (βh + θh_g) FOR ALL GROUPS
# =============================================================================
function plot_absolute_irf(irfs::Dict{Int, DataFrame}, pw_bh::DataFrame;
                             output_path::String = joinpath(OUTPUT_DIR, "irf_absolute.pdf"))

    fig = Figure(resolution = (1400, 900))
    ax_list = []

    for (idx, g) in enumerate(1:9)
        !haskey(irfs, g) && continue
        df  = irfs[g]
        row = ceil(Int, idx / 3)
        col = mod1(idx, 3)
        ax  = Axis(fig[row, col],
                   title  = get(OCC_LABELS, g, "Group $g"),
                   xlabel = "Horizon (months)",
                #    ylabel = "Log wage change",
                   ylabel = "log employment change",
                   titlesize = 11)
        # ylims!(ax, -0.03, 0.02) # income, median
        # ylims!(ax, -0.01, 0.01) # hours
        # CairoMakie.ylims!(ax, -0.02, 0.02) # hourly_rate, emp
        # ylims!(ax, -0.03, 0.03) # inequality
        ylims!(ax, -0.015, 0.015) # unemployment
        # ylims!(ax, -0.005, 0.005) # income share

        # 95% bootstrap CI band (outer, lighter)
        band!(ax, df.horizon, df.ci_lo95, df.ci_hi95;
              color = (GROUP_COLORS[idx], 0.15))
        # 90% bootstrap CI band (inner, darker)
        band!(ax, df.horizon, df.ci_lo90, df.ci_hi90;
              color = (GROUP_COLORS[idx], 0.30))
        # 68% bootstrap CI band (inner, darker)
        band!(ax, df.horizon, df.ci_lo68, df.ci_hi68;
              color = (GROUP_COLORS[idx], 0.50))
        # Point estimates
        lines!(ax, df.horizon, df.beta_abs;
               color = GROUP_COLORS[idx], linewidth = 2)
        # Zero line
        hlines!(ax, [0.0]; color = :black, linewidth = 0.8, linestyle = :dash)

        # Mark BH-significant points for g > 1
        if g > 1
            sig_pts = @subset(pw_bh, :occ_group .== g, :bh_reject .== true)
            if nrow(sig_pts) > 0
                sig_vals = innerjoin(
                            select(sig_pts, :horizon, :occ_group),
                            select(df, :horizon, :beta_abs),
                            on = :horizon)

                CairoMakie.scatter!(ax, sig_vals.horizon, sig_vals.beta_abs;
                         color = GROUP_COLORS[idx], marker = :star5, markersize = 8)
            end
        end

        push!(ax_list, ax)
    end

    # Shared legend note
    Label(fig[4, 1:3],
          "Shaded bands: 68%, 90% and 95% block-bootstrap confidence intervals (B=$N_BOOT, block=$BLOCK_SIZE months). ★ = BH-significant (q=0.05). Baseline: Managerial (Group 1).",
          fontsize = 9, tellwidth = false)

    save(output_path, fig)
    println("Saved: $output_path")
    return fig
end

# =============================================================================
# 2. PLOT DIFFERENTIAL IRF (θh_g) FOR GROUPS 2-9
# =============================================================================
function plot_theta_irf(irfs::Dict{Int, DataFrame}, pw_bh::DataFrame;
                          output_path::String = joinpath(OUTPUT_DIR, "irf_theta.pdf"))

    fig = Figure(resolution = (1400, 750))

    for (idx, g) in enumerate(2:9)
        !haskey(irfs, g) && continue
        df  = irfs[g]
        row = ceil(Int, idx / 4)
        col = mod1(idx, 4)
        ax  = Axis(fig[row, col],
                   title  = get(OCC_LABELS, g, "Group $g"),
                   xlabel = "Horizon (months)",
                   ylabel = "θ (relative to Managerial)",
                   titlesize = 10)
        # ylims!(ax, -0.03, 0.02)  # income, median
        # ylims!(ax, -0.01, 0.01)  # hours
        # CairoMakie.ylims!(ax, -0.02, 0.02)  # hourly_rate, emp
        # ylims!(ax, -0.03, 0.03) # inequality
        ylims!(ax, -0.015, 0.015) # unemployment
        # ylims!(ax, -0.005, 0.005) # income share

        band!(ax, df.horizon, df.ci_lo95_theta, df.ci_hi95_theta;
              color = (GROUP_COLORS[g], 0.15))
        band!(ax, df.horizon, df.ci_lo90_theta, df.ci_hi90_theta;
              color = (GROUP_COLORS[g], 0.30))
        band!(ax, df.horizon, df.ci_lo68_theta, df.ci_hi68_theta;
              color = (GROUP_COLORS[g], 0.50))
        
        lines!(ax, df.horizon, df.theta;
               color = GROUP_COLORS[g], linewidth = 2)
        hlines!(ax, [0.0]; color = :black, linewidth = 0.8, linestyle = :dash)

        # BH-significant points
        sig_pts = @subset(pw_bh, :occ_group .== g, :bh_reject .== true)
        if nrow(sig_pts) > 0
            sig_vals = innerjoin(
                    select(sig_pts, :horizon, :occ_group),
                    select(df, :horizon, :theta),
                    on = :horizon)
            CairoMakie.scatter!(ax, sig_vals.horizon, sig_vals.theta;
                     color = GROUP_COLORS[g], marker = :star5, markersize = 8)
        end
    end

    Label(fig[3, 1:4],
          "θ = differential effect relative to Managerial (Group 1). Bands: 68%/90%/95% block-bootstrap CIs. ★ = BH-significant (q=0.05).",
          fontsize = 9, tellwidth = false)

    save(output_path, fig)
    println("Saved: $output_path")
    return fig
end

# =============================================================================
# 3. SUMMARY PANEL: PERSISTENCE CLASSIFICATION
# =============================================================================
function plot_persistence_summary(sig_table::DataFrame;
                                   output_path::String = joinpath(OUTPUT_DIR, "persistence_summary.pdf"))

    labels  = sig_table.occ_label
    types   = sig_table.persistence_type
    p_short = sig_table.p_short
    p_long  = sig_table.p_long

    type_colors = Dict(
        "Persistent"          => :tomato,
        "Temporary"           => :steelblue,
        "Delayed_persistent"  => :seagreen,
        "Not_significant"     => :lightgray,
    )

    fig = Figure(resolution = (900, 500))
    ax  = Axis(fig[1, 1],
               yticks = (1:length(labels), labels),
               xlabel = "p-value",
               title  = "Persistence Classification by Occupational Group")

    for (i, row) in enumerate(eachrow(sig_table))
        c = get(type_colors, row.persistence_type, :gray)
        scatter!(ax, [row.p_short], [i]; color = c, marker = :circle, markersize = 12,
                 label = "Short-run")
        scatter!(ax, [row.p_long],  [i]; color = c, marker = :diamond, markersize = 12,
                 label = "Long-run")
        lines!(ax, [row.p_short, row.p_long], [i, i]; color = (c, 0.5))
    end

    vlines!(ax, [0.05]; color = :black, linestyle = :dash, linewidth = 1,
            label = "p = 0.05")

    Legend(fig[1, 2],
           [MarkerElement(color = v, marker = :circle, markersize = 12)
            for v in values(type_colors)],
           collect(keys(type_colors)),
           "Persistence type", framevisible = false)

    save(output_path, fig)
    println("Saved: $output_path")
    return fig
end

# =============================================================================
# MAIN
# =============================================================================
function run_plots(irfs, sig_table, pw_bh)
    println("\n=== Generating Plots ===")
    plot_absolute_irf(irfs, pw_bh)
    plot_theta_irf(irfs, pw_bh)
    plot_persistence_summary(sig_table)
    println("All plots saved to: $OUTPUT_DIR")
end