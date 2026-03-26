"""
plot_occ_irf_comparison.jl

Plot IRF comparison across 9 occupations for:
- log Income (EARNWEEK)
- log Hours (UHRSWORKT)  
- log Wage (EARNWEEK/UHRSWORKT)

Requires: results dict from run_labor_lp_analysis()
"""

using Plots, ColorSchemes

# 9 occupation colors — distinct enough to tell apart
const OCC_COLORS = [
    colorant"#E63946",   # 1 Managerial          — red
    colorant"#457B9D",   # 2 Professional        — steel blue
    colorant"#2A9D8F",   # 3 High-tech           — teal
    colorant"#E9C46A",   # 4 Sales               — yellow
    colorant"#F4A261",   # 5 Administrative      — orange
    colorant"#6A4C93",   # 6 Service             — purple
    colorant"#4CAF50",   # 7 Farming/construction— green
    colorant"#795548",   # 8 Precision production— brown
    colorant"#607D8B",   # 9 Machine operators   — blue-grey
]

const OCC_SHORT_LABELS = [
    "Managerial",
    "Professional",
    "High-tech",
    "Sales",
    "Admin support",
    "Service",
    "Farming/construct.",
    "Precision prod.",
    "Machine/transport",
]

"""
    plot_irf_comparison(results, outcome; nhorz, save_dir, show_ci)

Plot all 9 occupations' IRFs on one figure for a given outcome.

results  : Dict from run_labor_lp_analysis()
outcome  : :mean_income, :mean_hours, or :mean_wage
show_ci  : show 90% CI bands? (default false — too cluttered with 9 lines)
"""
function plot_irf_comparison(results::Dict,
                              outcome::Symbol;
                              nhorz::Int    = 8,
                              save_dir::String = joinpath(pwd(), "results", "irf_comparison"),
                              show_ci::Bool = false)

    outcome_labels = Dict(
        :mean_income => "log Income (EARNWEEK)",
        :mean_hours  => "log Hours (UHRSWORKT)",
        :mean_wage   => "log Wage (EARNWEEK/UHRSWORKT)",
    )

    mkpath(save_dir)

    pl = plot(size=(900, 500),
              title  = "Markup Shock → $(outcome_labels[outcome])\nby Occupation",
              xlabel = "Horizon (quarters)",
              ylabel = "Response (log points)",
              legend = :outertopright,
              grid   = true,
              framestyle = :box)

    hline!(pl, [0]; color=:black, ls=:dot, lw=1, label="")

    for occ_id in 1:9
        key = (occ_id, outcome)
        !haskey(results, key) && continue

        r    = results[key]
        horz = r.horizons
        coef = r.coef
        lo   = r.lo
        hi   = r.hi
        col  = OCC_COLORS[occ_id]
        lbl  = OCC_SHORT_LABELS[occ_id]

        if show_ci
            plot!(pl, horz, coef;
                  ribbon    = (coef .- lo, hi .- coef),
                  fillalpha = 0.08,
                  fillcolor = col,
                  color     = col,
                  lw        = 2,
                  label     = lbl)
        else
            plot!(pl, horz, coef;
                  color = col,
                  lw    = 2,
                  label = lbl)
        end

        # Mark the point estimate at each horizon
        scatter!(pl, horz, coef;
                 color      = col,
                 markersize = 4,
                 markerstrokewidth = 0,
                 label      = "")
    end

    fname = joinpath(save_dir, "irf_comparison_$(outcome).png")
    savefig(pl, fname)
    println("Saved: $fname")
    display(pl)
    return pl
end

"""
    plot_all_irf_comparisons(results; nhorz, save_dir, show_ci)

Plot comparison figures for all three outcomes and combine into one panel.
"""
function plot_all_irf_comparisons(results::Dict;
                                   nhorz::Int    = 8,
                                   save_dir::String = joinpath(pwd(), "results", "irf_comparison"),
                                   show_ci::Bool = false)

    outcomes = [:mean_income, :mean_hours, :mean_wage]
    plots    = []

    for outcome in outcomes
        pl = plot_irf_comparison(results, outcome;
                                  nhorz    = nhorz,
                                  save_dir = save_dir,
                                  show_ci  = show_ci)
        push!(plots, pl)
    end

    # Combined 3-panel figure
    fig = plot(plots...; layout=(3, 1), size=(900, 1200),
               plot_title="Markup Shock → Labor Outcomes by Occupation")
    fname = joinpath(save_dir, "irf_comparison_all.png")
    savefig(fig, fname)
    println("Combined figure saved: $fname")
    display(fig)

    return plots
end

"""
    plot_irf_heatmap(results, outcome; save_dir)

Heatmap: x = horizon, y = occupation, color = IRF magnitude.
Good for spotting which occupation × horizon combinations drive the Wald test results.
"""
function plot_irf_heatmap(results::Dict,
                           outcome::Symbol;
                           save_dir::String = joinpath(pwd(), "results", "irf_comparison"))

    outcome_labels = Dict(
        :mean_income => "log Income",
        :mean_hours  => "log Hours",
        :mean_wage   => "log Wage",
    )

    mkpath(save_dir)

    nhorz   = 8
    mat     = zeros(9, nhorz)

    for occ_id in 1:9
        key = (occ_id, outcome)
        !haskey(results, key) && continue
        mat[occ_id, :] = results[key].coef[1:nhorz]
    end

    # Symmetric color scale around zero
    lim = maximum(abs.(mat))

    pl = heatmap(1:nhorz, 1:9, mat;
                 color      = :RdBu,
                 clims      = (-lim, lim),
                 xlabel     = "Horizon (quarters)",
                 ylabel     = "Occupation",
                 title      = "IRF Heatmap: Markup Shock → $(outcome_labels[outcome])",
                 yticks     = (1:9, OCC_SHORT_LABELS),
                 size       = (700, 450),
                 colorbar_title = "Response")

    fname = joinpath(save_dir, "irf_heatmap_$(outcome).png")
    savefig(pl, fname)
    println("Saved: $fname")
    display(pl)
    return pl
end
