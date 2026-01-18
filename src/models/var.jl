using Plots, DataFrames, Statistics

"""
    plot_VAR_inputs(df::DataFrame; savepath="results/", tag="method")

Plot time series of variables that enter the VAR.

Variables plotted:
- ln_gdp_diff     : GDP growth
- pi_p            : inflation
- du              : change in unemployment rate
- markup_growth   : markup growth
- iL              : long-term interest rate

The function uses a common sample across variables to ensure consistency
with VAR estimation.

Arguments:
- df        : DataFrame containing the VAR variables
- savepath  : directory to save the plot
- tag       : string tag to distinguish specifications (e.g. "method1", "method2")
"""
function plot_VAR_inputs(df::DataFrame; savepath="results/", tag="method")

    vars = [:ln_gdp_diff, :pi_p, :du, :markup_growth, :iL]
    titles = [
        "GDP growth (Δ log GDP)",
        "Inflation (π)",
        "Change in unemployment rate (Δu)",
        "Markup growth",
        "Long-term interest rate"
    ]

    # --- extract data and enforce common sample ---
    X = convert(Matrix{Float64}, coalesce.(Matrix(df[:, vars]), NaN))
    X_clean, fo, lo = CommonSample(X)

    T = size(X_clean, 1)

    # --- build plot ---
    plt = plot(layout=(length(vars), 1),
               size=(900, 1100),
               link=:x)

    for i in 1:length(vars)
        plot!(
            plt[i],
            1:T,
            X_clean[:, i],
            title=titles[i],
            legend=false,
            linewidth=1.5
        )
        hline!(plt[i], [0.0], linestyle=:dash, color=:black)
    end

    xlabel!(plt[end], "Time")

    # --- save output ---
    folder_path = joinpath(pwd(), savepath, "var")
    mkpath(folder_path)

    savefile = joinpath(folder_path, "VAR_inputs_$(tag).png")
    savefig(plt, savefile)

    println("VAR input series plot saved to: ", savefile)
    display(plt)

    return nothing
end
