using LinearAlgebra, Statistics, Plots, DataFrames

"""
    estimate_VAR_SR(df::DataFrame; max_lags::Int=4, shock_col::Int=4, savepath::String="results/")

Estimate a reduced-form VAR and apply sign restrictions.

Arguments:
- df: DataFrame containing columns [:ln_gdp_diff, :pi_p, :u, :markup_growth, :iL]
- max_lags: maximum lag to consider for VAR
- shock_col: column index of the shock of interest for plotting IRFs (default 4, markup_growth)
- savepath: folder to save IRF plots (default "results/")

Returns a dictionary with:
- :VAR => VAR model object
- :VARopt => VAR options used for SR
- :SRout => Sign restriction output
"""
function estimate_VAR_SR(df::DataFrame; max_lags::Int=4, shock_col::Int=1, savepath::String="results/" , method_name::String="method")
    folder_path = joinpath(pwd(), savepath, method_name)
    mkpath(folder_path)
    savefile = joinpath(folder_path, "Combined_IRFs_to_Shock_$(method_name).png")
    
    # 1. Prepare VAR data
    svar_data = df[:, [:ln_gdp_diff, :pi_p, :du, :markup_growth, :iL]]
    X = convert(Matrix{Float64}, coalesce.(Matrix(svar_data), NaN))
    X_clean, fo, lo = CommonSample(X)

    # 2. Determine optimal lag
    optimal_lag, criteria = VARlag(X_clean, max_lags, 1)
    println("Optimal lag: ", optimal_lag)

    # 3. Estimate VAR
    VAR, VARopt = VARmodel(X_clean, optimal_lag, 1)

    # --- Model diagnostics ---
    println("\n--- Model Diagnostics ---")
    println("Maximum eigenvalue: ", VAR[:maxEig])
    if VAR[:maxEig] < 1
        println("VAR system is stable")
    else
        println("Warning: VAR system may be unstable")
    end

    for j in 1:size(X_clean, 2)
        println("Equation $j R²: ", VAR[Symbol("eq$j")][:rsqr])
    end
    println("--------------------------\n")

    # --- Sign restriction setup ---
    SIGN = [
    -1 0 0 0 0;   # ln_gdp_diff  ↓
    1 0 0 0 0;   # pi_p         ↑
    1 0 0 0 0;   # du        ↑     
    1 0 0 0 0;   # markup_growth↑
    0 0 0 0 0    # iL           unrestricted
]

    VARopt[:nsteps] = 20
    VARopt[:ndraws] = 1000
    VARopt[:sr_hor] = 4
    VARopt[:pctg] = 68
    VARopt[:sr_draw] = 100000
    VARopt[:sr_rot] = 1000
    VARopt[:sr_mod] = 0
    VARopt[:mult] = 100

    # Ensure covariance matrix is symmetric
    VAR[:sigma] = 0.5*(VAR[:sigma] + VAR[:sigma]')

    # 4. Run sign restriction identification
    SRout = SR(VAR, SIGN, VARopt)

    # 5. Plot impulse responses for selected shock
    variable_names = ["GDP Growth (ln_gdp_diff)", 
                    "Price Inflation (pi_p)", 
                    "Unemployment Rate Growth (du)", 
                    "Markup Growth (markup_growth)", 
                    "Long-term Interest Rate (iL)"]

    nsteps_actual = size(SRout[:IRmed], 1)

    plt = plot(layout=(3, 2), size=(800, 600))

    for i in 1:5
        plot!(plt[i], 1:nsteps_actual, SRout[:IRinf][:, i, shock_col],
            fillrange=SRout[:IRsup][:, i, shock_col],
            fillalpha=0.3,
            fillcolor=:lightblue,
            linealpha=0,
            label="68% CI")
        
        plot!(plt[i], 1:nsteps_actual, SRout[:IRmed][:, i, shock_col],
            label="Median",
            linewidth=2,
            color=:blue)

        plot!(plt[i], [0, nsteps_actual + 1], [0, 0],
            linestyle=:dash,
            color=:black,
            label="")
        
        title!(plt[i], "$(variable_names[i])")
        xlabel!(plt[i], "Horizon")
    end
    savefig(plt, savefile)
    println("IRF figure saved to: ", savefile)
    display(plt)
    # for i in 1:5
    #     println("$(variable_names[i]): ", round(SRout[:IRmed][1, i, shock_col], digits=4))
    # end

    # Check if sign restrictions are satisfied
    # Construct price markup level by integrating markup_growth (第4列)
    price_markup_response = cumsum(SRout[:IRmed][:, 4, 1])  # markup_growth is 4th variable

    println("Price inflation impact response: ", round(SRout[:IRmed][1, 2, 1], digits=4))  # pi_p is 2nd variable
    println("Price markup impact response: ", round(price_markup_response[1], digits=4))

    if sign(SRout[:IRmed][1, 2, 1]) == sign(price_markup_response[1])
        println("Positive comovement constraint satisfied")
        println("Both variables move in the same direction")
    else
        println("Positive comovement constraint NOT satisfied")
        println("Variables move in opposite directions")
    end

    return Dict(:VAR => VAR, :VARopt => VARopt, :SRout => SRout)
end
