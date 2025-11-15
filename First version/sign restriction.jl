using LinearAlgebra, Statistics, DelimitedFiles, Plots, CSV, Measures, StatsBase, Random, DataFrames, JLD2

# 1. Data Preparation Functions

function prepare_data()
    """
    Load and transform raw data into stationary series for VAR analysis
    Returns: Matrix of stationary variables, original level data, and dates
    """

    comprnfb_data = CSV.File("data/COMPRNFB.csv")
    ophnfb_data = CSV.File("data/OPHNFB.csv")
    pcepilfe_data = CSV.File("data/PCEPILFE.csv")
    fedfunds_data = CSV.File("data/FEDFUNDS.csv")
    unrate_data = CSV.File("data/UNRATE.csv")

    # Extract numerical values
    comprnfb = collect(comprnfb_data.COMPRNFB)
    ophnfb = collect(ophnfb_data.OPHNFB)
    pcepilfe = collect(pcepilfe_data.PCEPILFE)
    fedfunds = collect(fedfunds_data.FEDFUNDS)
    unrate = collect(unrate_data.UNRATE)

    dates = nothing
    dates_available = false

    # Try to extract dates
    try
        dates = collect(ophnfb_data.observation_date)
        dates_available = true
    catch
        dates = nothing
        dates_available = false
        println("Warning: No DATE column found, will use observation numbers")
    end

    # Verify all series have same length
    n_obs = minimum([length(comprnfb), length(ophnfb), length(pcepilfe), length(fedfunds), length(unrate)])
    println("Number of observations: $n_obs")

    # Calculate productivity growth: log difference (quarterly growth rate)
    prod_growth = log.(ophnfb[2:n_obs]) - log.(ophnfb[1:n_obs-1])

    # Calculate real wage growth: log difference  
    rw_growth = log.(comprnfb[2:n_obs]) - log.(comprnfb[1:n_obs-1])

    # Calculate price inflation: year-over-year percentage change
    pi_p = zeros(n_obs - 4)
    for i in 5:n_obs
        pi_p[i-4] = (pcepilfe[i] / pcepilfe[i-4] - 1.0) * 100.0
    end

    # Align all variables to the same time period
    start_index = 5
    end_index = n_obs

    # Extract aligned series
    prod_aligned = prod_growth[start_index-1:end]
    rw_aligned = rw_growth[start_index-1:end]
    pi_aligned = pi_p[start_index-4:end]
    fed_aligned = fedfunds[start_index:end_index]
    unrate_aligned = unrate[start_index:end_index]

    # Also keep level data for markup calculation
    ophnfb_aligned = ophnfb[start_index:end_index]
    comprnfb_aligned = comprnfb[start_index:end_index]

    # Extract aligned dates if available
    if dates_available
        dates_aligned = dates[start_index:end_index]
    else
        dates_aligned = collect(start_index:end_index)
    end

    # Final alignment
    final_length = minimum([length(prod_aligned), length(rw_aligned), length(pi_aligned),
        length(fed_aligned), length(unrate_aligned)])

    Y = hcat(prod_aligned[1:final_length],
        rw_aligned[1:final_length],
        pi_aligned[1:final_length],
        unrate_aligned[1:final_length],
        fed_aligned[1:final_length])

    # Calculate price markup using level data (log difference)
    price_markup = log.(ophnfb_aligned[1:final_length]) - log.(comprnfb_aligned[1:final_length])

    # Final dates
    final_dates = dates_aligned[1:final_length]

    println("Final data matrix size: $(size(Y))")
    println("Price markup series length: $(length(price_markup))")
    println("Variables: [Productivity Growth, Real Wage Growth, Inflation, Unemployment Rate, Fed Funds Rate]")

    return Y, price_markup, final_dates
end

# 2. VAR Estimation Functions  

function olsvarc(y::Matrix{Float64}, p::Int)
    """
    Estimate VAR(p) model with constant term using OLS
    """
    t, q = size(y)

    # Construct design matrix X = [constant, y_{t-1}, ..., y_{t-p}]
    X = ones(t - p, 1)

    for i in 1:p
        X = hcat(X, y[p-i+1:t-i, :])
    end

    Y = y[p+1:end, :]
    A = (X'X) \ (X'Y)
    Uhat = Y - X * A
    Σ = (Uhat'Uhat) / (t - p - q * p - 1)

    return A, Σ, Uhat, X
end

function companion_matrix(A::Matrix{Float64}, q::Int, p::Int)
    """
    Build companion form matrix from VAR coefficients
    """
    A_lags = A[2:end, :]
    F = zeros(q * p, q * p)
    F[1:q, 1:q*p] = transpose(A_lags)

    if p > 1
        F[q+1:end, 1:q*(p-1)] = I(q * (p - 1))
    end

    return F
end

# 3. QR Decomposition Based Sign Restriction

function qr_rotation_matrix(K::Int)
    """
    Generate random orthogonal matrix using QR decomposition
    Following Rubio-Ramirez, Waggoner, and Zha (2010)
    
    Input:
      K - dimension of matrix
    Output:
      Q - K×K orthogonal matrix from uniform distribution over O(K)
    """
    # Draw random matrix from N(0,1)
    W = randn(K, K)

    # QR decomposition
    Q, R = qr(W)
    Q = Matrix(Q)  # Convert to regular matrix

    # Normalize: ensure diagonal of R is positive
    for i in 1:K
        if R[i, i] < 0
            Q[:, i] = -Q[:, i]
        end
    end

    return Q
end

function check_sign_restriction(IRF_impact::Matrix{Float64}, price_markup_impact::Vector{Float64},
    inflation_idx::Int, shock_idx::Int)
    """
    Check if sign restriction is satisfied:
    Price markup shock generates positive comovement between markup and inflation
    
    Input:
      IRF_impact - q×q matrix of impact responses
      price_markup_impact - vector of markup responses (prod - real wage)
      inflation_idx - index of inflation variable
      shock_idx - which shock we're testing
    Output:
      true if restriction satisfied, false otherwise
    """
    inflation_response = IRF_impact[inflation_idx, shock_idx]
    markup_response = price_markup_impact[shock_idx]

    # Both should have same sign (positive comovement)
    return (inflation_response > 0 && markup_response > 0) ||
           (inflation_response < 0 && markup_response < 0)
end

function generate_sign_restriction_draws(A::Matrix{Float64}, Σ::Matrix{Float64}, p::Int, h::Int,
    n_draws::Int=1000, max_tries::Int=100000)
    """
    Generate multiple draws satisfying sign restrictions
    
    Input:
      A - VAR coefficient matrix
      Σ - residual covariance matrix  
      p - number of lags
      h - horizon for IRFs
      n_draws - number of accepted draws to generate
      max_tries - maximum attempts
    Output:
      IRF_draws - q×q×(h+1)×n_draws array
      Q_draws - q×q×n_draws array of rotation matrices
      n_total_draws - total attempts needed
    """
    q = size(Σ, 1)
    inflation_idx = 3

    # Cholesky decomposition
    P = cholesky(Σ).L

    # Build companion matrix
    F = companion_matrix(A, q, p)
    J = hcat(I(q), zeros(q, q * (p - 1)))

    # Storage
    IRF_draws = zeros(Float64, q, q, h + 1, n_draws)
    Q_draws = zeros(Float64, q, q, n_draws)

    accepted = 0
    total_tries = 0

    println("\n" * "="^70)
    println("Generating $n_draws draws satisfying sign restrictions...")
    println("="^70)

    while accepted < n_draws && total_tries < max_tries
        total_tries += 1

        # Generate random orthogonal matrix
        Q = qr_rotation_matrix(q)

        # Impact matrix
        B = P * Q
        IRF_impact = J * J' * B

        # Calculate markup impact response
        markup_impact = IRF_impact[1, :] - IRF_impact[2, :]

        # Check all shocks
        for shock_idx in 1:q
            if check_sign_restriction(IRF_impact, markup_impact, inflation_idx, shock_idx)
                # Reorder so markup shock is first
                Q_reordered = copy(Q)
                if shock_idx != 1
                    Q_reordered[:, [1, shock_idx]] = Q_reordered[:, [shock_idx, 1]]
                end

                # Compute full IRF
                B_final = P * Q_reordered
                Fpow = I(q * p)

                for i in 1:h+1
                    G = (J * Fpow * J') * B_final
                    IRF_draws[:, :, i, accepted+1] = G
                    Fpow = Fpow * F
                end

                Q_draws[:, :, accepted+1] = Q_reordered
                accepted += 1

                # Progress report
                if accepted % 100 == 0
                    println("Accepted: $accepted/$n_draws (tried $total_tries)")
                end

                break
            end
        end
    end

    if accepted < n_draws
        error("Only found $accepted draws after $max_tries attempts")
    end

    println("Successfully generated $n_draws draws (total attempts: $total_tries)")
    println("Acceptance rate: $(round(100*n_draws/total_tries, digits=2))%")

    return IRF_draws, Q_draws, total_tries
end

function compute_irf_statistics(IRF_draws::Array{Float64,4})
    """
    Compute median and confidence bands from IRF draws
    """
    q, _, h_plus_1, n_draws = size(IRF_draws)

    median_IRF = zeros(q, q, h_plus_1)
    lower_IRF = zeros(q, q, h_plus_1)
    upper_IRF = zeros(q, q, h_plus_1)

    for i in 1:q, j in 1:q, t in 1:h_plus_1
        draws = IRF_draws[i, j, t, :]
        median_IRF[i, j, t] = median(draws)
        lower_IRF[i, j, t] = quantile(draws, 0.16)
        upper_IRF[i, j, t] = quantile(draws, 0.84)
    end

    return median_IRF, lower_IRF, upper_IRF
end

function extract_structural_shocks(Uhat::Matrix{Float64}, Q_draws::Array{Float64,3},
    P::Matrix{Float64})
    """
    Extract structural shocks from VAR residuals
    ε = B^(-1) * u, where B = P * Q
    
    Input:
      Uhat - T×q matrix of VAR residuals
      Q_draws - q×q×n_draws rotation matrices
      P - q×q Cholesky factor
    Output:
      shocks_draws - n_draws×T×q array (each draw, time, shock)
    """
    T, q = size(Uhat)
    n_draws = size(Q_draws, 3)

    shocks_draws = zeros(n_draws, T, q)

    for draw in 1:n_draws
        B = P * Q_draws[:, :, draw]
        B_inv = inv(B)

        # ε_t = B^(-1) * u_t for each time period
        for t in 1:T
            shocks_draws[draw, t, :] = B_inv * Uhat[t, :]
        end
    end

    return shocks_draws
end

# 4. Results Saving Functions

struct SVARResults
    # Core: for subsequent analysis
    median_shocks::Matrix{Float64}        # T×q matrix
    shock_dates::Vector                   # T vector of dates
    median_B::Matrix{Float64}             # q×q median impact matrix

    # Distribution: for confidence intervals
    shock_draws::Array{Float64,3}         # n_draws×T×q
    lower_shocks::Matrix{Float64}         # T×q (16th percentile)
    upper_shocks::Matrix{Float64}         # T×q (84th percentile)

    # IRF results: for comparison
    median_IRF::Array{Float64,3}          # q×q×(h+1)
    lower_IRF::Array{Float64,3}           # q×q×(h+1)
    upper_IRF::Array{Float64,3}           # q×q×(h+1)
    IRF_draws::Array{Float64,4}           # q×q×(h+1)×n_draws

    # VAR estimation: for documentation
    A::Matrix{Float64}
    Σ::Matrix{Float64}
    var_order::Int

    # Sample Q matrices: for checking identification
    Q_sample::Array{Float64,3}            # q×q×100 (first 100)

    # Meta information
    variable_names::Vector{String}
    shock_names::Vector{String}
    n_draws::Int
    n_total_tries::Int
end

function save_results(results::SVARResults, output_dir::String="output")
    """
    Save results in both CSV and JLD2 formats
    """
    # Create output directory if it doesn't exist
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    println("\n" * "="^70)
    println("Saving results to $output_dir/")
    println("="^70)

    # 1. Save structural shocks to CSV
    T = size(results.median_shocks, 1)
    q = size(results.median_shocks, 2)

    shock_df = DataFrame(
        date=results.shock_dates,
        markup_shock_median=results.median_shocks[:, 1],
        markup_shock_lower=results.lower_shocks[:, 1],
        markup_shock_upper=results.upper_shocks[:, 1]
    )

    # Add other shocks if needed
    for i in 2:q
        shock_df[!, Symbol("shock$(i)_median")] = results.median_shocks[:, i]
    end

    CSV.write("$output_dir/structural_shocks.csv", shock_df)
    println("✓ Saved structural_shocks.csv")

    # 2. Save median IRF to CSV (long format)
    h = size(results.median_IRF, 3) - 1
    irf_data = []

    for i in 1:q, j in 1:q, t in 0:h
        push!(irf_data, (
            response=results.variable_names[i],
            shock=results.shock_names[j],
            horizon=t,
            median=results.median_IRF[i, j, t+1],
            lower=results.lower_IRF[i, j, t+1],
            upper=results.upper_IRF[i, j, t+1]
        ))
    end

    irf_df = DataFrame(irf_data)
    CSV.write("$output_dir/median_irf.csv", irf_df)
    println("✓ Saved median_irf.csv")

    # 3. Save complete results to JLD2
    @save "$output_dir/full_svar_results.jld2" results
    println("✓ Saved full_svar_results.jld2")

    println("\nFiles saved:")
    println("  1. structural_shocks.csv - for occupation analysis")
    println("  2. median_irf.csv - aggregate IRFs with confidence bands")
    println("  3. full_svar_results.jld2 - complete results (Julia format)")
end

# 5. Visualization Functions

function plot_irf_with_bands(median_IRF::Array{Float64,3}, lower_IRF::Array{Float64,3},
    upper_IRF::Array{Float64,3}, variable_names::Vector{String},
    shock_name::String="Price Markup Shock", shock_idx::Int=1)
    """
    Plot IRFs with confidence bands for a specific shock
    """
    q = size(median_IRF, 1)
    h = size(median_IRF, 3) - 1

    plots = []

    # Plot each variable's response
    for i in 1:q
        median_resp = median_IRF[i, shock_idx, :]
        lower_resp = lower_IRF[i, shock_idx, :]
        upper_resp = upper_IRF[i, shock_idx, :]

        p = plot(0:h, median_resp,
            linewidth=2.5,
            color=:blue,
            label="Median",
            title="$(variable_names[i])\nResponse to $shock_name",
            titlefontsize=10,
            xlabel="Horizon (Quarters)",
            ylabel="Response",
            grid=true,
            gridstyle=:dash,
            gridalpha=0.3,
            framestyle=:box,
            legend=:topright)

        # Add confidence band
        plot!(0:h, lower_resp,
            fillrange=upper_resp,
            fillalpha=0.2,
            color=:blue,
            label="68% CI",
            linewidth=0)

        plot!(0:h, lower_resp, color=:blue, linewidth=1, linestyle=:dash, label="")
        plot!(0:h, upper_resp, color=:blue, linewidth=1, linestyle=:dash, label="")

        # Zero line
        hline!([0], line=:dash, color=:black, alpha=0.5, label="")

        push!(plots, p)
    end

    # Also plot markup response (prod - real wage)
    markup_median = median_IRF[1, shock_idx, :] - median_IRF[2, shock_idx, :]
    markup_lower = lower_IRF[1, shock_idx, :] - lower_IRF[2, shock_idx, :]
    markup_upper = upper_IRF[1, shock_idx, :] - upper_IRF[2, shock_idx, :]

    p_markup = plot(0:h, markup_median,
        linewidth=2.5,
        color=:red,
        label="Median",
        title="Price Markup\nResponse to $shock_name",
        titlefontsize=10,
        xlabel="Horizon (Quarters)",
        ylabel="Response",
        grid=true,
        gridstyle=:dash,
        gridalpha=0.3,
        framestyle=:box,
        legend=:topright)

    plot!(0:h, markup_lower,
        fillrange=markup_upper,
        fillalpha=0.2,
        color=:red,
        label="68% CI",
        linewidth=0)

    plot!(0:h, markup_lower, color=:red, linewidth=1, linestyle=:dash, label="")
    plot!(0:h, markup_upper, color=:red, linewidth=1, linestyle=:dash, label="")
    hline!([0], line=:dash, color=:black, alpha=0.5, label="")

    push!(plots, p_markup)

    combined_plot = plot(plots...,
        layout=(3, 2),
        size=(1200, 900),
        left_margin=10mm,
        bottom_margin=10mm,
        top_margin=15mm)

    return combined_plot
end

# 6. Main Analysis Script

function main()
    println("="^70)
    println("SVAR Analysis with Sign Restriction (Price Markup Shock)")
    println("="^70)

    Random.seed!(123)  # For reproducibility

    # Set parameters
    p = 4          # VAR lag order
    h = 20         # IRF horizon (5 years for quarterly data)
    n_draws = 1000 # Number of accepted draws

    # Step 1: Prepare data
    println("\n" * "="^70)
    println("Step 1: Data Preparation")
    println("="^70)
    Y, price_markup, dates = prepare_data()

    # Step 2: Estimate VAR model
    println("\n" * "="^70)
    println("Step 2: VAR Estimation")
    println("="^70)
    A, Σ, Uhat, X = olsvarc(Y, p)
    println("VAR estimation completed:")
    println("  Coefficient matrix size: $(size(A))")
    println("  Sigma matrix size: $(size(Σ))")

    # Check VAR stability
    q = size(Σ, 1)
    eig_vals = eigen(companion_matrix(A, q, p)).values
    max_eig = maximum(abs.(eig_vals))
    println("  Maximum eigenvalue modulus: $(round(max_eig, digits=4))")
    if max_eig < 1.0
        println("VAR model is stable")
    else
        println("VAR model may be unstable")
    end

    # Step 3: Generate draws satisfying sign restrictions
    println("\n" * "="^70)
    println("Step 3: Sign Restriction Identification")
    println("="^70)
    IRF_draws, Q_draws, n_total_tries = generate_sign_restriction_draws(A, Σ, p, h, n_draws)

    # Step 4: Compute statistics
    println("\n" * "="^70)
    println("Step 4: Computing IRF Statistics")
    println("="^70)
    median_IRF, lower_IRF, upper_IRF = compute_irf_statistics(IRF_draws)
    println("Computed median and 68% confidence bands")

    # Step 5: Extract structural shocks
    println("\n" * "="^70)
    println("Step 5: Extracting Structural Shocks")
    println("="^70)
    P = cholesky(Σ).L
    shock_draws = extract_structural_shocks(Uhat, Q_draws, Matrix(P))

    # Compute median and percentiles
    T = size(Uhat, 1)
    median_shocks = zeros(T, q)
    lower_shocks = zeros(T, q)
    upper_shocks = zeros(T, q)

    for t in 1:T, i in 1:q
        draws = shock_draws[:, t, i]
        median_shocks[t, i] = median(draws)
        lower_shocks[t, i] = quantile(draws, 0.16)
        upper_shocks[t, i] = quantile(draws, 0.84)
    end

    println("Extracted structural shocks for T=$T periods")

    # Align dates with shocks (accounting for VAR lag)
    shock_dates = dates[p+1:end]

    # Compute median impact matrix
    median_B = P * median(Q_draws, dims=3)[:, :, 1]

    # Step 6: Create results structure
    variable_names = ["Productivity Growth", "Real Wage Growth", "Price Inflation",
        "Unemployment Rate", "Fed Funds Rate"]
    shock_names = ["Markup Shock", "Shock 2", "Shock 3", "Shock 4", "Shock 5"]

    results = SVARResults(
        median_shocks, shock_dates, median_B,
        shock_draws, lower_shocks, upper_shocks,
        median_IRF, lower_IRF, upper_IRF, IRF_draws,
        A, Σ, p,
        Q_draws[:, :, 1:min(100, n_draws)],
        variable_names, shock_names,
        n_draws, n_total_tries
    )

    # Step 7: Save results
    save_results(results)

    # Step 8: Plot IRFs
    println("\n" * "="^70)
    println("Step 6: Generating Plots")
    println("="^70)
    combined_plot = plot_irf_with_bands(median_IRF, lower_IRF, upper_IRF,
        variable_names, "Price Markup Shock", 1)

    savefig(combined_plot, "output/sign_restriction_irfs.png")
    println("IRF plot saved as 'output/sign_restriction_irfs.png'")

    # Verification
    println("\n" * "="^70)
    println("Verification: Impact Responses to Markup Shock")
    println("="^70)
    markup_impact = median_IRF[1, 1, 1] - median_IRF[2, 1, 1]
    inflation_impact = median_IRF[3, 1, 1]
    println("Price Markup (impact):  $(round(markup_impact, digits=6))")
    println("Inflation (impact):     $(round(inflation_impact, digits=6))")
    if (markup_impact > 0 && inflation_impact > 0) || (markup_impact < 0 && inflation_impact < 0)
        println("Sign restriction satisfied!")
    end

    return results
end

# Run the analysis
results = main()