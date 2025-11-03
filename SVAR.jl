using LinearAlgebra, Statistics, DelimitedFiles, Plots, CSV, Measures, StatsBase

# 1. Data Preparation Functions

function prepare_data()
    """
    Load and transform raw data into stationary series for VAR analysis
    Returns: Matrix of stationary variables [prod_growth, rw_growth, pi_p, fedfunds]
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
    # All growth rates start from observation 2, inflation from observation 5
    # So the common start point is observation 5
    start_index = 5
    end_index = n_obs

    # Extract aligned series
    prod_aligned = prod_growth[start_index-1:end]  # -1 because growth rates are shifted
    rw_aligned = rw_growth[start_index-1:end]
    pi_aligned = pi_p[start_index-4:end]  # -4 because inflation starts at i=5
    fed_aligned = fedfunds[start_index:end_index]
    unrate_aligned = unrate[start_index:end_index]

    # Final alignment - ensure all have same length
    final_length = minimum([length(prod_aligned), length(rw_aligned), length(pi_aligned),
        length(fed_aligned), length(unrate_aligned)])

    Y = hcat(prod_aligned[1:final_length],
        rw_aligned[1:final_length],
        pi_aligned[1:final_length],
        unrate_aligned[1:final_length],
        fed_aligned[1:final_length])

    println("Final data matrix size: $(size(Y))")
    println("Variables: [Productivity Growth, Real Wage Growth, Inflation, Unemployment Rate, Fed Funds Rate]")

    return Y
end

# 2. VAR Estimation Functions  

function olsvarc(y::Matrix{Float64}, p::Int)
    """
    Estimate VAR(p) model with constant term using OLS
    Input:
      y - T * q matrix of endogenous variables
      p - number of lags
    Output:
      A  - coefficient matrix (includes constant + lag coefficients)
      Σ  - residual covariance matrix 
      Uhat - matrix of residuals
      X    - design matrix
    """
    t, q = size(y)

    # Construct design matrix X = [constant, y_{t-1}, ..., y_{t-p}]
    X = ones(t - p, 1)  # constant term

    for i in 1:p
        X = hcat(X, y[p-i+1:t-i, :])
    end

    # Dependent variable: y_t (excluding first p observations)
    Y = y[p+1:end, :]

    # OLS estimation: A = (X'X)^{-1}X'Y
    A = (X'X) \ (X'Y)

    # Residuals and covariance matrix
    Uhat = Y - X * A
    Σ = (Uhat'Uhat) / (t - p - q * p - 1)  # degrees of freedom adjustment

    return A, Σ, Uhat, X
end

function companion_matrix(A::Matrix{Float64}, q::Int, p::Int)
    """
    Build companion form matrix from VAR coefficients
    Input:
      A - coefficient matrix from VAR (including constant)
      q - number of variables
      p - number of lags
    Output:
      F - companion matrix (qp * qp)
    """
    # Extract lag coefficients (exclude constant)
    A_lags = A[2:end, :]  # remove first row (constant)

    # Build companion matrix
    F = zeros(q * p, q * p)
    F[1:q, 1:q*p] = transpose(A_lags)

    if p > 1
        F[q+1:end, 1:q*(p-1)] = I(q * (p - 1))
    end

    return F
end


# 3. Impulse Response Functions (Cholesky identification)

function irf_cholesky(A::Matrix{Float64}, Σ::Matrix{Float64}, p::Int, h::Int)
    """
    Compute impulse responses using Cholesky identification
    Input:
      A - VAR coefficient matrix
      Σ - residual covariance matrix  
      p - number of lags
      h - horizon for IRFs
    Output:
      IRF - 5*5*(h+1) array of impulse responses
    """
    q = size(Σ, 1)

    # Build companion matrix
    F = companion_matrix(A, q, p)

    # Cholesky decomposition for identification
    P = cholesky(Σ).L  # Lower triangular matrix

    # Selection matrix
    J = hcat(I(q), zeros(q, q * (p - 1)))

    # Storage for IRFs
    IRF = zeros(Float64, q * q, h + 1)  # 2D数组: (q² × 时间)

    # Compute IRFs using companion form
    Fpow = I(q * p)
    for i in 1:h+1
        G = (J * Fpow * J') * P
        IRF[:, i] = vec(G)
        Fpow = Fpow * F
    end

    return IRF
end


# 4. Main Analysis Script

function main()
    println("Starting SVAR Analysis...")

    # Set parameters
    p = 4    # VAR lag order
    h = 20   # IRF horizon (5 years for quarterly data)

    # Step 1: Prepare data
    Y = prepare_data()

    # Step 2: Estimate VAR model
    A, Σ, Uhat, X = olsvarc(Y, p)
    println("VAR estimation completed:")
    println("  Coefficient matrix size: $(size(A))")
    println("  Sigma matrix size: $(size(Σ))")

    # 1. Check VAR stability
    q = size(Σ, 1)
    eig_vals = eigen(companion_matrix(A, q, p)).values
    max_eig = maximum(abs.(eig_vals))
    println("  Maximum eigenvalue modulus: $(round(max_eig, digits=4))")
    if max_eig < 1.0
        println("VAR model is stable")
    else
        println("VAR model may be unstable")
    end

    # 2. Check residual autocorrelation
    println("\nResidual autocorrelation (lag 1):")
    variable_names = ["Productivity Growth", "Real Wage Growth", "Price Inflation", "Unemployment Rate", "Fed Funds Rate"]
    for i in 1:q
        acf_res = autocor(Uhat[:, i])
        lag1_acf = round(acf_res[2], digits=4)
        println("  $(variable_names[i]): $lag1_acf")
    end

    # Step 3: Compute impulse responses
    IRF_flat = irf_cholesky(A, Σ, p, h)

    # Step 4: Basic plotting of IRFs
    variable_names = ["Productivity Growth", "Real Wage Growth", "Price Inflation", "Unemployment Rate", "Fed Funds Rate"]

    q = size(Σ, 1)
    IRF = reshape(IRF_flat, (q, q, h + 1))

    # Plot responses to a productivity shock (first variable)

    plots = []
    for i in 1:5
        p_temp = plot(0:h, IRF[i, 1, :],
            linewidth=2.5,
            legend=false,
            title="$(variable_names[i])\nResponse to Productivity Shock",
            titlefontsize=10,
            xlabel="Horizon (Quarters)",
            ylabel="Response",
            grid=true,
            gridstyle=:dash,
            gridalpha=0.3,
            framestyle=:box)

        # Add zero line for reference
        hline!([0], line=:dash, color=:black, alpha=0.7, label="")

        push!(plots, p_temp)
    end

    # Create the combined plot with better layout
    combined_plot = plot(plots...,
        layout=(3, 2),
        size=(1200, 800),
        left_margin=10mm,
        bottom_margin=10mm,
        top_margin=15mm)

    # Save with higher resolution
    savefig(combined_plot, "basic_irfs.png")
    println("Basic IRF plot saved as 'basic_irfs.png'")

    return A, Σ, IRF, Y
end

# Run the analysis
A, Σ, IRF, Y = main()


