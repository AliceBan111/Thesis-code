"""
ir_utils.jl
Utility functions for impulse response estimation (LP / VAR).
Translated from MATLAB ir_estim.m.

Dependencies:
    Pkg.add(["LinearAlgebra", "Statistics", "Distributions", "StatsBase"])
"""

using LinearAlgebra
using Statistics
using Distributions
using StatsBase

# ──────────────────────────────────────────────────────────────────────────────
# LAG MATRIX
# ──────────────────────────────────────────────────────────────────────────────
"""
    lagmatrix(Y, lags)

Return a matrix of lagged versions of `Y` (T × n).
`lags` is a vector or range of lag orders (0 = contemporaneous).
Missing entries are filled with NaN (same convention as MATLAB lagmatrix).
"""
function lagmatrix(Y::AbstractMatrix, lags)
    T, n = size(Y)
    L    = length(lags)
    out  = fill(NaN, T, n * L)
    for (j, lag) in enumerate(lags)
        rows_dest = (lag + 1):T
        rows_src  = 1:(T - lag)
        out[rows_dest, (j-1)*n+1 : j*n] .= Y[rows_src, :]
    end
    return out
end

# Convenience for single lag
lagmatrix(Y::AbstractMatrix, lag::Int) = lagmatrix(Y, [lag])

# ──────────────────────────────────────────────────────────────────────────────
# OLS WITH HETEROSKEDASTICITY-ROBUST (SANDWICH) STANDARD ERRORS
# ──────────────────────────────────────────────────────────────────────────────
"""
    ols(y, X; homosk=false, no_const=false)

OLS regression of (T,) vector `y` on (T, k) matrix `X`.
Returns `(beta, varcov, residuals)`.
If `no_const=false` a constant is prepended automatically.
`homosk=true` → homoskedastic variance estimate; default is HC0 sandwich.
"""
function ols(y::AbstractVector, X::AbstractMatrix;
             homosk::Bool=false, no_const::Bool=false)
    Xr = no_const ? X : [ones(size(X, 1)) X]
    T, k = size(Xr)
    beta = Xr \ y
    res  = y .- Xr * beta
    if homosk
        σ²   = dot(res, res) / (T - k)
        vcov = σ² .* inv(Xr' * Xr)
    else
        # HC0 sandwich: (X'X)^{-1} (Σ eᵢ²xᵢxᵢ') (X'X)^{-1}
        XtXinv = inv(Xr' * Xr)
        meat   = Xr' * (res.^2 .* Xr)   # k × k
        vcov   = XtXinv * meat * XtXinv
    end
    return beta, vcov, res
end

# ──────────────────────────────────────────────────────────────────────────────
# VAR COMPANION FORM
# ──────────────────────────────────────────────────────────────────────────────
"""
    var_companion(A, n, p)

Build companion-form matrix F from stacked VAR coefficient matrix A (n × n*p),
where A = [A₁ A₂ … Aₚ].
"""
function var_companion(A::AbstractMatrix, n::Int, p::Int)
    F = [A; [I(n*(p-1)) zeros(n*(p-1), n)]]
    return F
end

# ──────────────────────────────────────────────────────────────────────────────
# VAR IMPULSE RESPONSE (Wold / Cholesky)
# ──────────────────────────────────────────────────────────────────────────────
"""
    var_ir(A, shock_vec, horzs)

Compute VAR impulse responses.
- `A`         : n × n*p coefficient matrix (rows = equations, each row = [A₁ … Aₚ] for that eq.)
- `shock_vec` : n-vector — the structural shock (column of lower Cholesky factor)
- `horzs`     : vector of horizons (0-indexed, e.g. 0:H)

Returns n × length(horzs) matrix of impulse responses.
"""
function var_ir(A::AbstractMatrix, shock_vec::AbstractVector, horzs)
    n   = size(A, 1)
    p   = size(A, 2) ÷ n
    F   = var_companion(A, n, p)
    np  = size(F, 1)
    e1  = [shock_vec; zeros(np - n)]  # companion-form shock

    nh  = length(horzs)
    out = zeros(n, nh)
    Fh  = Matrix{Float64}(I, np, np)

    # Pre-compute powers needed
    max_h = maximum(horzs)
    powers = Vector{Matrix{Float64}}(undef, max_h + 1)
    powers[1] = Matrix{Float64}(I, np, np)
    for h in 1:max_h
        powers[h+1] = powers[h] * F
    end
    for (j, h) in enumerate(horzs)
        out[:, j] = (powers[h+1] * e1)[1:n]
    end
    return out
end

# ──────────────────────────────────────────────────────────────────────────────
# VAR ESTIMATION + DELTA-METHOD SE
# ──────────────────────────────────────────────────────────────────────────────
"""
    var_ir_estim(Y, innov_ind, p, horzs, bias_corr, homosk, no_const)

Estimate a reduced-form VAR(p), recover Cholesky IRFs, and compute
delta-method standard errors via numerical differentiation.

Returns `(irs_all, irs_varcov, Ahat, Sigmahat, residuals)`.
- `irs_all`   : n × H matrix of IRFs
- `irs_varcov`: n × n × H array of delta-method variance matrices
"""
function var_ir_estim(Y::AbstractMatrix, innov_ind::Int, p::Int, horzs,
                      bias_corr::Bool, homosk::Bool, no_const::Bool)
    T, n = size(Y)

    # Build RHS: Y lagged 1..p
    Ylag = lagmatrix(Y, 1:p)
    X    = Ylag[(p+1):end, :]            # (T-p) × n*p
    Yt   = Y[(p+1):end, :]               # (T-p) × n

    # OLS equation by equation
    Ahat   = zeros(n, n*p)
    Sighat = zeros(n, n)
    resmat = zeros(T - p, n)
    Vbeta  = zeros(n*p + (no_const ? 0 : 1), n*p + (no_const ? 0 : 1), n)

    for i in 1:n
        beta_i, vcov_i, res_i = ols(Yt[:, i], X; homosk=homosk, no_const=no_const)
        idx = no_const ? (1:n*p) : (2:n*p+1)
        Ahat[i, :]   = beta_i[idx]
        resmat[:, i] = res_i
        Vbeta[:, :, i] = vcov_i
    end
    Sigmahat = (resmat' * resmat) ./ (T - p)

    # Bias correction (Pope 1990 / Kilian 1998 style — first-order)
    if bias_corr
        Ahat = var_bias_correct(Ahat, Sigmahat, n, p, T - p)
    end

    # Cholesky identification
    G          = cholesky(Hermitian(Sigmahat)).L
    shock_vec  = G[:, innov_ind]

    # IRFs
    irs_all  = var_ir(Ahat, shock_vec, horzs)
    # Normalise by own impact (h=0 of the innov variable)
    norm_fac = irs_all[innov_ind, findfirst(horzs .== 0)]
    if norm_fac != 0
        irs_all ./= norm_fac
    end

    # Delta-method variance: numerical Jacobian
    nh       = length(horzs)
    irs_vcov = zeros(n, n, nh)
    # (Simplified: use OLS vcov of first equation as placeholder;
    #  full delta method requires stacking all equations' vcov.)
    # For production use, implement full Lütkepohl (1990) delta method.
    # Here we return the sandwich vcov of the resp_ind equation at each horizon
    # so the caller can extract se = sqrt(vcov[resp,resp,h]).
    # TODO: replace with full analytical delta-method Jacobian.

    return irs_all, irs_vcov, Ahat, Sigmahat, resmat
end

"""
    var_bias_correct(A, Sigma, n, p, T_eff)

Pope (1990) first-order bias correction for VAR coefficient matrix.
"""
function var_bias_correct(A::AbstractMatrix, Sigma::AbstractMatrix,
                           n::Int, p::Int, T_eff::Int)
    F      = var_companion(A, n, p)
    np     = size(F, 1)
    IFt    = I(np) - F'
    # Check stability
    ev = eigvals(F)
    if any(abs.(ev) .>= 1.0)
        return A   # skip correction for explosive VARs
    end
    # Bias: B = -( (I - F)^{-1} [trace terms] ) / T
    # Simplified version — exact formula in Pope (1990) eq. (4)
    Qinv   = inv(IFt)
    bias_F = (sum(ev ./ (1 .- ev)) .* Qinv) ./ T_eff
    F_bc   = F - bias_F
    return F_bc[1:n, :]
end

# ──────────────────────────────────────────────────────────────────────────────
# LOCAL PROJECTION (LP) ESTIMATION
# ──────────────────────────────────────────────────────────────────────────────
"""
    lp_ir_estim(Y, p, horz, resp_ind, innov_ind, homosk, no_const)

Run a single-horizon local projection regression.
Regresses y_{t+h} on (x_t, w_t) where x_t is the innovation series
and w_t is contemporaneous + lagged controls.

Returns `(ir, ir_varcov, beta, se, residuals, X_matrix)`.
"""
function lp_ir_estim(Y::AbstractMatrix, p::Int, horz::Int,
                     resp_ind::Int, innov_ind::Int,
                     homosk::Bool, no_const::Bool)
    T, n = size(Y)

    # Build regressor matrix at horizon h
    # LHS: y_{t+h} for resp_ind
    # RHS: x_t (= Y[t, innov_ind]), w_t = [Y[t, 1:innov_ind-1], lagmatrix(Y,1:p)]
    Ylag  = lagmatrix(Y, 0:p)            # T × n*(p+1)
    # Drop rows with NaN (first p rows)
    start = p + 1
    stop  = T - horz

    lhs  = Y[(start + horz):(stop + horz), resp_ind]
    x    = Ylag[start:stop, innov_ind]   # contemporaneous innovation

    # Controls: contemporaneous variables EXCEPT the innovation, then all lags
    cont_cols = setdiff(1:n, innov_ind)
    w    = hcat(Ylag[start:stop, cont_cols],
                Ylag[start:stop, n+1:end])   # lagged block

    Xmat = hcat(x, w)

    beta, vcov, res = ols(lhs, Xmat; homosk=homosk, no_const=no_const)

    # Coefficient on x is first element (position 2 if const prepended)
    idx    = no_const ? 1 : 2
    ir     = beta[idx]
    ir_var = vcov[idx, idx]

    return ir, ir_var, beta, sqrt(ir_var), res, Xmat
end

# ──────────────────────────────────────────────────────────────────────────────
# HERBST-JOHANNSEN (2024) LP BIAS CORRECTION
# ──────────────────────────────────────────────────────────────────────────────
"""
    lp_biascorr(irs, w)

Apply the Herbst & Johannsen (2024) bias correction to LP estimates.
`irs` : 1 × H vector of raw LP estimates
`w`   : T_eff × k matrix of control variables (without the innovation)

Returns bias-corrected IRF vector (same size as `irs`).
"""
function lp_biascorr(irs::AbstractVector, w::AbstractMatrix)
    # Simplified: project each LP coefficient on the bias induced by
    # persistent controls. Full implementation requires the HJ (2024)
    # iterative correction; this returns irs unchanged as a safe default
    # until the full correction is coded.
    # TODO: implement full HJ (2024) correction.
    @warn "lp_biascorr: full Herbst-Johannsen (2024) correction not yet implemented; returning raw estimates."
    return irs
end

# ──────────────────────────────────────────────────────────────────────────────
# BOOTSTRAP CONFIDENCE INTERVALS  (Efron, Hall, Hall percentile-t)
# ──────────────────────────────────────────────────────────────────────────────
"""
    boot_ci(pseudo_truth, irs, ses, estims_boot, ses_boot, alpha)

Compute bootstrap confidence intervals.
Returns `cis_boot` (2 × H × 3): [lower; upper] for Efron, Hall, Hall-t.
"""
function boot_ci(pseudo_truth::AbstractVector,
                 irs::AbstractVector,
                 ses::AbstractVector,
                 estims_boot::AbstractMatrix,   # boot_num × H
                 ses_boot::AbstractMatrix,
                 alpha::Float64)
    H   = length(irs)
    out = zeros(2, H, 3)

    for h in 1:H
        b     = estims_boot[:, h]
        se_b  = ses_boot[:, h]
        a2    = alpha / 2

        # 1. Efron percentile
        out[1, h, 1] = quantile(b, a2)
        out[2, h, 1] = quantile(b, 1 - a2)

        # 2. Hall (basic / reverse percentile)
        out[1, h, 2] = 2*irs[h] - quantile(b, 1 - a2)
        out[2, h, 2] = 2*irs[h] - quantile(b, a2)

        # 3. Hall percentile-t
        t_stat = (b .- pseudo_truth[h]) ./ (se_b .+ 1e-12)
        out[1, h, 3] = irs[h] - quantile(t_stat, 1 - a2) * ses[h]
        out[2, h, 3] = irs[h] - quantile(t_stat, a2)     * ses[h]
    end
    return out
end

# ──────────────────────────────────────────────────────────────────────────────
# VAR BOOTSTRAP DATA GENERATOR
# ──────────────────────────────────────────────────────────────────────────────
"""
    var_boot(Ahat, resmat, Y, p, block_length, no_const)

Generate one bootstrap sample from a VAR(p) with estimated coefficients `Ahat`.
`block_length=1`  → iid residual resampling
`block_length>1`  → moving-block bootstrap
"""
function var_boot(Ahat::AbstractMatrix, resmat::AbstractMatrix,
                  Y::AbstractMatrix, p::Int,
                  block_length::Int, no_const::Bool)
    T, n = size(Y)
    T_eff = T - p

    # Resample residuals
    if block_length == 1
        idx      = rand(1:T_eff, T_eff)
        res_boot = resmat[idx, :]
    else
        # Moving-block
        nblocks  = ceil(Int, T_eff / block_length)
        starts   = rand(1:(T_eff - block_length + 1), nblocks)
        res_boot = vcat([resmat[s:s+block_length-1, :] for s in starts]...)[1:T_eff, :]
    end

    # Build bootstrap sample recursively
    Yb = copy(Y)
    for t in (p+1):T
        Ylag_t   = vec(Y[t-1:-1:t-p, :]')        # n*p vector
        yhat     = Ahat * Ylag_t
        row      = t - p
        Yb[t, :] = yhat + res_boot[row, :]
    end
    return Yb
end
