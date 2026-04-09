"""
ir_estim.jl
Julia translation of the MATLAB ir_estim wrapper.

Estimates impulse responses via VAR or Local Projection (LP),
with optional shrinkage (BVAR / Smoothed LP) and bootstrap inference.

Usage:
    include("ir_utils.jl")
    include("ir_estim.jl")

    # Minimal example (LP, delta-method CIs)
    irs, ses, cis = ir_estim(Y, p, 0:12)

    # With bootstrap
    irs, ses, cis, cis_boot, ses_boot =
        ir_estim(Y, p, 0:12; estimator=:lp, bootstrap=:var, boot_num=500)
"""

include("ir_utils.jl")

# ──────────────────────────────────────────────────────────────────────────────
# KEYWORD DEFAULTS STRUCT
# ──────────────────────────────────────────────────────────────────────────────
"""
Options for smoothed LP cross-validation.
"""
Base.@kwdef struct OptsSLP
    lambda_range::Vector{Float64} = vcat(0.001:0.005:0.021,
                                         0.05:0.1:1.05,
                                         2.0:1.0:19.0,
                                         20.0:20.0:100.0,
                                         200.0:200.0:2000.0) |> collect
    cv_folds::Int        = 5
    irf_limit_order::Int = 2    # shrink towards polynomial of this order
    undersmooth::Bool    = false # multiply optimal lambda by 0.1?
end

"""
Options for Bayesian VAR.
"""
Base.@kwdef struct OptsBVAR
    random_walk::Bool = true   # Random walk prior (true) or white noise (false)
    ndraw::Int        = 500    # Posterior draws; 0 = posterior mean only
end

# ──────────────────────────────────────────────────────────────────────────────
# MAIN FUNCTION
# ──────────────────────────────────────────────────────────────────────────────
"""
    ir_estim(Y, p, horzs; kwargs...) -> (irs, ses[, cis[, cis_boot, ses_bootstrap]])

Estimate impulse responses at horizons `horzs` using LP or VAR.

# Arguments
- `Y`     : T × n data matrix
- `p`     : lag length
- `horzs` : vector/range of horizons (e.g. `0:12`)

# Keyword arguments
| Keyword          | Default    | Description                                          |
|------------------|------------|------------------------------------------------------|
| `resp_ind`       | 1          | Index of response variable                          |
| `innov_ind`      | 1          | Index of innovation (Cholesky ordering)              |
| `estimator`      | `:lp`      | `:lp` or `:var`                                     |
| `shrinkage`      | false      | Smoothed LP (if LP) or BVAR (if VAR)                |
| `se_homosk`      | false      | Homoskedastic SEs/bootstrap                         |
| `no_const`       | false      | Omit intercept                                      |
| `alpha`          | 0.05       | Significance level for CIs                          |
| `bias_corr`      | true       | Bias-correct estimates                              |
| `opts_slp`       | OptsSLP()  | Smoothed LP options                                 |
| `opts_bvar`      | OptsBVAR() | BVAR options                                        |
| `bootstrap`      | nothing    | `:var` for VAR bootstrap; `nothing` = delta method  |
| `boot_num`       | 1000       | Bootstrap replications                              |
| `boot_blocklength`| 1         | Block length (1 = iid, >1 = moving block)           |

# Returns
- `irs`           : 1 × H vector of point estimates
- `ses`           : 1 × H vector of standard errors
- `cis`           : 2 × H matrix [lower; upper] (delta method or BVAR credible)
- `cis_boot`      : 2 × H × 3 array (Efron / Hall / Hall-t bootstrap CIs)
- `ses_bootstrap` : 1 × H bootstrap standard deviations

# Example — markup shock on wages (Cholesky, LP, 20 horizons)
```julia
Y     = hcat(log_wages, unemployment, log_income, log_markup)  # T × 4
horzs = 0:20
innov = 4   # markup shock ordered last

irs, ses, cis = ir_estim(Y, 4, horzs;
                          resp_ind  = 1,     # wages
                          innov_ind = innov,
                          estimator = :lp,
                          bias_corr = true)
```
"""
function ir_estim(Y::AbstractMatrix, p::Int, horzs;
                  resp_ind::Int         = 1,
                  innov_ind::Int        = 1,
                  estimator::Symbol     = :lp,
                  shrinkage::Bool       = false,
                  se_homosk::Bool       = false,
                  no_const::Bool        = false,
                  alpha::Float64        = 0.05,
                  bias_corr::Bool       = true,
                  opts_slp::OptsSLP     = OptsSLP(),
                  opts_bvar::OptsBVAR   = OptsBVAR(),
                  bootstrap             = nothing,
                  boot_num::Int         = 1000,
                  boot_blocklength::Int = 1)

    horzs_vec = collect(horzs)
    nh        = length(horzs_vec)
    T, n      = size(Y)

    # Default critical values (normal)
    cvs      = fill(quantile(Normal(), 1 - alpha/2), nh)
    cis_boot = fill(NaN, 2, nh, 3)

    # ── Initialise outputs (may be overwritten below) ─────────────────────────
    irs = zeros(nh)
    ses = zeros(nh)
    cis = zeros(2, nh)

    # ═══════════════════════════════════════════════════════════════════════════
    # POINT ESTIMATES AND VARIANCE
    # ═══════════════════════════════════════════════════════════════════════════

    if estimator == :var && !shrinkage
        # ── VAR, no shrinkage ─────────────────────────────────────────────────
        irs_all, irs_vcov, Ahat, Sigmahat, resmat =
            var_ir_estim(Y, innov_ind, p, horzs_vec, bias_corr, se_homosk, no_const)

        irs = irs_all[resp_ind, :]
        ses = [sqrt(max(irs_vcov[resp_ind, resp_ind, h], 0.0)) for h in 1:nh]

    elseif estimator == :lp && !shrinkage
        # ── LP, no shrinkage ──────────────────────────────────────────────────
        betahat = Vector{Vector{Float64}}(undef, nh)
        resid   = Vector{Vector{Float64}}(undef, nh)
        Xmat    = Vector{Matrix{Float64}}(undef, nh)

        for (hi, h) in enumerate(horzs_vec)
            ir_h, irvar_h, betahat[hi], _, resid[hi], Xmat[hi] =
                lp_ir_estim(Y, p, h, resp_ind, innov_ind, se_homosk, no_const)
            irs[hi] = ir_h
            ses[hi] = sqrt(max(irvar_h, 0.0))
        end

        if bias_corr
            # Controls: contemporaneous (except innovation) + lags 1..p
            Ylag  = lagmatrix(Y, 0:p)
            Ylag  = Ylag[(p+1):end, :]
            wcols = setdiff(1:n, innov_ind)
            w     = hcat(Ylag[:, wcols], Ylag[:, n+1:end])
            irs   = lp_biascorr(irs, w)
        end

    elseif estimator == :var && shrinkage
        # ── BVAR (Minnesota / Litterman) ─────────────────────────────────────
        # NOTE: bvarGLP is not translated here — plug in your preferred
        # Julia BVAR package (e.g. BayesianVAR.jl or a custom implementation).
        error("""BVAR shrinkage not yet implemented in Julia.
              Suggested approach:
                1. Install / write a bvarGLP equivalent.
                2. Draw (beta_draws, sigma_draws) of shape (1+n*p, n, ndraw).
                3. Call var_ir() for each draw.
              See the in-line comments for the full algorithm.""")

    elseif estimator == :lp && shrinkage
        # ── Smoothed LP (leave-one-out CV) ────────────────────────────────────
        # NOTE: locproj_cv and locproj are not translated here.
        # They require a penalised regression (ridge-type) over the full
        # horizon sequence. Placeholder error below.
        error("""Smoothed LP not yet implemented in Julia.
              You need to port locproj_cv() and locproj() from the MATLAB code.
              See OptsSLP for the tuning parameters.""")
    else
        error("Unknown estimator/shrinkage combination.")
    end

    # ═══════════════════════════════════════════════════════════════════════════
    # DELTA-METHOD CONFIDENCE INTERVALS
    # ═══════════════════════════════════════════════════════════════════════════
    # (Skip for BVAR — credible interval is already stored in cis.)
    if !(estimator == :var && shrinkage)
        cis[1, :] = irs .- cvs .* ses
        cis[2, :] = irs .+ cvs .* ses
    end

    # ═══════════════════════════════════════════════════════════════════════════
    # BOOTSTRAP CONFIDENCE INTERVALS
    # ═══════════════════════════════════════════════════════════════════════════
    if !isnothing(bootstrap) && bootstrap == :var
        estims_boot = zeros(boot_num, nh)
        ses_boot    = zeros(boot_num, nh)

        # Pseudo-truth from VAR fit
        irs_var, _, Ahat_var, _, res_var =
            var_ir_estim(Y, innov_ind, p, horzs_vec, bias_corr, se_homosk, no_const)
        pseudo_truth = irs_var[resp_ind, :]

        Threads.@threads for b in 1:boot_num
            Yb = var_boot(Ahat_var, res_var, Y, p, boot_blocklength, no_const)
            irs_b, ses_b, _ = ir_estim(Yb, p, horzs_vec;
                                        resp_ind   = resp_ind,
                                        innov_ind  = innov_ind,
                                        estimator  = estimator,
                                        shrinkage  = shrinkage,
                                        se_homosk  = se_homosk,
                                        no_const   = no_const,
                                        alpha      = alpha,
                                        bias_corr  = bias_corr,
                                        opts_slp   = opts_slp,
                                        opts_bvar  = opts_bvar,
                                        bootstrap  = nothing)  # no nested bootstrap
            estims_boot[b, :] = irs_b
            ses_boot[b, :]    = ses_b
        end

        cis_boot      = boot_ci(pseudo_truth, irs, ses, estims_boot, ses_boot, alpha)
        ses_bootstrap = vec(std(estims_boot, dims=1))
        return irs, ses, cis, cis_boot, ses_bootstrap
    end

    return irs, ses, cis
end
