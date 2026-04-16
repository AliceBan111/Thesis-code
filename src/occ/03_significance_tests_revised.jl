# =============================================================================
# 03_significance_tests.jl
# Significance testing — fully bootstrap-based
#
# Design alignment with 02_lp_estimation.jl:
#   [1] No baseline group. All 9 groups estimated equally (shock_g1..shock_g9).
#   [2] Coefficient names follow the pattern "shock_g{g}" (from build_cols!).
#   [3] IRF DataFrames have columns: horizon, beta, ci_lo90, ci_hi90,
#       ci_lo68, ci_hi68  (no "theta" column — the point estimate is "beta").
#   [4] boot_store :: Dict{Int, Matrix{Float64}} maps horizon h → (n_boot × K)
#       bootstrap draw matrix, where K matches the full coef vector from
#       bootstrap_lp (shock interactions + lags + macro controls).
# =============================================================================

using CSV, DataFrames, DataFramesMeta
using LinearAlgebra, Statistics, Distributions
using Printf

# These constants must be available (defined in 02_lp_estimation.jl or main):
#   N_BOOT, BLOCK_SIZE, H_MAX, OCC_LABELS, get_output_dir

# =============================================================================
# 1. POINTWISE BOOTSTRAP p-VALUE
# =============================================================================
# Two-sided p-value for H0: beta_{g,h} = 0.
#
# Method: sign-flip fraction — the proportion of bootstrap draws on the
# opposite side of zero from the point estimate, doubled.
# Equivalent to the bootstrap analogue of a two-sided p-value.
#
# p = 2 * min(P(beta* <= 0 | beta_hat > 0),
#             P(beta* >= 0 | beta_hat < 0))

function bootstrap_pvalue_pointwise(beta_hat::Float64,
                                     boot_draws::Vector{Float64})::Float64
    n = length(boot_draws)
    n == 0 && return NaN
    if beta_hat >= 0
        p = 2.0 * mean(boot_draws .<= 0.0)
    else
        p = 2.0 * mean(boot_draws .>= 0.0)
    end
    return clamp(p, 0.0, 1.0)
end

# =============================================================================
# 2. JOINT WALD BOOTSTRAP p-VALUE ACROSS A HORIZON RANGE
# =============================================================================
# H0: beta_{g,h} = 0  for all h in horizon_range.
#
# Test statistic: W = beta_hat' * Sigma_hat^{-1} * beta_hat
#   where Sigma_hat is the bootstrap covariance of (beta_{h1}, ..., beta_{hH}).
#
# Null distribution: recentered bootstrap draws
#   W*_b = (beta*_b - beta_hat)' * Sigma_hat^{-1} * (beta*_b - beta_hat)
#
# p-value: fraction of W*_b exceeding W.
# This correctly accounts for cross-horizon correlation.

function bootstrap_wald_joint(beta_hat::Vector{Float64},
                               boot_beta::Matrix{Float64})::NamedTuple
    # boot_beta: (n_boot x H) matrix of bootstrap draws for beta at H horizons
    n_boot, H = size(boot_beta)
    H == 0 && return (W = NaN, p_value = NaN, H = 0)

    # Bootstrap covariance of beta across horizons
    Sigma = cov(boot_beta)

    # Ridge regularization to avoid singularity
    Sigma += 1e-10 * I

    Sigma_inv = try
        inv(Sigma)
    catch
        @warn "Sigma singular — using pseudo-inverse for Wald test"
        pinv(Sigma)
    end

    # Observed Wald statistic
    W_obs = dot(beta_hat, Sigma_inv * beta_hat)

    # Recentered bootstrap Wald statistics under the null
    W_boot = Vector{Float64}(undef, n_boot)
    for b in 1:n_boot
        delta     = boot_beta[b, :] .- beta_hat   # recenter at point estimate
        W_boot[b] = dot(delta, Sigma_inv * delta)
    end

    p = mean(W_boot .>= W_obs)
    return (W = W_obs, p_value = p, H = H)
end

# =============================================================================
# 3. EXTRACT BOOTSTRAP BETA DRAWS OVER A HORIZON RANGE
# =============================================================================
# Pulls the column for group g's shock interaction (shock_g{g}) from the
# boot_store for each requested horizon, returning an (n_boot x H) matrix.

function extract_boot_beta(boot_store::Dict{Int, Matrix{Float64}},
                            coef_names::Vector{Symbol},
                            g::Int,
                            horizon_range::AbstractVector{Int})::Matrix{Float64}

    col_name = Symbol("shock_g$(g)")
    ti       = findfirst(==(col_name), coef_names)

    if isnothing(ti)
        @warn "Coefficient $col_name not found in coef_names"
        return Matrix{Float64}(undef, 0, 0)
    end

    # Keep only horizons present in boot_store
    avail_h = [h for h in horizon_range if haskey(boot_store, h)]
    isempty(avail_h) && return Matrix{Float64}(undef, 0, 0)

    n_boot = size(first(values(boot_store)), 1)
    mat    = Matrix{Float64}(undef, n_boot, length(avail_h))

    for (col, h) in enumerate(avail_h)
        mat[:, col] = boot_store[h][:, ti]
    end

    return mat
end

# =============================================================================
# 4. BUILD FULL SIGNIFICANCE TABLE
# =============================================================================
# Runs joint Wald tests for all 9 groups over three horizon ranges:
#   short  = h in 0:11   (impact / short-run response)
#   long   = h in 24:36  (persistence test)
#   full   = h in 0:H_MAX

function build_significance_table(irfs::Dict{Int, DataFrame},
                                   boot_store::Dict{Int, Matrix{Float64}},
                                   coef_names::Vector{Symbol};
                                   short_range::UnitRange{Int} = 0:11,
                                   long_range::UnitRange{Int}  = 24:36)::DataFrame

    full_range = 0:H_MAX
    rows       = []

    for g in 1:9   # All 9 groups — no baseline excluded
        !haskey(irfs, g) && continue
        df = irfs[g]

        # ── point estimates of beta over each range ──────────────────────────
        get_beta = (rng) -> begin
            sub = subset(df, :horizon => h -> rng.start .<= h .<= rng.stop)
            nrow(sub) == 0 ? Float64[] : Float64.(sub.beta)
        end

        beta_full  = get_beta(full_range.start:full_range.stop)
        beta_short = get_beta(short_range)
        beta_long  = get_beta(long_range)

        # ── bootstrap draws for each range ──────────────────────────────────
        B_full  = extract_boot_beta(boot_store, coef_names, g, collect(full_range))
        B_short = extract_boot_beta(boot_store, coef_names, g, collect(short_range))
        B_long  = extract_boot_beta(boot_store, coef_names, g, collect(long_range))

        nan_result = (W = NaN, p_value = NaN, H = 0)

        w_full  = (size(B_full,  1) > 0 && size(B_full,  2) > 0) ?
                      bootstrap_wald_joint(beta_full,  B_full)  : nan_result
        w_short = (size(B_short, 1) > 0 && size(B_short, 2) > 0) ?
                      bootstrap_wald_joint(beta_short, B_short) : nan_result
        w_long  = (size(B_long,  1) > 0 && size(B_long,  2) > 0) ?
                      bootstrap_wald_joint(beta_long,  B_long)  : nan_result

        # ── peak effect (largest |beta| over full horizon) ───────────────────
        if !isempty(beta_full)
            full_sub = subset(df, :horizon => h -> h .<= H_MAX)
            peak_idx = argmax(abs.(Float64.(full_sub.beta)))
            peak_h   = full_sub.horizon[peak_idx]
            peak_b   = full_sub.beta[peak_idx]

            # Bootstrap SE at peak horizon
            col_name = Symbol("shock_g$(g)")
            ti       = findfirst(==(col_name), coef_names)
            if !isnothing(ti) && haskey(boot_store, peak_h)
                peak_boot_se = std(boot_store[peak_h][:, ti])
            else
                peak_boot_se = NaN
            end
        else
            peak_h = 0; peak_b = NaN; peak_boot_se = NaN
        end

        push!(rows, (
            occ_group        = g,
            occ_label        = get(OCC_LABELS, g, "group_$g"),
            # Full path Wald (h = 0:H_MAX)
            W_full           = w_full.W,
            H_full           = w_full.H,
            p_full           = w_full.p_value,
            # Short-run Wald (h = 0:11)
            W_short          = w_short.W,
            H_short          = w_short.H,
            p_short          = w_short.p_value,
            # Long-run Wald (h = 24:36, persistence test)
            W_long           = w_long.W,
            H_long           = w_long.H,
            p_long           = w_long.p_value,
            # Peak response
            peak_horizon     = peak_h,
            peak_beta        = peak_b,
            peak_boot_se     = peak_boot_se,
        ))
    end

    return DataFrame(rows)
end

# =============================================================================
# 5. BENJAMINI-HOCHBERG (FDR) CORRECTION
# =============================================================================
# Standard BH procedure applied to a vector of p-values.
# Returns a BitVector of rejection decisions (true = reject H0).

function bh_correction(p_values::Vector{Float64};
                         fdr_level::Float64 = 0.10)::BitVector
    n      = length(p_values)
    ord    = sortperm(p_values)
    rank   = invperm(ord)
    reject = falses(n)
    for i in 1:n
        k = rank[i]
        if p_values[i] <= (k / n) * fdr_level
            reject[i] = true
        end
    end
    return reject
end

function apply_bh!(sig_table::DataFrame; fdr_level::Float64 = 0.10)::DataFrame
    sig_table[!, :bh_reject_full]  = bh_correction(sig_table.p_full;  fdr_level)
    sig_table[!, :bh_reject_short] = bh_correction(sig_table.p_short; fdr_level)
    sig_table[!, :bh_reject_long]  = bh_correction(sig_table.p_long;  fdr_level)
    return sig_table
end

# =============================================================================
# 6. PERSISTENCE CLASSIFICATION
# =============================================================================
# Classifies each group's response pattern based on BH-corrected rejections:
#
#   Persistent         : significant at both short and long horizons
#   Temporary          : significant short-run only (effect dissipates)
#   Delayed_persistent : significant long-run only (builds over time)
#   Not_significant    : neither short nor long

function classify_persistence!(sig_table::DataFrame)::DataFrame
    sig_table[!, :persistence_type] = map(eachrow(sig_table)) do r
        if r.bh_reject_short && r.bh_reject_long
            "Persistent"
        elseif r.bh_reject_short && !r.bh_reject_long
            "Temporary"
        elseif !r.bh_reject_short && r.bh_reject_long
            "Delayed_persistent"
        else
            "Not_significant"
        end
    end
    return sig_table
end

# =============================================================================
# 7. POINTWISE BOOTSTRAP p-VALUES — for all (g, h) pairs
# =============================================================================
# Used to annotate IRF plots with significance stars (after BH correction).
# Iterates over all 9 groups and all available horizons.

function build_pointwise_table(irfs::Dict{Int, DataFrame},
                                boot_store::Dict{Int, Matrix{Float64}},
                                coef_names::Vector{Symbol};
                                fdr_level::Float64 = 0.05)::DataFrame

    all_rows = []

    for g in 1:9   # All 9 groups
        !haskey(irfs, g) && continue

        col_name = Symbol("shock_g$(g)")
        ti       = findfirst(==(col_name), coef_names)
        isnothing(ti) && continue

        df = irfs[g]

        for row in eachrow(df)
            h       = row.horizon
            beta_h  = Float64(row.beta)
            !haskey(boot_store, h) && continue

            boot_beta = Float64.(boot_store[h][:, ti])
            p         = bootstrap_pvalue_pointwise(beta_h, boot_beta)
            boot_se   = std(boot_beta)

            push!(all_rows, (
                occ_group = g,
                occ_label = get(OCC_LABELS, g, "group_$g"),
                horizon   = h,
                beta      = beta_h,
                boot_se   = boot_se,
                p_value   = p,
            ))
        end
    end

    pw = DataFrame(all_rows)
    pw[!, :bh_reject] = bh_correction(pw.p_value; fdr_level = fdr_level)
    return pw
end

# =============================================================================
# 8. MAIN RUNNER
# =============================================================================
function run_significance_tests(irfs::Dict{Int, DataFrame},
                                 boot_store::Dict{Int, Matrix{Float64}},
                                 coef_names::Vector{Symbol},
                                 variant::Symbol)

    println("\n=== Significance Tests (bootstrap-based) for variant: $variant ===")
    println("Groups: all 9 (no baseline exclusion)")
    println("Bootstrap draws: B = $N_BOOT, block size = $BLOCK_SIZE months")

    output_dir = get_output_dir(variant)

    # Step 1: Joint Wald tests for all groups
    sig_table = build_significance_table(irfs, boot_store, coef_names;
                    short_range = 0:11,
                    long_range  = 24:36)

    # Step 2: BH multiple testing correction
    apply_bh!(sig_table)

    # Step 3: Persistence classification
    classify_persistence!(sig_table)

    # Step 4: Pointwise table for IRF plot annotation
    pw_table = build_pointwise_table(irfs, boot_store, coef_names)

    # Step 5: Save outputs
    CSV.write(joinpath(output_dir, "significance_joint.csv"),     sig_table)
    CSV.write(joinpath(output_dir, "significance_pointwise.csv"), pw_table)
    println("Significance tables saved to: $output_dir")

    # ── Print summary ────────────────────────────────────────────────────────
    println("\nOccupational Group Summary (all 9 groups):")
    println("-"^100)
    @printf("%-5s  %-35s  %-22s  %-22s  %-22s\n",
            "Group", "Label", "Short (h=0-11)", "Long (h=24-36)", "Persistence")
    println("-"^100)

    for row in eachrow(sig_table)
        @printf("%-5d  %-35s  p=%.3f %-6s       p=%.3f %-6s       %-20s\n",
            row.occ_group,
            row.occ_label,
            row.p_short, row.bh_reject_short ? "[sig*]" : "[ns]  ",
            row.p_long,  row.bh_reject_long  ? "[sig*]" : "[ns]  ",
            row.persistence_type)
    end

    println("-"^100)
    println("* = significant after BH FDR correction (q = 0.10)")

    return sig_table, pw_table
end
