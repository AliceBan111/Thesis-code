# =============================================================================
# 03_significance_tests_unified.jl
# Significance testing — fully bootstrap-based (无基准组统一共享版 / 方案A)
# 测试目标: 检验各行业的绝对效应 (beta_abs) 是否显著异于 0
# =============================================================================

using CSV, DataFrames, DataFramesMeta
using LinearAlgebra, Statistics, Distributions
using Printf

# =============================================================================
# 1. POINTWISE BOOTSTRAP p-VALUES
# =============================================================================
# Two-sided p-value for H0: βh_g = 0
# Method: fraction of bootstrap draws that have the opposite sign to β̂,
# multiplied by 2 (sign-flip fraction).

function bootstrap_pvalue_pointwise(beta_hat::Float64,
                                     boot_draws::Vector{Float64})::Float64
    n = length(boot_draws)
    n == 0 && return NaN
    if beta_hat >= 0
        p = 2 * mean(boot_draws .<= 0)
    else
        p = 2 * mean(boot_draws .>= 0)
    end
    return clamp(p, 0.0, 1.0)
end

# =============================================================================
# 2. JOINT WALD BOOTSTRAP p-VALUE ACROSS A HORIZON RANGE
# =============================================================================
# H0: βh_g = 0 for all h in horizon_range
# Test statistic: W = β̂' Σ̂^{-1} β̂

function bootstrap_wald_joint(beta_hat::Vector{Float64},
                               boot_beta::Matrix{Float64})::NamedTuple
    n_boot, H = size(boot_beta)
    H == 0 && return (W = NaN, p_value = NaN, H = 0)

    Σ = cov(boot_beta)
    Σ += 1e-10 * I # Regularise

    Σ_inv = try
        inv(Σ)
    catch
        @warn "Σ singular — using pseudo-inverse for Wald test"
        pinv(Σ)
    end

    W_obs = dot(beta_hat, Σ_inv * beta_hat)

    W_boot = Vector{Float64}(undef, n_boot)
    for b in 1:n_boot
        δ = boot_beta[b, :] .- beta_hat   # recentered at point estimate
        W_boot[b] = dot(δ, Σ_inv * δ)
    end

    p = mean(W_boot .>= W_obs)
    return (W = W_obs, p_value = p, H = H)
end

# =============================================================================
# 3. EXTRACT BOOTSTRAP β DRAWS OVER A HORIZON RANGE
# =============================================================================
function extract_boot_beta(boot_store::Dict{Int, Matrix{Float64}},
                             coef_names::Vector{Symbol},
                             g::Int,
                             horizon_range::AbstractVector{Int})::Matrix{Float64}

    # 动态匹配无基准组对应的 shock_x_ind_g 系数
    col_name = Symbol("shock_x_ind", g)
    ti       = findfirst(==(col_name), coef_names)
    isnothing(ti) && return Matrix{Float64}(undef, 0, 0)

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
function build_significance_table(irfs::Dict{Int, DataFrame},
                                   boot_store::Dict{Int, Matrix{Float64}},
                                   coef_names::Vector{Symbol};
                                   short_range::UnitRange{Int} = 0:11,
                                   long_range::UnitRange{Int}  = 24:36)::DataFrame

    full_range = 0:H_MAX
    rows = []
    group_ids = sort(collect(keys(IND_LABELS)))

    for g in group_ids # 所有组别都需要单独检验
        !haskey(irfs, g) && continue
        df = irfs[g]

        # ── point estimates of beta_abs over each range ─────────────────────
        get_beta = (rng) -> begin
            sub = @subset(df, rng.start .<= :horizon .<= rng.stop)
            nrow(sub) == 0 ? Float64[] : sub.beta_abs
        end

        β_full  = get_beta(full_range.start:full_range.stop)
        β_short = get_beta(short_range)
        β_long  = get_beta(long_range)

        # ── bootstrap draws for each range ──────────────────────────────────
        B_full  = extract_boot_beta(boot_store, coef_names, g, collect(full_range))
        B_short = extract_boot_beta(boot_store, coef_names, g, collect(short_range))
        B_long  = extract_boot_beta(boot_store, coef_names, g, collect(long_range))

        # ── joint Wald p-values ──────────────────────────────────────────────
        w_full  = size(B_full,  2) > 0 ? bootstrap_wald_joint(β_full,  B_full)  : (W = NaN, p_value = NaN, H = 0)
        w_short = size(B_short, 2) > 0 ? bootstrap_wald_joint(β_short, B_short) : (W = NaN, p_value = NaN, H = 0)
        w_long  = size(B_long,  2) > 0 ? bootstrap_wald_joint(β_long,  B_long)  : (W = NaN, p_value = NaN, H = 0)

        # ── peak effect ──────────────────────────────────────────────────────
        if !isempty(β_full)
            full_sub  = @subset(df, :horizon .<= H_MAX)
            peak_idx  = argmax(abs.(full_sub.beta_abs))
            peak_h    = full_sub.horizon[peak_idx]
            peak_β    = full_sub.beta_abs[peak_idx]
            peak_bse  = full_sub.boot_se[peak_idx] # 使用绝对效应的 boot_se
        else
            peak_h = 0; peak_β = NaN; peak_bse = NaN
        end

        push!(rows, (
            ind_group    = g,
            ind_label    = get(IND_LABELS, g, "group_$g"),
            # Full path Wald
            W_full       = w_full.W,
            H_full       = w_full.H,
            p_full       = w_full.p_value,
            # Short-run Wald (h = 0-11)
            W_short      = w_short.W,
            H_short      = w_short.H,
            p_short      = w_short.p_value,
            # Long-run Wald (h = 24-36) — persistence test
            W_long       = w_long.W,
            H_long       = w_long.H,
            p_long       = w_long.p_value,
            # Peak
            peak_horizon = peak_h,
            peak_beta    = peak_β,
            peak_boot_se = peak_bse,
        ))
    end

    return DataFrame(rows)
end

# =============================================================================
# 5. BENJAMINI-HOCHBERG (FDR) CORRECTION
# =============================================================================
function bh_correction(p_values::Vector{Float64};
                         fdr_level::Float64 = 0.10)::Vector{Bool}
    n      = length(p_values)
    ord    = sortperm(p_values)
    rank   = invperm(ord)
    reject = falses(n)
    for i in 1:n
        k = rank[i]
        p_values[i] <= (k / n) * fdr_level && (reject[i] = true)
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
# 7. POINTWISE BOOTSTRAP p-VALUES  — for all (g, h) pairs
# =============================================================================
function build_pointwise_table(irfs::Dict{Int, DataFrame},
                                boot_store::Dict{Int, Matrix{Float64}},
                                coef_names::Vector{Symbol};
                                fdr_level::Float64 = 0.05)::DataFrame

    all_rows = []
    group_ids = sort(collect(keys(IND_LABELS)))

    for g in group_ids
        !haskey(irfs, g) && continue
        col_name = Symbol("shock_x_ind", g)
        ti       = findfirst(==(col_name), coef_names)
        isnothing(ti) && continue

        df = irfs[g]
        for row in eachrow(df)
            h       = row.horizon
            beta_h  = row.beta_abs  # 使用 beta_abs 而不是 theta
            !haskey(boot_store, h) && continue

            boot_β = boot_store[h][:, ti]
            p      = bootstrap_pvalue_pointwise(beta_h, boot_β)

            push!(all_rows, (
                ind_group = g,
                ind_label = get(IND_LABELS, g, "group_$g"),
                horizon   = h,
                beta_abs  = beta_h,
                boot_se   = std(boot_β),
                p_value   = p,
            ))
        end
    end

    pw = DataFrame(all_rows)
    pw[!, :bh_reject] = bh_correction(pw.p_value; fdr_level = fdr_level)
    return pw
end

# =============================================================================
# MAIN RUNNER FOR SIGNIFICANCE TESTS (动态适配版)
# =============================================================================
function run_significance_tests(irfs::Dict{Int, DataFrame}, 
                                boot_store::Dict{Int, Matrix{Float64}}, 
                                coef_names::Vector{Symbol},
                                variant::Symbol)

    println("\n=== Significance Tests (bootstrap-based) for $variant ===")

    # 1. 动态获取与 LP 估计时一致的输出路径
    output_dir = get_output_dir_ind(variant)

    # 2. 运行 Joint Wald tests
    sig_table = build_significance_table(irfs, boot_store, coef_names;
                    short_range = 0:11,
                    long_range  = 24:36)

    # 3. 运行 BH correction
    apply_bh!(sig_table)

    # 4. 运行 Persistence classification
    classify_persistence!(sig_table)

    # 5. 生成 Pointwise table for plots
    pw_bh = build_pointwise_table(irfs, boot_store, coef_names)

    # 6. 保存表格到对应的文件夹
    CSV.write(joinpath(output_dir, "significance_joint.csv"), sig_table)
    CSV.write(joinpath(output_dir, "significance_pointwise.csv"), pw_bh)
    println("Significance tables saved to $output_dir/")

    # ── print summary ────────────────────────────────────────────────────────
    println("\nIndustry Group Summary:")
    println("-"^90)
    for row in eachrow(sig_table)
        @printf("%-35s | Short p=%.3f (%s) | Long p=%.3f (%s) | %-20s\n",
            row.ind_label,
            row.p_short, row.bh_reject_short ? "sig*" : "ns  ",
            row.p_long,  row.bh_reject_long  ? "sig*" : "ns  ",
            row.persistence_type)
    end
    println("-"^90)
    println("* = significant after BH FDR correction (q = 0.05)")
    println("All p-values from block bootstrap (B = $N_BOOT, block = $BLOCK_SIZE months)")

    return sig_table, pw_bh
end