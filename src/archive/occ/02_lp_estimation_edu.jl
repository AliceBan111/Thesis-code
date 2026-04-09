# =============================================================================
# 02_lp_estimation.jl
# Local Projection estimation
#   - Point estimates: OLS with cell fixed effects (within transform)
#   - Standard errors: Driscoll-Kraay (for reference / diagnostics)
#   - Confidence intervals: block bootstrap (primary inference)
# =============================================================================

using CSV, DataFrames, DataFramesMeta
using LinearAlgebra, Statistics
using Dates, Random
using Printf

include("../data/data_prep_edu.jl")   # H_MAX, L_LAG, EDU_LABELS, OUTPUT_DIR

# Bootstrap settings
const N_BOOT       = 500    # number of bootstrap replications
const BLOCK_SIZE   = 3     # months per block (one full seasonal cycle)
const BOOT_SEED    = 42
const CI_LEVELS    = [0.90, 0.95, 0.68]
const H_MAX = 36
const DK_BW = 4

# =============================================================================
# 1. DRISCOLL-KRAAY STANDARD ERRORS  (reference / diagnostics only)
# =============================================================================
# Follows Hoechle (2007). Handles cross-sectional dependence and serial
# correlation up to bandwidth m. CIs use block bootstrap instead.

function driscoll_kraay_vcov(X::Matrix{Float64}, e::Vector{Float64},
                              t_idx::Vector{Int}; m::Int = 0)::Matrix{Float64}
    N_obs, K = size(X)
    T_vals   = sort(unique(t_idx))
    T        = length(T_vals)
    m        = m == 0 ? floor(Int, T^(1/4)) : m
    t_map    = Dict(v => i for (i, v) in enumerate(T_vals))

    H = zeros(K, T)
    for obs in 1:N_obs
        ti = t_map[t_idx[obs]]
        H[:, ti] .+= X[obs, :] .* e[obs]
    end

    S = zeros(K, K)
    for t in 1:T; S .+= H[:, t] * H[:, t]'; end
    S ./= T

    for l in 1:m
        Γl = zeros(K, K)
        for t in (l+1):T; Γl .+= H[:, t] * H[:, t-l]'; end
        Γl ./= T
        w   = 1.0 - l / (m + 1)
        S .+= w .* (Γl .+ Γl')
    end

    XtX_inv = inv(X' * X)
    return XtX_inv * (T .* S) * XtX_inv
end

# =============================================================================
# 2. WITHIN (DEMEANING) TRANSFORMATION  — fixed effects
# =============================================================================
function within_transform(df::DataFrame,
                           yvars::Vector{Symbol},
                           xvars::Vector{Symbol},
                           id_var::Symbol)::DataFrame
    df_out = copy(df)
    for var in vcat(yvars, xvars)
        gdf   = groupby(df_out, id_var)
        means = combine(gdf, var => mean => Symbol(var, "_mean"))
        df_out = leftjoin(df_out, means, on = id_var)
        df_out[!, var] = df_out[!, var] .- df_out[!, Symbol(var, "_mean")]
        select!(df_out, Not(Symbol(var, "_mean")))
    end
    return df_out
end

# =============================================================================
# 3. BUILD REGRESSION DATA FOR HORIZON h
# =============================================================================
# dep_var_{c,t} = log_rwage_{c,t+h} − log_rwage_{c,t−1}

function build_lp_data(panel::DataFrame, h::Int)::Union{DataFrame, Nothing}
    sort!(panel, [:edu_group, :date])
    gdf    = groupby(panel, :edu_group)
    chunks = DataFrame[]

    for g in gdf
        sub = copy(g)
        n   = nrow(sub)

        lead_wage = Vector{Union{Float64,Missing}}(missing, n)
        h < n && (lead_wage[1:n-h] = sub.log_rincome[h+1:n])

        lag_wage = Vector{Union{Float64,Missing}}(missing, n)
        lag_wage[2:n] = sub.log_rincome[1:n-1]

        sub[!, :dep_var]  = lead_wage .- lag_wage
        sub[!, :lag_wage] = lag_wage
        push!(chunks, sub)
    end

    df_h = vcat(chunks...)
    required = vcat([:dep_var, :shock, :log_oil_lag1, :ffr_lag1],
                    [Symbol("shock_lag", l) for l in 1:L_LAG])
    df_h = dropmissing(df_h, required)
    return nrow(df_h) == 0 ? nothing : df_h
end

# =============================================================================
# 4. BUILD COLUMN LISTS  — interaction dummies added in-place
# =============================================================================
function build_col_lists!(df_h::DataFrame)
    inter_cols = Symbol[]
    lag_cols   = [Symbol("shock_lag", l) for l in 1:L_LAG]
    macro_cols = [:log_oil_lag1, :ffr_lag1]
    cell_cols  = [:log_rwage, :age_mean, :female_share, :married_share, :hours_mean, :unemp_rate]
    x_cols     = vcat([:shock], lag_cols, macro_cols, cell_cols)
    return x_cols, inter_cols
end

# =============================================================================
# 5. OLS ON WITHIN-TRANSFORMED DATA
# =============================================================================
function ols_within(df_h::DataFrame, y_col::Symbol,
                     x_cols::Vector{Symbol}, panel::DataFrame)
    df_w = within_transform(df_h, [y_col], x_cols, :edu_group)
    df_w = dropmissing(df_w, vcat([y_col], x_cols))
    nrow(df_w) == 0 && return nothing

    Y = Float64.(df_w[!, y_col])
    X = hcat(ones(nrow(df_w)), Matrix{Float64}(df_w[!, x_cols]))
    β = (X' * X) \ (X' * Y)
    e = Y .- X * β

    all_dates = sort(unique(panel.date))
    date_map  = Dict(d => i for (i, d) in enumerate(all_dates))
    t_idx     = [date_map[d] for d in df_w.date]

    return (β = β, e = e, X = X, t_idx = t_idx, df_w = df_w)
end

# =============================================================================
# 6. BLOCK BOOTSTRAP  — primary inference
# =============================================================================
# We resample *time blocks* rather than individual observations, preserving:
#   (a) serial correlation in the shock sequence
#   (b) cross-sectional co-movement (all cells share the same shock at time t)
#
# Circular block bootstrap: blocks wrap around at the end of the sample.
# Block size = 12 months (one seasonal cycle) — large enough to capture
# autocorrelation in oil shocks, small enough to keep variance reasonable.
#
# Confidence intervals: percentile method
#   CI_α = [q(α/2, β*), q(1−α/2, β*)]
#
# Bootstrap draws for ALL coefficients are returned so that
# 03_significance_tests.jl can compute joint Wald p-values from the
# bootstrap distribution without any normality assumption.


function block_bootstrap_lp(panel::DataFrame, h::Int;
                              n_boot::Int     = N_BOOT,
                              block_size::Int = BLOCK_SIZE,
                              rng::AbstractRNG = Random.default_rng())
@time begin
    # ── baseline data ───────────────────────────────────────────────────────
    df_h = build_lp_data(panel, h)
    isnothing(df_h) && return nothing

    x_cols, _ = build_col_lists!(df_h)
    y_col     = :dep_var

    df_w = within_transform(df_h, [y_col], x_cols, :edu_group)
    df_w = dropmissing(df_w, vcat([y_col], x_cols))
    nrow(df_w) == 0 && return nothing

    Y_base = Float64.(df_w[!, y_col])
    X_base = hcat(ones(nrow(df_w)), Matrix{Float64}(df_w[!, x_cols]))
    β_base = (X_base' * X_base) \ (X_base' * Y_base)
    e_base = Y_base .- X_base * β_base

    K          = length(β_base)
    coef_names = vcat([:intercept], x_cols)

    # ── time structure ──────────────────────────────────────────────────────
    all_dates  = sort(unique(df_w.date))
    T          = length(all_dates)
    n_blocks   = ceil(Int, T / block_size)

    # Map each date to the row indices belonging to it
    date_to_rows = Dict{Date, Vector{Int}}()
    for (i, row) in enumerate(eachrow(df_h))
        push!(get!(date_to_rows, row.date, Int[]), i)
    end
end
    # ── bootstrap replications ──────────────────────────────────────────────
    boot_matrix = fill(NaN, n_boot, K)
@time begin
    for b in 1:n_boot
        # Circular block bootstrap: draw block starts, wrap around
        starts     = rand(rng, 1:T, n_blocks)
        boot_dates = Date[]
        for s in starts
            for k in 0:(block_size - 1)
                push!(boot_dates, all_dates[mod1(s + k, T)])
            end
        end
        boot_dates = boot_dates[1:T]

        # Collect row indices for this bootstrap sample
        new_rows = Int[]
        for d in boot_dates
            append!(new_rows, get(date_to_rows, d, Int[]))
        end
        isempty(new_rows) && continue

        Y_boot = Y_base[new_rows]
        X_boot = X_base[new_rows, :]

        β_boot = (X_boot' * X_boot) \ (X_boot' * Y_boot)
        any(isnan, β_boot) && continue

        boot_matrix[b, :] = β_boot
    end
end
@time begin

    # Keep only valid replications (no NaN rows)
    valid      = [!any(isnan.(boot_matrix[b, :])) for b in 1:n_boot]
    boot_valid = boot_matrix[valid, :]
    n_valid    = sum(valid)

    n_valid < 50 &&
        @warn "h=$h: only $n_valid valid bootstrap draws; CIs may be unreliable."

    # ── percentile CIs ──────────────────────────────────────────────────────
    ci_lo = Dict{Float64, Vector{Float64}}()
    ci_hi = Dict{Float64, Vector{Float64}}()
    for lvl in CI_LEVELS
        α = 1.0 - lvl
        ci_lo[lvl] = [quantile(boot_valid[:, k], α / 2)       for k in 1:K]
        ci_hi[lvl] = [quantile(boot_valid[:, k], 1.0 - α / 2) for k in 1:K]
    end
end

    all_dates_panel = sort(unique(panel.date))
    date_map        = Dict(d => i for (i, d) in enumerate(all_dates_panel))
    t_idx           = [date_map[d] for d in df_w.date]
    V               = driscoll_kraay_vcov(X_base, e_base, t_idx; m = DK_BW)
    dk_se           = sqrt.(diag(V))

    # ── results DataFrame ───────────────────────────────────────────────────
    res = DataFrame(
        horizon      = h,
        coef_name    = coef_names,
        beta         = β_base,
        dk_se        = dk_se,
        boot_se      = [std(boot_valid[:, k]) for k in 1:K],
        ci_lo95      = ci_lo[0.95],
        ci_hi95      = ci_hi[0.95],
        ci_lo90      = ci_lo[0.90],
        ci_hi90      = ci_hi[0.90],
        n_boot_valid = n_valid,
    )

    return (results = res, boot_draws = boot_valid, coef_names = coef_names)
end

# =============================================================================
# 7. RUN FULL LP  — all horizons h = 0 … H_MAX
# =============================================================================
# Returns:
#   results_df  : stacked coefficient + CI table
#   boot_store  : Dict{h => (n_boot × K) bootstrap draw matrix}
#   coef_names  : coefficient name vector (same for all h)

function run_full_lp(panel::DataFrame;
                      h_max::Int      = H_MAX,
                      n_boot::Int     = N_BOOT,
                      block_size::Int = BLOCK_SIZE,
                      seed::Int       = BOOT_SEED
                      )::Tuple{DataFrame, Dict{Int,Matrix{Float64}}, Vector{Symbol}}

    println("\nRunning LP  h = 0 … $h_max")
    println("Block bootstrap: $n_boot replications, block size = $block_size months")

    rng        = MersenneTwister(seed)
    all_res    = DataFrame[]
    boot_store = Dict{Int, Matrix{Float64}}()
    coef_names_ref = Symbol[]

    for h in 0:h_max
        print("\r  h = $h / $h_max  ")
        @time out = block_bootstrap_lp(panel, h;
                                  n_boot     = n_boot,
                                  block_size = block_size,
                                  rng        = rng)
        isnothing(out) && continue
        push!(all_res, out.results)
        boot_store[h]  = out.boot_draws
        coef_names_ref = out.coef_names
    end

    println("\nLP estimation complete.")
    return vcat(all_res...), boot_store, coef_names_ref
end

# =============================================================================
# 8. EXTRACT IRF PATHS
# =============================================================================
# 分别做LP
#
# Bootstrap CIs for beta_abs use (β_boot + θ_boot) directly,
# correctly accounting for their covariance — no sum-of-SEs approximation.

function extract_irf(results_df::DataFrame,
                      boot_store::Dict{Int, Matrix{Float64}},
                      coef_names::Vector{Symbol})::DataFrame

    horizons = sort(unique(results_df.horizon))

    # Coefficient index positions (intercept = 1)
    β_idx = findfirst(==(:shock), coef_names)
    β_idx === nothing && error(":shock not found in coef_names")
    rows = NamedTuple[]
for h in horizons
        sub = filter(r -> r.horizon == h && r.coef_name == :shock, results_df)
        nrow(sub) == 0 && continue
        β_h = sub.beta[1]

        if haskey(boot_store, h)
            B = boot_store[h]
            push!(rows, (
                horizon  = h,
                beta_abs = β_h,
                boot_se  = std(B[:, β_idx]),
                ci_lo95  = quantile(B[:, β_idx], 0.025),
                ci_hi95  = quantile(B[:, β_idx], 0.975),
                ci_lo90  = quantile(B[:, β_idx], 0.05),
                ci_hi90  = quantile(B[:, β_idx], 0.95),
                ci_lo68  = quantile(B[:, β_idx], 0.16),
                ci_hi68  = quantile(B[:, β_idx], 0.84),
            ))
        end
    end

    return DataFrame(rows)
end

# =============================================================================
# 9. MAIN
# =============================================================================
function run_estimation(panel::DataFrame)
    edu_groups = sort(unique(panel.edu_group))
    
    all_results  = Dict{Int, DataFrame}()
    all_boots    = Dict{Int, Dict{Int, Matrix{Float64}}}()
    all_coefnames = Symbol[]
    
    for g in edu_groups
        println("\n=== Education group $g ===")
        panel_g = filter(r -> r.edu_group == g, panel)
        
        results_g, boot_g, coef_names = run_full_lp(panel_g)
        
        all_results[g]   = results_g
        all_boots[g]     = boot_g
        all_coefnames    = coef_names
        
        CSV.write(joinpath(OUTPUT_DIR, "lp_coefficients_edu$(g).csv"), results_g)
    end
    
    # 每组提取IRF
    all_irfs = Dict{Int, DataFrame}()
    for g in edu_groups
        all_irfs[g] = extract_irf(all_results[g], all_boots[g], all_coefnames)
        CSV.write(joinpath(OUTPUT_DIR, "irf_edu$(g).csv"), all_irfs[g])
    end
    
    return all_results, all_boots, all_coefnames, all_irfs
end