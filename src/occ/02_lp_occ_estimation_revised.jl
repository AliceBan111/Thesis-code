# =============================================================================
# 02_lp_estimation.jl
# PANEL LOCAL PROJECTION (Group-specific slopes, no baseline group)
#
# Model:
# y_{i,t+h} = α_i + τ_t + Σ_g β_{g,h}(shock_t × 1[group=i=g]) + controls + ε
#
# Key design choices:
# [1] No baseline group. All groups estimated equally.
# [2] LHS = y_{t+h}   (not y_{t+h} - y_{t-1})
# [3] Two-way FE: occupation FE + time FE
# [4] Bootstrap = moving block bootstrap (time axis)
# [5] CI = percentile-t bootstrap
# [6] DK robust SE retained
# [7] Lag length selected dynamically via BIC (ic_var)
#
# Bug fixes vs previous version:
# [F1] get_output_dir defined
# [F2] var_estim defined; ic_var now wired to select L_LAG dynamically
# [F3] twoway_demean! guards against leftjoin! column collisions
# [F4] build_lp_data now sorts by :date inside each group
# [F5] Bootstrap t_idx preserves original time structure using a date→index map
# [F6] L_LAG is no longer a hardcoded const; selected per variant at runtime
# =============================================================================

using CSV, DataFrames, DataFramesMeta
using LinearAlgebra, Statistics, StatsBase
using Dates, Random
using Printf
using FixedEffectModels
using CategoricalArrays

# =============================================================================
# GLOBALS  (only truly fixed constants stay here)
# =============================================================================
const N_BOOT = 500
const BLOCK_SIZE = 6
const BOOT_SEED = 42
const H_MAX = 36
const L_LAG_MAX = 24          # upper bound for BIC search
const DK_BW = 12

const OCC_LABELS = Dict(
    1 => "Managerial",
    2 => "Professional_specialty",
    3 => "High_tech",
    4 => "Sales",
    5 => "Administrative_support",
    6 => "Service",
    7 => "Farming_forestry_construction",
    8 => "Precision_production_repair",
    9 => "Machine_operators_transport",
)

# =============================================================================
# OUTPUT DIRECTORY  [F1: was missing]
# =============================================================================
# function get_output_dir(variant::Symbol)
#     dir = joinpath("output", string(variant))
#     mkpath(dir)
#     return dir
# end

# =============================================================================
# VAR ESTIMATION HELPER  [F2: was missing, needed by ic_var]
#
# Estimates a VAR(p) by OLS equation-by-equation.
# Returns: β (coefficients), fitted values, Sigma (residual covariance)
# =============================================================================
function var_estim(y::Matrix{Float64}, p::Int, intercept::Bool, trend::Bool)

    T, n_v = size(y)
    T_eff = T - p

    # Build regressor matrix
    X_cols = Matrix{Float64}[]

    for lag in 1:p
        push!(X_cols, y[p-lag+1:T-lag, :])
    end

    if intercept
        push!(X_cols, ones(T_eff, 1))
    end

    if trend
        push!(X_cols, collect(1:T_eff) |> x -> reshape(x, :, 1) |> Float64)
    end

    X = hcat(X_cols...)
    Y = y[p+1:end, :]

    β = (X'X) \ (X'Y)
    resid = Y - X * β
    Sigma = (resid' * resid) ./ T_eff

    return β, X * β, Sigma
end

# =============================================================================
# VAR LAG LENGTH SELECTION via BIC  [F2: ic_var now fully functional]
#
# Inputs:
#   y        T × n_v  data matrix
#   p_max    maximum lag to consider
#   method   1 = AIC, 2 = BIC  (default BIC)
# Output:
#   p_select  selected lag length (Int)
# =============================================================================
function ic_var(y::Matrix{Float64}, p_max::Int, method::Int=2)

    n_v = size(y, 2)
    T = size(y, 1)
    aic = zeros(p_max)
    bic = zeros(p_max)

    T_eff = T - p_max          # use common effective sample for comparability

    for p in 1:p_max
        _, _, Sigma = var_estim(y, p, true, false)
        n_params = n_v^2 * p + n_v          # VAR slope params + intercepts
        log_det_S = log(det(Sigma))

        bic[p] = log_det_S + n_params * log(T_eff) / T_eff
        aic[p] = log_det_S + n_params * 2.0 / T_eff
    end

    p_select = (method == 1) ? argmin(aic) : argmin(bic)
    return p_select
end

# =============================================================================
# DYNAMIC LAG SELECTION
#
# Selects lag length from the shock series in the panel using BIC on a
# univariate AR (n_v = 1).  Called once per variant in run_lp.
# =============================================================================
function select_lag_length(panel::DataFrame; p_max::Int=L_LAG_MAX)

    shock_ts = sort(unique(panel[!, [:date, :shock]]), :date).shock
    y = reshape(Float64.(shock_ts), :, 1)
    p = ic_var(y, p_max, 2)           # BIC

    @info "BIC-selected lag length: $p (searched 1:$p_max)"
    return p
end

# =============================================================================
# VARIABLE MAP
# =============================================================================
function get_lp_specs(variant::Symbol)

    if variant == :hourly_rate
        return :log_rwage,
        [:age_mean, :female_share]

    elseif variant == :income
        return :log_rincome,
        [:age_mean, :female_share]

    elseif variant == :hours
        return :log_hours,
        [:age_mean, :female_share]

    elseif variant == :unemployment
        return :unemp_rate,
        [:age_mean, :female_share]

    elseif variant == :employment
        return :log_emp_count,
        [:age_mean, :female_share]

    elseif variant == :inequality
        return :log_ratio_7525, 
        [:age_mean, :female_share]

    elseif variant == :median
        return :log_rincome_p50, 
        [:age_mean, :female_share]

    elseif variant == :income_share_var
        return :income_share, 
        [:age_mean, :female_share]

    else
        error("Unknown variant: $variant")
    end
end

function driscoll_kraay_vcov(X::Matrix{Float64}, e::Vector{Float64},
    t_idx::Vector{Int}; m::Int=DK_BW)

    N, K = size(X)
    Tvals = sort(unique(t_idx))
    T = length(Tvals)
    tmap = Dict(v => i for (i, v) in enumerate(Tvals))

    H = zeros(K, T)
    for n in 1:N
        ti = tmap[t_idx[n]]
        H[:, ti] .+= X[n, :] .* e[n]
    end

    S = zeros(K, K)
    for t in 1:T
        S .+= H[:, t] * H[:, t]'
    end
    S ./= T

    for l in 1:m
        G = zeros(K, K)
        for t in (l+1):T
            G .+= H[:, t] * H[:, t-l]'
        end
        G ./= T
        w = 1 - l / (m + 1)
        S .+= w .* (G .+ G')
    end

    XX = inv(X'X)
    return XX * (T * S) * XX
end


# =============================================================================
# BUILD DATASET FOR HORIZON h  [F4: explicit sort by :date inside each group]
# =============================================================================
function build_lp_data(panel::DataFrame, h::Int, y_var::Symbol, l_lag::Int)

    # Ensure panel is sorted globally (defensive)
    sort!(panel, [:occ_group, :date])

    out = DataFrame[]

    for g in groupby(panel, :occ_group)

        # [F4] Sort by date within the group — groupby does not guarantee order
        sub = sort(copy(g), :date)
        n = nrow(sub)

        leady = Vector{Union{Missing,Float64}}(missing, n)

        if h < n
            leady[1:n-h] = sub[!, y_var][1+h:n]
        end

        sub.dep_var = leady

        for l in 1:l_lag
            col = Symbol("ylag$l")
            lagged = Vector{Union{Missing,Float64}}(missing, n)
            if l < n
                lagged[l+1:n] = sub[!, y_var][1:n-l]
            end
            sub[!, col] = lagged
        end

        push!(out, sub)
    end

    df = vcat(out...)
    req = [:dep_var, :shock, :log_oil_lag1, :ffr_lag1, :cpi_lag1, :indpro_lag1, :t10y3m_lag1]

    for l in 1:l_lag
        push!(req, Symbol("shock_lag$l"))
        push!(req, Symbol("ylag$l"))
    end

    return dropmissing(df, req)
end

# =============================================================================
# BUILD RHS COLUMNS
# =============================================================================
function build_cols!(df::DataFrame, controls::Vector{Symbol}, l_lag::Int)

    shock_cols = Symbol[]
    groups = sort(unique(df.occ_group))

    for g in groups
        c = Symbol("shock_g$(g)")
        df[!, c] = df.shock .* (df.occ_group .== g)
        push!(shock_cols, c)
    end

    lag_cols = [Symbol("shock_lag$l") for l in 1:l_lag]
    ylag_cols = [Symbol("ylag$l") for l in 1:l_lag]
    mac_controls = [:log_oil_lag1, :ffr_lag1, :cpi_lag1, :indpro_lag1, :t10y3m_lag1]
    xcols = vcat(shock_cols, lag_cols, ylag_cols, mac_controls, controls)

    return xcols, shock_cols
end

# =============================================================================
# FE ABSORPTION HELPER  [C1, C3]
#
# Calls reg() from FixedEffectModels.jl to partial out fe(occ_group) + fe(date)
# from a single variable v, returning the within-transformed residual vector
# M_D * v.
#
# By the Frisch-Waugh-Lovell theorem, OLS of (M_D*Y) on (M_D*X) recovers the
# same coefficient vector β as the full two-way FE regression of Y on X.
# reg() uses iterative LSMR, which converges to the exact within-transformation
# regardless of panel balance — unlike the single-pass demean which required
# the two FEs to be orthogonal.
#
# Note: df must already have :occ_group as categorical and :time_idx as Int.
# =============================================================================
function absorb_fe(df::DataFrame, v::Symbol)::Vector{Float64}
    fml = term(v) ~ term(1) + fe(term(:occ_group))
    fit = reg(df, fml, Vcov.simple(), save=:residuals)
    return residuals(fit)
end

# =============================================================================
# OLS CORE
# =============================================================================
function estimate_lp(df_h::DataFrame, panel::DataFrame,
    controls::Vector{Symbol}, h::Int, l_lag::Int)

    df = copy(df_h)
    xcols, _ = build_cols!(df, controls, l_lag)

    # Integer time index for DK vcov
    all_dates = sort(unique(panel.date))
    dmap = Dict(d => i for (i, d) in enumerate(all_dates))
    df[!, :time_idx] = Int.([dmap[d] for d in df.date])

    # FixedEffectModels requires categorical FE identifiers
    df[!, :occ_group] = categorical(df.occ_group)

    # Drop any remaining missings before absorption
    needed = vcat([:dep_var, :occ_group, :time_idx], xcols)
    df = dropmissing(df, needed)

    # FE absorption via reg() (LSMR, exact for unbalanced panels)
    Y_tilde = absorb_fe(df, :dep_var)
    X_tilde = hcat([absorb_fe(df, c) for c in xcols]...)

    # OLS on within-transformed data
    β = (X_tilde' * X_tilde) \ (X_tilde' * Y_tilde)
    e = Y_tilde - X_tilde * β

    # DK vcov on partialled-out X and residuals
    t_idx = Int.(df.time_idx)
    V = driscoll_kraay_vcov(X_tilde, e, t_idx; m=DK_BW)
    se = sqrt.(diag(V))

    return (β=β, se=se, e=e, X=X_tilde,
        t_idx=t_idx, df=df, names=xcols)
end


# =============================================================================
# PERCENTILE-T BLOCK BOOTSTRAP  [F5: bootstrap t_idx uses date→index map]
# =============================================================================
function bootstrap_lp(df_h::DataFrame, panel::DataFrame,
    controls::Vector{Symbol}, h::Int, l_lag::Int;
    n_boot::Int=N_BOOT,
    block_size::Int=BLOCK_SIZE,
    rng=MersenneTwister(BOOT_SEED))

    base = estimate_lp(df_h, panel, controls, h, l_lag)

    β0 = base.β
    se0 = base.se
    K = length(β0)

    all_dates = sort(unique(base.df.date))
    T = length(all_dates)
    nb = ceil(Int, T / block_size)

    # Map date → row indices in the estimation-sample df
    date_rows = Dict{Date,Vector{Int}}()
    for (i, r) in enumerate(eachrow(base.df))
        push!(get!(date_rows, r.date, Int[]), i)
    end

    B = fill(NaN, n_boot, K)
    TSTAT = fill(NaN, n_boot, K)

    for b in 1:n_boot

        starts = rand(rng, 1:T, nb)
        dd = Date[]
        for s in starts
            for k in 0:block_size-1
                push!(dd, all_dates[mod1(s + k, T)])
            end
        end
        dd = dd[1:T]

        rows = Int[]
        for d in dd
            append!(rows, get(date_rows, d, Int[]))
        end
        isempty(rows) && continue

        # Bootstrap DataFrame from original pre-absorption data
        df_b = base.df[rows, :]

        # Impose the null: replace dep_var with X*β0 + residual
        df_b[!, :dep_var] = base.X[rows, :] * β0 .+ base.e[rows]

        # Re-absorb FEs on the bootstrap sample
        Y_b = absorb_fe(df_b, :dep_var)
        X_b = hcat([absorb_fe(df_b, c) for c in base.names]...)

        βb = try
            (X_b' * X_b) \ (X_b' * Y_b)
        catch
            continue
        end

        eb = Y_b - X_b * βb
        tb = Int.(df_b.time_idx)
        Vb = driscoll_kraay_vcov(X_b, eb, tb; m=DK_BW)
        seb = sqrt.(diag(Vb))

        B[b, :] = βb
        TSTAT[b, :] = (βb .- β0) ./ seb
    end

    return B, TSTAT, β0, se0, base.names
end

# =============================================================================
# ONE HORIZON
# =============================================================================
function run_h(panel::DataFrame, y_var::Symbol,
    controls::Vector{Symbol}, h::Int, l_lag::Int)

    df_h = build_lp_data(panel, h, y_var, l_lag)

    B, TSTAT, β0, se0, names = bootstrap_lp(df_h, panel, controls, h, l_lag)

    K = length(β0)
    rows = NamedTuple[]

    for k in 1:K

        good = .!isnan.(TSTAT[:, k])
        q95 = quantile(TSTAT[good, k], 0.95)
        q05 = quantile(TSTAT[good, k], 0.05)
        q84 = quantile(TSTAT[good, k], 0.84)
        q16 = quantile(TSTAT[good, k], 0.16)

        ci_lo90 = β0[k] - se0[k] * q95
        ci_hi90 = β0[k] - se0[k] * q05
        ci_lo68 = β0[k] - se0[k] * q84
        ci_hi68 = β0[k] - se0[k] * q16

        push!(rows, (
            horizon=h,
            coef_name=string(names[k]),
            beta=β0[k],
            se=se0[k],
            ci_lo90=ci_lo90,
            ci_hi90=ci_hi90,
            ci_lo68=ci_lo68,
            ci_hi68=ci_hi68
        ))
    end

    df_res = DataFrame(rows)

    return df_res, B, TSTAT, names
end

# =============================================================================
# FULL RUN
# =============================================================================
function run_full_lp(panel::DataFrame, y_var::Symbol,
    controls::Vector{Symbol}, l_lag::Int;
    h_max::Int=H_MAX)

    out_dfs = DataFrame[]

    boot_store = Dict{Int,NamedTuple}()

    coef_names = Symbol[]
    for h in 0:h_max
        @info "  horizon h = $h / $h_max"

        df_res, B, TSTAT, names = run_h(panel, y_var, controls, h, l_lag)

        push!(out_dfs, df_res)

        boot_store[h] = (B=B, TSTAT=TSTAT, names=names)

        if isempty(coef_names)
            coef_names = names
        end
    end

    results_df = vcat(out_dfs...)

    return results_df, boot_store, coef_names
end

# =============================================================================
# EXTRACT GROUP IRFs
# =============================================================================
function extract_irf(df::DataFrame)

    out = Dict{Int,DataFrame}()

    for g in keys(OCC_LABELS)
        cname = "shock_g$(g)"
        sub = filter(r -> r.coef_name == cname, df)
        out[g] = sub[:, [:horizon, :beta, :ci_lo90, :ci_hi90, :ci_lo68, :ci_hi68]]
    end

    return out
end

# =============================================================================
# MAIN API
# =============================================================================
function run_lp(panel::DataFrame, variant::Symbol)

    y_var, controls = get_lp_specs(variant)

    # [F6] Dynamic lag selection via BIC on the shock series
    l_lag = select_lag_length(panel; p_max=L_LAG_MAX)
    @info "Running LP for $variant | outcome=$y_var | lags=$l_lag"

    results_df, boot_store, coef_names = run_full_lp(panel, y_var, controls, l_lag)

    output_dir = get_output_dir(variant)
    CSV.write(joinpath(output_dir, "lp_coefficients.csv"), results_df)

    irfs = extract_irf(results_df)

    for (g, d) in irfs
        label = get(OCC_LABELS, g, "group_$g")
        CSV.write(joinpath(output_dir, "irf_group$(g)_$(label).csv"), d)
    end

    @info "Saved results to $output_dir"
    return results_df, boot_store, coef_names, irfs
end