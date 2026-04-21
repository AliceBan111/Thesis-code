# =============================================================================
# 02_lp_occ1990_estimation.jl
# PANEL LOCAL PROJECTION — occ1990-level
#
# Model:
#   y_{i,t+h} = α_i + τ_t + Σ_c β_{c,h}(shock_t × 1[occ1990=c]) + controls + ε
#
# Output: ONE wide CSV
#   columns: variable | horizon | {occ1990_code_1} | {occ1990_code_2} | ...
#
# Design preserved from original:
#   - Two-way FE: occ1990 FE + time FE
#   - Driscoll-Kraay SE
#   - Moving block bootstrap
#   - BIC lag selection
# =============================================================================

using CSV, DataFrames, DataFramesMeta
using LinearAlgebra, Statistics, StatsBase
using Dates, Random
using Printf
using FixedEffectModels
using CategoricalArrays

# =============================================================================
# CONSTANTS
# =============================================================================
const N_BOOT      = 500
const BLOCK_SIZE  = 6
const BOOT_SEED   = 42
const H_MAX       = 36
const L_LAG_MAX   = 24
const DK_BW       = 12

safe_float(v) = Float64.(collect(skipmissing(v)))

# =============================================================================
# VARIABLE MAP
# =============================================================================
function get_lp_specs(variant::Symbol)
    if variant == :hourly_rate
        return :log_rwage,    [:age_mean, :female_share]
    elseif variant == :income
        return :log_rincome,  [:age_mean, :female_share]
    elseif variant == :hours
        return :log_hours,    [:age_mean, :female_share]
    elseif variant == :unemployment
        return :unemp_rate,   [:age_mean, :female_share]
    elseif variant == :employment
        return :log_emp_count,[:age_mean, :female_share]
    elseif variant == :inequality
        return :log_ratio_7525,[:age_mean, :female_share]
    elseif variant == :median
        return :log_rincome_p50,[:age_mean, :female_share]
    elseif variant == :income_share_var
        return :income_share, [:age_mean, :female_share]
    else
        error("Unknown variant: $variant")
    end
end

# =============================================================================
# VAR / BIC LAG SELECTION
# =============================================================================
function var_estim(y::Matrix{Float64}, p::Int, intercept::Bool, trend::Bool)
    T, n_v = size(y)
    T_eff  = T - p
    X_cols = Matrix{Float64}[]
    for lag in 1:p
        push!(X_cols, y[p-lag+1:T-lag, :])
    end
    intercept && push!(X_cols, ones(T_eff, 1))
    trend && push!(X_cols, reshape(Float64.(1:T_eff), :, 1))
    X = hcat(X_cols...)
    Y = y[p+1:end, :]
    β = (X'X) \ (X'Y)
    resid = Y - X * β
    Sigma = (resid' * resid) ./ T_eff
    return β, X * β, Sigma
end

function ic_var(y::Matrix{Float64}, p_max::Int, method::Int=2)
    n_v   = size(y, 2)
    T     = size(y, 1)
    T_eff = T - p_max
    aic   = zeros(p_max)
    bic   = zeros(p_max)

    for p in 1:p_max
        _, _, Sigma = var_estim(y, p, true, false)
        n_params = n_v^2 * p + n_v
        log_det_S = log(det(Sigma))
        bic[p] = log_det_S + n_params * log(T_eff) / T_eff
        aic[p] = log_det_S + n_params * 2.0 / T_eff
    end

    return (method == 1) ? argmin(aic) : argmin(bic)
end

function select_lag_length(panel::DataFrame; p_max::Int=L_LAG_MAX)
    shock_ts = sort(unique(panel[!, [:date, :shock]]), :date).shock
    y = reshape(Float64.(shock_ts), :, 1)
    p = ic_var(y, p_max, 2)
    @info "BIC-selected lag length: $p"
    return p
end

# =============================================================================
# DRISCOLL-KRAAY VCOV
# =============================================================================
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
# FE ABSORPTION
# =============================================================================
function absorb_fe(df::DataFrame, v::Symbol)::Vector{Float64}

    df2 = dropmissing(df, v)
    nrow(df2) == 0 && return Float64[]

    df2[!, v] = Float64.(collect(skipmissing(df2[!, v])))

    fml = term(v) ~ term(1) + fe(term(:occ1990))

    fit = reg(df2, fml, Vcov.simple(), save=:residuals)

    r = residuals(fit)

    return Float64.(r)
end

# =============================================================================
# BUILD LHS + LAGS
# =============================================================================
function build_lp_data(panel::DataFrame, h::Int, y_var::Symbol, l_lag::Int)

    sort!(panel, [:occ1990, :date])

    out = DataFrame[]

    for g in groupby(panel, :occ1990)

        sub = sort(copy(g), :date)
        n = nrow(sub)

        leady = Vector{Union{Missing,Float64}}(missing, n)
        h < n && (leady[1:n-h] = sub[!, y_var][1+h:n])

        sub.dep_var = leady

        for l in 1:l_lag
            col = Symbol("ylag$l")
            lagged = Vector{Union{Missing,Float64}}(missing, n)
            l < n && (lagged[l+1:n] = sub[!, y_var][1:n-l])
            sub[!, col] = lagged
        end

        push!(out, sub)
    end

    df = vcat(out...)

    req = [:dep_var, :shock, :log_oil_lag1, :ffr_lag1,
           :cpi_lag1, :indpro_lag1, :t10y3m_lag1]

    for l in 1:l_lag
        push!(req, Symbol("shock_lag$l"))
        push!(req, Symbol("ylag$l"))
    end

    return dropmissing(df, req)
end

# =============================================================================
# BUILD RHS INTERACTIONS
# =============================================================================
function build_cols!(df::DataFrame, controls::Vector{Symbol}, l_lag::Int)

    occ_codes = sort(unique(df.occ1990))
    shock_cols = Symbol[]

    shock = coalesce.(df.shock, 0.0)

    for c in occ_codes
        col = Symbol("shock_c$(c)")
        df[!, col] = shock .* (df.occ1990 .== c)
        push!(shock_cols, col)
    end

    lag_cols  = [Symbol("shock_lag$l") for l in 1:l_lag]
    ylag_cols = [Symbol("ylag$l") for l in 1:l_lag]

    mac_cols = [:log_oil_lag1, :ffr_lag1, :cpi_lag1,
                :indpro_lag1, :t10y3m_lag1]

    return vcat(shock_cols, lag_cols, ylag_cols, mac_cols, controls),
           shock_cols,
           occ_codes
end

# =============================================================================
# OLS CORE
# =============================================================================
function estimate_lp(df_h::DataFrame, panel::DataFrame,
                     controls::Vector{Symbol}, h::Int, l_lag::Int)

    df = copy(df_h)

    xcols, _, occ_codes = build_cols!(df, controls, l_lag)

    # time index
    all_dates = sort(unique(panel.date))
    dmap = Dict(d => i for (i, d) in enumerate(all_dates))
    df[!, :time_idx] = Int.([dmap[d] for d in df.date])

    # drop missing on required cols
    needed = unique(vcat([:dep_var, :occ1990, :time_idx], xcols, controls))
    needed = filter(c -> c in names(df), needed)
    df = dropmissing(df, needed)

    nrow(df) == 0 && error("empty sample after dropmissing at h=$h")

    # FE residualization
    Y_tilde = absorb_fe(df, :dep_var)
    length(Y_tilde) == 0 && return nothing

    X_list = Vector{Vector{Float64}}()

    for c in xcols
        xc = absorb_fe(df, c)
        length(xc) == 0 && return nothing
        push!(X_list, xc)
    end

    X_tilde = hcat(X_list...)

    # ================================
    # DROP BAD COLUMNS (CRITICAL FIX)
    # ================================
    keep = [std(X_tilde[:, j]) > 1e-10 for j in 1:size(X_tilde,2)]

    X_tilde = X_tilde[:, keep]
    xcols   = xcols[keep]

    # still need enough regressors
    size(X_tilde,2) == 0 && return nothing

    # ================================
    # OLS with singular protection
    # ================================
    XtX = X_tilde' * X_tilde

    β = try
        cholesky(XtX) \ (X_tilde' * Y_tilde)
    catch
        pinv(X_tilde) * Y_tilde
    end

    e = Y_tilde - X_tilde * β

    t_idx = Int.(df.time_idx)

    V = try
        driscoll_kraay_vcov(X_tilde, e, t_idx; m=DK_BW)
    catch
        I = Matrix{Float64}(I, length(β), length(β))
        I
    end

    se = sqrt.(diag(V))

    return (
        β = β,
        se = se,
        e = e,
        X = X_tilde,
        t_idx = t_idx,
        df = df,
        names = xcols,
        occ_codes = occ_codes
    )
end

# =============================================================================
# BLOCK BOOTSTRAP
# =============================================================================
function bootstrap_lp(df_h::DataFrame, panel::DataFrame,
                      controls::Vector{Symbol}, h::Int, l_lag::Int;
                      n_boot=500, block_size=6,
                      rng=MersenneTwister(42))

    base = estimate_lp(df_h, panel, controls, h, l_lag)

    β0 = base.β
    K  = length(β0)

    dates = sort(unique(base.df.date))
    T = length(dates)

    nb = ceil(Int, T / block_size)

    date_map = Dict(d => i for (i,d) in enumerate(dates))
    row_map = Dict{Date,Vector{Int}}()

    for (i,r) in enumerate(eachrow(base.df))
        push!(get!(row_map, r.date, Int[]), i)
    end

    B = fill(NaN, n_boot, K)

    for b in 1:n_boot

        starts = rand(rng, 1:T, nb)
        dd = Date[]

        for s in starts, k in 0:block_size-1
            push!(dd, dates[mod1(s+k, T)])
        end

        dd = dd[1:T]

        rows = reduce(vcat, [get(row_map, d, Int[]) for d in dd])

        isempty(rows) && continue

        dfb = base.df[rows, :]

        dfb[!, :dep_var] = base.X[rows, :] * β0 .+ base.e[rows]

        Yb = absorb_fe(dfb, :dep_var)

        Xb = hcat([absorb_fe(dfb, c) for c in base.names]...)

        βb = try (Xb'Xb) \ (Xb'Yb) catch; continue end

        B[b,:] = βb
    end

    return B, β0, base.names, base.occ_codes
end

# =============================================================================
# 单独 occ1990 code 列名
# =============================================================================
function occ_colname(code::Int)::Symbol
    return Symbol(string(code))
end

# =============================================================================
# ONE HORIZON
# =============================================================================
function run_h(panel::DataFrame, y_var::Symbol,
               controls::Vector{Symbol}, h::Int, l_lag::Int)

    df_h = build_lp_data(panel, h, y_var, l_lag)

    # skip empty / too small sample
    if nrow(df_h) < 100
        @warn "skip h=$h (too few obs)"
        return Dict(:horizon => h), Int[]
    end

    res = bootstrap_lp(df_h, panel, controls, h, l_lag)

    if res === nothing
        @warn "skip h=$h (singular)"
        return Dict(:horizon => h), Int[]
    end

    B, TSTAT, β0, se0, names, occ_codes = res

    row = Dict{Symbol,Any}(:horizon => h)

    for (k, nm) in enumerate(names)
        s = string(nm)
        if startswith(s, "shock_c")
            code = parse(Int, s[8:end])
            row[occ_colname(code)] = β0[k]
        end
    end

    return row, occ_codes
end

# =============================================================================
# FULL RUN
# =============================================================================
function run_full_lp(panel::DataFrame, y_var::Symbol,
                     controls::Vector{Symbol}, l_lag::Int;
                     h_max::Int=H_MAX)

    rows = Vector{Dict{Symbol,Any}}()
    occ_codes = Int[]

    for h in 0:h_max

        @info "  horizon h = $h / $h_max"

        row, codes = run_h(panel, y_var, controls, h, l_lag)

        push!(rows, row)

        isempty(occ_codes) && (occ_codes = sort(codes))
    end

    df = DataFrame(horizon = [r[:horizon] for r in rows])

    for c in occ_codes
        col = occ_colname(c)
        df[!, col] = [get(r, col, missing) for r in rows]
    end

    return df, occ_codes
end

# =============================================================================
# run_lp
# =============================================================================
function run_lp(panel::DataFrame, variant::Symbol)

    y_var, controls = get_lp_specs(variant)

    l_lag = select_lag_length(panel; p_max=L_LAG_MAX)

    @info "Running LP | variant=$variant | outcome=$y_var | lags=$l_lag"

    wide_df, occ_codes = run_full_lp(panel, y_var, controls, l_lag)

    insertcols!(wide_df, 1,
        :variable => fill(string(variant), nrow(wide_df)))

    @info "Done: variant=$variant | $(nrow(wide_df)) horizons × $(length(occ_codes)) occ codes"

    return wide_df
end

# =============================================================================
# COMBINED RUNNER
# =============================================================================
function run_all_variants(panels::Dict)

    all_dfs = DataFrame[]

    variant_order = [:hourly_rate, :hours, :income, :income_share_var,
                     :inequality, :median, :unemployment, :employment]

    for variant in variant_order
        haskey(panels, variant) || continue
        wide = run_lp(panels[variant], variant)
        push!(all_dfs, wide)
    end

    master = vcat(all_dfs..., cols=:union)

    base = joinpath(@__DIR__, "../..", "result", "occ1990")
    mkpath(base)

    master_path = joinpath(base, "lp_betas_all_variants.csv")

    CSV.write(master_path, master)

    @info "Master CSV saved → $master_path ($(nrow(master)) rows × $(ncol(master)) cols)"

    return master
end