# =============================================================================
# 02_lp_estimation_occ_ind.jl
# PANEL LOCAL PROJECTION: Occupation groups within merged Industry groups
#
# Model (run separately per merged ind_group):
#   y_{g,t+h} = α_g + Σ_g β_{g,h}(shock_t × 1[occ=g]) + macro_controls + ε
#
# Design choices:
#   [1] No baseline occ group — all groups estimated jointly via interactions
#   [2] LHS = y_{t+h} (level, not differenced)
#   [3] One-way FE: occ_group absorbed via manual group demean (NOT reg())
#       — reg() is unstable on small bootstrap samples; demean is exact & fast
#   [4] Macro controls: log_oil_lag1, ffr_lag1, cpi_lag1, indpro_lag1, t10y3m_lag1
#   [5] Moving block bootstrap (circular), percentile-t CI
#   [6] DK robust SE (bandwidth = DK_BW)
#   [7] Lag length selected via BIC on shock series (global, shared across industries)
#   [8] Skip ind_group if fewer than MIN_OCC_GROUPS occ groups present
#   [9] Bootstrap guard: skip draw if any occ_group drops out of sample
# =============================================================================

using CSV, DataFrames, DataFramesMeta
using LinearAlgebra, Statistics, StatsBase
using Dates, Random
using Printf
using CategoricalArrays

# =============================================================================
# CONSTANTS
# =============================================================================
const N_BOOT         = 500
const BLOCK_SIZE     = 6
const BOOT_SEED      = 42
const H_MAX          = 36
const L_LAG_MAX      = 24
const DK_BW          = 12
const MIN_OCC_GROUPS = 3

# =============================================================================
# INDUSTRY MERGE  (12-group → 4-group by oil intensity)
#   1 => Energy-intensive      : Mining (2), Transportation/Utilities (6)
#   2 => Manufacturing & Const : Construction (3), Mfg nondurable (4), durable (5)
#   3 => Trade                 : Wholesale (7), Retail (8)
#   4 => Services              : Finance (9), Business (10), Personal (11),
#                                Professional (12), Agriculture (1)
# =============================================================================
# const IND_MERGE = Dict(
#     1  => 4,
#     2  => 1,
#     3  => 2,
#     4  => 2,
#     5  => 2,
#     6  => 1,
#     7  => 3,
#     8  => 3,
#     9  => 4,
#     10 => 4,
#     11 => 4,
#     12 => 4,
# )

const IND_LABELS = Dict(
    1 => "Energy_intensive",
    2 => "Manufacturing_Construction",
    3 => "Trade",
    4 => "Services",
)

# const IND_LABELS = Dict(
#     1  => "Agriculture_forestry_fishing",
#     2  => "Mining",
#     3  => "Construction",
#     4  => "Manufacturing_nondurable",
#     5  => "Manufacturing_durable",
#     6  => "Transportation_utilities",
#     7  => "Wholesale_trade",
#     8  => "Retail_trade",
#     9  => "Finance_insurance_realestate",
#     10 => "Business_repair_services",
#     11 => "Personal_entertainment_services",
#     12 => "Professional_related_services",
# )

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
# OUTPUT DIRECTORY
# =============================================================================
function get_output_dir(variant::Symbol, ind_group::Int)::String
    label = get(IND_LABELS, ind_group, "industry_$(ind_group)")
    dir   = joinpath(@__DIR__, "../..", "result", "occ_ind", string(variant), label)
    mkpath(dir)
    return dir
end

# =============================================================================
# INDUSTRY MERGE
# =============================================================================
# function merge_industry_groups!(panel::DataFrame)::DataFrame
#     panel[!, :ind_merged] = [get(IND_MERGE, g, missing) for g in panel.ind_group]
#     panel = dropmissing(panel, :ind_merged)
#     panel[!, :ind_merged] = Int.(panel.ind_merged)
#     return panel
# end

# =============================================================================
# VARIABLE MAP
# =============================================================================
function get_lp_specs(variant::Symbol)
    specs = Dict(
        :hourly_rate      => (:log_rwage,      [:age_mean, :female_share]),
        :income           => (:log_rincome,     [:age_mean, :female_share]),
        :hours            => (:log_hours,       [:age_mean, :female_share]),
        :unemployment     => (:unemp_rate,      [:age_mean, :female_share]),
        :employment       => (:log_emp_count,   [:age_mean, :female_share]),
        :inequality       => (:log_ratio_7525,  [:age_mean, :female_share]),
        :median           => (:log_rincome_p50, [:age_mean, :female_share]),
        :income_share_var => (:income_share,    [:age_mean, :female_share]),
    )
    haskey(specs, variant) || error("Unknown variant: $variant")
    return specs[variant]
end

# =============================================================================
# VAR HELPERS FOR BIC LAG SELECTION
# =============================================================================
function var_estim(y::Matrix{Float64}, p::Int, intercept::Bool, trend::Bool)
    T, n_v = size(y)
    T_eff  = T - p
    X_cols = Matrix{Float64}[]
    for lag in 1:p
        push!(X_cols, y[p-lag+1:T-lag, :])
    end
    intercept && push!(X_cols, ones(T_eff, 1))
    trend     && push!(X_cols, reshape(Float64.(1:T_eff), :, 1))
    X     = hcat(X_cols...)
    Y     = y[p+1:end, :]
    β     = (X'X) \ (X'Y)
    resid = Y - X * β
    Sigma = (resid' * resid) ./ T_eff
    return β, X * β, Sigma
end

function ic_var(y::Matrix{Float64}, p_max::Int, method::Int=2)
    T = size(y, 1)
    n_v   = size(y, 2)
    aic   = zeros(p_max)
    bic   = zeros(p_max)

    T_eff = T - p_max

    for p in 1:p_max
        _, _, Sigma = var_estim(y, p, true, false)
        n_params = n_v^2 * p + n_v
        bic[p]   = log(det(Sigma)) + n_params * log(T_eff) / T_eff
        aic[p]   = log(det(Sigma)) + n_params * 2.0 / T_eff
    end
    return method == 1 ? argmin(aic) : argmin(bic)
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
                              t_idx::Vector{Int}; m::Int=DK_BW)::Matrix{Float64}
    N, K  = size(X)
    Tvals = sort(unique(t_idx))
    T     = length(Tvals)
    tmap  = Dict(v => i for (i, v) in enumerate(Tvals))

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
        w  = 1 - l / (m + 1)
        S .+= w .* (G .+ G')
    end

    XX = inv(X'X)
    return XX * (T * S) * XX
end

# =============================================================================
# FE ABSORPTION via manual group demean
#
# Replaces reg() — reg() produces missing residuals when a bootstrap sample
# has too few obs per group, poisoning the entire draw.
# Manual demean is exact, never produces missing, and is faster.
# =============================================================================
function demean_by_group(x::Vector{Float64}, group::Vector)::Vector{Float64}
    out = copy(x)
    for g in unique(group)
        idx = findall(group .== g)
        if length(idx) < 2
            out[idx] .= 0.0
        else
            out[idx] .= x[idx] .- mean(x[idx])   # singleton: demeaned contribution = 0
        end
    end
    return out
end

function absorb_fe(df::DataFrame, v::Symbol)::Vector{Float64}
    x = Float64.(df[!, v])
    return demean_by_group(x, Vector(df.occ_group))
end

# =============================================================================
# BUILD LP DATASET FOR HORIZON h
# =============================================================================
function build_lp_data(panel::DataFrame, h::Int, y_var::Symbol, l_lag::Int)
    sort!(panel, [:occ_group, :date])
    out = DataFrame[]

    for g in groupby(panel, :occ_group)
        sub = sort(copy(g), :date)
        n   = nrow(sub)

        leady = Vector{Union{Missing,Float64}}(missing, n)
        h < n && (leady[1:n-h] = sub[!, y_var][1+h:n])
        sub[!, :dep_var] = leady

        for l in 1:l_lag
            col    = Symbol("ylag$(l)")
            lagged = Vector{Union{Missing,Float64}}(missing, n)
            l < n  && (lagged[l+1:n] = sub[!, y_var][1:n-l])
            sub[!, col] = lagged
        end

        push!(out, sub)
    end

    df  = vcat(out...)
    req = [:dep_var, :shock, :log_oil_lag1, :ffr_lag1,
           :cpi_lag1, :indpro_lag1, :t10y3m_lag1]
    for l in 1:l_lag
        push!(req, Symbol("shock_lag$(l)"))
        push!(req, Symbol("ylag$(l)"))
    end
    return dropmissing(df, req)
end

# =============================================================================
# BUILD RHS COLUMNS  (shock × occ_group interactions, no baseline dropped)
# =============================================================================
function build_cols!(df::DataFrame, controls::Vector{Symbol}, l_lag::Int)
    shock_cols = Symbol[]
    groups     = sort(unique(df.occ_group))

    for g in groups
        c = Symbol("shock_g$(g)")
        df[!, c] = Float64.(df.shock .* (df.occ_group .== g))
        push!(shock_cols, c)
    end

    lag_cols     = [Symbol("shock_lag$(l)") for l in 1:l_lag]
    ylag_cols    = [Symbol("ylag$(l)")      for l in 1:l_lag]
    mac_controls = [:log_oil_lag1, :ffr_lag1, :cpi_lag1, :indpro_lag1, :t10y3m_lag1]
    xcols        = vcat(shock_cols, lag_cols, ylag_cols, mac_controls, controls)

    return xcols, shock_cols
end

function base_lp_cols(controls::Vector{Symbol}, l_lag::Int)::Vector{Symbol}
    lag_cols     = [Symbol("shock_lag$(l)") for l in 1:l_lag]
    ylag_cols    = [Symbol("ylag$(l)")      for l in 1:l_lag]
    mac_controls = [:log_oil_lag1, :ffr_lag1, :cpi_lag1, :indpro_lag1, :t10y3m_lag1]
    return vcat(lag_cols, ylag_cols, mac_controls, controls)
end

# =============================================================================
# OLS CORE
# =============================================================================
function estimate_lp(df_h::DataFrame, panel::DataFrame,
                     controls::Vector{Symbol}, l_lag::Int)

    df = copy(df_h)

    all_dates        = sort(unique(panel.date))
    dmap             = Dict(d => i for (i, d) in enumerate(all_dates))
    df[!, :time_idx] = Int.([dmap[d] for d in df.date])

    # Drop rows with missing controls/lags before constructing occupation
    # interactions. Otherwise occupations that vanish after dropmissing leave
    # all-zero shock_g columns after FE absorption and blank the whole horizon.
    base_needed = vcat([:dep_var, :shock, :occ_group, :time_idx],
                       base_lp_cols(controls, l_lag))
    df = dropmissing(df, base_needed)
    nrow(df) == 0 && return nothing

    xcols, _ = build_cols!(df, controls, l_lag)

    needed = vcat([:dep_var, :occ_group, :time_idx], xcols)
    df     = dropmissing(df, needed)
    nrow(df) == 0 && return nothing

    # ensure Float64 before demean
    df[!, :dep_var] = Float64.(df.dep_var)
    for c in xcols
        df[!, c] = Float64.(df[!, c])
    end

    Y_tilde = absorb_fe(df, :dep_var)
    X_tilde = hcat([absorb_fe(df, c) for c in xcols]...)

    # guard: drop degenerate columns after demean. This mainly happens when a
    # sparse occupation has too few usable rows after controls/lags are applied.
    col_norms = vec(sum(X_tilde .^ 2, dims=1))
    keep = isfinite.(col_norms) .& (col_norms .>= 1e-14)
    any(keep) || return nothing
    if !all(keep)
        X_tilde = X_tilde[:, keep]
        xcols = xcols[keep]
    end
    any(startswith.(string.(xcols), "shock_g")) || return nothing

    β = try
        (X_tilde' * X_tilde) \ (X_tilde' * Y_tilde)
    catch
        return nothing
    end

    e  = Y_tilde - X_tilde * β
    V  = driscoll_kraay_vcov(X_tilde, e, Int.(df.time_idx); m=DK_BW)
    se = sqrt.(abs.(diag(V)))

    return (β=β, se=se, e=e, X=X_tilde,
            t_idx=Int.(df.time_idx), df=df, names=xcols)
end

# =============================================================================
# PERCENTILE-T BLOCK BOOTSTRAP
# =============================================================================
function bootstrap_lp(df_h::DataFrame, panel::DataFrame,
                      controls::Vector{Symbol}, l_lag::Int;
                      n_boot::Int=N_BOOT,
                      block_size::Int=BLOCK_SIZE,
                      rng=MersenneTwister(BOOT_SEED))

    base = estimate_lp(df_h, panel, controls, l_lag)
    isnothing(base) && return nothing

    β0  = base.β
    se0 = base.se
    K   = length(β0)

    all_dates     = sort(unique(base.df.date))
    T             = length(all_dates)
    nb            = ceil(Int, T / block_size)
    expected_occs = sort(unique(base.df.occ_group))

    date_rows = Dict{Date,Vector{Int}}()
    for (i, r) in enumerate(eachrow(base.df))
        push!(get!(date_rows, r.date, Int[]), i)
    end

    B     = fill(NaN, n_boot, K)
    TSTAT = fill(NaN, n_boot, K)

    for b in 1:n_boot
        # circular block bootstrap
        starts = rand(rng, 1:T, nb)
        dd     = Date[]
        for s in starts, k in 0:block_size-1
            push!(dd, all_dates[mod1(s + k, T)])
        end
        dd   = dd[1:T]
        rows = reduce(vcat, [get(date_rows, d, Int[]) for d in dd])
        isempty(rows) && continue

        df_b = copy(base.df[rows, :])

        # [9] skip draw if any occ_group has dropped out
        sampled_occs = sort(unique(df_b.occ_group))
        sampled_occs != expected_occs && continue

        # impose the null
        df_b[!, :dep_var] = Float64.(base.X[rows, :] * β0 .+ base.e[rows])

        Y_b = absorb_fe(df_b, :dep_var)
        X_b = hcat([absorb_fe(df_b, c) for c in base.names]...)

        # guard: degenerate columns
        col_norms = vec(sum(X_b .^ 2, dims=1))
        any(col_norms .< 1e-14) && continue

        βb = try
            (X_b' * X_b) \ (X_b' * Y_b)
        catch
            continue
        end

        eb  = Y_b - X_b * βb
        Vb  = driscoll_kraay_vcov(X_b, eb, Int.(df_b.time_idx); m=DK_BW)
        seb = sqrt.(abs.(diag(Vb)))

        any(seb .< 1e-14) && continue

        B[b, :]     = βb
        TSTAT[b, :] = (βb .- β0) ./ seb
    end

    return B, TSTAT, β0, se0, base.names
end

# =============================================================================
# ONE HORIZON
# =============================================================================
function run_h(panel::DataFrame, y_var::Symbol,
               controls::Vector{Symbol}, h::Int, l_lag::Int)

    occ_groups = sort(unique(panel.occ_group))

    function blank_result()
        rows = NamedTuple[]
        for g in occ_groups
            push!(rows, (horizon=h, coef_name="shock_g$(g)",
                         beta=NaN, se=NaN,
                         ci_lo90=NaN, ci_hi90=NaN,
                         ci_lo68=NaN, ci_hi68=NaN))
        end
        return DataFrame(rows), fill(NaN,0,0), fill(NaN,0,0), Symbol[]
    end

    df_h = build_lp_data(panel, h, y_var, l_lag)
    nrow(df_h) == 0 && return blank_result()

    boot_out = try
        bootstrap_lp(df_h, panel, controls, l_lag)
    catch err
        @warn "Bootstrap failed at h=$h" exception=(err, catch_backtrace())
        nothing
    end
    isnothing(boot_out) && return blank_result()

    B, TSTAT, β0, se0, names = boot_out
    rows = NamedTuple[]

    for k in eachindex(β0)
        good   = .!isnan.(TSTAT[:, k])
        n_good = sum(good)

        if n_good < 50   # too few valid draws → unreliable CI
            push!(rows, (horizon=h, coef_name=string(names[k]),
                         beta=β0[k], se=se0[k],
                         ci_lo90=NaN, ci_hi90=NaN,
                         ci_lo68=NaN, ci_hi68=NaN))
            continue
        end

        q95 = quantile(TSTAT[good, k], 0.95)
        q05 = quantile(TSTAT[good, k], 0.05)
        q84 = quantile(TSTAT[good, k], 0.84)
        q16 = quantile(TSTAT[good, k], 0.16)

        push!(rows, (
            horizon   = h,
            coef_name = string(names[k]),
            beta      = β0[k],
            se        = se0[k],
            ci_lo90   = β0[k] - se0[k] * q95,
            ci_hi90   = β0[k] - se0[k] * q05,
            ci_lo68   = β0[k] - se0[k] * q84,
            ci_hi68   = β0[k] - se0[k] * q16,
        ))
    end

    return DataFrame(rows), B, TSTAT, names
end

# =============================================================================
# FULL LP LOOP OVER HORIZONS
# =============================================================================
function run_full_lp(panel::DataFrame, y_var::Symbol,
                     controls::Vector{Symbol}, l_lag::Int;
                     h_max::Int=H_MAX)

    out_dfs    = DataFrame[]
    boot_store = Dict{Int,NamedTuple}()
    coef_names = Symbol[]

    for h in 0:h_max
        @info "  horizon h=$h / $h_max"
        df_res, B, TSTAT, names = run_h(panel, y_var, controls, h, l_lag)
        push!(out_dfs, df_res)
        boot_store[h] = (B=B, TSTAT=TSTAT, names=names)
        isempty(coef_names) && (coef_names = names)
    end

    return vcat(out_dfs...), boot_store, coef_names
end

# =============================================================================
# EXTRACT IRF PER OCC GROUP
# =============================================================================
function extract_irf(df::DataFrame, occ_groups::Vector{Int})
    out = Dict{Int,DataFrame}()
    for g in occ_groups
        cname = "shock_g$(g)"
        sub   = filter(r -> r.coef_name == cname, df)
        out[g] = nrow(sub) == 0 ?
            DataFrame(horizon = collect(0:H_MAX),
                      beta    = fill(NaN, H_MAX+1),
                      ci_lo90 = fill(NaN, H_MAX+1),
                      ci_hi90 = fill(NaN, H_MAX+1),
                      ci_lo68 = fill(NaN, H_MAX+1),
                      ci_hi68 = fill(NaN, H_MAX+1)) :
            sub[:, [:horizon, :beta, :ci_lo90, :ci_hi90, :ci_lo68, :ci_hi68]]
    end
    return out
end

# =============================================================================
# MAIN API
# =============================================================================
function run_lp(panel::DataFrame, variant::Symbol)

    y_var, controls = get_lp_specs(variant)

    # panel = merge_industry_groups!(copy(panel))
    panel = copy(panel)
    l_lag = select_lag_length(panel; p_max=L_LAG_MAX)
    @info "variant=$variant | outcome=$y_var | lags=$l_lag"

    all_results = Dict{Int,DataFrame}()

    for ind in sort(unique(panel.ind_group))
        # panel_ind  = subset(panel, :ind_merged => x -> x .== ind)
        panel_ind  = subset(panel, :ind_group => x -> x .== ind)
        occ_groups = sort(unique(panel_ind.occ_group))
        ind_label  = get(IND_LABELS, ind, "industry_$(ind)")

        if length(occ_groups) < MIN_OCC_GROUPS
            @warn "Skipping $ind_label: $(length(occ_groups)) occ group(s) < $MIN_OCC_GROUPS"
            continue
        end

        @info "Industry: $ind_label | occ groups: $occ_groups"

        # results_df, boot_store, coef_names = try
        #     run_full_lp(panel_ind, y_var, controls, l_lag)
        # catch err
        #     @warn "LP failed for $ind_label" exception=(err, catch_backtrace())
        #     continue
        # end

        results_df, boot_store, coef_names = run_full_lp(panel_ind, y_var, controls, l_lag)
        

        irfs       = extract_irf(results_df, occ_groups)
        output_dir = get_output_dir(variant, ind)

        CSV.write(joinpath(output_dir, "lp_coefficients.csv"), results_df)
        for g in occ_groups
            occ_label = get(OCC_LABELS, g, "group_$(g)")
            CSV.write(
                joinpath(output_dir, "irf_occ$(g)_$(occ_label).csv"),
                irfs[g])
        end

        @info "Saved → $output_dir"
        all_results[ind] = results_df
    end

    return all_results
end

# =============================================================================
# CONVENIENCE
# =============================================================================
function run_all_variants(panel::DataFrame,
    variants::Vector{Symbol} = [
        :hourly_rate, :hours, :income,
        :income_share_var, :inequality, :median,
        :unemployment, :employment,
    ])
    for v in variants
        @info "===== Variant: $v ====="
        run_lp(panel, v)
    end
end
