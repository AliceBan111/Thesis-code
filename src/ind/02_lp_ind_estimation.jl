# =============================================================================
# lp_estimation_ind_v2.jl
# Industry-Level Local Projection — 方法对齐 occ 版本（文件2）
#
# 相对于文件4（02_lp_estimation_unified.jl）的核心改动：
#   [M1] FE absorption: 用 FixedEffectModels.jl reg() 替代手工 within_transform
#        (Frisch-Waugh-Lovell，精确处理非平衡面板)
#   [M2] Bootstrap: percentile-t bootstrap 替代 percentile bootstrap
#        (t_stat = (β* - β̂) / se*，更准确的 size 控制)
#   [M3] Lag length: BIC 动态选择替代固定 L_LAG=12
#   [M4] 行业组: 12组（删除 group 13）
#   [M5] 新增: 截面 OilShare 回归分析
#        - 因变量: 累积 IRF (h=1:36)
#        - 解释变量: oil_intensity (ind 层面) / oil_exposure (occ 层面)
#        - 报告 R², δ 系数, variance share
#
# 调用方式:
#   include("lp_estimation_ind_v2.jl")
#   results_df, boot_store, coef_names, irfs = run_lp(panel, :hourly_rate)
#   cross_results = run_cross_sectional(irfs, ind_intensity_df)
# =============================================================================

using CSV, DataFrames, DataFramesMeta
using LinearAlgebra, Statistics, StatsBase
using Dates, Random
using Printf
using FixedEffectModels
using CategoricalArrays
using Distributions

# =============================================================================
# 0. GLOBALS
# =============================================================================
const N_BOOT      = 500
const BLOCK_SIZE  = 6
const BOOT_SEED   = 42
const H_MAX       = 36
const L_LAG_MAX   = 24
const DK_BW       = 12

const IND_LABELS = Dict(
    1  => "Agriculture_forestry_fishing",
    2  => "Mining",
    3  => "Construction",
    4  => "Manufacturing_nondurable",
    5  => "Manufacturing_durable",
    6  => "Transportation_utilities",
    7  => "Wholesale_trade",
    8  => "Retail_trade",
    9  => "Finance_insurance_realestate",
    10 => "Business_repair_services",
    11 => "Personal_entertainment_services",
    12 => "Professional_related_services",
    # group 13 已删除
)

# =============================================================================
# 1. OUTPUT DIRECTORY
# =============================================================================
function get_output_dir(variant::Symbol)::String
    base = joinpath(@__DIR__, "../..", "result", "ind_v2")
    dir  = joinpath(base, string(variant))
    mkpath(dir)
    return dir
end

# =============================================================================
# 2. VARIABLE MAP
# =============================================================================
function get_lp_specs(variant::Symbol)
    if variant == :hourly_rate
        return :log_rwage,       [:age_mean, :female_share]
    elseif variant == :income
        return :log_rincome,     [:age_mean, :female_share]
    elseif variant == :hours
        return :log_hours,       [:age_mean, :female_share]
    elseif variant == :unemployment
        return :unemp_rate,      [:age_mean, :female_share]
    elseif variant == :employment
        return :log_emp_count,   [:age_mean, :female_share]
    elseif variant == :inequality
        return :log_ratio_7525,  [:age_mean, :female_share]
    elseif variant == :median
        return :log_rincome_p50, [:age_mean, :female_share]
    elseif variant == :income_share_var
        return :income_share, 
        [:age_mean, :female_share]
    else
        error("Unknown variant: $variant")
    end
end

function load_occ_irfs(variant::Symbol;
                        base_dir::String = joinpath(@__DIR__, "../..", "result", "occ"))::Dict{Int, DataFrame}

    fpath = joinpath(base_dir, "analysis", "merged_occ_irf_trajectories.csv")
    df    = CSV.read(fpath, DataFrame)

    # 筛选当前 variant 对应的 outcome
    outcome_str = variant == :income_share_var ? "income_share" : string(variant)
    df_sub = filter(r -> r.outcome == outcome_str, df)
    nrow(df_sub) == 0 && error("No rows found for outcome=$outcome_str in $fpath")
    

    # 列名 → group id 的映射
    group_col_map = Dict(
        1 => "group1_Managerial",
        2 => "group2_Professional_specialty",
        3 => "group3_High_tech",
        4 => "group4_Sales",
        5 => "group5_Administrative_support",
        6 => "group6_Service",
        7 => "group7_Farming_forestry_construction",
        8 => "group8_Precision_production_repair",
        9 => "group9_Machine_operators_transport",
    )

    occ_irfs = Dict{Int, DataFrame}()
    for (g, colname) in group_col_map
        col = Symbol(colname)
        col in propertynames(df_sub) || (@warn "Column $colname not found"; continue)
        occ_irfs[g] = DataFrame(
            horizon  = df_sub.horizon,
            beta     = Float64.(df_sub[!, col]),
            # 这个文件只有 beta，没有 CI——截面分析只需要 beta
            ci_lo90  = fill(NaN, nrow(df_sub)),
            ci_hi90  = fill(NaN, nrow(df_sub)),
            ci_lo68  = fill(NaN, nrow(df_sub)),
            ci_hi68  = fill(NaN, nrow(df_sub)),
        )
    end
    

    println("Loaded occ IRFs for outcome=$outcome_str: $(length(occ_irfs)) groups")
    return occ_irfs
end

# =============================================================================
# 3. VAR / BIC LAG SELECTION  [M3]
# =============================================================================
function var_estim(y::Matrix{Float64}, p::Int, intercept::Bool)
    T, n_v = size(y)
    T_eff  = T - p
    X_cols = [y[p-l+1:T-l, :] for l in 1:p]
    intercept && push!(X_cols, ones(T_eff, 1))
    X = hcat(X_cols...)
    Y = y[p+1:end, :]
    β = (X'X) \ (X'Y)
    resid = Y - X * β
    Sigma = (resid'resid) ./ T_eff
    return β, Sigma
end

function ic_var(y::Matrix{Float64}, p_max::Int, method::Int=2)
    n_v = size(y, 2)
    T   = size(y, 1)
    T_eff = T - p_max
    bic = zeros(p_max)
    aic = zeros(p_max)
    for p in 1:p_max
        _, Sigma = var_estim(y, p, true)
        n_params = n_v^2 * p + n_v
        bic[p] = log(det(Sigma)) + n_params * log(T_eff) / T_eff
        aic[p] = log(det(Sigma)) + n_params * 2.0 / T_eff
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
# 4. DRISCOLL-KRAAY VCOV
# =============================================================================
function driscoll_kraay_vcov(X::Matrix{Float64}, e::Vector{Float64},
                              t_idx::Vector{Int}; m::Int=DK_BW)::Matrix{Float64}
    N, K    = size(X)
    Tvals   = sort(unique(t_idx))
    T       = length(Tvals)
    tmap    = Dict(v => i for (i,v) in enumerate(Tvals))
    H       = zeros(K, T)
    for n in 1:N
        ti = tmap[t_idx[n]]
        H[:, ti] .+= X[n, :] .* e[n]
    end
    S = zeros(K, K)
    for t in 1:T; S .+= H[:,t] * H[:,t]'; end
    S ./= T
    for l in 1:m
        G = zeros(K, K)
        for t in (l+1):T; G .+= H[:,t] * H[:,t-l]'; end
        G ./= T
        w = 1 - l/(m+1)
        S .+= w .* (G .+ G')
    end
    XX = inv(X'X)
    return XX * (T * S) * XX
end

# =============================================================================
# 5. FE ABSORPTION VIA FIXEDEFFECTMODELS  [M1]
# =============================================================================
function absorb_fe(df::DataFrame, v::Symbol)::Vector{Float64}
    fml = term(v) ~ term(1) + fe(term(:ind_group))
    fit = reg(df, fml, Vcov.simple(), save=:residuals)
    return residuals(fit)
end

# =============================================================================
# 6. BUILD DATASET FOR HORIZON h
# =============================================================================
function build_lp_data(panel::DataFrame, h::Int, y_var::Symbol, l_lag::Int)
    sort!(panel, [:ind_group, :date])
    out = DataFrame[]
    for g in groupby(panel, :ind_group)
        sub = sort(copy(g), :date)
        n   = nrow(sub)
        leady = Vector{Union{Missing,Float64}}(missing, n)
        h < n && (leady[1:n-h] = sub[!, y_var][1+h:n])
        sub[!, :dep_var] = leady
        for l in 1:l_lag
            col    = Symbol("ylag$l")
            lagged = Vector{Union{Missing,Float64}}(missing, n)
            l < n  && (lagged[l+1:n] = sub[!, y_var][1:n-l])
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
# 7. BUILD RHS COLUMNS (group-specific shock interactions)
# =============================================================================
function build_cols!(df::DataFrame, controls::Vector{Symbol}, l_lag::Int)
    shock_cols = Symbol[]
    groups = sort(unique(df.ind_group))
    for g in groups
        c = Symbol("shock_g$(g)")
        df[!, c] = df.shock .* (df.ind_group .== g)
        push!(shock_cols, c)
    end
    lag_cols     = [Symbol("shock_lag$l") for l in 1:l_lag]
    ylag_cols    = [Symbol("ylag$l")      for l in 1:l_lag]
    mac_controls = [:log_oil_lag1, :ffr_lag1, :cpi_lag1, :indpro_lag1, :t10y3m_lag1]
    xcols = vcat(shock_cols, lag_cols, ylag_cols, mac_controls, controls)
    return xcols, shock_cols
end

# =============================================================================
# 8. OLS ON FE-ABSORBED DATA  [M1]
# =============================================================================
function estimate_lp(df_h::DataFrame, panel::DataFrame,
                     controls::Vector{Symbol}, h::Int, l_lag::Int)
    df = copy(df_h)
    xcols, _ = build_cols!(df, controls, l_lag)

    all_dates = sort(unique(panel.date))
    dmap      = Dict(d => i for (i,d) in enumerate(all_dates))
    df[!, :time_idx]  = Int.([dmap[d] for d in df.date])
    df[!, :ind_group] = categorical(df.ind_group)

    needed = vcat([:dep_var, :ind_group, :time_idx], xcols)
    df = dropmissing(df, needed)

    Y_tilde = absorb_fe(df, :dep_var)
    X_tilde = hcat([absorb_fe(df, c) for c in xcols]...)

    β  = (X_tilde'X_tilde) \ (X_tilde'Y_tilde)
    e  = Y_tilde - X_tilde * β
    V  = driscoll_kraay_vcov(X_tilde, e, Int.(df.time_idx); m=DK_BW)
    se = sqrt.(diag(V))

    return (β=β, se=se, e=e, X=X_tilde,
            t_idx=Int.(df.time_idx), df=df, names=xcols)
end

# =============================================================================
# 9. PERCENTILE-T BOOTSTRAP  [M2]
# =============================================================================
function bootstrap_lp(df_h::DataFrame, panel::DataFrame,
                      controls::Vector{Symbol}, h::Int, l_lag::Int;
                      n_boot::Int=N_BOOT,
                      block_size::Int=BLOCK_SIZE,
                      rng=MersenneTwister(BOOT_SEED))

    base     = estimate_lp(df_h, panel, controls, h, l_lag)
    β0, se0  = base.β, base.se
    K        = length(β0)

    all_dates = sort(unique(base.df.date))
    T  = length(all_dates)
    nb = ceil(Int, T / block_size)

    date_rows = Dict{Date,Vector{Int}}()
    for (i, r) in enumerate(eachrow(base.df))
        push!(get!(date_rows, r.date, Int[]), i)
    end

    B     = fill(NaN, n_boot, K)
    TSTAT = fill(NaN, n_boot, K)

    for b in 1:n_boot
        starts = rand(rng, 1:T, nb)
        dd = Date[]
        for s in starts, k in 0:block_size-1
            push!(dd, all_dates[mod1(s+k, T)])
        end
        dd = dd[1:T]
        rows = reduce(vcat, [get(date_rows, d, Int[]) for d in dd])
        isempty(rows) && continue

        df_b = base.df[rows, :]
        # Impose null: dep_var = X*β0 + residual
        df_b[!, :dep_var] = base.X[rows,:] * β0 .+ base.e[rows]

        Y_b = absorb_fe(df_b, :dep_var)
        X_b = hcat([absorb_fe(df_b, c) for c in base.names]...)

        βb = try; (X_b'X_b) \ (X_b'Y_b); catch; continue; end
        eb = Y_b - X_b * βb
        Vb = driscoll_kraay_vcov(X_b, eb, Int.(df_b.time_idx); m=DK_BW)
        seb = sqrt.(diag(Vb))

        B[b,:]     = βb
        TSTAT[b,:] = (βb .- β0) ./ seb
    end

    return B, TSTAT, β0, se0, base.names
end

# =============================================================================
# 10. ONE HORIZON
# =============================================================================
function run_h(panel::DataFrame, y_var::Symbol,
               controls::Vector{Symbol}, h::Int, l_lag::Int)
    df_h = build_lp_data(panel, h, y_var, l_lag)
    B, TSTAT, β0, se0, names = bootstrap_lp(df_h, panel, controls, h, l_lag)
    K    = length(β0)
    rows = NamedTuple[]
    for k in 1:K
        good  = .!isnan.(TSTAT[:,k])
        q95   = quantile(TSTAT[good,k], 0.95)
        q05   = quantile(TSTAT[good,k], 0.05)
        q84   = quantile(TSTAT[good,k], 0.84)
        q16   = quantile(TSTAT[good,k], 0.16)
        push!(rows, (
            horizon  = h,
            coef_name = string(names[k]),
            beta      = β0[k],
            se        = se0[k],
            ci_lo90   = β0[k] - se0[k]*q95,
            ci_hi90   = β0[k] - se0[k]*q05,
            ci_lo68   = β0[k] - se0[k]*q84,
            ci_hi68   = β0[k] - se0[k]*q16,
        ))
    end
    return DataFrame(rows), B, TSTAT, names
end

# =============================================================================
# 11. FULL LP RUN
# =============================================================================
function run_full_lp(panel::DataFrame, y_var::Symbol,
                     controls::Vector{Symbol}, l_lag::Int;
                     h_max::Int=H_MAX)
    out_dfs    = DataFrame[]
    boot_store = Dict{Int,NamedTuple}()
    coef_names = Symbol[]
    for h in 0:h_max
        @info "  h = $h / $h_max"
        df_res, B, TSTAT, names = run_h(panel, y_var, controls, h, l_lag)
        push!(out_dfs, df_res)
        boot_store[h] = (B=B, TSTAT=TSTAT, names=names)
        isempty(coef_names) && (coef_names = names)
    end
    return vcat(out_dfs...), boot_store, coef_names
end

# =============================================================================
# 12. EXTRACT GROUP IRFs
# =============================================================================
function extract_irf(df::DataFrame)
    out = Dict{Int,DataFrame}()
    for g in keys(IND_LABELS)
        cname = "shock_g$(g)"
        sub = filter(r -> r.coef_name == cname, df)
        out[g] = sub[:, [:horizon,:beta,:ci_lo90,:ci_hi90,:ci_lo68,:ci_hi68]]
    end
    return out
end

# =============================================================================
# 13. CROSS-SECTIONAL OIL INTENSITY ANALYSIS  [M5]
#
#   Step 1: 计算每个 ind_group 的累积 IRF (h=1:36)
#   Step 2: 截面回归 cumIRF_g = α + δ·OilIntensity_g + u_g
#   Step 3: 报告 R², δ, variance share
#
#   也接受 occ_group 层面的 oil_exposure（Bartik），结构相同
# =============================================================================
function run_cross_sectional(
        irfs::Dict{Int,DataFrame},
        intensity_df::DataFrame;     # 列: group_id, intensity_value
        group_col::Symbol  = :ind_group,
        intens_col::Symbol = :oil_intensity,
        h_range::UnitRange = 1:36,
        label::String      = "Industry",
        group_labels::Dict = IND_LABELS
    )::DataFrame

    println("\n", "="^50)
    println("Cross-sectional OilIntensity analysis ($(label))")
    println("="^50)

    # 累积 IRF
    cum_irf = Dict{Int,Float64}()
    for (g, df) in irfs
        sub = filter(r -> r.horizon in h_range, df)
        nrow(sub) == 0 && continue
        cum_irf[g] = sum(sub.beta)
    end

    # 合并 intensity
    groups  = sort(collect(keys(cum_irf)))
    irf_vec = Float64[]
    oi_vec  = Float64[]

    for g in groups
        row = filter(r -> r[group_col] == g, intensity_df)
        nrow(row) == 0 && continue
        push!(irf_vec, cum_irf[g])
        push!(oi_vec,  row[1, intens_col])
    end

    n = length(irf_vec)
    println("  Valid groups: $n")

    # OLS: cumIRF = α + δ·OilIntensity
    X_cs  = hcat(ones(n), oi_vec)
    β_cs  = (X_cs'X_cs) \ (X_cs'irf_vec)
    α_hat, δ_hat = β_cs

    resid = irf_vec .- X_cs * β_cs
    ss_tot = sum((irf_vec .- mean(irf_vec)).^2)
    ss_res = sum(resid.^2)
    r2     = 1 - ss_res / ss_tot

    s2      = ss_res / (n - 2)                        # MSE
    se_δ    = sqrt(s2 / sum((oi_vec .- mean(oi_vec)).^2))
    se_α    = sqrt(s2 * (1/n + mean(oi_vec)^2 / sum((oi_vec .- mean(oi_vec)).^2)))
    t_δ     = δ_hat / se_δ
    t_α     = α_hat / se_α
    p_δ     = 2 * ccdf(TDist(n - 2), abs(t_δ))
    p_α     = 2 * ccdf(TDist(n - 2), abs(t_α))

    # Variance share: Var(δ̂·OI) / Var(cumIRF)
    fitted_industry = δ_hat .* oi_vec
    var_share = var(fitted_industry) / var(irf_vec)

    println(@sprintf("  α = %.4f  (se=%.4f, t=%.3f, p=%.3f)", α_hat, se_α, t_α, p_α))
    println(@sprintf("  δ = %.4f  (se=%.4f, t=%.3f, p=%.3f)", δ_hat, se_δ, t_δ, p_δ))
    println(@sprintf("  R² = %.4f", r2))
    println(@sprintf("  Variance share = %.4f", var_share))

    # Print group-level details
    println("\n  Group-level results:")
    println(@sprintf("  %-5s  %-38s  %8s  %8s", "ID", "Label", "cumIRF", "OilInt"))
    for (i, g) in enumerate(groups)
        lbl = get(group_labels, g, "group_$g")
        println(@sprintf("  %-5d  %-38s  %8.4f  %8.4f",
                         g, lbl[1:min(38,length(lbl))], irf_vec[i], oi_vec[i]))
    end

    # Save results
    res_df = DataFrame(
        group_id      = groups,
        cum_irf       = irf_vec,
        oil_intensity = oi_vec,
        fitted        = X_cs * β_cs,
        residual      = resid,
    )
    res_df[!, :alpha]      .= α_hat
    res_df[!, :delta]      .= δ_hat
    res_df[!, :r_squared]  .= r2
    res_df[!, :var_share]  .= var_share

    return res_df
end

# =============================================================================
# 14. MAIN API
# =============================================================================
"""
    run_lp(panel, variant) -> (results_df, boot_store, coef_names, irfs)

主入口，对齐 occ 版本的 02_lp_estimation.jl 方法：
  - BIC lag selection
  - FixedEffectModels FE absorption
  - Percentile-t bootstrap
  - DK standard errors
"""
function run_lp(panel::DataFrame, variant::Symbol)
    y_var, controls = get_lp_specs(variant)

    # [M3] BIC lag selection
    l_lag = select_lag_length(panel; p_max=L_LAG_MAX)
    @info "Running LP for $variant | outcome=$y_var | lags=$l_lag"

    results_df, boot_store, coef_names = run_full_lp(panel, y_var, controls, l_lag)

    output_dir = get_output_dir(variant)
    CSV.write(joinpath(output_dir, "lp_coefficients.csv"), results_df)

    irfs = extract_irf(results_df)
    for (g, d) in irfs
        label = get(IND_LABELS, g, "group_$g")
        CSV.write(joinpath(output_dir, "irf_group$(g)_$(label).csv"), d)
    end

    @info "Saved results to $output_dir"
    return results_df, boot_store, coef_names, irfs
end

"""
    run_cross_sectional_ind(irfs, ind_intensity_df) -> DataFrame

在 ind 层面跑截面 OilIntensity 回归，报告 R² 和 variance share。
ind_intensity_df 来自 data_prep_ind_v2.jl 的 main() 第三个返回值。
"""
function run_cross_sectional_ind(irfs::Dict{Int,DataFrame},
                                  ind_intensity_df::DataFrame;
                                  h_range::UnitRange=1:36)::DataFrame
    # ind_intensity_df 列: ind_group, oil_intensity
    rename_df = rename(ind_intensity_df, :ind_group => :ind_group)
    return run_cross_sectional(irfs, rename_df;
                               group_col  = :ind_group,
                               intens_col = :oil_intensity,
                               h_range    = h_range,
                               label      = "Industry")
end

"""
    run_cross_sectional_occ_from_ind(occ_irfs, oilshare_df) -> DataFrame

在 occ 层面用 Bartik OilShare 跑截面回归。
occ_irfs      : 来自 occ 版本 run_lp 的 irfs（键为 occ_group）
oilshare_df   : 来自 data_prep_ind_v2.jl 的 main() 第二个返回值
"""
const OCC_LABELS_CS = Dict(
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

function run_cross_sectional_occ_from_ind(
        occ_irfs::Dict{Int,DataFrame},
        oilshare_df::DataFrame;
        h_range::UnitRange=1:36)::DataFrame

    df_renamed = rename(oilshare_df, :occ_group => :ind_group, :oil_exposure => :oil_intensity)
    return run_cross_sectional(occ_irfs, df_renamed;
                               group_col  = :ind_group,
                               intens_col = :oil_intensity,
                               h_range    = h_range,
                               label      = "Occupation (Bartik OilShare)",
                               group_labels = OCC_LABELS_CS)
end