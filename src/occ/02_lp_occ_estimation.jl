# =============================================================================
# 02_lp_estimation.jl
# 统一的 Local Projection 估计库 (合并版)
# 调用方式示例: run_lp(panel, :hourly_rate)
# =============================================================================

using CSV, DataFrames, DataFramesMeta
using LinearAlgebra, Statistics
using Dates, Random
using Printf

# =============================================================================
# 0. 全局配置与字典映射
# =============================================================================
# 默认 Bootstrap 与宏观配置 [cite: 2]
const N_BOOT       = 500    
const BLOCK_SIZE   = 3      
const BOOT_SEED    = 42
const CI_LEVELS    = [0.90, 0.95, 0.68]
const H_MAX        = 36
const L_LAG        = 12     # 请根据你的实际滞后阶数调整
const DK_BW        = 4

# 默认的职业分组标签 (如果在外部已经定义，这段可以删除或作备用)
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

# 核心映射字典：根据 variant 自动获取 y_var 和 控制变量
function get_lp_specs(variant::Symbol)
    if variant == :hourly_rate
        return :log_rwage, [:log_rincome, :age_mean, :female_share, :married_share, :hours_mean, :unemp_rate]
    elseif variant == :unemployment
        return :unemp_rate, [:log_rincome, :log_rwage, :age_mean, :female_share, :married_share, :hours_mean]
    elseif variant == :hours
        return :log_hours, [:log_rincome, :log_rwage, :age_mean, :female_share, :married_share, :unemp_rate]
    elseif variant == :income
        return :log_rincome, [:log_rwage, :age_mean, :female_share, :married_share, :hours_mean, :unemp_rate]
    elseif variant == :employment
        return :log_emp_count, [:log_rincome, :log_rwage, :age_mean, :female_share, :married_share, :hours_mean]
    elseif variant == :inequality
        return :log_ratio_7525, [:log_rincome, :log_rwage, :age_mean, :female_share, :married_share, :unemp_rate] 
    elseif variant == :median
        return :log_rincome_p50, [:log_rwage_p50, :age_mean, :female_share, :married_share, :hours_mean, :unemp_rate] 
    elseif variant == :income_share_var
        return :income_share, [:age_mean, :female_share, :married_share, :hours_mean, :unemp_rate]
    else
        error("LP Specs not defined for variant: $variant")
    end
end

# =============================================================================
# 1. DRISCOLL-KRAAY STANDARD ERRORS
# =============================================================================
function driscoll_kraay_vcov(X::Matrix{Float64}, e::Vector{Float64},
                              t_idx::Vector{Int}; m::Int = 0)::Matrix{Float64}
    N_obs, K = size(X)
    T_vals   = sort(unique(t_idx))
    T        = length(T_vals)
    m        = m == 0 ? floor(Int, T^(1/4)) : m # [cite: 4, 5]
    t_map    = Dict(v => i for (i, v) in enumerate(T_vals))

    H = zeros(K, T)
    for obs in 1:N_obs
        ti = t_map[t_idx[obs]]
        H[:, ti] .+= X[obs, :] .* e[obs]
    end

    S = zeros(K, K)
    for t in 1:T; S .+= H[:, t] * H[:, t]'; end # [cite: 6]
    S ./= T

    for l in 1:m
        Γl = zeros(K, K)
        for t in (l+1):T; Γl .+= H[:, t] * H[:, t-l]'; end # [cite: 7]
        Γl ./= T
        w   = 1.0 - l / (m + 1)
        S .+= w .* (Γl .+ Γl')
    end

    XtX_inv = inv(X' * X)
    return XtX_inv * (T .* S) * XtX_inv
end

# =============================================================================
# 2. WITHIN (DEMEANING) TRANSFORMATION
# =============================================================================
function within_transform(df::DataFrame,
                           yvars::Vector{Symbol},
                           xvars::Vector{Symbol},
                           id_var::Symbol)::DataFrame
    df_out = copy(df)
    for var in vcat(yvars, xvars) # [cite: 8]
        gdf   = groupby(df_out, id_var)
        means = combine(gdf, var => mean => Symbol(var, "_mean")) # [cite: 9]
        df_out = leftjoin(df_out, means, on = id_var)
        df_out[!, var] = df_out[!, var] .- df_out[!, Symbol(var, "_mean")]
        select!(df_out, Not(Symbol(var, "_mean")))
    end
    return df_out
end

# =============================================================================
# 3. BUILD REGRESSION DATA FOR HORIZON h (已动态化)
# =============================================================================
function build_lp_data(panel::DataFrame, h::Int, y_var::Symbol)::Union{DataFrame, Nothing}
    sort!(panel, [:occ_group, :date])
    gdf    = groupby(panel, :occ_group)
    chunks = DataFrame[]

    for g in gdf
        sub = copy(g)
        n   = nrow(sub)

        # 动态替换原本的 lead_hourlyrate
        lead_y = Vector{Union{Float64,Missing}}(missing, n)
        h < n && (lead_y[1:n-h] = sub[!, y_var][h+1:n])

        # 动态替换原本的 lag_hourlyrate
        lag_y = Vector{Union{Float64,Missing}}(missing, n)
        lag_y[2:n] = sub[!, y_var][1:n-1]

        sub[!, :dep_var] = lead_y .- lag_y
        sub[!, Symbol("lag_", y_var)] = lag_y
     
        push!(chunks, sub)
    end

    df_h = vcat(chunks...)
    required = vcat([:dep_var, :shock, :log_oil_lag1, :ffr_lag1, :cpi_lag1, :indpro_lag1],
                    [Symbol("shock_lag", l) for l in 1:L_LAG])
    df_h = dropmissing(df_h, required)
    return nrow(df_h) == 0 ? nothing : df_h # [cite: 11, 12]
end

# =============================================================================
# 4. BUILD COLUMN LISTS (已动态化)
# =============================================================================
function build_col_lists!(df_h::DataFrame, cell_cols::Vector{Symbol})
    inter_cols = Symbol[]
    for g in 2:9
        col = Symbol("shock_x_occ", g)
        df_h[!, col] = df_h.shock .* Float64.(df_h.occ_group .== g)
        push!(inter_cols, col)
    end
    lag_cols   = [Symbol("shock_lag", l) for l in 1:L_LAG]
    macro_cols = [:log_oil_lag1, :ffr_lag1, :cpi_lag1, :indpro_lag1]
    
    # 动态插入各个变体特有的 cell_cols 控制变量
    x_cols     = vcat([:shock], inter_cols, lag_cols, macro_cols, cell_cols) 
    return x_cols, inter_cols
end

# =============================================================================
# 5. OLS ON WITHIN-TRANSFORMED DATA
# =============================================================================
function ols_within(df_h::DataFrame, y_col::Symbol,
                     x_cols::Vector{Symbol}, panel::DataFrame)
    df_w = within_transform(df_h, [y_col], x_cols, :occ_group)
    df_w = dropmissing(df_w, vcat([y_col], x_cols))
    nrow(df_w) == 0 && return nothing

    Y = Float64.(df_w[!, y_col])
    X = hcat(ones(nrow(df_w)), Matrix{Float64}(df_w[!, x_cols]))
    β = (X' * X) \ (X' * Y) # [cite: 14]
    e = Y .- X * β

    all_dates = sort(unique(panel.date))
    date_map  = Dict(d => i for (i, d) in enumerate(all_dates))
    t_idx     = [date_map[d] for d in df_w.date]

    return (β = β, e = e, X = X, t_idx = t_idx, df_w = df_w)
end

# =============================================================================
# 6. BLOCK BOOTSTRAP (已动态化)
# =============================================================================
function block_bootstrap_lp(panel::DataFrame, h::Int, y_var::Symbol, cell_cols::Vector{Symbol};
                              n_boot::Int     = N_BOOT,
                              block_size::Int = BLOCK_SIZE,
                              rng::AbstractRNG = Random.default_rng()) # [cite: 18, 19]
@time begin
    # ── baseline data ───────────────────────────────────────────────────────
    df_h = build_lp_data(panel, h, y_var)
    isnothing(df_h) && return nothing

    x_cols, _ = build_col_lists!(df_h, cell_cols)
    y_col     = :dep_var

    df_w = within_transform(df_h, [y_col], x_cols, :occ_group)
    df_w = dropmissing(df_w, vcat([y_col], x_cols))
    nrow(df_w) == 0 && return nothing

    Y_base = Float64.(df_w[!, y_col])
    X_base = hcat(ones(nrow(df_w)), Matrix{Float64}(df_w[!, x_cols]))
    β_base = (X_base' * X_base) \ (X_base' * Y_base)
    e_base = Y_base .- X_base * β_base # [cite: 19, 20]

    K          = length(β_base)
    coef_names = vcat([:intercept], x_cols)

    # ── time structure ──────────────────────────────────────────────────────
    all_dates  = sort(unique(df_w.date))
    T          = length(all_dates)
    n_blocks   = ceil(Int, T / block_size)

    date_to_rows = Dict{Date, Vector{Int}}()
    for (i, row) in enumerate(eachrow(df_h))
        push!(get!(date_to_rows, row.date, Int[]), i) # [cite: 21]
    end
end
    # ── bootstrap replications ──────────────────────────────────────────────
    boot_matrix = fill(NaN, n_boot, K)
@time begin
    for b in 1:n_boot
        starts     = rand(rng, 1:T, n_blocks)
        boot_dates = Date[]
        for s in starts
            for k in 0:(block_size - 1)
                push!(boot_dates, all_dates[mod1(s + k, T)]) # [cite: 22]
            end
        end
        boot_dates = boot_dates[1:T]

        new_rows = Int[]
        for d in boot_dates
            append!(new_rows, get(date_to_rows, d, Int[])) # [cite: 22, 23]
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
    valid      = [!any(isnan.(boot_matrix[b, :])) for b in 1:n_boot] # [cite: 24]
    boot_valid = boot_matrix[valid, :]
    n_valid    = sum(valid)

    n_valid < 50 &&
        @warn "h=$h: only $n_valid valid bootstrap draws; CIs may be unreliable." # [cite: 24, 25]

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
    date_map        = Dict(d => i for (i, d) in enumerate(all_dates_panel)) # [cite: 26]
    t_idx           = [date_map[d] for d in df_w.date]
    V               = driscoll_kraay_vcov(X_base, e_base, t_idx; m = DK_BW)
    dk_se           = sqrt.(diag(V))

    # ── results DataFrame ───────────────────────────────────────────────────
    res = DataFrame(
        horizon      = h, # [cite: 27]
        coef_name    = coef_names,
        beta         = β_base,
        dk_se        = dk_se,
        boot_se      = [std(boot_valid[:, k]) for k in 1:K],
        ci_lo95      = ci_lo[0.95],
        ci_hi95      = ci_hi[0.95],
        ci_lo90      = ci_lo[0.90],
        ci_hi90      = ci_hi[0.90],
        ci_lo68      = ci_lo[0.68],
        ci_hi68      = ci_hi[0.68],
        n_boot_valid = n_valid,
    )

    return (results = res, boot_draws = boot_valid, coef_names = coef_names) # [cite: 28]
end

# =============================================================================
# 7. RUN FULL LP (已动态化)
# =============================================================================
function run_full_lp(panel::DataFrame, y_var::Symbol, cell_cols::Vector{Symbol};
                      h_max::Int      = H_MAX,
                      n_boot::Int     = N_BOOT,
                      block_size::Int = BLOCK_SIZE,
                      seed::Int       = BOOT_SEED
                      )::Tuple{DataFrame, Dict{Int,Matrix{Float64}}, Vector{Symbol}} # [cite: 30, 31]

    println("\nRunning LP  h = 0 … $h_max for Variable: $y_var")
    println("Block bootstrap: $n_boot replications, block size = $block_size months")

    rng        = MersenneTwister(seed)
    all_res    = DataFrame[]
    boot_store = Dict{Int, Matrix{Float64}}()
    coef_names_ref = Symbol[]

    for h in 0:h_max
        print("\r  h = $h / $h_max  ")
        @time out = block_bootstrap_lp(panel, h, y_var, cell_cols;
                                  n_boot     = n_boot,
                                  block_size = block_size,
                                  rng        = rng) # [cite: 32, 33]
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
function extract_irf(results_df::DataFrame,
                      boot_store::Dict{Int, Matrix{Float64}},
                      coef_names::Vector{Symbol})::Dict{Int, DataFrame} # [cite: 35]

    irfs     = Dict{Int, DataFrame}()
    horizons = sort(unique(results_df.horizon))

    β_idx = findfirst(==(:shock), coef_names)
    θ_idx = Dict(g => findfirst(==(Symbol("shock_x_occ", g)), coef_names) for g in 2:9) # [cite: 35, 36]

    # ── Group 1: baseline ───────────────────────────────────────────────────
    rows1 = []
    for h in horizons
        sub   = filter(r -> r.horizon == h && r.coef_name == :shock, results_df)
        nrow(sub) == 0 && continue
        β_h = sub.beta[1]

        if haskey(boot_store, h)
            bd = boot_store[h][:, β_idx] # [cite: 37]
            push!(rows1, (horizon  = h,
                          beta_abs = β_h,
                          boot_se  = std(bd),
                          ci_lo95  = quantile(bd, 0.025), # [cite: 38]
                          ci_hi95  = quantile(bd, 0.975),
                          ci_lo90  = quantile(bd, 0.05),
                          ci_hi90  = quantile(bd, 0.95),
                          ci_lo68  = quantile(bd, 0.16),
                          ci_hi68  = quantile(bd, 0.84),
                          theta            = 0.0,
                          boot_se_theta    = NaN,
                          ci_lo95_theta    = NaN,
                          ci_hi95_theta    = NaN,
                          ci_lo90_theta    = NaN,
                          ci_hi90_theta    = NaN,
                          ci_lo68_theta    = NaN,
                          ci_hi68_theta    = NaN)) # [cite: 40, 41, 42]
        end
    end
    irfs[1] = DataFrame(rows1)

    # ── Groups 2-9 ──────────────────────────────────────────────────────────
    for g in 2:9
        ti = get(θ_idx, g, nothing)
        isnothing(ti) && continue
        rows_g = []

        for h in horizons
            β_sub = filter(r -> r.horizon == h && r.coef_name == :shock, results_df) # [cite: 43]
            θ_sub = filter(r -> r.horizon == h && r.coef_name == Symbol("shock_x_occ", g), results_df)
            (nrow(β_sub) == 0 || nrow(θ_sub) == 0) && continue # [cite: 44]

            β_h   = β_sub.beta[1]
            θ_h   = θ_sub.beta[1]
            abs_h = β_h + θ_h

            if haskey(boot_store, h)
                B        = boot_store[h]
                abs_boot = B[:, β_idx] .+ B[:, ti] # [cite: 45]
                θ_boot   = B[:, ti]

                push!(rows_g, (
                    horizon          = h,
                    beta_abs         = abs_h, # [cite: 46]
                    boot_se          = std(abs_boot),
                    ci_lo95          = quantile(abs_boot, 0.025),
                    ci_hi95          = quantile(abs_boot, 0.975), # [cite: 47]
                    ci_lo90          = quantile(abs_boot, 0.05),
                    ci_hi90          = quantile(abs_boot, 0.95),
                    ci_lo68          = quantile(abs_boot, 0.16), # [cite: 48]
                    ci_hi68          = quantile(abs_boot, 0.84),
                    theta            = θ_h,
                    boot_se_theta    = std(θ_boot),
                    ci_lo95_theta    = quantile(θ_boot, 0.025), # [cite: 49]
                    ci_hi95_theta    = quantile(θ_boot, 0.975),
                    ci_lo90_theta    = quantile(θ_boot, 0.05),
                    ci_hi90_theta    = quantile(θ_boot, 0.95),
                    ci_lo68_theta    = quantile(θ_boot, 0.16), # [cite: 50]
                    ci_hi68_theta    = quantile(θ_boot, 0.84),
                ))
            end
        end
        irfs[g] = DataFrame(rows_g)
    end

    return irfs # [cite: 51]
end

# =============================================================================
# 9. MAIN RUNNER (最终封装给外部调用的极简入口)
# =============================================================================
function run_lp(panel::DataFrame, variant::Symbol)
    println("="^60)
    println("Running Local Projection for variant: $variant")
    println("="^60)

    # 1. 自动获取被解释变量与控制变量
    y_var, cell_cols = get_lp_specs(variant)
    println("Dependent Variable: $y_var")
    println("Controls: $cell_cols")

    # 2. 动态创建此变体的输出文件夹，类似你原本的机制
    output_dir = get_output_dir(variant)

    # 3. 运行本地投影与 Bootstrap 过程
    results_df, boot_store, coef_names = run_full_lp(panel, y_var, cell_cols)

    # 4. 保存回归系数主表
    CSV.write(joinpath(output_dir, "lp_coefficients.csv"), results_df)
    println("Coefficients saved to $output_dir/lp_coefficients.csv")

    # 5. 提取并保存各分组的 IRF 路径表
    irfs = extract_irf(results_df, boot_store, coef_names)
    for (g, df) in irfs
        label = get(OCC_LABELS, g, "group_$g")
        CSV.write(joinpath(output_dir, "irf_group$(g)_$(label).csv"), df)
    end
    println("IRF tables saved.")

    return results_df, boot_store, coef_names, irfs
end