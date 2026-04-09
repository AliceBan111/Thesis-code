# =============================================================================
# 02_lp_estimation_unified.jl
# 统一的 Local Projection 估计库 (无基准组合并版 / 方案A)
# 调用方式示例: run_lp(panel, :hourly_rate)
# =============================================================================

using CSV, DataFrames, DataFramesMeta
using LinearAlgebra, Statistics
using Dates, Random
using Printf

# =============================================================================
# 0. 全局配置与字典映射
# =============================================================================
# Bootstrap 与宏观配置
const N_BOOT       = 500    
const BLOCK_SIZE   = 3      
const BOOT_SEED    = 42
const CI_LEVELS    = [0.90, 0.95, 0.68]
const H_MAX        = 36
const L_LAG        = 12     
const DK_BW        = 4

# 默认的行业分组标签
const IND_LABELS = Dict(
    1  => "Agriculture_forestry_fishing",  2  => "Mining",
    3  => "Construction",                  4  => "Manufacturing_nondurable",
    5  => "Manufacturing_durable",         6  => "Transportation_utilities",
    7  => "Wholesale_trade",               8  => "Retail_trade",
    9  => "Finance_insurance_realestate",  10 => "Business_repair_services",
    11 => "Personal_entertainment_services", 12 => "Professional_related_services",
    13 => "Public_administration",
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
# 2. WITHIN (DEMEANING) TRANSFORMATION
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
function build_lp_data(panel::DataFrame, h::Int, y_var::Symbol)::Union{DataFrame, Nothing}
    sort!(panel, [:ind_group, :date])
    gdf    = groupby(panel, :ind_group)
    chunks = DataFrame[]

    for g in gdf
        sub = copy(g)
        n   = nrow(sub)

        lead_y = Vector{Union{Float64,Missing}}(missing, n)
        h < n && (lead_y[1:n-h] = sub[!, y_var][h+1:n])

        lag_y = Vector{Union{Float64,Missing}}(missing, n)
        lag_y[2:n] = sub[!, y_var][1:n-1]

        sub[!, :dep_var] = lead_y .- lag_y
        sub[!, Symbol("lag_", y_var)] = lag_y
      
        push!(chunks, sub)
    end

    df_h = vcat(chunks...)
    required = vcat([:dep_var, :shock, :log_oil_lag1, :ffr_lag1, :cpi_lag1, :indpro_lag1],
                    [Symbol("shock_lag", l) for l in 1:L_LAG])
    
    # 因为存在部分列可能在数据源中没有被使用，我们可以采用子集交集判断以防报错
    actual_required = intersect(required, propertynames(df_h))
    df_h = dropmissing(df_h, actual_required)
    
    return nrow(df_h) == 0 ? nothing : df_h 
end

# =============================================================================
# 4. BUILD COLUMN LISTS (方案A：完全无基准组)
# =============================================================================
function build_col_lists!(df_h::DataFrame, cell_cols::Vector{Symbol})
    inter_cols = Symbol[]
    # 【改动点1】遍历所有的组 (1 到 13)，为每一个组直接生成各自的 shock 交互项
    group_ids = sort(collect(keys(IND_LABELS)))
    for g in group_ids
        col = Symbol("shock_x_ind", g)
        df_h[!, col] = df_h.shock .* Float64.(df_h.ind_group .== g)
        push!(inter_cols, col)
    end
    lag_cols   = [Symbol("shock_lag", l) for l in 1:L_LAG]
    macro_cols = [:log_oil_lag1, :ffr_lag1, :cpi_lag1, :indpro_lag1] 
    macro_cols = intersect(macro_cols, propertynames(df_h)) # 防止某些文件没这些宏观变量报错
    
    # 【改动点2】不再拼接单个的 [:shock] 变量
    x_cols     = vcat(inter_cols, lag_cols, macro_cols, cell_cols) 
    return x_cols, inter_cols
end

# =============================================================================
# 5. OLS ON WITHIN-TRANSFORMED DATA
# =============================================================================
function ols_within(df_h::DataFrame, y_col::Symbol,
                     x_cols::Vector{Symbol}, panel::DataFrame)
    df_w = within_transform(df_h, [y_col], x_cols, :ind_group)
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
# 6. BLOCK BOOTSTRAP 
# =============================================================================
function block_bootstrap_lp(panel::DataFrame, h::Int, y_var::Symbol, cell_cols::Vector{Symbol};
                              n_boot::Int     = N_BOOT,
                              block_size::Int = BLOCK_SIZE,
                              rng::AbstractRNG = Random.default_rng()) 
    # ── baseline data ───────────────────────────────────────────────────────
    df_h = build_lp_data(panel, h, y_var)
    isnothing(df_h) && return nothing

    x_cols, _ = build_col_lists!(df_h, cell_cols)
    y_col     = :dep_var

    df_w = within_transform(df_h, [y_col], x_cols, :ind_group)
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

    date_to_rows = Dict{Date, Vector{Int}}()
    for (i, row) in enumerate(eachrow(df_h))
        push!(get!(date_to_rows, row.date, Int[]), i) 
    end

    # ── bootstrap replications ──────────────────────────────────────────────
    boot_matrix = fill(NaN, n_boot, K)
    for b in 1:n_boot
        starts     = rand(rng, 1:T, n_blocks)
        boot_dates = Date[]
        for s in starts
            for k in 0:(block_size - 1)
                push!(boot_dates, all_dates[mod1(s + k, T)]) 
            end
        end
        boot_dates = boot_dates[1:T]

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
        ci_lo68      = ci_lo[0.68],
        ci_hi68      = ci_hi[0.68],
        n_boot_valid = n_valid,
    )

    return (results = res, boot_draws = boot_valid, coef_names = coef_names) 
end

# =============================================================================
# 7. RUN FULL LP
# =============================================================================
function run_full_lp(panel::DataFrame, y_var::Symbol, cell_cols::Vector{Symbol};
                      h_max::Int      = H_MAX,
                      n_boot::Int     = N_BOOT,
                      block_size::Int = BLOCK_SIZE,
                      seed::Int       = BOOT_SEED
                      )::Tuple{DataFrame, Dict{Int,Matrix{Float64}}, Vector{Symbol}} 

    println("\nRunning LP  h = 0 … $h_max for Variable: $y_var")
    println("Block bootstrap: $n_boot replications, block size = $block_size months")

    rng        = MersenneTwister(seed)
    all_res    = DataFrame[]
    boot_store = Dict{Int, Matrix{Float64}}()
    coef_names_ref = Symbol[]

    for h in 0:h_max
        print("\r  h = $h / $h_max  ")
        out = block_bootstrap_lp(panel, h, y_var, cell_cols;
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
# 8. EXTRACT IRF PATHS (方案A：提取每个组绝对效应)
# =============================================================================
function extract_irf(results_df::DataFrame,
                      boot_store::Dict{Int, Matrix{Float64}},
                      coef_names::Vector{Symbol})::Dict{Int, DataFrame} 

    irfs     = Dict{Int, DataFrame}()
    horizons = sort(unique(results_df.horizon))
    group_ids = sort(collect(keys(IND_LABELS)))

    # 【改动点3】索引字典直接寻找各自的 shock_x_ind_g
    β_indices = Dict(g => findfirst(==(Symbol("shock_x_ind", g)), coef_names) for g in group_ids) 

    # 【改动点4】使用统一循环处理所有组（1-13），不再需要把 baseline 的项加在一起
    for g in group_ids
        ti = get(β_indices, g, nothing)
        isnothing(ti) && continue
        rows_g = []

        for h in horizons
            # 直接找 shock_x_ind_g 的系数
            coef_sub = filter(r -> r.horizon == h && r.coef_name == Symbol("shock_x_ind", g), results_df)
            nrow(coef_sub) == 0 && continue

            beta_h = coef_sub.beta[1]

            if haskey(boot_store, h)
                B = boot_store[h]
                beta_boot = B[:, ti]   # 直接取该行业的 bootstrap draws

                push!(rows_g, (
                    horizon          = h,
                    beta_abs         = beta_h,
                    boot_se          = std(beta_boot),
                    ci_lo95          = quantile(beta_boot, 0.025),
                    ci_hi95          = quantile(beta_boot, 0.975),
                    ci_lo90          = quantile(beta_boot, 0.05),
                    ci_hi90          = quantile(beta_boot, 0.95),
                    ci_lo68          = quantile(beta_boot, 0.16),
                    ci_hi68          = quantile(beta_boot, 0.84),
                    # 方案 A 下不再有相对差异，theta 字段设为 NaN
                    # theta            = NaN,
                    # boot_se_theta    = NaN,
                    # ci_lo95_theta    = NaN,
                    # ci_hi95_theta    = NaN,
                    # ci_lo90_theta    = NaN,
                    # ci_hi90_theta    = NaN,
                    # ci_lo68_theta    = NaN,
                    # ci_hi68_theta    = NaN,
                ))
            end
        end
        irfs[g] = DataFrame(rows_g)
    end

    return irfs 
end

# =============================================================================
# 9. MAIN RUNNER 
# =============================================================================
function run_lp(panel::DataFrame, variant::Symbol)
    println("="^60)
    println("Running Local Projection for variant: $variant")
    println("="^60)

    # 1. 自动获取被解释变量与控制变量
    y_var, cell_cols = get_lp_specs(variant)
    println("Dependent Variable: $y_var")
    println("Controls: $cell_cols")

    # 2. 动态创建此变体的输出文件夹
    output_dir = get_output_dir_ind(variant)

    # 3. 运行本地投影与 Bootstrap 过程
    results_df, boot_store, coef_names = run_full_lp(panel, y_var, cell_cols)

    # 4. 保存回归系数主表
    CSV.write(joinpath(output_dir, "lp_coefficients.csv"), results_df)
    println("Coefficients saved to $output_dir/lp_coefficients.csv")

    # 5. 提取并保存各分组的 IRF 路径表
    irfs = extract_irf(results_df, boot_store, coef_names)
    for (g, df) in irfs
        label = get(IND_LABELS, g, "group_$g")
        CSV.write(joinpath(output_dir, "irf_group$(g)_$(label).csv"), df)
    end
    println("IRF tables saved.")

    return results_df, boot_store, coef_names, irfs
end