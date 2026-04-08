using CSV, DataFrames, DataFramesMeta
using XLSX, Dates, Statistics
using LinearAlgebra, Random
using Printf

# =============================================================================
# 0. PATHS & SETTINGS
# =============================================================================
const DATA_DIR   = joinpath(@__DIR__, "../../data")
const OUTPUT_DIR = joinpath(@__DIR__, "../../result/macro/cpi")
mkpath(OUTPUT_DIR)

const DATE_START = Date(1983, 4, 1)
const DATE_END   = Date(2025, 6, 1)

const H_MAX      = 36
const L_LAG      = 12
const N_BOOT     = 500
const BLOCK_SIZE = 12
const BOOT_SEED  = 42
const CI_LEVELS  = [0.68, 0.90, 0.95]
const DK_BW      = 4

# =============================================================================
# 1. DATE PARSING UTILITY  (same as data_prep_occ_income.jl)
# =============================================================================
function parse_yearmonth(raw)::Union{Date, Nothing}
    raw isa Date && return raw
    s = strip(string(raw))
    m = match(r"^(\d{4})[Mm](\d{2})$", s)
    isnothing(m) && return nothing
    yr = parse(Int, m.captures[1])
    mo = parse(Int, m.captures[2])
    return Date(yr, mo, 1)
end

# =============================================================================
# 2. LOAD MACRO DATA
# =============================================================================
function load_cpi_macro()::DataFrame
    println("Loading CPI from VARdata.xlsx...")
    path = joinpath(DATA_DIR, "VARdata.xlsx")
    xf   = XLSX.readxlsx(path)
    sh   = xf["Monthly"]

    dates    = Vector{Date}()
    cpi_vals = Vector{Float64}()

    for row in XLSX.eachrow(sh)
        XLSX.row_number(row) == 1 && continue
        d_raw = row[1]
        c_raw = row[7]   # Column G = CPI
        (ismissing(d_raw) || ismissing(c_raw)) && continue
        d = parse_yearmonth(d_raw)
        isnothing(d) && continue
        push!(dates,    d)
        push!(cpi_vals, Float64(c_raw))
    end

    df = DataFrame(date = dates, cpi = cpi_vals)
    df[!, :log_cpi] = log.(df.cpi)
    df = @subset(df, DATE_START .<= :date .<= DATE_END)
    sort!(df, :date)
    println("  CPI rows: ", nrow(df))
    return df
end

function load_oil_shock_macro()::DataFrame
    println("Loading oil supply news shock...")
    path = joinpath(DATA_DIR, "oilSupplyNewsShocks_2025M06.xlsx")
    xf   = XLSX.readxlsx(path)
    sh   = xf["Monthly"]

    dates  = Vector{Date}()
    shocks = Vector{Float64}()

    for row in XLSX.eachrow(sh)
        XLSX.row_number(row) == 1 && continue
        d_raw = row[1]
        s_raw = row[3]
        (ismissing(d_raw) || ismissing(s_raw)) && continue
        d = parse_yearmonth(d_raw)
        isnothing(d) && continue
        push!(dates,  d)
        push!(shocks, Float64(s_raw))
    end

    df = DataFrame(date = dates, shock = shocks)
    df = @subset(df, DATE_START .<= :date .<= DATE_END)
    sort!(df, :date)
    println("  Shock rows: ", nrow(df))
    return df
end

function load_oil_price_macro()::DataFrame
    println("Loading oil price (WTI) from VARdata.xlsx...")
    path = joinpath(DATA_DIR, "VARdata.xlsx")
    xf   = XLSX.readxlsx(path)
    sh   = xf["Monthly"]

    dates  = Vector{Date}()
    prices = Vector{Float64}()

    for row in XLSX.eachrow(sh)
        XLSX.row_number(row) == 1 && continue
        d_raw = row[1]
        p_raw = row[2]
        (ismissing(d_raw) || ismissing(p_raw)) && continue
        d = parse_yearmonth(d_raw)
        isnothing(d) && continue
        push!(dates,  d)
        push!(prices, Float64(p_raw))
    end

    df = DataFrame(date = dates, oil_price = prices)
    df[!, :log_oil_price] = log.(df.oil_price)
    df = @subset(df, DATE_START .<= :date .<= DATE_END)
    sort!(df, :date)
    return df
end

function load_fred_macro(filename::String, colname::Symbol)::DataFrame
    println("Loading $filename...")
    path = joinpath(DATA_DIR, filename)
    df = CSV.read(path, DataFrame;
        types=Dict(1 => Date, 2 => Float32), dateformat="yyyy-m-d")
    rename!(df, first(names(df)) => :date, names(df)[2] => colname)
    df = filter(row -> DATE_START <= row.date <= DATE_END, df)
    sort!(df, :date)
    return df
end

# =============================================================================
# 3. HELPER: lag vector
# =============================================================================
function lag_vec(v::AbstractVector, k::Int)
    T = Union{eltype(v), Missing}
    out = Vector{T}(missing, length(v))
    out[k+1:end] = v[1:end-k]
    return out
end

# =============================================================================
# 4. BUILD MACRO PANEL
# =============================================================================
function build_cpi_panel()::DataFrame
    println("Building CPI macro panel...")

    cpi_df   = load_cpi_macro()
    shock_df = load_oil_shock_macro()
    oil_df   = load_oil_price_macro()
    ffr_df   = load_fred_macro("FEDFUNDS.csv", :fedfunds)

    panel = outerjoin(cpi_df, shock_df, on = :date)
    panel = outerjoin(panel, oil_df,    on = :date)
    panel = outerjoin(panel, ffr_df,    on = :date)
    sort!(panel, :date)

    # Shock lags
    for l in 1:L_LAG
        panel[!, Symbol("shock_lag", l)] = lag_vec(panel.shock, l)
    end

    # Control lags
    panel[!, :log_oil_lag1] = lag_vec(panel.log_oil_price, 1)
    panel[!, :ffr_lag1]     = lag_vec(panel.fedfunds, 1)
    panel[!, :log_cpi_lag1] = lag_vec(panel.log_cpi, 1)

    panel = dropmissing(panel, vcat([:log_cpi, :shock, :log_oil_lag1, :ffr_lag1],
                                    [Symbol("shock_lag", l) for l in 1:L_LAG]))
    sort!(panel, :date)
    println("  Panel rows: ", nrow(panel))
    return panel
end

# =============================================================================
# 5. BUILD LP DATA FOR HORIZON h
#    dep_var_t = log_cpi_{t+h} - log_cpi_{t-1}
# =============================================================================
function build_lp_data_cpi(panel::DataFrame, h::Int)::Union{DataFrame, Nothing}
    df = sort(copy(panel), :date)
    n  = nrow(df)

    lead_cpi = Vector{Union{Float64, Missing}}(missing, n)
    h < n && (lead_cpi[1:n-h] = df.log_cpi[h+1:n])

    lag_cpi = Vector{Union{Float64, Missing}}(missing, n)
    lag_cpi[2:n] = df.log_cpi[1:n-1]

    df[!, :dep_var]  = lead_cpi .- lag_cpi   # log_cpi_{t+h} - log_cpi_{t-1}
    df[!, :lag_cpi]  = lag_cpi

    required = vcat([:dep_var, :shock, :log_oil_lag1, :ffr_lag1, :lag_cpi],
                    [Symbol("shock_lag", l) for l in 1:L_LAG])
    df = dropmissing(df, required)
    return nrow(df) == 0 ? nothing : df
end

# =============================================================================
# 6. OLS (time-series, no FE needed for aggregate panel)
# =============================================================================
function build_x_cols()::Vector{Symbol}
    lag_cols   = [Symbol("shock_lag", l) for l in 1:L_LAG]
    macro_cols = [:log_oil_lag1, :ffr_lag1, :lag_cpi]
    return vcat([:shock], lag_cols, macro_cols)
end

# =============================================================================
# 7. DRISCOLL-KRAAY SE
# =============================================================================
function driscoll_kraay_vcov_ts(X::Matrix{Float64}, e::Vector{Float64};
                                 m::Int = DK_BW)::Matrix{Float64}
    T, K = size(X)
    S    = zeros(K, K)
    for t in 1:T; S .+= (X[t, :] .* e[t]) * (X[t, :] .* e[t])'; end
    S ./= T

    for l in 1:m
        Γl = zeros(K, K)
        for t in (l+1):T
            Γl .+= (X[t, :] .* e[t]) * (X[t-l, :] .* e[t-l])'
        end
        Γl ./= T
        w   = 1.0 - l / (m + 1)
        S .+= w .* (Γl .+ Γl')
    end

    XtX_inv = inv(X' * X)
    return XtX_inv * (T .* S) * XtX_inv
end

# =============================================================================
# 8. BLOCK BOOTSTRAP (time-series version — no cross-section)
# =============================================================================
function block_bootstrap_cpi(panel::DataFrame, h::Int;
                               n_boot::Int     = N_BOOT,
                               block_size::Int = BLOCK_SIZE,
                               rng::AbstractRNG = Random.default_rng())

    df_h = build_lp_data_cpi(panel, h)
    isnothing(df_h) && return nothing

    x_cols = build_x_cols()
    y_col  = :dep_var

    df_h = dropmissing(df_h, vcat([y_col], x_cols))
    nrow(df_h) == 0 && return nothing

    Y_base = Float64.(df_h[!, y_col])
    X_base = hcat(ones(nrow(df_h)), Matrix{Float64}(df_h[!, x_cols]))
    β_base = (X_base' * X_base) \ (X_base' * Y_base)
    e_base = Y_base .- X_base * β_base

    K          = length(β_base)
    coef_names = vcat([:intercept], x_cols)

    T        = nrow(df_h)
    n_blocks = ceil(Int, T / block_size)

    # ── bootstrap replications ──────────────────────────────────────────────
    boot_matrix = fill(NaN, n_boot, K)
    for b in 1:n_boot
        starts     = rand(rng, 1:T, n_blocks)
        boot_idx   = Int[]
        for s in starts
            for k in 0:(block_size - 1)
                push!(boot_idx, mod1(s + k, T))
            end
        end
        boot_idx = boot_idx[1:T]

        Y_boot = Y_base[boot_idx]
        X_boot = X_base[boot_idx, :]

        β_boot = try
            (X_boot' * X_boot) \ (X_boot' * Y_boot)
        catch
            fill(NaN, K)
        end
        any(isnan, β_boot) && continue
        boot_matrix[b, :] = β_boot
    end

    valid      = [!any(isnan.(boot_matrix[b, :])) for b in 1:n_boot]
    boot_valid = boot_matrix[valid, :]
    n_valid    = sum(valid)
    n_valid < 50 && @warn "h=$h: only $n_valid valid bootstrap draws."

    # ── percentile CIs ──────────────────────────────────────────────────────
    ci_lo = Dict{Float64, Vector{Float64}}()
    ci_hi = Dict{Float64, Vector{Float64}}()
    for lvl in CI_LEVELS
        α = 1.0 - lvl
        ci_lo[lvl] = [quantile(boot_valid[:, k], α / 2)       for k in 1:K]
        ci_hi[lvl] = [quantile(boot_valid[:, k], 1.0 - α / 2) for k in 1:K]
    end

    # ── DK standard errors ──────────────────────────────────────────────────
    V     = driscoll_kraay_vcov_ts(X_base, e_base; m = DK_BW)
    dk_se = sqrt.(diag(V))

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
# 9. RUN FULL LP
# =============================================================================
function run_full_lp_cpi(panel::DataFrame)::Tuple{DataFrame, DataFrame}
    println("\nRunning CPI LP  h = 0 … $H_MAX")
    println("Block bootstrap: $N_BOOT replications, block size = $BLOCK_SIZE months")

    rng      = MersenneTwister(BOOT_SEED)
    all_res  = DataFrame[]
    irf_rows = NamedTuple[]

    for h in 0:H_MAX
        print("\r  h = $h / $H_MAX  ")
        out = block_bootstrap_cpi(panel, h; n_boot = N_BOOT,
                                            block_size = BLOCK_SIZE, rng = rng)
        isnothing(out) && continue
        push!(all_res, out.results)

        # Extract IRF for the :shock coefficient
        β_idx  = findfirst(==(:shock), out.coef_names)
        β_h    = out.results.beta[β_idx]
        bd     = out.boot_draws[:, β_idx]

        push!(irf_rows, (
            horizon  = h,
            beta     = β_h,
            boot_se  = std(bd),
            ci_lo95  = quantile(bd, 0.025),
            ci_hi95  = quantile(bd, 0.975),
            ci_lo90  = quantile(bd, 0.05),
            ci_hi90  = quantile(bd, 0.95),
            ci_lo68  = quantile(bd, 0.16),
            ci_hi68  = quantile(bd, 0.84),
        ))
    end

    println("\nLP estimation complete.")
    results_df = vcat(all_res...)
    irf_df     = DataFrame(irf_rows)
    return results_df, irf_df
end

# =============================================================================
# 10. MAIN
# =============================================================================
function main()
    panel = build_cpi_panel()

    results_df, irf_df = run_full_lp_cpi(panel)

    CSV.write(joinpath(OUTPUT_DIR, "lp_cpi_coefficients.csv"), results_df)
    CSV.write(joinpath(OUTPUT_DIR, "lp_cpi_irf.csv"),          irf_df)
    println("Results saved to $(OUTPUT_DIR)")
    println("\nIRF (first 6 horizons):")
    println(first(irf_df, 6))

    return results_df, irf_df
end

results_df, irf_df = main()