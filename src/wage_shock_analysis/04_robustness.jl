# =============================================================================
# 04_robustness.jl
# Robustness checks:
#   R1. Alternative baseline group (Service = group 6)
#   R2. Alternative lag lengths (L = 6, 24)
#   R3. Alternative cell definition (OccGroup × State)
#   R4. Year fixed effects instead of no time FE
# =============================================================================

using CSV, DataFrames, DataFramesMeta
using LinearAlgebra, Statistics
using Printf

include("02_lp_estimation.jl")

# =============================================================================
# R1. ALTERNATIVE BASELINE GROUP
# =============================================================================
function run_lp_baseline_group(panel::DataFrame, baseline_g::Int;
                                 h_max::Int = H_MAX,
                                 dk_bandwidth::Int = 12)::DataFrame
    println("\nRobustness R1: Baseline group = $baseline_g ($(OCC_LABELS[baseline_g]))")

    # Relabel: drop baseline_g from interaction terms
    all_results = DataFrame[]

    for h in 0:h_max
        df_h = build_lp_data(panel, h)
        isnothing(df_h) && continue

        # Interaction dummies excluding new baseline
        occ_groups = setdiff(1:9, [baseline_g])
        inter_cols = Symbol[]
        for g in occ_groups
            col = Symbol("shock_x_occ", g)
            df_h[!, col] = df_h.shock .* Float64.(df_h.occ_group .== g)
            push!(inter_cols, col)
        end

        lag_cols   = [Symbol("shock_lag", l) for l in 1:L_LAG]
        macro_cols = [:log_oil_lag1, :ffr_lag1, :unrate_lag1]
        cell_cols  = [:age_mean, :female_share, :married_share, :union_share]
        x_cols     = vcat([:shock], inter_cols, lag_cols, macro_cols, cell_cols)
        y_col      = :dep_var

        df_w = within_transform(df_h, [y_col], x_cols, :cell_id)
        df_w = dropmissing(df_w, vcat([y_col], x_cols))
        nrow(df_w) == 0 && continue

        Y = Float64.(df_w[!, y_col])
        X = hcat(ones(nrow(df_w)), Matrix{Float64}(df_w[!, x_cols]))
        β = (X' * X) \ (X' * Y)
        e = Y .- X * β

        all_dates = sort(unique(panel.date))
        date_map  = Dict(d => i for (i, d) in enumerate(all_dates))
        t_idx     = [date_map[d] for d in df_w.date]

        V  = driscoll_kraay_vcov(X, e, t_idx; m = dk_bandwidth)
        se = sqrt.(diag(V))

        coef_names = vcat([:intercept], x_cols)
        res = DataFrame(
            horizon   = h,
            coef_name = coef_names,
            beta      = β,
            se        = se,
            baseline  = baseline_g,
        )
        push!(all_results, res)
    end

    results = vcat(all_results...)
    out = joinpath(OUTPUT_DIR, "robustness_baseline_g$(baseline_g).csv")
    CSV.write(out, results)
    println("  Saved: $out")
    return results
end

# =============================================================================
# R2. ALTERNATIVE LAG LENGTHS
# =============================================================================
function run_lp_lag_length(panel::DataFrame, L::Int;
                             h_max::Int = H_MAX,
                             dk_bandwidth::Int = 12)::DataFrame
    println("\nRobustness R2: Lag length L = $L")

    # Recompute shock lags in panel for new L
    panel_L = copy(panel)
    all_dates = sort(unique(panel_L.date))

    # Need to recompute macro lags if L differs from original L_LAG
    # Assume macro shock lags are stored as shock_lag1..shock_lag_LLAG
    # For L < L_LAG: simply use fewer; for L > L_LAG: need original shock series

    if L > L_LAG
        @warn "L=$L > L_LAG=$L_LAG; additional lags not precomputed. Using L=$L_LAG."
        L = L_LAG
    end

    all_results = DataFrame[]

    for h in 0:h_max
        df_h = build_lp_data(panel_L, h)
        isnothing(df_h) && continue

        occ_groups = 2:9
        inter_cols = Symbol[]
        for g in occ_groups
            col = Symbol("shock_x_occ", g)
            df_h[!, col] = df_h.shock .* Float64.(df_h.occ_group .== g)
            push!(inter_cols, col)
        end

        lag_cols   = [Symbol("shock_lag", l) for l in 1:L]
        macro_cols = [:log_oil_lag1, :ffr_lag1, :unrate_lag1]
        cell_cols  = [:age_mean, :female_share, :married_share, :union_share]
        x_cols     = vcat([:shock], inter_cols, lag_cols, macro_cols, cell_cols)
        y_col      = :dep_var

        df_w = within_transform(df_h, [y_col], x_cols, :cell_id)
        df_w = dropmissing(df_w, vcat([y_col], x_cols))
        nrow(df_w) == 0 && continue

        Y = Float64.(df_w[!, y_col])
        X = hcat(ones(nrow(df_w)), Matrix{Float64}(df_w[!, x_cols]))
        β = (X' * X) \ (X' * Y)
        e = Y .- X * β

        all_dates2 = sort(unique(panel.date))
        date_map   = Dict(d => i for (i, d) in enumerate(all_dates2))
        t_idx      = [date_map[d] for d in df_w.date]

        V  = driscoll_kraay_vcov(X, e, t_idx; m = dk_bandwidth)
        se = sqrt.(diag(V))

        coef_names = vcat([:intercept], x_cols)
        res = DataFrame(
            horizon   = h,
            coef_name = coef_names,
            beta      = β,
            se        = se,
            lag_L     = L,
        )
        push!(all_results, res)
    end

    results = vcat(all_results...)
    out = joinpath(OUTPUT_DIR, "robustness_lagL$(L).csv")
    CSV.write(out, results)
    println("  Saved: $out")
    return results
end

# =============================================================================
# R3. ALTERNATIVE CELL: OccGroup × State
# =============================================================================
function build_cell_panel_state(cps::DataFrame)::DataFrame
    println("\nRobustness R3: Rebuilding cell panel as OccGroup × State...")

    # Need STATEFIP in cps — must be loaded in 01_data_prep
    if !hasproperty(cps, :statefip)
        error("STATEFIP not found in CPS data. Add it in 01_data_prep.jl.")
    end

    cps[!, :cell_id] = string.(cps.occ_group, "_ST", cps.statefip)

    gdf = groupby(cps, [:cell_id, :occ_group, :statefip, :date, :year, :month])
    panel = combine(gdf,
        [:log_rwage, :earnwt] =>
            ((w, wt) -> sum(w .* wt) / sum(wt)) => :log_rwage,
        [:age, :earnwt]     =>
            ((x, wt) -> sum(x .* wt) / sum(wt)) => :age_mean,
        [:female, :earnwt]  =>
            ((x, wt) -> sum(x .* wt) / sum(wt)) => :female_share,
        [:married, :earnwt] =>
            ((x, wt) -> sum(x .* wt) / sum(wt)) => :married_share,
        [:union_d, :earnwt] =>
            ((x, wt) -> sum(x .* wt) / sum(wt)) => :union_share,
        :log_rwage => length => :n_obs,
    )

    panel = @subset(panel, :n_obs .>= 30)
    println("  OccGroup × State cells: ", length(unique(panel.cell_id)))
    return sort(panel, [:cell_id, :date])
end

# =============================================================================
# R4. YEAR FIXED EFFECTS
# =============================================================================
function run_lp_year_fe(panel::DataFrame;
                          h_max::Int = H_MAX,
                          dk_bandwidth::Int = 12)::DataFrame
    println("\nRobustness R4: Adding year fixed effects...")

    all_results = DataFrame[]

    for h in 0:h_max
        df_h = build_lp_data(panel, h)
        isnothing(df_h) && continue

        # Year dummies (drop first year as reference)
        years     = sort(unique(df_h.year))
        ref_year  = first(years)
        year_cols = Symbol[]
        for yr in years[2:end]
            col = Symbol("yr_", yr)
            df_h[!, col] = Float64.(df_h.year .== yr)
            push!(year_cols, col)
        end

        occ_groups = 2:9
        inter_cols = Symbol[]
        for g in occ_groups
            col = Symbol("shock_x_occ", g)
            df_h[!, col] = df_h.shock .* Float64.(df_h.occ_group .== g)
            push!(inter_cols, col)
        end

        lag_cols   = [Symbol("shock_lag", l) for l in 1:L_LAG]
        macro_cols = [:log_oil_lag1, :ffr_lag1, :unrate_lag1]
        cell_cols  = [:age_mean, :female_share, :married_share, :union_share]

        # NOTE: year FEs are NOT demeaned (they vary cross-sectionally within year)
        # Shock_t still identified via cross-sectional interaction variation
        x_cols = vcat([:shock], inter_cols, lag_cols, macro_cols, cell_cols, year_cols)
        y_col  = :dep_var

        df_w = within_transform(df_h, [y_col],
                                 vcat([:shock], inter_cols, lag_cols, macro_cols, cell_cols),
                                 :cell_id)
        # Year dummies added after demeaning (not demeaned themselves)
        for col in year_cols
            df_w[!, col] = df_h[!, col]
        end

        all_x = vcat([:shock], inter_cols, lag_cols, macro_cols, cell_cols, year_cols)
        df_w = dropmissing(df_w, vcat([y_col], all_x))
        nrow(df_w) == 0 && continue

        Y = Float64.(df_w[!, y_col])
        X = hcat(ones(nrow(df_w)), Matrix{Float64}(df_w[!, all_x]))
        β = (X' * X) \ (X' * Y)
        e = Y .- X * β

        all_dates = sort(unique(panel.date))
        date_map  = Dict(d => i for (i, d) in enumerate(all_dates))
        t_idx     = [date_map[d] for d in df_w.date]

        V  = driscoll_kraay_vcov(X, e, t_idx; m = dk_bandwidth)
        se = sqrt.(diag(V))

        coef_names = vcat([:intercept], all_x)
        res = DataFrame(
            horizon   = h,
            coef_name = coef_names,
            beta      = β,
            se        = se,
        )
        push!(all_results, res)
    end

    results = vcat(all_results...)
    out = joinpath(OUTPUT_DIR, "robustness_year_fe.csv")
    CSV.write(out, results)
    println("  Saved: $out")
    return results
end

# =============================================================================
# MAIN ROBUSTNESS RUNNER
# =============================================================================
function run_all_robustness(panel::DataFrame)
    println("\n=== Running All Robustness Checks ===")

    r1 = run_lp_baseline_group(panel, 6)          # Service as baseline
    r2a = run_lp_lag_length(panel, 6)             # L = 6
    r2b = run_lp_lag_length(panel, L_LAG)         # L = 12 (same as main, sanity check)
    r4 = run_lp_year_fe(panel)                    # Year FE

    # R3 requires raw CPS data with STATEFIP — called separately if available
    println("\nNote: R3 (OccGroup × State) requires re-running 01_data_prep.jl")
    println("      with STATEFIP included in CPS extraction.")

    println("\nAll robustness checks complete.")
    return (r1 = r1, r2a = r2a, r2b = r2b, r4 = r4)
end
