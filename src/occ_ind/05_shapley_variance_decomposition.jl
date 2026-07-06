# =============================================================================
# 05_variance_decomposition_occ_ind.jl
#
# Variance decomposition for occupation-within-industry IRFs
# with Shapley R2 decomposition
#
# Input structure:
#   ../../result/occ_ind_4/{outcome}/{industry}/irf_occ*.csv
#
# Each CSV contains:
#   horizon, beta, ci_lo90, ci_hi90, ci_lo68, ci_hi68
#
# Outputs:
#   ../../result/occ_ind_4/variance_decomposition/
#       occ_ind_cirf_panel.csv
#       variance_decomposition_r2.csv
#       variance_decomposition_long.csv
#       variance_decomposition_shapley_long.csv
# =============================================================================

using CSV
using DataFrames
using GLM
using StatsModels
using Statistics
using Printf
using CategoricalArrays

# =============================================================================
# SETTINGS
# =============================================================================

const BASE_DIR = joinpath(@__DIR__, "..", "..", "result", "occ_ind_4")
const OUTPUT_DIR = joinpath(BASE_DIR, "variance_decomposition")

mkpath(OUTPUT_DIR)

const OUTCOMES = [
    "employment",
    "hourly_rate",
    "hours",
    "income",
    "income_share_var",
    "inequality",
    "median",
    "unemployment",
]

const INDUSTRIES = [
    "Energy_intensive",
    "Manufacturing_Construction",
    "Services",
    "Trade",
]

const HORIZONS = [12, 24, 36]

# =============================================================================
# HELPERS
# =============================================================================

function parse_occ_from_filename(fp::String)
    fname = basename(fp)

    m = match(r"irf_occ(\d+)_(.+)\.csv$", fname)

    if m === nothing
        error("Cannot parse occupation from filename: $fname")
    end

    occ_id = parse(Int, m.captures[1])
    occ_name = String(m.captures[2])

    return occ_id, occ_name
end


function compute_cirf(irf::DataFrame, H::Int)
    if !("horizon" in names(irf)) || !("beta" in names(irf))
        error("IRF file must contain columns: horizon and beta")
    end

    sub = irf[irf.horizon .<= H, :]

    if nrow(sub) == 0
        return missing
    end

    return sum(skipmissing(sub.beta))
end


function compute_peak(irf::DataFrame, H::Int)
    if !("horizon" in names(irf)) || !("beta" in names(irf))
        error("IRF file must contain columns: horizon and beta")
    end

    sub = irf[irf.horizon .<= H, :]

    if nrow(sub) == 0
        return missing, missing
    end

    idx = argmax(sub.beta)

    return sub.beta[idx], sub.horizon[idx]
end


function safe_r2(df::DataFrame, formula)
    try
        fit = lm(formula, df)

        r2_val = r2(fit)
        adjr2_val = adjr2(fit)

        nobs_val = nobs(fit)
        dof_resid_val = dof_residual(fit)

        # df_model = number of estimated coefficients excluding intercept
        df_model_val = length(coef(fit)) - 1

        return r2_val, adjr2_val, nobs_val, df_model_val, dof_resid_val

    catch e
        @warn "Regression failed" formula=formula exception=e
        return missing, missing, missing, missing, missing
    end
end


function finite_rows(df::DataFrame, y::Symbol)
    return .!ismissing.(df[!, y]) .& isfinite.(df[!, y])
end


function shapley_r2_two_groups(r2_group1, r2_group2, r2_both)
    if ismissing(r2_group1) || ismissing(r2_group2) || ismissing(r2_both)
        return missing, missing, missing, missing
    end

    shapley_group1 = 0.5 * (r2_group1 + (r2_both - r2_group2))
    shapley_group2 = 0.5 * (r2_group2 + (r2_both - r2_group1))

    if r2_both == 0
        share_group1 = missing
        share_group2 = missing
    else
        share_group1 = shapley_group1 / r2_both
        share_group2 = shapley_group2 / r2_both
    end

    return shapley_group1, shapley_group2, share_group1, share_group2
end


# =============================================================================
# 1. BUILD OCCUPATION × INDUSTRY CIRF PANEL
# =============================================================================

rows = DataFrame()

for outcome in OUTCOMES
    outcome_dir = joinpath(BASE_DIR, outcome)

    if !isdir(outcome_dir)
        @warn "Outcome folder not found" outcome_dir
        continue
    end

    for industry in INDUSTRIES
        industry_dir = joinpath(outcome_dir, industry)

        if !isdir(industry_dir)
            @warn "Industry folder not found" industry_dir
            continue
        end

        files = filter(
            fp -> occursin(r"irf_occ\d+_.+\.csv$", basename(fp)),
            readdir(industry_dir; join=true)
        )

        sort!(files)

        if isempty(files)
            @warn "No IRF CSV files found" industry_dir
            continue
        end

        for fp in files
            occ_id, occ_name = parse_occ_from_filename(fp)
            irf = CSV.read(fp, DataFrame)

            row = Dict{Symbol, Any}()

            row[:outcome] = outcome
            row[:industry] = industry
            row[:occ_id] = occ_id
            row[:occupation] = occ_name
            row[:file] = fp

            for H in HORIZONS
                row[Symbol("CIRF_$H")] = compute_cirf(irf, H)

                peak_beta, peak_horizon = compute_peak(irf, H)
                row[Symbol("Peak_$H")] = peak_beta
                row[Symbol("PeakMonth_$H")] = peak_horizon
            end

            push!(rows, row; cols=:union)
        end
    end
end

if nrow(rows) == 0
    error("No valid IRF files were found. Please check BASE_DIR, OUTCOMES, INDUSTRIES, and filename format.")
end

sort!(rows, [:outcome, :industry, :occ_id])

cirf_panel_path = joinpath(OUTPUT_DIR, "occ_ind_cirf_panel.csv")
CSV.write(cirf_panel_path, rows)

println("Saved CIRF panel to: $cirf_panel_path")
println(first(rows, min(6, nrow(rows))))


# =============================================================================
# 2. VARIANCE DECOMPOSITION WITH SHAPLEY R2
# =============================================================================

vd_rows = DataFrame()

for outcome in OUTCOMES
    df_outcome = rows[rows.outcome .== outcome, :]

    if nrow(df_outcome) == 0
        continue
    end

    for H in HORIZONS
        y = Symbol("CIRF_$H")

        df_reg = dropmissing(df_outcome, [y, :industry, :occ_id])
        df_reg = df_reg[finite_rows(df_reg, y), :]

        if nrow(df_reg) == 0
            continue
        end

        # Treat industry and occupation as categorical variables
        df_reg.industry = categorical(df_reg.industry)
        df_reg.occ_id = categorical(df_reg.occ_id)

        # Model 1: Industry only
        f_ind = Term(y) ~ term(:industry)

        # Model 2: Occupation only
        f_occ = Term(y) ~ term(:occ_id)

        # Model 3: Industry + Occupation
        f_both = Term(y) ~ term(:industry) + term(:occ_id)

        r2_ind, adjr2_ind, n_ind, dfm_ind, dfr_ind = safe_r2(df_reg, f_ind)
        r2_occ, adjr2_occ, n_occ, dfm_occ, dfr_occ = safe_r2(df_reg, f_occ)
        r2_both, adjr2_both, n_both, dfm_both, dfr_both = safe_r2(df_reg, f_both)

        # Incremental R2
        inc_occ_given_ind = ismissing(r2_both) || ismissing(r2_ind) ? missing : r2_both - r2_ind
        inc_ind_given_occ = ismissing(r2_both) || ismissing(r2_occ) ? missing : r2_both - r2_occ

        # Incremental adjusted R2
        inc_adj_occ_given_ind = ismissing(adjr2_both) || ismissing(adjr2_ind) ? missing : adjr2_both - adjr2_ind
        inc_adj_ind_given_occ = ismissing(adjr2_both) || ismissing(adjr2_occ) ? missing : adjr2_both - adjr2_occ

        # Shapley R2 decomposition
        #
        # Industry Shapley R2:
        #   0.5 * [R2(industry) + R2(industry + occupation) - R2(occupation)]
        #
        # Occupation Shapley R2:
        #   0.5 * [R2(occupation) + R2(industry + occupation) - R2(industry)]
        #
        # These two components sum to R2(industry + occupation).
        shapley_r2_industry,
        shapley_r2_occupation,
        shapley_share_industry,
        shapley_share_occupation = shapley_r2_two_groups(r2_ind, r2_occ, r2_both)

        # Numerical check:
        # This should be close to zero when R2 values are not missing.
        shapley_check = (
            ismissing(shapley_r2_industry) ||
            ismissing(shapley_r2_occupation) ||
            ismissing(r2_both)
        ) ? missing : shapley_r2_industry + shapley_r2_occupation - r2_both

        push!(
            vd_rows,
            (
                outcome = outcome,
                H = H,
                n_obs = nrow(df_reg),

                R2_industry_only = r2_ind,
                R2_occupation_only = r2_occ,
                R2_industry_occupation = r2_both,

                AdjR2_industry_only = adjr2_ind,
                AdjR2_occupation_only = adjr2_occ,
                AdjR2_industry_occupation = adjr2_both,

                Incremental_R2_occupation_given_industry = inc_occ_given_ind,
                Incremental_R2_industry_given_occupation = inc_ind_given_occ,

                Incremental_AdjR2_occupation_given_industry = inc_adj_occ_given_ind,
                Incremental_AdjR2_industry_given_occupation = inc_adj_ind_given_occ,

                Shapley_R2_industry = shapley_r2_industry,
                Shapley_R2_occupation = shapley_r2_occupation,

                Shapley_share_industry = shapley_share_industry,
                Shapley_share_occupation = shapley_share_occupation,

                Shapley_check = shapley_check,

                df_model_industry_only = dfm_ind,
                df_model_occupation_only = dfm_occ,
                df_model_industry_occupation = dfm_both,

                df_resid_industry_only = dfr_ind,
                df_resid_occupation_only = dfr_occ,
                df_resid_industry_occupation = dfr_both,
            );
            cols = :union
        )
    end
end

vd_path = joinpath(OUTPUT_DIR, "variance_decomposition_r2.csv")
CSV.write(vd_path, vd_rows)

println("Saved variance decomposition table to: $vd_path")
println(vd_rows)


# =============================================================================
# 3. LONG FORMAT OUTPUT FOR MODEL R2 PLOTTING
# =============================================================================

vd_long = DataFrame(
    outcome = String[],
    H = Int[],
    model = String[],
    R2 = Union{Missing, Float64}[],
    AdjR2 = Union{Missing, Float64}[],
)

for row in eachrow(vd_rows)
    push!(
        vd_long,
        (
            outcome = row.outcome,
            H = row.H,
            model = "Industry only",
            R2 = row.R2_industry_only,
            AdjR2 = row.AdjR2_industry_only,
        )
    )

    push!(
        vd_long,
        (
            outcome = row.outcome,
            H = row.H,
            model = "Occupation only",
            R2 = row.R2_occupation_only,
            AdjR2 = row.AdjR2_occupation_only,
        )
    )

    push!(
        vd_long,
        (
            outcome = row.outcome,
            H = row.H,
            model = "Industry + Occupation",
            R2 = row.R2_industry_occupation,
            AdjR2 = row.AdjR2_industry_occupation,
        )
    )
end

vd_long_path = joinpath(OUTPUT_DIR, "variance_decomposition_long.csv")
CSV.write(vd_long_path, vd_long)

println("Saved long-format variance decomposition table to: $vd_long_path")


# =============================================================================
# 4. LONG FORMAT OUTPUT FOR SHAPLEY R2 PLOTTING
# =============================================================================

shapley_long = DataFrame(
    outcome = String[],
    H = Int[],
    component = String[],
    Shapley_R2 = Union{Missing, Float64}[],
    Shapley_share = Union{Missing, Float64}[],
)

for row in eachrow(vd_rows)
    push!(
        shapley_long,
        (
            outcome = row.outcome,
            H = row.H,
            component = "Industry",
            Shapley_R2 = row.Shapley_R2_industry,
            Shapley_share = row.Shapley_share_industry,
        )
    )

    push!(
        shapley_long,
        (
            outcome = row.outcome,
            H = row.H,
            component = "Occupation",
            Shapley_R2 = row.Shapley_R2_occupation,
            Shapley_share = row.Shapley_share_occupation,
        )
    )
end

shapley_long_path = joinpath(OUTPUT_DIR, "variance_decomposition_shapley_long.csv")
CSV.write(shapley_long_path, shapley_long)

println("Saved long-format Shapley decomposition table to: $shapley_long_path")


# =============================================================================
# 5. PRINT UNEMPLOYMENT RESULT ONLY
# =============================================================================

unemp = vd_rows[vd_rows.outcome .== "unemployment", :]

if nrow(unemp) > 0
    println()
    println("Variance decomposition for unemployment:")

    show(
        unemp[
            :,
            [
                :outcome,
                :H,
                :n_obs,

                :R2_industry_only,
                :R2_occupation_only,
                :R2_industry_occupation,

                :Incremental_R2_occupation_given_industry,
                :Incremental_R2_industry_given_occupation,

                :Shapley_R2_industry,
                :Shapley_R2_occupation,
                :Shapley_share_industry,
                :Shapley_share_occupation,

                :Shapley_check,
            ],
        ],
        allrows = true,
        allcols = true,
    )

    println()
end


# =============================================================================
# 6. PRINT SHAPLEY CHECK SUMMARY
# =============================================================================

println()
println("Shapley decomposition check:")

if "Shapley_check" in names(vd_rows)
    valid_check = skipmissing(vd_rows.Shapley_check)

    if isempty(collect(valid_check))
        println("No valid Shapley checks available.")
    else
        max_abs_check = maximum(abs.(skipmissing(vd_rows.Shapley_check)))
        println("Maximum absolute Shapley check error: $max_abs_check")
        println("This should be numerically close to 0.")
    end
end

println()
println("All outputs saved to: $OUTPUT_DIR")