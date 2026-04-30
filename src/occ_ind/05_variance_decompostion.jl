# =============================================================================
# 05_variance_decomposition_occ_ind.jl
#
# Variance decomposition for occupation-within-industry IRFs
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

sort!(rows, [:outcome, :industry, :occ_id])

cirf_panel_path = joinpath(OUTPUT_DIR, "occ_ind_cirf_panel.csv")
CSV.write(cirf_panel_path, rows)

println("Saved CIRF panel to: $cirf_panel_path")
println(first(rows, min(6, nrow(rows))))


# =============================================================================
# 2. VARIANCE DECOMPOSITION
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

        inc_occ_given_ind = ismissing(r2_both) || ismissing(r2_ind) ? missing : r2_both - r2_ind
        inc_ind_given_occ = ismissing(r2_both) || ismissing(r2_occ) ? missing : r2_both - r2_occ

        inc_adj_occ_given_ind = ismissing(adjr2_both) || ismissing(adjr2_ind) ? missing : adjr2_both - adjr2_ind
        inc_adj_ind_given_occ = ismissing(adjr2_both) || ismissing(adjr2_occ) ? missing : adjr2_both - adjr2_occ

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
# 3. LONG FORMAT OUTPUT FOR PLOTTING
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
# 4. PRINT UNEMPLOYMENT RESULT ONLY
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
            ],
        ],
        allrows = true,
        allcols = true,
    )

    println()
end
