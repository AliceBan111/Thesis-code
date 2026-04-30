using DataFrames
using CSV
using XLSX
using StatFiles
using Statistics
using GLM
using StatsModels
using CategoricalArrays

# =============================================================================
# Within-industry task regression only
#
# LHS:
#   cumulative_beta_{o,i} = sum_h beta_{o,i,h}
#
# Regression:
#   cumulative_beta_{o,i}
#       = industry FE_i
#       + θ abstract_o
#       + γ routine_o
#       + δ manual_o
#       + ε_{o,i}
#
# File structure:
#   ../../result/occ_ind_4/{outcome}/{industry}/irf_occ*.csv
# =============================================================================


# =============================================================================
# 1. SETTINGS
# =============================================================================

const BASE_IRF_DIR = "../../result/occ_ind_4"

const OUTCOME_LIST = [
    "employment",
    "hourly_rate",
    "income_share_var",
    "inequality",
    "median",
    "unemployment",
    "income",
    "hours"
]

const INDUSTRY_LIST = [
    "Energy_intensive",
    "Manufacturing_Construction",
    "Services",
    "Trade"
]

const OUTPUT_DIR = "../../result/occ_ind_4/task_regression"

const H_MIN = 1
const H_MAX = 36

const DTA_PATH  = "../../data/occ1990dd_task_alm.dta"
const XLSX_PATH = "../../result/mapping/mapping_done.xlsx"


# =============================================================================
# 2. OCCUPATION MAPS
# =============================================================================

const OCC_FILE_MAP = Dict(
    1 => "irf_occ1_Managerial.csv",
    2 => "irf_occ2_Professional_specialty.csv",
    3 => "irf_occ3_High_tech.csv",
    4 => "irf_occ4_Sales.csv",
    5 => "irf_occ5_Administrative_support.csv",
    6 => "irf_occ6_Service.csv",
    7 => "irf_occ7_Farming_forestry_construction.csv",
    8 => "irf_occ8_Precision_production_repair.csv",
    9 => "irf_occ9_Machine_operators_transport.csv"
)

const OCC_NAME_MAP = Dict(
    1 => "group1_Managerial",
    2 => "group2_Professional_specialty",
    3 => "group3_High_tech",
    4 => "group4_Sales",
    5 => "group5_Administrative_support",
    6 => "group6_Service",
    7 => "group7_Farming_forestry_construction",
    8 => "group8_Precision_production_repair",
    9 => "group9_Machine_operators_transport"
)

const OCC_SHORT_NAME_MAP = Dict(
    1 => "Managerial",
    2 => "Prof.",
    3 => "High Tech",
    4 => "Sales",
    5 => "Admin",
    6 => "Service",
    7 => "Constr./Farm",
    8 => "Prod./Repair",
    9 => "Mach./Transp."
)


# =============================================================================
# 3. PREPARE DORN / AUTOR TASK DATA
# =============================================================================

function prepare_task_data(
    dta_path::String = DTA_PATH,
    xlsx_path::String = XLSX_PATH
)
    isfile(dta_path) || error("Task DTA file not found: $dta_path")
    isfile(xlsx_path) || error("Mapping XLSX file not found: $xlsx_path")

    df_tasks = DataFrame(load(dta_path))

    select!(
        df_tasks,
        :occ1990dd,
        :task_abstract,
        :task_routine,
        :task_manual
    )

    df_mapping = DataFrame(XLSX.readtable(xlsx_path, "Sheet1"))

    # B column -> occ1990dd
    # F column -> group
    # H column -> weights
    select!(df_mapping, 2 => :occ1990dd, 6 => :group, 8 => :weights)

    df_mapping.group = convert.(Int, df_mapping.group)
    df_mapping.weights = convert.(Float64, df_mapping.weights)

    df_merged = innerjoin(df_tasks, df_mapping, on = :occ1990dd)
    dropmissing!(df_merged)

    calc_wmean(x, w) = sum(w) == 0 ? missing : sum(x .* w) / sum(w)

    df_y = combine(groupby(df_merged, :group),
        [:task_abstract, :weights] => calc_wmean => :task_abstract,
        [:task_routine,  :weights] => calc_wmean => :task_routine,
        [:task_manual,   :weights] => calc_wmean => :task_manual
    )

    df_y.occ_id = df_y.group
    df_y.Occupational_Group = [OCC_NAME_MAP[g] for g in df_y.group]
    df_y.Plot_Label = [OCC_SHORT_NAME_MAP[g] for g in df_y.group]

    select!(
        df_y,
        :occ_id,
        :Occupational_Group,
        :Plot_Label,
        :task_abstract,
        :task_routine,
        :task_manual
    )

    sort!(df_y, :occ_id)

    return df_y
end


# =============================================================================
# 4. READ ONE IRF FILE
# =============================================================================

function read_single_irf(
    outcome::String,
    industry::String,
    occ_id::Int;
    base_dir::String = BASE_IRF_DIR,
    h_min::Int = H_MIN,
    h_max::Int = H_MAX
)
    filename = OCC_FILE_MAP[occ_id]
    filepath = joinpath(base_dir, outcome, industry, filename)

    if !isfile(filepath)
        @warn "Missing file: $filepath"
        return nothing
    end

    df = CSV.read(filepath, DataFrame; missingstring = ["", "NA", "N/A"])

    required_cols = [:horizon, :beta]
    for c in required_cols
        c in propertynames(df) || error("Column $c not found in $filepath")
    end

    # Keep only finite beta values within horizon range.
    filter!(
        row -> !ismissing(row.horizon) &&
               !ismissing(row.beta) &&
               isfinite(Float64(row.beta)) &&
               row.horizon >= h_min &&
               row.horizon <= h_max,
        df
    )

    if nrow(df) == 0
        @warn "No valid finite beta observations in: $filepath"
        return nothing
    end

    cumulative_beta = sum(Float64.(df.beta))

    return DataFrame(
        outcome = outcome,
        industry = industry,
        occ_id = occ_id,
        Occupational_Group = OCC_NAME_MAP[occ_id],
        cumulative_beta = cumulative_beta,
        n_horizons_used = nrow(df)
    )
end


# =============================================================================
# 5. BUILD INDUSTRY × OCCUPATION DATASET FOR ONE OUTCOME
# =============================================================================

function build_within_industry_irf_dataset(
    outcome::String;
    base_dir::String = BASE_IRF_DIR,
    industries::Vector{String} = INDUSTRY_LIST
)
    rows = DataFrame[]

    for industry in industries
        for occ_id in sort(collect(keys(OCC_FILE_MAP)))
            tmp = read_single_irf(outcome, industry, occ_id; base_dir = base_dir)

            if tmp !== nothing
                push!(rows, tmp)
            end
        end
    end

    isempty(rows) && error("No IRF data loaded for outcome = $outcome")

    df = vcat(rows...)
    sort!(df, [:industry, :occ_id])

    return df
end


# =============================================================================
# 6. RUN THREE SEPARATE TASK REGRESSIONS WITH INDUSTRY FIXED EFFECTS
#
# For each task k ∈ {abstract, routine, manual}:
#
#   cumulative_beta_{o,i}
#       = α + β_k task^k_o + industry FE_i + ε_{o,i}
#
# =============================================================================

const TASK_SPEC_LIST = [
    ("abstract", :task_abstract),
    ("routine",  :task_routine),
    ("manual",   :task_manual)
]

function run_single_task_regression(df::DataFrame, task_var::Symbol)
    df_reg = copy(df)
    df_reg.industry = categorical(df_reg.industry)

    fml = Term(:cumulative_beta) ~ Term(task_var) + Term(:industry)

    model = lm(fml, df_reg)

    return model
end

function run_three_separate_task_regressions(df::DataFrame)
    models = Dict{String, Any}()

    for (task_name, task_var) in TASK_SPEC_LIST
        models[task_name] = run_single_task_regression(df, task_var)
    end

    return models
end

# =============================================================================
# 7. SAVE REGRESSION TABLES
# =============================================================================

function save_single_task_regression_table(
    model,
    outcome::String,
    task_name::String,
    task_var::Symbol,
    output_dir::String
)
    ct = coeftable(model)

    out = DataFrame(
        outcome = outcome,
        model = "single_task_industry_FE",
        task = task_name,
        task_var = String(task_var),
        term = coefnames(model),
        estimate = ct.cols[1],
        std_error = ct.cols[2],
        t_stat = ct.cols[3],
        p_value = ct.cols[4],
        r2 = r2(model),
        adj_r2 = adjr2(model),
        nobs = nobs(model)
    )

    path = joinpath(output_dir, "$(outcome)_$(task_name)_industry_FE_results.csv")
    CSV.write(path, out)

    return out
end

function save_three_separate_task_regression_tables(
    models::Dict{String, Any},
    outcome::String,
    output_dir::String
)
    rows = DataFrame[]

    for (task_name, task_var) in TASK_SPEC_LIST
        model = models[task_name]

        reg_table = save_single_task_regression_table(
            model,
            outcome,
            task_name,
            task_var,
            output_dir
        )

        push!(rows, reg_table)
    end

    out = vcat(rows...)

    path = joinpath(output_dir, "$(outcome)_three_separate_task_industry_FE_results.csv")
    CSV.write(path, out)

    return out
end


# =============================================================================
# 8. PROCESS ONE OUTCOME
# =============================================================================

function process_one_outcome(
    outcome::String,
    df_tasks::DataFrame;
    output_dir::String = OUTPUT_DIR
)
    println("\n============================================================")
    println("Outcome: $outcome")
    println("============================================================")

    mkpath(output_dir)

    df_irf = build_within_industry_irf_dataset(outcome)

    df = innerjoin(
        df_irf,
        df_tasks,
        on = [:occ_id, :Occupational_Group]
    )

    dropmissing!(
        df,
        [
            :cumulative_beta,
            :task_abstract,
            :task_routine,
            :task_manual,
            :industry
        ]
    )

    filter!(
        row -> isfinite(row.cumulative_beta) &&
               isfinite(row.task_abstract) &&
               isfinite(row.task_routine) &&
               isfinite(row.task_manual),
        df
    )

    if nrow(df) < 10
        @warn "Too few valid observations after cleaning for outcome = $outcome. N = $(nrow(df))"
        return nothing
    end

    df.industry = categorical(df.industry)

    input_path = joinpath(output_dir, "$(outcome)_within_industry_task_data.csv")
    CSV.write(input_path, df)
    println("Saved regression input: $input_path")
    println("N = ", nrow(df))

    models = run_three_separate_task_regressions(df)

    reg_table = save_three_separate_task_regression_tables(
        models,
        outcome,
        output_dir
    )

    println("\n--- Separate task regression results with industry FE ---")

    for (task_name, task_var) in TASK_SPEC_LIST
        model = models[task_name]

        println("\nTask: $task_name")
        println("Formula: cumulative_beta ~ $(String(task_var)) + industry")
        println(coeftable(model))
        println("R² = ", r2(model))
        println("Adj. R² = ", adjr2(model))
    end

    return Dict(
        "data" => df,
        "models" => models,
        "reg_table" => reg_table
    )
end


# =============================================================================
# 9. BUILD SUMMARY TABLE ACROSS OUTCOMES
# =============================================================================

function build_summary_table(output_dir::String = OUTPUT_DIR)
    rows = DataFrame[]

    for outcome in OUTCOME_LIST
        path = joinpath(output_dir, "$(outcome)_three_separate_task_industry_FE_results.csv")

        if !isfile(path)
            @warn "Missing regression result file: $path"
            continue
        end

        df = CSV.read(path, DataFrame)

        for (task_name, task_var) in TASK_SPEC_LIST
            sub = filter(
                row -> row.task == task_name && row.term == String(task_var),
                df
            )

            if nrow(sub) == 1
                push!(
                    rows,
                    DataFrame(
                        outcome = outcome,
                        model = "single_task_industry_FE",
                        task = task_name,
                        task_var = String(task_var),
                        estimate = sub.estimate[1],
                        std_error = sub.std_error[1],
                        t_stat = sub.t_stat[1],
                        p_value = sub.p_value[1],
                        r2 = sub.r2[1],
                        adj_r2 = sub.adj_r2[1],
                        nobs = sub.nobs[1]
                    )
                )
            else
                @warn "Task coefficient not found or duplicated" outcome task_name task_var n = nrow(sub)
            end
        end
    end

    isempty(rows) && error("No summary rows found.")

    summary = vcat(rows...)

    output_path = joinpath(output_dir, "summary_three_separate_tasks_with_industry_FE.csv")
    CSV.write(output_path, summary)

    println("\nSaved summary table: $output_path")

    return summary
end


# =============================================================================
# 10. MAIN
# =============================================================================

function main()
    cd(@__DIR__)

    mkpath(OUTPUT_DIR)

    println("Preparing Dorn/Autor task data...")
    df_tasks = prepare_task_data()

    task_path = joinpath(OUTPUT_DIR, "occupation_task_measures.csv")
    CSV.write(task_path, df_tasks)
    println("Saved task measures: $task_path")

    results = Dict{String, Any}()

    for outcome in OUTCOME_LIST
        try
            res = process_one_outcome(outcome, df_tasks)
            if res !== nothing
                results[outcome] = res
            end
        catch err
            @warn "Failed for outcome = $outcome" exception = err
        end
    end

    summary = build_summary_table()

    println("\n============================================================")
    println("Summary: separate regressions cumulative_beta ~ task + industry FE")
    println("============================================================")
    show(summary, allrows = true, allcols = true)

    println("\n\nAll done.")
end


main()