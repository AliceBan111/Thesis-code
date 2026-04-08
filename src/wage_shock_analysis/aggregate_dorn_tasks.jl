# =============================================================================
# 10_aggregate_dorn_tasks.jl
# Aggregate Autor-Dorn task measures to 9 broad occupation groups
# Scheme 1: Employment-weighted aggregation + Coverage check
# Assumes Dorn CSV and cps_00017.dat are in the same data directory
# =============================================================================

using CSV, DataFrames, DataFramesMeta, Statistics, Dates, StatFiles

# =============================================================================
# 0. PATHS & CONFIGURATION
# =============================================================================
const DATA_DIR   = joinpath(@__DIR__, "..", "..", "data")  # Adjust if data folder differs
const OUTPUT_DIR = joinpath(@__DIR__, "..", "..", "result", "occ", "decomposition")
mkpath(OUTPUT_DIR)

const CPS_FILE   = joinpath(DATA_DIR, "cps_00017.dat")
const DORN_FILE  = joinpath(DATA_DIR, "occ1990dd_task_alm.dta") 

const BASE_START = Date(1994, 1, 1)
const BASE_END   = Date(1996, 12, 1)

# =============================================================================
# 1. OCCUPATION CLASSIFICATION (Exact match to your specification)
# =============================================================================
function classify_occ1990(occ::Union{Integer,Missing})::Union{Int,Missing}
    ismissing(occ) && return missing
    occ in 3:37                                           && return 1
    occ in 43:200                                         && return 2
    occ in 203:235                                        && return 3
    occ in 243:283                                        && return 4
    occ in 303:389                                        && return 5
    occ in 405:469                                        && return 6
    (occ in 473:498 || occ in 558:599 || occ in 614:617) && return 7
    (occ in 503:549 || occ in 628:699)                   && return 8
    (occ in 703:799 || occ in 803:889)                   && return 9
    return missing
end

const OCC_LABELS = Dict(
    1 => "Managerial", 2 => "Professional_specialty", 3 => "High_tech",
    4 => "Sales", 5 => "Administrative_support", 6 => "Service",
    7 => "Farming_forestry_construction", 8 => "Precision_production_repair",
    9 => "Machine_operators_transport"
)

# =============================================================================
# 2. PARSE CPS FIXED-WIDTH FOR BASELINE EMPLOYMENT WEIGHTS
# =============================================================================
function parse_cps_baseline_weights(cps_path::String)::DataFrame
    println("Parsing CPS baseline weights (1994-1996)...")
    
    # Pre-allocate vectors (rough estimate: ~5M lines / ~45 chars)
    est_rows = 5_000_000
    year_v    = Vector{Int32}(undef, est_rows)
    month_v   = Vector{Int32}(undef, est_rows)
    empstat_v = Vector{Int32}(undef, est_rows)
    age_v     = Vector{Int32}(undef, est_rows)
    occ_v     = Vector{Int32}(undef, est_rows)
    classwkr_v= Vector{Int32}(undef, est_rows)
    uhrs_v    = Vector{Float32}(undef, est_rows)
    earnwt_v  = Vector{Float64}(undef, est_rows)
    
    n = 0
    open(cps_path, "r") do f
        for line in eachline(f)
            length(line) < 47 && continue
            
            yr   = parse(Int32, @view line[1:4])
            mo   = parse(Int32, @view line[5:6])
            d    = Date(yr, mo, 1)
            (d < BASE_START || d > BASE_END) && continue
            
            # Sample selection filters (match your earnings baseline)
            empstat = parse(Int32, @view line[25:26])
            (empstat ∉ (10, 12)) && continue  # At work / has job not at work
            
            age = parse(Int32, @view line[21:22])
            (age < 16 || age > 64) && continue
            
            occ = tryparse(Int32, strip(@view line[27:29]))
            isnothing(occ) && continue
            
            classwkr = parse(Int32, @view line[33:34])
            (classwkr ∉ (21,22,23,24,25,27,28)) && continue  # Wage/salary only
            
            uhrs = tryparse(Float32, strip(@view line[35:37]))
            isnothing(uhrs) || uhrs <= 0.0 || uhrs > 105.0 && continue
            
            earnwt_raw = tryparse(Float64, strip(@view line[38:47]))
            isnothing(earnwt_raw) && continue
            
            n += 1
            if n > length(year_v)
                new_cap = round(Int, length(year_v) * 1.5)
                foreach(v -> resize!(v, new_cap), (year_v, month_v, empstat_v, age_v, occ_v, classwkr_v, uhrs_v, earnwt_v))
            end
            
            year_v[n]    = yr
            month_v[n]   = mo
            empstat_v[n] = empstat
            age_v[n]     = age
            occ_v[n]     = occ
            classwkr_v[n]= classwkr
            uhrs_v[n]    = uhrs
            earnwt_v[n]  = earnwt_raw / 10_000.0  # 4 implied decimals
        end
    end
    
    cps_df = DataFrame(
        year=year_v[1:n], month=month_v[1:n], empstat=empstat_v[1:n],
        age=age_v[1:n], occ1990=occ_v[1:n], classwkr=classwkr_v[1:n],
        uhrsworkt=uhrs_v[1:n], earnwt=earnwt_v[1:n]
    )
    
    println("  Baseline sample size: $(nrow(cps_df))")
    return cps_df
end

# =============================================================================
# 3. LOAD DORN TASK DATA & CHECK COVERAGE
# =============================================================================
function load_and_aggregate_tasks(cps_df::DataFrame, dorn_path::String)
    println("Loading Dorn task measures...")
    dorn = DataFrame(load(dorn_path))

    rename!(dorn, Symbol.(replace.(names(dorn), r"[^a-zA-Z0-9_]" => "_")))

    hasproperty(dorn, :occ1990dd) && rename!(dorn, :occ1990dd => :occ1990_code)
    hasproperty(dorn, Symbol("task_a~t")) && rename!(dorn, Symbol("task_a~t") => :task_abstract)
    hasproperty(dorn, Symbol("task_r~e")) && rename!(dorn, Symbol("task_r~e") => :task_routine)
    hasproperty(dorn, Symbol("task_m~l")) && rename!(dorn, Symbol("task_m~l") => :task_manual)
    
    
    println("  Dorn records: $(nrow(dorn)), Unique occ codes: $(length(unique(dorn.occ1990_code)))")
    
    # Compute employment weight per fine occupation
    cps_weights = combine(groupby(cps_df, :occ1990), :earnwt => sum => :emp_weight)
    println("  CPS fine occupations with weight: $(nrow(cps_weights))")
    
    # Merge & compute coverage
    merged = leftjoin(cps_weights, dorn, on = :occ1990 => :occ1990_code)
    
    total_emp     = sum(merged.emp_weight)
    covered_mask  = .!ismissing.(merged.task_abstract)
    covered_emp   = sum(merged.emp_weight[covered_mask])
    coverage_rate = covered_emp / total_emp
    
    println("\n" * "="^60)
    println("COVERAGE DIAGNOSTICS")
    println("="^60)
    println(@sprintf("Total baseline employment   : %.2f million", total_emp/1e6))
    println(@sprintf("Covered by Dorn data        : %.2f million (%.1f%%)", covered_emp/1e6, coverage_rate*100))
    println(@sprintf("Uncovered occupations       : %d / %d", sum(.!covered_mask), nrow(merged)))
    
    # Coverage by broad group
    merged[!, :occ_group] = classify_occ1990.(merged.occ1990)
    merged = @subset(merged, .!ismissing.(:occ_group))
    
    grp_cov = combine(groupby(merged, :occ_group),
        :emp_weight => sum => :group_emp,
        :task_abstract => (t -> sum(.!ismissing.(t))) => :n_covered,
        :task_abstract => (t -> sum(ismissing.(t))) => :n_uncovered
    )
    grp_cov[!, :cov_rate] = grp_cov.group_emp ./ [sum(merged.emp_weight[merged.occ_group .== g]) for g in grp_cov.occ_group]
    
    println("\nCoverage by broad occupation group:")
    for r in eachrow(grp_cov)
        println(@sprintf("  %-30s | %5.1f%% covered | %3d occ covered | %3d uncovered", 
                         OCC_LABELS[r.occ_group], r.cov_rate*100, r.n_covered, r.n_uncovered))
    end
    
    if coverage_rate < 0.70
        @warn "Overall coverage <70%. Consider supplementing with O*NET or group-mean imputation."
    end
    
    # =============================================================================
    # 4. WEIGHTED AGGREGATION TO 9 BROAD GROUPS
    # =============================================================================
    println("\nAggregating task measures (employment-weighted)...")
    
    covered_df = @subset(merged, .!ismissing.(:task_abstract))
    
    broad_agg = combine(groupby(covered_df, :occ_group),
        [:task_abstract, :emp_weight] => ((v,w) -> sum(v.*w)/sum(w)) => :task_abstract,
        [:task_routine,  :emp_weight] => ((v,w) -> sum(v.*w)/sum(w)) => :task_routine,
        [:task_manual,   :emp_weight] => ((v,w) -> sum(v.*w)/sum(w)) => :task_manual,
        :emp_weight => sum => :total_employment
    )
    
    broad_agg[!, :occ_name] = [OCC_LABELS[g] for g in broad_agg.occ_group]
    select!(broad_agg, [:occ_group, :occ_name, :task_abstract, :task_routine, :task_manual, :total_employment])
    sort!(broad_agg, :occ_group)
    
    # Optional: Normalize so abstract+routine+manual ≈ 1 per group
    broad_agg[!, :task_sum] = broad_agg.task_abstract .+ broad_agg.task_routine .+ broad_agg.task_manual
    broad_agg[!, :task_abstract_norm] = broad_agg.task_abstract ./ broad_agg.task_sum
    broad_agg[!, :task_routine_norm]  = broad_agg.task_routine  ./ broad_agg.task_sum
    broad_agg[!, :task_manual_norm]   = broad_agg.task_manual   ./ broad_agg.task_sum
    
    println("\nAggregated task measures (9 broad groups):")
    show(broad_agg[:, Not(:total_employment)], allrows=true, allcols=true)
    
    # Save results
    CSV.write(joinpath(OUTPUT_DIR, "task_features_aggregated.csv"), broad_agg)
    println("\n✅ Saved to: $(joinpath(OUTPUT_DIR, "task_features_aggregated.csv"))")
    
    return broad_agg, coverage_rate
end

# =============================================================================
# 5. MAIN EXECUTION
# =============================================================================
function main()
    println("Starting Dorn task aggregation pipeline...")
    
    # 1. Parse CPS baseline weights
    cps_df = parse_cps_baseline_weights(CPS_FILE)
    
    # 2. Load Dorn, check coverage, aggregate
    features, coverage = load_and_aggregate_tasks(cps_df, DORN_FILE)
    
    println("\n" * "="^60)
    if coverage >= 0.85
        println("✅ Coverage adequate. Proceed with correlation analysis.")
    elseif coverage >= 0.70
        println("⚠️ Moderate coverage. Results are usable but note limitations in paper.")
    else
        println("❌ Low coverage. Strongly recommend O*NET supplementation or mean imputation.")
    end
    println("="^60)
    
    return features
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

main()