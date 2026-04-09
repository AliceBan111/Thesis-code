# =============================================================================
# 03_calc_baseline_shares.jl
# Calculate baseline occupation-industry employment shares for shift-share decomposition
# Following Bartik (1991) and Goldsmith-Pinkham, Sorkin, and Swift (2020)
# =============================================================================

using CSV, DataFrames, DataFramesMeta
using Dates, Statistics
using Printf

# =============================================================================
# 0. PATHS AND CONSTANTS
# =============================================================================
const DATA_DIR   = joinpath(@__DIR__, "../..", "data")
const OUTPUT_DIR = joinpath(@__DIR__, "../..", "result", "decomposition")
mkpath(OUTPUT_DIR)

const CPS_PATH = joinpath(DATA_DIR, "cps_00017.dat")

# Baseline period for share calculation (main specification)
const BASE_START = Date(1994, 1, 1)
const BASE_END   = Date(1996, 12, 1)

# Sample selection criteria (must match LP estimation sample)
const AGE_MIN = 16
const AGE_MAX = 64
const HOURS_MIN = 0.0
const HOURS_MAX = 105.0
const MIN_CELL_OBS = 50          # minimum observations per occ-ind cell
const MIN_CELL_WEIGHT = 1000.0   # minimum weighted employment per cell

# =============================================================================
# 1. DATE PARSING UTILITIES
# =============================================================================
"""
    parse_yearmonth(raw) -> Union{Date, Nothing}

Parse various date formats and return the first day of the month as Date.
Handles "1975M04", "1975m04", or Julia Date objects.
"""
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
# 2. OCCUPATION CLASSIFICATION (OCC1990 -> 9 major groups)
# =============================================================================
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

# =============================================================================
# 3. INDUSTRY CLASSIFICATION (IND1990 -> 13 major groups)
# =============================================================================
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
    13 => "Public_administration",
)

function classify_ind1990(ind::Union{Integer,Missing})::Union{Int,Missing}
    ismissing(ind) && return missing
    ind in 10:32   && return 1
    ind in 40:50   && return 2
    ind == 60      && return 3
    ind in 100:222 && return 4
    ind in 230:392 && return 5
    ind in 400:472 && return 6
    ind in 500:571 && return 7
    ind in 580:691 && return 8
    ind in 700:712 && return 9
    ind in 721:760 && return 10
    ind in 761:810 && return 11
    ind in 812:893 && return 12
    ind in 900:932 && return 13
    return missing
end

# =============================================================================
# 4. PARSE CPS FIXED-WIDTH FILE (same as 01_data_prep.jl)
# =============================================================================
function parse_cps(path::String)::DataFrame
    println("Parsing CPS fixed-width file: $path")
    println("File size: ", round(filesize(path)/1e9, digits=2), " GB")

    est_rows = max(1_000_000, round(Int, filesize(path) / 45))

    year_v     = Vector{Int32}(undef, est_rows)
    month_v    = Vector{Int32}(undef, est_rows)
    wtfinl_v   = Vector{Float32}(undef, est_rows)
    age_v      = Vector{Int32}(undef, est_rows)
    sex_v      = Vector{Int32}(undef, est_rows)
    marst_v    = Vector{Int32}(undef, est_rows)
    empstat_v  = Vector{Int32}(undef, est_rows)
    occ1990_v  = Vector{Int32}(undef, est_rows)
    ind1990_v  = Vector{Int32}(undef, est_rows)
    classwkr_v = Vector{Int32}(undef, est_rows)
    uhrsworkt_v = Vector{Float32}(undef, est_rows)
    earnwt_v   = Vector{Float32}(undef, est_rows)
    earnweek_v = Vector{Float32}(undef, est_rows)

    n = 0

    open(path, "r") do f
        for line in eachline(f)
            length(line) < 40 && continue

            year  = parse(Int32, @view line[1:4])
            month = parse(Int32, @view line[5:6])

            d = Date(year, month, 1)
            (d < Date(1983,4,1) || d > Date(2025,6,1)) && continue

            wtfinl_raw = tryparse(Float64, strip(@view line[7:20]))
            isnothing(wtfinl_raw) && continue

            earnwt_raw = tryparse(Float64, strip(@view line[38:47]))
            isnothing(earnwt_raw) && continue

            earnwk_str = strip(@view line[48:55])
            isempty(earnwk_str) && continue
            earnwk_raw = tryparse(Float64, earnwk_str)
            isnothing(earnwk_raw) && continue

            uhrsworkt_str = strip(@view line[35:37])
            isempty(uhrsworkt_str) && continue
            uhrsworkt_raw = tryparse(Float64, uhrsworkt_str)
            isnothing(uhrsworkt_raw) && continue

            n += 1
            if n > length(year_v)
                new_cap = round(Int, length(year_v) * 1.5)
                foreach(v -> resize!(v, new_cap),
                    (year_v, month_v, wtfinl_v, age_v, sex_v, marst_v,
                     empstat_v, occ1990_v, ind1990_v, classwkr_v,
                     uhrsworkt_v, earnwt_v, earnweek_v))
            end

            year_v[n]     = year
            month_v[n]    = month
            wtfinl_v[n]   = Float32(wtfinl_raw / 10_000.0)
            age_v[n]      = parse(Int32, @view line[21:22])
            sex_v[n]      = parse(Int32, @view line[23:23])
            marst_v[n]    = parse(Int32, @view line[24:24])
            empstat_v[n]  = parse(Int32, @view line[25:26])
            occ1990_v[n]  = parse(Int32, @view line[27:29])
            ind1990_v[n]  = parse(Int32, @view line[30:32])
            classwkr_v[n] = parse(Int32, @view line[33:34])
            uhrsworkt_v[n] = Float32(uhrsworkt_raw)
            earnwt_v[n]   = Float32(earnwt_raw / 10_000.0)
            earnweek_v[n] = Float32(earnwk_raw / 100.0)
        end
    end

    df = DataFrame(
        year      = year_v[1:n],
        month     = month_v[1:n],
        wtfinl    = wtfinl_v[1:n],
        age       = age_v[1:n],
        sex       = sex_v[1:n],
        marst     = marst_v[1:n],
        empstat   = empstat_v[1:n],
        occ1990   = occ1990_v[1:n],
        ind1990   = ind1990_v[1:n],
        classwkr  = classwkr_v[1:n],
        uhrsworkt = uhrsworkt_v[1:n],
        earnwt    = earnwt_v[1:n],
        earnweek  = earnweek_v[1:n],
    )

    println("  Raw rows parsed: ", nrow(df))
    return df
end

# =============================================================================
# 5. SAMPLE SELECTION FOR BASELINE SHARES
# =============================================================================
"""
    select_earnings_sample(df::DataFrame, date_start::Date, date_end::Date)

Apply the same sample selection as in LP estimation:
- Employment status: at work or has job not at work
- Age 16-64
- Valid earnings: earnweek > 0, not topcoded (9999.99)
- Apply topcode limits by year
- Valid weight: earnwt > 0
- Valid hours: 0 < uhrsworkt <= 105
- Remove self-employed (classwkr not in [21:28])
- Classify occupation and industry
- Date variable construction
"""
function select_earnings_sample(df::DataFrame, date_start::Date, date_end::Date)::DataFrame
    println("Selecting earnings sample for baseline period...")
    
    # Filter by date range
    df = @subset(df, date_start .<= Date.(:year, :month, 1) .<= date_end)
    
    # Employment: 10=at work, 12=has job not at work
    df = @subset(df, :empstat .∈ Ref([10, 12]))
    
    # Age restriction
    df = @subset(df, AGE_MIN .<= :age .<= AGE_MAX)
    
    # Valid earnings
    df = @subset(df, :earnweek .> 0, :earnweek .!= 9999.99)
    
    # Apply topcode limits by year (same as LP estimation)
    function topcode_limit(year::Integer)
        year <= 1988 ? 999.0f0 :
        year <= 1997 ? 1923.0f0 : 2884.0f0
    end
    df = @transform(df, :earnweek = min.(:earnweek, topcode_limit.(:year)))
    
    # Valid weight
    df = @subset(df, :earnwt .> 0)
    
    # Valid hours
    df = @subset(df, HOURS_MIN .< :uhrsworkt .<= HOURS_MAX)
    
    # Remove self-employed (keep only wage/salary workers)
    df = @subset(df, :classwkr .∈ Ref([21, 22, 23, 24, 25, 27, 28]))
    
    # Classify occupation and industry
    df[!, :occ_group] = map(classify_occ1990, df.occ1990)
    df[!, :ind_group] = map(classify_ind1990, df.ind1990)
    
    # Drop unclassified
    df = @subset(df, .!ismissing.(:occ_group), .!ismissing.(:ind_group))
    df[!, :occ_group] = convert(Vector{Int}, df.occ_group)
    df[!, :ind_group] = convert(Vector{Int}, df.ind_group)
    
    # Date variable
    df[!, :date] = Date.(df.year, df.month, 1)
    
    println("  Sample size after selection: ", nrow(df))
    return df
end

# =============================================================================
# 6. CALCULATE OCCUPATION-INDUSTRY SHARES
# =============================================================================
"""
    calc_occ_ind_shares(df::DataFrame)

Calculate weighted employment shares s_ij = emp_ij / emp_i for each occupation i.
Uses earnwt as sampling weight.
Returns a DataFrame with columns: occ_group, occ_name, ind_group, ind_name, share
"""
function calc_occ_ind_shares(df::DataFrame)::DataFrame
    println("Calculating occupation-industry employment shares...")
    
    # Step 1: Calculate weighted employment for each occ-ind cell
    cell_emp = combine(groupby(df, [:occ_group, :ind_group]),
        :earnwt => sum => :cell_weight,
        nrow => :n_obs
    )
    
    # Step 2: Calculate total employment for each occupation
    occ_emp = combine(groupby(df, :occ_group),
        :earnwt => sum => :occ_total_weight
    )
    
    # Step 3: Merge and calculate shares
    shares = leftjoin(cell_emp, occ_emp, on = :occ_group)
    shares[!, :share] = shares.cell_weight ./ shares.occ_total_weight
    
    # Step 4: Add labels
    shares[!, :occ_name] = [OCC_LABELS[o] for o in shares.occ_group]
    shares[!, :ind_name] = [IND_LABELS[i] for i in shares.ind_group]
    
    # Step 5: Reorder columns
    select!(shares, [:occ_group, :occ_name, :ind_group, :ind_name, :share, :cell_weight, :n_obs])
    
    # Step 6: Sort for readability
    sort!(shares, [:occ_group, :ind_group])
    
    println("  Total occ-ind cells: ", nrow(shares))
    return shares
end

# =============================================================================
# 7. VALIDATION AND DIAGNOSTICS
# =============================================================================
"""
    validate_shares(shares::DataFrame)

Perform validation checks:
1. Row sums: Σⱼ s_ij = 1 for each occupation i
2. Share range: 0 ≤ s_ij ≤ 1
3. Small cell flag: mark cells with n_obs < MIN_CELL_OBS or weight < MIN_CELL_WEIGHT
"""
function validate_shares(shares::DataFrame)
    println("\n" * "="^60)
    println("VALIDATION CHECKS")
    println("="^60)
    
    # Check 1: Row sums = 1
    row_sums = combine(groupby(shares, :occ_group), :share => sum => :row_sum)
    max_deviation = maximum(abs.(row_sums.row_sum .- 1.0))
    println("Max deviation from row sum = 1: $(round(max_deviation, digits=8))")
    @assert max_deviation < 1e-4 "⚠️ Row sums not equal to 1! Check weight calculation."
    
    # Check 2: Share range
    invalid_range = @subset(shares, :share .< 0 .|| :share .> 1)
    println("Cells with share outside [0,1]: $(nrow(invalid_range))")
    @assert nrow(invalid_range) == 0 "⚠️ Invalid share values detected."
    
    # Check 3: Small cells
    shares[!, :small_cell] = (shares.n_obs .< MIN_CELL_OBS) .| 
                             (shares.cell_weight .< MIN_CELL_WEIGHT)
    n_small = sum(shares.small_cell)
    println("Small cells (n_obs < $MIN_CELL_OBS or weight < $MIN_CELL_WEIGHT): $n_small / $(nrow(shares))")
    
    # Check 4: Summary by occupation
    println("\nSummary by occupation:")
    occ_summary = combine(groupby(shares, :occ_group),
        :share => (s -> count(s .> 0.05)) => :n_major_industries,  # industries with >5% share
        :share => maximum => :max_share,
        :small_cell => sum => :n_small_cells
    )
    occ_summary[!, :occ_name] = [OCC_LABELS[o] for o in occ_summary.occ_group]
    
    for row in eachrow(occ_summary)
        println("  $(row.occ_name): $(row.n_major_industries) major industries, " *
                "max share=$(round(row.max_share*100, digits=1))%, " *
                "$(row.n_small_cells) small cells")
    end
    
    return shares
end

# =============================================================================
# 8. ROBUSTNESS: ALTERNATIVE BASELINE PERIODS
# =============================================================================
"""
    calc_shares_robustness(cps_raw::DataFrame)

Calculate shares using alternative baseline periods for robustness checks.
Returns a Dict with results for each specification.
"""
function calc_shares_robustness(cps_raw::DataFrame, main_shares::DataFrame)
    println("\n" * "="^60)
    println("ROBUSTNESS: ALTERNATIVE BASELINE PERIODS")
    println("="^60)
    
    results = Dict{String, DataFrame}()
    
    # Specification 1: Extended baseline (1994-1998)
    println("\nCalculating shares for 1994-1998...")
    df_ext = select_earnings_sample(cps_raw, Date(1994,1,1), Date(1998,12,1))
    shares_ext = calc_occ_ind_shares(df_ext)
    shares_ext = validate_shares(shares_ext)
    results["1994-1998"] = shares_ext
    
    # Specification 2: Later baseline (1997-1999)
    println("\nCalculating shares for 1997-1999...")
    df_late = select_earnings_sample(cps_raw, Date(1997,1,1), Date(1999,12,1))
    shares_late = calc_occ_ind_shares(df_late)
    shares_late = validate_shares(shares_late)
    results["1997-1999"] = shares_late
    
    # Check correlation between specifications and MAIN baseline
    println("\nCorrelation of shares vs main baseline (1994-1996):")
    for (spec, df) in results
        # Merge main shares with current spec on occ_group + ind_group
        merged = innerjoin(
            select(main_shares, :occ_group, :ind_group, :share => :share_main),
            select(df, :occ_group, :ind_group, :share => :share_alt),
            on = [:occ_group, :ind_group]
        )
        if nrow(merged) > 0
            corr = cor(merged.share_main, merged.share_alt)
            println("  $(spec) vs main: $(round(corr, digits=4))")
        else
            println("  $(spec) vs main: N/A (no overlapping cells)")
        end
    end
    
    return results
end

# =============================================================================
# 9. OUTPUT FUNCTIONS
# =============================================================================
"""
    save_shares(shares::DataFrame, suffix::String)

Save shares to CSV with optional suffix for robustness specifications.
"""
function save_shares(shares::DataFrame, suffix::String = "")
    filename = suffix == "" ? "occ_ind_shares.csv" : "occ_ind_shares_$(suffix).csv"
    filepath = joinpath(OUTPUT_DIR, filename)
    
    # Prepare output: keep essential columns
    output = select(shares, :occ_group, :occ_name, :ind_group, :ind_name, :share, :n_obs)
    CSV.write(filepath, output)
    println("Saved: $filepath")
end

"""
    save_summary(shares::DataFrame, robustness_results::Dict)

Save a summary table comparing main and robustness specifications.
"""
function save_summary(shares::DataFrame, robustness_results::Dict)
    summary_rows = []
    
    for occ in 1:9
        main_sub = @subset(shares, :occ_group .== occ)
        row = (
            occ_group = occ,
            occ_name = OCC_LABELS[occ],
            n_industries = count(main_sub.share .> 0),
            top3_share = sum(sort(main_sub.share, rev=true)[1:min(3, nrow(main_sub))]),
            herfindahl = sum(main_sub.share .^ 2),  # concentration measure
        )
        push!(summary_rows, row)
    end
    
    summary_df = DataFrame(summary_rows)
    filepath = joinpath(OUTPUT_DIR, "occ_shares_summary.csv")
    CSV.write(filepath, summary_df)
    println("Saved summary: $filepath")
end

# =============================================================================
# 10. MAIN EXECUTION
# =============================================================================
function main()
    println("\n" * "="^60)
    println("BASELINE OCCUPATION-INDUSTRY SHARES CALCULATION")
    println("Baseline period: $(BASE_START) to $(BASE_END)")
    println("="^60 * "\n")
    
    # Step 1: Load CPS data
    cps_raw = parse_cps(CPS_PATH)
    
    # Step 2: Select baseline sample
    cps_base = select_earnings_sample(cps_raw, BASE_START, BASE_END)
    
    # Step 3: Calculate shares
    shares = calc_occ_ind_shares(cps_base)
    
    # Step 4: Validate
    shares = validate_shares(shares)
    
    # Step 5: Save main results
    save_shares(shares, "1994-1996")
    
    # Step 6: Robustness checks (FIXED: pass main_shares as argument)
    robustness = calc_shares_robustness(cps_raw, shares)  # ← 传入 shares
    for (spec, df) in robustness
        save_shares(df, spec)
    end
    
    # Step 7: Save summary
    save_summary(shares, robustness)
    
    println("\n" * "="^60)
    println("CALCULATION COMPLETE")
    println("Output directory: $OUTPUT_DIR")
    println("="^60)
    
    return shares, robustness
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

main()

