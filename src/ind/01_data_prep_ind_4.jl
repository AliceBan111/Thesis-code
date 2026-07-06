# =============================================================================
# data_prep_ind_4.jl
# Oil Supply News Shock × Industry Heterogeneity + OilShare Exposure
#
# Based on the dat column layout and parsing logic from data_prep_occ.jl (file 1)
# Rewritten on top of the industry grouping framework from data_prep_ind (file 3)
#
# Added:
#   - OilShare (Bartik exposure) construction:
#       Exposure_o = Σ_i Share_{o,i} × OilIntensity_i
#     where Share_{o,i} is the share of employment in industry group i within occupation group o
#     (This file runs at the industry level; OilShare is used later as a cross-sectional feature)
#   - 4-group industry analysis: identify by the 12-group crosswalk first, then merge to the narrow industry classification
#
# Usage:
#   include("data_prep_ind_4.jl")
#   panel = main(:hourly_rate)
#   # The OilShare vector (12 industry groups) is defined in OIL_INTENSITY
# =============================================================================

using HTTP, Downloads, SHA
using CSV, DataFrames, DataFramesMeta
using XLSX
using Dates, Statistics
using Printf

# =============================================================================
# 0. PATHS
# =============================================================================
const DATA_DIR = joinpath(@__DIR__, "../..", "data")

# ---- Use the same CPS file as the occ version (includes both occ1990 and ind1990 columns) ----------
# Column layout for file 1 (cps_00018.dat):
#   YEAR       1-4
#   MONTH      5-6
#   WTFINL     7-20  (4 implied decimals)
#   EARNWEEK2  21-28 (2 implied decimals)
#   AGE        29-30
#   SEX        31-31
#   MARST      32-32
#   EMPSTAT    33-34
#   OCC1990    35-37
#   IND1990    38-40
#   CLASSWKR   41-42
#   UHRSWORKT  43-45
#   AHRSWORKT  46-48
#   EARNWT     49-58  (4 implied decimals)
const CPS_FILE_ID = "1bz1vYmddvubhLfp_TQ2x83hdUZ3NAZfu"
const CPS_SHA256  = "67CFDB7EE29C81587F8C7F269328FB8FE43E58D9868C5F6AB98C5139575A5D4E"
const CPS_PATH    = joinpath(DATA_DIR, "cps_00018.dat")

const DATE_START = Date(1983, 4, 1)
const DATE_END   = Date(2025, 6, 1)
const L_LAG      = 12

# =============================================================================
# 1. OIL INTENSITY (1997 BEA I-O Table, 12 industry groups)
#    Source: user-provided, based on the ind1990 → 12-group crosswalk
#    Group 13 (Public_administration) has no matching BEA data and is dropped
# =============================================================================
const OIL_INTENSITY = Dict(
    1  => 0.0365,   # Agriculture_forestry_fishing
    2  => 0.8997,   # Mining
    3  => 0.0387,   # Construction
    4  => 0.0819,   # Manufacturing_nondurable 
    5  => 0.1789,   # Manufacturing_durable
    6  => 0.6012,   # Transportation_utilities
    7  => 0.0023,   # Wholesale_trade
    8  => 0.0069,   # Retail_trade
    9  => 0.0265,   # Finance_insurance_realestate
    10 => 0.0143,   # Business_repair_services
    11 => 0.0283,   # Personal_entertainment_services
    12 => 0.0210,   # Professional_related_services
)

# =============================================================================
# 2. INDUSTRY & OCCUPATION CLASSIFICATION
# =============================================================================
const IND_MERGE = Dict(
    1  => 4,  # Agriculture_forestry_fishing
    2  => 1,  # Mining
    3  => 2,  # Construction
    4  => 2,  # Manufacturing_nondurable
    5  => 2,  # Manufacturing_durable
    6  => 1,  # Transportation_utilities
    7  => 3,  # Wholesale_trade
    8  => 3,  # Retail_trade
    9  => 4,  # Finance_insurance_realestate
    10 => 4,  # Business_repair_services
    11 => 4,  # Personal_entertainment_services
    12 => 4,  # Professional_related_services
)

const IND_LABELS = Dict(
    1 => "Energy_intensive",
    2 => "Manufacturing_Construction",
    3 => "Trade",
    4 => "Services",
)

const IND_LABELS_12 = Dict(
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
)

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

function classify_ind1990_12(ind::Union{Integer,Missing})::Union{Int,Missing}
    ismissing(ind) && return missing
    ind in 10:32   && return 1
    (ind in 40:50 || ind in 200:201)   && return 2
    ind == 60      && return 3
    (ind in 100:199 || ind in 202:222) && return 4
    ind in 230:392 && return 5
    ind in 400:472 && return 6
    ind in 500:571 && return 7
    ind in 580:691 && return 8
    ind in 700:712 && return 9
    ind in 721:760 && return 10
    ind in 761:810 && return 11
    ind in 812:893 && return 12
    # ind in 900:932 → group 13, dropped; return missing
    return missing
end

function classify_ind1990(ind::Union{Integer,Missing})::Union{Int,Missing}
    g12 = classify_ind1990_12(ind)
    ismissing(g12) && return missing
    return get(IND_MERGE, g12, missing)
end

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
# 3. OUTPUT DIRECTORY
# =============================================================================
function get_output_dir(variant::Symbol)::String
    base = joinpath(@__DIR__, "../..", "result", "ind_4")
    dir  = joinpath(base, string(variant))
    mkpath(dir)
    return dir
end

# =============================================================================
# 4. DOWNLOAD CPS
# =============================================================================
function download_gdrive_large(file_id::String, dest::String;
                                expected_hash::Union{String,Nothing}=nothing)
    if isfile(dest)
        first_bytes = String(read(dest, min(20, filesize(dest))))
        if startswith(first_bytes, "<!") || startswith(first_bytes, "<h")
            println("Existing file is HTML (bad download) -- deleting and retrying.")
            rm(dest)
        else
            println("File already exists: $dest -- skipping download.")
            return
        end
    end
    println("Downloading CPS data from Google Drive...")
    base_url = "https://drive.google.com/uc?export=download&id=$(file_id)"
    resp     = HTTP.get(base_url)
    body     = String(resp.body)
    uuid_m   = match(r"name=\"uuid\" value=\"([^\"]+)\"", body)
    real_url = if !isnothing(uuid_m)
        "https://drive.usercontent.google.com/download?id=$(file_id)&export=download&confirm=t&uuid=$(uuid_m.captures[1])"
    else
        "https://drive.usercontent.google.com/download?id=$(file_id)&export=download&confirm=t"
    end
    Downloads.download(real_url, dest)
    println("Download complete.")
    if !isnothing(expected_hash)
        actual = bytes2hex(open(sha256, dest))
        lowercase(actual) != lowercase(expected_hash) &&
            error("Hash mismatch!\n  expected: $expected_hash\n  actual:   $actual")
        println("Hash verified ✓")
    end
end

# =============================================================================
# 5. DATE PARSING
# =============================================================================
function parse_yearmonth(raw)::Union{Date,Nothing}
    raw isa Date && return raw
    s = strip(string(raw))
    m = match(r"^(\d{4})[Mm](\d{2})$", s)
    isnothing(m) && return nothing
    yr = parse(Int, m.captures[1])
    mo = parse(Int, m.captures[2])
    return Date(yr, mo, 1)
end

# =============================================================================
# 6. PARSE CPS FIXED-WIDTH FILE
#    Use the column layout from file 1 (cps_00018.dat), reading both OCC1990 and IND1990
# =============================================================================
function parse_cps(path::String)::DataFrame
    println("Parsing CPS fixed-width file (occ+ind layout): $path")
    println("File size: ", round(filesize(path)/1e9, digits=2), " GB")

    est_rows = max(1_000_000, round(Int, filesize(path) / 48))

    year_v       = Vector{Int32}(undef, est_rows)
    month_v      = Vector{Int32}(undef, est_rows)
    wtfinl_v     = Vector{Float32}(undef, est_rows)
    earnweek2_v  = Vector{Float32}(undef, est_rows)
    age_v        = Vector{Int32}(undef, est_rows)
    sex_v        = Vector{Int32}(undef, est_rows)
    marst_v      = Vector{Int32}(undef, est_rows)
    empstat_v    = Vector{Int32}(undef, est_rows)
    occ1990_v    = Vector{Int32}(undef, est_rows)
    ind1990_v    = Vector{Int32}(undef, est_rows)
    classwkr_v   = Vector{Int32}(undef, est_rows)
    hours_v      = Vector{Float32}(undef, est_rows)
    earnwt_v     = Vector{Float32}(undef, est_rows)

    n = 0

    open(path, "r") do f
        for line in eachline(f)
            length(line) < 50 && continue

            year  = parse(Int32, @view line[1:4])
            month = parse(Int32, @view line[5:6])
            d = Date(year, month, 1)
            (d < DATE_START || d > DATE_END) && continue

            wtfinl_raw = tryparse(Float64, strip(@view line[7:20]))
            isnothing(wtfinl_raw) && continue

            earnwt_raw = tryparse(Float64, strip(@view line[49:58]))
            isnothing(earnwt_raw) && continue

            earnwk_str = strip(@view line[21:28])
            isempty(earnwk_str) && continue
            earnwk_raw = tryparse(Float64, earnwk_str)
            isnothing(earnwk_raw) && continue

            # Hours: pre-1994 use AHRSWORKT (46-48), post use UHRSWORKT (43-45)
            hrs_raw = nothing
            if year < 1994
                hrs_str = strip(@view line[46:48])
                !isempty(hrs_str) && (hrs_raw = tryparse(Float64, hrs_str))
            else
                hrs_str = strip(@view line[43:45])
                !isempty(hrs_str) && (hrs_raw = tryparse(Float64, hrs_str))
            end
            isnothing(hrs_raw) && continue

            n += 1
            if n > length(year_v)
                new_cap = round(Int, length(year_v) * 1.5)
                foreach(v -> resize!(v, new_cap),
                    (year_v, month_v, wtfinl_v, earnweek2_v, age_v, sex_v, marst_v,
                     empstat_v, occ1990_v, ind1990_v, classwkr_v, hours_v, earnwt_v))
            end

            year_v[n]      = year
            month_v[n]     = month
            wtfinl_v[n]    = Float32(wtfinl_raw / 10_000.0)
            earnweek2_v[n] = Float32(earnwk_raw / 100.0)
            age_v[n]       = parse(Int32, @view line[29:30])
            sex_v[n]       = parse(Int32, @view line[31:31])
            marst_v[n]     = parse(Int32, @view line[32:32])
            empstat_v[n]   = parse(Int32, @view line[33:34])
            occ1990_v[n]   = parse(Int32, @view line[35:37])
            ind1990_v[n]   = parse(Int32, @view line[38:40])
            classwkr_v[n]  = parse(Int32, @view line[41:42])
            hours_v[n]     = Float32(hrs_raw)
            earnwt_v[n]    = Float32(earnwt_raw / 10_000.0)
        end
    end

    df = DataFrame(
        year         = year_v[1:n],
        month        = month_v[1:n],
        wtfinl       = wtfinl_v[1:n],
        earnweek     = earnweek2_v[1:n],
        age          = age_v[1:n],
        sex          = sex_v[1:n],
        marst        = marst_v[1:n],
        empstat      = empstat_v[1:n],
        occ1990      = occ1990_v[1:n],
        ind1990      = ind1990_v[1:n],
        classwkr     = classwkr_v[1:n],
        hours_worked = hours_v[1:n],
        earnwt       = earnwt_v[1:n],
    )
    println("  Raw rows in sample period: ", nrow(df))
    return df
end

# =============================================================================
# 7. OILSHARE CONSTRUCTION (Bartik exposure)
#
#   Exposure_o = Σ_i Share_{o,i} × OilIntensity_i
#
#   Here o = occ_group (9 groups), i = ind_group (4 groups)
#   Share_{o,i} = weighted share of employment in industry group i within occupation group o
#
#   Returns DataFrame: occ_group, oil_exposure
#   Also returns ind_group-level oil_intensity for the cross-sectional analysis in the ind LP
# =============================================================================
function build_oilshare(cps_raw::DataFrame)::Tuple{DataFrame, DataFrame}
    println("Building OilShare (Bartik exposure)...")

    # Use only employed workers
    df = @subset(cps_raw, :empstat .∈ Ref([10, 12]))
    df = @subset(df, 15 .<= :age .<= 64)
    df = @subset(df, :wtfinl .> 0)
    df = @subset(df, :classwkr .∈ Ref([21, 22, 23, 24, 25, 27, 28]))

    df[!, :occ_group] = map(classify_occ1990, df.occ1990)
    df[!, :ind_group12] = map(classify_ind1990_12, df.ind1990)

    # Drop rows with missing occ or ind (including ind group 13)
    df = @subset(df, .!ismissing.(:occ_group), .!ismissing.(:ind_group12))
    df[!, :occ_group] = convert(Vector{Int}, df.occ_group)
    df[!, :ind_group12] = convert(Vector{Int}, df.ind_group12)
    df[!, :ind_group] = [IND_MERGE[g] for g in df.ind_group12]

    # 4-group OilIntensity: CPS employment-weighted average of underlying 12 groups.
    ind12_weights = combine(
        groupby(df, [:ind_group, :ind_group12]),
        :wtfinl => sum => :emp_weight
    )
    ind12_weights[!, :oil_intensity12] = [OIL_INTENSITY[g] for g in ind12_weights.ind_group12]
    ind_intensity_df = combine(
        groupby(ind12_weights, :ind_group),
        [:oil_intensity12, :emp_weight] =>
            ((oi, w) -> sum(oi .* w) / sum(w)) => :oil_intensity
    )
    sort!(ind_intensity_df, :ind_group)

    # Build a stable occ×ind employment matrix over the full sample (time-invariant, using full-sample weights)
    crosstab = combine(
        groupby(df, [:occ_group, :ind_group]),
        :wtfinl => sum => :emp_weight
    )

    # Total employment weight for each occ_group
    occ_total = combine(groupby(crosstab, :occ_group),
                        :emp_weight => sum => :total_weight)
    crosstab = leftjoin(crosstab, occ_total, on = :occ_group)
    crosstab[!, :share] = crosstab.emp_weight ./ crosstab.total_weight

    # Attach the merged 4-group OilIntensity
    crosstab = leftjoin(crosstab, ind_intensity_df, on = :ind_group)

    # Bartik exposure: Σ_i share_{o,i} × OilIntensity_i
    oilshare_df = combine(
        groupby(crosstab, :occ_group),
        [:share, :oil_intensity] =>
            ((s, oi) -> sum(s .* oi)) => :oil_exposure
    )
    sort!(oilshare_df, :occ_group)

    println("  OilShare by occupation group:")
    for row in eachrow(oilshare_df)
        label = get(OCC_LABELS, row.occ_group, "group_$(row.occ_group)")
        @printf("    Group %2d %-35s %.4f\n", row.occ_group, label, row.oil_exposure)
    end

    return oilshare_df, ind_intensity_df
end

# =============================================================================
# 8. SAMPLE SELECTION
# =============================================================================
function clean_cps_labor(df::DataFrame)::DataFrame
    println("Cleaning CPS (labor sample, ind_group)...")
    df = @subset(df, :empstat .∈ Ref([10, 12, 20, 21, 22]))
    df = @subset(df, 15 .<= :age .<= 64)
    df = @subset(df, :wtfinl .> 0)
    df = @subset(df, :classwkr .∈ Ref([21, 22, 23, 24, 25, 27, 28]))

    df[!, :ind_group] = map(classify_ind1990, df.ind1990)
    df = @subset(df, .!ismissing.(:ind_group))   # Automatically excludes group 13
    df[!, :ind_group] = convert(Vector{Int}, df.ind_group)
    df[!, :date]      = Date.(df.year, df.month, 1)
    println("  Labor sample size: ", nrow(df))
    return df
end

function clean_cps_earn(df::DataFrame)::DataFrame
    println("Cleaning CPS (earnings sample, ind_group)...")
    df = @subset(df, :empstat .∈ Ref([10, 12]))
    df = @subset(df, 15 .<= :age .<= 64)
    df = @subset(df, :earnweek .> 0, :earnweek .!= 999999.99)

    function topcode_limit(year::Integer)
        year <= 1988 ? 999.0 :
        year <= 1997 ? 1923.0 : 2884.61
    end
    df = @transform(df, :earnweek = min.(:earnweek, topcode_limit.(:year)))
    df = @subset(df, :earnwt .> 0)
    df = @subset(df, 0 .< :hours_worked .<= 105)
    df = @subset(df, :classwkr .∈ Ref([21, 22, 23, 24, 25, 27, 28]))

    df[!, :ind_group] = map(classify_ind1990, df.ind1990)
    df = @subset(df, .!ismissing.(:ind_group))
    df[!, :ind_group] = convert(Vector{Int}, df.ind_group)
    df[!, :female]    = Int.(df.sex .== 2)
    df[!, :married]   = Int.(df.marst .∈ Ref([1, 2]))
    df[!, :date]      = Date.(df.year, df.month, 1)
    println("  Earnings sample size: ", nrow(df))
    return df
end

# =============================================================================
# 9. BUILD UNEMPLOYMENT PANEL
# =============================================================================
function build_panel_unemp(cps::DataFrame)::DataFrame
    gdf = groupby(cps, [:ind_group, :date, :year, :month])
    panel = combine(gdf) do sdf
        emp = sdf.empstat
        wt  = sdf.wtfinl
        unemployed    = sum((emp .∈ Ref((20,21,22))) .* wt)
        laborforce    = sum((emp .∈ Ref((10,12,20,21,22))) .* wt)
        employed      = sum((emp .∈ Ref((10,12))) .* wt)
        unemp_rate    = laborforce == 0 ? missing : unemployed / laborforce
        log_emp_count = employed > 0 ? log(employed) : missing
        pop_weight    = sum(wt)
        (; unemp_rate, log_emp_count, pop_weight)
    end
    panel = @subset(panel, :pop_weight .>= 1000)
    return sort(panel, [:ind_group, :date])
end

# =============================================================================
# 10. LOAD CPI AND DEFLATE
# =============================================================================
function load_cpi()::DataFrame
    println("Loading CPI from VARdata.xlsx...")
    path = joinpath(DATA_DIR, "VARdata.xlsx")
    xf   = XLSX.readxlsx(path)
    sh   = xf["Monthly"]
    dates    = Vector{Date}()
    cpi_vals = Vector{Float64}()
    for row in XLSX.eachrow(sh)
        XLSX.row_number(row) == 1 && continue
        d_raw = row[1]; c_raw = row[7]
        (ismissing(d_raw) || ismissing(c_raw)) && continue
        d = parse_yearmonth(d_raw)
        isnothing(d) && continue
        push!(dates, d); push!(cpi_vals, Float64(c_raw))
    end
    cpi_df = DataFrame(date = dates, cpi = cpi_vals)
    cpi_df = @subset(cpi_df, DATE_START .<= :date .<= DATE_END)
    sort!(cpi_df, :date)
    println("  CPI rows: ", nrow(cpi_df))
    return cpi_df
end

function deflate_wages!(cps::DataFrame, cpi::DataFrame)::DataFrame
    cps = leftjoin(cps, select(cpi, :date, :cpi), on = :date)
    cps = @subset(cps, .!ismissing.(:cpi))
    cps[!, :real_earnweek] = cps.earnweek ./ cps.cpi
    cps[!, :log_rincome]   = log.(cps.real_earnweek)
    return cps
end

function add_log_rwage!(cps::DataFrame, variant::Symbol)::DataFrame
    cps[!, :log_rwage] = cps.log_rincome .- log.(cps.hours_worked)
    if variant == :hours
        cps[!, :log_hours] = log.(cps.hours_worked)
    end
    return cps
end

# =============================================================================
# 11. WEIGHTED QUANTILE
# =============================================================================
function weighted_quantile(x::AbstractVector, w::AbstractVector, q::Float64)::Float64
    idx  = sortperm(x)
    xs, ws = x[idx], w[idx]
    cumw = cumsum(ws); cumw ./= cumw[end]
    i = searchsortedfirst(cumw, q)
    return xs[clamp(i, 1, length(xs))]
end

# =============================================================================
# 12. BUILD EARNINGS PANEL (variant-specific)
# =============================================================================
function build_panel_earn(cps::DataFrame, variant::Symbol)::DataFrame
    gdf = groupby(cps, [:ind_group, :date, :year, :month])

    panel = if variant == :hourly_rate
        combine(gdf,
            [:log_rincome, :earnwt] => ((w,wt) -> sum(w.*wt)/sum(wt)) => :log_rincome,
            [:log_rwage,   :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :log_rwage,
            [:age,         :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :age_mean,
            [:female,      :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :female_share,
            [:married,     :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :married_share,
            [:hours_worked,:earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :hours_mean,
            :log_rincome => length => :n_obs,
        )
    elseif variant == :hours
        combine(gdf,
            [:log_rincome, :earnwt] => ((w,wt) -> sum(w.*wt)/sum(wt)) => :log_rincome,
            [:log_rwage,   :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :log_rwage,
            [:log_hours,   :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :log_hours,
            [:age,         :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :age_mean,
            [:female,      :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :female_share,
            [:married,     :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :married_share,
            :log_rincome => length => :n_obs,
        )
    elseif variant == :income
        combine(gdf,
            [:log_rincome, :earnwt] => ((w,wt) -> sum(w.*wt)/sum(wt)) => :log_rincome,
            [:log_rwage,   :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :log_rwage,
            [:age,         :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :age_mean,
            [:female,      :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :female_share,
            [:married,     :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :married_share,
            [:hours_worked,:earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :hours_mean,
            :log_rincome => length => :n_obs,
        )
    elseif variant == :income_share_var
        p = combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt))  => :log_rincome,
            [:log_rwage,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :log_rwage,
            [:age,         :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :age_mean,
            [:female,      :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :female_share,
            [:married,     :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :married_share,
            [:hours_worked,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :hours_mean,
            [:earnweek,    :earnwt] => ((e, wt) -> sum(e .* wt))             => :total_income_weekly,
            [:earnwt]               => sum                                    => :n_employed,
            :log_rincome => length => :n_obs,
        )
        p = @subset(p, :n_obs .>= 30)
        monthly_agg = combine(groupby(p, :date),
            :total_income_weekly => sum => :total_income_all,
            :n_employed          => sum => :total_emp_all)
        p = leftjoin(p, monthly_agg, on = :date)
        p[!, :income_share] = p.total_income_weekly ./ p.total_income_all
        check = combine(groupby(p, :date), :income_share => sum => :sum_share)
        @assert all(isapprox.(check.sum_share, 1.0, atol=1e-6))
        p

    elseif variant == :inequality
        p = combine(gdf,
            [:log_rincome, :earnwt] => ((w,wt) -> sum(w.*wt)/sum(wt))         => :log_rincome,
            [:log_rwage,   :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt))         => :log_rwage,
            [:log_rincome, :earnwt] => ((x,wt) -> weighted_quantile(x,wt,0.25)) => :log_rincome_p25,
            [:log_rincome, :earnwt] => ((x,wt) -> weighted_quantile(x,wt,0.75)) => :log_rincome_p75,
            [:age,         :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt))         => :age_mean,
            [:female,      :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt))         => :female_share,
            [:married,     :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt))         => :married_share,
            [:hours_worked,:earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt))         => :hours_mean,
            :log_rincome => length => :n_obs,
        )
        p = @subset(p, :n_obs .>= 30)
        p[!, :log_ratio_7525] = p.log_rincome_p75 .- p.log_rincome_p25
        p
    elseif variant == :median
        combine(gdf,
            [:log_rincome, :earnwt] => ((w,wt) -> weighted_quantile(w,wt,0.5)) => :log_rincome_p50,
            [:log_rwage,   :earnwt] => ((x,wt) -> weighted_quantile(x,wt,0.5)) => :log_rwage_p50,
            [:age,         :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt))          => :age_mean,
            [:female,      :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt))          => :female_share,
            [:married,     :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt))          => :married_share,
            [:hours_worked,:earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt))          => :hours_mean,
            :log_rincome => length => :n_obs,
        )
    elseif variant in (:unemployment, :employment)
        combine(gdf,
            [:log_rincome, :earnwt] => ((w,wt) -> sum(w.*wt)/sum(wt)) => :log_rincome,
            [:log_rwage,   :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :log_rwage,
            [:age,         :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :age_mean,
            [:female,      :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :female_share,
            [:married,     :earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :married_share,
            [:hours_worked,:earnwt] => ((x,wt) -> sum(x.*wt)/sum(wt)) => :hours_mean,
            :log_rincome => length => :n_obs,
        )
    else
        error("Unknown variant in build_panel_earn: $variant")
    end

    if !(variant in (:inequality,))
        panel = @subset(panel, :n_obs .>= 30)
    end
    println("  Panel observations: ", nrow(panel))
    return sort(panel, [:ind_group, :date])
end

# =============================================================================
# 13. LOAD MACRO DATA
# =============================================================================
function load_oil_shock()::DataFrame
    println("Loading oil supply news shock...")
    path = joinpath(DATA_DIR, "oilSupplyNewsShocks_2025M06.xlsx")
    xf   = XLSX.readxlsx(path)
    sh   = xf["Monthly"]
    dates = Vector{Date}(); shocks = Vector{Float64}()
    for row in XLSX.eachrow(sh)
        XLSX.row_number(row) == 1 && continue
        d_raw = row[1]; s_raw = row[3]
        (ismissing(d_raw) || ismissing(s_raw)) && continue
        d = parse_yearmonth(d_raw); isnothing(d) && continue
        push!(dates, d); push!(shocks, Float64(s_raw))
    end
    df = DataFrame(date=dates, shock=shocks)
    df = @subset(df, DATE_START .<= :date .<= DATE_END)
    sort!(df, :date); println("  Shock rows: ", nrow(df))
    return df
end

function load_oil_price()::DataFrame
    println("Loading oil price (WTI)...")
    path = joinpath(DATA_DIR, "VARdata.xlsx")
    xf   = XLSX.readxlsx(path); sh = xf["Monthly"]
    dates = Vector{Date}(); prices = Vector{Float64}()
    for row in XLSX.eachrow(sh)
        XLSX.row_number(row) == 1 && continue
        d_raw = row[1]; p_raw = row[2]
        (ismissing(d_raw) || ismissing(p_raw)) && continue
        d = parse_yearmonth(d_raw); isnothing(d) && continue
        push!(dates, d); push!(prices, Float64(p_raw))
    end
    df = DataFrame(date=dates, oil_price=prices)
    df[!, :log_oil_price] = log.(df.oil_price)
    df = @subset(df, DATE_START .<= :date .<= DATE_END)
    sort!(df, :date); return df
end

function load_fred(filename::String, colname::Symbol)::DataFrame
    path = joinpath(DATA_DIR, filename)
    df = CSV.read(path, DataFrame; types=Dict(1=>Date,2=>Float32), dateformat="yyyy-m-d")
    rename!(df, first(names(df)) => :date, names(df)[2] => colname)
    filter!(row -> DATE_START <= row.date <= DATE_END, df)
    sort!(df, :date); return df
end

# =============================================================================
# 14. BUILD MACRO PANEL
# =============================================================================
function lag_vec(v::AbstractVector, k::Int)
    T   = Union{eltype(v), Missing}
    out = Vector{T}(missing, length(v))
    out[k+1:end] = v[1:end-k]
    return out
end

function build_macro_panel(shock_df, oil_df, ffr_df, cpi_df, indpro_df, t10y3m_df)::DataFrame
    println("Building macro panel with lags...")
    macro_df = outerjoin(shock_df, oil_df, on=:date)
    macro_df = outerjoin(macro_df, ffr_df, on=:date)
    macro_df = outerjoin(macro_df, cpi_df, on=:date)
    macro_df = outerjoin(macro_df, indpro_df, on=:date)
    macro_df = outerjoin(macro_df, t10y3m_df, on=:date)
    sort!(macro_df, :date)

    for l in 1:L_LAG
        macro_df[!, Symbol("shock_lag", l)] = lag_vec(macro_df.shock, l)
    end
    macro_df[!, :log_oil_lag1] = lag_vec(macro_df.log_oil_price, 1)
    macro_df[!, :ffr_lag1]     = lag_vec(macro_df.fedfunds, 1)
    macro_df[!, :cpi_lag1]     = lag_vec(macro_df.cpi, 1)
    macro_df[!, :indpro_lag1]  = lag_vec(macro_df.indpro, 1)
    macro_df[!, :t10y3m_lag1]  = lag_vec(macro_df.t10y3m, 1)
    macro_df = dropmissing(macro_df)
    println("  Macro panel rows: ", nrow(macro_df))
    return macro_df
end

# =============================================================================
# 15. MAIN
# =============================================================================
"""
    main(variant::Symbol) -> (panel::DataFrame, oilshare_df::DataFrame, ind_intensity_df::DataFrame)

Build the industry-level panel and output the OilShare vector for cross-sectional analysis.

Returns:
  panel           : (ind_group × date) panel for LP estimation
  oilshare_df     : occupation-group-level Bartik oil exposure for occ-level cross-sectional analysis
  ind_intensity_df: industry-group-level oil intensity for ind-level cross-sectional analysis
"""
function main(variant::Symbol = :hourly_rate)
    println("="^60)
    println("Running variant: $variant  [ind_4]")
    println("="^60)

    download_gdrive_large(CPS_FILE_ID, CPS_PATH; expected_hash=CPS_SHA256)

    cps_raw = parse_cps(CPS_PATH)
    cpi_df  = load_cpi()

    # OilShare construction (full sample, time-invariant)
    oilshare_df, ind_intensity_df = build_oilshare(cps_raw)

    # Save OilShare for later analysis
    output_dir = get_output_dir(variant)
    CSV.write(joinpath(dirname(output_dir), "oilshare_by_occ.csv"),    oilshare_df)
    CSV.write(joinpath(dirname(output_dir), "oil_intensity_by_ind.csv"), ind_intensity_df)
    println("OilShare saved.")

    # Labor panel
    cps_labor = clean_cps_labor(cps_raw)
    panel_labor = build_panel_unemp(cps_labor)

    # Earnings panel
    cps_earn = clean_cps_earn(cps_raw)
    cps_earn = deflate_wages!(cps_earn, cpi_df)
    cps_earn = add_log_rwage!(cps_earn, variant)
    panel_earn = build_panel_earn(cps_earn, variant)

    # Merge
    panel_all = leftjoin(panel_labor, panel_earn, on=[:ind_group, :date, :year, :month])

    # Macro
    shock_df  = load_oil_shock()
    oil_df    = load_oil_price()
    ffr_df    = load_fred("FEDFUNDS.csv", :fedfunds)
    indpro_df = load_fred("INDPRO.csv",   :indpro)
    t10y3m_df = load_fred("T10Y3M.csv",   :t10y3m)
    macro_df  = build_macro_panel(shock_df, oil_df, ffr_df, cpi_df, indpro_df, t10y3m_df)

    panel = leftjoin(panel_all, macro_df, on=:date)
    panel = dropmissing(panel, [:shock, :log_oil_lag1, :ffr_lag1, :cpi_lag1,
                                 :indpro_lag1, :t10y3m_lag1])
    sort!(panel, [:ind_group, :date])

    println("\nVariant '$variant' complete — $(nrow(panel)) obs, $(length(unique(panel.ind_group))) ind groups")
    return panel, oilshare_df, ind_intensity_df
end

 
