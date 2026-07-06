# =============================================================================
# data_prep_occ.jl
# Oil Supply News Shock × Occupational Wage Heterogeneity
# Step 1: Load all data, construct cell-level panel

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

const CPS_FILE_ID = "1bz1vYmddvubhLfp_TQ2x83hdUZ3NAZfu"
const CPS_SHA256  = "67CFDB7EE29C81587F8C7F269328FB8FE43E58D9868C5F6AB98C5139575A5D4E"
const CPS_PATH    = joinpath(DATA_DIR, "cps_00018.dat")

# Sample period
const DATE_START = Date(1983, 4, 1)
const DATE_END   = Date(2025, 6, 1)

# LP horizons
const L_LAG = 12   # shock lags

function get_output_dir(variant::Symbol)::String
    base = joinpath(@__DIR__, "../..", "result", "occ")
    dir = if variant == :hourly_rate
        joinpath(base, "hourly_rate")
    elseif variant == :hours
        joinpath(base, "hours")
    elseif variant == :income
        joinpath(base, "income")
    elseif variant == :income_share_var
        joinpath(base, "income_share")   
    elseif variant == :inequality
        joinpath(base, "inequality")
    elseif variant == :median
        joinpath(base, "median")
    elseif variant == :unemployment
        joinpath(base, "unemployment")  
    elseif variant == :employment
        joinpath(base, "employment")
    else
        error("Unknown variant: $variant. Choose from :hourly_rate, :hours, :income, :inequality, :median, :unemployment")
    end
    mkpath(dir)
    return dir
end

# =============================================================================
# 1. DOWNLOAD CPS (if needed)
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
# 2. DATE PARSING UTILITIES
# =============================================================================

"""
    parse_yearmonth(raw) -> Union{Date, Nothing}

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
# 3. OCCUPATION CLASSIFICATION
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
# 4. INDUSTRY CLASSIFICATION (IND1990 -> 13 major groups)  [currently unused]
# =============================================================================
# const IND_LABELS = Dict(
#     1  => "Agriculture_forestry_fishing",  2  => "Mining",
#     3  => "Construction",                  4  => "Manufacturing_nondurable",
#     5  => "Manufacturing_durable",         6  => "Transportation_utilities",
#     7  => "Wholesale_trade",               8  => "Retail_trade",
#     9  => "Finance_insurance_realestate",  10 => "Business_repair_services",
#     11 => "Personal_entertainment_services", 12 => "Professional_related_services",
#     13 => "Public_administration",
# )

# function classify_ind1990(ind::Union{Integer,Missing})::Union{Int,Missing}
#     ismissing(ind) && return missing
#     ind in 10:32   && return 1;  ind in 40:50   && return 2
#     ind == 60      && return 3;  ind in 100:222 && return 4
#     ind in 230:392 && return 5;  ind in 400:472 && return 6
#     ind in 500:571 && return 7;  ind in 580:691 && return 8
#     ind in 700:712 && return 9;  ind in 721:760 && return 10
#     ind in 761:810 && return 11; ind in 812:893 && return 12
#     ind in 900:932 && return 13
#     return missing
# end

# =============================================================================
# 5. PARSE CPS FIXED-WIDTH FILE  
# =============================================================================
# Column layout from data dictionary:
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

function parse_cps(path::String)::DataFrame
    println("Parsing CPS fixed-width file: $path")
    println("File size: ", round(filesize(path)/1e9, digits=2), " GB")

    est_rows = max(1_000_000, round(Int, filesize(path) / 45))

    year_v      = Vector{Int32}(undef, est_rows)
    month_v     = Vector{Int32}(undef, est_rows)
    wtfinl_v    = Vector{Float32}(undef, est_rows)
    earnweek2_v  = Vector{Float32}(undef, est_rows)
    age_v       = Vector{Int32}(undef, est_rows)
    sex_v       = Vector{Int32}(undef, est_rows)
    marst_v     = Vector{Int32}(undef, est_rows)
    empstat_v   = Vector{Int32}(undef, est_rows)
    occ1990_v   = Vector{Int32}(undef, est_rows)
    classwkr_v  = Vector{Int32}(undef, est_rows)
    hours_worked_v = Vector{Float32}(undef, est_rows)
    earnwt_v    = Vector{Float32}(undef, est_rows)

    n = 0

    open(path, "r") do f
        for line in eachline(f)
            length(line) < 40 && continue

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

            hrs_raw = nothing
            if year < 1994
                hrs_str = strip(@view line[46:48])
                if !isempty(hrs_str)
                    hrs_raw = tryparse(Float64, hrs_str)
                end
            else
                hrs_str = strip(@view line[43:45])
                if !isempty(hrs_str)
                    hrs_raw = tryparse(Float64, hrs_str)
                end
            end

            isnothing(hrs_raw) && continue

            n += 1
            if n > length(year_v)
                new_cap = round(Int, length(year_v) * 1.5)
                foreach(v -> resize!(v, new_cap),
                    (year_v, month_v, wtfinl_v, earnweek2_v, age_v, sex_v, marst_v,
                     empstat_v, occ1990_v, classwkr_v, hours_worked_v, earnwt_v))
            end

            year_v[n]      = year
            month_v[n]     = month
            wtfinl_v[n]    = Float32(wtfinl_raw / 10_000.0)
            age_v[n]       = parse(Int32, @view line[29:30])
            sex_v[n]       = parse(Int32, @view line[31:31])
            marst_v[n]     = parse(Int32, @view line[32:32])
            empstat_v[n]   = parse(Int32, @view line[33:34])
            occ1990_v[n]   = parse(Int32, @view line[35:37])
            classwkr_v[n]  = parse(Int32, @view line[41:42])
            hours_worked_v[n] = Float32(hrs_raw)
            earnwt_v[n]    = Float32(earnwt_raw / 10_000.0)
            earnweek2_v[n]  = Float32(earnwk_raw / 100.0)
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
        classwkr  = classwkr_v[1:n],
        hours_worked = hours_worked_v[1:n],
        earnwt    = earnwt_v[1:n],
        earnweek  = earnweek2_v[1:n],
    )

    println("  Raw rows in sample period: ", nrow(df))
    return df
end

# =============================================================================
# 6. SAMPLE SELECTION & VARIABLE CONSTRUCTION
# =============================================================================
function clean_cps_labor(df::DataFrame)::DataFrame
    println("Cleaning CPS (labor sample)...")

    df = @subset(df, :empstat .∈ Ref([10, 12, 20, 21, 22]))
    df = @subset(df, 15 .<= :age .<= 64)
    df = @subset(df, :wtfinl .> 0)
    df = @subset(df, :classwkr .∈ Ref([21, 22, 23, 24, 25, 27, 28]))

    df[!, :occ_group] = map(classify_occ1990, df.occ1990)
    df = @subset(df, .!ismissing.(:occ_group))
    df[!, :occ_group] = convert(Vector{Int}, df.occ_group)

    df[!, :date] = Date.(df.year, df.month, 1)

    println(" Labor sample size:", nrow(df))
    return df
end

function build_panel_unemp(cps::DataFrame)::DataFrame
    gdf = groupby(cps, [:occ_group, :date, :year, :month])

    panel = 
        combine(gdf) do sdf
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
    return sort(panel, [:occ_group, :date])
end

function clean_cps_earn(df::DataFrame)::DataFrame
    println("Cleaning CPS (earnings sample)...")

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

    df[!, :occ_group] = map(classify_occ1990, df.occ1990)
    df = @subset(df, .!ismissing.(:occ_group))
    df[!, :occ_group] = convert(Vector{Int}, df.occ_group)

    df[!, :female]  = Int.(df.sex .== 2)
    df[!, :married] = Int.(df.marst .∈ Ref([1, 2]))

    df[!, :date] = Date.(df.year, df.month, 1)

    println("   Earnings sample size: ", nrow(df))
    return df
end

# =============================================================================
# 7. LOAD CPI AND DEFLATE WAGES
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
        d_raw = row[1]
        c_raw = row[7]   # Column G = CPI
        (ismissing(d_raw) || ismissing(c_raw)) && continue
        d = parse_yearmonth(d_raw)
        isnothing(d) && continue
        push!(dates,    d)
        push!(cpi_vals, Float64(c_raw))
    end

    cpi_df = DataFrame(date = dates, cpi = cpi_vals)
    cpi_df = @subset(cpi_df, DATE_START .<= :date .<= DATE_END)
    sort!(cpi_df, :date)
    println("  CPI rows: ", nrow(cpi_df))
    return cpi_df
end

function deflate_wages!(cps::DataFrame, cpi::DataFrame)::DataFrame
    println("Deflating wages with CPI...")
    cps = leftjoin(cps, select(cpi, :date, :cpi), on = :date)
    cps = @subset(cps, .!ismissing.(:cpi))
    cps[!, :real_earnweek] = cps.earnweek ./ cps.cpi
    cps[!, :log_rincome]   = log.(cps.real_earnweek)
    return cps
end

function add_log_rwage!(cps::DataFrame, variant::Symbol)::DataFrame
    println("Computing log real hourly wage...")
    cps[!, :log_rwage] = cps.log_rincome .- log.(cps.hours_worked)
    # :hours variant also needs log hours
    if variant == :hours
        cps[!, :log_hours] = log.(cps.hours_worked)
    end
    return cps
end

# =============================================================================
# 8. CONSTRUCT CELL-LEVEL PANEL  (variant-specific)
# =============================================================================

# 加权分位数（inequality / median 共用）
function weighted_quantile(x::AbstractVector, w::AbstractVector, q::Float64)::Float64
    idx    = sortperm(x)
    xs, ws = x[idx], w[idx]
    cumw   = cumsum(ws)
    cumw ./= cumw[end]
    i = searchsortedfirst(cumw, q)
    return xs[clamp(i, 1, length(xs))]
end

function build_panel_earn(cps::DataFrame, variant::Symbol)::DataFrame
    gdf = groupby(cps, [:occ_group, :date, :year, :month])

    panel = if variant == :hourly_rate
        combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt))  => :log_rincome,
            [:log_rwage,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :log_rwage,
            [:age,         :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :age_mean,
            [:female,      :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :female_share,
            [:married,     :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :married_share,
            [:hours_worked,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :hours_mean,
            :log_rincome => length => :n_obs,
        )

    elseif variant == :hours
        combine(gdf,
            [:log_rincome,   :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt))  => :log_rincome,
            [:log_rwage,     :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :log_rwage,
            [:log_hours, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :log_hours,
            [:age,           :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :age_mean,
            [:female,        :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :female_share,
            [:married,       :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :married_share,
            :log_rincome => length => :n_obs,
        )

    elseif variant == :income
        panel_tmp = combine(gdf,
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
        
    elseif variant == :income_share_var
        panel_tmp = combine(gdf,
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
        panel_tmp = @subset(panel_tmp, :n_obs .>= 30)
        monthly_agg = combine(groupby(panel_tmp, :date),
            :total_income_weekly => sum => :total_income_all,
            :n_employed          => sum => :total_emp_all)
        panel_tmp = leftjoin(panel_tmp, monthly_agg, on = :date)
        panel_tmp[!, :income_share] = panel_tmp.total_income_weekly ./ panel_tmp.total_income_all
        check = combine(groupby(panel_tmp, :date), :income_share => sum => :sum_share)
        @assert all(isapprox.(check.sum_share, 1.0, atol=1e-6))
        panel_tmp

    elseif variant == :inequality
        panel_tmp = combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt))              => :log_rincome,
            [:log_rwage,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))              => :log_rwage,
            [:log_rincome, :earnwt] => ((x, wt) -> weighted_quantile(x, wt, 0.25))      => :log_rincome_p25,
            [:log_rincome, :earnwt] => ((x, wt) -> weighted_quantile(x, wt, 0.75))      => :log_rincome_p75,
            [:age,         :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))              => :age_mean,
            [:female,      :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))              => :female_share,
            [:married,     :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))              => :married_share,
            [:hours_worked,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))              => :hours_mean,
            :log_rincome => length => :n_obs,
        )
        panel_tmp = @subset(panel_tmp, :n_obs .>= 30)
        panel_tmp[!, :log_ratio_7525] = panel_tmp.log_rincome_p75 .- panel_tmp.log_rincome_p25
        panel_tmp

    elseif variant == :median
        combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> weighted_quantile(w, wt, 0.5))  => :log_rincome_p50,
            [:log_rwage,   :earnwt] => ((x, wt) -> weighted_quantile(x, wt, 0.5))  => :log_rwage_p50,
            [:age,         :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))          => :age_mean,
            [:female,      :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))          => :female_share,
            [:married,     :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))          => :married_share,
            [:hours_worked,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))          => :hours_mean,
            :log_rincome => length => :n_obs,
        )

    elseif variant == :unemployment
        combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt))  => :log_rincome,
            [:log_rwage,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :log_rwage,
            [:age,         :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :age_mean,
            [:female,      :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :female_share,
            [:married,     :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :married_share,
            [:hours_worked,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :hours_mean,
            :log_rincome => length => :n_obs,
        )
    
    elseif variant == :employment
        combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt))  => :log_rincome,
            [:log_rwage,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :log_rwage,
            [:age,         :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :age_mean,
            [:female,      :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :female_share,
            [:married,     :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :married_share,
            [:hours_worked,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))  => :hours_mean,
            :log_rincome => length => :n_obs,
        )

    else
        error("Unknown variant in build_panel_earn: $variant")
    end

    if !(variant in (:income, :inequality))
        panel = @subset(panel, :n_obs .>= 30)
    end

    println("  Total observations: ", nrow(panel))
    println("  Unique cells: ", length(unique(panel.occ_group)))

    return sort(panel, [:occ_group, :date])
end

# =============================================================================
# 9. LOAD MACRO VARIABLES
# =============================================================================
function load_oil_shock()::DataFrame
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

function load_oil_price()::DataFrame
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
    println("  Oil price rows: ", nrow(df))
    return df
end

function load_fred(filename::String, colname::Symbol)::DataFrame
    println("Loading $filename...")
    path = joinpath(DATA_DIR, filename)

    df = CSV.read(path, DataFrame;
        types=Dict(1 => Date, 2 => Float32), dateformat="yyyy-m-d")

    rename!(df, first(names(df)) => :date, names(df)[2] => colname)
    df = filter(row -> DATE_START <= row.date && row.date <= DATE_END, df)
    sort!(df, :date)
    return df
end

# =============================================================================
# 10. BUILD MACRO PANEL (with lags)
# =============================================================================
function lag_vec(v::AbstractVector, k::Int)
    T   = Union{eltype(v), Missing}
    out = Vector{T}(missing, length(v))
    out[k+1:end] = v[1:end-k]
    return out
end

function build_macro_panel(shock_df, oil_df, ffr_df, cpi_df, indpro_df, t10y3m_df)::DataFrame
    println("Building macro panel with lags...")

    macro_df = outerjoin(shock_df, oil_df, on = :date)
    macro_df = outerjoin(macro_df, ffr_df, on = :date)
    macro_df = outerjoin(macro_df, cpi_df, on = :date)
    macro_df = outerjoin(macro_df, indpro_df, on = :date)
    macro_df = outerjoin(macro_df, t10y3m_df, on = :date)
    sort!(macro_df, :date)

    for l in 1:L_LAG
        macro_df[!, Symbol("shock_lag", l)] = lag_vec(macro_df.shock, l)
    end

    macro_df[!, :log_oil_lag1] = lag_vec(macro_df.log_oil_price, 1)
    macro_df[!, :ffr_lag1]     = lag_vec(macro_df.fedfunds, 1)
    macro_df[!, :cpi_lag1]     = lag_vec(macro_df.cpi, 1)
    macro_df[!, :indpro_lag1]     = lag_vec(macro_df.indpro, 1)
    macro_df[!, :t10y3m_lag1]     = lag_vec(macro_df.t10y3m, 1)
    macro_df = dropmissing(macro_df)
    println("  Macro panel rows after lag construction: ", nrow(macro_df))
    return macro_df
end

# =============================================================================
# 11. MAIN 
# =============================================================================
"""
    main(variant::Symbol) -> DataFrame

构建并返回指定变体的面板数据集。

`variant` options:
  :hourly_rate   — weighted mean log real hourly wage
  :hours         — weighted mean log hours worked
  :income        — mean weekly income
  :income_share_var
  :inequality    — within-group 75-25 log income ratio (log ratio p75/p25)
  :median        — weighted median log income and hourly wage
  :unemployment  — unemployment rate
  :employment
"""
function main(variant::Symbol = :hourly_rate)
    println("="^60)
    println("Running variant: $variant")
    println("="^60)

    output_dir = get_output_dir(variant)

    # Download CPS if needed
    download_gdrive_large(CPS_FILE_ID, CPS_PATH; expected_hash = CPS_SHA256)

    # Parse CPS
    cps_raw = parse_cps(CPS_PATH)
    cpi_df  = load_cpi()

    # Labor sample (for unemployment panel)
    cps_clean_labor = clean_cps_labor(cps_raw)

    # Earnings sample
    cps_clean_earn = clean_cps_earn(cps_raw)
    cps_clean_earn = deflate_wages!(cps_clean_earn, cpi_df)
    cps_clean_earn = add_log_rwage!(cps_clean_earn, variant)

    # Build panels
    panel_labor = build_panel_unemp(cps_clean_labor)
    panel_earn  = build_panel_earn(cps_clean_earn, variant)
    panel_all   = leftjoin(panel_labor, panel_earn, on = [:occ_group, :date, :year, :month])

    describe(DataFrame(unemp_rate = panel_all.unemp_rate))

    # Load macro data
    shock_df = load_oil_shock()
    oil_df   = load_oil_price()
    ffr_df   = load_fred("FEDFUNDS.csv", :fedfunds)
    indpro_df   = load_fred("INDPRO.csv", :indpro)
    t10y3m_df   = load_fred("T10Y3M.csv", :t10y3m)
    macro_df = build_macro_panel(shock_df, oil_df, ffr_df, cpi_df, indpro_df, t10y3m_df)

    # Merge panel with macro
    panel = leftjoin(panel_all, macro_df, on = :date)
    panel = dropmissing(panel, [:shock, :log_oil_lag1, :ffr_lag1, :cpi_lag1, :indpro_lag1, :t10y3m_lag1])
    sort!(panel, [:occ_group, :date])

    println("\nVariant '$variant' complete — $(nrow(panel)) observations, output_dir = $output_dir")
    return panel
end