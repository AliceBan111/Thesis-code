# =============================================================================
# 01_data_prep.jl
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
const DATA_DIR   = joinpath(@__DIR__, "../..", "data")
const OUTPUT_DIR = joinpath(@__DIR__, "../..", "result", "occ", "unemployment", "unemp")
# const OUTPUT_DIR = joinpath(@__DIR__, "../..", "result", "occ", "unemployment", "emp")
mkpath(OUTPUT_DIR)

const CPS_FILE_ID  = "1X1x9XCU5LGaxcnKrzhQBEb4O165OW4D1"
const CPS_SHA256   = "34C48660BDCDA4F2B6C22FB4D135D91CDCFBE137FF64D7100903E359C98B40B3"
const CPS_PATH     = joinpath(DATA_DIR, "cps_00017.dat")

# Sample period
const DATE_START = Date(1983, 4, 1)
const DATE_END   = Date(2025, 6, 1)

# LP horizons
const L_LAG = 12   # shock lags

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

解析多种日期格式，统一返回当月第一天的 Date：
  - "1975M04" / "1975m04"  (Excel 里的格式)
  - Julia Date 对象         (XLSX 有时直接返回 Date)
  - 其他格式返回 nothing
"""
function parse_yearmonth(raw)::Union{Date, Nothing}
    raw isa Date && return raw
    s = strip(string(raw))
    # 匹配 "1975M04" 或 "1975m04"
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
# 4. INDUSTRY CLASSIFICATION (IND1990 -> 13 major groups)
# =============================================================================
# const IND_LABELS = Dict(
#     1  => "Agriculture_forestry_fishing",
#     2  => "Mining",
#     3  => "Construction",
#     4  => "Manufacturing_nondurable",
#     5  => "Manufacturing_durable",
#     6  => "Transportation_utilities",
#     7  => "Wholesale_trade",
#     8  => "Retail_trade",
#     9  => "Finance_insurance_realestate",
#     10 => "Business_repair_services",
#     11 => "Personal_entertainment_services",
#     12 => "Professional_related_services",
#     13 => "Public_administration",
# )

# function classify_ind1990(ind::Union{Integer,Missing})::Union{Int,Missing}
#     ismissing(ind) && return missing
#     ind in 10:32   && return 1
#     ind in 40:50   && return 2
#     ind == 60      && return 3
#     ind in 100:222 && return 4
#     ind in 230:392 && return 5
#     ind in 400:472 && return 6
#     ind in 500:571 && return 7
#     ind in 580:691 && return 8
#     ind in 700:712 && return 9
#     ind in 721:760 && return 10
#     ind in 761:810 && return 11
#     ind in 812:893 && return 12
#     ind in 900:932 && return 13
#     return missing
# end

# =============================================================================
# 5. PARSE CPS FIXED-WIDTH FILE  （性能优化 + bug修复版）
# =============================================================================
# Column layout from data dictionary:
#   YEAR       1-4
#   MONTH      5-6
#   WTFINL     7-20  (4 implied decimals)
#   AGE        21-22
#   SEX        23-23
#   MARST      24-24
#   EMPSTAT    25-26
#   OCC1990    27-29
#   IND1990    30-32
#   CLASSWKR   33-34
#   UHRSWORKT  35-37
#   EARNWT     38-47  (4 implied decimals)
#   EARNWEEK   48-55  (2 implied decimals)

function parse_cps(path::String)::DataFrame
    println("Parsing CPS fixed-width file: $path")
    println("File size: ", round(filesize(path)/1e9, digits=2), " GB")

    # ── 预分配列数组，避免每行 push!(NamedTuple) 的大量内存分配 ──
    # 粗估行数（每行约 45 字符），后续按需扩容
    est_rows = max(1_000_000, round(Int, filesize(path) / 45))

    year_v     = Vector{Int32}(undef, est_rows)
    month_v    = Vector{Int32}(undef, est_rows)
    wtfinl_v   = Vector{Float32}(undef, est_rows)
    age_v      = Vector{Int32}(undef, est_rows)
    sex_v      = Vector{Int32}(undef, est_rows)
    marst_v    = Vector{Int32}(undef, est_rows)
    empstat_v  = Vector{Int32}(undef, est_rows)
    occ1990_v  = Vector{Int32}(undef, est_rows)
    classwkr_v  = Vector{Int32}(undef, est_rows)
    # ind1990_v  = Vector{Int32}(undef, est_rows)
    uhrsworkt_v  = Vector{Float32}(undef, est_rows)
    earnwt_v   = Vector{Float32}(undef, est_rows)
    earnweek_v = Vector{Float32}(undef, est_rows)

    n = 0   # 实际写入行数

    open(path, "r") do f
        for line in eachline(f)
            length(line) < 40 && continue

            # ── 用 @view 做零拷贝切片，避免 string(line[...]) 创建新对象 ──
            year  = parse(Int32, @view line[1:4])
            month = parse(Int32, @view line[5:6])

            # 提前过滤：不在样本期的行直接跳过（节省后续所有解析）
            d = Date(year, month, 1)
            (d < DATE_START || d > DATE_END) && continue

            # ── EARNWT&WTFINL：tryparse + missing 跳过 ──
            wtfinl_raw = tryparse(Float64, strip(@view line[7:20]))
            isnothing(wtfinl_raw) && continue

            earnwt_raw = tryparse(Float64, strip(@view line[38:47]))
            isnothing(earnwt_raw) && continue

            # ── EARNWEEK：修复原始 bug（parse 遇空串崩溃 → tryparse + 空串检查）──
            earnwk_str = strip(@view line[48:55])
            isempty(earnwk_str) && continue
            earnwk_raw = tryparse(Float64, earnwk_str)
            isnothing(earnwk_raw) && continue

            uhrsworkt_str = strip(@view line[35:37])
            isempty(uhrsworkt_str) && continue
            uhrsworkt_raw = tryparse(Float64, uhrsworkt_str)
            isnothing(uhrsworkt_raw) && continue

            # ── 动态扩容（1.5 倍）──
            n += 1
            if n > length(year_v)
                new_cap = round(Int, length(year_v) * 1.5)
                foreach(v -> resize!(v, new_cap),
                    (year_v, month_v, wtfinl_v, age_v, sex_v, marst_v,
                     empstat_v, occ1990_v, uhrsworkt_v, earnwt_v, earnweek_v, classwkr_v))
            end

            # ── 写入（剩余字段只在通过过滤后才解析）──
            year_v[n]     = year
            month_v[n]    = month
            wtfinl_v[n]   = Float32(wtfinl_raw / 10_000.0) 
            age_v[n]      = parse(Int32, @view line[21:22])
            sex_v[n]      = parse(Int32, @view line[23:23])
            marst_v[n]    = parse(Int32, @view line[24:24])
            empstat_v[n]  = parse(Int32, @view line[25:26])
            occ1990_v[n]  = parse(Int32, @view line[27:29])
            classwkr_v[n]  = parse(Int32, @view line[33:34])
            # ind1990_v[n]  = parse(Int32, @view line[16:18])
            uhrsworkt_v[n]  = Float32(uhrsworkt_raw)
            earnwt_v[n]   = Float32(earnwt_raw / 10_000.0)   # 4 implied decimals
            earnweek_v[n] = Float32(earnwk_raw / 100.0)      # 2 implied decimals
        end
    end

    # ── 截断到实际行数后构造 DataFrame（一次性分配，零额外拷贝）──
    df = DataFrame(
        year     = year_v[1:n],
        month    = month_v[1:n],
        wtfinl   = wtfinl_v[1:n],
        age      = age_v[1:n],
        sex      = sex_v[1:n],
        marst    = marst_v[1:n],
        empstat  = empstat_v[1:n],
        occ1990  = occ1990_v[1:n],
        classwkr = classwkr_v[1:n],
        # ind1990  = ind1990_v[1:n],
        uhrsworkt = uhrsworkt_v[1:n],
        earnwt   = earnwt_v[1:n],
        earnweek = earnweek_v[1:n],
    )

    println("  Raw rows in sample period: ", nrow(df))
    return df
end

# =============================================================================
# 6. SAMPLE SELECTION & VARIABLE CONSTRUCTION
# =============================================================================
function clean_cps_labor(df::DataFrame)::DataFrame
    println("Cleaning CPS (labor sample)...")

    # 6a. Employment: EMPSTAT 10=at work, 12=has job not at work last week
    df = @subset(df, :empstat .∈ Ref([10, 12, 20, 21, 22]))

    # 6b. Age 16-64
    df = @subset(df, 16 .<= :age .<= 64)

    # 6c. Valid weight
    df = @subset(df, :wtfinl .> 0)
    
    # 6d. Remove self-employed
    df = @subset(df, :classwkr .∈ Ref([21, 22, 23, 24, 25, 27, 28]))

    # 6e. Classify occupation and industry
    df[!, :occ_group] = map(classify_occ1990, df.occ1990)

    # Drop unclassified
    df = @subset(df, .!ismissing.(:occ_group))
    df[!, :occ_group] = convert(Vector{Int}, df.occ_group)
    #df[!, :ind_group] = convert(Vector{Int}, df.ind_group)

    # 6f. Date variable（用整数列拼，不依赖字符串格式）
    df[!, :date] = Date.(df.year, df.month, 1)

    println(" Labor sample size:", nrow(df))
    return df
end

function build_panel_unemp(cps::DataFrame)::DataFrame

    gdf = groupby(cps, [:occ_group, :date, :year, :month])

    panel = combine(gdf) do sdf
        emp = sdf.empstat
        wt  = sdf.wtfinl
        unemployed = sum((emp .∈ Ref((20,21,22))) .* wt)
        laborforce = sum((emp .∈ Ref((10,12,20,21,22))) .* wt)
        employed   = sum((emp .∈ Ref((10,12))) .* wt)
        unemp_rate = laborforce == 0 ? missing : unemployed / laborforce
        log_emp_count = employed > 0 ? log(employed) : missing
        pop_weight = sum(wt)
        (; unemp_rate, log_emp_count, pop_weight)
    end

    panel = @subset(panel, :pop_weight .>= 1000)

    return sort(panel, [:occ_group, :date])
end

function clean_cps_earn(df::DataFrame)::DataFrame
    println("Cleaning CPS (earnings sample)...")

    # 6a. Employment: EMPSTAT 10=at work, 12=has job not at work last week
    df = @subset(df, :empstat .∈ Ref([10, 12]))

    # 6b. Age 16-64
    df = @subset(df, 16 .<= :age .<= 64)

    # 6c. Valid earnings: EARNWEEK > 0, not topcoded
    df = @subset(df, :earnweek .> 0, :earnweek .!= 9999.99)

    function topcode_limit(year::Integer)
        year <= 1988 ? 999.0f0 :
        year <= 1997 ? 1923.0f0 : 2884.0f0
    end

    df = @transform(df, :earnweek = min.(:earnweek, topcode_limit.(:year)))

    # function is_not_topcoded(earnweek::Float32, year::Integer)::Bool
    #     if year <= 1988
    #         return earnweek < 999.0
    #     elseif year <= 1997
    #         return earnweek < 1923.0
    #     else  # year >= 1998
    #         return earnweek < 2884.0
    #     end
    # end
    
    # df = @subset(df, is_not_topcoded.(:earnweek, :year))

    # 6d. Valid weight
    df = @subset(df, :earnwt .> 0)

    # 6e. Valid hours
    df = @subset(df, 0 .< :uhrsworkt .<= 105)

    # 6f. Remove self-employed
    df = @subset(df, :classwkr .∈ Ref([21, 22, 23, 24, 25, 27, 28]))

    # 6g. Classify occupation and industry
    df[!, :occ_group] = map(classify_occ1990, df.occ1990)
    #df[!, :ind_group] = map(classify_ind1990, df.ind1990)

    # Drop unclassified
    df = @subset(df, .!ismissing.(:occ_group))
    df[!, :occ_group] = convert(Vector{Int}, df.occ_group)
    #df[!, :ind_group] = convert(Vector{Int}, df.ind_group)

    # 6h. Demographic dummies
    df[!, :female]  = Int.(df.sex .== 2)
    df[!, :married] = Int.(df.marst .∈ Ref([1, 2]))  

    # 6i. Date variable（用整数列拼，不依赖字符串格式）
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
        XLSX.row_number(row) == 1 && continue   # skip header
        d_raw = row[1]
        c_raw = row[7]   # Column G = CPI
        (ismissing(d_raw) || ismissing(c_raw)) && continue

        # ── 修复："1975M04" 格式解析 ──
        d = parse_yearmonth(d_raw)
        isnothing(d) && continue

        push!(dates,    d)
        push!(cpi_vals, Float64(c_raw))
    end

    cpi_df = DataFrame(date = dates, cpi = cpi_vals)
    cpi_df = @subset(cpi_df, DATE_START .<= :date .<= DATE_END)
    sort!(cpi_df, :date)

    # Normalize CPI to DATE_START = 1
    # base_row = @subset(cpi_df, :date .== DATE_START)
    # isempty(base_row) && error("Base CPI date $DATE_START not found in VARdata.xlsx")
    # base_cpi = first(base_row).cpi
    # cpi_df[!, :cpi_norm] = cpi_df.cpi ./ base_cpi

    println("  CPI rows: ", nrow(cpi_df))
    return cpi_df
end

function deflate_wages!(cps::DataFrame, cpi::DataFrame)::DataFrame
    println("Deflating wages with CPI...")
    # cps = leftjoin(cps, select(cpi, :date, :cpi_norm), on = :date)
    cps = leftjoin(cps, select(cpi, :date, :cpi), on = :date)
    cps = @subset(cps, .!ismissing.(:cpi))
    cps[!, :real_earnweek] = cps.earnweek ./ cps.cpi
    cps[!, :log_rincome]     = log.(cps.real_earnweek)
    return cps
end

function add_log_rwage!(cps::DataFrame)::DataFrame
    println("Computing log real hourly wage...")
    
    # log real hourly wage
    cps[!, :log_rwage] = cps.log_rincome .- log.(cps.uhrsworkt)
    
    return cps
end

# =============================================================================
# 8. CONSTRUCT CELL-LEVEL PANEL
# =============================================================================
function build_panel_earn(cps::DataFrame)::DataFrame
    # println("Building cell-level panel (OccGroup × IndGroup × YearMonth)...")

    # cps[!, :cell_id] = string.(cps.occ_group, "_", cps.ind_group)

    # gdf = groupby(cps, [:cell_id, :occ_group, :ind_group, :date, :year, :month])
    gdf = groupby(cps, [:occ_group, :date, :year, :month])

    panel = combine(gdf,
        [:log_rincome, :earnwt] =>
            ((w, wt) -> sum(w .* wt) / sum(wt)) => :log_rincome,
        [:log_rwage, :earnwt] =>
            ((x, wt) -> sum(x .* wt) / sum(wt)) => :log_rwage,
        [:age, :earnwt] =>
            ((x, wt) -> sum(x .* wt) / sum(wt)) => :age_mean,
        [:female, :earnwt] =>
            ((x, wt) -> sum(x .* wt) / sum(wt)) => :female_share,
        [:married, :earnwt] =>
            ((x, wt) -> sum(x .* wt) / sum(wt)) => :married_share,
        [:uhrsworkt, :earnwt] => 
            ((x, wt) -> sum(x .* wt) / sum(wt)) => :hours_mean,
        :log_rincome => length => :n_obs,
    )

    # Drop cells with fewer than 30 observations in a given month
    panel = @subset(panel, :n_obs .>= 30)

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
        d_raw = row[1]   # Column A = Date
        s_raw = row[3]   # Column C = Oil supply news shock
        (ismissing(d_raw) || ismissing(s_raw)) && continue

        # ── 修复："1975M04" 格式解析 ──
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
        d_raw = row[1]   # Column A
        p_raw = row[2]   # Column B = POIL
        (ismissing(d_raw) || ismissing(p_raw)) && continue

        # ── 修复："1975M04" 格式解析 ──
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
    # println("  Rows: ", nrow(df))
    return df
end

# function main()
#     ffr_df    = load_fred("FEDFUNDS.csv",  :fedfunds)
#     println("First 5 rows:")
#     println(first(ffr_df, 5))
# end


# main()

# =============================================================================
# 10. BUILD MACRO PANEL (with lags)
# =============================================================================

# Helper: lag a vector by k periods
function lag_vec(v::AbstractVector, k::Int)
    T = Union{eltype(v), Missing}
    out = Vector{T}(missing, length(v))
    out[k+1:end] = v[1:end-k]
    return out
end

function build_macro_panel(shock_df, oil_df, ffr_df)::DataFrame
    println("Building macro panel with lags...")

    macro_df = outerjoin(shock_df, oil_df,    on = :date)
    macro_df = outerjoin(macro_df, ffr_df,    on = :date)
    # macro_df = outerjoin(macro_df, unrate_df, on = :date)
    sort!(macro_df, :date)

    # Shock lags l = 1 … L_LAG
    for l in 1:L_LAG
        macro_df[!, Symbol("shock_lag", l)] = lag_vec(macro_df.shock, l)
    end

    # One-period lags of controls
    macro_df[!, :log_oil_lag1] = lag_vec(macro_df.log_oil_price, 1)
    macro_df[!, :ffr_lag1]     = lag_vec(macro_df.fedfunds, 1)
    # macro_df[!, :unrate_lag1]  = lag_vec(macro_df.unrate, 1)

    macro_df = dropmissing(macro_df)
    println("  Macro panel rows after lag construction: ", nrow(macro_df))
    return macro_df
end

# =============================================================================
# 11. MAIN
# =============================================================================
function main()
    # Download CPS if needed
    download_gdrive_large(CPS_FILE_ID, CPS_PATH; expected_hash = CPS_SHA256)

    # Load and clean CPS
    cps_raw   = parse_cps(CPS_PATH)
    cpi_df    = load_cpi()
    cps_clean_labor = clean_cps_labor(cps_raw)
    cps_clean_earn = clean_cps_earn(cps_raw)
    cps_clean_earn = deflate_wages!(cps_clean_earn, cpi_df)
    cps_clean_earn = add_log_rwage!(cps_clean_earn)

    # Build panel
    panel_labor = build_panel_unemp(cps_clean_labor)
    panel_earn = build_panel_earn(cps_clean_earn)
    panel_all = leftjoin(panel_labor, panel_earn, on = [:occ_group, :date, :year, :month])
    describe(DataFrame(unemp_rate = panel_all.unemp_rate))

    # Load macro data
    shock_df  = load_oil_shock()
    oil_df    = load_oil_price()
    ffr_df    = load_fred("FEDFUNDS.csv",  :fedfunds)
    # indpro_df = load_fred("INDPRO.csv",    :indpro)
    # unrate_df = load_fred("UNRATE.csv",    :unrate)
    # macro_df  = build_macro_panel(shock_df, oil_df, ffr_df, unrate_df)
    macro_df  = build_macro_panel(shock_df, oil_df, ffr_df)

    # Merge panel with macro
    panel = leftjoin(panel_all, macro_df, on = :date)
    panel = dropmissing(panel, [:shock, :log_oil_lag1, :ffr_lag1])
    sort!(panel, [:occ_group, :date])

    return panel
end



# function check_topcode_proportion(cps_raw::DataFrame)
#     println("\n" * "="^60)
#     println("Topcode Proportion by Occupation")
#     println("="^60)
    
#     # 应用 earnings sample 基础筛选
#     df = @subset(cps_raw, 
#         :empstat .∈ Ref([10, 12]),  # 正在工作
#         16 .<= :age .<= 64,
#         :earnweek .> 0,
#         :earnweek .!= 9999.99,
#         :earnwt .> 0,
#         0 .< :uhrsworkt .<= 105,
#         :classwkr .∈ Ref([21, 22, 23, 24, 25, 27, 28])
#     )
    
#     # 分类职业
#     df[!, :occ_group] = map(classify_occ1990, df.occ1990)
#     df = @subset(df, .!ismissing.(:occ_group))
#     df[!, :occ_group] = convert(Vector{Int}, df.occ_group)
    
#     # 判断是否为 topcode
#     function is_topcoded(earnweek::Float32, year::Integer)::Bool
#         if year <= 1988
#             return earnweek >= 999.0
#         elseif year <= 1997
#             return earnweek >= 1923.0
#         else
#             return earnweek >= 2884.0
#         end
#     end
    
#     df[!, :is_topcoded] = is_topcoded.(df.earnweek, df.year)
    
#     # 按职业统计
#     for occ in sort(unique(df.occ_group))
#         sub = @subset(df, :occ_group .== occ)
#         n_total = nrow(sub)
#         n_top = sum(sub.is_topcoded)
#         pct = n_top / n_total * 100
#         println("$(rpad(OCC_LABELS[occ], 25)): $n_top / $n_total = $(round(pct, digits=2))%")
#     end
    
#     # 总计
#     total_top = sum(df.is_topcoded)
#     total_n = nrow(df)
#     println("-"^60)
#     println("$(rpad("ALL OCCUPATIONS", 25)): $total_top / $total_n = $(round(total_top/total_n*100, digits=2))%")
# end

# # 运行检查
# cps_raw = parse_cps(CPS_PATH)
# check_topcode_proportion(cps_raw)

# =============================================================================
# CHECK TOPCODE YEAR DISTRIBUTION FOR FARMING/FORESTRY/CONSTRUCTION
# =============================================================================

# function check_farming_topcode_by_year(cps_raw::DataFrame)
#     println("\n" * "="^80)
#     println("Topcode Year Distribution: Farming/Forestry/Construction")
#     println("="^80)
    
#     # 定义 Farming/forestry/construction 的 occ1990 范围
#     farming_occ_codes = vcat(473:498, 558:599, 614:617)
    
#     # 筛选该职业的样本
#     df = @subset(cps_raw, 
#         :occ1990 .∈ Ref(farming_occ_codes),
#         :empstat .∈ Ref([10, 12]),  # 正在工作
#         16 .<= :age .<= 64,
#         :earnweek .> 0,
#         :earnweek .!= 9999.99,
#         :earnwt .> 0,
#         0 .< :uhrsworkt .<= 105
#     )
    
#     println("Total observations in Farming/forestry/construction: ", nrow(df))
#     println()
    
#     # 判断是否为 topcode
#     function is_topcoded(earnweek::Float32, year::Integer)::Bool
#         if year <= 1988
#             return earnweek >= 999.0
#         elseif year <= 1997
#             return earnweek >= 1923.0
#         else  # year >= 1998
#             return earnweek >= 2884.0
#         end
#     end
    
#     df[!, :is_topcoded] = is_topcoded.(df.earnweek, df.year)
    
#     # 按年份统计
#     results = DataFrame(
#         year = Int[],
#         n_total = Int[],
#         n_topcoded = Int[],
#         pct_topcoded = Float64[],
#         threshold = Float64[]
#     )
    
#     for year in sort(unique(df.year))
#         sub = @subset(df, :year .== year)
#         n_total = nrow(sub)
#         n_top = sum(sub.is_topcoded)
#         pct = n_top / n_total * 100
        
#         # 确定该年份的阈值
#         threshold = if year <= 1988
#             999.0
#         elseif year <= 1997
#             1923.0
#         else
#             2884.0
#         end
        
#         push!(results, (year, n_total, n_top, pct, threshold))
#     end
    
#     # 打印结果
#     println("Year | Total N | Topcoded N | Percentage | Threshold")

#     for row in eachrow(results)
#         println(@sprintf("%4d | %7d | %9d | %8.2f%% | %8.1f", 
#             row.year, row.n_total, row.n_topcoded, row.pct_topcoded, row.threshold))
#     end
    
#     # 总计
#     total_top = sum(results.n_topcoded)
#     total_n = sum(results.n_total)

#     println(@sprintf("ALL  | %7d | %9d | %8.2f%% |", 
#         total_n, total_top, total_top/total_n*100))
    
#     # 额外：按时期汇总
#     println("\n" * "="^80)
#     println("Summary by Topcode Regime:")
    
#     period1 = @subset(results, :year .<= 1988)
#     period2 = @subset(results, 1989 .<= :year .<= 1997)
#     period3 = @subset(results, :year .>= 1998)
    
#     if !isempty(period1)
#         println("1983-1988 (threshold = 999):   ", 
#             sum(period1.n_topcoded), " / ", sum(period1.n_total), 
#             " = ", @sprintf("%.2f%%", sum(period1.n_topcoded)/sum(period1.n_total)*100))
#     end
    
#     if !isempty(period2)
#         println("1989-1997 (threshold = 1923):  ", 
#             sum(period2.n_topcoded), " / ", sum(period2.n_total), 
#             " = ", @sprintf("%.2f%%", sum(period2.n_topcoded)/sum(period2.n_total)*100))
#     end
    
#     if !isempty(period3)
#         println("1998-2025 (threshold = 2884):  ", 
#             sum(period3.n_topcoded), " / ", sum(period3.n_total), 
#             " = ", @sprintf("%.2f%%", sum(period3.n_topcoded)/sum(period3.n_total)*100))
#     end
    
#     return results
# end

# # 运行检查（在 main() 里或单独运行）
# results = check_farming_topcode_by_year(cps_raw)

# # =============================================================================
# # CHECK INCOME DISTRIBUTION FOR FARMING/FORESTRY/CONSTRUCTION
# # =============================================================================

# function check_farming_income_distribution(cps_raw::DataFrame)
#     println("\n" * "="^80)
#     println("Income Distribution: Farming/Forestry/Construction")
#     println("="^80)
    
#     # 定义 Farming/forestry/construction 的 occ1990 范围
#     farming_occ_codes = vcat(473:498, 558:599, 614:617)
    
#     # 筛选该职业的样本
#     df = @subset(cps_raw, 
#         :occ1990 .∈ Ref(farming_occ_codes),
#         :empstat .∈ Ref([10, 12]),  # 正在工作
#         16 .<= :age .<= 64,
#         :earnweek .> 0,
#         :earnweek .!= 9999.99,
#         :earnwt .> 0,
#         0 .< :uhrsworkt .<= 105
#     )
    
#     println("Total observations: ", nrow(df))
#     println()
    
#     # 判断是否为 topcode
#     function is_topcoded(earnweek::Float32, year::Integer)::Bool
#         if year <= 1988
#             return earnweek >= 999.0
#         elseif year <= 1997
#             return earnweek >= 1923.0
#         else
#             return earnweek >= 2884.0
#         end
#     end
    
#     df[!, :is_topcoded] = is_topcoded.(df.earnweek, df.year)
    
#     # 按时期分组统计
#     df[!, :period] = ifelse.(df.year .<= 1988, "1983-1988",
#                        ifelse.(df.year .<= 1997, "1989-1997", "1998-2025"))
    
#     # 分位数函数（带权重）
#     function weighted_quantile(values::Vector{Float32}, weights::Vector{Float32}, q::Float64)::Float64
#         # 按值排序
#         idx = sortperm(values)
#         sorted_vals = values[idx]
#         sorted_weights = weights[idx]
        
#         # 计算累积权重
#         cum_weight = cumsum(sorted_weights)
#         total_weight = sum(sorted_weights)
#         target = total_weight * q
        
#         # 找到分位数位置
#         pos = findfirst(cum_weight .>= target)
#         if pos === nothing
#             return NaN
#         else
#             return Float64(sorted_vals[pos])
#         end
#     end
    
#     # 按时期输出统计
#     results = DataFrame()
    
#     for period in ["1983-1988", "1989-1997", "1998-2025", "ALL"]
#         if period == "ALL"
#             sub = df
#         else
#             sub = @subset(df, :period .== period)
#         end
        
#         if nrow(sub) == 0
#             continue
#         end
        
#         # 提取数据
#         earn = sub.earnweek
#         wt = sub.earnwt
        
#         # 基础统计
#         n_obs = nrow(sub)
#         n_top = sum(sub.is_topcoded)
#         pct_top = n_top / n_obs * 100
        
#         # 加权分位数
#         p50 = weighted_quantile(earn, wt, 0.50)
#         p75 = weighted_quantile(earn, wt, 0.75)
#         p90 = weighted_quantile(earn, wt, 0.90)
#         p95 = weighted_quantile(earn, wt, 0.95)
#         p99 = weighted_quantile(earn, wt, 0.99)
        
#         # 加权均值
#         wmean = sum(earn .* wt) / sum(wt)
        
#         # 存储结果
#         push!(results, (
#             period = period,
#             n_obs = n_obs,
#             n_topcoded = n_top,
#             pct_topcoded = pct_top,
#             mean = wmean,
#             p50 = p50,
#             p75 = p75,
#             p90 = p90,
#             p95 = p95,
#             p99 = p99
#         ))
#     end
    
#     # 打印结果表格
#     println("\n" * "-"^100)
#     println(@sprintf("%-12s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s",
#         "Period", "N", "Top N", "Top%", "Mean", "p50", "p75", "p90", "p95", "p99"))
#     println("-"^100)
    
#     for row in eachrow(results)
#         println(@sprintf("%-12s | %8d | %8d | %7.2f%% | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f | %8.2f",
#             row.period,
#             row.n_obs,
#             row.n_topcoded,
#             row.pct_topcoded,
#             row.mean,
#             row.p50,
#             row.p75,
#             row.p90,
#             row.p95,
#             row.p99))
#     end
    
#     println("-"^100)
    
#     # 额外：画出收入分布（如果安装了 Plots.jl）
#     # 可选可视化
#     println("\n" * "="^80)
#     println("Income Distribution Summary (Weighted):")
#     println("="^80)
    
#     # 检查是否有 topcode 影响
#     all_period = results[results.period .== "ALL", :]
#     if !isempty(all_period)
#         p99_val = all_period.p99[1]
#         println("\nNote: p99 = \$ $(round(p99_val, digits=2))")
        
#         # 判断 topcode 阈值是否影响
#         if p99_val >= 2884.0
#             println("⚠️  p99 exceeds 1998+ topcode threshold (\$2884) → significant topcoding issue")
#         elseif p99_val >= 1923.0
#             println("⚠️  p99 exceeds 1989-1997 topcode threshold (\$1923) → moderate topcoding issue")
#         elseif p99_val >= 999.0
#             println("⚠️  p99 exceeds 1983-1988 topcode threshold (\$999) → mild topcoding issue")
#         else
#             println("✓  p99 below all topcode thresholds → topcoding not a major concern")
#         end
#     end
    
#     return results
# end

# # 运行检查
# results = check_farming_income_distribution(cps_raw)


# # =============================================================================
# # CHECK WHICH SPECIFIC OCCUPATIONS WITHIN FARMING HAVE TOPCODED OBSERVATIONS
# # =============================================================================

# function check_farming_topcode_by_detailed_occ(cps_raw::DataFrame)
#     println("\n" * "="^80)
#     println("Topcoded Observations by Detailed Occupation (occ1990)")
#     println("Farming/Forestry/Construction Group")
#     println("="^80)
    
#     # 定义 Farming/forestry/construction 的 occ1990 范围
#     farming_occ_codes = vcat(473:498, 558:599, 614:617)
    
#     # 筛选该职业组的样本
#     df = @subset(cps_raw, 
#         :occ1990 .∈ Ref(farming_occ_codes),
#         :empstat .∈ Ref([10, 12]),  # 正在工作
#         16 .<= :age .<= 64,
#         :earnweek .> 0,
#         :earnweek .!= 9999.99,
#         :earnwt .> 0,
#         0 .< :uhrsworkt .<= 105
#     )
    
#     println("Total observations in group: ", nrow(df))
#     println()
    
#     # 判断是否为 topcode
#     function is_topcoded(earnweek::Float32, year::Integer)::Bool
#         if year <= 1988
#             return earnweek >= 999.0
#         elseif year <= 1997
#             return earnweek >= 1923.0
#         else
#             return earnweek >= 2884.0
#         end
#     end
    
#     df[!, :is_topcoded] = is_topcoded.(df.earnweek, df.year)
    
#     # 按 occ1990 详细分类统计
#     results = DataFrame(
#         occ1990 = Int[],
#         n_total = Int[],
#         n_topcoded = Int[],
#         pct_topcoded = Float64[],
#         mean_earn = Float64[],
#         median_earn = Float64[],
#         max_earn = Float64[]
#     )
    
#     for occ in sort(unique(df.occ1990))
#         sub = @subset(df, :occ1990 .== occ)
#         n_total = nrow(sub)
#         n_top = sum(sub.is_topcoded)
#         pct = n_top / n_total * 100
        
#         # 只有存在 topcode 或者样本量足够大的才显示详细信息
#         if n_top > 0 || n_total >= 100
#             mean_earn = mean(sub.earnweek)
#             median_earn = median(sub.earnweek)
#             max_earn = maximum(sub.earnweek)
            
#             push!(results, (occ, n_total, n_top, pct, mean_earn, median_earn, max_earn))
#         end
#     end
    
#     # 按 topcode 比例排序
#     results = sort(results, :pct_topcoded, rev=true)
    
#     # 打印结果
#     println("\nDetailed Occupations with Topcoding (sorted by topcode %):")
#     println("-"^100)
#     println(@sprintf("%-8s | %10s | %10s | %10s | %10s | %10s | %10s", 
#         "occ1990", "Total N", "Top N", "Top %", "Mean", "Median", "Max"))
#     println("-"^100)
    
#     for row in eachrow(results)
#         if row.n_topcoded > 0
#             println(@sprintf("%-8d | %10d | %10d | %9.2f%% | %10.2f | %10.2f | %10.2f",
#                 row.occ1990, row.n_total, row.n_topcoded, row.pct_topcoded, 
#                 row.mean_earn, row.median_earn, row.max_earn))
#         else
#             # 没有 topcode 的也显示，但用不同格式
#             println(@sprintf("%-8d | %10d | %10d | %9.2f%% | %10.2f | %10.2f | %10.2f",
#                 row.occ1990, row.n_total, row.n_topcoded, row.pct_topcoded, 
#                 row.mean_earn, row.median_earn, row.max_earn))
#         end
#     end
    
#     println("-"^100)
    
#     # 汇总统计
#     total_top = sum(results.n_topcoded)
#     total_n = sum(results.n_total)
#     n_occ_with_top = count(results.n_topcoded .> 0)
#     n_occ_total = nrow(results)
    
#     println("\nSummary:")
#     println("  Total occupations in group: $n_occ_total")
#     println("  Occupations with topcoding: $n_occ_with_top")
#     println("  Total topcoded observations: $total_top / $total_n = $(round(total_top/total_n*100, digits=2))%")
    
#     # 列出 top 5 高收入细分职业（按 mean）
#     println("\n" * "="^80)
#     println("Top 5 Highest-Paid Detailed Occupations (by mean earnings):")
#     println("-"^60)
#     top5_mean = sort(results, :mean_earn, rev=true)[1:min(5, nrow(results)), :]
#     for row in eachrow(top5_mean)
#         println(@sprintf("  occ1990=%4d: mean=\$%8.2f, median=\$%8.2f, topcode=%5.2f%%, n=%7d",
#             row.occ1990, row.mean_earn, row.median_earn, row.pct_topcoded, row.n_total))
#     end
    
#     # 列出 topcode 比例最高的 5 个职业
#     println("\n" * "="^80)
#     println("Top 5 Occupations with Highest Topcoding Proportion:")
#     println("-"^60)
#     top5_top = results[results.n_topcoded .> 0, :][1:min(5, sum(results.n_topcoded .> 0)), :]
#     for row in eachrow(top5_top)
#         println(@sprintf("  occ1990=%4d: topcode=%5.2f%% (%d/%d), mean=\$%8.2f",
#             row.occ1990, row.pct_topcoded, row.n_topcoded, row.n_total, row.mean_earn))
#     end
    
#     # 额外：如果有 occ1990 标签字典，可以添加职业名称
#     println("\n" * "="^80)
#     println("Note: occ1990 codes in this group represent:")
#     println("  - 473-498:  Farming occupations")
#     println("  - 558-599:  Forestry and fishing occupations") 
#     println("  - 614-617:  Construction trades")
    
#     return results
# end

# # 运行检查
# detailed_results = check_farming_topcode_by_detailed_occ(cps_raw)


# # =============================================================================
# # CHECK WHICH SPECIFIC OCCUPATIONS WITHIN FARMING HAVE TOPCODED OBSERVATIONS
# # =============================================================================

# function check_farming_topcode_by_detailed_occ(cps_raw::DataFrame)
#     println("\n" * "="^80)
#     println("Topcoded Observations by Detailed Occupation (occ1990)")
#     println("Farming/Forestry/Construction Group")
#     println("="^80)
    
#     # 定义 Farming/forestry/construction 的 occ1990 范围
#     farming_occ_codes = vcat(473:498, 558:599, 614:617)
    
#     # 筛选该职业组的样本
#     df = @subset(cps_raw, 
#         :occ1990 .∈ Ref(farming_occ_codes),
#         :empstat .∈ Ref([10, 12]),  # 正在工作
#         16 .<= :age .<= 64,
#         :earnweek .> 0,
#         :earnweek .!= 9999.99,
#         :earnwt .> 0,
#         0 .< :uhrsworkt .<= 105
#     )
    
#     println("Total observations in group: ", nrow(df))
#     println()
    
#     # 判断是否为 topcode
#     function is_topcoded(earnweek::Float32, year::Integer)::Bool
#         if year <= 1988
#             return earnweek >= 999.0
#         elseif year <= 1997
#             return earnweek >= 1923.0
#         else
#             return earnweek >= 2884.0
#         end
#     end
    
#     df[!, :is_topcoded] = is_topcoded.(df.earnweek, df.year)
    
#     # 按 occ1990 详细分类统计
#     results = DataFrame(
#         occ1990 = Int[],
#         n_total = Int[],
#         n_topcoded = Int[],
#         pct_topcoded = Float64[],
#         mean_earn = Float64[],
#         median_earn = Float64[],
#         max_earn = Float64[]
#     )
    
#     for occ in sort(unique(df.occ1990))
#         sub = @subset(df, :occ1990 .== occ)
#         n_total = nrow(sub)
#         n_top = sum(sub.is_topcoded)
#         pct = n_top / n_total * 100
#         mean_earn = mean(sub.earnweek)
#         median_earn = median(sub.earnweek)
#         max_earn = maximum(sub.earnweek)
        
#         push!(results, (occ, n_total, n_top, pct, mean_earn, median_earn, max_earn))
#     end
    
#     # 按 pct_topcoded 从大到小排序
#     results = sort(results, :pct_topcoded, rev=true)
    
#     # 保存 CSV
#     csv_path = joinpath(OUTPUT_DIR, "farming_topcode_by_occ.csv")
#     CSV.write(csv_path, results)
#     println("\nCSV saved to: $csv_path")
    
#     # 打印前20行到控制台
#     println("\nFirst 20 rows (sorted by topcode % descending):")
#     println("-"^90)
#     println(@sprintf("%-8s | %10s | %10s | %10s | %10s | %10s | %10s", 
#         "occ1990", "Total N", "Top N", "Top %", "Mean", "Median", "Max"))
#     println("-"^90)
    
#     for row in eachrow(results[1:min(20, nrow(results)), :])
#         println(@sprintf("%-8d | %10d | %10d | %9.2f%% | %10.2f | %10.2f | %10.2f",
#             row.occ1990, row.n_total, row.n_topcoded, row.pct_topcoded, 
#             row.mean_earn, row.median_earn, row.max_earn))
#     end
    
#     if nrow(results) > 20
#         println("... and $(nrow(results)-20) more rows (see CSV for full results)")
#     end
    
#     println("-"^90)
    
#     # 汇总统计
#     total_top = sum(results.n_topcoded)
#     total_n = sum(results.n_total)
#     n_occ_with_top = count(results.n_topcoded .> 0)
#     n_occ_total = nrow(results)
    
#     println("\nSummary:")
#     println("  Total occupations in group: $n_occ_total")
#     println("  Occupations with topcoding: $n_occ_with_top")
#     println("  Total topcoded observations: $total_top / $total_n = $(round(total_top/total_n*100, digits=2))%")
    
#     return results
# end

# # 运行检查
# detailed_results = check_farming_topcode_by_detailed_occ(cps_raw)