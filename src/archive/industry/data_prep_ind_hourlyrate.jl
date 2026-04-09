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
const OUTPUT_DIR = joinpath(@__DIR__, "../..", "result", "ind", "hourly_rate")
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
    ind1990_v  = Vector{Int32}(undef, est_rows)
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
                     empstat_v, occ1990_v, uhrsworkt_v, earnwt_v, earnweek_v, ind1990_v, classwkr_v))
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
            ind1990_v[n]  = parse(Int32, @view line[30:32])
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
        ind1990  = ind1990_v[1:n],
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
    df[!, :ind_group] = map(classify_ind1990, df.ind1990)

    # Drop unclassified
    df = @subset(df, .!ismissing.(:ind_group))
    df[!, :ind_group] = convert(Vector{Int}, df.ind_group)
    #df[!, :ind_group] = convert(Vector{Int}, df.ind_group)

    # 6f. Date variable（用整数列拼，不依赖字符串格式）
    df[!, :date] = Date.(df.year, df.month, 1)

    println(" Labor sample size:", nrow(df))
    return df
end

function build_panel_unemp(cps::DataFrame)::DataFrame

    gdf = groupby(cps, [:ind_group, :date, :year, :month])

    panel = combine(gdf) do sdf
        emp = sdf.empstat
        wt  = sdf.wtfinl
        unemployed = sum((emp .∈ Ref((20,21,22))) .* wt)
        laborforce = sum((emp .∈ Ref((10,12,20,21,22))) .* wt)
        unemp_rate = laborforce == 0 ? missing : unemployed / laborforce
        pop_weight = sum(wt)
        (; unemp_rate, pop_weight)
    end

    panel = @subset(panel, :pop_weight .>= 1000)

    return sort(panel, [:ind_group, :date])
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
    # df[!, :occ_group] = map(classify_occ1990, df.occ1990)
    df[!, :ind_group] = map(classify_ind1990, df.ind1990)

    # Drop unclassified
    df = @subset(df, .!ismissing.(:ind_group))
    # df[!, :occ_group] = convert(Vector{Int}, df.occ_group)
    df[!, :ind_group] = convert(Vector{Int}, df.ind_group)

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
    gdf = groupby(cps, [:ind_group, :date, :year, :month])

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
    println("  Unique cells: ", length(unique(panel.ind_group)))

    return sort(panel, [:ind_group, :date])
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
    # macro_df = outerjoin(macro_df, indpro_df, on = :date)
    # macro_df = outerjoin(macro_df, unrate_df, on = :date)
    sort!(macro_df, :date)

    # Shock lags l = 1 … L_LAG
    for l in 1:L_LAG
        macro_df[!, Symbol("shock_lag", l)] = lag_vec(macro_df.shock, l)
    end

    # One-period lags of controls
    macro_df[!, :log_oil_lag1] = lag_vec(macro_df.log_oil_price, 1)
    macro_df[!, :ffr_lag1]     = lag_vec(macro_df.fedfunds, 1)
    # macro_df[!, :indpro_lag1]     = lag_vec(macro_df.indpro, 1)
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
    panel_all = leftjoin(panel_labor, panel_earn, on = [:ind_group, :date, :year, :month])
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
    sort!(panel, [:ind_group, :date])

    return panel
end
