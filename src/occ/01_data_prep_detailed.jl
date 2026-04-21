# =============================================================================
# 01_data_prep_occ1990.jl
# Oil Supply News Shock × Occupational Wage Heterogeneity
# Step 1: Load all data, construct occ1990-level panel
#
# Cell = occ1990 code × month
# occ1990 codes: those that previously mapped to groups 1–9
#
# Usage:
#   include("01_data_prep_occ1990.jl")
#   panel = main(:hourly_rate)
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

const DATE_START = Date(1983, 4, 1)
const DATE_END   = Date(2025, 6, 1)

const L_LAG = 12

# All occ1990 codes that previously mapped to groups 1–9
const VALID_OCC1990 = Set(vcat(
    collect(3:37),
    collect(43:200),
    collect(203:235),
    collect(243:283),
    collect(303:389),
    collect(405:469),
    collect(473:498), collect(558:599), collect(614:617),
    collect(503:549), collect(628:699),
    collect(703:799), collect(803:889),
))

function get_output_dir(variant::Symbol)::String
    base = joinpath(@__DIR__, "../..", "result", "occ1990")
    dir  = joinpath(base, string(variant))
    mkpath(dir)
    return dir
end

# =============================================================================
# 1. DOWNLOAD CPS
# =============================================================================
function download_gdrive_large(file_id::String, dest::String;
                                expected_hash::Union{String,Nothing}=nothing)
    if isfile(dest)
        first_bytes = String(read(dest, min(20, filesize(dest))))
        if startswith(first_bytes, "<!") || startswith(first_bytes, "<h")
            println("Existing file is HTML — deleting and retrying.")
            rm(dest)
        else
            println("File already exists: $dest — skipping download.")
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
# 2. DATE PARSING
# =============================================================================
function parse_yearmonth(raw)::Union{Date,Nothing}
    raw isa Date && return raw
    s = strip(string(raw))
    m = match(r"^(\d{4})[Mm](\d{2})$", s)
    isnothing(m) && return nothing
    return Date(parse(Int, m.captures[1]), parse(Int, m.captures[2]), 1)
end

# =============================================================================
# 3. PARSE CPS FIXED-WIDTH
# =============================================================================
function parse_cps(path::String)::DataFrame
    println("Parsing CPS fixed-width file: $path")
    println("File size: ", round(filesize(path)/1e9, digits=2), " GB")

    est_rows = max(1_000_000, round(Int, filesize(path) / 45))

    year_v         = Vector{Int32}(undef, est_rows)
    month_v        = Vector{Int32}(undef, est_rows)
    wtfinl_v       = Vector{Float32}(undef, est_rows)
    earnweek2_v    = Vector{Float32}(undef, est_rows)
    age_v          = Vector{Int32}(undef, est_rows)
    sex_v          = Vector{Int32}(undef, est_rows)
    marst_v        = Vector{Int32}(undef, est_rows)
    empstat_v      = Vector{Int32}(undef, est_rows)
    occ1990_v      = Vector{Int32}(undef, est_rows)
    classwkr_v     = Vector{Int32}(undef, est_rows)
    hours_worked_v = Vector{Float32}(undef, est_rows)
    earnwt_v       = Vector{Float32}(undef, est_rows)

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
                !isempty(hrs_str) && (hrs_raw = tryparse(Float64, hrs_str))
            else
                hrs_str = strip(@view line[43:45])
                !isempty(hrs_str) && (hrs_raw = tryparse(Float64, hrs_str))
            end
            isnothing(hrs_raw) && continue

            occ_raw = tryparse(Int32, strip(@view line[35:37]))
            isnothing(occ_raw) && continue
            !(occ_raw in VALID_OCC1990) && continue

            n += 1
            if n > length(year_v)
                new_cap = round(Int, length(year_v) * 1.5)
                foreach(v -> resize!(v, new_cap),
                    (year_v, month_v, wtfinl_v, earnweek2_v, age_v, sex_v, marst_v,
                     empstat_v, occ1990_v, classwkr_v, hours_worked_v, earnwt_v))
            end

            year_v[n]         = year
            month_v[n]        = month
            wtfinl_v[n]       = Float32(wtfinl_raw / 10_000.0)
            age_v[n]          = parse(Int32, @view line[29:30])
            sex_v[n]          = parse(Int32, @view line[31:31])
            marst_v[n]        = parse(Int32, @view line[32:32])
            empstat_v[n]      = parse(Int32, @view line[33:34])
            occ1990_v[n]      = occ_raw
            classwkr_v[n]     = parse(Int32, @view line[41:42])
            hours_worked_v[n] = Float32(hrs_raw)
            earnwt_v[n]       = Float32(earnwt_raw / 10_000.0)
            earnweek2_v[n]    = Float32(earnwk_raw / 100.0)
        end
    end

    df = DataFrame(
        year         = year_v[1:n],
        month        = month_v[1:n],
        wtfinl       = wtfinl_v[1:n],
        age          = age_v[1:n],
        sex          = sex_v[1:n],
        marst        = marst_v[1:n],
        empstat      = empstat_v[1:n],
        occ1990      = occ1990_v[1:n],
        classwkr     = classwkr_v[1:n],
        hours_worked = hours_worked_v[1:n],
        earnwt       = earnwt_v[1:n],
        earnweek     = earnweek2_v[1:n],
    )
    println("  Raw rows in sample period: ", nrow(df))
    return df
end

# =============================================================================
# 4. SAMPLE SELECTION
# =============================================================================
function clean_cps_labor(df::DataFrame)::DataFrame
    println("Cleaning CPS (labor sample)...")
    df = @subset(df, :empstat .∈ Ref([10, 12, 20, 21, 22]))
    df = @subset(df, 15 .<= :age .<= 64)
    df = @subset(df, :wtfinl .> 0)
    df = @subset(df, :classwkr .∈ Ref([21, 22, 23, 24, 25, 27, 28]))
    df[!, :date] = Date.(df.year, df.month, 1)
    println("  Labor sample: ", nrow(df))
    return df
end

function build_panel_unemp(cps::DataFrame)::DataFrame
    gdf = groupby(cps, [:occ1990, :date, :year, :month])
    panel = combine(gdf) do sdf
        emp        = sdf.empstat
        wt         = sdf.wtfinl
        unemployed = sum((emp .∈ Ref((20,21,22))) .* wt)
        laborforce = sum((emp .∈ Ref((10,12,20,21,22))) .* wt)
        employed   = sum((emp .∈ Ref((10,12))) .* wt)
        unemp_rate    = laborforce == 0 ? missing : unemployed / laborforce
        log_emp_count = employed > 0 ? log(employed) : missing
        pop_weight    = sum(wt)
        (; unemp_rate, log_emp_count, pop_weight)
    end
    panel = @subset(panel, :pop_weight .>= 1000)
    return sort(panel, [:occ1990, :date])
end

function clean_cps_earn(df::DataFrame)::DataFrame
    println("Cleaning CPS (earnings sample)...")
    df = @subset(df, :empstat .∈ Ref([10, 12]))
    df = @subset(df, 15 .<= :age .<= 64)
    df = @subset(df, :earnweek .> 0, :earnweek .!= 999999.99)
    function topcode_limit(year::Integer)
        year <= 1988 ? 999.0 : year <= 1997 ? 1923.0 : 2884.61
    end
    df = @transform(df, :earnweek = min.(:earnweek, topcode_limit.(:year)))
    df = @subset(df, :earnwt .> 0)
    df = @subset(df, 0 .< :hours_worked .<= 105)
    df = @subset(df, :classwkr .∈ Ref([21, 22, 23, 24, 25, 27, 28]))
    df[!, :female]  = Int.(df.sex .== 2)
    df[!, :married] = Int.(df.marst .∈ Ref([1, 2]))
    df[!, :date]    = Date.(df.year, df.month, 1)
    println("  Earnings sample: ", nrow(df))
    return df
end

# =============================================================================
# 5. CPI + DEFLATE
# =============================================================================
function load_cpi()::DataFrame
    println("Loading CPI...")
    path = joinpath(DATA_DIR, "VARdata.xlsx")
    xf   = XLSX.readxlsx(path)
    sh   = xf["Monthly"]
    dates = Vector{Date}(); cpi_vals = Vector{Float64}()
    for row in XLSX.eachrow(sh)
        XLSX.row_number(row) == 1 && continue
        d_raw = row[1]; c_raw = row[7]
        (ismissing(d_raw) || ismissing(c_raw)) && continue
        d = parse_yearmonth(d_raw); isnothing(d) && continue
        push!(dates, d); push!(cpi_vals, Float64(c_raw))
    end
    df = DataFrame(date = dates, cpi = cpi_vals)
    df = @subset(df, DATE_START .<= :date .<= DATE_END)
    sort!(df, :date)
    println("  CPI rows: ", nrow(df))
    return df
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
    variant == :hours && (cps[!, :log_hours] = log.(cps.hours_worked))
    return cps
end

# =============================================================================
# 6. CELL-LEVEL PANEL (occ1990 × date)
# =============================================================================
function weighted_quantile(x::AbstractVector, w::AbstractVector, q::Float64)::Float64
    idx = sortperm(x)
    xs, ws = x[idx], w[idx]
    cumw = cumsum(ws); cumw ./= cumw[end]
    i = searchsortedfirst(cumw, q)
    return xs[clamp(i, 1, length(xs))]
end

function build_panel_earn(cps::DataFrame, variant::Symbol)::DataFrame
    gdf = groupby(cps, [:occ1990, :date, :year, :month])

    panel = if variant == :hourly_rate
        combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt)) => :log_rincome,
            [:log_rwage,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :log_rwage,
            [:age,         :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :age_mean,
            [:female,      :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :female_share,
            [:married,     :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :married_share,
            [:hours_worked,:earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :hours_mean,
            :log_rincome => length => :n_obs,
        )

    elseif variant == :hours
        combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt)) => :log_rincome,
            [:log_rwage,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :log_rwage,
            [:log_hours,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :log_hours,
            [:age,         :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :age_mean,
            [:female,      :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :female_share,
            [:married,     :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :married_share,
            :log_rincome => length => :n_obs,
        )

    elseif variant == :income
        combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt)) => :log_rincome,
            [:log_rwage,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :log_rwage,
            [:age,         :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :age_mean,
            [:female,      :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :female_share,
            [:married,     :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :married_share,
            [:hours_worked,:earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :hours_mean,
            :log_rincome => length => :n_obs,
        )

    elseif variant == :inequality
        p = combine(gdf,
            [:log_rincome, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))         => :log_rincome,
            [:log_rwage,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))         => :log_rwage,
            [:log_rincome, :earnwt] => ((x, wt) -> weighted_quantile(x, wt, 0.25)) => :log_rincome_p25,
            [:log_rincome, :earnwt] => ((x, wt) -> weighted_quantile(x, wt, 0.75)) => :log_rincome_p75,
            [:age,         :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))         => :age_mean,
            [:female,      :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))         => :female_share,
            [:married,     :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))         => :married_share,
            [:hours_worked,:earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))         => :hours_mean,
            :log_rincome => length => :n_obs,
        )
        p = @subset(p, :n_obs .>= 30)
        p[!, :log_ratio_7525] = p.log_rincome_p75 .- p.log_rincome_p25
        p

    elseif variant == :median
        combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> weighted_quantile(w, wt, 0.5)) => :log_rincome_p50,
            [:log_rwage,   :earnwt] => ((x, wt) -> weighted_quantile(x, wt, 0.5)) => :log_rwage_p50,
            [:age,         :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))        => :age_mean,
            [:female,      :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))        => :female_share,
            [:married,     :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))        => :married_share,
            [:hours_worked,:earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt))        => :hours_mean,
            :log_rincome => length => :n_obs,
        )

    elseif variant in (:unemployment, :employment)
        combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt)) => :log_rincome,
            [:log_rwage,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :log_rwage,
            [:age,         :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :age_mean,
            [:female,      :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :female_share,
            [:married,     :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :married_share,
            [:hours_worked,:earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :hours_mean,
            :log_rincome => length => :n_obs,
        )

    elseif variant == :income_share_var
        p = combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt)) => :log_rincome,
            [:log_rwage,   :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :log_rwage,
            [:age,         :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :age_mean,
            [:female,      :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :female_share,
            [:married,     :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :married_share,
            [:hours_worked,:earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :hours_mean,
            [:earnweek,    :earnwt] => ((e, wt) -> sum(e .* wt))            => :total_income_weekly,
            [:earnwt]               => sum                                   => :n_employed,
            :log_rincome => length => :n_obs,
        )
        p = @subset(p, :n_obs .>= 30)
        monthly_agg = combine(groupby(p, :date),
            :total_income_weekly => sum => :total_income_all)
        p = leftjoin(p, monthly_agg, on = :date)
        p[!, :income_share] = p.total_income_weekly ./ p.total_income_all
        p

    else
        error("Unknown variant: $variant")
    end

    if !(variant in (:inequality, :income_share_var))
        panel = @subset(panel, :n_obs .>= 30)
    end

    println("  Earn panel rows: ", nrow(panel),
            " | unique occ1990: ", length(unique(panel.occ1990)))
    return sort(panel, [:occ1990, :date])
end

# =============================================================================
# 7. MACRO DATA
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
    df = DataFrame(date = dates, shock = shocks)
    df = @subset(df, DATE_START .<= :date .<= DATE_END)
    sort!(df, :date)
    println("  Shock rows: ", nrow(df))
    return df
end

function load_oil_price()::DataFrame
    path = joinpath(DATA_DIR, "VARdata.xlsx")
    xf   = XLSX.readxlsx(path)
    sh   = xf["Monthly"]
    dates = Vector{Date}(); prices = Vector{Float64}()
    for row in XLSX.eachrow(sh)
        XLSX.row_number(row) == 1 && continue
        d_raw = row[1]; p_raw = row[2]
        (ismissing(d_raw) || ismissing(p_raw)) && continue
        d = parse_yearmonth(d_raw); isnothing(d) && continue
        push!(dates, d); push!(prices, Float64(p_raw))
    end
    df = DataFrame(date = dates, oil_price = prices)
    df[!, :log_oil_price] = log.(df.oil_price)
    df = @subset(df, DATE_START .<= :date .<= DATE_END)
    sort!(df, :date)
    return df
end

function load_fred(filename::String, colname::Symbol)::DataFrame
    path = joinpath(DATA_DIR, filename)
    df = CSV.read(path, DataFrame;
        types=Dict(1 => Date, 2 => Float32), dateformat="yyyy-m-d")
    rename!(df, first(names(df)) => :date, names(df)[2] => colname)
    df = filter(row -> DATE_START <= row.date <= DATE_END, df)
    sort!(df, :date)
    return df
end

function lag_vec(v::AbstractVector, k::Int)
    T   = Union{eltype(v), Missing}
    out = Vector{T}(missing, length(v))
    out[k+1:end] = v[1:end-k]
    return out
end

function build_macro_panel(shock_df, oil_df, ffr_df, cpi_df, indpro_df, t10y3m_df)::DataFrame
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
    macro_df[!, :indpro_lag1]  = lag_vec(macro_df.indpro, 1)
    macro_df[!, :t10y3m_lag1]  = lag_vec(macro_df.t10y3m, 1)
    macro_df = dropmissing(macro_df)
    println("  Macro panel rows: ", nrow(macro_df))
    return macro_df
end

# =============================================================================
# 8. MAIN
# =============================================================================
"""
    main(variant::Symbol) -> DataFrame

Build occ1990-level panel for LP estimation.
Cell = occ1990 code × month.

variant ∈ :hourly_rate | :hours | :income | :income_share_var |
           :inequality | :median | :unemployment | :employment
"""
function main(variant::Symbol = :hourly_rate)
    println("="^60)
    println("Running variant: $variant  [occ1990-level]")
    println("="^60)

    output_dir = get_output_dir(variant)
    download_gdrive_large(CPS_FILE_ID, CPS_PATH; expected_hash = CPS_SHA256)

    cps_raw = parse_cps(CPS_PATH)
    cpi_df  = load_cpi()

    cps_clean_labor = clean_cps_labor(cps_raw)
    panel_labor     = build_panel_unemp(cps_clean_labor)

    cps_clean_earn  = clean_cps_earn(cps_raw)
    cps_clean_earn  = deflate_wages!(cps_clean_earn, cpi_df)
    cps_clean_earn  = add_log_rwage!(cps_clean_earn, variant)
    panel_earn      = build_panel_earn(cps_clean_earn, variant)

    panel_all = leftjoin(panel_labor, panel_earn, on = [:occ1990, :date, :year, :month])

    shock_df   = load_oil_shock()
    oil_df     = load_oil_price()
    ffr_df     = load_fred("FEDFUNDS.csv", :fedfunds)
    indpro_df  = load_fred("INDPRO.csv", :indpro)
    t10y3m_df  = load_fred("T10Y3M.csv", :t10y3m)
    macro_df   = build_macro_panel(shock_df, oil_df, ffr_df, cpi_df, indpro_df, t10y3m_df)

    panel = leftjoin(panel_all, macro_df, on = :date)
    panel = dropmissing(panel, [:shock, :log_oil_lag1, :ffr_lag1, :cpi_lag1,
                                :indpro_lag1, :t10y3m_lag1])
    sort!(panel, [:occ1990, :date])

    println("\n✓ Variant '$variant' — $(nrow(panel)) obs | ",
            length(unique(panel.occ1990)), " occ1990 codes | dir=$output_dir")
    return panel
end