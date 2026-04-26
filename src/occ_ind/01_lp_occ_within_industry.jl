using HTTP
using Downloads
using SHA
using CSV
using DataFrames
using DataFramesMeta
using XLSX
using Dates
using Statistics
using Printf
using LinearAlgebra
using StatsBase
using Random
using FixedEffectModels
using CategoricalArrays
using CairoMakie

const DATA_DIR = joinpath(@__DIR__, "..", "..", "data")
const RESULT_DIR = normpath(joinpath(@__DIR__, "..", "..", "result", "occ_ind"))

const CPS_FILE_ID = "1bz1vYmddvubhLfp_TQ2x83hdUZ3NAZfu"
const CPS_SHA256 = "67CFDB7EE29C81587F8C7F269328FB8FE43E58D9868C5F6AB98C5139575A5D4E"
const CPS_PATH = joinpath(DATA_DIR, "cps_00018.dat")

const DATE_START = Date(1983, 4, 1)
const DATE_END = Date(2025, 6, 1)

const N_BOOT = 500
const BLOCK_SIZE = 6
const BOOT_SEED = 42
const H_MAX = 36
const L_LAG = 12
const L_LAG_MAX = 24
const DK_BW = 12

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

const OCC_COLORS = [
    :steelblue,
    :tomato,
    :seagreen,
    :darkorange,
    :mediumpurple,
    :saddlebrown,
    :hotpink,
    :teal,
    :goldenrod,
]

function variant_folder_name(variant::Symbol)::String
    return string(variant)
end

function get_variant_output_dir(variant::Symbol)::String
    dir = joinpath(RESULT_DIR, variant_folder_name(variant))
    mkpath(dir)
    return dir
end

function get_industry_output_dir(variant::Symbol, ind_group::Int)::String
    variant_dir = get_variant_output_dir(variant)
    ind_label = get(IND_LABELS, ind_group, "industry_$(ind_group)")
    dir = joinpath(variant_dir, ind_label)
    mkpath(dir)
    return dir
end

function get_plot_labels(variant::Symbol)
    labels = Dict(
        :hourly_rate      => "Response of Log Hourly Wage",
        :unemployment     => "Response of Unemployment Rate",
        :hours            => "Response of Log Hours Worked",
        :income           => "Response of Log Income",
        :employment       => "Response of Log Employment",
        :inequality       => "Response of Inequality",
        :median           => "Response of Median Log Income",
        :income_share_var => "Response of Income Share",
    )
    return get(labels, variant, "Response of $(string(variant))")
end

function parse_yearmonth(raw)::Union{Date, Nothing}
    raw isa Date && return raw
    s = strip(string(raw))
    m = match(r"^(\d{4})[Mm](\d{2})$", s)
    isnothing(m) && return nothing
    yr = parse(Int, m.captures[1])
    mo = parse(Int, m.captures[2])
    return Date(yr, mo, 1)
end

function classify_ind1990(ind::Union{Integer, Missing})::Union{Int, Missing}
    ismissing(ind) && return missing
    ind in 10:32 && return 1
    (ind in 40:50 || ind in 200:201) && return 2
    ind == 60 && return 3
    (ind in 100:199 || ind in 202:222) && return 4
    ind in 230:392 && return 5
    ind in 400:472 && return 6
    ind in 500:571 && return 7
    ind in 580:691 && return 8
    ind in 700:712 && return 9
    ind in 721:760 && return 10
    ind in 761:810 && return 11
    ind in 812:893 && return 12
    return missing
end

function classify_occ1990(occ::Union{Integer, Missing})::Union{Int, Missing}
    ismissing(occ) && return missing
    occ in 3:37 && return 1
    occ in 43:200 && return 2
    occ in 203:235 && return 3
    occ in 243:283 && return 4
    occ in 303:389 && return 5
    occ in 405:469 && return 6
    (occ in 473:498 || occ in 558:599 || occ in 614:617) && return 7
    (occ in 503:549 || occ in 628:699) && return 8
    (occ in 703:799 || occ in 803:889) && return 9
    return missing
end

function download_gdrive_large(file_id::String, dest::String;
    expected_hash::Union{String, Nothing}=nothing)

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
    resp = HTTP.get(base_url)
    body = String(resp.body)
    uuid_m = match(r"name=\"uuid\" value=\"([^\"]+)\"", body)
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
        println("Hash verified.")
    end
end

function parse_cps(path::String)::DataFrame
    println("Parsing CPS fixed-width file: $path")
    println("File size: ", round(filesize(path) / 1e9, digits=2), " GB")

    est_rows = max(1_000_000, round(Int, filesize(path) / 48))

    year_v = Vector{Int32}(undef, est_rows)
    month_v = Vector{Int32}(undef, est_rows)
    wtfinl_v = Vector{Float32}(undef, est_rows)
    earnweek_v = Vector{Float32}(undef, est_rows)
    age_v = Vector{Int32}(undef, est_rows)
    sex_v = Vector{Int32}(undef, est_rows)
    marst_v = Vector{Int32}(undef, est_rows)
    empstat_v = Vector{Int32}(undef, est_rows)
    occ1990_v = Vector{Int32}(undef, est_rows)
    ind1990_v = Vector{Int32}(undef, est_rows)
    classwkr_v = Vector{Int32}(undef, est_rows)
    hours_worked_v = Vector{Float32}(undef, est_rows)
    earnwt_v = Vector{Float32}(undef, est_rows)

    n = 0

    open(path, "r") do f
        for line in eachline(f)
            length(line) < 58 && continue

            year = parse(Int32, @view line[1:4])
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
                foreach(v -> resize!(v, new_cap), (
                    year_v, month_v, wtfinl_v, earnweek_v, age_v, sex_v,
                    marst_v, empstat_v, occ1990_v, ind1990_v, classwkr_v,
                    hours_worked_v, earnwt_v
                ))
            end

            year_v[n] = year
            month_v[n] = month
            wtfinl_v[n] = Float32(wtfinl_raw / 10_000.0)
            earnweek_v[n] = Float32(earnwk_raw / 100.0)
            age_v[n] = parse(Int32, @view line[29:30])
            sex_v[n] = parse(Int32, @view line[31:31])
            marst_v[n] = parse(Int32, @view line[32:32])
            empstat_v[n] = parse(Int32, @view line[33:34])
            occ1990_v[n] = parse(Int32, @view line[35:37])
            ind1990_v[n] = parse(Int32, @view line[38:40])
            classwkr_v[n] = parse(Int32, @view line[41:42])
            hours_worked_v[n] = Float32(hrs_raw)
            earnwt_v[n] = Float32(earnwt_raw / 10_000.0)
        end
    end

    return DataFrame(
        year = year_v[1:n],
        month = month_v[1:n],
        wtfinl = wtfinl_v[1:n],
        earnweek = earnweek_v[1:n],
        age = age_v[1:n],
        sex = sex_v[1:n],
        marst = marst_v[1:n],
        empstat = empstat_v[1:n],
        occ1990 = occ1990_v[1:n],
        ind1990 = ind1990_v[1:n],
        classwkr = classwkr_v[1:n],
        hours_worked = hours_worked_v[1:n],
        earnwt = earnwt_v[1:n],
    )
end

function clean_cps_labor(df::DataFrame)::DataFrame
    println("Cleaning CPS labor sample...")
    df = @subset(df, in.(:empstat, Ref([10, 12, 20, 21, 22])))
    df = @subset(df, 15 .<= :age .<= 64)
    df = @subset(df, :wtfinl .> 0)
    df = @subset(df, .!in.(:classwkr, Ref([21, 22, 23, 24, 25, 27, 28])))

    df[!, :occ_group] = map(classify_occ1990, df.occ1990)
    df[!, :ind_group] = map(classify_ind1990, df.ind1990)
    df = @subset(df, .!ismissing.(:occ_group), .!ismissing.(:ind_group))

    df[!, :occ_group] = convert(Vector{Int}, df.occ_group)
    df[!, :ind_group] = convert(Vector{Int}, df.ind_group)
    df[!, :date] = Date.(df.year, df.month, 1)
    return df
end

function clean_cps_earn(df::DataFrame)::DataFrame
    println("Cleaning CPS earnings sample...")
    df = @subset(df, in.(:empstat, Ref([10, 12])))
    df = @subset(df, 15 .<= :age .<= 64)
    df = @subset(df, :earnweek .> 0, :earnweek .!= 999999.99)

    function topcode_limit(year::Integer)
        year <= 1988 ? 999.0 :
        year <= 1997 ? 1923.0 : 2884.61
    end

    df = @transform(df, :earnweek = min.(:earnweek, topcode_limit.(:year)))
    df = @subset(df, :earnwt .> 0)
    df = @subset(df, 0 .< :hours_worked .<= 105)
    df = @subset(df, .!in.(:classwkr, Ref([21, 22, 23, 24, 25, 27, 28])))

    df[!, :occ_group] = map(classify_occ1990, df.occ1990)
    df[!, :ind_group] = map(classify_ind1990, df.ind1990)
    df = @subset(df, .!ismissing.(:occ_group), .!ismissing.(:ind_group))

    df[!, :occ_group] = convert(Vector{Int}, df.occ_group)
    df[!, :ind_group] = convert(Vector{Int}, df.ind_group)
    df[!, :female] = Int.(df.sex .== 2)
    df[!, :married] = Int.(in.(df.marst, Ref([1, 2])))
    df[!, :date] = Date.(df.year, df.month, 1)
    return df
end

function load_cpi()::DataFrame
    path = joinpath(DATA_DIR, "VARdata.xlsx")
    xf = XLSX.readxlsx(path)
    sh = xf["Monthly"]

    dates = Date[]
    cpi_vals = Float64[]

    for row in XLSX.eachrow(sh)
        XLSX.row_number(row) == 1 && continue
        d_raw = row[1]
        c_raw = row[7]
        (ismissing(d_raw) || ismissing(c_raw)) && continue
        d = parse_yearmonth(d_raw)
        isnothing(d) && continue
        push!(dates, d)
        push!(cpi_vals, Float64(c_raw))
    end

    cpi_df = DataFrame(date = dates, cpi = cpi_vals)
    cpi_df = @subset(cpi_df, DATE_START .<= :date .<= DATE_END)
    sort!(cpi_df, :date)
    return cpi_df
end

function deflate_wages!(cps::DataFrame, cpi::DataFrame)::DataFrame
    cps = leftjoin(cps, select(cpi, :date, :cpi), on=:date)
    cps = @subset(cps, .!ismissing.(:cpi))
    cps[!, :real_earnweek] = cps.earnweek ./ cps.cpi
    cps[!, :log_rincome] = log.(cps.real_earnweek)
    return cps
end

function add_log_rwage!(cps::DataFrame, variant::Symbol)::DataFrame
    cps[!, :log_rwage] = cps.log_rincome .- log.(cps.hours_worked)
    if variant == :hours
        cps[!, :log_hours] = log.(cps.hours_worked)
    end
    return cps
end

function weighted_quantile(x::AbstractVector, w::AbstractVector, q::Float64)::Float64
    idx = sortperm(x)
    xs = x[idx]
    ws = w[idx]
    cumw = cumsum(ws)
    cumw ./= cumw[end]
    i = searchsortedfirst(cumw, q)
    return xs[clamp(i, 1, length(xs))]
end

function build_panel_unemp(cps::DataFrame)::DataFrame
    gdf = groupby(cps, [:ind_group, :occ_group, :date, :year, :month])

    panel = combine(gdf) do sdf
        emp = sdf.empstat
        wt = sdf.wtfinl
        unemployed = sum(in.(emp, Ref((20, 21, 22))) .* wt)
        laborforce = sum(in.(emp, Ref((10, 12, 20, 21, 22))) .* wt)
        employed = sum(in.(emp, Ref((10, 12))) .* wt)
        unemp_rate = laborforce == 0 ? missing : unemployed / laborforce
        log_emp_count = employed > 0 ? log(employed) : missing
        pop_weight = sum(wt)
        (; unemp_rate, log_emp_count, pop_weight)
    end

    panel = @subset(panel, :pop_weight .>= 1000)
    return sort(panel, [:ind_group, :occ_group, :date])
end

function build_panel_earn(cps::DataFrame, variant::Symbol)::DataFrame
    gdf = groupby(cps, [:ind_group, :occ_group, :date, :year, :month])

    panel = if variant == :hourly_rate
        combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt)) => :log_rincome,
            [:log_rwage, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :log_rwage,
            [:age, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :age_mean,
            [:female, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :female_share,
            [:married, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :married_share,
            [:hours_worked, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :hours_mean,
            :log_rincome => length => :n_obs,
        )
    elseif variant == :hours
        combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt)) => :log_rincome,
            [:log_rwage, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :log_rwage,
            [:log_hours, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :log_hours,
            [:age, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :age_mean,
            [:female, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :female_share,
            [:married, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :married_share,
            :log_rincome => length => :n_obs,
        )
    elseif variant == :income
        combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt)) => :log_rincome,
            [:log_rwage, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :log_rwage,
            [:age, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :age_mean,
            [:female, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :female_share,
            [:married, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :married_share,
            [:hours_worked, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :hours_mean,
            :log_rincome => length => :n_obs,
        )
    elseif variant == :income_share_var
        p = combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt)) => :log_rincome,
            [:log_rwage, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :log_rwage,
            [:age, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :age_mean,
            [:female, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :female_share,
            [:married, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :married_share,
            [:hours_worked, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :hours_mean,
            [:earnweek, :earnwt] => ((e, wt) -> sum(e .* wt)) => :total_income_weekly,
            [:earnwt] => sum => :n_employed,
            :log_rincome => length => :n_obs,
        )
        p = @subset(p, :n_obs .>= 30)
        monthly_agg = combine(groupby(p, :date),
            :total_income_weekly => sum => :total_income_all,
            :n_employed => sum => :total_emp_all)
        p = leftjoin(p, monthly_agg, on=:date)
        p[!, :income_share] = p.total_income_weekly ./ p.total_income_all
        check = combine(groupby(p, :date), :income_share => sum => :sum_share)
        @assert all(isapprox.(check.sum_share, 1.0, atol=1e-6))
        p
    elseif variant == :inequality
        p = combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt)) => :log_rincome,
            [:log_rwage, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :log_rwage,
            [:log_rincome, :earnwt] => ((x, wt) -> weighted_quantile(x, wt, 0.25)) => :log_rincome_p25,
            [:log_rincome, :earnwt] => ((x, wt) -> weighted_quantile(x, wt, 0.75)) => :log_rincome_p75,
            [:age, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :age_mean,
            [:female, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :female_share,
            [:married, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :married_share,
            [:hours_worked, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :hours_mean,
            :log_rincome => length => :n_obs,
        )
        p = @subset(p, :n_obs .>= 30)
        p[!, :log_ratio_7525] = p.log_rincome_p75 .- p.log_rincome_p25
        p
    elseif variant == :median
        combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> weighted_quantile(w, wt, 0.5)) => :log_rincome_p50,
            [:log_rwage, :earnwt] => ((x, wt) -> weighted_quantile(x, wt, 0.5)) => :log_rwage_p50,
            [:age, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :age_mean,
            [:female, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :female_share,
            [:married, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :married_share,
            [:hours_worked, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :hours_mean,
            :log_rincome => length => :n_obs,
        )
    elseif variant == :unemployment || variant == :employment
        combine(gdf,
            [:log_rincome, :earnwt] => ((w, wt) -> sum(w .* wt) / sum(wt)) => :log_rincome,
            [:log_rwage, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :log_rwage,
            [:age, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :age_mean,
            [:female, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :female_share,
            [:married, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :married_share,
            [:hours_worked, :earnwt] => ((x, wt) -> sum(x .* wt) / sum(wt)) => :hours_mean,
            :log_rincome => length => :n_obs,
        )
    else
        error("Unknown variant in build_panel_earn: $variant")
    end

    if !(variant in (:income_share_var, :inequality))
        panel = @subset(panel, :n_obs .>= 30)
    end

    return sort(panel, [:ind_group, :occ_group, :date])
end

function load_oil_shock()::DataFrame
    path = joinpath(DATA_DIR, "oilSupplyNewsShocks_2025M06.xlsx")
    xf = XLSX.readxlsx(path)
    sh = xf["Monthly"]

    dates = Date[]
    shocks = Float64[]

    for row in XLSX.eachrow(sh)
        XLSX.row_number(row) == 1 && continue
        d_raw = row[1]
        s_raw = row[3]
        (ismissing(d_raw) || ismissing(s_raw)) && continue
        d = parse_yearmonth(d_raw)
        isnothing(d) && continue
        push!(dates, d)
        push!(shocks, Float64(s_raw))
    end

    df = DataFrame(date = dates, shock = shocks)
    df = @subset(df, DATE_START .<= :date .<= DATE_END)
    sort!(df, :date)
    return df
end

function load_oil_price()::DataFrame
    path = joinpath(DATA_DIR, "VARdata.xlsx")
    xf = XLSX.readxlsx(path)
    sh = xf["Monthly"]

    dates = Date[]
    prices = Float64[]

    for row in XLSX.eachrow(sh)
        XLSX.row_number(row) == 1 && continue
        d_raw = row[1]
        p_raw = row[2]
        (ismissing(d_raw) || ismissing(p_raw)) && continue
        d = parse_yearmonth(d_raw)
        isnothing(d) && continue
        push!(dates, d)
        push!(prices, Float64(p_raw))
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
        types=Dict(1 => Date, 2 => Float32),
        dateformat="yyyy-m-d")

    rename!(df, first(names(df)) => :date, names(df)[2] => colname)
    df = filter(row -> DATE_START <= row.date <= DATE_END, df)
    sort!(df, :date)
    return df
end

function lag_vec(v::AbstractVector, k::Int)
    T = Union{eltype(v), Missing}
    out = Vector{T}(missing, length(v))
    out[k+1:end] = v[1:end-k]
    return out
end

function build_macro_panel(shock_df, oil_df, ffr_df, cpi_df, indpro_df, t10y3m_df)::DataFrame
    macro_df = outerjoin(shock_df, oil_df, on=:date)
    macro_df = outerjoin(macro_df, ffr_df, on=:date)
    macro_df = outerjoin(macro_df, cpi_df, on=:date)
    macro_df = outerjoin(macro_df, indpro_df, on=:date)
    macro_df = outerjoin(macro_df, t10y3m_df, on=:date)
    sort!(macro_df, :date)

    for l in 1:L_LAG
        macro_df[!, Symbol("shock_lag$(l)")] = lag_vec(macro_df.shock, l)
    end
    macro_df[!, :log_oil_lag1] = lag_vec(macro_df.log_oil_price, 1)
    macro_df[!, :ffr_lag1] = lag_vec(macro_df.fedfunds, 1)
    macro_df[!, :cpi_lag1] = lag_vec(macro_df.cpi, 1)
    macro_df[!, :indpro_lag1] = lag_vec(macro_df.indpro, 1)
    macro_df[!, :t10y3m_lag1] = lag_vec(macro_df.t10y3m, 1)

    return dropmissing(macro_df)
end

function get_lp_specs(variant::Symbol)
    if variant == :hourly_rate
        return :log_rwage, [:age_mean, :female_share]
    elseif variant == :income
        return :log_rincome, [:age_mean, :female_share]
    elseif variant == :hours
        return :log_hours, [:age_mean, :female_share]
    elseif variant == :unemployment
        return :unemp_rate, [:age_mean, :female_share]
    elseif variant == :employment
        return :log_emp_count, [:age_mean, :female_share]
    elseif variant == :inequality
        return :log_ratio_7525, [:age_mean, :female_share]
    elseif variant == :median
        return :log_rincome_p50, [:age_mean, :female_share]
    elseif variant == :income_share_var
        return :income_share, [:age_mean, :female_share]
    else
        error("Unknown variant: $variant")
    end
end

function var_estim(y::Matrix{Float64}, p::Int, intercept::Bool)
    T, n_v = size(y)
    T_eff = T - p
    X_cols = [y[p-l+1:T-l, :] for l in 1:p]
    intercept && push!(X_cols, ones(T_eff, 1))
    X = hcat(X_cols...)
    Y = y[p+1:end, :]
    beta = (X'X) \ (X'Y)
    resid = Y - X * beta
    Sigma = (resid' * resid) ./ T_eff
    return beta, Sigma
end

function ic_var(y::Matrix{Float64}, p_max::Int, method::Int=2)
    n_v = size(y, 2)
    T = size(y, 1)
    T_eff = T - p_max
    bic = zeros(p_max)
    aic = zeros(p_max)
    for p in 1:p_max
        _, Sigma = var_estim(y, p, true)
        n_params = n_v^2 * p + n_v
        bic[p] = log(det(Sigma)) + n_params * log(T_eff) / T_eff
        aic[p] = log(det(Sigma)) + n_params * 2.0 / T_eff
    end
    return method == 1 ? argmin(aic) : argmin(bic)
end

function select_lag_length(panel::DataFrame; p_max::Int=L_LAG_MAX)
    shock_ts = sort(unique(panel[!, [:date, :shock]]), :date).shock
    y = reshape(Float64.(shock_ts), :, 1)
    p = ic_var(y, p_max, 2)
    @info "BIC-selected lag length: $p"
    return p
end

function driscoll_kraay_vcov(X::Matrix{Float64}, e::Vector{Float64},
    t_idx::Vector{Int}; m::Int=DK_BW)::Matrix{Float64}

    N, K = size(X)
    Tvals = sort(unique(t_idx))
    T = length(Tvals)
    tmap = Dict(v => i for (i, v) in enumerate(Tvals))
    H = zeros(K, T)

    for n in 1:N
        ti = tmap[t_idx[n]]
        H[:, ti] .+= X[n, :] .* e[n]
    end

    S = zeros(K, K)
    for t in 1:T
        S .+= H[:, t] * H[:, t]'
    end
    S ./= T

    for l in 1:m
        G = zeros(K, K)
        for t in (l+1):T
            G .+= H[:, t] * H[:, t-l]'
        end
        G ./= T
        w = 1 - l / (m + 1)
        S .+= w .* (G .+ G')
    end

    XX = inv(X'X)
    return XX * (T * S) * XX
end

function absorb_fe(df::DataFrame, v::Symbol)::Vector{Float64}
    fml = term(v) ~ term(1) + fe(term(:occ_group))
    fit = reg(df, fml, Vcov.simple(), save=:residuals)
    return residuals(fit)
end

function build_lp_data(panel::DataFrame, h::Int, y_var::Symbol, l_lag::Int)
    sort!(panel, [:occ_group, :date])
    out = DataFrame[]

    for g in groupby(panel, :occ_group)
        sub = sort(copy(g), :date)
        n = nrow(sub)
        leady = Vector{Union{Missing, Float64}}(missing, n)
        h < n && (leady[1:n-h] = sub[!, y_var][1+h:n])
        sub[!, :dep_var] = leady

        for l in 1:l_lag
            col = Symbol("ylag$(l)")
            lagged = Vector{Union{Missing, Float64}}(missing, n)
            l < n && (lagged[l+1:n] = sub[!, y_var][1:n-l])
            sub[!, col] = lagged
        end
        push!(out, sub)
    end

    df = vcat(out...)
    req = [:dep_var, :shock, :log_oil_lag1, :ffr_lag1, :cpi_lag1, :indpro_lag1, :t10y3m_lag1]
    for l in 1:l_lag
        push!(req, Symbol("shock_lag$(l)"))
        push!(req, Symbol("ylag$(l)"))
    end
    return dropmissing(df, req)
end

function build_cols!(df::DataFrame, controls::Vector{Symbol}, l_lag::Int)
    shock_cols = Symbol[]
    groups = sort(unique(df.occ_group))
    for g in groups
        c = Symbol("shock_g$(g)")
        df[!, c] = df.shock .* (df.occ_group .== g)
        push!(shock_cols, c)
    end

    lag_cols = [Symbol("shock_lag$(l)") for l in 1:l_lag]
    ylag_cols = [Symbol("ylag$(l)") for l in 1:l_lag]
    mac_controls = [:log_oil_lag1, :ffr_lag1, :cpi_lag1, :indpro_lag1, :t10y3m_lag1]
    xcols = vcat(shock_cols, lag_cols, ylag_cols, mac_controls, controls)
    return xcols, shock_cols
end

function estimate_lp(df_h::DataFrame, panel::DataFrame,
    controls::Vector{Symbol}, l_lag::Int)

    df = copy(df_h)
    xcols, _ = build_cols!(df, controls, l_lag)

    all_dates = sort(unique(panel.date))
    dmap = Dict(d => i for (i, d) in enumerate(all_dates))
    df[!, :time_idx] = Int.([dmap[d] for d in df.date])
    df[!, :occ_group] = categorical(df.occ_group)

    needed = vcat([:dep_var, :occ_group, :time_idx], xcols)
    df = dropmissing(df, needed)

    Y_tilde = absorb_fe(df, :dep_var)
    X_tilde = hcat([absorb_fe(df, c) for c in xcols]...)
    beta = (X_tilde'X_tilde) \ (X_tilde'Y_tilde)
    e = Y_tilde - X_tilde * beta
    V = driscoll_kraay_vcov(X_tilde, e, Int.(df.time_idx); m=DK_BW)
    se = sqrt.(diag(V))

    return (beta=beta, se=se, e=e, X=X_tilde, t_idx=Int.(df.time_idx), df=df, names=xcols)
end

function bootstrap_lp(df_h::DataFrame, panel::DataFrame,
    controls::Vector{Symbol}, l_lag::Int;
    n_boot::Int=N_BOOT,
    block_size::Int=BLOCK_SIZE,
    rng=MersenneTwister(BOOT_SEED))

    base = estimate_lp(df_h, panel, controls, l_lag)
    beta0 = base.beta
    se0 = base.se
    K = length(beta0)

    all_dates = sort(unique(base.df.date))
    T = length(all_dates)
    nb = ceil(Int, T / block_size)

    date_rows = Dict{Date, Vector{Int}}()
    for (i, r) in enumerate(eachrow(base.df))
        push!(get!(date_rows, r.date, Int[]), i)
    end

    B = fill(NaN, n_boot, K)
    TSTAT = fill(NaN, n_boot, K)

    for b in 1:n_boot
        starts = rand(rng, 1:T, nb)
        dd = Date[]
        for s in starts, k in 0:block_size-1
            push!(dd, all_dates[mod1(s + k, T)])
        end
        dd = dd[1:T]
        rows = reduce(vcat, [get(date_rows, d, Int[]) for d in dd])
        isempty(rows) && continue

        df_b = base.df[rows, :]
        df_b[!, :dep_var] = base.X[rows, :] * beta0 .+ base.e[rows]

        Y_b = absorb_fe(df_b, :dep_var)
        X_b = hcat([absorb_fe(df_b, c) for c in base.names]...)

        betab = try
            (X_b'X_b) \ (X_b'Y_b)
        catch
            continue
        end

        eb = Y_b - X_b * betab
        Vb = driscoll_kraay_vcov(X_b, eb, Int.(df_b.time_idx); m=DK_BW)
        seb = sqrt.(diag(Vb))

        B[b, :] = betab
        TSTAT[b, :] = (betab .- beta0) ./ seb
    end

    return B, TSTAT, beta0, se0, base.names
end

function run_h(panel::DataFrame, y_var::Symbol, controls::Vector{Symbol}, h::Int, l_lag::Int)
    df_h = build_lp_data(panel, h, y_var, l_lag)
    B, TSTAT, beta0, se0, names = bootstrap_lp(df_h, panel, controls, l_lag)

    rows = NamedTuple[]
    for k in eachindex(beta0)
        good = .!isnan.(TSTAT[:, k])
        if !any(good)
            ci_lo90 = NaN
            ci_hi90 = NaN
            ci_lo68 = NaN
            ci_hi68 = NaN
        else
            q95 = quantile(TSTAT[good, k], 0.95)
            q05 = quantile(TSTAT[good, k], 0.05)
            q84 = quantile(TSTAT[good, k], 0.84)
            q16 = quantile(TSTAT[good, k], 0.16)
            ci_lo90 = beta0[k] - se0[k] * q95
            ci_hi90 = beta0[k] - se0[k] * q05
            ci_lo68 = beta0[k] - se0[k] * q84
            ci_hi68 = beta0[k] - se0[k] * q16
        end

        push!(rows, (
            horizon = h,
            coef_name = string(names[k]),
            beta = beta0[k],
            se = se0[k],
            ci_lo90 = ci_lo90,
            ci_hi90 = ci_hi90,
            ci_lo68 = ci_lo68,
            ci_hi68 = ci_hi68,
        ))
    end

    return DataFrame(rows), B, TSTAT, names
end

function run_full_lp(panel::DataFrame, y_var::Symbol, controls::Vector{Symbol}, l_lag::Int;
    h_max::Int=H_MAX)

    out_dfs = DataFrame[]
    boot_store = Dict{Int, NamedTuple}()
    coef_names = Symbol[]

    for h in 0:h_max
        @info "  horizon h = $h / $h_max"
        df_res, B, TSTAT, names = run_h(panel, y_var, controls, h, l_lag)
        push!(out_dfs, df_res)
        boot_store[h] = (B=B, TSTAT=TSTAT, names=names)
        isempty(coef_names) && (coef_names = names)
    end

    return vcat(out_dfs...), boot_store, coef_names
end

function empty_irf_df()
    return DataFrame(
        horizon = collect(0:H_MAX),
        beta = Vector{Union{Missing, Float64}}(fill(missing, H_MAX + 1)),
        ci_lo90 = Vector{Union{Missing, Float64}}(fill(missing, H_MAX + 1)),
        ci_hi90 = Vector{Union{Missing, Float64}}(fill(missing, H_MAX + 1)),
        ci_lo68 = Vector{Union{Missing, Float64}}(fill(missing, H_MAX + 1)),
        ci_hi68 = Vector{Union{Missing, Float64}}(fill(missing, H_MAX + 1)),
    )
end

function extract_irf(df::DataFrame)
    out = Dict{Int, DataFrame}()
    for g in sort(collect(keys(OCC_LABELS)))
        cname = "shock_g$(g)"
        sub = filter(r -> r.coef_name == cname, df)
        out[g] = nrow(sub) == 0 ? empty_irf_df() :
            sub[:, [:horizon, :beta, :ci_lo90, :ci_hi90, :ci_lo68, :ci_hi68]]
    end
    return out
end

function save_irf_csvs(irfs::Dict{Int, DataFrame}, output_dir::String)
    for g in sort(collect(keys(OCC_LABELS)))
        label = get(OCC_LABELS, g, "group_$(g)")
        CSV.write(joinpath(output_dir, "irf_group$(g)_$(label).csv"), irfs[g])
    end
end

function get_global_ylims(irfs::Dict{Int, DataFrame})
    g_min = Inf
    g_max = -Inf

    for g in sort(collect(keys(OCC_LABELS)))
        df = irfs[g]
        nrow(df) == 0 && continue
        valid_lo = collect(skipmissing(df.ci_lo90))
        valid_hi = collect(skipmissing(df.ci_hi90))
        valid_beta = collect(skipmissing(df.beta))

        isempty(valid_lo) && isempty(valid_hi) && isempty(valid_beta) && continue

        if !isempty(valid_lo)
            g_min = min(g_min, minimum(valid_lo))
        end
        if !isempty(valid_hi)
            g_max = max(g_max, maximum(valid_hi))
        end
        if !isempty(valid_beta)
            g_min = min(g_min, minimum(valid_beta))
            g_max = max(g_max, maximum(valid_beta))
        end
    end

    if !isfinite(g_min) || !isfinite(g_max)
        return (-0.01, 0.01)
    end

    span = g_max - g_min
    padding = span == 0 ? 0.1 : span * 0.1
    return (g_min - padding, g_max + padding)
end

function plot_irf_grid(irfs::Dict{Int, DataFrame}, variant::Symbol, ind_group::Int, output_dir::String)
    fig = Figure(size=(1350, 1000), fontsize=12)
    ylims_all = get_global_ylims(irfs)
    ind_label = get(IND_LABELS, ind_group, "industry_$(ind_group)")
    title_text = "$(ind_label): $(get_plot_labels(variant))"
    Label(fig[0, :], title_text, fontsize=18, font=:bold)

    for (idx, g) in enumerate(sort(collect(keys(OCC_LABELS))))
        row_pos = ceil(Int, idx / 3)
        col_pos = mod1(idx, 3)
        label = get(OCC_LABELS, g, "group_$(g)")
        base_color = OCC_COLORS[g]

        ax = Axis(fig[row_pos, col_pos],
            title = "($(g)) $(label)",
            xlabel = "Horizon (months)",
            ylabel = col_pos == 1 ? "Beta" : "",
            xticks = 0:12:H_MAX,
        )
        ylims!(ax, ylims_all[1], ylims_all[2])
        hlines!(ax, [0.0], color=:black, linewidth=1, linestyle=:dash)

        df = irfs[g]
        valid_beta = collect(skipmissing(df.beta))
        if isempty(valid_beta)
            y_mid = (ylims_all[1] + ylims_all[2]) / 2
            text!(ax, [H_MAX / 2], [y_mid], text=["No data"], align=(:center, :center), color=:gray40)
            continue
        end

        plot_df = dropmissing(df, [:beta, :ci_lo90, :ci_hi90, :ci_lo68, :ci_hi68])
        horizons = Float64.(plot_df.horizon)
        betas = Float64.(plot_df.beta)
        lo90 = Float64.(plot_df.ci_lo90)
        hi90 = Float64.(plot_df.ci_hi90)
        lo68 = Float64.(plot_df.ci_lo68)
        hi68 = Float64.(plot_df.ci_hi68)

        band!(ax, horizons, lo90, hi90, color=(base_color, 0.25))
        band!(ax, horizons, lo68, hi68, color=(base_color, 0.45))
        lines!(ax, horizons, betas, color=base_color, linewidth=2.5)
    end

    output_path = joinpath(output_dir, "irf_occupations.pdf")
    save(output_path, fig)
    println("Saved IRF grid plot to: $output_path")
end

function build_occ_ind_panel(variant::Symbol)::DataFrame
    println("="^60)
    println("Building occ-within-industry panel for variant: $variant")
    println("="^60)

    download_gdrive_large(CPS_FILE_ID, CPS_PATH; expected_hash=CPS_SHA256)

    cps_raw = parse_cps(CPS_PATH)
    cpi_df = load_cpi()

    cps_labor = clean_cps_labor(cps_raw)
    cps_earn = clean_cps_earn(cps_raw)
    cps_earn = deflate_wages!(cps_earn, cpi_df)
    cps_earn = add_log_rwage!(cps_earn, variant)

    panel_labor = build_panel_unemp(cps_labor)
    panel_earn = build_panel_earn(cps_earn, variant)
    panel_all = leftjoin(panel_labor, panel_earn, on=[:ind_group, :occ_group, :date, :year, :month])

    shock_df = load_oil_shock()
    oil_df = load_oil_price()
    ffr_df = load_fred("FEDFUNDS.csv", :fedfunds)
    indpro_df = load_fred("INDPRO.csv", :indpro)
    t10y3m_df = load_fred("T10Y3M.csv", :t10y3m)
    macro_df = build_macro_panel(shock_df, oil_df, ffr_df, cpi_df, indpro_df, t10y3m_df)

    panel = leftjoin(panel_all, macro_df, on=:date)
    panel = dropmissing(panel, [:shock, :log_oil_lag1, :ffr_lag1, :cpi_lag1, :indpro_lag1, :t10y3m_lag1])
    sort!(panel, [:ind_group, :occ_group, :date])

    println("Panel complete: $(nrow(panel)) rows")
    return panel
end

function run_lp_for_industry(panel_ind::DataFrame, variant::Symbol, ind_group::Int, l_lag::Int)
    y_var, controls = get_lp_specs(variant)
    @info "Running LP for variant=$variant | industry=$(get(IND_LABELS, ind_group, ind_group)) | lags=$l_lag"

    results_df, boot_store, coef_names = run_full_lp(panel_ind, y_var, controls, l_lag)
    irfs = extract_irf(results_df)

    output_dir = get_industry_output_dir(variant, ind_group)
    CSV.write(joinpath(output_dir, "lp_coefficients.csv"), results_df)
    save_irf_csvs(irfs, output_dir)
    plot_irf_grid(irfs, variant, ind_group, output_dir)

    return results_df, boot_store, coef_names, irfs
end

function run_variant(variant::Symbol)
    panel = build_occ_ind_panel(variant)
    l_lag = select_lag_length(panel; p_max=L_LAG_MAX)

    for ind_group in sort(unique(panel.ind_group))
        panel_ind = subset(panel, :ind_group => x -> x .== ind_group)
        occ_count = length(unique(panel_ind.occ_group))
        if occ_count == 0
            @warn "Skipping industry $ind_group because no occupation groups remain."
            continue
        end
        run_lp_for_industry(panel_ind, variant, ind_group, l_lag)
    end
end

function run_all_variants(variants::Vector{Symbol}=Symbol[
    :hourly_rate,
    :hours,
    :income,
    :income_share_var,
    :inequality,
    :median,
    :unemployment,
    :employment,
])
    for variant in variants
        run_variant(variant)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_all_variants()
end
