using CSV, DataFrames, DataFramesMeta
using CairoMakie
using Dates, Statistics
using Printf

# =============================================================================
# 0. PATHS
# =============================================================================
const DATA_PATH = joinpath(@__DIR__, "data", "cps_00017.dat")
const PIE_DIR   = joinpath(@__DIR__, "result", "industry_pie")
mkpath(PIE_DIR)

# =============================================================================
# 1. LABELS & CLASSIFIERS
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

const IND_LABELS = Dict(
    1  => "Agriculture",
    2  => "Mining",
    3  => "Construction",
    4  => "Mfg_nondurable",
    5  => "Mfg_durable",
    6  => "Transport_util",
    7  => "Wholesale",
    8  => "Retail",
    9  => "Finance_RE",
    10 => "Business_svc",
    11 => "Personal_entmt",
    12 => "Professional_svc",
    13 => "Public_admin",
)

const IND_COLORS = [
    colorant"#4e79a7", colorant"#f28e2b", colorant"#e15759",
    colorant"#76b7b2", colorant"#59a14f", colorant"#edc948",
    colorant"#b07aa1", colorant"#ff9da7", colorant"#9c755f",
    colorant"#bab0ac", colorant"#d37295", colorant"#fabfd2",
    colorant"#8cd17d",
]

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
# 2. PARSE FIXED-WIDTH FILE
# 列位置（根据你的 data dictionary）：
#   YEAR      1-4
#   MONTH     5-6
#   WTFINL    7-20   (4 implied decimals)
#   AGE       21-22
#   OCC1990   27-29
#   IND1990   30-32
# =============================================================================
function parse_cps_minimal(path::String)::DataFrame
    println("Parsing: $path")
    println("File size: ", round(filesize(path)/1e9, digits=2), " GB")

    est_rows = max(1_000_000, round(Int, filesize(path) / 45))
    year_v   = Vector{Int32}(undef, est_rows)
    age_v    = Vector{Int32}(undef, est_rows)
    wtfinl_v = Vector{Float32}(undef, est_rows)
    classwkr_v  = Vector{Int32}(undef, est_rows)
    occ_v    = Vector{Int32}(undef, est_rows)
    ind_v    = Vector{Int32}(undef, est_rows)
    n = 0

    open(path, "r") do f
        for line in eachline(f)
            length(line) < 32 && continue

            year = parse(Int32, @view line[1:4])

            age_raw  = tryparse(Int32,   strip(@view line[21:22]))
            isnothing(age_raw) && continue
            (age_raw < 16 || age_raw > 64) && continue

            wtfinl_raw = tryparse(Float64, strip(@view line[7:20]))
            isnothing(wtfinl_raw) || wtfinl_raw <= 0 && continue
            
            clswkr_raw = tryparse(Int32, strip(@view line[33:34]))
            isnothing(clswkr_raw) && continue

            occ_raw = tryparse(Int32, strip(@view line[27:29]))
            isnothing(occ_raw) && continue

            ind_raw = tryparse(Int32, strip(@view line[30:32]))
            isnothing(ind_raw) && continue

            n += 1
            if n > length(year_v)
                new_cap = round(Int, length(year_v) * 1.5)
                foreach(v -> resize!(v, new_cap),
                    (year_v, age_v, wtfinl_v, occ_v, ind_v, classwkr_v))
            end

            year_v[n]   = year
            age_v[n]    = age_raw
            wtfinl_v[n] = Float32(wtfinl_raw / 10_000.0)
            occ_v[n]    = occ_raw
            ind_v[n]    = ind_raw
            classwkr_v[n] = clswkr_raw
        end
    end

    df = DataFrame(
        year   = year_v[1:n],
        age    = age_v[1:n],
        wtfinl = wtfinl_v[1:n],
        classwkr = classwkr_v[1:n],
        occ1990 = occ_v[1:n],
        ind1990 = ind_v[1:n],
    )
    println("  Rows parsed: ", nrow(df))
    return df
end

# =============================================================================
# 3. CLASSIFY & FILTER
# =============================================================================
function prepare(df::DataFrame)::DataFrame
    df[!, :occ_group] = map(classify_occ1990, df.occ1990)
    df[!, :ind_group] = map(classify_ind1990, df.ind1990)
    df = @subset(df,
        .!ismissing.(:occ_group),
        .!ismissing.(:ind_group),
        :wtfinl .> 0,
        :classwkr .∈ Ref([21, 22, 23, 24, 25, 27, 28]),
    )
    df[!, :occ_group] = convert(Vector{Int}, df.occ_group)
    df[!, :ind_group] = convert(Vector{Int}, df.ind_group)
    return df
end

# =============================================================================
# 4. PIE CHARTS
# =============================================================================
function plot_pies(df::DataFrame)
    for occ_g in 1:9
        occ_label = OCC_LABELS[occ_g]
        sub = @subset(df, :occ_group .== occ_g)
        nrow(sub) == 0 && continue

        # 加权占比
        gdf    = groupby(sub, :ind_group)
        counts = combine(gdf, :wtfinl => sum => :wt_sum)
        total  = sum(counts.wt_sum)
        counts[!, :share] = counts.wt_sum ./ total
        sort!(counts, :ind_group)
        counts = @subset(counts, 1 .<= :ind_group .<= 13)

        shares = counts.share
        labels = [get(IND_LABELS, g, "Other") for g in counts.ind_group]
        colors = [IND_COLORS[g]               for g in counts.ind_group]

        fig = Figure(size = (950, 580))
        ax  = Axis(fig[1, 1],
            title     = occ_label,
            titlesize = 15,
            aspect    = DataAspect(),
        )
        hidedecorations!(ax)
        hidespines!(ax)

        # 扇形
        angles = cumsum([0.0; shares .* 2π])
        for i in 1:length(shares)
            θ1, θ2 = angles[i], angles[i+1]
            Δθ = θ2 - θ1
            Δθ < 1e-4 && continue

            n_pts = max(3, round(Int, Δθ / 0.02))
            θs    = range(θ1, θ2, length = n_pts)
            xs    = vcat(0.0, cos.(θs))
            ys    = vcat(0.0, sin.(θs))

            poly!(ax, Point2f.(xs, ys);
                  color       = colors[i],
                  strokecolor = :white,
                  strokewidth = 1.5)

            # 百分比标注
            if shares[i] >= 0.03
                θ_mid = (θ1 + θ2) / 2
                r_txt = shares[i] > 0.12 ? 0.55 : 0.72
                text!(ax, @sprintf("%.1f%%", shares[i] * 100);
                      position = Point2f(r_txt * cos(θ_mid), r_txt * sin(θ_mid)),
                      align    = (:center, :center),
                      fontsize = shares[i] > 0.08 ? 11 : 9,
                      color    = :white,
                      font     = :bold)
            end
        end

        # Legend
        Legend(fig[1, 2],
            [PolyElement(color = colors[i]) for i in 1:length(shares)],
            [@sprintf("%s (%.1f%%)", labels[i], shares[i]*100) for i in 1:length(shares)];
            framevisible = false,
            labelsize    = 11,
            rowgap       = 4,
        )

        out = joinpath(PIE_DIR, "pie_occ$(occ_g)_$(occ_label).pdf")
        save(out, fig)
        println("Saved: $out")
    end
end

# =============================================================================
# 5. MAIN
# =============================================================================
df_raw  = parse_cps_minimal(DATA_PATH)
df_clean = prepare(df_raw)
plot_pies(df_clean)
