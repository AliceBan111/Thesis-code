# =============================================================================
# oilshare_no_mining.jl
#
# Recompute occupation-level Bartik oil exposure after excluding
# industry group 2 (Mining), then run the same cross-sectional regression
# used by run_cross_sectional_occ_from_ind.
#
# Usage:
#   include("src/ind/01_data_prep.jl")
#   include("src/ind/02_lp_ind_estimation.jl")
#   include("src/ind/06_oilshare_no_mining.jl")
#
#   oilshare_no_mining = save_oilshare_no_mining(:hourly_rate)
#   cs_no_mining = run_cross_sectional_occ_no_mining(occ_irfs, oilshare_no_mining)
# =============================================================================

using CSV, DataFrames, DataFramesMeta
using Printf

"""
    build_oilshare_no_mining(cps_raw) -> DataFrame

Construct occupation-level Bartik oil exposure after dropping industry
group 2 (Mining). Shares are re-normalized within each occupation group
over the remaining industries.
"""
function build_oilshare_no_mining(cps_raw::DataFrame)::DataFrame
    println("Building OilShare (Bartik exposure), excluding Mining...")

    df = @subset(cps_raw, :empstat .∈ Ref([10, 12]))
    df = @subset(df, 15 .<= :age .<= 64)
    df = @subset(df, :wtfinl .> 0)
    df = @subset(df, :classwkr .∈ Ref([21, 22, 23, 24, 25, 27, 28]))

    df[!, :occ_group] = map(classify_occ1990, df.occ1990)
    df[!, :ind_group] = map(classify_ind1990, df.ind1990)

    df = @subset(df, .!ismissing.(:occ_group), .!ismissing.(:ind_group))
    df[!, :occ_group] = convert(Vector{Int}, df.occ_group)
    df[!, :ind_group] = convert(Vector{Int}, df.ind_group)

    # Drop Mining before computing occupation-industry shares, so shares are
    # re-normalized over non-mining employment within each occupation group.
    df = @subset(df, :ind_group .!= 2)

    crosstab = combine(
        groupby(df, [:occ_group, :ind_group]),
        :wtfinl => sum => :emp_weight
    )

    occ_total = combine(
        groupby(crosstab, :occ_group),
        :emp_weight => sum => :total_weight
    )
    crosstab = leftjoin(crosstab, occ_total, on = :occ_group)
    crosstab[!, :share] = crosstab.emp_weight ./ crosstab.total_weight

    crosstab[!, :oil_intensity] = [
        get(OIL_INTENSITY, row.ind_group, 0.0) for row in eachrow(crosstab)
    ]

    oilshare_df = combine(
        groupby(crosstab, :occ_group),
        [:share, :oil_intensity] =>
            ((s, oi) -> sum(s .* oi)) => :oil_exposure
    )
    sort!(oilshare_df, :occ_group)

    println("  OilShare by occupation group, excluding Mining:")
    for row in eachrow(oilshare_df)
        label = get(OCC_LABELS, row.occ_group, "group_$(row.occ_group)")
        @printf("    Group %2d %-35s %.4f\n", row.occ_group, label, row.oil_exposure)
    end

    return oilshare_df
end

"""
    save_oilshare_no_mining(variant=:hourly_rate) -> DataFrame

Parse the CPS file, recompute occupation oil exposure without Mining,
and save it beside the original `oilshare_by_occ.csv`.
"""
function save_oilshare_no_mining(variant::Symbol = :hourly_rate)::DataFrame
    download_gdrive_large(CPS_FILE_ID, CPS_PATH; expected_hash = CPS_SHA256)
    cps_raw = parse_cps(CPS_PATH)

    oilshare_df = build_oilshare_no_mining(cps_raw)

    output_dir = get_output_dir(variant)
    out_path = joinpath(dirname(output_dir), "oil_share_by_occ_no mining.csv")
    CSV.write(out_path, oilshare_df)
    println("Saved no-mining OilShare to: $out_path")

    return oilshare_df
end

"""
    run_cross_sectional_occ_no_mining(occ_irfs, oilshare_no_mining_df; h_range=1:36)

Run the occupation-level cross-sectional regression with the no-mining
oil exposure. This intentionally reuses the same implementation as
`run_cross_sectional_occ_from_ind` in `02_lp_ind_estimation.jl`.
"""
function run_cross_sectional_occ_no_mining(
        occ_irfs::Dict{Int,DataFrame},
        oilshare_no_mining_df::DataFrame;
        h_range::UnitRange = 1:36)::DataFrame

    return run_cross_sectional_occ_from_ind(
        occ_irfs,
        oilshare_no_mining_df;
        h_range = h_range
    )
end

const VARIANTS_TO_RUN_NO_MINING = [
    :hourly_rate,
    :hours,
    :income,
    :income_share_var,
    :inequality,
    :median,
    :unemployment,
    :employment,
]

"""
    run_cross_sectional_occ_no_mining_all(occ_irfs_by_variant; h_range=1:36)

Run the no-mining OilShare cross-sectional regression for all variants.

`occ_irfs_by_variant` should be a Dict where each key is a variant Symbol
and each value is that variant's occupation-level IRF Dict.
"""
function run_cross_sectional_occ_no_mining_all(
        occ_irfs_by_variant::Dict{Symbol,Dict{Int,DataFrame}};
        h_range::UnitRange = 1:36,
        variants_to_run::Vector{Symbol} = VARIANTS_TO_RUN_NO_MINING
    )::Dict{Symbol,DataFrame}

    oilshare_no_mining = save_oilshare_no_mining(first(variants_to_run))
    results = Dict{Symbol,DataFrame}()

    for variant in variants_to_run
        if !haskey(occ_irfs_by_variant, variant)
            @warn "Skipping $variant: occ_irfs_by_variant has no entry for this variant"
            continue
        end

        println("\n", "="^60)
        println("No-mining OilShare cross-sectional regression: $variant")
        println("="^60)

        cs_result = run_cross_sectional_occ_no_mining(
            occ_irfs_by_variant[variant],
            oilshare_no_mining;
            h_range = h_range
        )

        output_dir = get_output_dir(variant)
        out_path = joinpath(output_dir, "cross_sectional_occ_oilshare_no_mining.csv")
        CSV.write(out_path, cs_result)
        println("Saved no-mining cross-sectional result to: $out_path")

        results[variant] = cs_result
    end

    return results
end

"""
    run_cross_sectional_occ_no_mining_all(; h_range=1:36)

Load occupation IRFs from
`result/occ/analysis/merged_occ_irf_trajectories.csv`, then run the
no-mining OilShare cross-sectional regression for all variants.
"""
function run_cross_sectional_occ_no_mining_all(;
        h_range::UnitRange = 1:36,
        variants_to_run::Vector{Symbol} = VARIANTS_TO_RUN_NO_MINING,
        occ_irf_base_dir::String = joinpath(@__DIR__, "../..", "result", "occ")
    )::Dict{Symbol,DataFrame}

    oilshare_no_mining = save_oilshare_no_mining(first(variants_to_run))
    results = Dict{Symbol,DataFrame}()

    for variant in variants_to_run
        println("\n", "="^60)
        println("No-mining OilShare cross-sectional regression: $variant")
        println("="^60)

        # Left-hand side: occupation IRF trajectories from
        # result/occ/analysis/merged_occ_irf_trajectories.csv.
        occ_irfs = load_occ_irfs(variant; base_dir = occ_irf_base_dir)

        cs_result = run_cross_sectional_occ_no_mining(
            occ_irfs,
            oilshare_no_mining;
            h_range = h_range
        )

        output_dir = get_output_dir(variant)
        out_path = joinpath(output_dir, "cross_sectional_occ_oilshare_no_mining.csv")
        CSV.write(out_path, cs_result)
        println("Saved no-mining cross-sectional result to: $out_path")

        results[variant] = cs_result
    end

    return results
end
