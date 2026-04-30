include("src/ind/01_data_prep.jl")
include("src/ind/02_lp_ind_estimation.jl")


cd(@__DIR__)

variants_to_run = [:hourly_rate, :hours, :income, :income_share_var, :inequality, :median, :unemployment, :employment]

for my_variant in variants_to_run
    panel, oilshare_df, ind_intensity_df = main(my_variant)

    results_df, boot_store, coef_names, irfs = run_lp(panel, my_variant)

    cross_results = run_cross_sectional_ind(irfs, ind_intensity_df)

    occ_irfs = load_occ_irfs(my_variant)
    cross_results_occ = run_cross_sectional_occ_from_ind(occ_irfs, oilshare_df)
end


using CSV, DataFrames
ind_intensity_df = CSV.read("result/ind_v2/oil_intensity_by_ind.csv", DataFrame)
irf_df = CSV.read("result/ind_v2/merged_irf_betas.csv", DataFrame)
ind_intensity_df_filtered = filter(
    r -> r.ind_group != 2,
    ind_intensity_df
)

for my_variant in variants_to_run

    println("===================================")
    println("Running: ", my_variant)
    println("===================================")

    # ----------------------------------
    # 当前 outcome
    # ----------------------------------
    df_sub = filter(
        r -> r.outcome == string(my_variant),
        irf_df
    )

    # ----------------------------------
    # 构造 irfs
    # ----------------------------------
    irfs = Dict{Int,DataFrame}()

    for g in 1:12
        irfs[g] = DataFrame(
            horizon = df_sub.horizon,
            beta    = Float64.(df_sub[!, Symbol(string(g))])
        )
    end

    # ----------------------------------
    # 剔除 group 2
    # ----------------------------------
    irfs_filtered = Dict(
        k => v for (k,v) in irfs if k != 2
    )

    # ----------------------------------
    # 跑回归
    # ----------------------------------
    cross_results = run_cross_sectional_ind(
        irfs_filtered,
        ind_intensity_df_filtered
    )
end