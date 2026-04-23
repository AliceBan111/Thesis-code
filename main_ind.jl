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
