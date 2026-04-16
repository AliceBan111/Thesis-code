include("src/occ/01_data_prep_revised.jl")
include("src/occ/02_lp_occ_estimation_revised.jl")
include("src/occ/03_significance_tests_revised.jl")
include("src/occ/04_plots_revised.jl")

cd(@__DIR__)

variants_to_run = [:hourly_rate, :hours, :income, :income_share_var, :inequality, :median, :unemployment, :employment]

for my_variant in variants_to_run
    panel = main(my_variant)

    results_df, boot_store, coef_names, irfs = run_lp(panel, my_variant)

    sig_table, pw_table = run_significance_tests(irfs, boot_store, coef_names, my_variant)

    run_visualization(irfs, my_variant, pw_table)
end