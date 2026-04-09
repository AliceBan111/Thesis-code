include("src/ind/01_data_prep.jl")
include("src/ind/02_lp_ind_estimation.jl")
include("src/ind/03_significance_tests.jl")
include("src/ind/04_plots.jl")

cd(@__DIR__)

variants_to_run = [:hourly_rate, :hours, :income, :income_share_var, :inequality, :median, :unemployment, :employment]

for my_variant in variants_to_run
    panel = main(my_variant)

    results_df, boot_store, coef_names, irfs = run_lp(panel, my_variant)

    sig_table, pw_bh = run_significance_tests(irfs, boot_store, coef_names, my_variant)

    run_plots(irfs, pw_bh, sig_table, my_variant)
end