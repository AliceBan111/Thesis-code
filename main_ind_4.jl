include("src/ind/01_data_prep_ind_4.jl")
include("src/ind/02_lp_ind_4_estimation.jl")
include("src/ind/03_significance_tests_ind_4.jl")
include("src/ind/04_plots_ind_4.jl")

cd(@__DIR__)

const VARIANTS_TO_RUN_IND4 = [
    :hourly_rate,
    :hours,
    :income,
    :income_share_var,
    :inequality,
    :median,
    :unemployment,
    :employment,
]

function run_ind_4(; variants = VARIANTS_TO_RUN_IND4, make_plots::Bool = true)
    for variant in variants
        println("===================================")
        println("Running ind_4: ", variant)
        println("===================================")

        panel, oilshare_df, ind_intensity_df = main(variant)
        results_df, boot_store, coef_names, irfs = run_lp(panel, variant)

        output_dir = get_output_dir(variant)
        cross_results = run_cross_sectional_ind(irfs, ind_intensity_df)
        CSV.write(joinpath(output_dir, "cross_sectional_ind_oil_intensity.csv"), cross_results)

        sig_table, pw_bh = run_significance_tests(irfs, boot_store, coef_names, variant)
        make_plots && run_plots(irfs, pw_bh, sig_table, variant)
    end

    include("src/ind/05_merge_irf_betas_ind_4.jl")
end

run_ind_4()
