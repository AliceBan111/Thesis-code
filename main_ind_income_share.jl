include("src/industry/data_prep_ind_income.jl")
include("src/industry/02_lp_estimation_ind_income_share.jl")
include("src/industry/03_significance_tests_ind.jl")
include("src/industry/05_plots_ind.jl")

cd(@__DIR__)


# ── Step 1: Data preparation ──────────────────────────────────
panel = main()  # from 01_data_prep.jl
# describe(DataFrame(unemp_rate = panel.unemp_rate))


# ── Step 2: LP estimation ─────────────────────────────────────
results_df, boot_store, coef_names, irfs = run_estimation(panel)    # occ

# ── Step 3: Significance tests ────────────────────────────────
sig_table, pw_bh = run_significance_tests(irfs, boot_store, coef_names)   # occ

# ── Step 4: Robustness checks ─────────────────────────────────
# run_all_robustness(panel)   # from 04_robustness.jl

# ── Step 5: Plots ─────────────────────────────────────────────
run_plots(irfs, sig_table, pw_bh)   # from 05_plots.jl