include("src/data/data_prep_occ.jl")
# include("src/data/data_prep_edu.jl")
include("src/wage_shock_analysis/02_lp_estimation_occ.jl")
# include("src/wage_shock_analysis/02_lp_estimation_edu.jl")
include("src/wage_shock_analysis/03_significance_tests_occ.jl")
#include("src/wage_shock_analysis/04_robustness.jl")
include("src/wage_shock_analysis/05_plots.jl")

cd(@__DIR__)


# ── Step 1: Data preparation ──────────────────────────────────
panel = main()   # from 01_data_prep.jl
describe(DataFrame(unemp_rate = panel.unemp_rate))

# ── Step 2: LP estimation ─────────────────────────────────────
results_df, boot_store, coef_names, irfs = run_estimation(panel)    # from 02_lp_estimation.jl


# ── Step 3: Significance tests ────────────────────────────────
sig_table, pw_bh = run_significance_tests(irfs, boot_store, coef_names)   # from 03_significance_tests.jl

# ── Step 4: Robustness checks ─────────────────────────────────
# run_all_robustness(panel)   # from 04_robustness.jl

# ── Step 5: Plots ─────────────────────────────────────────────
run_plots(irfs, sig_table, pw_bh)   # from 05_plots.jl
