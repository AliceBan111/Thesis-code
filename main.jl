# Include utility functions
include("src/ToolBox/load_var.jl")
include("src/data/build_macro_data.jl")
include("src/data/build_markup_data_method2.jl")
include("src/models/var_sr.jl")
include("src/wage_shock_analysis/extract_shocks.jl")  
include("src/wage_shock_analysis/wage_heterogeneity.jl")

using CSV, DataFrames, Statistics, Plots, Dates, Downloads, SHA

cd(@__DIR__)

# 1. Load macro data
start_date_method2 = Date("1987-01-01")
end_date_method2   = Date("2012-12-31")

macro_df_method2 = build_macro_data(start_date_method2, end_date_method2)

# 2. Build markup data
df_method2 = build_markup_method2(macro_df_method2)

# 3. Estimate VAR + Sign Restrictions
res_method2 = estimate_VAR_SR(df_method2; method_name="method2")

# 4. extract shocks
markup_shocks_m2 = extract_structural_shocks(res_method2, shock_index=1)
n_shocks = length(markup_shocks_m2)

# Get the corresponding dates for the shocks
shock_dates = last(df_method2.observation_date, n_shocks)

println("Shocks length: ", n_shocks)
println("Dates length: ", length(shock_dates))

# 5. Run heterogeneity analysis
df_results, irf_data = run_wage_heterogeneity_analysis(
    markup_shocks_m2, 
    shock_dates; 
    output_suffix="method2"
)