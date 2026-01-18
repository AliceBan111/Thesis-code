# Include utility functions
include("src/ToolBox/load_var.jl")
include("src/data/build_macro_data.jl")
include("src/data/build_markup_data_method1.jl")
include("src/data/build_markup_data_method2.jl")
include("src/data/build_markup_data_method2_no_hp.jl")
include("src/data/build_markup_data_method2_no_hp_latest.jl")
include("src/models/var_sr.jl")
include("src/models/var.jl")

using Dates, Plots

cd(@__DIR__)

# 1. Load macro data
start_date_method1 = Date("1964-01-01")
end_date_method1   = Date("2017-12-31")

start_date_method2 = Date("1987-01-01")
end_date_method2   = Date("2012-12-31")

start_date_method2_latest = Date("1997-01-01")
end_date_method2_latest   = Date("2023-12-31")

macro_df_method1 = build_macro_data(start_date_method1, end_date_method1)
macro_df_method2 = build_macro_data(start_date_method2, end_date_method2)
macro_df_method2_no_hp = build_macro_data(start_date_method2, end_date_method2)
macro_df_method2_latest = build_macro_data(start_date_method2_latest, end_date_method2_latest)

# 2. Build markup data
df_method1 = build_markup_method1(macro_df_method1)
df_method2 = build_markup_method2(macro_df_method2)
df_method2_no_hp = build_markup_method2_no_hp(macro_df_method2_no_hp)
df_method2_no_hp_latest = build_markup_method2_no_hp_latest(macro_df_method2_latest, "data/KLEMS_latest.xlsx")

# 2.5 Plot VAR inputs
plot_VAR_inputs(df_method1; tag="method1")
plot_VAR_inputs(df_method2; tag="method2")
plot_VAR_inputs(df_method2_no_hp; tag="method2_no_hp")
plot_VAR_inputs(df_method2_no_hp_latest; tag="method2_no_hp_latest")

# 3. Estimate VAR + Sign Restrictions
res_method1 = estimate_VAR_SR(df_method1; method_name="method1")
res_method2 = estimate_VAR_SR(df_method2; method_name="method2")
res_method2_no_hp = estimate_VAR_SR(df_method2_no_hp; method_name="method2_no_hp")
res_method2_no_hp_latest = estimate_VAR_SR(df_method2_no_hp_latest; method_name="method2_no_hp_latest")