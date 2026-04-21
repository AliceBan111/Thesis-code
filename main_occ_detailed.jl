cd(@__DIR__)
include("src/occ/01_data_prep_detailed.jl")
include("src/occ/02_lp_occ_estimation_detailed.jl")

variants = [:hourly_rate, :hours, :income, :income_share_var,
            :inequality, :median, :unemployment, :employment]

panels = Dict(v => main(v) for v in variants)
master = run_all_variants(panels)   