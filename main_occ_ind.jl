include("src/occ_ind/01_data_prep_occ_ind.jl")
include("src/occ_ind/02_lp_estimation_occ_ind.jl")
include("src/occ_ind/03_lp_plots.jl")

cd(@__DIR__)

variants_to_run = [
    :hourly_rate,
    :hours,
    :income,
    :income_share_var,
    :inequality,
    :median,
    :unemployment,
    :employment,
]

for my_variant in variants_to_run
    @info "===== Variant: $my_variant ====="

    # Step 1: build (occ × ind × date) panel
    panel = main(my_variant)

    # Step 2: run LP for all 4 merged industry groups
    # Returns Dict{Int, DataFrame}: ind_merged => lp_coefficients
    all_results = run_lp(panel, my_variant)

    # Step 3: plot IRF grids
    run_visualization(all_results, my_variant)

    @info "===== Done: $my_variant | industries with results: $(length(all_results)) ====="
end