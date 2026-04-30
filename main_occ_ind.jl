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
    # :employment,
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



include("src/occ_ind/01_data_prep_occ_ind_q.jl")
include("src/occ_ind/02_lp_estimation_occ_ind_q.jl")
include("src/occ_ind/03_lp_plots_q.jl")

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

using DataFrames, CSV

function check_energy_intensive_obs(panel::DataFrame;
                                    ind_col::Symbol = :ind_group,
                                    energy_code::Int = 1,
                                    save_path::Union{String,Nothing}=nothing)

    # Energy-intensive only
    df = subset(panel, ind_col => ByRow(==(energy_code)))

    # Count how many non-missing values each occupation-quarter has
    info_cols = [
        :unemp_rate,
        :log_emp_count,
        :pop_weight,
        :n_obs,
        :log_rincome,
        :log_rwage,
        :log_hours,
        :hours_mean,
        :age_mean,
        :female_share,
        :married_share,
        :log_rincome_p25,
        :log_rincome_p75,
        :log_ratio_7525,
        :log_rincome_p50,
        :log_rwage_p50,
        :income_share,
        :shock,
        :log_oil_lag1,
        :ffr_lag1,
        :cpi_lag1,
        :indpro_lag1,
        :t10y3m_lag1,
    ]

    info_cols = [c for c in info_cols if c in propertynames(df)]

    out = combine(groupby(df, [:date, :occ_group]),
        nrow => :rows,
        :pop_weight => (x -> sum(skipmissing(x))) => :pop_weight_sum,
        [c => (x -> count(!ismissing, x)) => Symbol("valid_", c) for c in info_cols]...
    )

    sort!(out, [:date, :occ_group])

    if !isnothing(save_path)
        CSV.write(save_path, out)
        println("Saved to: ", save_path)
    end

    return out
end

energy_check = check_energy_intensive_obs(
    panel;
    save_path = joinpath(get_output_dir(:hourly_rate), "energy_intensive_occ_quarter_valid_counts.csv")
)



function check_energy_sample_size(panel::DataFrame;
                                  ind_col::Symbol = :ind_group,
                                  energy_code::Int = 1)

    df = subset(panel, ind_col => ByRow(==(energy_code)))

    specs = Any[
        nrow => :rows,
    ]

    if :pop_weight in propertynames(df)
        push!(specs, :pop_weight => (x -> sum(skipmissing(x))) => :pop_weight_sum)
    end

    if :n_obs in propertynames(df)
        push!(specs, :n_obs => (x -> sum(skipmissing(x))) => :n_obs_sum)
    end

    if :unemp_rate in propertynames(df)
        push!(specs, :unemp_rate => (x -> count(!ismissing, x)) => :valid_unemp_rate)
    end

    if :log_emp_count in propertynames(df)
        push!(specs, :log_emp_count => (x -> count(!ismissing, x)) => :valid_log_emp_count)
    end

    out = combine(groupby(df, [:date, :occ_group]), specs...)
    sort!(out, [:date, :occ_group])
    return out
end

sample_check = check_energy_sample_size(panel)

CSV.write(
    joinpath(get_output_dir(:hourly_rate), "energy_intensive_occ_quarter_sample_size.csv"),
    sample_check
)