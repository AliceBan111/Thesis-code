# =============================================================================
# 05_decomposition_R2.jl
# Bartik-style Shift-Share Decomposition with main() Interface
# Supports multi-dimension batch processing (mean_income, unemp_rate, etc.)
# =============================================================================

using CSV, DataFrames, DataFramesMeta, Statistics, Plots, Printf

# =============================================================================
# CONSTANTS & LABELS
# =============================================================================
const CUM_HORIZONS = 0:11

const OCC_LABELS = Dict(
    1 => "Managerial", 2 => "Professional_specialty", 3 => "High_tech",
    4 => "Sales", 5 => "Administrative_support", 6 => "Service",
    7 => "Farming_forestry_construction", 8 => "Precision_production_repair",
    9 => "Machine_operators_transport",
)

const IND_LABELS = Dict(
    1  => "Agriculture_forestry_fishing", 2  => "Mining", 3  => "Construction",
    4  => "Manufacturing_nondurable", 5  => "Manufacturing_durable", 
    6  => "Transportation_utilities", 7  => "Wholesale_trade", 8  => "Retail_trade",
    9  => "Finance_insurance_realestate", 10 => "Business_repair_services",
    11 => "Personal_entertainment_services", 12 => "Professional_related_services",
    13 => "Public_administration",
)

# =============================================================================
# CORE HELPERS
# =============================================================================
function extract_cumulative_irf(dir_path::String, labels::Dict{Int, String}, horizon_range::UnitRange)::DataFrame
    results = []
    for (gid, name) in labels
        fname = "irf_group$(gid)_$(name).csv"
        fpath = joinpath(dir_path, fname)
        !isfile(fpath) && continue
        
        df = CSV.read(fpath, DataFrame)
        df_sub = @subset(df, :horizon .∈ Ref(horizon_range))
        nrow(df_sub) == 0 && continue
        
        push!(results, (
            group_id = gid, name = name,
            cum_resp = sum(df_sub.beta_abs),
            cum_se   = sqrt(sum(df_sub.boot_se .^ 2))
        ))
    end
    return DataFrame(results)
end

function calculate_predicted_responses(shares_df::DataFrame, ind_resp_df::DataFrame)::DataFrame
    merged = leftjoin(shares_df, ind_resp_df, on = [:ind_group => :group_id])
    pred_df = combine(groupby(merged, :occ_group),
        :occ_name          => first       => :occ_name,
        [:share, :cum_resp] => ((s, r) -> sum(s .* r))                    => :predicted_response,
        [:share, :cum_se]  => ((s, se) -> sqrt(sum((s .* se) .^ 2)))      => :pred_se
    )
    return sort!(pred_df, :occ_group)
end

function compute_decomposition_and_R2(occ_resp_df::DataFrame, pred_df::DataFrame, dimension::String, decomp_dir::String)
    rename!(occ_resp_df, :group_id => :occ_group, :name => :occ_name)
    comp_df = innerjoin(occ_resp_df, pred_df, on = :occ_group, makeunique = true)
    
    comp_df.residual    = comp_df.cum_resp .- comp_df.predicted_response
    comp_df.residual_se = sqrt.(comp_df.cum_se .^ 2 .+ comp_df.pred_se .^ 2)
    
    total_var = var(comp_df.cum_resp)
    residual_var = var(comp_df.residual)
    r2 = total_var < 1e-10 ? NaN : 1.0 - residual_var / total_var
    comp_df[!, :R2] .= r2
    
    out_csv = joinpath(decomp_dir, "decomposition_$(dimension).csv")
    out_png = joinpath(decomp_dir, "actual_vs_predicted_$(dimension).png")
    CSV.write(out_csv, comp_df)
    
    p = scatter(comp_df.predicted_response, comp_df.cum_resp,
        xlabel = "Predicted (Industry Composition)", ylabel = "Actual (Estimated)",
        title = @sprintf("Responses vs. Predicted [%s] (R² = %.3f)", dimension, r2),
        label = "", markersize = 7, color = :steelblue, legend = :topleft, grid = true, size = (850, 600))
    
    xlims = extrema([comp_df.predicted_response; comp_df.cum_resp])
    margin = 0.1 * (xlims[2] - xlims[1])
    xlims = (xlims[1] - margin, xlims[2] + margin)
    plot!(p, [xlims[1], xlims[2]], [xlims[1], xlims[2]], linestyle = :dash, color = :gray, linewidth = 2, label = "45° line")
    for i in 1:nrow(comp_df)
        annotate!(p, comp_df.predicted_response[i], comp_df.cum_resp[i], text(comp_df.occ_name[i], 8, :left))
    end
    savefig(p, out_png)
    
    println(@sprintf("  R² = %.4f | Total Var = %.6f | Resid Var = %.6f", r2, total_var, residual_var))
    for row in eachrow(comp_df)
        println(@sprintf("    %-28s | Act: %+.4f | Pred: %+.4f | Res: %+.4f", row.occ_name, row.cum_resp, row.predicted_response, row.residual))
    end
    println("  Saved: $out_csv & $out_png\n")
    
    return comp_df, r2
end

function run_decomposition_for_dimension(dim::String, occ_dir::String, ind_dir::String, decomp_dir::String)
    shares_df = CSV.read(joinpath(decomp_dir, "occ_ind_shares_1994-1996.csv"), DataFrame)
    @assert nrow(shares_df) > 0 "Shares file missing in $decomp_dir"
    
    ind_resp = extract_cumulative_irf(ind_dir, IND_LABELS, CUM_HORIZONS)
    @assert nrow(ind_resp) == 13 "Expected 13 industries for $dim"
    
    occ_resp = extract_cumulative_irf(occ_dir, OCC_LABELS, CUM_HORIZONS)
    @assert nrow(occ_resp) == 9 "Expected 9 occupations for $dim"
    
    pred_resp = calculate_predicted_responses(shares_df, ind_resp)
    return compute_decomposition_and_R2(occ_resp, pred_resp, dim, decomp_dir)
end

# =============================================================================
# MAIN INTERFACE
# =============================================================================
"""
    main(; config=nothing, base_dir=joinpath(@__DIR__, "..", ".."))

Run Bartik-style decomposition for one or multiple dimensions.
Returns Dict{dimension => (comp_df, r2)}.
"""
function main(; config=nothing, base_dir=joinpath(@__DIR__, "..", ".."))
    decomp_dir = joinpath(base_dir, "result", "decomposition")
    mkpath(decomp_dir)
    
    if isnothing(config)
        config = Dict(
            "mean_income" => (
                occ_dir = joinpath(base_dir, "result", "occ", "income", "income_mean"),
                ind_dir = joinpath(base_dir, "result", "ind", "income")
            ),
            "unemp_rate" => (
                occ_dir = joinpath(base_dir, "result", "occ", "unemployment", "unemp"),
                ind_dir = joinpath(base_dir, "result", "ind", "unemployment")
            ),
            "hourly_wage" => (
                occ_dir = joinpath(base_dir, "result", "occ", "hourly_rate"),
                ind_dir = joinpath(base_dir, "result", "ind", "hourly_rate")
            ),
            "hours" => (
                occ_dir = joinpath(base_dir, "result", "occ", "hours"),
                ind_dir = joinpath(base_dir, "result", "ind", "hours")
            ),
            "inequality" => (
                occ_dir = joinpath(base_dir, "result", "occ", "inequality"),
                ind_dir = joinpath(base_dir, "result", "ind", "inequality")
            ),
            "income_share" => (
                occ_dir = joinpath(base_dir, "result", "occ", "income", "income_share"),
                ind_dir = joinpath(base_dir, "result", "ind", "income_share")
            ),
             "median" => (
                occ_dir = joinpath(base_dir, "result", "occ", "median"),
                ind_dir = joinpath(base_dir, "result", "ind", "median")
            ),
            "employment" => (
                occ_dir = joinpath(base_dir, "result", "occ",  "unemployment", "emp"),
                ind_dir = joinpath(base_dir, "result", "ind", "employment")
            )
        )
    end
    
    results = Dict{String, Tuple{DataFrame, Float64}}()
    println("🚀 Starting decomposition pipeline...")
    for (dim, paths) in config
        println("\n📦 Processing: $dim")
        try
            comp_df, r2 = run_decomposition_for_dimension(dim, paths.occ_dir, paths.ind_dir, decomp_dir)
            results[dim] = (comp_df, r2)
        catch e
            @warn "❌ Failed dimension '$dim': $e"
        end
    end
    println("\n✅ Pipeline complete. Processed $(length(results)) dimension(s).")
    return results
end

# =============================================================================
# ENTRY POINT
# =============================================================================
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

main()