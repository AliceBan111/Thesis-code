# =============================================================================
# 08_pattern_matrix_heatmap.jl
# Step 1 & 2: Build occupational response pattern matrix & visualize via heatmap
# Extracts cumulative responses (0-11 months) across 8 dimensions, 
# standardizes by dimension, and generates a publication-ready heatmap.
# =============================================================================

using CSV, DataFrames, DataFramesMeta, Statistics, Plots, Printf

# =============================================================================
# 0. PATHS & CONSTANTS
# =============================================================================
const BASE_DIR   = joinpath(@__DIR__, "..", "..", "result", "occ")
const DECOMP_DIR = joinpath(BASE_DIR, "decomposition")
mkpath(DECOMP_DIR)

# Cumulative horizon window (matches decomposition logic: 1-year effect)
const CUM_HORIZONS = 0:11

# Occupation labels (must match file naming)
const OCC_LABELS = Dict(
    1 => "Managerial",
    2 => "Professional_specialty",
    3 => "High_tech",
    4 => "Sales",
    5 => "Administrative_support",
    6 => "Service",
    7 => "Farming_forestry_construction",
    8 => "Precision_production_repair",
    9 => "Machine_operators_transport",
)

# Dimension directories (paths relative to project root)
const DIM_PATHS = Dict(
    "hourly_wage"   => joinpath(BASE_DIR, "hourly_rate"),
    "hours"         => joinpath(BASE_DIR, "hours"),
    "mean_income"   => joinpath(BASE_DIR, "income", "income_mean"),
    "median_income" => joinpath(BASE_DIR, "median"),
    "income_share"  => joinpath(BASE_DIR, "income", "income_share"),
    "ratio_75_25"   => joinpath(BASE_DIR, "inequality"), 
    "unemp_rate"    => joinpath(BASE_DIR, "unemployment", "unemp"),
    "emp_number"    => joinpath(BASE_DIR, "unemployment", "emp"),
)

# =============================================================================
# 1. EXTRACT CUMULATIVE RESPONSES PER DIMENSION
# =============================================================================
"""
    extract_cumulative_for_dimension(dir_path, horizon_range) -> Vector{Float64}
Reads 9 IRF CSVs, sums beta_abs over horizon_range, returns 9-element vector.
"""
function extract_cumulative_for_dimension(dir_path::String, horizon_range::UnitRange)::Vector{Float64}
    cum_vals = Float64[]
    for occ in 1:9
        fname = "irf_group$(occ)_$(OCC_LABELS[occ]).csv"
        fpath = joinpath(dir_path, fname)
        
        if !isfile(fpath)
            @warn "File not found: $fpath -- skipping occupation $occ"
            push!(cum_vals, NaN)
            continue
        end
        
        df = CSV.read(fpath, DataFrame)
        df_sub = @subset(df, :horizon .∈ Ref(horizon_range))
        
        if nrow(df_sub) == 0
            @warn "No horizons in range for $fname"
            push!(cum_vals, NaN)
            continue
        end
        
        push!(cum_vals, sum(df_sub.beta_abs))
    end
    return cum_vals
end

# =============================================================================
# 2. BUILD & STANDARDIZE PATTERN MATRIX
# =============================================================================
function build_and_standardize_matrix()::DataFrame
    println("Building occupational response pattern matrix...")
    
    pattern_df = DataFrame(
        occ_group = 1:9,
        occ_name  = [OCC_LABELS[i] for i in 1:9]
    )
    
    for (dim_name, dim_dir) in DIM_PATHS
        println("  Loading $dim_name from $dim_dir")
        pattern_df[!, dim_name] = extract_cumulative_for_dimension(dim_dir, CUM_HORIZONS)
    end
    
    std_df = copy(pattern_df)
    
    for col in names(pattern_df)
        col in ["occ_group", "occ_name"] && continue
        
        vals = Float64.(pattern_df[!, col])
        μ = mean(vals)
        σ = std(vals)
        
        std_df[!, col * "_std"] = σ > 1e-10 ? (vals .- μ) ./ σ : zeros(length(vals))
    end
    
    CSV.write(joinpath(DECOMP_DIR, "pattern_matrix_raw.csv"), pattern_df)
    CSV.write(joinpath(DECOMP_DIR, "pattern_matrix_standardized.csv"), std_df)
    println("Saved pattern matrices to $DECOMP_DIR")
    
    return std_df
end

# =============================================================================
# 3. HEATMAP VISUALIZATION
# =============================================================================
const OCC_SHORT_LABELS = Dict(
    1 => "Managerial",
    2 => "Professional",
    3 => "High-tech",
    4 => "Sales",
    5 => "Admin Support",
    6 => "Service",
    7 => "Farming/Constr.",
    8 => "Precision Prod.",
    9 => "Machine/Trans.",
)

function plot_response_heatmap(std_df::DataFrame)
    println("Generating standardized response heatmap...")
    
    std_cols = [Symbol(c, "_std") for c in keys(DIM_PATHS)]
    mat = Matrix(std_df[:, std_cols])'

    x_labels = [OCC_SHORT_LABELS[i] for i in 1:9]
    y_labels = [replace(String(c), "_std" => "") for c in std_cols]
    
    p = heatmap(mat,
        xticks = (1:9, x_labels),
        yticks = (1:length(std_cols), y_labels),
        xrotation = 30,          # 倾斜角度
        xtickfontsize = 9,
        ytickfontsize = 10,
        size = (1200, 700),
        left_margin  = 10Plots.mm,
        bottom_margin = 20Plots.mm,  # 给倾斜标签留空间
        xlabel = "Occupation",
        ylabel = "Response Dimension",
        title  = "Standardized Cumulative Responses (Horizon 0-11)",
        color  = :diverging_bkr_55_10_c35_n256,
        grid   = false,
        framestyle = :box
    )
    
    for i in 1:size(mat, 1), j in 1:size(mat, 2)
        val = mat[i, j]
        txt_color = abs(val) > 1.5 ? :white : :black
        annotate!(p, j, i, text(@sprintf("%.2f", val), 7, txt_color, :center))
    end
    
    hline!(p, collect(1.5:1.0:(length(std_cols)-0.5)), color=:gray, linewidth=0.5, label="")
    
    plot_path = joinpath(DECOMP_DIR, "response_pattern_heatmap.png")
    savefig(p, plot_path)
    println("Heatmap saved: $plot_path")
    
    return p
end

# =============================================================================
# 4. MAIN
# =============================================================================
function main()
    std_df = build_and_standardize_matrix()
    plot_response_heatmap(std_df)
    
    println("\n" * "="^60)
    println("MATRIX PREVIEW (Standardized)")
    println("="^60)
    # Show only standardized columns
    std_cols = [Symbol(c, "_std") for c in keys(DIM_PATHS)]
    show(std_df[:, [:occ_name; std_cols]], allrows=true, allcols=true)
    println("\n" * "="^60)
    println("Interpretation guide:")
    println("• Blue cells (<0): below-average response for that dimension")
    println("• Red cells  (>0): above-average response for that dimension")
    println("• Look for vertical clustering (similar occupations) and")
    println("  horizontal clustering (similar dimensions).")
    println("="^60)
end

# Run script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

main()