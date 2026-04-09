# =============================================================================
# 11_feature_analysis.jl
# Step 1 & 2: Spearman Correlation Analysis & Visualization
# Analyzes the relationship between decomposition residuals and task features.
# Designed for small sample (N=9 broad occupation groups).
# =============================================================================

using CSV, DataFrames, DataFramesMeta, Statistics, Plots, Printf, StatsBase

# =============================================================================
# 0. PATHS & CONFIGURATION
# =============================================================================
const BASE_DIR   = joinpath(@__DIR__, "..", "..")
const DECOMP_DIR = joinpath(BASE_DIR, "result", "occ", "decomposition")

const RESIDUAL_FILE = joinpath(BASE_DIR, "result", "decomposition", "decomposition_results.csv")
const TASK_FILE     = joinpath(DECOMP_DIR, "task_features_aggregated.csv")
const OUTPUT_DIR    = DECOMP_DIR  # Save plots and tables in the same directory

# =============================================================================
# 1. LOAD & MERGE DATA
# =============================================================================
function load_and_prepare_data()::DataFrame
    println("Loading decomposition results and task features...")
    
    # Load residuals (actual - predicted)
    if !isfile(RESIDUAL_FILE)
        error("File not found: $RESIDUAL_FILE. Please run decomposition script first.")
    end
    decomp = CSV.read(RESIDUAL_FILE, DataFrame)
    
    # Load task features (Dorn aggregated)
    if !isfile(TASK_FILE)
        error("File not found: $TASK_FILE. Please run task aggregation script first.")
    end
    tasks = CSV.read(TASK_FILE, DataFrame)
    
    # Check column names (handle potential variations)
    if !("residual" in names(decomp))
        # Try to find the column
        residual_col = findfirst(contains("resid"), names(decomp))
        if isnothing(residual_col)
            error("Could not find 'residual' column in decomposition results.")
        end
        rename!(decomp, residual_col => :residual)
    end
    
    # Merge on occ_group (inner join to ensure only valid matches)
    merged = innerjoin(decomp, tasks, on = :occ_group, makeunique = true)
    
    # Sort for consistent ordering
    sort!(merged, :occ_group)
    
    println("  Merged data: $(nrow(merged)) occupations")
    println("  Columns: $(names(merged))")
    
    return merged
end

# =============================================================================
# 2. SPEARMAN CORRELATION ANALYSIS (Step 1)
# =============================================================================
function compute_spearman_correlations(df::DataFrame)::DataFrame
    println("\n" * "="^60)
    println("STEP 1: SPEARMAN RANK CORRELATIONS (N=$(nrow(df)))")
    println("="^60)
    
    task_cols = [:task_abstract, :task_routine, :task_manual]
    results = []
    
    for col in task_cols
        mask = .!ismissing.(df.residual) .& .!ismissing.(df[!, col])
        x = df.residual[mask]
        y = df[!, col][mask]
        # Compute Spearman correlation (robust to outliers, monotonic relationships)
        r = cor(tiedrank(x), tiedrank(y))
        
        # Store results
        push!(results, (
            feature = String(col),
            spearman_rho = r,
            interpretation = r > 0.5 ? "Strong Positive" : 
                             r < -0.5 ? "Strong Negative" :
                             abs(r) > 0.3 ? "Moderate" : "Weak"
        ))
        
        # Print formatted output
        sign_str = r >= 0 ? "+" : ""
        println(@sprintf("  %-20s : ρ = %s%.3f  (%s)", col, sign_str, r, 
                         r > 0.5 ? "Strong Positive" : 
                         r < -0.5 ? "Strong Negative" : "Moderate/Weak"))
    end

    println("\n" * "="^60)
    println("LINEAR REGRESSION R² (Simple OLS)")
    println("="^60)
    
    for col in task_cols
        # 过滤缺失值
        mask = .!ismissing.(df.residual) .& .!ismissing.(df[!, col])
        x = df[!, col][mask]
        y = df.residual[mask]
        
        # Pearson r 和 R²
        r_pearson = cor(x, y)
        r_squared = r_pearson^2
        
        # 简单 OLS: y = β₀ + β₁*x
        β₁ = cov(x, y) / var(x)
        β₀ = mean(y) - β₁ * mean(x)
        
        println(@sprintf("  %-20s : R² = %.3f  (β = %+.4f)", 
                         col, r_squared, β₁))
    end
    
    results_df = DataFrame(results)
    
    # Save to CSV
    CSV.write(joinpath(OUTPUT_DIR, "spearman_correlations.csv"), results_df)
    println("\n✅ Saved correlations to: $(joinpath(OUTPUT_DIR, "spearman_correlations.csv"))")
    
    return results_df
end

# =============================================================================
# 3. VISUALIZATION: SCATTER PLOTS (Step 2)
# =============================================================================
function plot_task_residual_scatters(df::DataFrame)
    println("\n" * "="^60)
    println("STEP 2: VISUALIZING TASK-RESIDUAL RELATIONSHIPS")
    println("="^60)
    
    task_cols = [:task_abstract, :task_routine, :task_manual]
    
    # Standardize task features for better plot scaling (optional but recommended)
    # We create a copy to avoid modifying original data
    plot_df = copy(df)
    for col in task_cols
        μ = mean(plot_df[!, col])
        σ = std(plot_df[!, col])
        plot_df[!, Symbol(col, "_std")] = σ > 0 ? (plot_df[!, col] .- μ) ./ σ : zeros(nrow(plot_df))
    end
    
    # Prepare subplots
    plots = []
    
    for (i, col) in enumerate(task_cols)
        std_col = Symbol(col, "_std")
        
        # Spearman trend line (non-parametric fit)
        # Sort by X for clean line plotting
        sorted_idx = sortperm(plot_df[!, std_col])
        x_vals = plot_df[!, std_col][sorted_idx]
        y_vals = plot_df.residual[sorted_idx]
        
        p = scatter(plot_df[!, std_col], plot_df.residual,
            xlabel = "$(replace(String(col), "_" => " ")) (Standardized)",
            ylabel = "Residual (Actual - Predicted)",
            title = "$(replace(String(col), "_" => " ")) vs. Residual",
            label = "", markersize = 8, color = :steelblue, 
            legend = :topleft, grid = true, size = (450, 350),
            framestyle = :box
        )
        
        # Add trend line
        plot!(p, x_vals, y_vals, color = :red, linewidth = 2, label = "Spearman Trend", linestyle = :solid)
        
        # Annotate each occupation
        for j in 1:nrow(plot_df)
            # Offset text slightly to avoid overlapping with points
            annotate!(p, plot_df[!, std_col][j], plot_df.residual[j],
                      text(plot_df.occ_name[j], 7, :left))
        end
        
        # Add zero line reference
        hline!(p, [0], color = :gray, linewidth = 1, linestyle = :dot, label = "")
        
        push!(plots, p)
    end
    
    # Combine plots into one figure
    combined = plot(plots..., layout = (1, 3), size = (1400, 400))
    
    # Save plot
    plot_path = joinpath(OUTPUT_DIR, "task_residual_scatters.png")
    savefig(combined, plot_path)
    println("✅ Saved scatter plots to: $plot_path")
    
    return combined
end

# =============================================================================
# 4. MAIN EXECUTION
# =============================================================================
function main()
    println("Starting Feature-Residual Analysis...")
    
    # Load data
    df = load_and_prepare_data()
    
    # Step 1: Spearman Correlations
    corr_results = compute_spearman_correlations(df)
    
    # Step 2: Visualization
    plot_task_residual_scatters(df)
    
    println("\n" * "="^60)
    println("ANALYSIS COMPLETE")
    println("="^60)
    println("Check the generated CSV and PNG files in: $OUTPUT_DIR")
    println("Interpretation guide:")
    println("  • ρ > 0.5 or ρ < -0.5 : Strong association (focus on these)")
    println("  • ρ near 0 : No clear monotonic relationship")
    println("  • Scatter plots show which specific occupations drive the correlation")
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

main()