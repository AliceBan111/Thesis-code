using CSV
using DataFrames
using Statistics
using Clustering
using CairoMakie
using Printf

const BASE_DIR = normpath(joinpath(@__DIR__, "..", "..", "result", "occ"))
const OUTPUT_DIR = joinpath(BASE_DIR, "analysis")

const IRF_FILES = [
    ("employment", joinpath(BASE_DIR, "wide_employment.csv")),
    ("hourly_rate", joinpath(BASE_DIR, "wide_hourly_rate.csv")),
    ("hours", joinpath(BASE_DIR, "wide_hours.csv")),
    ("income", joinpath(BASE_DIR, "wide_income.csv")),
    ("income_share", joinpath(BASE_DIR, "wide_income_share.csv")),
    ("inequality", joinpath(BASE_DIR, "wide_inequality.csv")),
    ("median", joinpath(BASE_DIR, "wide_median.csv")),
    ("unemployment", joinpath(BASE_DIR, "wide_unemployment.csv")),
]

const OCCUPATION_LABELS = Dict(
    :group1_Managerial => "Managerial",
    :group2_Professional_specialty => "Professional specialty",
    :group3_High_tech => "High tech",
    :group4_Sales => "Sales",
    :group5_Administrative_support => "Administrative support",
    :group6_Service => "Service",
    :group7_Farming_forestry_construction => "Farming / forestry / construction",
    :group8_Precision_production_repair => "Precision production / repair",
    :group9_Machine_operators_transport => "Machine operators / transport",
)

"""
    load_irf_table(path, outcome; horizon_range=nothing)

Read one wide IRF CSV, optionally filter horizons, and append the outcome name.
"""
function load_irf_table(path::String, outcome::String; horizon_range=nothing)
    df = CSV.read(path, DataFrame)
    sort!(df, :horizon)

    if horizon_range !== nothing
        df = filter(:horizon => in(horizon_range), df)
    end

    insertcols!(df, 1, :outcome => fill(outcome, nrow(df)))
    return df
end

"""
    build_merged_irf_df(; horizon_range=nothing)

Stack the outcome-specific wide IRF tables into one long DataFrame with columns:
`outcome`, `horizon`, and the 9 occupation series.
"""
function build_standardized_merged_irf_df(; horizon_range=nothing)
    dfs = DataFrame[]
    occupation_cols = Symbol[]

    for (outcome, path) in IRF_FILES
        isfile(path) || error("IRF table not found: $path")

        # load one outcome table
        df = load_irf_table(path, outcome; horizon_range=horizon_range)

        # identify occupation columns
        current_occ_cols = names(df, Not([:outcome, :horizon]))

        if isempty(occupation_cols)
            occupation_cols = Symbol.(current_occ_cols)
        elseif Symbol.(current_occ_cols) != occupation_cols
            error("Occupation columns in $path do not match previous IRF tables.")
        end

        # --------------------------------------------------
        # Z-score standardization within this outcome only
        # each occupation column separately
        # --------------------------------------------------
        # for col in occupation_cols
        #     x = Float64.(df[!, col])

        #     μ = mean(x)
        #     σ = std(x)

        #     # avoid divide-by-zero if constant series
        #     if σ < 1e-12
        #         df[!, col] .= 0.0
        #     else
        #         df[!, col] = (x .- μ) ./ σ
        #     end
        # end

        push!(dfs, df)
    end

    merged_df = vcat(dfs..., cols=:setequal)
    sort!(merged_df, [:outcome, :horizon])

    return merged_df, occupation_cols
end

"""
    occupation_matrix(df, occupation_cols)

Convert the 9 occupation columns into a dense Float64 matrix.
Rows are outcome-horizon combinations; columns are occupations.
"""
function occupation_matrix(df::DataFrame, occupation_cols::Vector{Symbol})
    mat = Matrix(select(df, occupation_cols))
    any(ismissing, mat) && error("The merged IRF DataFrame contains missing values.")
    return Float64.(mat)
end

occupation_labels(occupation_cols::Vector{Symbol}) = [
    get(OCCUPATION_LABELS, col, replace(String(col), "_" => " ")) for col in occupation_cols
]

"""
    correlation_matrix_df(corr_mat, occupation_cols)

Turn the correlation matrix into a DataFrame for CSV export.
"""
function correlation_matrix_df(corr_mat::AbstractMatrix, occupation_cols::Vector{Symbol})
    out = DataFrame(occupation = String.(occupation_cols))
    for (j, col) in enumerate(occupation_cols)
        out[!, col] = corr_mat[:, j]
    end
    return out
end

"""
    plot_correlation_heatmap(corr_mat, occupation_cols, output_stem)

Plot the occupation correlation matrix as a labeled heatmap.
"""
function plot_correlation_heatmap(
    corr_mat::AbstractMatrix,
    occupation_cols::Vector{Symbol},
    output_stem::String,
)
    labels = occupation_labels(occupation_cols)
    n = length(labels)

    fig = Figure(size=(1150, 1000))
    ax = Axis(
        fig[1, 1],
        title="Correlation of Occupation IRF Trajectories",
        xlabel="Occupation",
        ylabel="Occupation",
        xticks=(1:n, labels),
        yticks=(1:n, labels),
        xticklabelrotation=pi / 4,
        yreversed=true,
        aspect=DataAspect(),
    )

    hm = heatmap!(
        ax,
        1:n,
        1:n,
        corr_mat;
        colormap=:RdBu,
        colorrange=(-1, 1),
    )

    for i in 1:n, j in 1:n
        val = corr_mat[i, j]
        text_color = abs(val) >= 0.60 ? :white : :black
        text!(
            ax,
            j,
            i;
            text=@sprintf("%.2f", val),
            align=(:center, :center),
            color=text_color,
            fontsize=12,
        )
    end

    Colorbar(fig[1, 2], hm, label="Correlation")

    save(output_stem * ".png", fig)
    save(output_stem * ".pdf", fig)
    return fig
end

"""
    correlation_distance(corr_mat)

Convert a correlation matrix to a valid distance matrix for hierarchical clustering.
"""
function correlation_distance(corr_mat::AbstractMatrix)
    dist_mat = clamp.(1 .- corr_mat, 0.0, 2.0)
    for i in axes(dist_mat, 1)
        dist_mat[i, i] = 0.0
    end
    return dist_mat
end

"""
    reference_coordinates(ref, leaf_x, cluster_x, heights)

Return the x/y location associated with a dendrogram node reference.
Negative refs denote leaves; positive refs denote previous merges.
"""
function reference_coordinates(
    ref::Int,
    leaf_x::Vector{Float64},
    cluster_x::Vector{Float64},
    heights::AbstractVector,
)
    if ref < 0
        leaf_idx = -ref
        return leaf_x[leaf_idx], 0.0
    end

    return cluster_x[ref], heights[ref]
end

"""
    plot_dendrogram(hc, occupation_cols, output_stem)

Draw a dendrogram from a `Clustering.Hclust` object using CairoMakie.
"""
function plot_dendrogram(hc::Hclust, occupation_cols::Vector{Symbol}, output_stem::String)
    labels = occupation_labels(occupation_cols)
    n = length(labels)

    leaf_x = zeros(Float64, n)
    for (position, leaf_idx) in enumerate(hc.order)
        leaf_x[leaf_idx] = position
    end

    ordered_labels = labels[hc.order]
    cluster_x = zeros(Float64, n - 1)

    fig = Figure(size=(1150, 750))
    ax = Axis(
        fig[1, 1],
        title="Hierarchical Clustering of Occupation IRF Trajectories",
        xlabel="Occupation",
        ylabel="Correlation distance (1 - correlation)",
        xticks=(1:n, ordered_labels),
        xticklabelrotation=pi / 4,
    )

    for merge_idx in 1:(n - 1)
        left_ref = hc.merges[merge_idx, 1]
        right_ref = hc.merges[merge_idx, 2]

        x_left, y_left = reference_coordinates(left_ref, leaf_x, cluster_x, hc.heights)
        x_right, y_right = reference_coordinates(right_ref, leaf_x, cluster_x, hc.heights)
        merge_height = hc.heights[merge_idx]

        lines!(ax, [x_left, x_left], [y_left, merge_height], color=:black, linewidth=2)
        lines!(ax, [x_right, x_right], [y_right, merge_height], color=:black, linewidth=2)
        lines!(ax, [x_left, x_right], [merge_height, merge_height], color=:black, linewidth=2)

        cluster_x[merge_idx] = (x_left + x_right) / 2
    end

    xlims!(ax, 0.5, n + 0.5)
    ylims!(ax, 0.0, maximum(hc.heights) * 1.05)

    save(output_stem * ".png", fig)
    save(output_stem * ".pdf", fig)
    return fig
end

"""
    clustering_order_df(hc, occupation_cols)

Save the dendrogram leaf order for quick inspection.
"""
function clustering_order_df(hc::Hclust, occupation_cols::Vector{Symbol})
    labels = occupation_labels(occupation_cols)
    return DataFrame(
        order=1:length(hc.order),
        occupation_raw=String.(occupation_cols[hc.order]),
        occupation_label=labels[hc.order],
    )
end

# function analyze_extremes_and_turning_points(df::DataFrame, occupation_cols::Vector{Symbol}, output_path::String)
#     results = DataFrame(
#         outcome = String[],
#         occupation = String[],
#         max_abs_shock = Float64[],
#         max_shock_horizon = Int[],
#         turning_point_horizons = String[]
#     )

#     for sub_df in groupby(df, :outcome)
#         outcome_name = first(sub_df.outcome)
#         horizons = sub_df.horizon

#         for col in occupation_cols
#             y = sub_df[!, col]
            
#             # 1. 极值定位：最大绝对值对应的 Horizon 及其大小
#             max_abs_idx = argmax(abs.(y))
#             max_abs_shock = y[max_abs_idx]
#             max_shock_horizon = horizons[max_abs_idx]

#             # 2. 拐点定位：寻找响应曲线中发生反转的 Horizon
#             # 通过判断前后一阶差分 (斜率) 是否改变符号来识别局部极值点
#             tps = Int[]
#             for i in 2:(length(y)-1)
#                 dy_prev = y[i] - y[i-1]
#                 dy_next = y[i+1] - y[i]
                
#                 # 若斜率发生正负反转且数值有明显波动，标记为拐点
#                 if dy_prev * dy_next < 0
#                     push!(tps, horizons[i])
#                 end
#             end
            
#             # 将该职业的信息记录入表
#             push!(results, (
#                 outcome_name,
#                 get(OCCUPATION_LABELS, col, String(col)),
#                 max_abs_shock,
#                 max_shock_horizon,
#                 isempty(tps) ? "None" : join(tps, ", ")
#             ))
#         end
#     end

#     CSV.write(output_path, results)
#     println("Saved Extremes & Turning Points to: $output_path")
#     return results
# end

"""
    main(; horizon_range=nothing, linkage=:average)

Run the full occupation IRF correlation and clustering analysis.
"""
function main(; horizon_range=nothing, linkage::Symbol=:average)
    mkpath(OUTPUT_DIR)

    merged_df, occupation_cols = build_standardized_merged_irf_df(; horizon_range=horizon_range)
    merged_path = joinpath(OUTPUT_DIR, "merged_occ_irf_trajectories.csv")
    CSV.write(merged_path, merged_df)

    traj_mat = occupation_matrix(merged_df, occupation_cols)
    corr_mat = cor(traj_mat)
    corr_path = joinpath(OUTPUT_DIR, "occupation_irf_correlation_matrix.csv")
    CSV.write(corr_path, correlation_matrix_df(corr_mat, occupation_cols))

    corr_heatmap_stem = joinpath(OUTPUT_DIR, "occupation_irf_correlation_heatmap")
    plot_correlation_heatmap(corr_mat, occupation_cols, corr_heatmap_stem)

    dist_mat = correlation_distance(corr_mat)
    hc = hclust(dist_mat; linkage=linkage)

    dendrogram_stem = joinpath(OUTPUT_DIR, "occupation_irf_dendrogram")
    plot_dendrogram(hc, occupation_cols, dendrogram_stem)

    cluster_order_path = joinpath(OUTPUT_DIR, "occupation_irf_cluster_order.csv")
    CSV.write(cluster_order_path, clustering_order_df(hc, occupation_cols))

    # extremes_path = joinpath(OUTPUT_DIR, "occupation_irf_extremes_turning_points.csv")
    # analyze_extremes_and_turning_points(merged_df, occupation_cols, extremes_path)

    println("Saved merged IRF table to: $merged_path")
    println("Saved correlation matrix to: $corr_path")
    println("Saved correlation heatmap to: $(corr_heatmap_stem).png/.pdf")
    println("Saved dendrogram to: $(dendrogram_stem).png/.pdf")
    println("Saved clustering order to: $cluster_order_path")

    return (
        merged_df=merged_df,
        corr_mat=corr_mat,
        distance_mat=dist_mat,
        hclust=hc,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

main()