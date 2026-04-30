# =============================================================================
# 06_variance_decomposition_table.jl
#
# Render variance_decomposition_r2.csv as a publication-style table image.
#
# Input:
#   result/occ_ind_4/variance_decomposition/variance_decomposition_r2.csv
#
# Output:
#   result/occ_ind_4/variance_decomposition/variance_decomposition_r2_table.png
#   result/occ_ind_4/variance_decomposition/variance_decomposition_r2_table.pdf
# =============================================================================

using CSV
using DataFrames
using CairoMakie
using Printf

const BASE_DIR = joinpath(@__DIR__, "..", "..", "result", "occ_ind_4")
const OUTPUT_DIR = joinpath(BASE_DIR, "variance_decomposition")
const INPUT_CSV = joinpath(OUTPUT_DIR, "variance_decomposition_r2.csv")

const OUTCOME_LABELS = Dict(
    "employment" => "Employment",
    "hourly_rate" => "Hourly rate",
    "hours" => "Hours",
    "income" => "Income",
    "income_share_var" => "Income share",
    "inequality" => "Inequality",
    "median" => "Median income",
    "unemployment" => "Unemployment",
)

const TABLE_COLUMNS = [
    :outcome,
    :H,
    :R2_industry_only,
    :AdjR2_industry_only,
    :R2_occupation_only,
    :AdjR2_occupation_only,
    :R2_industry_occupation,
    :AdjR2_industry_occupation,
]

const COLUMN_LABELS = Dict(
    :outcome => "Outcome",
    :H => "H",
    :R2_industry_only => "R2: ind.",
    :AdjR2_industry_only => "Adj. R2: ind.",
    :R2_occupation_only => "R2: occ.",
    :AdjR2_occupation_only => "Adj. R2: occ.",
    :R2_industry_occupation => "R2: both",
    :AdjR2_industry_occupation => "Adj. R2: both",
)

const COLUMN_WIDTHS = Dict(
    :outcome => 160.0,
    :H => 55.0,
    :R2_industry_only => 95.0,
    :AdjR2_industry_only => 105.0,
    :R2_occupation_only => 95.0,
    :AdjR2_occupation_only => 105.0,
    :R2_industry_occupation => 95.0,
    :AdjR2_industry_occupation => 105.0,
)

const R2_COLUMNS = [
    :R2_industry_only,
    :R2_occupation_only,
    :R2_industry_occupation,
]

const HIGHLIGHT_BY_R2 = Dict(
    :R2_industry_only => :R2_industry_only,
    :AdjR2_industry_only => :R2_industry_only,
    :R2_occupation_only => :R2_occupation_only,
    :AdjR2_occupation_only => :R2_occupation_only,
    :R2_industry_occupation => :R2_industry_occupation,
    :AdjR2_industry_occupation => :R2_industry_occupation,
)

function format_cell(x)
    if ismissing(x)
        return ""
    elseif x isa AbstractFloat
        return isnan(x) ? "NaN" : @sprintf("%.3f", x)
    elseif x isa Integer
        return string(x)
    else
        return string(x)
    end
end

function format_outcome(x)
    return get(OUTCOME_LABELS, String(x), replace(String(x), "_" => " "))
end

function build_display_table(df::DataFrame)::DataFrame
    table = select(df, TABLE_COLUMNS)
    table[!, :outcome] = format_outcome.(table.outcome)
    sort!(table, [:outcome, :H])
    return table
end

function draw_table_image(table::DataFrame, output_path::String;
                          title::String = "Variance decomposition R2")
    n_rows = nrow(table)
    row_h = 30.0
    header_h = 34.0
    title_h = 52.0
    foot_h = 34.0
    left_pad = 20.0
    right_pad = 20.0

    col_widths = [COLUMN_WIDTHS[c] for c in TABLE_COLUMNS]
    table_w = sum(col_widths)
    fig_w = Int(round(left_pad + table_w + right_pad))
    fig_h = Int(round(title_h + header_h + n_rows * row_h + foot_h))

    fig = Figure(size = (fig_w, fig_h), backgroundcolor = :white)
    ax = Axis(fig[1, 1])
    hidedecorations!(ax)
    hidespines!(ax)
    xlims!(ax, 0, fig_w)
    ylims!(ax, 0, fig_h)

    x0 = left_pad
    y_top = fig_h - title_h

    text!(ax, fig_w / 2, fig_h - 24;
          text = title, align = (:center, :center),
          fontsize = 20, font = :bold, color = :black)

    x_edges = [x0]
    for w in col_widths
        push!(x_edges, last(x_edges) + w)
    end

    header_bottom = y_top - header_h
    table_bottom = header_bottom - n_rows * row_h

    band_color = RGBf(0.965, 0.972, 0.982)
    header_color = RGBf(0.88, 0.91, 0.95)
    grid_color = RGBf(0.72, 0.76, 0.80)

    poly!(ax, Rect2f(x0, header_bottom, table_w, header_h);
          color = header_color, strokewidth = 0)

    for r in 1:n_rows
        y = header_bottom - r * row_h
        if isodd(r)
            poly!(ax, Rect2f(x0, y, table_w, row_h);
                  color = band_color, strokewidth = 0)
        end
    end

    for x in x_edges
        lines!(ax, [x, x], [table_bottom, y_top]; color = grid_color, linewidth = 1)
    end
    for y in vcat([y_top, header_bottom], [header_bottom - r * row_h for r in 1:n_rows])
        lines!(ax, [x0, x0 + table_w], [y, y]; color = grid_color, linewidth = 1)
    end

    for (j, col) in enumerate(TABLE_COLUMNS)
        xc = (x_edges[j] + x_edges[j + 1]) / 2
        text!(ax, xc, header_bottom + header_h / 2;
              text = COLUMN_LABELS[col], align = (:center, :center),
              fontsize = 11, font = :bold, color = :black)
    end

    max_r2_by_row = map(eachrow(table)) do row
        vals = [
            row[c] for c in R2_COLUMNS
            if !ismissing(row[c]) && row[c] isa Number && isfinite(row[c])
        ]
        isempty(vals) ? missing : maximum(vals)
    end

    for (i, row) in enumerate(eachrow(table))
        yc = header_bottom - (i - 0.5) * row_h
        for (j, col) in enumerate(TABLE_COLUMNS)
            xc = (x_edges[j] + x_edges[j + 1]) / 2
            value = col == :outcome ? row[col] : format_cell(row[col])
            align = col == :outcome ? (:left, :center) : (:center, :center)
            xpos = col == :outcome ? x_edges[j] + 8 : xc
            highlight_col = get(HIGHLIGHT_BY_R2, col, nothing)
            is_best_r2 = highlight_col !== nothing &&
                         !ismissing(max_r2_by_row[i]) &&
                         row[highlight_col] == max_r2_by_row[i]

            if is_best_r2
                cell_pad = 5.0
                poly!(ax,
                      Rect2f(x_edges[j] + cell_pad,
                             yc - row_h / 2 + cell_pad,
                             col_widths[j] - 2cell_pad,
                             row_h - 2cell_pad);
                      color = RGBf(1.0, 0.91, 0.58),
                      strokewidth = 0)
            end

            text!(ax, xpos, yc;
                  text = value, align = align,
                  fontsize = 10,
                  font = is_best_r2 ? :bold : :regular,
                  color = is_best_r2 ? RGBf(0.10, 0.10, 0.10) : :black)
        end
    end

    text!(ax, x0, 18;
          text = "Source: variance_decomposition_r2.csv. Values rounded to three decimals. Highlight = model with largest R2 within row.",
          align = (:left, :center), fontsize = 9, color = RGBf(0.25, 0.25, 0.25))

    save(output_path, fig, px_per_unit = 2)
    return output_path
end

function main(; input_csv::String = INPUT_CSV,
                output_dir::String = OUTPUT_DIR)
    isfile(input_csv) || error("Input CSV not found: $input_csv")
    mkpath(output_dir)

    df = CSV.read(input_csv, DataFrame)
    table = build_display_table(df)

    png_path = joinpath(output_dir, "variance_decomposition_r2_table.png")
    pdf_path = joinpath(output_dir, "variance_decomposition_r2_table.pdf")

    draw_table_image(table, png_path)
    draw_table_image(table, pdf_path)

    println("Saved table image:")
    println("  ", png_path)
    println("  ", pdf_path)

    return table
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
