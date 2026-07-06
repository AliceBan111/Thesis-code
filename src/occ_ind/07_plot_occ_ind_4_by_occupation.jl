using CairoMakie
using CSV
using DataFrames
using Printf

const RESULT_ROOT = joinpath(@__DIR__, "..", "..", "result", "occ_ind_4")
const OUT_ROOT = joinpath(RESULT_ROOT, "by_occupation")

const EXCLUDED_RESULT_DIRS = Set([
    "candidate",
    "task_regression",
    "task_regression_total",
    "variance_decomposition",
    "by_occupation",
])

const OCC_LABELS_BY_FILE = Dict(
    "irf_occ1_Managerial.csv" => (1, "Managerial"),
    "irf_occ2_Professional_specialty.csv" => (2, "Professional specialty"),
    "irf_occ3_High_tech.csv" => (3, "High tech"),
    "irf_occ4_Sales.csv" => (4, "Sales"),
    "irf_occ5_Administrative_support.csv" => (5, "Administrative support"),
    "irf_occ6_Service.csv" => (6, "Service"),
    "irf_occ7_Farming_forestry_construction.csv" => (7, "Farming/forestry/construction"),
    "irf_occ8_Precision_production_repair.csv" => (8, "Precision production/repair"),
    "irf_occ9_Machine_operators_transport.csv" => (9, "Machine operators/transport"),
)

const INDUSTRY_ORDER = [
    "Energy_intensive",
    "Manufacturing_Construction",
    "Trade",
    "Services",
]

const INDUSTRY_DISPLAY = Dict(
    "Energy_intensive" => "Energy intensive",
    "Manufacturing_Construction" => "Mfg. + construction",
    "Trade" => "Trade",
    "Services" => "Services",
)

const INDUSTRY_COLORS = Dict(
    "Energy_intensive" => :darkorange,
    "Manufacturing_Construction" => :steelblue,
    "Trade" => :seagreen,
    "Services" => :mediumpurple,
)

const OUTCOME_LABELS = Dict(
    "hourly_rate" => "Log hourly wage",
    "hours" => "Log hours worked",
    "income" => "Log weekly income",
    "income_share_var" => "Income share",
    "inequality" => "Log income ratio (P75/P25)",
    "median" => "Median log income",
    "unemployment" => "Unemployment rate",
    "employment" => "Log employment",
)

finite_values(x) = Float64[Float64(v) for v in skipmissing(x) if isfinite(Float64(v))]

function outcome_dirs()
    dirs = filter(readdir(RESULT_ROOT; join=false)) do name
        path = joinpath(RESULT_ROOT, name)
        isdir(path) && !(name in EXCLUDED_RESULT_DIRS)
    end
    return sort(dirs)
end

function read_outcome_irfs(outcome::String)::DataFrame
    rows = DataFrame()

    for industry in INDUSTRY_ORDER
        industry_dir = joinpath(RESULT_ROOT, outcome, industry)
        isdir(industry_dir) || continue

        for file in sort(readdir(industry_dir; join=false))
            haskey(OCC_LABELS_BY_FILE, file) || continue
            occ_id, occ_label = OCC_LABELS_BY_FILE[file]
            df = CSV.read(joinpath(industry_dir, file), DataFrame)
            df.outcome .= outcome
            df.industry .= industry
            df.industry_label .= get(INDUSTRY_DISPLAY, industry, industry)
            df.occ_group .= occ_id
            df.occupation .= occ_label
            append!(rows, df; cols=:union)
        end
    end

    return rows
end

function global_ylims(df::DataFrame)
    vals = Float64[]
    for col in (:beta, :ci_lo90, :ci_hi90)
        col in propertynames(df) || continue
        append!(vals, finite_values(df[!, col]))
    end

    isempty(vals) && return (-0.01, 0.01)
    y_min = minimum(vals)
    y_max = maximum(vals)
    span = y_max - y_min
    pad = span == 0 ? max(abs(y_min), 0.01) * 0.2 : span * 0.08
    return (y_min - pad, y_max + pad)
end

function plot_outcome_by_occupation(df::DataFrame, outcome::String; shared_y::Bool=false)
    out_dir = joinpath(OUT_ROOT, outcome)
    mkpath(out_dir)

    sort!(df, [:occ_group, :industry, :horizon])

    title = get(OUTCOME_LABELS, outcome, outcome)
    fig = Figure(size=(1450, 1050), fontsize=13)
    scale_note = shared_y ? "shared y-axis" : "occupation-specific y-axis"
    Label(fig[0, 1:3], "IRFs by occupation across industries: $(title) ($(scale_note))";
          fontsize=20, font=:bold, tellwidth=false)

    shared_ylims = global_ylims(df)

    for occ in 1:9
        row = ceil(Int, occ / 3)
        col = mod1(occ, 3)
        occ_df = filter(:occ_group => ==(occ), df)
        occ_label = isempty(occ_df.occupation) ? "Occupation $(occ)" : first(occ_df.occupation)

        ax = Axis(fig[row, col];
            title="($(occ)) $(occ_label)",
            xlabel="Horizon (months)",
            ylabel=col == 1 ? "Response of $(title)" : "",
            xticks=0:12:36,
        )
        panel_ylims = shared_y ? shared_ylims : global_ylims(occ_df)
        ylims!(ax, panel_ylims[1], panel_ylims[2])
        hlines!(ax, [0.0]; color=:gray35, linewidth=1, linestyle=:dash)

        for industry in INDUSTRY_ORDER
            sub = filter(:industry => ==(industry), occ_df)
            nrow(sub) == 0 && continue
            sort!(sub, :horizon)
            line_df = dropmissing(sub, [:horizon, :beta])
            nrow(line_df) == 0 && continue
            keep = isfinite.(Float64.(line_df.horizon)) .& isfinite.(Float64.(line_df.beta))
            line_df = line_df[keep, :]
            nrow(line_df) == 0 && continue

            lines!(ax, Float64.(line_df.horizon), Float64.(line_df.beta);
                   color=get(INDUSTRY_COLORS, industry, :black),
                   linewidth=2.6,
                   label=get(INDUSTRY_DISPLAY, industry, industry))
        end
    end

    labels = [get(INDUSTRY_DISPLAY, ind, ind) for ind in INDUSTRY_ORDER]
    elements = [LineElement(color=get(INDUSTRY_COLORS, ind, :black), linewidth=3)
                for ind in INDUSTRY_ORDER]
    Legend(fig[4, 1:3], elements, labels; orientation=:horizontal, framevisible=false,
           tellheight=true, tellwidth=false)

    stem = shared_y ? "irf_by_occupation_shared_y" : "irf_by_occupation"
    for ext in ("png", "pdf")
        save(joinpath(out_dir, "$(stem).$ext"), fig)
    end
    println("Saved $(outcome): ", joinpath(out_dir, "$(stem).pdf"))
end

function main()
    mkpath(OUT_ROOT)
    combined = DataFrame()

    for outcome in outcome_dirs()
        df = read_outcome_irfs(outcome)
        if nrow(df) == 0
            @warn "No IRF CSVs found for outcome=$outcome"
            continue
        end
        mkpath(joinpath(OUT_ROOT, outcome))
        CSV.write(joinpath(OUT_ROOT, outcome, "irf_by_occupation_long.csv"), df)
        plot_outcome_by_occupation(df, outcome; shared_y=false)
        plot_outcome_by_occupation(df, outcome; shared_y=true)
        append!(combined, df; cols=:union)
    end

    if nrow(combined) > 0
        CSV.write(joinpath(OUT_ROOT, "occ_ind_4_irf_by_occupation_long.csv"), combined)
        println("Saved combined long file: ", joinpath(OUT_ROOT, "occ_ind_4_irf_by_occupation_long.csv"))
    end
end

main()
