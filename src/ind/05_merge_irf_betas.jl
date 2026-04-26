using CSV
using DataFrames

const RESULT_DIR = normpath(joinpath(@__DIR__, "..", "..", "result", "ind_v2"))
const OUTPUT_FILE = joinpath(RESULT_DIR, "merged_irf_betas.csv")
const HORIZONS = 0:36
const GROUPS = 1:12

function list_outcomes(result_dir::AbstractString)
    names = readdir(result_dir)
    dirs = filter(name -> isdir(joinpath(result_dir, name)), names)
    sort!(dirs)
    return dirs
end

function load_group_beta_map(outcome_dir::AbstractString)
    beta_map = Dict{Tuple{Int, Int}, Union{Missing, Float64}}()

    for file in readdir(outcome_dir)
        match_obj = match(r"^irf_group(\d+)_.*\.csv$", file)
        isnothing(match_obj) && continue

        group = parse(Int, match_obj.captures[1])
        group in GROUPS || continue

        df = CSV.read(joinpath(outcome_dir, file), DataFrame)
        for row in eachrow(df)
            horizon = Int(row.horizon)
            horizon in HORIZONS || continue
            beta_map[(group, horizon)] = ismissing(row.beta) ? missing : Float64(row.beta)
        end
    end

    return beta_map
end

function build_merged_table(result_dir::AbstractString)
    rows = NamedTuple[]

    for outcome in list_outcomes(result_dir)
        outcome_dir = joinpath(result_dir, outcome)
        beta_map = load_group_beta_map(outcome_dir)

        for horizon in HORIZONS
            row = (
                outcome = outcome,
                horizon = horizon,
                Symbol("1") => get(beta_map, (1, horizon), missing),
                Symbol("2") => get(beta_map, (2, horizon), missing),
                Symbol("3") => get(beta_map, (3, horizon), missing),
                Symbol("4") => get(beta_map, (4, horizon), missing),
                Symbol("5") => get(beta_map, (5, horizon), missing),
                Symbol("6") => get(beta_map, (6, horizon), missing),
                Symbol("7") => get(beta_map, (7, horizon), missing),
                Symbol("8") => get(beta_map, (8, horizon), missing),
                Symbol("9") => get(beta_map, (9, horizon), missing),
                Symbol("10") => get(beta_map, (10, horizon), missing),
                Symbol("11") => get(beta_map, (11, horizon), missing),
                Symbol("12") => get(beta_map, (12, horizon), missing),
            )
            push!(rows, row)
        end
    end

    return DataFrame(rows)
end

function main()
    merged = build_merged_table(RESULT_DIR)
    CSV.write(OUTPUT_FILE, merged)
    println("Merged file saved to: $(OUTPUT_FILE)")
end

main()
