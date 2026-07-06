using CSV
using DataFrames

const RESULT_DIR = normpath(joinpath(@__DIR__, "..", "..", "result", "ind_4"))
const OUTPUT_FILE = joinpath(RESULT_DIR, "merged_irf_betas.csv")
const HORIZONS = 0:36
const GROUPS = 1:4

function list_outcomes(result_dir::AbstractString)
    !isdir(result_dir) && return String[]
    dirs = filter(name -> isdir(joinpath(result_dir, name)), readdir(result_dir))
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
            push!(rows, (
                outcome = outcome,
                horizon = horizon,
                Symbol("1") => get(beta_map, (1, horizon), missing),
                Symbol("2") => get(beta_map, (2, horizon), missing),
                Symbol("3") => get(beta_map, (3, horizon), missing),
                Symbol("4") => get(beta_map, (4, horizon), missing),
            ))
        end
    end

    return DataFrame(rows)
end

function main()
    mkpath(RESULT_DIR)
    merged = build_merged_table(RESULT_DIR)
    if nrow(merged) == 0
        println("No ind_4 outcome folders found under $(RESULT_DIR); merged file not written.")
        return merged
    end
    CSV.write(OUTPUT_FILE, merged)
    println("Merged file saved to: $(OUTPUT_FILE)")
    return merged
end

main()
