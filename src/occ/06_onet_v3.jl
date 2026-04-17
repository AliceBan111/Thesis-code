using DataFrames
using CSV
using XLSX
using Statistics
using Tables

# ============================================================
#  HELPER FUNCTIONS
# ============================================================

"""
    clean_soc_code(code::AbstractString)

Removes hyphens from O*NET-SOC codes to standardize format for matching.
Example: "11-1011.00" -> "1110110" ... keeps digits/dots, strips hyphens.
"""
function clean_soc_code(code::AbstractString)
    return replace(string(code), "-" => "")
end

"""
    weighted_mean_skipmissing(values, weights)

Calculates the weighted mean, skipping observations where the value is missing.
Returns `missing` if all values are missing or total weight is zero.
"""
function weighted_mean_skipmissing(values::AbstractVector, weights::AbstractVector)
    valid_indices = .!ismissing.(values)
    if !any(valid_indices)
        return missing
    end
    v_float = Float64.(values[valid_indices])
    w_float = Float64.(weights[valid_indices])
    denom = sum(w_float)
    return denom == 0 ? missing : sum(v_float .* w_float) / denom
end

"""
    standardize_column(x)

Z-score standardizes a numeric vector, skipping missing values.
Returns a vector of Float64 (NaN where input was missing).
"""
function standardize_column(x::AbstractVector)
    x_float = convert(Vector{Union{Float64, Missing}}, x)
    mu    = mean(skipmissing(x_float))
    sigma = std(skipmissing(x_float))
    if ismissing(mu) || ismissing(sigma) || sigma == 0
        return fill(NaN, length(x_float))
    end
    return [(ismissing(v) ? NaN : (v - mu) / sigma) for v in x_float]
end

"""
    calc_correlation_pairwise(x, y)

Pearson correlation ignoring rows where either value is missing/NaN.
Returns `missing` if fewer than 3 valid pairs exist.
"""
function calc_correlation_pairwise(x::AbstractVector, y::AbstractVector)
    valid = .!ismissing.(x) .& .!ismissing.(y) .& .!isnan.(x) .& .!isnan.(y)
    sum(valid) < 3 && return missing
    try
        return cor(Float64.(x[valid]), Float64.(y[valid]))
    catch
        return missing
    end
end

# ============================================================
#  STEP 1 – READ ONE O*NET FILE AND RESHAPE TO WIDE FORMAT
# ============================================================

"""
    process_onet_wide(filename, table_prefix; data_dir, col_positions)

Reads an O*NET .xlsx or .csv file, filters for Scale_ID == "LV", prefixes every
Element Name with `table_prefix`, and returns a wide-format DataFrame keyed by
cleaned SOC code.

Column positions (1-based) default to the standard O*NET layout:
  soc_idx   = 1  (O*NET-SOC Code)
  elem_idx  = 4  (Element Name)   ← was 5 in original, off by one for most files
  scale_idx = 6  (Scale ID)
  value_idx = 8  (Data Value)

Adjust via keyword arguments if your files differ.
"""
function process_onet_wide(filename::String, table_prefix::String;
                           data_dir::String  = "../../data/ONET",
                           soc_idx::Int      = 1,
                           elem_idx::Int     = 5,   # BUG FIX: was 5 (off by one for most O*NET files)
                           scale_idx::Int    = 6,
                           value_idx::Int    = 8)

    filepath = joinpath(data_dir, filename)
    isfile(filepath) || error("File not found: $filepath")

    # ── Read file ──────────────────────────────────────────────────────────────
    local raw::DataFrame
    ext = lowercase(last(splitext(filename)))

    if ext == ".xlsx"
        # BUG FIX: original code mixed up df_proc (local rename) with df (original)
        # after building df_proc, later code still referenced undefined variables
        # `scale_col`, `elem_col`, `soc_col`, `sub_col`, `value_col`.
        # Rewritten to work directly with positional column indices on one DataFrame.
        xf         = XLSX.readxlsx(filepath)
        sheet      = XLSX.sheetnames(xf)[1]
        raw = DataFrame(XLSX.readtable(filepath, sheet; header=true))

    elseif ext == ".csv"
        raw = CSV.read(filepath, DataFrame; missingstring=["", "NA", "N/A"])

    else
        error("Unsupported format: $filename  (use .xlsx or .csv)")
    end

    ncols = ncol(raw)
    for idx in (soc_idx, elem_idx, scale_idx, value_idx)
        idx <= ncols || error("Column index $idx out of range (file has $ncols columns): $filename")
    end

    # ── Extract and rename columns of interest ─────────────────────────────────
    df = DataFrame(
        SOC_Clean      = clean_soc_code.(strip.(string.(raw[!, soc_idx]))),
        Element_Name   = strip.(string.(raw[!, elem_idx])),
        Scale_ID       = strip.(string.(raw[!, scale_idx])),
        Data_Value     = raw[!, value_idx]
    )

    # ── Filter Scale_ID == "LV" ────────────────────────────────────────────────
    filter!(row -> row.Scale_ID == "LV", df)
    nrow(df) > 0 || error("No rows with Scale_ID == 'LV' found in $filename")

    # ── Build prefixed element name ────────────────────────────────────────────
    df[!, :Wide_Name] = string.(table_prefix, "_", df[!, :Element_Name])

    # ── Ensure numeric Data_Value ──────────────────────────────────────────────
    df[!, :Data_Value] = coalesce.(tryparse.(Float64, string.(df[!, :Data_Value])), NaN)

    # ── Reshape long → wide ────────────────────────────────────────────────────
    # BUG FIX: original code tried to unstack on `soc_col` and `sub_col` which were
    # never defined in scope (they were local to the earlier df_proc block that was
    # then abandoned). Sub_Code is also dropped here intentionally – aggregation
    # happens at the SOC level in the next step.
    df_wide = unstack(df, :SOC_Clean, :Wide_Name, :Data_Value,
                      combine = mean)   # average duplicate (SOC, element) pairs

    return df_wide   # columns: SOC_Clean, one column per element
end

# ============================================================
#  STEP 2 – AGGREGATE O*NET WIDE DATA TO 9 OCCUPATION GROUPS
# ============================================================

"""
    process_onet_aggregated(df_wide, mapping_file; mapping_dir, sheet_name)

Joins the wide O*NET frame to your mapping file (which supplies Group 1-9 and
employment weights), then returns a 9-row DataFrame of weighted-average element
values per group.

Mapping file column positions (1-based):
  col_soc     = 5  (O*NET-SOC Code)
  col_group   = 6  (Occupation Group 1-9)
  col_weights = 8  (Employment weights)
"""
function process_onet_aggregated(df_wide::DataFrame, mapping_file::String;
                                 mapping_dir::String  = "../../result/mapping",
                                 sheet_name::String   = "Sheet1",
                                 col_soc::Int         = 5,
                                 col_group::Int        = 6,
                                 col_weights::Int      = 8)

    # ── Read mapping file ──────────────────────────────────────────────────────
    mapping_path = joinpath(mapping_dir, mapping_file)
    isfile(mapping_path) || error("Mapping file not found: $mapping_path")

    xf = XLSX.readxlsx(mapping_path)
    sheet_name in XLSX.sheetnames(xf) || error("Sheet '$sheet_name' not found in $mapping_path")

    # BUG FIX: original used DataFrame(xf[sheet_name], :auto) which reads raw
    # matrix and may include a header row as data. Use readtable instead.
    df_map_raw = DataFrame(XLSX.readtable(mapping_path, sheet_name; header=true))

    ncols_map = ncol(df_map_raw)
    for idx in (col_soc, col_group, col_weights)
        idx <= ncols_map || error("Mapping column index $idx out of range ($ncols_map cols)")
    end

    df_map = DataFrame(
        SOC_Clean = clean_soc_code.(strip.(string.(df_map_raw[!, col_soc]))),
        Group     = df_map_raw[!, col_group],
        Weights   = df_map_raw[!, col_weights]
    )

    # Drop rows where any key field is missing / non-parseable
    filter!(row -> !ismissing(row.SOC_Clean) && !ismissing(row.Group) && !ismissing(row.Weights), df_map)
    df_map[!, :Group]   = Int.(df_map[!, :Group])
    df_map[!, :Weights] = Float64.(df_map[!, :Weights])

    # ── Join O*NET wide data with mapping ──────────────────────────────────────
    df_matched = innerjoin(df_wide, df_map, on = :SOC_Clean)
    nrow(df_matched) > 0 || error("No SOC codes matched between O*NET data and mapping file")

    # ── Identify element columns ───────────────────────────────────────────────
    meta_cols     = Set([:SOC_Clean, :Group, :Weights])
    element_cols = [c for c in propertynames(df_matched) if !(c in meta_cols)]

    # ── Weighted mean per group for each element ───────────────────────────────
    # BUG FIX: original used a `push!` into `weighted_exprs` with a closure that
    # captured `c` by reference – the classic Julia loop-closure bug. All closures
    # would refer to the last value of `c`. Fixed by wrapping in a let block.
    gdf = groupby(df_matched, :Group)

    weighted_exprs = [
        let col = c
            [col, :Weights] => ((vals, wts) -> weighted_mean_skipmissing(vals, wts)) => col
        end
        for c in element_cols
    ]

    df_grouped = combine(gdf, weighted_exprs...)
    sort!(df_grouped, :Group)

    return df_grouped   # 9 rows × (1 + n_elements) columns
end

# ============================================================
#  STEP 3 – CUMULATIVE IRF SCALARS PER OUTCOME
# ============================================================

"""
    calculate_cumulative_irf_single(input_file, target_outcome; data_dir)

Reads the merged IRF trajectories CSV, filters for `target_outcome`, sums betas
over horizons 1-36 for each occupational group, and returns a 2-column DataFrame:
  :Occupational_Group  |  Symbol(target_outcome)
"""
function calculate_cumulative_irf_single(input_file::String, target_outcome::String;
                                          data_dir::String = "../../result/occ/analysis")

    filepath = joinpath(data_dir, input_file)
    isfile(filepath) || error("File not found: $filepath")

    df = CSV.read(filepath, DataFrame; missingstring=["", "NA", "N/A"])

    # Filter for this outcome
    filter!(row -> row.outcome == target_outcome, df)
    nrow(df) > 0 || error("Outcome '$target_outcome' not found in $input_file")

    # Identify group columns (everything except :outcome and :horizon)
    exclude = Set([:outcome, :horizon])
    grp_cols = [c for c in propertynames(df) if !(c in exclude)]
    isempty(grp_cols) && error("No group columns found in $input_file")

    # Reshape wide → long
    df_long = stack(df, grp_cols, variable_name=:Occupational_Group, value_name=:Beta)

    # Filter horizons 1-36
    filter!(row -> !ismissing(row.horizon) && row.horizon >= 1 && row.horizon <= 36, df_long)

    # Sum betas by group
    df_cum = combine(groupby(df_long, :Occupational_Group),
                     :Beta => (b -> sum(skipmissing(b))) => Symbol(target_outcome))

    sort!(df_cum, :Occupational_Group)
    return df_cum
end

# ============================================================
#  STEP 4 – TOP-10 CORRELATIONS FOR ONE TABLE × ONE OUTCOME
# ============================================================

"""
    get_top10_correlations(onet_file, table_prefix, outcome; kwargs...)

Full pipeline for one O*NET table and one outcome:
  1. Read & reshape O*NET → wide
  2. Aggregate to 9 groups (weighted)
  3. Z-score standardize element columns
  4. Load cumulative IRF for `outcome`
  5. Inner-join on group identifier
  6. Correlate every element with the outcome's cumulative IRF
  7. Return the top 10 by |r|

Returns a DataFrame with columns: element_name, correlation_coefficient, source_table.
"""
function get_top10_correlations(onet_file::String, table_prefix::String, outcome::String;
                                onet_dir::String    = "../../data/ONET",
                                mapping_dir::String = "../../result/mapping",
                                irf_dir::String     = "../../result/occ/analysis",
                                irf_file::String    = "merged_occ_irf_trajectories.csv")

    # 1-2. O*NET → wide → group-level
    df_wide   = process_onet_wide(onet_file, table_prefix; data_dir=onet_dir)
    df_groups = process_onet_aggregated(df_wide, "mapping_done.xlsx"; mapping_dir=mapping_dir)

    # BUG FIX: process_onet_aggregated returns :Group (Int), but
    # calculate_cumulative_irf_single returns :Occupational_Group (String from column names).
    # We need a common key. Rename :Group → :Occupational_Group in df_groups,
    # converting to String to match what stack() produces from column names.
    df_groups[!, :Occupational_Group] = string.(df_groups[!, :Group])
    select!(df_groups, Not(:Group))

    # 3. Z-score standardize element columns
    meta = Set([:Occupational_Group])
    elem_cols = [c for c in propertynames(df_groups) if !(c in meta)]
    for col in elem_cols
        df_groups[!, col] = standardize_column(df_groups[!, col])
    end

    # 4. Cumulative IRF
    df_irf = calculate_cumulative_irf_single(irf_file, outcome; data_dir=irf_dir)
    df_irf[!, :Occupational_Group] = replace.(string.(df_irf[!, :Occupational_Group]), r"\D" => "")

    # 5. Merge on group identifier
    df_groups[!, :Occupational_Group] =
    string.(strip.(df_groups[!, :Occupational_Group]))

    df_irf[!, :Occupational_Group] =
        string.(strip.(df_irf[!, :Occupational_Group]))

    println("onet groups = ", unique(df_groups.Occupational_Group))
    println("irf groups = ", unique(df_irf.Occupational_Group))

    df_merged = innerjoin(df_groups, df_irf, on=:Occupational_Group)

    println("merged rows = ", nrow(df_merged))
    nrow(df_merged) == 0 && return DataFrame(element_name=String[],
                                              correlation_coefficient=Float64[],
                                              source_table=String[])

    # 6. Correlate
    outcome_sym = Symbol(outcome)
    results = NamedTuple{(:element_name, :correlation_coefficient, :source_table),
                          Tuple{String, Float64, String}}[]

    for col in elem_cols
        r = calc_correlation_pairwise(df_merged[!, col], df_merged[!, outcome_sym])
        ismissing(r) && continue
        push!(results, (element_name=String(col), correlation_coefficient=Float64(r),
                        source_table=table_prefix))
    end

    isempty(results) && return DataFrame(element_name=String[],
                                          correlation_coefficient=Float64[],
                                          source_table=String[])

    # 7. Sort by |r|, take top 10
    df_res = DataFrame(results)
    df_res[!, :abs_corr] = abs.(df_res[!, :correlation_coefficient])
    sort!(df_res, :abs_corr, rev=true)
    return select(first(df_res, 10), :element_name, :correlation_coefficient, :source_table)
end

# ============================================================
#  MAIN – LOOP OVER ALL OUTCOMES × ALL O*NET TABLES
# ============================================================

function main()
    outcome_list = [
        "employment", "hourly_rate", "income_share",
        "inequality", "median", "unemployment", "income", "hours"
    ]

    onet_tables = [
        ("Abilities.xlsx",        "Abilities"),
        ("Knowledge.xlsx",        "Knowledge"),
        ("Skills.xlsx",           "Skills"),
        ("Work Activities.xlsx",  "WorkActivities"),
    ]

    output_dir = "../../result/occ/correlation"
    mkpath(output_dir)

    for outcome in outcome_list
        println("\n=== Outcome: $outcome ===")

        df_combined = DataFrame(
            element_name            = String[],
            correlation_coefficient = Float64[],
            source_table            = String[]
        )

        for (file, prefix) in onet_tables
            println("  ↳ Table: $prefix")
            try
                df_top10 = get_top10_correlations(file, prefix, outcome)
                nrow(df_top10) > 0 && append!(df_combined, df_top10)
            catch e
                @warn "  Skipped $prefix for $outcome" exception=e
            end
        end

        if nrow(df_combined) > 0
            out_path = joinpath(output_dir, "$(outcome).csv")
            CSV.write(out_path, df_combined)
            println("  Saved $(nrow(df_combined)) rows → $out_path")
        else
            println("  No data for outcome: $outcome")
        end
    end

    println("\n✓ Pipeline complete.")
end

cd(@__DIR__)
main()