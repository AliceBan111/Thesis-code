using DataFrames, XLSX, Statistics, Interpolations, Dates

"""
    build_markup_method2(df::DataFrame)
    non HP-filtered markup (Bils et al. style)
    markup = log (gross output / (energy compensation + materials compensation + service compensation))

Adds two columns to df:
- markup_level: interpolated to quarterly
- markup_growth: first difference of markup

Requires: data/KLEMS.xlsx
"""
function build_markup_method2_no_hp_latest(df::DataFrame, file_path::String)
    df = copy(df)
    
    # 1. load data
    markup_data = XLSX.readxlsx(file_path)
    
    function get_long_df(sheet_name, value_col_name)
        sh = markup_data[sheet_name]
        data = sh[:]

        header = [Symbol("Industry"); Symbol.(data[2, 2:end])]
        body = data[3:end, :]
        df_inner = DataFrame(body, header)

        filter!(row -> !ismissing(row.Industry), df_inner)

        df_long = stack(df_inner, Not(:Industry), variable_name=:year, value_name=value_col_name)
        df_long.year = parse.(Int, string.(df_long.year))
        df_long[!, value_col_name] = Float64.(df_long[!, value_col_name])
        
        return df_long
    end

    df_e  = get_long_df("Energy Compensation", :energy)
    df_m  = get_long_df("Materials Compensation", :material)
    df_s  = get_long_df("Service Compensation", :service)
    df_go = get_long_df("Gross Output", :go)
    df_va = get_long_df("Value Added", :va)

    merged = innerjoin(df_e, df_m, on=[:Industry, :year])
    merged = innerjoin(merged, df_s, on=[:Industry, :year])
    merged = innerjoin(merged, df_go, on=[:Industry, :year])
    merged = innerjoin(merged, df_va, on=[:Industry, :year])

    # 2. calculate markup
    merged.markup_ind = log.(merged.go ./ (merged.energy .+ merged.material .+ merged.service))

    agg_markup = combine(groupby(merged, :year)) do sdf
        valid = .!isnan.(sdf.markup_ind) .& .!isinf.(sdf.markup_ind)
        
        sub_va = sdf.va[valid]
        sub_mu = sdf.markup_ind[valid]

        weighted_val = sum(sub_mu .* (sub_va ./ sum(sub_va)))
        
        return (annual_markup = weighted_val,) 
    end

    sort!(agg_markup, :year)

    # 3. interpolate to quarterly
    values = agg_markup.annual_markup
    itp = interpolate(values, BSpline(Cubic(Line(OnGrid()))))
    
    n_quarters = length(values) * 4
    quarterly_indices = range(1, stop=length(values), length=n_quarters)
    interpolated_values = [itp(x) for x in quarterly_indices]

    quarterly_markup_pool = DataFrame(
        "observation_date" => Date[],
        "markup_level" => Float64[]
    )

    current_idx = 1
    for yr in agg_markup.year
        y = Int(yr)
        for m in [1, 4, 7, 10]
            if current_idx <= length(interpolated_values)
                push!(quarterly_markup_pool, (
                    observation_date = Date(y, m, 1),
                    markup_level = interpolated_values[current_idx]
                ))
                current_idx += 1
            end
        end
    end

    # 4. join to input DataFrame
    df_result = leftjoin(df, quarterly_markup_pool, on="observation_date")
    sort!(df_result, "observation_date")
    
    # 5. compute quarterly growth
    markup_vals = df_result.markup_level
    growth = Vector{Union{Missing, Float64}}(missing, nrow(df_result))
    
    for i in 2:nrow(df_result)
        if !ismissing(markup_vals[i]) && !ismissing(markup_vals[i-1])
            growth[i] = markup_vals[i] - markup_vals[i-1]
        end
    end
    
    df_result.markup_growth = growth

    return df_result
end
