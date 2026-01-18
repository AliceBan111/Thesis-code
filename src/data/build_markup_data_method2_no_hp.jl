using DataFrames, XLSX, Dates, Interpolations, ReadStatTables

"""
    build_markup_method2(df::DataFrame)
    no HP-filtered markup (Bils et al. style)

Adds two columns to df:
- markup_level: interpolated to quarterly
- markup_growth: first difference of markup

Requires: data/Hklems_combined_updatedMarkAug2015.dta
"""
function build_markup_method2_no_hp(df::DataFrame)
    # 1. load data
    markup_data = DataFrame(readstat("data/klems_combined_updatedMarkAug2015.dta"))

    # 2. aggregate by year
    markup_data.s_m = markup_data.ninterm ./ markup_data.noutput

    agg_markup = combine(
        groupby(dropmissing(markup_data, :s_m), :year),
        :s_m => mean => :s_m_avg
    )
    sort!(agg_markup, :year)

    years = agg_markup.year
    values = -log.(agg_markup.s_m_avg)

    # 3. interpolate to quarterly
    itp = interpolate(values, BSpline(Cubic(Line(OnGrid()))))
    quarterly_indices = range(1, stop=length(values), length=length(values)*4)
    interpolated_values = [itp(x) for x in quarterly_indices]

    quarterly_markup_pool = DataFrame(
        observation_date = Date[],
        markup_level = Float64[]
    )

    current_idx = 1
    for row in eachrow(agg_markup)
        y = Int(row.year)
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
    df_result = leftjoin(df, quarterly_markup_pool, on=:observation_date)

    # 5. compute quarterly growth
    df_result.markup_growth = [missing; diff(df_result.markup_level)]

    return df_result
end
