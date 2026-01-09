using DataFrames, XLSX, Dates, Interpolations, ReadStatTables

"""
    build_markup_method2(df::DataFrame)
    HP-filtered markup (Bils et al. style)

Adds two columns to df:
- markup: hp_l_MWint87 interpolated to quarterly
- markup_growth: first difference of markup

Requires: data/HPfiltered_AllIndustries_energy.dta
"""
function build_markup_method2(df::DataFrame)
    # 1. load data
    markup_data = DataFrame(readstat("data/HPfiltered_AllIndustries_energy.dta"))

    # 2. aggregate by year
    agg_markup = combine(
        groupby(dropmissing(markup_data, :hp_l_MWint87), :year),
        :hp_l_MWint87 => mean => :hp_l_MWint87_agg
    )

    years = agg_markup.year
    values = agg_markup.hp_l_MWint87_agg

    # 3. interpolate to quarterly
    itp = interpolate(values, BSpline(Cubic(Line(OnGrid()))))
    quarterly_indices = range(1, stop=length(values), length=length(values)*4)
    interpolated_values = [itp(x) for x in quarterly_indices]

    quarterly_markup = DataFrame(
        observation_date = Date[],
        hp_l_MWint87 = Float64[]
    )

    current_idx = 1
    for row in eachrow(agg_markup)
        y = Int(row.year)
        for m in [1, 4, 7, 10]
            if current_idx <= length(interpolated_values)
                push!(quarterly_markup, (
                    observation_date = Date(y, m, 1),
                    hp_l_MWint87 = interpolated_values[current_idx]
                ))
                current_idx += 1
            end
        end
    end

    sort!(quarterly_markup, :observation_date)

    # 4. join to input DataFrame
    df = copy(df)
    df = leftjoin(df, quarterly_markup, on=:observation_date)

    # 5. compute quarterly growth
    df.markup = df.hp_l_MWint87
    df.markup_growth = [missing; diff(df.markup)]

    return df
end
