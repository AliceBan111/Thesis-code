using DataFrames, XLSX, Dates, Interpolations, ReadStatTables

"""
    build_markup_method1(df::DataFrame)
    inverse of labor share (Nekarda & Ramey (2020))

Adds two columns to df:
- markup: -log(Labor / Value-added)
- markup_growth: first difference of markup

Requires: data/markup.xlsx (Sheet1)
"""
function build_markup_method1(df::DataFrame)
    # 1. load data
    markup_data = DataFrame(XLSX.readtable("data/markup.xlsx", "Sheet1", infer_eltypes=true))

    labor_row = markup_data[markup_data.Measure .== "Labor compensation", :]
    va_row = markup_data[markup_data.Measure .== "Value-added output", :]

    # 2. find out year column
    cols = names(markup_data)
    time_cols = cols[5:end]   

    selected_cols = filter(c -> begin
            year = parse(Int, split(c, " ")[1])
            (year ≥ 1964) && (year ≤ 2017)
        end,
        time_cols
    )

    # 3. markup
    labor = collect(labor_row[1, selected_cols])
    value_added = collect(va_row[1, selected_cols])

    markup = -log.(labor ./ value_added)
    markup_growth = [missing; diff(markup)]

    # 4.add to df
    df = copy(df)   
    df.markup = markup
    df.markup_growth = markup_growth

    return df
end
