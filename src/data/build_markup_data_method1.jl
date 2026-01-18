using DataFrames, XLSX, Dates, Interpolations, ReadStatTables

"""
    build_markup_method1(df::DataFrame)
    inverse of labor share (Nekarda & Ramey (2020))

Adds two columns to df:
- markup_level: -log(Labor / Value-added)
- markup_growth: first difference of markup

Requires: data/markup.xlsx (Sheet1)
"""
function build_markup_method1(df::DataFrame)
    # 1. obtain df time range
    start_date = minimum(df.observation_date)
    end_date = maximum(df.observation_date)

    # 2. load data
    markup_data = DataFrame(XLSX.readtable(
        joinpath(@__DIR__, "../../data/markup.xlsx"), "Sheet1"; 
        infer_eltypes = true
    ))
    
    labor_row = markup_data[markup_data.Measure .== "Labor compensation", :]
    va_row = markup_data[markup_data.Measure .== "Value-added output", :]

    # 3. data conversion
    function parse_q_date(s)
        parts = split(s, " ")
        y = parse(Int, parts[1])
        q = parts[2]
        m = q == "Q1" ? 1 : q == "Q2" ? 4 : q == "Q3" ? 7 : 10
        return Date(y, m, 1)
    end

    # 4. filter columns according to the time range
    all_cols = names(markup_data)[5:end]
   
    selected_cols = filter(c -> begin
            dt = parse_q_date(c)
            return dt >= (start_date - Month(3)) && dt <= end_date
        end, all_cols)

    # 5. markup calculation
    labor = Float64.(collect(Vector(labor_row[1, selected_cols])))
    value_added = Float64.(collect(Vector(va_row[1, selected_cols])))
    
    markup_vec = -log.(labor ./ value_added)

    # 6. construct temp
    temp_markup_df = DataFrame(
        observation_date = parse_q_date.(selected_cols),
        markup_level = markup_vec,
        markup_growth = [missing; diff(markup_vec)]
    )

    # 7. using leftjoin to align the date
    res_df = leftjoin(df, temp_markup_df, on = :observation_date)

    return res_df
end