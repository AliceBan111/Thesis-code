using CSV, DataFrames, Dates

"""
    build_macro_data(start_date::Date, end_date::Date)

Returns a DataFrame with:
- ln_gdp_diff: GDP growth (log difference)
- pi_p: inflation (log difference)
- u: unemployment rate
- du: unemployment change
- iL: 10-year government bond
"""
function build_macro_data(start_date::Date, end_date::Date)
    df = DataFrame(observation_date = Date[])

    # 1. GDP
    gdp_data = CSV.read("data/GDPC1.csv", DataFrame)
    df.observation_date = gdp_data[(gdp_data.observation_date .>= start_date) .& 
                                   (gdp_data.observation_date .<= end_date), :observation_date]
    df.ln_gdp = log.(gdp_data[gdp_data.observation_date .∈ Set(df.observation_date), :GDPC1])
    df.ln_gdp_diff = [missing; diff(df.ln_gdp)]

    # 2. Inflation
    inflation_data = CSV.read("data/GDPDEF.csv", DataFrame)
    df.pi_p = [missing; diff(log.(inflation_data[inflation_data.observation_date .∈ Set(df.observation_date), :GDPDEF]))]

    # 3. Unemployment
    u_data = CSV.read("data/UNRATE.csv", DataFrame)
    df.u = u_data[(u_data.observation_date .>= start_date) .& 
                  (u_data.observation_date .<= end_date), :UNRATE]
    df.du = [missing; diff(df.u)]

    # 4. 10-year government bond
    iL_data = CSV.read("data/GS10.csv", DataFrame)
    df.iL = iL_data[(iL_data.observation_date .>= start_date) .&
                    (iL_data.observation_date .<= end_date), :GS10]


    return df[2:end, :]
end
