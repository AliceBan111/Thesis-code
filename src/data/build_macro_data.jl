cd(@__DIR__)
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
    # 1. Create a DataFrame with observation dates
    gdp_data = CSV.read(
    joinpath(@__DIR__, "../../data/GDPC1.csv"),
    DataFrame
)
    df = gdp_data[(gdp_data.observation_date .>= start_date) .& 
                  (gdp_data.observation_date .<= end_date), [:observation_date]]

    # 2. GDP
    gdp_filtered = gdp_data[(gdp_data.observation_date .>= start_date) .& 
                            (gdp_data.observation_date .<= end_date), :]
    df.ln_gdp = log.(gdp_filtered.GDPC1)
    df.ln_gdp_diff = [missing; diff(df.ln_gdp)]

    # 3. Inflation
    inflation_data = CSV.read(
    joinpath(@__DIR__, "../../data/GDPDEF.csv"),
    DataFrame
)
    infl_filtered = inflation_data[(inflation_data.observation_date .>= start_date) .& 
                                   (inflation_data.observation_date .<= end_date), :]
    df.pi_p = [missing; diff(log.(infl_filtered.GDPDEF))]

    # 4. Unemployment
    u_data = CSV.read(
    joinpath(@__DIR__, "../../data/UNRATE.csv"),
    DataFrame
)
    u_filtered = u_data[(u_data.observation_date .>= start_date) .& 
                        (u_data.observation_date .<= end_date), :]
    df.u = u_filtered.UNRATE
    df.du = [missing; diff(df.u)]

    # 5. 10-year government bond
    iL_data = CSV.read(
    joinpath(@__DIR__, "../../data/GS10.csv"),
    DataFrame
)
    iL_filtered = iL_data[(iL_data.observation_date .>= start_date) .&
                          (iL_data.observation_date .<= end_date), :]
    df.iL = iL_filtered.GS10

    # Remove first row due to diff → missing
    return df[2:end, :]
end
