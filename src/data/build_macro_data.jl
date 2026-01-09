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
    lookback_date = start_date - Month(3)

    # 1. Create a DataFrame with observation dates
    gdp_data = CSV.read(
    joinpath(@__DIR__, "../../data/GDPC1.csv"),
    DataFrame
)
    # 2. GDP
    gdp_filtered = gdp_data[(gdp_data.observation_date .>= lookback_date) .& 
                            (gdp_data.observation_date .<= end_date), :]
    full_ln_gdp = log.(gdp_filtered.GDPC1)
    full_ln_gdp_diff = [missing; diff(full_ln_gdp)]

    temp_gdp = DataFrame(
        observation_date = gdp_filtered.observation_date,
        ln_gdp = full_ln_gdp,
        ln_gdp_diff = full_ln_gdp_diff
    )

    df = temp_gdp[temp_gdp.observation_date .>= start_date, :]

    # 3. Inflation
    inflation_data = CSV.read(
    joinpath(@__DIR__, "../../data/GDPDEF.csv"),
    DataFrame
)
    infl_filtered = inflation_data[(inflation_data.observation_date .>= lookback_date) .& 
                                   (inflation_data.observation_date .<= end_date), :]
    
    infl_log_diff = [missing; diff(log.(infl_filtered.GDPDEF))]
    temp_infl = DataFrame(
        observation_date = infl_filtered.observation_date,
        pi_p = infl_log_diff
    )

    df = leftjoin(df, temp_infl, on = :observation_date)

    # 4. Unemployment
    u_data = CSV.read(
    joinpath(@__DIR__, "../../data/UNRATE.csv"),
    DataFrame
)
    u_filtered = u_data[(u_data.observation_date .>= lookback_date) .& 
                        (u_data.observation_date .<= end_date), :]
    
    temp_u = DataFrame(
        observation_date = u_filtered.observation_date,
        u = u_filtered.UNRATE,
        du = [missing; diff(u_filtered.UNRATE)]
    )
    df = leftjoin(df, temp_u, on = :observation_date)

    # 5. 10-year government bond
    iL_data = CSV.read(
    joinpath(@__DIR__, "../../data/GS10.csv"),
    DataFrame
)
    iL_filtered = iL_data[(iL_data.observation_date .>= start_date) .&
                          (iL_data.observation_date .<= end_date), :]
    
    temp_iL = iL_filtered[:, [:observation_date, :GS10]]
    rename!(temp_iL, :GS10 => :iL)
    
    df = leftjoin(df, temp_iL, on = :observation_date)

    return df
end
