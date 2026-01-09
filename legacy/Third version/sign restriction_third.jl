# inverse of labor share (Nekarda & Ramey (2020))
cd(@__DIR__)

include("ToolBox/load_var.jl")

using CSV, DataFrames, Dates, LinearAlgebra, Statistics, Random, Distributions, Plots, SpecialFunctions, XLSX, HypothesisTests

start_date = Date("1964-01-01")
end_date = Date("2017-12-31")

# 1. Create dataframe

# log_real gdp
gdp_data = CSV.read("data/GDPC1.csv", DataFrame)
main_df = gdp_data[(gdp_data.observation_date.>=start_date).&(gdp_data.observation_date.<=end_date), [:observation_date]]
main_df.ln_gdp = log.(gdp_data[gdp_data.observation_date.∈[Set(main_df.observation_date)], :GDPC1])
main_df.ln_gdp_diff = [missing; diff(main_df.ln_gdp)]

# log difference inflation data
inflation_data = CSV.read("data/GDPDEF.csv", DataFrame)
main_df.ln_inflation = log.(inflation_data[inflation_data.observation_date.∈[Set(main_df.observation_date)], :GDPDEF])
main_df.pi_p = [missing; diff(main_df.ln_inflation)]

# Unemployment rate
u_data = CSV.read("data/UNRATE.csv", DataFrame)
u_filtered = u_data[(u_data.observation_date.>=start_date).&(u_data.observation_date.<=end_date), :]
main_df.u = u_filtered.UNRATE

#markup
markup_data = DataFrame(XLSX.readtable("data/markup.xlsx", "Sheet1", infer_eltypes=true))

labor_row = markup_data[markup_data.Measure.=="Labor compensation", :]
va_row = markup_data[markup_data.Measure.=="Value-added output", :]

cols = names(markup_data)

time_cols = cols[5:end]

selected_cols = filter(c ->
        begin
            year = parse(Int, split(c, " ")[1])
            (year ≥ 1964) && (year ≤ 2017)
        end,
    time_cols
)

labor = collect(labor_row[1, selected_cols])
value_added = collect(va_row[1, selected_cols])

markup = -log.(labor ./ value_added)

main_df.markup = markup
main_df.markup_growth = [missing; diff(markup)]

# 10-year government bond
iL_data = CSV.read("data/GS10.csv", DataFrame)
iL_filtered = iL_data[(iL_data.observation_date.>=start_date).&(iL_data.observation_date.<=end_date), :]
main_df.iL = iL_filtered.GS10

final_df = main_df[2:end, :]

svar_data = final_df[:, [:ln_gdp_diff, :pi_p, :u, :markup_growth, :iL]]

# 2. Estimate reduced-form VAR
X = Matrix(svar_data)
X_float = convert(Matrix{Float64}, coalesce.(X, NaN))
X_clean, fo, lo = CommonSample(X_float)

max_lags = 8
optimal_lag, criteria = VARlag(X_clean, max_lags, 1)
println("Optimal lag: ", optimal_lag)

VAR, VARopt = VARmodel(X_clean, optimal_lag, 1)  # const=1

# Display estimation results
println("VAR coefficient matrices:")
display(VAR[:F])
println("Residual covariance matrix:")
display(VAR[:sigma])

# Model diagnostics
println("Maximum eigenvalue: ", VAR[:maxEig])
if VAR[:maxEig] < 1
    println("VAR system is stable")
else
    println("VAR system may be unstable")
end

# Check R-squared for each equation
for j in 1:size(X_clean, 2)
    println("Equation $j R²: ", VAR[Symbol("eq$j")][:rsqr])
end

# 3. Sign Restriction

# Define sign restriction matrix

SIGN = [
    -1 0 0 0 0;   # ln_gdp_diff  ↓
    1 0 0 0 0;   # pi_p         ↑
    1 0 0 0 0;   # u            
    1 0 0 0 0;   # markup_growth↑
    0 0 0 0 0    # iL           unrestricted
]


# Set sign restriction parameters using dictionary syntax
VARopt[:nsteps] = 20
VARopt[:ndraws] = 1000
VARopt[:sr_hor] = 4
VARopt[:pctg] = 68
VARopt[:sr_draw] = 100000
VARopt[:sr_rot] = 1000
VARopt[:sr_mod] = 0
VARopt[:mult] = 100

# Force the covariance matrix to be symmetric (acceptable due to numerical precision)
VAR[:sigma] = 0.5 * (VAR[:sigma] + VAR[:sigma]')

# Re-check the matrix properties
println("Is covariance matrix symmetric? ", issymmetric(VAR[:sigma]))
println("Is covariance matrix positive definite? ", isposdef(VAR[:sigma]))


# Run sign restriction identification
SRout = SR(VAR, SIGN, VARopt)

# Display results
println("Number of accepted rotations: ", size(SRout[:Ball], 3))

# Extract median impulse responses for price markup shock
SRout[:IRmed] = SRout[:IRmed]
SRout[:IRinf] = SRout[:IRinf]
SRout[:IRsup] = SRout[:IRsup]

# Plot impulse responses for price markup shock (first shock)
variable_names = ["GDP Growth (ln_gdp_diff)", 
                  "Price Inflation (pi_p)", 
                  "Unemployment Rate (u)", 
                  "Markup Growth (markup_growth)", 
                  "Long-term Interest Rate (iL)"]

# Plot all variables' responses to price markup shock
nsteps_actual = size(SRout[:IRmed], 1)

plt = plot(layout=(3, 2), size=(800, 600))
for i in 1:5
    plot!(plt[i], 1:nsteps_actual, SRout[:IRinf][:, i, 1],
        fillrange=SRout[:IRsup][:, i, 1],
        fillalpha=0.3,
        fillcolor=:lightblue,
        linealpha=0,
        label="68% CI")
    plot!(plt[i], 1:nsteps_actual, SRout[:IRmed][:, i, 1],
        label="Median",
        linewidth=2,
        color=:blue)
    plot!(plt[i], [0, nsteps_actual + 1], [0, 0],
        linestyle=:dash,
        color=:black,
        label="")
    title!(plt[i], "$(variable_names[i])")
    xlabel!(plt[i], "Horizon")
end
display(plt)


for i in 1:5
    println("$(variable_names[i]): ", round(SRout[:IRmed][1, i, 1], digits=4))
end

# Check if sign restrictions are satisfied
# Construct price markup level by integrating markup_growth
price_markup_response = cumsum(SRout[:IRmed][:, 4, 1])  # markup_growth is 4th variable

println("Price inflation impact response: ", round(SRout[:IRmed][1, 2, 1], digits=4))  # pi_p is 2nd variable
println("Price markup impact response: ", round(price_markup_response[1], digits=4))

if sign(SRout[:IRmed][1, 2, 1]) == sign(price_markup_response[1])
    println("Positive comovement constraint satisfied")
    println("Both variables move in the same direction")
else
    println("Positive comovement constraint NOT satisfied")
    println("Variables move in opposite directions")
end
