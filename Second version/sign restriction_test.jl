cd(@__DIR__)

include("ToolBox/load_var.jl")

using CSV, DataFrames, Dates, LinearAlgebra, Statistics, Random, Distributions, Plots, SpecialFunctions

start_date = Date("1964-01-01")
end_date = Date("2017-12-31")

# 1. Create dataframe

# log_real gdp
gdp_data = CSV.read("data/GDPC1.csv", DataFrame)
main_df = gdp_data[(gdp_data.observation_date.>=start_date).&(gdp_data.observation_date.<=end_date), [:observation_date]]
main_df.ln_gdp = log.(gdp_data[gdp_data.observation_date.∈[Set(main_df.observation_date)], :GDPC1])

# Log_hours of all persons in the nonfarm business sector
hours_data = CSV.read("data/HOABS.csv", DataFrame)
hours_filtered = hours_data[(hours_data.observation_date.>=start_date).&(hours_data.observation_date.<=end_date), :]
main_df.ln_hours = log.(hours_filtered.HOABS)

# log difference wage data
wage_data = CSV.read("data/AHETPI.csv", DataFrame)
main_df.ln_wage = log.(wage_data[wage_data.observation_date.∈[Set(main_df.observation_date)], :AHETPI])
main_df.pi_w = [missing; diff(main_df.ln_wage)]

# log difference inflation data
inflation_data = CSV.read("data/GDPDEF.csv", DataFrame)
main_df.ln_inflation = log.(inflation_data[inflation_data.observation_date.∈[Set(main_df.observation_date)], :GDPDEF])
main_df.pi_p = [missing; diff(main_df.ln_inflation)]

# Unemployment rate
u_data = CSV.read("data/UNRATE.csv", DataFrame)
u_filtered = u_data[(u_data.observation_date.>=start_date).&(u_data.observation_date.<=end_date), :]
main_df.u = u_filtered.UNRATE

# 10-year government bond
iL_data = CSV.read("data/GS10.csv", DataFrame)
iL_filtered = iL_data[(iL_data.observation_date.>=start_date).&(iL_data.observation_date.<=end_date), :]
main_df.iL = iL_filtered.GS10

main_df.d_lprod = [missing; diff(main_df.ln_gdp - main_df.ln_hours)]

final_df = main_df[2:end, :]

svar_data = final_df[:, [:d_lprod, :pi_w, :pi_p, :u, :iL]]

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

SIGN = [0  0  0  0  0;   # d_lprod: labor productivity growth (unrestricted)
        0  0  0  0  0;   # pi_w: wage inflation (unrestricted)
        1  0  0  0  0;   # pi_p: price inflation (positive for price markup shock)
        0  0  0  0  0;   # u: unemployment rate (unrestricted)
        0  0  0  0  0]   # iL: long-term interest rate (unrestricted)

# Set sign restriction parameters using dictionary syntax
VARopt[:nsteps] = 20        
VARopt[:ndraws] = 1000      
VARopt[:sr_hor] = 20         
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
variable_names = ["Labor Productivity Growth", "Wage Inflation", "Price Inflation", 
                  "Unemployment Rate", "Long-term Interest Rate"]

# Plot all variables' responses to price markup shock
nsteps_actual = size(SRout[:IRmed], 1)

plt = plot(layout = (3, 2), size = (800, 600))
for i in 1:5
    plot!(plt[i], 1:nsteps_actual, SRout[:IRinf][:, i, 1], 
          fillrange = SRout[:IRsup][:, i, 1], 
          fillalpha = 0.3, 
          fillcolor = :lightblue, 
          linealpha = 0,  
          label = "68% CI")
    plot!(plt[i], 1:nsteps_actual, SRout[:IRmed][:, i, 1], 
          label = "Median", 
          linewidth = 2, 
          color = :blue)
    plot!(plt[i], [0, nsteps_actual+1], [0, 0], 
          linestyle = :dash, 
          color = :black, 
          label = "")
    title!(plt[i], "$(variable_names[i])")
    xlabel!(plt[i], "Horizon")
end
display(plt)


# Display impact responses
for i in 1:5
    println("$(variable_names[i]): ", round(SRout[:IRmed][1, i, 1], digits=4))
end

# Check if sign restrictions are satisfied
price_markup_response = construct_price_markup_irf(SRout[:IRmed], 1, size(SRout[:IRmed], 1))

println("Price inflation impact response: ", round(SRout[:IRmed][1, 3, 1], digits=4))
println("Price markup impact response: ", round(price_markup_response[1], digits=4))

if sign(SRout[:IRmed][1, 3, 1]) == sign(price_markup_response[1])
    println("Positive comovement constraint satisfied")
    println("Both variables move in the same direction")
else
    println("Positive comovement constraint NOT satisfied")
    println("Variables move in opposite directions")
end

