include("C:/Me/thesis/Thesis code/Forth version/sign restriction_forth.jl")
cd(@__DIR__)
using CSV, DataFrames, Statistics, Plots, Statistics, Dates

# extract price markup shock series

resid_matrix  = VAR[:resid]
println(size(resid_matrix ))

n_accepted = size(SRout[:Ball], 3)
median_idx = div(n_accepted, 2)
B_median = SRout[:Ball][:, :, median_idx]

println(size(B_median))

# structual shock = residual * B^(-1)'
structural_shocks = resid_matrix  * inv(B_median)'

println(size(structural_shocks))

# extract markup shock
markup_shocks = structural_shocks[:, 1]

#statistics of markup shock
println("  Average: $(round(mean(markup_shocks), digits=6))")
println("  Standard Deviation: $(round(std(markup_shocks), digits=4))")
println("  Maximum: $(round(maximum(markup_shocks), digits=4))")
println("  Minimum: $(round(minimum(markup_shocks), digits=4))")

# date align

p = optimal_lag
shock_dates = final_df.observation_date[(p+1):end]

shock_df = DataFrame(
    date = shock_dates,
    markup_shock = markup_shocks
)

shock_df[!, :quarter] = Dates.year.(shock_df.date) .+ (Dates.quarterofyear.(shock_df.date) .- 1) ./ 4

println("\nFirst 5 rows:")
println(first(shock_df, 5))
println("\nLast 5 rows:")
println(last(shock_df, 5))

# wage

df_working = DataFrame(
    YEAR = Int[],
    MONTH = Int[],
    EARNWT = Int[],
    EMPSTAT = Int[],
    OCC = Int[],
    EARNWEEK = Int[]
)

open("cps_00006.dat", "r") do io
    for line in eachline(io)
        empstat = parse(Int, strip(line[7:8]))
        if empstat == 10 || empstat == 12
            occ = parse(Int, strip(line[9:11]))
            earn = parse(Int, strip(line[22:29]))
            earnwt = parse(Int, strip(line[12:21]))
            if occ != 0 && earn != 999999 && earnwt > 0
                push!(df_working, (
                    parse(Int, strip(line[1:4])),
                    parse(Int, strip(line[5:6])),
                    earnwt,   
                    empstat,
                    occ,
                    earn
                ))
            end
        end
    end
end

df_working.EARNWEEK .= df_working.EARNWEEK ./ 100
df_working.EARNWT .= df_working.EARNWT ./ 10000

function white_blue_occ1990(occ::Int)
    if occ in vcat(473:498, 503:549, 558:599, 614:617, 628:699, 703:799, 803:889)
        return 0 # blue collar
    elseif occ in vcat(3:37, 43:200, 203:235, 243:283, 303:389, 405:469)
        return 1 # white collar
    else
        return missing
    end
end

df_working[!, :collar] = [white_blue_occ1990(occ) for occ in df_working.OCC]
df_working = dropmissing(df_working, :collar)

df_working[!, :quarter] = df_working.YEAR .+ ((df_working.MONTH .- 1) .÷ 3) ./ 4

df_q = combine(
    groupby(df_working, [:quarter, :collar]),
    [:EARNWEEK, :EARNWT] => ((w, wt) -> sum(w .* wt) / sum(wt)) => :avg_wage
)

df_n = df_q[df_q.collar .== 1, :]  # white collar
df_y = df_q[df_q.collar .== 0, :]  # blue collar

df_n[:, :dlog_wage] = [missing; diff(log.(df_n.avg_wage))]
df_y[:, :dlog_wage] = [missing; diff(log.(df_y.avg_wage))]


# combine shocks and wage

df_n_shock = leftjoin(df_n, shock_df[:, [:quarter, :markup_shock]], on=:quarter)
df_y_shock = leftjoin(df_y, shock_df[:, [:quarter, :markup_shock]], on=:quarter)


df_plot = innerjoin(
    df_n_shock[:, [:quarter, :dlog_wage, :markup_shock]],
    df_y_shock[:, [:quarter, :dlog_wage]],
    on=:quarter,
    makeunique=true
)
rename!(df_plot, :dlog_wage => :white_wage, :dlog_wage_1 => :blue_wage)
df_plot = dropmissing(df_plot)

# plot
plt1 = plot(size=(1200, 600), layout=(1,1))

plot!(plt1, df_plot.quarter, df_plot.white_wage, 
    label="White-collar wage growth", color=:blue, linewidth=2, alpha=0.8)
plot!(plt1, df_plot.quarter, df_plot.blue_wage, 
    label="Blue-collar wage growth", color=:green, linewidth=2, alpha=0.8)
ylabel!(plt1, "Quarterly log wage change", guidefontsize=10)

plt1_twin = twinx(plt1)
plot!(plt1_twin, df_plot.quarter, df_plot.markup_shock, 
    label="Historical markup shock", color=:red, linewidth=2.5, linestyle=:dash, alpha=0.7)
hline!(plt1_twin, [0], color=:gray, linestyle=:dot, linewidth=1, label="")
ylabel!(plt1_twin, "Markup shock size", guidefontsize=10)

xlabel!(plt1, "Quarter", guidefontsize=10)
title!(plt1, "Wage Growth and Historical Price Markup Shocks", titlefontsize=12)

display(plt1)

# correlation analysis

cor_white = cor(df_plot.white_wage, df_plot.markup_shock)
cor_blue = cor(df_plot.blue_wage, df_plot.markup_shock)

# lag correlation analysis
println("Lag | White-collar | Blue-collar")
println("----|--------------|------------")
for lag in 0:8
    if lag < nrow(df_plot)
        cor_white_lag = cor(df_plot.markup_shock[1:end-lag], df_plot.white_wage[lag+1:end])
        cor_blue_lag = cor(df_plot.markup_shock[1:end-lag], df_plot.blue_wage[lag+1:end])
        println("  $lag |    $(round(cor_white_lag, digits=4))    |   $(round(cor_blue_lag, digits=4))")
    end
end

# scatter plot
plt2 = scatter(df_plot.markup_shock, df_plot.white_wage,
    label="White-collar", color=:blue, alpha=0.6, markersize=4,
    xlabel="Markup shock", ylabel="Wage growth",
    title="Wage Growth vs Markup Shocks", size=(1000, 500))
    
scatter!(plt2[1], df_plot.markup_shock, df_plot.blue_wage,
    label="Blue-collar", color=:green, alpha=0.6, markersize=4)

display(plt2)

# wage IRF
max_lag = 4
n_lags = max_lag + 1

white_response = zeros(n_lags)
blue_response = zeros(n_lags)

# wage correlation
for lag in 0:max_lag
    if lag < nrow(df_plot)
        valid_idx = 1:(nrow(df_plot)-lag)
        white_response[lag+1] = cor(df_plot.markup_shock[valid_idx], df_plot.white_wage[lag .+ valid_idx])
        blue_response[lag+1] = cor(df_plot.markup_shock[valid_idx], df_plot.blue_wage[lag .+ valid_idx])
    end
end

plt_lag = plot(size=(1000, 500), title="Wage Response to Price Markup Shock", xlabel="Quarters after shock", ylabel="Correlation")

plot!(0:max_lag, white_response,
    label="White-collar", color=:blue, lw=3, marker=:circle)
plot!(0:max_lag, blue_response,
    label="Blue-collar", color=:green, lw=3, marker=:diamond)
hline!([0], linestyle=:dash, color=:black, label="")

display(plt_lag)

println("Lag | White-collar | Blue-collar")
println("----|--------------|------------")
for lag in 0:max_lag
    println("  $lag |    $(round(white_response[lag+1], digits=4))    |   $(round(blue_response[lag+1], digits=4))")
end