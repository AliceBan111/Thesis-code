using DataFrames, XLSX, Statistics, Dates, LinearAlgebra, Plots

file_path = "data/KLEMS_latest.xlsx"  # ← 改成你的路径

# ============================================================
# 1. Load data
# ============================================================
markup_data = XLSX.readxlsx(file_path)

function get_long_df(markup_data, sheet_name, value_col_name)
    sh = markup_data[sheet_name]
    data = sh[:]
    header = [Symbol("Industry"); Symbol.(data[2, 2:end])]
    body = data[3:end, :]
    df_inner = DataFrame(body, header)
    filter!(row -> !ismissing(row.Industry), df_inner)
    df_long = stack(df_inner, Not(:Industry), variable_name=:year, value_name=value_col_name)
    df_long.year = parse.(Int, string.(df_long.year))
    df_long[!, value_col_name] = Float64.(df_long[!, value_col_name])
    return df_long
end

df_e  = get_long_df(markup_data, "Energy Compensation",    :energy)
df_m  = get_long_df(markup_data, "Materials Compensation", :material)
df_s  = get_long_df(markup_data, "Service Compensation",   :service)
df_go = get_long_df(markup_data, "Gross Output",           :go)
df_va = get_long_df(markup_data, "Value Added",            :va)

merged = innerjoin(df_e, df_m,  on=[:Industry, :year])
merged = innerjoin(merged, df_s,  on=[:Industry, :year])
merged = innerjoin(merged, df_go, on=[:Industry, :year])
merged = innerjoin(merged, df_va, on=[:Industry, :year])

# ============================================================
# 2. Sanity check
# ============================================================
merged.implied_go = merged.va .+ merged.energy .+ merged.material .+ merged.service
merged.go_ratio   = merged.go ./ merged.implied_go
println("\n===== Sanity Check: GO / (VA + Intermediates) =====")
println("均值: ",   round(mean(skipmissing(merged.go_ratio)), digits=4))
println("中位数: ", round(median(skipmissing(merged.go_ratio)), digits=4))
println("应接近1.0")

# ============================================================
# 3. Calculate markup
# ============================================================
merged.markup_all = merged.go ./ (merged.energy .+ merged.material .+ merged.service)
merged.markup_me  = merged.go ./ (merged.energy .+ merged.material)
merged.s_m        = (merged.energy .+ merged.material .+ merged.service) ./ merged.go

# 固定权重
va_avg        = combine(groupby(merged, :Industry), :va => mean => :va_mean)
va_avg.weight = va_avg.va_mean ./ sum(va_avg.va_mean)
merged        = leftjoin(merged, va_avg[:, [:Industry, :weight]], on=:Industry)

# 加权聚合
function weighted_agg(merged, markup_col)
    result = combine(groupby(merged, :year)) do sdf
        valid = .!isnan.(sdf[!, markup_col]) .& 
                .!isinf.(sdf[!, markup_col]) .& 
                .!ismissing.(sdf[!, markup_col])
        return (annual_markup = sum(sdf[valid, markup_col] .* sdf.weight[valid]),)
    end
    sort!(result, :year)
    return result
end

agg_all   = weighted_agg(merged, :markup_all)
agg_me    = weighted_agg(merged, :markup_me)
agg_share = combine(groupby(merged, :year)) do sdf
    valid = .!isnan.(sdf.s_m) .& .!isinf.(sdf.s_m)
    return (share = sum(sdf.s_m[valid] .* sdf.weight[valid]) * 100,)
end
sort!(agg_share, :year)

# ============================================================
# 4. 统计摘要
# ============================================================
println("\n===== Markup序列统计 =====")
println("--- All intermediates ---")
println("均值:   ", round(mean(agg_all.annual_markup), digits=4))
println("标准差: ", round(std(agg_all.annual_markup),  digits=4))
println("极差:   ", round(maximum(agg_all.annual_markup) - minimum(agg_all.annual_markup), digits=4))
println("--- Materials + Energy only ---")
println("均值:   ", round(mean(agg_me.annual_markup), digits=4))
println("标准差: ", round(std(agg_me.annual_markup),  digits=4))
println("极差:   ", round(maximum(agg_me.annual_markup) - minimum(agg_me.annual_markup), digits=4))
println("\n===== Intermediate Share =====")
println("均值: ",   round(mean(agg_share.share), digits=2), "%  (论文参考: 44-48%)")

# ============================================================
# 5. 画图
# ============================================================
p1 = plot(agg_all.year, agg_all.annual_markup,
          label="All intermediates",
          color=:blue, linewidth=2,
          xlabel="Year", ylabel="Log Markup",
          title="Aggregate Price Markup")
plot!(p1, agg_me.year, agg_me.annual_markup,
          label="Materials + Energy only",
          color=:red, linewidth=2, linestyle=:dash)

p2 = plot(agg_share.year, agg_share.share,
          label="Intermediate share (%)",
          color=:green, linewidth=2,
          xlabel="Year", ylabel="Share (%)",
          title="Aggregate Intermediate Share (论文Figure 4参考: 44-48%)")
hline!(p2, [44.0, 48.0],
           label="论文参考范围",
           color=:gray, linestyle=:dot)

display(plot(p1, p2, layout=(2,1), size=(900, 700)))


p1 = plot(agg_all.year, agg_all.annual_markup,
          label="All intermediates",
          color=:blue, linewidth=2,
          xlabel="Year", ylabel="Log Markup",
          title="Aggregate Price Markup (All Intermediates)")

p2 = plot(agg_me.year, agg_me.annual_markup,
          label="Materials + Energy only",
          color=:red, linewidth=2, linestyle=:dash,
          xlabel="Year", ylabel="Log Markup",
          title="Aggregate Price Markup (Materials + Energy)")

display(plot(p1, p2, layout=(2,1), size=(900, 700)))