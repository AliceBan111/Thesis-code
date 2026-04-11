using DataFrames
using StatFiles
using FileIO
using XLSX

# 1. 路径配置
input_dir  = "../../data"
output_dir = "../../result/mapping"
mkpath(output_dir)
cd(@__DIR__)


file_1990 = joinpath(input_dir, "occ1990_occ1990dd.dta")
file_2010 = joinpath(input_dir, "occ2010_occ1990dd.dta")
output_xlsx = joinpath(output_dir, "mapping.xlsx")

df_1990 = DataFrame(load(file_1990))
df_2010 = DataFrame(load(file_2010))

# 2. 避免列名冲突，提前重命名
rename!(df_1990, :occ => :occ1990)
rename!(df_2010, :occ => :occ2010)

# 3. 聚合 occ2010：将同一个 occ1990dd 对应的多个 occ2010 合并到一行
df_2010_agg = combine(groupby(df_2010, :occ1990dd),
    :occ2010 => (x -> begin
        # 去重 + 过滤缺失值 + 转为字符串
        vals = string.(unique(filter(!ismissing, x)))
        # 若无有效值则保留 missing，否则用 "; " 拼接
        isempty(vals) ? missing : join(vals, "; ")
    end) => :occ2010
)

# 4. 左连接
df_mapping = leftjoin(df_1990, df_2010_agg, on = :occ1990dd)

# 调整列顺序（仅保留核心三列）
select!(df_mapping, :occ1990, :occ1990dd, :occ2010)

# 5. 导出为 Excel

XLSX.writetable(output_xlsx, collect(eachcol(df_mapping)), names(df_mapping))

println("Saved file: $output_xlsx")
