using DataFrames, XLSX, Statistics, Interpolations, Dates, LinearAlgebra

"""
    build_markup_method2(df::DataFrame)
    non HP-filtered markup (Bils et al. style)
    markup = log (gross output / (energy compensation + materials compensation + service compensation))

Adds two columns to df:
- markup_level: interpolated to quarterly
- markup_growth: first difference of markup

Requires: data/KLEMS.xlsx
"""
function build_markup_method2_no_hp_latest(df::DataFrame, file_path::String)
    df = copy(df)
    
    # 1. load data
    markup_data = XLSX.readxlsx(file_path)
    
    function get_long_df(sheet_name, value_col_name)
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

    df_e  = get_long_df("Energy Compensation", :energy)
    df_m  = get_long_df("Materials Compensation", :material)
    df_s  = get_long_df("Service Compensation", :service)
    df_go = get_long_df("Gross Output", :go)
    df_va = get_long_df("Value Added", :va)

    merged = innerjoin(df_e, df_m, on=[:Industry, :year])
    merged = innerjoin(merged, df_s, on=[:Industry, :year])
    merged = innerjoin(merged, df_go, on=[:Industry, :year])
    merged = innerjoin(merged, df_va, on=[:Industry, :year])

    # 2. calculate markup
    merged.markup_ind = log.(merged.go ./ (merged.energy .+ merged.material .+ merged.service))

    va_avg        = combine(groupby(merged, :Industry), :va => mean => :va_mean)
    va_avg.weight = va_avg.va_mean ./ sum(va_avg.va_mean)
    merged        = leftjoin(merged, va_avg[:, [:Industry, :weight]], on=:Industry)

    agg_markup = combine(groupby(merged, :year)) do sdf
        valid = .!isnan.(sdf.markup_ind) .& .!isinf.(sdf.markup_ind) .& .!ismissing.(sdf.markup_ind)
        return (annual_markup = sum(sdf.markup_ind[valid] .* sdf.weight[valid]),)
    end

    sort!(agg_markup, :year)

    # 3. Fernandez method interpolation to quarterly
    # --------------------------------------------------
    # Fernandez (1981): p = q + D'(DD')^{-1}(y - Dq)
    # where q is a random-walk prior (first-difference smoother),
    # D is the aggregation matrix (sum of 4 quarters = annual),
    # y is the annual observed series.
    # --------------------------------------------------
    function fernandez_interpolate(y::Vector{Float64}, m::Int=4)
        n = length(y)
        N = n * m

        # --- 1. First-difference matrix D1: (N-1) x N ---
        D1 = zeros(N-1, N)
        for i in 1:N-1
            D1[i, i]   =  1.0
            D1[i, i+1] = -1.0
        end

        # --- 2. Prior covariance Σ = (D1'D1)^+ ---
        AtA   = D1' * D1
        Sigma = pinv(AtA)   # N x N

        # --- 3. Aggregation matrix C: n x N ---
        C = zeros(n, N)
        for i in 1:n
            C[i, (i-1)*m+1 : i*m] .= 1.0
        end

        # --- 4. High-frequency indicator vector x (constant term) ---
        # Without external indicator, use x = ones(N)
        # Corresponding low-frequency aggregation: X = C * x = m * ones(n)
        x = ones(N)       # N x 1
        X = C * x         # n x 1, each element = m = 4

        # --- 5. Estimate β (GLS) ---
        # β̂ = (X' (CΣC')^{-1} X)^{-1}  X' (CΣC')^{-1} y
        CSigmaC     = C * Sigma * C'          # n x n
        CSigmaC_inv = pinv(CSigmaC)

        beta_num = (X' * CSigmaC_inv * X)   # scalar
        beta_den = (X' * CSigmaC_inv * y)   # scalar
        beta_hat = beta_den / beta_num       # scalar

        # --- 6. Full Fernandez formula ---
        # p̂ = β̂·x + Σ C' (CΣC')^{-1} (y - C·β̂·x)
        trend     = beta_hat .* x                           # N x 1, high-frequency trend
        residuals = y .- C * trend                          # n x 1, annual residuals
        p_hat     = trend .+ Sigma * C' * (CSigmaC_inv * residuals)

        return p_hat

        # --- 7. Verify annual aggregation constraint ---
        reconstructed = C * p_hat
        max_error = maximum(abs.(reconstructed .- y))
        if max_error > 1e-6
            @warn "Annual aggregation constraint not satisfied, max error = $max_error"
        else
            @info "Annual aggregation constraint satisfied, max error = $max_error"
        end
    end

    y_annual = Float64.(agg_markup.annual_markup)
    interpolated_values = fernandez_interpolate(y_annual, 4)

    # Build quarterly date pool
    quarterly_markup_pool = DataFrame(
        "observation_date" => Date[],
        "markup_level" => Float64[]
    )

    current_idx = 1
    for yr in agg_markup.year
        y = Int(yr)
        for m in [1, 4, 7, 10]
            if current_idx <= length(interpolated_values)
                push!(quarterly_markup_pool, (
                    observation_date = Date(y, m, 1),
                    markup_level = interpolated_values[current_idx]
                ))
                current_idx += 1
            end
        end
    end

    # 4. join to input DataFrame
    df_result = leftjoin(df, quarterly_markup_pool, on="observation_date")
    sort!(df_result, "observation_date")
    
    # 5. compute quarterly growth
    markup_vals = df_result.markup_level
    growth = Vector{Union{Missing, Float64}}(missing, nrow(df_result))
    
    for i in 2:nrow(df_result)
        if !ismissing(markup_vals[i]) && !ismissing(markup_vals[i-1])
            growth[i] = markup_vals[i] - markup_vals[i-1]
        end
    end
    
    df_result.markup_growth = growth

    return df_result
end