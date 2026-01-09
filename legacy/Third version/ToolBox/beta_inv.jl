using SpecialFunctions

function beta_inv(p, a, b)
    # input validation
    if a <= 0 || b <= 0
        error("beta_inv: parameter a or b is nonpositive")
    end
    if any(abs.(2 .* p .- 1) .> 1)
        error("beta_inv: p must satisfy 0 ≤ p ≤ 1")
    end

    # --- ensure p is array (avoid tuple errors) ---
    scalar_flag = false
    if isa(p, Number)
        p = [p]           # convert scalar → vector
        scalar_flag = true
    else
        p = collect(p)    # ensure it's Array, not Tuple
    end

    # --- initial guess, same shape as p ---
    x = fill(a / (a + b), length(p))
    dx = ones(length(p))

    while any(abs.(dx) .> 256 * eps() .* max.(x, 1))
        # beta_inc returns (value, error), take only the value [1]
        cdf_val = [beta_inc(a, b, xi)[1] for xi in x]  
        pdf_val = beta_pdf.(x, a, b)
        dx = (cdf_val .- p) ./ pdf_val

        x .= x .- dx

        # protect boundaries
        x .= max.(x, 0.0)
        x .= min.(x, 1.0)
    end

    return scalar_flag ? x[1] : x
end