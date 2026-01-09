using SpecialFunctions

function beta_pdf(x, a, b)
    # PURPOSE: pdf of the beta(a,b) distribution
    # x: scalar, vector, or matrix
    # a, b: scalar beta parameters
    
    # --- input check ---
    if a <= 0 || b <= 0
        error("Parameter a or b is nonpositive")
    end

    # --- ensure beta(a,b) is a Float64 ---
    B = beta(a, b)
    B = B isa Tuple ? B[1] : B

    # --- calculate pdf ---
    pdf = x .^ (a - 1) .* (1 .- x) .^ (b - 1) ./ B

    # --- set out-of-range values to 0 ---
    pdf = ifelse.((x .< 0) .| (x .> 1), 0.0, pdf)

    return pdf
end