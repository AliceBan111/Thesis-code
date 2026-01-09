using SpecialFunctions

function tdis_prb(x::Union{Float64, Vector{Float64}}, n::Real)
    # =======================================================================
    # calculates t-probabilities for elements in x-vector
    # =======================================================================
    
    n = Float64(n) 
    
    if n <= 0
        error("dof is negative or zero in tdis_prb")
    end

    scalar_flag = (x isa Float64)
    x_vec = scalar_flag ? [x] : x
    
    x2 = n ./ (n .+ x_vec.^2)

    x2 = clamp.(x2, 1e-12, 1 - 1e-12)

    tmp = 1.0 .- 0.5 .* [beta_inc(0.5*n, 0.5, xi)[1] for xi in x2]
    y = 2 .* (1 .- tmp)
    
    return scalar_flag ? y[1] : y
end