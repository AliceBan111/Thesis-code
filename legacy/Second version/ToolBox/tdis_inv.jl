function tdis_inv(p::Union{Float64, Vector{Float64}}, a::Real)
    a = Float64(a)
    s = p .< 0.5
    p = p .+ (1 .- 2 .* p) .* s
    p = 1 .- (2 .* (1 .- p))
    x = beta_inv(p, 1/2, a/2)
    x = x .* a ./ (1 .- x)
    x = (1 .- 2 .* s) .* sqrt.(x)
    return x
end