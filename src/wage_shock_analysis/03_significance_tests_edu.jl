# =============================================================================
# 03_significance_tests_edu.jl
# Bootstrap significance tests for LP impulse responses
#   - Per-horizon pointwise p-values for each education group
#   - Pairwise group difference tests at each horizon
# =============================================================================

include("02_lp_estimation_edu.jl")

const EDU_GROUPS = [1, 2, 3]
const PAIRS      = [(1,2), (1,3), (2,3)]

# =============================================================================
# 1. POINTWISE SIGNIFICANCE  — each group, each horizon
# =============================================================================
# H0: β_h = 0
# p-value: bootstrap percentile, two-sided
# p = 2 * min(P(β* > 0), P(β* < 0))
#
# Note: bootstrap draws are already centered around β̂, not zero.
# For percentile p-value we use the shifted draws: β* − β̂
# so that H0 is evaluated at zero.

function pointwise_significance(all_results::Dict{Int, DataFrame},
                                 all_boots::Dict{Int, Dict{Int, Matrix{Float64}}},
                                 coef_names::Vector{Symbol};
                                 h_max::Int = H_MAX)::DataFrame

    β_idx = findfirst(==(:shock), coef_names)
    β_idx === nothing && error(":shock not found in coef_names")

    rows = NamedTuple[]

    for g in EDU_GROUPS
        results_g = all_results[g]
        boots_g   = all_boots[g]

        for h in 0:h_max
            haskey(boots_g, h) || continue

            sub = filter(r -> r.horizon == h && r.coef_name == :shock, results_g)
            nrow(sub) == 0 && continue

            β̂     = sub.beta[1]
            draws  = boots_g[h][:, β_idx]

            # Shift draws to be centered under H0: β = 0
            shifted = draws .- mean(draws) 

            p_two  = 2 * min(mean(shifted .> β̂), mean(shifted .< β̂))
            p_two  = clamp(p_two, 0.0, 1.0)

            push!(rows, (
                edu_group = g,
                horizon   = h,
                beta      = β̂,
                boot_se   = std(draws),
                p_value   = p_two,
                sig_90    = p_two < 0.10,
                sig_95    = p_two < 0.05,
                sig_68    = p_two < 0.32,
            ))
        end
    end

    return DataFrame(rows)
end

# =============================================================================
# 2. PAIRWISE GROUP DIFFERENCE TESTS  — each pair, each horizon
# =============================================================================
# H0: β_h^(g1) − β_h^(g2) = 0
#
# Since the two groups are estimated independently, their bootstrap draws
# are independent — direct subtraction gives the correct null distribution.
#
# Shifted difference: d* = (β*_g1 − β*_g2) − (β̂_g1 − β̂_g2)
# p-value: 2 * min(P(d* > d̂), P(d* < d̂))  where d̂ = β̂_g1 − β̂_g2

function pairwise_difference_tests(all_results::Dict{Int, DataFrame},
                                    all_boots::Dict{Int, Dict{Int, Matrix{Float64}}},
                                    coef_names::Vector{Symbol};
                                    h_max::Int = H_MAX)::DataFrame

    β_idx = findfirst(==(:shock), coef_names)
    β_idx === nothing && error(":shock not found in coef_names")

    rows = NamedTuple[]

    for (g1, g2) in PAIRS
        for h in 0:h_max
            haskey(all_boots[g1], h) || continue
            haskey(all_boots[g2], h) || continue

            sub1 = filter(r -> r.horizon == h && r.coef_name == :shock, all_results[g1])
            sub2 = filter(r -> r.horizon == h && r.coef_name == :shock, all_results[g2])
            (nrow(sub1) == 0 || nrow(sub2) == 0) && continue

            β̂1 = sub1.beta[1]
            β̂2 = sub2.beta[1]
            d̂  = β̂1 - β̂2

            draws1 = all_boots[g1][h][:, β_idx]
            draws2 = all_boots[g2][h][:, β_idx]

            # Match lengths in case n_valid differs slightly between groups
            n      = min(length(draws1), length(draws2))
            diff   = draws1[1:n] .- draws2[1:n]

            # Shift to center under H0
            shifted = diff .- mean(diff)

            p_two = 2 * min(mean(shifted .> d̂), mean(shifted .< d̂))
            p_two = clamp(p_two, 0.0, 1.0)

            push!(rows, (
                g1      = g1,
                g2      = g2,
                horizon = h,
                diff    = d̂,
                boot_se = std(diff),
                ci_lo95 = quantile(diff, 0.025),
                ci_hi95 = quantile(diff, 0.975),
                ci_lo90 = quantile(diff, 0.05),
                ci_hi90 = quantile(diff, 0.95),
                p_value = p_two,
                sig_90  = p_two < 0.10,
                sig_95  = p_two < 0.05,
                sig_68  = p_two < 0.32,
            ))
        end
    end

    return DataFrame(rows)
end

# =============================================================================
# 3. MAIN
# =============================================================================
function run_significance_tests(all_results::Dict{Int, DataFrame},
                                 all_boots::Dict{Int, Dict{Int, Matrix{Float64}}},
                                 all_coefnames::Vector{Symbol})

    println("\nRunning pointwise significance tests...")
    pw = pointwise_significance(all_results, all_boots, all_coefnames)
    println("  Done. $(nrow(pw)) group-horizon cells.")

    println("Running pairwise difference tests...")
    pd = pairwise_difference_tests(all_results, all_boots, all_coefnames)
    println("  Done. $(nrow(pd)) pair-horizon cells.")

    return pw, pd
end