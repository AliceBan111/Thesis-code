function compute_dynamic_multiplier(r_wage, r_unemp; nhorz=8)

    shock_idx = 1  # markup_shock is the first xname

    # Extract wage IRF coefficients
    betas_wage  = [r_wage.B[shock_idx,  1, h] for h in 1:nhorz]

    # Extract unemployment IRF coefficients
    betas_unemp = [r_unemp.B[shock_idx, 1, h] for h in 1:nhorz]

    # Cumulative sums
    cum_wage  = cumsum(betas_wage)
    cum_unemp = cumsum(betas_unemp)

    # Dynamic multiplier at each horizon H
    multipliers = [abs(cum_unemp[h]) > 1e-10 ? cum_wage[h] / cum_unemp[h] : missing for h in 1:nhorz]

    # Summary table
    println("\nDynamic Multiplier Φ_w(h) = Σwage_IRF / Σunemp_IRF")
    println("H  |  ΣWage IRF     |  ΣUnemp IRF    |  Φ_w(H)")
    for h in 1:nhorz
        println("$h  |  $(round(cum_wage[h],  sigdigits=4))  |  $(round(cum_unemp[h], sigdigits=4))  |  $(round(multipliers[h], sigdigits=4))")
    end

    return multipliers, cum_wage, cum_unemp
end

