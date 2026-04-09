function run_wald_tests(df_agg::DataFrame,
                         shock_df::DataFrame,
                         macro_q::DataFrame;
                         nlag::Int  = 2,
                         nhorz::Int = 8,
                         controls   = (:ln_gdp_diff, :pi_p, :Interest))

    outcomes = [:mean_income, :mean_hours, :mean_wage, :unemp_rate]
    outcome_labels = Dict(
        :mean_income => "log Income",
        :mean_hours  => "log Hours",
        :mean_wage   => "log Wage",
        :unemp_rate  => "Unemployment",
    )

    # ── Build stacked panel (原始列，去掉 moving_avg) ───────────────
    df_panel = DataFrame()
    for occ_id in 1:9
        df_occ = filter(r -> r.occ_cat == occ_id, df_agg)
        sort!(df_occ, :quarter)
        dropmissing!(df_occ)
        df_panel = vcat(df_panel, df_occ)
    end

    # ── Merge with shock and macro data ─────────────────────────────
    df_panel = innerjoin(df_panel, shock_df; on=:quarter)
    df_panel = innerjoin(df_panel, macro_q;  on=:quarter)
    dropmissing!(df_panel)
    sort!(df_panel, [:occ_cat, :quarter])

    # ── Joint test table ─────────────────────────────────────────────
    W = 72
    println("\n" * "="^W)
    println("JOINT WALD TEST: H₀: IRF equal across all 9 occupations")
    println("="^W)
    @printf("%-14s  %8s  %6s  %8s\n", "Outcome", "W stat", "df", "p-value")
    println("-"^W)

    wald_results = Dict()
    for outcome in outcomes
        res = wald_test_joint_irfs(df_panel, outcome;
                                   nlag=nlag, nhorz=nhorz, controls=controls)
        wald_results[outcome] = res
        sig = res.p_value < 0.01 ? "***" :
              res.p_value < 0.05 ? "**"  :
              res.p_value < 0.10 ? "*"   : ""
        @printf("%-14s  %8.3f  %6d  %8.4f %s\n",
                outcome_labels[outcome], res.W_stat, res.df, res.p_value, sig)
    end
    println("-"^W)
    println("*** p<0.01  ** p<0.05  * p<0.10")

    # ── Per-horizon table ─────────────────────────────────────────────
    println("\nPer-horizon p-values (H₀: β equal across 9 occs at horizon h):")
    @printf("%-6s  %-12s  %-12s  %-12s  %-12s\n",
            "Horizon", "log Income", "log Hours", "log Wage", "Unemployment")
    println("-"^W)
    for hi in 1:nhorz
        @printf("h=%-4d  %12.4f  %12.4f  %12.4f  %12.4f\n",
                hi,
                wald_results[:mean_income].p_by_h[hi], 
                wald_results[:mean_hours].p_by_h[hi],
                wald_results[:mean_wage].p_by_h[hi],
                wald_results[:unemp_rate].p_by_h[hi])
    end

    return wald_results
end