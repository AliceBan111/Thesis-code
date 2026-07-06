# Labor Adjustment Sequencing: Thesis Interpretation

The empirical exercise compares the early response of hours with the later response
of unemployment across occupation-industry cells. A negative hours response is
interpreted as adjustment on the intensive margin: firms reduce labor input by
cutting hours before reducing headcount. A positive unemployment response at
longer horizons is interpreted as adjustment on the extensive margin.

The baseline specification uses the cumulative hours response over horizons 0-12
and the cumulative unemployment response over horizons 24-36. In this
specification, the Pearson correlation between early hours and later unemployment
is -0.460 (p = 0.005). Because the
hours measure is signed in levels, a negative correlation means that cells with
more negative early hours responses tend to have larger positive unemployment
responses later. Equivalently, if early hours declines are multiplied by -1, the
sequencing hypothesis predicts a positive association between the intensity of
the early hours decline and subsequent unemployment.

The baseline OLS regression of later unemployment on early hours gives a slope of
-0.198 (heteroskedasticity-robust p = 0.004). A negative coefficient is
consistent with sequencing: lower early hours are associated with higher later
unemployment. Industry-fixed-effect specifications ask whether the same pattern
holds after comparing occupation cells within broad industries.

For sign consistency, 31 out of
35 cells (88.6%)
show both an early decline in hours and a later increase in unemployment. Using
pointwise 90 percent confidence intervals, 29
cells (82.9%) show both a
significantly negative hours response in the early window and a significantly
positive unemployment response in the later window.

For cumulative windows, significance is classified using the pointwise confidence
intervals within the relevant window: a cumulative hours decline is flagged as
significant when the cumulative response is negative and at least one horizon in
the window has a 90 percent confidence interval entirely below zero; the analogous
rule is used for positive unemployment. This is a transparent descriptive rule,
not a joint test of the cumulative impulse response, because the covariance
matrix across horizons is not available in the CSV files.

Overall, evidence supporting the intensive-before-extensive adjustment mechanism
would consist of: (i) negative correlations between early hours and later
unemployment, (ii) negative OLS coefficients when unemployment is regressed on
early hours, (iii) positive coefficients when unemployment is regressed on the
early-hours-decline measure, and (iv) a large share of cells with hours < 0 and
unemployment > 0. Robustness windows test whether this pattern is specific to
one horizon choice or appears across nearby definitions of early and later
adjustment.
