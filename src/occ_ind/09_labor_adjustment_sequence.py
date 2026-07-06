"""
Labor adjustment sequencing analysis for occupation-industry IRFs.

This script tests whether oil-news supply shocks first reduce hours
(intensive margin) and later increase unemployment (extensive margin)
across occupation-industry cells.

Default input:
    result/occ_ind_4

Default output:
    result/occ_ind_4/labor_adjustment_sequence
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats


REQUIRED_COLUMNS = ["horizon", "beta", "ci_lo90", "ci_hi90", "ci_lo68", "ci_hi68"]

OCC_LABELS = {
    1: "Managerial",
    2: "Professional specialty",
    3: "High tech",
    4: "Sales",
    5: "Administrative support",
    6: "Service",
    7: "Farming/forestry/construction",
    8: "Precision production/repair",
    9: "Machine operators/transport",
}

INDUSTRY_LABELS = {
    "Energy_intensive": "Energy intensive",
    "Manufacturing_Construction": "Mfg. + construction",
    "Services": "Services",
    "Trade": "Trade",
}

INDUSTRY_MARKERS = {
    "Energy_intensive": "o",
    "Manufacturing_Construction": "s",
    "Services": "D",
    "Trade": "^",
}


@dataclass(frozen=True)
class SequenceSpec:
    name: str
    hours_kind: str
    hours_start: int
    hours_end: int
    unemployment_kind: str
    unemployment_start: int
    unemployment_end: int
    description: str


SPECS = [
    SequenceSpec(
        name="A_cumulative_h0_12_u24_36",
        hours_kind="cumulative",
        hours_start=0,
        hours_end=12,
        unemployment_kind="cumulative",
        unemployment_start=24,
        unemployment_end=36,
        description="Cumulative hours response at H=0-12 and cumulative unemployment response at H=24-36.",
    ),
    SequenceSpec(
        name="B_endpoint_h12_u36",
        hours_kind="endpoint",
        hours_start=12,
        hours_end=12,
        unemployment_kind="endpoint",
        unemployment_start=36,
        unemployment_end=36,
        description="Endpoint hours response at H=12 and endpoint unemployment response at H=36.",
    ),
    SequenceSpec(
        name="R1_cumulative_h0_6_u18_30",
        hours_kind="cumulative",
        hours_start=0,
        hours_end=6,
        unemployment_kind="cumulative",
        unemployment_start=18,
        unemployment_end=30,
        description="Robustness: shorter early hours window and earlier later-unemployment window.",
    ),
    SequenceSpec(
        name="R2_cumulative_h0_18_u18_36",
        hours_kind="cumulative",
        hours_start=0,
        hours_end=18,
        unemployment_kind="cumulative",
        unemployment_start=18,
        unemployment_end=36,
        description="Robustness: broader early and later windows.",
    ),
    SequenceSpec(
        name="R3_cumulative_h6_12_u30_36",
        hours_kind="cumulative",
        hours_start=6,
        hours_end=12,
        unemployment_kind="cumulative",
        unemployment_start=30,
        unemployment_end=36,
        description="Robustness: delayed early-hours window and late unemployment window.",
    ),
    SequenceSpec(
        name="R4_endpoint_h6_u24",
        hours_kind="endpoint",
        hours_start=6,
        hours_end=6,
        unemployment_kind="endpoint",
        unemployment_start=24,
        unemployment_end=24,
        description="Robustness: earlier endpoint comparison.",
    ),
]


def pretty_label(value: str) -> str:
    return INDUSTRY_LABELS.get(value, value.replace("_", " "))


def parse_irf_path(path: Path, root: Path) -> dict[str, object] | None:
    rel = path.relative_to(root)
    if len(rel.parts) < 3:
        return None

    outcome, industry, filename = rel.parts[0], rel.parts[1], rel.parts[-1]
    match = re.match(r"irf_occ(?P<occ_id>\d+)_(?P<occupation>.+)\.csv$", filename)
    if not match:
        return None

    occ_id = int(match.group("occ_id"))
    occupation = match.group("occupation").replace("_", " ")

    return {
        "outcome": outcome,
        "industry": industry,
        "industry_label": pretty_label(industry),
        "occ_id": occ_id,
        "occupation": OCC_LABELS.get(occ_id, occupation),
        "source_file": str(path),
    }


def load_irf_panel(root: Path) -> pd.DataFrame:
    rows = []
    for path in sorted(root.rglob("*.csv")):
        parsed = parse_irf_path(path, root)
        if parsed is None:
            continue

        df = pd.read_csv(path)
        missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
        if missing:
            raise ValueError(f"{path} is missing required columns: {missing}")

        df = df[REQUIRED_COLUMNS].copy()
        for col in REQUIRED_COLUMNS:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        df = df.dropna(subset=["horizon", "beta"])
        df["horizon"] = df["horizon"].astype(int)

        for key, value in parsed.items():
            df[key] = value
        rows.append(df)

    if not rows:
        raise FileNotFoundError(f"No IRF occupation CSV files found under {root}")

    panel = pd.concat(rows, ignore_index=True)
    panel["cell_id"] = (
        panel["industry"].astype(str)
        + "__occ"
        + panel["occ_id"].astype(str)
        + "_"
        + panel["occupation"].str.replace(r"\W+", "_", regex=True).str.strip("_")
    )
    return panel.sort_values(["outcome", "industry", "occ_id", "horizon"]).reset_index(drop=True)


def summarize_response(df: pd.DataFrame, outcome: str, kind: str, start: int, end: int) -> pd.DataFrame:
    sub = df[(df["outcome"] == outcome) & (df["horizon"].between(start, end))].copy()
    if sub.empty:
        raise ValueError(f"No rows for outcome={outcome}, H={start}-{end}")

    id_cols = ["industry", "industry_label", "occ_id", "occupation", "cell_id"]
    grouped = sub.groupby(id_cols, dropna=False)

    if kind == "cumulative":
        out = grouped.agg(
            value=("beta", "sum"),
            mean_beta=("beta", "mean"),
            n_horizons=("horizon", "count"),
            min_horizon=("horizon", "min"),
            max_horizon=("horizon", "max"),
            sig_neg90_any=("ci_hi90", lambda x: bool((x < 0).any())),
            sig_pos90_any=("ci_lo90", lambda x: bool((x > 0).any())),
            sig_neg68_any=("ci_hi68", lambda x: bool((x < 0).any())),
            sig_pos68_any=("ci_lo68", lambda x: bool((x > 0).any())),
            sig_neg90_horizons=("ci_hi90", lambda x: int((x < 0).sum())),
            sig_pos90_horizons=("ci_lo90", lambda x: int((x > 0).sum())),
        ).reset_index()
    elif kind == "endpoint":
        sub = sub[sub["horizon"] == start].copy()
        out = sub[id_cols + REQUIRED_COLUMNS].rename(columns={"beta": "value"})
        out["mean_beta"] = out["value"]
        out["n_horizons"] = 1
        out["min_horizon"] = start
        out["max_horizon"] = end
        out["sig_neg90_any"] = out["ci_hi90"] < 0
        out["sig_pos90_any"] = out["ci_lo90"] > 0
        out["sig_neg68_any"] = out["ci_hi68"] < 0
        out["sig_pos68_any"] = out["ci_lo68"] > 0
        out["sig_neg90_horizons"] = out["sig_neg90_any"].astype(int)
        out["sig_pos90_horizons"] = out["sig_pos90_any"].astype(int)
    else:
        raise ValueError(f"Unknown summary kind: {kind}")

    out["outcome"] = outcome
    out["measure_kind"] = kind
    out["window_start"] = start
    out["window_end"] = end
    return out


def build_sequence_measures(panel: pd.DataFrame, specs: list[SequenceSpec]) -> pd.DataFrame:
    rows = []
    for spec in specs:
        hours = summarize_response(
            panel, "hours", spec.hours_kind, spec.hours_start, spec.hours_end
        )
        unemp = summarize_response(
            panel,
            "unemployment",
            spec.unemployment_kind,
            spec.unemployment_start,
            spec.unemployment_end,
        )

        keep = ["industry", "industry_label", "occ_id", "occupation", "cell_id"]
        merged = hours.merge(unemp, on=keep, suffixes=("_hours", "_unemployment"))
        merged["spec"] = spec.name
        merged["spec_description"] = spec.description
        merged["early_hours"] = merged["value_hours"]
        merged["early_hours_decline"] = -merged["early_hours"]
        merged["later_unemployment"] = merged["value_unemployment"]
        merged["hours_negative"] = merged["early_hours"] < 0
        merged["unemployment_positive"] = merged["later_unemployment"] > 0
        merged["sequence_sign_consistent"] = (
            merged["hours_negative"] & merged["unemployment_positive"]
        )
        merged["hours_sig_neg90"] = (
            merged["hours_negative"] & merged["sig_neg90_any_hours"]
        )
        merged["unemployment_sig_pos90"] = (
            merged["unemployment_positive"] & merged["sig_pos90_any_unemployment"]
        )
        merged["sequence_sig90_consistent"] = (
            merged["hours_sig_neg90"] & merged["unemployment_sig_pos90"]
        )
        merged["hours_sig_neg68"] = (
            merged["hours_negative"] & merged["sig_neg68_any_hours"]
        )
        merged["unemployment_sig_pos68"] = (
            merged["unemployment_positive"] & merged["sig_pos68_any_unemployment"]
        )
        merged["sequence_sig68_consistent"] = (
            merged["hours_sig_neg68"] & merged["unemployment_sig_pos68"]
        )
        rows.append(merged)

    return pd.concat(rows, ignore_index=True)


def write_coverage_diagnostics(panel: pd.DataFrame, output_dir: Path) -> None:
    cell_cols = ["industry", "industry_label", "occ_id", "occupation", "cell_id"]
    cells = panel[cell_cols].drop_duplicates()
    outcomes = sorted(panel["outcome"].unique())

    coverage = []
    for outcome in outcomes:
        present = set(panel.loc[panel["outcome"] == outcome, "cell_id"])
        coverage.append(
            {
                "outcome": outcome,
                "n_cells": len(present),
                "missing_relative_to_all_cells": len(set(cells["cell_id"]) - present),
            }
        )
    pd.DataFrame(coverage).to_csv(output_dir / "coverage_by_outcome.csv", index=False)

    required = {"hours", "unemployment"}
    if required.issubset(outcomes):
        hours_cells = set(panel.loc[panel["outcome"] == "hours", "cell_id"])
        unemployment_cells = set(panel.loc[panel["outcome"] == "unemployment", "cell_id"])
        diagnostic_rows = []
        for cell_id in sorted(hours_cells | unemployment_cells):
            row = cells[cells["cell_id"] == cell_id].iloc[0].to_dict()
            row["has_hours"] = cell_id in hours_cells
            row["has_unemployment"] = cell_id in unemployment_cells
            row["included_in_sequence_tests"] = row["has_hours"] and row["has_unemployment"]
            diagnostic_rows.append(row)
        pd.DataFrame(diagnostic_rows).to_csv(
            output_dir / "coverage_hours_unemployment_cells.csv", index=False
        )


def pearson_spearman(df: pd.DataFrame) -> pd.Series:
    if len(df) < 3:
        return pd.Series(
            {
                "pearson_r": np.nan,
                "pearson_p": np.nan,
                "spearman_r": np.nan,
                "spearman_p": np.nan,
            }
        )
    pearson = stats.pearsonr(df["early_hours"], df["later_unemployment"])
    spearman = stats.spearmanr(df["early_hours"], df["later_unemployment"])
    return pd.Series(
        {
            "pearson_r": pearson.statistic,
            "pearson_p": pearson.pvalue,
            "spearman_r": spearman.statistic,
            "spearman_p": spearman.pvalue,
        }
    )


def run_regressions(sequence: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    coef_tables = []

    for spec, sub in sequence.groupby("spec", sort=False):
        data = sub.replace([np.inf, -np.inf], np.nan).dropna(
            subset=["early_hours", "early_hours_decline", "later_unemployment"]
        )
        if len(data) < 3:
            continue

        models = {
            "ols_raw_hours": "later_unemployment ~ early_hours",
            "ols_decline_intensity": "later_unemployment ~ early_hours_decline",
            "ols_decline_plus_industry_fe": "later_unemployment ~ early_hours_decline + C(industry)",
        }

        for model_name, formula in models.items():
            model = smf.ols(formula, data=data).fit(cov_type="HC1")
            rows.append(
                {
                    "spec": spec,
                    "model": model_name,
                    "formula": formula,
                    "n": int(model.nobs),
                    "r2": model.rsquared,
                    "adj_r2": model.rsquared_adj,
                    "f_pvalue": model.f_pvalue,
                }
            )
            coef = model.summary2().tables[1].reset_index().rename(columns={"index": "term"})
            coef.insert(0, "model", model_name)
            coef.insert(0, "spec", spec)
            coef_tables.append(coef)

    model_summary = pd.DataFrame(rows)
    coefficient_table = pd.concat(coef_tables, ignore_index=True) if coef_tables else pd.DataFrame()
    return model_summary, coefficient_table


def make_tables(sequence: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    corr_rows = []
    for spec, sub in sequence.groupby("spec", sort=False):
        row = pearson_spearman(sub).to_dict()
        row["spec"] = spec
        corr_rows.append(row)
    corr = pd.DataFrame(corr_rows)[
        ["spec", "pearson_r", "pearson_p", "spearman_r", "spearman_p"]
    ]

    sign = (
        sequence.groupby("spec", sort=False)
        .agg(
            n_cells=("cell_id", "nunique"),
            hours_negative_n=("hours_negative", "sum"),
            unemployment_positive_n=("unemployment_positive", "sum"),
            sequence_sign_consistent_n=("sequence_sign_consistent", "sum"),
            hours_sig_neg90_n=("hours_sig_neg90", "sum"),
            unemployment_sig_pos90_n=("unemployment_sig_pos90", "sum"),
            sequence_sig90_consistent_n=("sequence_sig90_consistent", "sum"),
            hours_sig_neg68_n=("hours_sig_neg68", "sum"),
            unemployment_sig_pos68_n=("unemployment_sig_pos68", "sum"),
            sequence_sig68_consistent_n=("sequence_sig68_consistent", "sum"),
        )
        .reset_index()
    )

    for col in [c for c in sign.columns if c.endswith("_n")]:
        share_col = col.removesuffix("_n") + "_share"
        sign[share_col] = sign[col] / sign["n_cells"]

    return corr, sign


def save_publication_tables(
    corr: pd.DataFrame,
    sign: pd.DataFrame,
    model_summary: pd.DataFrame,
    coefficients: pd.DataFrame,
    output_dir: Path,
) -> None:
    corr.to_csv(output_dir / "table_01_correlations.csv", index=False)
    sign.to_csv(output_dir / "table_02_sign_consistency_and_significance.csv", index=False)
    model_summary.to_csv(output_dir / "table_03_ols_model_summary.csv", index=False)
    coefficients.to_csv(output_dir / "table_04_ols_coefficients.csv", index=False)

    with pd.ExcelWriter(output_dir / "labor_adjustment_sequence_tables.xlsx") as writer:
        corr.to_excel(writer, sheet_name="correlations", index=False)
        sign.to_excel(writer, sheet_name="sign_consistency", index=False)
        model_summary.to_excel(writer, sheet_name="ols_summary", index=False)
        coefficients.to_excel(writer, sheet_name="ols_coefficients", index=False)

    latex_specs = [
        (corr, "table_01_correlations.tex", "Correlation between early hours and later unemployment responses"),
        (sign, "table_02_sign_consistency_and_significance.tex", "Sign consistency and pointwise significance counts"),
        (model_summary, "table_03_ols_model_summary.tex", "OLS model summary"),
    ]
    for df, filename, caption in latex_specs:
        write_simple_latex_table(df, output_dir / filename, caption)


def latex_escape(value: object) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, float):
        text = f"{value:.4f}"
    else:
        text = str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def write_simple_latex_table(df: pd.DataFrame, path: Path, caption: str) -> None:
    """Write a plain LaTeX table without pandas Styler/Jinja2."""
    col_spec = "l" * len(df.columns)
    lines = [
        r"\begin{table}",
        r"\centering",
        rf"\caption{{{latex_escape(caption)}}}",
        rf"\begin{{tabular}}{{{col_spec}}}",
        r"\hline",
        " & ".join(latex_escape(col) for col in df.columns) + r" \\",
        r"\hline",
    ]
    for _, row in df.iterrows():
        lines.append(" & ".join(latex_escape(value) for value in row) + r" \\")
    lines.extend([r"\hline", r"\end{tabular}", r"\end{table}", ""])
    path.write_text("\n".join(lines))


def add_regression_line(ax, data: pd.DataFrame) -> None:
    model = smf.ols("later_unemployment ~ early_hours", data=data).fit()
    x_grid = np.linspace(data["early_hours"].min(), data["early_hours"].max(), 150)
    pred = model.get_prediction(pd.DataFrame({"early_hours": x_grid})).summary_frame(alpha=0.10)
    ax.plot(x_grid, pred["mean"], color="black", linewidth=1.4, label="OLS fit")
    ax.fill_between(
        x_grid,
        pred["mean_ci_lower"].to_numpy(),
        pred["mean_ci_upper"].to_numpy(),
        color="black",
        alpha=0.12,
        linewidth=0,
    )


def plot_scatter(
    data: pd.DataFrame,
    spec: str,
    output_dir: Path,
    label_by: str = "occupation",
) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 6.2), constrained_layout=True)

    for industry, sub in data.groupby("industry"):
        ax.scatter(
            sub["early_hours"],
            sub["later_unemployment"],
            s=62,
            alpha=0.86,
            marker=INDUSTRY_MARKERS.get(industry, "o"),
            label=pretty_label(industry),
        )

    add_regression_line(ax, data)
    ax.axvline(0, color="0.25", linewidth=0.8)
    ax.axhline(0, color="0.25", linewidth=0.8)

    x_span = data["early_hours"].max() - data["early_hours"].min()
    y_span = data["later_unemployment"].max() - data["later_unemployment"].min()
    x_pad = 0.012 * (x_span if not np.isclose(x_span, 0) else 1.0)
    y_pad = 0.012 * (y_span if not np.isclose(y_span, 0) else 1.0)

    for _, row in data.iterrows():
        label = row["occupation"] if label_by == "occupation" else pretty_label(row["industry"])
        ax.text(
            row["early_hours"] + x_pad,
            row["later_unemployment"] + y_pad,
            label,
            fontsize=6.8,
            alpha=0.78,
        )

    ax.set_title(spec)
    ax.set_xlabel("Early hours response")
    ax.set_ylabel("Later unemployment response")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(output_dir / f"scatter_{spec}_{label_by}.png", dpi=300)
    fig.savefig(output_dir / f"scatter_{spec}_{label_by}.pdf")
    plt.close(fig)


def make_plots(sequence: pd.DataFrame, output_dir: Path) -> None:
    plot_dir = output_dir / "figures"
    plot_dir.mkdir(parents=True, exist_ok=True)

    for spec, data in sequence.groupby("spec", sort=False):
        data = data.replace([np.inf, -np.inf], np.nan).dropna(
            subset=["early_hours", "later_unemployment"]
        )
        if data.empty:
            continue
        plot_scatter(data, spec, plot_dir, label_by="occupation")
        plot_scatter(data, spec, plot_dir, label_by="industry")


def write_interpretation(
    output_dir: Path,
    corr: pd.DataFrame,
    sign: pd.DataFrame,
    coefficients: pd.DataFrame,
) -> None:
    main_spec = "A_cumulative_h0_12_u24_36"
    main_corr = corr[corr["spec"] == main_spec].iloc[0]
    main_sign = sign[sign["spec"] == main_spec].iloc[0]

    coef = coefficients[
        (coefficients["spec"] == main_spec)
        & (coefficients["model"] == "ols_raw_hours")
        & (coefficients["term"] == "early_hours")
    ]
    slope = float(coef["Coef."].iloc[0]) if not coef.empty else np.nan
    pval = float(coef["P>|z|"].iloc[0]) if not coef.empty else np.nan

    text = f"""# Labor Adjustment Sequencing: Thesis Interpretation

The empirical exercise compares the early response of hours with the later response
of unemployment across occupation-industry cells. A negative hours response is
interpreted as adjustment on the intensive margin: firms reduce labor input by
cutting hours before reducing headcount. A positive unemployment response at
longer horizons is interpreted as adjustment on the extensive margin.

The baseline specification uses the cumulative hours response over horizons 0-12
and the cumulative unemployment response over horizons 24-36. In this
specification, the Pearson correlation between early hours and later unemployment
is {main_corr['pearson_r']:.3f} (p = {main_corr['pearson_p']:.3f}). Because the
hours measure is signed in levels, a negative correlation means that cells with
more negative early hours responses tend to have larger positive unemployment
responses later. Equivalently, if early hours declines are multiplied by -1, the
sequencing hypothesis predicts a positive association between the intensity of
the early hours decline and subsequent unemployment.

The baseline OLS regression of later unemployment on early hours gives a slope of
{slope:.3f} (heteroskedasticity-robust p = {pval:.3f}). A negative coefficient is
consistent with sequencing: lower early hours are associated with higher later
unemployment. Industry-fixed-effect specifications ask whether the same pattern
holds after comparing occupation cells within broad industries.

For sign consistency, {int(main_sign['sequence_sign_consistent_n'])} out of
{int(main_sign['n_cells'])} cells ({main_sign['sequence_sign_consistent_share']:.1%})
show both an early decline in hours and a later increase in unemployment. Using
pointwise 90 percent confidence intervals, {int(main_sign['sequence_sig90_consistent_n'])}
cells ({main_sign['sequence_sig90_consistent_share']:.1%}) show both a
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
"""
    (output_dir / "thesis_interpretation.md").write_text(text)


def first_horizon_where(df: pd.DataFrame, condition_col: str) -> float:
    hits = df.loc[df[condition_col], "horizon"].sort_values()
    if hits.empty:
        return np.nan
    return float(hits.iloc[0])


def first_horizon_two_consecutive(df: pd.DataFrame, condition_col: str) -> float:
    horizons = set(df.loc[df[condition_col], "horizon"].astype(int))
    for horizon in sorted(horizons):
        if horizon + 1 in horizons:
            return float(horizon)
    return np.nan


def positive_part_sum(df: pd.DataFrame, start: int, end: int, sign: int) -> float:
    sub = df[df["horizon"].between(start, end)]
    if sub.empty:
        return np.nan
    return float(np.maximum(sign * sub["beta"].to_numpy(), 0).sum())


def signed_sum(df: pd.DataFrame, start: int, end: int) -> float:
    sub = df[df["horizon"].between(start, end)]
    if sub.empty:
        return np.nan
    return float(sub["beta"].sum())


def make_ratio(numerator: float, denominator: float) -> float:
    if pd.isna(numerator) or pd.isna(denominator) or np.isclose(denominator, 0):
        return np.nan
    return numerator / denominator


def build_timing_cell_measures(
    panel: pd.DataFrame,
    full_start: int = 0,
    full_end: int = 36,
    early_start: int = 0,
    early_end: int = 12,
    later_start: int = 24,
    later_end: int = 36,
) -> pd.DataFrame:
    timing_panel = panel[panel["outcome"].isin(["hours", "unemployment"])].copy()
    cell_cols = ["industry", "industry_label", "occ_id", "occupation", "cell_id"]

    hours_cells = set(timing_panel.loc[timing_panel["outcome"] == "hours", "cell_id"])
    unemp_cells = set(timing_panel.loc[timing_panel["outcome"] == "unemployment", "cell_id"])
    cells = timing_panel[cell_cols].drop_duplicates("cell_id").set_index("cell_id")

    rows = []
    for cell_id in sorted(hours_cells & unemp_cells):
        hours = timing_panel[
            (timing_panel["cell_id"] == cell_id) & (timing_panel["outcome"] == "hours")
        ].sort_values("horizon")
        unemp = timing_panel[
            (timing_panel["cell_id"] == cell_id) & (timing_panel["outcome"] == "unemployment")
        ].sort_values("horizon")

        hours = hours.assign(sig_neg90=hours["ci_hi90"] < 0)
        unemp = unemp.assign(sig_pos90=unemp["ci_lo90"] > 0)

        hours_early_decline = -signed_sum(hours, early_start, early_end)
        early_unemployment = signed_sum(unemp, early_start, early_end)
        later_unemployment = signed_sum(unemp, later_start, later_end)

        hours_early_decline_pos = positive_part_sum(hours, early_start, early_end, sign=-1)
        hours_full_decline_pos = positive_part_sum(hours, full_start, full_end, sign=-1)
        unemp_early_rise_pos = positive_part_sum(unemp, early_start, early_end, sign=1)
        unemp_full_rise_pos = positive_part_sum(unemp, full_start, full_end, sign=1)

        hours_peak_horizon = float(hours.loc[hours["beta"].idxmin(), "horizon"])
        unemp_peak_horizon = float(unemp.loc[unemp["beta"].idxmax(), "horizon"])

        row = cells.loc[cell_id].to_dict()
        row["cell_id"] = cell_id
        row.update(
            {
                "hours_first_sig_neg90": first_horizon_where(hours, "sig_neg90"),
                "unemp_first_sig_pos90": first_horizon_where(unemp, "sig_pos90"),
                "hours_first_sig_neg90_2consec": first_horizon_two_consecutive(hours, "sig_neg90"),
                "unemp_first_sig_pos90_2consec": first_horizon_two_consecutive(unemp, "sig_pos90"),
                "hours_peak_horizon": hours_peak_horizon,
                "unemp_peak_horizon": unemp_peak_horizon,
                "peak_gap": unemp_peak_horizon - hours_peak_horizon,
                "hours_early_share": make_ratio(hours_early_decline_pos, hours_full_decline_pos),
                "unemp_early_share": make_ratio(unemp_early_rise_pos, unemp_full_rise_pos),
                "early_hours_decline_0_12": hours_early_decline,
                "early_unemployment_0_12": early_unemployment,
                "later_unemployment_24_36": later_unemployment,
                "hours_early_decline_pos_0_12": hours_early_decline_pos,
                "hours_full_decline_pos_0_36": hours_full_decline_pos,
                "unemp_early_rise_pos_0_12": unemp_early_rise_pos,
                "unemp_full_rise_pos_0_36": unemp_full_rise_pos,
            }
        )
        rows.append(row)

    timing = pd.DataFrame(rows)
    if timing.empty:
        return timing

    timing["onset_gap_90"] = timing["unemp_first_sig_pos90"] - timing["hours_first_sig_neg90"]
    timing["onset_gap_90_2consec"] = (
        timing["unemp_first_sig_pos90_2consec"] - timing["hours_first_sig_neg90_2consec"]
    )
    timing["hours_before_unemp_90"] = np.where(
        timing["onset_gap_90"].notna(), timing["onset_gap_90"] > 0, np.nan
    )
    timing["hours_before_unemp_90_2consec"] = np.where(
        timing["onset_gap_90_2consec"].notna(), timing["onset_gap_90_2consec"] > 0, np.nan
    )
    timing["frontload_diff"] = timing["hours_early_share"] - timing["unemp_early_share"]
    timing["hours_peak_before_unemp"] = timing["peak_gap"] > 0

    return timing


def onset_summary_row(timing: pd.DataFrame, label: str, hours_col: str, unemp_col: str) -> dict[str, object]:
    both = timing.dropna(subset=[hours_col, unemp_col]).copy()
    gap = both[unemp_col] - both[hours_col]
    n_both = len(both)

    return {
        "definition": label,
        "n_cells": len(timing),
        "hours_onset_observed_n": int(timing[hours_col].notna().sum()),
        "unemp_onset_observed_n": int(timing[unemp_col].notna().sum()),
        "both_onsets_observed_n": n_both,
        "hours_earlier_n": int((gap > 0).sum()),
        "unemp_earlier_n": int((gap < 0).sum()),
        "same_horizon_n": int((gap == 0).sum()),
        "hours_earlier_share_among_both": (gap > 0).mean() if n_both else np.nan,
        "unemp_earlier_share_among_both": (gap < 0).mean() if n_both else np.nan,
        "same_horizon_share_among_both": (gap == 0).mean() if n_both else np.nan,
        "median_onset_gap": float(gap.median()) if n_both else np.nan,
        "mean_onset_gap": float(gap.mean()) if n_both else np.nan,
    }


def sign_test_greater(values: pd.Series) -> tuple[int, int, float]:
    valid = values.dropna()
    nonzero = valid[~np.isclose(valid, 0)]
    n = len(nonzero)
    if n == 0:
        return 0, 0, np.nan
    successes = int((nonzero > 0).sum())
    return successes, n, stats.binomtest(successes, n=n, p=0.5, alternative="greater").pvalue


def wilcoxon_greater(values: pd.Series) -> float:
    valid = values.dropna()
    nonzero = valid[~np.isclose(valid, 0)]
    if len(nonzero) == 0:
        return np.nan
    try:
        return stats.wilcoxon(nonzero, alternative="greater").pvalue
    except ValueError:
        return np.nan


def make_timing_summaries(timing: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    onset_summary = pd.DataFrame(
        [
            onset_summary_row(
                timing,
                "first_significant_horizon",
                "hours_first_sig_neg90",
                "unemp_first_sig_pos90",
            ),
            onset_summary_row(
                timing,
                "first_two_consecutive_significant_horizons",
                "hours_first_sig_neg90_2consec",
                "unemp_first_sig_pos90_2consec",
            ),
        ]
    )

    rows = []
    peak_gap = timing["peak_gap"].dropna()
    peak_successes, peak_n, peak_sign_p = sign_test_greater(peak_gap)
    rows.append(
        {
            "metric": "peak_gap",
            "definition": "unemp_peak_horizon - hours_peak_horizon",
            "n_valid": len(peak_gap),
            "positive_n": int((peak_gap > 0).sum()),
            "negative_n": int((peak_gap < 0).sum()),
            "zero_n": int((peak_gap == 0).sum()),
            "positive_share": (peak_gap > 0).mean() if len(peak_gap) else np.nan,
            "median": float(peak_gap.median()) if len(peak_gap) else np.nan,
            "mean": float(peak_gap.mean()) if len(peak_gap) else np.nan,
            "sign_test_positive_n": peak_successes,
            "sign_test_nonzero_n": peak_n,
            "sign_test_greater_p": peak_sign_p,
            "wilcoxon_greater_p": wilcoxon_greater(peak_gap),
        }
    )

    frontload_diff = timing["frontload_diff"].dropna()
    frontload_successes, frontload_n, frontload_sign_p = sign_test_greater(frontload_diff)
    rows.append(
        {
            "metric": "frontload_diff",
            "definition": "hours_early_share - unemp_early_share",
            "n_valid": len(frontload_diff),
            "positive_n": int((frontload_diff > 0).sum()),
            "negative_n": int((frontload_diff < 0).sum()),
            "zero_n": int((frontload_diff == 0).sum()),
            "positive_share": (frontload_diff > 0).mean() if len(frontload_diff) else np.nan,
            "median": float(frontload_diff.median()) if len(frontload_diff) else np.nan,
            "mean": float(frontload_diff.mean()) if len(frontload_diff) else np.nan,
            "sign_test_positive_n": frontload_successes,
            "sign_test_nonzero_n": frontload_n,
            "sign_test_greater_p": frontload_sign_p,
            "wilcoxon_greater_p": wilcoxon_greater(frontload_diff),
        }
    )

    return onset_summary, pd.DataFrame(rows)


def run_timing_regressions(timing: pd.DataFrame) -> pd.DataFrame:
    data = timing.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["later_unemployment_24_36", "early_hours_decline_0_12", "early_unemployment_0_12"]
    )
    if len(data) < 3:
        return pd.DataFrame()

    rows = []
    models = {
        "baseline_controls_early_unemp": (
            "later_unemployment_24_36 ~ early_hours_decline_0_12 + early_unemployment_0_12"
        ),
        "industry_fe_controls_early_unemp": (
            "later_unemployment_24_36 ~ early_hours_decline_0_12 + early_unemployment_0_12 + C(industry)"
        ),
    }
    for model_name, formula in models.items():
        model = smf.ols(formula, data=data).fit(cov_type="HC1")
        coef = model.summary2().tables[1].reset_index().rename(columns={"index": "term"})
        coef.insert(0, "model", model_name)
        coef.insert(1, "formula", formula)
        coef.insert(2, "n", int(model.nobs))
        coef.insert(3, "r2", model.rsquared)
        coef.insert(4, "adj_r2", model.rsquared_adj)
        rows.append(coef)

    return pd.concat(rows, ignore_index=True)


def format_count_share(row: pd.Series, count_col: str, denom_col: str) -> str:
    denom = row[denom_col]
    if pd.isna(denom) or denom == 0:
        return "0/0"
    count = int(row[count_col])
    return f"{count}/{int(denom)} ({count / denom:.1%})"


def write_timing_interpretation(
    output_dir: Path,
    onset_summary: pd.DataFrame,
    dynamic_summary: pd.DataFrame,
    regressions: pd.DataFrame,
) -> None:
    one_hit = onset_summary[onset_summary["definition"] == "first_significant_horizon"].iloc[0]
    two_consec = onset_summary[
        onset_summary["definition"] == "first_two_consecutive_significant_horizons"
    ].iloc[0]
    peak = dynamic_summary[dynamic_summary["metric"] == "peak_gap"].iloc[0]
    frontload = dynamic_summary[dynamic_summary["metric"] == "frontload_diff"].iloc[0]

    reg_line = ""
    reg = regressions[
        (regressions["model"] == "baseline_controls_early_unemp")
        & (regressions["term"] == "early_hours_decline_0_12")
    ]
    if not reg.empty:
        reg_line = (
            "In the regression controlling for early unemployment, the coefficient on "
            f"early hours decline is {float(reg['Coef.'].iloc[0]):.3f} "
            f"(robust p = {float(reg['P>|z|'].iloc[0]):.3f})."
        )

    text = f"""# Timing Tests: Hours and Unemployment

These tests ask whether the hours response appears earlier or is more front-loaded
than the unemployment response across occupation-industry cells. They are
descriptive timing evidence. They do not by themselves identify causal firm
margin-switching, because unemployment can also respond early and both margins
may reflect common exposure to the oil-news shock.

## Significance Onset

Using the first 90 percent pointwise-significant horizon, hours becomes
significantly negative before unemployment becomes significantly positive in
{format_count_share(one_hit, 'hours_earlier_n', 'both_onsets_observed_n')} cells
with both onsets observed. The median onset gap, defined as unemployment onset
minus hours onset, is {one_hit['median_onset_gap']:.2f} horizons.

Using the stricter two-consecutive-horizon rule, hours becomes significantly
negative before unemployment becomes significantly positive in
{format_count_share(two_consec, 'hours_earlier_n', 'both_onsets_observed_n')} cells
with both onsets observed. The median onset gap is
{two_consec['median_onset_gap']:.2f} horizons.

## Peak Timing and Front Loading

The peak-timing measure defines peak_gap as the unemployment peak horizon minus
the hours trough horizon. A positive value means the hours trough occurs earlier.
The median peak gap is {peak['median']:.2f} horizons; hours peaks earlier in
{format_count_share(peak, 'positive_n', 'n_valid')} valid cells.

The front-loading measure compares the share of each response that occurs during
horizons 0-12. A positive frontload_diff means the hours decline is more
front-loaded than the unemployment rise. The median frontload_diff is
{frontload['median']:.3f}; it is positive in
{format_count_share(frontload, 'positive_n', 'n_valid')} valid cells
(sign-test p = {frontload['sign_test_greater_p']:.3f};
Wilcoxon p = {frontload['wilcoxon_greater_p']:.3f}).

## Controlling for Early Unemployment

The timing regression tests whether early hours declines predict later
unemployment after controlling for the unemployment response that has already
appeared in horizons 0-12. {reg_line}

The interpretation should therefore be stated carefully: evidence supports an
earlier or more front-loaded hours margin only where these timing measures are
positive. It should not be written as if unemployment is zero in the early
horizons.
"""
    (output_dir / "timing_interpretation.md").write_text(text)


def run_timing_tests(panel: pd.DataFrame, output_dir: Path) -> None:
    timing_dir = output_dir / "timing_tests"
    timing_dir.mkdir(parents=True, exist_ok=True)

    timing = build_timing_cell_measures(panel)
    timing.to_csv(timing_dir / "timing_cell_measures.csv", index=False)

    onset_summary, dynamic_summary = make_timing_summaries(timing)
    regressions = run_timing_regressions(timing)

    onset_summary.to_csv(timing_dir / "timing_onset_summary.csv", index=False)
    dynamic_summary.to_csv(timing_dir / "timing_dynamic_summary.csv", index=False)
    regressions.to_csv(timing_dir / "timing_regressions.csv", index=False)
    write_timing_interpretation(timing_dir, onset_summary, dynamic_summary, regressions)


def run_workflow(input_dir: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    panel = load_irf_panel(input_dir)
    panel.to_csv(output_dir / "irf_panel_long.csv", index=False)
    write_coverage_diagnostics(panel, output_dir)

    sequence = build_sequence_measures(panel, SPECS)
    sequence.to_csv(output_dir / "cell_sequence_measures.csv", index=False)

    corr, sign = make_tables(sequence)
    model_summary, coefficients = run_regressions(sequence)
    save_publication_tables(corr, sign, model_summary, coefficients, output_dir)
    make_plots(sequence, output_dir)
    write_interpretation(output_dir, corr, sign, coefficients)
    run_timing_tests(panel, output_dir)

    spec_doc = pd.DataFrame([spec.__dict__ for spec in SPECS])
    spec_doc.to_csv(output_dir / "sequence_specifications.csv", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("result/occ_ind_4"),
        help="Root directory containing outcome/industry/irf_occ*.csv files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("result/occ_ind_4/labor_adjustment_sequence"),
        help="Directory where tables, figures, and interpretation text are saved.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_workflow(args.input_dir, args.output_dir)


if __name__ == "__main__":
    main()
