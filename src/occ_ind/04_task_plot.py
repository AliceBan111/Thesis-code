import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf


# =============================================================================
# SETTINGS
# =============================================================================

BASE_DIR = "../../result/occ_ind_4/task_regression"
OUTPUT_DIR = os.path.join(BASE_DIR, "python_plots")
os.makedirs(OUTPUT_DIR, exist_ok=True)

OUTCOME_LIST = [
    "employment",
    "hourly_rate",
    "hours",
    "income_share_var",
    "income",
    "inequality",
    "median",
    "unemployment",
]

TASKS = [
    ("task_abstract", "Abstract task"),
    ("task_routine", "Routine task"),
    ("task_manual", "Manual task"),
]

INDUSTRY_LABELS = {
    "Energy_intensive": "Energy intensive",
    "Manufacturing_Construction": "Mfg. + Construction",
    "Services": "Services",
    "Trade": "Trade",
}

INDUSTRY_MARKERS = {
    "Energy_intensive": "o",
    "Manufacturing_Construction": "s",
    "Services": "D",
    "Trade": "^",
}

OCC_LABELS = {
    1: "Managerial",
    2: "Prof.",
    3: "High Tech",
    4: "Sales",
    5: "Admin",
    6: "Service",
    7: "Constr./Farm",
    8: "Prod./Repair",
    9: "Mach./Transp.",
}


# =============================================================================
# DATA CLEANING
# =============================================================================

def load_outcome_data(outcome: str) -> pd.DataFrame:
    path = os.path.join(BASE_DIR, f"{outcome}_within_industry_task_data.csv")

    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing file: {path}")

    df = pd.read_csv(path)

    required = [
        "industry",
        "occ_id",
        "cumulative_beta",
        "task_abstract",
        "task_routine",
        "task_manual",
    ]

    for col in required:
        if col not in df.columns:
            raise ValueError(f"{col} not found in {path}")

    for col in ["cumulative_beta", "task_abstract", "task_routine", "task_manual"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["occ_id"] = pd.to_numeric(df["occ_id"], errors="coerce").astype("Int64")

    if "Plot_Label" not in df.columns:
        df["Plot_Label"] = df["occ_id"].map(OCC_LABELS)
    else:
        df["Plot_Label"] = df["Plot_Label"].fillna(df["occ_id"].map(OCC_LABELS))

    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(
        subset=[
            "industry",
            "occ_id",
            "cumulative_beta",
            "task_abstract",
            "task_routine",
            "task_manual",
        ]
    )

    df["industry"] = df["industry"].astype(str)

    return df


# =============================================================================
# REGRESSION / PARTIALIZATION
# =============================================================================

def run_main_model(df: pd.DataFrame):
    formula = (
        "cumulative_beta ~ task_abstract + task_routine "
        "+ task_manual + C(industry)"
    )
    return smf.ols(formula, data=df).fit()


def residualize(df: pd.DataFrame, lhs: str, rhs_terms: list[str]) -> np.ndarray:
    formula = lhs + " ~ " + " + ".join(rhs_terms)
    model = smf.ols(formula, data=df).fit()
    return model.resid.values


def make_partial_data(df: pd.DataFrame, task_col: str) -> pd.DataFrame:
    other_tasks = [t for t, _ in TASKS if t != task_col]
    controls = other_tasks + ["C(industry)"]

    out = df.copy()
    out["partial_beta"] = residualize(out, "cumulative_beta", controls)
    out["partial_task"] = residualize(out, task_col, controls)

    return out


def fit_line_with_ci(x: np.ndarray, y: np.ndarray):
    tmp = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()

    model = smf.ols("y ~ x", data=tmp).fit()

    x_min, x_max = tmp["x"].min(), tmp["x"].max()

    if np.isclose(x_min, x_max):
        x_min -= 0.5
        x_max += 0.5

    x_grid = np.linspace(x_min, x_max, 150)
    pred_df = pd.DataFrame({"x": x_grid})
    pred = model.get_prediction(pred_df).summary_frame(alpha=0.10)

    return model, x_grid, pred


# =============================================================================
# PLOTTING
# =============================================================================

def add_industry_points(ax, df: pd.DataFrame, x_col: str, y_col: str):
    for industry, sub in df.groupby("industry"):
        ax.scatter(
            sub[x_col],
            sub[y_col],
            marker=INDUSTRY_MARKERS.get(industry, "o"),
            s=55,
            alpha=0.85,
            label=INDUSTRY_LABELS.get(industry, industry),
        )


def add_selected_labels(ax, df: pd.DataFrame, x_col: str, y_col: str):
    """
    To avoid clutter, label only a few interpretable occupation groups.
    Change selected_occ_ids if needed.
    """
    selected_occ_ids = {1, 3, 4, 5, 8, 9}

    x_span = df[x_col].max() - df[x_col].min()
    y_span = df[y_col].max() - df[y_col].min()

    if np.isclose(x_span, 0):
        x_span = 1.0
    if np.isclose(y_span, 0):
        y_span = 1.0

    for _, row in df.iterrows():
        occ_id = int(row["occ_id"])
        if occ_id not in selected_occ_ids:
            continue

        label = row.get("Plot_Label", OCC_LABELS.get(occ_id, str(occ_id)))

        ax.text(
            row[x_col] + 0.012 * x_span,
            row[y_col] + 0.012 * y_span,
            label,
            fontsize=7,
            alpha=0.85,
        )


def plot_raw_three_tasks(df: pd.DataFrame, outcome: str):
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(8.5, 12),
        constrained_layout=True,
    )

    for ax, (task_col, task_label) in zip(axes, TASKS):
        add_industry_points(ax, df, task_col, "cumulative_beta")
        add_selected_labels(ax, df, task_col, "cumulative_beta")

        ax.axhline(0, linewidth=0.8)
        ax.set_title(task_label)
        ax.set_xlabel(task_label)
        ax.set_ylabel("Cumulative IRF beta")
        ax.grid(True, linestyle=":", linewidth=0.6, alpha=0.7)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 1.02),
    )

    fig.suptitle(f"Raw task relationships by industry: {outcome}", y=1.06)

    pdf_path = os.path.join(OUTPUT_DIR, f"{outcome}_raw_three_tasks_by_industry.pdf")
    png_path = os.path.join(OUTPUT_DIR, f"{outcome}_raw_three_tasks_by_industry.png")

    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return pdf_path


def plot_partial_three_tasks(df: pd.DataFrame, outcome: str):
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(8.5, 12),
        constrained_layout=True,
    )

    main_model = run_main_model(df)

    for ax, (task_col, task_label) in zip(axes, TASKS):
        pdf = make_partial_data(df, task_col)

        line_model, x_grid, pred = fit_line_with_ci(
            pdf["partial_task"].values,
            pdf["partial_beta"].values,
        )

        ax.plot(x_grid, pred["mean"].values, linewidth=2)
        ax.fill_between(
            x_grid,
            pred["mean_ci_lower"].values,
            pred["mean_ci_upper"].values,
            alpha=0.20,
        )

        add_industry_points(ax, pdf, "partial_task", "partial_beta")
        add_selected_labels(ax, pdf, "partial_task", "partial_beta")

        coef = main_model.params.get(task_col, np.nan)
        pval = main_model.pvalues.get(task_col, np.nan)

        if np.isfinite(pval):
            pstr = "<0.001" if pval < 0.001 else f"{pval:.3f}"
        else:
            pstr = "NA"

        ax.text(
            0.02,
            0.95,
            f"coef = {coef:.4f}, p = {pstr}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
        )

        ax.axhline(0, linewidth=0.8)
        ax.axvline(0, linewidth=0.8)
        ax.set_title(task_label)
        ax.set_xlabel(f"Partial {task_label}")
        ax.set_ylabel("Partial cumulative beta")
        ax.grid(True, linestyle=":", linewidth=0.6, alpha=0.7)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 1.02),
    )

    fig.suptitle(
        f"Partial task relationships controlling for industry FE: {outcome}",
        y=1.06,
    )

    pdf_path = os.path.join(OUTPUT_DIR, f"{outcome}_partial_three_tasks_industry_FE.pdf")
    png_path = os.path.join(OUTPUT_DIR, f"{outcome}_partial_three_tasks_industry_FE.png")

    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return pdf_path


# =============================================================================
# MAIN
# =============================================================================

def main():
    summary_rows = []

    for outcome in OUTCOME_LIST:
        print(f"\nProcessing outcome: {outcome}")

        try:
            df = load_outcome_data(outcome)
        except Exception as e:
            print(f"[skip] {outcome}: {e}")
            continue

        print(f"N = {len(df)}")

        if len(df) < 10:
            print(f"[skip] too few observations: {outcome}")
            continue

        model = run_main_model(df)

        for term in ["task_abstract", "task_routine", "task_manual"]:
            summary_rows.append(
                {
                    "outcome": outcome,
                    "term": term,
                    "estimate": model.params.get(term, np.nan),
                    "std_error": model.bse.get(term, np.nan),
                    "t_stat": model.tvalues.get(term, np.nan),
                    "p_value": model.pvalues.get(term, np.nan),
                    "r2": model.rsquared,
                    "adj_r2": model.rsquared_adj,
                    "nobs": model.nobs,
                }
            )

        # Save full statsmodels text summary
        summary_txt = os.path.join(OUTPUT_DIR, f"{outcome}_statsmodels_summary.txt")
        with open(summary_txt, "w", encoding="utf-8") as f:
            f.write(str(model.summary()))

        try:
            raw_path = plot_raw_three_tasks(df, outcome)
            partial_path = plot_partial_three_tasks(df, outcome)
            print(f"Saved raw plot:     {raw_path}")
            print(f"Saved partial plot: {partial_path}")
        except Exception as e:
            print(f"[plot failed] {outcome}: {e}")

    summary = pd.DataFrame(summary_rows)
    summary_path = os.path.join(
        OUTPUT_DIR,
        "summary_three_tasks_with_industry_FE_python.csv",
    )
    summary.to_csv(summary_path, index=False)

    print("\nSaved Python summary:")
    print(summary_path)
    print(summary)


if __name__ == "__main__":
    main()