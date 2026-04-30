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

def run_single_task_model(df: pd.DataFrame, task_col: str):
    formula = f"cumulative_beta ~ {task_col} + C(industry)"
    return smf.ols(formula, data=df).fit()


def residualize(df: pd.DataFrame, lhs: str, rhs_terms: list[str]) -> np.ndarray:
    formula = lhs + " ~ " + " + ".join(rhs_terms)
    model = smf.ols(formula, data=df).fit()
    return model.resid.values


def make_partial_data(df: pd.DataFrame, task_col: str) -> pd.DataFrame:
    controls = ["C(industry)"]

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

FIGSIZE = (10.5, 13.5)
POINT_SIZE = 58
LABEL_FONTSIZE = 6.5
TITLE_FONTSIZE = 11
AXIS_LABEL_FONTSIZE = 9
LEGEND_FONTSIZE = 9

# Whether to label every occupation-industry point.
# If False, only selected occupation groups will be labelled.
LABEL_ALL_POINTS = True

SELECTED_OCC_IDS = {1, 2, 3, 4, 5, 6, 7, 8, 9}


def add_industry_points(ax, df: pd.DataFrame, x_col: str, y_col: str):
    """
    Scatter points by industry. Marker identifies industry.
    """
    for industry, sub in df.groupby("industry"):
        ax.scatter(
            sub[x_col],
            sub[y_col],
            marker=INDUSTRY_MARKERS.get(industry, "o"),
            s=POINT_SIZE,
            alpha=0.85,
            edgecolor="black",
            linewidth=0.35,
            label=INDUSTRY_LABELS.get(industry, industry),
            zorder=3,
        )


def get_label_offset(occ_id: int, industry: str):
    """
    Deterministic label offsets in screen points.

    Purpose:
    - keep labels away from markers;
    - separate labels belonging to different industries;
    - reduce overlap for observations with similar x-y coordinates.
    """

    base_offsets = {
        "Energy_intensive": (7, 7),
        "Manufacturing_Construction": (-7, 7),
        "Services": (7, -9),
        "Trade": (-7, -9),
    }

    dx, dy = base_offsets.get(industry, (7, 7))

    # Small occupation-specific jitter to reduce label overlap.
    dx += ((occ_id % 3) - 1) * 4
    dy += (((occ_id + 1) % 3) - 1) * 3

    return dx, dy


def add_occupation_labels(ax, df: pd.DataFrame, x_col: str, y_col: str):
    """
    Label occupation points.

    This version labels all occupations by default.
    Labels are offset by industry and occupation ID to reduce overlap.
    """

    for _, row in df.iterrows():
        occ_id = int(row["occ_id"])

        if (not LABEL_ALL_POINTS) and (occ_id not in SELECTED_OCC_IDS):
            continue

        industry = str(row["industry"])
        label = row.get("Plot_Label", OCC_LABELS.get(occ_id, str(occ_id)))

        dx, dy = get_label_offset(occ_id, industry)

        ha = "left" if dx >= 0 else "right"
        va = "bottom" if dy >= 0 else "top"

        ax.annotate(
            label,
            xy=(row[x_col], row[y_col]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=LABEL_FONTSIZE,
            alpha=0.90,
            ha=ha,
            va=va,
            zorder=4,
            clip_on=True,
            bbox=dict(
                boxstyle="round,pad=0.12",
                facecolor="white",
                edgecolor="none",
                alpha=0.55,
            ),
        )


def add_axis_padding(ax, df: pd.DataFrame, x_col: str, y_col: str):
    """
    Add explicit x/y padding so labels near plot boundaries are less likely clipped.
    """

    x_min = df[x_col].min()
    x_max = df[x_col].max()
    y_min = df[y_col].min()
    y_max = df[y_col].max()

    x_span = x_max - x_min
    y_span = y_max - y_min

    if np.isclose(x_span, 0):
        x_span = 1.0
    if np.isclose(y_span, 0):
        y_span = 1.0

    ax.set_xlim(x_min - 0.12 * x_span, x_max + 0.16 * x_span)
    ax.set_ylim(y_min - 0.18 * y_span, y_max + 0.18 * y_span)


def get_unique_legend_handles_labels(ax):
    """
    Remove duplicate legend entries.
    """

    handles, labels = ax.get_legend_handles_labels()

    unique = {}
    for h, l in zip(handles, labels):
        if l not in unique:
            unique[l] = h

    return list(unique.values()), list(unique.keys())


def format_axis(ax):
    ax.axhline(0, linewidth=0.8, color="black", alpha=0.75, zorder=1)
    ax.grid(True, linestyle=":", linewidth=0.6, alpha=0.7, zorder=0)
    ax.tick_params(axis="both", labelsize=8)


def add_figure_legend(fig, ax):
    """
    Put legend below the whole figure, not above the first subplot.
    This avoids covering subplot titles.
    """

    handles, labels = get_unique_legend_handles_labels(ax)

    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=4,
        frameon=False,
        fontsize=LEGEND_FONTSIZE,
        bbox_to_anchor=(0.5, 0.015),
    )


def plot_raw_three_tasks(df: pd.DataFrame, outcome: str):
    fig, axes = plt.subplots(
        3,
        1,
        figsize=FIGSIZE,
        constrained_layout=False,
    )

    for ax, (task_col, task_label) in zip(axes, TASKS):
        add_industry_points(ax, df, task_col, "cumulative_beta")
        add_occupation_labels(ax, df, task_col, "cumulative_beta")
        add_axis_padding(ax, df, task_col, "cumulative_beta")

        format_axis(ax)

        ax.set_title(task_label, fontsize=TITLE_FONTSIZE, pad=10)
        ax.set_xlabel(task_label, fontsize=AXIS_LABEL_FONTSIZE)
        ax.set_ylabel("Cumulative IRF beta", fontsize=AXIS_LABEL_FONTSIZE)

    fig.suptitle(
        f"Raw task relationships by industry: {outcome}",
        fontsize=13,
        y=0.975,
    )

    add_figure_legend(fig, axes[0])

    fig.subplots_adjust(
        top=0.93,
        bottom=0.08,
        left=0.10,
        right=0.97,
        hspace=0.48,
    )

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
        figsize=FIGSIZE,
        constrained_layout=False,
    )

    for ax, (task_col, task_label) in zip(axes, TASKS):
        pdf = make_partial_data(df, task_col)

        line_model, x_grid, pred = fit_line_with_ci(
            pdf["partial_task"].values,
            pdf["partial_beta"].values,
        )

        ax.plot(
            x_grid,
            pred["mean"].values,
            linewidth=2,
            zorder=2,
        )

        ax.fill_between(
            x_grid,
            pred["mean_ci_lower"].values,
            pred["mean_ci_upper"].values,
            alpha=0.20,
            zorder=1,
        )

        add_industry_points(ax, pdf, "partial_task", "partial_beta")
        add_occupation_labels(ax, pdf, "partial_task", "partial_beta")
        add_axis_padding(ax, pdf, "partial_task", "partial_beta")

        model = run_single_task_model(df, task_col)

        coef = model.params.get(task_col, np.nan)
        pval = model.pvalues.get(task_col, np.nan)

        if np.isfinite(pval):
            pstr = "<0.001" if pval < 0.001 else f"{pval:.3f}"
        else:
            pstr = "NA"

        ax.text(
            0.015,
            0.965,
            f"coef = {coef:.4f}, p = {pstr}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox=dict(
                boxstyle="round,pad=0.25",
                facecolor="white",
                edgecolor="none",
                alpha=0.75,
            ),
            zorder=5,
        )

        ax.axvline(0, linewidth=0.8, color="black", alpha=0.75, zorder=1)
        format_axis(ax)

        ax.set_title(task_label, fontsize=TITLE_FONTSIZE, pad=10)
        ax.set_xlabel(
            f"Within-industry residualized {task_label}",
            fontsize=AXIS_LABEL_FONTSIZE,
        )
        ax.set_ylabel(
            "Within-industry residualized cumulative beta",
            fontsize=AXIS_LABEL_FONTSIZE,
        )

    fig.suptitle(
        f"Separate task relationships controlling for industry FE: {outcome}",
        fontsize=13,
        y=0.975,
    )

    add_figure_legend(fig, axes[0])

    fig.subplots_adjust(
        top=0.93,
        bottom=0.08,
        left=0.10,
        right=0.97,
        hspace=0.48,
    )

    pdf_path = os.path.join(
        OUTPUT_DIR,
        f"{outcome}_partial_separate_tasks_industry_FE.pdf",
    )
    png_path = os.path.join(
        OUTPUT_DIR,
        f"{outcome}_partial_separate_tasks_industry_FE.png",
    )

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

        for task_col, task_label in TASKS:
            model = run_single_task_model(df, task_col)

            summary_rows.append(
                {
                    "outcome": outcome,
                    "model": "single_task_industry_FE",
                    "task": task_label,
                    "term": task_col,
                    "estimate": model.params.get(task_col, np.nan),
                    "std_error": model.bse.get(task_col, np.nan),
                    "t_stat": model.tvalues.get(task_col, np.nan),
                    "p_value": model.pvalues.get(task_col, np.nan),
                    "r2": model.rsquared,
                    "adj_r2": model.rsquared_adj,
                    "nobs": model.nobs,
                }
            )

            summary_txt = os.path.join(
                OUTPUT_DIR,
                f"{outcome}_{task_col}_single_task_industry_FE_statsmodels_summary.txt"
            )               

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