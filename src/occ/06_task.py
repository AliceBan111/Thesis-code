import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf


# =============================================================================
# Occupation-level task correlation plots
#
# Regression direction:
#
#   Cumulative IRF beta_o = alpha + beta * Task_o + error_o
#
# Plot:
#   x-axis = task measure
#   y-axis = cumulative IRF beta
#
# Inputs:
#   ../../result/occ/analysis/merged_occ_irf_trajectories.csv
#   ../../data/occ1990dd_task_alm.dta
#   ../../result/mapping/mapping_done.xlsx
#
# Outputs:
#   ../../result/occ/correlation/task_python/{outcome}_all_tasks.pdf
#   ../../result/occ/correlation/task_python/{outcome}_all_tasks.png
# =============================================================================


# =============================================================================
# 1. SETTINGS
# =============================================================================

IRF_FILE = "../../result/occ/analysis/merged_occ_irf_trajectories.csv"
DTA_PATH = "../../data/occ1990dd_task_alm.dta"
XLSX_PATH = "../../result/mapping/mapping_done.xlsx"

OUTPUT_DIR = "../../result/occ/correlation/task_python"
os.makedirs(OUTPUT_DIR, exist_ok=True)

OUTCOME_LIST = [
    "employment",
    "hourly_rate",
    "income_share",
    "inequality",
    "median",
    "unemployment",
    "income",
    "hours",
]

TASK_LIST = [
    ("task_abstract", "Abstract Task Weight"),
    ("task_routine", "Routine Task Weight"),
    ("task_manual", "Manual Task Weight"),
]

GROUP_MAP = {
    1: "group1_Managerial",
    2: "group2_Professional_specialty",
    3: "group3_High_tech",
    4: "group4_Sales",
    5: "group5_Administrative_support",
    6: "group6_Service",
    7: "group7_Farming_forestry_construction",
    8: "group8_Precision_production_repair",
    9: "group9_Machine_operators_transport",
}

PLOT_NAME_MAP = {
    "group1_Managerial": "Managerial",
    "group2_Professional_specialty": "Prof. specialty",
    "group3_High_tech": "High Tech",
    "group4_Sales": "Sales",
    "group5_Administrative_support": "Admin",
    "group6_Service": "Service",
    "group7_Farming_forestry_construction": "Constr., Extract., Farm",
    "group8_Precision_production_repair": "Prod, Repair",
    "group9_Machine_operators_transport": "Machinists, Transp.",
}


# =============================================================================
# 2. CALCULATE OCCUPATION-LEVEL CUMULATIVE IRF
# =============================================================================

def calculate_cumulative_irf_single(
    input_file: str,
    target_outcome: str,
    standardize: bool = True,
    h_min: int = 1,
    h_max: int = 36,
) -> pd.DataFrame:
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"File not found: {input_file}")

    df = pd.read_csv(input_file, na_values=["", "NA", "N/A"])

    if "outcome" not in df.columns or "horizon" not in df.columns:
        raise ValueError("IRF file must contain columns: outcome, horizon")

    df = df[df["outcome"] == target_outcome].copy()

    if df.empty:
        raise ValueError(f"Outcome '{target_outcome}' not found in {input_file}")

    group_cols = [c for c in df.columns if c not in ["outcome", "horizon"]]

    if not group_cols:
        raise ValueError(f"No occupation group columns found in {input_file}")

    df_long = df.melt(
        id_vars=["outcome", "horizon"],
        value_vars=group_cols,
        var_name="Occupational_Group",
        value_name="Beta",
    )

    df_long["horizon"] = pd.to_numeric(df_long["horizon"], errors="coerce")
    df_long["Beta"] = pd.to_numeric(df_long["Beta"], errors="coerce")

    df_long = df_long.replace([np.inf, -np.inf], np.nan)
    df_long = df_long.dropna(subset=["horizon", "Beta"])

    df_long = df_long[
        (df_long["horizon"] >= h_min)
        & (df_long["horizon"] <= h_max)
    ].copy()

    if df_long.empty:
        raise ValueError(f"No valid beta observations for outcome {target_outcome}")

    df_cum = (
        df_long.groupby("Occupational_Group", as_index=False)["Beta"]
        .sum()
        .rename(columns={"Beta": target_outcome})
    )

    if standardize:
        vals = df_cum[target_outcome].astype(float)
        mu = vals.mean()
        sigma = vals.std(ddof=1)

        if sigma == 0 or not np.isfinite(sigma):
            df_cum[target_outcome] = 0.0
        else:
            df_cum[target_outcome] = (vals - mu) / sigma

    df_cum = df_cum.sort_values("Occupational_Group").reset_index(drop=True)

    return df_cum


# =============================================================================
# 3. PREPARE DORN / AUTOR TASK DATA
# =============================================================================

def prepare_task_data(
    dta_path: str = DTA_PATH,
    xlsx_path: str = XLSX_PATH,
) -> pd.DataFrame:
    if not os.path.exists(dta_path):
        raise FileNotFoundError(f"Task DTA file not found: {dta_path}")

    if not os.path.exists(xlsx_path):
        raise FileNotFoundError(f"Mapping Excel file not found: {xlsx_path}")

    df_tasks = pd.read_stata(dta_path)

    task_cols = ["occ1990dd", "task_abstract", "task_routine", "task_manual"]
    for col in task_cols:
        if col not in df_tasks.columns:
            raise ValueError(f"{col} not found in task data")

    df_tasks = df_tasks[task_cols].copy()

    df_mapping = pd.read_excel(xlsx_path, sheet_name="Sheet1")

    # Julia used column positions:
    # B column (2) -> occ1990dd
    # F column (6) -> group
    # H column (8) -> weights
    #
    # Python is 0-indexed:
    # B -> index 1
    # F -> index 5
    # H -> index 7
    if df_mapping.shape[1] < 8:
        raise ValueError("Mapping Excel must have at least 8 columns.")

    df_mapping = df_mapping.iloc[:, [1, 5, 7]].copy()
    df_mapping.columns = ["occ1990dd", "group", "weights"]

    df_mapping["group"] = pd.to_numeric(df_mapping["group"], errors="coerce")
    df_mapping["weights"] = pd.to_numeric(df_mapping["weights"], errors="coerce")
    df_mapping["occ1990dd"] = pd.to_numeric(df_mapping["occ1990dd"], errors="coerce")

    df_mapping = df_mapping.dropna(subset=["occ1990dd", "group", "weights"]).copy()
    df_mapping["group"] = df_mapping["group"].astype(int)

    df_tasks["occ1990dd"] = pd.to_numeric(df_tasks["occ1990dd"], errors="coerce")

    df_merged = pd.merge(df_tasks, df_mapping, on="occ1990dd", how="inner")

    df_merged = df_merged.replace([np.inf, -np.inf], np.nan)
    df_merged = df_merged.dropna(
        subset=["task_abstract", "task_routine", "task_manual", "weights", "group"]
    )

    def weighted_mean(g: pd.DataFrame, col: str) -> float:
        w = g["weights"].to_numpy(dtype=float)
        x = g[col].to_numpy(dtype=float)

        if np.sum(w) == 0:
            return np.nan

        return float(np.sum(x * w) / np.sum(w))

    rows = []
    for group, g in df_merged.groupby("group"):
        rows.append(
            {
                "group": int(group),
                "task_abstract": weighted_mean(g, "task_abstract"),
                "task_routine": weighted_mean(g, "task_routine"),
                "task_manual": weighted_mean(g, "task_manual"),
            }
        )

    df_y = pd.DataFrame(rows)

    df_y["Occupational_Group"] = df_y["group"].map(GROUP_MAP)

    df_y = df_y.dropna(subset=["Occupational_Group"]).copy()

    df_y = df_y[
        [
            "Occupational_Group",
            "task_abstract",
            "task_routine",
            "task_manual",
        ]
    ].copy()

    return df_y


# =============================================================================
# 4. PLOT FUNCTION
# =============================================================================

def plot_correlation_with_ci(
    df_merged: pd.DataFrame,
    x_col: str,
    y_col: str,
    x_label: str,
    y_label: str,
    ax,
    show_labels: bool = True,
):
    df_plot = df_merged.copy()

    df_plot["Plot_Label"] = df_plot["Occupational_Group"].map(PLOT_NAME_MAP)
    df_plot["Plot_Label"] = df_plot["Plot_Label"].fillna(df_plot["Occupational_Group"])

    df_plot[x_col] = pd.to_numeric(df_plot[x_col], errors="coerce")
    df_plot[y_col] = pd.to_numeric(df_plot[y_col], errors="coerce")

    df_plot = df_plot.replace([np.inf, -np.inf], np.nan)
    df_plot = df_plot.dropna(subset=[x_col, y_col])

    if len(df_plot) < 3:
        raise ValueError(f"Too few observations for regression: {len(df_plot)}")

    # OLS: y ~ x
    model = smf.ols(f"{y_col} ~ {x_col}", data=df_plot).fit()

    x_min, x_max = df_plot[x_col].min(), df_plot[x_col].max()
    y_min, y_max = df_plot[y_col].min(), df_plot[y_col].max()

    x_span = x_max - x_min
    y_span = y_max - y_min

    if np.isclose(x_span, 0):
        x_span = 1.0

    if np.isclose(y_span, 0):
        y_span = 1.0

    x_pad = x_span * 0.12
    y_pad = y_span * 0.12

    x_grid = np.linspace(x_min - x_pad, x_max + x_pad, 150)
    pred_df = pd.DataFrame({x_col: x_grid})
    pred = model.get_prediction(pred_df).summary_frame(alpha=0.10)

    ax.plot(x_grid, pred["mean"], linewidth=2)
    ax.fill_between(
        x_grid,
        pred["mean_ci_lower"],
        pred["mean_ci_upper"],
        alpha=0.20,
    )

    ax.scatter(
        df_plot[x_col],
        df_plot[y_col],
        s=45,
        facecolors="white",
        edgecolors="red",
        linewidths=1.4,
        zorder=3,
    )

    if show_labels:
        for _, row in df_plot.iterrows():
            ax.text(
                row[x_col] + 0.015 * x_span,
                row[y_col] + 0.020 * y_span,
                row["Plot_Label"],
                fontsize=8,
                ha="left",
                va="bottom",
            )

    coef = model.params[x_col]
    pval = model.pvalues[x_col]
    pstr = "<0.001" if pval < 0.001 else f"{pval:.3f}"

    ax.text(
        0.02,
        0.95,
        f"coef = {coef:.3f}, p = {pstr}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    ax.set_xlim(x_min - x_pad, x_max + x_pad * 2.2)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)

    ax.grid(False)

    return model


# =============================================================================
# 5. MAIN
# =============================================================================

def main():
    print("Preparing task data...")
    df_tasks = prepare_task_data()

    summary_rows = []

    for outcome in OUTCOME_LIST:
        print(f"\n--- Processing outcome: {outcome} ---")

        try:
            df_x = calculate_cumulative_irf_single(
                IRF_FILE,
                outcome,
                standardize=True,
                h_min=1,
                h_max=36,
            )
        except Exception as e:
            print(f"[skip] {outcome}: {e}")
            continue

        df_merged = pd.merge(df_x, df_tasks, on="Occupational_Group", how="inner")

        print(f"Number of occupation groups: {len(df_merged)}")

        if len(df_merged) < 3:
            print(f"[skip] too few occupation groups for {outcome}")
            continue

        fig, axes = plt.subplots(
            3,
            1,
            figsize=(7.5, 12),
            constrained_layout=True,
        )

        for ax, (task_col, task_label) in zip(axes, TASK_LIST):
            model = plot_correlation_with_ci(
                df_merged,
                x_col=task_col,
                y_col=outcome,
                x_label=task_label,
                y_label=f"Cumulative IRF Beta ({outcome})",
                ax=ax,
                show_labels=True,
            )

            summary_rows.append(
                {
                    "outcome": outcome,
                    "task": task_col,
                    "estimate": model.params[task_col],
                    "std_error": model.bse[task_col],
                    "t_stat": model.tvalues[task_col],
                    "p_value": model.pvalues[task_col],
                    "r2": model.rsquared,
                    "adj_r2": model.rsquared_adj,
                    "nobs": model.nobs,
                }
            )

        fig.suptitle(
            f"Task Characteristics vs Cumulative IRF: {outcome}",
            y=1.02,
            fontsize=12,
        )

        pdf_path = os.path.join(OUTPUT_DIR, f"{outcome}_all_tasks.pdf")
        png_path = os.path.join(OUTPUT_DIR, f"{outcome}_all_tasks.png")

        fig.savefig(pdf_path, bbox_inches="tight")
        fig.savefig(png_path, dpi=300, bbox_inches="tight")

        plt.close(fig)

        print(f"Saved: {pdf_path}")

    summary = pd.DataFrame(summary_rows)
    summary_path = os.path.join(OUTPUT_DIR, "occupation_task_correlation_summary.csv")
    summary.to_csv(summary_path, index=False)

    print("\nSaved summary:")
    print(summary_path)
    print(summary)


if __name__ == "__main__":
    main()