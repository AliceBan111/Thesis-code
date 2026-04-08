import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

# ── paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR  = os.getcwd()
IRF_PATH    = os.path.join(SCRIPT_DIR, "../../result/macro/indpro/lp_indpro_irf.csv")
OUTPUT_PATH = os.path.join(SCRIPT_DIR, "../../result/macro/indpro/lp_indpro_irf.pdf")

# ── style ─────────────────────────────────────────────────────────────────────
PLT_COLOR  = "#185FA5"   # blue-600
CI95_COLOR = "#B5D4F4"   # blue-100
CI90_COLOR = "#85B7EB"   # blue-200
CI68_COLOR = "#378ADD"   # blue-400

def plot_irf(irf_path: str, output_path: str) -> None:
    df = pd.read_csv(irf_path)
    df = df.sort_values("horizon").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8, 4.2), dpi=150)

    h = df["horizon"].values

    # shaded confidence bands (widest first so narrower bands render on top)
    ax.fill_between(h, df["ci_lo95"], df["ci_hi95"],
                    color=CI95_COLOR, alpha=0.55, linewidth=0, label="95% CI")
    ax.fill_between(h, df["ci_lo90"], df["ci_hi90"],
                    color=CI90_COLOR, alpha=0.60, linewidth=0, label="90% CI")
    ax.fill_between(h, df["ci_lo68"], df["ci_hi68"],
                    color=CI68_COLOR, alpha=0.55, linewidth=0, label="68% CI")

    # point estimate
    ax.plot(h, df["beta"], color=PLT_COLOR, linewidth=1.8, label="Point estimate", zorder=5)

    # zero line
    ax.axhline(0, color="#444441", linewidth=0.7, linestyle="--", alpha=0.6)

    # axes
    ax.set_xlabel("Horizon (months)", fontsize=11)
    ax.set_ylabel("Response of log Industrial Production\n(cumulative from t−1)", fontsize=11)
    ax.set_title("IRF: Oil Supply News Shock → Industrial Production", fontsize=13, fontweight="medium")
    ax.set_xlim(h[0], h[-1])

    # y-axis: convert log difference to approximate % change label
    yticks = ax.get_yticks()
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{v*100:.1f}%" for v in yticks], fontsize=9)
    ax.set_xticks(np.arange(0, h[-1]+1, 6))

    ax.tick_params(axis="both", labelsize=10)
    ax.spines[["top","right"]].set_visible(False)
    ax.spines[["left","bottom"]].set_linewidth(0.6)
    ax.grid(axis="y", linewidth=0.4, linestyle=":", color="#D3D1C7", alpha=0.8)

    # legend
    handles = [
        plt.Line2D([0],[0], color=PLT_COLOR, linewidth=1.8),
        mpatches.Patch(color=CI68_COLOR, alpha=0.55),
        mpatches.Patch(color=CI90_COLOR, alpha=0.60),
        mpatches.Patch(color=CI95_COLOR, alpha=0.55),
    ]
    labels = ["Point estimate", "68% CI", "90% CI", "95% CI"]
    ax.legend(handles, labels, fontsize=9, frameon=False, loc="upper left",
              ncol=2, handlelength=1.2, columnspacing=1.0)

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    print(f"Saved: {output_path}")
    plt.show()


if __name__ == "__main__":
    plot_irf(IRF_PATH, OUTPUT_PATH)