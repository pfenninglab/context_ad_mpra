#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analyze Nguyen cortical controls in the Pseudobarcode R1R2 MAD-score table.

This script:
1. Loads 20260625_MAD_Neuron_Pseudobarcode_control_R1R2.csv
2. Defines:
   - Negative control: control == TRUE or ID contains Nguyen_neg
   - Nguyen cortical: ID contains Nguyen_cortical
   - Mouse Nguyen: Nguyen_cortical IDs containing _mm_
   - Human Nguyen: Nguyen_cortical IDs containing _hs_
3. Performs one-sided Mann-Whitney U tests:
   - All Nguyen cortical > Negative control
   - Mouse Nguyen cortical > Negative control
4. Applies Benjamini-Hochberg FDR correction across these two planned tests
5. Saves:
   - test summary CSV
   - plot data CSV
   - PNG/SVG/PDF figure with human blue and mouse orange
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


# -----------------------------------------------------------------------------
# Input / output
# -----------------------------------------------------------------------------

input_csv = "20260625_MAD_Neuron_Pseudobarcode_control_R1R2.csv"
out_prefix = "mad_pseudobarcode_R1R2_nguyen_one_sided_fdr"


# -----------------------------------------------------------------------------
# Plot settings for Illustrator-editable vector output
# -----------------------------------------------------------------------------

plt.rcParams["svg.fonttype"] = "none"   # keep SVG text editable
plt.rcParams["pdf.fonttype"] = 42       # embed TrueType fonts in PDF
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["figure.dpi"] = 300
plt.rcParams["axes.linewidth"] = 0.8


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

def parse_bool(x):
    """Robustly parse boolean-like values from CSV."""
    return str(x).strip().lower() in {"true", "t", "1", "yes", "y"}


def classify_group(row):
    """Classify each sequence into Negative control, Nguyen cortical, or Other."""
    seq_id = str(row["id"]).lower()
    is_control = parse_bool(row.get("control", False))

    if is_control or "nguyen_neg" in seq_id:
        return "Negative control"
    if "nguyen_cortical" in seq_id:
        return "Nguyen cortical"
    return "Other"


def classify_species(seq_id):
    """Classify Nguyen cortical sequence species from the ID."""
    seq_id = str(seq_id)
    if "_mm_" in seq_id:
        return "Mouse"
    if "_hs_" in seq_id:
        return "Human"
    return "Unknown"


def benjamini_hochberg(pvalues):
    """Benjamini-Hochberg FDR correction. Returns q-values in original order."""
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)

    order = np.argsort(pvalues)
    ranked_p = pvalues[order]

    q_ranked = ranked_p * n / np.arange(1, n + 1)

    # Enforce monotonicity from largest rank to smallest rank
    q_ranked = np.minimum.accumulate(q_ranked[::-1])[::-1]
    q_ranked = np.minimum(q_ranked, 1.0)

    qvalues = np.empty(n, dtype=float)
    qvalues[order] = q_ranked

    return qvalues


def one_sided_mwu(group_values, control_values):
    """
    One-sided Mann-Whitney U test.

    Alternative hypothesis:
        group_values > control_values
    """
    x = np.asarray(group_values, dtype=float)
    y = np.asarray(control_values, dtype=float)

    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]

    mwu = stats.mannwhitneyu(x, y, alternative="greater")
    welch = stats.ttest_ind(x, y, alternative="greater", equal_var=False)

    # Effect size: Cliff's delta and P(x > y)
    gt = sum(a > b for a in x for b in y)
    lt = sum(a < b for a in x for b in y)
    n_pairs = len(x) * len(y)

    cliffs_delta = (gt - lt) / n_pairs
    p_x_gt_y = gt / n_pairs

    return {
        "n_group": len(x),
        "n_control": len(y),
        "median_group": float(np.median(x)),
        "median_control": float(np.median(y)),
        "delta_median": float(np.median(x) - np.median(y)),
        "mannwhitney_U": float(mwu.statistic),
        "mannwhitney_one_sided_p": float(mwu.pvalue),
        "welch_t": float(welch.statistic),
        "welch_one_sided_p": float(welch.pvalue),
        "cliffs_delta": float(cliffs_delta),
        "p_group_gt_control": float(p_x_gt_y),
    }


# -----------------------------------------------------------------------------
# Load and annotate data
# -----------------------------------------------------------------------------

df = pd.read_csv(input_csv)

# Usually the sequence ID is in an unnamed first column
id_col = "Unnamed: 0" if "Unnamed: 0" in df.columns else df.columns[0]
df = df.rename(columns={id_col: "id"})

df["group"] = df.apply(classify_group, axis=1)
df["species"] = df["id"].apply(classify_species)

plot_df = df[df["group"].isin(["Negative control", "Nguyen cortical"])].copy()

negative = plot_df.loc[
    plot_df["group"] == "Negative control",
    "mad.score"
].dropna()

nguyen_all = plot_df.loc[
    plot_df["group"] == "Nguyen cortical",
    "mad.score"
].dropna()

nguyen_mouse = plot_df.loc[
    (plot_df["group"] == "Nguyen cortical") &
    (plot_df["species"] == "Mouse"),
    "mad.score"
].dropna()

nguyen_human = plot_df.loc[
    (plot_df["group"] == "Nguyen cortical") &
    (plot_df["species"] == "Human"),
    "mad.score"
].dropna()


# -----------------------------------------------------------------------------
# Planned one-sided tests
# -----------------------------------------------------------------------------

test_rows = []

comparisons = {
    "All Nguyen cortical vs Negative control": nguyen_all,
    "Mouse Nguyen cortical vs Negative control": nguyen_mouse,
}

for comparison_name, values in comparisons.items():
    res = one_sided_mwu(values, negative)
    res["comparison"] = comparison_name
    test_rows.append(res)

test_summary = pd.DataFrame(test_rows)

# BH-FDR across the two planned one-sided MWU tests
test_summary["bh_fdr_q"] = benjamini_hochberg(
    test_summary["mannwhitney_one_sided_p"].values
)

# Reorder columns
test_summary = test_summary[
    [
        "comparison",
        "n_group",
        "n_control",
        "median_group",
        "median_control",
        "delta_median",
        "mannwhitney_U",
        "mannwhitney_one_sided_p",
        "bh_fdr_q",
        "welch_t",
        "welch_one_sided_p",
        "cliffs_delta",
        "p_group_gt_control",
    ]
]


# -----------------------------------------------------------------------------
# Save data tables
# -----------------------------------------------------------------------------

plot_df[["id", "group", "species", "mad.score", "control"]].to_csv(
    f"{out_prefix}_plot_data.csv",
    index=False,
)

test_summary.to_csv(
    f"{out_prefix}_test_summary.csv",
    index=False,
)


# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(5.7, 5.0), constrained_layout=True)

pos_neg = 0
pos_nguyen = 1

rng = np.random.default_rng(7)

negative_values = negative.to_numpy()
mouse_values = nguyen_mouse.to_numpy()
human_values = nguyen_human.to_numpy()
nguyen_all_values = nguyen_all.to_numpy()

# Negative control: gray
neg_jitter = rng.normal(0, 0.06, size=len(negative_values))
ax.scatter(
    np.full(len(negative_values), pos_neg) + neg_jitter,
    negative_values,
    s=42,
    alpha=0.9,
    color="0.45",
    label=f"Negative control (n={len(negative_values)})",
    zorder=3,
)

# Mouse Nguyen: orange
mouse_jitter = rng.normal(-0.06, 0.035, size=len(mouse_values))
ax.scatter(
    np.full(len(mouse_values), pos_nguyen) + mouse_jitter,
    mouse_values,
    s=48,
    alpha=0.9,
    color="tab:orange",
    label=f"Nguyen cortical mouse (n={len(mouse_values)})",
    zorder=4,
)

# Human Nguyen: blue
human_jitter = rng.normal(0.06, 0.035, size=len(human_values))
ax.scatter(
    np.full(len(human_values), pos_nguyen) + human_jitter,
    human_values,
    s=48,
    alpha=0.9,
    color="tab:blue",
    label=f"Nguyen cortical human (n={len(human_values)})",
    zorder=4,
)

# Boxplots: control vs all Nguyen cortical
ax.boxplot(
    [negative_values, nguyen_all_values],
    positions=[pos_neg, pos_nguyen],
    widths=0.5,
    patch_artist=False,
    showfliers=False,
    medianprops={"linewidth": 1.4},
    boxprops={"linewidth": 1.0},
    whiskerprops={"linewidth": 1.0},
    capprops={"linewidth": 1.0},
)

ax.axhline(0, linestyle="--", linewidth=0.8)

# Annotation
all_row = test_summary[
    test_summary["comparison"] == "All Nguyen cortical vs Negative control"
].iloc[0]

mouse_row = test_summary[
    test_summary["comparison"] == "Mouse Nguyen cortical vs Negative control"
].iloc[0]

y_all = plot_df["mad.score"].to_numpy()
y_min = float(np.nanmin(y_all))
y_max = float(np.nanmax(y_all))
y_range = y_max - y_min if y_max > y_min else 1.0

y1 = y_max + y_range * 0.12
h = y_range * 0.04

ax.plot(
    [pos_neg, pos_neg, pos_nguyen, pos_nguyen],
    [y1, y1 + h, y1 + h, y1],
    linewidth=1.0,
)

ax.text(
    0.5,
    y1 + h + y_range * 0.02,
    f"All Nguyen vs control\nMWU p={all_row['mannwhitney_one_sided_p']:.3g}, FDR q={all_row['bh_fdr_q']:.3g}",
    ha="center",
    va="bottom",
    fontsize=8.5,
)

ax.text(
    0.98,
    0.04,
    f"Mouse vs control:\n"
    f"MWU p={mouse_row['mannwhitney_one_sided_p']:.3g}, FDR q={mouse_row['bh_fdr_q']:.3g}\n"
    f"median {mouse_row['median_group']:.3f} vs {mouse_row['median_control']:.3f}\n"
    f"n={int(mouse_row['n_group'])} vs {int(mouse_row['n_control'])}",
    transform=ax.transAxes,
    ha="right",
    va="bottom",
    fontsize=8.5,
)

ax.set_xticks([pos_neg, pos_nguyen])
ax.set_xticklabels(["Negative control", "Nguyen cortical"])
ax.set_ylabel("MAD score")
ax.set_title("Pseudobarcode R1R2", fontsize=11)
ax.legend(frameon=False, fontsize=8, loc="upper left")

ax.set_ylim(y_min - y_range * 0.08, y1 + h + y_range * 0.18)

fig.savefig(f"{out_prefix}.png", bbox_inches="tight", dpi=300)
fig.savefig(f"{out_prefix}_editable.svg", bbox_inches="tight")
fig.savefig(f"{out_prefix}_editable.pdf", bbox_inches="tight")

plt.close(fig)


# -----------------------------------------------------------------------------
# Print summary
# -----------------------------------------------------------------------------

print(test_summary.to_string(index=False))
