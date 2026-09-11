# -*- coding: utf-8 -*-
"""Draw SFig.9 for the PC-CFRL-Bench supplement.

The figure summarizes query-disjoint few-shot target-adaptation protocols from
the 2026-06-15 fairness/active-learning runs. It keeps the existing manuscript
figure contract: white background, Times New Roman, thin boxed axes, fixed
palette, and PDF plus 600 dpi PNG export.
"""

from __future__ import annotations

from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "remote_extra" / "fewshot_20260615"


def load_palette() -> dict[str, str]:
    with open(HERE / "figure_palette_fixed.json", "r", encoding="utf-8") as fh:
        return json.load(fh)


PALETTE = load_palette()
PROTOCOL_STYLE = {
    "natural_frac": ("natural random", PALETTE["pc_cfrl"], "o"),
    "stratified_frac_previous": ("stratified previous", PALETTE["extratrees"], "s"),
    "best_active": ("best active", PALETTE["scalar_only"], "D"),
}
BENCH_LABEL = {"bindingdb_raw": "BindingDB", "chembl_raw": "ChEMBL"}


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Times New Roman",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.edgecolor": PALETTE["border"],
            "axes.linewidth": 0.45,
            "axes.labelsize": 10.5,
            "xtick.labelsize": 9.2,
            "ytick.labelsize": 9.2,
            "legend.fontsize": 8.0,
            "lines.linewidth": 0.95,
            "xtick.major.width": 0.45,
            "ytick.major.width": 0.45,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "axes.grid": False,
        }
    )


def boxed_axis(ax: mpl.axes.Axes) -> None:
    for side in ("left", "right", "top", "bottom"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(0.45)
        ax.spines[side].set_color(PALETTE["border"])
    ax.tick_params(width=0.45, colors=PALETTE["border"])
    ax.grid(axis="y", color=PALETTE["grid"], linewidth=0.35)


def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(
        -0.13,
        1.045,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        fontweight="bold",
        color=PALETTE["border"],
    )


def plot_metric(ax: mpl.axes.Axes, df: pd.DataFrame, benchmark: str, metric: str) -> None:
    sub = df[df["benchmark"].eq(benchmark)].copy()
    for protocol, (name, color, marker) in PROTOCOL_STYLE.items():
        cur = sub[sub["protocol"].eq(protocol)].sort_values("support_percent")
        ax.plot(
            cur["support_percent"],
            cur[metric],
            color=color,
            marker=marker,
            markersize=3.8,
            markeredgewidth=0.45,
            markeredgecolor=PALETTE["border"],
            label=name,
        )
    ax.set_xlim(3, 52)
    ax.set_xticks([5, 10, 20, 30, 50])
    ax.set_xlabel("Support labels (%)")
    ax.set_ylabel("ROC-AUC" if metric == "roc_auc_mean" else "MCC")
    boxed_axis(ax)


def main() -> None:
    apply_style()
    path = DATA_DIR / "tkde_fewshot_protocol_fairness_active_final_20260615_stratified_natural_active_compare.csv"
    df = pd.read_csv(path)
    df.to_csv(HERE / "SFig9_source_data.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.4), constrained_layout=False)
    panels = [
        ("bindingdb_raw", "roc_auc_mean"),
        ("chembl_raw", "roc_auc_mean"),
        ("bindingdb_raw", "mcc_mean"),
        ("chembl_raw", "mcc_mean"),
    ]
    labels = ["(a)", "(b)", "(c)", "(d)"]
    for ax, (benchmark, metric), label in zip(axes.flat, panels, labels):
        plot_metric(ax, df, benchmark, metric)
        panel_label(ax, label)
        if metric == "roc_auc_mean":
            ax.set_ylim(0.55, 0.83)
        else:
            ax.set_ylim(0.10, 0.46)
        ax.text(
            0.96,
            0.08,
            BENCH_LABEL[benchmark],
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8.5,
            color=PALETTE["border"],
        )

    handles, labels_ = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels_,
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.005),
        handlelength=1.7,
        columnspacing=1.8,
    )
    fig.subplots_adjust(left=0.085, right=0.985, top=0.965, bottom=0.15, wspace=0.23, hspace=0.28)

    pad = 25 / 72
    fig.savefig(HERE / "SFig9_fewshot_label_efficiency.pdf", bbox_inches="tight", pad_inches=pad)
    fig.savefig(HERE / "SFig9_fewshot_label_efficiency.png", dpi=600, bbox_inches="tight", pad_inches=pad)
    plt.close(fig)


if __name__ == "__main__":
    main()
