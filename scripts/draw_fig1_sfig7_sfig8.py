# -*- coding: utf-8 -*-
"""Draw Fig.1, SFig.7, and SFig.8 for the PC-CFRL TKDE package.

The figures follow the existing manuscript figure contract:
- white background, Times New Roman, thin boxed axes/borders;
- one fixed colorblind-safe palette across main and supplementary figures;
- legends placed below the plotting area;
- PDF plus 600 dpi PNG exports, with source CSV files saved next to figures.
"""

from __future__ import annotations

from pathlib import Path
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd


PALETTE = {
    "pc_cfrl": "#0072B2",
    "ecfp": "#4D4D4D",
    "extratrees": "#009E73",
    "scalar_only": "#E69F00",
    "pair_only": "#56B4E9",
    "negative": "#D55E00",
    "neutral": "#B3B3B3",
    "light_neutral": "#E8E8E8",
    "border": "#1A1A1A",
    "grid": "#EAEAEA",
    "missing": "#F5F5F5",
    "white": "#FFFFFF",
}

CATEGORY_COLORS = {
    "similarity_counts": PALETTE["pair_only"],
    "descriptor_absdiff": PALETTE["pc_cfrl"],
    "descriptor_sum": PALETTE["scalar_only"],
    "diff_bits": PALETTE["extratrees"],
}

METHOD_COLORS = {
    "ecfp_absdiff_xgb": PALETTE["ecfp"],
    "rich_diff_scalars": PALETTE["pc_cfrl"],
}

SPLIT_LABELS = {
    "target_family": "target-cluster",
    "family_scaffold_purged": "target-cluster+scaffold",
    "bindingdb_to_chembl_pair_scaffold_purged": "BindingDB->ChEMBL\npair+scaffold",
    "chembl_to_bindingdb_pair_scaffold_purged": "ChEMBL->BindingDB\npair+scaffold",
    "chembl_future_after_2005_scaffold_source_purged": "future>2005\nscaffold+source",
    "chembl_future_after_2010_scaffold_source_purged": "future>2010\nscaffold+source",
    "chembl_future_after_2015_scaffold_source_purged": "future>2015\nscaffold+source",
}


def out_dir() -> Path:
    return Path(__file__).resolve().parent


def repo_root() -> Path:
    return out_dir().parent


def table_dir() -> Path:
    return repo_root() / "PC-CFRL_TKDE_trans_framework_20260607" / "data_remote" / "tkde_experiment_tables_2026_05_20"


def remote_extra_dir() -> Path:
    return out_dir() / "remote_extra"


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
            "lines.linewidth": 0.85,
            "xtick.major.width": 0.45,
            "ytick.major.width": 0.45,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "axes.grid": False,
        }
    )


def save_figure(fig: mpl.figure.Figure, stem: str) -> None:
    target = out_dir() / stem
    pad = 25 / 72
    fig.savefig(target.with_suffix(".pdf"), bbox_inches="tight", pad_inches=pad)
    fig.savefig(target.with_suffix(".png"), dpi=600, bbox_inches="tight", pad_inches=pad)
    plt.close(fig)


def set_box(ax: mpl.axes.Axes) -> None:
    for side in ("left", "right", "top", "bottom"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(0.45)
        ax.spines[side].set_color(PALETTE["border"])
    ax.tick_params(width=0.45, colors=PALETTE["border"])


def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        fontweight="bold",
        color=PALETTE["border"],
    )


def mean_std(value) -> tuple[float, float]:
    if pd.isna(value):
        return np.nan, np.nan
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value), np.nan
    match = re.match(r"\s*([+-]?\d+(?:\.\d+)?)\s*\+/-\s*([+-]?\d+(?:\.\d+)?)", str(value))
    if match:
        return float(match.group(1)), float(match.group(2))
    try:
        return float(value), np.nan
    except ValueError:
        return np.nan, np.nan


def draw_box(ax, xy, wh, text, face, edge=None, fontsize=9.0, weight="normal", radius=0.025):
    edge = edge or PALETTE["border"]
    box = FancyBboxPatch(
        xy,
        wh[0],
        wh[1],
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        linewidth=0.45,
        edgecolor=edge,
        facecolor=face,
        mutation_aspect=1,
    )
    ax.add_patch(box)
    ax.text(
        xy[0] + wh[0] / 2,
        xy[1] + wh[1] / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight=weight,
        color=PALETTE["border"],
        linespacing=1.15,
    )
    return box


def draw_arrow(ax, start, end, color=PALETTE["neutral"], style="solid", rad=0.0):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=8,
        lw=0.55,
        color=color,
        linestyle=style,
        connectionstyle=f"arc3,rad={rad}",
        shrinkA=2,
        shrinkB=2,
    )
    ax.add_patch(arrow)


def draw_fig1() -> None:
    fig, ax = plt.subplots(figsize=(7.15, 3.15))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    stages = [
        {
            "x": 0.035,
            "w": 0.205,
            "accent": PALETTE["ecfp"],
            "label": "Raw sources",
            "items": ["ACNet", "MoleculeACE", "ChEMBL", "BindingDB"],
        },
        {
            "x": 0.285,
            "w": 0.205,
            "accent": PALETTE["negative"],
            "label": "Leakage contract",
            "items": ["target cluster", "scaffold / ligand", "source / temporal", "overlap audit"],
        },
        {
            "x": 0.535,
            "w": 0.205,
            "accent": PALETTE["pc_cfrl"],
            "label": "Pair features",
            "items": ["ECFP absdiff", "similarity / counts", "descriptor Δ and Σ", "difference bits"],
        },
        {
            "x": 0.785,
            "w": 0.18,
            "accent": PALETTE["extratrees"],
            "label": "Evidence",
            "items": ["hard-OOD", "cross-db\n/ temporal", "reliability", "claim boundary"],
        },
    ]

    nodes = []
    source_rows = []
    for stage in stages:
        x, w = stage["x"], stage["w"]
        ax.add_patch(
            Rectangle(
                (x, 0.14),
                w,
                0.72,
                facecolor=PALETTE["white"],
                edgecolor=PALETTE["light_neutral"],
                linewidth=0.45,
            )
        )
        ax.add_patch(Rectangle((x, 0.14), 0.008, 0.72, facecolor=stage["accent"], edgecolor="none"))
        ax.text(
            x + 0.018,
            0.80,
            stage["label"],
            ha="left",
            va="center",
            fontsize=8.2,
            fontweight="bold",
            color=PALETTE["border"],
        )
        for idx, item in enumerate(stage["items"]):
            y = 0.665 - idx * 0.125
            face = "#F7FBFD" if stage["accent"] == PALETTE["pc_cfrl"] else PALETTE["missing"]
            node = draw_box(ax, (x + 0.026, y - 0.034), (w - 0.052, 0.068), item, face=face, fontsize=7.65)
            nodes.append(node)
            source_rows.append({"stage": stage["label"], "item": item})

    draw_arrow(ax, (0.242, 0.50), (0.282, 0.50), color=PALETTE["ecfp"])
    draw_arrow(ax, (0.492, 0.50), (0.532, 0.50), color=PALETTE["negative"])
    draw_arrow(ax, (0.742, 0.50), (0.782, 0.50), color=PALETTE["pc_cfrl"])

    draw_box(
        ax,
        (0.555, 0.168),
        (0.165, 0.052),
        "PC-CFRL + controls",
        face="#EAF4FB",
        edge=PALETTE["pc_cfrl"],
        fontsize=7.5,
        weight="bold",
    )

    ax.text(
        0.50,
        0.075,
        "All main claims are read from frozen manifests and paired seed-fold evaluations; negative controls are retained to bound the scope.",
        ha="center",
        va="center",
        fontsize=7.4,
        color=PALETTE["border"],
    )

    pd.DataFrame(source_rows).to_csv(out_dir() / "Fig1_source_schema.csv", index=False)
    save_figure(fig, "Fig1_overview_workflow")


def draw_sfig7() -> None:
    gain = pd.read_csv(remote_extra_dir() / "acnet_explanation_stability_seed0_4_n10_category_summary.csv")
    shap = pd.read_csv(remote_extra_dir() / "acnet_treeshap_stability_seed0_4_n10_category_summary.csv")

    gain = gain.assign(source="total_gain", share=gain["gain_share_mean"], std=gain["gain_share_std"])
    shap = shap.assign(source="treeshap", share=shap["contribution_share_mean"], std=shap["contribution_share_std"])
    cols = ["source", "split_mode", "feature_category", "category_label", "share", "std", "n_folds"]
    src = pd.concat([gain[cols], shap[cols]], ignore_index=True)
    src.to_csv(out_dir() / "SFig7_source_data.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(7.15, 2.65), sharey=True)
    panels = [("total_gain", axes[0], "Total-gain share (%)"), ("treeshap", axes[1], "TreeSHAP contribution share (%)")]
    split_order = ["target_family", "family_scaffold_purged"]
    cat_order = ["similarity_counts", "descriptor_absdiff", "descriptor_sum", "diff_bits"]

    for idx, (source, ax, xlabel) in enumerate(panels):
        panel = src[src["source"] == source]
        left = np.zeros(len(split_order))
        y = np.arange(len(split_order))
        for cat in cat_order:
            vals = []
            for split in split_order:
                row = panel[(panel["split_mode"] == split) & (panel["feature_category"] == cat)].iloc[0]
                vals.append(row["share"] * 100)
            ax.barh(
                y,
                vals,
                left=left,
                color=CATEGORY_COLORS[cat],
                edgecolor=PALETTE["border"],
                linewidth=0.28,
                height=0.46,
            )
            for yi, v, lft in zip(y, vals, left):
                if v >= 8:
                    ax.text(lft + v / 2, yi, f"{v:.1f}", ha="center", va="center", fontsize=7.1, color=PALETTE["border"])
            left += np.array(vals)
        ax.set_xlim(0, 100)
        ax.set_yticks(y)
        ax.set_yticklabels([SPLIT_LABELS[s] for s in split_order])
        ax.invert_yaxis()
        ax.set_xlabel(xlabel)
        ax.xaxis.grid(True, color=PALETTE["grid"], linewidth=0.35)
        set_box(ax)
        panel_label(ax, f"({chr(ord('a') + idx)})")

    handles = [
        Rectangle((0, 0), 1, 1, facecolor=CATEGORY_COLORS[c], edgecolor=PALETTE["border"], linewidth=0.28)
        for c in cat_order
    ]
    labels = ["similarity/counts", "descriptor absdiff", "descriptor sums", "difference bits"]
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.005), ncol=4, frameon=False)
    fig.subplots_adjust(left=0.16, right=0.985, top=0.91, bottom=0.26, wspace=0.12)
    save_figure(fig, "SFig7_explanation_stability")


def draw_sfig8() -> None:
    df = pd.read_csv(table_dir() / "table_05_crossdb_prospective_topk.csv")
    split_order = [
        "bindingdb_to_chembl_pair_scaffold_purged",
        "chembl_to_bindingdb_pair_scaffold_purged",
        "chembl_future_after_2005_scaffold_source_purged",
        "chembl_future_after_2010_scaffold_source_purged",
        "chembl_future_after_2015_scaffold_source_purged",
    ]
    variants = ["ecfp_absdiff_xgb", "rich_diff_scalars"]
    metrics = [
        ("bedroc20", "BEDROC20"),
        ("ef_5pct", "EF@5%"),
        ("precision_at_50", "P@50"),
    ]

    records = []
    for split in split_order:
        for variant in variants:
            row = df[(df["split_mode"] == split) & (df["variant"] == variant)].iloc[0]
            for metric, label in metrics:
                mean, std = mean_std(row[metric])
                records.append(
                    {
                        "split_mode": split,
                        "split_label": SPLIT_LABELS[split],
                        "variant": variant,
                        "metric": metric,
                        "metric_label": label,
                        "mean": mean,
                        "std": std,
                    }
                )
    src = pd.DataFrame(records)
    src.to_csv(out_dir() / "SFig8_source_data.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(7.15, 3.05), sharey=True)
    y = np.arange(len(split_order))
    for idx, ((metric, xlabel), ax) in enumerate(zip(metrics, axes)):
        panel = src[src["metric"] == metric]
        ecfp = panel[panel["variant"] == "ecfp_absdiff_xgb"].set_index("split_mode").loc[split_order]
        pc = panel[panel["variant"] == "rich_diff_scalars"].set_index("split_mode").loc[split_order]
        for yi, split in enumerate(split_order):
            x0 = ecfp.loc[split, "mean"]
            x1 = pc.loc[split, "mean"]
            ax.plot([x0, x1], [yi, yi], color=PALETTE["light_neutral"], linewidth=1.15, zorder=1)
        ax.errorbar(
            ecfp["mean"],
            y,
            xerr=ecfp["std"],
            fmt="o",
            ms=3.2,
            color=PALETTE["ecfp"],
            ecolor=PALETTE["ecfp"],
            elinewidth=0.55,
            capsize=1.8,
            markerfacecolor=PALETTE["white"],
            markeredgewidth=0.55,
            zorder=3,
        )
        ax.errorbar(
            pc["mean"],
            y,
            xerr=pc["std"],
            fmt="o",
            ms=3.2,
            color=PALETTE["pc_cfrl"],
            ecolor=PALETTE["pc_cfrl"],
            elinewidth=0.55,
            capsize=1.8,
            markerfacecolor=PALETTE["pc_cfrl"],
            markeredgewidth=0.55,
            zorder=4,
        )
        ax.set_xlabel(xlabel)
        ax.set_yticks(y)
        ax.set_yticklabels([SPLIT_LABELS[s] for s in split_order])
        ax.invert_yaxis()
        ax.xaxis.grid(True, color=PALETTE["grid"], linewidth=0.35)
        ax.margins(x=0.12)
        set_box(ax)
        panel_label(ax, f"({chr(ord('a') + idx)})")

    handles = [
        Line2D([0], [0], color=PALETTE["ecfp"], marker="o", markerfacecolor=PALETTE["white"], lw=0, label="ECFP-XGB"),
        Line2D([0], [0], color=PALETTE["pc_cfrl"], marker="o", markerfacecolor=PALETTE["pc_cfrl"], lw=0, label="PC-CFRL"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.005), ncol=2, frameon=False)
    fig.subplots_adjust(left=0.215, right=0.985, top=0.90, bottom=0.255, wspace=0.13)
    save_figure(fig, "SFig8_topk_ranking_utility")


def main() -> None:
    apply_style()
    draw_fig1()
    draw_sfig7()
    draw_sfig8()


if __name__ == "__main__":
    main()
