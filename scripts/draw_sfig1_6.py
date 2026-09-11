# -*- coding: utf-8 -*-
"""Draw all supplementary figures for the PC-CFRL TKDE manuscript.

Figure contracts:
- SFig.1: split construction and leakage predicates.
- SFig.2: raw-data quality and fold-level curation distributions.
- SFig.3: full method-by-split performance heatmaps.
- SFig.4: calibration and risk-coverage diagnostics.
- SFig.5: threshold/noise/censoring sensitivity.
- SFig.6: RDKit molecular-pair case gallery.

Backend: Python/matplotlib. RDKit is used only for 2D molecule cards.
Export: white background, Times New Roman, PDF + 600 dpi PNG, fixed palette.
"""

from __future__ import annotations

from pathlib import Path
import json
import math
import re
from textwrap import shorten

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch
import numpy as np
import pandas as pd
from PIL import Image

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdFMCS


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

METHOD_STYLE = {
    "ecfp_absdiff_xgb": ("ECFP-XGB", PALETTE["ecfp"], "o"),
    "rich_diff_scalars": ("PC-CFRL", PALETTE["pc_cfrl"], "o"),
    "rich_pair_only": ("pair-only", PALETTE["pair_only"], "s"),
    "rich_scalars_only": ("scalar-only", PALETTE["scalar_only"], "D"),
    "extratrees_scalars": ("ExtraTrees", PALETTE["extratrees"], "^"),
}

SPLIT_NAME = {
    "target_family": "target-cluster",
    "family_scaffold_purged": "target-cluster+scaffold",
    "family_scaffold_source_purged": "target-cluster+scaffold+source",
    "temporal_scaffold_source_purged": "temporal+source",
    "family_exact_ligand_purged": "target-cluster+ligand",
    "lohi_ligand_scaffold_source_purged": "LoHi ligand",
    "lohi_pair_scaffold_source_purged": "LoHi pair",
    "bindingdb_to_chembl_scaffold_purged": "BDB->ChEMBL scaf.",
    "chembl_to_bindingdb_scaffold_purged": "ChEMBL->BDB scaf.",
    "bindingdb_to_chembl_exact_ligand_purged": "BDB->ChEMBL lig.",
    "chembl_to_bindingdb_exact_ligand_purged": "ChEMBL->BDB lig.",
    "bindingdb_to_chembl_pair_scaffold_purged": "BDB->ChEMBL pair",
    "chembl_to_bindingdb_pair_scaffold_purged": "ChEMBL->BDB pair",
    "bindingdb_to_chembl_matched_targets": "BDB->ChEMBL target",
    "chembl_to_bindingdb_matched_targets": "ChEMBL->BDB target",
    "chembl_future_after_2005_scaffold_source_purged": "future>2005 scaf.",
    "chembl_future_after_2010_scaffold_source_purged": "future>2010 scaf.",
    "chembl_future_after_2015_scaffold_source_purged": "future>2015 scaf.",
    "chembl_future_after_2005_future_only": "future>2005",
    "chembl_future_after_2010_future_only": "future>2010",
    "chembl_future_after_2015_future_only": "future>2015",
    "temporal_forward": "temporal",
}

BENCH_NAME = {
    "acnet": "ACNet",
    "chembl_raw": "ChEMBL",
    "bindingdb_raw": "BindingDB",
    "moleculeace": "MoleculeACE",
    "ACNet": "ACNet",
    "ChEMBL raw": "ChEMBL",
    "BindingDB raw": "BindingDB",
    "MoleculeACE": "MoleculeACE",
    "LoHi raw": "LoHi",
    "Noise raw": "Noise",
    "Remaining P0-P1": "P0-P1",
}

SHORT_BENCH_NAME = {
    "acnet": "ACN",
    "chembl_raw": "ChE",
    "bindingdb_raw": "BDB",
    "moleculeace": "ACE",
}

SHORT_SPLIT_NAME = {
    "target_family": "t-clust.",
    "family_scaffold_purged": "t-clust.+scaf.",
    "family_scaffold_source_purged": "t-clust.+scaf.+src",
    "temporal_scaffold_source_purged": "tmp",
}


def data_root() -> Path:
    here = Path(__file__).resolve()
    candidates = [
        here.parent.parent / "PC-CFRL_TKDE_trans_framework_20260607" / "data_remote",
        Path(r"C:\codex_tmp\paper15_tkde_remote"),
    ]
    for cand in candidates:
        if (cand / "tkde_experiment_tables_2026_05_20").exists():
            return cand
    raise FileNotFoundError("Cannot locate data_remote/tkde_experiment_tables_2026_05_20")


def table_dir() -> Path:
    return data_root() / "tkde_experiment_tables_2026_05_20"


def extra_dir() -> Path:
    root = data_root()
    if (root / "extra_csv").exists():
        return root / "extra_csv"
    return Path(r"C:\codex_tmp\paper15_tkde_remote\extra_csv")


def out_dir() -> Path:
    out = Path(__file__).resolve().parent
    out.mkdir(parents=True, exist_ok=True)
    return out


def read_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(table_dir() / name)


def read_extra(name: str) -> pd.DataFrame:
    return pd.read_csv(extra_dir() / name)


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
            "xtick.labelsize": 9.0,
            "ytick.labelsize": 9.0,
            "legend.fontsize": 8.0,
            "lines.linewidth": 0.85,
            "xtick.major.width": 0.45,
            "ytick.major.width": 0.45,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "axes.grid": False,
        }
    )


def save_figure(fig: mpl.figure.Figure, stem: str, pdf_from_png: bool = False) -> None:
    target = out_dir() / stem
    pad = 25 / 72
    pdf_path = target.with_suffix(".pdf")
    png_path = target.with_suffix(".png")
    if pdf_from_png:
        fig.savefig(png_path, dpi=600, bbox_inches="tight", pad_inches=pad)
        with Image.open(png_path) as img:
            if img.mode == "RGBA":
                white = Image.new("RGB", img.size, "white")
                white.paste(img, mask=img.getchannel("A"))
                img = white
            else:
                img = img.convert("RGB")
            img.save(pdf_path, "PDF", resolution=600.0)
    else:
        fig.savefig(pdf_path, bbox_inches="tight", pad_inches=pad)
        fig.savefig(png_path, dpi=600, bbox_inches="tight", pad_inches=pad)
    plt.close(fig)


def set_box(ax: mpl.axes.Axes) -> None:
    for side in ("left", "right", "top", "bottom"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(0.45)
        ax.spines[side].set_color(PALETTE["border"])
    ax.tick_params(width=0.45, colors=PALETTE["border"])


def remove_box(ax: mpl.axes.Axes) -> None:
    for side in ("left", "right", "top", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)


def panel_label(ax: mpl.axes.Axes, label: str, x: float = -0.08, y: float = 1.035) -> None:
    ax.text(
        x,
        y,
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
    match = re.match(r"\s*([+-]?\d+(?:\.\d+)?)", str(value))
    if match:
        return float(match.group(1)), np.nan
    return np.nan, np.nan


def parse_series(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        if col in out.columns:
            vals = out[col].apply(lambda x: mean_std(x)[0])
            out[col] = pd.to_numeric(vals, errors="coerce")
    return out


def task_label(row: pd.Series) -> str:
    bench = BENCH_NAME.get(str(row.get("benchmark", "")), str(row.get("benchmark", "")))
    split = SPLIT_NAME.get(str(row.get("split_mode", "")), str(row.get("split_mode", "")))
    condition = str(row.get("condition", ""))
    if condition and condition not in {"default", "nan", "None"}:
        condition = {
            "equal_relation_rebuilt": "exact relations only",
        }.get(condition, condition.replace("_", " "))
        return f"{bench}\n{condition}"
    return f"{bench}\n{split}"


def draw_arrow(ax: mpl.axes.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=8,
        lw=0.65,
        color=PALETTE["border"],
        shrinkA=3,
        shrinkB=3,
    )
    ax.add_patch(arrow)


def flow_box(
    ax: mpl.axes.Axes,
    xy: tuple[float, float],
    wh: tuple[float, float],
    text: str,
    face: str,
    edge: str,
) -> None:
    box = FancyBboxPatch(
        xy,
        wh[0],
        wh[1],
        boxstyle="round,pad=0.012,rounding_size=0.018",
        facecolor=face,
        edgecolor=edge,
        linewidth=0.6,
    )
    ax.add_patch(box)
    ax.text(
        xy[0] + wh[0] / 2,
        xy[1] + wh[1] / 2,
        text,
        ha="center",
        va="center",
        fontsize=8.5,
        color=PALETTE["border"],
    )


def draw_sfig1() -> None:
    audit = read_csv("table_01_split_leakage_audit.csv")
    audit = parse_series(
        audit,
        [
            "train_retention_after_purge",
            "target_overlap_rate",
            "exact_ligand_overlap_rate",
            "scaffold_overlap_rate",
            "exact_pair_overlap_rate",
            "pair_scaffold_overlap_rate",
            "document_source_overlap_rate",
            "assay_source_overlap_rate",
        ],
    )
    fig, ax_mat = plt.subplots(figsize=(6.7, 4.65))
    plt.subplots_adjust(left=0.31, right=0.88, top=0.96, bottom=0.18)
    cols = [
        ("target_overlap_rate", "target"),
        ("exact_ligand_overlap_rate", "ligand"),
        ("scaffold_overlap_rate", "scaffold"),
        ("exact_pair_overlap_rate", "pair"),
        ("document_source_overlap_rate", "doc."),
        ("assay_source_overlap_rate", "assay"),
    ]
    ordered = audit.copy()
    ordered["row_label"] = ordered.apply(
        lambda r: f"{BENCH_NAME.get(r['benchmark'], r['benchmark'])} | {SPLIT_NAME.get(r['split_mode'], r['split_mode'])}",
        axis=1,
    )
    mat = ordered[[c[0] for c in cols]].fillna(0.0).to_numpy(dtype=float)
    cmap = LinearSegmentedColormap.from_list("leakage", ["#FFFFFF", PALETTE["negative"]])
    im = ax_mat.imshow(mat, aspect="auto", cmap=cmap, norm=Normalize(vmin=0, vmax=max(1.0, np.nanmax(mat))))
    ax_mat.set_xticks(np.arange(len(cols)))
    ax_mat.set_xticklabels([c[1] for c in cols], rotation=35, ha="right")
    ax_mat.set_yticks(np.arange(len(ordered)))
    ax_mat.set_yticklabels(ordered["row_label"].tolist(), fontsize=7.9)
    ax_mat.set_xticks(np.arange(-0.5, len(cols), 1), minor=True)
    ax_mat.set_yticks(np.arange(-0.5, len(ordered), 1), minor=True)
    ax_mat.grid(which="minor", color=PALETTE["grid"], linewidth=0.45)
    ax_mat.tick_params(which="minor", bottom=False, left=False)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            value = mat[i, j]
            label = "0" if abs(value) < 0.0005 else f"{value:.2f}"
            ax_mat.text(j, i, label, ha="center", va="center", fontsize=6.8, color=PALETTE["border"])
    set_box(ax_mat)
    cbar = fig.colorbar(im, ax=ax_mat, fraction=0.045, pad=0.018)
    cbar.set_label("overlap rate", fontsize=8.5)
    cbar.ax.tick_params(labelsize=7.6, width=0.4, length=2)
    ordered.to_csv(out_dir() / "SFig1_source_data.csv", index=False)
    save_figure(fig, "SFig1_split_construction_leakage")


def grouped_barh(ax: mpl.axes.Axes, df: pd.DataFrame, categories: list[str], xscale: str = "linear") -> None:
    y = np.arange(len(categories))
    offsets = [-0.17, 0.17]
    log_floor = 1e2
    if xscale == "log":
        ax.set_xscale("log")
        ax.set_xlim(log_floor, 1e6)
    for offset, (_, row), color in zip(offsets, df.iterrows(), [PALETTE["pc_cfrl"], PALETTE["extratrees"]]):
        values: list[float] = []
        for cat in categories:
            val = float(row.get(cat, np.nan))
            values.append(val if val > 0 else np.nan)
        values_array = np.asarray(values, dtype=float)
        if xscale == "log":
            valid = np.isfinite(values_array) & (values_array >= log_floor)
            widths = np.where(valid, values_array - log_floor, np.nan)
            ax.barh(
                y + offset,
                widths,
                left=log_floor,
                height=0.28,
                color=color,
                alpha=0.82,
                edgecolor=PALETTE["border"],
                linewidth=0.25,
                zorder=3,
            )
            for yi, is_valid in zip(y + offset, valid):
                if not is_valid:
                    ax.text(log_floor * 1.12, yi, "NA", color=color, fontsize=7.2, ha="left", va="center", zorder=4)
        else:
            ax.barh(
                y + offset,
                values_array,
                height=0.28,
                color=color,
                alpha=0.82,
                edgecolor=PALETTE["border"],
                linewidth=0.25,
                zorder=3,
            )
    ax.set_yticks(y)
    ax.set_yticklabels([c.replace("_", " ") for c in categories])
    ax.grid(axis="x", color=PALETTE["grid"], linewidth=0.35, zorder=0)
    set_box(ax)


def dot_rate_panel(ax: mpl.axes.Axes, df: pd.DataFrame, columns: list[tuple[str, str]], xlabel: str) -> None:
    y = np.arange(len(columns))
    for _, row in df.iterrows():
        dataset = str(row["benchmark"])
        color = PALETTE["pc_cfrl"] if dataset == "chembl_raw" else PALETTE["extratrees"]
        marker = "o" if dataset == "chembl_raw" else "s"
        values = [float(row[col]) if pd.notna(row[col]) else np.nan for col, _ in columns]
        ax.plot(values, y, color=color, marker=marker, ms=4.1, lw=0.85, label=BENCH_NAME.get(dataset, dataset))
    ax.set_yticks(y)
    ax.set_yticklabels([label for _, label in columns])
    ax.set_xlabel(xlabel)
    ax.set_xlim(left=0)
    ax.grid(axis="x", color=PALETTE["grid"], linewidth=0.35)
    set_box(ax)


def draw_sfig2() -> None:
    quality = read_extra("tkde_raw_data_quality_summary.csv")
    manifest = read_extra("tkde_benchmark_split_manifest.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.45, 5.95))
    axes = axes.ravel()
    plt.subplots_adjust(left=0.14, right=0.985, top=0.955, bottom=0.145, wspace=0.43, hspace=0.47)

    scale_cols = ["raw_rows", "pair_rows", "target_ligand_groups", "n_source_ids", "n_assay_ids"]
    grouped_barh(axes[0], quality, scale_cols, xscale="log")
    axes[0].set_xlabel("count (log scale)")
    panel_label(axes[0], "(a)")

    rate_cols = [
        ("multi_measurement_group_rate", "multi-measure"),
        ("high_conflict_std_gt_0_5_rate", "std>0.5"),
        ("high_conflict_range_gt_1_rate", "range>1"),
        ("censored_relation_rate", "censored"),
        ("pair_near_cliff_boundary_0_8_1_2_rate", "near boundary"),
        ("pair_positive_rate", "positive pair"),
    ]
    dot_rate_panel(axes[1], quality, rate_cols, "rate")
    panel_label(axes[1], "(b)")

    var_cols = [
        ("pactivity_std_mean", "std mean"),
        ("pactivity_std_p90", "std p90"),
        ("pactivity_range_p90", "range p90"),
    ]
    dot_rate_panel(axes[2], quality, var_cols, "pActivity units")
    panel_label(axes[2], "(c)")

    m = manifest.copy()
    m["split_group"] = np.where(m["split_mode"].str.contains("purged"), "purged", "target-cluster")
    m.loc[m["split_mode"].str.contains("temporal"), "split_group"] = "temporal"
    groups = ["target-cluster", "purged", "temporal"]
    data = [m.loc[m["split_group"] == g, "train_retention_after_purge"].dropna().to_numpy(dtype=float) for g in groups]
    box = axes[3].boxplot(
        data,
        tick_labels=groups,
        patch_artist=True,
        showfliers=False,
        widths=0.52,
        boxprops={"edgecolor": PALETTE["border"], "linewidth": 0.55},
        medianprops={"color": PALETTE["pc_cfrl"], "linewidth": 0.9},
        whiskerprops={"color": PALETTE["border"], "linewidth": 0.55},
        capprops={"color": PALETTE["border"], "linewidth": 0.55},
    )
    for patch, color in zip(box["boxes"], [PALETTE["light_neutral"], PALETTE["pair_only"], PALETTE["scalar_only"]]):
        patch.set_facecolor(color)
        patch.set_alpha(0.45)
    rng = np.random.default_rng(2026)
    for idx, vals in enumerate(data, start=1):
        axes[3].scatter(
            idx + rng.uniform(-0.09, 0.09, size=len(vals)),
            vals,
            s=11,
            color=PALETTE["ecfp"],
            alpha=0.42,
            linewidth=0,
        )
    axes[3].set_ylabel("train retention")
    axes[3].set_ylim(0, 1.04)
    axes[3].grid(axis="y", color=PALETTE["grid"], linewidth=0.35)
    set_box(axes[3])
    panel_label(axes[3], "(d)")

    handles = [
        Line2D([0], [0], color=PALETTE["pc_cfrl"], marker="o", lw=0.85, label="ChEMBL"),
        Line2D([0], [0], color=PALETTE["extratrees"], marker="s", lw=0.85, label="BindingDB"),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.54, 0.035),
        ncol=2,
        frameon=True,
        framealpha=0.92,
        facecolor="white",
        edgecolor="none",
        borderpad=0.35,
    )
    quality.to_csv(out_dir() / "SFig2_source_data_quality.csv", index=False)
    manifest.to_csv(out_dir() / "SFig2_source_data_manifest.csv", index=False)
    save_figure(fig, "SFig2_raw_data_quality")


def build_delta_heatmaps() -> tuple[list[str], list[str], dict[str, pd.DataFrame], float]:
    lead = read_extra("tkde_leaderboard_snapshot.csv")
    lead = lead.copy()
    lead["task"] = lead.apply(task_label, axis=1)
    sort_cols = ["benchmark", "source", "condition", "split_mode"]
    lead = lead.sort_values(sort_cols)
    tasks = lead.drop_duplicates("task")["task"].tolist()
    variants = [v for v in METHOD_STYLE if v in set(lead["variant"])]
    metrics = ["roc_auc_mean", "pr_auc_mean", "mcc_mean"]
    matrices: dict[str, pd.DataFrame] = {}
    all_values: list[float] = []
    for metric in metrics:
        pivot = lead.pivot_table(index="variant", columns="task", values=metric, aggfunc="mean")
        base = pivot.loc["ecfp_absdiff_xgb"] if "ecfp_absdiff_xgb" in pivot.index else 0
        delta = pivot.subtract(base, axis=1)
        delta = delta.reindex(index=variants, columns=tasks)
        matrices[metric] = delta
        vals = delta.to_numpy(dtype=float).ravel()
        all_values.extend(vals[np.isfinite(vals)].tolist())
    max_abs = float(np.nanpercentile(np.abs(all_values), 98)) if all_values else 0.1
    max_abs = max(0.06, min(0.18, max_abs))
    return variants, tasks, matrices, max_abs


def draw_sfig3() -> None:
    variants, tasks, matrices, max_abs = build_delta_heatmaps()
    metric_labels = [
        ("roc_auc_mean", "ROC-AUC delta"),
        ("pr_auc_mean", "PR-AUC delta"),
        ("mcc_mean", "MCC delta"),
    ]
    cmap = LinearSegmentedColormap.from_list("signed_pc", [PALETTE["negative"], "#FFFFFF", PALETTE["pc_cfrl"]])
    cmap.set_bad(PALETTE["missing"])
    fig = plt.figure(figsize=(10.8, 6.9))
    gs = fig.add_gridspec(3, 2, width_ratios=[1, 0.025], left=0.115, right=0.955, top=0.955, bottom=0.265, hspace=0.34, wspace=0.035)
    axes = []
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs)
    im = None
    for row, (metric, ylabel) in enumerate(metric_labels):
        ax = fig.add_subplot(gs[row, 0])
        mat = matrices[metric]
        im = ax.imshow(mat.to_numpy(dtype=float), aspect="auto", cmap=cmap, norm=norm)
        ax.set_yticks(np.arange(len(variants)))
        ax.set_yticklabels([METHOD_STYLE[v][0] for v in variants])
        ax.set_xticks(np.arange(len(tasks)))
        if row == 2:
            ax.set_xticklabels(tasks, rotation=90, ha="center", fontsize=6.7)
        else:
            ax.set_xticklabels([])
        ax.set_ylabel(ylabel)
        ax.set_xticks(np.arange(-0.5, len(tasks), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(variants), 1), minor=True)
        ax.grid(which="minor", color=PALETTE["grid"], linewidth=0.35)
        ax.tick_params(which="minor", bottom=False, left=False)
        set_box(ax)
        panel_label(ax, f"({chr(97 + row)})", x=-0.105, y=1.025)
        axes.append(ax)
    cax = fig.add_subplot(gs[:, 1])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("delta vs ECFP-XGB", fontsize=8.5)
    cbar.ax.tick_params(labelsize=7.6, width=0.4, length=2)
    source = []
    for metric, mat in matrices.items():
        long = mat.reset_index().melt(id_vars="variant", var_name="task", value_name="delta")
        long["metric"] = metric
        source.append(long)
    pd.concat(source, ignore_index=True).to_csv(out_dir() / "SFig3_source_data.csv", index=False)
    save_figure(fig, "SFig3_full_method_split_heatmap")


def aggregate_calibration(df: pd.DataFrame, experiment: str, variant: str) -> pd.DataFrame:
    sub = df[(df["experiment"] == experiment) & (df["variant"] == variant) & (df["n"] > 0)].copy()
    if sub.empty:
        return sub
    return (
        sub.groupby("bin", as_index=False)
        .agg(mean_prob=("mean_prob", "mean"), observed_rate=("observed_rate", "mean"), n=("n", "sum"))
        .sort_values("bin")
    )


def aggregate_risk(df: pd.DataFrame, experiment: str, variant: str) -> pd.DataFrame:
    sub = df[(df["experiment"] == experiment) & (df["variant"] == variant)].copy()
    if sub.empty:
        return sub
    sub["coverage_round"] = pd.to_numeric(sub["coverage"], errors="coerce").round(2)
    return (
        sub.groupby("coverage_round", as_index=False)
        .agg(risk=("risk", "mean"), positive_lift=("positive_lift", "mean"))
        .rename(columns={"coverage_round": "coverage"})
        .sort_values("coverage")
    )


def draw_sfig4() -> None:
    cal = read_csv("table_14a_p0_p1_calibration_bins.csv")
    risk = read_csv("table_14b_p0_p1_risk_coverage.csv")
    variants = ["ecfp_absdiff_xgb", "rich_diff_scalars", "rich_scalars_only", "extratrees_scalars"]
    fig, axes = plt.subplots(2, 2, figsize=(7.45, 5.65))
    axes = axes.ravel()
    plt.subplots_adjust(left=0.105, right=0.985, top=0.955, bottom=0.16, wspace=0.30, hspace=0.42)

    experiments = [("cross_database_transfer", "cross-database"), ("prospective_temporal", "historical release")]
    source_records = []
    for col, (experiment, _) in enumerate(experiments):
        ax = axes[col]
        ax.plot([0, 1], [0, 1], color=PALETTE["neutral"], lw=0.55, ls=(0, (2, 2)), zorder=0)
        for variant in variants:
            label, color, marker = METHOD_STYLE[variant]
            curve = aggregate_calibration(cal, experiment, variant)
            if curve.empty:
                continue
            ax.plot(curve["mean_prob"], curve["observed_rate"], color=color, marker=marker, ms=3.5, lw=0.85, label=label)
            tmp = curve.copy()
            tmp["experiment"] = experiment
            tmp["variant"] = variant
            tmp["panel"] = f"calibration_{experiment}"
            source_records.append(tmp)
        ax.set_xlabel("mean predicted probability")
        ax.set_ylabel("observed positive rate")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.grid(color=PALETTE["grid"], linewidth=0.35)
        set_box(ax)
        panel_label(ax, f"({chr(97 + col)})")

        ax_risk = axes[2 + col]
        for variant in variants:
            label, color, marker = METHOD_STYLE[variant]
            curve = aggregate_risk(risk, experiment, variant)
            if curve.empty:
                continue
            ax_risk.plot(curve["coverage"], curve["risk"], color=color, marker=marker, ms=3.5, lw=0.85, label=label)
            tmp = curve.copy()
            tmp["experiment"] = experiment
            tmp["variant"] = variant
            tmp["panel"] = f"risk_{experiment}"
            source_records.append(tmp)
        ax_risk.set_xlabel("retained coverage")
        ax_risk.set_ylabel("risk")
        ax_risk.set_xlim(0.08, 1.02)
        ax_risk.grid(color=PALETTE["grid"], linewidth=0.35)
        set_box(ax_risk)
        panel_label(ax_risk, f"({chr(99 + col)})")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.54, 0.035),
        ncol=4,
        frameon=True,
        framealpha=0.92,
        facecolor="white",
        edgecolor="none",
        borderpad=0.35,
    )
    if source_records:
        pd.concat(source_records, ignore_index=True, sort=False).to_csv(out_dir() / "SFig4_source_data.csv", index=False)
    save_figure(fig, "SFig4_calibration_risk_coverage")


def draw_condition_panel(
    ax: mpl.axes.Axes,
    df: pd.DataFrame,
    source: str,
    metric: str,
    show_xlabel: bool,
) -> None:
    order = ["all_pairs", "high_margin", "multi_source", "high_similarity", "equal_relation_rebuilt"]
    labels = ["all", "high\nmargin", "multi\nsource", "high\nsim.", "exact\nrelations"]
    sub = df[df["source"] == source].copy()
    for variant in ["ecfp_absdiff_xgb", "rich_diff_scalars", "rich_scalars_only"]:
        label, color, marker = METHOD_STYLE[variant]
        rows = []
        for cond in order:
            val = sub.loc[(sub["condition"] == cond) & (sub["variant"] == variant), metric]
            rows.append(float(val.iloc[0]) if len(val) else np.nan)
        ax.plot(np.arange(len(order)), rows, color=color, marker=marker, ms=3.7, lw=0.85, label=label)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(labels if show_xlabel else [])
    ax.set_ylabel(metric.replace("_", "-").upper() if metric == "mcc" else metric.replace("_", "-"))
    ax.grid(axis="y", color=PALETTE["grid"], linewidth=0.35)
    set_box(ax)


def draw_sfig5() -> None:
    sens = read_csv("table_07_noise_censoring_sensitivity.csv")
    sens = parse_series(sens, ["roc_auc", "pr_auc", "mcc", "balanced_accuracy", "f1", "brier", "ece_10"])
    fig, axes = plt.subplots(2, 3, figsize=(7.45, 5.25))
    plt.subplots_adjust(left=0.09, right=0.985, top=0.955, bottom=0.17, wspace=0.34, hspace=0.34)
    configs = [
        ("chembl_raw", "roc_auc"),
        ("chembl_raw", "pr_auc"),
        ("chembl_raw", "mcc"),
        ("bindingdb_raw", "roc_auc"),
        ("bindingdb_raw", "pr_auc"),
        ("bindingdb_raw", "mcc"),
    ]
    for idx, (ax, (source, metric)) in enumerate(zip(axes.ravel(), configs)):
        draw_condition_panel(ax, sens, source, metric, show_xlabel=idx >= 3)
        panel_label(ax, f"({chr(97 + idx)})", x=-0.12)
        if idx in (0, 3):
            ax.text(
                -0.34,
                0.5,
                BENCH_NAME.get(source, source),
                transform=ax.transAxes,
                rotation=90,
                ha="center",
                va="center",
                fontsize=9.4,
                color=PALETTE["border"],
            )
    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.54, 0.035),
        ncol=3,
        frameon=True,
        framealpha=0.92,
        facecolor="white",
        edgecolor="none",
        borderpad=0.35,
    )
    sens.to_csv(out_dir() / "SFig5_source_data.csv", index=False)
    save_figure(fig, "SFig5_sensitivity_curves")


def mol_image(smiles: str, highlight_atoms: list[int] | None = None, size: tuple[int, int] = (360, 260)) -> Image.Image:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(
        mol,
        size=size,
        kekulize=True,
        wedgeBonds=True,
        fitImage=True,
        highlightAtoms=highlight_atoms or [],
        highlightColor=(86 / 255, 180 / 255, 233 / 255),
    )
    return img.convert("RGBA")


def matched_atoms(smiles1: str, smiles2: str) -> tuple[list[int], list[int]]:
    m1 = Chem.MolFromSmiles(str(smiles1))
    m2 = Chem.MolFromSmiles(str(smiles2))
    if m1 is None or m2 is None:
        return [], []
    res = rdFMCS.FindMCS([m1, m2], timeout=2, ringMatchesRingOnly=True, completeRingsOnly=True)
    if not res.smartsString:
        return [], []
    patt = Chem.MolFromSmarts(res.smartsString)
    if patt is None:
        return [], []
    return list(m1.GetSubstructMatch(patt)), list(m2.GetSubstructMatch(patt))


def select_sfig6_cases() -> pd.DataFrame:
    df = read_extra("tkde_failure_case_studies.csv")
    df = df[df["variant"] == "rich_diff_scalars"].copy()
    choices = [
        "true_positive_high_confidence",
        "false_positive_high_confidence",
        "false_negative_low_score",
        "true_negative_low_score",
    ]
    rows = []
    for case_type in choices:
        sub = df[df["case_type"] == case_type].copy()
        sub["complexity"] = sub.apply(
            lambda r: (Chem.MolFromSmiles(str(r["smiles1"])).GetNumHeavyAtoms() if Chem.MolFromSmiles(str(r["smiles1"])) else 999)
            + (Chem.MolFromSmiles(str(r["smiles2"])).GetNumHeavyAtoms() if Chem.MolFromSmiles(str(r["smiles2"])) else 999),
            axis=1,
        )
        if "false_negative" in case_type or "true_negative" in case_type:
            sub = sub.sort_values(["source", "complexity", "probability"], ascending=[True, True, True])
        else:
            sub = sub.sort_values(["source", "complexity", "probability"], ascending=[True, True, False])
        rows.extend(sub.head(2).to_dict("records"))
    return pd.DataFrame(rows)


def draw_case_card(ax: mpl.axes.Axes, row: pd.Series, label: str) -> None:
    h1, h2 = matched_atoms(str(row["smiles1"]), str(row["smiles2"]))
    img1 = mol_image(str(row["smiles1"]), h1, size=(420, 250))
    img2 = mol_image(str(row["smiles2"]), h2, size=(420, 250))
    ax.imshow(np.asarray(img1), extent=(0.025, 0.475, 0.42, 0.95), aspect="auto")
    ax.imshow(np.asarray(img2), extent=(0.525, 0.975, 0.42, 0.95), aspect="auto")
    ax.annotate(
        "",
        xy=(0.515, 0.685),
        xytext=(0.485, 0.685),
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="->", lw=0.65, color=PALETTE["border"]),
    )
    case = str(row["case_type"]).replace("_", " ")
    target = shorten(str(row["target_pref_name"]), width=34, placeholder="...")
    metric = (
        f"{case} | p={float(row['probability']):.2f}  y={int(row['label'])}  "
        f"D={float(row['activity_delta']):.2f}  T={float(row['tanimoto']):.2f}"
    )
    ax.text(0.03, 0.135, target, transform=ax.transAxes, ha="left", va="bottom", fontsize=7.9)
    ax.text(0.03, 0.055, shorten(metric, width=76, placeholder="..."), transform=ax.transAxes, ha="left", va="bottom", fontsize=7.35)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    remove_box(ax)
    panel_label(ax, label, x=-0.015, y=0.97)


def draw_sfig6() -> None:
    cases = select_sfig6_cases()
    fig, axes = plt.subplots(4, 2, figsize=(7.45, 8.45))
    axes = axes.ravel()
    plt.subplots_adjust(left=0.045, right=0.985, top=0.975, bottom=0.095, wspace=0.11, hspace=0.20)
    for ax, (_, row), idx in zip(axes, cases.iterrows(), range(len(cases))):
        draw_case_card(ax, row, f"({chr(97 + idx)})")
    legend_handles = [Patch(facecolor=PALETTE["pair_only"], edgecolor="none", alpha=0.55, label="matched core")]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.54, 0.035),
        ncol=1,
        frameon=True,
        framealpha=0.92,
        facecolor="white",
        edgecolor="none",
        borderpad=0.35,
    )
    cases.to_csv(out_dir() / "SFig6_source_data.csv", index=False)
    save_figure(fig, "SFig6_rdkit_case_gallery", pdf_from_png=True)


def verify_main_palette() -> None:
    palette_path = out_dir() / "figure_palette_fixed.json"
    if palette_path.exists():
        recorded = json.loads(palette_path.read_text(encoding="utf-8"))
        if recorded != PALETTE:
            raise ValueError("Existing main-figure palette differs from supplementary palette.")
    (out_dir() / "figure_palette_fixed.json").write_text(json.dumps(PALETTE, indent=2), encoding="utf-8")


def main() -> None:
    apply_style()
    verify_main_palette()
    draw_sfig1()
    draw_sfig2()
    draw_sfig3()
    draw_sfig4()
    draw_sfig5()
    draw_sfig6()
    print(f"Wrote SFig.1--SFig.6 to {out_dir()}")


if __name__ == "__main__":
    main()
