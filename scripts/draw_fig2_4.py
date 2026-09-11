# -*- coding: utf-8 -*-
"""Draw Fig.2--Fig.4 for the PC-CFRL-Bench manuscript.

Figure contract, following the local nature-figure workflow:
- Fig.2 conclusion: target-cluster-balanced splitting sharply reduces pair,
  ligand, and scaffold overlap without materially changing class balance.
- Fig.3 conclusion: PC-CFRL is consistently stronger than ECFP under hard-OOD
  settings, while remaining visually separable from pair-only/scalar controls.
- Fig.4 conclusion: PC-CFRL improves over ECFP under cross-database and
  prospective temporal shifts; paired dumbbells show the deployment delta.

Backend: Python/matplotlib only.
Export: white background, Times New Roman, PDF + 600 dpi PNG, fixed palette.
"""

from __future__ import annotations

from pathlib import Path
import json
import math
import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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

METHOD_STYLE = {
    "ECFP-XGB": {"color": PALETTE["ecfp"], "marker": "o"},
    "PC-CFRL": {"color": PALETTE["pc_cfrl"], "marker": "o"},
    "pair-only": {"color": PALETTE["pair_only"], "marker": "s"},
    "scalar-only": {"color": PALETTE["scalar_only"], "marker": "D"},
}

NAME_MAP = {
    "ecfp_absdiff_xgb": "ECFP-XGB",
    "rich_diff_scalars": "PC-CFRL",
    "pccfrl_scalar_xgb": "PC-CFRL",
    "rich_pair_only": "pair-only",
    "rich_scalars_only": "scalar-only",
    "extratrees_scalars": "ExtraTrees",
    "family_scaffold_purged": "family+scaffold",
    "family_scaffold_source_purged": "source+scaffold",
    "temporal_scaffold_source_purged": "temporal+source",
    "bindingdb_to_chembl_scaffold_purged": "BDB->ChEMBL\nscaffold",
    "chembl_to_bindingdb_scaffold_purged": "ChEMBL->BDB\nscaffold",
    "chembl_future_after_2005_scaffold_source_purged": "ChEMBL\nfuture>2005",
    "chembl_future_after_2010_scaffold_source_purged": "ChEMBL\nfuture>2010",
    "chembl_future_after_2015_scaffold_source_purged": "ChEMBL\nfuture>2015",
    "release_new_evidence": "ChEMBL36-37\nnew evidence",
    "release_scaffold_source_purged": "ChEMBL36-37\nsource purged",
    "release_cold_ligand_both": "ChEMBL36-37\ncold ligand",
    "chembl_raw": "ChEMBL",
    "bindingdb_raw": "BindingDB",
}


def project_root() -> Path:
    here = Path(__file__).resolve()
    configured = os.environ.get("PC_CFRL_DATA_ROOT")
    candidates = ([Path(configured).expanduser()] if configured else []) + [
        here.parent.parent / "data_remote",
    ]
    for cand in candidates:
        if (cand / "tkde_experiment_tables_2026_05_20").exists():
            return cand
    raise FileNotFoundError(
        "Set PC_CFRL_DATA_ROOT to the reconstructed experiment-data directory."
    )


def table_dir() -> Path:
    return project_root() / "tkde_experiment_tables_2026_05_20"


def extra_dir() -> Path:
    return project_root() / "extra_csv"


def out_dir() -> Path:
    out = Path(__file__).resolve().parent
    out.mkdir(parents=True, exist_ok=True)
    return out


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
        -0.10,
        1.045,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        fontweight="bold",
        color=PALETTE["border"],
    )


def parse_mean_std(value) -> tuple[float, float]:
    if pd.isna(value):
        return np.nan, np.nan
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value), np.nan
    text = str(value)
    match = re.match(r"\s*([+-]?\d+(?:\.\d+)?)\s*\+/-\s*([+-]?\d+(?:\.\d+)?)", text)
    if match:
        return float(match.group(1)), float(match.group(2))
    match = re.match(r"\s*([+-]?\d+(?:\.\d+)?)", text)
    if match:
        return float(match.group(1)), np.nan
    return np.nan, np.nan


def read_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(table_dir() / name)


def read_extra(name: str) -> pd.DataFrame:
    return pd.read_csv(extra_dir() / name)


def results_1_4_dir() -> Path:
    root = project_root()
    candidates = [
        root / "remote_1_4_experiments_20260611" / "results",
        root / "results",
    ]
    for cand in candidates:
        if cand.exists():
            return cand
    raise FileNotFoundError("Could not locate remote_1_4_experiments_20260611/results")


def read_results_1_4(name: str) -> pd.DataFrame:
    return pd.read_csv(results_1_4_dir() / name)


def results_dir() -> Path:
    root = project_root()
    candidates = [
        root / "results",
        root / "15PC-CFRL" / "results",
    ]
    for cand in candidates:
        if cand.exists():
            return cand
    raise FileNotFoundError("Could not locate 15PC-CFRL/results")


def read_result(name: str) -> pd.DataFrame:
    return pd.read_csv(results_dir() / name)


def read_fold_leakage_audit() -> pd.DataFrame:
    root = project_root()
    candidates = [
        root / "results" / "leakage_and_split_audit.csv",
        root / "15PC-CFRL" / "results" / "leakage_and_split_audit.csv",
        root / "extra_csv" / "leakage_and_split_audit.csv",
    ]
    for path in candidates:
        if path.exists():
            return pd.read_csv(path)
    raise FileNotFoundError("Could not locate leakage_and_split_audit.csv")


def leakage_box_panel(
    ax: mpl.axes.Axes,
    audit: pd.DataFrame,
    metric: str,
    ylabel: str,
    label: str,
    ylim: tuple[float, float],
    panel_index: int,
) -> None:
    modes = ["random", "target_cluster_balanced"]
    edge_colors = [PALETTE["ecfp"], PALETTE["pc_cfrl"]]
    fill_colors = [PALETTE["light_neutral"], PALETTE["pair_only"]]

    for pos, (mode, edge, fill) in enumerate(zip(modes, edge_colors, fill_colors)):
        values = audit.loc[audit["fold_mode"] == mode, metric].dropna().to_numpy(dtype=float)
        box = ax.boxplot(
            [values],
            positions=[pos],
            widths=0.43,
            patch_artist=True,
            showfliers=False,
            whis=1.5,
            boxprops={"facecolor": fill, "edgecolor": edge, "linewidth": 0.75, "alpha": 0.32},
            medianprops={"color": edge, "linewidth": 1.15},
            whiskerprops={"color": edge, "linewidth": 0.65},
            capprops={"color": edge, "linewidth": 0.65},
        )
        box["boxes"][0].set_zorder(2)

        rng = np.random.default_rng(2026 + panel_index * 10 + pos)
        jitter = rng.uniform(-0.115, 0.115, size=len(values))
        ax.scatter(
            pos + jitter,
            values,
            s=15,
            color=edge,
            alpha=0.72,
            edgecolor="white",
            linewidth=0.28,
            zorder=3,
        )
        ax.scatter(
            pos,
            float(np.mean(values)),
            marker="D",
            s=29,
            color=PALETTE["scalar_only"],
            edgecolor=PALETTE["border"],
            linewidth=0.45,
            zorder=4,
        )

    counts = audit.groupby("fold_mode")[metric].count()
    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        [
            f"Random\nn={int(counts['random'])}",
            f"Target-cluster\nbalanced, n={int(counts['target_cluster_balanced'])}",
        ]
    )
    ax.set_ylabel(ylabel)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(*ylim)
    ax.yaxis.grid(True, color=PALETTE["grid"], linewidth=0.48)
    ax.set_axisbelow(True)
    set_box(ax)
    panel_label(ax, label)


def draw_fig2() -> None:
    audit = read_fold_leakage_audit()
    audit = audit[audit["fold_mode"].isin(["random", "target_cluster_balanced"])].copy()

    panels = [
        ("pair_overlap_rate", "Pair overlap rate", "(a)", (0.0, 1.04)),
        ("any_ligand_overlap_rate", "Any-ligand overlap rate", "(b)", (0.0, 1.04)),
        ("pair_scaffold_overlap_rate", "Pair-scaffold overlap rate", "(c)", (0.0, 1.04)),
        ("test_positive_rate", "Positive-pair rate", "(d)", (0.0, 0.28)),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.55))
    axes = axes.ravel()
    plt.subplots_adjust(left=0.115, right=0.975, top=0.955, bottom=0.20, wspace=0.28, hspace=0.42)
    for idx, (ax, (metric, ylabel, label, ylim)) in enumerate(zip(axes, panels)):
        leakage_box_panel(ax, audit, metric, ylabel, label, ylim, idx)

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markersize=5,
            markerfacecolor=PALETTE["ecfp"],
            markeredgecolor="white",
            markeredgewidth=0.35,
            label="Random split",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markersize=5,
            markerfacecolor=PALETTE["pc_cfrl"],
            markeredgecolor="white",
            markeredgewidth=0.35,
            label="Target-cluster-balanced split",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            linestyle="none",
            markersize=4.8,
            markerfacecolor=PALETTE["scalar_only"],
            markeredgecolor=PALETTE["border"],
            markeredgewidth=0.45,
            label="Mean",
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.54, 0.035),
        ncol=3,
        frameon=True,
        framealpha=0.92,
        facecolor="white",
        edgecolor="none",
        handletextpad=0.45,
        columnspacing=1.4,
        borderpad=0.35,
    )
    fig.patch.set_facecolor("white")

    source_rows = []
    panel_names = ["fig2a_pair_overlap", "fig2b_ligand_overlap", "fig2c_scaffold_overlap", "fig2d_positive_rate"]
    for panel_name, (metric, ylabel, _, _) in zip(panel_names, panels):
        for _, row in audit.iterrows():
            source_rows.append(
                {
                    "panel": panel_name,
                    "metric": metric,
                    "display_metric": ylabel,
                    "seed": row["seed"],
                    "fold": row["fold"],
                    "fold_mode": row["fold_mode"],
                    "value": row[metric],
                    "source_artifact": "results/leakage_and_split_audit.csv",
                    "audit_scope": "ACNet random vs target-cluster-balanced fold leakage audit",
                }
            )
    pd.DataFrame(source_rows).to_csv(out_dir() / "Fig2_source_data.csv", index=False)
    save_figure(fig, "Fig2_data_quality_leakage_matrix")


def paired_delta_ci(delta: float, std: float, n: int) -> tuple[float, float]:
    if not np.isfinite(delta):
        return np.nan, np.nan
    if not np.isfinite(std) or n <= 1:
        return delta, delta
    # n=25 is the fixed seed-fold count for the main paired-delta audits.
    tcrit = 2.0639 if n == 25 else 1.96
    half_width = tcrit * std / math.sqrt(n)
    return delta - half_width, delta + half_width


def collect_hard_ood_effects() -> pd.DataFrame:
    rows = []

    def add_delta_file(
        filename: str,
        dataset: str,
        split_mode: str,
        split_label: str,
        variants: list[str],
    ) -> None:
        df = read_result(filename)
        sub = df[
            (df["split_mode"] == split_mode)
            & (df["baseline"] == "ecfp_absdiff_xgb")
            & (df["variant"].isin(variants))
            & (df["metric"].isin(["roc_auc", "pr_auc", "mcc"]))
        ]
        for _, row in sub.iterrows():
            n = int(float(row["n_pairs"]))
            delta = float(row["delta_mean"])
            ci_low, ci_high = paired_delta_ci(delta, float(row["delta_std"]), n)
            rows.append(
                {
                    "dataset": dataset,
                    "split": split_label,
                    "group": f"{dataset}\n{split_label}",
                    "variant": NAME_MAP.get(row["variant"], row["variant"]),
                    "baseline": "ECFP-XGB",
                    "metric": row["metric"],
                    "n": n,
                    "delta": delta,
                    "ci95_low": ci_low,
                    "ci95_high": ci_high,
                    "win_rate": float(row.get("win_rate_favorable", np.nan)),
                    "p_value": float(row.get("wilcoxon_p", np.nan)),
                    "ci_method": "paired t 95% CI",
                }
            )

    def add_audit_rows(
        experiment: str,
        dataset: str,
        split_mode: str,
        split_label: str,
        variants: list[str],
    ) -> None:
        df = read_result("tkde_statistical_effects_audit.csv")
        sub = df[
            (df["experiment"] == experiment)
            & (df["split_mode"] == split_mode)
            & (df["baseline"] == "ecfp_absdiff_xgb")
            & (df["variant"].isin(variants))
            & (df["metric"].isin(["roc_auc", "pr_auc", "mcc"]))
        ]
        for _, row in sub.iterrows():
            rows.append(
                {
                    "dataset": dataset,
                    "split": split_label,
                    "group": f"{dataset}\n{split_label}",
                    "variant": NAME_MAP.get(row["variant"], row["variant"]),
                    "baseline": "ECFP-XGB",
                    "metric": row["metric"],
                    "n": int(float(row["n"])),
                    "delta": float(row["mean_delta_favorable"]),
                    "ci95_low": float(row["ci95_low"]),
                    "ci95_high": float(row["ci95_high"]),
                    "win_rate": float(row.get("win_rate", np.nan)),
                    "p_value": float(row.get("holm_p", np.nan)),
                    "ci_method": "bootstrap 95% CI",
                }
            )

    add_delta_file(
        "acnet_target_family_best_joint_seed0_4_allfold_n10_paired_deltas.csv",
        "ACNet",
        "family_scaffold_purged",
        "family+scaffold",
        ["rich_pair_only", "rich_diff_scalars"],
    )
    add_delta_file(
        "moleculeace_target_family_pairs_seed0_4_n120_paired_deltas.csv",
        "MoleculeACE",
        "family_scaffold_purged",
        "family+scaffold",
        ["rich_scalars_only", "rich_diff_scalars"],
    )
    add_audit_rows(
        "raw_chembl",
        "ChEMBL",
        "family_scaffold_source_purged",
        "source+scaffold",
        ["rich_scalars_only", "rich_diff_scalars"],
    )
    add_audit_rows(
        "raw_chembl",
        "ChEMBL",
        "temporal_scaffold_source_purged",
        "temporal+source",
        ["rich_scalars_only", "rich_diff_scalars"],
    )
    add_audit_rows(
        "raw_bindingdb",
        "BindingDB",
        "family_scaffold_source_purged",
        "source+scaffold",
        ["rich_scalars_only", "rich_diff_scalars"],
    )
    df = pd.DataFrame(rows)
    group_order = {
        "ACNet\nfamily+scaffold": 0,
        "MoleculeACE\nfamily+scaffold": 1,
        "ChEMBL\nsource+scaffold": 2,
        "ChEMBL\ntemporal+source": 3,
        "BindingDB\nsource+scaffold": 4,
    }
    metric_order = {"roc_auc": 0, "pr_auc": 1, "mcc": 2}
    method_order = {"pair-only": 0, "scalar-only": 1, "PC-CFRL": 2}
    df["method_order"] = df["variant"].map(method_order)
    df["group_order"] = df["group"].map(group_order)
    df["metric_order"] = df["metric"].map(metric_order)
    return df.sort_values(["group_order", "metric_order", "method_order"]).drop(
        columns=["group_order", "metric_order", "method_order"]
    )


def paired_delta_forest_panel(
    ax: mpl.axes.Axes,
    df: pd.DataFrame,
    metric: str,
    groups: list[str],
    variants: list[str],
    label: str,
    xlim: tuple[float, float],
) -> None:
    sub = df[df["metric"] == metric]
    ybase = np.arange(len(groups))[::-1]
    offsets = {"pair-only": -0.16, "scalar-only": 0.0, "PC-CFRL": 0.16}
    for variant in variants:
        vdf = sub[sub["variant"] == variant]
        xs, ys, xerr_low, xerr_high = [], [], [], []
        for i, group in enumerate(groups):
            row = vdf[vdf["group"] == group]
            if len(row):
                delta = float(row.iloc[0]["delta"])
                xs.append(delta)
                xerr_low.append(delta - float(row.iloc[0]["ci95_low"]))
                xerr_high.append(float(row.iloc[0]["ci95_high"]) - delta)
            else:
                xs.append(np.nan)
                xerr_low.append(np.nan)
                xerr_high.append(np.nan)
            ys.append(ybase[i] + offsets[variant])
        style = METHOD_STYLE[variant]
        ax.errorbar(
            xs,
            ys,
            xerr=np.vstack([xerr_low, xerr_high]),
            fmt=style["marker"],
            ms=4.1,
            color=style["color"],
            ecolor=style["color"],
            elinewidth=0.60,
            capsize=1.8,
            capthick=0.60,
            markerfacecolor=style["color"],
            markeredgecolor="white",
            markeredgewidth=0.35,
            label=variant,
            zorder=3,
        )
    ax.set_yticks(ybase)
    ax.set_yticklabels(groups)
    ax.axvline(0, color=PALETTE["neutral"], linewidth=0.55, linestyle=(0, (2.0, 2.0)), zorder=1)
    ax.grid(axis="x", color=PALETTE["grid"], linewidth=0.35)
    set_box(ax)
    panel_label(ax, label)
    labels = {"roc_auc": r"$\Delta$ROC-AUC", "pr_auc": r"$\Delta$PR-AUC", "mcc": r"$\Delta$MCC"}
    ax.set_xlabel(f"{labels[metric]} vs. ECFP-XGB")
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.7, len(groups) - 0.3)
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, _pos: "0" if abs(x) < 0.0005 else f"{x:+.2f}"))


def draw_fig3() -> None:
    df = collect_hard_ood_effects()
    groups = [
        "ACNet\nfamily+scaffold",
        "MoleculeACE\nfamily+scaffold",
        "ChEMBL\nsource+scaffold",
        "ChEMBL\ntemporal+source",
        "BindingDB\nsource+scaffold",
    ]
    variants = ["pair-only", "scalar-only", "PC-CFRL"]
    ci_bounds = df[["ci95_low", "ci95_high"]].to_numpy(dtype=float)
    xlow = min(-0.02, float(np.nanmin(ci_bounds)) - 0.012)
    xhigh = max(0.16, float(np.nanmax(ci_bounds)) + 0.012)
    fig, axes = plt.subplots(1, 3, figsize=(7.25, 3.55), sharey=True)
    plt.subplots_adjust(left=0.205, right=0.985, top=0.94, bottom=0.30, wspace=0.14)
    for ax, metric, label in zip(axes, ["roc_auc", "pr_auc", "mcc"], ["(a)", "(b)", "(c)"]):
        paired_delta_forest_panel(ax, df, metric, groups, variants, label, (xlow, xhigh))
    handles = [
        Line2D(
            [0],
            [0],
            color=METHOD_STYLE[name]["color"],
            marker=METHOD_STYLE[name]["marker"],
            linestyle="None",
            markersize=4.6,
            markerfacecolor=METHOD_STYLE[name]["color"],
            markeredgecolor="white",
            markeredgewidth=0.35,
        )
        for name in variants
    ]
    fig.legend(
        handles,
        variants,
        loc="lower center",
        bbox_to_anchor=(0.60, 0.035),
        ncol=3,
        frameon=True,
        framealpha=0.92,
        facecolor="white",
        edgecolor="none",
        handlelength=1.1,
        borderpad=0.35,
        labelspacing=0.3,
    )
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
    df.to_csv(out_dir() / "Fig3_source_data.csv", index=False)
    save_figure(fig, "Fig3_hard_ood_forest")


def collect_transfer() -> pd.DataFrame:
    df = read_csv("table_05_crossdb_prospective_topk.csv")
    keep_splits = [
        "bindingdb_to_chembl_scaffold_purged",
        "chembl_to_bindingdb_scaffold_purged",
    ]
    keep_vars = ["ecfp_absdiff_xgb", "rich_diff_scalars"]
    sub = df[df["split_mode"].isin(keep_splits) & df["variant"].isin(keep_vars)].copy()
    rows = []
    for _, row in sub.iterrows():
        for metric in ["roc_auc", "pr_auc", "mcc", "ef_5pct"]:
            mean, std = parse_mean_std(row[metric])
            rows.append(
                {
                    "split": NAME_MAP.get(row["split_mode"], row["split_mode"]),
                    "variant": NAME_MAP.get(row["variant"], row["variant"]),
                    "metric": metric,
                    "mean": mean,
                    "std": std,
                }
            )
    release = read_results_1_4("chembl36_to_37_release_compact_summary.csv")
    release_splits = [
        "release_new_evidence",
        "release_scaffold_source_purged",
        "release_cold_ligand_both",
    ]
    release_vars = ["ecfp_absdiff_xgb", "pccfrl_scalar_xgb"]
    rel = release[release["split_mode"].isin(release_splits) & release["variant"].isin(release_vars)].copy()
    for _, row in rel.iterrows():
        for metric in ["roc_auc", "pr_auc", "mcc"]:
            rows.append(
                {
                    "split": NAME_MAP.get(row["split_mode"], row["split_mode"]),
                    "variant": NAME_MAP.get(row["variant"], row["variant"]),
                    "metric": metric,
                    "mean": float(row[f"{metric}_mean"]),
                    "std": float(row[f"{metric}_std"]),
                }
            )
        rows.append(
            {
                "split": NAME_MAP.get(row["split_mode"], row["split_mode"]),
                "variant": NAME_MAP.get(row["variant"], row["variant"]),
                "metric": "ef_5pct",
                "mean": np.nan,
                "std": np.nan,
            }
        )
    return pd.DataFrame(rows)


def dumbbell_panel(
    ax: mpl.axes.Axes,
    df: pd.DataFrame,
    metric: str,
    splits: list[str],
    label: str,
) -> None:
    sub = df[df["metric"] == metric]
    y = np.arange(len(splits))[::-1]
    xs_ecfp, xs_pc, err_ecfp, err_pc = [], [], [], []
    reference_lines = {"roc_auc": 0.5, "mcc": 0.0, "ef_5pct": 1.0}
    for split in splits:
        e = sub[(sub["split"] == split) & (sub["variant"] == "ECFP-XGB")].iloc[0]
        p = sub[(sub["split"] == split) & (sub["variant"] == "PC-CFRL")].iloc[0]
        xs_ecfp.append(e["mean"])
        xs_pc.append(p["mean"])
        err_ecfp.append(e["std"])
        err_pc.append(p["std"])
    xs_ecfp_arr = np.array(xs_ecfp, dtype=float)
    xs_pc_arr = np.array(xs_pc, dtype=float)
    err_ecfp_arr = np.array(err_ecfp, dtype=float)
    err_pc_arr = np.array(err_pc, dtype=float)
    finite_ecfp = np.isfinite(xs_ecfp_arr)
    finite_pc = np.isfinite(xs_pc_arr)
    finite_pair = finite_ecfp & finite_pc
    if metric in reference_lines:
        ax.axvline(
            reference_lines[metric],
            color=PALETTE["neutral"],
            linewidth=0.48,
            linestyle=(0, (2.0, 2.0)),
            alpha=0.80,
            zorder=0,
        )
    for yi, x0, x1 in zip(y[finite_pair], xs_ecfp_arr[finite_pair], xs_pc_arr[finite_pair]):
        ax.plot([x0, x1], [yi, yi], color=PALETTE["light_neutral"], linewidth=1.15, zorder=1)
    ax.errorbar(
        xs_ecfp_arr[finite_ecfp],
        y[finite_ecfp],
        xerr=err_ecfp_arr[finite_ecfp],
        fmt="o",
        ms=3.8,
        color=PALETTE["ecfp"],
        ecolor=PALETTE["ecfp"],
        elinewidth=0.50,
        capsize=1.4,
        markerfacecolor=PALETTE["ecfp"],
        markeredgecolor="white",
        markeredgewidth=0.35,
        label="ECFP-XGB",
        zorder=3,
    )
    ax.errorbar(
        xs_pc_arr[finite_pc],
        y[finite_pc],
        xerr=err_pc_arr[finite_pc],
        fmt="o",
        ms=4.2,
        color=PALETTE["pc_cfrl"],
        ecolor=PALETTE["pc_cfrl"],
        elinewidth=0.55,
        capsize=1.4,
        markerfacecolor=PALETTE["pc_cfrl"],
        markeredgecolor="white",
        markeredgewidth=0.35,
        label="PC-CFRL",
        zorder=4,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(splits)
    ax.grid(axis="x", color=PALETTE["grid"], linewidth=0.35)
    set_box(ax)
    panel_label(ax, label)
    labels = {"roc_auc": "ROC-AUC", "pr_auc": "PR-AUC", "mcc": "MCC", "ef_5pct": "EF@5%"}
    ax.set_xlabel(labels[metric])
    values = np.array(xs_ecfp + xs_pc + ([reference_lines[metric]] if metric in reference_lines else []), dtype=float)
    values = values[np.isfinite(values)]
    low, high = np.nanmin(values), np.nanmax(values)
    span = high - low if high > low else 1.0
    ax.set_xlim(low - span * 0.22 - 0.005, high + span * 0.22 + 0.005)
    ax.set_ylim(-0.7, len(splits) - 0.3)


def draw_fig4() -> None:
    df = collect_transfer()
    splits = [
        "BDB->ChEMBL\nscaffold",
        "ChEMBL->BDB\nscaffold",
        "ChEMBL36-37\nnew evidence",
        "ChEMBL36-37\nsource purged",
        "ChEMBL36-37\ncold ligand",
    ]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.35), sharey=True)
    axes = axes.ravel()
    plt.subplots_adjust(left=0.20, right=0.97, top=0.95, bottom=0.19, wspace=0.18, hspace=0.34)
    for ax, metric, label in zip(axes, ["roc_auc", "pr_auc", "mcc", "ef_5pct"], ["(a)", "(b)", "(c)", "(d)"]):
        dumbbell_panel(ax, df, metric, splits, label)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.58, 0.035),
        ncol=2,
        frameon=True,
        framealpha=0.92,
        facecolor="white",
        edgecolor="none",
        handlelength=1.2,
        borderpad=0.35,
        labelspacing=0.3,
    )
    for ax in [axes[1], axes[3]]:
        ax.tick_params(labelleft=False)
    df.to_csv(out_dir() / "Fig4_source_data.csv", index=False)
    save_figure(fig, "Fig4_crossdb_prospective_dumbbell")


def main() -> None:
    apply_style()
    (out_dir() / "figure_palette_fixed.json").write_text(json.dumps(PALETTE, indent=2), encoding="utf-8")
    draw_fig2()
    draw_fig3()
    draw_fig4()
    print(f"Wrote figures to {out_dir()}")


if __name__ == "__main__":
    main()
