# -*- coding: utf-8 -*-
"""Draw Fig.5--Fig.7 for the PC-CFRL TKDE manuscript.

Figure contracts following the local nature-figure workflow:
- Fig.5 conclusion: reliability improves when PC-CFRL is evaluated under
  selective and conformal prediction; all panels are line-curve summaries.
- Fig.6 conclusion: robustness gains are not uniform; a signed heatmap exposes
  where PC-CFRL is strong, competitive, or limited.
- Fig.7 conclusion: real RDKit molecular-pair cases ground the benchmark in
  chemical examples without implying a docking/3D mechanism.

Backend: Python/matplotlib only. RDKit is used only for 2D molecule rendering.
Export: white background, Times New Roman, PDF + 600 dpi PNG, fixed palette.
"""

from __future__ import annotations

from io import BytesIO
from pathlib import Path
import json
import math
import re
from textwrap import shorten

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Patch
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

METHOD = {
    "ecfp_absdiff_xgb": ("ECFP-XGB", PALETTE["ecfp"], "o"),
    "rich_diff_scalars": ("PC-CFRL", PALETTE["pc_cfrl"], "o"),
    "pccfrl_scalar_xgb": ("PC-CFRL", PALETTE["pc_cfrl"], "o"),
}

NAME_MAP = {
    "chembl_raw": "ChEMBL",
    "bindingdb_raw": "BindingDB",
    "family_scaffold_source_purged": "target-cluster+scaffold+source",
    "temporal_scaffold_source_purged": "temporal+source",
    "temporal_forward": "temporal",
    "target_family": "target cluster",
    "ligand_cold": "cold ligand",
    "pair_scaffold_group": "pair scaffold",
    "document_source_group": "doc source",
    "random": "random",
    "seen_like_ge_0p8": "seen-like",
    "mild_0p6_0p8": "mild",
    "moderate_0p4_0p6": "moderate",
    "severe_lt_0p4": "severe",
    "lohi_ligand_scaffold_source_purged": "LoHi ligand",
    "lohi_pair_scaffold_source_purged": "LoHi pair",
    "all_pairs": "all",
    "high_margin": "high margin",
    "multi_source": "multi-source",
    "high_similarity": "high sim.",
    "equal_relation_rebuilt": "exact relations only",
    "raw_chembl": "ChEMBL",
    "raw_bindingdb": "BindingDB",
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
    raise FileNotFoundError("Cannot locate remote table package.")


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


def read_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(table_dir() / name)


def read_extra(name: str) -> pd.DataFrame:
    return pd.read_csv(extra_dir() / name)


def results_1_4_dir() -> Path:
    root = data_root()
    candidates = [
        root / "remote_1_4_experiments_20260611" / "results",
        Path(r"C:\codex_tmp\paper15_tkde_remote\remote_1_4_experiments_20260611\results"),
    ]
    for cand in candidates:
        if cand.exists():
            return cand
    raise FileNotFoundError("Cannot locate remote_1_4_experiments_20260611/results.")


def results_5_7_dir() -> Path:
    root = data_root()
    candidates = [
        root / "remote_5_7_experiments_20260611" / "results",
        Path(r"C:\codex_tmp\paper15_tkde_remote\remote_5_7_experiments_20260611\results"),
    ]
    for cand in candidates:
        if cand.exists():
            return cand
    raise FileNotFoundError("Cannot locate remote_5_7_experiments_20260611/results.")


def read_results_1_4(name: str) -> pd.DataFrame:
    return pd.read_csv(results_1_4_dir() / name)


def read_results_5_7(name: str) -> pd.DataFrame:
    return pd.read_csv(results_5_7_dir() / name)


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


def aggregate_curve(df: pd.DataFrame, x_col: str, y_col: str, variant: str, group_filter: pd.Series) -> pd.DataFrame:
    sub = df[group_filter & (df["variant"] == variant)].copy()
    rows = []
    for _, row in sub.iterrows():
        y, _ = mean_std(row[y_col])
        x = float(row[x_col])
        rows.append({"x": x, "y": y})
    tmp = pd.DataFrame(rows)
    if tmp.empty:
        return tmp
    return tmp.groupby("x")["y"].agg(["mean", "min", "max"]).reset_index()


def line_with_band(ax: mpl.axes.Axes, curve: pd.DataFrame, label: str, color: str, marker: str) -> None:
    curve = curve.sort_values("x")
    ax.plot(curve["x"], curve["mean"], color=color, marker=marker, ms=4.0, lw=0.95, label=label, zorder=3)
    if {"min", "max"}.issubset(curve.columns):
        ax.fill_between(curve["x"].to_numpy(), curve["min"].to_numpy(), curve["max"].to_numpy(), color=color, alpha=0.12, lw=0, zorder=1)


def draw_fig5() -> None:
    missing = read_results_5_7("tkde_experiments_5_7_exp7_missing_source_trend.csv")
    risk = read_results_5_7("tkde_experiments_5_7_exp7_risk_coverage_overview.csv")
    variants = ["ecfp_absdiff_xgb", "pccfrl_scalar_xgb"]

    panels = [
        (
            missing,
            "missing_source_rate",
            "roc_auc_mean",
            missing["benchmark"].isin(["chembl_raw", "bindingdb_raw"]),
            "Missing source rate",
            "ROC-AUC",
        ),
        (
            missing,
            "missing_source_rate",
            "ece_10_mean",
            missing["benchmark"].isin(["chembl_raw", "bindingdb_raw"]),
            "Missing source rate",
            "ECE",
        ),
        (
            risk,
            "coverage",
            "selective_accuracy_mean",
            risk["benchmark"].isin(["chembl_raw", "bindingdb_raw"]) & (risk["missing_source_rate"].astype(float) == 0.0),
            "Retained coverage",
            "Retained-set accuracy",
        ),
        (
            risk,
            "coverage",
            "selective_accuracy_mean",
            risk["benchmark"].isin(["chembl_raw", "bindingdb_raw"]) & (risk["missing_source_rate"].astype(float) == 1.0),
            "Retained coverage",
            "Retained-set accuracy",
        ),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.35))
    axes = axes.ravel()
    plt.subplots_adjust(left=0.12, right=0.98, top=0.95, bottom=0.19, wspace=0.28, hspace=0.38)

    source_records = []
    for ax, (df, x_col, y_col, filt, xlabel, ylabel), label in zip(axes, panels, ["(a)", "(b)", "(c)", "(d)"]):
        for variant in variants:
            m_label, color, marker = METHOD[variant]
            curve = aggregate_curve(df, x_col, y_col, variant, filt)
            line_with_band(ax, curve, m_label, color, marker)
            if not curve.empty:
                c = curve.copy()
                c["panel"] = label
                c["metric"] = y_col
                c["variant"] = m_label
                source_records.append(c)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(color=PALETTE["grid"], linewidth=0.35)
        set_box(ax)
        panel_label(ax, label)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.55, 0.035),
        ncol=2,
        frameon=True,
        framealpha=0.92,
        facecolor="white",
        edgecolor="none",
        handlelength=1.4,
        borderpad=0.35,
    )
    pd.concat(source_records, ignore_index=True).to_csv(out_dir() / "Fig5_source_data.csv", index=False)
    save_figure(fig, "Fig5_reliability_curves")


def signed_heatmap(
    ax: mpl.axes.Axes,
    mat: pd.DataFrame,
    label: str,
    cbar_label: str,
    limit: float | None = None,
) -> None:
    arr = mat.to_numpy(dtype=float)
    if limit is None:
        limit = float(np.nanmax(np.abs(arr))) if np.isfinite(arr).any() else 1.0
    cmap = LinearSegmentedColormap.from_list("signed_pc", [PALETTE["negative"], "#FFFFFF", PALETTE["pc_cfrl"]])
    cmap.set_bad(PALETTE["missing"])
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-limit, vmax=limit)
    im = ax.imshow(np.ma.masked_invalid(arr), cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_xticklabels(mat.columns, rotation=35, ha="right", rotation_mode="anchor")
    ax.set_yticklabels(mat.index)
    ax.set_xticks(np.arange(-0.5, mat.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, mat.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.6)
    ax.tick_params(which="minor", bottom=False, left=False)
    set_box(ax)
    panel_label(ax, label)
    bbox = ax.get_position()
    cax = ax.figure.add_axes([bbox.x1 + 0.010, bbox.y0, 0.011, bbox.height])
    cb = plt.colorbar(im, cax=cax, orientation="vertical")
    cb.outline.set_linewidth(0.35)
    cb.ax.tick_params(labelsize=7.5, width=0.35, length=2, pad=1.2)
    cb.set_label(cbar_label, fontsize=7.8, labelpad=3, rotation=90)


def draw_fig6() -> None:
    shift = read_results_5_7("tkde_experiments_5_7_exp6_delta_vs_ecfp.csv")
    severity = read_results_5_7("tkde_experiments_5_7_exp6_severity_overview.csv")
    multitask = read_results_5_7("tkde_experiments_5_7_exp5_delta_vs_ecfp.csv")
    missing = read_results_5_7("tkde_experiments_5_7_exp7_missing_source_trend.csv")

    def num(value) -> float:
        return float(value) if not pd.isna(value) and str(value) != "" else np.nan

    def label_pair(bench: str, split: str) -> str:
        return f"{NAME_MAP.get(bench, bench)}\n{NAME_MAP.get(split, split)}"

    def summary_delta(df: pd.DataFrame, bench: str, split: str, variant: str, metric: str) -> float:
        pc = df[(df["benchmark"] == bench) & (df["split_mode"] == split) & (df["variant"] == variant)]
        base = df[(df["benchmark"] == bench) & (df["split_mode"] == split) & (df["variant"] == "ecfp_absdiff_xgb")]
        if pc.empty or base.empty:
            return np.nan
        return num(pc.iloc[0][f"{metric}_mean"]) - num(base.iloc[0][f"{metric}_mean"])

    shift_rows = [
        ("chembl_raw", "family_scaffold_source_purged"),
        ("chembl_raw", "temporal_forward"),
        ("chembl_raw", "ligand_cold"),
        ("bindingdb_raw", "family_scaffold_source_purged"),
        ("bindingdb_raw", "ligand_cold"),
        ("bindingdb_raw", "target_family"),
    ]
    mat_a = pd.DataFrame(
        [
            [
                summary_delta(shift, bench, split, "pccfrl_scalar_xgb", "roc_auc"),
                summary_delta(shift, bench, split, "pccfrl_scalar_xgb", "pr_auc"),
                summary_delta(shift, bench, split, "pccfrl_scalar_xgb", "mcc"),
            ]
            for bench, split in shift_rows
        ],
        index=[label_pair(bench, split) for bench, split in shift_rows],
        columns=[r"$\Delta$ROC", r"$\Delta$PR", r"$\Delta$MCC"],
    )

    severity_rows = [
        ("chembl_raw", "family_scaffold_source_purged"),
        ("chembl_raw", "temporal_forward"),
        ("bindingdb_raw", "family_scaffold_source_purged"),
        ("bindingdb_raw", "ligand_cold"),
    ]
    severity_cols = ["seen_like_ge_0p8", "mild_0p6_0p8", "moderate_0p4_0p6", "severe_lt_0p4"]
    sev_values = []
    for bench, split in severity_rows:
        vals = []
        for sev in severity_cols:
            pc = severity[
                (severity["benchmark"] == bench)
                & (severity["split_mode"] == split)
                & (severity["severity_bin"] == sev)
                & (severity["variant"] == "pccfrl_scalar_xgb")
            ]
            base = severity[
                (severity["benchmark"] == bench)
                & (severity["split_mode"] == split)
                & (severity["severity_bin"] == sev)
                & (severity["variant"] == "ecfp_absdiff_xgb")
            ]
            vals.append(num(pc.iloc[0]["roc_auc_mean"]) - num(base.iloc[0]["roc_auc_mean"]) if len(pc) and len(base) else np.nan)
        sev_values.append(vals)
    mat_b = pd.DataFrame(
        sev_values,
        index=[label_pair(bench, split) for bench, split in severity_rows],
        columns=[NAME_MAP.get(c, c) for c in severity_cols],
    )

    task_rows = [
        ("chembl_raw", "family_scaffold_source_purged"),
        ("chembl_raw", "target_family"),
        ("bindingdb_raw", "family_scaffold_source_purged"),
        ("bindingdb_raw", "target_family"),
    ]
    task_values = []
    for bench, split in task_rows:
        reg_pc = multitask[
            (multitask["experiment"] == "regression_delta")
            & (multitask["benchmark"] == bench)
            & (multitask["split_mode"] == split)
            & (multitask["variant"] == "pccfrl_scalar_xgb")
        ]
        reg_base = multitask[
            (multitask["experiment"] == "regression_delta")
            & (multitask["benchmark"] == bench)
            & (multitask["split_mode"] == split)
            & (multitask["variant"] == "ecfp_absdiff_xgb")
        ]
        dir_pc = multitask[
            (multitask["experiment"] == "direction_sign")
            & (multitask["benchmark"] == bench)
            & (multitask["split_mode"] == split)
            & (multitask["variant"] == "directed_pccfrl_xgb")
        ]
        dir_base = multitask[
            (multitask["experiment"] == "direction_sign")
            & (multitask["benchmark"] == bench)
            & (multitask["split_mode"] == split)
            & (multitask["variant"] == "directed_ecfp_xgb")
        ]
        mae_gain = num(reg_base.iloc[0]["mae_mean"]) - num(reg_pc.iloc[0]["mae_mean"]) if len(reg_pc) and len(reg_base) else np.nan
        spearman_gain = num(reg_pc.iloc[0]["spearman_mean"]) - num(reg_base.iloc[0]["spearman_mean"]) if len(reg_pc) and len(reg_base) else np.nan
        dir_gain = num(dir_pc.iloc[0]["roc_auc_mean"]) - num(dir_base.iloc[0]["roc_auc_mean"]) if len(dir_pc) and len(dir_base) else np.nan
        task_values.append([mae_gain, spearman_gain, dir_gain])
    mat_c = pd.DataFrame(
        task_values,
        index=[label_pair(bench, split) for bench, split in task_rows],
        columns=["MAE gain", "Spearman gain", "Dir. ROC gain"],
    )

    missing_rows = [
        ("chembl_raw", 0.0),
        ("chembl_raw", 0.5),
        ("chembl_raw", 1.0),
        ("bindingdb_raw", 0.0),
        ("bindingdb_raw", 0.5),
        ("bindingdb_raw", 1.0),
    ]
    miss_values = []
    for bench, rate in missing_rows:
        pc = missing[
            (missing["benchmark"] == bench)
            & (missing["missing_source_rate"].astype(float) == rate)
            & (missing["variant"] == "pccfrl_scalar_xgb")
        ]
        base = missing[
            (missing["benchmark"] == bench)
            & (missing["missing_source_rate"].astype(float) == rate)
            & (missing["variant"] == "ecfp_absdiff_xgb")
        ]
        miss_values.append(
            [
                num(pc.iloc[0]["roc_auc_mean"]) - num(base.iloc[0]["roc_auc_mean"]) if len(pc) and len(base) else np.nan,
                num(pc.iloc[0]["pr_auc_mean"]) - num(base.iloc[0]["pr_auc_mean"]) if len(pc) and len(base) else np.nan,
                num(base.iloc[0]["ece_10_mean"]) - num(pc.iloc[0]["ece_10_mean"]) if len(pc) and len(base) else np.nan,
            ]
        )
    mat_d = pd.DataFrame(
        miss_values,
        index=[f"{NAME_MAP.get(bench, bench)}\nmissing={rate:g}" for bench, rate in missing_rows],
        columns=[r"$\Delta$ROC", r"$\Delta$PR", "ECE gain"],
    )

    fig, axes = plt.subplots(2, 2, figsize=(7.45, 6.8))
    axes = axes.ravel()
    plt.subplots_adjust(left=0.15, right=0.93, top=0.96, bottom=0.12, wspace=0.70, hspace=0.42)
    signed_heatmap(axes[0], mat_a, "(a)", "delta vs ECFP", limit=0.12)
    signed_heatmap(axes[1], mat_b, "(b)", "ROC delta vs ECFP", limit=0.16)
    signed_heatmap(axes[2], mat_c, "(c)", "favorable delta", limit=0.08)
    signed_heatmap(axes[3], mat_d, "(d)", "delta vs ECFP", limit=0.12)
    pd.concat(
        {"fig6a_multisplit": mat_a, "fig6b_severity": mat_b, "fig6c_multitask": mat_c, "fig6d_missing_source": mat_d},
        names=["panel", "row"],
    ).to_csv(out_dir() / "Fig6_source_data.csv")
    save_figure(fig, "Fig6_robustness_boundary_heatmap")


def mol_image(smiles: str, highlight_atoms: list[int] | None = None, size: tuple[int, int] = (520, 260)) -> Image.Image:
    mol = Chem.MolFromSmiles(smiles)
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
    m1 = Chem.MolFromSmiles(smiles1)
    m2 = Chem.MolFromSmiles(smiles2)
    if m1 is None or m2 is None:
        return [], []
    res = rdFMCS.FindMCS([m1, m2], timeout=2, ringMatchesRingOnly=True, completeRingsOnly=True)
    if not res.smartsString:
        return [], []
    patt = Chem.MolFromSmarts(res.smartsString)
    if patt is None:
        return [], []
    match1 = list(m1.GetSubstructMatch(patt))
    match2 = list(m2.GetSubstructMatch(patt))
    return match1, match2


def select_cases() -> pd.DataFrame:
    df = read_extra("tkde_failure_case_studies.csv")
    choices = [
        "true_positive_high_confidence",
        "false_positive_high_confidence",
        "false_negative_low_score",
        "true_negative_low_score",
    ]
    rows = []
    for case_type in choices:
        sub = df[(df["case_type"] == case_type) & (df["variant"] == "rich_diff_scalars")].copy()
        scored = []
        for _, row in sub.iterrows():
            m1 = Chem.MolFromSmiles(str(row["smiles1"]))
            m2 = Chem.MolFromSmiles(str(row["smiles2"]))
            if m1 is None or m2 is None:
                continue
            total_atoms = m1.GetNumHeavyAtoms() + m2.GetNumHeavyAtoms()
            prob = float(row["probability"])
            prob_rank = prob if ("false_negative" in case_type or "true_negative" in case_type) else -prob
            scored.append((total_atoms, prob_rank, row))
        if scored:
            rows.append(sorted(scored, key=lambda x: (x[0], x[1]))[0][2])
    return pd.DataFrame(rows)


def draw_case_card(ax: mpl.axes.Axes, row: pd.Series, label: str) -> None:
    h1, h2 = matched_atoms(str(row["smiles1"]), str(row["smiles2"]))
    img1 = mol_image(str(row["smiles1"]), h1)
    img2 = mol_image(str(row["smiles2"]), h2)
    ax.imshow(np.asarray(img1), extent=(0.03, 0.47, 0.43, 0.93), aspect="auto")
    ax.imshow(np.asarray(img2), extent=(0.53, 0.97, 0.43, 0.93), aspect="auto")
    ax.annotate(
        "",
        xy=(0.515, 0.68),
        xytext=(0.485, 0.68),
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="->", lw=0.65, color=PALETTE["border"]),
    )
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    metric = (
        f"p={float(row['probability']):.2f}  y={int(row['label'])}  "
        f"Delta={float(row['activity_delta']):.2f}  T={float(row['tanimoto']):.2f}"
    )
    target = shorten(str(row["target_pref_name"]), width=34, placeholder="...")
    ax.text(0.03, 0.145, target, transform=ax.transAxes, ha="left", va="bottom", fontsize=8.1)
    ax.text(0.03, 0.055, metric, transform=ax.transAxes, ha="left", va="bottom", fontsize=7.9)
    ax.set_xticks([])
    ax.set_yticks([])
    remove_box(ax)
    panel_label(ax, label)


def draw_fig7() -> None:
    cases = select_cases()
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.9))
    axes = axes.ravel()
    plt.subplots_adjust(left=0.05, right=0.98, top=0.96, bottom=0.16, wspace=0.14, hspace=0.24)
    for ax, (_, row), label in zip(axes, cases.iterrows(), ["(a)", "(b)", "(c)", "(d)"]):
        draw_case_card(ax, row, label)
    legend_handles = [Patch(facecolor=PALETTE["pair_only"], edgecolor="none", alpha=0.55, label="matched core")]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.53, 0.035),
        ncol=1,
        frameon=True,
        framealpha=0.92,
        facecolor="white",
        edgecolor="none",
        borderpad=0.35,
    )
    cases.to_csv(out_dir() / "Fig7_source_data.csv", index=False)
    save_figure(fig, "Fig7_rdkit_case_cards", pdf_from_png=True)


def main() -> None:
    apply_style()
    (out_dir() / "figure_palette_fixed.json").write_text(json.dumps(PALETTE, indent=2), encoding="utf-8")
    draw_fig5()
    draw_fig6()
    draw_fig7()
    print(f"Wrote Fig.5--Fig.7 to {out_dir()}")


if __name__ == "__main__":
    main()
