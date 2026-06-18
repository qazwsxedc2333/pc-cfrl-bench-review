from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PALETTE = {
    "pair_prior": "#999999",
    "old": "#D55E00",
    "rich_pair": "#0072B2",
    "rich_target": "#009E73",
    "shuffled": "#CC79A7",
    "neural": "#E69F00",
}


def save(fig: plt.Figure, outdir: Path, name: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(outdir / f"{name}.png", dpi=300)
    fig.savefig(outdir / f"{name}.pdf")
    plt.close(fig)


def load_csv(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p)


def plot_threshold(summary: pd.DataFrame, outdir: Path) -> None:
    if summary.empty:
        return
    threshold_col = "threshold" if "threshold" in summary else "cliff_threshold"
    keep = {
        "pair_mean_cliff_score": ("Pair prior", PALETTE["pair_prior"]),
        "classifier_blend_rich_pair_only": ("Rich pair", PALETTE["rich_pair"]),
        "classifier_blend_rich_interaction_saprot": ("Rich target", PALETTE["rich_target"]),
        "classifier_blend_rich_interaction_saprot_shuffled": ("Shuffled target", PALETTE["shuffled"]),
    }
    fig, ax = plt.subplots(figsize=(5.4, 3.2))
    for variant, (label, color) in keep.items():
        sub = summary[summary["variant"] == variant].sort_values(threshold_col)
        if sub.empty:
            continue
        ax.errorbar(
            sub[threshold_col],
            sub["cliff_auc_mean"],
            yerr=sub.get("cliff_auc_std"),
            marker="o",
            linewidth=2.0,
            capsize=3,
            label=label,
            color=color,
        )
    ax.set_xlabel("Cliff threshold")
    ax.set_ylabel("Cliff AUC")
    ax.set_title("Threshold sensitivity")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, ncol=2, fontsize=8)
    save(fig, outdir, "fig_rich_pair_threshold_sensitivity")


def plot_ood(summary: pd.DataFrame, outdir: Path) -> None:
    if summary.empty:
        return
    variants = [
        ("classifier_cliff_score_rich_pair_only", "Rich pair", PALETTE["rich_pair"]),
        ("classifier_cliff_score_rich_interaction_saprot", "Rich target", PALETTE["rich_target"]),
        ("classifier_cliff_score_rich_interaction_saprot_shuffled", "Shuffled target", PALETTE["shuffled"]),
    ]
    splits = ["cold_pair", "pair_scaffold", "ligand_component"]
    fig, ax = plt.subplots(figsize=(6.0, 3.3))
    x = np.arange(len(splits))
    width = 0.24
    for idx, (variant, label, color) in enumerate(variants):
        vals = []
        errs = []
        for split in splits:
            row = summary[(summary["split_mode"] == split) & (summary["variant"] == variant)]
            vals.append(float(row["cliff_auc_mean"].iloc[0]) if not row.empty else np.nan)
            errs.append(float(row["cliff_auc_std"].iloc[0]) if not row.empty and "cliff_auc_std" in row else 0.0)
        ax.bar(x + (idx - 1) * width, vals, width, yerr=errs, capsize=2, label=label, color=color)
    ax.set_xticks(x)
    ax.set_xticklabels(["Cold pair", "Pair scaffold", "Ligand component"])
    ax.set_ylabel("Cliff AUC")
    ax.set_ylim(0.45, 0.90)
    ax.set_title("Chemical OOD")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    save(fig, outdir, "fig_rich_pair_chemical_ood")


def plot_acnet(summary: pd.DataFrame, baseline: pd.DataFrame, neural: pd.DataFrame, outdir: Path) -> None:
    if summary.empty:
        return
    rows = []
    rich_map = {
        "rich_pair_only": ("Rich pair", PALETTE["rich_pair"]),
        "rich_seq_interaction": ("Rich target", PALETTE["rich_target"]),
        "rich_seq_interaction_shuffled": ("Shuffled target", PALETTE["shuffled"]),
    }
    for _, row in summary.iterrows():
        label, color = rich_map.get(row["variant"], (row["variant"], "#333333"))
        rows.append({"split": row["split_mode"], "method": label, "roc": row["roc_auc_mean"], "pr": row["pr_auc_mean"], "color": color})
    if not baseline.empty:
        for _, row in baseline[baseline["variant"] == "ecfp_absdiff_xgb"].iterrows():
            rows.append({"split": row["split_mode"], "method": "ECFP-XGB", "roc": row["roc_auc_mean"], "pr": row["pr_auc_mean"], "color": PALETTE["old"]})
    if not neural.empty:
        for variant, label in [("siamese_ecfp", "Siamese ECFP"), ("target_seq_film", "Seq-FiLM")]:
            for _, row in neural[neural["variant"] == variant].iterrows():
                rows.append({"split": row["split_mode"], "method": label, "roc": row["roc_auc_mean"], "pr": row["pr_auc_mean"], "color": PALETTE["neural"]})
    data = pd.DataFrame(rows)
    if data.empty:
        return
    for metric, ylabel, name in [("roc", "ROC-AUC", "fig_acnet_external_roc"), ("pr", "PR-AUC", "fig_acnet_external_pr")]:
        methods = list(dict.fromkeys(data["method"]))
        splits = ["random", "target"]
        fig, ax = plt.subplots(figsize=(7.0, 3.4))
        width = 0.8 / max(len(methods), 1)
        x = np.arange(len(splits))
        for idx, method in enumerate(methods):
            vals = []
            color = data[data["method"] == method]["color"].iloc[0]
            for split in splits:
                sub = data[(data["split"] == split) & (data["method"] == method)]
                vals.append(float(sub[metric].iloc[0]) if not sub.empty else np.nan)
            ax.bar(x - 0.4 + width / 2 + idx * width, vals, width, label=method, color=color)
        ax.set_xticks(x)
        ax.set_xticklabels(["Random", "Target-DG"])
        ax.set_ylabel(ylabel)
        ax.set_title(f"ACNet external {ylabel}")
        ax.grid(axis="y", alpha=0.25)
        ax.legend(frameon=False, fontsize=7, ncol=3)
        save(fig, outdir, name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="figures")
    parser.add_argument("--threshold", default="results/rich_pair_threshold_sensitivity_summary.csv")
    parser.add_argument("--ood", default="results/rich_pair_chemical_ood_summary.csv")
    parser.add_argument("--acnet", default="results/acnet_rich_pair_external_screened_seed0_summary.csv")
    parser.add_argument("--acnet-baseline", default="results/acnet_external_ecfp_xgb_summary.csv")
    parser.add_argument("--acnet-neural", default="results/acnet_sequence_neural_baselines_summary.csv")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    plot_threshold(load_csv(args.threshold), outdir)
    plot_ood(load_csv(args.ood), outdir)
    plot_acnet(load_csv(args.acnet), load_csv(args.acnet_baseline), load_csv(args.acnet_neural), outdir)
    print(f"Wrote figures to {outdir}")


if __name__ == "__main__":
    main()
