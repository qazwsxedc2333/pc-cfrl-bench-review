from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


SOURCES = [
    ("ACNet", "results/acnet_target_family_best_joint_seed0_4_allfold_n10_summary.csv"),
    ("MoleculeACE", "results/moleculeace_target_family_pairs_seed0_4_n120_summary.csv"),
    ("ChEMBL raw", "results/chembl_raw_source_temporal_pairs_seed0_4_n40_summary.csv"),
    ("BindingDB raw", "results/bindingdb_raw_source_pairs_seed0_4_n80_summary.csv"),
    ("LoHi raw", "results/tkde_lohi_leakage_splits_seed0_2_summary.csv"),
    ("Noise raw", "results/tkde_noise_sensitivity_seed0_2_summary.csv"),
    ("Remaining P0-P1", "results/tkde_remaining_p0_p1_seed0_2_summary.csv"),
]


def read_summary(label: str, path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p)
    if "benchmark" not in df.columns:
        df["benchmark"] = label
    else:
        df["benchmark"] = df["benchmark"].fillna(label)
    if "source" in df.columns:
        df["source"] = df["source"].fillna(label)
    else:
        df["source"] = label
    if "condition" not in df.columns:
        df["condition"] = "default"
    keep = [
        "benchmark",
        "source",
        "condition",
        "split_mode",
        "variant",
        "roc_auc_mean",
        "roc_auc_std",
        "pr_auc_mean",
        "pr_auc_std",
        "mcc_mean",
        "mcc_std",
        "roc_auc_count",
        "n_train_mean",
        "n_test_mean",
        "train_retention_after_purge_mean",
        "target_overlap_rate_mean",
        "exact_ligand_overlap_rate_mean",
        "scaffold_overlap_rate_mean",
        "document_source_overlap_rate_mean",
        "ef_1pct_mean",
        "ef_5pct_mean",
        "bedroc20_mean",
        "precision_at_50_mean",
    ]
    return df[[c for c in keep if c in df.columns]].copy()


def write_schema(out: Path, snapshot_path: str) -> None:
    lines = [
        "# TKDE Leaderboard Package Schema",
        "",
        f"Snapshot file: `{snapshot_path}`",
        "",
        "## Required Columns",
        "",
        "| column | meaning |",
        "|---|---|",
        "| benchmark/source | benchmark block and raw source |",
        "| condition | default or robustness condition |",
        "| split_mode | split protocol name |",
        "| variant | model or feature variant |",
        "| roc_auc_mean/std, pr_auc_mean/std, mcc_mean/std | primary metrics over seed/fold runs |",
        "| ef_1pct_mean, ef_5pct_mean, bedroc20_mean, precision_at_50_mean | drug-discovery ranking utility metrics when available |",
        "| *_overlap_rate_mean | leakage audit fields |",
        "",
        "## Reporting Rule",
        "",
        "Every leaderboard entry must cite its split protocol and overlap audit. Hard-OOD claims should prioritize source/scaffold-purged and frozen TKDE protocol rows.",
    ]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snapshot", default="results/tkde_leaderboard_snapshot.csv")
    parser.add_argument("--schema", default="docs/tkde_leaderboard_schema.md")
    args = parser.parse_args()
    frames = [read_summary(label, path) for label, path in SOURCES]
    frames = [df for df in frames if not df.empty]
    if not frames:
        raise RuntimeError("No summary files found")
    snapshot = pd.concat(frames, ignore_index=True)
    Path(args.snapshot).parent.mkdir(parents=True, exist_ok=True)
    snapshot.to_csv(args.snapshot, index=False)
    write_schema(Path(args.schema), args.snapshot)
    print({"rows": int(len(snapshot)), "snapshot": args.snapshot, "schema": args.schema})


if __name__ == "__main__":
    main()
