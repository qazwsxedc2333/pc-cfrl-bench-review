from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

import run_chembl_raw_source_temporal_pairs as raw


def parse_csv(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()]


def pair_scaffold_folds(records: pd.DataFrame, seed: int, n_folds: int) -> list[np.ndarray]:
    stats = records.groupby("pair_scaffold_key")["label"].agg(["size", "sum"]).reset_index()
    stats = stats.sample(frac=1.0, random_state=seed).sort_values(["sum", "size"], ascending=False)
    folds: list[list[str]] = [[] for _ in range(n_folds)]
    fold_rows = np.zeros(n_folds, dtype=np.int64)
    fold_pos = np.zeros(n_folds, dtype=np.int64)
    for row in stats.itertuples(index=False):
        idx = int(np.argmin(fold_rows + 3 * fold_pos))
        folds[idx].append(str(row.pair_scaffold_key))
        fold_rows[idx] += int(row.size)
        fold_pos[idx] += int(row.sum)
    key = records["pair_scaffold_key"].astype(str)
    return [key.isin(fold).to_numpy(dtype=bool) for fold in folds if fold]


def ligand_block_folds(records: pd.DataFrame, seed: int, n_folds: int) -> list[np.ndarray]:
    ligand_to_scaffold: dict[str, str] = {}
    for row in records.itertuples(index=False):
        ligand_to_scaffold[str(row.smiles1)] = str(row.scaffold1)
        ligand_to_scaffold[str(row.smiles2)] = str(row.scaffold2)
    scaffold_stats = []
    for scaffold in sorted(set(ligand_to_scaffold.values())):
        mask = records["scaffold1"].astype(str).eq(scaffold) | records["scaffold2"].astype(str).eq(scaffold)
        scaffold_stats.append({"scaffold": scaffold, "size": int(mask.sum()), "sum": int(records.loc[mask, "label"].sum())})
    stats = pd.DataFrame(scaffold_stats).sample(frac=1.0, random_state=seed).sort_values(["sum", "size"], ascending=False)
    folds: list[list[str]] = [[] for _ in range(n_folds)]
    fold_rows = np.zeros(n_folds, dtype=np.int64)
    fold_pos = np.zeros(n_folds, dtype=np.int64)
    for row in stats.itertuples(index=False):
        idx = int(np.argmin(fold_rows + 3 * fold_pos))
        folds[idx].append(str(row.scaffold))
        fold_rows[idx] += int(row.size)
        fold_pos[idx] += int(row.sum)
    s1 = records["scaffold1"].astype(str)
    s2 = records["scaffold2"].astype(str)
    return [(s1.isin(fold) | s2.isin(fold)).to_numpy(dtype=bool) for fold in folds if fold]


def purge(records: pd.DataFrame, train_mask: np.ndarray, test_mask: np.ndarray, source_purge: bool) -> tuple[np.ndarray, float]:
    n_base = int(train_mask.sum())
    test_scaffolds = set(records.loc[test_mask, "scaffold1"].astype(str)) | set(records.loc[test_mask, "scaffold2"].astype(str))
    keep = ~(
        records["scaffold1"].astype(str).isin(test_scaffolds).to_numpy()
        | records["scaffold2"].astype(str).isin(test_scaffolds).to_numpy()
    )
    train_mask = train_mask & keep
    if source_purge:
        train_mask = raw.source_overlap_mask(records, train_mask, test_mask)
    return train_mask, float(train_mask.sum() / max(n_base, 1))


def build_x(diff_bits: sparse.csr_matrix, scalars: np.ndarray, train_mask: np.ndarray, variant: str):
    return raw.build_x(diff_bits, scalars, train_mask, variant)


def run_source(args: argparse.Namespace, source: str, pair_csv: str) -> pd.DataFrame:
    records = pd.read_csv(pair_csv)
    diff_bits, scalars = raw.build_feature_matrices(records)
    variants = parse_csv(args.variants)
    seeds = [int(x) for x in parse_csv(args.seeds)]
    rows = []
    for seed in seeds:
        fold_sets = {
            "lohi_pair_scaffold_source_purged": pair_scaffold_folds(records, seed, args.folds),
            "lohi_ligand_scaffold_source_purged": ligand_block_folds(records, seed, args.folds),
        }
        for split_mode, masks in fold_sets.items():
            for fold, test_mask in enumerate(masks[: args.folds]):
                base_train = ~test_mask
                train_mask, retention = purge(records, base_train, test_mask, "source" in split_mode)
                if train_mask.sum() < args.min_train_rows or test_mask.sum() < args.min_test_rows:
                    continue
                y_train = records.loc[train_mask, "label"].to_numpy(dtype=np.int32)
                y_test = records.loc[test_mask, "label"].to_numpy(dtype=np.int32)
                if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
                    continue
                for variant in variants:
                    x = build_x(diff_bits, scalars, train_mask, variant)
                    model = raw.make_model(seed, y_train, args.model_threads, args.n_estimators, args.max_depth)
                    model.fit(x[train_mask], y_train)
                    prob = model.predict_proba(x[test_mask])[:, 1].astype(np.float32)
                    row = {
                        "source": source,
                        "split_mode": split_mode,
                        "variant": variant,
                        "seed": seed,
                        "fold": fold,
                        "n_train": int(train_mask.sum()),
                        "n_test": int(test_mask.sum()),
                        "train_retention_after_purge": float(retention),
                        "test_positive_rate": float(y_test.mean()),
                    }
                    row.update(raw.audit_overlaps(records, train_mask, test_mask))
                    row.update(raw.row_metrics(y_test, prob))
                    rows.append(row)
                    print(f"[{source} {split_mode} seed={seed} fold={fold} {variant}] roc={row['roc_auc']:.4f} pr={row['pr_auc']:.4f}", flush=True)
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in raw.METRICS + ["n_train", "n_test", "train_retention_after_purge", "target_overlap_rate", "exact_ligand_overlap_rate", "scaffold_overlap_rate", "document_source_overlap_rate"] if c in df]
    out = df.groupby(["source", "split_mode", "variant"], as_index=False)[cols].agg(["mean", "std", "count"])
    out.columns = ["_".join([x for x in col if x]).strip("_") for col in out.columns.to_flat_index()]
    return out.reset_index(drop=True)


def paired_deltas(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (source, split), sub in df.groupby(["source", "split_mode"], sort=False):
        base = sub[sub["variant"] == "ecfp_absdiff_xgb"].set_index(["seed", "fold"])
        for variant in sorted(sub["variant"].astype(str).unique()):
            if variant == "ecfp_absdiff_xgb":
                continue
            cur = sub[sub["variant"] == variant].set_index(["seed", "fold"])
            idx = cur.index.intersection(base.index)
            if len(idx) == 0:
                continue
            for metric in raw.METRICS:
                if metric not in cur or metric not in base:
                    continue
                diff = cur.loc[idx, metric].astype(float) - base.loc[idx, metric].astype(float)
                diff = diff.replace([np.inf, -np.inf], np.nan).dropna()
                if diff.empty:
                    continue
                rows.append({
                    "source": source,
                    "split_mode": split,
                    "variant": variant,
                    "baseline": "ecfp_absdiff_xgb",
                    "metric": metric,
                    "n_pairs": int(len(diff)),
                    "delta_mean": float(diff.mean()),
                    "wins": int((diff > 0).sum()) if raw.HIGHER_IS_BETTER[metric] else int((diff < 0).sum()),
                })
    return pd.DataFrame(rows)


def write_report(summary: pd.DataFrame, deltas: pd.DataFrame, out: Path) -> None:
    compact_cols = ["source", "split_mode", "variant", "roc_auc_mean", "pr_auc_mean", "mcc_mean", "n_test_mean", "scaffold_overlap_rate_mean", "document_source_overlap_rate_mean"]
    focus = deltas[deltas["metric"].isin(["roc_auc", "pr_auc", "mcc"])].copy()
    lines = [
        "# LoHi / DataSAIL-Style Leakage-Aware Raw Splits",
        "",
        "This experiment adds external-style chemical block splits. Test folds are built from held-out pair-scaffold or ligand-scaffold blocks, then train rows sharing test scaffolds and publication/assay sources are removed.",
        "",
        "## Summary",
        "",
        raw.markdown_table(summary[[c for c in compact_cols if c in summary.columns]]),
        "",
        "## Paired Deltas vs ECFP",
        "",
        raw.markdown_table(focus),
        "",
        "## Reading",
        "",
        "- This is a leakage-aware sanity check inspired by LoHi/DataSAIL-style chemical separation, not a claim that the exact DataSAIL optimizer is used.",
        "- Use it to show that the raw conclusions are not tied only to our target-family split definition.",
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sources", default="chembl_raw,bindingdb_raw")
    parser.add_argument("--chembl-pairs", default="data/chembl_raw_moleculeace_target_pairs.csv")
    parser.add_argument("--bindingdb-pairs", default="data/bindingdb_raw_moleculeace_target_pairs.csv")
    parser.add_argument("--variants", default="ecfp_absdiff_xgb,rich_diff_scalars,rich_scalars_only")
    parser.add_argument("--seeds", default="0,1,2")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--min-train-rows", type=int, default=500)
    parser.add_argument("--min-test-rows", type=int, default=200)
    parser.add_argument("--n-estimators", type=int, default=40)
    parser.add_argument("--max-depth", type=int, default=4)
    parser.add_argument("--model-threads", type=int, default=2)
    parser.add_argument("--output-prefix", default="results/tkde_lohi_leakage_splits_seed0_2")
    args = parser.parse_args()
    paths = {"chembl_raw": args.chembl_pairs, "bindingdb_raw": args.bindingdb_pairs}
    frames = [run_source(args, src, paths[src]) for src in parse_csv(args.sources)]
    results = pd.concat(frames, ignore_index=True)
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    results.to_json(f"{prefix}.jsonl", orient="records", lines=True)
    summary = summarize(results)
    deltas = paired_deltas(results)
    summary.to_csv(f"{prefix}_summary.csv", index=False)
    deltas.to_csv(f"{prefix}_paired_deltas.csv", index=False)
    write_report(summary, deltas, Path(f"{prefix}_report.md"))
    print(json.dumps({"rows": int(len(results)), "output_prefix": str(prefix)}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
