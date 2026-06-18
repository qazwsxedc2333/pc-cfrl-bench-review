from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import binomtest, wilcoxon
from sklearn.ensemble import ExtraTreesClassifier, HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

import run_chembl_raw_source_temporal_pairs as raw


HIGHER_IS_BETTER = {
    "roc_auc": True,
    "pr_auc": True,
    "balanced_accuracy": True,
    "f1": True,
    "mcc": True,
    "brier": False,
    "ece_10": False,
}


def parse_csv_arg(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()]


def fit_predict_variant(
    variant: str,
    diff_bits: sparse.csr_matrix,
    scalars_raw: np.ndarray,
    records: pd.DataFrame,
    train_mask: np.ndarray,
    test_mask: np.ndarray,
    seed: int,
    threads: int,
    n_estimators: int,
) -> tuple[np.ndarray, float, float]:
    y_train = records.loc[train_mask, "label"].to_numpy(dtype=np.int32)
    train_start = time.perf_counter()
    predict_seconds = 0.0

    if variant == "sgd_logloss_diff":
        model = SGDClassifier(
            loss="log_loss",
            alpha=1e-5,
            class_weight="balanced",
            max_iter=1500,
            tol=1e-4,
            random_state=seed,
            n_jobs=threads,
        )
        model.fit(diff_bits[train_mask], y_train)
        train_seconds = time.perf_counter() - train_start
        pred_start = time.perf_counter()
        prob = model.predict_proba(diff_bits[test_mask])[:, 1].astype(np.float32)
        predict_seconds = time.perf_counter() - pred_start
        return prob, train_seconds, predict_seconds

    scaler = StandardScaler().fit(scalars_raw[train_mask])
    scalars = scaler.transform(scalars_raw).astype(np.float32)

    if variant == "histgb_scalars":
        model = HistGradientBoostingClassifier(
            max_iter=n_estimators,
            learning_rate=0.04,
            max_leaf_nodes=31,
            l2_regularization=0.1,
            random_state=seed,
        )
        model.fit(scalars[train_mask], y_train)
        train_seconds = time.perf_counter() - train_start
        pred_start = time.perf_counter()
        prob = model.predict_proba(scalars[test_mask])[:, 1].astype(np.float32)
        predict_seconds = time.perf_counter() - pred_start
        return prob, train_seconds, predict_seconds

    if variant in {"rf_scalars", "extratrees_scalars"}:
        cls = RandomForestClassifier if variant == "rf_scalars" else ExtraTreesClassifier
        model = cls(
            n_estimators=n_estimators,
            max_depth=None,
            min_samples_leaf=3,
            class_weight="balanced_subsample" if variant == "rf_scalars" else "balanced",
            n_jobs=threads,
            random_state=seed,
        )
        model.fit(scalars[train_mask], y_train)
        train_seconds = time.perf_counter() - train_start
        pred_start = time.perf_counter()
        prob = model.predict_proba(scalars[test_mask])[:, 1].astype(np.float32)
        predict_seconds = time.perf_counter() - pred_start
        return prob, train_seconds, predict_seconds

    if variant == "knn_sar_diff":
        nn = NearestNeighbors(n_neighbors=min(15, int(train_mask.sum())), metric="cosine", n_jobs=threads)
        nn.fit(diff_bits[train_mask])
        train_seconds = time.perf_counter() - train_start
        pred_start = time.perf_counter()
        _, idx = nn.kneighbors(diff_bits[test_mask], return_distance=True)
        prob = y_train[idx].mean(axis=1).astype(np.float32)
        predict_seconds = time.perf_counter() - pred_start
        return prob, train_seconds, predict_seconds

    if variant == "mmp_rule_pair_scaffold":
        train = records.loc[train_mask, ["pair_scaffold_key", "label"]].copy()
        prior = float(train["label"].mean())
        table = train.groupby("pair_scaffold_key")["label"].mean().to_dict()
        prob = records.loc[test_mask, "pair_scaffold_key"].astype(str).map(table).fillna(prior).to_numpy(dtype=np.float32)
        train_seconds = time.perf_counter() - train_start
        return prob, train_seconds, predict_seconds

    raise ValueError(variant)


def paired_deltas(df: pd.DataFrame, baseline: str) -> pd.DataFrame:
    rows = []
    for split in sorted(df["split_mode"].astype(str).unique()):
        sub = df[df["split_mode"] == split]
        base = sub[sub["variant"] == baseline].set_index(["seed", "fold"])
        for variant in sorted(sub["variant"].astype(str).unique()):
            if variant == baseline:
                continue
            cur = sub[sub["variant"] == variant].set_index(["seed", "fold"])
            idx = cur.index.intersection(base.index)
            if len(idx) == 0:
                continue
            for metric, higher in HIGHER_IS_BETTER.items():
                if metric not in cur or metric not in base:
                    continue
                diff = cur.loc[idx, metric].astype(float) - base.loc[idx, metric].astype(float)
                diff = diff.replace([np.inf, -np.inf], np.nan).dropna()
                if diff.empty:
                    continue
                favorable = diff > 0 if higher else diff < 0
                wins = int(favorable.sum())
                try:
                    wilcoxon_p = float(wilcoxon(diff).pvalue) if (diff != 0).any() else 1.0
                except ValueError:
                    wilcoxon_p = float("nan")
                rows.append(
                    {
                        "split_mode": split,
                        "variant": variant,
                        "baseline": baseline,
                        "metric": metric,
                        "n_pairs": int(len(diff)),
                        "delta_mean": float(diff.mean()),
                        "delta_std": float(diff.std(ddof=1)) if len(diff) > 1 else 0.0,
                        "wins_favorable": wins,
                        "sign_p_favorable": float(binomtest(wins, len(diff), 0.5, alternative="greater").pvalue),
                        "wilcoxon_p": wilcoxon_p,
                    }
                )
    return pd.DataFrame(rows)


def write_report(results: pd.DataFrame, summary: pd.DataFrame, deltas: pd.DataFrame, out: Path) -> None:
    compact_cols = [
        "benchmark",
        "split_mode",
        "variant",
        "roc_auc_mean",
        "roc_auc_std",
        "pr_auc_mean",
        "pr_auc_std",
        "train_seconds_mean",
        "predict_seconds_mean",
        "target_overlap_rate_mean",
        "exact_ligand_overlap_rate_mean",
        "scaffold_overlap_rate_mean",
        "document_source_overlap_rate_mean",
    ]
    compact = summary[[c for c in compact_cols if c in summary.columns]]
    delta_focus = deltas[deltas["metric"].isin(["roc_auc", "pr_auc", "mcc", "f1"])].copy()
    lines = [
        "# TKDE Non-Neural Baseline Stress Test",
        "",
        "## Compact Summary",
        "",
        raw.markdown_table(compact),
        "",
        "## Paired Deltas",
        "",
        raw.markdown_table(delta_focus),
        "",
        "## Interpretation",
        "",
        "- These baselines answer whether the result is only an artifact of comparing against one weak ECFP-XGB baseline.",
        "- `knn_sar_diff` and `mmp_rule_pair_scaffold` are intentionally simple SAR/rule baselines; they are useful even when they underperform because they bound nearest-neighbor memorization.",
        "- Tree/scalar baselines test whether descriptor-level shifts alone explain the raw database results.",
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-csv", required=True)
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--split-modes", default="target_family,family_scaffold_source_purged")
    parser.add_argument("--variants", default="sgd_logloss_diff,histgb_scalars,rf_scalars,extratrees_scalars,knn_sar_diff,mmp_rule_pair_scaffold")
    parser.add_argument("--seeds", default="0,1,2,3,4")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--target-family-clusters", type=int, default=10)
    parser.add_argument("--sequence-hash-features", type=int, default=128)
    parser.add_argument("--n-estimators", type=int, default=80)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    records = pd.read_csv(args.pair_csv)
    diff_bits, scalars_raw = raw.build_feature_matrices(records)
    seeds = [int(x) for x in parse_csv_arg(args.seeds)]
    split_modes = parse_csv_arg(args.split_modes)
    variants = parse_csv_arg(args.variants)
    rows = []
    for seed in seeds:
        masks, target_to_family = raw.family_folds(
            records,
            seed,
            args.folds,
            args.sequence_hash_features,
            args.target_family_clusters,
        )
        family_arr = records["target"].astype(str).map(target_to_family).to_numpy(dtype=np.int32)
        for fold, test_mask in enumerate(masks[: args.folds]):
            base_train_mask = ~test_mask
            for split_mode in split_modes:
                train_mask, retention = raw.purge_train_mask(records, base_train_mask, test_mask, split_mode)
                y_train = records.loc[train_mask, "label"].to_numpy(dtype=np.int32)
                y_test = records.loc[test_mask, "label"].to_numpy(dtype=np.int32)
                if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
                    continue
                train_families = set(family_arr[train_mask])
                test_families = set(family_arr[test_mask])
                for variant in variants:
                    prob, train_seconds, predict_seconds = fit_predict_variant(
                        variant,
                        diff_bits,
                        scalars_raw,
                        records,
                        train_mask,
                        test_mask,
                        seed,
                        args.threads,
                        args.n_estimators,
                    )
                    row = {
                        "benchmark": args.benchmark,
                        "split_mode": split_mode,
                        "variant": variant,
                        "seed": seed,
                        "fold": fold,
                        "n_train": int(train_mask.sum()),
                        "n_test": int(test_mask.sum()),
                        "train_retention_after_purge": float(retention),
                        "train_positive_rate": float(y_train.mean()),
                        "test_positive_rate": float(y_test.mean()),
                        "n_train_targets": int(records.loc[train_mask, "target"].nunique()),
                        "n_test_targets": int(records.loc[test_mask, "target"].nunique()),
                        "n_train_families": int(len(train_families)),
                        "n_test_families": int(len(test_families)),
                        "family_overlap_rate": float(len(train_families & test_families) / max(len(test_families), 1)),
                        "train_seconds": float(train_seconds),
                        "predict_seconds": float(predict_seconds),
                    }
                    row.update(raw.audit_overlaps(records, train_mask, test_mask))
                    row.update(raw.row_metrics(y_test, prob))
                    rows.append(row)
                    print(
                        f"[{args.benchmark} {split_mode} seed={seed} fold={fold} {variant}] "
                        f"roc={row['roc_auc']:.4f} pr={row['pr_auc']:.4f} train_s={train_seconds:.2f}",
                        flush=True,
                    )
    results = pd.DataFrame(rows)
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    results.to_json(f"{prefix}.jsonl", orient="records", lines=True)
    metric_cols = [c for c in list(raw.METRICS) + ["train_seconds", "predict_seconds", "n_train", "n_test", "train_retention_after_purge", "target_overlap_rate", "family_overlap_rate", "exact_ligand_overlap_rate", "scaffold_overlap_rate", "document_source_overlap_rate"] if c in results.columns]
    summary = results.groupby(["benchmark", "split_mode", "variant"], as_index=False)[metric_cols].agg(["mean", "std", "count"])
    summary.columns = ["_".join([x for x in col if x]).strip("_") for col in summary.columns.to_flat_index()]
    summary = summary.reset_index(drop=True)
    deltas = paired_deltas(results, baseline=variants[0])
    summary.to_csv(f"{prefix}_summary.csv", index=False)
    deltas.to_csv(f"{prefix}_paired_deltas.csv", index=False)
    write_report(results, summary, deltas, Path(f"{prefix}_analysis.md"))
    print(json.dumps({"rows": int(len(results)), "output_prefix": str(prefix)}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
