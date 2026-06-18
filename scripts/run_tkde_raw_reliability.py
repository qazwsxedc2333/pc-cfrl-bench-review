from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon
from sklearn.model_selection import StratifiedShuffleSplit

import run_chembl_raw_source_temporal_pairs as raw


def parse_csv(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()]


def selected_metrics(y: np.ndarray, prob: np.ndarray, coverage: float) -> dict[str, float | int]:
    k = max(1, int(np.ceil(float(coverage) * len(y))))
    idx = np.argsort(-prob)[:k]
    yy = y[idx]
    precision = float(yy.mean())
    recall = float(yy.sum() / max(float(y.sum()), 1.0))
    return {
        "coverage": float(k / len(y)),
        "coverage_target": float(coverage),
        "precision": precision,
        "recall": recall,
        "f1": float(2 * precision * recall / max(precision + recall, 1e-12)),
        "positive_precision_lift": float(precision / max(float(y.mean()), 1e-12)),
        "avg_score": float(prob[idx].mean()),
    }


def class_probs(prob: np.ndarray) -> np.ndarray:
    return np.vstack([1.0 - prob, prob]).T.astype(np.float32)


def conformal_pvals(y_cal: np.ndarray, prob_cal: np.ndarray, prob_test: np.ndarray) -> np.ndarray:
    cal = class_probs(prob_cal)
    test = class_probs(prob_test)
    cal_scores = 1.0 - cal[np.arange(len(y_cal)), y_cal.astype(int)]
    out = np.zeros_like(test, dtype=np.float32)
    denom = float(len(cal_scores) + 1)
    for klass in (0, 1):
        scores = 1.0 - test[:, klass]
        out[:, klass] = (np.sum(cal_scores[:, None] >= scores[None, :], axis=0) + 1.0) / denom
    return out


def conformal_metrics(y: np.ndarray, pvals: np.ndarray, alpha: float) -> dict[str, float | int]:
    pred_sets = pvals > float(alpha)
    sizes = pred_sets.sum(axis=1)
    contains = pred_sets[np.arange(len(y)), y.astype(int)]
    singleton = sizes == 1
    positive_singleton = singleton & pred_sets[:, 1]
    pos_y = y[positive_singleton]
    prec = float(pos_y.mean()) if len(pos_y) else float("nan")
    rec = float(pos_y.sum() / max(float((y == 1).sum()), 1.0)) if len(pos_y) else 0.0
    return {
        "alpha": float(alpha),
        "coverage": float(contains.mean()),
        "coverage_gap": float(contains.mean() - (1.0 - alpha)),
        "avg_set_size": float(sizes.mean()),
        "singleton_rate": float(singleton.mean()),
        "positive_singleton_precision": prec,
        "positive_singleton_recall": rec,
        "positive_singleton_f1": float(2 * prec * rec / max(prec + rec, 1e-12)) if len(pos_y) else 0.0,
        "empty_rate": float((sizes == 0).mean()),
    }


def fit_calib_split(train_mask: np.ndarray, y_train: np.ndarray, seed: int, fraction: float) -> tuple[np.ndarray, np.ndarray]:
    idx = np.flatnonzero(train_mask)
    splitter = StratifiedShuffleSplit(n_splits=1, test_size=fraction, random_state=seed)
    fit_rel, cal_rel = next(splitter.split(idx, y_train))
    fit = np.zeros_like(train_mask, dtype=bool)
    cal = np.zeros_like(train_mask, dtype=bool)
    fit[idx[fit_rel]] = True
    cal[idx[cal_rel]] = True
    return fit, cal


def run_source(args: argparse.Namespace, source: str, pair_csv: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    records = pd.read_csv(pair_csv)
    diff_bits, scalars = raw.build_feature_matrices(records)
    seeds = [int(x) for x in parse_csv(args.seeds)]
    variants = parse_csv(args.variants)
    coverages = [float(x) for x in parse_csv(args.coverages)]
    alphas = [float(x) for x in parse_csv(args.alphas)]
    selective_rows = []
    conformal_rows = []
    for seed in seeds:
        masks, target_to_family = raw.family_folds(records, seed, args.folds, args.sequence_hash_features, args.target_family_clusters)
        for fold in range(args.folds):
            split_modes = ["family_scaffold_source_purged"]
            if source == "chembl_raw":
                split_modes.append("temporal_scaffold_source_purged")
            for split_mode in split_modes:
                if split_mode.startswith("temporal"):
                    base_train, test_mask, cutoff = raw.temporal_masks(records, fold, args.folds)
                else:
                    test_mask = masks[fold]
                    base_train = ~test_mask
                    cutoff = -1
                train_mask, retention = raw.purge_train_mask(records, base_train, test_mask, split_mode)
                if train_mask.sum() < args.min_train_rows or test_mask.sum() < args.min_test_rows:
                    continue
                y_train = records.loc[train_mask, "label"].to_numpy(dtype=np.int32)
                y_test = records.loc[test_mask, "label"].to_numpy(dtype=np.int32)
                if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
                    continue
                audit = raw.audit_overlaps(records, train_mask, test_mask)
                for variant in variants:
                    fit_mask, cal_mask = fit_calib_split(train_mask, y_train, 8000 + seed * 101 + fold, args.calib_fraction)
                    y_fit = records.loc[fit_mask, "label"].to_numpy(dtype=np.int32)
                    y_cal = records.loc[cal_mask, "label"].to_numpy(dtype=np.int32)
                    x = raw.build_x(diff_bits, scalars, fit_mask, variant)
                    model = raw.make_model(seed, y_fit, args.model_threads, args.n_estimators, args.max_depth)
                    model.fit(x[fit_mask], y_fit)
                    prob_cal = model.predict_proba(x[cal_mask])[:, 1].astype(np.float32)
                    prob_test = model.predict_proba(x[test_mask])[:, 1].astype(np.float32)
                    base = {
                        "source": source,
                        "split_mode": split_mode,
                        "variant": variant,
                        "seed": seed,
                        "fold": fold,
                        "temporal_cutoff_year": int(cutoff),
                        "n_train": int(train_mask.sum()),
                        "n_calib": int(cal_mask.sum()),
                        "n_test": int(test_mask.sum()),
                        "train_retention_after_purge": float(retention),
                        "test_positive_rate": float(y_test.mean()),
                        **audit,
                    }
                    for coverage in coverages:
                        row = dict(base)
                        row.update(selected_metrics(y_test, prob_test, coverage))
                        selective_rows.append(row)
                    pvals = conformal_pvals(y_cal, prob_cal, prob_test)
                    for alpha in alphas:
                        row = dict(base)
                        row.update(conformal_metrics(y_test, pvals, alpha))
                        conformal_rows.append(row)
                    print(f"[{source} {split_mode} seed={seed} fold={fold} {variant}] reliability_done", flush=True)
    return pd.DataFrame(selective_rows), pd.DataFrame(conformal_rows)


def summarize(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    metric_cols = [c for c in df.columns if c not in set(group_cols + ["seed", "fold", "temporal_cutoff_year"])]
    metric_cols = [c for c in metric_cols if pd.api.types.is_numeric_dtype(df[c])]
    out = df.groupby(group_cols, as_index=False)[metric_cols].agg(["mean", "std", "count"])
    out.columns = ["_".join([x for x in col if x]).strip("_") for col in out.columns.to_flat_index()]
    return out.reset_index(drop=True)


def paired(df: pd.DataFrame, setting_cols: list[str], metrics: list[str]) -> pd.DataFrame:
    rows = []
    for keys, sub in df.groupby(setting_cols, sort=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        base = sub[sub["variant"] == "ecfp_absdiff_xgb"].set_index(["seed", "fold"])
        for variant in sorted(sub["variant"].astype(str).unique()):
            if variant == "ecfp_absdiff_xgb":
                continue
            cur = sub[sub["variant"] == variant].set_index(["seed", "fold"])
            idx = cur.index.intersection(base.index)
            if len(idx) == 0:
                continue
            for metric in metrics:
                if metric not in cur or metric not in base:
                    continue
                diff = cur.loc[idx, metric].astype(float) - base.loc[idx, metric].astype(float)
                diff = diff.replace([np.inf, -np.inf], np.nan).dropna()
                if diff.empty:
                    continue
                row = dict(zip(setting_cols, keys))
                row.update({"variant": variant, "baseline": "ecfp_absdiff_xgb", "metric": metric, "n_pairs": int(len(diff)), "delta_mean": float(diff.mean()), "wins": int((diff > 0).sum())})
                try:
                    row["wilcoxon_p"] = float(wilcoxon(diff).pvalue) if (diff != 0).any() else 1.0
                except ValueError:
                    row["wilcoxon_p"] = float("nan")
                row["sign_p"] = float(binomtest(int((diff > 0).sum()), len(diff), 0.5, alternative="greater").pvalue)
                rows.append(row)
    return pd.DataFrame(rows)


def write_report(sel_sum: pd.DataFrame, conf_sum: pd.DataFrame, sel_delta: pd.DataFrame, conf_delta: pd.DataFrame, out: Path) -> None:
    lines = [
        "# Raw External Selective And Conformal Reliability",
        "",
        "## Selective Summary",
        "",
        raw.markdown_table(sel_sum),
        "",
        "## Selective Paired Deltas",
        "",
        raw.markdown_table(sel_delta),
        "",
        "## Conformal Summary",
        "",
        raw.markdown_table(conf_sum),
        "",
        "## Conformal Paired Deltas",
        "",
        raw.markdown_table(conf_delta),
        "",
        "## Reading",
        "",
        "- This extends reliability evidence from ACNet to raw ChEMBL/BindingDB source-purged settings.",
        "- Use selective precision/lift for triage utility and conformal coverage/set size for uncertainty-aware deployment wording.",
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sources", default="chembl_raw,bindingdb_raw")
    parser.add_argument("--chembl-pairs", default="data/chembl_raw_moleculeace_target_pairs.csv")
    parser.add_argument("--bindingdb-pairs", default="data/bindingdb_raw_moleculeace_target_pairs.csv")
    parser.add_argument("--variants", default="ecfp_absdiff_xgb,rich_diff_scalars")
    parser.add_argument("--seeds", default="0,1,2")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--coverages", default="0.4,0.2,0.1")
    parser.add_argument("--alphas", default="0.1,0.2")
    parser.add_argument("--calib-fraction", type=float, default=0.2)
    parser.add_argument("--target-family-clusters", type=int, default=10)
    parser.add_argument("--sequence-hash-features", type=int, default=128)
    parser.add_argument("--min-train-rows", type=int, default=500)
    parser.add_argument("--min-test-rows", type=int, default=200)
    parser.add_argument("--n-estimators", type=int, default=40)
    parser.add_argument("--max-depth", type=int, default=4)
    parser.add_argument("--model-threads", type=int, default=2)
    parser.add_argument("--output-prefix", default="results/tkde_raw_reliability_seed0_2")
    args = parser.parse_args()
    paths = {"chembl_raw": args.chembl_pairs, "bindingdb_raw": args.bindingdb_pairs}
    sel_frames, conf_frames = [], []
    for source in parse_csv(args.sources):
        sel, conf = run_source(args, source, paths[source])
        sel_frames.append(sel)
        conf_frames.append(conf)
    selective = pd.concat(sel_frames, ignore_index=True)
    conformal = pd.concat(conf_frames, ignore_index=True)
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    selective.to_csv(f"{prefix}_selective_rows.csv", index=False)
    conformal.to_csv(f"{prefix}_conformal_rows.csv", index=False)
    sel_sum = summarize(selective, ["source", "split_mode", "variant", "coverage_target"])
    conf_sum = summarize(conformal, ["source", "split_mode", "variant", "alpha"])
    sel_delta = paired(selective, ["source", "split_mode", "coverage_target"], ["precision", "positive_precision_lift", "recall", "f1"])
    conf_delta = paired(conformal, ["source", "split_mode", "alpha"], ["coverage", "avg_set_size", "singleton_rate", "positive_singleton_precision", "positive_singleton_f1"])
    sel_sum.to_csv(f"{prefix}_selective_summary.csv", index=False)
    conf_sum.to_csv(f"{prefix}_conformal_summary.csv", index=False)
    sel_delta.to_csv(f"{prefix}_selective_paired_deltas.csv", index=False)
    conf_delta.to_csv(f"{prefix}_conformal_paired_deltas.csv", index=False)
    write_report(sel_sum, conf_sum, sel_delta, conf_delta, Path(f"{prefix}_report.md"))
    print(json.dumps({"selective_rows": int(len(selective)), "conformal_rows": int(len(conformal)), "output_prefix": str(prefix)}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
