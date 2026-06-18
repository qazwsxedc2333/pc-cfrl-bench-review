from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

import run_chembl_raw_source_temporal_pairs as raw
import run_literature_guided_multirf_fusion as multirf
import run_tkde_latent_anchor_portfolio as latent_anchor
import run_tkde_source_anchor_portfolio as source_anchor
import run_tnnls_multi_anchor_local_evidence as anchor


def parse_csv(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()]


def parse_float_csv(value: str) -> list[float]:
    return [float(x) for x in parse_csv(value)]


def parse_int_csv(value: str) -> list[int]:
    out: list[int] = []
    for x in parse_csv(value):
        if x.lower() in {"none", "null", "na", "nan", "-"}:
            continue
        out.append(int(x))
    return out


def choose_random_n(idx: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    n = max(1, min(int(n), max(len(idx) - 2, 1)))
    return rng.choice(idx, size=n, replace=False)


def masks_from_chosen(n_records: int, test_mask: np.ndarray, chosen: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    support_mask = np.zeros(n_records, dtype=bool)
    support_mask[np.asarray(chosen, dtype=np.int64)] = True
    query_mask = test_mask & ~support_mask
    return support_mask, query_mask


def make_stratified_frac_support(
    records: pd.DataFrame,
    test_mask: np.ndarray,
    seed: int,
    fold: int,
    frac: float,
) -> tuple[np.ndarray, np.ndarray, dict]:
    rng = np.random.default_rng(seed * 1009 + fold * 9173 + int(round(frac * 1000)) + 11)
    support_mask = np.zeros(len(records), dtype=bool)
    test_idx = np.flatnonzero(test_mask)
    y = records.loc[test_mask, "label"].to_numpy(dtype=np.int32)
    for lab in [0, 1]:
        local = test_idx[np.flatnonzero(y == lab)]
        if len(local) == 0:
            continue
        n = max(1, int(round(len(local) * frac)))
        n = min(n, max(len(local) - 1, 1))
        support_mask[rng.choice(local, size=n, replace=False)] = True
    query_mask = test_mask & ~support_mask
    return support_mask, query_mask, {"fewshot_support_protocol": "stratified_frac", "fewshot_requested_frac": float(frac)}


def make_natural_frac_support(
    records: pd.DataFrame,
    test_mask: np.ndarray,
    seed: int,
    fold: int,
    frac: float,
) -> tuple[np.ndarray, np.ndarray, dict]:
    rng = np.random.default_rng(seed * 1009 + fold * 9173 + int(round(frac * 1000)) + 23)
    test_idx = np.flatnonzero(test_mask)
    n = max(2, int(round(len(test_idx) * frac)))
    chosen = choose_random_n(test_idx, n, rng)
    support_mask, query_mask = masks_from_chosen(len(records), test_mask, chosen)
    return support_mask, query_mask, {"fewshot_support_protocol": "natural_frac", "fewshot_requested_frac": float(frac)}


def make_fixedk_balanced_support(
    records: pd.DataFrame,
    test_mask: np.ndarray,
    seed: int,
    fold: int,
    k_per_class: int,
) -> tuple[np.ndarray, np.ndarray, dict]:
    rng = np.random.default_rng(seed * 1009 + fold * 9173 + int(k_per_class) * 37 + 31)
    support_mask = np.zeros(len(records), dtype=bool)
    test_idx = np.flatnonzero(test_mask)
    y = records.loc[test_mask, "label"].to_numpy(dtype=np.int32)
    for lab in [0, 1]:
        local = test_idx[np.flatnonzero(y == lab)]
        if len(local) == 0:
            continue
        n = min(int(k_per_class), max(len(local) - 1, 1))
        support_mask[rng.choice(local, size=n, replace=False)] = True
    query_mask = test_mask & ~support_mask
    return support_mask, query_mask, {"fewshot_support_protocol": "fixedk_balanced", "fewshot_support_k_per_class": int(k_per_class)}


def active_uncertainty_diverse_support(
    records: pd.DataFrame,
    scalars_raw: np.ndarray,
    test_mask: np.ndarray,
    zero_score: np.ndarray,
    seed: int,
    fold: int,
    frac: float,
) -> tuple[np.ndarray, np.ndarray, dict]:
    test_idx = np.flatnonzero(test_mask)
    n = max(2, int(round(len(test_idx) * frac)))
    n = min(n, max(len(test_idx) - 2, 1))
    rng = np.random.default_rng(seed * 1009 + fold * 9173 + int(round(frac * 1000)) + 43)
    uncertainty = np.abs(np.asarray(zero_score, dtype=np.float32) - 0.5)
    order = np.argsort(uncertainty, kind="mergesort")
    n_uncertain = min(n, max(1, int(round(n * 0.7))))
    chosen_local = list(order[:n_uncertain])
    if len(chosen_local) < n:
        remaining_local = np.setdiff1d(np.arange(len(test_idx)), np.asarray(chosen_local, dtype=np.int64), assume_unique=False)
        if len(remaining_local) > 0:
            x = scalars_raw[test_idx[remaining_local]]
            x = StandardScaler().fit_transform(x).astype(np.float32)
            center = x.mean(axis=0, keepdims=True)
            dist = np.linalg.norm(x - center, axis=1)
            # Add a tiny deterministic jitter to avoid tie-order artifacts.
            jitter = rng.normal(0.0, 1.0e-8, size=len(dist))
            diverse_order = np.argsort(-(dist + jitter), kind="mergesort")
            chosen_local.extend(remaining_local[diverse_order[: n - len(chosen_local)]].tolist())
    chosen = test_idx[np.asarray(chosen_local[:n], dtype=np.int64)]
    support_mask, query_mask = masks_from_chosen(len(records), test_mask, chosen)
    return support_mask, query_mask, {
        "fewshot_support_protocol": "active_uncertainty_diverse_frac",
        "fewshot_requested_frac": float(frac),
        "active_uncertainty_share": 0.7,
    }


def active_prediction_balanced_support(
    records: pd.DataFrame,
    test_mask: np.ndarray,
    zero_score: np.ndarray,
    seed: int,
    fold: int,
    frac: float,
) -> tuple[np.ndarray, np.ndarray, dict]:
    test_idx = np.flatnonzero(test_mask)
    n = max(2, int(round(len(test_idx) * frac)))
    n = min(n, max(len(test_idx) - 2, 1))
    low_n = n // 2
    high_n = n - low_n
    rng = np.random.default_rng(seed * 1009 + fold * 9173 + int(round(frac * 1000)) + 59)
    score = np.asarray(zero_score, dtype=np.float32)
    jitter = rng.normal(0.0, 1.0e-8, size=len(score))
    low_order = np.argsort(score + jitter, kind="mergesort")
    high_order = np.argsort(-(score + jitter), kind="mergesort")
    chosen_local: list[int] = []
    chosen_set: set[int] = set()
    for local in low_order:
        if len(chosen_local) >= low_n:
            break
        chosen_local.append(int(local))
        chosen_set.add(int(local))
    for local in high_order:
        if len(chosen_local) >= n:
            break
        if int(local) not in chosen_set:
            chosen_local.append(int(local))
            chosen_set.add(int(local))
    if len(chosen_local) < n:
        remaining = np.asarray([i for i in range(len(test_idx)) if i not in chosen_set], dtype=np.int64)
        fill = rng.choice(remaining, size=n - len(chosen_local), replace=False)
        chosen_local.extend([int(x) for x in fill])
    chosen = test_idx[np.asarray(chosen_local[:n], dtype=np.int64)]
    support_mask, query_mask = masks_from_chosen(len(records), test_mask, chosen)
    return support_mask, query_mask, {
        "fewshot_support_protocol": "active_prediction_balanced_frac",
        "fewshot_requested_frac": float(frac),
        "active_low_score_share": float(low_n / max(n, 1)),
    }


def support_knn_prob(
    scalars_raw: np.ndarray,
    y_all: np.ndarray,
    support_mask: np.ndarray,
    query_mask: np.ndarray,
    k: int,
) -> np.ndarray | None:
    support_idx = np.flatnonzero(support_mask)
    query_idx = np.flatnonzero(query_mask)
    if len(support_idx) < 2 or len(query_idx) == 0:
        return None
    scaler = StandardScaler().fit(scalars_raw[support_idx])
    x_support = scaler.transform(scalars_raw[support_idx]).astype(np.float32)
    x_query = scaler.transform(scalars_raw[query_idx]).astype(np.float32)
    n_neighbors = max(1, min(int(k), len(support_idx)))
    nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
    nn.fit(x_support)
    dist, ind = nn.kneighbors(x_query, return_distance=True)
    weights = 1.0 / (dist + 1.0e-6)
    labels = y_all[support_idx][ind]
    prob = (weights * labels).sum(axis=1) / np.maximum(weights.sum(axis=1), 1.0e-8)
    return source_anchor.clip_prob(prob.astype(np.float32))


def support_extratrees_prob(
    scalars_raw: np.ndarray,
    y_all: np.ndarray,
    support_mask: np.ndarray,
    query_mask: np.ndarray,
    seed: int,
    threads: int,
    n_estimators: int,
) -> np.ndarray | None:
    support_idx = np.flatnonzero(support_mask)
    query_idx = np.flatnonzero(query_mask)
    if len(support_idx) < 8 or len(np.unique(y_all[support_idx])) < 2:
        return None
    scaler = StandardScaler().fit(scalars_raw[support_idx])
    x_support = scaler.transform(scalars_raw[support_idx]).astype(np.float32)
    x_query = scaler.transform(scalars_raw[query_idx]).astype(np.float32)
    model = ExtraTreesClassifier(
        n_estimators=max(60, min(int(n_estimators), 160)),
        min_samples_leaf=1,
        class_weight="balanced",
        n_jobs=threads,
        random_state=seed * 1291 + int(support_idx[0]),
    )
    model.fit(x_support, y_all[support_idx])
    return source_anchor.clip_prob(model.predict_proba(x_query)[:, 1].astype(np.float32))


def support_logreg_prob(
    scalars_raw: np.ndarray,
    y_all: np.ndarray,
    support_mask: np.ndarray,
    query_mask: np.ndarray,
    seed: int,
) -> np.ndarray | None:
    support_idx = np.flatnonzero(support_mask)
    query_idx = np.flatnonzero(query_mask)
    if len(support_idx) < 8 or len(np.unique(y_all[support_idx])) < 2:
        return None
    scaler = StandardScaler().fit(scalars_raw[support_idx])
    x_support = scaler.transform(scalars_raw[support_idx]).astype(np.float32)
    x_query = scaler.transform(scalars_raw[query_idx]).astype(np.float32)
    model = LogisticRegression(
        C=0.5,
        class_weight="balanced",
        max_iter=500,
        solver="liblinear",
        random_state=seed * 1777 + int(support_idx[0]),
    )
    model.fit(x_support, y_all[support_idx])
    return source_anchor.clip_prob(model.predict_proba(x_query)[:, 1].astype(np.float32))


def add_row(rows, records, family_arr, train_mask, query_mask, source, split_mode, variant, seed, fold, score, prob, seconds, extra=None):
    source_anchor.add_row(rows, records, family_arr, train_mask, query_mask, source, split_mode, variant, seed, fold, score, prob, seconds, extra or {})


def support_extra(records: pd.DataFrame, y_all: np.ndarray, support_mask: np.ndarray, query_mask: np.ndarray, protocol_extra: dict) -> dict:
    support_y = y_all[support_mask]
    query_y = y_all[query_mask]
    out = {
        **protocol_extra,
        "fewshot_support_rows": int(support_mask.sum()),
        "fewshot_query_rows": int(query_mask.sum()),
        "fewshot_support_frac_actual": float(support_mask.sum() / max((support_mask | query_mask).sum(), 1)),
        "fewshot_support_positive_rate": float(support_y.mean()) if len(support_y) else 0.0,
        "fewshot_query_positive_rate": float(query_y.mean()) if len(query_y) else 0.0,
        "fewshot_support_label_entropy": float(raw.entropy_binary(float(support_y.mean()))) if hasattr(raw, "entropy_binary") and len(support_y) else 0.0,
    }
    if "fewshot_requested_frac" not in out:
        out["fewshot_requested_frac"] = float(out["fewshot_support_frac_actual"])
    return out


def make_support_sets(
    args: argparse.Namespace,
    records: pd.DataFrame,
    scalars_raw: np.ndarray,
    heldout_mask: np.ndarray,
    zero_score: np.ndarray,
    seed: int,
    fold: int,
) -> list[tuple[str, np.ndarray, np.ndarray, dict]]:
    sets: list[tuple[str, np.ndarray, np.ndarray, dict]] = []
    for protocol in parse_csv(args.fraction_protocols):
        for frac in parse_float_csv(args.support_fracs):
            if protocol == "natural":
                support_mask, query_mask, extra = make_natural_frac_support(records, heldout_mask, seed, fold, frac)
            elif protocol == "stratified":
                support_mask, query_mask, extra = make_stratified_frac_support(records, heldout_mask, seed, fold, frac)
            elif protocol == "active_uncertainty_diverse":
                support_mask, query_mask, extra = active_uncertainty_diverse_support(records, scalars_raw, heldout_mask, zero_score, seed, fold, frac)
            elif protocol == "active_prediction_balanced":
                support_mask, query_mask, extra = active_prediction_balanced_support(records, heldout_mask, zero_score, seed, fold, frac)
            else:
                raise ValueError(f"Unknown fraction protocol: {protocol}")
            sets.append((f"{extra['fewshot_support_protocol']}_f{frac:g}", support_mask, query_mask, extra))
    for k_per_class in parse_int_csv(args.fixed_ks):
        support_mask, query_mask, extra = make_fixedk_balanced_support(records, heldout_mask, seed, fold, k_per_class)
        sets.append((f"fixedk_balanced_k{k_per_class}", support_mask, query_mask, extra))
    return sets


def run_source(args: argparse.Namespace, source: str) -> pd.DataFrame:
    pair_csv = {
        "chembl_raw": "data/chembl_raw_moleculeace_target_pairs.csv",
        "bindingdb_raw": "data/bindingdb_raw_moleculeace_target_pairs.csv",
    }[source]
    records = pd.read_csv(pair_csv)
    diff_bits, scalars_raw = raw.build_feature_matrices(records)
    seq_raw = anchor.target_sequence_matrix(records, args.sequence_hash_features)
    latent_spaces = latent_anchor.make_latent_spaces(records, seq_raw, args)
    y_all = records["label"].to_numpy(dtype=np.int32)
    rows = []
    for seed in parse_int_csv(args.seeds):
        masks, target_to_family = raw.family_folds(records, seed, args.folds, args.sequence_hash_features, args.target_family_clusters)
        family_arr = records["target"].astype(str).map(target_to_family).to_numpy(dtype=np.int32)
        for fold, heldout_mask in enumerate(masks[: args.folds]):
            base_train_mask, _ = raw.purge_train_mask(records, ~heldout_mask, heldout_mask, args.split_mode)
            if len(np.unique(y_all[base_train_mask])) < 2 or len(np.unique(y_all[heldout_mask])) < 2:
                continue
            start = time.perf_counter()
            zero_tree = multirf.fit_tree_fusions(scalars_raw, records, base_train_mask, heldout_mask, seed, args.threads, args.n_estimators)
            zero_robust_rank, _ = anchor.robust_predictions(records, diff_bits, scalars_raw, family_arr, base_train_mask, heldout_mask, seed, args.threads, args.n_estimators)
            zero_elapsed = time.perf_counter() - start
            zero_tree_prob = source_anchor.clip_prob(zero_tree["tree_prob_fusion_scalars"][0])
            zero_tree_rank = source_anchor.rank01(zero_tree["multirf_rank_fusion_scalars"][0])
            zero_score = (0.5 * zero_tree_rank + 0.5 * np.asarray(zero_robust_rank, dtype=np.float32)).astype(np.float32)
            zero_extra = {
                "fewshot_support_protocol": "zero_source_only",
                "fewshot_requested_frac": 0.0,
                "fewshot_support_rows": 0,
                "fewshot_query_rows": int(heldout_mask.sum()),
                "fewshot_support_frac_actual": 0.0,
                "fewshot_support_positive_rate": 0.0,
                "fewshot_query_positive_rate": float(y_all[heldout_mask].mean()),
            }
            add_row(rows, records, family_arr, base_train_mask, heldout_mask, source, f"{args.split_mode}_support_protocols", "zero_source_only_base_rank", seed, fold, zero_score, zero_tree_prob, zero_elapsed, zero_extra)

            support_sets = make_support_sets(args, records, scalars_raw, heldout_mask, zero_score, seed, fold)
            for protocol_name, support_mask, query_mask, protocol_extra in support_sets:
                if len(np.unique(y_all[query_mask])) < 2:
                    print(f"[SKIP {source} seed={seed} fold={fold} {protocol_name}] query single-class", flush=True)
                    continue
                train_mask = base_train_mask | support_mask
                if len(np.unique(y_all[train_mask])) < 2:
                    continue
                extra = support_extra(records, y_all, support_mask, query_mask, protocol_extra)

                support_knn = support_knn_prob(scalars_raw, y_all, support_mask, query_mask, args.support_knn_k)
                if support_knn is not None:
                    add_row(rows, records, family_arr, train_mask, query_mask, source, f"{args.split_mode}_support_protocols", f"{protocol_name}_support_only_knn_k{args.support_knn_k}", seed, fold, source_anchor.rank01(support_knn), support_knn, 0.0, extra)

                support_et = support_extratrees_prob(scalars_raw, y_all, support_mask, query_mask, seed, args.threads, args.n_estimators)
                if support_et is not None:
                    add_row(rows, records, family_arr, train_mask, query_mask, source, f"{args.split_mode}_support_protocols", f"{protocol_name}_support_only_extratrees", seed, fold, source_anchor.rank01(support_et), support_et, 0.0, extra)

                support_lr = support_logreg_prob(scalars_raw, y_all, support_mask, query_mask, seed)
                if support_lr is not None:
                    add_row(rows, records, family_arr, train_mask, query_mask, source, f"{args.split_mode}_support_protocols", f"{protocol_name}_support_only_logreg", seed, fold, source_anchor.rank01(support_lr), support_lr, 0.0, extra)

                start = time.perf_counter()
                tree_preds = multirf.fit_tree_fusions(scalars_raw, records, train_mask, query_mask, seed, args.threads, args.n_estimators)
                robust_rank, _ = anchor.robust_predictions(records, diff_bits, scalars_raw, family_arr, train_mask, query_mask, seed, args.threads, args.n_estimators)
                elapsed = time.perf_counter() - start
                tree_prob = source_anchor.clip_prob(tree_preds["tree_prob_fusion_scalars"][0])
                tree_rank = source_anchor.rank01(tree_preds["multirf_rank_fusion_scalars"][0])
                base_rank = (0.5 * tree_rank + 0.5 * np.asarray(robust_rank, dtype=np.float32)).astype(np.float32)
                add_row(rows, records, family_arr, train_mask, query_mask, source, f"{args.split_mode}_support_protocols", f"{protocol_name}_source_plus_support_base_rank", seed, fold, base_rank, tree_prob, elapsed, extra)

                if support_knn is not None:
                    knn_rank = source_anchor.rank01(support_knn)
                    for blend in parse_float_csv(args.knn_blend_weights):
                        blend_score = ((1.0 - blend) * base_rank + blend * knn_rank).astype(np.float32)
                        blend_prob = source_anchor.clip_prob((1.0 - blend) * tree_prob + blend * support_knn)
                        add_row(rows, records, family_arr, train_mask, query_mask, source, f"{args.split_mode}_support_protocols", f"{protocol_name}_base_plus_support_knn_b{blend:g}", seed, fold, blend_score, blend_prob, 0.0, extra)

                for mode in parse_csv(args.anchor_modes):
                    for k in parse_int_csv(args.anchor_ks):
                        try:
                            a_prob, support, a_extra = latent_anchor.latent_anchor_probability(
                                mode, diff_bits, latent_spaces, y_all, train_mask, query_mask, k, args
                            )
                        except Exception as exc:
                            print(f"[WARN fair anchor {mode}] {type(exc).__name__}: {exc}", flush=True)
                            continue
                        a_rank = source_anchor.rank01(a_prob)
                        for wa in parse_float_csv(args.anchor_weights):
                            score = ((1.0 - wa) * base_rank + wa * a_rank).astype(np.float32)
                            add_row(
                                rows,
                                records,
                                family_arr,
                                train_mask,
                                query_mask,
                                source,
                                f"{args.split_mode}_support_protocols",
                                f"{protocol_name}_latent_anchor_{mode}_k{k}_a{wa:g}",
                                seed,
                                fold,
                                score,
                                tree_prob,
                                0.0,
                                {**extra, **a_extra, "anchor_support_mean": float(np.mean(support))},
                            )
                print(f"[{source} seed={seed} fold={fold} {protocol_name}] support={support_mask.sum()} query={query_mask.sum()}", flush=True)
    return pd.DataFrame(rows)


def summarize(results: pd.DataFrame, prefix: Path) -> None:
    results.to_json(f"{prefix}.jsonl", orient="records", lines=True)
    metric_cols = [
        c for c in list(raw.METRICS) + [
            "train_seconds", "predict_seconds", "n_train", "n_test", "fewshot_requested_frac",
            "fewshot_support_frac_actual", "fewshot_support_rows", "fewshot_query_rows",
            "fewshot_support_positive_rate", "fewshot_query_positive_rate", "fewshot_support_k_per_class",
            "anchor_support_mean", "latent_anchor_coverage_mean", "latent_anchor_n_spaces",
            "target_overlap_rate", "family_overlap_rate", "exact_ligand_overlap_rate",
            "scaffold_overlap_rate", "document_source_overlap_rate", "assay_source_overlap_rate",
        ] if c in results.columns
    ]
    summary = results.groupby(["benchmark", "split_mode", "fewshot_support_protocol", "variant"], as_index=False)[metric_cols].agg(["mean", "std", "count"])
    summary.columns = ["_".join([x for x in col if x]).strip("_") for col in summary.columns.to_flat_index()]
    summary = summary.reset_index(drop=True)
    summary.to_csv(f"{prefix}_summary.csv", index=False)
    raw.paired_deltas(results).to_csv(f"{prefix}_paired_deltas.csv", index=False)
    lines = ["# Few-Shot Support Fairness, Baselines, and Active Label Efficiency", ""]
    for benchmark in sorted(summary.benchmark.unique()):
        sub = summary[summary.benchmark == benchmark].sort_values(["roc_auc_mean", "pr_auc_mean", "mcc_mean"], ascending=False)
        lines.append(f"## {benchmark}")
        keep = [
            "benchmark", "fewshot_support_protocol", "variant", "roc_auc_mean", "roc_auc_std",
            "pr_auc_mean", "pr_auc_std", "mcc_mean", "mcc_std", "fewshot_requested_frac_mean",
            "fewshot_support_frac_actual_mean", "fewshot_support_rows_mean", "fewshot_query_rows_mean",
            "roc_auc_count",
        ]
        lines.append(raw.markdown_table(sub[[c for c in keep if c in sub.columns]].head(80)))
        lines.append("")
    Path(f"{prefix}_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sources", default="bindingdb_raw,chembl_raw")
    parser.add_argument("--split-mode", default="family_scaffold_source_purged")
    parser.add_argument("--seeds", default="0,1,2,3,4")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--target-family-clusters", type=int, default=10)
    parser.add_argument("--sequence-hash-features", type=int, default=128)
    parser.add_argument("--n-estimators", type=int, default=160)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--embedding-dir", default="data/pretrained_embeddings")
    parser.add_argument("--embedding-names", default="molformer,chemberta")
    parser.add_argument("--embedding-pair-modes", default="absdiff")
    parser.add_argument("--include-target-latent", action="store_true")
    parser.add_argument("--target-latent-weight", type=float, default=0.35)
    parser.add_argument("--anchor-modes", default="diffbit+molformer_absdiff+chemberta_absdiff")
    parser.add_argument("--anchor-ks", default="16")
    parser.add_argument("--anchor-shrink", type=float, default=1.0)
    parser.add_argument("--anchor-bandwidth", default="median")
    parser.add_argument("--anchor-weights", default="0.35,0.5")
    parser.add_argument("--fraction-protocols", default="natural,active_uncertainty_diverse")
    parser.add_argument("--support-fracs", default="0.05,0.1,0.2,0.3,0.5")
    parser.add_argument("--fixed-ks", default="16,32,64")
    parser.add_argument("--support-knn-k", type=int, default=7)
    parser.add_argument("--knn-blend-weights", default="0.25")
    parser.add_argument("--output-prefix", default="results/tkde_fewshot_fairness_active_baselines")
    args = parser.parse_args()

    frames = [run_source(args, source) for source in parse_csv(args.sources)]
    results = pd.concat(frames, ignore_index=True)
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    summarize(results, prefix)
    print(json.dumps({"rows": int(len(results)), "output_prefix": str(prefix)}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()

