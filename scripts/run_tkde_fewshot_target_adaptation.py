from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

import run_chembl_raw_source_temporal_pairs as raw
import run_literature_guided_multirf_fusion as multirf
import run_tkde_latent_anchor_portfolio as latent_anchor
import run_tkde_source_anchor_portfolio as source_anchor
import run_tnnls_multi_anchor_local_evidence as anchor


def parse_csv(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()]


def make_support_query_masks(records: pd.DataFrame, test_mask: np.ndarray, seed: int, fold: int, frac: float) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed * 1009 + fold * 9173 + int(round(frac * 1000)))
    support_mask = np.zeros(len(records), dtype=bool)
    test_idx = np.flatnonzero(test_mask)
    y = records.loc[test_mask, "label"].to_numpy(dtype=np.int32)
    local_by_label = {lab: test_idx[np.flatnonzero(y == lab)] for lab in [0, 1]}
    for lab, idx in local_by_label.items():
        if len(idx) == 0:
            continue
        n = max(1, int(round(len(idx) * float(frac))))
        n = min(n, max(len(idx) - 1, 1))
        chosen = rng.choice(idx, size=n, replace=False)
        support_mask[chosen] = True
    query_mask = test_mask & ~support_mask
    return support_mask, query_mask


def add_row(rows, records, family_arr, train_mask, query_mask, source, split_mode, variant, seed, fold, score, prob, seconds, extra=None):
    source_anchor.add_row(rows, records, family_arr, train_mask, query_mask, source, split_mode, variant, seed, fold, score, prob, seconds, extra or {})


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
    for seed in [int(x) for x in parse_csv(args.seeds)]:
        masks, target_to_family = raw.family_folds(records, seed, args.folds, args.sequence_hash_features, args.target_family_clusters)
        family_arr = records["target"].astype(str).map(target_to_family).to_numpy(dtype=np.int32)
        for fold, heldout_mask in enumerate(masks[: args.folds]):
            base_train_mask, _ = raw.purge_train_mask(records, ~heldout_mask, heldout_mask, args.split_mode)
            for frac in [float(x) for x in parse_csv(args.support_fracs)]:
                support_mask, query_mask = make_support_query_masks(records, heldout_mask, seed, fold, frac)
                train_mask = base_train_mask | support_mask
                if len(np.unique(y_all[train_mask])) < 2 or len(np.unique(y_all[query_mask])) < 2:
                    continue
                extra = {
                    "fewshot_support_frac": float(frac),
                    "fewshot_support_rows": int(support_mask.sum()),
                    "fewshot_query_rows": int(query_mask.sum()),
                    "fewshot_support_positive_rate": float(y_all[support_mask].mean()) if support_mask.any() else 0.0,
                }
                start = time.perf_counter()
                tree_preds = multirf.fit_tree_fusions(scalars_raw, records, train_mask, query_mask, seed, args.threads, args.n_estimators)
                robust_rank, _ = anchor.robust_predictions(records, diff_bits, scalars_raw, family_arr, train_mask, query_mask, seed, args.threads, args.n_estimators)
                elapsed = time.perf_counter() - start
                tree_prob = source_anchor.clip_prob(tree_preds["tree_prob_fusion_scalars"][0])
                tree_rank = source_anchor.rank01(tree_preds["multirf_rank_fusion_scalars"][0])
                base_rank = (0.5 * tree_rank + 0.5 * np.asarray(robust_rank, dtype=np.float32)).astype(np.float32)
                add_row(rows, records, family_arr, train_mask, query_mask, source, f"{args.split_mode}_fewshot", f"fewshot_base_rank_f{frac:g}", seed, fold, base_rank, tree_prob, elapsed, extra)
                for mode in parse_csv(args.anchor_modes):
                    for k in [int(x) for x in parse_csv(args.anchor_ks)]:
                        try:
                            a_prob, support, a_extra = latent_anchor.latent_anchor_probability(
                                mode, diff_bits, latent_spaces, y_all, train_mask, query_mask, k, args
                            )
                        except Exception as exc:
                            print(f"[WARN fewshot anchor {mode}] {type(exc).__name__}: {exc}", flush=True)
                            continue
                        a_rank = source_anchor.rank01(a_prob)
                        for wa in [float(x) for x in parse_csv(args.anchor_weights)]:
                            score = ((1.0 - wa) * base_rank + wa * a_rank).astype(np.float32)
                            add_row(
                                rows,
                                records,
                                family_arr,
                                train_mask,
                                query_mask,
                                source,
                                f"{args.split_mode}_fewshot",
                                f"fewshot_latent_anchor_{mode}_k{k}_a{wa:g}_f{frac:g}",
                                seed,
                                fold,
                                score,
                                tree_prob,
                                0.0,
                                {**extra, **a_extra, "anchor_support_mean": float(np.mean(support))},
                            )
                print(f"[{source} seed={seed} fold={fold} frac={frac:g}] support={support_mask.sum()} query={query_mask.sum()}", flush=True)
    return pd.DataFrame(rows)


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
    parser.add_argument("--anchor-modes", default="diffbit+molformer_absdiff,diffbit+molformer_absdiff+chemberta_absdiff")
    parser.add_argument("--anchor-ks", default="16,64")
    parser.add_argument("--anchor-shrink", type=float, default=1.0)
    parser.add_argument("--anchor-bandwidth", default="median")
    parser.add_argument("--anchor-weights", default="0.2,0.35,0.5")
    parser.add_argument("--support-fracs", default="0.05,0.1,0.2")
    parser.add_argument("--output-prefix", default="results/tkde_fewshot_target_adaptation")
    args = parser.parse_args()

    frames = [run_source(args, source) for source in parse_csv(args.sources)]
    results = pd.concat(frames, ignore_index=True)
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    results.to_json(f"{prefix}.jsonl", orient="records", lines=True)
    metric_cols = [
        c for c in list(raw.METRICS) + [
            "train_seconds", "predict_seconds", "n_train", "n_test", "fewshot_support_frac",
            "fewshot_support_rows", "fewshot_query_rows", "fewshot_support_positive_rate",
            "anchor_support_mean", "latent_anchor_coverage_mean", "latent_anchor_n_spaces",
            "target_overlap_rate", "family_overlap_rate", "exact_ligand_overlap_rate",
            "scaffold_overlap_rate", "document_source_overlap_rate", "assay_source_overlap_rate",
        ] if c in results.columns
    ]
    summary = results.groupby(["benchmark", "split_mode", "variant"], as_index=False)[metric_cols].agg(["mean", "std", "count"])
    summary.columns = ["_".join([x for x in col if x]).strip("_") for col in summary.columns.to_flat_index()]
    summary = summary.reset_index(drop=True)
    summary.to_csv(f"{prefix}_summary.csv", index=False)
    raw.paired_deltas(results).to_csv(f"{prefix}_paired_deltas.csv", index=False)
    Path(f"{prefix}_report.md").write_text(
        "# Few-Shot Target Adaptation\n\n"
        + raw.markdown_table(summary.sort_values(["benchmark", "roc_auc_mean", "pr_auc_mean"], ascending=[True, False, False]).head(120))
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"rows": int(len(results)), "output_prefix": str(prefix)}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()

