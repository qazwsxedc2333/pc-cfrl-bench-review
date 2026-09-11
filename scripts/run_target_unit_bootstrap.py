from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from run_chembl_raw_source_temporal_pairs import (
    build_feature_matrices,
    build_x,
    family_folds,
    make_model,
    purge_train_mask,
)


def parse_ints(value: str) -> list[int]:
    return [int(item.strip()) for item in value.split(",") if item.strip()]


def generate_oof_predictions(args: argparse.Namespace, records: pd.DataFrame) -> pd.DataFrame:
    diff_bits, scalars = build_feature_matrices(records)
    prediction_rows: list[pd.DataFrame] = []

    for seed in parse_ints(args.seeds):
        masks, target_to_cluster = family_folds(
            records,
            seed,
            args.folds,
            args.sequence_hash_features,
            args.target_clusters,
        )
        target_cluster = records["target"].astype(str).map(target_to_cluster).to_numpy(dtype=np.int32)

        for fold, test_mask in enumerate(masks):
            base_train_mask = ~test_mask
            train_mask, _ = purge_train_mask(
                records,
                base_train_mask,
                test_mask,
                "family_scaffold_source_purged",
            )
            y_train = records.loc[train_mask, "label"].to_numpy(dtype=np.int32)
            y_test = records.loc[test_mask, "label"].to_numpy(dtype=np.int32)
            if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
                continue

            scores: dict[str, np.ndarray] = {}
            for variant, column in (
                ("ecfp_absdiff_xgb", "ecfp_score"),
                ("rich_diff_scalars", "pc_cfrl_full_score"),
            ):
                matrix = build_x(diff_bits, scalars, train_mask, variant)
                model = make_model(seed, y_train, args.threads, args.n_estimators, args.max_depth)
                model.fit(matrix[train_mask], y_train)
                scores[column] = model.predict_proba(matrix[test_mask])[:, 1].astype(np.float32)

            test_rows = records.loc[test_mask]
            biological_target = (
                test_rows["protein_target"].astype(str)
                if "protein_target" in records.columns
                else test_rows["target"].astype(str)
            )
            prediction_rows.append(
                pd.DataFrame(
                    {
                        "dataset": args.dataset,
                        "seed": seed,
                        "fold": fold,
                        "row_index": np.flatnonzero(test_mask),
                        "target": test_rows["target"].astype(str).to_numpy(),
                        "biological_target": biological_target.to_numpy(),
                        "target_cluster": target_cluster[test_mask],
                        "label": y_test,
                        "ecfp_score": scores["ecfp_score"],
                        "pc_cfrl_full_score": scores["pc_cfrl_full_score"],
                    }
                )
            )

    if not prediction_rows:
        raise RuntimeError("No valid out-of-fold predictions were generated")
    return pd.concat(prediction_rows, ignore_index=True)


def cluster_bootstrap(
    averaged: pd.DataFrame,
    unit_column: str,
    n_bootstrap: int,
    seed: int,
) -> dict[str, float | int | str]:
    units = averaged[unit_column].astype(str).unique()
    unit_indices = {
        unit: np.flatnonzero(averaged[unit_column].astype(str).to_numpy() == unit)
        for unit in units
    }
    y = averaged["label"].to_numpy(dtype=np.int32)
    ecfp = averaged["ecfp_score"].to_numpy(dtype=float)
    full = averaged["pc_cfrl_full_score"].to_numpy(dtype=float)
    observed_ecfp = float(roc_auc_score(y, ecfp))
    observed_full = float(roc_auc_score(y, full))

    rng = np.random.default_rng(seed)
    deltas: list[float] = []
    for _ in range(n_bootstrap):
        sampled_units = rng.choice(units, size=len(units), replace=True)
        sampled_indices = np.concatenate([unit_indices[str(unit)] for unit in sampled_units])
        sampled_y = y[sampled_indices]
        if len(np.unique(sampled_y)) < 2:
            continue
        deltas.append(
            float(
                roc_auc_score(sampled_y, full[sampled_indices])
                - roc_auc_score(sampled_y, ecfp[sampled_indices])
            )
        )

    if not deltas:
        raise RuntimeError(f"No valid bootstrap replicates for {unit_column}")
    delta_array = np.asarray(deltas, dtype=float)
    return {
        "dataset": str(averaged["dataset"].iloc[0]),
        "resampling_unit": unit_column,
        "n_units": int(len(units)),
        "n_rows": int(len(averaged)),
        "n_bootstrap": int(len(delta_array)),
        "ecfp_roc_auc": observed_ecfp,
        "pc_cfrl_full_roc_auc": observed_full,
        "delta_roc_auc": observed_full - observed_ecfp,
        "ci_2_5": float(np.quantile(delta_array, 0.025)),
        "ci_97_5": float(np.quantile(delta_array, 0.975)),
        "bootstrap_positive_rate": float(np.mean(delta_array > 0)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate paired OOF predictions and bootstrap ROC-AUC differences over biological units."
    )
    parser.add_argument("--dataset", choices=("chembl", "bindingdb"), required=True)
    parser.add_argument("--pair-csv", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--seeds", default="0,1,2,3,4")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--target-clusters", type=int, default=10)
    parser.add_argument("--sequence-hash-features", type=int, default=128)
    parser.add_argument("--n-estimators", type=int, default=80)
    parser.add_argument("--max-depth", type=int, default=4)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260911)
    args = parser.parse_args()

    records = pd.read_csv(args.pair_csv)
    predictions = generate_oof_predictions(args, records)
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(prefix.with_name(prefix.name + "_oof_predictions.csv"), index=False)

    averaged = (
        predictions.groupby(
            ["dataset", "row_index", "target", "biological_target", "target_cluster", "label"],
            as_index=False,
        )[["ecfp_score", "pc_cfrl_full_score"]]
        .mean()
    )
    summaries = [
        cluster_bootstrap(averaged, "biological_target", args.n_bootstrap, args.bootstrap_seed),
        cluster_bootstrap(averaged, "target_cluster", args.n_bootstrap, args.bootstrap_seed + 1),
    ]
    summary_path = prefix.with_name(prefix.name + "_summary.csv")
    pd.DataFrame(summaries).to_csv(summary_path, index=False)
    manifest = {
        "dataset": args.dataset,
        "pair_csv": args.pair_csv,
        "split_mode": "family_scaffold_source_purged",
        "variants": ["ecfp_absdiff_xgb", "rich_diff_scalars"],
        "seeds": parse_ints(args.seeds),
        "folds": args.folds,
        "target_clusters": args.target_clusters,
        "n_bootstrap": args.n_bootstrap,
        "prediction_file": str(prefix.with_name(prefix.name + "_oof_predictions.csv")),
        "summary_file": str(summary_path),
    }
    prefix.with_name(prefix.name + "_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
