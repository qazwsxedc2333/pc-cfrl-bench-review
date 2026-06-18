from __future__ import annotations

import argparse
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

import run_chembl_raw_source_temporal_pairs as raw


RECORDS: pd.DataFrame | None = None
FEATURES: dict[str, np.ndarray | sparse.csr_matrix] = {}
SCALARS: np.ndarray | None = None


def parse_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_sources(value: str) -> list[tuple[str, str]]:
    out = []
    for item in parse_csv(value):
        name, path = item.split("=", 1)
        out.append((name.strip(), path.strip()))
    return out


def load_embedding_cache(path: Path) -> tuple[dict[str, int], np.ndarray]:
    data = np.load(path, allow_pickle=True)
    smiles = data["smiles"].astype(str)
    matrix = data["embeddings"].astype(np.float32)
    return {smile: idx for idx, smile in enumerate(smiles)}, matrix


def build_embedding_pair_features(records: pd.DataFrame, cache_path: Path) -> np.ndarray:
    lookup, embeddings = load_embedding_cache(cache_path)
    left = np.vstack([embeddings[lookup[str(smile)]] for smile in records["smiles1"]])
    right = np.vstack([embeddings[lookup[str(smile)]] for smile in records["smiles2"]])
    return np.hstack([np.abs(left - right), left * right]).astype(np.float32)


def prepare_features(records: pd.DataFrame, embedding_dir: Path, embeddings: list[str]) -> None:
    global FEATURES, SCALARS
    FEATURES = {}
    diff_bits, scalars = raw.build_feature_matrices(records)
    FEATURES["ecfp_absdiff_xgb"] = diff_bits
    SCALARS = scalars
    for alias in embeddings:
        path = embedding_dir / f"{alias}.npz"
        if path.exists():
            FEATURES[f"{alias}_sym_xgb"] = build_embedding_pair_features(records, path)


def random_fold_masks(records: pd.DataFrame, seed: int, folds: int) -> list[np.ndarray]:
    splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    y = records["label"].to_numpy(dtype=np.int32)
    masks = []
    for _, test_idx in splitter.split(np.zeros(len(records)), y):
        mask = np.zeros(len(records), dtype=bool)
        mask[test_idx] = True
        masks.append(mask)
    return masks


def split_masks(records: pd.DataFrame, split_mode: str, seed: int, folds: int, clusters: int) -> list[np.ndarray]:
    if split_mode == "random":
        return random_fold_masks(records, seed, folds)
    masks, _ = raw.family_folds(records, seed, folds, 128, clusters)
    return masks


def make_model(variant: str, seed: int, y_train: np.ndarray):
    positives = float(y_train.sum())
    negatives = float(len(y_train) - positives)
    weight = negatives / positives if positives > 0 else 1.0
    if variant in {"rf_scalars", "extratrees_scalars"}:
        cls = RandomForestClassifier if variant == "rf_scalars" else ExtraTreesClassifier
        return cls(
            n_estimators=300,
            max_features="sqrt",
            min_samples_leaf=2,
            class_weight="balanced",
            n_jobs=1,
            random_state=seed,
        )
    return XGBClassifier(
        n_estimators=120,
        max_depth=4,
        learning_rate=0.04,
        subsample=0.85,
        colsample_bytree=0.8,
        reg_lambda=4.0,
        objective="binary:logistic",
        eval_metric="aucpr",
        tree_method="hist",
        device="cpu",
        n_jobs=1,
        random_state=seed,
        scale_pos_weight=weight,
    )


def build_x(variant: str, train_mask: np.ndarray):
    assert SCALARS is not None
    if variant == "ecfp_absdiff_xgb" or variant.endswith("_sym_xgb"):
        return FEATURES[variant]
    scaler = StandardScaler().fit(SCALARS[train_mask])
    scalar_matrix = scaler.transform(SCALARS).astype(np.float32)
    if variant in {"rf_scalars", "extratrees_scalars", "pccfrl_scalar_xgb"}:
        return scalar_matrix
    if variant == "pccfrl_xgb":
        return sparse.hstack(
            [FEATURES["ecfp_absdiff_xgb"], sparse.csr_matrix(scalar_matrix)],
            format="csr",
        )
    raise ValueError(variant)


def metrics(y: np.ndarray, prob: np.ndarray) -> dict[str, float]:
    pred = (prob >= 0.5).astype(np.int32)
    return {
        "roc_auc": float(roc_auc_score(y, prob)) if len(np.unique(y)) == 2 else float("nan"),
        "pr_auc": float(average_precision_score(y, prob)) if len(np.unique(y)) == 2 else float("nan"),
        "mcc": float(matthews_corrcoef(y, pred))
        if len(np.unique(y)) > 1 and len(np.unique(pred)) > 1
        else 0.0,
    }


def run_task(task: dict) -> tuple[dict, list[dict]]:
    assert RECORDS is not None
    records = RECORDS
    split_mode = task["split_mode"]
    masks = split_masks(records, split_mode, task["seed"], task["folds"], task["clusters"])
    test_mask = masks[task["fold"]]
    base_train = ~test_mask
    if split_mode == "random":
        train_mask, retention = base_train, 1.0
    elif split_mode == "target_family":
        train_mask, retention = base_train, 1.0
    else:
        train_mask, retention = raw.purge_train_mask(records, base_train, test_mask, split_mode)
    y_train = records.loc[train_mask, "label"].to_numpy(dtype=np.int32)
    y_test = records.loc[test_mask, "label"].to_numpy(dtype=np.int32)
    if (
        train_mask.sum() < task["min_train"]
        or test_mask.sum() < task["min_test"]
        or len(np.unique(y_train)) < 2
        or len(np.unique(y_test)) < 2
    ):
        return {}, []
    x = build_x(task["variant"], train_mask)
    model = make_model(task["variant"], task["seed"], y_train)
    model.fit(x[train_mask], y_train)
    prob = model.predict_proba(x[test_mask])[:, 1].astype(np.float32)
    row = {
        "benchmark": task["source"],
        "split_mode": split_mode,
        "variant": task["variant"],
        "seed": task["seed"],
        "fold": task["fold"],
        "n_train": int(train_mask.sum()),
        "n_test": int(test_mask.sum()),
        "train_retention_after_purge": float(retention),
        **raw.audit_overlaps(records, train_mask, test_mask),
        **metrics(y_test, prob),
    }
    target_rows = []
    test_indices = np.flatnonzero(test_mask)
    test_targets = records.iloc[test_indices]["target"].astype(str).to_numpy()
    for target in sorted(set(test_targets)):
        local = test_targets == target
        if local.sum() < task["min_target_rows"] or len(np.unique(y_test[local])) < 2:
            continue
        target_rows.append(
            {
                "benchmark": task["source"],
                "split_mode": split_mode,
                "variant": task["variant"],
                "seed": task["seed"],
                "fold": task["fold"],
                "target": target,
                "n_rows": int(local.sum()),
                "positive_rate": float(y_test[local].mean()),
                **metrics(y_test[local], prob[local]),
            }
        )
    return row, target_rows


def run_source(args: argparse.Namespace, source: str, pair_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    global RECORDS
    RECORDS = pd.read_csv(pair_path)
    prepare_features(RECORDS, Path(args.embedding_dir), parse_csv(args.embeddings))
    variants = parse_csv(args.variants)
    variants += [name for name in FEATURES if name.endswith("_sym_xgb")]
    variants = list(dict.fromkeys(variants))
    tasks = []
    for split_mode in parse_csv(args.split_modes):
        for seed in [int(x) for x in parse_csv(args.seeds)]:
            n_masks = len(split_masks(RECORDS, split_mode, seed, args.folds, args.target_family_clusters))
            for fold in range(n_masks):
                for variant in variants:
                    tasks.append(
                        {
                            "source": source,
                            "split_mode": split_mode,
                            "seed": seed,
                            "fold": fold,
                            "folds": args.folds,
                            "clusters": args.target_family_clusters,
                            "variant": variant,
                            "min_train": args.min_train_rows,
                            "min_test": args.min_test_rows,
                            "min_target_rows": args.min_target_rows,
                        }
                    )
    fold_rows: list[dict] = []
    target_rows: list[dict] = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run_task, task) for task in tasks]
        for future in as_completed(futures):
            row, targets = future.result()
            if not row:
                continue
            fold_rows.append(row)
            target_rows.extend(targets)
            print(
                f"[{source} {row['split_mode']} {row['variant']} seed={row['seed']} fold={row['fold']}] "
                f"roc={row['roc_auc']:.4f} pr={row['pr_auc']:.4f}",
                flush=True,
            )
    return pd.DataFrame(fold_rows), pd.DataFrame(target_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sources",
        default=(
            "chembl_raw=data/chembl_raw_moleculeace_target_pairs.csv,"
            "bindingdb_raw=data/bindingdb_raw_moleculeace_target_pairs.csv"
        ),
    )
    parser.add_argument("--embedding-dir", default="data/pretrained_embeddings")
    parser.add_argument("--embeddings", default="chemberta,molformer")
    parser.add_argument(
        "--variants",
        default="ecfp_absdiff_xgb,pccfrl_xgb,pccfrl_scalar_xgb,rf_scalars,extratrees_scalars",
    )
    parser.add_argument(
        "--split-modes",
        default="random,target_family,family_scaffold_source_purged",
    )
    parser.add_argument("--seeds", default="0,1,2,3,4")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--target-family-clusters", type=int, default=10)
    parser.add_argument("--workers", type=int, default=64)
    parser.add_argument("--min-train-rows", type=int, default=300)
    parser.add_argument("--min-test-rows", type=int, default=100)
    parser.add_argument("--min-target-rows", type=int, default=30)
    parser.add_argument("--output-prefix", default="results/tkde_unified_baseline_matrix")
    args = parser.parse_args()

    all_folds = []
    all_targets = []
    for source, path in parse_sources(args.sources):
        folds, targets = run_source(args, source, path)
        all_folds.append(folds)
        all_targets.append(targets)
    fold_df = pd.concat(all_folds, ignore_index=True)
    target_df = pd.concat(all_targets, ignore_index=True)
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fold_df.to_csv(f"{prefix}_fold_metrics.csv", index=False)
    target_df.to_csv(f"{prefix}_target_metrics.csv", index=False)
    summary = (
        fold_df.groupby(["benchmark", "split_mode", "variant"], as_index=False)[
            ["roc_auc", "pr_auc", "mcc", "n_train", "n_test", "train_retention_after_purge"]
        ]
        .agg(["mean", "std", "count"])
    )
    summary.columns = ["_".join([x for x in col if x]) for col in summary.columns.to_flat_index()]
    summary.reset_index(drop=True).to_csv(f"{prefix}_summary.csv", index=False)
    print(json.dumps({"fold_rows": len(fold_df), "target_rows": len(target_df)}), flush=True)


if __name__ == "__main__":
    main()
