from __future__ import annotations

import argparse
import hashlib
import json
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy import sparse
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import ExtraTreesClassifier, ExtraTreesRegressor, RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    matthews_corrcoef,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier, XGBRegressor

import run_chembl_raw_source_temporal_pairs as raw
import run_tkde_unified_baseline_matrix as unified


RECORDS: pd.DataFrame | None = None
ORIGINAL_RECORDS: pd.DataFrame | None = None
SYMMETRIC_FEATURES: dict[str, np.ndarray | sparse.csr_matrix] = {}
DIRECTED_FEATURES: dict[str, np.ndarray | sparse.csr_matrix] = {}
SCALARS: np.ndarray | None = None
DIRECTED_SCALARS: np.ndarray | None = None
FP_CACHE: dict[str, object] = {}
FP_ARRAY_CACHE: dict[str, np.ndarray] = {}
LIGAND_INDEX: dict[str, int] = {}
LIGAND_FP_MATRIX: sparse.csr_matrix | None = None
LIGAND_FP_POPCOUNT: np.ndarray | None = None
LIGAND_SIM_MATRIX: sparse.csr_matrix | None = None


def parse_csv(value: str) -> list[str]:
    return [item.strip() for item in str(value).split(",") if item.strip()]


def parse_sources(value: str) -> list[tuple[str, str, str | None]]:
    out = []
    for item in parse_csv(value):
        parts = item.split("=")
        if len(parts) == 2:
            out.append((parts[0].strip(), parts[1].strip(), None))
        elif len(parts) == 3:
            out.append((parts[0].strip(), parts[1].strip(), parts[2].strip()))
        else:
            raise ValueError(f"Bad source spec: {item}")
    return out


def stable_int(text: str, seed: int = 0) -> int:
    digest = hashlib.blake2b(f"{seed}:{text}".encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "little")


def fp_obj(smiles: str):
    if smiles not in FP_CACHE:
        mol = Chem.MolFromSmiles(str(smiles))
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        FP_CACHE[smiles] = gen.GetFingerprint(mol) if mol is not None else None
    return FP_CACHE[smiles]


def fp_array(smiles: str) -> np.ndarray:
    if smiles not in FP_ARRAY_CACHE:
        arr = np.zeros((2048,), dtype=np.float32)
        fp = fp_obj(smiles)
        if fp is not None:
            DataStructs.ConvertToNumpyArray(fp, arr)
        FP_ARRAY_CACHE[smiles] = arr
    return FP_ARRAY_CACHE[smiles]


def make_activity_lookup(source: str, activity_path: str) -> pd.Series:
    activity = pd.read_csv(activity_path)
    if source.startswith("bindingdb"):
        activity["pair_target"] = activity["uniprot_accession"].astype(str) + "::" + activity["standard_type"].astype(str)
    else:
        activity["pair_target"] = activity["target_chembl_id"].astype(str)
    activity["canonical_smiles"] = activity["canonical_smiles"].astype(str)
    activity["pchembl_value"] = pd.to_numeric(activity["pchembl_value"], errors="coerce")
    activity = activity[np.isfinite(activity["pchembl_value"])]
    lookup = activity.groupby(["pair_target", "canonical_smiles"], sort=False)["pchembl_value"].median()
    return lookup


def augment_pairs(source: str, pair_path: str, activity_path: str | None) -> pd.DataFrame:
    records = pd.read_csv(pair_path)
    if activity_path is None:
        records["pchembl1"] = np.nan
        records["pchembl2"] = np.nan
        records["signed_delta"] = np.nan
        records["direction_label"] = np.nan
        return records
    lookup = make_activity_lookup(source, activity_path)
    key1 = pd.MultiIndex.from_arrays([records["target"].astype(str), records["smiles1"].astype(str)])
    key2 = pd.MultiIndex.from_arrays([records["target"].astype(str), records["smiles2"].astype(str)])
    records["pchembl1"] = lookup.reindex(key1).to_numpy()
    records["pchembl2"] = lookup.reindex(key2).to_numpy()
    records["signed_delta"] = records["pchembl2"] - records["pchembl1"]
    records["direction_label"] = (records["signed_delta"] > 0).astype(float)
    return records


def build_directed_feature_matrices(records: pd.DataFrame) -> tuple[sparse.csr_matrix, np.ndarray]:
    desc_cache: dict[str, np.ndarray] = {}

    def desc(smiles: str) -> np.ndarray:
        if smiles not in desc_cache:
            desc_cache[smiles] = raw.molecular_descriptors(smiles)
        return desc_cache[smiles]

    indices: list[int] = []
    indptr = [0]
    data: list[float] = []
    scalar_rows = []
    for row in records.itertuples(index=False):
        fa = fp_array(str(row.smiles1))
        fb = fp_array(str(row.smiles2))
        diff = fb - fa
        nz = np.flatnonzero(diff)
        indices.extend(nz.tolist())
        data.extend(diff[nz].astype(float).tolist())
        indptr.append(len(indices))
        common = float(np.count_nonzero((fa > 0) & (fb > 0)))
        union = float(np.count_nonzero((fa > 0) | (fb > 0)))
        tanimoto = common / union if union > 0 else 0.0
        pop_a = float(np.count_nonzero(fa))
        pop_b = float(np.count_nonzero(fb))
        da = desc(str(row.smiles1))
        db = desc(str(row.smiles2))
        scalar_rows.append(
            np.concatenate(
                [
                    np.asarray([tanimoto, pop_a, pop_b, pop_b - pop_a, common, union - common], dtype=np.float32),
                    db - da,
                    da + db,
                ]
            )
        )
    signed_bits = sparse.csr_matrix(
        (np.asarray(data, dtype=np.float32), np.asarray(indices, dtype=np.int32), np.asarray(indptr, dtype=np.int64)),
        shape=(len(records), 2048),
        dtype=np.float32,
    )
    return signed_bits, np.vstack(scalar_rows).astype(np.float32)


def prepare_features(records: pd.DataFrame, with_ligand_similarity: bool = False) -> None:
    global SYMMETRIC_FEATURES, DIRECTED_FEATURES, SCALARS, DIRECTED_SCALARS
    global LIGAND_INDEX, LIGAND_FP_MATRIX, LIGAND_FP_POPCOUNT, LIGAND_SIM_MATRIX
    SYMMETRIC_FEATURES = {}
    DIRECTED_FEATURES = {}
    LIGAND_INDEX = {}
    LIGAND_FP_MATRIX = None
    LIGAND_FP_POPCOUNT = None
    LIGAND_SIM_MATRIX = None
    diff_bits, scalars = raw.build_feature_matrices(records)
    signed_bits, directed_scalars = build_directed_feature_matrices(records)
    SYMMETRIC_FEATURES["ecfp_absdiff_xgb"] = diff_bits
    DIRECTED_FEATURES["directed_ecfp_xgb"] = signed_bits
    SCALARS = scalars
    DIRECTED_SCALARS = directed_scalars
    if not with_ligand_similarity:
        return

    ligands = pd.unique(pd.concat([records["smiles1"].astype(str), records["smiles2"].astype(str)], ignore_index=True))
    LIGAND_INDEX = {smiles: idx for idx, smiles in enumerate(ligands)}
    indices: list[int] = []
    indptr = [0]
    data: list[float] = []
    popcounts = np.zeros(len(ligands), dtype=np.float32)
    for smiles in ligands:
        bits = np.flatnonzero(fp_array(str(smiles)) > 0)
        indices.extend(bits.astype(np.int32).tolist())
        data.extend([1.0] * len(bits))
        popcounts[LIGAND_INDEX[str(smiles)]] = float(len(bits))
        indptr.append(len(indices))
    LIGAND_FP_MATRIX = sparse.csr_matrix(
        (np.asarray(data, dtype=np.float32), np.asarray(indices, dtype=np.int32), np.asarray(indptr, dtype=np.int64)),
        shape=(len(ligands), 2048),
        dtype=np.float32,
    )
    LIGAND_FP_POPCOUNT = popcounts
    sim = LIGAND_FP_MATRIX.dot(LIGAND_FP_MATRIX.T).tocsr().astype(np.float32)
    row_ids = np.repeat(np.arange(sim.shape[0], dtype=np.int32), np.diff(sim.indptr))
    denom = LIGAND_FP_POPCOUNT[row_ids] + LIGAND_FP_POPCOUNT[sim.indices] - sim.data
    valid = denom > 0
    sim.data[valid] = sim.data[valid] / denom[valid]
    sim.data[~valid] = 0.0
    sim.eliminate_zeros()
    LIGAND_SIM_MATRIX = sim


def random_masks(records: pd.DataFrame, seed: int, folds: int) -> list[np.ndarray]:
    y = records["label"].to_numpy(dtype=np.int32)
    splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    masks = []
    for _, test_idx in splitter.split(np.zeros(len(records)), y):
        mask = np.zeros(len(records), dtype=bool)
        mask[test_idx] = True
        masks.append(mask)
    return masks


def group_hash_masks(records: pd.DataFrame, group_col: str, seed: int, folds: int) -> list[np.ndarray]:
    groups = records[group_col].fillna("MISSING").astype(str)
    fold_id = groups.map(lambda value: stable_int(value, seed) % folds).to_numpy()
    return [fold_id == fold for fold in range(folds)]


def split_masks(records: pd.DataFrame, split_mode: str, seed: int, folds: int, clusters: int) -> list[tuple[np.ndarray, np.ndarray, int]]:
    out = []
    if split_mode == "random":
        masks = random_masks(records, seed, folds)
        for mask in masks:
            out.append((~mask, mask, -1))
    elif split_mode in {"target_family", "family_scaffold_source_purged"}:
        masks, _ = raw.family_folds(records, seed, folds, 128, clusters)
        for mask in masks:
            base_train = ~mask
            train, _ = raw.purge_train_mask(records, base_train, mask, split_mode)
            out.append((train, mask, -1))
    elif split_mode == "pair_scaffold_group":
        masks = group_hash_masks(records, "pair_scaffold_key", seed, folds)
        for mask in masks:
            out.append((~mask, mask, -1))
    elif split_mode == "document_source_group":
        masks = group_hash_masks(records, "document_source_key", seed, folds)
        for mask in masks:
            out.append((~mask, mask, -1))
    elif split_mode == "ligand_cold":
        for mask in random_masks(records, seed, folds):
            base_train = ~mask
            test_ligands = set(records.loc[mask, "smiles1"].astype(str)) | set(records.loc[mask, "smiles2"].astype(str))
            keep = ~(
                records["smiles1"].astype(str).isin(test_ligands).to_numpy()
                | records["smiles2"].astype(str).isin(test_ligands).to_numpy()
            )
            out.append((base_train & keep, mask, -1))
    elif split_mode == "temporal_forward":
        for fold in range(folds):
            try:
                train, test, cutoff = raw.temporal_masks(records, fold, folds)
            except Exception:
                continue
            out.append((train, test, cutoff))
    else:
        raise ValueError(split_mode)
    return out


def build_symmetric_x(variant: str, train_mask: np.ndarray):
    assert SCALARS is not None
    if variant == "ecfp_absdiff_xgb":
        return SYMMETRIC_FEATURES["ecfp_absdiff_xgb"]
    scaler = StandardScaler().fit(SCALARS[train_mask])
    scalars = scaler.transform(SCALARS).astype(np.float32)
    if variant in {"pccfrl_scalar_xgb", "rf_scalars", "extratrees_scalars"}:
        return scalars
    if variant == "pccfrl_xgb":
        return sparse.hstack([SYMMETRIC_FEATURES["ecfp_absdiff_xgb"], sparse.csr_matrix(scalars)], format="csr")
    raise ValueError(variant)


def build_directed_x(variant: str, train_mask: np.ndarray):
    assert DIRECTED_SCALARS is not None
    if variant == "directed_ecfp_xgb":
        return DIRECTED_FEATURES["directed_ecfp_xgb"]
    scaler = StandardScaler().fit(DIRECTED_SCALARS[train_mask])
    scalars = scaler.transform(DIRECTED_SCALARS).astype(np.float32)
    if variant in {"directed_scalars_xgb", "directed_rf_scalars", "directed_extratrees_scalars"}:
        return scalars
    if variant == "directed_pccfrl_xgb":
        return sparse.hstack([DIRECTED_FEATURES["directed_ecfp_xgb"], sparse.csr_matrix(scalars)], format="csr")
    raise ValueError(variant)


def classifier(variant: str, seed: int, y_train: np.ndarray):
    if variant in {"rf_scalars", "directed_rf_scalars"}:
        return RandomForestClassifier(n_estimators=300, max_features="sqrt", min_samples_leaf=2, class_weight="balanced", n_jobs=1, random_state=seed)
    if variant in {"extratrees_scalars", "directed_extratrees_scalars"}:
        return ExtraTreesClassifier(n_estimators=300, max_features="sqrt", min_samples_leaf=2, class_weight="balanced", n_jobs=1, random_state=seed)
    positives = float(y_train.sum())
    negatives = float(len(y_train) - positives)
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
        scale_pos_weight=(negatives / positives) if positives > 0 else 1.0,
    )


def regressor(variant: str, seed: int):
    if variant == "rf_scalars":
        return RandomForestRegressor(n_estimators=300, max_features="sqrt", min_samples_leaf=2, n_jobs=1, random_state=seed)
    if variant == "extratrees_scalars":
        return ExtraTreesRegressor(n_estimators=300, max_features="sqrt", min_samples_leaf=2, n_jobs=1, random_state=seed)
    return XGBRegressor(
        n_estimators=140,
        max_depth=4,
        learning_rate=0.04,
        subsample=0.85,
        colsample_bytree=0.8,
        reg_lambda=4.0,
        objective="reg:squarederror",
        tree_method="hist",
        device="cpu",
        n_jobs=1,
        random_state=seed,
    )


def binary_metrics(y: np.ndarray, prob: np.ndarray) -> dict[str, float]:
    pred = (prob >= 0.5).astype(np.int32)
    return {
        "roc_auc": float(roc_auc_score(y, prob)) if len(np.unique(y)) == 2 else float("nan"),
        "pr_auc": float(average_precision_score(y, prob)) if len(np.unique(y)) == 2 else float("nan"),
        "accuracy": float(accuracy_score(y, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)) if len(np.unique(y)) == 2 else float("nan"),
        "mcc": float(matthews_corrcoef(y, pred)) if len(np.unique(pred)) > 1 and len(np.unique(y)) > 1 else 0.0,
        "brier": float(brier_score_loss(y, prob)),
        "ece_10": expected_calibration_error(y, prob, 10),
    }


def regression_metrics(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    return {
        "mae": float(mean_absolute_error(y, pred)),
        "rmse": float(math.sqrt(mean_squared_error(y, pred))),
        "r2": float(r2_score(y, pred)) if len(y) > 1 else float("nan"),
        "pearson": float(pearsonr(y, pred).statistic) if len(np.unique(y)) > 1 and len(np.unique(pred)) > 1 else float("nan"),
        "spearman": float(spearmanr(y, pred).statistic) if len(np.unique(y)) > 1 and len(np.unique(pred)) > 1 else float("nan"),
    }


def expected_calibration_error(y: np.ndarray, prob: np.ndarray, bins: int = 10) -> float:
    edges = np.linspace(0, 1, bins + 1)
    ece = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (prob >= lo) & (prob <= hi) if hi == 1 else (prob >= lo) & (prob < hi)
        if mask.any():
            ece += float(mask.mean()) * abs(float(prob[mask].mean()) - float(y[mask].mean()))
    return float(ece)


def source_tokens(value: str) -> set[str]:
    return raw.source_tokens(str(value))


def corrupt_sources(records: pd.DataFrame, rate: float, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    out = records.copy()
    for col in ["document_source_key", "assay_source_key"]:
        if col not in out.columns:
            out[col] = ""
        out[col] = out[col].fillna("").astype(str)
    if rate <= 0:
        return out
    mask = rng.random(len(out)) < rate
    out.loc[mask, "document_source_key"] = ""
    out.loc[mask, "assay_source_key"] = ""
    return out


def run_multitask(args: argparse.Namespace, source: str, pair_path: str, activity_path: str | None) -> tuple[pd.DataFrame, pd.DataFrame]:
    global RECORDS
    RECORDS = augment_pairs(source, pair_path, activity_path).dropna(subset=["activity_delta"]).reset_index(drop=True)
    prepare_features(RECORDS)
    regression_variants = parse_csv(args.regression_variants)
    direction_variants = parse_csv(args.direction_variants)
    rows = []
    target_rows = []

    tasks = []
    for split_mode in parse_csv(args.multitask_split_modes):
        for seed in [int(x) for x in parse_csv(args.multitask_seeds)]:
            for fold, (train_mask, test_mask, cutoff) in enumerate(split_masks(RECORDS, split_mode, seed, args.folds, args.target_family_clusters)):
                tasks.append((split_mode, seed, fold, cutoff, train_mask, test_mask))

    def one_task(task):
        split_mode, seed, fold, cutoff, train_mask, test_mask = task
        local_rows = []
        local_targets = []
        y_train_reg = RECORDS.loc[train_mask, "activity_delta"].to_numpy(dtype=np.float32)
        y_test_reg = RECORDS.loc[test_mask, "activity_delta"].to_numpy(dtype=np.float32)
        if train_mask.sum() < args.min_train_rows or test_mask.sum() < args.min_test_rows:
            return local_rows, local_targets
        for variant in regression_variants:
            x = build_symmetric_x(variant, train_mask)
            model = regressor(variant, seed)
            model.fit(x[train_mask], y_train_reg)
            pred = model.predict(x[test_mask]).astype(np.float32)
            row = {
                "experiment": "regression_delta",
                "benchmark": source,
                "split_mode": split_mode,
                "variant": variant,
                "seed": seed,
                "fold": fold,
                "temporal_cutoff_year": cutoff,
                "n_train": int(train_mask.sum()),
                "n_test": int(test_mask.sum()),
                **raw.audit_overlaps(RECORDS, train_mask, test_mask),
                **regression_metrics(y_test_reg, pred),
            }
            local_rows.append(row)
        direction_records = RECORDS[np.abs(RECORDS["signed_delta"].fillna(0).to_numpy(dtype=np.float32)) >= args.direction_min_abs_delta]
        if len(direction_records) < args.min_test_rows:
            return local_rows, local_targets
        dir_idx = direction_records.index.to_numpy()
        dir_train = train_mask[dir_idx]
        dir_test = test_mask[dir_idx]
        if dir_train.sum() < args.min_train_rows or dir_test.sum() < args.min_test_rows:
            return local_rows, local_targets
        y_train_dir = direction_records.loc[dir_train, "direction_label"].to_numpy(dtype=np.int32)
        y_test_dir = direction_records.loc[dir_test, "direction_label"].to_numpy(dtype=np.int32)
        if len(np.unique(y_train_dir)) < 2 or len(np.unique(y_test_dir)) < 2:
            return local_rows, local_targets
        full_dir_train = np.zeros(len(RECORDS), dtype=bool)
        full_dir_test = np.zeros(len(RECORDS), dtype=bool)
        full_dir_train[dir_idx[dir_train]] = True
        full_dir_test[dir_idx[dir_test]] = True
        for variant in direction_variants:
            x = build_directed_x(variant, full_dir_train)
            model = classifier(variant, seed, y_train_dir)
            model.fit(x[full_dir_train], y_train_dir)
            prob = model.predict_proba(x[full_dir_test])[:, 1].astype(np.float32)
            row = {
                "experiment": "direction_sign",
                "benchmark": source,
                "split_mode": split_mode,
                "variant": variant,
                "seed": seed,
                "fold": fold,
                "temporal_cutoff_year": cutoff,
                "n_train": int(full_dir_train.sum()),
                "n_test": int(full_dir_test.sum()),
                **raw.audit_overlaps(RECORDS, full_dir_train, full_dir_test),
                **binary_metrics(y_test_dir, prob),
            }
            local_rows.append(row)
        return local_rows, local_targets

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(one_task, task) for task in tasks]
        for future in as_completed(futures):
            got_rows, got_targets = future.result()
            rows.extend(got_rows)
            target_rows.extend(got_targets)
            for row in got_rows:
                key = "roc_auc" if row["experiment"] == "direction_sign" else "mae"
                print(f"[5 {row['benchmark']} {row['experiment']} {row['split_mode']} {row['variant']} seed={row['seed']} fold={row['fold']}] {key}={row[key]:.4f}", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(target_rows)


def nearest_ligand_bins(records: pd.DataFrame, train_mask: np.ndarray, test_mask: np.ndarray, sample_size: int, seed: int) -> pd.Series:
    train_smiles = sorted(set(records.loc[train_mask, "smiles1"].astype(str)) | set(records.loc[train_mask, "smiles2"].astype(str)))
    if sample_size > 0 and len(train_smiles) > sample_size:
        train_smiles = sorted(train_smiles, key=lambda smiles: stable_int(smiles, seed))[:sample_size]
    train_fps = [fp_obj(s) for s in train_smiles]
    train_fps = [fp for fp in train_fps if fp is not None]
    test_frame = records.loc[test_mask, ["smiles1", "smiles2"]].astype(str)
    test_smiles = pd.unique(pd.concat([test_frame["smiles1"], test_frame["smiles2"]], ignore_index=True))
    nearest_map: dict[str, float] = {}

    for smiles in test_smiles:
        fp = fp_obj(str(smiles))
        if fp is None or not train_fps:
            nearest_map[str(smiles)] = 0.0
        else:
            nearest_map[str(smiles)] = float(max(DataStructs.BulkTanimotoSimilarity(fp, train_fps)))

    vals = []
    for row in test_frame.itertuples(index=False):
        vals.append(min(nearest_map.get(str(row.smiles1), 0.0), nearest_map.get(str(row.smiles2), 0.0)))
    nn = pd.Series(vals, index=records.index[test_mask], name="ligand_nn_min_tanimoto")
    return nn


def severity_bin(value: float) -> str:
    if value >= 0.8:
        return "seen_like_ge_0p8"
    if value >= 0.6:
        return "mild_0p6_0p8"
    if value >= 0.4:
        return "moderate_0p4_0p6"
    return "severe_lt_0p4"


def run_shift(args: argparse.Namespace, source: str, pair_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    global RECORDS
    RECORDS = pd.read_csv(pair_path).reset_index(drop=True)
    prepare_features(RECORDS)
    variants = parse_csv(args.shift_variants)
    fold_rows = []
    severity_rows = []
    tasks = []
    for split_mode in parse_csv(args.shift_split_modes):
        for seed in [int(x) for x in parse_csv(args.shift_seeds)]:
            for fold, (train_mask, test_mask, cutoff) in enumerate(split_masks(RECORDS, split_mode, seed, args.folds, args.target_family_clusters)):
                tasks.append((split_mode, seed, fold, cutoff, train_mask, test_mask))

    def one_task(task):
        split_mode, seed, fold, cutoff, train_mask, test_mask = task
        local_folds = []
        local_bins = []
        if train_mask.sum() < args.min_train_rows or test_mask.sum() < args.min_test_rows:
            return local_folds, local_bins
        y_train = RECORDS.loc[train_mask, "label"].to_numpy(dtype=np.int32)
        y_test = RECORDS.loc[test_mask, "label"].to_numpy(dtype=np.int32)
        if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
            return local_folds, local_bins
        nn = nearest_ligand_bins(RECORDS, train_mask, test_mask, args.nearest_train_sample_size, seed + 1009 * (fold + 1))
        bins = nn.map(severity_bin)
        for variant in variants:
            x = build_symmetric_x(variant, train_mask)
            model = classifier(variant, seed, y_train)
            model.fit(x[train_mask], y_train)
            prob = model.predict_proba(x[test_mask])[:, 1].astype(np.float32)
            row = {
                "benchmark": source,
                "split_mode": split_mode,
                "variant": variant,
                "seed": seed,
                "fold": fold,
                "temporal_cutoff_year": cutoff,
                "nearest_train_sample_size": int(args.nearest_train_sample_size),
                "n_train": int(train_mask.sum()),
                "n_test": int(test_mask.sum()),
                **raw.audit_overlaps(RECORDS, train_mask, test_mask),
                **binary_metrics(y_test, prob),
            }
            local_folds.append(row)
            test_index = np.flatnonzero(test_mask)
            for bin_name in ["seen_like_ge_0p8", "mild_0p6_0p8", "moderate_0p4_0p6", "severe_lt_0p4"]:
                local = bins.to_numpy() == bin_name
                if local.sum() < args.min_bin_rows or len(np.unique(y_test[local])) < 2:
                    continue
                local_bins.append(
                    {
                        "benchmark": source,
                        "split_mode": split_mode,
                        "variant": variant,
                        "seed": seed,
                        "fold": fold,
                        "severity_bin": bin_name,
                        "nearest_train_sample_size": int(args.nearest_train_sample_size),
                        "ligand_nn_min_tanimoto_mean": float(nn.to_numpy()[local].mean()),
                        "n_rows": int(local.sum()),
                        **binary_metrics(y_test[local], prob[local]),
                    }
                )
        return local_folds, local_bins

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(one_task, task) for task in tasks]
        for future in as_completed(futures):
            got_folds, got_bins = future.result()
            fold_rows.extend(got_folds)
            severity_rows.extend(got_bins)
            for row in got_folds:
                print(f"[6 {row['benchmark']} {row['split_mode']} {row['variant']} seed={row['seed']} fold={row['fold']}] roc={row['roc_auc']:.4f}", flush=True)
    return pd.DataFrame(fold_rows), pd.DataFrame(severity_rows)


def selective_rows(y: np.ndarray, prob: np.ndarray, coverages: list[float]) -> list[dict]:
    pred = (prob >= 0.5).astype(np.int32)
    confidence = np.maximum(prob, 1 - prob)
    order = np.argsort(-confidence)
    rows = []
    for coverage in coverages:
        n = max(1, int(math.ceil(len(y) * coverage)))
        idx = order[:n]
        rows.append(
            {
                "coverage": float(coverage),
                "n_selected": int(n),
                "selective_error": float((pred[idx] != y[idx]).mean()),
                "selective_accuracy": float((pred[idx] == y[idx]).mean()),
                "mean_confidence": float(confidence[idx].mean()),
            }
        )
    return rows


def run_reliability(args: argparse.Namespace, source: str, pair_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    global RECORDS, ORIGINAL_RECORDS
    original = pd.read_csv(pair_path).reset_index(drop=True)
    ORIGINAL_RECORDS = original
    RECORDS = original
    prepare_features(original)
    variants = parse_csv(args.reliability_variants)
    fold_rows = []
    rc_rows = []
    tasks = []
    for missing_rate in [float(x) for x in parse_csv(args.missing_source_rates)]:
        corrupted = corrupt_sources(original, missing_rate, args.missing_seed)
        for seed in [int(x) for x in parse_csv(args.reliability_seeds)]:
            for fold, (train_mask, test_mask, cutoff) in enumerate(split_masks(corrupted, "family_scaffold_source_purged", seed, args.folds, args.target_family_clusters)):
                tasks.append((missing_rate, seed, fold, cutoff, corrupted, train_mask, test_mask))

    def one_task(task):
        missing_rate, seed, fold, cutoff, corrupted, train_mask, test_mask = task
        local_rows = []
        local_rc = []
        if train_mask.sum() < args.min_train_rows or test_mask.sum() < args.min_test_rows:
            return local_rows, local_rc
        y_train = original.loc[train_mask, "label"].to_numpy(dtype=np.int32)
        y_test = original.loc[test_mask, "label"].to_numpy(dtype=np.int32)
        if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
            return local_rows, local_rc
        observed_audit = raw.audit_overlaps(corrupted, train_mask, test_mask)
        true_audit = raw.audit_overlaps(original, train_mask, test_mask)
        for variant in variants:
            x = build_symmetric_x(variant, train_mask)
            model = classifier(variant, seed, y_train)
            model.fit(x[train_mask], y_train)
            prob = model.predict_proba(x[test_mask])[:, 1].astype(np.float32)
            row = {
                "benchmark": source,
                "missing_source_rate": missing_rate,
                "split_mode": "family_scaffold_source_purged_observed_missing",
                "variant": variant,
                "seed": seed,
                "fold": fold,
                "n_train": int(train_mask.sum()),
                "n_test": int(test_mask.sum()),
                "observed_document_source_overlap_rate": observed_audit["document_source_overlap_rate"],
                "observed_assay_source_overlap_rate": observed_audit["assay_source_overlap_rate"],
                "true_document_source_overlap_rate": true_audit["document_source_overlap_rate"],
                "true_assay_source_overlap_rate": true_audit["assay_source_overlap_rate"],
                "hidden_document_source_overlap_rate": true_audit["document_source_overlap_rate"] - observed_audit["document_source_overlap_rate"],
                "hidden_assay_source_overlap_rate": true_audit["assay_source_overlap_rate"] - observed_audit["assay_source_overlap_rate"],
                **{f"true_{k}": v for k, v in true_audit.items() if k not in {"document_source_overlap_rate", "assay_source_overlap_rate"}},
                **binary_metrics(y_test, prob),
            }
            local_rows.append(row)
            for rc in selective_rows(y_test, prob, [float(x) for x in parse_csv(args.coverages)]):
                rc.update(
                    {
                        "benchmark": source,
                        "missing_source_rate": missing_rate,
                        "variant": variant,
                        "seed": seed,
                        "fold": fold,
                    }
                )
                local_rc.append(rc)
        return local_rows, local_rc

    # Source corruption changes only the observed split/purge metadata, not chemistry.
    # Keep feature matrices fixed to the original row order so threaded tasks cannot
    # overwrite shared feature state while estimating missing-source sensitivity.
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(one_task, task) for task in tasks]
        for future in as_completed(futures):
            got_rows, got_rc = future.result()
            fold_rows.extend(got_rows)
            rc_rows.extend(got_rc)
            for row in got_rows:
                print(f"[7 {row['benchmark']} miss={row['missing_source_rate']} {row['variant']} seed={row['seed']} fold={row['fold']}] roc={row['roc_auc']:.4f} ece={row['ece_10']:.4f}", flush=True)
    return pd.DataFrame(fold_rows), pd.DataFrame(rc_rows)


def write_summary(prefix: Path, name: str, df: pd.DataFrame, group_cols: list[str], metric_cols: list[str]) -> None:
    if df.empty:
        return
    summary = df.groupby(group_cols, as_index=False)[metric_cols].agg(["mean", "std", "count"])
    summary.columns = ["_".join([part for part in col if part]) for col in summary.columns.to_flat_index()]
    summary.reset_index(drop=True).to_csv(f"{prefix}_{name}_summary.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sources",
        default=(
            "chembl_raw=data/chembl_raw_moleculeace_target_pairs.csv=data/chembl_raw_moleculeace_target_activities.csv,"
            "bindingdb_raw=data/bindingdb_raw_moleculeace_target_pairs.csv=data/bindingdb_raw_moleculeace_target_activities.csv"
        ),
    )
    parser.add_argument("--experiments", default="5,6,7")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--target-family-clusters", type=int, default=10)
    parser.add_argument("--workers", type=int, default=24)
    parser.add_argument("--min-train-rows", type=int, default=300)
    parser.add_argument("--min-test-rows", type=int, default=100)
    parser.add_argument("--min-bin-rows", type=int, default=100)
    parser.add_argument("--output-prefix", default="results/tkde_experiments_5_7")
    parser.add_argument("--regression-variants", default="ecfp_absdiff_xgb,pccfrl_xgb,pccfrl_scalar_xgb,rf_scalars,extratrees_scalars")
    parser.add_argument("--direction-variants", default="directed_ecfp_xgb,directed_pccfrl_xgb,directed_scalars_xgb,directed_rf_scalars,directed_extratrees_scalars")
    parser.add_argument("--multitask-split-modes", default="random,target_family,family_scaffold_source_purged")
    parser.add_argument("--multitask-seeds", default="0,1,2,3,4")
    parser.add_argument("--direction-min-abs-delta", type=float, default=0.3)
    parser.add_argument("--shift-variants", default="ecfp_absdiff_xgb,pccfrl_scalar_xgb,rf_scalars")
    parser.add_argument("--shift-split-modes", default="random,target_family,family_scaffold_source_purged,pair_scaffold_group,document_source_group,ligand_cold,temporal_forward")
    parser.add_argument("--shift-seeds", default="0,1,2")
    parser.add_argument("--nearest-train-sample-size", type=int, default=4096)
    parser.add_argument("--reliability-variants", default="ecfp_absdiff_xgb,pccfrl_scalar_xgb,rf_scalars")
    parser.add_argument("--reliability-seeds", default="0,1,2,3,4")
    parser.add_argument("--missing-source-rates", default="0,0.25,0.5,0.75,1.0")
    parser.add_argument("--missing-seed", type=int, default=20260611)
    parser.add_argument("--coverages", default="0.5,0.7,0.9,1.0")
    args = parser.parse_args()

    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    experiments = set(parse_csv(args.experiments))
    all_manifest = {}
    for source, pair_path, activity_path in parse_sources(args.sources):
        print(json.dumps({"source": source, "pair_path": pair_path, "activity_path": activity_path}), flush=True)
        if "5" in experiments:
            rows, _ = run_multitask(args, source, pair_path, activity_path)
            rows.to_csv(f"{prefix}_5_multitask_{source}.csv", index=False)
            write_summary(
                prefix,
                f"5_multitask_{source}",
                rows,
                ["experiment", "benchmark", "split_mode", "variant"],
                ["roc_auc", "pr_auc", "mcc", "mae", "rmse", "r2", "pearson", "spearman", "n_train", "n_test"],
            )
            all_manifest[f"5_{source}_rows"] = int(len(rows))
        if "6" in experiments:
            folds, bins = run_shift(args, source, pair_path)
            folds.to_csv(f"{prefix}_6_shift_folds_{source}.csv", index=False)
            bins.to_csv(f"{prefix}_6_shift_severity_bins_{source}.csv", index=False)
            write_summary(
                prefix,
                f"6_shift_folds_{source}",
                folds,
                ["benchmark", "split_mode", "variant"],
                ["roc_auc", "pr_auc", "mcc", "brier", "ece_10", "n_train", "n_test", "target_overlap_rate", "exact_ligand_overlap_rate", "scaffold_overlap_rate", "document_source_overlap_rate", "assay_source_overlap_rate"],
            )
            write_summary(
                prefix,
                f"6_shift_bins_{source}",
                bins,
                ["benchmark", "split_mode", "variant", "severity_bin"],
                ["roc_auc", "pr_auc", "mcc", "brier", "ece_10", "n_rows", "ligand_nn_min_tanimoto_mean"],
            )
            all_manifest[f"6_{source}_fold_rows"] = int(len(folds))
            all_manifest[f"6_{source}_bin_rows"] = int(len(bins))
        if "7" in experiments:
            folds, rc = run_reliability(args, source, pair_path)
            folds.to_csv(f"{prefix}_7_missing_source_reliability_{source}.csv", index=False)
            rc.to_csv(f"{prefix}_7_risk_coverage_{source}.csv", index=False)
            write_summary(
                prefix,
                f"7_missing_source_reliability_{source}",
                folds,
                ["benchmark", "missing_source_rate", "variant"],
                ["roc_auc", "pr_auc", "mcc", "brier", "ece_10", "observed_document_source_overlap_rate", "true_document_source_overlap_rate", "hidden_document_source_overlap_rate", "observed_assay_source_overlap_rate", "true_assay_source_overlap_rate", "hidden_assay_source_overlap_rate", "n_train", "n_test"],
            )
            write_summary(
                prefix,
                f"7_risk_coverage_{source}",
                rc,
                ["benchmark", "missing_source_rate", "variant", "coverage"],
                ["selective_error", "selective_accuracy", "mean_confidence", "n_selected"],
            )
            all_manifest[f"7_{source}_fold_rows"] = int(len(folds))
            all_manifest[f"7_{source}_rc_rows"] = int(len(rc))
    Path(f"{prefix}_manifest.json").write_text(json.dumps(all_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(all_manifest, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
