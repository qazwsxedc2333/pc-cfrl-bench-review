from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    matthews_corrcoef,
    roc_auc_score,
)
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier


AMINO = "ACDEFGHIKLMNPQRSTVWY"
RECORDS: pd.DataFrame | None = None
DIFF_BITS_X: sparse.csr_matrix | None = None
SCALARS_RAW_X: np.ndarray | None = None
FOLD_CACHE: dict[tuple[int, int, int], tuple[list[np.ndarray], dict[str, int]]] = {}


def parse_csv_arg(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()]


def stable_bucket(text: str, n_hash: int) -> int:
    digest = hashlib.blake2b(text.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "little") % n_hash


def sequence_features(sequence: str, n_hash: int) -> np.ndarray:
    seq = "".join(aa for aa in str(sequence).upper() if aa in AMINO)
    x = np.zeros(20 + n_hash, dtype=np.float32)
    if not seq:
        return x
    for aa in seq:
        x[AMINO.index(aa)] += 1.0
    x[:20] /= max(float(len(seq)), 1.0)
    for k in (2, 3):
        if len(seq) < k:
            continue
        denom = max(float(len(seq) - k + 1), 1.0)
        for idx in range(len(seq) - k + 1):
            token = f"{k}:{seq[idx:idx+k]}"
            x[20 + stable_bucket(token, n_hash)] += 1.0 / denom
    return x


def parse_sequence(fasta: str) -> tuple[str, str]:
    text = str(fasta)
    if "\n" not in text:
        return text, ""
    header, seq = text.split("\n", 1)
    return header.strip(), "".join(seq.split())


def canonical_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return str(smiles)
    return Chem.MolToSmiles(mol, canonical=True)


def scaffold_key(smiles: str) -> str:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return f"BAD:{smiles}"
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    if scaffold:
        return scaffold
    return f"NOSCAFFOLD:{Chem.MolToSmiles(mol, canonical=True)}"


def make_pair_key(a: str, b: str) -> str:
    ca = canonical_smiles(a)
    cb = canonical_smiles(b)
    return "||".join(sorted([ca, cb]))


def molecular_descriptors(smiles: str) -> np.ndarray:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return np.zeros(10, dtype=np.float32)
    return np.asarray(
        [
            Descriptors.MolWt(mol),
            Crippen.MolLogP(mol),
            rdMolDescriptors.CalcTPSA(mol),
            Lipinski.NumHDonors(mol),
            Lipinski.NumHAcceptors(mol),
            Lipinski.NumRotatableBonds(mol),
            rdMolDescriptors.CalcNumRings(mol),
            rdMolDescriptors.CalcNumAromaticRings(mol),
            rdMolDescriptors.CalcFractionCSP3(mol),
            mol.GetNumHeavyAtoms(),
        ],
        dtype=np.float32,
    )


def load_moleculeace_pairs(
    input_csv: str,
    tanimoto_threshold: float,
    cliff_delta: float,
    smooth_delta: float,
    neg_pos_ratio: float,
    downsample_seed: int,
) -> pd.DataFrame:
    df = pd.read_csv(input_csv)
    required = {"Dataset", "SMILES", "y [pEC50/pKi]", "fasta"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing columns in {input_csv}: {missing}")

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp_cache: dict[str, object] = {}
    scaffold_cache: dict[str, str] = {}
    desc_cache: dict[str, np.ndarray] = {}
    rows: list[dict] = []
    rng = np.random.default_rng(downsample_seed)

    def fp(smiles: str):
        if smiles not in fp_cache:
            mol = Chem.MolFromSmiles(str(smiles))
            fp_cache[smiles] = gen.GetFingerprint(mol) if mol is not None else None
        return fp_cache[smiles]

    def scaffold(smiles: str) -> str:
        if smiles not in scaffold_cache:
            scaffold_cache[smiles] = scaffold_key(smiles)
        return scaffold_cache[smiles]

    def desc(smiles: str) -> np.ndarray:
        if smiles not in desc_cache:
            desc_cache[smiles] = molecular_descriptors(smiles)
        return desc_cache[smiles]

    for dataset, group in df.groupby("Dataset", sort=True):
        local = group.reset_index(drop=True).copy()
        smiles = local["SMILES"].astype(str).tolist()
        y = local["y [pEC50/pKi]"].to_numpy(dtype=np.float32)
        fps = [fp(smi) for smi in smiles]
        header, seq = parse_sequence(str(local["fasta"].iloc[0]))
        pos_rows: list[dict] = []
        neg_rows: list[dict] = []
        for i in range(len(local) - 1):
            fp_i = fps[i]
            if fp_i is None:
                continue
            sims = DataStructs.BulkTanimotoSimilarity(fp_i, fps[i + 1 :])
            for j, sim in enumerate(sims, start=i + 1):
                if sim < tanimoto_threshold:
                    continue
                delta = abs(float(y[i] - y[j]))
                if delta >= cliff_delta:
                    label = 1
                elif delta <= smooth_delta:
                    label = 0
                else:
                    continue
                a = smiles[i]
                b = smiles[j]
                sa = scaffold(a)
                sb = scaffold(b)
                fa = np.asarray(fps[i], dtype=object)
                fb = np.asarray(fps[j], dtype=object)
                # Keep descriptor calls here so invalid molecules are caught during construction.
                _ = desc(a)
                _ = desc(b)
                row = {
                    "target": str(dataset),
                    "smiles1": a,
                    "smiles2": b,
                    "label": int(label),
                    "activity_delta": float(delta),
                    "tanimoto": float(sim),
                    "target_header": header,
                    "target_sequence": seq,
                    "scaffold1": sa,
                    "scaffold2": sb,
                    "pair_key": make_pair_key(a, b),
                    "pair_scaffold_key": "||".join(sorted([sa, sb])),
                }
                # Prevent linters from treating fp arrays above as unused in older RDKit wrappers.
                if fa is None or fb is None:
                    continue
                (pos_rows if label else neg_rows).append(row)
        if neg_pos_ratio > 0 and pos_rows:
            max_neg = int(np.ceil(len(pos_rows) * neg_pos_ratio))
            if len(neg_rows) > max_neg:
                keep = rng.choice(len(neg_rows), size=max_neg, replace=False)
                neg_rows = [neg_rows[int(i)] for i in keep]
        rows.extend(pos_rows)
        rows.extend(neg_rows)

    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("No MoleculeACE pair rows constructed")
    out = out.sample(frac=1.0, random_state=downsample_seed).reset_index(drop=True)
    return out


def build_feature_matrices(records: pd.DataFrame) -> tuple[sparse.csr_matrix, np.ndarray]:
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp_cache: dict[str, np.ndarray] = {}
    desc_cache: dict[str, np.ndarray] = {}

    def fp_np(smiles: str) -> np.ndarray:
        if smiles not in fp_cache:
            arr = np.zeros((2048,), dtype=np.float32)
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is not None:
                bv = gen.GetFingerprint(mol)
                DataStructs.ConvertToNumpyArray(bv, arr)
            fp_cache[smiles] = arr
        return fp_cache[smiles]

    def desc(smiles: str) -> np.ndarray:
        if smiles not in desc_cache:
            desc_cache[smiles] = molecular_descriptors(smiles)
        return desc_cache[smiles]

    indices: list[int] = []
    indptr = [0]
    data: list[float] = []
    scalars: list[np.ndarray] = []
    for row in records.itertuples(index=False):
        fa = fp_np(str(row.smiles1))
        fb = fp_np(str(row.smiles2))
        diff_idx = np.flatnonzero(fa != fb)
        common_idx = np.flatnonzero((fa > 0) & (fb > 0))
        indices.extend(diff_idx.tolist())
        data.extend([1.0] * len(diff_idx))
        indptr.append(len(indices))

        union_count = float(len(diff_idx) + len(common_idx))
        common_count = float(len(common_idx))
        tanimoto = common_count / union_count if union_count > 0 else 0.0
        pop_a = float(np.count_nonzero(fa))
        pop_b = float(np.count_nonzero(fb))
        da = desc(str(row.smiles1))
        db = desc(str(row.smiles2))
        scalars.append(
            np.concatenate(
                [
                    np.asarray(
                        [tanimoto, pop_a, pop_b, abs(pop_a - pop_b), common_count, float(len(diff_idx))],
                        dtype=np.float32,
                    ),
                    np.abs(da - db),
                    da + db,
                ]
            )
        )
    diff_bits = sparse.csr_matrix(
        (np.asarray(data, dtype=np.float32), np.asarray(indices, dtype=np.int32), np.asarray(indptr, dtype=np.int64)),
        shape=(len(records), 2048),
        dtype=np.float32,
    )
    return diff_bits, np.vstack(scalars).astype(np.float32)


def target_family_folds(records: pd.DataFrame, seed: int, family_seed: int, n_folds: int, n_families: int, n_hash: int):
    key = (seed, family_seed, n_families)
    if key in FOLD_CACHE:
        return FOLD_CACHE[key]
    targets = sorted(records["target"].astype(str).unique())
    seq_by_target = records.groupby("target")["target_sequence"].first().to_dict()
    x = np.vstack([sequence_features(str(seq_by_target.get(target, "")), n_hash) for target in targets]).astype(np.float32)
    x = StandardScaler().fit_transform(x).astype(np.float32)
    n_clusters = max(n_folds, min(n_families, len(targets)))
    labels = KMeans(n_clusters=n_clusters, random_state=family_seed, n_init=20).fit_predict(x)
    target_to_family = {target: int(label) for target, label in zip(targets, labels)}

    stats = records.groupby("target")["label"].agg(["size", "sum"]).reset_index()
    stats["family"] = stats["target"].astype(str).map(target_to_family)
    family_stats = stats.groupby("family")[["size", "sum"]].sum().reset_index()
    family_stats = family_stats.sample(frac=1.0, random_state=seed).sort_values(["sum", "size"], ascending=False)

    folds: list[list[int]] = [[] for _ in range(n_folds)]
    fold_rows = np.zeros(n_folds, dtype=np.int64)
    fold_pos = np.zeros(n_folds, dtype=np.int64)
    for row in family_stats.itertuples(index=False):
        score = fold_rows + 3 * fold_pos
        fold_idx = int(np.argmin(score))
        folds[fold_idx].append(int(row.family))
        fold_rows[fold_idx] += int(row.size)
        fold_pos[fold_idx] += int(row.sum)
    row_family = records["target"].astype(str).map(target_to_family).to_numpy(dtype=np.int32)
    masks = [np.isin(row_family, fold).astype(bool) for fold in folds if fold]
    FOLD_CACHE[key] = (masks, target_to_family)
    return masks, target_to_family


def purge_training_mask(base_train_mask: np.ndarray, test_mask: np.ndarray, split_mode: str) -> tuple[np.ndarray, float]:
    assert RECORDS is not None
    train_mask = base_train_mask.copy()
    base_n = int(base_train_mask.sum())
    if split_mode == "target_family":
        return train_mask, 1.0
    if split_mode == "family_scaffold_purged":
        test_scaffolds = set(RECORDS.loc[test_mask, "scaffold1"].astype(str)) | set(RECORDS.loc[test_mask, "scaffold2"].astype(str))
        keep = ~(
            RECORDS["scaffold1"].astype(str).isin(test_scaffolds).to_numpy()
            | RECORDS["scaffold2"].astype(str).isin(test_scaffolds).to_numpy()
        )
        train_mask &= keep
    elif split_mode == "family_exact_ligand_purged":
        test_ligands = set(RECORDS.loc[test_mask, "smiles1"].astype(str)) | set(RECORDS.loc[test_mask, "smiles2"].astype(str))
        keep = ~(
            RECORDS["smiles1"].astype(str).isin(test_ligands).to_numpy()
            | RECORDS["smiles2"].astype(str).isin(test_ligands).to_numpy()
        )
        train_mask &= keep
    else:
        raise ValueError(split_mode)
    return train_mask, float(train_mask.sum() / max(base_n, 1))


def build_features(variant: str, train_mask: np.ndarray):
    assert DIFF_BITS_X is not None and SCALARS_RAW_X is not None
    if variant == "ecfp_absdiff_xgb" or variant == "rich_diff_bits_only":
        return DIFF_BITS_X
    scaler = StandardScaler().fit(SCALARS_RAW_X[train_mask])
    scalars = sparse.csr_matrix(scaler.transform(SCALARS_RAW_X).astype(np.float32))
    if variant == "rich_scalars_only":
        return scalars
    if variant == "rich_diff_scalars":
        return sparse.hstack([DIFF_BITS_X, scalars], format="csr")
    raise ValueError(variant)


def safe_auc(y: np.ndarray, score: np.ndarray) -> float:
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, score))


def safe_ap(y: np.ndarray, score: np.ndarray) -> float:
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(average_precision_score(y, score))


def ece_score(y: np.ndarray, prob: np.ndarray, n_bins: int = 10) -> float:
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    for left, right in zip(edges[:-1], edges[1:]):
        mask = (prob >= left) & (prob <= right) if right == 1.0 else (prob >= left) & (prob < right)
        if np.any(mask):
            ece += float(mask.mean()) * abs(float(prob[mask].mean()) - float(y[mask].mean()))
    return float(ece)


def row_metrics(y: np.ndarray, prob: np.ndarray) -> dict[str, float]:
    pred = (prob >= 0.5).astype(np.int32)
    return {
        "roc_auc": safe_auc(y, prob),
        "pr_auc": safe_ap(y, prob),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "f1": float(f1_score(y, pred, zero_division=0)),
        "mcc": float(matthews_corrcoef(y, pred)) if len(np.unique(pred)) > 1 else 0.0,
        "brier": float(brier_score_loss(y, prob)) if len(np.unique(y)) > 1 else float("nan"),
        "ece_10": ece_score(y, prob),
        "score_mean": float(np.mean(prob)),
        "score_std": float(np.std(prob)),
    }


def audit_overlaps(records: pd.DataFrame, train_mask: np.ndarray, test_mask: np.ndarray, target_to_family: dict[str, int]) -> dict[str, float]:
    train_targets = set(records.loc[train_mask, "target"].astype(str))
    train_families = {target_to_family[str(t)] for t in train_targets}
    train_ligands = set(records.loc[train_mask, "smiles1"].astype(str)) | set(records.loc[train_mask, "smiles2"].astype(str))
    train_scaffolds = set(records.loc[train_mask, "scaffold1"].astype(str)) | set(records.loc[train_mask, "scaffold2"].astype(str))
    train_pairs = set(records.loc[train_mask, "pair_key"].astype(str))
    train_pair_scaffolds = set(records.loc[train_mask, "pair_scaffold_key"].astype(str))
    test_targets = records.loc[test_mask, "target"].astype(str)
    test_families = test_targets.map(target_to_family)
    exact_ligand_overlap = (
        records.loc[test_mask, "smiles1"].astype(str).isin(train_ligands).to_numpy()
        | records.loc[test_mask, "smiles2"].astype(str).isin(train_ligands).to_numpy()
    )
    scaffold_overlap = (
        records.loc[test_mask, "scaffold1"].astype(str).isin(train_scaffolds).to_numpy()
        | records.loc[test_mask, "scaffold2"].astype(str).isin(train_scaffolds).to_numpy()
    )
    return {
        "target_overlap_rate": float(test_targets.isin(train_targets).mean()),
        "family_overlap_rate": float(test_families.isin(train_families).mean()),
        "exact_ligand_overlap_rate": float(exact_ligand_overlap.mean()),
        "scaffold_overlap_rate": float(scaffold_overlap.mean()),
        "exact_pair_overlap_rate": float(records.loc[test_mask, "pair_key"].astype(str).isin(train_pairs).mean()),
        "pair_scaffold_overlap_rate": float(records.loc[test_mask, "pair_scaffold_key"].astype(str).isin(train_pair_scaffolds).mean()),
        "n_train_targets": int(records.loc[train_mask, "target"].nunique()),
        "n_test_targets": int(records.loc[test_mask, "target"].nunique()),
        "n_train_families": int(pd.Series(records.loc[train_mask, "target"].astype(str).map(target_to_family)).nunique()),
        "n_test_families": int(pd.Series(records.loc[test_mask, "target"].astype(str).map(target_to_family)).nunique()),
    }


def make_model(seed: int, y_train: np.ndarray, threads: int, n_estimators: int, max_depth: int) -> XGBClassifier:
    positives = float(y_train.sum())
    negatives = float(len(y_train) - positives)
    return XGBClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        learning_rate=0.05,
        subsample=0.85,
        colsample_bytree=0.8,
        reg_lambda=4.0,
        objective="binary:logistic",
        eval_metric="aucpr",
        tree_method="hist",
        device="cpu",
        n_jobs=threads,
        random_state=seed,
        scale_pos_weight=(negatives / positives) if positives > 0 else 1.0,
    )


def run_task(task: dict) -> dict:
    assert RECORDS is not None
    seed = int(task["seed"])
    fold = int(task["fold"])
    split_mode = str(task["split_mode"])
    variant = str(task["variant"])
    masks, target_to_family = target_family_folds(
        RECORDS,
        seed=seed,
        family_seed=int(task["family_seed"]),
        n_folds=int(task["folds"]),
        n_families=int(task["n_families"]),
        n_hash=int(task["sequence_hash_features"]),
    )
    test_mask = masks[fold]
    base_train_mask = ~test_mask
    train_mask, retention = purge_training_mask(base_train_mask, test_mask, split_mode)
    y_train = RECORDS.loc[train_mask, "label"].to_numpy(dtype=np.int32)
    y_test = RECORDS.loc[test_mask, "label"].to_numpy(dtype=np.int32)
    row = {
        "benchmark": "MoleculeACE",
        "split_mode": split_mode,
        "variant": variant,
        "seed": seed,
        "fold": fold,
        "n_train": int(train_mask.sum()),
        "n_test": int(test_mask.sum()),
        "train_positive_rate": float(y_train.mean()) if len(y_train) else float("nan"),
        "test_positive_rate": float(y_test.mean()) if len(y_test) else float("nan"),
        "train_retention_after_purge": retention,
        "n_families": int(task["n_families"]),
        "tanimoto_threshold": float(task["tanimoto_threshold"]),
        "cliff_delta": float(task["cliff_delta"]),
        "smooth_delta": float(task["smooth_delta"]),
        "neg_pos_ratio": float(task["neg_pos_ratio"]),
    }
    row.update(audit_overlaps(RECORDS, train_mask, test_mask, target_to_family))
    if len(y_train) < 20 or len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
        row["error"] = "insufficient_class_support"
        return row
    x = build_features(variant, train_mask)
    model = make_model(seed, y_train, int(task["model_threads"]), int(task["n_estimators"]), int(task["max_depth"]))
    model.fit(x[train_mask], y_train)
    prob = model.predict_proba(x[test_mask])[:, 1].astype(np.float32)
    row.update(row_metrics(y_test, prob))
    return row


def read_done_keys(output: Path) -> set[tuple[str, int, int, str]]:
    done: set[tuple[str, int, int, str]] = set()
    if not output.exists():
        return done
    with output.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            done.add((str(row.get("split_mode")), int(row.get("seed")), int(row.get("fold")), str(row.get("variant"))))
    return done


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="external/MTPNet/data/MoleculeACE.csv")
    parser.add_argument("--variants", default="ecfp_absdiff_xgb,rich_diff_bits_only,rich_scalars_only,rich_diff_scalars")
    parser.add_argument("--split-modes", default="target_family,family_scaffold_purged")
    parser.add_argument("--seeds", default="0,1,2,3,4")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--fold-filter", default="")
    parser.add_argument("--family-seed", type=int, default=31337)
    parser.add_argument("--n-families", type=int, default=10)
    parser.add_argument("--sequence-hash-features", type=int, default=64)
    parser.add_argument("--tanimoto-threshold", type=float, default=0.7)
    parser.add_argument("--cliff-delta", type=float, default=1.0)
    parser.add_argument("--smooth-delta", type=float, default=0.3)
    parser.add_argument("--neg-pos-ratio", type=float, default=5.0)
    parser.add_argument("--downsample-seed", type=int, default=2026)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--model-threads", type=int, default=4)
    parser.add_argument("--n-estimators", type=int, default=60)
    parser.add_argument("--max-depth", type=int, default=4)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--output", default="results/moleculeace_target_family_pairs_seed0_4.jsonl")
    args = parser.parse_args()

    started = time.time()
    global RECORDS, DIFF_BITS_X, SCALARS_RAW_X
    RECORDS = load_moleculeace_pairs(
        args.input,
        tanimoto_threshold=args.tanimoto_threshold,
        cliff_delta=args.cliff_delta,
        smooth_delta=args.smooth_delta,
        neg_pos_ratio=args.neg_pos_ratio,
        downsample_seed=args.downsample_seed,
    )
    DIFF_BITS_X, SCALARS_RAW_X = build_feature_matrices(RECORDS)

    variants = parse_csv_arg(args.variants)
    split_modes = parse_csv_arg(args.split_modes)
    seeds = [int(x) for x in parse_csv_arg(args.seeds)]
    fold_filter = {int(x) for x in parse_csv_arg(args.fold_filter)} if args.fold_filter else None
    tasks = []
    for split_mode in split_modes:
        for seed in seeds:
            masks, _ = target_family_folds(
                RECORDS,
                seed=seed,
                family_seed=args.family_seed,
                n_folds=args.folds,
                n_families=args.n_families,
                n_hash=args.sequence_hash_features,
            )
            for fold in range(len(masks)):
                if fold_filter is not None and fold not in fold_filter:
                    continue
                for variant in variants:
                    tasks.append(
                        {
                            "split_mode": split_mode,
                            "seed": seed,
                            "fold": fold,
                            "folds": args.folds,
                            "variant": variant,
                            "family_seed": args.family_seed,
                            "n_families": args.n_families,
                            "sequence_hash_features": args.sequence_hash_features,
                            "tanimoto_threshold": args.tanimoto_threshold,
                            "cliff_delta": args.cliff_delta,
                            "smooth_delta": args.smooth_delta,
                            "neg_pos_ratio": args.neg_pos_ratio,
                            "model_threads": args.model_threads,
                            "n_estimators": args.n_estimators,
                            "max_depth": args.max_depth,
                        }
                    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if args.resume:
        done = read_done_keys(output)
        tasks = [t for t in tasks if (t["split_mode"], t["seed"], t["fold"], t["variant"]) not in done]
    print(
        f"pairs={len(RECORDS)} targets={RECORDS['target'].nunique()} pos_rate={RECORDS['label'].mean():.4f} "
        f"features=diff_bits:{DIFF_BITS_X.shape} scalars:{SCALARS_RAW_X.shape} tasks={len(tasks)} "
        f"elapsed_build={time.time() - started:.1f}s output={output}",
        flush=True,
    )
    with output.open("a", encoding="utf-8") as handle:
        if args.workers <= 1:
            iterator = (run_task(task) for task in tasks)
            pool = None
        else:
            ctx = mp.get_context("fork")
            pool = ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx)
            iterator = (future.result() for future in as_completed([pool.submit(run_task, task) for task in tasks]))
        try:
            for row in iterator:
                handle.write(json.dumps(row, sort_keys=True) + "\n")
                handle.flush()
                if "error" in row:
                    status = f"ERROR={row['error']}"
                else:
                    status = f"roc={row['roc_auc']:.4f} pr={row['pr_auc']:.4f}"
                print(f"[{row['split_mode']} seed={row['seed']} fold={row['fold']} {row['variant']}] {status}", flush=True)
        finally:
            if pool is not None:
                pool.shutdown(wait=True)


if __name__ == "__main__":
    main()
