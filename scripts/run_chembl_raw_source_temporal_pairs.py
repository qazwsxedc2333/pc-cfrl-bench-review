from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy import sparse
from scipy.stats import binomtest, wilcoxon
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
CH_EMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
METRICS = ["roc_auc", "pr_auc", "balanced_accuracy", "f1", "mcc", "brier", "ece_10"]
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


def chembl_get_json(path: str, params: dict | None = None, retries: int = 4) -> dict:
    url = f"{CH_EMBL_BASE}/{path}"
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            response = requests.get(url, params=params, timeout=45)
            response.raise_for_status()
            return response.json()
        except Exception as exc:
            last_error = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"ChEMBL request failed: {url} params={params} error={last_error}")


def dataset_targets_from_moleculeace(path: str, max_targets: int, activity_types: set[str]) -> pd.DataFrame:
    df = pd.read_csv(path)
    rows = []
    dataset_names = [str(x) for x in df["Dataset"].dropna().unique() if str(x).startswith("CHEMBL")]
    for name in sorted(dataset_names):
        parts = name.split("_")
        if len(parts) < 2:
            continue
        target_id = parts[0]
        activity_type = parts[1]
        if activity_type not in activity_types:
            continue
        rows.append({"dataset": name, "target_chembl_id": target_id, "standard_type": activity_type})
    out = pd.DataFrame(rows).drop_duplicates(["target_chembl_id", "standard_type"])
    if max_targets > 0:
        out = out.head(max_targets)
    return out.reset_index(drop=True)


def target_metadata(target_chembl_id: str) -> dict:
    data = chembl_get_json(f"target/{target_chembl_id}.json")
    sequence = ""
    accession = ""
    component_id = None
    for comp in data.get("target_components", []) or []:
        if str(comp.get("component_type", "")).upper() == "PROTEIN":
            accession = str(comp.get("accession") or "")
            component_id = comp.get("component_id")
            break
    if component_id is not None:
        comp_data = chembl_get_json(f"target_component/{component_id}.json")
        sequence = str(comp_data.get("sequence") or "")
    return {
        "target_chembl_id": target_chembl_id,
        "target_pref_name": str(data.get("pref_name") or target_chembl_id),
        "target_type": str(data.get("target_type") or ""),
        "organism": str(data.get("organism") or ""),
        "component_accession": accession,
        "target_sequence": sequence,
    }


def canonical_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def clean_pchembl(row: dict) -> float | None:
    pchembl = row.get("pchembl_value")
    try:
        if pchembl is not None and str(pchembl).strip():
            val = float(pchembl)
            if math.isfinite(val):
                return val
    except Exception:
        pass
    try:
        standard_value = float(row.get("standard_value"))
        standard_units = str(row.get("standard_units") or "").lower()
        if standard_units in {"nm", "nanomolar"} and standard_value > 0:
            return float(-math.log10(standard_value * 1e-9))
    except Exception:
        return None
    return None


def fetch_target_activities(
    target_chembl_id: str,
    standard_type: str,
    per_target_limit: int,
    page_size: int,
) -> list[dict]:
    rows = []
    offset = 0
    while len(rows) < per_target_limit:
        params = {
            "target_chembl_id": target_chembl_id,
            "standard_type": standard_type,
            "standard_units": "nM",
            "limit": page_size,
            "offset": offset,
        }
        data = chembl_get_json("activity.json", params=params)
        activities = data.get("activities", []) or []
        if not activities:
            break
        rows.extend(activities)
        if len(activities) < page_size:
            break
        offset += page_size
        time.sleep(0.15)
    return rows[:per_target_limit]


def fetch_or_load_raw_chembl(args: argparse.Namespace) -> pd.DataFrame:
    out = Path(args.raw_output)
    if args.reuse_raw and out.exists():
        return pd.read_csv(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    target_specs = dataset_targets_from_moleculeace(
        args.moleculeace_csv,
        args.max_targets,
        set(parse_csv_arg(args.activity_types)),
    )
    rows: list[dict] = []
    for spec in target_specs.itertuples(index=False):
        meta = target_metadata(str(spec.target_chembl_id))
        activities = fetch_target_activities(
            str(spec.target_chembl_id),
            str(spec.standard_type),
            int(args.per_target_limit),
            int(args.page_size),
        )
        kept = 0
        for act in activities:
            smiles = canonical_smiles(str(act.get("canonical_smiles") or ""))
            if not smiles:
                continue
            pchembl = clean_pchembl(act)
            if pchembl is None:
                continue
            relation = str(act.get("standard_relation") or "").strip()
            if args.require_equal_relation and relation not in {"=", ""}:
                continue
            year = act.get("document_year")
            try:
                year_val = int(year) if year is not None and not pd.isna(year) else np.nan
            except Exception:
                year_val = np.nan
            row = dict(meta)
            row.update(
                {
                    "dataset_source": str(spec.dataset),
                    "standard_type": str(spec.standard_type),
                    "activity_id": act.get("activity_id"),
                    "molecule_chembl_id": str(act.get("molecule_chembl_id") or ""),
                    "parent_molecule_chembl_id": str(act.get("parent_molecule_chembl_id") or act.get("molecule_chembl_id") or ""),
                    "canonical_smiles": smiles,
                    "pchembl_value": float(pchembl),
                    "standard_relation": relation,
                    "standard_value": act.get("standard_value"),
                    "standard_units": act.get("standard_units"),
                    "assay_chembl_id": str(act.get("assay_chembl_id") or ""),
                    "document_chembl_id": str(act.get("document_chembl_id") or ""),
                    "document_journal": str(act.get("document_journal") or ""),
                    "document_year": year_val,
                    "potential_duplicate": act.get("potential_duplicate"),
                }
            )
            rows.append(row)
            kept += 1
        print(f"fetched target={spec.target_chembl_id} type={spec.standard_type} raw={len(activities)} kept={kept}", flush=True)
    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("No ChEMBL activities fetched")
    df.to_csv(out, index=False)
    return df


def scaffold_key(smiles: str) -> str:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return f"BAD:{smiles}"
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    if scaffold:
        return scaffold
    return f"NOSCAFFOLD:{Chem.MolToSmiles(mol, canonical=True)}"


def make_pair_key(a: str, b: str) -> str:
    return "||".join(sorted([canonical_smiles(a), canonical_smiles(b)]))


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


def source_tokens(text: str) -> set[str]:
    if not isinstance(text, str) or not text:
        return set()
    return {x for x in text.split("|") if x}


def aggregate_activities(raw: pd.DataFrame, max_mols_per_target: int, seed: int) -> pd.DataFrame:
    raw = raw.copy()
    raw = raw[np.isfinite(raw["pchembl_value"].astype(float))]
    raw["document_year"] = pd.to_numeric(raw["document_year"], errors="coerce")
    group_cols = ["target_chembl_id", "canonical_smiles"]
    rows = []
    rng = np.random.default_rng(seed)
    for (target, smiles), group in raw.groupby(group_cols, sort=True):
        first = group.iloc[0]
        years = group["document_year"].dropna().astype(int).tolist()
        docs = sorted({str(x) for x in group["document_chembl_id"] if str(x)})
        assays = sorted({str(x) for x in group["assay_chembl_id"] if str(x)})
        rows.append(
            {
                "target": str(target),
                "target_chembl_id": str(target),
                "target_pref_name": str(first.get("target_pref_name", target)),
                "target_sequence": str(first.get("target_sequence", "")),
                "component_accession": str(first.get("component_accession", "")),
                "smiles": str(smiles),
                "pchembl": float(group["pchembl_value"].median()),
                "activity_count": int(len(group)),
                "document_year_median": float(np.median(years)) if years else np.nan,
                "document_year_min": int(min(years)) if years else -1,
                "document_year_max": int(max(years)) if years else -1,
                "document_ids": "|".join(docs[:50]),
                "assay_ids": "|".join(assays[:50]),
                "n_documents": int(len(docs)),
                "n_assays": int(len(assays)),
            }
        )
    agg = pd.DataFrame(rows)
    if max_mols_per_target > 0:
        kept = []
        for _, group in agg.groupby("target", sort=True):
            if len(group) > max_mols_per_target:
                idx = rng.choice(group.index.to_numpy(), size=max_mols_per_target, replace=False)
                kept.append(agg.loc[idx])
            else:
                kept.append(group)
        agg = pd.concat(kept, ignore_index=True)
    return agg


def build_pair_records(args: argparse.Namespace, raw: pd.DataFrame) -> pd.DataFrame:
    agg = aggregate_activities(raw, int(args.max_mols_per_target), int(args.downsample_seed))
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp_cache: dict[str, object] = {}
    scaffold_cache: dict[str, str] = {}
    rows: list[dict] = []
    rng = np.random.default_rng(int(args.downsample_seed))

    def fp(smiles: str):
        if smiles not in fp_cache:
            mol = Chem.MolFromSmiles(str(smiles))
            fp_cache[smiles] = gen.GetFingerprint(mol) if mol is not None else None
        return fp_cache[smiles]

    def scaffold(smiles: str) -> str:
        if smiles not in scaffold_cache:
            scaffold_cache[smiles] = scaffold_key(smiles)
        return scaffold_cache[smiles]

    for target, group in agg.groupby("target", sort=True):
        group = group.reset_index(drop=True)
        if len(group) < int(args.min_mols_per_target):
            continue
        smiles = group["smiles"].astype(str).tolist()
        y = group["pchembl"].to_numpy(dtype=np.float32)
        fps = [fp(s) for s in smiles]
        pos_rows: list[dict] = []
        neg_rows: list[dict] = []
        for i in range(len(group) - 1):
            fp_i = fps[i]
            if fp_i is None:
                continue
            sims = DataStructs.BulkTanimotoSimilarity(fp_i, fps[i + 1 :])
            for j, sim in enumerate(sims, start=i + 1):
                if sim < float(args.tanimoto_threshold):
                    continue
                delta = abs(float(y[i] - y[j]))
                if delta >= float(args.cliff_delta):
                    label = 1
                elif delta <= float(args.smooth_delta):
                    label = 0
                else:
                    continue
                a = smiles[i]
                b = smiles[j]
                sa = scaffold(a)
                sb = scaffold(b)
                doc_tokens = sorted(source_tokens(str(group.loc[i, "document_ids"])) | source_tokens(str(group.loc[j, "document_ids"])))
                assay_tokens = sorted(source_tokens(str(group.loc[i, "assay_ids"])) | source_tokens(str(group.loc[j, "assay_ids"])))
                years = [
                    int(group.loc[i, "document_year_max"]),
                    int(group.loc[j, "document_year_max"]),
                ]
                years = [year for year in years if year > 0]
                row = {
                    "target": str(target),
                    "target_pref_name": str(group.loc[i, "target_pref_name"]),
                    "component_accession": str(group.loc[i, "component_accession"]),
                    "target_sequence": str(group.loc[i, "target_sequence"]),
                    "smiles1": a,
                    "smiles2": b,
                    "label": int(label),
                    "activity_delta": float(delta),
                    "tanimoto": float(sim),
                    "scaffold1": sa,
                    "scaffold2": sb,
                    "pair_key": make_pair_key(a, b),
                    "pair_scaffold_key": "||".join(sorted([sa, sb])),
                    "document_source_key": "|".join(doc_tokens[:100]),
                    "assay_source_key": "|".join(assay_tokens[:100]),
                    "pair_year": int(max(years)) if years else -1,
                    "n_document_sources": int(len(doc_tokens)),
                    "n_assay_sources": int(len(assay_tokens)),
                }
                (pos_rows if label else neg_rows).append(row)
        if float(args.neg_pos_ratio) > 0 and pos_rows:
            max_neg = int(np.ceil(len(pos_rows) * float(args.neg_pos_ratio)))
            if len(neg_rows) > max_neg:
                keep = rng.choice(len(neg_rows), size=max_neg, replace=False)
                neg_rows = [neg_rows[int(k)] for k in keep]
        rows.extend(pos_rows)
        rows.extend(neg_rows)
        print(f"pairs target={target} mols={len(group)} pos={len(pos_rows)} neg={len(neg_rows)}", flush=True)
    pairs = pd.DataFrame(rows)
    if pairs.empty:
        raise RuntimeError("No ChEMBL pair rows constructed")
    pairs = pairs.sample(frac=1.0, random_state=int(args.downsample_seed)).reset_index(drop=True)
    Path(args.pair_output).parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(args.pair_output, index=False)
    return pairs


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
                    np.asarray([tanimoto, pop_a, pop_b, abs(pop_a - pop_b), common_count, float(len(diff_idx))], dtype=np.float32),
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


def family_folds(records: pd.DataFrame, seed: int, n_folds: int, n_hash: int, n_families: int) -> tuple[list[np.ndarray], dict[str, int]]:
    targets = np.array(sorted(records["target"].astype(str).unique())).astype(str)
    seq_by_target = records.groupby("target")["target_sequence"].first().to_dict()
    x = np.vstack([sequence_features(str(seq_by_target.get(target, "")), n_hash) for target in targets]).astype(np.float32)
    x = StandardScaler().fit_transform(x).astype(np.float32)
    n_clusters = max(n_folds, min(n_families, len(targets)))
    labels = KMeans(n_clusters=n_clusters, random_state=31337, n_init=20).fit_predict(x)
    target_to_family = {target: int(label) for target, label in zip(targets, labels)}
    stats = records.groupby("target")["label"].agg(["size", "sum"]).reset_index()
    stats["family"] = stats["target"].astype(str).map(target_to_family)
    family_stats = stats.groupby("family")[["size", "sum"]].sum().reset_index()
    family_stats = family_stats.sample(frac=1.0, random_state=seed).sort_values(["sum", "size"], ascending=False)
    folds: list[list[int]] = [[] for _ in range(n_folds)]
    fold_rows = np.zeros(n_folds, dtype=np.int64)
    fold_pos = np.zeros(n_folds, dtype=np.int64)
    for row in family_stats.itertuples(index=False):
        idx = int(np.argmin(fold_rows + 3 * fold_pos))
        folds[idx].append(int(row.family))
        fold_rows[idx] += int(row.size)
        fold_pos[idx] += int(row.sum)
    family_arr = records["target"].astype(str).map(target_to_family).to_numpy(dtype=np.int32)
    masks = [np.isin(family_arr, fold).astype(bool) for fold in folds if fold]
    return masks, target_to_family


def source_overlap_mask(records: pd.DataFrame, train_mask: np.ndarray, test_mask: np.ndarray) -> np.ndarray:
    test_docs: set[str] = set()
    test_assays: set[str] = set()
    for value in records.loc[test_mask, "document_source_key"].astype(str):
        test_docs.update(source_tokens(value))
    for value in records.loc[test_mask, "assay_source_key"].astype(str):
        test_assays.update(source_tokens(value))
    keep = np.ones(len(records), dtype=bool)
    if test_docs:
        keep &= ~records["document_source_key"].astype(str).map(lambda x: bool(source_tokens(x) & test_docs)).to_numpy()
    if test_assays:
        keep &= ~records["assay_source_key"].astype(str).map(lambda x: bool(source_tokens(x) & test_assays)).to_numpy()
    return train_mask & keep


def purge_train_mask(records: pd.DataFrame, base_train_mask: np.ndarray, test_mask: np.ndarray, split_mode: str) -> tuple[np.ndarray, float]:
    train_mask = base_train_mask.copy()
    n_base = int(base_train_mask.sum())
    if split_mode in {"target_family", "temporal_forward"}:
        return train_mask, 1.0
    if "scaffold" in split_mode:
        test_scaffolds = set(records.loc[test_mask, "scaffold1"].astype(str)) | set(records.loc[test_mask, "scaffold2"].astype(str))
        keep = ~(
            records["scaffold1"].astype(str).isin(test_scaffolds).to_numpy()
            | records["scaffold2"].astype(str).isin(test_scaffolds).to_numpy()
        )
        train_mask &= keep
    if "source" in split_mode:
        train_mask = source_overlap_mask(records, train_mask, test_mask)
    return train_mask, float(train_mask.sum() / max(n_base, 1))


def temporal_masks(records: pd.DataFrame, fold: int, n_folds: int) -> tuple[np.ndarray, np.ndarray, int]:
    valid = records["pair_year"].astype(int).to_numpy() > 0
    years = np.array(sorted(pd.unique(records.loc[valid, "pair_year"].astype(int))))
    if len(years) < 4:
        raise RuntimeError("Not enough valid years for temporal split")
    quantiles = np.linspace(0.55, 0.85, n_folds)
    cutoff = int(np.quantile(years, quantiles[min(fold, len(quantiles) - 1)]))
    train_mask = valid & (records["pair_year"].astype(int).to_numpy() <= cutoff)
    test_mask = valid & (records["pair_year"].astype(int).to_numpy() > cutoff)
    return train_mask, test_mask, cutoff


def expected_calibration_error(y: np.ndarray, prob: np.ndarray, n_bins: int = 10) -> float:
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    for left, right in zip(edges[:-1], edges[1:]):
        mask = (prob >= left) & (prob <= right) if right == 1.0 else (prob >= left) & (prob < right)
        if np.any(mask):
            ece += float(mask.mean()) * abs(float(prob[mask].mean()) - float(y[mask].mean()))
    return float(ece)


def row_metrics(y: np.ndarray, prob: np.ndarray) -> dict[str, float]:
    pred = (prob >= 0.5).astype(np.int32)
    out: dict[str, float] = {}
    out["roc_auc"] = float(roc_auc_score(y, prob)) if len(np.unique(y)) == 2 else float("nan")
    out["pr_auc"] = float(average_precision_score(y, prob)) if len(np.unique(y)) == 2 else float("nan")
    out["balanced_accuracy"] = float(balanced_accuracy_score(y, pred)) if len(np.unique(y)) == 2 else float("nan")
    out["f1"] = float(f1_score(y, pred, zero_division=0))
    out["mcc"] = float(matthews_corrcoef(y, pred)) if len(np.unique(pred)) > 1 and len(np.unique(y)) > 1 else 0.0
    out["brier"] = float(brier_score_loss(y, prob))
    out["ece_10"] = expected_calibration_error(y, prob, 10)
    return out


def make_model(seed: int, y_train: np.ndarray, threads: int, n_estimators: int, max_depth: int) -> XGBClassifier:
    positives = float(y_train.sum())
    negatives = float(len(y_train) - positives)
    return XGBClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        learning_rate=0.04,
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


def audit_overlaps(records: pd.DataFrame, train_mask: np.ndarray, test_mask: np.ndarray) -> dict[str, float]:
    train_targets = set(records.loc[train_mask, "target"].astype(str))
    train_ligands = set(records.loc[train_mask, "smiles1"].astype(str)) | set(records.loc[train_mask, "smiles2"].astype(str))
    train_scaffolds = set(records.loc[train_mask, "scaffold1"].astype(str)) | set(records.loc[train_mask, "scaffold2"].astype(str))
    train_pairs = set(records.loc[train_mask, "pair_key"].astype(str))
    train_pair_scaffolds = set(records.loc[train_mask, "pair_scaffold_key"].astype(str))
    train_docs: set[str] = set()
    train_assays: set[str] = set()
    for value in records.loc[train_mask, "document_source_key"].astype(str):
        train_docs.update(source_tokens(value))
    for value in records.loc[train_mask, "assay_source_key"].astype(str):
        train_assays.update(source_tokens(value))
    doc_overlap = records.loc[test_mask, "document_source_key"].astype(str).map(lambda x: bool(source_tokens(x) & train_docs)).to_numpy()
    assay_overlap = records.loc[test_mask, "assay_source_key"].astype(str).map(lambda x: bool(source_tokens(x) & train_assays)).to_numpy()
    ligand_overlap = (
        records.loc[test_mask, "smiles1"].astype(str).isin(train_ligands).to_numpy()
        | records.loc[test_mask, "smiles2"].astype(str).isin(train_ligands).to_numpy()
    )
    scaffold_overlap = (
        records.loc[test_mask, "scaffold1"].astype(str).isin(train_scaffolds).to_numpy()
        | records.loc[test_mask, "scaffold2"].astype(str).isin(train_scaffolds).to_numpy()
    )
    train_years = records.loc[train_mask & (records["pair_year"].astype(int) > 0), "pair_year"].astype(int)
    test_years = records.loc[test_mask & (records["pair_year"].astype(int) > 0), "pair_year"].astype(int)
    return {
        "target_overlap_rate": float(records.loc[test_mask, "target"].astype(str).isin(train_targets).mean()),
        "exact_ligand_overlap_rate": float(ligand_overlap.mean()) if len(ligand_overlap) else float("nan"),
        "scaffold_overlap_rate": float(scaffold_overlap.mean()) if len(scaffold_overlap) else float("nan"),
        "exact_pair_overlap_rate": float(records.loc[test_mask, "pair_key"].astype(str).isin(train_pairs).mean()),
        "pair_scaffold_overlap_rate": float(records.loc[test_mask, "pair_scaffold_key"].astype(str).isin(train_pair_scaffolds).mean()),
        "document_source_overlap_rate": float(doc_overlap.mean()) if len(doc_overlap) else float("nan"),
        "assay_source_overlap_rate": float(assay_overlap.mean()) if len(assay_overlap) else float("nan"),
        "train_max_year": int(train_years.max()) if len(train_years) else -1,
        "test_min_year": int(test_years.min()) if len(test_years) else -1,
    }


def build_x(diff_bits: sparse.csr_matrix, scalars: np.ndarray, train_mask: np.ndarray, variant: str):
    if variant == "ecfp_absdiff_xgb":
        return diff_bits
    scaler = StandardScaler().fit(scalars[train_mask])
    scalar_sparse = sparse.csr_matrix(scaler.transform(scalars).astype(np.float32))
    if variant == "rich_scalars_only":
        return scalar_sparse
    if variant == "rich_diff_scalars":
        return sparse.hstack([diff_bits, scalar_sparse], format="csr")
    raise ValueError(variant)


def run_tasks(args: argparse.Namespace, records: pd.DataFrame) -> pd.DataFrame:
    diff_bits, scalars = build_feature_matrices(records)
    variants = parse_csv_arg(args.variants)
    split_modes = parse_csv_arg(args.split_modes)
    seeds = [int(x) for x in parse_csv_arg(args.seeds)]
    rows = []
    for seed in seeds:
        masks, target_to_family = family_folds(records, seed, int(args.folds), int(args.sequence_hash_features), int(args.target_family_clusters))
        family_arr = records["target"].astype(str).map(target_to_family).to_numpy(dtype=np.int32)
        for fold in range(int(args.folds)):
            for split_mode in split_modes:
                if split_mode.startswith("temporal"):
                    base_train_mask, test_mask, cutoff = temporal_masks(records, fold, int(args.folds))
                else:
                    test_mask = masks[fold]
                    base_train_mask = ~test_mask
                    cutoff = -1
                train_mask, retention = purge_train_mask(records, base_train_mask, test_mask, split_mode)
                if train_mask.sum() < int(args.min_train_rows) or test_mask.sum() < int(args.min_test_rows):
                    print(f"skip split={split_mode} seed={seed} fold={fold} train={train_mask.sum()} test={test_mask.sum()}", flush=True)
                    continue
                y_train = records.loc[train_mask, "label"].to_numpy(dtype=np.int32)
                y_test = records.loc[test_mask, "label"].to_numpy(dtype=np.int32)
                if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
                    print(f"skip one-class split={split_mode} seed={seed} fold={fold}", flush=True)
                    continue
                for variant in variants:
                    x = build_x(diff_bits, scalars, train_mask, variant)
                    model = make_model(seed, y_train, int(args.model_threads), int(args.n_estimators), int(args.max_depth))
                    model.fit(x[train_mask], y_train)
                    prob = model.predict_proba(x[test_mask])[:, 1].astype(np.float32)
                    train_families = set(family_arr[train_mask])
                    test_families = set(family_arr[test_mask])
                    row = {
                        "benchmark": "ChEMBL_raw_MoleculeACE_targets",
                        "split_mode": split_mode,
                        "variant": variant,
                        "seed": seed,
                        "fold": fold,
                        "temporal_cutoff_year": int(cutoff),
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
                    }
                    row.update(audit_overlaps(records, train_mask, test_mask))
                    row.update(row_metrics(y_test, prob))
                    rows.append(row)
                    print(
                        f"[{split_mode} seed={seed} fold={fold} {variant}] "
                        f"roc={row['roc_auc']:.4f} pr={row['pr_auc']:.4f} "
                        f"train={row['n_train']} test={row['n_test']} "
                        f"docov={row['document_source_overlap_rate']:.3f} assayov={row['assay_source_overlap_rate']:.3f}",
                        flush=True,
                    )
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    audit_cols = [
        "n_train",
        "n_test",
        "train_retention_after_purge",
        "train_positive_rate",
        "test_positive_rate",
        "target_overlap_rate",
        "family_overlap_rate",
        "exact_ligand_overlap_rate",
        "scaffold_overlap_rate",
        "exact_pair_overlap_rate",
        "pair_scaffold_overlap_rate",
        "document_source_overlap_rate",
        "assay_source_overlap_rate",
        "train_max_year",
        "test_min_year",
    ]
    cols = [c for c in METRICS + audit_cols if c in df.columns]
    out = df.groupby(["split_mode", "variant"], as_index=False)[cols].agg(["mean", "std", "count"])
    out.columns = ["_".join([x for x in col if x]).strip("_") for col in out.columns.to_flat_index()]
    return out.reset_index(drop=True)


def paired_deltas(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for split in sorted(df["split_mode"].astype(str).unique()):
        sub = df[df["split_mode"] == split]
        variants = sorted(sub["variant"].astype(str).unique())
        for baseline in variants:
            base = sub[sub["variant"] == baseline].set_index(["seed", "fold"])
            for variant in variants:
                if variant == baseline:
                    continue
                cur = sub[sub["variant"] == variant].set_index(["seed", "fold"])
                idx = cur.index.intersection(base.index)
                if len(idx) == 0:
                    continue
                for metric in [m for m in METRICS if m in cur.columns and m in base.columns]:
                    diff = cur.loc[idx, metric].astype(float) - base.loc[idx, metric].astype(float)
                    diff = diff.replace([np.inf, -np.inf], np.nan).dropna()
                    if diff.empty:
                        continue
                    favorable = diff > 0 if HIGHER_IS_BETTER[metric] else diff < 0
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


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_float_dtype(out[col]):
            out[col] = out[col].map(lambda x: "" if pd.isna(x) else f"{float(x):.6g}")
    return out.to_markdown(index=False)


def write_report(raw: pd.DataFrame, pairs: pd.DataFrame, results: pd.DataFrame, summary: pd.DataFrame, deltas: pd.DataFrame, out: Path) -> None:
    compact_cols = [
        "split_mode",
        "variant",
        "roc_auc_mean",
        "roc_auc_std",
        "pr_auc_mean",
        "pr_auc_std",
        "target_overlap_rate_mean",
        "family_overlap_rate_mean",
        "exact_ligand_overlap_rate_mean",
        "scaffold_overlap_rate_mean",
        "document_source_overlap_rate_mean",
        "assay_source_overlap_rate_mean",
        "train_retention_after_purge_mean",
        "roc_auc_count",
    ]
    compact = summary[[c for c in compact_cols if c in summary.columns]]
    delta_focus = deltas[(deltas["baseline"] == "ecfp_absdiff_xgb") & (deltas["metric"].isin(["roc_auc", "pr_auc", "f1", "mcc"]))].copy()
    if not delta_focus.empty:
        delta_focus = delta_focus.sort_values(["split_mode", "metric", "delta_mean"], ascending=[True, True, False])
    lines = [
        "# Raw ChEMBL Source/Temporal Target-Family Pair-Cliff Benchmark",
        "",
        "## Scope",
        "",
        "This benchmark queries ChEMBL raw activity records for MoleculeACE target IDs, keeps assay/document/year metadata, constructs high-similarity pair-cliff labels, and evaluates source/temporal leakage controls.",
        "",
        "## Dataset",
        "",
        f"- raw ChEMBL activity rows: {len(raw):,}",
        f"- constructed pair rows: {len(pairs):,}",
        f"- targets: {pairs['target'].nunique():,}",
        f"- positive pair rate: {pairs['label'].mean():.4f}",
        f"- year range: {int(pairs.loc[pairs['pair_year'] > 0, 'pair_year'].min())} - {int(pairs.loc[pairs['pair_year'] > 0, 'pair_year'].max())}",
        "",
        "## Compact Summary",
        "",
        markdown_table(compact),
        "",
        "## Paired Deltas vs ECFP",
        "",
        markdown_table(delta_focus),
        "",
        "## Interpretation",
        "",
        "- `family_scaffold_source_purged` is the hardest target-family split because target family, scaffold, and assay/document sources are all purged from training.",
        "- `temporal_forward` measures forward-in-time generalization and should be interpreted separately from target-family OOD.",
        "- Use this as a raw-data external evidence block. If gains are weaker than ACNet/MoleculeACE, report it honestly as a deployment-shift stress test.",
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--moleculeace-csv", default="external/MTPNet/data/MoleculeACE.csv")
    parser.add_argument("--raw-output", default="data/chembl_raw_moleculeace_target_activities.csv")
    parser.add_argument("--pair-output", default="data/chembl_raw_moleculeace_target_pairs.csv")
    parser.add_argument("--reuse-raw", action="store_true")
    parser.add_argument("--reuse-pairs", action="store_true")
    parser.add_argument("--max-targets", type=int, default=30)
    parser.add_argument("--activity-types", default="Ki,IC50,Kd")
    parser.add_argument("--per-target-limit", type=int, default=1200)
    parser.add_argument("--page-size", type=int, default=500)
    parser.add_argument("--require-equal-relation", action="store_true")
    parser.add_argument("--min-mols-per-target", type=int, default=40)
    parser.add_argument("--max-mols-per-target", type=int, default=650)
    parser.add_argument("--tanimoto-threshold", type=float, default=0.7)
    parser.add_argument("--cliff-delta", type=float, default=1.0)
    parser.add_argument("--smooth-delta", type=float, default=0.3)
    parser.add_argument("--neg-pos-ratio", type=float, default=5.0)
    parser.add_argument("--downsample-seed", type=int, default=2026)
    parser.add_argument("--variants", default="ecfp_absdiff_xgb,rich_scalars_only,rich_diff_scalars")
    parser.add_argument("--split-modes", default="target_family,family_scaffold_source_purged,temporal_forward,temporal_scaffold_source_purged")
    parser.add_argument("--seeds", default="0,1,2,3,4")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--target-family-clusters", type=int, default=10)
    parser.add_argument("--sequence-hash-features", type=int, default=128)
    parser.add_argument("--min-train-rows", type=int, default=500)
    parser.add_argument("--min-test-rows", type=int, default=200)
    parser.add_argument("--n-estimators", type=int, default=80)
    parser.add_argument("--max-depth", type=int, default=4)
    parser.add_argument("--model-threads", type=int, default=2)
    parser.add_argument("--output-prefix", default="results/chembl_raw_source_temporal_pairs_seed0_4_n80")
    args = parser.parse_args()

    raw = fetch_or_load_raw_chembl(args)
    if args.reuse_pairs and Path(args.pair_output).exists():
        pairs = pd.read_csv(args.pair_output)
    else:
        pairs = build_pair_records(args, raw)
    results = run_tasks(args, pairs)
    if results.empty:
        raise RuntimeError("No experiment rows produced")
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    results.to_json(f"{prefix}.jsonl", orient="records", lines=True)
    summary = summarize(results)
    deltas = paired_deltas(results)
    summary.to_csv(f"{prefix}_summary.csv", index=False)
    deltas.to_csv(f"{prefix}_paired_deltas.csv", index=False)
    write_report(raw, pairs, results, summary, deltas, Path(f"{prefix}_analysis.md"))
    manifest = {
        "raw_rows": int(len(raw)),
        "pair_rows": int(len(pairs)),
        "targets": int(pairs["target"].nunique()),
        "positive_rate": float(pairs["label"].mean()),
        "output_prefix": str(prefix),
    }
    Path(f"{prefix}_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
