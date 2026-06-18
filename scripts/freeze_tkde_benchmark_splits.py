from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_acnet_external_baseline as acnet
import run_acnet_rich_pair_external as rich
import run_acnet_target_family_ablation_joint as joint
import run_chembl_raw_source_temporal_pairs as chembl
import run_moleculeace_target_family_pairs as molace
from tkde_benchmark_evaluator import sha256_file


def parse_csv_arg(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()]


def row_hash(values: list[object]) -> str:
    text = "||".join(str(v) for v in values)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:20]


def write_npz(path: Path, train_idx: np.ndarray, test_idx: np.ndarray) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, train_idx=train_idx.astype(np.int32), test_idx=test_idx.astype(np.int32))
    return sha256_file(path)


def save_dataset_rows(path: Path, row_ids: list[str]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame({"row_index": np.arange(len(row_ids), dtype=np.int32), "row_id": row_ids})
    df.to_csv(path, index=False)
    return sha256_file(path)


def freeze_raw_pair_dataset(
    benchmark: str,
    pair_csv: str,
    split_modes: list[str],
    seeds: list[int],
    folds: int,
    out_dir: Path,
    sequence_hash_features: int,
    target_family_clusters: int,
) -> list[dict]:
    records = pd.read_csv(pair_csv)
    row_ids = [
        row_hash([r.target, r.smiles1, r.smiles2, r.label, getattr(r, "activity_delta", "")])
        for r in records.itertuples(index=False)
    ]
    row_file = out_dir / f"{benchmark}_rows.csv"
    row_hash_value = save_dataset_rows(row_file, row_ids)
    rows = []
    for seed in seeds:
        masks, target_to_family = chembl.family_folds(
            records,
            seed,
            folds,
            sequence_hash_features,
            target_family_clusters,
        )
        for fold, test_mask in enumerate(masks[:folds]):
            for split_mode in split_modes:
                if split_mode.startswith("temporal"):
                    if "pair_year" not in records.columns or (pd.to_numeric(records["pair_year"], errors="coerce") > 0).sum() == 0:
                        continue
                    base_train_mask, test_mask_temporal, cutoff = chembl.temporal_masks(records, fold, folds)
                    train_mask, retention = chembl.purge_train_mask(records, base_train_mask, test_mask_temporal, split_mode)
                    cur_test = test_mask_temporal
                else:
                    base_train_mask = ~test_mask
                    train_mask, retention = chembl.purge_train_mask(records, base_train_mask, test_mask, split_mode)
                    cur_test = test_mask
                    cutoff = -1
                train_idx = np.flatnonzero(train_mask)
                test_idx = np.flatnonzero(cur_test)
                split_file = out_dir / f"{benchmark}_{split_mode}_seed{seed}_fold{fold}.npz"
                split_hash = write_npz(split_file, train_idx, test_idx)
                audit = chembl.audit_overlaps(records, train_mask, cur_test)
                rows.append(
                    {
                        "benchmark": benchmark,
                        "source_path": pair_csv,
                        "source_sha256": sha256_file(pair_csv),
                        "row_file": str(row_file),
                        "row_file_sha256": row_hash_value,
                        "split_mode": split_mode,
                        "seed": seed,
                        "fold": fold,
                        "split_file": str(split_file),
                        "split_file_sha256": split_hash,
                        "n_rows": int(len(records)),
                        "n_train": int(len(train_idx)),
                        "n_test": int(len(test_idx)),
                        "train_retention_after_purge": float(retention),
                        "temporal_cutoff_year": int(cutoff),
                        **audit,
                    }
                )
    return rows


def factorize_pair_values(left: np.ndarray, right: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = pd.Series(np.concatenate([left.astype(str), right.astype(str)]), dtype="string")
    codes, _ = pd.factorize(values, sort=False)
    n = len(left)
    return codes[:n].astype(np.int32), codes[n:].astype(np.int32)


def overlap_rate(test_values: np.ndarray, train_values: np.ndarray) -> float:
    if len(test_values) == 0:
        return 0.0
    train_unique = np.unique(train_values)
    return float(np.isin(test_values, train_unique).mean())


def pair_overlap_rate(
    test_left: np.ndarray,
    test_right: np.ndarray,
    train_left: np.ndarray,
    train_right: np.ndarray,
) -> float:
    if len(test_left) == 0:
        return 0.0
    train_unique = np.unique(np.concatenate([train_left, train_right]))
    return float((np.isin(test_left, train_unique) | np.isin(test_right, train_unique)).mean())


def acnet_fast_audit(
    records: pd.DataFrame,
    train_mask: np.ndarray,
    test_mask: np.ndarray,
    codes: dict[str, np.ndarray],
) -> dict[str, float]:
    return {
        "target_overlap_rate": overlap_rate(codes["target"][test_mask], codes["target"][train_mask]),
        "exact_ligand_overlap_rate": pair_overlap_rate(
            codes["ligand_a"][test_mask],
            codes["ligand_b"][test_mask],
            codes["ligand_a"][train_mask],
            codes["ligand_b"][train_mask],
        ),
        "scaffold_overlap_rate": pair_overlap_rate(
            codes["scaffold_a"][test_mask],
            codes["scaffold_b"][test_mask],
            codes["scaffold_a"][train_mask],
            codes["scaffold_b"][train_mask],
        ),
        "exact_pair_overlap_rate": overlap_rate(codes["pair"][test_mask], codes["pair"][train_mask]),
        "pair_scaffold_overlap_rate": overlap_rate(codes["pair_scaffold"][test_mask], codes["pair_scaffold"][train_mask]),
    }


def acnet_fast_purge(
    base_train_mask: np.ndarray,
    test_mask: np.ndarray,
    split_mode: str,
    codes: dict[str, np.ndarray],
) -> tuple[np.ndarray, float]:
    train_mask = base_train_mask.copy()
    n_base_train = int(base_train_mask.sum())
    if split_mode == "target_family":
        return train_mask, 1.0
    if split_mode != "family_scaffold_purged":
        raise ValueError(split_mode)
    test_scaffolds = np.unique(np.concatenate([codes["scaffold_a"][test_mask], codes["scaffold_b"][test_mask]]))
    keep = ~(np.isin(codes["scaffold_a"], test_scaffolds) | np.isin(codes["scaffold_b"], test_scaffolds))
    train_mask &= keep
    return train_mask, float(train_mask.sum() / max(n_base_train, 1))


def setup_acnet(args: argparse.Namespace) -> tuple[pd.DataFrame, dict[int, list[np.ndarray]], dict[str, np.ndarray]]:
    records = acnet.load_acnet_json(args.acnet_input)
    records = acnet.downsample_negatives(records, args.neg_pos_ratio, args.downsample_seed)
    target_emb = joint.load_target_embeddings(args.target_embeddings)
    scaffold_a, scaffold_b, pair_key, pair_scaffold_key = joint.make_scaffold_arrays(records)
    ligand_a, ligand_b = factorize_pair_values(
        records["smiles1"].astype(str).to_numpy(),
        records["smiles2"].astype(str).to_numpy(),
    )
    scaffold_a_code, scaffold_b_code = factorize_pair_values(scaffold_a, scaffold_b)
    codes = {
        "target": pd.factorize(records["target"].astype(str), sort=False)[0].astype(np.int32),
        "ligand_a": ligand_a,
        "ligand_b": ligand_b,
        "scaffold_a": scaffold_a_code,
        "scaffold_b": scaffold_b_code,
        "pair": pd.factorize(pd.Series(pair_key, dtype="string"), sort=False)[0].astype(np.int32),
        "pair_scaffold": pd.factorize(pd.Series(pair_scaffold_key, dtype="string"), sort=False)[0].astype(np.int32),
    }
    fold_masks = {}
    for seed in [int(x) for x in parse_csv_arg(args.seeds)]:
        masks, mapping = joint.target_family_folds(
            records,
            target_emb,
            seed,
            args.family_cluster_seed,
            args.folds,
            args.target_family_clusters,
        )
        fold_masks[seed] = masks
    return records, fold_masks, codes


def freeze_acnet(args: argparse.Namespace, out_dir: Path) -> list[dict]:
    seeds = [int(x) for x in parse_csv_arg(args.seeds)]
    split_modes = parse_csv_arg(args.acnet_split_modes)
    records, fold_masks, codes = setup_acnet(args)
    row_ids = [
        row_hash([r.target, r.smiles1, r.smiles2, r.label])
        for r in records.itertuples(index=False)
    ]
    row_file = out_dir / "acnet_rows.csv"
    row_hash_value = save_dataset_rows(row_file, row_ids)
    rows = []
    for seed in seeds:
        for fold, test_mask in enumerate(fold_masks[seed]):
            base_train_mask = ~test_mask
            for split_mode in split_modes:
                train_mask, retention = acnet_fast_purge(base_train_mask, test_mask, split_mode, codes)
                train_idx = np.flatnonzero(train_mask)
                test_idx = np.flatnonzero(test_mask)
                split_file = out_dir / f"acnet_{split_mode}_seed{seed}_fold{fold}.npz"
                split_hash = write_npz(split_file, train_idx, test_idx)
                rows.append(
                    {
                        "benchmark": "acnet",
                        "source_path": args.acnet_input,
                        "source_sha256": sha256_file(args.acnet_input),
                        "row_file": str(row_file),
                        "row_file_sha256": row_hash_value,
                        "split_mode": split_mode,
                        "seed": seed,
                        "fold": fold,
                        "split_file": str(split_file),
                        "split_file_sha256": split_hash,
                        "n_rows": int(len(records)),
                        "n_train": int(len(train_idx)),
                        "n_test": int(len(test_idx)),
                        "train_retention_after_purge": float(retention),
                        **acnet_fast_audit(records, train_mask, test_mask, codes),
                    }
                )
    return rows


def freeze_moleculeace(args: argparse.Namespace, out_dir: Path) -> list[dict]:
    seeds = [int(x) for x in parse_csv_arg(args.seeds)]
    records = molace.load_moleculeace_pairs(
        args.moleculeace_csv,
        args.tanimoto_threshold,
        args.cliff_delta,
        args.smooth_delta,
        args.neg_pos_ratio,
        args.downsample_seed,
    )
    molace.RECORDS = records
    row_ids = [
        row_hash([r.target, r.smiles1, r.smiles2, r.label, r.activity_delta])
        for r in records.itertuples(index=False)
    ]
    row_file = out_dir / "moleculeace_rows.csv"
    row_hash_value = save_dataset_rows(row_file, row_ids)
    rows = []
    for seed in seeds:
        masks, mapping = molace.target_family_folds(
            records,
            seed,
            args.family_cluster_seed,
            args.folds,
            args.target_family_clusters,
            args.sequence_hash_features,
        )
        for fold, test_mask in enumerate(masks[: args.folds]):
            base_train_mask = ~test_mask
            for split_mode in parse_csv_arg(args.moleculeace_split_modes):
                train_mask, retention = molace.purge_training_mask(base_train_mask, test_mask, split_mode)
                train_idx = np.flatnonzero(train_mask)
                test_idx = np.flatnonzero(test_mask)
                split_file = out_dir / f"moleculeace_{split_mode}_seed{seed}_fold{fold}.npz"
                split_hash = write_npz(split_file, train_idx, test_idx)
                rows.append(
                    {
                        "benchmark": "moleculeace",
                        "source_path": args.moleculeace_csv,
                        "source_sha256": sha256_file(args.moleculeace_csv),
                        "row_file": str(row_file),
                        "row_file_sha256": row_hash_value,
                        "split_mode": split_mode,
                        "seed": seed,
                        "fold": fold,
                        "split_file": str(split_file),
                        "split_file_sha256": split_hash,
                        "n_rows": int(len(records)),
                        "n_train": int(len(train_idx)),
                        "n_test": int(len(test_idx)),
                        "train_retention_after_purge": float(retention),
                        **molace.audit_overlaps(records, train_mask, test_mask, mapping),
                    }
                )
    return rows


def write_report(manifest: pd.DataFrame, out: Path) -> None:
    compact = manifest.groupby(["benchmark", "split_mode"], as_index=False).agg(
        n_splits=("split_file", "count"),
        n_rows=("n_rows", "mean"),
        n_train=("n_train", "mean"),
        n_test=("n_test", "mean"),
        train_retention_after_purge=("train_retention_after_purge", "mean"),
        target_overlap_rate=("target_overlap_rate", "mean"),
        exact_ligand_overlap_rate=("exact_ligand_overlap_rate", "mean"),
        scaffold_overlap_rate=("scaffold_overlap_rate", "mean"),
    )
    lines = [
        "# TKDE Benchmark Split Freeze Report",
        "",
        "This report freezes benchmark rows and train/test index files for the main hard-OOD settings. Split files are stored as compressed NumPy archives with `train_idx` and `test_idx` arrays; row identity files carry stable row hashes.",
        "",
        compact.to_markdown(index=False),
        "",
        "## Integrity Notes",
        "",
        "- `*_rows.csv` files bind each row index to a stable hash derived from target, molecule pair, label, and activity delta when available.",
        "- Manifest rows include SHA-256 checksums for source files, row files, and split files.",
        "- Split audit columns are included so downstream users can verify target, ligand, scaffold, and source overlap properties.",
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="splits/tkde_protocol")
    parser.add_argument("--manifest", default="results/tkde_benchmark_split_manifest.csv")
    parser.add_argument("--report", default="results/tkde_benchmark_split_freeze_report.md")
    parser.add_argument("--benchmarks", default="acnet,moleculeace,chembl_raw,bindingdb_raw")
    parser.add_argument("--seeds", default="0,1,2,3,4")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--family-cluster-seed", type=int, default=31337)
    parser.add_argument("--target-family-clusters", type=int, default=10)
    parser.add_argument("--sequence-hash-features", type=int, default=128)
    parser.add_argument("--neg-pos-ratio", type=float, default=5.0)
    parser.add_argument("--downsample-seed", type=int, default=2026)
    parser.add_argument("--pair-proj-dim", type=int, default=32)
    parser.add_argument("--target-embeddings", default="external/ACNet_data/target_esm2_t6_embeddings_all.csv")
    parser.add_argument("--acnet-input", default="external/ACNet/ACNet/ACComponents/ACDataset/data_files/generated_datasets/MMP_AC.json")
    parser.add_argument("--acnet-split-modes", default="target_family,family_scaffold_purged")
    parser.add_argument("--moleculeace-csv", default="external/MTPNet/data/MoleculeACE.csv")
    parser.add_argument("--moleculeace-split-modes", default="target_family,family_scaffold_purged")
    parser.add_argument("--chembl-pairs", default="data/chembl_raw_moleculeace_target_pairs.csv")
    parser.add_argument("--chembl-split-modes", default="target_family,family_scaffold_source_purged,temporal_scaffold_source_purged")
    parser.add_argument("--bindingdb-pairs", default="data/bindingdb_raw_moleculeace_target_pairs.csv")
    parser.add_argument("--bindingdb-split-modes", default="target_family,family_scaffold_source_purged")
    parser.add_argument("--tanimoto-threshold", type=float, default=0.7)
    parser.add_argument("--cliff-delta", type=float, default=1.0)
    parser.add_argument("--smooth-delta", type=float, default=0.3)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    benchmarks = set(parse_csv_arg(args.benchmarks))
    seeds = [int(x) for x in parse_csv_arg(args.seeds)]
    rows: list[dict] = []
    if "acnet" in benchmarks:
        rows.extend(freeze_acnet(args, out_dir))
    if "moleculeace" in benchmarks:
        rows.extend(freeze_moleculeace(args, out_dir))
    if "chembl_raw" in benchmarks:
        rows.extend(
            freeze_raw_pair_dataset(
                "chembl_raw",
                args.chembl_pairs,
                parse_csv_arg(args.chembl_split_modes),
                seeds,
                args.folds,
                out_dir,
                args.sequence_hash_features,
                args.target_family_clusters,
            )
        )
    if "bindingdb_raw" in benchmarks:
        rows.extend(
            freeze_raw_pair_dataset(
                "bindingdb_raw",
                args.bindingdb_pairs,
                parse_csv_arg(args.bindingdb_split_modes),
                seeds,
                args.folds,
                out_dir,
                args.sequence_hash_features,
                args.target_family_clusters,
            )
        )
    manifest = pd.DataFrame(rows)
    Path(args.manifest).parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.manifest, index=False)
    write_report(manifest, Path(args.report))
    print(json.dumps({"manifest": args.manifest, "rows": int(len(manifest)), "output_dir": str(out_dir)}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
