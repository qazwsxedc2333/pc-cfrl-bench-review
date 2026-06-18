from __future__ import annotations

import argparse
import json
import math
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from run_chembl_raw_source_temporal_pairs import (
    audit_overlaps,
    build_feature_matrices,
    build_x,
    canonical_smiles,
    dataset_targets_from_moleculeace,
    family_folds,
    make_model,
    make_pair_key,
    markdown_table,
    paired_deltas,
    parse_csv_arg,
    purge_train_mask,
    row_metrics,
    scaffold_key,
    source_tokens,
    summarize,
    target_metadata,
)


BINDINGDB_REST = "https://bindingdb.org/rest/getLigandsByUniprots"


def clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def parse_bindingdb_affinity(value: object) -> tuple[float | None, str]:
    text = clean_text(value)
    if not text:
        return None, ""
    relation = ""
    if text[0] in "<>~=":
        relation = text[0]
    match = re.search(r"([0-9]+(?:\.[0-9]+)?(?:[eE][-+]?[0-9]+)?)", text.replace(",", ""))
    if not match:
        return None, relation
    try:
        val = float(match.group(1))
    except ValueError:
        return None, relation
    if not math.isfinite(val) or val <= 0:
        return None, relation
    return val, relation


def load_target_specs(args: argparse.Namespace) -> pd.DataFrame:
    metadata = Path(args.target_metadata_csv)
    if metadata.exists():
        df = pd.read_csv(metadata)
        cols = [
            "target_chembl_id",
            "target_pref_name",
            "component_accession",
            "target_sequence",
        ]
        out = df[[c for c in cols if c in df.columns]].drop_duplicates("target_chembl_id")
    else:
        specs = dataset_targets_from_moleculeace(
            args.moleculeace_csv,
            int(args.max_targets),
            set(parse_csv_arg(args.activity_types)),
        )
        rows = []
        for target in specs["target_chembl_id"].drop_duplicates():
            meta = target_metadata(str(target))
            rows.append(
                {
                    "target_chembl_id": str(target),
                    "target_pref_name": str(meta.get("target_pref_name", target)),
                    "component_accession": str(meta.get("component_accession", "")),
                    "target_sequence": str(meta.get("target_sequence", "")),
                }
            )
        out = pd.DataFrame(rows)
    out = out.rename(columns={"component_accession": "uniprot_accession"})
    out["uniprot_accession"] = out["uniprot_accession"].astype(str).str.strip()
    out = out[out["uniprot_accession"].ne("")]
    out = out.drop_duplicates("uniprot_accession").sort_values("target_chembl_id")
    if int(args.max_targets) > 0:
        out = out.head(int(args.max_targets))
    return out.reset_index(drop=True)


def bindingdb_get(uniprot: str, cutoff_nm: float, retries: int = 4) -> list[dict]:
    params = {
        "uniprot": uniprot,
        "cutoff": str(float(cutoff_nm)),
        "response": "application/json",
    }
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            response = requests.get(BINDINGDB_REST, params=params, timeout=90)
            response.raise_for_status()
            data = response.json()
            return data.get("getLindsByUniprotsResponse", {}).get("affinities", []) or []
        except Exception as exc:
            last_error = exc
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"BindingDB request failed for {uniprot}: {last_error}")


def fetch_or_load_raw_bindingdb(args: argparse.Namespace) -> pd.DataFrame:
    out = Path(args.raw_output)
    if args.reuse_raw and out.exists():
        return pd.read_csv(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    targets = load_target_specs(args)
    allowed_types = set(parse_csv_arg(args.activity_types))
    rows: list[dict] = []
    for target in targets.itertuples(index=False):
        uniprot = str(target.uniprot_accession)
        records = bindingdb_get(uniprot, float(args.affinity_cutoff_nm))
        kept = 0
        for rec in records:
            affinity_type = clean_text(rec.get("affinity_type"))
            if affinity_type not in allowed_types:
                continue
            affinity_nm, relation = parse_bindingdb_affinity(rec.get("affinity"))
            if affinity_nm is None:
                continue
            if args.require_equal_relation and relation in {"<", ">"}:
                continue
            smiles = canonical_smiles(clean_text(rec.get("smile")))
            if not smiles:
                continue
            doi = clean_text(rec.get("doi"))
            pmid = clean_text(rec.get("pmid"))
            source_id = doi or (f"PMID:{pmid}" if pmid else "")
            rows.append(
                {
                    "target_chembl_id": clean_text(getattr(target, "target_chembl_id", "")),
                    "target_pref_name": clean_text(getattr(target, "target_pref_name", uniprot)),
                    "uniprot_accession": uniprot,
                    "target_sequence": clean_text(getattr(target, "target_sequence", "")),
                    "standard_type": affinity_type,
                    "bindingdb_query": clean_text(rec.get("query")),
                    "bindingdb_monomerid": clean_text(rec.get("monomerid")),
                    "canonical_smiles": smiles,
                    "pchembl_value": float(-math.log10(affinity_nm * 1e-9)),
                    "standard_relation": relation,
                    "standard_value": float(affinity_nm),
                    "standard_units": "nM",
                    "doi": doi,
                    "pmid": pmid,
                    "source_id": source_id,
                }
            )
            kept += 1
        print(f"fetched BindingDB target={uniprot} raw={len(records)} kept={kept}", flush=True)
        time.sleep(float(args.request_sleep))
    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("No BindingDB activities fetched")
    df.to_csv(out, index=False)
    return df


def aggregate_activities(raw: pd.DataFrame, max_mols_per_target: int, seed: int) -> pd.DataFrame:
    raw = raw.copy()
    raw["pchembl_value"] = pd.to_numeric(raw["pchembl_value"], errors="coerce")
    raw = raw[np.isfinite(raw["pchembl_value"])]
    raw["target"] = raw["uniprot_accession"].astype(str) + "::" + raw["standard_type"].astype(str)
    rows = []
    rng = np.random.default_rng(seed)
    for (target, smiles), group in raw.groupby(["target", "canonical_smiles"], sort=True):
        first = group.iloc[0]
        sources = sorted({clean_text(x) for x in group.get("source_id", pd.Series(dtype=str)) if clean_text(x)})
        rows.append(
            {
                "target": str(target),
                "protein_target": str(first.get("uniprot_accession", "")),
                "target_pref_name": clean_text(first.get("target_pref_name", target)),
                "target_sequence": clean_text(first.get("target_sequence", "")),
                "smiles": str(smiles),
                "pchembl": float(group["pchembl_value"].median()),
                "activity_count": int(len(group)),
                "document_ids": "|".join(sources[:50]),
                "assay_ids": "",
                "n_documents": int(len(sources)),
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
                row = {
                    "target": str(target),
                    "protein_target": str(group.loc[i, "protein_target"]),
                    "target_pref_name": str(group.loc[i, "target_pref_name"]),
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
                    "assay_source_key": "",
                    "pair_year": -1,
                    "n_document_sources": int(len(doc_tokens)),
                    "n_assay_sources": 0,
                }
                (pos_rows if label else neg_rows).append(row)
        if float(args.neg_pos_ratio) > 0 and pos_rows:
            max_neg = int(np.ceil(len(pos_rows) * float(args.neg_pos_ratio)))
            if len(neg_rows) > max_neg:
                keep = rng.choice(len(neg_rows), size=max_neg, replace=False)
                neg_rows = [neg_rows[int(k)] for k in keep]
        rows.extend(pos_rows)
        rows.extend(neg_rows)
        print(f"pairs BindingDB target={target} mols={len(group)} pos={len(pos_rows)} neg={len(neg_rows)}", flush=True)
    pairs = pd.DataFrame(rows)
    if pairs.empty:
        raise RuntimeError("No BindingDB pair rows constructed")
    pairs = pairs.sample(frac=1.0, random_state=int(args.downsample_seed)).reset_index(drop=True)
    Path(args.pair_output).parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(args.pair_output, index=False)
    return pairs


def run_tasks(args: argparse.Namespace, records: pd.DataFrame) -> pd.DataFrame:
    diff_bits, scalars = build_feature_matrices(records)
    variants = parse_csv_arg(args.variants)
    split_modes = parse_csv_arg(args.split_modes)
    seeds = [int(x) for x in parse_csv_arg(args.seeds)]
    rows = []
    for seed in seeds:
        masks, target_to_family = family_folds(
            records,
            seed,
            int(args.folds),
            int(args.sequence_hash_features),
            int(args.target_family_clusters),
        )
        family_arr = records["target"].astype(str).map(target_to_family).to_numpy(dtype=np.int32)
        for fold in range(min(int(args.folds), len(masks))):
            for split_mode in split_modes:
                if split_mode.startswith("temporal"):
                    continue
                test_mask = masks[fold]
                base_train_mask = ~test_mask
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
                        "benchmark": "BindingDB_raw_MoleculeACE_targets",
                        "split_mode": split_mode,
                        "variant": variant,
                        "seed": seed,
                        "fold": fold,
                        "temporal_cutoff_year": -1,
                        "n_train": int(train_mask.sum()),
                        "n_test": int(test_mask.sum()),
                        "train_retention_after_purge": float(retention),
                        "train_positive_rate": float(y_train.mean()),
                        "test_positive_rate": float(y_test.mean()),
                        "n_train_targets": int(records.loc[train_mask, "target"].nunique()),
                        "n_test_targets": int(records.loc[test_mask, "target"].nunique()),
                        "n_train_protein_targets": int(records.loc[train_mask, "protein_target"].nunique()) if "protein_target" in records else -1,
                        "n_test_protein_targets": int(records.loc[test_mask, "protein_target"].nunique()) if "protein_target" in records else -1,
                        "n_train_families": int(len(train_families)),
                        "n_test_families": int(len(test_families)),
                        "family_overlap_rate": float(len(train_families & test_families) / max(len(test_families), 1)),
                    }
                    row.update(audit_overlaps(records, train_mask, test_mask))
                    row.update(row_metrics(y_test, prob))
                    rows.append(row)
                    print(
                        f"[BindingDB {split_mode} seed={seed} fold={fold} {variant}] "
                        f"roc={row['roc_auc']:.4f} pr={row['pr_auc']:.4f} "
                        f"train={row['n_train']} test={row['n_test']} "
                        f"docov={row['document_source_overlap_rate']:.3f}",
                        flush=True,
                    )
    return pd.DataFrame(rows)


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
        "train_retention_after_purge_mean",
        "roc_auc_count",
    ]
    compact = summary[[c for c in compact_cols if c in summary.columns]]
    delta_focus = deltas[(deltas["baseline"] == "ecfp_absdiff_xgb") & (deltas["metric"].isin(["roc_auc", "pr_auc", "f1", "mcc"]))].copy()
    if not delta_focus.empty:
        delta_focus = delta_focus.sort_values(["split_mode", "metric", "delta_mean"], ascending=[True, True, False])
    lines = [
        "# Raw BindingDB Source-Purged Target-Family Pair-Cliff Benchmark",
        "",
        "## Scope",
        "",
        "This benchmark queries BindingDB raw affinity records for MoleculeACE/ChEMBL protein targets through UniProt accessions, keeps PMID/DOI source identifiers, constructs high-similarity pair-cliff labels, and evaluates target-family plus scaffold/source purging.",
        "",
        "## Dataset",
        "",
        f"- raw BindingDB activity rows: {len(raw):,}",
        f"- constructed pair rows: {len(pairs):,}",
        f"- protein targets: {pairs['protein_target'].nunique():,}",
        f"- target/type tasks: {pairs['target'].nunique():,}",
        f"- positive pair rate: {pairs['label'].mean():.4f}",
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
        "- This is the second raw-database external check, independent from the raw ChEMBL API collection.",
        "- BindingDB exposes DOI/PMID source metadata but not a stable assay/year field in this endpoint, so the strict variant purges publication source and chemical scaffold rather than temporal source.",
        "- Treat this as reproducibility and database-transfer evidence; ChEMBL remains the stronger block for explicit temporal testing.",
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--moleculeace-csv", default="external/MTPNet/data/MoleculeACE.csv")
    parser.add_argument("--target-metadata-csv", default="data/chembl_raw_moleculeace_target_activities.csv")
    parser.add_argument("--raw-output", default="data/bindingdb_raw_moleculeace_target_activities.csv")
    parser.add_argument("--pair-output", default="data/bindingdb_raw_moleculeace_target_pairs.csv")
    parser.add_argument("--reuse-raw", action="store_true")
    parser.add_argument("--reuse-pairs", action="store_true")
    parser.add_argument("--max-targets", type=int, default=30)
    parser.add_argument("--activity-types", default="Ki,IC50,Kd")
    parser.add_argument("--affinity-cutoff-nm", type=float, default=100000.0)
    parser.add_argument("--request-sleep", type=float, default=0.2)
    parser.add_argument("--require-equal-relation", action="store_true")
    parser.add_argument("--min-mols-per-target", type=int, default=40)
    parser.add_argument("--max-mols-per-target", type=int, default=650)
    parser.add_argument("--tanimoto-threshold", type=float, default=0.7)
    parser.add_argument("--cliff-delta", type=float, default=1.0)
    parser.add_argument("--smooth-delta", type=float, default=0.3)
    parser.add_argument("--neg-pos-ratio", type=float, default=5.0)
    parser.add_argument("--downsample-seed", type=int, default=2026)
    parser.add_argument("--variants", default="ecfp_absdiff_xgb,rich_scalars_only,rich_diff_scalars")
    parser.add_argument("--split-modes", default="target_family,family_scaffold_source_purged")
    parser.add_argument("--seeds", default="0,1,2,3,4")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--target-family-clusters", type=int, default=10)
    parser.add_argument("--sequence-hash-features", type=int, default=128)
    parser.add_argument("--min-train-rows", type=int, default=500)
    parser.add_argument("--min-test-rows", type=int, default=200)
    parser.add_argument("--n-estimators", type=int, default=80)
    parser.add_argument("--max-depth", type=int, default=4)
    parser.add_argument("--model-threads", type=int, default=2)
    parser.add_argument("--output-prefix", default="results/bindingdb_raw_source_pairs_seed0_4_n80")
    args = parser.parse_args()

    raw = fetch_or_load_raw_bindingdb(args)
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
        "protein_targets": int(pairs["protein_target"].nunique()),
        "target_type_tasks": int(pairs["target"].nunique()),
        "positive_rate": float(pairs["label"].mean()),
        "output_prefix": str(prefix),
        "bindingdb_rest": BINDINGDB_REST,
    }
    Path(f"{prefix}_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
