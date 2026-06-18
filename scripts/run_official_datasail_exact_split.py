from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from datasail.sail import datasail

import run_chembl_raw_source_temporal_pairs as raw


def parse_csv(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()]


def stable_ligand_ids(records: pd.DataFrame) -> tuple[dict[str, str], dict[str, str]]:
    ligands = sorted(set(records["smiles1"].astype(str)) | set(records["smiles2"].astype(str)))
    id_to_smiles = {f"m{idx:06d}": smiles for idx, smiles in enumerate(ligands)}
    smiles_to_id = {smiles: mid for mid, smiles in id_to_smiles.items()}
    return smiles_to_id, id_to_smiles


def pair_tuple(row: pd.Series, smiles_to_id: dict[str, str]) -> tuple[str, str]:
    a = smiles_to_id[str(row["smiles1"])]
    b = smiles_to_id[str(row["smiles2"])]
    return tuple(sorted((a, b)))


def run_datasail_split(
    records: pd.DataFrame,
    source: str,
    args: argparse.Namespace,
) -> tuple[list[dict[tuple[str, str], str]], pd.DataFrame]:
    smiles_to_id, id_to_smiles = stable_ligand_ids(records)
    pairs = records.apply(lambda row: pair_tuple(row, smiles_to_id), axis=1)
    interactions = sorted(set(pairs.tolist()))
    e_data = dict(id_to_smiles)
    f_data = dict(id_to_smiles)
    names = [f"fold{i}" for i in range(int(args.folds))]
    split_sizes = [1.0 / int(args.folds)] * int(args.folds)
    output_dir = Path(args.datasail_output_dir) / source
    output_dir.mkdir(parents=True, exist_ok=True)
    _, _, inter_splits = datasail(
        techniques=["C2"],
        inter=interactions,
        output=str(output_dir),
        max_sec=int(args.datasail_max_sec),
        verbose=args.datasail_verbose,
        splits=split_sizes,
        names=names,
        runs=int(args.datasail_runs),
        solver="SCIP",
        cache=bool(args.datasail_cache),
        e_type="M",
        e_data=e_data,
        e_sim="ecfp",
        e_clusters=int(args.datasail_clusters),
        f_type="M",
        f_data=f_data,
        f_sim="ecfp",
        f_clusters=int(args.datasail_clusters),
        threads=int(args.datasail_threads),
    )
    assignments = inter_splits["C2"]
    assignment_rows = []
    for run_idx, mapping in enumerate(assignments):
        for (e_id, f_id), split in mapping.items():
            assignment_rows.append(
                {
                    "source": source,
                    "datasail_run": run_idx,
                    "e_id": e_id,
                    "f_id": f_id,
                    "smiles1": id_to_smiles[e_id],
                    "smiles2": id_to_smiles[f_id],
                    "split": split,
                }
            )
    return assignments, pd.DataFrame(assignment_rows)


def source_purge(records: pd.DataFrame, train_mask: np.ndarray, test_mask: np.ndarray, enabled: bool) -> tuple[np.ndarray, float]:
    if not enabled:
        return train_mask, 1.0
    n_base = int(train_mask.sum())
    train_mask = raw.source_overlap_mask(records, train_mask, test_mask)
    return train_mask, float(train_mask.sum() / max(n_base, 1))


def normalize_source_fields(records: pd.DataFrame) -> pd.DataFrame:
    """Keep missing source identifiers from becoming the literal token 'nan'."""
    records = records.copy()
    for col in ["document_source_key", "assay_source_key"]:
        if col in records.columns:
            records[col] = (
                records[col]
                .fillna("")
                .astype(str)
                .replace({"nan": "", "NaN": "", "None": "", "<NA>": ""})
            )
    return records


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    metrics = raw.METRICS + [
        "n_train",
        "n_test",
        "train_retention_after_source_purge",
        "train_positive_rate",
        "test_positive_rate",
        "target_overlap_rate",
        "exact_ligand_overlap_rate",
        "scaffold_overlap_rate",
        "exact_pair_overlap_rate",
        "pair_scaffold_overlap_rate",
        "document_source_overlap_rate",
        "assay_source_overlap_rate",
    ]
    cols = [c for c in metrics if c in df.columns]
    out = df.groupby(["source", "split_mode", "variant"], as_index=False)[cols].agg(["mean", "std", "count"])
    out.columns = ["_".join([x for x in col if x]).strip("_") for col in out.columns.to_flat_index()]
    return out.reset_index(drop=True)


def paired_deltas(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (source, split_mode), sub in df.groupby(["source", "split_mode"], sort=False):
        base = sub[sub["variant"].eq("ecfp_absdiff_xgb")].set_index(["datasail_run", "fold"])
        for variant in sorted(sub["variant"].astype(str).unique()):
            if variant == "ecfp_absdiff_xgb":
                continue
            cur = sub[sub["variant"].eq(variant)].set_index(["datasail_run", "fold"])
            idx = cur.index.intersection(base.index)
            if len(idx) == 0:
                continue
            for metric in raw.METRICS:
                diff = cur.loc[idx, metric].astype(float) - base.loc[idx, metric].astype(float)
                diff = diff.replace([np.inf, -np.inf], np.nan).dropna()
                if diff.empty:
                    continue
                favorable = diff > 0 if raw.HIGHER_IS_BETTER[metric] else diff < 0
                rows.append(
                    {
                        "source": source,
                        "split_mode": split_mode,
                        "variant": variant,
                        "baseline": "ecfp_absdiff_xgb",
                        "metric": metric,
                        "n_pairs": int(len(diff)),
                        "delta_mean": float(diff.mean()),
                        "delta_std": float(diff.std(ddof=1)) if len(diff) > 1 else 0.0,
                        "wins_favorable": int(favorable.sum()),
                    }
                )
    return pd.DataFrame(rows)


def run_source(args: argparse.Namespace, source: str, pair_csv: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    records = normalize_source_fields(pd.read_csv(pair_csv).reset_index(drop=True))
    if int(args.max_rows) > 0 and len(records) > int(args.max_rows):
        records = records.sample(n=int(args.max_rows), random_state=int(args.sample_seed)).reset_index(drop=True)
    smiles_to_id, _ = stable_ligand_ids(records)
    records["_datasail_pair"] = records.apply(lambda row: pair_tuple(row, smiles_to_id), axis=1)
    assignments, assignment_df = run_datasail_split(records, source, args)
    diff_bits, scalars = raw.build_feature_matrices(records)
    variants = parse_csv(args.variants)
    rows = []
    for run_idx, mapping in enumerate(assignments):
        pair_split = records["_datasail_pair"].map(mapping).astype(str)
        valid_fold_names = [f"fold{i}" for i in range(int(args.folds))]
        selected_mask = pair_split.isin(valid_fold_names).to_numpy(dtype=bool)
        for fold_name in valid_fold_names:
            test_mask = pair_split.eq(fold_name).to_numpy(dtype=bool)
            base_train = selected_mask & ~test_mask
            train_mask, retention = source_purge(records, base_train, test_mask, not bool(args.no_source_purge))
            if train_mask.sum() < int(args.min_train_rows) or test_mask.sum() < int(args.min_test_rows):
                print(
                    f"skip source={source} run={run_idx} fold={fold_name} train={train_mask.sum()} test={test_mask.sum()}",
                    flush=True,
                )
                continue
            y_train = records.loc[train_mask, "label"].to_numpy(dtype=np.int32)
            y_test = records.loc[test_mask, "label"].to_numpy(dtype=np.int32)
            if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
                print(f"skip one-class source={source} run={run_idx} fold={fold_name}", flush=True)
                continue
            fold_id = int(str(fold_name).replace("fold", "")) if str(fold_name).startswith("fold") else len(rows)
            for variant in variants:
                x = raw.build_x(diff_bits, scalars, train_mask, variant)
                model = raw.make_model(
                    seed=int(args.model_seed) + 97 * run_idx + fold_id,
                    y_train=y_train,
                    threads=int(args.model_threads),
                    n_estimators=int(args.n_estimators),
                    max_depth=int(args.max_depth),
                )
                model.fit(x[train_mask], y_train)
                prob = model.predict_proba(x[test_mask])[:, 1].astype(np.float32)
                row = {
                    "source": source,
                    "split_mode": "official_datasail_c2_ecfp_source_purged"
                    if not bool(args.no_source_purge)
                    else "official_datasail_c2_ecfp",
                    "variant": variant,
                    "datasail_run": int(run_idx),
                    "fold": fold_id,
                    "fold_name": str(fold_name),
                    "n_train": int(train_mask.sum()),
                    "n_test": int(test_mask.sum()),
                    "train_retention_after_source_purge": float(retention),
                    "train_positive_rate": float(y_train.mean()),
                    "test_positive_rate": float(y_test.mean()),
                }
                row.update(raw.audit_overlaps(records, train_mask, test_mask))
                row.update(raw.row_metrics(y_test, prob))
                rows.append(row)
                print(
                    f"[DataSAIL {source} run={run_idx} fold={fold_name} {variant}] "
                    f"roc={row['roc_auc']:.4f} pr={row['pr_auc']:.4f} train={row['n_train']} test={row['n_test']} "
                    f"ligov={row['exact_ligand_overlap_rate']:.3f} scafov={row['scaffold_overlap_rate']:.3f}",
                    flush=True,
                )
    return pd.DataFrame(rows), assignment_df


def write_report(summary: pd.DataFrame, deltas: pd.DataFrame, assignment_df: pd.DataFrame, out: Path) -> None:
    focus_cols = [
        "source",
        "split_mode",
        "variant",
        "roc_auc_mean",
        "roc_auc_std",
        "pr_auc_mean",
        "pr_auc_std",
        "mcc_mean",
        "mcc_std",
        "n_test_mean",
        "exact_ligand_overlap_rate_mean",
        "scaffold_overlap_rate_mean",
        "document_source_overlap_rate_mean",
        "assay_source_overlap_rate_mean",
        "roc_auc_count",
    ]
    delta_focus = deltas[deltas["metric"].isin(["roc_auc", "pr_auc", "mcc"])].copy()
    lines = [
        "# Official DataSAIL Exact Split Audit",
        "",
        "This follow-up uses the official DataSAIL 1.3.0 C2 optimizer with ECFP molecular similarity. DataSAIL assigns unique unordered molecular-pair interactions to fold names; the training side is then source-purged using the same document/assay overlap rule as the raw-source benchmark.",
        "",
        "## Assignment Scale",
        "",
        raw.markdown_table(
            assignment_df.groupby(["source", "datasail_run", "split"], as_index=False).size().rename(columns={"size": "n_unique_interactions"})
        ),
        "",
        "## Model Summary",
        "",
        raw.markdown_table(summary[[c for c in focus_cols if c in summary.columns]]),
        "",
        "## Paired Deltas vs ECFP-XGB",
        "",
        raw.markdown_table(delta_focus),
        "",
        "## Manuscript Boundary",
        "",
        "- This closes the environment-level blocker: official DataSAIL is installed and executed in a separate Python 3.12 environment.",
        "- Because DataSAIL C2 optimizes molecule-side cluster assignment rather than the paper's manually purged target/source contract, these rows should be reported as an additional official-split audit, not as a replacement for the frozen main split contract.",
        "- Missing document/assay source identifiers are normalized to empty source sets before source purging, so an absent assay key is not treated as a shared leakage token.",
        "- If a source/fold combination still has too few selected or source-purged rows, it is skipped and the feasibility audit should be reported with the metric table.",
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sources", default="chembl_raw,bindingdb_raw")
    parser.add_argument("--chembl-pairs", default="data/chembl_raw_moleculeace_target_pairs.csv")
    parser.add_argument("--bindingdb-pairs", default="data/bindingdb_raw_moleculeace_target_pairs.csv")
    parser.add_argument("--variants", default="ecfp_absdiff_xgb,rich_diff_scalars,rich_scalars_only")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--datasail-runs", type=int, default=1)
    parser.add_argument("--datasail-clusters", type=int, default=40)
    parser.add_argument("--datasail-max-sec", type=int, default=120)
    parser.add_argument("--datasail-threads", type=int, default=8)
    parser.add_argument("--datasail-cache", action="store_true")
    parser.add_argument("--datasail-verbose", default="W")
    parser.add_argument("--datasail-output-dir", default="results/official_datasail_outputs")
    parser.add_argument("--no-source-purge", action="store_true")
    parser.add_argument("--max-rows", type=int, default=0)
    parser.add_argument("--sample-seed", type=int, default=20260611)
    parser.add_argument("--min-train-rows", type=int, default=500)
    parser.add_argument("--min-test-rows", type=int, default=200)
    parser.add_argument("--model-seed", type=int, default=2026)
    parser.add_argument("--model-threads", type=int, default=4)
    parser.add_argument("--n-estimators", type=int, default=60)
    parser.add_argument("--max-depth", type=int, default=4)
    parser.add_argument("--output-prefix", default="results/official_datasail_exact_split")
    args = parser.parse_args()

    paths = {
        "chembl_raw": args.chembl_pairs,
        "bindingdb_raw": args.bindingdb_pairs,
    }
    frames = []
    assignment_frames = []
    for source in parse_csv(args.sources):
        frame, assignment = run_source(args, source, paths[source])
        if not frame.empty:
            frames.append(frame)
        if not assignment.empty:
            assignment_frames.append(assignment)
    if not frames:
        raise RuntimeError("No official DataSAIL experiment rows were produced")
    results = pd.concat(frames, ignore_index=True)
    assignments = pd.concat(assignment_frames, ignore_index=True) if assignment_frames else pd.DataFrame()
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    results.to_json(f"{prefix}.jsonl", orient="records", lines=True)
    results.to_csv(f"{prefix}_fold_metrics.csv", index=False)
    assignments.to_csv(f"{prefix}_assignments.csv", index=False)
    summary = summarize(results)
    deltas = paired_deltas(results)
    summary.to_csv(f"{prefix}.csv", index=False)
    deltas.to_csv(f"{prefix}_paired_deltas.csv", index=False)
    write_report(summary, deltas, assignments, Path(f"{prefix}_report.md"))
    print(
        json.dumps(
            {
                "rows": int(len(results)),
                "assignments": int(len(assignments)),
                "output_prefix": str(prefix),
                "sources": parse_csv(args.sources),
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
