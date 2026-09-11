from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from run_chembl_raw_source_temporal_pairs import family_folds


def cluster_map(records: pd.DataFrame, n_hash: int, n_clusters: int) -> dict[str, int]:
    _, mapping = family_folds(records, 0, 5, n_hash, n_clusters)
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export the exact target-to-cluster memberships used by the raw-source benchmarks."
    )
    parser.add_argument("--chembl-pairs", required=True)
    parser.add_argument("--chembl-activities", required=True)
    parser.add_argument("--bindingdb-pairs", required=True)
    parser.add_argument("--bindingdb-activities", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sequence-hash-features", type=int, default=128)
    parser.add_argument("--target-clusters", type=int, default=10)
    args = parser.parse_args()

    chembl_pairs = pd.read_csv(args.chembl_pairs)
    binding_pairs = pd.read_csv(args.bindingdb_pairs)
    chembl_activities = pd.read_csv(args.chembl_activities)
    binding_activities = pd.read_csv(args.bindingdb_activities)

    chembl_clusters = cluster_map(
        chembl_pairs, args.sequence_hash_features, args.target_clusters
    )
    binding_task_clusters = cluster_map(
        binding_pairs, args.sequence_hash_features, args.target_clusters
    )

    chembl_meta = (
        chembl_activities[
            [
                "target_chembl_id",
                "component_accession",
                "target_pref_name",
                "target_sequence",
            ]
        ]
        .drop_duplicates("target_chembl_id")
        .rename(columns={"component_accession": "uniprot_accession"})
    )
    chembl_meta["chembl_cluster"] = chembl_meta["target_chembl_id"].astype(str).map(
        chembl_clusters
    )
    chembl_meta["sequence_length"] = chembl_meta["target_sequence"].astype(str).str.len()

    binding_task = binding_pairs[["target", "protein_target"]].drop_duplicates()
    binding_task["bindingdb_cluster"] = binding_task["target"].astype(str).map(
        binding_task_clusters
    )
    if binding_task.groupby("protein_target")["bindingdb_cluster"].nunique().max() != 1:
        raise RuntimeError("A BindingDB biological target was assigned to multiple clusters")
    binding_cluster = (
        binding_task.drop_duplicates("protein_target")
        .set_index("protein_target")["bindingdb_cluster"]
        .to_dict()
    )
    binding_id = (
        binding_activities[["target_chembl_id", "uniprot_accession"]]
        .drop_duplicates()
        .set_index("target_chembl_id")["uniprot_accession"]
        .to_dict()
    )
    chembl_meta["uniprot_accession"] = chembl_meta.apply(
        lambda row: str(row["uniprot_accession"])
        if str(row["uniprot_accession"]).strip()
        else str(binding_id.get(str(row["target_chembl_id"]), "")),
        axis=1,
    )
    chembl_meta["bindingdb_cluster"] = chembl_meta["uniprot_accession"].map(
        binding_cluster
    )

    output = chembl_meta[
        [
            "target_chembl_id",
            "uniprot_accession",
            "target_pref_name",
            "sequence_length",
            "chembl_cluster",
            "bindingdb_cluster",
        ]
    ].sort_values("target_chembl_id")
    if len(output) != 23 or output.isna().any().any():
        raise RuntimeError("Expected complete memberships for 23 biological targets")
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(path, index=False)
    print(output.to_string(index=False))


if __name__ == "__main__":
    main()
