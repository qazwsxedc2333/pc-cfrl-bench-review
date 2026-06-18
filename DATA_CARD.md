# Data Card

## Data Sources

- ACNet activity-cliff pairs and target metadata.
- MoleculeACE target-family activity-cliff tasks.
- Raw ChEMBL activity records reconstructed into pairwise activity-cliff rows.
- Raw BindingDB affinity records reconstructed into pairwise activity-cliff rows.

## Split Fields

The benchmark contract audits target, target family, exact ligand, scaffold, exact pair, pair scaffold, document source, assay source, and temporal fields when available. Raw ChEMBL includes document and assay source identifiers; raw BindingDB includes document/DOI/PMID-like source identifiers but has no measured assay-source field in the current table package.

## Missing Values

Missing document or assay source identifiers are treated as empty source sets for source-purge logic. They are not converted to literal `nan` tokens.

## Main Evidence Tables

- `tables/table_01_split_leakage_audit.csv`
- `tables/table_02_acnet_hard_ood_main.csv`
- `tables/table_03_raw_source_temporal_main.csv`
- `tables/table_04_lohi_datasail_style.csv`
- `tables/table_05_crossdb_prospective_topk.csv`
- `tables/table_06a_selective_prediction.csv`
- `tables/table_06b_conformal_prediction.csv`
- `tables/table_07_noise_censoring_sensitivity.csv`
- `tables/table_08_non_neural_baselines.csv`
- `tables/table_16_official_datasail_exact_split.csv`
- `tables/table_17_official_datasail_multicluster_sensitivity.csv`

## Benchmark Toolkit Files

- `case_studies/tkde_leaderboard_snapshot.csv`
- `toolkit/pc_cfrl_bench_schema.json`
- `toolkit/run_contract_audit.py`
- `toolkit/generated/contract_leaderboard_deltas.csv`
- `toolkit/generated/contract_audit_summary.json`

## Known Boundaries

ChEMBL hard splits expose negative or weak rows under some LoHi/DataSAIL-style settings. These rows are retained intentionally to bound the claim. The artifact is designed for reproducibility and auditing, not for hiding unfavorable configurations.
