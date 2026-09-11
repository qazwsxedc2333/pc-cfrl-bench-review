# Data Card And Reproducibility Package

## Repository State

- git head: `45eb6ba`
- git status: follow-up manuscript/documentation edits are present after the consistency audit
- python: `Python 3.13.13`
- DataSAIL audit: separate `.env_datasail_py312` was created because official `datasail` requires Python `<3.13`; DataSAIL 1.3.0 was installed there and used for the official C2/ECFP exact-split audit.

## Core Data Sources

- ACNet generated MMP activity-cliff dataset: `external/ACNet/.../MMP_AC.json`.
- ACNet target metadata and UniProt/ESM2 embeddings: `external/ACNet_data/`.
- MoleculeACE ChEMBL-derived target datasets: `external/MTPNet/data/MoleculeACE.csv`.
- Raw ChEMBL activity records queried by `scripts/run_chembl_raw_source_temporal_pairs.py` from the ChEMBL REST API.
- Raw BindingDB affinity records queried by `scripts/run_bindingdb_raw_source_pairs.py` from the BindingDB REST API.

## Split Integrity Fields

All hard-OOD reports include, when available:

- target overlap rate;
- sequence-derived target-cluster overlap rate;
- exact ligand overlap rate;
- Bemis-Murcko scaffold overlap rate;
- exact pair and pair-scaffold overlap;
- document/assay source overlap for raw ChEMBL;
- DOI/PMID source overlap for raw BindingDB;
- train/test year ranges for temporal splits.

## Main Result Files

- `results/acnet_target_family_best_joint_seed0_4_allfold_n10_summary.csv`: 9 rows x 68 columns
- `results/moleculeace_target_family_pairs_seed0_4_n120_summary.csv`: 6 rows x 68 columns
- `results/acnet_selective_prediction_seed0_4_n10_summary.csv`: 24 rows x 33 columns
- `results/acnet_explanation_stability_seed0_4_n10_category_summary.csv`: 8 rows x 9 columns
- `results/acnet_treeshap_stability_seed0_4_n10_category_summary.csv`: 8 rows x 6 columns
- `results/acnet_conformal_prediction_seed0_4_n10_summary.csv`: 12 rows x 93 columns
- `results/chembl_raw_source_temporal_pairs_seed0_4_n40_summary.csv`: 12 rows x 68 columns
- `results/bindingdb_raw_source_pairs_seed0_4_n80_summary.csv`: 6 rows x 68 columns
- `results/tkde_benchmark_split_manifest.csv`: 225 rows x 29 columns
- `results/tkde_raw_data_quality_summary.csv`: 2 rows x 20 columns
- `results/tkde_non_neural_baselines_chembl_seed0_4_summary.csv`: 12 rows x 54 columns
- `results/tkde_non_neural_baselines_bindingdb_seed0_4_summary.csv`: 12 rows x 54 columns
- `results/raw_threshold_sensitivity_seed0_2_summary.csv`: 30 rows x 74 columns
- `results/tkde_stratified_raw_eval_seed0_2_fold_summary.csv`: 4 rows x 21 columns
- `results/tkde_runtime_scalability.csv`: 6 rows x 17 columns
- `results/tkde_lohi_leakage_splits_seed0_2_summary.csv`: 12 rows x 45 columns
- `results/tkde_raw_reliability_seed0_2_selective_summary.csv`: 18 rows x 64 columns
- `results/tkde_raw_reliability_seed0_2_conformal_summary.csv`: 12 rows x 70 columns
- `results/tkde_noise_sensitivity_seed0_2_summary.csv`: 30 rows x 70 columns
- `results/tkde_failure_case_studies.csv`: 80 rows x 24 columns
- `results/tkde_statistical_effects_audit.csv`: 200 rows x 15 columns
- `results/tkde_leaderboard_snapshot.csv`: 131 rows x 23 columns
- `results/tkde_remaining_p0_p1_seed0_2_summary.csv`: 56 rows x 135 columns
- `results/tkde_remaining_p0_p1_seed0_2_risk_coverage.csv`: 1008 rows x 12 columns
- `results/tkde_remaining_p0_p1_seed0_2_target_failure_taxonomy.csv`: 2856 rows x 36 columns
- `results/official_datasail_exact_split.csv`: official DataSAIL C2/ECFP source-purged summary rows
- `results/official_datasail_exact_split_fold_metrics.csv`: official DataSAIL C2/ECFP fold metrics
- `results/official_datasail_exact_split_paired_deltas.csv`: paired deltas for official DataSAIL C2/ECFP
- `results/official_datasail_source_purge_feasibility.csv`: source-purge feasibility for DataSAIL-selected folds
- `results/official_datasail_multicluster_sensitivity_report.md`: clusters=20/40/80 sensitivity and solver-boundary audit

## One-Command Reproduction

Run `configs/reproduce_trans_top_journal_experiments.sh` from the repository root. The script uses `${PY:-./.env/bin/python}` and writes `_reproduce` result files so it does not overwrite the curated outputs.

The reproduction script covers 23 blocks: core ACNet, neural/fusion controls, MoleculeACE, reliability, explanations, raw ChEMBL, raw BindingDB, split freezing, raw-data auditing, non-neural stress tests, threshold sensitivity, stratified robustness, runtime profiling, LoHi/DataSAIL-style hard splits, raw selective/conformal reliability, noise/censoring sensitivity, qualitative case studies, statistical effect-size auditing, leaderboard packaging, remaining P0-P1 cross-database/prospective/top-k utility experiments, and final TKDE/Transactions reporting. Official DataSAIL exact-split rows are generated separately in the Python 3.12 DataSAIL environment and mirrored into the table package.

## Key Environment Packages

```text
numpy==2.4.4
pandas==3.0.2
pandas-flavor==0.8.1
rdkit==2025.9.3
scikit-fingerprints==2.0.0
scikit-learn==1.8.0
scipy==1.17.1
torch==2.11.0
xgboost==3.2.0
```

## Reporting Notes

- Report all paired tests over seed/fold matched runs.
- Use the internal `family_scaffold_purged` mode for the ACNet target-cluster+scaffold/cold-ligand analysis.
- Use the internal `family_scaffold_source_purged` and `temporal_scaffold_source_purged` modes for raw ChEMBL target-cluster+scaffold+source and same-target temporal+scaffold+source evidence.
- Use BindingDB source-purged results as a database-transfer check because this endpoint has DOI/PMID source identifiers but no stable assay/year fields.
- Use `results/tkde_benchmark_split_manifest.csv` as the frozen benchmark contract for reviewer-visible train/test indices.
- Report non-neural and threshold-sensitivity results as stress tests; they prevent the story from depending on one weak baseline or one arbitrary cliff threshold.
- Keep negative target/assay pretraining results as limitations rather than hiding them.
