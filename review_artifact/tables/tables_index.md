# Experiment Table Index

The CSV files are the frozen, machine-readable evidence records. LaTeX exports are retained where the original table-building workflow generated them.

| Table ID | Rows | Evidence |
|---|---:|---|
| `table_01_split_leakage_audit` | 9 | Frozen split leakage audit |
| `table_02_acnet_hard_ood_main` | 9 | ACNet hard-OOD comparison |
| `table_03_raw_source_temporal_main` | 18 | Raw ChEMBL and BindingDB source/temporal comparison |
| `table_04_lohi_datasail_style` | 12 | LoHi and DataSAIL-style splitter stress tests |
| `table_05_crossdb_prospective_topk` | 56 | Cross-database, release-prospective, and ranking utility results |
| `table_06a_selective_prediction` | 42 | Score-thresholded retained-set reliability |
| `table_06b_conformal_prediction` | 24 | Conformal prediction reliability |
| `table_07_noise_censoring_sensitivity` | 30 | Noise, censoring, and exact-relations-only sensitivity |
| `table_08_non_neural_baselines` | 24 | Non-neural learner controls |
| `table_09a_moleculeace_external` | 6 | MoleculeACE external target-cluster controls |
| `table_09b_mtpnet_controls` | 9 | MTPNet split and cold-target controls |
| `table_10a_acnet_neural_controls` | 8 | ACNet neural and shuffled-target controls |
| `table_10b_target_assay_pretraining_controls` | 8 | Target/assay pretraining controls |
| `table_10c_few_pretrained_representations` | 32 | Few-shot pretrained-representation controls |
| `internal_graph_gcn_baseline` | 2 | End-to-end graph-pair GCN stress control |
| `table_11_statistical_effects` | 200 | Paired effect-size and multiplicity audit |
| `table_12_runtime_scalability` | 6 | Runtime and scalability profile |
| `table_13_target_failure_taxonomy` | 2,856 | Target-level failure taxonomy |
| `table_14a_p0_p1_calibration_bins` | 1,680 | Cross-database/prospective calibration bins |
| `table_14b_p0_p1_risk_coverage` | 1,008 | Cross-database/prospective risk-coverage records |
| `table_15_experiment_completion_matrix` | 20 | Experiment coverage matrix |
| `table_16_official_datasail_exact_split` | variable | Official DataSAIL exact-split audit |
| `table_17_official_datasail_multicluster_sensitivity` | variable | DataSAIL cluster-budget feasibility audit |
| `table_18_biological_unit_bootstrap` | 8 | Biological-target and target-cluster bootstrap sensitivity |
| `table_19_target_cluster_membership` | 23 | Complete sequence-derived target-cluster membership |
| `table_20_target_level_distribution` | 23 | Strict-contract target-level ROC-AUC differences in both sources |
| `table_21_target_distribution_summary` | 4 | Biological-target effect distribution summary |

See `RESULT_MANIFEST.csv` at the repository root for manuscript links, inputs, preprocessing, configurations, commands, intermediate records, and final outputs.
