#!/usr/bin/env bash
set -euo pipefail

export PYTHONPATH=scripts:src
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-2}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-2}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-2}
PY=${PY:-./.env/bin/python}

echo "[1/24] ACNet joint target-family/cold-ligand OOD"
$PY scripts/run_acnet_target_family_ablation_joint.py \
  --variants ecfp_absdiff_xgb,rich_pair_only,rich_diff_scalars \
  --split-modes target_family,family_exact_ligand_purged,family_scaffold_purged \
  --seeds 0,1,2,3,4 --folds 5 --n-estimators 10 --workers 4 --model-threads 2 \
  --output results/acnet_target_family_best_joint_seed0_4_allfold_n10_reproduce.jsonl

echo "[2/24] ACNet neural/fusion controls"
$PY scripts/run_acnet_target_family_neural_controls.py \
  --split-modes target_family,family_scaffold_purged \
  --variants ecfp_mlp,rich_diff_scalars_mlp,target_esm_concat_mlp,target_esm_concat_mlp_shuffled \
  --seeds 0,1,2,3,4 --folds 5 --epochs 12 --workers 4 \
  --output results/acnet_target_family_neural_controls_full25_reproduce.jsonl

echo "[3/24] MoleculeACE external target-family pairs"
$PY scripts/run_moleculeace_target_family_pairs.py \
  --variants ecfp_absdiff_xgb,rich_scalars_only,rich_diff_scalars \
  --split-modes target_family,family_scaffold_purged \
  --seeds 0,1,2,3,4 --folds 5 --n-estimators 120 --workers 4 --model-threads 2 \
  --output results/moleculeace_target_family_pairs_seed0_4_n120_reproduce.jsonl

echo "[4/24] ACNet selective prediction"
$PY scripts/run_acnet_selective_prediction.py \
  --split-modes target_family,family_scaffold_purged \
  --variants ecfp_absdiff_xgb,rich_diff_scalars \
  --seeds 0,1,2,3,4 --folds 5 --n-estimators 10 --workers 4 --model-threads 2 \
  --output-prefix results/acnet_selective_prediction_seed0_4_n10_reproduce

echo "[5/24] ACNet conformal prediction"
$PY scripts/run_acnet_conformal_prediction.py \
  --split-modes target_family,family_scaffold_purged \
  --variants ecfp_absdiff_xgb,rich_diff_scalars \
  --seeds 0,1,2,3,4 --folds 5 --n-estimators 10 --workers 4 --model-threads 2 \
  --output-prefix results/acnet_conformal_prediction_seed0_4_n10_reproduce

echo "[6/24] ACNet TreeSHAP-style explanation stability"
$PY scripts/run_acnet_treeshap_stability.py \
  --split-modes target_family,family_scaffold_purged \
  --seeds 0,1,2,3,4 --folds 5 --n-estimators 10 --workers 4 --model-threads 2 \
  --output-prefix results/acnet_treeshap_stability_seed0_4_n10_reproduce

echo "[7/24] Raw ChEMBL source/temporal benchmark"
$PY scripts/run_chembl_raw_source_temporal_pairs.py \
  --max-targets 30 --per-target-limit 1000 --max-mols-per-target 550 \
  --split-modes target_family,family_scaffold_source_purged,temporal_forward,temporal_scaffold_source_purged \
  --seeds 0,1,2,3,4 --folds 5 --n-estimators 40 --model-threads 2 \
  --output-prefix results/chembl_raw_source_temporal_pairs_seed0_4_n40_reproduce

echo "[8/24] Raw BindingDB source-purged benchmark"
$PY scripts/run_bindingdb_raw_source_pairs.py \
  --max-targets 30 --min-mols-per-target 40 --max-mols-per-target 650 \
  --split-modes target_family,family_scaffold_source_purged \
  --seeds 0,1,2,3,4 --folds 5 --n-estimators 80 --model-threads 2 \
  --output-prefix results/bindingdb_raw_source_pairs_seed0_4_n80_reproduce

echo "[9/24] TKDE benchmark split freeze"
$PY scripts/freeze_tkde_benchmark_splits.py \
  --benchmarks acnet,moleculeace,chembl_raw,bindingdb_raw \
  --seeds 0,1,2,3,4 --folds 5 \
  --output-dir splits/tkde_protocol_reproduce \
  --manifest results/tkde_benchmark_split_manifest_reproduce.csv \
  --report results/tkde_benchmark_split_freeze_report_reproduce.md

echo "[10/24] Raw data quality audit"
$PY scripts/analyze_tkde_raw_data_quality.py \
  --output-prefix results/tkde_raw_data_quality_reproduce

echo "[11/24] ChEMBL non-neural baseline stress test"
$PY scripts/run_tkde_non_neural_baselines.py \
  --pair-csv data/chembl_raw_moleculeace_target_pairs.csv --benchmark chembl_raw \
  --split-modes target_family,family_scaffold_source_purged \
  --seeds 0,1,2,3,4 --folds 5 --n-estimators 80 --threads 2 \
  --output-prefix results/tkde_non_neural_baselines_chembl_seed0_4_reproduce

echo "[12/24] BindingDB non-neural baseline stress test"
$PY scripts/run_tkde_non_neural_baselines.py \
  --pair-csv data/bindingdb_raw_moleculeace_target_pairs.csv --benchmark bindingdb_raw \
  --split-modes target_family,family_scaffold_source_purged \
  --seeds 0,1,2,3,4 --folds 5 --n-estimators 80 --threads 2 \
  --output-prefix results/tkde_non_neural_baselines_bindingdb_seed0_4_reproduce

echo "[13/24] Unified raw-source representation-control matrix"
$PY scripts/run_tkde_unified_baseline_matrix.py \
  --sources chembl_raw=data/chembl_raw_moleculeace_target_pairs.csv,bindingdb_raw=data/bindingdb_raw_moleculeace_target_pairs.csv \
  --embedding-dir data/pretrained_embeddings --embeddings chemberta,molformer \
  --variants ecfp_absdiff_xgb,pccfrl_xgb,pccfrl_scalar_xgb,rf_scalars,extratrees_scalars \
  --split-modes random,target_family,family_scaffold_source_purged \
  --seeds 0,1,2,3,4 --folds 5 --target-family-clusters 10 --workers 4 \
  --output-prefix results/tkde_unified_baseline_matrix_reproduce

echo "[14/24] Raw threshold sensitivity"
$PY scripts/run_raw_threshold_sensitivity.py \
  --sources chembl_raw,bindingdb_raw --seeds 0,1,2 --folds 5 \
  --n-estimators 30 --model-threads 2 --max-mols-per-target 500 \
  --pair-prefix data/raw_threshold_sensitivity_reproduce \
  --output-prefix results/raw_threshold_sensitivity_seed0_2_reproduce

echo "[15/24] Stratified raw robustness"
$PY scripts/run_tkde_stratified_raw_eval.py \
  --sources chembl_raw,bindingdb_raw --seeds 0,1,2 --folds 5 \
  --n-estimators 30 --model-threads 2 \
  --output-prefix results/tkde_stratified_raw_eval_seed0_2_reproduce

echo "[16/24] Runtime scalability profile"
$PY scripts/profile_tkde_runtime_scalability.py \
  --sources chembl_raw,bindingdb_raw --fractions 0.25,0.5,1.0 \
  --n-estimators 30 --model-threads 2 \
  --output results/tkde_runtime_scalability_reproduce.csv \
  --report results/tkde_runtime_scalability_reproduce_report.md

echo "[17/24] LoHi/DataSAIL-style chemical block splits"
$PY scripts/run_tkde_lohi_leakage_splits.py \
  --sources chembl_raw,bindingdb_raw \
  --variants ecfp_absdiff_xgb,rich_diff_scalars,rich_scalars_only \
  --seeds 0,1,2 --folds 5 --n-estimators 40 --model-threads 2 \
  --output-prefix results/tkde_lohi_leakage_splits_seed0_2_reproduce

echo "[18/24] Raw external selective and conformal reliability"
$PY scripts/run_tkde_raw_reliability.py \
  --sources chembl_raw,bindingdb_raw \
  --variants ecfp_absdiff_xgb,rich_diff_scalars \
  --seeds 0,1,2 --folds 5 --n-estimators 40 --model-threads 2 \
  --output-prefix results/tkde_raw_reliability_seed0_2_reproduce

echo "[19/24] Raw label-noise and censoring sensitivity"
$PY scripts/run_tkde_noise_sensitivity.py \
  --sources chembl_raw,bindingdb_raw \
  --conditions all_pairs,high_margin,multi_source,high_similarity,equal_relation_rebuilt \
  --variants ecfp_absdiff_xgb,rich_diff_scalars,rich_scalars_only \
  --seeds 0,1,2 --folds 5 --n-estimators 30 --model-threads 2 \
  --output-prefix results/tkde_noise_sensitivity_seed0_2_reproduce

echo "[20/24] Raw success/failure case studies"
$PY scripts/run_tkde_failure_case_studies.py \
  --sources chembl_raw,bindingdb_raw \
  --seeds 0,1,2 --folds 5 --top-k 10 --n-estimators 40 --model-threads 2 \
  --csv-out results/tkde_failure_case_studies_reproduce.csv \
  --report results/tkde_failure_case_studies_reproduce.md

echo "[21/24] Statistical effect-size audit"
$PY scripts/run_tkde_statistical_effects.py \
  --n-boot 5000 --seed 20260518 \
  --csv-out results/tkde_statistical_effects_audit_reproduce.csv \
  --report results/tkde_statistical_effects_audit_reproduce.md

echo "[22/24] TKDE leaderboard package"
$PY scripts/export_tkde_leaderboard_package.py \
  --snapshot results/tkde_leaderboard_snapshot_reproduce.csv \
  --schema docs/tkde_leaderboard_schema_reproduce.md

echo "[23/24] Remaining P0-P1 cross-database, prospective, and ranking utility"
$PY scripts/run_tkde_remaining_p0_p1_experiments.py \
  --seeds 0,1,2 \
  --crossdb-modes matched_targets,exact_ligand_purged,pair_scaffold_purged,scaffold_purged \
  --prospective-modes future_only,scaffold_source_purged \
  --cutoff-years 2005,2010,2015 \
  --variants ecfp_absdiff_xgb,rich_diff_scalars,rich_scalars_only,extratrees_scalars \
  --n-estimators 40 --model-threads 8 \
  --output-prefix results/tkde_remaining_p0_p1_seed0_2_reproduce

echo "[24/24] TKDE / Transactions gap-closure and P0-P1-P2 reports"
$PY scripts/write_tkde_transactions_gap_closure_report.py \
  --output results/tkde_transactions_gap_closure_report_reproduce.md
$PY scripts/write_tkde_p0_p1_p2_completion_report.py \
  --output results/tkde_p0_p1_p2_completion_report_reproduce.md
