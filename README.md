# PC-CFRL-Bench

This is the lightweight research repository for PC-CFRL-Bench, a leakage-audited information contract for pairwise molecular activity-cliff benchmarks under distribution shift.

The repository is intentionally table-first and compact. It contains the code needed to inspect benchmark contracts, reproduce contract-level summaries from frozen CSV tables, and run small smoke checks. It does not include full raw ChEMBL, BindingDB, ACNet, or MoleculeACE releases.

## What Is Included

- `src/pccfrl/`: core utility code for fingerprints and metrics.
- `scripts/`: selected experiment, audit, table-building, and plotting scripts from the frozen experiment workspace.
- `review_artifact/tables/`: CSV and LaTeX table sources used for the reported evidence.
- `review_artifact/case_studies/`: leaderboard snapshot, split-audit summary, case-study rows, and SAR provenance notes.
- `review_artifact/official_datasail/`: official DataSAIL audit outputs used as splitter stress tests.
- `review_artifact/toolkit/`: lightweight contract-audit utility.
- `review_artifact/figures/pdf/` and `review_artifact/figures/source/`: exported figures and their frozen numerical source records.
- `data/examples/`: small CSV samples only; full public data should be obtained from the original sources.
- `DATA_MANIFEST.csv` and `SOURCE_SNAPSHOT_MANIFEST.csv`: source access records, retrieval dates, release information, and content hashes.
- `RESULT_MANIFEST.csv`: figure- and table-level inputs, preprocessing, fixed configuration, commands, intermediates, and outputs.
- `SOFTWARE_ENVIRONMENT.csv`: exact package versions for the main benchmark and DataSAIL audit environments.

## Quick Check

Install a minimal Python environment:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Run the lightweight audit:

```bash
python reproduce_quick.py
```

Expected output:

```text
contracts_with_both_variants: 39
favorable_roc_auc: 38
```

The quick check regenerates:

```text
review_artifact/toolkit/generated/contract_leaderboard_deltas.csv
review_artifact/toolkit/generated/contract_audit_summary.json
```

The frozen 39-row set includes same-contract LoHi boundary rows. Official DataSAIL assignments, stronger-input controls, label-assisted adaptation, reliability analyses, and auxiliary tasks are reported separately and are not counted in the 39-row summary.

The headline random-versus-audited diagnostic is independently inspectable in `review_artifact/tables/raw_audit_15_summary.csv` and its 120 run-level rows. The representation-selection audit is exposed as a complete 42-row matrix plus two 12-row summaries. Across the 12 source--audited-split--metric comparisons, the selected representation changes in 6 cases, and 85 of 252 pairwise representation orderings reverse.

## Full Experiment Reproduction

Full reproduction requires downloading the public source datasets and rebuilding pair rows and split definitions. The selected scripts in `scripts/` document the execution path used for the main benchmark components:

- public activity-cliff anchors: `run_moleculeace_target_family_pairs.py`
- raw ChEMBL/BindingDB reconstruction: `run_chembl_raw_source_temporal_pairs.py`, `run_bindingdb_raw_source_pairs.py`
- splitter stress tests: `run_official_datasail_exact_split.py`, `run_tkde_lohi_leakage_splits.py`
- reliability and auxiliary checks: `run_tkde_raw_reliability.py`, `run_tkde_experiments_5_7.py`
- table and artifact export: `build_tkde_experiment_tables.py`, `export_tkde_leaderboard_package.py`
- biological-unit uncertainty: `run_target_unit_bootstrap.py`
- exact target-cluster membership: `export_target_cluster_membership.py`
- publication figures: `draw_fig2_4.py`, `draw_fig5_7.py`, `draw_sfig1_6.py`, `draw_fig1_sfig7_sfig8.py`, and `draw_sfig9_fewshot_label_efficiency.py`

The repository omits large raw data, pretrained embedding caches, conda environments, logs, and generated intermediate model files to remain lightweight. Frozen numerical figure sources remain available for direct inspection. Full figure regeneration follows the commands in `RESULT_MANIFEST.csv` after reconstructing the source table package. `REPRODUCIBILITY.md` distinguishes the quick contract audit from full data reconstruction and model reruns. The frozen ChEMBL reconstruction is a capped snapshot: retrieval stops at exhaustion or 1,000 activity rows per retained target, and pair construction then keeps at most 550 ligands per target.

Set `PC_CFRL_DATA_ROOT` to the reconstructed experiment-data directory before running publication plotting scripts. Neural and transformer controls additionally require the exact PyTorch and Transformers versions recorded in `SOFTWARE_ENVIRONMENT.csv`; these optional heavy dependencies are not installed by the lightweight quick-check requirements.

## Data Sources

The benchmark uses public molecular resources including ACNet, MoleculeACE, ChEMBL, and BindingDB. This repository includes only small examples and frozen result tables. Users should follow the original resource licenses and download terms when rebuilding full data.

## Repository Scope

This repository contains no private credentials, local machine paths, or restricted source data. Full source records must be retrieved from their original public providers under the corresponding licenses and access terms.
