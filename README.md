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
- `review_artifact/figures/pdf/`: exported figure PDFs for visual inspection.
- `data/examples/`: small CSV samples only; full public data should be obtained from the original sources.
- `DATA_MANIFEST.csv` and `RESULT_MANIFEST.csv`: source-data and manuscript-evidence maps.

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

## Full Experiment Reproduction

Full reproduction requires downloading the public source datasets and rebuilding pair rows and split definitions. The selected scripts in `scripts/` document the execution path used for the main benchmark components:

- public activity-cliff anchors: `run_moleculeace_target_family_pairs.py`
- raw ChEMBL/BindingDB reconstruction: `run_chembl_raw_source_temporal_pairs.py`, `run_bindingdb_raw_source_pairs.py`
- splitter stress tests: `run_official_datasail_exact_split.py`, `run_tkde_lohi_leakage_splits.py`
- reliability and auxiliary checks: `run_tkde_raw_reliability.py`, `run_tkde_experiments_5_7.py`
- table and artifact export: `build_tkde_experiment_tables.py`, `export_tkde_leaderboard_package.py`
- biological-unit uncertainty: `run_target_unit_bootstrap.py`

The repository omits large raw data, pretrained embedding caches, conda environments, logs, and generated intermediate model files to remain lightweight. `REPRODUCIBILITY.md` distinguishes the quick contract audit from full data reconstruction and model reruns.

## Data Sources

The benchmark uses public molecular resources including ACNet, MoleculeACE, ChEMBL, and BindingDB. This repository includes only small examples and frozen result tables. Users should follow the original resource licenses and download terms when rebuilding full data.

## Repository Scope

This repository contains no private credentials, local machine paths, or restricted source data. Full source records must be retrieved from their original public providers under the corresponding licenses and access terms.
