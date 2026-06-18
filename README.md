# PC-CFRL-Bench Review Artifact

This is an anonymized, lightweight review repository for PC-CFRL-Bench, a leakage-audited information contract for pairwise molecular prediction benchmarks under distribution shift.

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

## Full Experiment Reproduction

Full reproduction requires downloading the public source datasets and rebuilding pair rows and split definitions. The selected scripts in `scripts/` document the execution path used for the main benchmark components:

- public activity-cliff anchors: `run_moleculeace_target_family_pairs.py`
- raw ChEMBL/BindingDB reconstruction: `run_chembl_raw_source_temporal_pairs.py`, `run_bindingdb_raw_source_pairs.py`
- splitter stress tests: `run_official_datasail_exact_split.py`, `run_tkde_lohi_leakage_splits.py`
- reliability and auxiliary checks: `run_tkde_raw_reliability.py`, `run_tkde_experiments_5_7.py`
- table and artifact export: `build_tkde_experiment_tables.py`, `export_tkde_leaderboard_package.py`

The repository omits large raw data, pretrained embedding caches, conda environments, logs, and generated intermediate model files to keep the review artifact lightweight.

## Data Sources

The benchmark uses public molecular resources including ACNet, MoleculeACE, ChEMBL, and BindingDB. This repository includes only small examples and frozen result tables. Users should follow the original resource licenses and download terms when rebuilding full data.

## Anonymous Review Note

This repository is prepared for anonymous peer review. It intentionally avoids author names, affiliations, local machine paths, remote server addresses, private credentials, and commit history from the experiment workspace.

