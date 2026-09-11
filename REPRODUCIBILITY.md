# Reproducibility Guide

PC-CFRL-Bench separates benchmark inspection from full model reruns so that the central contract claim can be checked without downloading large molecular resources.

## Level 1: Contract Audit

Create the environment and run:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python reproduce_quick.py
```

The audit reads the frozen leaderboard and verifies that 39 comparisons satisfy the same-contract rule and that 38 have a positive mean ROC-AUC difference. LoHi rows are included; DataSAIL assignments and stronger-input, adaptation, reliability, and auxiliary-task rows are excluded.

## Level 2: Data Reconstruction

`DATA_MANIFEST.csv` records each public source, access route, local role, and redistribution status. The raw-source scripts reconstruct activity rows and pair rows from ChEMBL and BindingDB. Frozen pair CSVs are intentionally not redistributed in this lightweight repository.

## Level 3: Model Reruns

The principal entry points are:

```bash
python scripts/run_chembl_raw_source_temporal_pairs.py --help
python scripts/run_bindingdb_raw_source_pairs.py --help
python scripts/run_moleculeace_target_family_pairs.py --help
python scripts/run_official_datasail_exact_split.py --help
python scripts/run_tkde_lohi_leakage_splits.py --help
```

`configs/reproduce_trans_top_journal_experiments.sh` documents the complete experiment sequence. Commands write new outputs rather than replacing the frozen tables.

## Biological-Unit Bootstrap

After reconstructing the raw pair CSVs, run:

```bash
python scripts/run_target_unit_bootstrap.py \
  --dataset chembl \
  --pair-csv data/chembl_raw_moleculeace_target_pairs.csv \
  --output-prefix results/chembl_target_unit_bootstrap

python scripts/run_target_unit_bootstrap.py \
  --dataset bindingdb \
  --pair-csv data/bindingdb_raw_moleculeace_target_pairs.csv \
  --output-prefix results/bindingdb_target_unit_bootstrap
```

The script generates paired out-of-fold scores under the family, scaffold, and source-purged contract, averages each pair across five seeds, and recomputes the ROC-AUC difference in 2,000 bootstrap samples of biological targets and sequence-derived target clusters.

## Evidence Map

`RESULT_MANIFEST.csv` links each central manuscript item to its source table, script, comparison protocol, and inspection output. SHA-256 checksums for the distributed files are listed in `checksums_sha256.txt`.
