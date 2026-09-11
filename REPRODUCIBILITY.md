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

`DATA_MANIFEST.csv` records each public source, access route, local role, and redistribution status. `SOURCE_SNAPSHOT_MANIFEST.csv` binds the cached ChEMBL 36 and BindingDB inputs and derived pair tables to their retrieval date and SHA-256 hashes. ChEMBL 36 was the production release when the unversioned web service was queried on 18 May 2026; ChEMBL 37 was released publicly on 29 May 2026 and is used only by the separately identified historical-release replay. The BindingDB endpoint did not expose a release identifier, so the retrieval record, query specification, and content hash jointly identify that snapshot. The ChEMBL retrieval queries 23 retained MoleculeACE target/endpoint combinations, all Ki in the cached reconstruction, with pagination to exhaustion or a hard cap of 1,000 API activity rows per target; nine targets reach the cap. Pair construction then retains at most 550 ligands per target. The raw-source scripts reconstruct activity rows and pair rows from these declared snapshots. Frozen pair CSVs are intentionally not redistributed in this lightweight repository.

## Leakage and Representation Audits

`review_artifact/tables/raw_audit_15_fold_metrics.csv` contains the matched diagnostic rows for seeds 0--2 and folds 0--4. `raw_audit_15_summary.csv` is the eight-row summary used for the random-to-audited score-drop claim. This dedicated audit configuration is distinct from the canonical 25-run reference evaluation.

`review_artifact/tables/representation_selection_42.csv` contains two sources by three split modes by seven representations. `representation_selection_changes_12.csv` compares the random-split winner with the audited-split winner in two sources by two audited splits by three metrics; 6/12 selections change. `representation_ranking_reversal_12.csv` evaluates all 21 pairwise orderings in each cell; 85/252 orderings reverse, yielding 33.73%.

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

Official DataSAIL v1.3.0 audits use C2/ECFP, molecular entities on both sides, SCIP, five equal target folds, eight threads, and one optimizer run per configuration. Time limits are 120 seconds for 20/40 clusters and 300 seconds for 80 clusters. The 80-cluster status means that no accepted assignment was found within this budget; it is not a proof of mathematical infeasibility. LoHi-style pair partitions group rows by paired Bemis--Murcko scaffold pattern, while ligand partitions group individual ligands by Bemis--Murcko scaffold. Seeds 0--2 are greedily balanced across five folds by positive and total row counts before document-source purging.

## Biological-Unit Bootstrap

After reconstructing the raw pair CSVs, run:

```bash
python scripts/run_target_unit_bootstrap.py \
  --dataset chembl \
  --pair-csv data/chembl_raw_moleculeace_target_pairs.csv \
  --output-prefix results/chembl_target_unit_bootstrap \
  --split-mode family_scaffold_source_purged

python scripts/run_target_unit_bootstrap.py \
  --dataset bindingdb \
  --pair-csv data/bindingdb_raw_moleculeace_target_pairs.csv \
  --output-prefix results/bindingdb_target_unit_bootstrap \
  --split-mode family_scaffold_source_purged
```

These commands implement the strict target-cluster+scaffold+source contract. Repeat each command with `--split-mode target_family` and a distinct output prefix for the target-cluster-only contract. The command-line values retain the experiment code's original identifiers; manuscript-facing names are recorded in `RESULT_MANIFEST.csv`. The script generates paired out-of-fold scores, averages each pair across five seeds, reports every evaluable biological target, and recomputes the ROC-AUC difference in 2,000 bootstrap samples of biological targets and sequence-derived target clusters. `scripts/export_target_cluster_membership.py` exports the exact 23-target membership under both raw-source reconstructions.

## Figure Reproduction

`review_artifact/figures/source/` contains the frozen numerical record behind each data figure, and `review_artifact/figures/pdf/` contains the corresponding publication output. The plotting scripts reconstruct those source records from the full experiment table package and then render the figures. Because the lightweight repository does not redistribute the large raw/intermediate package, run the plotting commands only after Level 2 reconstruction. `RESULT_MANIFEST.csv` gives the input records, preprocessing, fixed configuration, exact command, intermediate source CSV, and final PDF for every main and supplementary figure.

Point the plotting scripts to the reconstructed package before running them:

```bash
export PC_CFRL_DATA_ROOT=/path/to/reconstructed/experiment-data
python scripts/draw_fig2_4.py
python scripts/draw_fig5_7.py
python scripts/draw_sfig1_6.py
```

## Software Environments

`SOFTWARE_ENVIRONMENT.csv` lists exact versions for the main benchmark environment and the separate DataSAIL environment. These records correspond to the frozen reruns summarized in the manuscript and supplement.

## Evidence Map

`RESULT_MANIFEST.csv` is the figure- and table-level execution map. For every manuscript evidence block it identifies inputs, preprocessing, fixed configuration, executable command, intermediate records, and final distributed output. SHA-256 checksums for the distributed files are listed in `checksums_sha256.txt`.
