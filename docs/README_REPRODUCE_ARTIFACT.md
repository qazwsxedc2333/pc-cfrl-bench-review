# PC-CFRL-Bench Anonymous Reproducibility Artifact

This package is an anonymized evidence bundle for the PC-CFRL Transactions manuscript draft. It is designed to let reviewers inspect the frozen split contract, table-generating CSVs, leaderboard snapshots, official DataSAIL audits, toolkit scripts, and manuscript source snapshots without depending on local figure assets.

## Contents

- `tables/`: publication-table CSV/TEX files mirrored from the frozen result package.
- `official_datasail/`: official DataSAIL C2/ECFP summaries, fold metrics, paired deltas, and feasibility audits.
- `reports/`: data card, completion reports, official DataSAIL reports, and solver-boundary logs.
- `case_studies/`: DOI/PubMed-anchored SAR case diagnostics used by Supplementary Table S23.
- `scripts/`: table and official DataSAIL audit scripts used for the final manuscript package.
- `toolkit/`: lightweight PC-CFRL-Bench schema and audit script for regenerating same-contract leaderboard deltas.
- `manuscript/`: the single-file `main.tex` and `supplementary.tex` snapshots used to compile the PDFs.
- `图片/`: PDF figure assets referenced by the manuscript snapshots.
- `checksums_sha256.txt`: SHA-256 hashes for files in this artifact.

## One-Command Table Inspection

The manuscript tables are already materialized. To inspect the current official DataSAIL evidence, open:

```text
official_datasail/official_datasail_exact_split.csv
official_datasail/official_datasail_exact_split_paired_deltas.csv
official_datasail/official_datasail_source_purge_feasibility.csv
tables/table_16_official_datasail_exact_split.csv
tables/table_17_official_datasail_multicluster_sensitivity.csv
case_studies/tkde_sar_literature_case_table.csv
toolkit/pc_cfrl_bench_schema.json
toolkit/run_contract_audit.py
toolkit/generated/contract_leaderboard_deltas.csv
toolkit/generated/contract_audit_summary.json
```

The full training workflow is intended for release through an anonymized public archive or review-only repository. This package contains the frozen table-level artifact used to verify the manuscript during review; public repository and archival DOI information will be inserted after de-anonymization.

## Environment Notes

The main experiment environment used Python 3.13 for the frozen table package. Official DataSAIL requires Python `<3.13`, so a separate Python 3.12 environment was created for the DataSAIL exact-split audit. Missing document/assay source identifiers are normalized to empty source sets before source-purge filtering.

## Claim Boundary

This artifact supports a benchmark-and-reliability contribution. It should not be interpreted as a claim of universal neural SOTA or as experimental biochemical validation. ChEMBL official DataSAIL rows are retained as boundary evidence, and the 80-cluster official DataSAIL configuration is retained as a solver-boundary audit. The toolkit-generated leaderboard deltas are derived from frozen CSV rows and do not constitute additional training runs.
