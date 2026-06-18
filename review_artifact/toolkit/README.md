# PC-CFRL-Bench Toolkit

This folder contains a lightweight review artifact for inspecting the PC-CFRL benchmark contract without rerunning the full training workflow. The toolkit is intentionally table-first: it reads frozen CSV files from the anonymized artifact and regenerates contract-level audit summaries.

## Contents

- `pc_cfrl_bench_schema.json`: schema snapshot for leaderboard rows, split-audit rows, and contract-delta rows.
- `run_contract_audit.py`: computes same-contract deltas between PC-CFRL-full (`rich_diff_scalars`) and ECFP-XGB (`ecfp_absdiff_xgb`) from the frozen leaderboard snapshot.

## Quick Start

From the artifact root:

```text
python toolkit/run_contract_audit.py \
  --leaderboard case_studies/tkde_leaderboard_snapshot.csv \
  --out-dir toolkit/generated
```

The command writes:

```text
toolkit/generated/contract_leaderboard_deltas.csv
toolkit/generated/contract_audit_summary.json
```

## Interpretation

The generated deltas are not a new experiment. They are a leaderboard-style view of the same frozen results used by the manuscript. A favorable contrast means that PC-CFRL-full has higher ROC-AUC than ECFP-XGB within the same benchmark, condition, and split contract. Stronger-information controls, target/assay-profile rows, support-only rows, and active-acquisition rows are intentionally reported elsewhere because they answer different information-protocol questions.

## Claim Boundary

PC-CFRL-Bench is an evaluation contract and audit package. It does not claim universal neural SOTA, wet-lab validation, or superiority over methods using additional biological context. Its purpose is to make train-test shortcuts, split semantics, leaderboard rows, and reproducibility files inspectable before a model comparison is interpreted.
