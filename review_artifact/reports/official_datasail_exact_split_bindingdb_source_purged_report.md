# Official DataSAIL Exact Split Audit

This follow-up uses the official DataSAIL 1.3.0 C2 optimizer with ECFP molecular similarity. DataSAIL assigns unique unordered molecular-pair interactions to fold names; the training side is then source-purged using the same document/assay overlap rule as the raw-source benchmark.

## Assignment Scale

| source        |   datasail_run | split        |   n_unique_interactions |
|:--------------|---------------:|:-------------|------------------------:|
| bindingdb_raw |              0 | fold0        |                    1682 |
| bindingdb_raw |              0 | fold1        |                     645 |
| bindingdb_raw |              0 | fold2        |                      23 |
| bindingdb_raw |              0 | fold3        |                    1118 |
| bindingdb_raw |              0 | fold4        |                     299 |
| bindingdb_raw |              0 | not selected |                   18244 |

## Model Summary

| source        | split_mode                              | variant           |   roc_auc_mean |   roc_auc_std |   pr_auc_mean |   pr_auc_std |   mcc_mean |   mcc_std |   n_test_mean |   exact_ligand_overlap_rate_mean |   scaffold_overlap_rate_mean |   document_source_overlap_rate_mean |   assay_source_overlap_rate_mean |   roc_auc_count |
|:--------------|:----------------------------------------|:------------------|---------------:|--------------:|--------------:|-------------:|-----------:|----------:|--------------:|---------------------------------:|-----------------------------:|------------------------------------:|---------------------------------:|----------------:|
| bindingdb_raw | official_datasail_c2_ecfp_source_purged | ecfp_absdiff_xgb  |       0.462669 |     0.0478536 |      0.311705 |     0.137257 | -0.0702223 | 0.0945002 |           952 |                                0 |                   0.00343845 |                                   0 |                                0 |               4 |
| bindingdb_raw | official_datasail_c2_ecfp_source_purged | rich_diff_scalars |       0.609585 |     0.0448624 |      0.428919 |     0.120137 |  0.154027  | 0.0714145 |           952 |                                0 |                   0.00343845 |                                   0 |                                0 |               4 |
| bindingdb_raw | official_datasail_c2_ecfp_source_purged | rich_scalars_only |       0.601298 |     0.0653736 |      0.41593  |     0.116161 |  0.126623  | 0.0725991 |           952 |                                0 |                   0.00343845 |                                   0 |                                0 |               4 |

## Paired Deltas vs ECFP-XGB

| source        | split_mode                              | variant           | baseline         | metric   |   n_pairs |   delta_mean |   delta_std |   wins_favorable |
|:--------------|:----------------------------------------|:------------------|:-----------------|:---------|----------:|-------------:|------------:|-----------------:|
| bindingdb_raw | official_datasail_c2_ecfp_source_purged | rich_diff_scalars | ecfp_absdiff_xgb | roc_auc  |         4 |     0.146916 |   0.0665054 |                4 |
| bindingdb_raw | official_datasail_c2_ecfp_source_purged | rich_diff_scalars | ecfp_absdiff_xgb | pr_auc   |         4 |     0.117214 |   0.0448327 |                4 |
| bindingdb_raw | official_datasail_c2_ecfp_source_purged | rich_diff_scalars | ecfp_absdiff_xgb | mcc      |         4 |     0.224249 |   0.11862   |                4 |
| bindingdb_raw | official_datasail_c2_ecfp_source_purged | rich_scalars_only | ecfp_absdiff_xgb | roc_auc  |         4 |     0.138629 |   0.0643283 |                4 |
| bindingdb_raw | official_datasail_c2_ecfp_source_purged | rich_scalars_only | ecfp_absdiff_xgb | pr_auc   |         4 |     0.104225 |   0.0378072 |                4 |
| bindingdb_raw | official_datasail_c2_ecfp_source_purged | rich_scalars_only | ecfp_absdiff_xgb | mcc      |         4 |     0.196845 |   0.0873013 |                4 |

## Manuscript Boundary

- This closes the environment-level blocker: official DataSAIL is installed and executed in a separate Python 3.12 environment.
- Because DataSAIL C2 optimizes molecule-side cluster assignment rather than the paper's manually purged target/source contract, these rows should be reported as an additional official-split audit, not as a replacement for the frozen main split contract.
- Missing document/assay source identifiers are normalized to empty source sets before source purging, so an absent assay key is not treated as a shared leakage token.
- If a source/fold combination still has too few selected or source-purged rows, it is skipped and the feasibility audit should be reported with the metric table.
