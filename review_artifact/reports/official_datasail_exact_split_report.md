# Official DataSAIL Exact Split Audit

This follow-up uses official DataSAIL 1.3.0 C2 with ECFP molecular similarity, then applies the paper source-purge rule. Missing document/assay source identifiers are normalized to empty source sets so absent assay IDs are not treated as shared leakage tokens. DataSAIL not-selected interactions are excluded from model training and testing.

## Summary
| Source        | Split                                    | Variant             |   Folds |   Train |   Test | ROC-AUC         | PR-AUC          | MCC              | Delta ROC   | Delta PR   | Delta MCC   |   Ligand ov. |   Scaffold ov. |   Doc ov. |   Assay ov. |
|:--------------|:-----------------------------------------|:--------------------|--------:|--------:|-------:|:----------------|:----------------|:-----------------|:------------|:-----------|:------------|-------------:|---------------:|----------:|------------:|
| ChEMBL raw    | Official DataSAIL C2/ECFP + source purge | ECFP abs-diff XGB   |       5 |    1786 |    462 | 0.598 +/- 0.039 | 0.476 +/- 0.101 | 0.168 +/- 0.079  | --          | --         | --          |            0 |          0.015 |         0 |           0 |
| ChEMBL raw    | Official DataSAIL C2/ECFP + source purge | PC-CFRL rich-diff   |       5 |    1786 |    462 | 0.582 +/- 0.114 | 0.480 +/- 0.132 | 0.071 +/- 0.188  | -0.016      | +0.003     | -0.098      |            0 |          0.015 |         0 |           0 |
| ChEMBL raw    | Official DataSAIL C2/ECFP + source purge | PC-CFRL scalar-only |       5 |    1786 |    462 | 0.569 +/- 0.085 | 0.461 +/- 0.169 | 0.084 +/- 0.165  | -0.029      | -0.015     | -0.084      |            0 |          0.015 |         0 |           0 |
| BindingDB raw | Official DataSAIL C2/ECFP + source purge | ECFP abs-diff XGB   |       4 |    2722 |    952 | 0.463 +/- 0.048 | 0.312 +/- 0.137 | -0.070 +/- 0.095 | --          | --         | --          |            0 |          0.003 |         0 |           0 |
| BindingDB raw | Official DataSAIL C2/ECFP + source purge | PC-CFRL rich-diff   |       4 |    2722 |    952 | 0.610 +/- 0.045 | 0.429 +/- 0.120 | 0.154 +/- 0.071  | +0.147      | +0.117     | +0.224      |            0 |          0.003 |         0 |           0 |
| BindingDB raw | Official DataSAIL C2/ECFP + source purge | PC-CFRL scalar-only |       4 |    2722 |    952 | 0.601 +/- 0.065 | 0.416 +/- 0.116 | 0.127 +/- 0.073  | +0.139      | +0.104     | +0.197      |            0 |          0.003 |         0 |           0 |

## Feasibility
| source        | fold   |   selected_train_before_source_purge |   selected_train_after_source_purge |   selected_test | source_purge_feasible   |
|:--------------|:-------|-------------------------------------:|------------------------------------:|----------------:|:------------------------|
| chembl_raw    | fold0  |                                 1779 |                                1731 |             530 | True                    |
| chembl_raw    | fold1  |                                 1373 |                                1196 |             936 | True                    |
| chembl_raw    | fold2  |                                 2239 |                                2203 |              70 | True                    |
| chembl_raw    | fold3  |                                 1668 |                                1643 |             641 | True                    |
| chembl_raw    | fold4  |                                 2177 |                                2158 |             132 | True                    |
| bindingdb_raw | fold0  |                                 2124 |                                1699 |            1707 | True                    |
| bindingdb_raw | fold1  |                                 3183 |                                3110 |             648 | True                    |
| bindingdb_raw | fold2  |                                 3808 |                                3792 |              23 | False                   |
| bindingdb_raw | fold3  |                                 2681 |                                2575 |            1150 | True                    |
| bindingdb_raw | fold4  |                                 3528 |                                3504 |             303 | True                    |

## Boundary
- ChEMBL source-purged official DataSAIL is retained as a negative/boundary audit: PC-CFRL scalar variants do not beat ECFP on mean ROC/MCC, though rich-diff is close on PR-AUC.
- BindingDB source-purged official DataSAIL remains strongly positive: rich-diff improves ROC-AUC, PR-AUC, and MCC over ECFP across the four evaluable folds; fold2 is skipped because the DataSAIL-selected test fold has only 23 interactions.
- These rows close the previous official-DataSAIL environment blocker but remain supplementary diagnostics, not replacements for the frozen target/source/temporal split contract used for the main claims.
