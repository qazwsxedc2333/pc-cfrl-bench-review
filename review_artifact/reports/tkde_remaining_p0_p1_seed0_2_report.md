# Remaining TKDE / Transactions P0-P1 Experiments

## Official DataSAIL Attempt

- A separate Python 3.12 DataSAIL environment was used during the follow-up consistency audit.
- The repository `.env` uses Python 3.13, while the official `datasail` package declares `Requires-Python >=3.9,<3.13`; therefore the official audit is isolated in `.env_datasail_py312`.
- Official DataSAIL 1.3.0 C2/ECFP exact-split results are now included in the frozen manuscript package as Supplementary Tables S21--S22 and associated CSV artifacts.
- The runs below cover cross-database transfer, prospective temporal validation, top-k utility, calibration bins, risk-coverage, and target-level failure taxonomy.

## Main Summary

| experiment              | split_mode                                      | variant            |   roc_auc_mean |   pr_auc_mean |    mcc_mean |   ef_1pct_mean |   ef_5pct_mean |   bedroc20_mean |   precision_at_50_mean |   n_train_mean |   n_test_mean |   train_retention_after_purge_mean |   scaffold_overlap_rate_mean |   exact_ligand_overlap_rate_mean |
|:------------------------|:------------------------------------------------|:-------------------|---------------:|--------------:|------------:|---------------:|---------------:|----------------:|-----------------------:|---------------:|--------------:|-----------------------------------:|-----------------------------:|---------------------------------:|
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | ecfp_absdiff_xgb   |       0.55323  |      0.426702 |  0.0607835  |       1.45001  |       1.4662   |        0.507003 |               0.626667 |          19845 |         11885 |                           0.882628 |                    0.492554  |                        0         |
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | extratrees_scalars |       0.642765 |      0.522864 |  0.191046   |       2.17869  |       1.72823  |        0.635713 |               0.853333 |          19845 |         11885 |                           0.882628 |                    0.492554  |                        0         |
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | rich_diff_scalars  |       0.629016 |      0.506984 |  0.175835   |       1.74443  |       1.63696  |        0.604771 |               0.673333 |          19845 |         11885 |                           0.882628 |                    0.492554  |                        0         |
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | rich_scalars_only  |       0.629931 |      0.508372 |  0.182261   |       1.78123  |       1.67377  |        0.609598 |               0.733333 |          19845 |         11885 |                           0.882628 |                    0.492554  |                        0         |
| cross_database_transfer | bindingdb_to_chembl_matched_targets             | ecfp_absdiff_xgb   |       0.584004 |      0.447393 |  0.115651   |       1.72971  |       1.51625  |        0.540013 |               0.713333 |          22484 |         11885 |                           1        |                    0.745562  |                        0.351283  |
| cross_database_transfer | bindingdb_to_chembl_matched_targets             | extratrees_scalars |       0.723234 |      0.616815 |  0.309165   |       2.35534  |       2.12128  |        0.759675 |               0.906667 |          22484 |         11885 |                           1        |                    0.745562  |                        0.351283  |
| cross_database_transfer | bindingdb_to_chembl_matched_targets             | rich_diff_scalars  |       0.661051 |      0.531765 |  0.215984   |       1.78123  |       1.71793  |        0.633393 |               0.673333 |          22484 |         11885 |                           1        |                    0.745562  |                        0.351283  |
| cross_database_transfer | bindingdb_to_chembl_matched_targets             | rich_scalars_only  |       0.657619 |      0.527857 |  0.204968   |       1.79595  |       1.72823  |        0.629938 |               0.573333 |          22484 |         11885 |                           1        |                    0.745562  |                        0.351283  |
| cross_database_transfer | bindingdb_to_chembl_pair_scaffold_purged        | ecfp_absdiff_xgb   |       0.558284 |      0.427325 |  0.0662947  |       1.42057  |       1.45148  |        0.506885 |               0.426667 |          19307 |         11885 |                           0.8587   |                    0.402945  |                        0.11199   |
| cross_database_transfer | bindingdb_to_chembl_pair_scaffold_purged        | extratrees_scalars |       0.644502 |      0.51242  |  0.191448   |       2.04621  |       1.62372  |        0.603304 |               0.706667 |          19307 |         11885 |                           0.8587   |                    0.402945  |                        0.11199   |
| cross_database_transfer | bindingdb_to_chembl_pair_scaffold_purged        | rich_diff_scalars  |       0.640539 |      0.50802  |  0.195635   |       1.77387  |       1.57219  |        0.591325 |               0.613333 |          19307 |         11885 |                           0.8587   |                    0.402945  |                        0.11199   |
| cross_database_transfer | bindingdb_to_chembl_pair_scaffold_purged        | rich_scalars_only  |       0.639262 |      0.509806 |  0.187225   |       1.71498  |       1.63844  |        0.601152 |               0.626667 |          19307 |         11885 |                           0.8587   |                    0.402945  |                        0.11199   |
| cross_database_transfer | bindingdb_to_chembl_scaffold_purged             | ecfp_absdiff_xgb   |       0.540888 |      0.417402 |  0.0476586  |       1.32488  |       1.38229  |        0.486489 |               0.406667 |          18079 |         11885 |                           0.804083 |                    0         |                        0         |
| cross_database_transfer | bindingdb_to_chembl_scaffold_purged             | extratrees_scalars |       0.616483 |      0.492803 |  0.161902   |       1.95788  |       1.60311  |        0.585102 |               0.746667 |          18079 |         11885 |                           0.804083 |                    0         |                        0         |
| cross_database_transfer | bindingdb_to_chembl_scaffold_purged             | rich_diff_scalars  |       0.623328 |      0.498374 |  0.174511   |       1.72971  |       1.59133  |        0.584304 |               0.726667 |          18079 |         11885 |                           0.804083 |                    0         |                        0         |
| cross_database_transfer | bindingdb_to_chembl_scaffold_purged             | rich_scalars_only  |       0.619815 |      0.494834 |  0.166318   |       1.64138  |       1.60311  |        0.583783 |               0.62     |          18079 |         11885 |                           0.804083 |                    0         |                        0         |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | ecfp_absdiff_xgb   |       0.541676 |      0.345401 |  0.0429841  |       1.72896  |       1.35664  |        0.407508 |               0.733333 |           7710 |         22484 |                           0.648717 |                    0.144369  |                        0         |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | extratrees_scalars |       0.598133 |      0.411708 |  0.121346   |       2.4158   |       1.76117  |        0.515606 |               0.793333 |           7710 |         22484 |                           0.648717 |                    0.144369  |                        0         |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | rich_diff_scalars  |       0.592556 |      0.395481 |  0.106244   |       1.9658   |       1.63801  |        0.476824 |               0.82     |           7710 |         22484 |                           0.648717 |                    0.144369  |                        0         |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | rich_scalars_only  |       0.597233 |      0.398853 |  0.118662   |       2.07001  |       1.6579   |        0.482639 |               0.7      |           7710 |         22484 |                           0.648717 |                    0.144369  |                        0         |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | ecfp_absdiff_xgb   |       0.550252 |      0.357009 |  0.0535164  |       1.93738  |       1.47222  |        0.436619 |               0.74     |          11885 |         22484 |                           1        |                    0.195917  |                        0.117372  |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | extratrees_scalars |       0.632817 |      0.462944 |  0.169071   |       2.75686  |       2.16948  |        0.613079 |               0.833333 |          11885 |         22484 |                           1        |                    0.195917  |                        0.117372  |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | rich_diff_scalars  |       0.604947 |      0.40638  |  0.137494   |       2.28317  |       1.48264  |        0.475155 |               0.793333 |          11885 |         22484 |                           1        |                    0.195917  |                        0.117372  |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | rich_scalars_only  |       0.607686 |      0.411936 |  0.147337   |       2.48212  |       1.55464  |        0.485033 |               0.793333 |          11885 |         22484 |                           1        |                    0.195917  |                        0.117372  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | ecfp_absdiff_xgb   |       0.532171 |      0.338496 |  0.0409101  |       1.67211  |       1.20032  |        0.371376 |               0.7      |           5933 |         22484 |                           0.499201 |                    0.143702  |                        0.063245  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | extratrees_scalars |       0.588533 |      0.394198 |  0.105626   |       1.8758   |       1.5859   |        0.473656 |               0.78     |           5933 |         22484 |                           0.499201 |                    0.143702  |                        0.063245  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | rich_diff_scalars  |       0.58468  |      0.389687 |  0.103109   |       1.99896  |       1.55274  |        0.465901 |               0.713333 |           5933 |         22484 |                           0.499201 |                    0.143702  |                        0.063245  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | rich_scalars_only  |       0.59766  |      0.400778 |  0.128027   |       1.79527  |       1.62759  |        0.485301 |               0.646667 |           5933 |         22484 |                           0.499201 |                    0.143702  |                        0.063245  |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | ecfp_absdiff_xgb   |       0.525583 |      0.326024 |  0.0397825  |       1.01843  |       0.979584 |        0.327985 |               0.426667 |           3024 |         22484 |                           0.254438 |                    0         |                        0         |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | extratrees_scalars |       0.558754 |      0.366563 |  0.0806363  |       1.56317  |       1.48264  |        0.435311 |               0.493333 |           3024 |         22484 |                           0.254438 |                    0         |                        0         |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | rich_diff_scalars  |       0.584553 |      0.388731 |  0.11516    |       1.54895  |       1.58022  |        0.465681 |               0.52     |           3024 |         22484 |                           0.254438 |                    0         |                        0         |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | rich_scalars_only  |       0.578782 |      0.385222 |  0.105706   |       1.56317  |       1.57927  |        0.466187 |               0.48     |           3024 |         22484 |                           0.254438 |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2005_future_only            | ecfp_absdiff_xgb   |       0.519614 |      0.300058 |  0.0474595  |       1.00817  |       0.79579  |        0.258845 |               0.286667 |           6975 |          4503 |                           1        |                    0.0952698 |                        0.044859  |
| prospective_temporal    | chembl_future_after_2005_future_only            | extratrees_scalars |       0.620184 |      0.43131  |  0.175371   |       2.70486  |       1.99198  |        0.553588 |               0.813333 |           6975 |          4503 |                           1        |                    0.0952698 |                        0.044859  |
| prospective_temporal    | chembl_future_after_2005_future_only            | rich_diff_scalars  |       0.63039  |      0.422305 |  0.16051    |       2.1393   |       2.00199  |        0.532481 |               0.626667 |           6975 |          4503 |                           1        |                    0.0952698 |                        0.044859  |
| prospective_temporal    | chembl_future_after_2005_future_only            | rich_scalars_only  |       0.624968 |      0.416812 |  0.16189    |       2.23766  |       1.98197  |        0.523068 |               0.66     |           6975 |          4503 |                           1        |                    0.0952698 |                        0.044859  |
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | ecfp_absdiff_xgb   |       0.508848 |      0.298755 | -0.0266868  |       0.762278 |       0.670666 |        0.253837 |               0.213333 |           5561 |          4503 |                           0.797276 |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | extratrees_scalars |       0.596428 |      0.368116 |  0.12962    |       1.32784  |       1.39138  |        0.410246 |               0.386667 |           5561 |          4503 |                           0.797276 |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | rich_diff_scalars  |       0.638134 |      0.389075 |  0.13751    |       1.49997  |       1.43643  |        0.417953 |               0.433333 |           5561 |          4503 |                           0.797276 |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | rich_scalars_only  |       0.62639  |      0.375653 |  0.110172   |       1.52456  |       1.40139  |        0.404803 |               0.426667 |           5561 |          4503 |                           0.797276 |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2010_future_only            | ecfp_absdiff_xgb   |       0.497342 |      0.284736 |  0.0238299  |       2.23476  |       1.30961  |        0.315079 |               0.506667 |           9244 |          2234 |                           1        |                    0.0425246 |                        0.0170098 |
| prospective_temporal    | chembl_future_after_2010_future_only            | extratrees_scalars |       0.579446 |      0.328188 |  0.11318    |       0.926609 |       1.43274  |        0.37125  |               0.313333 |           9244 |          2234 |                           1        |                    0.0425246 |                        0.0170098 |
| prospective_temporal    | chembl_future_after_2010_future_only            | rich_diff_scalars  |       0.611429 |      0.346908 |  0.107077   |       1.41717  |       1.48871  |        0.385948 |               0.386667 |           9244 |          2234 |                           1        |                    0.0425246 |                        0.0170098 |
| prospective_temporal    | chembl_future_after_2010_future_only            | rich_scalars_only  |       0.624324 |      0.355776 |  0.125857   |       1.09013  |       1.57825  |        0.397101 |               0.393333 |           9244 |          2234 |                           1        |                    0.0425246 |                        0.0170098 |
| prospective_temporal    | chembl_future_after_2010_scaffold_source_purged | ecfp_absdiff_xgb   |       0.486263 |      0.282532 | -0.00201578 |       2.23476  |       1.35439  |        0.317497 |               0.5      |           8891 |          2234 |                           0.961813 |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2010_scaffold_source_purged | extratrees_scalars |       0.591521 |      0.332252 |  0.107761   |       0.708583 |       1.39916  |        0.363391 |               0.273333 |           8891 |          2234 |                           0.961813 |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2010_scaffold_source_purged | rich_diff_scalars  |       0.59557  |      0.340865 |  0.10453    |       1.63519  |       1.57825  |        0.393455 |               0.406667 |           8891 |          2234 |                           0.961813 |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2010_scaffold_source_purged | rich_scalars_only  |       0.587631 |      0.334927 |  0.107443   |       1.03562  |       1.56706  |        0.391094 |               0.38     |           8891 |          2234 |                           0.961813 |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2015_future_only            | ecfp_absdiff_xgb   |       0.553132 |      0.264061 |  0.0103781  |       6.02381  |       2.77721  |        0.392297 |               0.46     |          10374 |          1104 |                           1        |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2015_future_only            | extratrees_scalars |       0.552709 |      0.188874 |  0.0554977  |       1.09524  |       1.56463  |        0.225828 |               0.24     |          10374 |          1104 |                           1        |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2015_future_only            | rich_diff_scalars  |       0.651696 |      0.302058 |  0.227734   |       3.46825  |       3.44218  |        0.447366 |               0.513333 |          10374 |          1104 |                           1        |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2015_future_only            | rich_scalars_only  |       0.695952 |      0.336694 |  0.232366   |       2.92063  |       3.2466   |        0.466986 |               0.486667 |          10374 |          1104 |                           1        |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2015_scaffold_source_purged | ecfp_absdiff_xgb   |       0.553132 |      0.264061 |  0.0103781  |       6.02381  |       2.77721  |        0.392297 |               0.46     |          10374 |          1104 |                           1        |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2015_scaffold_source_purged | extratrees_scalars |       0.552709 |      0.188874 |  0.0554977  |       1.09524  |       1.56463  |        0.225828 |               0.24     |          10374 |          1104 |                           1        |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2015_scaffold_source_purged | rich_diff_scalars  |       0.651696 |      0.302058 |  0.227734   |       3.46825  |       3.44218  |        0.447366 |               0.513333 |          10374 |          1104 |                           1        |                    0         |                        0         |
| prospective_temporal    | chembl_future_after_2015_scaffold_source_purged | rich_scalars_only  |       0.695952 |      0.336694 |  0.232366   |       2.92063  |       3.2466   |        0.466986 |               0.486667 |          10374 |          1104 |                           1        |                    0         |                        0         |

## Paired Deltas vs ECFP

| split_mode                               | variant            | baseline         | metric          |   n |   delta_mean |   wins |
|:-----------------------------------------|:-------------------|:-----------------|:----------------|----:|-------------:|-------:|
| chembl_to_bindingdb_matched_targets      | extratrees_scalars | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0825642 |      3 |
| chembl_to_bindingdb_matched_targets      | extratrees_scalars | ecfp_absdiff_xgb | pr_auc          |   3 |    0.105935  |      3 |
| chembl_to_bindingdb_matched_targets      | extratrees_scalars | ecfp_absdiff_xgb | mcc             |   3 |    0.115555  |      3 |
| chembl_to_bindingdb_matched_targets      | extratrees_scalars | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.819478  |      3 |
| chembl_to_bindingdb_matched_targets      | extratrees_scalars | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.697266  |      3 |
| chembl_to_bindingdb_matched_targets      | extratrees_scalars | ecfp_absdiff_xgb | bedroc20        |   3 |    0.17646   |      3 |
| chembl_to_bindingdb_matched_targets      | extratrees_scalars | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.0933333 |      3 |
| chembl_to_bindingdb_matched_targets      | rich_diff_scalars  | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0546949 |      3 |
| chembl_to_bindingdb_matched_targets      | rich_diff_scalars  | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0493706 |      3 |
| chembl_to_bindingdb_matched_targets      | rich_diff_scalars  | ecfp_absdiff_xgb | mcc             |   3 |    0.0839777 |      3 |
| chembl_to_bindingdb_matched_targets      | rich_diff_scalars  | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.345791  |      3 |
| chembl_to_bindingdb_matched_targets      | rich_diff_scalars  | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.0104211 |      2 |
| chembl_to_bindingdb_matched_targets      | rich_diff_scalars  | ecfp_absdiff_xgb | bedroc20        |   3 |    0.0385361 |      3 |
| chembl_to_bindingdb_matched_targets      | rich_diff_scalars  | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.0533333 |      2 |
| chembl_to_bindingdb_matched_targets      | rich_scalars_only  | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0574339 |      3 |
| chembl_to_bindingdb_matched_targets      | rich_scalars_only  | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0549271 |      3 |
| chembl_to_bindingdb_matched_targets      | rich_scalars_only  | ecfp_absdiff_xgb | mcc             |   3 |    0.0938205 |      3 |
| chembl_to_bindingdb_matched_targets      | rich_scalars_only  | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.544739  |      3 |
| chembl_to_bindingdb_matched_targets      | rich_scalars_only  | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.0824214 |      3 |
| chembl_to_bindingdb_matched_targets      | rich_scalars_only  | ecfp_absdiff_xgb | bedroc20        |   3 |    0.0484139 |      3 |
| chembl_to_bindingdb_matched_targets      | rich_scalars_only  | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.0533333 |      2 |
| chembl_to_bindingdb_exact_ligand_purged  | extratrees_scalars | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0564569 |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | extratrees_scalars | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0663064 |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | extratrees_scalars | ecfp_absdiff_xgb | mcc             |   3 |    0.0783624 |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | extratrees_scalars | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.686845  |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | extratrees_scalars | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.404528  |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | extratrees_scalars | ecfp_absdiff_xgb | bedroc20        |   3 |    0.108098  |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | extratrees_scalars | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.06      |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | rich_diff_scalars  | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0508798 |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | rich_diff_scalars  | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0500797 |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | rich_diff_scalars  | ecfp_absdiff_xgb | mcc             |   3 |    0.0632599 |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | rich_diff_scalars  | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.236843  |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | rich_diff_scalars  | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.28137   |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | rich_diff_scalars  | ecfp_absdiff_xgb | bedroc20        |   3 |    0.0693164 |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | rich_diff_scalars  | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.0866667 |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | rich_scalars_only  | ecfp_absdiff_xgb | roc_auc         |   3 |    0.055556  |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | rich_scalars_only  | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0534511 |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | rich_scalars_only  | ecfp_absdiff_xgb | mcc             |   3 |    0.0756781 |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | rich_scalars_only  | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.341054  |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | rich_scalars_only  | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.301265  |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | rich_scalars_only  | ecfp_absdiff_xgb | bedroc20        |   3 |    0.0751308 |      3 |
| chembl_to_bindingdb_exact_ligand_purged  | rich_scalars_only  | ecfp_absdiff_xgb | precision_at_50 |   3 |   -0.0333333 |      0 |
| chembl_to_bindingdb_pair_scaffold_purged | extratrees_scalars | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0563613 |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | extratrees_scalars | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0557021 |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | extratrees_scalars | ecfp_absdiff_xgb | mcc             |   3 |    0.0647161 |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | extratrees_scalars | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.203685  |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | extratrees_scalars | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.385581  |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | extratrees_scalars | ecfp_absdiff_xgb | bedroc20        |   3 |    0.10228   |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | extratrees_scalars | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.08      |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | rich_diff_scalars  | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0525085 |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | rich_diff_scalars  | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0511909 |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | rich_diff_scalars  | ecfp_absdiff_xgb | mcc             |   3 |    0.0621993 |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | rich_diff_scalars  | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.326844  |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | rich_diff_scalars  | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.352423  |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | rich_diff_scalars  | ecfp_absdiff_xgb | bedroc20        |   3 |    0.0945244 |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | rich_diff_scalars  | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.0133333 |      2 |
| chembl_to_bindingdb_pair_scaffold_purged | rich_scalars_only  | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0654888 |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | rich_scalars_only  | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0622814 |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | rich_scalars_only  | ecfp_absdiff_xgb | mcc             |   3 |    0.0871168 |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | rich_scalars_only  | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.123158  |      2 |
| chembl_to_bindingdb_pair_scaffold_purged | rich_scalars_only  | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.427265  |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | rich_scalars_only  | ecfp_absdiff_xgb | bedroc20        |   3 |    0.113924  |      3 |
| chembl_to_bindingdb_pair_scaffold_purged | rich_scalars_only  | ecfp_absdiff_xgb | precision_at_50 |   3 |   -0.0533333 |      0 |
| chembl_to_bindingdb_scaffold_purged      | extratrees_scalars | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0331717 |      3 |
| chembl_to_bindingdb_scaffold_purged      | extratrees_scalars | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0405384 |      3 |
| chembl_to_bindingdb_scaffold_purged      | extratrees_scalars | ecfp_absdiff_xgb | mcc             |   3 |    0.0408538 |      3 |
| chembl_to_bindingdb_scaffold_purged      | extratrees_scalars | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.544739  |      3 |
| chembl_to_bindingdb_scaffold_purged      | extratrees_scalars | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.503055  |      3 |
| chembl_to_bindingdb_scaffold_purged      | extratrees_scalars | ecfp_absdiff_xgb | bedroc20        |   3 |    0.107326  |      3 |
| chembl_to_bindingdb_scaffold_purged      | extratrees_scalars | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.0666667 |      3 |
| chembl_to_bindingdb_scaffold_purged      | rich_diff_scalars  | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0589705 |      3 |
| chembl_to_bindingdb_scaffold_purged      | rich_diff_scalars  | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0627067 |      3 |
| chembl_to_bindingdb_scaffold_purged      | rich_diff_scalars  | ecfp_absdiff_xgb | mcc             |   3 |    0.0753771 |      3 |
| chembl_to_bindingdb_scaffold_purged      | rich_diff_scalars  | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.530529  |      3 |
| chembl_to_bindingdb_scaffold_purged      | rich_diff_scalars  | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.600634  |      3 |
| chembl_to_bindingdb_scaffold_purged      | rich_diff_scalars  | ecfp_absdiff_xgb | bedroc20        |   3 |    0.137696  |      3 |
| chembl_to_bindingdb_scaffold_purged      | rich_diff_scalars  | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.0933333 |      2 |
| chembl_to_bindingdb_scaffold_purged      | rich_scalars_only  | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0531996 |      3 |
| chembl_to_bindingdb_scaffold_purged      | rich_scalars_only  | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0591977 |      3 |
| chembl_to_bindingdb_scaffold_purged      | rich_scalars_only  | ecfp_absdiff_xgb | mcc             |   3 |    0.065924  |      3 |
| chembl_to_bindingdb_scaffold_purged      | rich_scalars_only  | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.544739  |      3 |
| chembl_to_bindingdb_scaffold_purged      | rich_scalars_only  | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.599687  |      3 |
| chembl_to_bindingdb_scaffold_purged      | rich_scalars_only  | ecfp_absdiff_xgb | bedroc20        |   3 |    0.138202  |      3 |
| chembl_to_bindingdb_scaffold_purged      | rich_scalars_only  | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.0533333 |      1 |
| bindingdb_to_chembl_matched_targets      | extratrees_scalars | ecfp_absdiff_xgb | roc_auc         |   3 |    0.13923   |      3 |
| bindingdb_to_chembl_matched_targets      | extratrees_scalars | ecfp_absdiff_xgb | pr_auc          |   3 |    0.169421  |      3 |
| bindingdb_to_chembl_matched_targets      | extratrees_scalars | ecfp_absdiff_xgb | mcc             |   3 |    0.193513  |      3 |
| bindingdb_to_chembl_matched_targets      | extratrees_scalars | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.625638  |      3 |
| bindingdb_to_chembl_matched_targets      | extratrees_scalars | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.605029  |      3 |
| bindingdb_to_chembl_matched_targets      | extratrees_scalars | ecfp_absdiff_xgb | bedroc20        |   3 |    0.219662  |      3 |
| bindingdb_to_chembl_matched_targets      | extratrees_scalars | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.193333  |      3 |
| bindingdb_to_chembl_matched_targets      | rich_diff_scalars  | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0770462 |      3 |
| bindingdb_to_chembl_matched_targets      | rich_diff_scalars  | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0843716 |      3 |
| bindingdb_to_chembl_matched_targets      | rich_diff_scalars  | ecfp_absdiff_xgb | mcc             |   3 |    0.100332  |      3 |
| bindingdb_to_chembl_matched_targets      | rich_diff_scalars  | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.0515232 |      1 |
| bindingdb_to_chembl_matched_targets      | rich_diff_scalars  | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.201676  |      3 |
| bindingdb_to_chembl_matched_targets      | rich_diff_scalars  | ecfp_absdiff_xgb | bedroc20        |   3 |    0.0933803 |      3 |
| bindingdb_to_chembl_matched_targets      | rich_diff_scalars  | ecfp_absdiff_xgb | precision_at_50 |   3 |   -0.04      |      1 |
| bindingdb_to_chembl_matched_targets      | rich_scalars_only  | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0736143 |      3 |
| bindingdb_to_chembl_matched_targets      | rich_scalars_only  | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0804637 |      3 |
| bindingdb_to_chembl_matched_targets      | rich_scalars_only  | ecfp_absdiff_xgb | mcc             |   3 |    0.0893162 |      3 |
| bindingdb_to_chembl_matched_targets      | rich_scalars_only  | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.0662441 |      1 |
| bindingdb_to_chembl_matched_targets      | rich_scalars_only  | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.211981  |      3 |
| bindingdb_to_chembl_matched_targets      | rich_scalars_only  | ecfp_absdiff_xgb | bedroc20        |   3 |    0.0899245 |      3 |
| bindingdb_to_chembl_matched_targets      | rich_scalars_only  | ecfp_absdiff_xgb | precision_at_50 |   3 |   -0.14      |      0 |
| bindingdb_to_chembl_exact_ligand_purged  | extratrees_scalars | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0895354 |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | extratrees_scalars | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0961619 |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | extratrees_scalars | ecfp_absdiff_xgb | mcc             |   3 |    0.130263  |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | extratrees_scalars | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.728685  |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | extratrees_scalars | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.262032  |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | extratrees_scalars | ecfp_absdiff_xgb | bedroc20        |   3 |    0.12871   |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | extratrees_scalars | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.226667  |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | rich_diff_scalars  | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0757866 |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | rich_diff_scalars  | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0802819 |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | rich_diff_scalars  | ecfp_absdiff_xgb | mcc             |   3 |    0.115052  |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | rich_diff_scalars  | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.294418  |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | rich_diff_scalars  | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.170762  |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | rich_diff_scalars  | ecfp_absdiff_xgb | bedroc20        |   3 |    0.0977674 |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | rich_diff_scalars  | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.0466667 |      2 |
| bindingdb_to_chembl_exact_ligand_purged  | rich_scalars_only  | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0767016 |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | rich_scalars_only  | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0816704 |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | rich_scalars_only  | ecfp_absdiff_xgb | mcc             |   3 |    0.121478  |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | rich_scalars_only  | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.33122   |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | rich_scalars_only  | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.207565  |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | rich_scalars_only  | ecfp_absdiff_xgb | bedroc20        |   3 |    0.102594  |      3 |
| bindingdb_to_chembl_exact_ligand_purged  | rich_scalars_only  | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.106667  |      2 |
| bindingdb_to_chembl_pair_scaffold_purged | extratrees_scalars | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0862176 |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | extratrees_scalars | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0850948 |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | extratrees_scalars | ecfp_absdiff_xgb | mcc             |   3 |    0.125153  |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | extratrees_scalars | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.625638  |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | extratrees_scalars | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.172235  |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | extratrees_scalars | ecfp_absdiff_xgb | bedroc20        |   3 |    0.096419  |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | extratrees_scalars | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.28      |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | rich_diff_scalars  | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0822545 |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | rich_diff_scalars  | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0806951 |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | rich_diff_scalars  | ecfp_absdiff_xgb | mcc             |   3 |    0.12934   |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | rich_diff_scalars  | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.353302  |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | rich_diff_scalars  | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.120711  |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | rich_diff_scalars  | ecfp_absdiff_xgb | bedroc20        |   3 |    0.0844403 |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | rich_diff_scalars  | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.186667  |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | rich_scalars_only  | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0809777 |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | rich_scalars_only  | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0824809 |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | rich_scalars_only  | ecfp_absdiff_xgb | mcc             |   3 |    0.12093   |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | rich_scalars_only  | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.294418  |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | rich_scalars_only  | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.186955  |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | rich_scalars_only  | ecfp_absdiff_xgb | bedroc20        |   3 |    0.094267  |      3 |
| bindingdb_to_chembl_pair_scaffold_purged | rich_scalars_only  | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.2       |      3 |
| bindingdb_to_chembl_scaffold_purged      | extratrees_scalars | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0755952 |      3 |
| bindingdb_to_chembl_scaffold_purged      | extratrees_scalars | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0754012 |      3 |
| bindingdb_to_chembl_scaffold_purged      | extratrees_scalars | ecfp_absdiff_xgb | mcc             |   3 |    0.114243  |      3 |
| bindingdb_to_chembl_scaffold_purged      | extratrees_scalars | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.632999  |      3 |
| bindingdb_to_chembl_scaffold_purged      | extratrees_scalars | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.220814  |      3 |
| bindingdb_to_chembl_scaffold_purged      | extratrees_scalars | ecfp_absdiff_xgb | bedroc20        |   3 |    0.0986128 |      3 |
| bindingdb_to_chembl_scaffold_purged      | extratrees_scalars | ecfp_absdiff_xgb | precision_at_50 |   3 |    0.34      |      3 |
| bindingdb_to_chembl_scaffold_purged      | rich_diff_scalars  | ecfp_absdiff_xgb | roc_auc         |   3 |    0.0824394 |      3 |
| bindingdb_to_chembl_scaffold_purged      | rich_diff_scalars  | ecfp_absdiff_xgb | pr_auc          |   3 |    0.0809723 |      3 |
| bindingdb_to_chembl_scaffold_purged      | rich_diff_scalars  | ecfp_absdiff_xgb | mcc             |   3 |    0.126852  |      3 |
| bindingdb_to_chembl_scaffold_purged      | rich_diff_scalars  | ecfp_absdiff_xgb | ef_1pct         |   3 |    0.404825  |      3 |
| bindingdb_to_chembl_scaffold_purged      | rich_diff_scalars  | ecfp_absdiff_xgb | ef_5pct         |   3 |    0.209037  |      3 |
| bindingdb_to_chembl_scaffold_purged      | rich_diff_scalars  | ecfp_absdiff_xgb | bedroc20        |   3 |    0.0978145 |      3 |

## Target-Level Failure Taxonomy Preview

| experiment              | split_mode                                      | variant          | target_accession   | target_display                        |   n_test |   positive_rate |   roc_auc |   pr_auc |   ef_5pct |   bedroc20 |
|:------------------------|:------------------------------------------------|:-----------------|:-------------------|:--------------------------------------|---------:|----------------:|----------:|---------:|----------:|-----------:|
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | ecfp_absdiff_xgb | P35367             | Histamine H1 receptor                 |      396 |        0.161616 |  0.440865 | 0.142228 |  0        |  0.0212593 |
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | ecfp_absdiff_xgb | P35367             | Histamine H1 receptor                 |      396 |        0.161616 |  0.455572 | 0.142864 |  0        |  0.0256535 |
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | ecfp_absdiff_xgb | P35367             | Histamine H1 receptor                 |      396 |        0.161616 |  0.45289  | 0.146206 |  0.309375 |  0.0402326 |
| prospective_temporal    | chembl_future_after_2005_future_only            | ecfp_absdiff_xgb | P35367             | Histamine H1 receptor                 |      396 |        0.161616 |  0.45489  | 0.175374 |  0.61875  |  0.111336  |
| prospective_temporal    | chembl_future_after_2005_future_only            | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1128 |        0.166667 |  0.527982 | 0.178641 |  0.736842 |  0.168847  |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | ecfp_absdiff_xgb | P49841             | Glycogen synthase kinase-3 beta       |     1515 |        0.184158 |  0.497241 | 0.179539 |  0.500141 |  0.139451  |
| prospective_temporal    | chembl_future_after_2005_future_only            | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1128 |        0.166667 |  0.533284 | 0.179719 |  0.947368 |  0.155631  |
| prospective_temporal    | chembl_future_after_2005_future_only            | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1128 |        0.166667 |  0.525942 | 0.179952 |  0.526316 |  0.157876  |
| prospective_temporal    | chembl_future_after_2005_future_only            | ecfp_absdiff_xgb | P35367             | Histamine H1 receptor                 |      396 |        0.161616 |  0.471974 | 0.180153 |  0.61875  |  0.11809   |
| prospective_temporal    | chembl_future_after_2005_future_only            | ecfp_absdiff_xgb | P35367             | Histamine H1 receptor                 |      396 |        0.161616 |  0.474068 | 0.181648 |  0.61875  |  0.114382  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | ecfp_absdiff_xgb | P49841             | Glycogen synthase kinase-3 beta       |     1515 |        0.184158 |  0.445533 | 0.181996 |  0.857385 |  0.182726  |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | ecfp_absdiff_xgb | P49841             | Glycogen synthase kinase-3 beta       |     1515 |        0.184158 |  0.500254 | 0.183905 |  0.57159  |  0.150772  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | ecfp_absdiff_xgb | O60674             | Tyrosine-protein kinase JAK2          |     2270 |        0.177974 |  0.532898 | 0.185589 |  0.985756 |  0.176196  |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | ecfp_absdiff_xgb | P49841             | Glycogen synthase kinase-3 beta       |     1515 |        0.184158 |  0.506452 | 0.186482 |  0.785937 |  0.16611   |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | ecfp_absdiff_xgb | O60674             | Tyrosine-protein kinase JAK2          |     2270 |        0.177974 |  0.541568 | 0.193227 |  1.2322   |  0.214718  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | ecfp_absdiff_xgb | O60674             | Tyrosine-protein kinase JAK2          |     2270 |        0.177974 |  0.543581 | 0.194271 |  1.18291  |  0.220268  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | ecfp_absdiff_xgb | P49841             | Glycogen synthase kinase-3 beta       |     1515 |        0.184158 |  0.476629 | 0.194631 |  0.928834 |  0.200145  |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | ecfp_absdiff_xgb | O60674             | Tyrosine-protein kinase JAK2          |     2270 |        0.177974 |  0.565927 | 0.196163 |  0.69003  |  0.145846  |
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1128 |        0.166667 |  0.557532 | 0.198072 |  1.78947  |  0.233774  |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | ecfp_absdiff_xgb | O60674             | Tyrosine-protein kinase JAK2          |     2270 |        0.177974 |  0.574813 | 0.198861 |  0.936469 |  0.180933  |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | ecfp_absdiff_xgb | O60674             | Tyrosine-protein kinase JAK2          |     2270 |        0.177974 |  0.56892  | 0.199986 |  0.985756 |  0.170199  |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | ecfp_absdiff_xgb | O60674             | Tyrosine-protein kinase JAK2          |     2270 |        0.177974 |  0.539069 | 0.200386 |  1.38006  |  0.225602  |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | ecfp_absdiff_xgb | O60674             | Tyrosine-protein kinase JAK2          |     2270 |        0.177974 |  0.545156 | 0.201923 |  1.2322   |  0.216954  |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | ecfp_absdiff_xgb | O60674             | Tyrosine-protein kinase JAK2          |     2270 |        0.177974 |  0.546252 | 0.204181 |  1.13362  |  0.205953  |
| cross_database_transfer | bindingdb_to_chembl_scaffold_purged             | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1128 |        0.166667 |  0.548715 | 0.206571 |  1.68421  |  0.276871  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | ecfp_absdiff_xgb | P49841             | Glycogen synthase kinase-3 beta       |     1515 |        0.184158 |  0.473737 | 0.206756 |  1.14318  |  0.236132  |
| prospective_temporal    | chembl_future_after_2010_scaffold_source_purged | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1128 |        0.166667 |  0.493108 | 0.209632 |  2.31579  |  0.300101  |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | ecfp_absdiff_xgb | P49841             | Glycogen synthase kinase-3 beta       |     1515 |        0.184158 |  0.53523  | 0.210458 |  1.78622  |  0.313929  |
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1128 |        0.166667 |  0.582792 | 0.214049 |  1.78947  |  0.240306  |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | ecfp_absdiff_xgb | P49841             | Glycogen synthase kinase-3 beta       |     1515 |        0.184158 |  0.52742  | 0.216118 |  2.07201  |  0.336553  |
| cross_database_transfer | bindingdb_to_chembl_scaffold_purged             | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1128 |        0.166667 |  0.573353 | 0.216281 |  2.10526  |  0.273318  |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | ecfp_absdiff_xgb | O60674             | Tyrosine-protein kinase JAK2          |     2270 |        0.177974 |  0.542017 | 0.216831 |  1.47863  |  0.278771  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | ecfp_absdiff_xgb | Q01959             | Sodium-dependent dopamine transporter |      528 |        0.223485 |  0.486668 | 0.217541 |  0.828625 |  0.193766  |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | ecfp_absdiff_xgb | P49841             | Glycogen synthase kinase-3 beta       |     1515 |        0.184158 |  0.498744 | 0.217774 |  2.00057  |  0.305286  |
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1128 |        0.166667 |  0.601457 | 0.220391 |  1.78947  |  0.257989  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | ecfp_absdiff_xgb | Q01959             | Sodium-dependent dopamine transporter |      528 |        0.223485 |  0.484363 | 0.220646 |  1.3258   |  0.242106  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | ecfp_absdiff_xgb | Q01959             | Sodium-dependent dopamine transporter |      528 |        0.223485 |  0.48828  | 0.221267 |  1.16008  |  0.249245  |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | ecfp_absdiff_xgb | Q01959             | Sodium-dependent dopamine transporter |      528 |        0.223485 |  0.511875 | 0.225868 |  0.497175 |  0.128366  |
| cross_database_transfer | bindingdb_to_chembl_scaffold_purged             | ecfp_absdiff_xgb | P35367             | Histamine H1 receptor                 |      594 |        0.213805 |  0.486208 | 0.226106 |  1.24724  |  0.227313  |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | ecfp_absdiff_xgb | P49841             | Glycogen synthase kinase-3 beta       |     1515 |        0.184158 |  0.518706 | 0.226403 |  2.07201  |  0.332371  |
| prospective_temporal    | chembl_future_after_2010_future_only            | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1128 |        0.166667 |  0.554555 | 0.226833 |  2.52632  |  0.337634  |
| cross_database_transfer | bindingdb_to_chembl_pair_scaffold_purged        | ecfp_absdiff_xgb | P35367             | Histamine H1 receptor                 |      594 |        0.213805 |  0.500152 | 0.228331 |  1.40315  |  0.243039  |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | ecfp_absdiff_xgb | P49841             | Glycogen synthase kinase-3 beta       |     1515 |        0.184158 |  0.570171 | 0.22901  |  2.07201  |  0.322983  |
| cross_database_transfer | bindingdb_to_chembl_matched_targets             | ecfp_absdiff_xgb | P11309             | Serine/threonine-protein kinase pim-1 |      120 |        0.166667 |  0.615    | 0.23018  |  0        |  0.18178   |
| cross_database_transfer | bindingdb_to_chembl_pair_scaffold_purged        | ecfp_absdiff_xgb | P00519             | Tyrosine-protein kinase ABL1          |      174 |        0.189655 |  0.537825 | 0.230233 |  1.75758  |  0.311412  |
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1128 |        0.166667 |  0.623837 | 0.230952 |  1.47368  |  0.266209  |
| cross_database_transfer | bindingdb_to_chembl_scaffold_purged             | ecfp_absdiff_xgb | P35367             | Histamine H1 receptor                 |      594 |        0.213805 |  0.491519 | 0.231054 |  1.71496  |  0.266747  |
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1128 |        0.166667 |  0.62745  | 0.231835 |  1.57895  |  0.270751  |
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1128 |        0.166667 |  0.61224  | 0.233372 |  1.78947  |  0.285564  |
| cross_database_transfer | bindingdb_to_chembl_pair_scaffold_purged        | ecfp_absdiff_xgb | P35367             | Histamine H1 receptor                 |      594 |        0.213805 |  0.496636 | 0.233488 |  1.24724  |  0.251812  |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | ecfp_absdiff_xgb | Q01959             | Sodium-dependent dopamine transporter |      528 |        0.223485 |  0.532979 | 0.2335   |  0.99435  |  0.218558  |
| prospective_temporal    | chembl_future_after_2010_scaffold_source_purged | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1128 |        0.166667 |  0.525427 | 0.23439  |  2.73684  |  0.359514  |
| cross_database_transfer | bindingdb_to_chembl_scaffold_purged             | ecfp_absdiff_xgb | P35367             | Histamine H1 receptor                 |      594 |        0.213805 |  0.481284 | 0.234724 |  0.935433 |  0.24573   |
| cross_database_transfer | bindingdb_to_chembl_matched_targets             | ecfp_absdiff_xgb | P11309             | Serine/threonine-protein kinase pim-1 |      120 |        0.166667 |  0.63925  | 0.234917 |  0        |  0.140119  |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | ecfp_absdiff_xgb | O60674             | Tyrosine-protein kinase JAK2          |     2270 |        0.177974 |  0.564964 | 0.235284 |  2.07009  |  0.316137  |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | ecfp_absdiff_xgb | O60674             | Tyrosine-protein kinase JAK2          |     2270 |        0.177974 |  0.551832 | 0.236872 |  2.07009  |  0.337072  |
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | ecfp_absdiff_xgb | P35367             | Histamine H1 receptor                 |      594 |        0.213805 |  0.493315 | 0.237081 |  1.09134  |  0.248705  |
| prospective_temporal    | chembl_future_after_2010_future_only            | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1128 |        0.166667 |  0.529838 | 0.237964 |  2.63158  |  0.358271  |
| cross_database_transfer | bindingdb_to_chembl_matched_targets             | ecfp_absdiff_xgb | P00519             | Tyrosine-protein kinase ABL1          |      174 |        0.189655 |  0.443047 | 0.238724 |  2.34343  |  0.373085  |
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | ecfp_absdiff_xgb | Q9HAZ1             | Dual specificity protein kinase CLK4  |       91 |        0.252747 |  0.436701 | 0.239231 |  0        |  0.06898   |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | ecfp_absdiff_xgb | Q01959             | Sodium-dependent dopamine transporter |      528 |        0.223485 |  0.542724 | 0.241321 |  1.3258   |  0.22189   |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | ecfp_absdiff_xgb | P49841             | Glycogen synthase kinase-3 beta       |     1515 |        0.184158 |  0.538652 | 0.24179  |  2.35781  |  0.374614  |
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | ecfp_absdiff_xgb | Q01959             | Sodium-dependent dopamine transporter |      283 |        0.257951 |  0.420483 | 0.242336 |  1.29224  |  0.318723  |
| cross_database_transfer | bindingdb_to_chembl_pair_scaffold_purged        | ecfp_absdiff_xgb | P35367             | Histamine H1 receptor                 |      594 |        0.213805 |  0.504527 | 0.242376 |  1.24724  |  0.269833  |
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | ecfp_absdiff_xgb | P35367             | Histamine H1 receptor                 |      594 |        0.213805 |  0.504923 | 0.242811 |  1.09134  |  0.250252  |
| cross_database_transfer | bindingdb_to_chembl_matched_targets             | ecfp_absdiff_xgb | P35367             | Histamine H1 receptor                 |      594 |        0.213805 |  0.498626 | 0.243104 |  1.40315  |  0.319917  |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | ecfp_absdiff_xgb | Q01959             | Sodium-dependent dopamine transporter |      528 |        0.223485 |  0.538539 | 0.244926 |  1.49153  |  0.3097    |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | ecfp_absdiff_xgb | Q01959             | Sodium-dependent dopamine transporter |      528 |        0.223485 |  0.534208 | 0.245021 |  1.49153  |  0.294239  |
| cross_database_transfer | bindingdb_to_chembl_pair_scaffold_purged        | ecfp_absdiff_xgb | P11309             | Serine/threonine-protein kinase pim-1 |      120 |        0.166667 |  0.61     | 0.245158 |  1        |  0.271778  |
| cross_database_transfer | bindingdb_to_chembl_matched_targets             | ecfp_absdiff_xgb | P11309             | Serine/threonine-protein kinase pim-1 |      120 |        0.166667 |  0.632    | 0.245905 |  0        |  0.185838  |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | ecfp_absdiff_xgb | Q01959             | Sodium-dependent dopamine transporter |      528 |        0.223485 |  0.541794 | 0.246202 |  1.16008  |  0.314734  |
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | ecfp_absdiff_xgb | Q9HAZ1             | Dual specificity protein kinase CLK4  |       91 |        0.252747 |  0.423274 | 0.247153 |  0        |  0.0969956 |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | ecfp_absdiff_xgb | Q01959             | Sodium-dependent dopamine transporter |      528 |        0.223485 |  0.541629 | 0.247258 |  1.16008  |  0.319286  |
| cross_database_transfer | bindingdb_to_chembl_matched_targets             | ecfp_absdiff_xgb | O60674             | Tyrosine-protein kinase JAK2          |      153 |        0.169935 |  0.610236 | 0.247982 |  2.20673  |  0.321489  |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | ecfp_absdiff_xgb | Q01959             | Sodium-dependent dopamine transporter |      528 |        0.223485 |  0.53149  | 0.248839 |  1.49153  |  0.321488  |
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | ecfp_absdiff_xgb | Q01959             | Sodium-dependent dopamine transporter |      283 |        0.257951 |  0.437019 | 0.249235 |  1.29224  |  0.351446  |
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | ecfp_absdiff_xgb | P11309             | Serine/threonine-protein kinase pim-1 |      120 |        0.166667 |  0.57125  | 0.249511 |  2        |  0.358406  |
| prospective_temporal    | chembl_future_after_2015_future_only            | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1098 |        0.150273 |  0.536399 | 0.250043 |  2.41983  |  0.370531  |
| prospective_temporal    | chembl_future_after_2015_scaffold_source_purged | ecfp_absdiff_xgb | P23458             | Tyrosine-protein kinase JAK1          |     1098 |        0.150273 |  0.536399 | 0.250043 |  2.41983  |  0.370531  |
| cross_database_transfer | bindingdb_to_chembl_scaffold_purged             | ecfp_absdiff_xgb | P11309             | Serine/threonine-protein kinase pim-1 |      120 |        0.166667 |  0.59275  | 0.252052 |  1        |  0.28704   |

## Reading

- Cross-database transfer is the strongest remaining external-generalization test because train and test databases are disjoint while target accessions are matched.
- Prospective temporal validation tests whether older ChEMBL records support newer records; the source/scaffold-purged mode is the strictest cell.
- EF/BEDROC/AP@K add drug-discovery ranking utility beyond ROC-AUC and PR-AUC.
- Calibration bins and risk-coverage rows are written separately for publication-ready reliability plots.
