# TKDE / Transactions P0-P1-P2 Completion Report

Date: 2026-05-18

## Completion Matrix

| priority   | experiment                                                         | status   |
|:-----------|:-------------------------------------------------------------------|:---------|
| P0         | LoHi/DataSAIL-style chemical block split                           | done     |
| P0         | Published pretrained representation controls                       | done     |
| P0         | Raw external selective/conformal reliability                       | done     |
| P0         | Cross-database and prospective temporal validation                 | done     |
| P0         | Bootstrap CI, effect size, Holm correction                         | done     |
| P1         | Raw failure/success medicinal chemistry cases                      | done     |
| P1         | Label-noise/censoring sensitivity                                  | done     |
| P1         | Top-k enrichment, risk-coverage, and target-level failure taxonomy | done     |
| P2         | Leaderboard snapshot and schema                                    | done     |

## P0 LoHi/DataSAIL-Style Split

| source        | split_mode                         | variant           |   roc_auc_mean |   pr_auc_mean |   mcc_mean |   scaffold_overlap_rate_mean |   document_source_overlap_rate_mean |
|:--------------|:-----------------------------------|:------------------|---------------:|--------------:|-----------:|-----------------------------:|------------------------------------:|
| bindingdb_raw | lohi_ligand_scaffold_source_purged | ecfp_absdiff_xgb  |       0.517817 |      0.332828 |  0.024348  |                            0 |                                   0 |
| bindingdb_raw | lohi_ligand_scaffold_source_purged | rich_diff_scalars |       0.561178 |      0.371333 |  0.0845165 |                            0 |                                   0 |
| bindingdb_raw | lohi_ligand_scaffold_source_purged | rich_scalars_only |       0.556927 |      0.367278 |  0.0793842 |                            0 |                                   0 |
| bindingdb_raw | lohi_pair_scaffold_source_purged   | ecfp_absdiff_xgb  |       0.517697 |      0.327168 |  0.0308271 |                            0 |                                   0 |
| bindingdb_raw | lohi_pair_scaffold_source_purged   | rich_diff_scalars |       0.55398  |      0.360284 |  0.0738055 |                            0 |                                   0 |
| bindingdb_raw | lohi_pair_scaffold_source_purged   | rich_scalars_only |       0.551501 |      0.357334 |  0.0769539 |                            0 |                                   0 |
| chembl_raw    | lohi_ligand_scaffold_source_purged | ecfp_absdiff_xgb  |       0.531267 |      0.410229 |  0.0349133 |                            0 |                                   0 |
| chembl_raw    | lohi_ligand_scaffold_source_purged | rich_diff_scalars |       0.545287 |      0.437142 |  0.0700262 |                            0 |                                   0 |
| chembl_raw    | lohi_ligand_scaffold_source_purged | rich_scalars_only |       0.550749 |      0.44251  |  0.0700464 |                            0 |                                   0 |
| chembl_raw    | lohi_pair_scaffold_source_purged   | ecfp_absdiff_xgb  |       0.547116 |      0.423277 |  0.0546931 |                            0 |                                   0 |
| chembl_raw    | lohi_pair_scaffold_source_purged   | rich_diff_scalars |       0.542789 |      0.43277  |  0.0483985 |                            0 |                                   0 |
| chembl_raw    | lohi_pair_scaffold_source_purged   | rich_scalars_only |       0.538442 |      0.428139 |  0.0407034 |                            0 |                                   0 |

## P0 Raw Reliability

Selective prediction:

| source        | split_mode                      | variant           |   coverage_target |   precision_mean |   positive_precision_lift_mean |   recall_mean |   f1_mean |
|:--------------|:--------------------------------|:------------------|------------------:|-----------------:|-------------------------------:|--------------:|----------:|
| bindingdb_raw | family_scaffold_source_purged   | ecfp_absdiff_xgb  |               0.1 |         0.343584 |                       1.04036  |      0.104111 |  0.15858  |
| bindingdb_raw | family_scaffold_source_purged   | ecfp_absdiff_xgb  |               0.2 |         0.351778 |                       1.078    |      0.215709 |  0.264751 |
| bindingdb_raw | family_scaffold_source_purged   | ecfp_absdiff_xgb  |               0.4 |         0.337765 |                       1.04469  |      0.418055 |  0.369737 |
| bindingdb_raw | family_scaffold_source_purged   | rich_diff_scalars |               0.1 |         0.435103 |                       1.3813   |      0.138243 |  0.208043 |
| bindingdb_raw | family_scaffold_source_purged   | rich_diff_scalars |               0.2 |         0.418918 |                       1.32758  |      0.26565  |  0.321684 |
| bindingdb_raw | family_scaffold_source_purged   | rich_diff_scalars |               0.4 |         0.384491 |                       1.2056   |      0.482444 |  0.423344 |
| chembl_raw    | family_scaffold_source_purged   | ecfp_absdiff_xgb  |               0.1 |         0.425737 |                       1.1066   |      0.1109   |  0.17434  |
| chembl_raw    | family_scaffold_source_purged   | ecfp_absdiff_xgb  |               0.2 |         0.403025 |                       1.04651  |      0.209478 |  0.272468 |
| chembl_raw    | family_scaffold_source_purged   | ecfp_absdiff_xgb  |               0.4 |         0.409196 |                       1.06281  |      0.425314 |  0.412093 |
| chembl_raw    | family_scaffold_source_purged   | rich_diff_scalars |               0.1 |         0.506415 |                       1.38457  |      0.138751 |  0.215359 |
| chembl_raw    | family_scaffold_source_purged   | rich_diff_scalars |               0.2 |         0.482761 |                       1.30123  |      0.260466 |  0.333671 |
| chembl_raw    | family_scaffold_source_purged   | rich_diff_scalars |               0.4 |         0.454081 |                       1.19612  |      0.478656 |  0.459861 |
| chembl_raw    | temporal_scaffold_source_purged | ecfp_absdiff_xgb  |               0.1 |         0.302393 |                       1.16377  |      0.11656  |  0.167125 |
| chembl_raw    | temporal_scaffold_source_purged | ecfp_absdiff_xgb  |               0.2 |         0.263539 |                       0.993866 |      0.198963 |  0.224989 |
| chembl_raw    | temporal_scaffold_source_purged | ecfp_absdiff_xgb  |               0.4 |         0.270078 |                       1.00785  |      0.403322 |  0.321387 |
| chembl_raw    | temporal_scaffold_source_purged | rich_diff_scalars |               0.1 |         0.399008 |                       1.54467  |      0.154701 |  0.22146  |
| chembl_raw    | temporal_scaffold_source_purged | rich_diff_scalars |               0.2 |         0.390652 |                       1.50337  |      0.300975 |  0.337336 |
| chembl_raw    | temporal_scaffold_source_purged | rich_diff_scalars |               0.4 |         0.348363 |                       1.31804  |      0.527448 |  0.416658 |

Conformal prediction:

| source        | split_mode                      | variant           |   alpha |   coverage_mean |   coverage_gap_mean |   avg_set_size_mean |   singleton_rate_mean |   positive_singleton_precision_mean |
|:--------------|:--------------------------------|:------------------|--------:|----------------:|--------------------:|--------------------:|----------------------:|------------------------------------:|
| bindingdb_raw | family_scaffold_source_purged   | ecfp_absdiff_xgb  |     0.1 |        0.856986 |          -0.0430136 |             1.72802 |              0.271976 |                            0.353735 |
| bindingdb_raw | family_scaffold_source_purged   | ecfp_absdiff_xgb  |     0.2 |        0.783086 |          -0.0169141 |             1.56419 |              0.43581  |                            0.359697 |
| bindingdb_raw | family_scaffold_source_purged   | rich_diff_scalars |     0.1 |        0.861805 |          -0.0381945 |             1.6914  |              0.308603 |                            0.42718  |
| bindingdb_raw | family_scaffold_source_purged   | rich_diff_scalars |     0.2 |        0.721859 |          -0.0781414 |             1.3849  |              0.615103 |                            0.395901 |
| chembl_raw    | family_scaffold_source_purged   | ecfp_absdiff_xgb  |     0.1 |        0.873603 |          -0.0263968 |             1.71025 |              0.289753 |                            0.417352 |
| chembl_raw    | family_scaffold_source_purged   | ecfp_absdiff_xgb  |     0.2 |        0.82761  |           0.0276096 |             1.60572 |              0.394283 |                            0.409539 |
| chembl_raw    | family_scaffold_source_purged   | rich_diff_scalars |     0.1 |        0.88468  |          -0.0153204 |             1.70331 |              0.296688 |                            0.512433 |
| chembl_raw    | family_scaffold_source_purged   | rich_diff_scalars |     0.2 |        0.740641 |          -0.0593589 |             1.35376 |              0.646245 |                            0.482707 |
| chembl_raw    | temporal_scaffold_source_purged | ecfp_absdiff_xgb  |     0.1 |        0.827342 |          -0.0726576 |             1.63785 |              0.362152 |                            0.270246 |
| chembl_raw    | temporal_scaffold_source_purged | ecfp_absdiff_xgb  |     0.2 |        0.729203 |          -0.0707973 |             1.44591 |              0.554088 |                            0.263544 |
| chembl_raw    | temporal_scaffold_source_purged | rich_diff_scalars |     0.1 |        0.925175 |           0.0251755 |             1.64357 |              0.356434 |                            0.406686 |
| chembl_raw    | temporal_scaffold_source_purged | rich_diff_scalars |     0.2 |        0.828674 |           0.028674  |             1.29151 |              0.708493 |                            0.377255 |

## P0/P1 Noise And Censoring Sensitivity

| source        | condition              | split_mode                    | variant           |   roc_auc_mean |   pr_auc_mean |    mcc_mean |   n_test_mean |
|:--------------|:-----------------------|:------------------------------|:------------------|---------------:|--------------:|------------:|--------------:|
| chembl_raw    | all_pairs              | family_scaffold_source_purged | ecfp_absdiff_xgb  |       0.53247  |      0.413368 |  0.0424569  |        2377   |
| chembl_raw    | all_pairs              | family_scaffold_source_purged | rich_diff_scalars |       0.570722 |      0.455945 |  0.113882   |        2377   |
| chembl_raw    | all_pairs              | family_scaffold_source_purged | rich_scalars_only |       0.577772 |      0.462323 |  0.125992   |        2377   |
| chembl_raw    | high_margin            | family_scaffold_source_purged | ecfp_absdiff_xgb  |       0.531453 |      0.411828 |  0.0392392  |        1778.8 |
| chembl_raw    | high_margin            | family_scaffold_source_purged | rich_diff_scalars |       0.56662  |      0.468371 |  0.0915655  |        1778.8 |
| chembl_raw    | high_margin            | family_scaffold_source_purged | rich_scalars_only |       0.579039 |      0.474592 |  0.124577   |        1778.8 |
| chembl_raw    | multi_source           | family_scaffold_source_purged | ecfp_absdiff_xgb  |       0.523618 |      0.591512 | -0.00363115 |         411.6 |
| chembl_raw    | multi_source           | family_scaffold_source_purged | rich_diff_scalars |       0.594268 |      0.651249 |  0.135057   |         411.6 |
| chembl_raw    | multi_source           | family_scaffold_source_purged | rich_scalars_only |       0.587308 |      0.647837 |  0.114172   |         411.6 |
| chembl_raw    | high_similarity        | family_scaffold_source_purged | ecfp_absdiff_xgb  |       0.533487 |      0.367914 | -0.0106892  |         706.2 |
| chembl_raw    | high_similarity        | family_scaffold_source_purged | rich_diff_scalars |       0.55874  |      0.421701 |  0.111957   |         706.2 |
| chembl_raw    | high_similarity        | family_scaffold_source_purged | rich_scalars_only |       0.566893 |      0.418947 |  0.104217   |         706.2 |
| chembl_raw    | equal_relation_rebuilt | family_scaffold_source_purged | ecfp_absdiff_xgb  |       0.526026 |      0.411384 |  0.0750484  |        2041.4 |
| chembl_raw    | equal_relation_rebuilt | family_scaffold_source_purged | rich_diff_scalars |       0.601157 |      0.48037  |  0.150869   |        2041.4 |
| chembl_raw    | equal_relation_rebuilt | family_scaffold_source_purged | rich_scalars_only |       0.591264 |      0.471139 |  0.139998   |        2041.4 |
| bindingdb_raw | all_pairs              | family_scaffold_source_purged | ecfp_absdiff_xgb  |       0.51303  |      0.333908 |  0.0273278  |        4496.8 |
| bindingdb_raw | all_pairs              | family_scaffold_source_purged | rich_diff_scalars |       0.572667 |      0.384719 |  0.0977772  |        4496.8 |
| bindingdb_raw | all_pairs              | family_scaffold_source_purged | rich_scalars_only |       0.569329 |      0.383898 |  0.0927697  |        4496.8 |
| bindingdb_raw | high_margin            | family_scaffold_source_purged | ecfp_absdiff_xgb  |       0.515642 |      0.318795 |  0.0476886  |        3404.4 |
| bindingdb_raw | high_margin            | family_scaffold_source_purged | rich_diff_scalars |       0.566848 |      0.361789 |  0.0885822  |        3404.4 |
| bindingdb_raw | high_margin            | family_scaffold_source_purged | rich_scalars_only |       0.573033 |      0.368011 |  0.0976138  |        3404.4 |
| bindingdb_raw | multi_source           | family_scaffold_source_purged | ecfp_absdiff_xgb  |       0.517542 |      0.354684 |  0.0129939  |        1223   |
| bindingdb_raw | multi_source           | family_scaffold_source_purged | rich_diff_scalars |       0.533267 |      0.373012 |  0.0177483  |        1223   |
| bindingdb_raw | multi_source           | family_scaffold_source_purged | rich_scalars_only |       0.535893 |      0.369884 |  0.0315772  |        1223   |
| bindingdb_raw | high_similarity        | family_scaffold_source_purged | ecfp_absdiff_xgb  |       0.52289  |      0.28504  |  0.0260535  |        1165.6 |
| bindingdb_raw | high_similarity        | family_scaffold_source_purged | rich_diff_scalars |       0.523671 |      0.295788 |  0.0260227  |        1165.6 |
| bindingdb_raw | high_similarity        | family_scaffold_source_purged | rich_scalars_only |       0.533198 |      0.299649 |  0.0278816  |        1165.6 |
| bindingdb_raw | equal_relation_rebuilt | family_scaffold_source_purged | ecfp_absdiff_xgb  |       0.506862 |      0.357894 |  0.023098   |        2875.6 |
| bindingdb_raw | equal_relation_rebuilt | family_scaffold_source_purged | rich_diff_scalars |       0.559349 |      0.416838 |  0.0834494  |        2875.6 |
| bindingdb_raw | equal_relation_rebuilt | family_scaffold_source_purged | rich_scalars_only |       0.559767 |      0.41638  |  0.0782804  |        2875.6 |

## P0/P1 Cross-Database, Prospective, And Ranking Utility

| experiment              | split_mode                                      | variant            |   roc_auc_mean |   pr_auc_mean |    mcc_mean |   ef_1pct_mean |   ef_5pct_mean |   bedroc20_mean |   n_train_mean |   n_test_mean |   scaffold_overlap_rate_mean |
|:------------------------|:------------------------------------------------|:-------------------|---------------:|--------------:|------------:|---------------:|---------------:|----------------:|---------------:|--------------:|-----------------------------:|
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | ecfp_absdiff_xgb   |       0.55323  |      0.426702 |  0.0607835  |       1.45001  |       1.4662   |        0.507003 |          19845 |         11885 |                    0.492554  |
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | extratrees_scalars |       0.642765 |      0.522864 |  0.191046   |       2.17869  |       1.72823  |        0.635713 |          19845 |         11885 |                    0.492554  |
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | rich_diff_scalars  |       0.629016 |      0.506984 |  0.175835   |       1.74443  |       1.63696  |        0.604771 |          19845 |         11885 |                    0.492554  |
| cross_database_transfer | bindingdb_to_chembl_exact_ligand_purged         | rich_scalars_only  |       0.629931 |      0.508372 |  0.182261   |       1.78123  |       1.67377  |        0.609598 |          19845 |         11885 |                    0.492554  |
| cross_database_transfer | bindingdb_to_chembl_matched_targets             | ecfp_absdiff_xgb   |       0.584004 |      0.447393 |  0.115651   |       1.72971  |       1.51625  |        0.540013 |          22484 |         11885 |                    0.745562  |
| cross_database_transfer | bindingdb_to_chembl_matched_targets             | extratrees_scalars |       0.723234 |      0.616815 |  0.309165   |       2.35534  |       2.12128  |        0.759675 |          22484 |         11885 |                    0.745562  |
| cross_database_transfer | bindingdb_to_chembl_matched_targets             | rich_diff_scalars  |       0.661051 |      0.531765 |  0.215984   |       1.78123  |       1.71793  |        0.633393 |          22484 |         11885 |                    0.745562  |
| cross_database_transfer | bindingdb_to_chembl_matched_targets             | rich_scalars_only  |       0.657619 |      0.527857 |  0.204968   |       1.79595  |       1.72823  |        0.629938 |          22484 |         11885 |                    0.745562  |
| cross_database_transfer | bindingdb_to_chembl_pair_scaffold_purged        | ecfp_absdiff_xgb   |       0.558284 |      0.427325 |  0.0662947  |       1.42057  |       1.45148  |        0.506885 |          19307 |         11885 |                    0.402945  |
| cross_database_transfer | bindingdb_to_chembl_pair_scaffold_purged        | extratrees_scalars |       0.644502 |      0.51242  |  0.191448   |       2.04621  |       1.62372  |        0.603304 |          19307 |         11885 |                    0.402945  |
| cross_database_transfer | bindingdb_to_chembl_pair_scaffold_purged        | rich_diff_scalars  |       0.640539 |      0.50802  |  0.195635   |       1.77387  |       1.57219  |        0.591325 |          19307 |         11885 |                    0.402945  |
| cross_database_transfer | bindingdb_to_chembl_pair_scaffold_purged        | rich_scalars_only  |       0.639262 |      0.509806 |  0.187225   |       1.71498  |       1.63844  |        0.601152 |          19307 |         11885 |                    0.402945  |
| cross_database_transfer | bindingdb_to_chembl_scaffold_purged             | ecfp_absdiff_xgb   |       0.540888 |      0.417402 |  0.0476586  |       1.32488  |       1.38229  |        0.486489 |          18079 |         11885 |                    0         |
| cross_database_transfer | bindingdb_to_chembl_scaffold_purged             | extratrees_scalars |       0.616483 |      0.492803 |  0.161902   |       1.95788  |       1.60311  |        0.585102 |          18079 |         11885 |                    0         |
| cross_database_transfer | bindingdb_to_chembl_scaffold_purged             | rich_diff_scalars  |       0.623328 |      0.498374 |  0.174511   |       1.72971  |       1.59133  |        0.584304 |          18079 |         11885 |                    0         |
| cross_database_transfer | bindingdb_to_chembl_scaffold_purged             | rich_scalars_only  |       0.619815 |      0.494834 |  0.166318   |       1.64138  |       1.60311  |        0.583783 |          18079 |         11885 |                    0         |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | ecfp_absdiff_xgb   |       0.541676 |      0.345401 |  0.0429841  |       1.72896  |       1.35664  |        0.407508 |           7710 |         22484 |                    0.144369  |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | extratrees_scalars |       0.598133 |      0.411708 |  0.121346   |       2.4158   |       1.76117  |        0.515606 |           7710 |         22484 |                    0.144369  |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | rich_diff_scalars  |       0.592556 |      0.395481 |  0.106244   |       1.9658   |       1.63801  |        0.476824 |           7710 |         22484 |                    0.144369  |
| cross_database_transfer | chembl_to_bindingdb_exact_ligand_purged         | rich_scalars_only  |       0.597233 |      0.398853 |  0.118662   |       2.07001  |       1.6579   |        0.482639 |           7710 |         22484 |                    0.144369  |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | ecfp_absdiff_xgb   |       0.550252 |      0.357009 |  0.0535164  |       1.93738  |       1.47222  |        0.436619 |          11885 |         22484 |                    0.195917  |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | extratrees_scalars |       0.632817 |      0.462944 |  0.169071   |       2.75686  |       2.16948  |        0.613079 |          11885 |         22484 |                    0.195917  |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | rich_diff_scalars  |       0.604947 |      0.40638  |  0.137494   |       2.28317  |       1.48264  |        0.475155 |          11885 |         22484 |                    0.195917  |
| cross_database_transfer | chembl_to_bindingdb_matched_targets             | rich_scalars_only  |       0.607686 |      0.411936 |  0.147337   |       2.48212  |       1.55464  |        0.485033 |          11885 |         22484 |                    0.195917  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | ecfp_absdiff_xgb   |       0.532171 |      0.338496 |  0.0409101  |       1.67211  |       1.20032  |        0.371376 |           5933 |         22484 |                    0.143702  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | extratrees_scalars |       0.588533 |      0.394198 |  0.105626   |       1.8758   |       1.5859   |        0.473656 |           5933 |         22484 |                    0.143702  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | rich_diff_scalars  |       0.58468  |      0.389687 |  0.103109   |       1.99896  |       1.55274  |        0.465901 |           5933 |         22484 |                    0.143702  |
| cross_database_transfer | chembl_to_bindingdb_pair_scaffold_purged        | rich_scalars_only  |       0.59766  |      0.400778 |  0.128027   |       1.79527  |       1.62759  |        0.485301 |           5933 |         22484 |                    0.143702  |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | ecfp_absdiff_xgb   |       0.525583 |      0.326024 |  0.0397825  |       1.01843  |       0.979584 |        0.327985 |           3024 |         22484 |                    0         |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | extratrees_scalars |       0.558754 |      0.366563 |  0.0806363  |       1.56317  |       1.48264  |        0.435311 |           3024 |         22484 |                    0         |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | rich_diff_scalars  |       0.584553 |      0.388731 |  0.11516    |       1.54895  |       1.58022  |        0.465681 |           3024 |         22484 |                    0         |
| cross_database_transfer | chembl_to_bindingdb_scaffold_purged             | rich_scalars_only  |       0.578782 |      0.385222 |  0.105706   |       1.56317  |       1.57927  |        0.466187 |           3024 |         22484 |                    0         |
| prospective_temporal    | chembl_future_after_2005_future_only            | ecfp_absdiff_xgb   |       0.519614 |      0.300058 |  0.0474595  |       1.00817  |       0.79579  |        0.258845 |           6975 |          4503 |                    0.0952698 |
| prospective_temporal    | chembl_future_after_2005_future_only            | extratrees_scalars |       0.620184 |      0.43131  |  0.175371   |       2.70486  |       1.99198  |        0.553588 |           6975 |          4503 |                    0.0952698 |
| prospective_temporal    | chembl_future_after_2005_future_only            | rich_diff_scalars  |       0.63039  |      0.422305 |  0.16051    |       2.1393   |       2.00199  |        0.532481 |           6975 |          4503 |                    0.0952698 |
| prospective_temporal    | chembl_future_after_2005_future_only            | rich_scalars_only  |       0.624968 |      0.416812 |  0.16189    |       2.23766  |       1.98197  |        0.523068 |           6975 |          4503 |                    0.0952698 |
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | ecfp_absdiff_xgb   |       0.508848 |      0.298755 | -0.0266868  |       0.762278 |       0.670666 |        0.253837 |           5561 |          4503 |                    0         |
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | extratrees_scalars |       0.596428 |      0.368116 |  0.12962    |       1.32784  |       1.39138  |        0.410246 |           5561 |          4503 |                    0         |
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | rich_diff_scalars  |       0.638134 |      0.389075 |  0.13751    |       1.49997  |       1.43643  |        0.417953 |           5561 |          4503 |                    0         |
| prospective_temporal    | chembl_future_after_2005_scaffold_source_purged | rich_scalars_only  |       0.62639  |      0.375653 |  0.110172   |       1.52456  |       1.40139  |        0.404803 |           5561 |          4503 |                    0         |
| prospective_temporal    | chembl_future_after_2010_future_only            | ecfp_absdiff_xgb   |       0.497342 |      0.284736 |  0.0238299  |       2.23476  |       1.30961  |        0.315079 |           9244 |          2234 |                    0.0425246 |
| prospective_temporal    | chembl_future_after_2010_future_only            | extratrees_scalars |       0.579446 |      0.328188 |  0.11318    |       0.926609 |       1.43274  |        0.37125  |           9244 |          2234 |                    0.0425246 |
| prospective_temporal    | chembl_future_after_2010_future_only            | rich_diff_scalars  |       0.611429 |      0.346908 |  0.107077   |       1.41717  |       1.48871  |        0.385948 |           9244 |          2234 |                    0.0425246 |
| prospective_temporal    | chembl_future_after_2010_future_only            | rich_scalars_only  |       0.624324 |      0.355776 |  0.125857   |       1.09013  |       1.57825  |        0.397101 |           9244 |          2234 |                    0.0425246 |
| prospective_temporal    | chembl_future_after_2010_scaffold_source_purged | ecfp_absdiff_xgb   |       0.486263 |      0.282532 | -0.00201578 |       2.23476  |       1.35439  |        0.317497 |           8891 |          2234 |                    0         |
| prospective_temporal    | chembl_future_after_2010_scaffold_source_purged | extratrees_scalars |       0.591521 |      0.332252 |  0.107761   |       0.708583 |       1.39916  |        0.363391 |           8891 |          2234 |                    0         |
| prospective_temporal    | chembl_future_after_2010_scaffold_source_purged | rich_diff_scalars  |       0.59557  |      0.340865 |  0.10453    |       1.63519  |       1.57825  |        0.393455 |           8891 |          2234 |                    0         |
| prospective_temporal    | chembl_future_after_2010_scaffold_source_purged | rich_scalars_only  |       0.587631 |      0.334927 |  0.107443   |       1.03562  |       1.56706  |        0.391094 |           8891 |          2234 |                    0         |
| prospective_temporal    | chembl_future_after_2015_future_only            | ecfp_absdiff_xgb   |       0.553132 |      0.264061 |  0.0103781  |       6.02381  |       2.77721  |        0.392297 |          10374 |          1104 |                    0         |
| prospective_temporal    | chembl_future_after_2015_future_only            | extratrees_scalars |       0.552709 |      0.188874 |  0.0554977  |       1.09524  |       1.56463  |        0.225828 |          10374 |          1104 |                    0         |
| prospective_temporal    | chembl_future_after_2015_future_only            | rich_diff_scalars  |       0.651696 |      0.302058 |  0.227734   |       3.46825  |       3.44218  |        0.447366 |          10374 |          1104 |                    0         |
| prospective_temporal    | chembl_future_after_2015_future_only            | rich_scalars_only  |       0.695952 |      0.336694 |  0.232366   |       2.92063  |       3.2466   |        0.466986 |          10374 |          1104 |                    0         |
| prospective_temporal    | chembl_future_after_2015_scaffold_source_purged | ecfp_absdiff_xgb   |       0.553132 |      0.264061 |  0.0103781  |       6.02381  |       2.77721  |        0.392297 |          10374 |          1104 |                    0         |
| prospective_temporal    | chembl_future_after_2015_scaffold_source_purged | extratrees_scalars |       0.552709 |      0.188874 |  0.0554977  |       1.09524  |       1.56463  |        0.225828 |          10374 |          1104 |                    0         |
| prospective_temporal    | chembl_future_after_2015_scaffold_source_purged | rich_diff_scalars  |       0.651696 |      0.302058 |  0.227734   |       3.46825  |       3.44218  |        0.447366 |          10374 |          1104 |                    0         |
| prospective_temporal    | chembl_future_after_2015_scaffold_source_purged | rich_scalars_only  |       0.695952 |      0.336694 |  0.232366   |       2.92063  |       3.2466   |        0.466986 |          10374 |          1104 |                    0         |

## P0 Statistical Effect Audit

| experiment    | variant           | baseline         | metric   |   n |   mean_delta_favorable |     ci95_low |   ci95_high |   win_rate |   cliffs_delta |   wilcoxon_p | split_mode                         | source        | condition       |      holm_p |
|:--------------|:------------------|:-----------------|:---------|----:|-----------------------:|-------------:|------------:|-----------:|---------------:|-------------:|:-----------------------------------|:--------------|:----------------|------------:|
| raw_chembl    | rich_diff_scalars | ecfp_absdiff_xgb | roc_auc  |  25 |             0.0325339  |  0.0130739   |   0.0538827 |   0.64     |      0.28      |  0.0550705   | target_family                      | nan           | nan             | 1           |
| raw_chembl    | rich_diff_scalars | ecfp_absdiff_xgb | pr_auc   |  25 |             0.0495803  |  0.0241149   |   0.0738538 |   0.8      |      0.6       |  0.00181568  | target_family                      | nan           | nan             | 0.203356    |
| raw_chembl    | rich_diff_scalars | ecfp_absdiff_xgb | mcc      |  25 |             0.055305   |  0.016303    |   0.09714   |   0.64     |      0.28      |  0.0318078   | target_family                      | nan           | nan             | 1           |
| raw_chembl    | rich_scalars_only | ecfp_absdiff_xgb | roc_auc  |  25 |             0.0365652  |  0.0126304   |   0.0596602 |   0.64     |      0.28      |  0.0318078   | target_family                      | nan           | nan             | 1           |
| raw_chembl    | rich_scalars_only | ecfp_absdiff_xgb | pr_auc   |  25 |             0.050978   |  0.0306754   |   0.0732574 |   0.88     |      0.76      |  1.82986e-05 | target_family                      | nan           | nan             | 0.00303757  |
| raw_chembl    | rich_scalars_only | ecfp_absdiff_xgb | mcc      |  25 |             0.0760667  |  0.0361793   |   0.116198  |   0.68     |      0.36      |  0.00557882  | target_family                      | nan           | nan             | 0.591354    |
| raw_chembl    | rich_diff_scalars | ecfp_absdiff_xgb | roc_auc  |  25 |             0.037352   |  0.0172954   |   0.0571258 |   0.64     |      0.28      |  0.0080688   | family_scaffold_source_purged      | nan           | nan             | 0.814949    |
| raw_chembl    | rich_diff_scalars | ecfp_absdiff_xgb | pr_auc   |  25 |             0.0397899  |  0.0223974   |   0.0588555 |   0.76     |      0.52      |  0.000556409 | family_scaffold_source_purged      | nan           | nan             | 0.0712204   |
| raw_chembl    | rich_diff_scalars | ecfp_absdiff_xgb | mcc      |  25 |             0.0680282  |  0.0308492   |   0.106472  |   0.68     |      0.36      |  0.00882232  | family_scaffold_source_purged      | nan           | nan             | 0.864588    |
| raw_chembl    | rich_scalars_only | ecfp_absdiff_xgb | roc_auc  |  25 |             0.0450222  |  0.0236026   |   0.0665943 |   0.64     |      0.28      |  0.00250793  | family_scaffold_source_purged      | nan           | nan             | 0.275872    |
| raw_chembl    | rich_scalars_only | ecfp_absdiff_xgb | pr_auc   |  25 |             0.0467234  |  0.02957     |   0.0642998 |   1        |      1         |  5.96046e-08 | family_scaffold_source_purged      | nan           | nan             | 1.19209e-05 |
| raw_chembl    | rich_scalars_only | ecfp_absdiff_xgb | mcc      |  25 |             0.0896826  |  0.0556954   |   0.123673  |   0.92     |      0.84      |  2.563e-06   | family_scaffold_source_purged      | nan           | nan             | 0.000428021 |
| raw_chembl    | rich_diff_scalars | ecfp_absdiff_xgb | roc_auc  |  25 |             0.110448   |  0.0911187   |   0.128638  |   1        |      1         |  5.96046e-08 | temporal_forward                   | nan           | nan             | 1.19209e-05 |
| raw_chembl    | rich_diff_scalars | ecfp_absdiff_xgb | pr_auc   |  25 |             0.0738535  |  0.0582657   |   0.0908792 |   1        |      1         |  5.96046e-08 | temporal_forward                   | nan           | nan             | 1.19209e-05 |
| raw_chembl    | rich_diff_scalars | ecfp_absdiff_xgb | mcc      |  25 |             0.130471   |  0.102989    |   0.157824  |   0.96     |      0.92      |  1.19209e-07 | temporal_forward                   | nan           | nan             | 2.03848e-05 |
| raw_chembl    | rich_scalars_only | ecfp_absdiff_xgb | roc_auc  |  25 |             0.118818   |  0.0972694   |   0.142055  |   1        |      1         |  5.96046e-08 | temporal_forward                   | nan           | nan             | 1.19209e-05 |
| raw_chembl    | rich_scalars_only | ecfp_absdiff_xgb | pr_auc   |  25 |             0.0851612  |  0.0666149   |   0.103376  |   0.96     |      0.92      |  1.19209e-07 | temporal_forward                   | nan           | nan             | 2.03848e-05 |
| raw_chembl    | rich_scalars_only | ecfp_absdiff_xgb | mcc      |  25 |             0.159414   |  0.132691    |   0.187257  |   1        |      1         |  5.96046e-08 | temporal_forward                   | nan           | nan             | 1.19209e-05 |
| raw_chembl    | rich_diff_scalars | ecfp_absdiff_xgb | roc_auc  |  25 |             0.110366   |  0.089585    |   0.129262  |   1        |      1         |  5.96046e-08 | temporal_scaffold_source_purged    | nan           | nan             | 1.19209e-05 |
| raw_chembl    | rich_diff_scalars | ecfp_absdiff_xgb | pr_auc   |  25 |             0.070128   |  0.0564949   |   0.0825042 |   1        |      1         |  5.96046e-08 | temporal_scaffold_source_purged    | nan           | nan             | 1.19209e-05 |
| raw_chembl    | rich_diff_scalars | ecfp_absdiff_xgb | mcc      |  25 |             0.104989   |  0.0883413   |   0.120725  |   1        |      1         |  5.96046e-08 | temporal_scaffold_source_purged    | nan           | nan             | 1.19209e-05 |
| raw_chembl    | rich_scalars_only | ecfp_absdiff_xgb | roc_auc  |  25 |             0.107592   |  0.0803288   |   0.134616  |   0.88     |      0.76      |  1.13249e-06 | temporal_scaffold_source_purged    | nan           | nan             | 0.000190258 |
| raw_chembl    | rich_scalars_only | ecfp_absdiff_xgb | pr_auc   |  25 |             0.0691278  |  0.0522256   |   0.0846329 |   0.92     |      0.84      |  4.17233e-07 | temporal_scaffold_source_purged    | nan           | nan             | 7.05123e-05 |
| raw_chembl    | rich_scalars_only | ecfp_absdiff_xgb | mcc      |  25 |             0.127986   |  0.10802     |   0.149883  |   1        |      1         |  5.96046e-08 | temporal_scaffold_source_purged    | nan           | nan             | 1.19209e-05 |
| raw_bindingdb | rich_diff_scalars | ecfp_absdiff_xgb | roc_auc  |  25 |             0.0663405  |  0.0596966   |   0.0732734 |   1        |      1         |  5.96046e-08 | target_family                      | nan           | nan             | 1.19209e-05 |
| raw_bindingdb | rich_diff_scalars | ecfp_absdiff_xgb | pr_auc   |  25 |             0.0601892  |  0.0523213   |   0.0690556 |   1        |      1         |  5.96046e-08 | target_family                      | nan           | nan             | 1.19209e-05 |
| raw_bindingdb | rich_diff_scalars | ecfp_absdiff_xgb | mcc      |  25 |             0.063402   |  0.0548359   |   0.0712004 |   1        |      1         |  5.96046e-08 | target_family                      | nan           | nan             | 1.19209e-05 |
| raw_bindingdb | rich_scalars_only | ecfp_absdiff_xgb | roc_auc  |  25 |             0.0679191  |  0.0609404   |   0.0748264 |   1        |      1         |  5.96046e-08 | target_family                      | nan           | nan             | 1.19209e-05 |
| raw_bindingdb | rich_scalars_only | ecfp_absdiff_xgb | pr_auc   |  25 |             0.0652203  |  0.0558129   |   0.0749695 |   1        |      1         |  5.96046e-08 | target_family                      | nan           | nan             | 1.19209e-05 |
| raw_bindingdb | rich_scalars_only | ecfp_absdiff_xgb | mcc      |  25 |             0.0698222  |  0.0605934   |   0.0785691 |   1        |      1         |  5.96046e-08 | target_family                      | nan           | nan             | 1.19209e-05 |
| raw_bindingdb | rich_diff_scalars | ecfp_absdiff_xgb | roc_auc  |  25 |             0.0553475  |  0.0430574   |   0.0684146 |   1        |      1         |  5.96046e-08 | family_scaffold_source_purged      | nan           | nan             | 1.19209e-05 |
| raw_bindingdb | rich_diff_scalars | ecfp_absdiff_xgb | pr_auc   |  25 |             0.0484917  |  0.0351933   |   0.0634942 |   1        |      1         |  5.96046e-08 | family_scaffold_source_purged      | nan           | nan             | 1.19209e-05 |
| raw_bindingdb | rich_diff_scalars | ecfp_absdiff_xgb | mcc      |  25 |             0.0580622  |  0.0469824   |   0.0690094 |   1        |      1         |  5.96046e-08 | family_scaffold_source_purged      | nan           | nan             | 1.19209e-05 |
| raw_bindingdb | rich_scalars_only | ecfp_absdiff_xgb | roc_auc  |  25 |             0.0535828  |  0.0429319   |   0.0642928 |   1        |      1         |  5.96046e-08 | family_scaffold_source_purged      | nan           | nan             | 1.19209e-05 |
| raw_bindingdb | rich_scalars_only | ecfp_absdiff_xgb | pr_auc   |  25 |             0.0484455  |  0.0360369   |   0.0618403 |   1        |      1         |  5.96046e-08 | family_scaffold_source_purged      | nan           | nan             | 1.19209e-05 |
| raw_bindingdb | rich_scalars_only | ecfp_absdiff_xgb | mcc      |  25 |             0.0551397  |  0.0466021   |   0.0629464 |   1        |      1         |  5.96046e-08 | family_scaffold_source_purged      | nan           | nan             | 1.19209e-05 |
| lohi          | rich_diff_scalars | ecfp_absdiff_xgb | roc_auc  |  15 |            -0.00432685 | -0.0330021   |   0.0235337 |   0.533333 |      0.0666667 |  1           | lohi_pair_scaffold_source_purged   | chembl_raw    | nan             | 1           |
| lohi          | rich_diff_scalars | ecfp_absdiff_xgb | pr_auc   |  15 |             0.0094926  | -0.0137943   |   0.0324559 |   0.533333 |      0.0666667 |  0.421204    | lohi_pair_scaffold_source_purged   | chembl_raw    | nan             | 1           |
| lohi          | rich_diff_scalars | ecfp_absdiff_xgb | mcc      |  15 |            -0.00629455 | -0.0599567   |   0.0448392 |   0.6      |      0.2       |  0.934082    | lohi_pair_scaffold_source_purged   | chembl_raw    | nan             | 1           |
| lohi          | rich_scalars_only | ecfp_absdiff_xgb | roc_auc  |  15 |            -0.00867408 | -0.0338971   |   0.0142102 |   0.533333 |      0.0666667 |  0.846924    | lohi_pair_scaffold_source_purged   | chembl_raw    | nan             | 1           |
| lohi          | rich_scalars_only | ecfp_absdiff_xgb | pr_auc   |  15 |             0.00486191 | -0.01484     |   0.0236489 |   0.6      |      0.2       |  0.599487    | lohi_pair_scaffold_source_purged   | chembl_raw    | nan             | 1           |
| lohi          | rich_scalars_only | ecfp_absdiff_xgb | mcc      |  15 |            -0.0139897  | -0.0602669   |   0.0284021 |   0.6      |      0.2       |  0.934082    | lohi_pair_scaffold_source_purged   | chembl_raw    | nan             | 1           |
| lohi          | rich_diff_scalars | ecfp_absdiff_xgb | roc_auc  |  15 |             0.0140204  | -0.0119599   |   0.0413237 |   0.466667 |     -0.0666667 |  0.524475    | lohi_ligand_scaffold_source_purged | chembl_raw    | nan             | 1           |
| lohi          | rich_diff_scalars | ecfp_absdiff_xgb | pr_auc   |  15 |             0.0269133  |  0.008192    |   0.046598  |   0.866667 |      0.733333  |  0.0124512   | lohi_ligand_scaffold_source_purged | chembl_raw    | nan             | 1           |
| lohi          | rich_diff_scalars | ecfp_absdiff_xgb | mcc      |  15 |             0.0351129  |  0.000626055 |   0.0737824 |   0.6      |      0.2       |  0.168823    | lohi_ligand_scaffold_source_purged | chembl_raw    | nan             | 1           |
| lohi          | rich_scalars_only | ecfp_absdiff_xgb | roc_auc  |  15 |             0.0194825  | -0.00697204  |   0.045649  |   0.533333 |      0.0666667 |  0.276855    | lohi_ligand_scaffold_source_purged | chembl_raw    | nan             | 1           |
| lohi          | rich_scalars_only | ecfp_absdiff_xgb | pr_auc   |  15 |             0.0322809  |  0.0119887   |   0.0528223 |   0.8      |      0.6       |  0.0180664   | lohi_ligand_scaffold_source_purged | chembl_raw    | nan             | 1           |
| lohi          | rich_scalars_only | ecfp_absdiff_xgb | mcc      |  15 |             0.035133   | -0.00231806  |   0.072333  |   0.666667 |      0.333333  |  0.151428    | lohi_ligand_scaffold_source_purged | chembl_raw    | nan             | 1           |
| lohi          | rich_diff_scalars | ecfp_absdiff_xgb | roc_auc  |  15 |             0.0362834  |  0.025161    |   0.0477077 |   0.866667 |      0.733333  |  0.000305176 | lohi_pair_scaffold_source_purged   | bindingdb_raw | nan             | 0.0415039   |
| lohi          | rich_diff_scalars | ecfp_absdiff_xgb | pr_auc   |  15 |             0.0331158  |  0.020406    |   0.0459226 |   0.8      |      0.6       |  0.000854492 | lohi_pair_scaffold_source_purged   | bindingdb_raw | nan             | 0.105103    |
| lohi          | rich_diff_scalars | ecfp_absdiff_xgb | mcc      |  15 |             0.0429784  |  0.0248042   |   0.0622859 |   0.8      |      0.6       |  0.00152588  | lohi_pair_scaffold_source_purged   | bindingdb_raw | nan             | 0.178528    |
| lohi          | rich_scalars_only | ecfp_absdiff_xgb | roc_auc  |  15 |             0.033804   |  0.0215274   |   0.047006  |   0.933333 |      0.866667  |  0.000183105 | lohi_pair_scaffold_source_purged   | bindingdb_raw | nan             | 0.0258179   |
| lohi          | rich_scalars_only | ecfp_absdiff_xgb | pr_auc   |  15 |             0.030166   |  0.0172541   |   0.0425703 |   0.8      |      0.6       |  0.00152588  | lohi_pair_scaffold_source_purged   | bindingdb_raw | nan             | 0.178528    |
| lohi          | rich_scalars_only | ecfp_absdiff_xgb | mcc      |  15 |             0.0461268  |  0.0279948   |   0.0657435 |   0.866667 |      0.733333  |  0.000305176 | lohi_pair_scaffold_source_purged   | bindingdb_raw | nan             | 0.0415039   |
| lohi          | rich_diff_scalars | ecfp_absdiff_xgb | roc_auc  |  15 |             0.0433609  |  0.0323762   |   0.0535814 |   1        |      1         |  6.10352e-05 | lohi_ligand_scaffold_source_purged | bindingdb_raw | nan             | 0.0100708   |
| lohi          | rich_diff_scalars | ecfp_absdiff_xgb | pr_auc   |  15 |             0.0385044  |  0.0305556   |   0.0460963 |   1        |      1         |  6.10352e-05 | lohi_ligand_scaffold_source_purged | bindingdb_raw | nan             | 0.0100708   |
| lohi          | rich_diff_scalars | ecfp_absdiff_xgb | mcc      |  15 |             0.0601684  |  0.0414184   |   0.078778  |   0.933333 |      0.866667  |  0.000183105 | lohi_ligand_scaffold_source_purged | bindingdb_raw | nan             | 0.0258179   |
| lohi          | rich_scalars_only | ecfp_absdiff_xgb | roc_auc  |  15 |             0.0391099  |  0.028984    |   0.0491903 |   0.933333 |      0.866667  |  0.00012207  | lohi_ligand_scaffold_source_purged | bindingdb_raw | nan             | 0.0178223   |
| lohi          | rich_scalars_only | ecfp_absdiff_xgb | pr_auc   |  15 |             0.0344499  |  0.0265546   |   0.0425659 |   1        |      1         |  6.10352e-05 | lohi_ligand_scaffold_source_purged | bindingdb_raw | nan             | 0.0100708   |
| lohi          | rich_scalars_only | ecfp_absdiff_xgb | mcc      |  15 |             0.0550361  |  0.0375357   |   0.0710918 |   0.933333 |      0.866667  |  0.000183105 | lohi_ligand_scaffold_source_purged | bindingdb_raw | nan             | 0.0258179   |
| noise         | rich_diff_scalars | ecfp_absdiff_xgb | roc_auc  |  15 |             0.0382522  |  0.0103115   |   0.066337  |   0.6      |      0.2       |  0.0946045   | family_scaffold_source_purged      | chembl_raw    | all_pairs       | 1           |
| noise         | rich_diff_scalars | ecfp_absdiff_xgb | pr_auc   |  15 |             0.0425775  |  0.0178002   |   0.0707921 |   0.733333 |      0.466667  |  0.0150757   | family_scaffold_source_purged      | chembl_raw    | all_pairs       | 1           |
| noise         | rich_diff_scalars | ecfp_absdiff_xgb | mcc      |  15 |             0.0714249  |  0.0224521   |   0.124475  |   0.733333 |      0.466667  |  0.0353394   | family_scaffold_source_purged      | chembl_raw    | all_pairs       | 1           |
| noise         | rich_scalars_only | ecfp_absdiff_xgb | roc_auc  |  15 |             0.0453019  |  0.0130834   |   0.0770692 |   0.666667 |      0.333333  |  0.0353394   | family_scaffold_source_purged      | chembl_raw    | all_pairs       | 1           |
| noise         | rich_scalars_only | ecfp_absdiff_xgb | pr_auc   |  15 |             0.0489552  |  0.0236381   |   0.0755248 |   0.866667 |      0.733333  |  0.000854492 | family_scaffold_source_purged      | chembl_raw    | all_pairs       | 0.105103    |
| noise         | rich_scalars_only | ecfp_absdiff_xgb | mcc      |  15 |             0.0835347  |  0.0328326   |   0.135114  |   0.733333 |      0.466667  |  0.0124512   | family_scaffold_source_purged      | chembl_raw    | all_pairs       | 1           |
| noise         | rich_diff_scalars | ecfp_absdiff_xgb | roc_auc  |  15 |             0.0351667  | -0.0050049   |   0.076604  |   0.466667 |     -0.0666667 |  0.359131    | family_scaffold_source_purged      | chembl_raw    | high_margin     | 1           |
| noise         | rich_diff_scalars | ecfp_absdiff_xgb | pr_auc   |  15 |             0.0565433  |  0.00995803  |   0.102359  |   0.6      |      0.2       |  0.151428    | family_scaffold_source_purged      | chembl_raw    | high_margin     | 1           |
| noise         | rich_diff_scalars | ecfp_absdiff_xgb | mcc      |  15 |             0.0523264  | -0.00623932  |   0.112702  |   0.533333 |      0.0666667 |  0.302795    | family_scaffold_source_purged      | chembl_raw    | high_margin     | 1           |
| noise         | rich_scalars_only | ecfp_absdiff_xgb | roc_auc  |  15 |             0.0475858  |  0.00185687  |   0.0949584 |   0.533333 |      0.0666667 |  0.276855    | family_scaffold_source_purged      | chembl_raw    | high_margin     | 1           |
| noise         | rich_scalars_only | ecfp_absdiff_xgb | pr_auc   |  15 |             0.0627634  |  0.0170778   |   0.109437  |   0.666667 |      0.333333  |  0.072998    | family_scaffold_source_purged      | chembl_raw    | high_margin     | 1           |
| noise         | rich_scalars_only | ecfp_absdiff_xgb | mcc      |  15 |             0.0853378  |  0.0212013   |   0.154879  |   0.666667 |      0.333333  |  0.168823    | family_scaffold_source_purged      | chembl_raw    | high_margin     | 1           |
| noise         | rich_diff_scalars | ecfp_absdiff_xgb | roc_auc  |  15 |             0.0706492  |  0.0348858   |   0.105456  |   0.8      |      0.6       |  0.00262451  | family_scaffold_source_purged      | chembl_raw    | multi_source    | 0.286072    |
| noise         | rich_diff_scalars | ecfp_absdiff_xgb | pr_auc   |  15 |             0.0597368  |  0.0385954   |   0.0798175 |   0.866667 |      0.733333  |  0.000305176 | family_scaffold_source_purged      | chembl_raw    | multi_source    | 0.0415039   |
| noise         | rich_diff_scalars | ecfp_absdiff_xgb | mcc      |  15 |             0.138689   |  0.0876032   |   0.187155  |   0.866667 |      0.733333  |  0.000610352 | family_scaffold_source_purged      | chembl_raw    | multi_source    | 0.0775146   |
| noise         | rich_scalars_only | ecfp_absdiff_xgb | roc_auc  |  15 |             0.0636896  |  0.0262997   |   0.100563  |   0.8      |      0.6       |  0.00671387  | family_scaffold_source_purged      | chembl_raw    | multi_source    | 0.691528    |
| noise         | rich_scalars_only | ecfp_absdiff_xgb | pr_auc   |  15 |             0.0563254  |  0.0366836   |   0.0768614 |   0.933333 |      0.866667  |  0.00012207  | family_scaffold_source_purged      | chembl_raw    | multi_source    | 0.0178223   |
| noise         | rich_scalars_only | ecfp_absdiff_xgb | mcc      |  15 |             0.117804   |  0.066176    |   0.170472  |   0.933333 |      0.866667  |  0.000305176 | family_scaffold_source_purged      | chembl_raw    | multi_source    | 0.0415039   |
| noise         | rich_diff_scalars | ecfp_absdiff_xgb | roc_auc  |  15 |             0.0252532  | -0.00574863  |   0.0562723 |   0.6      |      0.2       |  0.187622    | family_scaffold_source_purged      | chembl_raw    | high_similarity | 1           |
| noise         | rich_diff_scalars | ecfp_absdiff_xgb | pr_auc   |  15 |             0.0537872  |  0.0308097   |   0.0737889 |   0.866667 |      0.733333  |  0.000610352 | family_scaffold_source_purged      | chembl_raw    | high_similarity | 0.0775146   |

## P1 Case Studies

- selected cases: 80
- file: `results/tkde_failure_case_studies.md`

## P2 Leaderboard Package

- leaderboard rows: 131
- schema: `docs/tkde_leaderboard_schema.md`

## Manuscript Reading

- The added P0 block directly addresses the strongest Transactions objections: split definition, pretrained baselines, raw reliability, and statistical multiplicity.
- The added P1 block turns raw-database noise from a vulnerability into an explicit robustness and limitation analysis.
- The new remaining-gap block adds cross-database transfer, prospective temporal validation, top-k enrichment, risk-coverage, calibration-bin data, and target-level failure taxonomy.
- The added P2 block makes the work look like a reusable benchmark/protocol package, which is the cleanest TKDE positioning without a new neural architecture.
