# Neural Cliff Baselines

## Mean/Std Metrics

| benchmark        | split_mode              | variant        |   roc_auc_mean |   roc_auc_std |   pr_auc_mean |   pr_auc_std |   balanced_accuracy_mean |   f1_mean |   mcc_mean |   abs_delta_spearman_mean |   pair_overlap_rate_mean |   any_ligand_overlap_rate_mean |
|:-----------------|:------------------------|:---------------|---------------:|--------------:|--------------:|-------------:|-------------------------:|----------:|-----------:|--------------------------:|-------------------------:|-------------------------------:|
| PC-CFRL-internal | target_cluster_balanced | graph_gcn_pair |        0.48918 |       0.06783 |       0.17254 |      0.08871 |                  0.50452 |   0.2565  |    0.0073  |                   0.00797 |                  0.11507 |                        0.12117 |
| PC-CFRL-internal | ligand_component        | graph_gcn_pair |        0.52505 |       0.08648 |       0.20898 |      0.10576 |                  0.50892 |   0.26598 |    0.02272 |                   0.05903 |                  0       |                        0       |

## Paired Deltas

No paired deltas available.

## Interpretation Template

Use target-conditioned neural variants as positive evidence only when they beat both ligand-only/siamese neural baselines and shuffled-target controls on the same split. Otherwise they are SOTA-style controls or negative evidence.
