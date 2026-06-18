# Model Card

## Model Family

PC-CFRL is reported as a simple counterfactual pair representation coupled to XGBoost or non-neural tree controls. It is not presented as a new large neural architecture.

## Inputs

Each instance is a ligand pair within a target context. Features include absolute-difference ECFP fingerprints, scalar descriptor differences, pair similarity/activity-gap fields, and target/source fields available under the audited split contract.

## Outputs

The model outputs an activity-cliff probability. Evaluation reports ROC-AUC, PR-AUC, MCC, BEDROC20, enrichment factors, selective-prediction precision, conformal singleton precision, calibration diagnostics, paired effects, and runtime.

## Intended Use

The intended use is benchmark evaluation for leakage-controlled pairwise activity-cliff prediction. The model should be treated as a strong reproducible baseline and audit instrument within PC-CFRL-Bench, not as a clinically validated screening tool.

## Claim Boundary

- No biochemical wet-lab validation is included.
- RDKit case figures and SAR tables are qualitative diagnostics, not mechanistic proof.
- Official DataSAIL C2/ECFP at 80 clusters failed to return assignments in the current solver/time setting.
- ChEMBL official DataSAIL rows remain a boundary case where PC-CFRL is not uniformly better than ECFP.
