# Official DataSAIL Exact Split Attempt

Date: 2026-06-11

- Main experiment environment: `.env` uses Python 3.13.13.
- Official DataSAIL package range: Python >=3.9,<3.13.
- Follow-up action: created `.env_datasail_py312` with Python 3.12.13 in a separate compute environment to avoid contaminating the frozen main environment.
- Outcome: official DataSAIL 1.3.0 was installed in the separate Python 3.12 environment and used for the C2/ECFP exact-split audit.
- Manuscript implication: official DataSAIL exact-split evidence is now reported as an additional audit in Supplementary Tables S21--S22, while the main text still uses the broader wording `LoHi/DataSAIL-style` for the original stress-test block.
