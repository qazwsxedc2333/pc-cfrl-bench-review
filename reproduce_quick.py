from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
LEADERBOARD = ROOT / "review_artifact" / "case_studies" / "tkde_leaderboard_snapshot.csv"
OUT_DIR = ROOT / "review_artifact" / "toolkit" / "generated"
AUDIT_SCRIPT = ROOT / "review_artifact" / "toolkit" / "run_contract_audit.py"


def main() -> None:
    if not LEADERBOARD.exists():
        raise SystemExit(f"Missing leaderboard snapshot: {LEADERBOARD}")
    subprocess.run(
        [
            sys.executable,
            str(AUDIT_SCRIPT),
            "--leaderboard",
            str(LEADERBOARD),
            "--out-dir",
            str(OUT_DIR),
        ],
        check=True,
    )
    summary_path = OUT_DIR / "contract_audit_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    print(f"contracts_with_both_variants: {summary['n_contracts_with_both_variants']}")
    print(f"favorable_roc_auc: {summary['n_favorable_roc_auc']}")
    print(f"mean_delta_roc_auc: {summary['mean_delta_roc_auc']}")
    print(f"summary: {summary_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

