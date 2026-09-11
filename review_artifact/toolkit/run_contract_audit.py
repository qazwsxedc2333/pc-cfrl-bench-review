"""Generate PC-CFRL-Bench contract-level leaderboard deltas.

This script reads the frozen leaderboard snapshot distributed with the
anonymous artifact and compares PC-CFRL-full (`rich_diff_scalars`) with
ECFP-XGB (`ecfp_absdiff_xgb`) inside the same benchmark/condition/split
contract. It does not train models or change any manuscript result.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


KEYS = ["benchmark", "source", "condition", "split_mode"]
BASELINE = "ecfp_absdiff_xgb"
CANDIDATE = "rich_diff_scalars"


def _float(row: dict[str, str], key: str) -> float | None:
    value = row.get(key, "")
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def build_deltas(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str], dict[str, dict[str, str]]] = {}
    for row in rows:
        key = tuple(row.get(k, "") for k in KEYS)
        grouped.setdefault(key, {})[row.get("variant", "")] = row

    deltas: list[dict[str, object]] = []
    for key, variants in sorted(grouped.items()):
        if BASELINE not in variants or CANDIDATE not in variants:
            continue
        base = variants[BASELINE]
        cand = variants[CANDIDATE]
        out: dict[str, object] = dict(zip(KEYS, key))
        out.update(
            {
                "baseline_variant": BASELINE,
                "candidate_variant": CANDIDATE,
                "baseline_roc_auc": _float(base, "roc_auc_mean"),
                "candidate_roc_auc": _float(cand, "roc_auc_mean"),
                "delta_roc_auc": None,
                "baseline_pr_auc": _float(base, "pr_auc_mean"),
                "candidate_pr_auc": _float(cand, "pr_auc_mean"),
                "delta_pr_auc": None,
                "baseline_mcc": _float(base, "mcc_mean"),
                "candidate_mcc": _float(cand, "mcc_mean"),
                "delta_mcc": None,
                "n_train_mean": _float(cand, "n_train_mean"),
                "n_test_mean": _float(cand, "n_test_mean"),
                "target_overlap_rate_mean": _float(cand, "target_overlap_rate_mean"),
                "exact_ligand_overlap_rate_mean": _float(cand, "exact_ligand_overlap_rate_mean"),
                "scaffold_overlap_rate_mean": _float(cand, "scaffold_overlap_rate_mean"),
                "document_source_overlap_rate_mean": _float(cand, "document_source_overlap_rate_mean"),
            }
        )
        for metric in ["roc_auc", "pr_auc", "mcc"]:
            b = out[f"baseline_{metric}"]
            c = out[f"candidate_{metric}"]
            if b is not None and c is not None:
                out[f"delta_{metric}"] = round(float(c) - float(b), 10)
        deltas.append(out)
    return deltas


def write_csv(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        with path.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize(deltas: list[dict[str, object]]) -> dict[str, object]:
    favorable = [d for d in deltas if (d.get("delta_roc_auc") is not None and float(d["delta_roc_auc"]) > 0)]
    roc_deltas = [float(d["delta_roc_auc"]) for d in deltas if d.get("delta_roc_auc") is not None]
    return {
        "baseline_variant": BASELINE,
        "candidate_variant": CANDIDATE,
        "n_contracts_with_both_variants": len(deltas),
        "n_favorable_roc_auc": len(favorable),
        "favorable_roc_auc_rate": round(len(favorable) / len(deltas), 6) if deltas else None,
        "mean_delta_roc_auc": round(sum(roc_deltas) / len(roc_deltas), 6) if roc_deltas else None,
        "min_delta_roc_auc": round(min(roc_deltas), 6) if roc_deltas else None,
        "max_delta_roc_auc": round(max(roc_deltas), 6) if roc_deltas else None,
        "note": "Computed from frozen leaderboard rows; not a new training run.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--leaderboard", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    rows = load_rows(args.leaderboard)
    deltas = build_deltas(rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(deltas, args.out_dir / "contract_leaderboard_deltas.csv")
    summary_path = args.out_dir / "contract_audit_summary.json"
    with summary_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(summarize(deltas), indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
