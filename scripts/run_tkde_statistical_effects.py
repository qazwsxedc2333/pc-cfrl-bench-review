from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


SPECS = [
    ("raw_chembl", "results/chembl_raw_source_temporal_pairs_seed0_4_n40.jsonl", ["split_mode"]),
    ("raw_bindingdb", "results/bindingdb_raw_source_pairs_seed0_4_n80.jsonl", ["split_mode"]),
    ("lohi", "results/tkde_lohi_leakage_splits_seed0_2.jsonl", ["source", "split_mode"]),
    ("noise", "results/tkde_noise_sensitivity_seed0_2.jsonl", ["source", "condition", "split_mode"]),
]


def load_jsonl(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    return pd.DataFrame(json.loads(line) for line in p.read_text(encoding="utf-8").splitlines() if line.strip())


def bootstrap_ci(values: np.ndarray, n_boot: int, seed: int) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    means = [float(rng.choice(values, size=len(values), replace=True).mean()) for _ in range(n_boot)]
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def cliffs_delta(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return float("nan")
    return float(((values > 0).sum() - (values < 0).sum()) / len(values))


def holm(pvals: list[float]) -> list[float]:
    order = np.argsort([1.0 if not np.isfinite(p) else p for p in pvals])
    adjusted = np.ones(len(pvals), dtype=float)
    running = 0.0
    m = len(pvals)
    for rank, idx in enumerate(order):
        p = pvals[idx]
        val = 1.0 if not np.isfinite(p) else min(1.0, (m - rank) * p)
        running = max(running, val)
        adjusted[idx] = running
    return adjusted.tolist()


def analyze_frame(name: str, df: pd.DataFrame, group_cols: list[str], n_boot: int, seed: int) -> list[dict]:
    if df.empty or "ecfp_absdiff_xgb" not in set(df.get("variant", [])):
        return []
    rows = []
    for keys, sub in df.groupby(group_cols, sort=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        base = sub[sub["variant"] == "ecfp_absdiff_xgb"].set_index(["seed", "fold"])
        for variant in sorted(sub["variant"].astype(str).unique()):
            if variant == "ecfp_absdiff_xgb":
                continue
            cur = sub[sub["variant"] == variant].set_index(["seed", "fold"])
            idx = cur.index.intersection(base.index)
            if len(idx) == 0:
                continue
            for metric in ["roc_auc", "pr_auc", "mcc", "brier", "ece_10"]:
                if metric not in cur or metric not in base:
                    continue
                diff = cur.loc[idx, metric].astype(float) - base.loc[idx, metric].astype(float)
                if metric in {"brier", "ece_10"}:
                    diff = -diff
                diff = diff.replace([np.inf, -np.inf], np.nan).dropna().to_numpy(dtype=float)
                if len(diff) == 0:
                    continue
                try:
                    p = float(wilcoxon(diff).pvalue) if np.any(diff != 0) else 1.0
                except ValueError:
                    p = float("nan")
                lo, hi = bootstrap_ci(diff, n_boot, seed + len(rows))
                row = {"experiment": name, "variant": variant, "baseline": "ecfp_absdiff_xgb", "metric": metric, "n": int(len(diff)), "mean_delta_favorable": float(diff.mean()), "ci95_low": lo, "ci95_high": hi, "win_rate": float((diff > 0).mean()), "cliffs_delta": cliffs_delta(diff), "wilcoxon_p": p}
                row.update(dict(zip(group_cols, keys)))
                rows.append(row)
    return rows


def write_report(df: pd.DataFrame, out: Path) -> None:
    focus = df[df["metric"].isin(["roc_auc", "pr_auc", "mcc"])].copy()
    strong = focus[(focus["ci95_low"] > 0) & (focus["holm_p"] < 0.05)]
    lines = [
        "# TKDE Statistical Effect-Size Audit",
        "",
        "Intervals are bootstrap 95% CIs over seed/fold paired deltas. Brier/ECE are sign-flipped so positive means favorable. P-values are Wilcoxon signed-rank with Holm correction.",
        "",
        "## Main Table",
        "",
        focus.round(6).to_markdown(index=False),
        "",
        "## Strong Positive Cells",
        "",
        strong.round(6).to_markdown(index=False) if not strong.empty else "_No Holm-significant positive cells._",
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-boot", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260518)
    parser.add_argument("--csv-out", default="results/tkde_statistical_effects_audit.csv")
    parser.add_argument("--report", default="results/tkde_statistical_effects_audit.md")
    args = parser.parse_args()
    rows = []
    for i, (name, path, groups) in enumerate(SPECS):
        rows.extend(analyze_frame(name, load_jsonl(path), groups, args.n_boot, args.seed + i * 1000))
    out = pd.DataFrame(rows)
    if out.empty:
        raise RuntimeError("No statistical rows produced")
    out["holm_p"] = holm(out["wilcoxon_p"].tolist())
    Path(args.csv_out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.csv_out, index=False)
    write_report(out, Path(args.report))
    print(json.dumps({"rows": int(len(out)), "csv": args.csv_out, "report": args.report}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
