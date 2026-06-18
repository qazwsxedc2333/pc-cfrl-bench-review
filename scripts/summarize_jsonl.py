from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("jsonl")
    parser.add_argument("--output", default="")
    args = parser.parse_args()

    rows = []
    with Path(args.jsonl).open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                rows.append(json.loads(line))
    df = pd.DataFrame(rows)
    ok = df[~df.get("rmse", pd.Series(index=df.index)).isna()].copy()
    if ok.empty:
        print("No successful rows.")
        print(df)
        return

    metrics = [
        "rmse",
        "mae",
        "pearson",
        "spearman",
        "ci",
        "pairwise_rank_acc",
        "delta_spearman",
        "cliff_mae",
        "cliff_auc",
        "cliff_pr_auc",
    ]
    summary = ok.groupby("model")[metrics].mean(numeric_only=True).sort_values("rmse")
    counts = ok.groupby("model")["dataset"].nunique().rename("datasets")
    summary = counts.to_frame().join(summary)
    print(summary.round(4).to_string())

    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out)
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()

