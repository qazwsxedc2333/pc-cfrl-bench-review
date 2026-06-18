from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    matthews_corrcoef,
    roc_auc_score,
)


METRIC_ORDER = [
    "roc_auc",
    "pr_auc",
    "balanced_accuracy",
    "f1",
    "mcc",
    "brier",
    "ece_10",
]


def sha256_file(path: str | Path, chunk_size: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def expected_calibration_error(y: np.ndarray, prob: np.ndarray, n_bins: int = 10) -> float:
    y = np.asarray(y, dtype=np.float32)
    prob = np.asarray(prob, dtype=np.float32)
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    for left, right in zip(edges[:-1], edges[1:]):
        mask = (prob >= left) & (prob <= right) if right == 1.0 else (prob >= left) & (prob < right)
        if np.any(mask):
            ece += float(mask.mean()) * abs(float(prob[mask].mean()) - float(y[mask].mean()))
    return float(ece)


def evaluate_binary(y_true: np.ndarray, prob_pos: np.ndarray, threshold: float = 0.5) -> dict[str, float]:
    y = np.asarray(y_true, dtype=np.int32)
    prob = np.asarray(prob_pos, dtype=np.float32)
    pred = (prob >= threshold).astype(np.int32)
    return {
        "roc_auc": float(roc_auc_score(y, prob)) if len(np.unique(y)) == 2 else float("nan"),
        "pr_auc": float(average_precision_score(y, prob)) if len(np.unique(y)) == 2 else float("nan"),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)) if len(np.unique(y)) == 2 else float("nan"),
        "f1": float(f1_score(y, pred, zero_division=0)),
        "mcc": float(matthews_corrcoef(y, pred)) if len(np.unique(y)) > 1 and len(np.unique(pred)) > 1 else 0.0,
        "brier": float(brier_score_loss(y, prob)),
        "ece_10": expected_calibration_error(y, prob, 10),
    }


def summarize_results(df: pd.DataFrame, group_cols: list[str] | None = None) -> pd.DataFrame:
    if group_cols is None:
        group_cols = ["benchmark", "split_mode", "variant"]
    metric_cols = [c for c in METRIC_ORDER if c in df.columns]
    audit_cols = [
        c
        for c in [
            "n_train",
            "n_test",
            "train_retention_after_purge",
            "target_overlap_rate",
            "family_overlap_rate",
            "exact_ligand_overlap_rate",
            "scaffold_overlap_rate",
            "document_source_overlap_rate",
            "assay_source_overlap_rate",
            "train_seconds",
            "predict_seconds",
        ]
        if c in df.columns
    ]
    out = df.groupby(group_cols, as_index=False)[metric_cols + audit_cols].agg(["mean", "std", "count"])
    out.columns = ["_".join([x for x in col if x]).strip("_") for col in out.columns.to_flat_index()]
    return out.reset_index(drop=True)


def markdown_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if df.empty:
        return "_No rows._"
    out = df.copy()
    if max_rows is not None:
        out = out.head(max_rows)
    for col in out.columns:
        if pd.api.types.is_float_dtype(out[col]):
            out[col] = out[col].map(lambda x: "" if pd.isna(x) else f"{float(x):.6g}")
    return out.to_markdown(index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", required=True, help="CSV/JSONL with y_true and prob_pos columns.")
    parser.add_argument("--output", default="")
    args = parser.parse_args()

    path = Path(args.predictions)
    if path.suffix == ".jsonl":
        df = pd.read_json(path, lines=True)
    else:
        df = pd.read_csv(path)
    required = {"y_true", "prob_pos"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    metrics = evaluate_binary(df["y_true"].to_numpy(), df["prob_pos"].to_numpy())
    text = json.dumps(metrics, indent=2, sort_keys=True) + "\n"
    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
