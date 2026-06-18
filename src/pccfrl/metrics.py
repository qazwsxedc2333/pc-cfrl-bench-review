from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import average_precision_score, mean_absolute_error, mean_squared_error, roc_auc_score

from pccfrl.fingerprints import bitvectors_for_similarity, tanimoto


def _safe_corr(fn, y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if len(y_true) < 2 or np.std(y_true) == 0 or np.std(y_pred) == 0:
        return float("nan")
    return float(fn(y_true, y_pred).statistic)


def concordance_index(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    total = 0
    concordant = 0.0
    n = len(y_true)
    for i in range(n):
        for j in range(i + 1, n):
            dy = y_true[i] - y_true[j]
            if dy == 0:
                continue
            total += 1
            dp = y_pred[i] - y_pred[j]
            if dp == 0:
                concordant += 0.5
            elif np.sign(dy) == np.sign(dp):
                concordant += 1.0
    return float(concordant / total) if total else float("nan")


def regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    return {
        "rmse": float(mean_squared_error(y_true, y_pred) ** 0.5),
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "pearson": _safe_corr(pearsonr, y_true, y_pred),
        "spearman": _safe_corr(spearmanr, y_true, y_pred),
        "ci": concordance_index(y_true, y_pred),
    }


@dataclass
class CliffMetricConfig:
    tanimoto_threshold: float = 0.7
    cliff_delta: float = 1.0
    smooth_delta: float = 0.3


def cliff_pair_metrics(
    smiles: list[str],
    y_true: np.ndarray,
    y_pred: np.ndarray,
    config: CliffMetricConfig = CliffMetricConfig(),
) -> dict[str, float | int]:
    fps = bitvectors_for_similarity(smiles)
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)

    true_deltas: list[float] = []
    pred_deltas: list[float] = []
    cliff_rank_hits = 0
    cliff_rank_total = 0
    binary_labels: list[int] = []
    binary_scores: list[float] = []

    n = len(smiles)
    high_sim_pairs = 0
    for i in range(n):
        for j in range(i + 1, n):
            sim = tanimoto(fps[i], fps[j])
            if sim < config.tanimoto_threshold:
                continue
            high_sim_pairs += 1
            dy = y_true[i] - y_true[j]
            dp = y_pred[i] - y_pred[j]
            abs_dy = abs(dy)
            abs_dp = abs(dp)

            if abs_dy >= config.cliff_delta:
                true_deltas.append(abs_dy)
                pred_deltas.append(abs_dp)
                cliff_rank_total += 1
                if np.sign(dy) == np.sign(dp):
                    cliff_rank_hits += 1
                binary_labels.append(1)
                binary_scores.append(abs_dp)
            elif abs_dy <= config.smooth_delta:
                binary_labels.append(0)
                binary_scores.append(abs_dp)

    out: dict[str, float | int] = {
        "high_sim_pairs": high_sim_pairs,
        "cliff_pairs": cliff_rank_total,
        "pairwise_rank_acc": float(cliff_rank_hits / cliff_rank_total) if cliff_rank_total else float("nan"),
        "cliff_mae": float(np.mean(np.abs(np.asarray(true_deltas) - np.asarray(pred_deltas)))) if true_deltas else float("nan"),
        "delta_spearman": _safe_corr(spearmanr, np.asarray(true_deltas), np.asarray(pred_deltas)) if len(true_deltas) >= 2 else float("nan"),
        "cliff_auc": float("nan"),
        "cliff_pr_auc": float("nan"),
    }
    if len(set(binary_labels)) == 2:
        out["cliff_auc"] = float(roc_auc_score(binary_labels, binary_scores))
        out["cliff_pr_auc"] = float(average_precision_score(binary_labels, binary_scores))
    out.update({f"cliff_cfg_{k}": v for k, v in asdict(config).items()})
    return out

