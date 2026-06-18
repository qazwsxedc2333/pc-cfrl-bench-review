from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

CLASS_METRICS = [
    "roc_auc",
    "pr_auc",
    "mcc",
    "balanced_accuracy",
    "f1",
    "brier",
    "ece_10",
]

RANKING_METRICS = [
    "bedroc20",
    "ef_1pct",
    "ef_5pct",
    "precision_at_50",
    "precision_at_100",
]

REGRESSION_METRICS = [
    "rmse",
    "mae",
    "pearson",
    "spearman",
    "ci",
    "pairwise_rank_acc",
    "delta_spearman",
    "cliff_auc",
    "cliff_pr_auc",
    "cliff_mol_rmse",
]


@dataclass
class TableRecord:
    table_id: str
    title: str
    rows: int
    csv_path: str
    tex_path: str
    sources: str


def repo_path(path: str | Path) -> Path:
    return ROOT / Path(path)


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT)).replace("\\", "/")


def read_csv(path: str) -> pd.DataFrame:
    p = repo_path(path)
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p)


def read_jsonl(path: str) -> pd.DataFrame:
    p = repo_path(path)
    if not p.exists():
        return pd.DataFrame()
    return pd.read_json(p, lines=True)


def existing_sources(paths: Iterable[str]) -> str:
    return "; ".join(p for p in paths if repo_path(p).exists())


def fmt_float(value: object, digits: int = 4) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        return str(value)
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def fmt_int(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        return str(value)
    try:
        return f"{int(round(float(value)))}"
    except (TypeError, ValueError):
        return str(value)


def fmt_mean_std(row: pd.Series, metric: str, digits: int = 4) -> str:
    mean_col = f"{metric}_mean"
    std_col = f"{metric}_std"
    mean = row.get(mean_col, row.get(metric, np.nan))
    if pd.isna(mean):
        return ""
    std = row.get(std_col, np.nan)
    if pd.isna(std):
        return fmt_float(mean, digits)
    return f"{float(mean):.{digits}f} +/- {float(std):.{digits}f}"


def first_count(row: pd.Series, metrics: Iterable[str]) -> str:
    for metric in metrics:
        count = row.get(f"{metric}_count", np.nan)
        if not pd.isna(count):
            return fmt_int(count)
    return ""


def add_delta_columns(
    df: pd.DataFrame,
    group_cols: list[str],
    metrics: list[str],
    baseline_variant: str = "ecfp_absdiff_xgb",
) -> pd.DataFrame:
    if df.empty or "variant" not in df.columns:
        return df.copy()
    out = df.copy()
    for metric in metrics:
        mean_col = f"{metric}_mean"
        if mean_col not in out.columns:
            continue
        base = out[out["variant"].eq(baseline_variant)][group_cols + [mean_col]].copy()
        if base.empty:
            continue
        base = base.rename(columns={mean_col: f"baseline_{metric}_mean"})
        out = out.merge(base, on=group_cols, how="left")
        out[f"delta_{metric}_vs_ecfp"] = out[mean_col] - out[f"baseline_{metric}_mean"]
        out = out.drop(columns=[f"baseline_{metric}_mean"])
    return out


def class_table(
    df: pd.DataFrame,
    id_cols: list[str],
    group_cols: list[str] | None = None,
    metrics: list[str] | None = None,
    extra_mean_cols: list[str] | None = None,
    add_deltas: bool = True,
) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    metrics = metrics or CLASS_METRICS
    extra_mean_cols = extra_mean_cols or []
    raw = df.copy()
    if add_deltas and group_cols:
        raw = add_delta_columns(raw, group_cols, ["roc_auc", "pr_auc", "mcc"])
    rows: list[dict[str, object]] = []
    for _, row in raw.iterrows():
        item: dict[str, object] = {}
        for col in id_cols:
            if col in raw.columns:
                item[col] = row.get(col, "")
        item["n_runs"] = first_count(row, metrics)
        for col in extra_mean_cols:
            if f"{col}_mean" in raw.columns or col in raw.columns:
                digits = 0 if col.startswith("n_") else 4
                item[col] = fmt_mean_std(row, col, digits=digits)
        for metric in metrics:
            if f"{metric}_mean" in raw.columns or metric in raw.columns:
                item[metric] = fmt_mean_std(row, metric)
        for metric in ["roc_auc", "pr_auc", "mcc"]:
            delta = row.get(f"delta_{metric}_vs_ecfp", np.nan)
            if not pd.isna(delta):
                item[f"delta_{metric}_vs_ecfp"] = fmt_float(delta)
        rows.append(item)
    return pd.DataFrame(rows)


def regression_summary(df: pd.DataFrame, group_cols: list[str], metrics: list[str]) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    present = [m for m in metrics if m in df.columns]
    if not present:
        return pd.DataFrame()
    grouped = df.groupby(group_cols, dropna=False)[present].agg(["mean", "std", "count"]).reset_index()
    grouped.columns = [
        "_".join(str(x) for x in col if x) if isinstance(col, tuple) else str(col)
        for col in grouped.columns
    ]
    return grouped


def regression_table(df: pd.DataFrame, id_cols: list[str], metrics: list[str] | None = None) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    metrics = metrics or REGRESSION_METRICS
    rows: list[dict[str, object]] = []
    for _, row in df.iterrows():
        item: dict[str, object] = {}
        for col in id_cols:
            if col in df.columns:
                item[col] = row.get(col, "")
        item["n_runs"] = first_count(row, metrics)
        for metric in metrics:
            if f"{metric}_mean" in df.columns:
                item[metric] = fmt_mean_std(row, metric)
        rows.append(item)
    return pd.DataFrame(rows)


def compact_numeric(df: pd.DataFrame, digits: int = 4) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_float_dtype(out[col]):
            out[col] = out[col].round(digits)
    return out


def save_table(
    records: list[TableRecord],
    out_dir: Path,
    table_id: str,
    title: str,
    df: pd.DataFrame,
    sources: list[str],
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    df = df.replace([np.inf, -np.inf], np.nan)
    csv_path = out_dir / f"{table_id}.csv"
    tex_path = out_dir / f"{table_id}.tex"
    df.to_csv(csv_path, index=False)
    if df.empty:
        tex = f"% {title}\n% Empty table: no source rows found.\n"
    else:
        tex = df.to_latex(
            index=False,
            longtable=True,
            escape=True,
            caption=title,
            label=f"tab:{table_id}",
        )
    tex_path.write_text(tex, encoding="utf-8")
    records.append(
        TableRecord(
            table_id=table_id,
            title=title,
            rows=int(len(df)),
            csv_path=rel(csv_path),
            tex_path=rel(tex_path),
            sources=existing_sources(sources),
        )
    )


def build_split_audit(records: list[TableRecord], out_dir: Path) -> None:
    source = "results/tkde_benchmark_split_manifest.csv"
    df = read_csv(source)
    if df.empty:
        save_table(records, out_dir, "table_01_split_leakage_audit", "Frozen split leakage audit", df, [source])
        return
    metrics = [
        "n_rows",
        "n_train",
        "n_test",
        "train_retention_after_purge",
        "target_overlap_rate",
        "exact_ligand_overlap_rate",
        "scaffold_overlap_rate",
        "exact_pair_overlap_rate",
        "pair_scaffold_overlap_rate",
        "family_overlap_rate",
        "document_source_overlap_rate",
        "assay_source_overlap_rate",
        "train_max_year",
        "test_min_year",
    ]
    present = [m for m in metrics if m in df.columns]
    grouped = (
        df.groupby(["benchmark", "split_mode"], dropna=False)
        .agg(
            n_splits=("split_file", "count"),
            n_seeds=("seed", "nunique"),
            n_folds=("fold", "nunique"),
            **{m: (m, "mean") for m in present},
        )
        .reset_index()
    )
    save_table(
        records,
        out_dir,
        "table_01_split_leakage_audit",
        "Frozen split leakage audit",
        compact_numeric(grouped),
        [source],
    )


def build_acnet_main(records: list[TableRecord], out_dir: Path) -> None:
    source = "results/acnet_target_family_best_joint_seed0_4_allfold_n10_summary.csv"
    df = read_csv(source)
    table = class_table(
        df,
        id_cols=["split_mode", "variant"],
        group_cols=["split_mode"],
        extra_mean_cols=[
            "n_train",
            "n_test",
            "train_retention_after_purge",
            "exact_ligand_overlap_rate",
            "scaffold_overlap_rate",
            "family_overlap_rate",
        ],
    )
    save_table(records, out_dir, "table_02_acnet_hard_ood_main", "ACNet hard-OOD main results", table, [source])


def build_raw_main(records: list[TableRecord], out_dir: Path) -> None:
    sources = [
        ("chembl_raw", "results/chembl_raw_source_temporal_pairs_seed0_4_n40_summary.csv"),
        ("bindingdb_raw", "results/bindingdb_raw_source_pairs_seed0_4_n80_summary.csv"),
    ]
    frames = []
    for label, path in sources:
        df = read_csv(path)
        if not df.empty:
            df.insert(0, "source", label)
            frames.append(df)
    raw = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    table = class_table(
        raw,
        id_cols=["source", "split_mode", "variant"],
        group_cols=["source", "split_mode"],
        extra_mean_cols=[
            "n_train",
            "n_test",
            "train_retention_after_purge",
            "target_overlap_rate",
            "exact_ligand_overlap_rate",
            "scaffold_overlap_rate",
            "document_source_overlap_rate",
            "assay_source_overlap_rate",
        ],
    )
    save_table(
        records,
        out_dir,
        "table_03_raw_source_temporal_main",
        "Raw ChEMBL and BindingDB source/temporal results",
        table,
        [p for _, p in sources],
    )


def build_lohi(records: list[TableRecord], out_dir: Path) -> None:
    source = "results/tkde_lohi_leakage_splits_seed0_2_summary.csv"
    df = read_csv(source)
    table = class_table(
        df,
        id_cols=["source", "split_mode", "variant"],
        group_cols=["source", "split_mode"],
        extra_mean_cols=[
            "n_train",
            "n_test",
            "train_retention_after_purge",
            "target_overlap_rate",
            "exact_ligand_overlap_rate",
            "scaffold_overlap_rate",
            "document_source_overlap_rate",
        ],
    )
    save_table(records, out_dir, "table_04_lohi_datasail_style", "LoHi/DataSAIL-style leakage-aware splits", table, [source])


def build_remaining_p0_p1(records: list[TableRecord], out_dir: Path) -> None:
    source = "results/tkde_remaining_p0_p1_seed0_2_summary.csv"
    df = read_csv(source)
    table = class_table(
        df,
        id_cols=["experiment", "split_mode", "variant"],
        group_cols=["experiment", "split_mode"],
        metrics=CLASS_METRICS + RANKING_METRICS,
        extra_mean_cols=[
            "n_train",
            "n_test",
            "train_retention_after_purge",
            "target_overlap_rate",
            "exact_ligand_overlap_rate",
            "scaffold_overlap_rate",
            "document_source_overlap_rate",
            "assay_source_overlap_rate",
            "cutoff_year",
        ],
    )
    save_table(
        records,
        out_dir,
        "table_05_crossdb_prospective_topk",
        "Cross-database, prospective, and top-k utility results",
        table,
        [source],
    )


def build_reliability(records: list[TableRecord], out_dir: Path) -> None:
    selective_sources = [
        "results/tkde_raw_reliability_seed0_2_selective_summary.csv",
        "results/acnet_selective_prediction_seed0_4_n10_summary.csv",
    ]
    frames = []
    for path in selective_sources:
        df = read_csv(path)
        if df.empty:
            continue
        df = df.copy()
        if "source" not in df.columns:
            df.insert(0, "source", "acnet")
        frames.append(df)
    selective = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    selective_table = class_table(
        selective,
        id_cols=["source", "split_mode", "variant", "coverage_target"],
        metrics=["coverage", "precision", "recall", "f1", "positive_precision_lift", "n_selected"],
        extra_mean_cols=["n_test", "train_retention_after_purge", "exact_ligand_overlap_rate", "scaffold_overlap_rate"],
        add_deltas=False,
    )
    save_table(
        records,
        out_dir,
        "table_06a_selective_prediction",
        "Selective prediction reliability",
        selective_table,
        selective_sources,
    )

    conformal_sources = [
        "results/tkde_raw_reliability_seed0_2_conformal_summary.csv",
        "results/acnet_conformal_prediction_seed0_4_n10_summary.csv",
    ]
    frames = []
    for path in conformal_sources:
        df = read_csv(path)
        if df.empty:
            continue
        df = df.copy()
        if "source" not in df.columns:
            df.insert(0, "source", "acnet")
        frames.append(df)
    conformal = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    conformal_table = class_table(
        conformal,
        id_cols=["source", "split_mode", "variant", "alpha"],
        metrics=[
            "coverage",
            "target_coverage",
            "coverage_gap",
            "avg_set_size",
            "positive_singleton_precision",
            "positive_singleton_recall",
            "positive_singleton_f1",
            "singleton_accuracy",
        ],
        extra_mean_cols=["n_test", "train_retention_after_purge", "exact_ligand_overlap_rate", "scaffold_overlap_rate"],
        add_deltas=False,
    )
    save_table(
        records,
        out_dir,
        "table_06b_conformal_prediction",
        "Conformal prediction reliability",
        conformal_table,
        conformal_sources,
    )


def build_noise(records: list[TableRecord], out_dir: Path) -> None:
    source = "results/tkde_noise_sensitivity_seed0_2_summary.csv"
    df = read_csv(source)
    table = class_table(
        df,
        id_cols=["source", "condition", "split_mode", "variant"],
        group_cols=["source", "condition", "split_mode"],
        extra_mean_cols=[
            "n_train",
            "n_test",
            "train_retention_after_purge",
            "test_positive_rate",
            "exact_ligand_overlap_rate",
            "scaffold_overlap_rate",
            "document_source_overlap_rate",
            "assay_source_overlap_rate",
        ],
    )
    save_table(records, out_dir, "table_07_noise_censoring_sensitivity", "Raw label-noise and censoring sensitivity", table, [source])


def build_non_neural(records: list[TableRecord], out_dir: Path) -> None:
    sources = [
        "results/tkde_non_neural_baselines_chembl_seed0_4_summary.csv",
        "results/tkde_non_neural_baselines_bindingdb_seed0_4_summary.csv",
    ]
    frames = [read_csv(path) for path in sources]
    frames = [df for df in frames if not df.empty]
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    table = class_table(
        df,
        id_cols=["benchmark", "split_mode", "variant"],
        group_cols=["benchmark", "split_mode"],
        extra_mean_cols=[
            "n_train",
            "n_test",
            "train_seconds",
            "predict_seconds",
            "train_retention_after_purge",
            "exact_ligand_overlap_rate",
            "scaffold_overlap_rate",
            "document_source_overlap_rate",
        ],
    )
    save_table(records, out_dir, "table_08_non_neural_baselines", "Non-neural baseline stress tests", table, sources)


def build_external_controls(records: list[TableRecord], out_dir: Path) -> None:
    moleculeace_source = "results/moleculeace_target_family_pairs_seed0_4_n120_summary.csv"
    mol = read_csv(moleculeace_source)
    mol_table = class_table(
        mol,
        id_cols=["split_mode", "variant"],
        group_cols=["split_mode"],
        extra_mean_cols=["n_train", "n_test", "target_overlap_rate", "exact_ligand_overlap_rate", "scaffold_overlap_rate"],
    )
    save_table(records, out_dir, "table_09a_moleculeace_external", "MoleculeACE external target-family controls", mol_table, [moleculeace_source])

    mtp_split_source = "results/mtpnet_csv_xgb_seed0_4.jsonl"
    mtp_cold_sources = [
        "results/mtpnet_cold_target_xgb_seed0_4.jsonl",
        "results/mtpnet_cold_target_saprot_xgb_seed0_4.jsonl",
    ]
    split = read_jsonl(mtp_split_source)
    if not split.empty:
        split = split.copy()
        split.insert(0, "protocol", "mtpnet_csv_split")
    cold_frames = []
    for path in mtp_cold_sources:
        df = read_jsonl(path)
        if not df.empty:
            df = df.copy()
            df.insert(0, "protocol", "cold_target")
            cold_frames.append(df)
    frames = []
    if not split.empty:
        frames.append(split)
    frames.extend(cold_frames)
    mtp_raw = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    mtp_summary = regression_summary(mtp_raw, ["protocol", "variant"], REGRESSION_METRICS + ["n_train", "n_test", "seconds"])
    mtp_table = regression_table(
        mtp_summary,
        id_cols=["protocol", "variant"],
        metrics=REGRESSION_METRICS + ["n_train", "n_test", "seconds"],
    )
    save_table(records, out_dir, "table_09b_mtpnet_controls", "MTPNet split and cold-target controls", mtp_table, [mtp_split_source] + mtp_cold_sources)


def build_neural_and_pretraining(records: list[TableRecord], out_dir: Path) -> None:
    neural_source = "results/acnet_target_family_neural_controls_full25_summary.csv"
    neural = read_csv(neural_source)
    neural_table = class_table(
        neural,
        id_cols=["split_mode", "variant"],
        group_cols=["split_mode"],
        extra_mean_cols=["n_train", "n_test", "train_retention_after_purge", "exact_ligand_overlap_rate", "scaffold_overlap_rate"],
    )
    save_table(records, out_dir, "table_10a_acnet_neural_controls", "ACNet neural and shuffled target controls", neural_table, [neural_source])

    tac_source = "results/acnet_target_assay_strengthening_current_no_shufpairs_summary.csv"
    tac = read_csv(tac_source)
    tac_table = class_table(
        tac,
        id_cols=["split_mode", "variant"],
        group_cols=["split_mode"],
        metrics=["roc_auc", "pr_auc", "mcc", "balanced_accuracy"],
    )
    save_table(records, out_dir, "table_10b_target_assay_pretraining_controls", "Target/assay pretraining controls", tac_table, [tac_source])

    few_source = "results/acnet_few_pretrained_representation_summary.csv"
    few = read_csv(few_source)
    few_table = class_table(
        few,
        id_cols=[
            "benchmark",
            "subset",
            "split_mode",
            "variant",
            "feature_family",
            "representation",
            "pair_feature_mode",
        ],
        group_cols=["benchmark", "subset", "split_mode", "feature_family", "representation", "pair_feature_mode"],
        metrics=["roc_auc", "pr_auc", "mcc", "balanced_accuracy", "f1", "label_score_spearman"],
        extra_mean_cols=["test_positive_rate", "target_overlap_rate"],
    )
    save_table(records, out_dir, "table_10c_few_pretrained_representations", "Few-shot pretrained representation controls", few_table, [few_source])


def build_statistics(records: list[TableRecord], out_dir: Path) -> None:
    source = "results/tkde_statistical_effects_audit.csv"
    df = read_csv(source)
    if not df.empty:
        keep = [
            "experiment",
            "source",
            "condition",
            "split_mode",
            "variant",
            "baseline",
            "metric",
            "n",
            "mean_delta_favorable",
            "ci95_low",
            "ci95_high",
            "win_rate",
            "cliffs_delta",
            "wilcoxon_p",
            "holm_p",
        ]
        df = compact_numeric(df[[c for c in keep if c in df.columns]])
    save_table(records, out_dir, "table_11_statistical_effects", "Paired statistical effect-size audit", df, [source])


def build_runtime(records: list[TableRecord], out_dir: Path) -> None:
    source = "results/tkde_runtime_scalability.csv"
    df = compact_numeric(read_csv(source))
    save_table(records, out_dir, "table_12_runtime_scalability", "Runtime and scalability profile", df, [source])


def build_failure_and_diagnostics(records: list[TableRecord], out_dir: Path) -> None:
    tax_source = "results/tkde_remaining_p0_p1_seed0_2_target_failure_taxonomy.csv"
    tax = compact_numeric(read_csv(tax_source))
    save_table(records, out_dir, "table_13_target_failure_taxonomy", "Target-level failure taxonomy", tax, [tax_source])

    cal_source = "results/tkde_remaining_p0_p1_seed0_2_calibration_bins.csv"
    cal = compact_numeric(read_csv(cal_source))
    save_table(records, out_dir, "table_14a_p0_p1_calibration_bins", "Cross-database/prospective calibration bins", cal, [cal_source])

    risk_source = "results/tkde_remaining_p0_p1_seed0_2_risk_coverage.csv"
    risk = compact_numeric(read_csv(risk_source))
    save_table(records, out_dir, "table_14b_p0_p1_risk_coverage", "Cross-database/prospective risk-coverage curves", risk, [risk_source])


def build_completion_matrix(records: list[TableRecord], out_dir: Path) -> None:
    expected = [
        ("P0", "Frozen leakage-audited benchmark splits", "results/tkde_benchmark_split_manifest.csv", "done"),
        ("P0", "ACNet target-family/cold-ligand hard OOD", "results/acnet_target_family_best_joint_seed0_4_allfold_n10_summary.csv", "done"),
        ("P0", "MoleculeACE external target-family controls", "results/moleculeace_target_family_pairs_seed0_4_n120_summary.csv", "done"),
        ("P0", "MTPNet same-source split and cold-target controls", "results/mtpnet_cold_target_analysis.md", "done"),
        ("P0", "Raw ChEMBL source/temporal benchmark", "results/chembl_raw_source_temporal_pairs_seed0_4_n40_summary.csv", "done"),
        ("P0", "Raw BindingDB source-purged benchmark", "results/bindingdb_raw_source_pairs_seed0_4_n80_summary.csv", "done"),
        ("P0", "LoHi/DataSAIL-style leakage-aware split", "results/tkde_lohi_leakage_splits_seed0_2_summary.csv", "done"),
        ("P0", "Raw external selective/conformal reliability", "results/tkde_raw_reliability_seed0_2_conformal_summary.csv", "done"),
        ("P0", "Cross-database transfer", "results/tkde_remaining_p0_p1_seed0_2_summary.csv", "done"),
        ("P0", "Prospective temporal validation", "results/tkde_remaining_p0_p1_seed0_2_summary.csv", "done"),
        ("P0", "Top-k enrichment and BEDROC utility", "results/tkde_remaining_p0_p1_seed0_2_summary.csv", "done"),
        ("P0", "Statistical effect sizes and multiplicity", "results/tkde_statistical_effects_audit.csv", "done"),
        ("P1", "Noise/censoring sensitivity", "results/tkde_noise_sensitivity_seed0_2_summary.csv", "done"),
        ("P1", "Non-neural baseline stress tests", "results/tkde_non_neural_baselines_chembl_seed0_4_summary.csv", "done"),
        ("P1", "Neural and shuffled target controls", "results/acnet_target_family_neural_controls_full25_summary.csv", "done"),
        ("P1", "Target/assay pretraining controls", "results/acnet_target_assay_strengthening_current_no_shufpairs_summary.csv", "done"),
        ("P1", "Success/failure case studies", "results/tkde_failure_case_studies.csv", "done"),
        ("P2", "Leaderboard snapshot and schema", "results/tkde_leaderboard_snapshot.csv", "done"),
        ("P2", "Runtime and scalability", "results/tkde_runtime_scalability.csv", "done"),
        ("P2", "Official DataSAIL optimizer exact split", "results/official_datasail_exact_split.csv", "blocked_env"),
    ]
    rows = []
    for priority, experiment, path, planned_status in expected:
        exists = repo_path(path).exists()
        status = "done" if exists and planned_status != "blocked_env" else planned_status
        detail = ""
        if planned_status == "blocked_env":
            detail = "official DataSAIL requires Python >=3.9,<3.13; a separate Python 3.12 environment was created for follow-up testing, but the package install did not complete within the audit run, so exact DataSAIL is not included in this frozen table package"
        rows.append(
            {
                "priority": priority,
                "experiment": experiment,
                "status": status,
                "artifact": path,
                "artifact_exists": bool(exists),
                "detail": detail,
            }
        )
    table = pd.DataFrame(rows)
    save_table(records, out_dir, "table_15_experiment_completion_matrix", "P0/P1/P2 experiment completion matrix", table, [row[2] for row in expected])


def write_index(records: list[TableRecord], out_dir: Path) -> None:
    df = pd.DataFrame([record.__dict__ for record in records])
    df.to_csv(out_dir / "tables_index.csv", index=False)
    lines = [
        "# TKDE Experiment Tables Package",
        "",
        "This package contains table artifacts only: compact CSV files and matching LaTeX longtables.",
        "",
        "| table_id | rows | csv | tex | sources |",
        "|---|---:|---|---|---|",
    ]
    for record in records:
        lines.append(
            f"| {record.table_id} | {record.rows} | `{record.csv_path}` | `{record.tex_path}` | {record.sources} |"
        )
    (out_dir / "tables_index.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        default="results/tkde_experiment_tables_2026_05_20",
        help="Directory for compact CSV and LaTeX table artifacts.",
    )
    args = parser.parse_args()

    out_dir = repo_path(args.output_dir)
    records: list[TableRecord] = []
    build_split_audit(records, out_dir)
    build_acnet_main(records, out_dir)
    build_raw_main(records, out_dir)
    build_lohi(records, out_dir)
    build_remaining_p0_p1(records, out_dir)
    build_reliability(records, out_dir)
    build_noise(records, out_dir)
    build_non_neural(records, out_dir)
    build_external_controls(records, out_dir)
    build_neural_and_pretraining(records, out_dir)
    build_statistics(records, out_dir)
    build_runtime(records, out_dir)
    build_failure_and_diagnostics(records, out_dir)
    build_completion_matrix(records, out_dir)
    write_index(records, out_dir)
    print({"output_dir": rel(out_dir), "tables": len(records), "rows": sum(r.rows for r in records)})


if __name__ == "__main__":
    main()
