from __future__ import annotations

from pathlib import Path

import pandas as pd


BASE = Path("results")
TABLE_DIR = BASE / "tkde_experiment_tables_2026_05_20"


def latex_escape(value: object) -> str:
    text = str(value)
    return (
        text.replace("\\", "\\textbackslash{}")
        .replace("_", "\\_")
        .replace("%", "\\%")
        .replace("&", "\\&")
    )


def fmt_pm(mean: float, std: float) -> str:
    return f"{mean:.3f} +/- {std:.3f}"


def main() -> None:
    chem = pd.read_csv(BASE / "official_datasail_exact_split_chembl.csv")
    bind = pd.read_csv(BASE / "official_datasail_exact_split_bindingdb_source_purged.csv")
    combined = pd.concat([chem, bind], ignore_index=True)
    combined.to_csv(BASE / "official_datasail_exact_split.csv", index=False)

    fold = pd.concat(
        [
            pd.read_csv(BASE / "official_datasail_exact_split_chembl_fold_metrics.csv"),
            pd.read_csv(BASE / "official_datasail_exact_split_bindingdb_source_purged_fold_metrics.csv"),
        ],
        ignore_index=True,
    )
    fold.to_csv(BASE / "official_datasail_exact_split_fold_metrics.csv", index=False)

    deltas = pd.concat(
        [
            pd.read_csv(BASE / "official_datasail_exact_split_chembl_paired_deltas.csv"),
            pd.read_csv(BASE / "official_datasail_exact_split_bindingdb_source_purged_paired_deltas.csv"),
        ],
        ignore_index=True,
    )
    deltas.to_csv(BASE / "official_datasail_exact_split_paired_deltas.csv", index=False)

    assignments = pd.concat(
        [
            pd.read_csv(BASE / "official_datasail_exact_split_chembl_assignments.csv"),
            pd.read_csv(BASE / "official_datasail_exact_split_bindingdb_source_purged_assignments.csv"),
        ],
        ignore_index=True,
    )
    assignments.to_csv(BASE / "official_datasail_exact_split_assignments.csv", index=False)

    variant_map = {
        "ecfp_absdiff_xgb": "ECFP abs-diff XGB",
        "rich_diff_scalars": "PC-CFRL rich-diff",
        "rich_scalars_only": "PC-CFRL scalar-only",
    }
    source_map = {"chembl_raw": "ChEMBL raw", "bindingdb_raw": "BindingDB raw"}
    split_map = {
        "official_datasail_c2_ecfp_source_purged": "Official DataSAIL C2/ECFP + source purge"
    }

    rows: list[dict[str, object]] = []
    for (_source, _split), sub in combined.groupby(["source", "split_mode"], sort=False):
        base_row = sub[sub["variant"].eq("ecfp_absdiff_xgb")]
        base_roc = float(base_row["roc_auc_mean"].iloc[0])
        base_pr = float(base_row["pr_auc_mean"].iloc[0])
        base_mcc = float(base_row["mcc_mean"].iloc[0])
        for _, r in sub.iterrows():
            is_base = r["variant"] == "ecfp_absdiff_xgb"
            rows.append(
                {
                    "Source": source_map.get(r["source"], r["source"]),
                    "Split": split_map.get(r["split_mode"], r["split_mode"]),
                    "Variant": variant_map.get(r["variant"], r["variant"]),
                    "Folds": int(r.get("roc_auc_count", r.get("n_test_count", 0))),
                    "Train": f"{r['n_train_mean']:.0f}",
                    "Test": f"{r['n_test_mean']:.0f}",
                    "ROC-AUC": fmt_pm(float(r["roc_auc_mean"]), float(r["roc_auc_std"])),
                    "PR-AUC": fmt_pm(float(r["pr_auc_mean"]), float(r["pr_auc_std"])),
                    "MCC": fmt_pm(float(r["mcc_mean"]), float(r["mcc_std"])),
                    "Delta ROC": "--" if is_base else f"{float(r['roc_auc_mean']) - base_roc:+.3f}",
                    "Delta PR": "--" if is_base else f"{float(r['pr_auc_mean']) - base_pr:+.3f}",
                    "Delta MCC": "--" if is_base else f"{float(r['mcc_mean']) - base_mcc:+.3f}",
                    "Ligand ov.": f"{float(r['exact_ligand_overlap_rate_mean']):.3f}",
                    "Scaffold ov.": f"{float(r['scaffold_overlap_rate_mean']):.3f}",
                    "Doc ov.": f"{float(r['document_source_overlap_rate_mean']):.3f}",
                    "Assay ov.": f"{float(r['assay_source_overlap_rate_mean']):.3f}",
                }
            )

    compact = pd.DataFrame(rows)
    TABLE_DIR.mkdir(exist_ok=True)
    compact.to_csv(TABLE_DIR / "table_16_official_datasail_exact_split.csv", index=False)

    caption = (
        "Official DataSAIL C2/ECFP exact split audit with the paper source-purge rule. "
        "Missing document/assay identifiers are treated as empty source sets; BindingDB fold2 "
        "is omitted by the prespecified minimum-test threshold."
    )
    cols = list(compact.columns)
    tex = [
        "\\begin{table*}[!t]",
        "\\centering",
        "\\scriptsize",
        "\\caption{" + caption + "}",
        "\\label{tab:supp_official_datasail}",
        "\\begin{tabular}{lllrrrrrrrrrrrrr}",
        "\\toprule",
        " & ".join(cols).replace("Delta", "$\\Delta$") + " \\\\",
        "\\midrule",
    ]
    for _, r in compact.iterrows():
        vals = [latex_escape(r[c]).replace("+/-", "$\\pm$") for c in cols]
        tex.append(" & ".join(vals) + " \\\\")
    tex.extend(["\\bottomrule", "\\end{tabular}", "\\end{table*}"])
    (TABLE_DIR / "table_16_official_datasail_exact_split.tex").write_text("\n".join(tex) + "\n")

    matrix_csv = TABLE_DIR / "table_15_experiment_completion_matrix.csv"
    if matrix_csv.exists():
        matrix = pd.read_csv(matrix_csv)
        mask = matrix["experiment"].astype(str).eq("Official DataSAIL optimizer exact split")
        if mask.any():
            matrix.loc[mask, "status"] = "done"
            matrix.loc[mask, "artifact"] = "results/official_datasail_exact_split.csv"
            matrix.loc[mask, "artifact_exists"] = True
            if "detail" in matrix.columns:
                matrix.loc[
                    mask,
                    "detail",
                ] = (
                    "official DataSAIL 1.3.0 C2/ECFP source-purged split audit completed in "
                    "Python 3.12; missing source identifiers are normalized to empty source sets"
                )
            matrix.loc[
                mask,
                "notes",
            ] = (
                "official DataSAIL 1.3.0 C2/ECFP source-purged split audit completed in "
                "Python 3.12; not-selected interactions are excluded from model train/test "
                "and BindingDB fold2 is skipped by the minimum-test threshold"
            )
        matrix.to_csv(matrix_csv, index=False)

        tex2 = [
            "\\begin{table*}[!t]",
            "\\centering",
            "\\scriptsize",
            "\\caption{STable 20: Experiment completion matrix.}",
            "\\label{tab:supp_completion_matrix}",
            "\\begin{tabular}{llllll}",
            "\\toprule",
            "Priority & Experiment & Status & Artifact & Artifact exists & Notes \\\\",
            "\\midrule",
        ]
        for _, r in matrix.iterrows():
            vals = [latex_escape(r.get(c, "")) for c in ["priority", "experiment", "status", "artifact", "artifact_exists", "notes"]]
            tex2.append(" & ".join(vals) + " \\\\")
        tex2.extend(["\\bottomrule", "\\end{tabular}", "\\end{table*}"])
        (TABLE_DIR / "table_15_experiment_completion_matrix.tex").write_text("\n".join(tex2) + "\n")

    index_csv = TABLE_DIR / "tables_index.csv"
    if index_csv.exists():
        index = pd.read_csv(index_csv)
        row = {
            "table_id": "table_16_official_datasail_exact_split",
            "title": "Official DataSAIL exact split audit",
            "csv": "table_16_official_datasail_exact_split.csv",
            "tex": "table_16_official_datasail_exact_split.tex",
            "source": "results/official_datasail_exact_split.csv",
        }
        index = index[index["table_id"].astype(str) != row["table_id"]]
        index = pd.concat([index, pd.DataFrame([row])], ignore_index=True)
        index.to_csv(index_csv, index=False)
        lines = ["# TKDE experiment table index", ""]
        for _, r in index.iterrows():
            lines.append(
                f"- {r['table_id']}: {r.get('title', '')} ({r.get('csv', '')}, {r.get('tex', '')})"
            )
        (TABLE_DIR / "tables_index.md").write_text("\n".join(lines) + "\n")

    report = [
        "# Official DataSAIL Exact Split Audit",
        "",
        (
            "This follow-up uses official DataSAIL 1.3.0 C2 with ECFP molecular similarity, "
            "then applies the paper source-purge rule. Missing document/assay source identifiers "
            "are normalized to empty source sets so absent assay IDs are not treated as shared "
            "leakage tokens. DataSAIL not-selected interactions are excluded from model training "
            "and testing."
        ),
        "",
        "## Summary",
        compact.to_markdown(index=False),
        "",
        "## Feasibility",
        pd.read_csv(BASE / "official_datasail_source_purge_feasibility.csv").to_markdown(index=False),
        "",
        "## Boundary",
        (
            "- ChEMBL source-purged official DataSAIL is retained as a negative/boundary audit: "
            "PC-CFRL scalar variants do not beat ECFP on mean ROC/MCC, though rich-diff is close "
            "on PR-AUC."
        ),
        (
            "- BindingDB source-purged official DataSAIL remains strongly positive: rich-diff "
            "improves ROC-AUC, PR-AUC, and MCC over ECFP across the four evaluable folds; fold2 "
            "is skipped because the DataSAIL-selected test fold has only 23 interactions."
        ),
        (
            "- These rows close the previous official-DataSAIL environment blocker but remain "
            "supplementary diagnostics, not replacements for the frozen target/source/temporal "
            "split contract used for the main claims."
        ),
    ]
    (BASE / "official_datasail_exact_split_report.md").write_text("\n".join(report) + "\n")

    print(compact.to_string(index=False))
    print(
        {
            "combined_rows": len(combined),
            "fold_rows": len(fold),
            "delta_rows": len(deltas),
            "assignments": len(assignments),
        }
    )


if __name__ == "__main__":
    main()
