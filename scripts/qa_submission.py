from __future__ import annotations

import csv
import hashlib
import re
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ERRORS: list[str] = []


def require(condition: bool, message: str) -> None:
    if not condition:
        ERRORS.append(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_checksums() -> int:
    checksum_path = ROOT / "checksums_sha256.txt"
    rows = 0
    for line in checksum_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        require(match is not None, f"Malformed checksum line: {line}")
        if match is None:
            continue
        expected, relative = match.groups()
        path = ROOT / relative
        require(path.is_file(), f"Checksum target is missing: {relative}")
        if path.is_file():
            require(sha256(path) == expected, f"Checksum mismatch: {relative}")
        rows += 1
    return rows


def verify_file_manifest() -> int:
    with (ROOT / "MANIFEST.csv").open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    declared: set[str] = set()
    for row in rows:
        relative = row["path"]
        path = ROOT / relative
        require(relative not in declared, f"Duplicate manifest path: {relative}")
        require(path.is_file(), f"Manifest target is missing: {relative}")
        if path.is_file():
            require(path.stat().st_size == int(row["bytes"]), f"Manifest size mismatch: {relative}")
        declared.add(relative)
    return len(rows)


def verify_result_manifest() -> int:
    with (ROOT / "RESULT_MANIFEST.csv").open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        result_id = row["result_id"]
        for relative in (item.strip() for item in row["final_outputs"].split(";")):
            if not relative or relative == "none" or "*" in relative or "<" in relative:
                continue
            require((ROOT / relative).exists(), f"Missing final output for {result_id}: {relative}")
        for relative in re.findall(r"python\s+([^;\s]+\.py)", row["command"]):
            require((ROOT / relative).is_file(), f"Missing command entry point for {result_id}: {relative}")
    return len(rows)


def verify_control_contract() -> tuple[int, int]:
    table_path = ROOT / "review_artifact" / "tables" / "representation_selection_42.csv"
    with table_path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))

    expected_sources = {"ChEMBL", "BindingDB"}
    expected_splits = {"random", "target-cluster", "target-cluster+scaffold+source"}
    expected_representations = {
        "ECFP-XGB",
        "PC-CFRL-full",
        "PC-CFRL-scalar",
        "RF-scalars",
        "ExtraTrees-scalars",
        "ChemBERTa-sym-XGB",
        "MolFormer-sym-XGB",
    }
    require(len(rows) == 42, f"Representation matrix has {len(rows)} rows, expected 42")
    require({row["source"] for row in rows} == expected_sources, "Unexpected representation-matrix sources")
    require({row["split"] for row in rows} == expected_splits, "Unexpected representation-matrix split labels")
    require(
        {row["representation"] for row in rows} == expected_representations,
        "Unexpected representation-matrix methods",
    )
    counts = Counter((row["source"], row["split"]) for row in rows)
    require(
        all(counts[(source, split)] == 7 for source in expected_sources for split in expected_splits),
        "Each source/split cell must contain seven representations",
    )
    require(all(int(row["runs"]) == 25 for row in rows), "Every representation row must contain 25 runs")

    config = (ROOT / "configs" / "reproduce_trans_top_journal_experiments.sh").read_text(encoding="utf-8")
    unified = (ROOT / "scripts" / "run_tkde_unified_baseline_matrix.py").read_text(encoding="utf-8")
    raw = (ROOT / "scripts" / "run_chembl_raw_source_temporal_pairs.py").read_text(encoding="utf-8")
    require("scripts/run_tkde_unified_baseline_matrix.py" in config, "Unified control sweep is absent from the full sequence")
    require(
        "--split-modes random,target_family,family_scaffold_source_purged" in config,
        "Unified control sweep does not request the strict internal contract",
    )
    require(
        'default="random,target_family,family_scaffold_source_purged"' in unified,
        "Unified control script default does not include the strict internal contract",
    )
    require(
        "raw.family_folds" in unified and "raw.purge_train_mask" in unified,
        "Unified controls do not use target-cluster folds followed by purge logic",
    )
    require(
        'if "scaffold" in split_mode' in raw and 'if "source" in split_mode' in raw,
        "Strict split implementation does not expose both scaffold and source purges",
    )
    strict_rows = sum(row["split"] == "target-cluster+scaffold+source" for row in rows)
    require(strict_rows == 14, f"Strict display contract has {strict_rows} rows, expected 14")
    return len(rows), strict_rows


def main() -> None:
    checksum_rows = verify_checksums()
    manifest_rows = verify_file_manifest()
    result_rows = verify_result_manifest()
    matrix_rows, strict_rows = verify_control_contract()
    if ERRORS:
        for error in ERRORS:
            print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
    print("submission_qa: PASS")
    print(f"checksums_verified: {checksum_rows}")
    print(f"manifest_entries_verified: {manifest_rows}")
    print(f"result_manifest_rows_verified: {result_rows}")
    print(f"representation_rows_verified: {matrix_rows}")
    print(f"strict_contract_rows_verified: {strict_rows}")
    print("contract_mapping: family_scaffold_source_purged -> target-cluster+scaffold+source")


if __name__ == "__main__":
    main()
