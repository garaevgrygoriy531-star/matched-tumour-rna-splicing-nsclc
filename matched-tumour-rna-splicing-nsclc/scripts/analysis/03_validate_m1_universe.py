#!/usr/bin/env python3
"""Streaming validation for the v2.3 M1 RNA universe."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pandas as pd


V23_ROOT = Path(__file__).resolve().parents[1]
TABLE = V23_ROOT / "results" / "tables" / "V23_RARE_D2_NULL_UNIVERSE.tsv.gz"
AUDIT = V23_ROOT / "results" / "audit"


def main() -> int:
    status_path = AUDIT / "V23_M1_STATUS.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    cols = [
        "histology", "variant_uid", "gene", "junction_id", "cluster_id", "patient",
        "carrier_status", "target_count", "cluster_total", "junction_usage",
        "functional_potential_status", "R2_FUNCTIONALLY_ELIGIBLE",
    ]
    rows = 0
    pair_keys: set[tuple[str, ...]] = set()
    histology_counts: dict[str, int] = {}
    carrier_counts: dict[str, int] = {}
    function_counts: dict[str, int] = {}
    eligibility_counts: dict[str, int] = {}
    target_gt_total = 0
    target_negative = 0
    usage_outside = 0
    target_na = 0
    min_total = None
    max_total = None

    for chunk in pd.read_csv(TABLE, sep="\t", compression="gzip", usecols=cols, chunksize=250_000, low_memory=False):
        rows += len(chunk)
        pair_keys.update(map(tuple, chunk[["histology", "variant_uid", "gene", "junction_id", "cluster_id"]].itertuples(index=False, name=None)))
        for value, count in chunk["histology"].value_counts(dropna=False).items():
            histology_counts[str(value)] = histology_counts.get(str(value), 0) + int(count)
        for value, count in chunk["carrier_status"].value_counts(dropna=False).items():
            carrier_counts[str(value)] = carrier_counts.get(str(value), 0) + int(count)
        for value, count in chunk["functional_potential_status"].value_counts(dropna=False).items():
            function_counts[str(value)] = function_counts.get(str(value), 0) + int(count)
        for value, count in chunk["R2_FUNCTIONALLY_ELIGIBLE"].value_counts(dropna=False).items():
            eligibility_counts[str(value)] = eligibility_counts.get(str(value), 0) + int(count)
        target_na += int(chunk["target_count"].isna().sum())
        target_gt_total += int((chunk["target_count"] > chunk["cluster_total"]).sum())
        target_negative += int((chunk["target_count"] < 0).sum())
        usage_outside += int(((chunk["junction_usage"] < 0) | (chunk["junction_usage"] > 1)).sum())
        chunk_min = float(chunk["cluster_total"].min())
        chunk_max = float(chunk["cluster_total"].max())
        min_total = chunk_min if min_total is None else min(min_total, chunk_min)
        max_total = chunk_max if max_total is None else max(max_total, chunk_max)

    checks = [
        ("TABLE_EXISTS", True, TABLE.is_file()),
        ("ROWS_MATCH_STATUS", int(status["sample_level_rows"]), rows),
        ("PAIR_COUNT_MATCH_STATUS", int(status["raw_rare_d2_candidate_pairs"]), len(pair_keys)),
        ("TARGET_COUNT_NOT_NULL", 0, target_na),
        ("TARGET_COUNT_LE_CLUSTER_TOTAL", 0, target_gt_total),
        ("TARGET_COUNT_NONNEGATIVE", 0, target_negative),
        ("JUNCTION_USAGE_IN_0_1", 0, usage_outside),
        ("FUNCTIONAL_STATUS_HOLD_ONLY", ["HOLD_MISSING_GENCODE50_MANE_SEQUENCE_CONTEXT"], sorted(function_counts)),
        ("FUNCTIONAL_ELIGIBILITY_HOLD_ONLY", ["HOLD"], sorted(eligibility_counts)),
        ("MIN_CLUSTER_TOTAL", 10, min_total),
    ]
    rows_out = []
    for check_id, expected, observed in checks:
        if isinstance(expected, float) or isinstance(observed, float):
            passed = abs(float(expected) - float(observed)) < 1e-9
        else:
            passed = expected == observed
        rows_out.append({"check_id": check_id, "expected": json.dumps(expected, ensure_ascii=False), "observed": json.dumps(observed, ensure_ascii=False), "status": "PASS" if passed else "FAIL"})
    validation_pass = all(row["status"] == "PASS" for row in rows_out)
    with (AUDIT / "V23_M1_VALIDATION.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check_id", "expected", "observed", "status"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)

    status.update({
        "validation_status": "PASS" if validation_pass else "FAIL",
        "validation_checks_failed": sum(row["status"] == "FAIL" for row in rows_out),
        "validated_rows": rows,
        "validated_candidate_pairs": len(pair_keys),
        "histology_row_counts": histology_counts,
        "carrier_status_counts": carrier_counts,
        "functional_status_counts": function_counts,
        "eligibility_counts": eligibility_counts,
        "min_cluster_total": min_total,
        "max_cluster_total": max_total,
    })
    status_path.write_text(json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(status, indent=2, ensure_ascii=False))
    return 0 if validation_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
