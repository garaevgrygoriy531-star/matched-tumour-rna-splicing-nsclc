#!/usr/bin/env python3
"""Freeze the pre-existing GTEx lung context for the post-M2/M3 frozen R2 set."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


EXPECTED_N = 9
KEY = ("variant_uid", "gene", "junction_id", "cluster_id")
EXPECTED_CLASSES = {
    "NORMAL_ABSENT": 1,
    "NORMAL_RARE": 2,
    "NORMAL_LOW": 0,
    "PREEXISTING_NORMAL_ROUTE": 6,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"missing header: {path}")
        return list(reader.fieldnames), list(reader)


def require_columns(columns: list[str], required: set[str], label: str) -> None:
    missing = sorted(required - set(columns))
    if missing:
        raise RuntimeError(f"{label} missing required columns: {missing}")


def row_key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row[field].strip() for field in KEY)


def route_class(detected: int, total: int, prevalence: float) -> str:
    if total <= 0 or detected < 0 or detected > total or not math.isfinite(prevalence):
        raise RuntimeError("invalid GTEx detected/total/prevalence values")
    expected = detected / total
    if not math.isclose(prevalence, expected, rel_tol=0.0, abs_tol=5e-15):
        raise RuntimeError(f"prevalence is inconsistent with detected/total: {prevalence} != {expected}")
    if detected == 0:
        if prevalence != 0.0:
            raise RuntimeError("zero detected samples must have zero prevalence")
        return "NORMAL_ABSENT"
    if prevalence <= 0.01:
        return "NORMAL_RARE"
    if prevalence <= 0.05:
        return "NORMAL_LOW"
    return "PREEXISTING_NORMAL_ROUTE"


def interpretation(route: str) -> str:
    if route in {"NORMAL_ABSENT", "NORMAL_RARE"}:
        return "CRYPTIC_OR_RARE_ROUTE_ACTIVATION"
    if route in {"NORMAL_LOW", "PREEXISTING_NORMAL_ROUTE"}:
        return "OVER_UTILIZATION_OF_PRE_EXISTING_NORMAL_SPLICE_ROUTE"
    raise RuntimeError(f"unexpected route class: {route}")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def self_test() -> int:
    assert route_class(0, 100, 0.0) == "NORMAL_ABSENT"
    assert route_class(1, 100, 0.01) == "NORMAL_RARE"
    assert route_class(2, 100, 0.02) == "NORMAL_LOW"
    assert route_class(5, 100, 0.05) == "NORMAL_LOW"
    assert route_class(6, 100, 0.06) == "PREEXISTING_NORMAL_ROUTE"
    print(json.dumps({"status": "SELF_TEST_PASS", "tests": ["zero_absent", "one_percent_inclusive_rare", "five_percent_inclusive_low", "above_five_preexisting"]}))
    return 0


def run(args: argparse.Namespace) -> int:
    for path in (args.frozen_keys, args.gtex_v22):
        if not path.is_file():
            raise RuntimeError(f"missing input: {path}")
    if sha256(args.frozen_keys) != args.expected_frozen_sha256.lower():
        raise RuntimeError("frozen9 SHA256 mismatch")
    if sha256(args.gtex_v22) != args.expected_gtex_sha256.lower():
        raise RuntimeError("GTEx v2.2 SHA256 mismatch")

    frozen_cols, frozen = read_tsv(args.frozen_keys)
    gtex_cols, gtex = read_tsv(args.gtex_v22)
    require_columns(frozen_cols, set(KEY) | {"histology", "selected_transcript"}, "frozen9")
    require_columns(gtex_cols, set(KEY) | {
        "R_class", "event_class", "gtex_lung_detected_samples", "gtex_lung_total_samples",
        "gtex_lung_prevalence", "gtex_lung_evaluable_cluster_samples_ge10",
        "gtex_lung_median_usage_evaluable", "gtex_lung_q90_usage_evaluable",
        "gtex_lung_q99_usage_evaluable", "normal_usage_status", "interpretation_guard",
    }, "GTEx v2.2")
    if len(frozen) != EXPECTED_N or len({row_key(row) for row in frozen}) != EXPECTED_N:
        raise RuntimeError("frozen key file is not exactly 9 unique events")
    gtex_index: dict[tuple[str, ...], dict[str, str]] = {}
    for row in gtex:
        key = row_key(row)
        if key in gtex_index:
            raise RuntimeError(f"duplicate GTEx key: {key}")
        gtex_index[key] = row

    out_rows: list[dict[str, object]] = []
    for frozen_row in frozen:
        key = row_key(frozen_row)
        if key not in gtex_index:
            raise RuntimeError(f"frozen event missing from GTEx v2.2: {key}")
        source = gtex_index[key]
        if source["normal_usage_status"] != "PASS":
            raise RuntimeError(f"GTEx source row did not pass: {key}")
        detected_float = float(source["gtex_lung_detected_samples"])
        total_float = float(source["gtex_lung_total_samples"])
        if not detected_float.is_integer() or not total_float.is_integer():
            raise RuntimeError(f"non-integer GTEx counts: {key}")
        detected, total = int(detected_float), int(total_float)
        prevalence = float(source["gtex_lung_prevalence"])
        cls = route_class(detected, total, prevalence)
        out_rows.append({
            "histology": frozen_row["histology"],
            "variant_uid": frozen_row["variant_uid"],
            "gene": frozen_row["gene"],
            "junction_id": frozen_row["junction_id"],
            "cluster_id": frozen_row["cluster_id"],
            "selected_transcript": frozen_row["selected_transcript"],
            "phase0_R_class": source["R_class"],
            "phase0_event_class": source["event_class"],
            "gtex_lung_detected_samples": detected,
            "gtex_lung_total_samples": total,
            "gtex_lung_prevalence": format(prevalence, ".17g"),
            "gtex_lung_evaluable_cluster_samples_ge10": source["gtex_lung_evaluable_cluster_samples_ge10"],
            "gtex_lung_median_usage_evaluable": source["gtex_lung_median_usage_evaluable"],
            "gtex_lung_q90_usage_evaluable": source["gtex_lung_q90_usage_evaluable"],
            "gtex_lung_q99_usage_evaluable": source["gtex_lung_q99_usage_evaluable"],
            "gtex_route_class": cls,
            "biological_route_language": interpretation(cls),
            "interpretation_guard": "GTEX_ABSENCE_OR_LOW_PREVALENCE_IS_NOT_PROOF_OF_TUMOR_SPECIFICITY",
            "step11_status": "PASS_REUSED_FROZEN_V22_GTEX_CONTEXT",
        })

    counts = Counter(str(row["gtex_route_class"]) for row in out_rows)
    if dict(counts) != {key: value for key, value in EXPECTED_CLASSES.items() if value > 0}:
        raise RuntimeError(f"unexpected Step11 class counts: {dict(counts)}")

    fields = list(out_rows[0])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output, out_rows, fields)
    audit = {
        "module": "STEP11_REUSE_GTEX_CONTEXT",
        "status": "PASS",
        "row_count": len(out_rows),
        "key_fields": list(KEY),
        "class_rule": {
            "NORMAL_ABSENT": "detected_samples == 0 and prevalence == 0",
            "NORMAL_RARE": "0 < prevalence <= 0.01",
            "NORMAL_LOW": "0.01 < prevalence <= 0.05",
            "PREEXISTING_NORMAL_ROUTE": "prevalence > 0.05",
        },
        "class_counts": EXPECTED_CLASSES,
        "classification_uses": ["gtex_lung_detected_samples", "gtex_lung_total_samples", "gtex_lung_prevalence"],
        "usage_fields_are_descriptive_only": True,
        "interpretation_boundary": "GTEx absence or rarity is normal-context annotation, not proof of tumor specificity",
        "inputs": {
            str(args.frozen_keys): sha256(args.frozen_keys),
            str(args.gtex_v22): sha256(args.gtex_v22),
        },
        "script_sha256": sha256(Path(__file__).resolve()),
        "output_sha256": sha256(args.output),
    }
    args.audit.write_text(json.dumps(audit, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    manifest_rows = [
        {"role": "script", "path": str(Path(__file__).resolve()), "sha256": audit["script_sha256"]},
        {"role": "input_frozen9", "path": str(args.frozen_keys), "sha256": sha256(args.frozen_keys)},
        {"role": "input_gtex_v22", "path": str(args.gtex_v22), "sha256": sha256(args.gtex_v22)},
        {"role": "output", "path": str(args.output), "sha256": sha256(args.output)},
        {"role": "audit", "path": str(args.audit), "sha256": sha256(args.audit)},
    ]
    write_tsv(args.manifest, manifest_rows, ["role", "path", "sha256"])
    print(json.dumps({"status": "PASS", "rows": len(out_rows), "class_counts": EXPECTED_CLASSES, "output_sha256": sha256(args.output)}, sort_keys=True))
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--self-test", action="store_true")
    result.add_argument("--frozen-keys", type=Path)
    result.add_argument("--gtex-v22", type=Path)
    result.add_argument("--expected-frozen-sha256")
    result.add_argument("--expected-gtex-sha256")
    result.add_argument("--output", type=Path)
    result.add_argument("--audit", type=Path)
    result.add_argument("--manifest", type=Path)
    return result


def main() -> int:
    args = parser().parse_args()
    if args.self_test:
        return self_test()
    required = (args.frozen_keys, args.gtex_v22, args.expected_frozen_sha256, args.expected_gtex_sha256, args.output, args.audit, args.manifest)
    if any(value is None for value in required):
        raise RuntimeError("formal execution requires all input, expected-hash, and output arguments")
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
