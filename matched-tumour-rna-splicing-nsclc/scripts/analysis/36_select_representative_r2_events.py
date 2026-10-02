#!/usr/bin/env python3
"""Deterministically select representative frozen R2 events for main figures."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


EXPECTED_INPUT_SHA256 = "41dd8701f60183dafc9f9543016dd17a47c99583e78289e45ff275140a108cd9"
EXPECTED_SELECTED = {
    "chr19:45420410:C:T",  # synonymous + RNA-added
    "chr19:1220373:G:T",   # splice-region + cryptic
    "chr13:48307406:G:T",  # missense + RNA-added
    "chr2:224481892:C:A",  # fully concordant comparator
}
FULL_KEY = ("histology", "variant_uid", "gene", "junction_id", "cluster_id", "selected_transcript")


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
            raise RuntimeError("missing evidence-matrix header")
        return list(reader.fieldnames), list(reader)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def eligible(row: dict[str, str]) -> tuple[bool, str]:
    checks = [
        (row["phase0_high_confidence"].lower() == "true", "NOT_PHASE0_HIGH_CONFIDENCE"),
        (row["phase0_R_class"].startswith("R2"), "NOT_R2"),
        (row["matrix_status"] == "PASS_NO_COMPOSITE_SCORE", "FINAL_MATRIX_NOT_PASS"),
        (row["m3_reciprocal_class"] in {"RECIPROCAL_STRONG", "RECIPROCAL_DIRECTIONAL"}, "NO_CLEAR_RECIPROCAL_EVIDENCE"),
        (row["altered_splice_type"] not in {"", "NA", "UNRESOLVED"}, "ALTERED_SPLICE_TYPE_UNRESOLVED"),
        (row["frame_consequence_rebuilt"] == "FRAMESHIFT_PTC", "FRAME_CONSEQUENCE_UNRESOLVED"),
        (row["nmd_geometry_status"] == "NMD_COMPATIBLE", "NMD_GEOMETRY_UNRESOLVED"),
        (row["protein_consequence_scope"] == "PREDICTED_ALTERED_ORF_IF_TRANSLATED", "PROTEIN_SCOPE_INVALID"),
    ]
    try:
        m2_p = float(row["m2_event_empirical_p_one_sided"])
        checks.append((math.isfinite(m2_p) and 0.0 <= m2_p <= 1.0, "M2_SUPPORT_INVALID"))
    except ValueError:
        checks.append((False, "M2_SUPPORT_INVALID"))
    failed = [reason for passed, reason in checks if not passed]
    return not failed, "PASS" if not failed else "|".join(failed)


def priority(row: dict[str, str]) -> tuple[int, int, str, str]:
    reciprocal = 0 if row["m3_reciprocal_class"] == "RECIPROCAL_STRONG" else 1
    legacy = 0 if row["legacy_sensitivity_class"] == "LEGACY_CONCORDANT" else 1
    return reciprocal, legacy, row["variant_uid"], row["junction_id"]


def choose(candidates: list[dict[str, str]], used: set[str], predicate, slot: str) -> dict[str, str]:
    available = [row for row in candidates if row["variant_uid"] not in used and predicate(row)]
    if not available:
        raise RuntimeError(f"no eligible candidate for slot {slot}")
    selected = sorted(available, key=priority)[0]
    used.add(selected["variant_uid"])
    return selected


def run(args: argparse.Namespace) -> int:
    if not args.input.is_file():
        raise RuntimeError(f"missing input: {args.input}")
    if sha256(args.input) != EXPECTED_INPUT_SHA256:
        raise RuntimeError("final evidence matrix SHA256 mismatch")
    columns, rows = read_tsv(args.input)
    required = set(FULL_KEY) | {
        "phase0_R_class", "phase0_high_confidence", "D2_subclass", "m2_event_empirical_p_one_sided",
        "m3_reciprocal_class", "altered_splice_type", "frame_consequence_rebuilt",
        "nmd_geometry_status", "protein_consequence_scope", "predictor_consensus_d50",
        "legacy_sensitivity_class", "matrix_status",
    }
    missing = sorted(required - set(columns))
    if missing:
        raise RuntimeError(f"input missing columns: {missing}")
    if len(rows) != 9 or len({tuple(row[field] for field in FULL_KEY) for row in rows}) != 9:
        raise RuntimeError("input is not exactly nine unique frozen events")

    eligibility = {row["variant_uid"]: eligible(row) for row in rows}
    candidates = [row for row in rows if eligibility[row["variant_uid"]][0]]
    used: set[str] = set()
    selections: list[tuple[str, dict[str, str]]] = []
    selections.append(("SYNONYMOUS_RNA_ADDED", choose(
        candidates, used,
        lambda row: row["D2_subclass"] == "synonymous" and row["predictor_consensus_d50"] == "RNA_ADDS_GEOMETRY_INFORMATION",
        "SYNONYMOUS_RNA_ADDED",
    )))
    selections.append(("SPLICE_REGION_CRYPTIC", choose(
        candidates, used,
        lambda row: row["D2_subclass"] == "splice-region" and row["altered_splice_type"].startswith("CRYPTIC_"),
        "SPLICE_REGION_CRYPTIC",
    )))
    selections.append(("MISSENSE_RNA_ADDED", choose(
        candidates, used,
        lambda row: row["D2_subclass"] == "missense" and row["predictor_consensus_d50"] == "RNA_ADDS_GEOMETRY_INFORMATION",
        "MISSENSE_RNA_ADDED",
    )))
    selections.append(("FULLY_CONCORDANT_COMPARATOR", choose(
        candidates, used,
        lambda row: row["predictor_consensus_d50"] == "PREDICTOR_FULLY_CONCORDANT",
        "FULLY_CONCORDANT_COMPARATOR",
    )))
    if used != EXPECTED_SELECTED:
        raise RuntimeError(f"selection replay mismatch: {sorted(used)}")

    slot_by_variant = {row["variant_uid"]: (order, slot) for order, (slot, row) in enumerate(selections, start=1)}
    nonselection = {
        "chr12:25227342:T:A": "NO_RECIPROCAL_NOT_ELIGIBLE",
        "chr19:10489648:C:A": "REDUNDANT_FULLY_CONCORDANT_DIRECTIONAL_AFTER_CUL3",
        "chr6:43443117:G:T": "REDUNDANT_MISSENSE_RNA_ADDS_DIRECTIONAL_AFTER_RB1",
        "chr17:7673570:T:C": "REDUNDANT_MISSENSE_RNA_ADDS_LEGACY_DIFFERENT_AFTER_RB1",
        "chr17:7676379:T:A": "REDUNDANT_SPLICE_REGION_AND_FULLY_CONCORDANT_NO_ADDITIONAL_SUBTYPE",
    }
    out: list[dict[str, object]] = []
    for row in rows:
        selected = row["variant_uid"] in used
        order, slot = slot_by_variant.get(row["variant_uid"], ("NA", "NOT_SELECTED"))
        is_eligible, eligibility_reason = eligibility[row["variant_uid"]]
        reason = f"SELECTED_{slot}" if selected else nonselection.get(row["variant_uid"], eligibility_reason)
        out.append({
            **{field: row[field] for field in FULL_KEY},
            "selected_for_main_figure": str(selected).upper(),
            "selection_slot": slot,
            "figure_order": order,
            "selection_eligible": str(is_eligible).upper(),
            "eligibility_reason": eligibility_reason,
            "D2_subclass": row["D2_subclass"],
            "altered_splice_type": row["altered_splice_type"],
            "m3_reciprocal_class": row["m3_reciprocal_class"],
            "predictor_consensus_d50": row["predictor_consensus_d50"],
            "legacy_sensitivity_class": row["legacy_sensitivity_class"],
            "selection_or_nonselection_reason": reason,
            "tie_break_fields": "m3_class_priority|legacy_concordance_priority|variant_uid|junction_id",
            "used_gene_identity_for_selection": "FALSE",
            "used_m2_p_for_ranking": "FALSE",
            "used_effect_size_for_ranking": "FALSE",
            "composite_score_present": "FALSE",
            "selection_rule_version": "V23_STEP13_RULE_1.0",
        })
    out.sort(key=lambda row: (0 if row["selected_for_main_figure"] == "TRUE" else 1, int(row["figure_order"]) if row["figure_order"] != "NA" else 999, str(row["variant_uid"])))
    selected_rows = [row for row in out if row["selected_for_main_figure"] == "TRUE"]
    if len(selected_rows) != 4 or not (3 <= len(selected_rows) <= 5):
        raise RuntimeError("selected event count is outside the planned 3-5 range")
    if {row["D2_subclass"] for row in selected_rows} != {"synonymous", "missense", "splice-region"}:
        raise RuntimeError("selected events do not cover all required D2 subclasses")
    if not any(str(row["altered_splice_type"]).startswith("CRYPTIC_") for row in selected_rows):
        raise RuntimeError("selected events lack a cryptic splice subtype")
    if not any(row["predictor_consensus_d50"] == "RNA_ADDS_GEOMETRY_INFORMATION" for row in selected_rows):
        raise RuntimeError("selected events lack an RNA-added example")
    if not any(row["predictor_consensus_d50"] == "PREDICTOR_FULLY_CONCORDANT" for row in selected_rows):
        raise RuntimeError("selected events lack a fully concordant comparator")

    write_tsv(args.output, out, list(out[0]))
    audit = {
        "module": "STEP13_REPRESENTATIVE_EVENT_SELECTION",
        "status": "PASS",
        "candidate_count": len(rows),
        "selected_count": len(selected_rows),
        "selected_variant_uids_in_figure_order": [row["variant_uid"] for row in selected_rows],
        "selection_slots_in_figure_order": [row["selection_slot"] for row in selected_rows],
        "selection_scope": "MECHANISM_DISPLAY_ONLY_NOT_RECLASSIFICATION",
        "used_gene_identity_for_selection": False,
        "used_m2_p_for_ranking": False,
        "used_effect_size_for_ranking": False,
        "composite_score_present": False,
        "r3_rows": 0,
        "tie_break": ["m3 class priority", "legacy concordance priority", "variant_uid", "junction_id"],
        "interpretation_guards": [
            "STATISTICAL_LINK_NOT_MOLECULARLY_PHASED",
            "NMD_COMPATIBLE_NOT_MEASURED_DECAY",
            "PREDICTED_ALTERED_ORF_IF_TRANSLATED",
            "DOMAIN_UNRESOLVED",
            "GTEX_ABSENT_OR_RARE_NOT_TUMOR_SPECIFIC",
            "RNA_ADDS_GEOMETRY_INFORMATION_NOT_PREDICTOR_FAILURE",
            "REGTOOLS_RECOVERY_NOT_FUNCTIONAL_VALIDATION",
        ],
        "input_sha256": sha256(args.input),
        "script_sha256": sha256(Path(__file__).resolve()),
        "output_sha256": sha256(args.output),
    }
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    args.audit.write_text(json.dumps(audit, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    manifest = [
        {"role": "script", "path": str(Path(__file__).resolve()), "sha256": audit["script_sha256"]},
        {"role": "input", "path": str(args.input), "sha256": sha256(args.input)},
        {"role": "output", "path": str(args.output), "sha256": sha256(args.output)},
        {"role": "audit", "path": str(args.audit), "sha256": sha256(args.audit)},
    ]
    write_tsv(args.manifest, manifest, ["role", "path", "sha256"])
    print(json.dumps({"status": "PASS", "selected": [row["variant_uid"] for row in selected_rows], "output_sha256": sha256(args.output)}, sort_keys=True))
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--input", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--audit", type=Path, required=True)
    result.add_argument("--manifest", type=Path, required=True)
    return result


if __name__ == "__main__":
    raise SystemExit(run(parser().parse_args()))
