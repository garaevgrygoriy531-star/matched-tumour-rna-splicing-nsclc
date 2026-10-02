#!/usr/bin/env python3
"""Build the final nine-event R2 evidence matrix without composite scoring."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


EXPECTED_N = 9
FULL_KEY = ("histology", "variant_uid", "gene", "junction_id", "cluster_id", "selected_transcript")
EVENT_KEY = ("variant_uid", "gene", "junction_id")
EXPECTED_SHA256 = {
    "frozen": "a587b8ecf8c6c286b6b285f2cab99041ed5cbb5b68d565c3762200c6bc844c81",
    "phase0": "66e2204a29c973259e6e216bf774530d52453f9884f700f99a101753deab2d97",
    "m2_event": "f415d30c460cf8fc3ad93c0a93e7894105d158aef92bdc24ddf7a5454010c90d",
    "m2_global": "3dec99ee61548b89ece911c621a82aae1be679f2966790d5ff274ae87184661a",
    "m2s": "5b5b35a415aa97ecb26330d7728336eb78a4d2cbccccae98c6fa77b0c2c59498",
    "m3": "a4af130e32c018d77394e7494fc230d4b520ac806e0a2071510c249aaed7395c",
    "predictor": "079e8f78a76b822843e1e8fde21acc88a6b87d14d1a5a8edaaa6a057cbdb4b8a",
    "gtex": "86098b6466f47c61f3e4e4c0c06412610fefe1749d64150f34ab1bd7ad89ea88",
    "protein": "f8c18933d4cb86414217373a4fd4914dfa2ffc7149c14c89610c014ba619994d",
    "m3_global_sidecar": "564aae4599783edbe4a788b59717fea81782b248e2d3f98aaeb9dced14ab2e5e",
    "m6_global_sidecar": "055ac6f0c791da90be03fcf906f983a1af06019c93082e78d3d613d8058e67eb",
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
            raise RuntimeError(f"missing TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def require_columns(columns: list[str], required: set[str], label: str) -> None:
    missing = sorted(required - set(columns))
    if missing:
        raise RuntimeError(f"{label} missing columns: {missing}")


def index_unique(rows: list[dict[str, str]], fields: tuple[str, ...], label: str) -> dict[tuple[str, ...], dict[str, str]]:
    result: dict[tuple[str, ...], dict[str, str]] = {}
    for row in rows:
        key = tuple(row[field].strip() for field in fields)
        if key in result:
            raise RuntimeError(f"duplicate {label} key: {key}")
        result[key] = row
    return result


def frozen_key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row[field].strip() for field in FULL_KEY)


def phase0_key(row: dict[str, str]) -> tuple[str, ...]:
    return (
        row["histologies"].strip(), row["variant_uid"].strip(), row["gene"].strip(),
        row["junction_id"].strip(), row["cluster_id"].strip(), row["selected_transcript"].strip(),
    )


def d2_subclass(consequence: str) -> str:
    value = consequence.lower()
    if "synonymous_variant" in value:
        return "synonymous"
    if "missense_variant" in value:
        return "missense"
    if any(token in value for token in ("splice_region", "splice_donor_region", "polypyrimidine_tract")):
        return "splice-region"
    if any(token in value for token in ("inframe_insertion", "inframe_deletion", "protein_altering_variant")):
        return "inframe/exon-edge"
    raise RuntimeError(f"unsupported frozen D2 consequence: {consequence}")


def one_row(rows: list[dict[str, str]], label: str) -> dict[str, str]:
    if len(rows) != 1:
        raise RuntimeError(f"{label} must contain exactly one row")
    return rows[0]


def require_equal(label: str, *values: str) -> None:
    normalized = {value.strip() for value in values}
    if len(normalized) != 1:
        raise RuntimeError(f"cross-module mismatch for {label}: {values}")


def run(args: argparse.Namespace) -> int:
    inputs = {
        "frozen": args.frozen,
        "phase0": args.phase0,
        "m2_event": args.m2_event,
        "m2_global": args.m2_global,
        "m2s": args.m2s,
        "m3": args.m3,
        "predictor": args.predictor,
        "gtex": args.gtex,
        "protein": args.protein,
        "m3_global_sidecar": args.m3_global_sidecar,
        "m6_global_sidecar": args.m6_global_sidecar,
    }
    for label, path in inputs.items():
        if not path.is_file():
            raise RuntimeError(f"missing {label} input: {path}")
        actual = sha256(path)
        if actual != EXPECTED_SHA256[label]:
            raise RuntimeError(f"{label} SHA256 mismatch: {actual} != {EXPECTED_SHA256[label]}")

    frozen_cols, frozen = read_tsv(args.frozen)
    phase0_cols, phase0_rows = read_tsv(args.phase0)
    m2_cols, m2_rows = read_tsv(args.m2_event)
    _, m2_global_rows = read_tsv(args.m2_global)
    _, m2s_rows = read_tsv(args.m2s)
    m3_cols, m3_rows = read_tsv(args.m3)
    predictor_cols, predictor_rows = read_tsv(args.predictor)
    gtex_cols, gtex_rows = read_tsv(args.gtex)
    protein_cols, protein_rows = read_tsv(args.protein)

    require_columns(frozen_cols, set(FULL_KEY), "frozen9")
    require_columns(phase0_cols, {
        "histologies", "variant_uid", "gene", "junction_id", "cluster_id", "selected_transcript",
        "patients", "DNA_class", "primary_consequence", "R_class_final", "event_class_final",
        "high_confidence_final", "regtools_match_class", "regtools_known_recovered",
        "regtools_same_histology_exact_event", "regtools_exact_event_min_adjusted_p",
        "v26_annotation_status", "legacy_version_sensitive", "legacy_sensitivity_class",
        "post_legacy_event_class",
    }, "phase0")
    require_columns(m2_cols, set(FULL_KEY) | {
        "event_empirical_p_one_sided", "event_bh_q_within_frozen_9", "event_empirical_p_basis",
        "q_scope", "q_interpretation",
    }, "M2 event")
    require_columns(m3_cols, set(FULL_KEY) | {"classification", "aberrant_delta", "canonical_delta", "m3_status"}, "M3")
    require_columns(predictor_cols, set(EVENT_KEY) | {
        "selected_transcript", "spliceai_d50_status", "spliceai_d50_score",
        "spliceai_d50_relevant_score", "spliceai_d50_geometry_class", "pangolin_d50_status",
        "pangolin_d50_score", "pangolin_d50_signed_delta", "pangolin_d50_geometry_class",
        "consensus_d50_class", "interpretation_boundary",
    }, "predictor")
    require_columns(gtex_cols, set(FULL_KEY) | {
        "gtex_lung_prevalence", "gtex_route_class", "biological_route_language",
        "interpretation_guard", "step11_status",
    }, "GTEx")
    require_columns(protein_cols, set(EVENT_KEY) | {
        "histology", "transcript_id", "altered_splice_type", "frame_consequence_rebuilt",
        "nmd_geometry_status", "predicted_length_fraction", "exact_reference_residue_fraction_retained",
        "ptc_position_aa_1based", "domain_status", "domain_reason", "consequence_scope",
        "protein_expression_status", "row_status",
    }, "protein")
    if len(frozen) != EXPECTED_N or len({frozen_key(row) for row in frozen}) != EXPECTED_N:
        raise RuntimeError("frozen input is not exactly nine unique events")

    phase0_index = index_unique(phase0_rows, ("histologies", "variant_uid", "gene", "junction_id", "cluster_id", "selected_transcript"), "phase0")
    m2_index = index_unique(m2_rows, FULL_KEY, "M2")
    m3_index = index_unique(m3_rows, FULL_KEY, "M3")
    predictor_index = index_unique(predictor_rows, EVENT_KEY, "predictor")
    gtex_index = index_unique(gtex_rows, FULL_KEY, "GTEx")
    protein_index = index_unique(protein_rows, ("histology", "variant_uid", "gene", "junction_id"), "protein")
    m2_global = one_row(m2_global_rows, "M2 global")
    m2s = one_row(m2s_rows, "M2-S")
    if m2s["ROBUST"] != "YES":
        raise RuntimeError("M2-S robustness gate did not pass")

    output: list[dict[str, object]] = []
    for skeleton in frozen:
        full = frozen_key(skeleton)
        event = tuple(skeleton[field] for field in EVENT_KEY)
        protein_key = (skeleton["histology"],) + event
        try:
            p0 = phase0_index[full]
            m2 = m2_index[full]
            m3 = m3_index[full]
            pred = predictor_index[event]
            gtex = gtex_index[full]
            protein = protein_index[protein_key]
        except KeyError as exc:
            raise RuntimeError(f"MISSING_JOIN for frozen event {full}: {exc}") from exc
        require_equal("selected transcript", skeleton["selected_transcript"], m2["selected_transcript"], m3["selected_transcript"], pred["selected_transcript"], gtex["selected_transcript"], protein["transcript_id"])
        require_equal("phase0 R class", p0["R_class_final"], gtex["phase0_R_class"], protein["phase0_R_class"])
        require_equal("phase0 event class", p0["event_class_final"], gtex["phase0_event_class"], protein["phase0_event_class"])
        if p0["DNA_class"] != "D2" or p0["high_confidence_final"].lower() != "true":
            raise RuntimeError(f"non-HC D2 event entered final matrix: {full}")
        if m3["m3_status"] != "RECIPROCAL_DIRECTION_CONCORDANT_DESCRIPTIVE" or gtex["step11_status"] != "PASS_REUSED_FROZEN_V22_GTEX_CONTEXT" or protein["row_status"] != "PASS":
            raise RuntimeError(f"upstream module status failed for {full}")
        if pred["spliceai_d50_status"] != "SCORED" or pred["pangolin_d50_status"] != "SCORED":
            raise RuntimeError(f"unscored predictor event: {full}")
        relevant_score = pred["spliceai_d50_relevant_score"].strip()
        output.append({
            "histology": skeleton["histology"],
            "variant_uid": skeleton["variant_uid"],
            "gene": skeleton["gene"],
            "junction_id": skeleton["junction_id"],
            "cluster_id": skeleton["cluster_id"],
            "selected_transcript": skeleton["selected_transcript"],
            "patient": p0["patients"],
            "phase0_R_class": p0["R_class_final"],
            "phase0_event_class": p0["event_class_final"],
            "phase0_high_confidence": p0["high_confidence_final"],
            "DNA_class": p0["DNA_class"],
            "D2_subclass": d2_subclass(p0["primary_consequence"]),
            "primary_DNA_consequence": p0["primary_consequence"],
            "m2_event_empirical_p_one_sided": m2["event_empirical_p_one_sided"],
            "m2_event_bh_q_within_frozen9": m2["event_bh_q_within_frozen_9"],
            "m2_event_p_basis": m2["event_empirical_p_basis"],
            "m2_event_scope": m2["q_scope"],
            "m2_event_interpretation": m2["q_interpretation"],
            "context_m2_global_empirical_p_one_sided": m2_global["empirical_p_one_sided"],
            "context_m2_global_observed_r2": m2_global["observed_frozen_r2_count"],
            "context_m2_global_null_mean": m2_global["null_mean"],
            "context_m2s_empirical_p_one_sided": m2s["empirical_p"],
            "context_m2s_robust": m2s["ROBUST"],
            "context_m2s_null_mean": m2s["null_mean"],
            "context_m2s_observed_expected_ratio": m2s["observed_expected_ratio"],
            "m3_reciprocal_class": m3["classification"],
            "m3_aberrant_delta": m3["aberrant_delta"],
            "m3_canonical_delta": m3["canonical_delta"],
            "m3_status": m3["m3_status"],
            "spliceai_d50_score": pred["spliceai_d50_score"],
            "spliceai_d50_relevant_score": relevant_score if relevant_score else "NA",
            "spliceai_d50_relevant_score_status": "AVAILABLE" if relevant_score else "NO_RELEVANT_GEOMETRY_SCORE",
            "spliceai_d50_geometry_class": pred["spliceai_d50_geometry_class"],
            "pangolin_d50_score": pred["pangolin_d50_score"],
            "pangolin_d50_signed_delta": pred["pangolin_d50_signed_delta"],
            "pangolin_d50_geometry_class": pred["pangolin_d50_geometry_class"],
            "predictor_consensus_d50": pred["consensus_d50_class"],
            "gtex_lung_prevalence": gtex["gtex_lung_prevalence"],
            "gtex_route_class": gtex["gtex_route_class"],
            "gtex_biological_route_language": gtex["biological_route_language"],
            "gtex_interpretation_guard": gtex["interpretation_guard"],
            "altered_splice_type": protein["altered_splice_type"],
            "frame_consequence_rebuilt": protein["frame_consequence_rebuilt"],
            "ptc_position_aa_1based": protein["ptc_position_aa_1based"],
            "nmd_geometry_status": protein["nmd_geometry_status"],
            "predicted_altered_orf_length_fraction": protein["predicted_length_fraction"],
            "exact_reference_residue_fraction_retained": protein["exact_reference_residue_fraction_retained"],
            "domain_status": protein["domain_status"],
            "domain_reason": protein["domain_reason"],
            "protein_consequence_scope": protein["consequence_scope"],
            "protein_expression_status": protein["protein_expression_status"],
            "regtools_match_class": p0["regtools_match_class"],
            "regtools_known_recovered": p0["regtools_known_recovered"],
            "regtools_same_histology_exact_event": p0["regtools_same_histology_exact_event"],
            "regtools_exact_event_min_adjusted_p": p0["regtools_exact_event_min_adjusted_p"],
            "v26_annotation_status": p0["v26_annotation_status"],
            "legacy_version_sensitive": p0["legacy_version_sensitive"],
            "legacy_sensitivity_class": p0["legacy_sensitivity_class"],
            "post_legacy_event_class": p0["post_legacy_event_class"],
            "matrix_interpretation_guards": "M2_EVENT_P_CONDITIONAL_WITHIN_FROZEN9|NMD_COMPATIBLE_NOT_MEASURED_DECAY|PROTEIN_CONSEQUENCE_IF_TRANSLATED|DOMAIN_UNRESOLVED_NOT_NO_EFFECT|GTEX_ABSENT_OR_RARE_NOT_TUMOR_SPECIFIC|PREDICTOR_LOW_NOT_PREDICTOR_FAILURE|REGTOOLS_RECOVERY_NOT_FUNCTIONAL_VALIDATION|V26_NOVEL_IS_ANNOTATION_VERSION_STATUS",
            "matrix_status": "PASS_NO_COMPOSITE_SCORE",
        })

    if len(output) != EXPECTED_N or len({tuple(str(row[field]) for field in FULL_KEY) for row in output}) != EXPECTED_N:
        raise RuntimeError("final output is not exactly nine unique frozen events")
    for context_field in (
        "context_m2_global_empirical_p_one_sided", "context_m2_global_observed_r2",
        "context_m2_global_null_mean", "context_m2s_empirical_p_one_sided",
        "context_m2s_robust", "context_m2s_null_mean", "context_m2s_observed_expected_ratio",
    ):
        if len({str(row[context_field]) for row in output}) != 1:
            raise RuntimeError(f"context field is not constant: {context_field}")
    fields = list(output[0])
    write_tsv(args.output, output, fields)

    audit = {
        "module": "STEP12_FINAL_R2_EVIDENCE_MATRIX",
        "status": "PASS",
        "row_count": len(output),
        "left_table": "V23_POST_M2_M3_FROZEN9_KEYS.tsv",
        "canonical_key": list(FULL_KEY),
        "join_contract": "restrict each module to frozen9 and validate one-to-one; MISSING_JOIN is fatal",
        "r3_rows": 0,
        "composite_score_present": False,
        "d500_fields_in_primary_matrix": False,
        "m6_pair_rows_joined": False,
        "global_sidecars_not_event_attributes": [str(args.m3_global_sidecar), str(args.m6_global_sidecar)],
        "protein_fraction_semantics": {
            "predicted_altered_orf_length_fraction": "altered ORF length divided by reference protein length, if translated",
            "exact_reference_residue_fraction_retained": "reference residues retained by exact codon provenance",
        },
        "input_sha256": {label: sha256(path) for label, path in inputs.items()},
        "script_sha256": sha256(Path(__file__).resolve()),
        "output_sha256": sha256(args.output),
    }
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    args.audit.write_text(json.dumps(audit, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    manifest_rows = [
        {"role": "script", "path": str(Path(__file__).resolve()), "sha256": audit["script_sha256"]},
        *({"role": f"input_{label}", "path": str(path), "sha256": sha256(path)} for label, path in inputs.items()),
        {"role": "output", "path": str(args.output), "sha256": sha256(args.output)},
        {"role": "audit", "path": str(args.audit), "sha256": sha256(args.audit)},
    ]
    write_tsv(args.manifest, manifest_rows, ["role", "path", "sha256"])
    print(json.dumps({"status": "PASS", "rows": len(output), "output_sha256": sha256(args.output)}, sort_keys=True))
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    for name in EXPECTED_SHA256:
        result.add_argument("--" + name.replace("_", "-"), dest=name, type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--audit", type=Path, required=True)
    result.add_argument("--manifest", type=Path, required=True)
    return result


if __name__ == "__main__":
    raise SystemExit(run(parser().parse_args()))
