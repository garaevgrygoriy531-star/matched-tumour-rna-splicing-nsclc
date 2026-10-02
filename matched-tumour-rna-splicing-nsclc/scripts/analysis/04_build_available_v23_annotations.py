#!/usr/bin/env python3
"""Materialize v2.3 annotations that are supported by frozen v2.2 assets.

No v2.2 field is overwritten.  Modules requiring missing reference exon,
expression, or predictor assets are represented as explicit HOLD/SKIPPED
rows rather than fabricated negative results.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import pandas as pd


V23_ROOT = Path(__file__).resolve().parents[1]
V22_ROOT = Path(os.environ.get("V22_ROOT", str(V23_ROOT.parent / "v2.2")))
V22_FREEZE = V22_ROOT / "final_freeze" / "v2.2"
TABLES = V23_ROOT / "results" / "tables"
AUDIT = V23_ROOT / "results" / "audit"
STATS = V23_ROOT / "results" / "stats"


def write_df(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def main() -> int:
    event = pd.read_csv(V22_FREEZE / "A_PHASE0_EVENT_TABLE_V22.tsv", sep="\t", dtype=str).fillna("")
    r2 = pd.read_csv(V22_FREEZE / "A_PHASE0_R2_UNMASKED_V22.tsv", sep="\t", dtype=str).fillna("")
    r3 = pd.read_csv(V22_FREEZE / "A_PHASE0_R3_RESCUE_V22.tsv", sep="\t", dtype=str).fillna("")

    # 79-event transition atlas: all fields are newly prefixed v23_* or
    # explicitly named phase0_* snapshots, preserving the frozen columns.
    atlas_cols = [
        "variant_uid", "gene", "junction_id", "histologies", "patients", "DNA_class",
        "primary_consequence", "event_class_final", "R_class_final", "high_confidence_final",
        "altered_splice_type", "frame_consequence", "nmd_status", "functional_consensus",
        "reference_protein_length_aa", "altered_protein_length_aa", "domain_consequence",
        "gtex_lung_context_class", "gtex_lung_prevalence", "regtools_match_class",
    ]
    atlas = event[atlas_cols].copy()
    atlas = atlas.rename(columns={
        "event_class_final": "phase0_event_class",
        "R_class_final": "phase0_R_class",
        "high_confidence_final": "phase0_high_confidence",
        "domain_consequence": "phase0_domain_consequence",
        "gtex_lung_context_class": "phase0_gtex_route_class",
    })
    atlas["v23_transition"] = (
        atlas["primary_consequence"].astype(str) + " -> "
        + atlas["altered_splice_type"].astype(str) + " -> "
        + atlas["frame_consequence"].astype(str) + " -> "
        + atlas["nmd_status"].astype(str) + " -> predicted_protein_length"
    )
    atlas["v23_status"] = "V22_FROZEN_TRANSITION_REPACKAGED_NO_CLASSIFICATION_CHANGE"
    write_df(atlas, TABLES / "V23_FUNCTIONAL_TRANSITION_ATLAS.tsv")

    route = r2[[
        "variant_uid", "gene", "junction_id", "histologies", "patients", "event_class_final",
        "R_class_final", "gtex_lung_detected_samples", "gtex_lung_total_samples",
        "gtex_lung_prevalence", "gtex_lung_context_class", "gtex_lung_q99_usage_evaluable",
    ]].copy()
    route = route.rename(columns={"event_class_final": "phase0_event_class", "R_class_final": "phase0_R_class"})

    def route_mode(value: str) -> str:
        if value in {"NOT_DETECTED_IN_GTEX_LUNG", "RARE_IN_GTEX_LUNG_LE_1PCT"}:
            return "MODE_A_CRYPTIC_OR_RARE_ROUTE_ACTIVATION"
        if value in {"LOW_PREVALENCE_IN_GTEX_LUNG_LE_5PCT", "COMMON_IN_GTEX_LUNG_GT_5PCT"}:
            return "MODE_B_PREEXISTING_ROUTE_OVERUTILIZATION"
        return "UNRESOLVED"

    route["v23_route_mode"] = route["gtex_lung_context_class"].map(route_mode)
    route["v23_interpretation_boundary"] = "GTEx_absence_or_low_prevalence_is_not_tumor_specific_proof"
    route["v23_status"] = "REUSED_V22_GTEX_LUNG_CONTEXT"
    write_df(route, TABLES / "V23_GTEX_ROUTE_CONTEXT.tsv")

    # Exact transcript-protein/domain cross-reference is not present in the
    # local asset set.  Preserve the frozen unresolved state explicitly.
    domain_frames = []
    for label, frame in (("R2", r2), ("R3_TIER_B", r3)):
        out = frame[[
            "variant_uid", "gene", "junction_id", "histologies", "patients", "R_class_final",
            "selected_transcript", "reference_protein_length_aa", "altered_protein_length_aa",
            "domain_consequence", "frame_consequence", "nmd_status",
        ]].copy()
        out = out.rename(columns={"R_class_final": "phase0_R_class", "domain_consequence": "phase0_domain_consequence"})
        out["v23_postphase0_refinement"] = "DOMAIN_UNRESOLVED"
        out["domain_status"] = "UNRESOLVED"
        out["domain_reason"] = "NO_EXACT_TRANSCRIPT_PROTEIN_DOMAIN_CROSSREF_IN_LOCAL_ASSETS"
        out["v23_source_group"] = label
        domain_frames.append(out)
    write_df(pd.concat(domain_frames, ignore_index=True), TABLES / "V23_DOMAIN_REFINEMENT.tsv")
    r3_out = domain_frames[1].copy()
    r3_out["v23_postphase0_R3_status"] = "R3_TIER_B_REFINEMENT_UNRESOLVED"
    write_df(r3_out, TABLES / "V23_R3_POSTPHASE0_REFINEMENT.tsv")

    nmd = r2[["variant_uid", "gene", "junction_id", "histologies", "patients", "R_class_final", "nmd_status"]].copy()
    nmd = nmd.rename(columns={"R_class_final": "phase0_R_class"})
    nmd["v23_nmd_expression_status"] = "SKIPPED_NO_MATCHED_EXPRESSION"
    nmd["v23_nmd_interpretation"] = "NMD_consistency_not_assessed; no matched TCGA expression matrix present"
    write_df(nmd, TABLES / "V23_NMD_EXPRESSION_CONSISTENCY.tsv")

    # Explicit module holds: these are not empty result tables and are not
    # used as statistical evidence.
    hold_rows = [
        {"module": "M2_EMPIRICAL_NULL", "status": "HOLD_MISSING_INDEPENDENT_FUNCTIONAL_POTENTIAL", "reason": "GENCODE50/MANE/GRCh38 sequence assets unavailable locally; M2 cannot be run on outcome-selected proxies."},
        {"module": "M3_RECIPROCAL_SPLICE", "status": "HOLD_MISSING_REFERENCE_EXON_CHAIN", "reason": "Frozen output contains altered chain but not complete reference exon chain needed to predefine canonical competitors."},
        {"module": "M4_M5_PREDICTORS", "status": "HOLD_PREDICTOR_ENVIRONMENT_AND_REFERENCE_ASSETS", "reason": "SpliceAI/Pangolin frozen environments/models and GRCh38 reference are not available locally."},
        {"module": "M6_MATCHED_CONTROLS", "status": "HOLD_UPSTREAM_PREDICTOR_FEATURES", "reason": "Control matching requires frozen predictor geometry features."},
        {"module": "M8_NMD_EXPRESSION", "status": "SKIPPED_NO_MATCHED_EXPRESSION", "reason": "No matched TCGA expression matrix was found in v2.2 assets."},
        {"module": "M9_DOMAIN", "status": "HOLD_EXACT_CROSSREF_MISSING", "reason": "No exact transcript-protein-domain cross-reference asset is available locally."},
        {"module": "M11_CLINICAL_CONTEXT", "status": "DEFERRED_CORE_MODULES_INCOMPLETE", "reason": "v2.3 requires M2-M10 before low-priority clinical context."},
    ]
    with (AUDIT / "V23_MODULE_STATUS.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["module", "status", "reason"], delimiter="\t")
        writer.writeheader()
        writer.writerows(hold_rows)

    (AUDIT / "V23_AVAILABLE_ANNOTATIONS_STATUS.json").write_text(json.dumps({
        "transition_atlas_rows": int(len(atlas)),
        "gtex_r2_rows": int(len(route)),
        "domain_r2_rows": int(len(r2)),
        "domain_r3_rows": int(len(r3)),
        "modules": hold_rows,
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"transition_atlas_rows": len(atlas), "gtex_r2_rows": len(route), "domain_r2_rows": len(r2), "domain_r3_rows": len(r3), "status": "AVAILABLE_ANNOTATIONS_WRITTEN"}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
