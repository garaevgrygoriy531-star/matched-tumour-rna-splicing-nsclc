#!/usr/bin/env python3
"""Assemble v2.2 Phase-0 endpoints, QC table, and prespecified gate decision."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import pandas as pd


KEYS = ["variant_uid", "gene", "junction_id"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def norm_histology(value: object) -> str:
    return str(value).replace("TCGA-", "").upper()


def joined(values) -> str:
    clean = sorted({str(value) for value in values if str(value) not in {"", "None", "nan"}})
    return ";".join(clean) if clean else "NA"


def wilson(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float, float]:
    if total <= 0:
        return float("nan"), float("nan"), float("nan")
    p = successes / total
    denom = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denom
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denom
    low = 0.0 if successes == 0 else max(0.0, center - half)
    high = 1.0 if successes == total else min(1.0, center + half)
    return p, low, high


def consequence_group(value: object) -> str:
    text = str(value)
    if text == "synonymous_variant":
        return "synonymous"
    if text == "missense_variant":
        return "missense"
    if text in {"splice_region_variant", "splice_donor_region_variant", "splice_polypyrimidine_tract_variant"}:
        return "splice-region"
    if text in {"inframe_insertion", "inframe_deletion", "protein_altering_variant"}:
        return "inframe/exon-edge"
    return "other"


def read_audits(specifications: list[str]) -> tuple[list[dict], dict[str, dict]]:
    rows = []
    loaded = {}
    for item in specifications:
        label, raw_path = item.split("=", 1)
        path = Path(raw_path)
        payload = json.loads(path.read_text(encoding="utf-8"))
        loaded[label] = payload
        status = str(payload.get("status", "MISSING_STATUS"))
        rows.append({
            "metric": f"audit_{label}",
            "value": status,
            "status": "PASS" if status == "PASS" else "FAIL",
            "note": str(path),
        })
    return rows, loaded


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--evidence", required=True, type=Path)
    ap.add_argument("--reconstruction", required=True, type=Path)
    ap.add_argument("--recurrent", required=True, type=Path)
    ap.add_argument("--gtex", required=True, type=Path)
    ap.add_argument("--regtools", required=True, type=Path)
    ap.add_argument("--legacy", required=True, type=Path)
    ap.add_argument("--qc-audit", action="append", default=[])
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    evidence_all = pd.read_csv(args.evidence, sep="\t", compression="gzip", low_memory=False)
    evidence_all["histology"] = evidence_all["histology"].map(norm_histology)
    evidence_all["high_confidence_discordance"] = evidence_all["high_confidence_discordance"].map(parse_bool)
    event_seed_rows = []
    stable_columns = [
        "cluster_id", "DNA_class", "primary_consequence", "selected_transcript", "junction_class",
        "link_rule", "event_class", "event_class_reason", "functional_consensus", "altered_splice_type",
        "frame_consequence", "nmd_status", "domain_consequence", "variant_retention", "R_class", "R_reason",
        "high_confidence_discordance", "bh_fdr_within_histology", "robust_z", "support_type",
    ]
    for key, group in evidence_all.groupby(KEYS, sort=True, dropna=False):
        first = group.iloc[0]
        for column in ["DNA_class", "primary_consequence", "selected_transcript", "event_class", "R_class", "high_confidence_discordance"]:
            if group[column].dropna().astype(str).nunique() > 1:
                raise RuntimeError(f"event-level inconsistency for {key}: {column}")
        row = {column: first.get(column) for column in stable_columns}
        row.update({"variant_uid": key[0], "gene": key[1], "junction_id": key[2]})
        row["histologies"] = joined(group["histology"])
        row["patients"] = joined(group["patient"])
        row["patient_evidence_count"] = int(group["patient"].nunique())
        event_seed_rows.append(row)
    events = pd.DataFrame(event_seed_rows)

    reconstruction = pd.read_csv(args.reconstruction, sep="\t", compression="gzip", low_memory=False)
    reconstruction = reconstruction.drop_duplicates(KEYS)
    reconstruction_cols = KEYS + [
        "altered_exon_chain", "reference_protein_length_aa", "altered_protein_length_aa",
        "v26_annotation_status", "v50_annotation_status", "current_transcripts_supporting",
        "legacy_version_sensitive", "current_exact_coding_transcripts", "cis_status", "HGVSc", "HGVSp",
    ]
    events = events.merge(reconstruction[reconstruction_cols], on=KEYS, how="left", validate="one_to_one")

    gtex = pd.read_csv(args.gtex, sep="\t", low_memory=False)
    gtex_cols = KEYS + [
        "gtex_lung_detected_samples", "gtex_lung_total_samples", "gtex_lung_prevalence",
        "gtex_lung_evaluable_cluster_samples_ge10", "gtex_lung_median_usage_evaluable",
        "gtex_lung_q90_usage_evaluable", "gtex_lung_q99_usage_evaluable", "gtex_lung_context_class",
        "normal_usage_status",
    ]
    events = events.merge(gtex[gtex_cols], on=KEYS, how="left", validate="one_to_one")

    regtools = pd.read_csv(args.regtools, sep="\t", low_memory=False)
    regtools_cols = KEYS + [
        "regtools_match_class", "regtools_known_recovered", "regtools_same_histology_exact_event",
        "regtools_exact_event_cohorts", "regtools_exact_event_anchors", "regtools_exact_event_tags",
        "regtools_exact_event_min_adjusted_p",
    ]
    events = events.merge(regtools[regtools_cols], on=KEYS, how="left", validate="one_to_one")

    legacy = pd.read_csv(args.legacy, sep="\t", low_memory=False)
    legacy["automatic_downgrade_to_E4"] = legacy["automatic_downgrade_to_E4"].map(parse_bool)
    legacy_cols = KEYS + [
        "legacy_sensitivity_class", "legacy_selected_transcript_full_id", "legacy_selected_altered_splice_type",
        "legacy_selected_frame_consequence", "legacy_selected_nmd_status", "legacy_resolved_functional_categories",
        "automatic_downgrade_to_E4", "post_legacy_event_class",
    ]
    events = events.merge(legacy[legacy_cols], on=KEYS, how="left", validate="one_to_one")
    events["automatic_downgrade_to_E4"] = events["automatic_downgrade_to_E4"].fillna(False).map(parse_bool)
    events["event_class_final"] = events["event_class"]
    events.loc[events["automatic_downgrade_to_E4"], "event_class_final"] = "E4_LEGACY_FUNCTIONAL_CONFLICT"
    events["R_class_final"] = events["R_class"]
    events.loc[events["automatic_downgrade_to_E4"], "R_class_final"] = "R4_LEGACY_FUNCTIONAL_CONFLICT"
    events["high_confidence_final"] = events["high_confidence_discordance"] & ~events["automatic_downgrade_to_E4"]
    events["novelty_label"] = events["regtools_known_recovered"].map(parse_bool).map(
        {True: "KNOWN_RECOVERED_NOT_NOVEL", False: "NOT_EXACTLY_REPORTED_BY_REGTOOLS_NOT_PROOF_OF_NOVELTY"}
    )
    events = events.sort_values(["high_confidence_final", "R_class_final", "gene", "variant_uid"], ascending=[False, True, True, True])

    recurrent = pd.read_csv(args.recurrent, sep="\t", compression="gzip", low_memory=False)
    recurrent["histology"] = recurrent["histology"].map(norm_histology)
    evaluable = recurrent[
        recurrent["DNA_class"].isin(["D1", "D2"])
        & recurrent["n_carriers_evaluable"].fillna(0).ge(1)
    ].copy()
    evaluable = evaluable.sort_values("n_carriers_evaluable", ascending=False).drop_duplicates(["histology", "variant_uid"])
    evaluable["consequence_group"] = evaluable["primary_consequence"].map(consequence_group)

    high_hist = evidence_all[evidence_all["high_confidence_discordance"]].copy()
    high_hist = high_hist.drop_duplicates(["histology"] + KEYS)
    conflict_keys = set(map(tuple, legacy.loc[legacy["automatic_downgrade_to_E4"], KEYS].to_numpy()))
    high_hist = high_hist[~high_hist[KEYS].apply(tuple, axis=1).isin(conflict_keys)]

    endpoint_rows = []
    for histology in ["LUAD", "LUSC"]:
        denominator = int(evaluable[evaluable["histology"].eq(histology)].shape[0])
        numerator = int(high_hist[high_hist["histology"].eq(histology)].shape[0])
        rate, low, high_ci = wilson(numerator, denominator)
        endpoint_rows.append({
            "histology": histology, "numerator_high_confidence_events": numerator,
            "denominator_rna_evaluable_d1_d2_variant_histology_instances": denominator,
            "rate": rate, "wilson95_low": low, "wilson95_high": high_ci,
        })
    pooled_numerator = int(len(high_hist))
    pooled_denominator = int(len(evaluable))
    rate, low, high_ci = wilson(pooled_numerator, pooled_denominator)
    endpoint_rows.append({
        "histology": "COMBINED_DESCRIPTIVE_POOLED", "numerator_high_confidence_events": pooled_numerator,
        "denominator_rna_evaluable_d1_d2_variant_histology_instances": pooled_denominator,
        "rate": rate, "wilson95_low": low, "wilson95_high": high_ci,
    })
    endpoint = pd.DataFrame(endpoint_rows)

    r2 = events[events["R_class_final"].astype(str).str.startswith("R2") & events["high_confidence_final"]].copy()
    r2["D2_subclass"] = r2["primary_consequence"].map(consequence_group)
    r2_endpoint_rows = []
    for histology in ["LUAD", "LUSC", "COMBINED_DESCRIPTIVE_POOLED"]:
        denom_scope = evaluable if histology.startswith("COMBINED") else evaluable[evaluable["histology"].eq(histology)]
        num_scope = high_hist[high_hist["R_class"].astype(str).str.startswith("R2")]
        if not histology.startswith("COMBINED"):
            num_scope = num_scope[num_scope["histology"].eq(histology)]
        for category in ["synonymous", "missense", "splice-region", "inframe/exon-edge"]:
            denominator = int((denom_scope["DNA_class"].eq("D2") & denom_scope["consequence_group"].eq(category)).sum())
            numerator = int((num_scope["primary_consequence"].map(consequence_group) == category).sum())
            p, lo, hi = wilson(numerator, denominator)
            r2_endpoint_rows.append({
                "histology": histology, "D2_subclass": category, "R2_events": numerator,
                "rna_evaluable_D2_variants": denominator, "proportion": p,
                "wilson95_low": lo, "wilson95_high": hi,
            })
    r2_endpoint = pd.DataFrame(r2_endpoint_rows)

    r3 = events[events["R_class_final"].astype(str).str.startswith("R3")].copy()
    r3a = r3[r3["R_class_final"].astype(str).str.contains("TIER_A", case=False, na=False)]
    d1_total = int(evaluable["DNA_class"].eq("D1").sum())
    r3_rate, r3_low, r3_high = wilson(int(len(r3a)), d1_total)

    qc_rows, audits = read_audits(args.qc_audit)
    event_counts = events["event_class_final"].astype(str).str.extract(r"^(E[1-5])")[0].value_counts().to_dict()
    r_counts = events["R_class_final"].astype(str).str.extract(r"^(R[0-5])")[0].value_counts().to_dict()
    interpretable = sum(event_counts.get(label, 0) for label in ["E1", "E2", "E3"])
    interpretable_fraction = interpretable / len(events)
    high_events = events[events["high_confidence_final"]]
    high_genes = int(high_events["gene"].nunique())
    noncanonical = int((~high_events["primary_consequence"].isin(["splice_acceptor_variant", "splice_donor_variant"])).sum())
    patient_counter = Counter()
    for patients in high_events["patients"]:
        for patient in str(patients).split(";"):
            if patient and patient != "NA":
                patient_counter[patient] += 1
    max_patient_events = max(patient_counter.values(), default=0)
    max_patient_fraction = max_patient_events / len(high_events) if len(high_events) else float("nan")
    no_single_sample_driver = len(patient_counter) > 1 and max_patient_fraction < 0.5
    gtex_eligible = events[events["R_class_final"].astype(str).str.startswith(("R1", "R2", "R3", "R4"))]
    gtex_evaluable = int(gtex_eligible["normal_usage_status"].eq("PASS").sum())
    audit_pass = all(row["status"] == "PASS" for row in qc_rows)
    required_qc = [
        ("event_count", len(events), "PASS" if len(events) > 0 else "FAIL", "Strict68 event count is derived after panel restriction; must be non-zero"),
        ("interpretable_E1_E2_E3_fraction", interpretable_fraction, "PASS" if interpretable_fraction >= 0.80 else "FAIL", "Adequate event-level interpretability threshold >=0.80"),
        ("high_confidence_events", len(high_events), "PASS" if len(high_events) >= 8 else "FAIL", "GO_MAIN_STRONG threshold >=8"),
        ("high_confidence_genes", high_genes, "PASS" if high_genes >= 4 else "FAIL", "GO_MAIN_STRONG threshold >=4"),
        ("high_confidence_noncanonical_consequence_events", noncanonical, "PASS" if noncanonical >= 3 else "FAIL", "Not primary splice_donor/acceptor consequence; threshold >=3"),
        ("high_confidence_R2_events", len(r2), "PASS" if len(r2) >= 3 else "FAIL", "GO_MAIN_STRONG R2 threshold >=3"),
        ("unique_high_confidence_patients", len(patient_counter), "PASS" if no_single_sample_driver else "FAIL", "Single hypermutated-sample concentration guard"),
        ("max_high_confidence_events_per_patient", max_patient_events, "PASS" if no_single_sample_driver else "FAIL", f"Fraction={max_patient_fraction:.6g}"),
        ("legacy_functional_conflicts", int(events["automatic_downgrade_to_E4"].sum()), "PASS" if not events["automatic_downgrade_to_E4"].any() else "WARN", "Conflicts are automatically downgraded to E4"),
        ("GTEx_normal_usage_evaluable_R1_R4_events", gtex_evaluable, "PASS" if gtex_evaluable == len(gtex_eligible) else "FAIL", "GTEx is required for all Strict68 R1-R4 events; R5 unresolved events are outside this context step. Absence is not tumor-specific proof"),
        ("RegTools_high_confidence_exact_known_recovered", int(high_events["regtools_known_recovered"].map(parse_bool).sum()), "PASS", "Technical overlap only; recovered events are not novel"),
    ]
    qc_rows.extend({"metric": m, "value": v, "status": s, "note": n} for m, v, s, n in required_qc)
    qc = pd.DataFrame(qc_rows)
    qc_system_pass = audit_pass and not (qc["status"] == "FAIL").any()

    gate = {
        "high_confidence_events_ge8": len(high_events) >= 8,
        "genes_ge4": high_genes >= 4,
        "noncanonical_events_ge3": noncanonical >= 3,
        "r2_or_r3_rule": len(r2) >= 3 or (len(r3a) >= 1 and len(r2) >= 3),
        "not_single_hypermutated_sample": no_single_sample_driver,
        "qc_no_systematic_error": qc_system_pass,
        "event_interpretation_sufficient": interpretable_fraction >= 0.80,
    }
    if all(gate.values()):
        verdict = "GO_MAIN_STRONG"
    elif 4 <= len(high_events) <= 7 and high_genes >= 3 and qc_system_pass:
        verdict = "GO_MAIN_CONDITIONAL"
    else:
        verdict = "STOP_MAIN"

    event_path = args.outdir / "A_PHASE0_EVENT_TABLE_V22.tsv"
    r2_path = args.outdir / "A_PHASE0_R2_UNMASKED_V22.tsv"
    r3_path = args.outdir / "A_PHASE0_R3_RESCUE_V22.tsv"
    qc_path = args.outdir / "A_PHASE0_QC_V22.tsv"
    endpoint_path = args.outdir / "A_PHASE0_PRIMARY_ENDPOINTS_V22.tsv"
    r2_endpoint_path = args.outdir / "A_PHASE0_R2_PROPORTIONS_V22.tsv"
    decision_path = args.outdir / "A_PHASE0_DECISION_V22.md"
    events.to_csv(event_path, sep="\t", index=False, na_rep="NA")
    r2.to_csv(r2_path, sep="\t", index=False, na_rep="NA")
    r3.to_csv(r3_path, sep="\t", index=False, na_rep="NA")
    qc.to_csv(qc_path, sep="\t", index=False, na_rep="NA")
    endpoint.to_csv(endpoint_path, sep="\t", index=False, na_rep="NA")
    r2_endpoint.to_csv(r2_endpoint_path, sep="\t", index=False, na_rep="NA")

    endpoint_lines = [
        f"- {row.histology}: {int(row.numerator_high_confidence_events)}/{int(row.denominator_rna_evaluable_d1_d2_variant_histology_instances)} "
        f"= {100*row.rate:.3f}% (Wilson 95% CI {100*row.wilson95_low:.3f}–{100*row.wilson95_high:.3f}%)"
        for row in endpoint.itertuples(index=False)
    ]
    checklist = [f"- {'PASS' if passed else 'FAIL'} — {name}" for name, passed in gate.items()]
    decision_text = f"""# A v2.2 Phase-0 Decision

## Verdict

`{verdict}`

This is the prespecified Phase-0 gate result. Per protocol, execution stops here for human review; CPTAC, PDC, protein validation, and new cohorts are not started automatically.

## Primary endpoint: DNA–RNA functional discordance

{chr(10).join(endpoint_lines)}

The combined denominator is the pooled count of histology-specific RNA-evaluable D1+D2 variant instances (`cluster_total >=10` in at least one prespecified candidate pair); it is descriptive and does not merge LUAD/LUSC contexts.

## Event summary

- Supported unique events: {len(events)}
- E1/E2/E3/E4/E5: {event_counts.get('E1',0)}/{event_counts.get('E2',0)}/{event_counts.get('E3',0)}/{event_counts.get('E4',0)}/{event_counts.get('E5',0)}
- R1/R2/R3/R4/R5: {r_counts.get('R1',0)}/{r_counts.get('R2',0)}/{r_counts.get('R3',0)}/{r_counts.get('R4',0)}/{r_counts.get('R5',0)}
- High-confidence discordant events: {len(high_events)} across {high_genes} genes and {len(patient_counter)} patients
- Non-obvious canonical splice-site consequences: {noncanonical}
- R3 Tier A rescue: {len(r3a)}/{d1_total} RNA-evaluable D1 variant instances = {100*r3_rate:.3f}% (Wilson 95% CI {100*r3_low:.3f}–{100*r3_high:.3f}%)

## Context and sensitivity

- RegTools exact known/recovered among high-confidence events: {int(high_events['regtools_known_recovered'].map(parse_bool).sum())}/{len(high_events)}. These are not novel discoveries.
- Legacy v26 functional conflicts: {int(events['automatic_downgrade_to_E4'].sum())}; automatic E4 downgrades: {int(events['automatic_downgrade_to_E4'].sum())}.
- GTEx zero/rarity is used only as normal-lung context and is not proof of tumor specificity.

## Gate checklist

{chr(10).join(checklist)}

## Scientific boundary

The event associations are statistical DNA–RNA links and are not molecular phasing proof. Final biological claims require manual review and orthogonal validation.
"""
    decision_path.write_text(decision_text, encoding="utf-8")

    output_paths = [event_path, r2_path, r3_path, qc_path, endpoint_path, r2_endpoint_path, decision_path]
    audit_path = args.outdir / "A_PHASE0_FINALIZATION_AUDIT.json"
    report = {
        "status": "PASS" if verdict in {"GO_MAIN_STRONG", "GO_MAIN_CONDITIONAL", "STOP_MAIN"} and qc_system_pass else "FAIL",
        "verdict": verdict,
        "gate": gate,
        "primary_endpoint": endpoint.to_dict("records"),
        "rna_evaluable_d1_d2_variant_histology_instances": pooled_denominator,
        "rna_evaluable_unique_genomic_variants_across_histologies": int(evaluable["variant_uid"].nunique()),
        "event_class_counts": {str(k): int(v) for k, v in event_counts.items()},
        "R_class_counts": {str(k): int(v) for k, v in r_counts.items()},
        "high_confidence_events": int(len(high_events)),
        "high_confidence_genes": high_genes,
        "high_confidence_patients": int(len(patient_counter)),
        "R3_tier_A": int(len(r3a)),
        "input_sha256": {path.name: sha256(path) for path in [args.evidence, args.reconstruction, args.recurrent, args.gtex, args.regtools, args.legacy]},
        "output_sha256": {path.name: sha256(path) for path in output_paths},
        "stop_after_phase0_for_human_review": True,
        "automatic_CPTAC_PDC_protein_new_cohorts_forbidden": True,
    }
    audit_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
