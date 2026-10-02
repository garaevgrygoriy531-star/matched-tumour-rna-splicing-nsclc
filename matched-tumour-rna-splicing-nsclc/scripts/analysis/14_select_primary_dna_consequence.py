#!/usr/bin/env python3
"""Select deterministic primary DNA consequences without severity-driven transcript picking."""

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path


D1 = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant"}
D2 = {
    "missense_variant", "synonymous_variant", "splice_region_variant", "inframe_insertion",
    "inframe_deletion", "protein_altering_variant", "coding_sequence_variant",
    "splice_donor_5th_base_variant", "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
}


def value_present(value: str) -> bool:
    return value not in {"", "-", "."}


def tier(row: dict[str, str]) -> tuple[int, str]:
    if value_present(row.get("MANE_SELECT", "")):
        return 1, "MANE_SELECT"
    if value_present(row.get("MANE_PLUS_CLINICAL", "")):
        return 2, "MANE_PLUS_CLINICAL"
    if row.get("CANONICAL", "") == "YES":
        return 3, "VEP_CANONICAL"
    return 4, "NO_PREFERRED_TRANSCRIPT"


def d_class(consequence: str) -> str:
    terms = set(consequence.split(","))
    if terms & D1:
        return "D1"
    if terms & D2:
        return "D2"
    return "D3"


def read_vep(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        reader = None
        for line in handle:
            if line.startswith("##") or not line.strip():
                continue
            if line.startswith("#variant_uid"):
                header = line.rstrip("\n").lstrip("#").split("\t")
                reader = header
                continue
            if line.startswith("#"):
                continue
            if reader is None:
                raise RuntimeError("VEP header not found")
            values = line.rstrip("\n").split("\t")
            yield dict(zip(reader, values))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vep", required=True, type=Path)
    ap.add_argument("--mc3-panel-rows", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    intended: dict[str, set[str]] = defaultdict(set)
    with gzip.open(args.mc3_panel_rows, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            intended[row["variant_uid"]].add(row["gene"])
    by_pair: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    all_rows_by_uid: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_vep(args.vep):
        uid = row["variant_uid"]
        all_rows_by_uid[uid].append(row)
        symbol = row.get("SYMBOL", "")
        if symbol in intended.get(uid, set()):
            by_pair[(uid, symbol)].append(row)

    primary_path = args.outdir / "DNA_VARIANT_PRIMARY_CONSEQUENCE.tsv"
    conflicts_path = args.outdir / "DNA_VARIANT_TRANSCRIPT_CONFLICTS.tsv"
    audit_path = args.outdir / "DNA_VARIANT_PRIMARY_CONSEQUENCE_AUDIT.json"
    primary_fields = [
        "variant_uid", "gene", "selection_tier", "selection_rule", "selected_transcript",
        "selected_transcript_type", "primary_consequence", "all_consequences_on_selected_row",
        "DNA_class", "IMPACT", "MANE_SELECT", "MANE_PLUS_CLINICAL", "CANONICAL",
        "HGVSc", "HGVSp", "Protein_position", "Amino_acids", "Codons", "EXON", "INTRON",
        "VARIANT_CLASS", "selected_tier_row_count", "selected_tier_conflict", "status",
    ]
    conflict_fields = [
        "variant_uid", "gene", "selection_tier", "selection_rule", "Feature", "Consequence",
        "HGVSc", "HGVSp", "MANE_SELECT", "MANE_PLUS_CLINICAL", "CANONICAL", "conflict_reason",
    ]
    expected_pairs = {(uid, gene) for uid, genes in intended.items() for gene in genes}
    resolved = 0
    conflict_pairs = 0
    no_preferred_pairs = 0
    class_counts = defaultdict(int)

    with primary_path.open("w", encoding="utf-8", newline="") as primary_handle, \
         conflicts_path.open("w", encoding="utf-8", newline="") as conflict_handle:
        primary_writer = csv.DictWriter(primary_handle, fieldnames=primary_fields, delimiter="\t", lineterminator="\n")
        conflict_writer = csv.DictWriter(conflict_handle, fieldnames=conflict_fields, delimiter="\t", lineterminator="\n")
        primary_writer.writeheader()
        conflict_writer.writeheader()
        for uid, gene in sorted(expected_pairs):
            rows = by_pair.get((uid, gene), [])
            if not rows:
                # Preserve unresolved variants rather than selecting another gene/transcript by severity.
                primary_writer.writerow({
                    "variant_uid": uid, "gene": gene, "selection_tier": 9,
                    "selection_rule": "UNRESOLVED_PANEL_GENE", "DNA_class": "D3",
                    "status": "UNRESOLVED_NO_SYMBOL_MATCH",
                })
                class_counts["D3"] += 1
                continue
            best_tier = min(tier(row)[0] for row in rows)
            candidates = [row for row in rows if tier(row)[0] == best_tier]
            candidates.sort(key=lambda row: (row.get("Feature", ""), row.get("Consequence", ""), row.get("HGVSc", "")))
            chosen = candidates[0]
            selection_rule = tier(chosen)[1]
            signatures = {
                (row.get("Consequence", ""), row.get("HGVSc", ""), row.get("HGVSp", ""))
                for row in candidates
            }
            has_conflict = len(signatures) > 1
            if has_conflict:
                conflict_pairs += 1
                for row in candidates:
                    conflict_writer.writerow({
                        "variant_uid": uid, "gene": gene, "selection_tier": best_tier,
                        "selection_rule": selection_rule, "Feature": row.get("Feature", ""),
                        "Consequence": row.get("Consequence", ""), "HGVSc": row.get("HGVSc", ""),
                        "HGVSp": row.get("HGVSp", ""), "MANE_SELECT": row.get("MANE_SELECT", ""),
                        "MANE_PLUS_CLINICAL": row.get("MANE_PLUS_CLINICAL", ""),
                        "CANONICAL": row.get("CANONICAL", ""),
                        "conflict_reason": "MULTIPLE_PREFERRED_ROWS_WITH_DIFFERENT_CONSEQUENCES",
                    })
            consequence = chosen.get("Consequence", "")
            primary_term = consequence.split(",")[0] if consequence else ""
            if best_tier == 4:
                # A lexical transcript is recorded only for provenance. It is not
                # allowed to create a D1/D2 event without a preferred transcript.
                klass = "D3"
                no_preferred_pairs += 1
                row_status = "UNRESOLVED_NO_PREFERRED_TRANSCRIPT"
            else:
                klass = d_class(consequence)
                row_status = "PASS"
            class_counts[klass] += 1
            resolved += 1
            primary_writer.writerow({
                "variant_uid": uid, "gene": gene, "selection_tier": best_tier,
                "selection_rule": selection_rule, "selected_transcript": chosen.get("Feature", ""),
                "selected_transcript_type": chosen.get("Feature_type", ""),
                "primary_consequence": primary_term,
                "all_consequences_on_selected_row": consequence, "DNA_class": klass,
                "IMPACT": chosen.get("IMPACT", ""), "MANE_SELECT": chosen.get("MANE_SELECT", ""),
                "MANE_PLUS_CLINICAL": chosen.get("MANE_PLUS_CLINICAL", ""),
                "CANONICAL": chosen.get("CANONICAL", ""), "HGVSc": chosen.get("HGVSc", ""),
                "HGVSp": chosen.get("HGVSp", ""), "Protein_position": chosen.get("Protein_position", ""),
                "Amino_acids": chosen.get("Amino_acids", ""), "Codons": chosen.get("Codons", ""),
                "EXON": chosen.get("EXON", ""), "INTRON": chosen.get("INTRON", ""),
                "VARIANT_CLASS": chosen.get("VARIANT_CLASS", ""),
                "selected_tier_row_count": len(candidates),
                "selected_tier_conflict": str(has_conflict), "status": row_status,
            })

    audit = {
        "status": "PASS" if resolved == len(expected_pairs) else "HOLD",
        "expected_variant_gene_pairs": len(expected_pairs),
        "resolved_variant_gene_pairs": resolved,
        "unresolved_variant_gene_pairs": len(expected_pairs) - resolved,
        "preferred_tier_conflict_pairs": conflict_pairs,
        "no_preferred_transcript_pairs": no_preferred_pairs,
        "DNA_class_counts": dict(sorted(class_counts.items())),
        "selection_hierarchy": ["MANE_SELECT", "MANE_PLUS_CLINICAL", "VEP_CANONICAL", "NO_PREFERRED_TRANSCRIPT"],
        "deterministic_tie_break": "Feature transcript ID lexicographic; never consequence severity",
    }
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0 if audit["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
