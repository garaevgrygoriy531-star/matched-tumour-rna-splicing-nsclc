#!/usr/bin/env python3
"""Build the complete RNA-evaluable rare-D2 candidate/sample universe.

This is an M1 preparation step.  It deliberately does not infer functional
consequence from RNA-positive events.  Without the GENCODE50/MANE transcript
and sequence context, functional-potential columns remain explicit HOLD
values and M2 must not be run.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
from pathlib import Path

import pandas as pd


V23_ROOT = Path(__file__).resolve().parents[1]
V22_ROOT = Path(os.environ.get("V22_ROOT", str(V23_ROOT.parent / "v2.2")))
TABLES = V23_ROOT / "results" / "tables"
AUDIT = V23_ROOT / "results" / "audit"

PRIMARY_LINK_RULES = {
    "AFFECTED_CANONICAL_OR_NORMAL_FLANKING",
    "EXON_SKIP",
    "CRYPTIC_OR_ALTERNATIVE_DONOR_WITHIN_100",
    "CRYPTIC_OR_ALTERNATIVE_ACCEPTOR_WITHIN_100",
}
VALID_GEOMETRY = {"EXON_OVERLAP", "NEAR_EXON_BOUNDARY"}


def main() -> int:
    link_path = V22_ROOT / "06_variant_junction_link" / "v2.2_full74" / "FULL74_VARIANT_JUNCTION_CANDIDATE_LINKS.tsv.gz"
    carrier_path = V22_ROOT / "06_variant_junction_link" / "v2.2_full74" / "FULL74_VARIANT_CARRIERS.tsv.gz"
    totals_path = V22_ROOT / "05_junction_universe" / "v2.2_full74" / "FULL74_CLUSTER_TOTALS.tsv.gz"
    sparse_path = V22_ROOT / "05_junction_universe" / "v2.2_full74" / "FULL74_JUNCTION_COUNTS_SPARSE.tsv.gz"
    for path in (link_path, carrier_path, totals_path, sparse_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    links = pd.read_csv(link_path, sep="\t", compression="gzip", dtype=str, low_memory=False).fillna("")
    carriers = pd.read_csv(carrier_path, sep="\t", compression="gzip", dtype=str, low_memory=False).fillna("")
    totals = pd.read_csv(totals_path, sep="\t", compression="gzip", dtype={"cluster_total": "float64"}, low_memory=False)
    sparse = pd.read_csv(sparse_path, sep="\t", compression="gzip", dtype={"count": "float64"}, low_memory=False)

    carriers["histology"] = carriers["histology"].astype(str).str.replace("TCGA-", "", regex=False)
    carriers["patient"] = carriers["patient"].astype(str)
    totals["histology"] = totals["histology"].astype(str)
    totals["patient"] = totals["patient"].astype(str)
    sparse["histology"] = sparse["histology"].astype(str)
    sparse["patient"] = sparse["patient"].astype(str)

    link_filter = links["DNA_class"].eq("D2") & links["link_rule"].isin(PRIMARY_LINK_RULES) & links["variant_geometry_status"].isin(VALID_GEOMETRY)
    d2_links = links.loc[link_filter].copy()

    carrier_counts = (
        carriers.groupby(["histology", "variant_uid", "gene"], as_index=False)["patient"]
        .nunique()
        .rename(columns={"patient": "n_carriers_total"})
    )
    # A pair is tested separately in each histology; do not manufacture the
    # other histology when the variant has no carrier there.
    candidate_pairs = d2_links.merge(carrier_counts, on=["variant_uid", "gene"], how="inner")
    candidate_pairs = candidate_pairs[candidate_pairs["n_carriers_total"].isin([1, 2])].copy()

    eval_totals = totals[totals["cluster_total"] >= 10].copy()
    eval_counts = (
        eval_totals.groupby(["histology", "gene", "cluster_id"], as_index=False)["patient"]
        .nunique()
        .rename(columns={"patient": "n_rna_evaluable_patients"})
    )
    candidate_pairs = candidate_pairs.merge(eval_counts, on=["histology", "gene", "cluster_id"], how="left")
    candidate_pairs["n_rna_evaluable_patients"] = candidate_pairs["n_rna_evaluable_patients"].fillna(0).astype(int)

    # The original v2.2 rare implementation evaluates non-carriers among
    # cluster_total >= 10.  Subtract observed carrier patients in that same
    # cluster; carriers absent from the cluster remain non-carrier in the
    # observed RNA matrix, exactly as in v2.2.
    carrier_key = carriers[["histology", "patient", "variant_uid", "gene"]].drop_duplicates()
    carrier_key["carrier_status"] = "TRUE_CARRIER"
    carrier_key = carrier_key.merge(
        candidate_pairs[["histology", "variant_uid", "gene", "cluster_id"]].drop_duplicates(),
        on=["histology", "variant_uid", "gene"],
        how="inner",
    )
    carrier_key = carrier_key.merge(
        eval_totals[["histology", "gene", "cluster_id", "patient"]].drop_duplicates(),
        on=["histology", "gene", "cluster_id", "patient"],
        how="inner",
    )
    carrier_in_cluster = (
        carrier_key.groupby(["histology", "variant_uid", "gene", "cluster_id"], as_index=False)["patient"]
        .nunique()
        .rename(columns={"patient": "n_carriers_evaluable"})
    )
    candidate_pairs = candidate_pairs.merge(
        carrier_in_cluster,
        on=["histology", "variant_uid", "gene", "cluster_id"],
        how="left",
    )
    candidate_pairs["n_carriers_evaluable"] = candidate_pairs["n_carriers_evaluable"].fillna(0).astype(int)
    candidate_pairs["n_noncarriers_evaluable"] = candidate_pairs["n_rna_evaluable_patients"] - candidate_pairs["n_carriers_evaluable"]
    candidate_pairs = candidate_pairs[candidate_pairs["n_noncarriers_evaluable"] >= 20].copy()
    candidate_pairs = candidate_pairs.drop_duplicates(["histology", "variant_uid", "gene", "junction_id", "cluster_id"])

    # Sparse input is already nonzero-only.  Aggregate defensively in case a
    # future producer emits duplicate rows for the same junction/patient.
    target = (
        sparse.groupby(["histology", "junction_id", "patient"], as_index=False)["count"]
        .sum()
        .rename(columns={"count": "target_count"})
    )
    pair_columns = [
        "histology", "variant_uid", "gene", "junction_id", "cluster_id", "DNA_class",
        "primary_consequence", "selected_transcript", "variant_geometry_status", "link_rule",
        "n_carriers_total", "n_carriers_evaluable", "n_rna_evaluable_patients", "n_noncarriers_evaluable",
    ]
    candidate_pairs = candidate_pairs[pair_columns].copy()

    output = TABLES / "V23_RARE_D2_NULL_UNIVERSE.tsv.gz"
    output.parent.mkdir(parents=True, exist_ok=True)
    first = True
    sample_rows = 0
    pair_count = len(candidate_pairs)
    first = True
    output.unlink(missing_ok=True)
    with gzip.open(output, "wt", encoding="utf-8", newline="") as handle:
        for start in range(0, pair_count, 500):
            block = candidate_pairs.iloc[start : start + 500].copy()
            block = block.merge(
                eval_totals[["histology", "gene", "cluster_id", "patient", "cluster_total"]],
                on=["histology", "gene", "cluster_id"], how="inner"
            )
            block = block.merge(target, on=["histology", "junction_id", "patient"], how="left")
            block["target_count"] = block["target_count"].fillna(0.0)
            block["cluster_total"] = block["cluster_total"].astype(float)
            block["junction_usage"] = block["target_count"] / block["cluster_total"]
            block = block.merge(
                carrier_key[["histology", "variant_uid", "gene", "cluster_id", "patient", "carrier_status"]].drop_duplicates(),
                on=["histology", "variant_uid", "gene", "cluster_id", "patient"], how="left"
            )
            block["carrier_status"] = block["carrier_status"].fillna("NONCARRIER")
            block["coverage_evaluable"] = True
            block["functional_potential_status"] = "HOLD_MISSING_GENCODE50_MANE_SEQUENCE_CONTEXT"
            for col in ["event_class", "altered_splice_type", "frame_consequence", "nmd_status", "functional_consensus"]:
                block[col] = ""
            block["R2_FUNCTIONALLY_ELIGIBLE"] = "HOLD"
            columns = [
                "histology", "variant_uid", "gene", "junction_id", "cluster_id", "DNA_class",
                "primary_consequence", "selected_transcript", "variant_geometry_status", "link_rule",
                "patient", "carrier_status", "target_count", "cluster_total", "junction_usage",
                "n_carriers_total", "n_carriers_evaluable", "n_rna_evaluable_patients", "n_noncarriers_evaluable",
                "coverage_evaluable", "functional_potential_status", "event_class", "altered_splice_type",
                "frame_consequence", "nmd_status", "functional_consensus", "R2_FUNCTIONALLY_ELIGIBLE",
            ]
            block[columns].to_csv(handle, sep="\t", index=False, header=first, float_format="%.10g")
            first = False
            sample_rows += len(block)

    status = {
        "m1_status": "HOLD_MISSING_GENCODE50_MANE_SEQUENCE_CONTEXT",
        "raw_rare_d2_candidate_pairs": int(pair_count),
        "sample_level_rows": int(sample_rows),
        "coverage_rule": "cluster_total >= 10; non-carrier evaluable count >= 20",
        "candidate_source_rows_before_filters": int(len(links)),
        "d2_primary_geometry_link_rows": int(len(d2_links)),
        "output": str(output),
        "functional_potential_is_independent": False,
        "m2_permutation_allowed": False,
        "hold_reason": "Local workspace lacks the GENCODE50/MANE transcript exon/CDS and GRCh38 sequence assets required to compute event_class, frame_consequence and nmd_status independently of carrier RNA evidence."
    }
    (AUDIT / "V23_M1_STATUS.json").write_text(json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(status, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
