#!/usr/bin/env python3
"""Create deterministic TCGA LUAD/LUSC primary-tumor RNA patient maps."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


META_COLUMNS = [
    "rail_id",
    "study",
    "tcga.tcga_barcode",
    "tcga.gdc_cases.project.project_id",
    "tcga.cgc_sample_sample_type",
    "tcga.gdc_cases.samples.sample_type",
    "recount_qc.star.number_of_input_reads",
    "recount_qc.star.uniquely_mapped_reads_%",
    "recount_qc.star.uniquely_mapped_reads_number",
    "recount_qc.junction_coverage",
]
PATIENT_RE = re.compile(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", re.I)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_local_depth(path: Path, panel: set[str]) -> pd.Series:
    pieces = []
    for chunk in pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        usecols=["query_gene", "rail_id", "count"],
        dtype={"query_gene": "string", "rail_id": "string"},
        chunksize=1_000_000,
    ):
        chunk = chunk[chunk["query_gene"].isin(panel)]
        chunk["count"] = pd.to_numeric(chunk["count"], errors="coerce").fillna(0)
        pieces.append(chunk.groupby("rail_id", sort=False)["count"].sum())
    if not pieces:
        return pd.Series(dtype=float)
    merged = pd.concat(pieces)
    return merged.groupby(level=0).sum()


def primary_tumor_mask(frame: pd.DataFrame) -> pd.Series:
    type_a = frame["tcga.cgc_sample_sample_type"].fillna("").str.casefold()
    type_b = frame["tcga.gdc_cases.samples.sample_type"].fillna("").str.casefold()
    barcode = frame["rna_barcode"].fillna("")
    sample_code_01 = barcode.str.slice(13, 15).eq("01")
    return (type_a.eq("primary tumor") | type_b.eq("primary tumor")) & sample_code_01


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True, type=Path)
    ap.add_argument("--gdc-map", required=True, type=Path)
    ap.add_argument("--panel-crosswalk", required=True, type=Path)
    ap.add_argument("--luad-counts", required=True, type=Path)
    ap.add_argument("--lusc-counts", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--min-input-reads", type=int, default=10_000_000)
    ap.add_argument("--min-unique-mapped-pct", type=float, default=50.0)
    ap.add_argument("--min-junction-coverage", type=int, default=1_000_000)
    args = ap.parse_args()

    panel_df = pd.read_csv(args.panel_crosswalk, sep="\t", dtype=str).fillna("")
    panel = set(panel_df["gene_symbol"])
    if len(panel) != 74:
        raise ValueError(f"expected 74 frozen panel genes, found {len(panel)}")

    metadata = pd.read_csv(
        args.metadata,
        sep="\t",
        usecols=META_COLUMNS,
        dtype=str,
        low_memory=False,
    ).fillna("")
    metadata = metadata[metadata["study"].isin(["LUAD", "LUSC"])].copy()
    metadata["rail_id"] = metadata["rail_id"].astype(str)
    metadata["rna_barcode"] = metadata["tcga.tcga_barcode"].astype(str).str.upper()
    metadata["patient"] = metadata["rna_barcode"].str.extract(f"({PATIENT_RE.pattern})", expand=False).fillna("")
    metadata = metadata[metadata["patient"].ne("") & primary_tumor_mask(metadata)].copy()

    numeric_map = {
        "input_reads": "recount_qc.star.number_of_input_reads",
        "unique_mapped_pct": "recount_qc.star.uniquely_mapped_reads_%",
        "unique_mapped_reads": "recount_qc.star.uniquely_mapped_reads_number",
        "global_junction_coverage": "recount_qc.junction_coverage",
    }
    for new, old in numeric_map.items():
        metadata[new] = pd.to_numeric(metadata[old], errors="coerce")

    metadata["qc_pass"] = (
        metadata["input_reads"].ge(args.min_input_reads)
        & metadata["unique_mapped_pct"].ge(args.min_unique_mapped_pct)
        & metadata["global_junction_coverage"].ge(args.min_junction_coverage)
    )
    metadata["qc_failure_reason"] = ""
    metadata.loc[metadata["input_reads"].lt(args.min_input_reads) | metadata["input_reads"].isna(), "qc_failure_reason"] += "LOW_OR_MISSING_INPUT_READS;"
    metadata.loc[metadata["unique_mapped_pct"].lt(args.min_unique_mapped_pct) | metadata["unique_mapped_pct"].isna(), "qc_failure_reason"] += "LOW_OR_MISSING_UNIQUE_MAP_PCT;"
    metadata.loc[metadata["global_junction_coverage"].lt(args.min_junction_coverage) | metadata["global_junction_coverage"].isna(), "qc_failure_reason"] += "LOW_OR_MISSING_JUNCTION_COVERAGE;"
    metadata["qc_failure_reason"] = metadata["qc_failure_reason"].str.rstrip(";")

    gdc = pd.read_csv(args.gdc_map, sep="\t", dtype=str).fillna("")
    dna_cases = {
        histology: set(gdc.loc[gdc["project_id"].eq(f"TCGA-{histology}"), "case_submitter_id"])
        for histology in ("LUAD", "LUSC")
    }

    depths = {
        "LUAD": load_local_depth(args.luad_counts, panel),
        "LUSC": load_local_depth(args.lusc_counts, panel),
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    audit = {
        "status": "PASS",
        "parameters": {
            "primary_tumor_required": True,
            "min_input_reads": args.min_input_reads,
            "min_unique_mapped_pct": args.min_unique_mapped_pct,
            "min_junction_coverage": args.min_junction_coverage,
            "selection_order": ["QC_PASS", "PRIMARY74_LOCAL_JUNCTION_DEPTH_DESC", "BARCODE_ASC", "RAIL_ID_ASC"],
        },
        "histologies": {},
        "source_sha256": {
            "metadata": sha256(args.metadata),
            "gdc_map": sha256(args.gdc_map),
            "panel_crosswalk": sha256(args.panel_crosswalk),
            "luad_counts": sha256(args.luad_counts),
            "lusc_counts": sha256(args.lusc_counts),
        },
    }

    for histology in ("LUAD", "LUSC"):
        candidates = metadata[metadata["study"].eq(histology)].copy()
        candidates["local_primary74_junction_depth"] = candidates["rail_id"].map(depths[histology]).fillna(0)
        candidates["dna_case_present"] = candidates["patient"].isin(dna_cases[histology])
        candidate_counts = candidates.groupby("patient")["rail_id"].size()
        qc_counts = candidates[candidates["qc_pass"]].groupby("patient")["rail_id"].size()
        eligible = candidates[candidates["qc_pass"] & candidates["dna_case_present"]].copy()
        eligible = eligible.sort_values(
            ["patient", "local_primary74_junction_depth", "rna_barcode", "rail_id"],
            ascending=[True, False, True, True],
            kind="mergesort",
        )
        selected = eligible.drop_duplicates("patient", keep="first").copy()
        selected["primary_tumor_rna_candidate_count"] = selected["patient"].map(candidate_counts).astype(int)
        selected["qc_pass_candidate_count"] = selected["patient"].map(qc_counts).fillna(0).astype(int)
        selected["selection_reason"] = "QC_PASS;MAX_PRIMARY74_LOCAL_JUNCTION_DEPTH;BARCODE_LEXICOGRAPHIC_TIEBREAK"
        output_cols = [
            "study", "patient", "rail_id", "rna_barcode", "tcga.gdc_cases.project.project_id",
            "tcga.cgc_sample_sample_type", "tcga.gdc_cases.samples.sample_type", "qc_pass",
            "qc_failure_reason", "input_reads", "unique_mapped_pct", "unique_mapped_reads",
            "global_junction_coverage", "local_primary74_junction_depth", "dna_case_present",
            "primary_tumor_rna_candidate_count", "qc_pass_candidate_count", "selection_reason",
        ]
        selected = selected[output_cols].rename(columns={"study": "histology"})
        out_path = args.outdir / f"TCGA_{histology}_PATIENT_MAP_V22.tsv"
        selected.to_csv(out_path, sep="\t", index=False)
        audit["histologies"][histology] = {
            "primary_tumor_rna_candidates": int(len(candidates)),
            "candidate_patients": int(candidates["patient"].nunique()),
            "qc_pass_candidates": int(candidates["qc_pass"].sum()),
            "dna_matched_qc_pass_candidates": int(len(eligible)),
            "selected_unique_patients": int(len(selected)),
            "selected_unique_rail_ids": int(selected["rail_id"].nunique()),
            "selected_duplicate_patients": int(selected["patient"].duplicated().sum()),
            "selected_duplicate_rail_ids": int(selected["rail_id"].duplicated().sum()),
            "output": str(out_path),
            "output_sha256": sha256(out_path),
        }
        if selected.empty or selected["patient"].duplicated().any() or selected["rail_id"].duplicated().any():
            audit["status"] = "FAIL"

    audit_path = args.outdir / "TCGA_DNA_RNA_MAPPING_AUDIT_V22.json"
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0 if audit["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
