#!/usr/bin/env python3
"""Independent v2.3 functional-potential reconstruction for all eligible pairs.

Candidate eligibility is defined from the complete D2 candidate-link and
carrier/coverage universe.  Functional reconstruction then uses only variant
and junction geometry plus GENCODE50/MANE/GRCh38 assets; RNA target counts are
not used to select the functional class.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pysam


V23_ROOT = Path(__file__).resolve().parents[1]
V22_ROOT = Path(os.environ.get("V22_ROOT", str(V23_ROOT.parent / "v2.2")))
V22_SCRIPT = V22_ROOT / "scripts_v22" / "20_reconstruct_pilot_events.py"
PRIMARY_LINK_RULES = {
    "AFFECTED_CANONICAL_OR_NORMAL_FLANKING",
    "EXON_SKIP",
    "CRYPTIC_OR_ALTERNATIVE_DONOR_WITHIN_100",
    "CRYPTIC_OR_ALTERNATIVE_ACCEPTOR_WITHIN_100",
}
VALID_GEOMETRY = {"EXON_OVERLAP", "NEAR_EXON_BOUNDARY"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_v22_module():
    spec = importlib.util.spec_from_file_location("v22_reconstruction", V22_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import v2.2 reconstruction module: {V22_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_candidate_pairs(links: pd.DataFrame, carriers: pd.DataFrame, totals: pd.DataFrame) -> pd.DataFrame:
    carriers = carriers.copy()
    carriers["histology"] = carriers["histology"].astype(str).str.replace("TCGA-", "", regex=False)
    carriers["patient"] = carriers["patient"].astype(str)
    totals = totals.copy()
    totals["histology"] = totals["histology"].astype(str)
    totals["patient"] = totals["patient"].astype(str)
    d2 = links[
        links["DNA_class"].eq("D2")
        & links["link_rule"].isin(PRIMARY_LINK_RULES)
        & links["variant_geometry_status"].isin(VALID_GEOMETRY)
    ].copy()
    carrier_counts = (
        carriers.groupby(["histology", "variant_uid", "gene"], as_index=False)["patient"]
        .nunique()
        .rename(columns={"patient": "n_carriers_total"})
    )
    candidate = d2.merge(carrier_counts, on=["variant_uid", "gene"], how="inner")
    candidate = candidate[candidate["n_carriers_total"].isin([1, 2])].copy()
    eval_totals = totals[totals["cluster_total"] >= 10].copy()
    eval_counts = (
        eval_totals.groupby(["histology", "gene", "cluster_id"], as_index=False)["patient"]
        .nunique()
        .rename(columns={"patient": "n_rna_evaluable_patients"})
    )
    candidate = candidate.merge(eval_counts, on=["histology", "gene", "cluster_id"], how="left")
    candidate["n_rna_evaluable_patients"] = candidate["n_rna_evaluable_patients"].fillna(0).astype(int)
    carrier_key = carriers[["histology", "patient", "variant_uid", "gene"]].drop_duplicates()
    carrier_key = carrier_key.merge(
        candidate[["histology", "variant_uid", "gene", "cluster_id"]].drop_duplicates(),
        on=["histology", "variant_uid", "gene"], how="inner",
    )
    carrier_key = carrier_key.merge(
        eval_totals[["histology", "gene", "cluster_id", "patient"]].drop_duplicates(),
        on=["histology", "gene", "cluster_id", "patient"], how="inner",
    )
    carrier_in_cluster = (
        carrier_key.groupby(["histology", "variant_uid", "gene", "cluster_id"], as_index=False)["patient"]
        .nunique()
        .rename(columns={"patient": "n_carriers_evaluable"})
    )
    candidate = candidate.merge(carrier_in_cluster, on=["histology", "variant_uid", "gene", "cluster_id"], how="left")
    candidate["n_carriers_evaluable"] = candidate["n_carriers_evaluable"].fillna(0).astype(int)
    candidate["n_noncarriers_evaluable"] = candidate["n_rna_evaluable_patients"] - candidate["n_carriers_evaluable"]
    candidate = candidate[candidate["n_noncarriers_evaluable"] >= 20].copy()
    candidate = candidate.drop_duplicates(["histology", "variant_uid", "gene", "junction_id", "cluster_id"])
    candidate["candidate_histologies"] = candidate.groupby(["variant_uid", "gene", "junction_id"])["histology"].transform(lambda values: ",".join(sorted(set(values))))
    # Keep the exact histology-specific M1 universe from the supplied inputs.  The
    # geometry/function is identical across histologies, but retaining this
    # key prevents a later merge from silently changing the null universe.
    return candidate.drop_duplicates(["histology", "variant_uid", "gene", "junction_id", "cluster_id"]).copy()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--links", required=True, type=Path)
    ap.add_argument("--carriers", required=True, type=Path)
    ap.add_argument("--cluster-totals", required=True, type=Path)
    ap.add_argument("--junction-annotation", required=True, type=Path)
    ap.add_argument("--primary-consequence", required=True, type=Path)
    ap.add_argument("--transcripts", required=True, type=Path)
    ap.add_argument("--exons", required=True, type=Path)
    ap.add_argument("--cds", required=True, type=Path)
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--max-pairs", type=int, default=0, help="Smoke-test cap; 0 means all eligible histology-specific pairs")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    v22 = load_v22_module()

    links = pd.read_csv(args.links, sep="\t", compression="gzip", dtype=str, low_memory=False).fillna("")
    carriers = pd.read_csv(args.carriers, sep="\t", compression="gzip", low_memory=False).fillna("")
    totals = pd.read_csv(args.cluster_totals, sep="\t", compression="gzip", low_memory=False)
    candidate = build_candidate_pairs(links, carriers, totals)
    if args.max_pairs > 0:
        candidate = candidate.head(args.max_pairs).copy()
    genes = sorted(candidate["gene"].unique())

    junction = pd.read_csv(args.junction_annotation, sep="\t", compression="gzip", low_memory=False)
    junction = junction[junction["gene"].isin(genes) & junction["junction_id"].isin(set(candidate["junction_id"]))].copy()
    duplicate = junction.duplicated(["gene", "junction_id"], keep=False)
    if duplicate.any():
        raise RuntimeError("ambiguous junction annotation for candidate functional reconstruction")
    annotation_lookup = junction.set_index(["gene", "junction_id"]).to_dict("index")
    missing_annotation = sorted(set(zip(candidate["gene"], candidate["junction_id"])) - set(annotation_lookup))
    if missing_annotation:
        raise RuntimeError(f"candidate junction annotations missing: {missing_annotation[:10]}")

    tx_columns = ["chrom", "strand", "gene_name", "transcript_id", "transcript_version", "transcript_type", "tag"]
    transcripts = pq.read_table(args.transcripts, columns=tx_columns, filters=[("gene_name", "in", genes)]).to_pandas()
    exons = pq.read_table(args.exons, filters=[("gene_name", "in", genes)]).to_pandas()
    cds = pq.read_table(args.cds, filters=[("gene_name", "in", genes)]).to_pandas()
    fasta = pysam.FastaFile(str(args.fasta))
    models = {}
    for transcript_id, group in transcripts.groupby("transcript_id", sort=False):
        metadata = group.iloc[0].to_dict()
        tx_exons = exons[exons["transcript_id"] == transcript_id]
        tx_cds = cds[cds["transcript_id"] == transcript_id]
        if not tx_exons.empty:
            models[str(transcript_id)] = v22.TranscriptModel(metadata, tx_exons, tx_cds, fasta)

    primary = pd.read_csv(args.primary_consequence, sep="\t", low_memory=False)
    primary_cols = ["variant_uid", "gene", "HGVSc", "HGVSp", "Protein_position", "Amino_acids", "Codons", "EXON", "INTRON", "VARIANT_CLASS", "MANE_SELECT", "MANE_PLUS_CLINICAL"]
    primary = primary[primary_cols].drop_duplicates(["variant_uid", "gene"])
    candidates = candidate.merge(primary, on=["variant_uid", "gene"], how="left", validate="many_to_one")

    event_rows = []
    transcript_rows = []
    unresolved_selected = []
    for seed in candidates.to_dict("records"):
        annotation = annotation_lookup[(seed["gene"], seed["junction_id"])]
        selected_id = v22.strip_version(seed["selected_transcript"])
        selected_model = models.get(selected_id)
        gene_models = [model for model in models.values() if model.gene == seed["gene"]]
        mane_models = [model for model in gene_models if model.is_mane and model.cds]
        supporting_ids = [v22.strip_version(value) for value in str(annotation.get("current_transcripts_supporting", "") or "").split("|") if value]
        supporting_models = [models[tx] for tx in supporting_ids if tx in models]
        supporting_coding = [model for model in supporting_models if model.cds]
        evaluation_ids = []
        for model in ([selected_model] if selected_model is not None else []) + sorted(mane_models, key=lambda item: (not item.is_mane_select, item.transcript_id)) + supporting_coding:
            if model.transcript_id not in evaluation_ids:
                evaluation_ids.append(model.transcript_id)
        results_by_tx = {}
        for tx_id in evaluation_ids:
            model = models[tx_id]
            result = model.apply_junction(int(annotation["intron_start"]), int(annotation["intron_end"]))
            result.update({
                "variant_uid": seed["variant_uid"], "gene": seed["gene"], "junction_id": seed["junction_id"],
                "transcript_id": model.transcript_id, "transcript_full_id": model.full_id,
                "transcript_type": model.transcript_type, "is_mane": model.is_mane,
                "is_mane_select": model.is_mane_select, "is_selected_transcript": model.transcript_id == selected_id,
                "is_exact_supporting_transcript": model.transcript_id in supporting_ids,
            })
            results_by_tx[tx_id] = result
            transcript_rows.append(result)
        selected_result = results_by_tx.get(selected_id)
        if selected_result is None:
            unresolved_selected.append(seed["selected_transcript"])
        mane_results = [results_by_tx[model.transcript_id] for model in mane_models if model.transcript_id in results_by_tx]
        event_class, event_reason, functional_consensus = v22.classify_event(seed["junction_class"], selected_result, mane_results, len(supporting_coding))
        primary_result = selected_result or {"status": "UNRESOLVED_SELECTED_TRANSCRIPT_MISSING"}
        variant_retention = "VARIANT_RETENTION_UNRESOLVED"
        if primary_result.get("altered_exon_chain"):
            _chrom, variant_pos, variant_ref, _variant_alt = v22.parse_variant_uid(seed["variant_uid"])
            variant_positions = set(range(variant_pos, variant_pos + len(variant_ref)))
            altered_intervals = []
            for token in str(primary_result["altered_exon_chain"]).split("|"):
                match = re.match(r"^chr[^:]+:(\d+)-(\d+):[+-]$", token)
                if match:
                    altered_intervals.append((int(match.group(1)), int(match.group(2))))
            retained = any(any(start <= position <= end for start, end in altered_intervals) for position in variant_positions)
            variant_retention = "VARIANT_RETAINED_IN_ALTERED_TRANSCRIPT" if retained else "VARIANT_EXCLUDED_FROM_ALTERED_TRANSCRIPT"
        row = {
            "histology": seed["histology"], "variant_uid": seed["variant_uid"], "gene": seed["gene"], "junction_id": seed["junction_id"],
            "cluster_id": seed["cluster_id"], "candidate_histologies": seed["candidate_histologies"],
            "n_carriers_total": seed["n_carriers_total"], "n_noncarriers_evaluable": seed["n_noncarriers_evaluable"],
            "DNA_class": seed["DNA_class"], "primary_consequence": seed["primary_consequence"],
            "selected_transcript": seed["selected_transcript"], "junction_class": seed["junction_class"],
            "link_rule": seed["link_rule"], "variant_geometry_status": seed["variant_geometry_status"],
            "event_class": event_class, "event_class_reason": event_reason,
            "functional_consensus": bool(functional_consensus),
            "altered_splice_type": primary_result.get("altered_splice_type", "UNRESOLVED"),
            "altered_exon_chain": primary_result.get("altered_exon_chain", ""),
            "frame_consequence": primary_result.get("frame_consequence", "FRAME_UNRESOLVED"),
            "nmd_status": primary_result.get("nmd_status", "NMD_UNRESOLVED"),
            "reference_protein_length_aa": primary_result.get("reference_protein_length_aa", np.nan),
            "altered_protein_length_aa": primary_result.get("altered_protein_length_aa", np.nan),
            "variant_retention": variant_retention,
            "v23_functional_potential_status": "PASS_INDEPENDENT_GEOMETRY_RECONSTRUCTION",
            "R2_FUNCTIONALLY_ELIGIBLE": bool(event_class.startswith(("E1", "E2", "E3")) and functional_consensus),
            "mane_transcripts_evaluated": len(mane_results),
            "current_exact_coding_transcripts": len(supporting_coding),
            "v23_domain_status": "UNRESOLVED_EXACT_DOMAIN_CROSSREF_NOT_RUN",
            "cis_status": "STATISTICAL_LINK_NOT_MOLECULARLY_PHASED",
        }
        row.update({key: seed.get(key, "") for key in ["HGVSc", "HGVSp", "Protein_position", "Amino_acids", "Codons", "EXON", "INTRON", "VARIANT_CLASS", "MANE_SELECT", "MANE_PLUS_CLINICAL"]})
        event_rows.append(row)

    events = pd.DataFrame(event_rows)
    transcript = pd.DataFrame(transcript_rows)
    compression = {"method": "gzip", "compresslevel": 6, "mtime": 0}
    event_path = args.outdir / "V23_FULL_D2_FUNCTIONAL_POTENTIAL.tsv.gz"
    transcript_path = args.outdir / "V23_FULL_D2_TRANSCRIPT_RECONSTRUCTION.tsv.gz"
    audit_path = args.outdir / "V23_FULL_D2_FUNCTIONAL_POTENTIAL_AUDIT.json"
    events.to_csv(event_path, sep="\t", index=False, compression=compression)
    transcript.to_csv(transcript_path, sep="\t", index=False, compression=compression)
    audit = {
        "status": "PASS" if len(events) == len(candidates) and not unresolved_selected else "PASS_WITH_UNRESOLVED_SELECTED_TRANSCRIPTS",
        "candidate_pairs": int(len(candidates)),
        "functional_rows": int(len(events)),
        "transcript_rows": int(len(transcript)),
        "event_class_counts": {str(k): int(v) for k, v in events["event_class"].value_counts().to_dict().items()},
        "r2_functionally_eligible": int(events["R2_FUNCTIONALLY_ELIGIBLE"].sum()),
        "unresolved_selected_transcript_count": len(unresolved_selected),
        "unresolved_selected_transcripts": sorted(set(unresolved_selected))[:50],
        "functional_independence": "variant_geometry_plus_GENCODE50_MANE_GRCh38_only; no RNA target_count used",
        "source_sha256": {key: sha256(path) for key, path in {"links": args.links, "carriers": args.carriers, "cluster_totals": args.cluster_totals, "junction_annotation": args.junction_annotation, "primary_consequence": args.primary_consequence, "transcripts": args.transcripts, "exons": args.exons, "cds": args.cds, "fasta": args.fasta}.items()},
        "output_sha256": {event_path.name: sha256(event_path), transcript_path.name: sha256(transcript_path)},
    }
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
