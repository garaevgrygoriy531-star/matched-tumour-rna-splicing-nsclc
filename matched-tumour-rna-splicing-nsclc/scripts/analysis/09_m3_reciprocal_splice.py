#!/usr/bin/env python3
"""M3 reciprocal splice analysis for the frozen nine R2 events.

Reference competitors are generated from the same GENCODE50 selected-
transcript exon chains used by v2.2/v2.3 functional reconstruction.  TCGA
usage is then evaluated from the frozen v2.2 sparse junction infrastructure
within the candidate's own cluster and coverage universe.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pysam


MAX_CRYPTIC_EXTENSION = 100


def strip_version(value: str) -> str:
    return str(value or "").split(".", 1)[0]


def parse_junction(value: str) -> tuple[str, int, int, str]:
    match = re.match(r"^(chr[^:]+):(\d+)-(\d+):([+-])$", str(value))
    if not match:
        raise ValueError(f"invalid junction_id: {value}")
    return match.group(1), int(match.group(2)), int(match.group(3)), match.group(4)


def import_reconstruction(path: Path):
    spec = importlib.util.spec_from_file_location("v22_reconstruction", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import v2.2 reconstruction implementation: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bool_value(value) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def reference_junction_id(model, left_index: int) -> str:
    left = model.exons[left_index]
    right = model.exons[left_index + 1]
    intron_start = min(int(left["genomic_end"]), int(right["genomic_end"])) + 1
    intron_end = max(int(left["genomic_start"]), int(right["genomic_start"])) - 1
    return f"{model.chrom}:{intron_start}-{intron_end}:{model.strand}"


def reference_chain(model) -> str:
    return "|".join(reference_junction_id(model, index) for index in range(len(model.exons) - 1))


def choose_flanking_exons(model, intron_start: int, intron_end: int) -> tuple[int, int] | None:
    if model.strand == "+":
        upstream_base = int(intron_start) - 1
        downstream_base = int(intron_end) + 1
    else:
        upstream_base = -(int(intron_end) + 1)
        downstream_base = -(int(intron_start) - 1)
    upstream = []
    downstream = []
    for index, exon in enumerate(model.exons):
        if exon["oriented_start"] <= upstream_base <= exon["oriented_end"] + MAX_CRYPTIC_EXTENSION:
            extension = max(0, upstream_base - exon["oriented_end"])
            truncation = max(0, exon["oriented_end"] - upstream_base)
            upstream.append((extension + truncation, index))
        if exon["oriented_start"] - MAX_CRYPTIC_EXTENSION <= downstream_base <= exon["oriented_end"]:
            extension = max(0, exon["oriented_start"] - downstream_base)
            truncation = max(0, downstream_base - exon["oriented_start"])
            downstream.append((extension + truncation, index))
    pairs = [(up[0] + down[0], up[1], down[1]) for up in upstream for down in downstream if up[1] <= down[1]]
    if not pairs:
        return None
    minimum = min(item[0] for item in pairs)
    best = {(item[1], item[2]) for item in pairs if item[0] == minimum}
    if len(best) != 1:
        return None
    return next(iter(best))


def stats(values: np.ndarray) -> dict[str, float | int]:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return {"n": 0, "mean": np.nan, "median": np.nan, "n_nonzero": 0, "sum": 0.0}
    return {
        "n": int(values.size),
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "n_nonzero": int(np.sum(values > 0)),
        "sum": float(np.sum(values)),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--r2-freeze", required=True, type=Path)
    ap.add_argument("--functional", required=True, type=Path)
    ap.add_argument("--carriers", required=True, type=Path)
    ap.add_argument("--sparse-counts", required=True, type=Path)
    ap.add_argument("--cluster-totals", required=True, type=Path)
    ap.add_argument("--transcripts", required=True, type=Path)
    ap.add_argument("--exons", required=True, type=Path)
    ap.add_argument("--cds", required=True, type=Path)
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--reconstruction-implementation", required=True, type=Path)
    ap.add_argument("--reference-sha256", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    r2 = pd.read_csv(args.r2_freeze, sep="\t", dtype=str, low_memory=False).fillna("")
    if len(r2) != 9:
        raise ValueError(f"frozen R2 positive control must contain 9 rows, found {len(r2)}")
    functional = pd.read_csv(args.functional, sep="\t", compression="gzip", low_memory=False).fillna("")
    key = ["variant_uid", "gene", "junction_id"]
    func_unique = functional.sort_values(["gene", "variant_uid", "junction_id", "histology"]).drop_duplicates(key)
    events = r2.merge(
        func_unique[key + ["histology", "cluster_id", "selected_transcript", "event_class", "altered_splice_type", "altered_exon_chain", "functional_consensus"]],
        on=key, how="left", suffixes=("_v22", "_v23"), validate="one_to_one",
    )
    if events["selected_transcript_v23"].eq("").any():
        raise RuntimeError("one or more frozen R2 events is missing from the independent v2.3 functional table")

    genes = sorted(events["gene"].unique())
    transcripts = pq.read_table(
        args.transcripts,
        columns=["chrom", "strand", "gene_name", "transcript_id", "transcript_version", "transcript_type", "tag"],
        filters=[("gene_name", "in", genes)],
    ).to_pandas()
    exons = pq.read_table(args.exons, filters=[("gene_name", "in", genes)]).to_pandas()
    cds = pq.read_table(args.cds, filters=[("gene_name", "in", genes)]).to_pandas()
    fasta = pysam.FastaFile(str(args.fasta))
    reconstruction = import_reconstruction(args.reconstruction_implementation)
    models = {}
    for transcript_id, group in transcripts.groupby("transcript_id", sort=False):
        metadata = group.iloc[0].to_dict()
        tx_exons = exons[exons["transcript_id"] == transcript_id]
        tx_cds = cds[cds["transcript_id"] == transcript_id]
        if not tx_exons.empty:
            models[str(transcript_id)] = reconstruction.TranscriptModel(metadata, tx_exons, tx_cds, fasta)

    carrier_table = pd.read_csv(args.carriers, sep="\t", compression="gzip", dtype=str).fillna("")
    carrier_sets = {}
    for row in carrier_table.itertuples(index=False):
        hist = str(row.histology).replace("TCGA-", "")
        carrier_sets.setdefault((hist, str(row.variant_uid), str(row.gene)), set()).add(str(row.patient))

    sparse = pd.read_csv(args.sparse_counts, sep="\t", compression="gzip", low_memory=False)
    sparse["histology"] = sparse["histology"].astype(str).str.replace("TCGA-", "", regex=False)
    totals = pd.read_csv(args.cluster_totals, sep="\t", compression="gzip", low_memory=False)
    totals["histology"] = totals["histology"].astype(str).str.replace("TCGA-", "", regex=False)
    count_key = ["histology", "patient", "gene", "junction_id", "cluster_id"]
    if sparse.duplicated(count_key).any():
        raise ValueError("v2.2 sparse counts are not unique at the expected M3 lookup key")
    count_lookup = {
        tuple(row[key] for key in count_key): float(row["count"])
        for row in sparse.to_dict("records")
    }
    total_key = ["histology", "patient", "gene", "cluster_id"]
    if totals.duplicated(total_key).any():
        raise ValueError("v2.2 cluster totals are not unique at the expected M3 lookup key")
    total_lookup = {
        tuple(row[key] for key in total_key): float(row["cluster_total"])
        for row in totals.to_dict("records")
    }
    cluster_junctions = {}
    for row in sparse.to_dict("records"):
        cluster_junctions.setdefault((row["histology"], row["gene"], row["cluster_id"], row["junction_id"]), True)

    chain_rows = []
    stat_rows = []
    summary_rows = []
    model_failures = []
    for row in events.to_dict("records"):
        histology = str(row["histologies"]).split("|")[0].replace("TCGA-", "")
        gene = str(row["gene"])
        selected_id = strip_version(row["selected_transcript_v23"])
        model = models.get(selected_id)
        if model is None:
            model_failures.append({"gene": gene, "transcript": selected_id, "reason": "MISSING_GENCODE50_MODEL"})
            continue
        chrom, intron_start, intron_end, strand = parse_junction(row["junction_id"])
        if chrom != model.chrom or strand != model.strand:
            model_failures.append({"gene": gene, "transcript": selected_id, "reason": "MODEL_JUNCTION_CHROM_OR_STRAND_MISMATCH"})
            continue
        flanks = choose_flanking_exons(model, intron_start, intron_end)
        if flanks is None:
            model_failures.append({"gene": gene, "transcript": selected_id, "reason": "AMBIGUOUS_OR_MISSING_REFERENCE_FLANKS"})
            continue
        up_index, down_index = flanks
        reciprocal_ids = [reference_junction_id(model, index) for index in range(up_index, down_index)]
        reference_ids = [reference_junction_id(model, index) for index in range(len(model.exons) - 1)]
        cluster_id = str(row["cluster_id_v23"] or row["cluster_id"])
        altered_id = str(row["junction_id"])
        chain_rows.append({
            "histology": histology,
            "variant_uid": row["variant_uid"],
            "gene": gene,
            "junction_id": altered_id,
            "cluster_id": cluster_id,
            "selected_transcript": selected_id,
            "event_class": row["event_class_v23"],
            "altered_splice_type": row["altered_splice_type_v23"],
            "altered_exon_chain": row["altered_exon_chain_v23"],
            "reference_exon_chain": "|".join(reference_ids),
            "reciprocal_reference_junctions": "|".join(reciprocal_ids),
            "reciprocal_reference_junction_count": len(reciprocal_ids),
            "reference_flanking_exon_indices": f"{up_index}:{down_index}",
            "reference_asset_release": "GENCODE50_index_v2.2 / GRCh38 primary assembly",
        })

        candidate_patients = [
            (patient, total)
            for (h, patient, g, c), total in total_lookup.items()
            if h == histology and g == gene and c == cluster_id and total >= 10
        ]
        carriers = carrier_sets.get((histology, str(row["variant_uid"]), gene), set())
        carrier_patients = [patient for patient, _total in candidate_patients if patient in carriers]
        noncarrier_patients = [patient for patient, _total in candidate_patients if patient not in carriers]
        junctions = [("ALTERED", altered_id)] + [(f"REFERENCE_{index + 1}", junction_id) for index, junction_id in enumerate(reciprocal_ids)]
        event_stat_map = {}
        for role, junction_id in junctions:
            present_in_cluster = (histology, gene, cluster_id, junction_id) in cluster_junctions
            carrier_values = np.array([
                count_lookup.get((histology, patient, gene, junction_id, cluster_id), 0.0) / total_lookup[(histology, patient, gene, cluster_id)]
                for patient in carrier_patients
            ], dtype=float)
            noncarrier_values = np.array([
                count_lookup.get((histology, patient, gene, junction_id, cluster_id), 0.0) / total_lookup[(histology, patient, gene, cluster_id)]
                for patient in noncarrier_patients
            ], dtype=float)
            carrier_stats = stats(carrier_values)
            noncarrier_stats = stats(noncarrier_values)
            event_stat_map[role] = {
                "carrier_mean": carrier_stats["mean"],
                "noncarrier_mean": noncarrier_stats["mean"],
                "present_in_cluster": present_in_cluster,
            }
            stat_rows.append({
                **chain_rows[-1],
                "junction_role": role,
                "measured_junction_id": junction_id,
                "present_in_candidate_cluster": present_in_cluster,
                "n_carrier_evaluable": carrier_stats["n"],
                "n_noncarrier_evaluable": noncarrier_stats["n"],
                "carrier_mean_usage": carrier_stats["mean"],
                "noncarrier_mean_usage": noncarrier_stats["mean"],
                "carrier_median_usage": carrier_stats["median"],
                "noncarrier_median_usage": noncarrier_stats["median"],
                "carrier_n_nonzero": carrier_stats["n_nonzero"],
                "noncarrier_n_nonzero": noncarrier_stats["n_nonzero"],
                "carrier_total_junction_count": carrier_stats["sum"],
                "noncarrier_total_junction_count": noncarrier_stats["sum"],
                "carrier_status": "PASS" if carrier_stats["n"] else "NO_EVALUABLE_CARRIERS",
                "coverage_rule": "cluster_total >= 10",
            })

        altered = event_stat_map["ALTERED"]
        reference_stats = [event_stat_map[role] for role, _jid in junctions if role != "ALTERED"]
        directional = [
            altered["carrier_mean"] > altered["noncarrier_mean"]
            and ref["carrier_mean"] < ref["noncarrier_mean"]
            for ref in reference_stats
            if ref["present_in_cluster"]
        ]
        if not reciprocal_ids:
            status = "NOT_EVALUABLE_NO_REFERENCE_JUNCTION"
        elif not any(ref["present_in_cluster"] for ref in reference_stats):
            status = "NOT_EVALUABLE_REFERENCE_NOT_IN_CLUSTER"
        elif any(directional):
            status = "RECIPROCAL_DIRECTION_CONCORDANT_DESCRIPTIVE"
        elif altered["carrier_mean"] > altered["noncarrier_mean"]:
            status = "ALTERED_ENRICHED_NO_RECIPROCAL_DEPLETION"
        else:
            status = "NO_RECIPROCAL_DIRECTIONAL_SUPPORT"
        summary_rows.append({
            **chain_rows[-1],
            "n_carriers_total": len(carriers),
            "n_carrier_evaluable": len(carrier_patients),
            "n_noncarrier_evaluable": len(noncarrier_patients),
            "altered_carrier_mean_usage": altered["carrier_mean"],
            "altered_noncarrier_mean_usage": altered["noncarrier_mean"],
            "altered_usage_delta_carrier_minus_noncarrier": altered["carrier_mean"] - altered["noncarrier_mean"],
            "reference_junctions_present_in_cluster": sum(ref["present_in_cluster"] for ref in reference_stats),
            "reciprocal_direction_concordant_count": sum(directional),
            "m3_status": status,
            "m3_interpretation": "descriptive reciprocal splice evidence; not a causal or clinical endpoint",
        })

    if model_failures:
        raise RuntimeError(json.dumps({"m3_status": "HOLD_MODEL_FAILURE", "model_failures": model_failures}, ensure_ascii=False))
    if len(chain_rows) != 9 or len(summary_rows) != 9:
        raise RuntimeError(f"M3 expected 9 event summaries, got chains={len(chain_rows)} summaries={len(summary_rows)}")

    compression = {"method": "gzip", "mtime": 0}
    chain_path = args.outdir / "V23_M3_REFERENCE_EXON_CHAIN_MAP.tsv.gz"
    stat_path = args.outdir / "V23_M3_RECIPROCAL_JUNCTION_STATS.tsv.gz"
    summary_path = args.outdir / "V23_M3_RECIPROCAL_EVENT_SUMMARY.tsv"
    pd.DataFrame(chain_rows).to_csv(chain_path, sep="\t", index=False, compression=compression)
    pd.DataFrame(stat_rows).to_csv(stat_path, sep="\t", index=False, compression=compression)
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)

    reference_audit = args.reference_sha256.read_text(encoding="utf-8")
    audit = {
        "status": "PASS",
        "scope": "frozen_9_R2",
        "events": 9,
        "reference_release": "GENCODE50_index_v2.2 / GRCh38 primary assembly",
        "reference_asset_sha256_manifest": str(args.reference_sha256),
        "reference_asset_sha256_manifest_sha256": hashlib.sha256(reference_audit.encode()).hexdigest(),
        "reconstruction_implementation": str(args.reconstruction_implementation),
        "reconstruction_implementation_sha256": hashlib.sha256(args.reconstruction_implementation.read_bytes()).hexdigest(),
        "same_reference_exon_chain_for_all_events": True,
        "coverage_rule": "candidate cluster_total >= 10",
        "reciprocal_definition": "canonical GENCODE50 junctions spanning the selected transcript flanking exon interval of the altered junction",
        "functional_fields_used_to_define_reciprocal_junctions": False,
        "clinical_or_survival_analysis_run": False,
        "m3_status_counts": pd.Series([row["m3_status"] for row in summary_rows]).value_counts().to_dict(),
        "outputs": [str(chain_path), str(stat_path), str(summary_path)],
    }
    (args.outdir / "V23_M3_AUDIT.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
