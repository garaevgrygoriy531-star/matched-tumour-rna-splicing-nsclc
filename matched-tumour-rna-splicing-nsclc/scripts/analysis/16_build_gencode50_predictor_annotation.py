#!/usr/bin/env python3
"""Build outcome-blind GENCODE50 predictor annotation assets.

The builder consumes only a GENCODE GTF and a FASTA index.  It selects the
single ``Ensembl_canonical`` transcript per gene, validates the parent/exon
chain, and writes deterministic SpliceAI custom annotation, canonical GTF,
and transcript mapping assets.  No RNA, variant, M2, or M3 result is read.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


GTf_FEATURES = {"gene", "transcript", "exon"}
SELECTION_TAG = "Ensembl_canonical"
SELECTION_RULE = "Ensembl_canonical_tag"
CUSTOM_NAME = "V23_GENCODE50_SPLICEAI_CUSTOM.txt"
CANONICAL_GTF = "V23_GENCODE50_ENSEMBL_CANONICAL.gtf.gz"
TRANSCRIPT_TABLE = "V23_GENCODE50_CANONICAL_TRANSCRIPTS.tsv"
SHA_TABLE = "V23_M4_M5_ANNOTATION_SHA256.tsv"
AUDIT_NAME = "V23_M4_M5_ANNOTATION_BUILD_AUDIT.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def input_file_metadata(path: Path) -> dict[str, Any]:
    return {"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)}


def load_fai(path: Path) -> tuple[list[str], dict[str, int]]:
    order: list[str] = []
    lengths: dict[str, int] = {}
    with path.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            text = line.rstrip("\r\n")
            if not text:
                continue
            fields = text.split("\t")
            if len(fields) < 2 or not fields[0]:
                raise ValueError(f"invalid FASTA .fai line {line_number}: {text!r}")
            name = fields[0]
            if name in lengths:
                raise ValueError(f"duplicate FASTA contig in .fai: {name}")
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"invalid FASTA contig length on .fai line {line_number}") from exc
            if length < 1:
                raise ValueError(f"non-positive FASTA contig length on .fai line {line_number}")
            order.append(name)
            lengths[name] = length
    if not order:
        raise ValueError("FASTA .fai has no contigs")
    return order, lengths


def parse_attributes(text: str) -> dict[str, list[str]]:
    """Parse standard GTF key/value attributes while retaining repeated keys."""
    attributes: dict[str, list[str]] = defaultdict(list)
    for field in text.strip().rstrip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        match = re.match(r"^(\S+)\s+\"((?:[^\"]|\\.)*)\"$", field)
        if match is None:
            match = re.match(r"^(\S+)\s+(.+)$", field)
        if match is None:
            raise ValueError(f"unparseable GTF attribute field: {field!r}")
        key, value = match.group(1), match.group(2)
        value = value.replace('\\"', '"').replace('\\\\', '\\')
        attributes[key].append(value)
    return dict(attributes)


def first_attr(attributes: dict[str, list[str]], key: str, default: str = "") -> str:
    values = attributes.get(key, [])
    return str(values[0]) if values else default


def all_attr(attributes: dict[str, list[str]], key: str) -> list[str]:
    return [str(value) for value in attributes.get(key, [])]


def strip_version(transcript_id: str) -> tuple[str, str]:
    value = str(transcript_id)
    match = re.match(r"^(.*)\.([0-9]+)$", value)
    if match:
        return match.group(1), match.group(2)
    return value, ""


def iter_gtf(path: Path) -> Iterable[dict[str, Any]]:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            raw = line.rstrip("\r\n")
            if not raw or raw.startswith("#"):
                continue
            fields = raw.split("\t")
            if len(fields) != 9:
                raise ValueError(f"GTF line {line_number} does not have 9 tab-separated fields")
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"non-integer GTF coordinates on line {line_number}") from exc
            if start < 1 or end < 1 or start > end:
                raise ValueError(f"invalid GTF coordinates on line {line_number}: {start}-{end}")
            yield {
                "line_number": line_number,
                "raw": raw,
                "chrom": fields[0],
                "source": fields[1],
                "feature": fields[2],
                "start": start,
                "end": end,
                "score": fields[5],
                "strand": fields[6],
                "frame": fields[7],
                "attributes": parse_attributes(fields[8]),
            }


def sort_key(record: dict[str, Any], chrom_rank: dict[str, int]) -> tuple[Any, ...]:
    feature_rank = {"gene": 0, "transcript": 1, "exon": 2}
    return (
        chrom_rank[record["chrom"]],
        int(record["start"]),
        str(record.get("transcript_id", "")),
        feature_rank.get(record["feature"], 9),
        int(record["end"]),
        int(record.get("line_number", 0)),
    )


def require_attr(attributes: dict[str, list[str]], key: str, feature: str, line_number: int) -> str:
    values = attributes.get(key, [])
    if len(values) != 1 or not values[0]:
        raise ValueError(
            f"{feature} line {line_number} requires exactly one non-empty {key!r} attribute; "
            f"observed {values!r}"
        )
    return str(values[0])


def deterministic_gzip_text(path: Path, lines: Iterable[str]) -> None:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(fileobj=raw_handle, mode="wb", filename="", mtime=0) as gzip_handle:
            with io.TextIOWrapper(gzip_handle, encoding="utf-8", newline="\n") as text_handle:
                for line in lines:
                    text_handle.write(line.rstrip("\r\n") + "\n")


def write_tsv(path: Path, columns: list[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(column, "")) for column in columns) + "\n")


def build_assets(gtf: Path, fai: Path, outdir: Path, audit_dir: Path) -> dict[str, Any]:
    fai_order, fai_lengths = load_fai(fai)
    chrom_rank = {chrom: index for index, chrom in enumerate(fai_order)}
    retained_features = Counter()
    dropped_chrom = Counter()
    ignored_features = Counter()
    missing_attributes = Counter()
    selection_drops = Counter()
    transcript_records: dict[str, dict[str, Any]] = {}
    transcript_genes: dict[str, set[str]] = defaultdict(set)
    canonical_by_gene: dict[str, set[str]] = defaultdict(set)
    gene_records: dict[str, dict[str, Any]] = {}
    gene_name_by_id: dict[str, str] = {}
    pass1_records = 0

    # Pass 1: collect only gene/transcript metadata and determine the frozen
    # canonical transcript set. Exon lines are intentionally not retained yet.
    for record in iter_gtf(gtf):
        pass1_records += 1
        feature = record["feature"]
        if feature not in GTf_FEATURES:
            ignored_features[feature] += 1
            continue
        chrom = record["chrom"]
        if chrom not in chrom_rank:
            dropped_chrom[feature] += 1
            if feature == "transcript" and SELECTION_TAG in all_attr(record["attributes"], "tag"):
                selection_drops["canonical_transcript_contig_not_in_fai"] += 1
            continue
        retained_features[feature] += 1
        attributes = record["attributes"]
        gene_id = require_attr(attributes, "gene_id", feature, record["line_number"])
        if feature == "gene":
            gene_name = first_attr(attributes, "gene_name")
            if not gene_name:
                missing_attributes["gene.gene_name"] += 1
            current = dict(record)
            current.update({"gene_id": gene_id, "gene_name": gene_name})
            previous = gene_records.get(gene_id)
            if previous is not None and previous["raw"] != current["raw"]:
                raise ValueError(f"conflicting duplicate gene records for {gene_id}")
            gene_records[gene_id] = current
            gene_name_by_id[gene_id] = gene_name
            continue
        if feature == "transcript":
            raw_transcript_id = require_attr(attributes, "transcript_id", feature, record["line_number"])
            base_id, derived_version = strip_version(raw_transcript_id)
            transcript_version = first_attr(attributes, "transcript_version", derived_version)
            transcript_type = first_attr(attributes, "transcript_type")
            if not transcript_type:
                missing_attributes["transcript.transcript_type"] += 1
            gene_name = first_attr(attributes, "gene_name", gene_name_by_id.get(gene_id, ""))
            current_tags = all_attr(attributes, "tag")
            current = dict(record)
            current.update({
                "gene_id": gene_id,
                "gene_name": gene_name,
                "transcript_id": base_id,
                "gtf_transcript_id": raw_transcript_id,
                "transcript_version": transcript_version,
                "transcript_type": transcript_type,
                "tags": current_tags,
            })
            transcript_genes[base_id].add(gene_id)
            previous = transcript_records.get(base_id)
            if previous is not None:
                invariant = ("gene_id", "chrom", "strand", "start", "end", "transcript_type")
                if any(previous.get(field) != current.get(field) for field in invariant):
                    raise ValueError(f"conflicting duplicate transcript records for {base_id}")
                raise ValueError(f"duplicate transcript_id/#NAME record for {base_id}")
            else:
                transcript_records[base_id] = current
            if SELECTION_TAG in current["tags"]:
                canonical_by_gene[gene_id].add(base_id)
            continue
        if feature == "exon":
            require_attr(attributes, "transcript_id", feature, record["line_number"])

    ambiguous = {gene: sorted(ids) for gene, ids in canonical_by_gene.items() if len(ids) > 1}
    if ambiguous:
        examples = dict(list(ambiguous.items())[:20])
        raise ValueError("more than one distinct Ensembl_canonical transcript per gene: " + json.dumps(examples, ensure_ascii=False))
    transcript_gene_conflicts = {tx: sorted(genes) for tx, genes in transcript_genes.items() if len(genes) > 1}
    if transcript_gene_conflicts:
        examples = dict(list(transcript_gene_conflicts.items())[:20])
        raise ValueError("one transcript_id/#NAME belongs to multiple genes: " + json.dumps(examples, ensure_ascii=False))
    selected_ids = {next(iter(ids)) for ids in canonical_by_gene.values() if ids}
    if not selected_ids:
        raise ValueError("no Ensembl_canonical transcripts found on FASTA contigs")
    if len(selected_ids) != len({str(item) for item in selected_ids}):
        raise ValueError("transcript_id/#NAME is not unique")
    selected_gene_ids = {transcript_records[tx]["gene_id"] for tx in selected_ids if tx in transcript_records}
    missing_gene_parents = sorted(gene for gene in selected_gene_ids if gene not in gene_records)
    if missing_gene_parents:
        raise ValueError("selected transcript gene parents are missing: " + ", ".join(missing_gene_parents[:20]))
    parent_mismatches = []
    for tx in selected_ids:
        transcript = transcript_records[tx]
        gene = gene_records[transcript["gene_id"]]
        if (
            gene["chrom"] != transcript["chrom"]
            or gene["strand"] != transcript["strand"]
            or not (gene["start"] <= transcript["start"] <= transcript["end"] <= gene["end"])
        ):
            parent_mismatches.append(tx)
    if parent_mismatches:
        raise ValueError("selected transcript/gene parent chrom/strand/bounds mismatch: " + ", ".join(parent_mismatches[:20]))

    invalid_selected_coordinates = []
    for tx in selected_ids:
        transcript = transcript_records[tx]
        gene = gene_records[transcript["gene_id"]]
        if (
            transcript["strand"] not in {"+", "-"}
            or transcript["end"] > fai_lengths[transcript["chrom"]]
            or gene["strand"] not in {"+", "-"}
            or gene["end"] > fai_lengths[gene["chrom"]]
        ):
            invalid_selected_coordinates.append(tx)
    if invalid_selected_coordinates:
        raise ValueError(
            "selected transcripts have invalid strand or exceed FASTA contig bounds: "
            + ", ".join(sorted(invalid_selected_coordinates)[:20])
        )

    # Pass 2: collect only exons belonging to the frozen selected transcript
    # set and retain the original GTF records for deterministic canonical GTF.
    exons_by_transcript: dict[str, list[dict[str, Any]]] = defaultdict(list)
    exon_coords: dict[str, set[tuple[int, int]]] = defaultdict(set)
    exon_ids: dict[str, set[str]] = defaultdict(set)
    selected_transcript_lines: dict[str, dict[str, Any]] = {}
    pass2_records = 0
    selected_missing_or_inconsistent = []
    for record in iter_gtf(gtf):
        pass2_records += 1
        feature = record["feature"]
        if feature not in GTf_FEATURES:
            continue
        attributes = record["attributes"]
        if record["chrom"] not in chrom_rank:
            if feature == "exon":
                raw_id = first_attr(attributes, "transcript_id")
                base_id, _version = strip_version(raw_id)
                if base_id in selected_ids:
                    selected_missing_or_inconsistent.append({
                        "line_number": record["line_number"],
                        "transcript_id": base_id,
                        "reason": "selected_exon_chrom_not_in_fai",
                    })
            continue
        if feature == "transcript":
            raw_id = first_attr(attributes, "transcript_id")
            base_id, _version = strip_version(raw_id)
            if base_id not in selected_ids:
                continue
            selected_transcript_lines[base_id] = record
            selected_record = transcript_records[base_id]
            if record["raw"] != selected_record["raw"]:
                raise ValueError(f"selected transcript changed between GTF passes for {base_id}")
        elif feature == "exon":
            raw_id = require_attr(attributes, "transcript_id", feature, record["line_number"])
            gene_id = require_attr(attributes, "gene_id", feature, record["line_number"])
            base_id, _version = strip_version(raw_id)
            if base_id not in selected_ids:
                continue
            exon_id = require_attr(attributes, "exon_id", feature, record["line_number"])
            selected = transcript_records[base_id]
            if gene_id != selected["gene_id"] or record["chrom"] != selected["chrom"] or record["strand"] != selected["strand"]:
                selected_missing_or_inconsistent.append({
                    "line_number": record["line_number"],
                    "transcript_id": base_id,
                    "reason": "exon_parent_chrom_or_strand_mismatch",
                })
            if not (selected["start"] <= record["start"] <= record["end"] <= selected["end"]):
                selected_missing_or_inconsistent.append({
                    "line_number": record["line_number"],
                    "transcript_id": base_id,
                    "reason": "exon_outside_transcript_bounds",
                })
            coord = (int(record["start"]), int(record["end"]))
            if coord in exon_coords[base_id]:
                raise ValueError(f"duplicate exon coordinate for selected transcript {base_id}: {coord}")
            if exon_id in exon_ids[base_id]:
                raise ValueError(f"duplicate exon_id for selected transcript {base_id}: {exon_id}")
            exon_coords[base_id].add(coord)
            exon_ids[base_id].add(exon_id)
            exon = dict(record)
            exon["transcript_id"] = base_id
            exons_by_transcript[base_id].append(exon)

    if selected_missing_or_inconsistent:
        raise ValueError("selected transcript exon validation failed: " + json.dumps(selected_missing_or_inconsistent[:20], ensure_ascii=False))
    missing_transcript_lines = sorted(selected_ids - set(selected_transcript_lines))
    if missing_transcript_lines:
        raise ValueError("selected transcript lines missing in second GTF pass: " + ", ".join(missing_transcript_lines[:20]))
    missing_exons = sorted(tx for tx in selected_ids if not exons_by_transcript.get(tx))
    if missing_exons:
        raise ValueError("selected transcripts without exons: " + ", ".join(missing_exons[:20]))

    selected_records: list[dict[str, Any]] = []
    transcript_rows: list[dict[str, Any]] = []
    custom_rows: list[dict[str, Any]] = []
    for tx in selected_ids:
        record = transcript_records[tx]
        exons = sorted(exons_by_transcript[tx], key=lambda item: (int(item["start"]), int(item["end"]), int(item["line_number"])))
        exon_starts = ",".join(str(int(exon["start"]) - 1) for exon in exons)
        exon_ends = ",".join(str(int(exon["end"])) for exon in exons)
        metadata = {
            "spliceai_name": tx,
            "gene_id": record["gene_id"],
            "gene_name": record.get("gene_name") or gene_records[record["gene_id"]].get("gene_name", ""),
            "transcript_id": tx,
            "gtf_transcript_id": record["gtf_transcript_id"],
            "transcript_version": record["transcript_version"],
            "transcript_type": record["transcript_type"],
            "chrom": record["chrom"],
            "strand": record["strand"],
            "tx_start_1based": int(record["start"]),
            "tx_end_1based": int(record["end"]),
            "n_exons": len(exons),
            "selection_rule": SELECTION_RULE,
        }
        transcript_rows.append(metadata)
        custom_rows.append({
            "name": tx,
            "chrom": record["chrom"],
            "strand": record["strand"],
            "tx_start": int(record["start"]) - 1,
            "tx_end": int(record["end"]),
            "exon_start": exon_starts,
            "exon_end": exon_ends,
        })
        tx_line = dict(selected_transcript_lines[tx])
        tx_line["transcript_id"] = tx
        selected_records.append(tx_line)
        for exon in exons:
            selected_records.append(exon)
    for gene_id in selected_gene_ids:
        selected_records.append(gene_records[gene_id])

    transcript_rows.sort(key=lambda row: (chrom_rank[row["chrom"]], int(row["tx_start_1based"]), str(row["transcript_id"])))
    custom_rows.sort(key=lambda row: (chrom_rank[row["chrom"]], int(row["tx_start"]), str(row["name"])))
    selected_records.sort(key=lambda row: sort_key(row, chrom_rank))

    outdir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)
    custom_path = outdir / CUSTOM_NAME
    canonical_gtf_path = outdir / CANONICAL_GTF
    transcript_path = outdir / TRANSCRIPT_TABLE
    deterministic_gzip_text(
        canonical_gtf_path,
        (record["raw"] for record in selected_records),
    )
    with custom_path.open("wt", encoding="utf-8", newline="") as handle:
        handle.write("#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n")
        for row in custom_rows:
            handle.write("\t".join(str(row[key]) for key in ["name", "chrom", "strand", "tx_start", "tx_end", "exon_start", "exon_end"]) + "\n")
    transcript_columns = [
        "spliceai_name", "gene_id", "gene_name", "transcript_id", "gtf_transcript_id", "transcript_version",
        "transcript_type", "chrom", "strand", "tx_start_1based", "tx_end_1based", "n_exons", "selection_rule",
    ]
    write_tsv(transcript_path, transcript_columns, transcript_rows)

    output_paths = [custom_path, canonical_gtf_path, transcript_path]
    script_path = Path(__file__).resolve()
    sha_rows = [
        {
            "role": "script",
            "name": script_path.name,
            "path": str(script_path),
            "size_bytes": script_path.stat().st_size,
            "sha256": sha256_file(script_path),
        }
    ]
    sha_rows.extend(
        {
            "role": "input",
            "name": path.name,
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in (gtf, fai)
    )
    sha_rows.extend(
        {
            "role": "output",
            "name": path.name,
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in output_paths
    )
    sha_path = audit_dir / SHA_TABLE
    write_tsv(sha_path, ["role", "name", "path", "size_bytes", "sha256"], sha_rows)

    audit = {
        "status": "PASS",
        "gencode_release": "GENCODE50",
        "selection_tag": SELECTION_TAG,
        "selection_rule": SELECTION_RULE,
        "outcome_blind_contract": {
            "outcome_fields_used": [],
            "variant_fields_used": [],
            "allowed_inputs": ["GENCODE release 50 GTF", "reference FASTA .fai"],
        },
        "coordinate_contract": {
            "transcript_start": "GTF 1-based inclusive start minus 1 -> SpliceAI 0-based TX_START",
            "transcript_end": "GTF 1-based inclusive end unchanged -> SpliceAI TX_END",
            "exon_start": "GTF 1-based inclusive start minus 1 -> SpliceAI EXON_START",
            "exon_end": "GTF 1-based inclusive end unchanged -> SpliceAI EXON_END",
            "exon_order": "ascending genomic start, independent of strand",
        },
        "inputs": {"gtf": input_file_metadata(gtf), "fasta_fai": input_file_metadata(fai)},
        "fai_contigs": {"count": len(fai_order), "order": fai_order},
        "counts": {
            "gtf_records_seen_pass1": pass1_records,
            "gtf_records_seen_pass2": pass2_records,
            "retained_features_pass1": dict(retained_features),
            "dropped_chrom_by_feature": dict(dropped_chrom),
            "ignored_features": dict(ignored_features),
            "genes_retained": len(gene_records),
            "transcripts_retained": len(transcript_records),
            "selected_genes": len(selected_gene_ids),
            "selected_transcripts": len(selected_ids),
            "selected_exons": sum(len(value) for value in exons_by_transcript.values()),
        },
        "drops": {
            "records_on_contigs_not_in_fai": dict(dropped_chrom),
            "selection": {
                **dict(selection_drops),
                "retained_genes_without_ensembl_canonical_transcript": len(gene_records) - len(selected_gene_ids),
            },
        },
        "missing_attributes": dict(missing_attributes),
        "validation_flags": {
            "one_ensembl_canonical_transcript_per_selected_gene": True,
            "transcript_id_unique": len(transcript_records) == len(set(transcript_records)),
            "spliceai_name_unique": len(selected_ids) == len(set(selected_ids)),
            "selected_transcripts_have_exons": True,
            "exons_within_transcript_and_consistent_chrom_strand": True,
            "exon_coordinates_and_ids_unique_within_transcript": True,
            "gene_parents_complete": True,
            "selected_coordinates_within_fasta_contig_bounds": True,
            "only_fai_contigs_retained": True,
            "outcome_blind_inputs_only": True,
        },
        "outputs": {
            "spliceai_custom": {"path": str(custom_path), "size_bytes": custom_path.stat().st_size, "sha256": sha256_file(custom_path)},
            "canonical_gtf": {"path": str(canonical_gtf_path), "size_bytes": canonical_gtf_path.stat().st_size, "sha256": sha256_file(canonical_gtf_path)},
            "canonical_transcripts": {"path": str(transcript_path), "size_bytes": transcript_path.stat().st_size, "sha256": sha256_file(transcript_path)},
            "sha256_tsv": {"path": str(sha_path), "size_bytes": sha_path.stat().st_size, "sha256": sha256_file(sha_path)},
        },
        "script_sha256": sha256_file(script_path),
    }
    audit_path = audit_dir / AUDIT_NAME
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return audit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--fasta-fai", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--audit-dir", required=True, type=Path)
    args = parser.parse_args()
    for path in (args.gtf, args.fasta_fai):
        if not path.is_file():
            raise FileNotFoundError(path)
    audit = build_assets(args.gtf, args.fasta_fai, args.outdir, args.audit_dir)
    print(json.dumps({"status": audit["status"], "counts": audit["counts"], "script_sha256": audit["script_sha256"]}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
