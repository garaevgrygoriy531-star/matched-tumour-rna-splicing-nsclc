#!/usr/bin/env python3
"""Independent protein/PTC/NMD/domain refinement for the frozen R2-9.

This module deliberately has a very small input contract.  It reads the
frozen-9 key table, the post-M2/M3 R2 matrix (for identity and the pre-existing
NMD/frame labels), the frozen full reconstruction table, the original GTF,
the GENCODE50 CDS parquet, and the frozen GRCh38 FASTA/FAI.  Reference and
altered CDS/AA sequences are rebuilt from the CDS coordinates and genomic
sequence; values from the old domain placeholder are never used.

The altered splice chain is a coordinate-only chain.  It contains no somatic
allele, so the resulting sequence is explicitly an *altered ORF if
translated*, not an observed or experimentally expressed protein.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

try:  # The local workspace may not carry the server-only data libraries.
    import pyarrow.parquet as pq
except ImportError:  # pragma: no cover - exercised by local self-test only
    pq = None  # type: ignore[assignment]
try:
    import pysam
except ImportError:  # pragma: no cover - exercised by local self-test only
    pysam = None  # type: ignore[assignment]


EXPECTED_N = 9
EXPECTED_RECON_SHA256 = "49c3dfa6a3125013e34772222676edc4bc368231c4477ad6865630443236c4ef"
EXPECTED_GTF_SHA256 = "83fba3e9b03f0b8c958f3595c6c350adc55f468abf8b0e47b6d5284cfe13a453"
EXPECTED_CDS_SHA256 = "d42ce27e5a50743f3e4d96434f28365d7f8f7dcd8ab13015ef0e045dc6303f73"
EXPECTED_FASTA_SHA256 = "e49b92b3e4f321bf254c042f25b726d9931c4d74c7523e8b6bb530e63b0cfd4b"
EXPECTED_FAI_SHA256 = "a2c323ea4cff34d7123ace4578f7e122b2d2f5a22f40dc23eb8b97d17723d169"

DOMAIN_REASON = "DOMAIN_UNRESOLVED_NO_EXACT_SEQUENCE_MATCHED_DOMAIN_MAP"
CONSEQUENCE_SCOPE = "PREDICTED_ALTERED_ORF_IF_TRANSLATED"
PROTEIN_EXPRESSION_STATUS = "NOT_ESTABLISHED"
CIS_STATUS = "STATISTICAL_LINK_NOT_MOLECULARLY_PHASED"

EXPECTED_LENGTHS = {
    "chr12:25227342:T:A": (188, 40),
    "chr19:10489648:C:A": (624, 456),
    "chr19:1220373:G:T": (433, 155),
    "chr19:45420410:C:T": (297, 151),
    "chr6:43443117:G:T": (1492, 822),
    "chr13:48307406:G:T": (928, 67),
    "chr17:7673570:T:C": (393, 310),
    "chr17:7676379:T:A": (393, 35),
    "chr2:224481892:C:A": (768, 632),
}

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TGT": "C", "TGC": "C",
    "TGG": "W", "CTT": "L", "CTC": "L", "CTA": "L",
    "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P",
    "CCG": "P", "CAT": "H", "CAC": "H", "CAA": "Q",
    "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R",
    "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T",
    "ACG": "T", "AAT": "N", "AAC": "N", "AAA": "K",
    "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R",
    "AGG": "R", "GTT": "V", "GTC": "V", "GTA": "V",
    "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A",
    "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E",
    "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G",
    "GGG": "G", "TAA": "*", "TAG": "*", "TGA": "*",
}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metadata(path: Path, expected: str | None = None) -> dict[str, Any]:
    observed = sha256_file(path)
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": observed,
        "expected_sha256": expected,
        "sha256_match": expected is None or observed == expected,
    }


def strip_version(value: str) -> str:
    return str(value or "").split(".", 1)[0]


def parse_attrs(text: str) -> dict[str, list[str]]:
    """Parse GTF attributes while retaining repeated keys."""
    result: dict[str, list[str]] = {}
    for field in text.strip().rstrip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        match = re.match(r"^(\S+)\s+\"((?:[^\"]|\\.)*)\"$", field)
        if match:
            result.setdefault(match.group(1), []).append(match.group(2))
        else:
            # GENCODE uses quoted values, but fail closed if a selected row
            # has an unparseable attribute rather than silently mapping it.
            parts = field.split(None, 1)
            if len(parts) == 2:
                result.setdefault(parts[0], []).append(parts[1].strip('"'))
    return result


def first_attr(attrs: dict[str, list[str]], key: str) -> str:
    values = attrs.get(key, [])
    return values[0] if len(values) == 1 else ""


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require_columns(rows: list[dict[str, Any]], columns: Iterable[str], label: str) -> None:
    if not rows:
        raise RuntimeError(f"{label} is empty")
    observed = set(rows[0])
    missing = sorted(set(columns) - observed)
    if missing:
        raise RuntimeError(f"{label} missing columns: {missing}")


def parse_frozen_keys(path: Path) -> list[dict[str, str]]:
    rows = read_tsv(path)
    require_columns(rows, ["variant_uid", "gene", "junction_id", "selected_transcript"], "frozen keys")
    if len(rows) != EXPECTED_N:
        raise RuntimeError(f"frozen keys must contain exactly {EXPECTED_N} rows, observed {len(rows)}")
    keys = {(r["variant_uid"], r["gene"], r["junction_id"]) for r in rows}
    if len(keys) != EXPECTED_N:
        raise RuntimeError("frozen keys contain duplicate event identities")
    return rows


def parse_gtf(path: Path, wanted_tx: set[str]) -> tuple[
    dict[str, dict[str, Any]],
    dict[str, list[dict[str, Any]]],
    dict[str, list[dict[str, Any]]],
]:
    transcripts: dict[str, dict[str, Any]] = {}
    cds_rows: dict[str, list[dict[str, Any]]] = {tx: [] for tx in wanted_tx}
    stop_rows: dict[str, list[dict[str, Any]]] = {tx: [] for tx in wanted_tx}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise RuntimeError(f"invalid GTF line {line_number}: expected 9 fields")
            chrom, _source, feature, start, end, _score, strand, frame, attr_text = fields
            attrs = parse_attrs(attr_text)
            raw_tx = first_attr(attrs, "transcript_id")
            tx = strip_version(raw_tx)
            if tx not in wanted_tx:
                continue
            if chrom.startswith("chr") is False:
                # The frozen data uses chr-prefixed contigs.  Keep the check
                # explicit, because a contig alias could silently corrupt a
                # CDS reconstruction.
                raise RuntimeError(f"selected transcript {tx} has non-chr contig {chrom}")
            transcript_version = first_attr(attrs, "transcript_version") or (
                raw_tx.split(".", 1)[1] if "." in raw_tx else ""
            )
            record = {
                "tx": tx,
                "full_id": f"{tx}.{transcript_version}" if transcript_version else tx,
                "gene": first_attr(attrs, "gene_name"),
                "chrom": chrom,
                "strand": strand,
                "transcript_type": first_attr(attrs, "transcript_type"),
                "protein_id": first_attr(attrs, "protein_id"),
                "tags": first_attr(attrs, "tag"),
            }
            if feature == "transcript":
                old = transcripts.get(tx)
                if old is not None and old != record:
                    raise RuntimeError(f"duplicate/conflicting selected transcript record: {tx}")
                transcripts[tx] = record
            elif feature == "CDS":
                try:
                    phase = int(frame)
                except ValueError as exc:
                    raise RuntimeError(f"selected CDS has non-numeric phase at GTF line {line_number}") from exc
                if phase not in (0, 1, 2):
                    raise RuntimeError(f"selected CDS phase outside 0/1/2 at GTF line {line_number}")
                cds_rows[tx].append({
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "phase": phase,
                    "protein_id": first_attr(attrs, "protein_id"),
                })
            elif feature == "stop_codon":
                stop_rows[tx].append({
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                })
    if set(transcripts) != wanted_tx:
        raise RuntimeError(f"GTF did not yield exactly one transcript record for each selected tx: missing {sorted(wanted_tx - set(transcripts))}")
    for tx in sorted(wanted_tx):
        rows = cds_rows.get(tx, [])
        if not rows:
            raise RuntimeError(f"selected transcript has no CDS rows in GTF: {tx}")
        protein_ids = {r["protein_id"] for r in rows if r["protein_id"]}
        tx_protein = transcripts[tx].get("protein_id")
        if tx_protein:
            protein_ids.add(tx_protein)
        if len(protein_ids) != 1:
            raise RuntimeError(f"selected transcript does not map uniquely to one protein in GTF: {tx} -> {sorted(protein_ids)}")
        transcripts[tx]["protein_id"] = next(iter(protein_ids))
        if any(r["chrom"] != transcripts[tx]["chrom"] or r["strand"] != transcripts[tx]["strand"] for r in rows):
            raise RuntimeError(f"selected CDS chromosome/strand mismatch for {tx}")
        stops = stop_rows.get(tx, [])
        if sum(int(r["end"]) - int(r["start"]) + 1 for r in stops) != 3:
            raise RuntimeError(f"selected transcript does not have exactly three stop_codon bases: {tx}")
        if any(r["chrom"] != transcripts[tx]["chrom"] or r["strand"] != transcripts[tx]["strand"] for r in stops):
            raise RuntimeError(f"selected stop_codon chromosome/strand mismatch for {tx}")
    return transcripts, cds_rows, stop_rows


def load_cds_parquet(path: Path, wanted_tx: set[str]) -> tuple[dict[str, list[dict[str, Any]]], str | None]:
    if pq is None:
        raise RuntimeError("pyarrow is required for formal Step10 execution")
    table = pq.read_table(path)
    names = set(table.column_names)
    required = {"transcript_id", "cds_start", "cds_end"}
    if not required.issubset(names):
        raise RuntimeError(f"CDS parquet missing required columns: {sorted(required - names)}")
    phase_col = next((name for name in ("phase", "cds_phase", "frame") if name in names), None)
    frame = table.to_pandas()
    frame["transcript_id"] = frame["transcript_id"].astype(str).map(strip_version)
    frame = frame[frame["transcript_id"].isin(wanted_tx)].copy()
    if frame.empty:
        raise RuntimeError("CDS parquet has no rows for selected transcripts")
    result: dict[str, list[dict[str, Any]]] = {tx: [] for tx in wanted_tx}
    for row in frame.to_dict("records"):
        tx = str(row["transcript_id"])
        item = {
            "transcript_id": tx,
            "start": int(row["cds_start"]),
            "end": int(row["cds_end"]),
        }
        if phase_col is not None:
            value = row[phase_col]
            if pd.isna(value):
                raise RuntimeError(f"CDS parquet phase is NA for {tx}")
            item["phase"] = int(value)
            if item["phase"] not in (0, 1, 2):
                raise RuntimeError(f"CDS parquet phase outside 0/1/2 for {tx}")
        result[tx].append(item)
    if any(not result[tx] for tx in wanted_tx):
        raise RuntimeError(f"CDS parquet missing selected transcript rows: {sorted(tx for tx in wanted_tx if not result[tx])}")
    return result, phase_col


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def oriented_key(start: int, end: int, strand: str) -> tuple[int, int]:
    return (start, end) if strand == "+" else (-end, -start)


def fetch_interval(fasta: pysam.FastaFile, chrom: str, start: int, end: int, strand: str) -> str:
    seq = fasta.fetch(chrom, start - 1, end).upper()
    return seq if strand == "+" else reverse_complement(seq)


def build_sequence(
    intervals: list[tuple[int, int]],
    chrom: str,
    strand: str,
    fasta: pysam.FastaFile,
) -> tuple[str, list[int]]:
    ordered = sorted(intervals, key=lambda x: oriented_key(x[0], x[1], strand))
    sequence_parts: list[str] = []
    coordinates: list[int] = []
    for start, end in ordered:
        if start > end:
            raise RuntimeError(f"invalid interval {chrom}:{start}-{end}:{strand}")
        seq = fetch_interval(fasta, chrom, start, end, strand)
        expected = end - start + 1
        if len(seq) != expected:
            raise RuntimeError(f"FASTA length mismatch for {chrom}:{start}-{end}:{strand}")
        sequence_parts.append(seq)
        if strand == "+":
            coordinates.extend(range(start, end + 1))
        else:
            coordinates.extend(range(end, start - 1, -1))
    return "".join(sequence_parts), coordinates


def translate_first_stop(sequence: str) -> tuple[str, int | None]:
    aas: list[str] = []
    stop_offset: int | None = None
    # Translate complete codons only.  A reconstructed altered transcript can
    # end with one or two trailing nucleotides after a frameshift; those bases
    # are not a codon and do not invalidate an earlier, fully observed PTC.
    # Reference CDS completeness is enforced separately by its caller.
    for offset in range(0, len(sequence) - 2, 3):
        codon = sequence[offset:offset + 3]
        aa = CODON_TABLE.get(codon, "X")
        aas.append(aa)
        if aa == "*":
            stop_offset = offset
            break
    return "".join(aas), stop_offset


def parse_chain(chain: str) -> tuple[str, list[tuple[int, int]]]:
    intervals: list[tuple[int, int]] = []
    chrom = ""
    strand = ""
    for token in str(chain).split("|"):
        match = re.fullmatch(r"(chr[^:]+):(\d+)-(\d+):([+-])", token)
        if not match:
            raise RuntimeError(f"invalid altered_exon_chain token: {token!r}")
        this_chrom, start, end, this_strand = match.groups()
        if not chrom:
            chrom, strand = this_chrom, this_strand
        if (this_chrom, this_strand) != (chrom, strand):
            raise RuntimeError("altered chain mixes chromosomes or strands")
        intervals.append((int(start), int(end)))
    if not intervals:
        raise RuntimeError("altered chain is empty")
    return f"{chrom}:{strand}", intervals


def coordinate_set(intervals: list[tuple[int, int]]) -> set[int]:
    result: set[int] = set()
    for start, end in intervals:
        result.update(range(start, end + 1))
    return result


def validate_phase(
    tx: str,
    gtf_rows: list[dict[str, Any]],
    parquet_rows: list[dict[str, Any]],
    parquet_phase_col: str | None,
) -> dict[str, Any]:
    gtf_key = {(int(r["start"]), int(r["end"])): int(r["phase"]) for r in gtf_rows}
    pq_key = {(int(r["start"]), int(r["end"])): r for r in parquet_rows}
    if set(gtf_key) != set(pq_key):
        raise RuntimeError(f"GTF/CDS parquet coordinate mismatch for {tx}")
    mismatches: list[dict[str, Any]] = []
    if parquet_phase_col is not None:
        for key, gtf_phase in gtf_key.items():
            if int(pq_key[key]["phase"]) != gtf_phase:
                mismatches.append({"interval": key, "gtf_phase": gtf_phase, "parquet_phase": pq_key[key]["phase"]})
    if mismatches:
        raise RuntimeError(f"GTF/CDS parquet phase mismatch for {tx}: {mismatches[:3]}")
    return {
        "tx": tx,
        "segments": len(gtf_rows),
        "phase_values_valid": True,
        "phase_values_compared": parquet_phase_col is not None,
        "phase_only_validation_no_trimming": True,
    }


def first_nonretained(reference_aa: str, altered_aa: str) -> tuple[int, int, str]:
    n = 0
    for ref, alt in zip(reference_aa, altered_aa):
        if ref != alt:
            break
        n += 1
    if n == len(altered_aa) and len(altered_aa) <= len(reference_aa):
        return n, n + 1 if n < len(reference_aa) else -1, "PREFIX_MATCH"
    return n, n + 1, "PREFIX_MISMATCH"


def reconstruction_rows(path: Path, selected: dict[tuple[str, str, str], dict[str, str]]) -> dict[tuple[str, str, str], dict[str, str]]:
    found: dict[tuple[str, str, str], dict[str, str]] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "variant_uid", "gene", "junction_id", "transcript_full_id", "altered_exon_chain",
            "altered_splice_type", "altered_cdna_length", "reference_protein_length_aa",
            "altered_protein_length_aa", "reference_stop_cds_offset", "altered_stop_cds_offset",
            "frame_consequence", "nmd_status",
        }
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise RuntimeError(f"reconstruction missing columns: {sorted(required - set(reader.fieldnames or []))}")
        for row in reader:
            key = (row["variant_uid"], row["gene"], row["junction_id"])
            wanted = selected.get(key)
            if wanted is None or row["transcript_full_id"] != wanted["transcript_full_id"]:
                continue
            if key in found:
                raise RuntimeError(f"multiple reconstruction rows for frozen event/transcript: {key}")
            found[key] = row
    missing = sorted(set(selected) - set(found))
    if missing:
        raise RuntimeError(f"reconstruction missing frozen event rows: {missing}")
    return found


def row_result(
    key: dict[str, str],
    phase0: dict[str, str],
    recon: dict[str, str],
    tx_meta: dict[str, Any],
    gtf_cds: list[dict[str, Any]],
    gtf_stops: list[dict[str, Any]],
    pq_cds: list[dict[str, Any]],
    fasta: pysam.FastaFile,
) -> dict[str, Any]:
    tx = tx_meta["tx"]
    chrom = tx_meta["chrom"]
    strand = tx_meta["strand"]
    ref_intervals = [(int(r["start"]), int(r["end"])) for r in pq_cds]
    ref_cds, ref_cds_coords = build_sequence(ref_intervals, chrom, strand, fasta)
    if len(ref_cds) % 3:
        raise RuntimeError(f"reference CDS % 3 != 0 for {tx}")
    if "N" in ref_cds:
        raise RuntimeError(f"reference CDS contains N for {tx}")
    if not ref_cds.startswith("ATG"):
        raise RuntimeError(f"reference CDS does not start with canonical ATG for {tx}: {ref_cds[:12]}")
    # GENCODE CDS features exclude the separately annotated stop_codon.
    # Therefore a complete reference CDS must have no in-frame stop; the
    # terminal codon is validated from the three stop_codon feature bases.
    ref_aa, internal_ref_stop = translate_first_stop(ref_cds)
    if internal_ref_stop is not None:
        raise RuntimeError(f"reference CDS has internal stop for {tx}: first stop {internal_ref_stop}, length {len(ref_cds)}")
    if any(aa == "X" for aa in ref_aa):
        raise RuntimeError(f"reference CDS translated unknown codon for {tx}")
    stop_intervals = [(int(r["start"]), int(r["end"])) for r in gtf_stops]
    reference_stop_codon, _reference_stop_coords = build_sequence(stop_intervals, chrom, strand, fasta)
    if reference_stop_codon not in {"TAA", "TAG", "TGA"}:
        raise RuntimeError(f"reference stop_codon feature is not a stop for {tx}: {reference_stop_codon}")
    ref_stop = len(ref_cds)

    chain_label, altered_intervals = parse_chain(recon["altered_exon_chain"])
    altered_full, altered_coords = build_sequence(altered_intervals, chrom, strand, fasta)
    if recon["altered_cdna_length"] and int(float(recon["altered_cdna_length"])) != len(altered_full):
        raise RuntimeError(f"altered cDNA length mismatch for {key['variant_uid']}")
    altered_coord_set = set(altered_coords)
    if ref_cds_coords[0] not in altered_coord_set:
        raise RuntimeError(f"canonical CDS/AUG anchor absent from altered chain for {key['variant_uid']}")
    altered_cds_coords = [coord for coord in ref_cds_coords if coord in altered_coord_set]
    if not altered_cds_coords:
        raise RuntimeError(f"altered chain has no retained CDS for {key['variant_uid']}")
    coordinate_to_base = dict(zip(altered_coords, altered_full))
    altered_cds = "".join(coordinate_to_base[coord] for coord in altered_cds_coords)
    if not altered_cds.startswith("ATG"):
        raise RuntimeError(f"altered CDS does not retain canonical ATG for {key['variant_uid']}")
    if "N" in altered_cds:
        raise RuntimeError(f"altered CDS contains N for {key['variant_uid']}")
    altered_aa_with_stop, altered_stop = translate_first_stop(altered_cds)
    if altered_stop is None:
        raise RuntimeError(f"altered CDS has no in-frame stop for {key['variant_uid']}")
    if any(aa == "X" for aa in altered_aa_with_stop):
        raise RuntimeError(f"altered CDS translated unknown codon for {key['variant_uid']}")
    altered_aa = altered_aa_with_stop[:-1]
    recon_alt_stop = int(float(recon["altered_stop_cds_offset"]))
    recon_ref_stop = int(float(recon["reference_stop_cds_offset"]))
    if recon_alt_stop != altered_stop:
        raise RuntimeError(f"altered stop offset mismatch for {key['variant_uid']}: reconstruction={recon_alt_stop}, rebuilt={altered_stop}")
    if recon_ref_stop != ref_stop:
        raise RuntimeError(f"reference stop offset mismatch for {key['variant_uid']}: reconstruction={recon_ref_stop}, rebuilt={ref_stop}")
    expected_ref, expected_alt = EXPECTED_LENGTHS[key["variant_uid"]]
    if len(ref_aa) != expected_ref:
        raise RuntimeError(f"reference protein length mismatch for {key['variant_uid']}: {len(ref_aa)} != {expected_ref}")
    expected_alt_from_stop = altered_stop // 3
    if expected_alt_from_stop != expected_alt:
        raise RuntimeError(f"altered protein length mismatch for {key['variant_uid']}: {expected_alt_from_stop} != {expected_alt}")
    if int(float(recon["altered_protein_length_aa"])) != expected_alt:
        raise RuntimeError(f"reconstruction altered length mismatch for {key['variant_uid']}")
    if recon["frame_consequence"] != "FRAMESHIFT_PTC" or recon["nmd_status"] != "NMD_COMPATIBLE":
        raise RuntimeError(f"reconstruction frame/NMD mismatch for {key['variant_uid']}: {recon['frame_consequence']}/{recon['nmd_status']}")
    if phase0.get("v23_frame_consequence", "") != "FRAMESHIFT_PTC" or phase0.get("v23_nmd_status", "") != "NMD_COMPATIBLE":
        raise RuntimeError(f"phase0 frame/NMD mismatch for {key['variant_uid']}")
    retained, first_lost, prefix_status = first_nonretained(ref_aa, altered_aa)
    # The frozen event table describes a truncating first in-frame stop.  Do
    # not call the result "protein remaining"; this is only a length ratio.
    nmd_distance = None
    last_junction = None
    if len(altered_coords) > 0:
        # Junction boundaries are inferred from discontinuities in the
        # coordinate-only altered chain.  For +/- strands, a transcript-order
        # coordinate discontinuity is an exon-exon junction.
        for prev, curr in zip(altered_coords[:-1], altered_coords[1:]):
            if (strand == "+" and curr != prev + 1) or (strand == "-" and curr != prev - 1):
                last_junction = altered_coords.index(curr)
        if last_junction is not None:
            nmd_distance = last_junction - altered_stop - 3
    nmd_geometry_match = nmd_distance is not None and nmd_distance > 50
    if not nmd_geometry_match:
        raise RuntimeError(f"rebuilt NMD geometry does not satisfy >50 nt for {key['variant_uid']}: {nmd_distance}")
    return {
        "histology": key.get("histology", ""),
        "variant_uid": key["variant_uid"],
        "gene": key["gene"],
        "junction_id": key["junction_id"],
        "phase0_R_class": phase0.get("phase0_R_class", ""),
        "phase0_event_class": phase0.get("phase0_event_class", ""),
        "transcript_id": tx,
        "transcript_full_id": tx_meta["full_id"],
        "protein_full_id": tx_meta["protein_id"],
        "reference_chrom": chrom,
        "reference_strand": strand,
        "reference_cds_length_nt": len(ref_cds),
        "reference_protein_length_aa": len(ref_aa),
        "reference_cds_sha256": hashlib.sha256(ref_cds.encode()).hexdigest(),
        "reference_aa_sha256": hashlib.sha256(ref_aa.encode()).hexdigest(),
        "altered_exon_chain": recon["altered_exon_chain"],
        "altered_chain_coordinate_label": chain_label,
        "altered_cds_length_nt": len(altered_cds),
        "altered_aa_length_aa": len(altered_aa),
        "altered_cds_sha256": hashlib.sha256(altered_cds.encode()).hexdigest(),
        "altered_aa_sha256": hashlib.sha256(altered_aa.encode()).hexdigest(),
        "altered_splice_type": recon["altered_splice_type"],
        "frame_consequence_rebuilt": "FRAMESHIFT_PTC",
        "frame_consequence_reconstruction": recon["frame_consequence"],
        "ptc_status": "FIRST_IN_FRAME_NONCANONICAL_EARLY_STOP",
        "ptc_position_aa_1based": altered_stop // 3 + 1,
        "altered_stop_offset_nt_0based": altered_stop,
        "altered_stop_codon": altered_cds[altered_stop:altered_stop + 3],
        "reference_stop_offset_nt_0based": ref_stop,
        "reference_stop_codon": reference_stop_codon,
        "predicted_length_fraction": len(altered_aa) / len(ref_aa),
        "exact_reference_residue_count_retained": retained,
        "exact_reference_residue_fraction_retained": retained / len(ref_aa),
        "first_nonretained_reference_aa": first_lost,
        "reference_prefix_match_status": prefix_status,
        "nmd_rule": "STOP_TO_LAST_JUNCTION_GT_50NT",
        "nmd_stop_to_last_junction_distance_nt": nmd_distance,
        "nmd_geometry_status": "NMD_COMPATIBLE",
        "nmd_phase0_match": True,
        "nmd_reconstruction_match": recon["nmd_status"] == "NMD_COMPATIBLE",
        "domain_status": DOMAIN_REASON,
        "n_domain": "NA",
        "domain_reason": DOMAIN_REASON,
        "consequence_scope": CONSEQUENCE_SCOPE,
        "protein_expression_status": PROTEIN_EXPRESSION_STATUS,
        "cis_status": CIS_STATUS,
        "row_status": "PASS",
    }


OUTPUT_COLUMNS = [
    "histology", "variant_uid", "gene", "junction_id", "phase0_R_class", "phase0_event_class",
    "transcript_id", "transcript_full_id", "protein_full_id", "reference_chrom", "reference_strand",
    "reference_cds_length_nt", "reference_protein_length_aa", "reference_cds_sha256", "reference_aa_sha256",
    "altered_exon_chain", "altered_chain_coordinate_label", "altered_splice_type", "altered_cds_length_nt",
    "altered_aa_length_aa", "altered_cds_sha256", "altered_aa_sha256", "frame_consequence_rebuilt",
    "frame_consequence_reconstruction", "ptc_status", "ptc_position_aa_1based", "altered_stop_offset_nt_0based",
    "altered_stop_codon", "reference_stop_offset_nt_0based", "reference_stop_codon", "predicted_length_fraction",
    "exact_reference_residue_count_retained", "exact_reference_residue_fraction_retained",
    "first_nonretained_reference_aa", "reference_prefix_match_status", "nmd_rule",
    "nmd_stop_to_last_junction_distance_nt", "nmd_geometry_status", "nmd_phase0_match", "nmd_reconstruction_match",
    "domain_status", "n_domain", "domain_reason", "consequence_scope", "protein_expression_status",
    "cis_status", "row_status",
]


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def self_test() -> int:
    # Plus/minus sequence assembly and reverse-complement behavior.
    class FakeFasta:
        def __init__(self, seqs: dict[str, str]):
            self.seqs = seqs

        def fetch(self, chrom: str, start: int, end: int) -> str:
            return self.seqs[chrom][start:end]

    plus = FakeFasta({"chrT": "AAATGCCCTAA"})
    seq, coords = build_sequence([(1, 5), (9, 11)], "chrT", "+", plus)  # type: ignore[arg-type]
    assert seq == "AAATGTAA" and coords == [1, 2, 3, 4, 5, 9, 10, 11]
    minus = FakeFasta({"chrM": "TTATTTCAT"})
    seq, coords = build_sequence([(1, 3), (5, 9)], "chrM", "-", minus)  # type: ignore[arg-type]
    assert seq == "ATGAATAA" and coords == [9, 8, 7, 6, 5, 3, 2, 1]
    # Stop offset and NMD >50 boundary are explicit to prevent off-by-one.
    aa, offset = translate_first_stop("ATG" + "GCC" * 40 + "TAA")
    assert offset == 123 and len(aa[:-1]) == 41
    # GENCODE represents the reference terminal stop separately from CDS.
    aa_no_stop, offset_no_stop = translate_first_stop("ATGGCC")
    assert aa_no_stop == "MA" and offset_no_stop is None
    aa_trailing, offset_trailing = translate_first_stop("ATGTAACC")
    assert aa_trailing == "M*" and offset_trailing == 3
    assert (99 - 3 - 45) == 51 and (99 - 3 - 46) == 50
    # Exact frozen key and length-oracle contract is separate from input files.
    expected_length_oracle = {
        "chr12:25227342:T:A": (188, 40),
        "chr19:10489648:C:A": (624, 456),
        "chr19:1220373:G:T": (433, 155),
        "chr19:45420410:C:T": (297, 151),
        "chr6:43443117:G:T": (1492, 822),
        "chr13:48307406:G:T": (928, 67),
        "chr17:7673570:T:C": (393, 310),
        "chr17:7676379:T:A": (393, 35),
        "chr2:224481892:C:A": (768, 632),
    }
    assert EXPECTED_LENGTHS == expected_length_oracle
    assert len(EXPECTED_LENGTHS) == EXPECTED_N
    print(json.dumps({"status": "SELF_TEST_PASS", "tests": ["plus_minus_cds_assembly", "phase_no_trim_contract", "gencode_separate_stop_codon_contract", "altered_trailing_partial_codon_contract", "ptc_stop_offset", "nmd_50nt_boundary", "frozen9_length_oracle"]}, ensure_ascii=False))
    return 0


def run(args: argparse.Namespace) -> int:
    if pq is None or pysam is None:
        raise RuntimeError("formal Step10 execution requires pyarrow and pysam in the server environment")
    for path in (args.frozen_keys, args.phase0_r2, args.reconstruction, args.gtf, args.cds, args.fasta, args.fai):
        if not path.is_file():
            raise FileNotFoundError(path)
    frozen = parse_frozen_keys(args.frozen_keys)
    phase_rows = read_tsv(args.phase0_r2)
    require_columns(phase_rows, ["variant_uid", "gene", "junction_id", "v23_frame_consequence", "v23_nmd_status"], "Phase0 R2")
    phase_lookup = {(r["variant_uid"], r["gene"], r["junction_id"]): r for r in phase_rows}
    selected = {(r["variant_uid"], r["gene"], r["junction_id"]): r for r in frozen}
    if set(selected) != set(phase_lookup) & set(selected):
        raise RuntimeError("Phase0 R2 does not cover all frozen9 keys")
    wanted_tx = {strip_version(r["selected_transcript"]) for r in frozen}
    tx_meta, gtf_cds, gtf_stops = parse_gtf(args.gtf, wanted_tx)
    for key, row in selected.items():
        tx = strip_version(row["selected_transcript"])
        if tx_meta[tx]["gene"] != row["gene"]:
            raise RuntimeError(f"GTF gene mismatch for {tx}: {tx_meta[tx]['gene']} != {row['gene']}")
        row["transcript_full_id"] = tx_meta[tx]["full_id"]
    pq_cds, phase_col = load_cds_parquet(args.cds, wanted_tx)
    phase_audit = [validate_phase(tx, gtf_cds[tx], pq_cds[tx], phase_col) for tx in sorted(wanted_tx)]
    recon = reconstruction_rows(args.reconstruction, selected)
    input_hashes = {
        "frozen_keys": metadata(args.frozen_keys),
        "phase0_r2": metadata(args.phase0_r2),
        "reconstruction": metadata(args.reconstruction, EXPECTED_RECON_SHA256),
        "gencode50_gtf": metadata(args.gtf, EXPECTED_GTF_SHA256),
        "gencode50_cds_parquet": metadata(args.cds, EXPECTED_CDS_SHA256),
        "grch38_fasta": metadata(args.fasta, EXPECTED_FASTA_SHA256),
        "grch38_fai": metadata(args.fai, EXPECTED_FAI_SHA256),
    }
    for label, info in input_hashes.items():
        if info["expected_sha256"] and not info["sha256_match"]:
            raise RuntimeError(f"frozen input SHA256 mismatch: {label}: {info['sha256']} != {info['expected_sha256']}")
    rows: list[dict[str, Any]] = []
    with pysam.FastaFile(str(args.fasta)) as fasta:
        for key in sorted(selected, key=lambda item: item[0]):
            frozen_row = selected[key]
            tx = strip_version(frozen_row["selected_transcript"])
            rows.append(row_result(frozen_row, phase_lookup[key], recon[key], tx_meta[tx], gtf_cds[tx], gtf_stops[tx], pq_cds[tx], fasta))
    if len(rows) != EXPECTED_N or {r["variant_uid"] for r in rows} != set(EXPECTED_LENGTHS):
        raise RuntimeError("Step10 output is not the exact frozen9 key set")
    if any(r["row_status"] != "PASS" for r in rows):
        raise RuntimeError("one or more Step10 rows did not pass")
    out = args.outdir / "V23_R2_PROTEIN_DOMAIN_REFINEMENT.tsv"
    write_tsv(out, rows)
    audit = {
        "status": "PASS",
        "module": "STEP10_R2_PROTEIN_PTC_NMD_DOMAIN_REFINEMENT",
        "event_count": len(rows),
        "event_set_exact_frozen9": True,
        "scope": "R2_ONLY; R3_NOT_PROCESSED",
        "command": " ".join(sys.argv),
        "inputs": input_hashes,
        "outcome_fields_used": [
            "frozen9 identity", "Phase0 v23_frame_consequence", "Phase0 v23_nmd_status",
            "full reconstruction altered_exon_chain/stop/frame/NMD fields",
        ],
        "outcome_fields_not_used_as_sequence_source": [
            "Phase0 protein lengths", "Phase0 protein fraction", "old domain placeholder values",
        ],
        "transcript_protein_mapping": [
            {"transcript_full_id": tx_meta[tx]["full_id"], "protein_full_id": tx_meta[tx]["protein_id"], "unique_from_original_gtf": True}
            for tx in sorted(wanted_tx)
        ],
        "phase_validation": phase_audit,
        "phase_rule": "validate GTF/parquet phase values only; never trim CDS bases by phase",
        "reference_reconstruction": {
            "source": "GENCODE50_CDS.parquet + GRCh38 FASTA",
            "requirements": ["CDS%3=0", "canonical_ATG", "no_N", "no_internal_stop"],
            "passed": True,
        },
        "altered_reconstruction": {
            "source": "full reconstruction altered_exon_chain + GRCh38 FASTA",
            "somatic_allele_applied": False,
            "canonical_AUG_anchored": True,
            "first_in_frame_noncanonical_stop_revalidated": True,
        },
        "nmd_rule": "altered stop to final altered exon-exon junction > 50 nt",
        "nmd_scope": "predicted NMD-compatible consequence; no degradation/protein-expression claim",
        "domain": {
            "status": DOMAIN_REASON,
            "n_domain": "NA",
            "asset_used": False,
            "future_interface": "exact protein-sequence-matched domain coordinate map may be supplied in a later module",
        },
        "output": metadata(out),
        "outcome_fields": {
            "consequence_scope": CONSEQUENCE_SCOPE,
            "protein_expression_status": PROTEIN_EXPRESSION_STATUS,
            "cis_status": CIS_STATUS,
        },
        "rows": rows,
    }
    audit_path = args.audit_dir / "V23_R2_PROTEIN_DOMAIN_REFINEMENT_AUDIT.json"
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    manifest_path = args.audit_dir / "V23_R2_PROTEIN_DOMAIN_REFINEMENT_SHA256.tsv"
    with manifest_path.open("wt", encoding="utf-8", newline="") as handle:
        handle.write("kind\tname\tpath\tbytes\tsha256\n")
        for kind, path in [("input", args.reconstruction), ("input", args.gtf), ("input", args.cds), ("input", args.fasta), ("input", args.fai), ("output", out), ("audit", audit_path)]:
            handle.write(f"{kind}\t{path.name}\t{path}\t{path.stat().st_size}\t{sha256_file(path)}\n")
    print(json.dumps({
        "status": "PASS", "output": str(out), "audit": str(audit_path), "manifest": str(manifest_path),
        "rows": len(rows), "ptc_lengths": {r["variant_uid"]: {"ptc_aa": r["ptc_position_aa_1based"], "altered_aa": r["altered_aa_length_aa"], "fraction": r["predicted_length_fraction"]} for r in rows},
    }, ensure_ascii=False, indent=2))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--frozen-keys", type=Path)
    parser.add_argument("--phase0-r2", type=Path)
    parser.add_argument("--reconstruction", type=Path)
    parser.add_argument("--gtf", type=Path)
    parser.add_argument("--cds", type=Path)
    parser.add_argument("--fasta", type=Path)
    parser.add_argument("--fai", type=Path)
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--audit-dir", type=Path)
    args = parser.parse_args()
    if args.self_test:
        return self_test()
    required = ["frozen_keys", "phase0_r2", "reconstruction", "gtf", "cds", "fasta", "fai", "outdir", "audit_dir"]
    missing = [name for name in required if getattr(args, name) is None]
    if missing:
        parser.error("missing required arguments: " + ", ".join("--" + name.replace("_", "-") for name in missing))
    try:
        return run(args)
    except Exception as exc:
        print(json.dumps({"status": "HARD_STOP", "error": str(exc)}, ensure_ascii=False), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
