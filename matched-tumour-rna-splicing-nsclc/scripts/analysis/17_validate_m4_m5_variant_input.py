#!/usr/bin/env python3
"""Outcome-blind technical validation of the frozen M4/M5 variant input.

Only the frozen VCF, frozen variant-universe TSV, frozen input JSON, and the
reference FASTA/.fai are consumed.  No association or carrier outcome is read.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import tempfile
from pathlib import Path
from typing import Any


EXPECTED_N = 2374
EXPECTED_VCF_SHA256 = "0274386b83094adb5b9e5796568aa29b7fe6b124eb516f633d475072b34a7f99"
EXPECTED_TSV_SHA256 = "55d73120ffaa7fc5f7e14e5673a47edb8e5433a854143bbc0943f8d0b4a8583e"
EXPECTED_FREEZE_JSON_SHA256 = "d73f0aef34bc1014ea92c4c0ae84da640ffbf9b1c1bff8d6cadd126ec33b2668"
OUTPUT_TSV = "V23_M4_M5_VARIANT_REF_VALIDATION.tsv"
OUTPUT_AUDIT = "V23_M4_M5_VARIANT_REF_VALIDATION_AUDIT.json"
ALLELE_RE = re.compile(r"^[ACGT]+$")
CHROM_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
EXPECTED_TSV_COLUMNS = {
    "predictor_variant_id", "variant_uid", "chrom", "pos_1based", "ref", "alt",
    "outcome_blind_input",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metadata(path: Path, sha256: str | None = None) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256 or sha256_file(path),
    }


def expected_id(index: int) -> str:
    return f"V23VAR{index:05d}"


def parse_vcf(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    header_seen = False
    with path.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            text = line.rstrip("\r\n")
            if not text or text.startswith("##"):
                continue
            if text.startswith("#CHROM\t"):
                if header_seen:
                    raise ValueError("VCF has more than one #CHROM header")
                columns = text.split("\t")
                if columns != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]:
                    raise ValueError("VCF must contain exactly the eight frozen site-only columns")
                header_seen = True
                continue
            if text.startswith("#"):
                raise ValueError(f"unexpected VCF header line {line_number}")
            if not header_seen:
                raise ValueError("VCF data appeared before the #CHROM header")
            fields = text.split("\t")
            if len(fields) != 8:
                raise ValueError(f"VCF line {line_number} does not have exactly eight columns")
            chrom, pos_text, variant_id, ref, alt, _qual, filt, _info = fields
            try:
                pos = int(pos_text)
            except ValueError as exc:
                raise ValueError(f"VCF line {line_number} has non-integer POS") from exc
            record_number = len(records) + 1
            if variant_id != expected_id(record_number):
                raise ValueError(
                    f"VCF ID sequence failure at record {record_number}: {variant_id!r}"
                )
            if not CHROM_RE.fullmatch(chrom):
                raise ValueError(f"VCF line {line_number} has invalid CHROM {chrom!r}")
            if pos < 1:
                raise ValueError(f"VCF line {line_number} has POS < 1")
            if not ALLELE_RE.fullmatch(ref) or not ALLELE_RE.fullmatch(alt):
                raise ValueError(
                    f"VCF line {line_number} is not an uppercase simple A/C/G/T biallelic record"
                )
            if ref == alt:
                raise ValueError(f"VCF line {line_number} has REF == ALT")
            if filt != "PASS":
                raise ValueError(f"VCF line {line_number} FILTER is not PASS")
            records.append({
                "predictor_variant_id": variant_id,
                "variant_uid": f"{chrom}:{pos}:{ref}:{alt}",
                "chrom": chrom,
                "pos_1based": pos,
                "ref": ref,
                "alt": alt,
            })
    if not header_seen:
        raise ValueError("VCF is missing the #CHROM header")
    if len(records) != EXPECTED_N:
        raise ValueError(f"VCF record count is {len(records)}, expected {EXPECTED_N}")
    ids = [row["predictor_variant_id"] for row in records]
    variants = [(row["chrom"], row["pos_1based"], row["ref"], row["alt"]) for row in records]
    if len(ids) != len(set(ids)):
        raise ValueError("VCF IDs are not unique")
    if len(variants) != len(set(variants)):
        raise ValueError("VCF variant records are not unique")
    return records


def parse_variant_universe(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("variant-universe TSV has no header")
        missing = EXPECTED_TSV_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(f"variant-universe TSV missing columns: {sorted(missing)}")
        rows = [dict(row) for row in reader]
    if len(rows) != EXPECTED_N:
        raise ValueError(f"variant-universe row count is {len(rows)}, expected {EXPECTED_N}")
    ids: list[str] = []
    uids: list[str] = []
    for index, row in enumerate(rows, 1):
        variant_id = row["predictor_variant_id"]
        if variant_id != expected_id(index):
            raise ValueError(f"TSV ID sequence failure at row {index}: {variant_id!r}")
        if row["outcome_blind_input"] != "TRUE":
            raise ValueError(f"TSV row {index} outcome_blind_input is not exactly TRUE")
        try:
            pos = int(row["pos_1based"])
        except ValueError as exc:
            raise ValueError(f"TSV row {index} has non-integer pos_1based") from exc
        expected_uid = f"{row['chrom']}:{pos}:{row['ref']}:{row['alt']}"
        if row["variant_uid"] != expected_uid:
            raise ValueError(f"TSV row {index} variant_uid is inconsistent with its allele fields")
        ids.append(variant_id)
        uids.append(row["variant_uid"])
    if len(ids) != len(set(ids)) or len(uids) != len(set(uids)):
        raise ValueError("variant-universe TSV IDs or variant_uid values are not unique")
    return rows


def parse_fai(path: Path) -> dict[str, tuple[int, int, int, int]]:
    index: dict[str, tuple[int, int, int, int]] = {}
    with path.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            text = line.rstrip("\r\n")
            if not text:
                continue
            fields = text.split("\t")
            if len(fields) < 5:
                raise ValueError(f"FASTA .fai line {line_number} has fewer than five columns")
            chrom = fields[0]
            if chrom in index:
                raise ValueError(f"duplicate contig in FASTA .fai: {chrom}")
            try:
                length, offset, line_bases, line_width = map(int, fields[1:5])
            except ValueError as exc:
                raise ValueError(f"non-integer FASTA .fai fields on line {line_number}") from exc
            if length < 1 or offset < 0 or line_bases < 1 or line_width < line_bases:
                raise ValueError(f"invalid FASTA .fai geometry on line {line_number}")
            index[chrom] = (length, offset, line_bases, line_width)
    if not index:
        raise ValueError("FASTA .fai is empty")
    return index


def fetch_reference(
    fasta_handle: Any,
    fai_entry: tuple[int, int, int, int],
    pos_1based: int,
    length: int,
) -> str:
    _contig_length, offset, line_bases, line_width = fai_entry
    zero_based = pos_1based - 1
    byte_offset = offset + (zero_based // line_bases) * line_width + (zero_based % line_bases)
    fasta_handle.seek(byte_offset)
    sequence = bytearray()
    while len(sequence) < length:
        block = fasta_handle.read(max(64, length - len(sequence) + 4))
        if not block:
            break
        sequence.extend(base for base in block if base not in (10, 13))
    try:
        return bytes(sequence[:length]).decode("ascii")
    except UnicodeDecodeError as exc:
        raise ValueError("reference FASTA contains non-ASCII sequence bytes") from exc


def safe_commit(temp_path: Path, target_path: Path) -> str:
    new_hash = sha256_file(temp_path)
    if target_path.exists():
        old_hash = sha256_file(target_path)
        if old_hash != new_hash:
            raise FileExistsError(
                f"refusing to replace existing target with different hash: {target_path}; "
                f"existing={old_hash}, new={new_hash}"
            )
        temp_path.unlink()
        return old_hash
    os.replace(temp_path, target_path)
    return new_hash


def write_validation_table(path: Path, rows: list[dict[str, Any]]) -> None:
    columns = [
        "predictor_variant_id", "variant_uid", "chrom", "pos_1based", "ref", "alt",
        "fasta_ref", "vcf_pass", "tsv_identity_match", "contig_present",
        "coordinate_in_bounds", "ref_match_case_sensitive", "validation_status",
    ]
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def validate(args: argparse.Namespace) -> dict[str, Any]:
    paths = [args.vcf, args.variant_universe, args.fasta, args.fasta_fai, args.freeze_json]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)

    vcf_hash = sha256_file(args.vcf)
    tsv_hash = sha256_file(args.variant_universe)
    freeze_hash = sha256_file(args.freeze_json)
    if vcf_hash != EXPECTED_VCF_SHA256:
        raise ValueError(f"frozen VCF SHA256 mismatch: {vcf_hash}")
    if tsv_hash != EXPECTED_TSV_SHA256:
        raise ValueError(f"frozen variant-universe TSV SHA256 mismatch: {tsv_hash}")
    if freeze_hash != EXPECTED_FREEZE_JSON_SHA256:
        raise ValueError(f"frozen input JSON SHA256 mismatch: {freeze_hash}")

    freeze = json.loads(args.freeze_json.read_text(encoding="utf-8"))
    if freeze.get("unique_variants") != EXPECTED_N:
        raise ValueError("freeze JSON unique_variants is not 2374")
    expected_fasta_hash = freeze.get("reference_fasta_sha256_expected")
    if not isinstance(expected_fasta_hash, str) or not re.fullmatch(r"[0-9a-f]{64}", expected_fasta_hash):
        raise ValueError("freeze JSON has no valid reference_fasta_sha256_expected")
    fasta_hash = sha256_file(args.fasta)
    if fasta_hash != expected_fasta_hash:
        raise ValueError(f"reference FASTA SHA256 mismatch: {fasta_hash}")
    fai_hash = sha256_file(args.fasta_fai)

    vcf_rows = parse_vcf(args.vcf)
    tsv_rows = parse_variant_universe(args.variant_universe)
    fai = parse_fai(args.fasta_fai)
    validation_rows: list[dict[str, Any]] = []
    with args.fasta.open("rb") as fasta_handle:
        for index, (vcf_row, tsv_row) in enumerate(zip(vcf_rows, tsv_rows, strict=True), 1):
            identity = {
                "predictor_variant_id": tsv_row["predictor_variant_id"],
                "variant_uid": tsv_row["variant_uid"],
                "chrom": tsv_row["chrom"],
                "pos_1based": int(tsv_row["pos_1based"]),
                "ref": tsv_row["ref"],
                "alt": tsv_row["alt"],
            }
            if identity != vcf_row:
                raise ValueError(f"VCF/TSV identity mismatch at record {index}")
            chrom = vcf_row["chrom"]
            if chrom not in fai:
                raise ValueError(f"variant {vcf_row['predictor_variant_id']} contig absent from .fai: {chrom}")
            contig_length = fai[chrom][0]
            end = vcf_row["pos_1based"] + len(vcf_row["ref"]) - 1
            if end > contig_length:
                raise ValueError(f"variant {vcf_row['predictor_variant_id']} exceeds contig bounds")
            fasta_ref = fetch_reference(
                fasta_handle, fai[chrom], vcf_row["pos_1based"], len(vcf_row["ref"])
            )
            if fasta_ref != vcf_row["ref"]:
                raise ValueError(
                    f"FASTA REF mismatch for {vcf_row['predictor_variant_id']}: "
                    f"VCF={vcf_row['ref']}, FASTA={fasta_ref}"
                )
            validation_rows.append({
                **vcf_row,
                "fasta_ref": fasta_ref,
                "vcf_pass": "TRUE",
                "tsv_identity_match": "TRUE",
                "contig_present": "TRUE",
                "coordinate_in_bounds": "TRUE",
                "ref_match_case_sensitive": "TRUE",
                "validation_status": "PASS",
            })

    if len(validation_rows) != EXPECTED_N:
        raise AssertionError("internal validation row-count failure")
    args.audit_dir.mkdir(parents=True, exist_ok=True)
    script_path = Path(__file__).resolve()
    output_path = args.audit_dir / OUTPUT_TSV
    audit_path = args.audit_dir / OUTPUT_AUDIT
    with tempfile.NamedTemporaryFile(
        mode="w", delete=False, dir=args.audit_dir, prefix=OUTPUT_TSV + ".", suffix=".tmp"
    ) as handle:
        temp_output = Path(handle.name)
    try:
        write_validation_table(temp_output, validation_rows)
        output_hash = safe_commit(temp_output, output_path)
    finally:
        if temp_output.exists():
            temp_output.unlink()

    audit = {
        "status": "PASS",
        "analysis_step": "M4_M5_OUTCOME_BLIND_VARIANT_REF_VALIDATION",
        "contracts": {
            "expected_variants": EXPECTED_N,
            "variant_ids": "V23VAR00001..V23VAR02374, unique and consecutive",
            "vcf": "site-only, biallelic uppercase simple A/C/G/T alleles, FILTER=PASS",
            "identity": "VCF and frozen TSV ID/UID/chrom/pos/ref/alt one-to-one exact",
            "reference": "case-sensitive exact FASTA REF match on .fai contigs and in bounds",
            "outcome_fields_used": [],
            "forbidden_inputs_read": [],
        },
        "counts": {
            "vcf_records": len(vcf_rows),
            "variant_universe_rows": len(tsv_rows),
            "one_to_one_identity_matches": len(validation_rows),
            "fai_contig_present": len(validation_rows),
            "coordinates_in_bounds": len(validation_rows),
            "fasta_ref_case_sensitive_matches": len(validation_rows),
            "pass_rows": len(validation_rows),
            "fail_rows": 0,
        },
        "inputs": {
            "frozen_vcf": metadata(args.vcf, vcf_hash),
            "frozen_variant_universe": metadata(args.variant_universe, tsv_hash),
            "frozen_input_json": metadata(args.freeze_json, freeze_hash),
            "reference_fasta": metadata(args.fasta, fasta_hash),
            "reference_fasta_fai": metadata(args.fasta_fai, fai_hash),
        },
        "script": metadata(script_path),
        "output": {
            "path": str(output_path.resolve()),
            "size_bytes": output_path.stat().st_size,
            "sha256": output_hash,
        },
    }
    audit_text = json.dumps(audit, indent=2, ensure_ascii=False, sort_keys=True) + "\n"
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", delete=False,
        dir=args.audit_dir, prefix=OUTPUT_AUDIT + ".", suffix=".tmp"
    ) as handle:
        handle.write(audit_text)
        temp_audit = Path(handle.name)
    try:
        safe_commit(temp_audit, audit_path)
    finally:
        if temp_audit.exists():
            temp_audit.unlink()
    return audit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--variant-universe", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--fasta-fai", required=True, type=Path)
    parser.add_argument("--freeze-json", required=True, type=Path)
    parser.add_argument("--audit-dir", required=True, type=Path)
    args = parser.parse_args()
    audit = validate(args)
    print(json.dumps({"status": audit["status"], "counts": audit["counts"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
