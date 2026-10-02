#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


EXPECTED_VARIANTS = 2374


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def split_vcf(path: Path) -> tuple[list[str], list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    chrom_index = next((i for i, line in enumerate(lines) if line.startswith("#CHROM\t")), None)
    if chrom_index is None:
        raise ValueError(f"#CHROM missing: {path}")
    return lines[: chrom_index + 1], lines[chrom_index + 1 :]


def keys(lines: list[str]) -> list[tuple[str, str, str, str, str]]:
    result = []
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            raise ValueError(f"invalid VCF record: {line!r}")
        result.append(tuple(fields[:5]))
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--transport-vcf", required=True, type=Path)
    parser.add_argument("--input-shard-dir", required=True, type=Path)
    parser.add_argument("--prediction-shard-dir", required=True, type=Path)
    parser.add_argument("--output-vcf", required=True, type=Path)
    parser.add_argument("--audit-json", required=True, type=Path)
    parser.add_argument("--predictor", required=True, choices=("SpliceAI", "Pangolin"))
    parser.add_argument("--distance", required=True, type=int, choices=(50, 500))
    parser.add_argument("--shards", type=int, default=8)
    args = parser.parse_args()
    for target in (args.output_vcf, args.audit_json):
        if target.exists():
            raise FileExistsError(f"refusing to overwrite: {target}")

    _, transport_records = split_vcf(args.transport_vcf)
    expected_global_keys = keys(transport_records)
    if len(expected_global_keys) != EXPECTED_VARIANTS:
        raise ValueError("transport variant count mismatch")

    common_header: list[str] | None = None
    merged_records: list[str] = []
    shard_audits: list[dict[str, object]] = []
    annotated = 0
    info_token = f"{args.predictor}="
    for index in range(args.shards):
        name = f"shard_{index:03d}.vcf"
        input_path = args.input_shard_dir / name
        prediction_path = args.prediction_shard_dir / name
        for path in (input_path, prediction_path):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)
        _, input_records = split_vcf(input_path)
        header, prediction_records = split_vcf(prediction_path)
        if common_header is None:
            common_header = header
        elif header != common_header:
            raise ValueError(f"prediction VCF headers differ at shard {index}")
        if keys(prediction_records) != keys(input_records):
            raise ValueError(f"prediction/input identity or order mismatch at shard {index}")
        annotated += sum(info_token in line.split("\t")[7] for line in prediction_records)
        merged_records.extend(prediction_records)
        shard_audits.append(
            {
                "index": index,
                "n_variants": len(prediction_records),
                "input_sha256": sha256_file(input_path),
                "prediction_sha256": sha256_file(prediction_path),
            }
        )
    if common_header is None or keys(merged_records) != expected_global_keys:
        raise ValueError("merged prediction identity/order differs from transport VCF")
    if len(merged_records) != EXPECTED_VARIANTS:
        raise ValueError("merged prediction count mismatch")
    if not any(line.startswith(f"##INFO=<ID={args.predictor},") for line in common_header):
        raise ValueError(f"{args.predictor} INFO definition missing")

    args.output_vcf.parent.mkdir(parents=True, exist_ok=True)
    args.audit_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_vcf.write_text("".join(common_header + merged_records), encoding="utf-8", newline="")
    audit = {
        "status": "PASS",
        "predictor": args.predictor,
        "distance": args.distance,
        "unmasked": True,
        "n_variants": len(merged_records),
        "n_records_with_predictor_info": annotated,
        "identity_and_order_match_transport": True,
        "transport_vcf": {"path": str(args.transport_vcf), "sha256": sha256_file(args.transport_vcf)},
        "shards": shard_audits,
        "output_vcf": {"path": str(args.output_vcf), "size_bytes": args.output_vcf.stat().st_size, "sha256": sha256_file(args.output_vcf)},
        "script": {"path": str(Path(__file__).resolve()), "sha256": sha256_file(Path(__file__).resolve())},
        "outcome_fields_used": [],
    }
    args.audit_json.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "predictor": args.predictor, "distance": args.distance, "variants": len(merged_records), "annotated": annotated, "output_sha256": audit["output_vcf"]["sha256"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
