#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


EXPECTED_VCF_SHA256 = "0274386b83094adb5b9e5796568aa29b7fe6b124eb516f633d475072b34a7f99"
EXPECTED_FAI_SHA256 = "a2c323ea4cff34d7123ace4578f7e122b2d2f5a22f40dc23eb8b97d17723d169"
EXPECTED_VARIANTS = 2374


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--shards", type=int, default=8)
    args = parser.parse_args()
    if args.shards < 1 or args.shards > 16:
        raise ValueError("--shards must be within 1..16")

    root = args.root.resolve()
    source = root / "v2.3/results/inputs/m4_m5/V23_M4_M5_ALL_UNIQUE_VARIANTS.vcf"
    fai = root / "01_data_acquisition/annotations/GRCh38.primary_assembly.genome.fa.fai"
    output_dir = root / "v2.3/results/inputs/m4_m5/predictor_transport"
    transport = output_dir / "V23_M4_M5_ALL_UNIQUE_VARIANTS_WITH_CONTIG.vcf"
    manifest = output_dir / "V23_M4_M5_TRANSPORT_SHARDS.tsv"
    audit_dir = root / "v2.3/results/audit/m4_m5_phase1"
    audit_path = audit_dir / "V23_M4_M5_TRANSPORT_VCF_AUDIT.json"
    for path in (source, fai):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    if sha256_file(source) != EXPECTED_VCF_SHA256:
        raise ValueError("source VCF SHA256 mismatch")
    if sha256_file(fai) != EXPECTED_FAI_SHA256:
        raise ValueError("FAI SHA256 mismatch")

    shard_paths = [output_dir / f"shard_{index:03d}.vcf" for index in range(args.shards)]
    for target in (transport, manifest, audit_path, *shard_paths):
        if target.exists():
            raise FileExistsError(f"refusing to overwrite: {target}")

    source_lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    chrom_index = next((i for i, line in enumerate(source_lines) if line.startswith("#CHROM\t")), None)
    if chrom_index is None:
        raise ValueError("VCF #CHROM header missing")
    meta = source_lines[:chrom_index]
    column_header = source_lines[chrom_index]
    records = source_lines[chrom_index + 1 :]
    if len(records) != EXPECTED_VARIANTS or any(not line.endswith("\n") for line in records):
        raise ValueError(f"unexpected VCF record structure/count: {len(records)}")
    ids = [line.split("\t", 3)[2] for line in records]
    if len(set(ids)) != EXPECTED_VARIANTS or ids != [f"V23VAR{i:05d}" for i in range(1, EXPECTED_VARIANTS + 1)]:
        raise ValueError("variant IDs are not unique sequential V23VAR00001..V23VAR02374")

    contigs: list[str] = []
    contig_names: set[str] = set()
    for raw in fai.read_text(encoding="utf-8").splitlines():
        fields = raw.split("\t")
        if len(fields) < 2 or not fields[1].isdigit() or fields[0] in contig_names:
            raise ValueError(f"invalid/duplicate FAI row: {raw!r}")
        contig_names.add(fields[0])
        contigs.append(f"##contig=<ID={fields[0]},length={fields[1]}>\n")
    record_contigs = {line.split("\t", 1)[0] for line in records}
    if not record_contigs.issubset(contig_names):
        raise ValueError(f"VCF contigs absent from FAI: {sorted(record_contigs - contig_names)}")
    if any(line.startswith("##contig=") for line in meta):
        raise ValueError("source VCF unexpectedly already contains contig metadata")

    output_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)
    transport_header = meta + contigs + [column_header]
    transport.write_text("".join(transport_header + records), encoding="utf-8", newline="")
    if transport.read_text(encoding="utf-8").splitlines(keepends=True)[-(len(records)) :] != records:
        raise ValueError("transport VCF data lines differ from frozen source")

    base, remainder = divmod(EXPECTED_VARIANTS, args.shards)
    shard_rows: list[dict[str, object]] = []
    start = 0
    reconstructed: list[str] = []
    for index, path in enumerate(shard_paths):
        count = base + (1 if index < remainder else 0)
        subset = records[start : start + count]
        path.write_text("".join(transport_header + subset), encoding="utf-8", newline="")
        reconstructed.extend(subset)
        shard_rows.append(
            {
                "shard_index": index,
                "path": str(path),
                "n_variants": count,
                "first_id": ids[start],
                "last_id": ids[start + count - 1],
                "sha256": sha256_file(path),
            }
        )
        start += count
    if reconstructed != records or start != EXPECTED_VARIANTS:
        raise ValueError("shards do not reconstruct frozen VCF records exactly")

    with manifest.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(shard_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(shard_rows)

    script = Path(__file__).resolve()
    audit = {
        "status": "PASS",
        "source_vcf": {"path": str(source), "sha256": sha256_file(source), "n_variants": len(records)},
        "fai": {"path": str(fai), "sha256": sha256_file(fai), "n_contigs": len(contigs)},
        "transport_vcf": {"path": str(transport), "size_bytes": transport.stat().st_size, "sha256": sha256_file(transport)},
        "transformation": "insert FAI-derived ##contig headers immediately before #CHROM; variant data lines byte-identical",
        "variant_data_lines_byte_identical": True,
        "variant_ids_unique_and_sequential": True,
        "shards": shard_rows,
        "shards_reconstruct_variant_data_lines_exactly": True,
        "manifest": {"path": str(manifest), "sha256": sha256_file(manifest)},
        "script": {"path": str(script), "sha256": sha256_file(script)},
        "outcome_fields_used": [],
    }
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "variants": len(records), "contigs": len(contigs), "shards": args.shards, "transport_sha256": sha256_file(transport)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
