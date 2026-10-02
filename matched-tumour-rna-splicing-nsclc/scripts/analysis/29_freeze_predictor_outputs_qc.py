#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from collections import Counter
from pathlib import Path


EXPECTED_VARIANTS = 2374
EXPECTED_HASHES = {
    "universe": "55d73120ffaa7fc5f7e14e5673a47edb8e5433a854143bbc0943f8d0b4a8583e",
    "mapping": "a82b30c2385e5cfb0afe7cb641b8cbf72613139a4220d6d11d037e61335cb863",
    "spliceai_d50": "5041fde57acc5fda0a4e6d934d1caedeab70021996f19f7b18d01a901bf11c48",
    "spliceai_d500": "2efb3bef9e9710190160658e6a4795c3bbbc6c8c226fb0fbf5a5bca2aa2a1bb6",
    "pangolin_d50": "5db203e433549080d3c1627f0e968c5310098a8b73473dd485d57d39f5e445bd",
    "pangolin_d500": "289731b99d8d3dbf1e5705e39139b9205bb35cadf693221198215caa4842889d",
}
SPLICEAI_SCORE_NAMES = ("DS_AG", "DS_AL", "DS_DG", "DS_DL")
SPLICEAI_OFFSET_NAMES = ("DP_AG", "DP_AL", "DP_DG", "DP_DL")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_info(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    if text == ".":
        return result
    for item in text.split(";"):
        key, separator, value = item.partition("=")
        if not separator:
            value = "TRUE"
        if not key or key in result:
            raise ValueError(f"invalid/duplicate INFO item: {item!r}")
        result[key] = value
    return result


def read_vcf(path: Path) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"invalid VCF record in {path}: {line!r}")
            records.append(
                {
                    "chrom": fields[0],
                    "pos": int(fields[1]),
                    "id": fields[2],
                    "ref": fields[3],
                    "alt": fields[4],
                    "info": parse_info(fields[7]),
                }
            )
    return records


def record_key(record: dict[str, object]) -> tuple[object, ...]:
    return tuple(record[name] for name in ("chrom", "pos", "id", "ref", "alt"))


def finite_float(text: str, label: str) -> float:
    value = float(text)
    if not math.isfinite(value):
        raise ValueError(f"non-finite {label}: {text!r}")
    return value


def format_number(value: float | int | None) -> str:
    if value is None:
        return ""
    if isinstance(value, int):
        return str(value)
    return format(value, ".17g")


def write_tsv_gz(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=columns, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)


def parse_spliceai_bundles(raw: str) -> list[dict[str, object]]:
    bundles: list[dict[str, object]] = []
    if not raw:
        return bundles
    for bundle_text in raw.split(","):
        fields = bundle_text.split("|")
        if len(fields) != 10:
            raise ValueError(f"invalid SpliceAI bundle: {bundle_text!r}")
        scores = [finite_float(fields[index], SPLICEAI_SCORE_NAMES[index - 2]) for index in range(2, 6)]
        offsets = [int(fields[index]) for index in range(6, 10)]
        if any(score < 0 or score > 1 for score in scores):
            raise ValueError(f"SpliceAI score outside [0,1]: {bundle_text!r}")
        bundles.append(
            {
                "raw": bundle_text,
                "allele": fields[0],
                "spliceai_name": fields[1],
                "scores": scores,
                "offsets": offsets,
            }
        )
    return bundles


def parse_pangolin_bundles(raw: str) -> list[dict[str, object]]:
    bundles: list[dict[str, object]] = []
    if not raw:
        return bundles
    for bundle_text in raw.split(","):
        fields = bundle_text.split("|")
        if len(fields) != 4 or not fields[3].startswith("Warnings:"):
            raise ValueError(f"invalid Pangolin bundle: {bundle_text!r}")
        calls = []
        for call_kind, token in zip(("GAIN", "LOSS"), fields[1:3]):
            offset_text, separator, score_text = token.partition(":")
            if not separator:
                raise ValueError(f"invalid Pangolin call: {token!r}")
            score = finite_float(score_text, "Pangolin score")
            if abs(score) > 1.0000001:
                raise ValueError(f"Pangolin score outside [-1,1]: {score!r}")
            calls.append({"call_kind": call_kind, "offset": int(offset_text), "score": score})
        bundles.append({"raw": bundle_text, "gene_id": fields[0], "calls": calls, "warnings": fields[3]})
    return bundles


def choose_spliceai(
    raw: str,
    universe: dict[str, str],
    mapping_by_spliceai: dict[str, dict[str, str]],
) -> dict[str, object]:
    all_bundles = parse_spliceai_bundles(raw)
    target = []
    for bundle in all_bundles:
        metadata = mapping_by_spliceai.get(str(bundle["spliceai_name"]))
        if metadata is None:
            raise ValueError(f"SpliceAI transcript absent from frozen mapping: {bundle['spliceai_name']}")
        if (
            metadata["gene_name"] == universe["gene"]
            and metadata["chrom"] == universe["chrom"]
            and int(metadata["tx_start_1based"]) <= int(universe["pos_1based"]) <= int(metadata["tx_end_1based"])
        ):
            target.append((bundle, metadata))
    if not target:
        return {"status": "TARGET_NOT_ANNOTATED", "all_raw": raw, "target_raw": "", "target_count": 0}
    target.sort(key=lambda item: (-max(item[0]["scores"]), str(item[0]["spliceai_name"])))
    bundle, metadata = target[0]
    scores = list(bundle["scores"])
    offsets = list(bundle["offsets"])
    max_index = max(range(4), key=lambda index: (scores[index], -index))
    result: dict[str, object] = {
        "status": "SCORED",
        "all_raw": raw,
        "target_raw": ",".join(str(item[0]["raw"]) for item in target),
        "target_count": len(target),
        "target_transcript": bundle["spliceai_name"],
        "target_gene_id": metadata["gene_id"],
        "strand": metadata["strand"],
        "ds_max": scores[max_index],
        "ds_max_type": SPLICEAI_SCORE_NAMES[max_index],
        "ds_max_offset": offsets[max_index],
        "ds_max_coordinate": int(universe["pos_1based"]) + offsets[max_index],
    }
    for name, value in zip(SPLICEAI_SCORE_NAMES, scores):
        result[name] = value
    for name, value in zip(SPLICEAI_OFFSET_NAMES, offsets):
        result[name] = value
        result[name.replace("DP_", "COORD_")] = int(universe["pos_1based"]) + value
    return result


def choose_pangolin(
    raw: str,
    universe: dict[str, str],
    mapping_by_gene_id: dict[str, dict[str, str]],
) -> dict[str, object]:
    all_bundles = parse_pangolin_bundles(raw)
    target = []
    for bundle in all_bundles:
        metadata = mapping_by_gene_id.get(str(bundle["gene_id"]))
        if metadata is None:
            raise ValueError(f"Pangolin gene absent from frozen mapping: {bundle['gene_id']}")
        if (
            metadata["gene_name"] == universe["gene"]
            and metadata["chrom"] == universe["chrom"]
            and int(metadata["tx_start_1based"]) <= int(universe["pos_1based"]) <= int(metadata["tx_end_1based"])
        ):
            target.append((bundle, metadata))
    if not target:
        return {"status": "TARGET_NOT_ANNOTATED", "all_raw": raw, "target_raw": "", "target_count": 0}
    calls = []
    for bundle, metadata in target:
        for call in bundle["calls"]:
            calls.append((call, bundle, metadata))
    calls.sort(key=lambda item: (-abs(float(item[0]["score"])), 0 if item[0]["call_kind"] == "GAIN" else 1, str(item[1]["gene_id"]), int(item[0]["offset"])))
    call, bundle, metadata = calls[0]
    return {
        "status": "SCORED",
        "all_raw": raw,
        "target_raw": ",".join(str(item[0]["raw"]) for item in target),
        "target_count": len(target),
        "target_gene_id": bundle["gene_id"],
        "target_transcript": metadata["spliceai_name"],
        "strand": metadata["strand"],
        "max_abs_delta": abs(float(call["score"])),
        "signed_delta": float(call["score"]),
        "relative_offset": int(call["offset"]),
        "predicted_coordinate": int(universe["pos_1based"]) + int(call["offset"]),
        "direction": call["call_kind"],
        "site_type": "GENERIC_SPLICE_SITE_USAGE",
        "site_type_provided_by_model": False,
        "warnings": bundle["warnings"],
    }


def base_output(universe: dict[str, str]) -> dict[str, object]:
    return {
        "predictor_variant_id": universe["predictor_variant_id"],
        "variant_uid": universe["variant_uid"],
        "chrom": universe["chrom"],
        "pos_1based": universe["pos_1based"],
        "ref": universe["ref"],
        "alt": universe["alt"],
        "gene": universe["gene"],
        "histologies": universe["histologies"],
        "candidate_pair_count": universe["candidate_pair_count"],
        "allele_class": universe["allele_class"],
        "outcome_blind_input": universe["outcome_blind_input"],
    }


def add_prefixed(output: dict[str, object], selected: dict[str, object], prefix: str, numeric_names: set[str]) -> None:
    for key, value in selected.items():
        name = f"{prefix}{key}" if prefix else key
        if key in numeric_names:
            output[name] = format_number(value)  # type: ignore[arg-type]
        elif isinstance(value, bool):
            output[name] = str(value).upper()
        else:
            output[name] = value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    input_dir = root / "v2.3/results/inputs/m4_m5"
    raw_dir = root / "v2.3/results/raw_predictors/m4_m5_phase1"
    tables_dir = root / "v2.3/results/tables"
    audit_dir = root / "v2.3/results/audit"
    phase_audit_dir = audit_dir / "m4_m5_phase1"
    paths = {
        "universe": input_dir / "V23_M4_M5_VARIANT_UNIVERSE.tsv",
        "mapping": input_dir / "V23_GENCODE50_CANONICAL_TRANSCRIPTS.tsv",
        "spliceai_d50": raw_dir / "V23_SPLICEAI_D50_RAW.vcf",
        "spliceai_d500": raw_dir / "V23_SPLICEAI_D500_RAW.vcf",
        "pangolin_d50": raw_dir / "V23_PANGOLIN_D50_RAW.vcf",
        "pangolin_d500": raw_dir / "V23_PANGOLIN_D500_RAW.vcf",
    }
    validation_audit = phase_audit_dir / "V23_M4_M5_VARIANT_REF_VALIDATION_AUDIT.json"
    spliceai_table = tables_dir / "V23_M4_SPLICEAI_ALL_VARIANTS.tsv.gz"
    pangolin_table = tables_dir / "V23_M5_PANGOLIN_ALL_VARIANTS.tsv.gz"
    qc_path = audit_dir / "V23_PREDICTOR_QC.tsv"
    audit_path = phase_audit_dir / "V23_PREDICTOR_RAW_FREEZE_AUDIT.json"
    sha_path = phase_audit_dir / "V23_PREDICTOR_RAW_FREEZE_SHA256.tsv"
    for path in (*paths.values(), validation_audit):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    for role, expected in EXPECTED_HASHES.items():
        observed = sha256_file(paths[role])
        if observed != expected:
            raise ValueError(f"{role} SHA256 mismatch: {observed}")
    validation = json.loads(validation_audit.read_text(encoding="utf-8"))
    if validation.get("status") != "PASS":
        raise ValueError("variant REF validation audit is not PASS")
    for output in (spliceai_table, pangolin_table, qc_path, audit_path, sha_path):
        if output.exists():
            raise FileExistsError(f"refusing to overwrite: {output}")

    universe_rows = read_tsv(paths["universe"])
    mapping_rows = read_tsv(paths["mapping"])
    if len(universe_rows) != EXPECTED_VARIANTS:
        raise ValueError("variant universe count mismatch")
    expected_ids = [f"V23VAR{i:05d}" for i in range(1, EXPECTED_VARIANTS + 1)]
    if [row["predictor_variant_id"] for row in universe_rows] != expected_ids:
        raise ValueError("variant universe IDs/order mismatch")
    mapping_by_spliceai = {row["spliceai_name"]: row for row in mapping_rows}
    mapping_by_gene_id = {row["gene_id"]: row for row in mapping_rows}
    if len(mapping_by_spliceai) != len(mapping_rows) or len(mapping_by_gene_id) != len(mapping_rows):
        raise ValueError("canonical mapping predictor identifiers are not unique")

    raw_records = {role: read_vcf(paths[role]) for role in ("spliceai_d50", "spliceai_d500", "pangolin_d50", "pangolin_d500")}
    universe_keys = [
        (row["chrom"], int(row["pos_1based"]), row["predictor_variant_id"], row["ref"], row["alt"])
        for row in universe_rows
    ]
    for role, records in raw_records.items():
        if len(records) != EXPECTED_VARIANTS or [record_key(record) for record in records] != universe_keys:
            raise ValueError(f"{role} record identity/order mismatch")
        for universe, record in zip(universe_rows, records):
            info = record["info"]
            if info.get("GENE") != universe["gene"] or info.get("HISTOLOGIES") != universe["histologies"] or info.get("PAIR_COUNT") != universe["candidate_pair_count"]:
                raise ValueError(f"{role} frozen input INFO mismatch at {universe['predictor_variant_id']}")

    spliceai_rows: list[dict[str, object]] = []
    pangolin_rows: list[dict[str, object]] = []
    spliceai_numeric = set(SPLICEAI_SCORE_NAMES) | set(SPLICEAI_OFFSET_NAMES) | {name.replace("DP_", "COORD_") for name in SPLICEAI_OFFSET_NAMES} | {"target_count", "ds_max", "ds_max_offset", "ds_max_coordinate"}
    pangolin_numeric = {"target_count", "max_abs_delta", "signed_delta", "relative_offset", "predicted_coordinate"}
    for index, universe in enumerate(universe_rows):
        splice_output = base_output(universe)
        s50 = choose_spliceai(str(raw_records["spliceai_d50"][index]["info"].get("SpliceAI", "")), universe, mapping_by_spliceai)
        s500 = choose_spliceai(str(raw_records["spliceai_d500"][index]["info"].get("SpliceAI", "")), universe, mapping_by_spliceai)
        add_prefixed(splice_output, s50, "", spliceai_numeric)
        add_prefixed(splice_output, s500, "D500_", spliceai_numeric)
        spliceai_rows.append(splice_output)

        pangolin_output = base_output(universe)
        p50 = choose_pangolin(str(raw_records["pangolin_d50"][index]["info"].get("Pangolin", "")), universe, mapping_by_gene_id)
        p500 = choose_pangolin(str(raw_records["pangolin_d500"][index]["info"].get("Pangolin", "")), universe, mapping_by_gene_id)
        add_prefixed(pangolin_output, p50, "", pangolin_numeric)
        add_prefixed(pangolin_output, p500, "D500_", pangolin_numeric)
        pangolin_rows.append(pangolin_output)

    splice_status = Counter(str(row["status"]) for row in spliceai_rows)
    splice500_status = Counter(str(row["D500_status"]) for row in spliceai_rows)
    pangolin_status = Counter(str(row["status"]) for row in pangolin_rows)
    pangolin500_status = Counter(str(row["D500_status"]) for row in pangolin_rows)
    allowed_statuses = {"SCORED", "TARGET_NOT_ANNOTATED"}
    for label, counter in (
        ("spliceai_d50", splice_status), ("spliceai_d500", splice500_status),
        ("pangolin_d50", pangolin_status), ("pangolin_d500", pangolin500_status),
    ):
        if sum(counter.values()) != EXPECTED_VARIANTS or not set(counter).issubset(allowed_statuses):
            raise ValueError(f"unaccounted target predictor status for {label}: {counter}")

    common_columns = list(base_output(universe_rows[0]))
    splice_columns = common_columns + list(dict.fromkeys(key for row in spliceai_rows for key in row if key not in common_columns))
    pangolin_columns = common_columns + list(dict.fromkeys(key for row in pangolin_rows for key in row if key not in common_columns))
    tables_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)
    phase_audit_dir.mkdir(parents=True, exist_ok=True)
    write_tsv_gz(spliceai_table, spliceai_rows, splice_columns)
    write_tsv_gz(pangolin_table, pangolin_rows, pangolin_columns)

    s_monotonic_eligible = [row for row in spliceai_rows if row["status"] == "SCORED" and row["D500_status"] == "SCORED"]
    p_monotonic_eligible = [row for row in pangolin_rows if row["status"] == "SCORED" and row["D500_status"] == "SCORED"]
    s_monotonic_violations = sum(float(row["D500_ds_max"]) + 1e-9 < float(row["ds_max"]) for row in s_monotonic_eligible)
    p_monotonic_violations = sum(float(row["D500_max_abs_delta"]) + 1e-7 < float(row["max_abs_delta"]) for row in p_monotonic_eligible)
    qc_rows = [
        ("ALL", "RAW_RECORD_COUNT", "PASS", "2374 each", "2374 each", "Four raw VCFs"),
        ("ALL", "IDENTITY_AND_ORDER", "PASS", "exact", "exact", "CHROM/POS/ID/REF/ALT versus frozen universe"),
        ("ALL", "FROZEN_INPUT_INFO", "PASS", "exact", "exact", "GENE/HISTOLOGIES/PAIR_COUNT preserved"),
        ("SpliceAI_D50", "TARGET_STATUS_ACCOUNTED", "PASS", json.dumps(dict(splice_status), sort_keys=True), str(EXPECTED_VARIANTS), "Mapped ENST to frozen candidate gene; explicit TARGET_NOT_ANNOTATED retained"),
        ("SpliceAI_D500", "TARGET_STATUS_ACCOUNTED", "PASS", json.dumps(dict(splice500_status), sort_keys=True), str(EXPECTED_VARIANTS), "Mapped ENST to frozen candidate gene; explicit TARGET_NOT_ANNOTATED retained"),
        ("Pangolin_D50", "TARGET_STATUS_ACCOUNTED", "PASS", json.dumps(dict(pangolin_status), sort_keys=True), str(EXPECTED_VARIANTS), "Mapped ENSG to frozen candidate gene; explicit TARGET_NOT_ANNOTATED retained"),
        ("Pangolin_D500", "TARGET_STATUS_ACCOUNTED", "PASS", json.dumps(dict(pangolin500_status), sort_keys=True), str(EXPECTED_VARIANTS), "Mapped ENSG to frozen candidate gene; explicit TARGET_NOT_ANNOTATED retained"),
        ("SpliceAI", "D500_GE_D50_DS_MAX", "PASS" if s_monotonic_violations == 0 else "ADVISORY", str(s_monotonic_violations), "0", "Sensitivity window should contain primary window"),
        ("Pangolin", "D500_GE_D50_MAX_ABS_DELTA", "PASS" if p_monotonic_violations == 0 else "ADVISORY", str(p_monotonic_violations), "0", "Sensitivity window should contain primary window"),
        ("Pangolin", "SITE_TYPE_SEMANTICS", "PASS", "GENERIC_SPLICE_SITE_USAGE", "no fabricated donor/acceptor", "Pangolin 1.0.2 does not emit donor/acceptor class"),
        ("ALL", "OUTCOME_BLIND", "PASS", "[]", "[]", "No M2/M3/frozen9 inputs read"),
    ]
    with qc_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["predictor", "check", "status", "observed", "expected", "detail"])
        writer.writerows(qc_rows)

    script = Path(__file__).resolve()
    audit = {
        "status": "PASS",
        "n_variants": EXPECTED_VARIANTS,
        "outcome_blind": True,
        "outcome_fields_used": [],
        "inputs": {role: {"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)} for role, path in paths.items()},
        "reference_validation_audit": {"path": str(validation_audit), "sha256": sha256_file(validation_audit), "status": "PASS"},
        "mapping": {
            "SpliceAI": "raw ENST spliceai_name -> frozen canonical mapping -> candidate gene symbol/chrom/bounds",
            "Pangolin": "raw versioned ENSG gene_id -> frozen canonical mapping -> candidate gene symbol/chrom/bounds",
        },
        "primary": {"distance": 50, "mask": False, "continuous_score_primary": True},
        "sensitivity": {"distance": 500, "mask": False},
        "target_status_counts": {
            "spliceai_d50": dict(splice_status), "spliceai_d500": dict(splice500_status),
            "pangolin_d50": dict(pangolin_status), "pangolin_d500": dict(pangolin500_status),
        },
        "advisories": {
            "spliceai_d500_lt_d50_ds_max": s_monotonic_violations,
            "pangolin_d500_lt_d50_max_abs_delta": p_monotonic_violations,
        },
        "pangolin_site_type_semantics": "GENERIC_SPLICE_SITE_USAGE; donor/acceptor class is not emitted by Pangolin 1.0.2 and was not inferred",
        "outputs": {
            "spliceai_table": {"path": str(spliceai_table), "size_bytes": spliceai_table.stat().st_size, "sha256": sha256_file(spliceai_table)},
            "pangolin_table": {"path": str(pangolin_table), "size_bytes": pangolin_table.stat().st_size, "sha256": sha256_file(pangolin_table)},
            "qc": {"path": str(qc_path), "size_bytes": qc_path.stat().st_size, "sha256": sha256_file(qc_path)},
        },
        "script": {"path": str(script), "sha256": sha256_file(script)},
    }
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    sha_rows = [("script", script), *[("input", path) for path in paths.values()], ("input_audit", validation_audit), ("output", spliceai_table), ("output", pangolin_table), ("output", qc_path), ("output", audit_path)]
    with sha_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["role", "name", "path", "size_bytes", "sha256"])
        for role, path in sha_rows:
            writer.writerow([role, path.name, str(path), path.stat().st_size, sha256_file(path)])
    print(json.dumps({"status": "PASS", "variants": EXPECTED_VARIANTS, "spliceai_table_sha256": sha256_file(spliceai_table), "pangolin_table_sha256": sha256_file(pangolin_table), "qc_sha256": sha256_file(qc_path), "advisories": audit["advisories"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
