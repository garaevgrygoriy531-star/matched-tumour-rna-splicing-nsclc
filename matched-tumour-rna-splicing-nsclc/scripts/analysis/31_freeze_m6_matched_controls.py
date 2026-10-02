#!/usr/bin/env python3
"""Freeze the outcome-blind M6 matched-control set.

This module intentionally reads only the frozen denominator, Phase-0 event
exclusions, M4/M5 variant universe, v2.2 boundary links, the v2.3 coverage
universe, and the frozen nine case keys.  It does not read predictor scores,
the Step-8 geometry table, target counts/usage, robust-z values, or M2/M3
effect tables.

The executable path is deterministic and single-threaded.  The matching
algorithm is defined at ``histology x variant_uid`` level; controls are
globally without replacement by ``variant_uid`` and the mandatory covariates
are never relaxed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

import numpy as np


SEED = "V23_M6_MATCH_SEED_20260823"
ALLOWED_HISTOLOGIES = ("LUAD", "LUSC")
ALLOWED_ALLELES = {"SNV", "INS", "DEL"}
BOUNDARY_BINS = ("0", "1_2", "3_8", "GE9")
CONSEQUENCE_MAP = {
    "synonymous_variant": "synonymous",
    "missense_variant": "missense",
    "splice_region_variant": "splice-region",
    "splice_donor_region_variant": "splice-region",
    "splice_polypyrimidine_tract_variant": "splice-region",
    "inframe_insertion": "inframe/exon-edge",
    "inframe_deletion": "inframe/exon-edge",
    "protein_altering_variant": "inframe/exon-edge",
}
REQUIRED_COVERAGE_COLUMNS = {
    "histology",
    "variant_uid",
    "patient",
    "cluster_id",
    "carrier_status",
    "cluster_total",
}


class FreezeError(RuntimeError):
    """A frozen-design or input-integrity violation."""


def norm_histology(value: object) -> str:
    text = str(value).strip().replace("TCGA-", "").upper()
    return text


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_size(path: Path) -> int:
    return int(path.stat().st_size)


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def iter_tsv(path: Path) -> Iterator[dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise FreezeError(f"missing header: {path}")
        for row in reader:
            yield {str(k): ("" if v is None else str(v)) for k, v in row.items()}


def require_columns(path: Path, columns: set[str]) -> None:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise FreezeError(f"empty input: {path}") from exc
    missing = sorted(columns.difference(header))
    if missing:
        raise FreezeError(f"{path} missing columns: {missing}")


def as_int(value: object, name: str) -> int:
    text = str(value).strip()
    if text == "":
        raise FreezeError(f"blank integer in {name}")
    try:
        number = float(text)
    except ValueError as exc:
        raise FreezeError(f"non-numeric integer in {name}: {text!r}") from exc
    if not math.isfinite(number) or number != int(number):
        raise FreezeError(f"invalid integer in {name}: {text!r}")
    return int(number)


def as_float(value: object, name: str) -> float:
    text = str(value).strip()
    if text == "":
        raise FreezeError(f"blank float in {name}")
    try:
        number = float(text)
    except ValueError as exc:
        raise FreezeError(f"non-numeric float in {name}: {text!r}") from exc
    if not math.isfinite(number):
        raise FreezeError(f"non-finite float in {name}: {text!r}")
    return number


def consequence_group(value: object) -> str:
    """Exactly reuse the v2.2 Phase-0 mapping."""

    text = str(value)
    return CONSEQUENCE_MAP.get(text, "other")


def boundary_bin(distance: int) -> str:
    if distance < 0:
        raise FreezeError(f"negative variant_boundary_distance: {distance}")
    if distance == 0:
        return "0"
    if distance <= 2:
        return "1_2"
    if distance <= 8:
        return "3_8"
    return "GE9"


def stable_json(payload: object) -> str:
    return json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n"


def deterministic_gzip_write(path: Path, header: Sequence[str], rows: Sequence[Mapping[str, object]]) -> None:
    """Write a byte-deterministic UTF-8 TSV gzip (mtime=0, blank filename)."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    payload_header = "\t".join(str(value) for value in header) + "\n"
    with temp.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            zipped.write(payload_header.encode("utf-8"))
            for row in rows:
                values = []
                for column in header:
                    value = row.get(column, "")
                    if value is None:
                        value = ""
                    if isinstance(value, bool):
                        value = "TRUE" if value else "FALSE"
                    elif isinstance(value, float):
                        if not math.isfinite(value):
                            raise FreezeError(f"non-finite output value {column}: {value}")
                        value = format(value, ".17g")
                    values.append(str(value))
                zipped.write(("\t".join(values) + "\n").encode("utf-8"))
    os.replace(temp, path)


def median_log1p(values: Sequence[float]) -> float:
    if not values:
        raise FreezeError("cannot take median of empty values")
    return float(statistics.median(values))


def numpy_linear_quantiles(values: Sequence[float]) -> tuple[float, float, float]:
    if not values:
        raise FreezeError("cannot compute control cutpoints from an empty pool")
    array = np.asarray(sorted(float(v) for v in values), dtype=np.float64)
    try:
        quantiles = np.quantile(array, [0.25, 0.50, 0.75], method="linear")
    except TypeError:  # NumPy <1.22 compatibility.
        quantiles = np.quantile(array, [0.25, 0.50, 0.75], interpolation="linear")
    return tuple(float(v) for v in quantiles)


def coverage_bin(value: float, cutpoints: Sequence[float]) -> tuple[int, str]:
    if not math.isfinite(value) or len(cutpoints) != 3:
        raise FreezeError(f"invalid coverage/cutpoints: {value!r}, {cutpoints!r}")
    index = int(np.searchsorted(np.asarray(cutpoints, dtype=np.float64), value, side="right"))
    return index, f"Q{index + 1}"


def match_seed(case_key: str, control_key: str) -> str:
    text = f"{SEED}|{case_key}|{control_key}".encode("utf-8")
    return hashlib.sha256(text).hexdigest()


def collapse_denominator(rows: Sequence[Mapping[str, str]]) -> list[dict[str, object]]:
    """Apply the frozen v2.2 denominator rule and deterministic collapse."""

    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    selected_rows = 0
    for row in rows:
        if str(row.get("DNA_class", "")).strip() != "D2":
            continue
        value = as_int(row.get("n_carriers_evaluable", ""), "n_carriers_evaluable")
        if value < 1:
            continue
        selected_rows += 1
        subclass = consequence_group(row.get("primary_consequence", ""))
        if subclass == "other":
            raise FreezeError(
                "unknown D2 consequence_group in denominator: "
                f"{row.get('primary_consequence', '')!r}"
            )
        hist = norm_histology(row.get("histology", ""))
        variant = str(row.get("variant_uid", "")).strip()
        if hist not in ALLOWED_HISTOLOGIES or not variant:
            raise FreezeError(f"invalid denominator histology/variant: {hist!r}, {variant!r}")
        grouped[(hist, variant)].append(
            {
                "histology": hist,
                "variant_uid": variant,
                "gene": str(row.get("gene", "")).strip(),
                "junction_id": str(row.get("junction_id", "")).strip(),
                "cluster_id": str(row.get("cluster_id", "")).strip(),
                "DNA_class": "D2",
                "primary_consequence": str(row.get("primary_consequence", "")).strip(),
                "D2_subclass": subclass,
                "n_carriers_evaluable": value,
            }
        )
    collapsed: list[dict[str, object]] = []
    for key in sorted(grouped):
        candidates = sorted(
            grouped[key],
            key=lambda item: (
                -int(item["n_carriers_evaluable"]),
                str(item["junction_id"]),
                str(item["cluster_id"]),
            ),
        )
        winner = dict(candidates[0])
        collapsed.append(winner)
    if selected_rows == 0 or not collapsed:
        raise FreezeError("no D2 evaluable denominator rows")
    return collapsed


def load_frozen_cases(path: Path) -> list[dict[str, object]]:
    require_columns(path, {"histology", "variant_uid", "gene", "junction_id", "cluster_id", "selected_transcript"})
    cases: list[dict[str, object]] = []
    keys: set[tuple[str, str]] = set()
    for row in iter_tsv(path):
        hist = norm_histology(row["histology"])
        variant = row["variant_uid"].strip()
        if hist not in ALLOWED_HISTOLOGIES or not variant:
            raise FreezeError(f"invalid frozen case key: {hist!r}, {variant!r}")
        key = (hist, variant)
        if key in keys:
            raise FreezeError(f"duplicate frozen case key: {key}")
        keys.add(key)
        cases.append(
            {
                "histology": hist,
                "variant_uid": variant,
                "gene": row["gene"].strip(),
                "junction_id": row["junction_id"].strip(),
                "cluster_id": row["cluster_id"].strip(),
                "selected_transcript": row["selected_transcript"].strip(),
                "case_key": f"{hist}|{variant}",
            }
        )
    if len(cases) != 9:
        raise FreezeError(f"expected exactly 9 frozen cases, observed {len(cases)}")
    return cases


def load_universe(path: Path) -> dict[str, dict[str, object]]:
    require_columns(path, {"variant_uid", "gene", "histologies", "allele_class"})
    universe: dict[str, dict[str, object]] = {}
    for row in iter_tsv(path):
        variant = row["variant_uid"].strip()
        if not variant:
            raise FreezeError(f"blank variant_uid in M4/M5 universe: {path}")
        if variant in universe:
            raise FreezeError(f"duplicate variant_uid in M4/M5 universe: {variant}")
        allele = row["allele_class"].strip()
        if allele not in ALLOWED_ALLELES:
            raise FreezeError(f"invalid universe allele_class for {variant}: {allele!r}")
        histologies = frozenset(
            norm_histology(value)
            for value in row["histologies"].split(",")
            if value.strip()
        )
        if not histologies or not histologies.issubset(set(ALLOWED_HISTOLOGIES)):
            raise FreezeError(f"invalid universe histologies for {variant}: {row['histologies']!r}")
        universe[variant] = {
            "variant_uid": variant,
            "gene": row["gene"].strip(),
            "allele_class": allele,
            "histologies": histologies,
        }
    if not universe:
        raise FreezeError("empty M4/M5 variant universe")
    return universe


def load_supported_variant_uids(path: Path) -> set[str]:
    require_columns(path, {"variant_uid", "R_class_final"})
    supported: set[str] = set()
    for row in iter_tsv(path):
        variant = row["variant_uid"].strip()
        if not variant:
            raise FreezeError(f"blank variant_uid in Phase-0 event table: {path}")
        r_class = row["R_class_final"].strip()
        if r_class.startswith(("R2", "R4", "R5")):
            supported.add(variant)
    return supported


def load_boundary_distances(path: Path) -> dict[str, int]:
    require_columns(path, {"variant_uid", "variant_boundary_distance"})
    distances: dict[str, int] = {}
    seen_counts: Counter[str] = Counter()
    for row in iter_tsv(path):
        variant = row["variant_uid"].strip()
        distance = as_int(row["variant_boundary_distance"], "variant_boundary_distance")
        if variant in distances and distances[variant] != distance:
            raise FreezeError(
                f"variant_boundary_distance is not unique for {variant}: "
                f"{distances[variant]} vs {distance}"
            )
        distances[variant] = distance
        seen_counts[variant] += 1
    if not distances:
        raise FreezeError("empty boundary-link table")
    return distances


def load_coverage(path: Path, relevant_variants: set[str]) -> tuple[dict[tuple[str, str], float], dict[str, int]]:
    """Read TRUE_CARRIER rows and return median log1p cluster totals.

    Duplicate rows are reduced at histology x variant x patient x cluster_id
    by the median of log1p(cluster_total).  The variant-level covariate is the
    median of those deduplicated patient/cluster values.
    """

    require_columns(path, REQUIRED_COVERAGE_COLUMNS)
    dedup_values: dict[tuple[str, str, str, str], list[float]] = defaultdict(list)
    counters: Counter[str] = Counter()
    for row in iter_tsv(path):
        counters["rows_total"] += 1
        if row["carrier_status"].strip() != "TRUE_CARRIER":
            continue
        counters["true_carrier_rows"] += 1
        variant = row["variant_uid"].strip()
        if variant not in relevant_variants:
            counters["true_carrier_rows_outside_relevant_variants"] += 1
            continue
        hist = norm_histology(row["histology"])
        if hist not in ALLOWED_HISTOLOGIES:
            raise FreezeError(f"invalid coverage histology: {hist!r}")
        patient = row["patient"].strip()
        cluster = row["cluster_id"].strip()
        if not patient or not cluster:
            raise FreezeError(f"blank patient/cluster in TRUE_CARRIER coverage row: {row}")
        cluster_total = as_float(row["cluster_total"], "cluster_total")
        if cluster_total < 0:
            raise FreezeError(f"negative cluster_total in TRUE_CARRIER row: {cluster_total}")
        dedup_values[(hist, variant, patient, cluster)].append(math.log1p(cluster_total))
    variant_values: dict[tuple[str, str], list[float]] = defaultdict(list)
    for key, values in dedup_values.items():
        variant_values[key[:2]].append(median_log1p(values))
    variant_medians = {key: median_log1p(values) for key, values in variant_values.items()}
    counters["dedup_patient_cluster_keys"] = len(dedup_values)
    counters["variant_hist_coverage_values"] = len(variant_medians)
    return variant_medians, dict(counters)


def add_case_covariates(
    cases: Sequence[dict[str, object]],
    universe: Mapping[str, Mapping[str, object]],
    boundary_distances: Mapping[str, int],
    coverage: Mapping[tuple[str, str], float],
    cutpoints: Mapping[str, tuple[float, float, float]],
) -> None:
    for case in cases:
        variant = str(case["variant_uid"])
        hist = str(case["histology"])
        if variant not in universe:
            raise FreezeError(f"frozen case absent from M4/M5 universe: {variant}")
        if variant not in boundary_distances:
            raise FreezeError(f"frozen case absent from boundary links: {variant}")
        coverage_key = (hist, variant)
        if coverage_key not in coverage:
            raise FreezeError(f"frozen case has no TRUE_CARRIER coverage: {coverage_key}")
        distance = boundary_distances[variant]
        cov = float(coverage[coverage_key])
        cov_index, cov_bin = coverage_bin(cov, cutpoints[hist])
        case.update(
            {
                "boundary_distance": distance,
                "boundary_bin": boundary_bin(distance),
                "allele_class": universe[variant]["allele_class"],
                "coverage_log1p_median": cov,
                "coverage_bin_index": cov_index,
                "coverage_bin": cov_bin,
            }
        )


def annotate_cases_from_denominator(
    cases: Sequence[dict[str, object]], denominator: Sequence[dict[str, object]]
) -> None:
    """Attach the collapsed denominator subclass/count to each fixed case."""

    denominator_by_key = {
        (str(row["histology"]), str(row["variant_uid"])): row for row in denominator
    }
    for case in cases:
        key = (str(case["histology"]), str(case["variant_uid"]))
        if key not in denominator_by_key:
            raise FreezeError(f"frozen case absent from collapsed D2 denominator: {key}")
        row = denominator_by_key[key]
        if str(case["gene"]) != str(row["gene"]):
            raise FreezeError(f"frozen case gene differs from denominator: {key}")
        case.update(
            {
                "D2_subclass": row["D2_subclass"],
                "n_carriers_evaluable": row["n_carriers_evaluable"],
                "DNA_class": row["DNA_class"],
            }
        )


def make_controls(
    denominator: Sequence[dict[str, object]],
    supported_uids: set[str],
    frozen_uids: set[str],
    universe: Mapping[str, Mapping[str, object]],
    boundary_distances: Mapping[str, int],
    coverage: Mapping[tuple[str, str], float],
) -> tuple[list[dict[str, object]], dict[str, int]]:
    """Apply global exclusions and construct fully covariated controls."""

    counters: Counter[str] = Counter()
    counters["denominator_collapsed_hist_variant"] = len(denominator)
    controls: list[dict[str, object]] = []
    for item in denominator:
        variant = str(item["variant_uid"])
        hist = str(item["histology"])
        if variant in supported_uids:
            counters["excluded_global_phase0_supported"] += 1
            continue
        if variant in frozen_uids:
            counters["excluded_frozen9_case_variant"] += 1
            continue
        if variant not in universe:
            counters["excluded_absent_m4_m5_universe"] += 1
            continue
        if variant not in boundary_distances:
            raise FreezeError(f"control absent from boundary links: {variant}")
        universe_row = universe[variant]
        if str(item["gene"]) != str(universe_row["gene"]):
            raise FreezeError(f"gene mismatch between denominator and universe: {variant}")
        coverage_key = (hist, variant)
        if coverage_key not in coverage:
            counters["excluded_missing_true_carrier_coverage"] += 1
            continue
        distance = int(boundary_distances[variant])
        cov = float(coverage[coverage_key])
        controls.append(
            {
                **item,
                "gene": str(universe_row["gene"]),
                "allele_class": str(universe_row["allele_class"]),
                "boundary_distance": distance,
                "boundary_bin": boundary_bin(distance),
                "coverage_log1p_median": cov,
                "control_key": f"{hist}|{variant}",
                "proxy_status": "R0_PROXY_NO_SUPPORTED_PHASE0_EVENT",
            }
        )
    counters["eligible_control_pool"] = len(controls)
    return controls, dict(counters)


def candidate_sort_key(case: Mapping[str, object], control: Mapping[str, object]) -> tuple[object, ...]:
    case_gene = str(case["gene"])
    control_gene = str(control["gene"])
    coverage_diff = abs(float(control["coverage_log1p_median"]) - float(case["coverage_log1p_median"]))
    boundary_diff = abs(int(control["boundary_distance"]) - int(case["boundary_distance"]))
    case_n = 0 if case["n_carriers_evaluable"] == "NA" else int(case["n_carriers_evaluable"])
    n_diff = abs(int(control["n_carriers_evaluable"]) - case_n)
    case_key = str(case["case_key"])
    control_key = str(control["control_key"])
    return (
        0 if case_gene == control_gene else 1,
        coverage_diff,
        boundary_diff,
        n_diff,
        match_seed(case_key, control_key),
        control_key,
    )


def mandatory_match(case: Mapping[str, object], control: Mapping[str, object]) -> bool:
    return (
        str(case["histology"]) == str(control["histology"])
        and str(case["D2_subclass"]) == str(control["D2_subclass"])
        and str(case["boundary_bin"]) == str(control["boundary_bin"])
        and str(case["allele_class"]) == str(control["allele_class"])
        and str(case["coverage_bin"]) == str(control["coverage_bin"])
    )


def match_controls(
    cases: Sequence[dict[str, object]], controls: Sequence[dict[str, object]],
) -> tuple[list[dict[str, object]], dict[str, dict[str, object]], dict[str, object]]:
    """Greedy frozen allocation, constrained by global no replacement."""

    initial_candidates: dict[str, list[dict[str, object]]] = {}
    for case in cases:
        key = str(case["case_key"])
        initial_candidates[key] = [control for control in controls if mandatory_match(case, control)]
    allocation_order = sorted(
        cases,
        key=lambda case: (len(initial_candidates[str(case["case_key"])]), str(case["case_key"])),
    )
    used_variant_uids: set[str] = set()
    rows: list[dict[str, object]] = []
    case_audit: dict[str, dict[str, object]] = {}
    for allocation_rank, case in enumerate(allocation_order, start=1):
        case_key = str(case["case_key"])
        remaining = [
            control for control in initial_candidates[case_key]
            if str(control["variant_uid"]) not in used_variant_uids
        ]
        ordered = sorted(remaining, key=lambda control: candidate_sort_key(case, control))
        selected = ordered[: min(10, len(ordered))]
        case_audit[case_key] = {
            "allocation_rank": allocation_rank,
            "initial_eligible_candidate_count": len(initial_candidates[case_key]),
            "eligible_remaining_at_assignment": len(remaining),
            "selected_count": len(selected),
            "status": "PASS" if selected else "UNMATCHABLE",
        }
        for rank, control in enumerate(selected, start=1):
            control_uid = str(control["variant_uid"])
            if control_uid in used_variant_uids:
                raise FreezeError(f"global control reuse detected during assignment: {control_uid}")
            used_variant_uids.add(control_uid)
            signed_cov_diff = float(control["coverage_log1p_median"]) - float(case["coverage_log1p_median"])
            signed_boundary_diff = int(control["boundary_distance"]) - int(case["boundary_distance"])
            signed_n_diff = int(control["n_carriers_evaluable"]) - int(case["n_carriers_evaluable"])
            rows.append(
                {
                    "case_key": case_key,
                    "case_histology": case["histology"],
                    "case_variant_uid": case["variant_uid"],
                    "case_gene": case["gene"],
                    "case_junction_id": case["junction_id"],
                    "case_cluster_id": case["cluster_id"],
                    "case_D2_subclass": case["D2_subclass"],
                    "case_n_carriers_evaluable": case["n_carriers_evaluable"],
                    "case_allele_class": case["allele_class"],
                    "case_boundary_distance": case["boundary_distance"],
                    "case_boundary_bin": case["boundary_bin"],
                    "case_coverage_log1p_median": case["coverage_log1p_median"],
                    "case_coverage_bin_index": case["coverage_bin_index"],
                    "case_coverage_bin": case["coverage_bin"],
                    "control_key": control["control_key"],
                    "control_histology": control["histology"],
                    "control_variant_uid": control["variant_uid"],
                    "control_gene": control["gene"],
                    "control_junction_id": control["junction_id"],
                    "control_cluster_id": control["cluster_id"],
                    "control_D2_subclass": control["D2_subclass"],
                    "control_allele_class": control["allele_class"],
                    "control_n_carriers_evaluable": control["n_carriers_evaluable"],
                    "control_boundary_distance": control["boundary_distance"],
                    "control_boundary_bin": control["boundary_bin"],
                    "control_coverage_log1p_median": control["coverage_log1p_median"],
                    "control_coverage_bin_index": control["coverage_bin_index"],
                    "control_coverage_bin": control["coverage_bin"],
                    "control_proxy_status": control["proxy_status"],
                    "same_gene": str(case["gene"]) == str(control["gene"]),
                    "coverage_log1p_diff_control_minus_case": signed_cov_diff,
                    "abs_coverage_log1p_diff": abs(signed_cov_diff),
                    "boundary_distance_diff_control_minus_case": signed_boundary_diff,
                    "abs_boundary_distance_diff": abs(signed_boundary_diff),
                    "n_carriers_evaluable_diff_control_minus_case": signed_n_diff,
                    "abs_n_carriers_evaluable_diff": abs(signed_n_diff),
                    "match_rank": rank,
                    "case_allocation_rank": allocation_rank,
                    "match_seed_sha256": match_seed(case_key, str(control["control_key"])),
                    "match_rule": "M6_EXACT_MANDATORY_NO_RELAXATION",
                }
            )
    rows.sort(key=lambda row: (str(row["case_key"]), int(row["match_rank"])))
    all_control_ids = [str(row["control_variant_uid"]) for row in rows]
    no_reuse = len(all_control_ids) == len(set(all_control_ids))
    exact = all(
        str(row["case_histology"]) == str(row["control_histology"])
        and str(row["case_D2_subclass"]) == str(row["control_D2_subclass"])
        and str(row["case_boundary_bin"]) == str(row["control_boundary_bin"])
        and str(row["case_allele_class"]) == str(row["control_allele_class"])
        and str(row["case_coverage_bin"]) == str(row["control_coverage_bin"])
        for row in rows
    )
    global_audit = {
        "global_control_variant_uid_no_replacement": no_reuse,
        "mandatory_exact_matching_pass": exact,
        "no_relaxation": True,
        "selected_rows": len(rows),
        "selected_unique_control_variant_uids": len(set(all_control_ids)),
    }
    if not no_reuse or not exact:
        raise FreezeError(f"matching QC failed: {global_audit}")
    return rows, case_audit, global_audit


OUTPUT_HEADER = (
    "case_key",
    "case_histology",
    "case_variant_uid",
    "case_gene",
    "case_junction_id",
    "case_cluster_id",
    "case_D2_subclass",
    "case_n_carriers_evaluable",
    "case_allele_class",
    "case_boundary_distance",
    "case_boundary_bin",
    "case_coverage_log1p_median",
    "case_coverage_bin_index",
    "case_coverage_bin",
    "control_key",
    "control_histology",
    "control_variant_uid",
    "control_gene",
    "control_junction_id",
    "control_cluster_id",
    "control_D2_subclass",
    "control_allele_class",
    "control_n_carriers_evaluable",
    "control_boundary_distance",
    "control_boundary_bin",
    "control_coverage_log1p_median",
    "control_coverage_bin_index",
    "control_coverage_bin",
    "control_proxy_status",
    "same_gene",
    "coverage_log1p_diff_control_minus_case",
    "abs_coverage_log1p_diff",
    "boundary_distance_diff_control_minus_case",
    "abs_boundary_distance_diff",
    "n_carriers_evaluable_diff_control_minus_case",
    "abs_n_carriers_evaluable_diff",
    "match_rank",
    "case_allocation_rank",
    "match_seed_sha256",
    "match_rule",
)


def run_self_test() -> dict[str, object]:
    """Small synthetic tests for the frozen mechanics, without project inputs."""

    denominator = collapse_denominator(
        [
            {"histology": "LUAD", "variant_uid": "v1", "gene": "G", "junction_id": "j2", "cluster_id": "c2", "DNA_class": "D2", "primary_consequence": "missense_variant", "n_carriers_evaluable": "1"},
            {"histology": "LUAD", "variant_uid": "v1", "gene": "G", "junction_id": "j1", "cluster_id": "c1", "DNA_class": "D2", "primary_consequence": "missense_variant", "n_carriers_evaluable": "2"},
            {"histology": "LUAD", "variant_uid": "v2", "gene": "H", "junction_id": "j1", "cluster_id": "c1", "DNA_class": "D2", "primary_consequence": "synonymous_variant", "n_carriers_evaluable": "1"},
        ]
    )
    assert len(denominator) == 2 and denominator[0]["variant_uid"] == "v1" and denominator[0]["junction_id"] == "j1"
    universe = {
        "v1": {"gene": "G", "allele_class": "SNV"},
        "v2": {"gene": "H", "allele_class": "SNV"},
        "v3": {"gene": "G", "allele_class": "SNV"},
    }
    coverage = {("LUAD", "v1"): 1.0, ("LUAD", "v2"): 2.0, ("LUAD", "v3"): 3.0}
    controls, _ = make_controls(denominator, {"v2"}, {"v9"}, universe, {"v1": 0, "v2": 2, "v3": 3}, coverage)
    assert len(controls) == 1 and controls[0]["variant_uid"] == "v1"
    cuts = {"LUAD": numpy_linear_quantiles([1.0, 1.0, 2.0, 3.0])}
    assert coverage_bin(cuts["LUAD"][0], cuts["LUAD"]) [0] >= 1
    case = {"histology": "LUAD", "variant_uid": "case", "gene": "G", "junction_id": "j", "cluster_id": "c", "case_key": "LUAD|case", "D2_subclass": "missense", "allele_class": "SNV", "boundary_distance": 0, "boundary_bin": "0", "coverage_log1p_median": 1.0, "coverage_bin_index": 0, "coverage_bin": "Q1", "n_carriers_evaluable": 1}
    controls_for_match = []
    for index in range(3):
        controls_for_match.append({"histology": "LUAD", "variant_uid": f"ctrl{index}", "gene": "G", "junction_id": "j", "cluster_id": "c", "D2_subclass": "missense", "allele_class": "SNV", "n_carriers_evaluable": 1, "boundary_distance": 0, "boundary_bin": "0", "coverage_log1p_median": 1.0, "coverage_bin_index": 0, "coverage_bin": "Q1", "control_key": f"LUAD|ctrl{index}", "proxy_status": "R0_PROXY_NO_SUPPORTED_PHASE0_EVENT"})
    matched, audit, global_audit = match_controls([case, {**case, "case_key": "LUAD|case2", "variant_uid": "case2"}], controls_for_match)
    assert len(matched) == 3 and global_audit["global_control_variant_uid_no_replacement"] and all(v["selected_count"] >= 0 for v in audit.values())
    impossible = {**case, "case_key": "LUAD|impossible", "coverage_bin": "Q4"}
    matched_impossible, audit_impossible, _ = match_controls([impossible], controls_for_match)
    assert not matched_impossible and audit_impossible["LUAD|impossible"]["status"] == "UNMATCHABLE"
    return {
        "status": "PASS",
        "collapse": "PASS",
        "exclusions": "PASS",
        "quantile_ties_and_searchsorted_right": "PASS",
        "global_no_replacement": "PASS",
        "no_relaxation": "PASS",
    }


def project_paths(project_root: Path) -> dict[str, Path]:
    # Local exports keep v2.2 under ``<project>/v2.2/final_freeze/v2.2``;
    # the server keeps the same frozen files under ``<project>/final_freeze/v2.2``.
    v22_candidates = (
        project_root.parent / "v2.2" / "final_freeze" / "v2.2",
        project_root.parent / "final_freeze" / "v2.2",
    )
    v22 = next((candidate for candidate in v22_candidates if (candidate / "FULL74_RECURRENT_MODELS.tsv.gz").is_file()), v22_candidates[0])
    return {
        "recurrent": v22 / "FULL74_RECURRENT_MODELS.tsv.gz",
        "phase0": v22 / "A_PHASE0_EVENT_TABLE_V22.tsv",
        "links": v22 / "FULL74_VARIANT_JUNCTION_CANDIDATE_LINKS.tsv.gz",
        "frozen9": project_root / "results" / "audit" / "V23_POST_M2_M3_FROZEN9_KEYS.tsv",
        "universe": project_root / "results" / "inputs" / "m4_m5" / "V23_M4_M5_VARIANT_UNIVERSE.tsv",
        "coverage": project_root / "results" / "tables" / "V23_RARE_D2_NULL_UNIVERSE.tsv.gz",
        "table": project_root / "results" / "tables" / "V23_M6_MATCHED_CONTROLS.tsv.gz",
        "audit": project_root / "results" / "audit" / "m4_m5_phase1" / "V23_M6_MATCH_FREEZE_AUDIT.json",
        "manifest": project_root / "results" / "audit" / "m4_m5_phase1" / "V23_M6_MATCH_FREEZE_SHA256.tsv",
    }


def check_inputs(paths: Mapping[str, Path]) -> None:
    for label in ("recurrent", "phase0", "links", "frozen9", "universe", "coverage"):
        path = paths[label]
        if not path.is_file():
            raise FreezeError(f"missing required input {label}: {path}")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=None)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--allow-overwrite", action="store_true", help="allow replacement of existing formal outputs")
    args = parser.parse_args(argv)
    if args.self_test:
        result = run_self_test()
        print(stable_json(result), end="")
        return 0

    project_root = (args.project_root or Path(__file__).resolve().parents[1]).resolve()
    paths = project_paths(project_root)
    check_inputs(paths)
    if not args.allow_overwrite and any(paths[label].exists() for label in ("table", "audit", "manifest")):
        existing = [str(paths[label]) for label in ("table", "audit", "manifest") if paths[label].exists()]
        raise FreezeError(f"formal output already exists; use --allow-overwrite only for an explicitly authorized rerun: {existing}")

    # Hash every actual input and this implementation before processing.
    input_labels = ("recurrent", "phase0", "links", "frozen9", "universe", "coverage")
    input_sha = {label: sha256_file(paths[label]) for label in input_labels}
    script_path = Path(__file__).resolve()
    input_sha["script_31_freeze_m6_matched_controls.py"] = sha256_file(script_path)

    cases = load_frozen_cases(paths["frozen9"])
    frozen_uids = {str(case["variant_uid"]) for case in cases}
    universe = load_universe(paths["universe"])
    supported = load_supported_variant_uids(paths["phase0"])
    boundary_distances = load_boundary_distances(paths["links"])
    recurrent_rows = list(iter_tsv(paths["recurrent"]))
    denominator = collapse_denominator(recurrent_rows)
    denominator_hist_counts = dict(sorted(Counter(str(row["histology"]) for row in denominator).items()))
    if len(denominator) < 1:
        raise FreezeError("frozen D2 denominator is empty")
    annotate_cases_from_denominator(cases, denominator)
    relevant_variants = {str(row["variant_uid"]) for row in denominator} | frozen_uids
    coverage, coverage_counts = load_coverage(paths["coverage"], relevant_variants)
    controls, pool_counts = make_controls(denominator, supported, frozen_uids, universe, boundary_distances, coverage)
    if not 1 <= len(controls) <= len(denominator):
        raise FreezeError(f"predictor-qualified control pool count conflict: observed {len(controls)}")

    # Cutpoints are computed from eligible controls only and then applied to cases.
    controls_by_hist = defaultdict(list)
    for control in controls:
        controls_by_hist[str(control["histology"])].append(float(control["coverage_log1p_median"]))
    cutpoints: dict[str, tuple[float, float, float]] = {}
    for hist in ALLOWED_HISTOLOGIES:
        cutpoints[hist] = numpy_linear_quantiles(controls_by_hist[hist])
        for control in controls:
            if str(control["histology"]) == hist:
                index, label = coverage_bin(float(control["coverage_log1p_median"]), cutpoints[hist])
                control["coverage_bin_index"] = index
                control["coverage_bin"] = label
    add_case_covariates(cases, universe, boundary_distances, coverage, cutpoints)
    match_rows, case_audit, global_audit = match_controls(cases, controls)

    # A deterministic output is written once; the formal path refuses silent overwrite.
    deterministic_gzip_write(paths["table"], OUTPUT_HEADER, match_rows)
    output_sha = {"V23_M6_MATCHED_CONTROLS.tsv.gz": sha256_file(paths["table"])}
    output_sizes = {"V23_M6_MATCHED_CONTROLS.tsv.gz": file_size(paths["table"])}
    audit = {
        "status": "PASS",
        "analysis": "M6_MATCHED_CONTROLS_PHASE1_OUTCOME_BLIND_FREEZE",
        "seed": SEED,
        "design": {
            "analysis_unit": "histology x variant_uid",
            "case_count": len(cases),
            "case_keys_fixed_from": str(paths["frozen9"]),
            "case_selection": "fixed_frozen9_no_reselection",
            "denominator_rule": "v2.2 FULL74_RECURRENT_MODELS: DNA_class=D2 and n_carriers_evaluable>=1; deterministic collapse by n_carriers_evaluable desc, junction_id asc, cluster_id asc",
            "consequence_group_mapping": "v2.2 30_finalize_phase0.py exact mapping",
            "global_supported_exclusion": "all variant_uid with Phase-0 R_class_final starting R2/R4/R5; frozen9 also excluded globally",
            "control_universe_rule": "variant_uid must exist in V23_M4_M5_VARIANT_UNIVERSE.tsv",
            "allele_rule": "exact allele_class in SNV/INS/DEL",
            "boundary_bins": list(BOUNDARY_BINS),
            "coverage_rule": "TRUE_CARRIER only; deduplicate histology x variant_uid x patient x cluster_id using median(log1p(cluster_total)), then median across deduplicated patient/cluster values",
            "coverage_cutpoints": "numpy linear q25/q50/q75 computed separately on eligible LUAD/LUSC control pools; case uses control cutpoints; numpy.searchsorted(side=right)",
            "mandatory_covariates": ["histology", "D2_subclass", "boundary_bin", "allele_class", "coverage_bin"],
            "ranking": "same gene, abs continuous logcoverage difference, abs raw boundary-distance difference, abs n_carriers_evaluable difference, SHA256(seed|case_key|control_key)",
            "global_control_no_replacement": True,
            "per_case_max_controls": 10,
            "no_relaxation": True,
        },
        "outcome_fields_used": [],
        "predictor_fields_used": [],
        "forbidden_fields_not_used": ["target_count", "junction_usage", "usage", "robust_z", "Step8_geometry", "M2_effect", "M3_effect"],
        "input_sha256": input_sha,
        "counts": {
            "denominator_collapsed_hist_variant": len(denominator),
            "denominator_histology_counts": denominator_hist_counts,
            "global_phase0_supported_variant_uids": len(supported),
            "frozen9_variant_uids": len(frozen_uids),
            "m4_m5_universe_variant_uids": len(universe),
            "boundary_link_variant_uids": len(boundary_distances),
            "pool_steps": pool_counts,
        },
        "coverage": {
            **coverage_counts,
            "eligible_control_histology_counts": dict(sorted(Counter(str(row["histology"]) for row in controls).items())),
            "control_cutpoints_numpy_linear": {hist: list(cutpoints[hist]) for hist in ALLOWED_HISTOLOGIES},
            "control_bin_counts": {
                hist: dict(sorted(Counter(str(row["coverage_bin"]) for row in controls if str(row["histology"]) == hist).items()))
                for hist in ALLOWED_HISTOLOGIES
            },
        },
        "cases": {key: case_audit[key] for key in sorted(case_audit)},
        "global_matching_qc": global_audit,
        "outputs": {"V23_M6_MATCHED_CONTROLS.tsv.gz": {"rows": len(match_rows), "sha256": output_sha["V23_M6_MATCHED_CONTROLS.tsv.gz"], "size_bytes": output_sizes["V23_M6_MATCHED_CONTROLS.tsv.gz"]}},
    }
    paths["audit"].parent.mkdir(parents=True, exist_ok=True)
    paths["audit"].write_text(stable_json(audit), encoding="utf-8")
    manifest_targets: list[tuple[str, str, Path]] = [
        ("input", label, paths[label]) for label in input_labels
    ]
    manifest_targets.extend(
        [
            ("input", "script_31_freeze_m6_matched_controls.py", script_path),
            ("output", "V23_M6_MATCHED_CONTROLS.tsv.gz", paths["table"]),
            ("output", "V23_M6_MATCH_FREEZE_AUDIT.json", paths["audit"]),
        ]
    )
    manifest_rows = []
    for kind, label, path in manifest_targets:
        manifest_rows.append({"kind": kind, "label": label, "path": str(path), "sha256": sha256_file(path), "size_bytes": file_size(path)})
    manifest_rows.sort(key=lambda row: (row["kind"], row["label"], row["path"]))
    paths["manifest"].write_text(
        "kind\tlabel\tpath\tsha256\tsize_bytes\n" + "".join(
            f"{row['kind']}\t{row['label']}\t{row['path']}\t{row['sha256']}\t{row['size_bytes']}\n" for row in manifest_rows
        ),
        encoding="utf-8",
    )
    print(stable_json({"status": "PASS", "table": str(paths["table"]), "rows": len(match_rows), "case_selected_counts": {key: case_audit[key]["selected_count"] for key in sorted(case_audit)}, "global_matching_qc": global_audit, "table_sha256": output_sha["V23_M6_MATCHED_CONTROLS.tsv.gz"], "audit_sha256": sha256_file(paths["audit"]), "manifest_sha256": sha256_file(paths["manifest"])}), end="")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FreezeError as exc:
        print(f"M6_FREEZE_FAIL: {exc}", file=sys.stderr)
        raise SystemExit(2)
