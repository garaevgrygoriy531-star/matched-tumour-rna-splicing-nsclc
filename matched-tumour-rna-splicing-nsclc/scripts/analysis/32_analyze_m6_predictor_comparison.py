#!/usr/bin/env python3
"""Step 9, phase 2: continuous predictor comparison on frozen M6 sets.

The M6 control table is already frozen and outcome-blind.  This script only
joins the four frozen predictor score columns to those 9 disjoint matched
sets and performs the pre-specified matched-set calculations.  It deliberately
does not read the Step-8 geometry table, RNA usage/target fields, M2/M3
effects, or any other outcome field.

The exact permutation test treats each matched set as a small label
permutation block: one of the K+1 units is selected as pseudo-case and is
compared with the remaining K units.  Per-set Fraction distributions are
convolved with a Counter, so the reported p values are exact rational
numerator/denominator values (no pseudocount).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import os
import statistics
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


EXPECTED_MATCH_SHA256 = "cc7b8f783fbd908b317d074f1f9aeb56b0b18239bfb89ad8a9fc74a795017096"
EXPECTED_MATCH_ROWS = 74
EXPECTED_CASE_COUNT = 9
EXPECTED_CASE_K = {
    "LUAD|chr12:25227342:T:A": 10,
    "LUAD|chr19:10489648:C:A": 10,
    "LUAD|chr19:1220373:G:T": 3,
    "LUAD|chr19:45420410:C:T": 10,
    "LUAD|chr6:43443117:G:T": 10,
    "LUSC|chr13:48307406:G:T": 10,
    "LUSC|chr17:7673570:T:C": 10,
    "LUSC|chr17:7676379:T:A": 1,
    "LUSC|chr2:224481892:C:A": 10,
}
EXPECTED_STEP7 = {
    "results/tables/V23_M4_SPLICEAI_ALL_VARIANTS.tsv.gz": "525e896f358ee1d6c7a2e83ae54207c353c31a0e6861968af4dfb420358f68f6",
    "results/tables/V23_M5_PANGOLIN_ALL_VARIANTS.tsv.gz": "1c9f75d93c48e9d67f031ceeb5c571bf3ca2c916fd3e6fd8a3936ee31797dc2d",
    "results/audit/m4_m5_phase1/V23_PREDICTOR_RAW_FREEZE_AUDIT.json": "cea4ff166808319d8b2998d629a904a37b6f6874a50adab72fb613d6922447d9",
}
EXPECTED_MATCH_AUDIT_FIELDS = (
    "outcome_fields_used",
    "predictor_fields_used",
    "design",
    "global_matching_qc",
    "status",
)
THRESHOLDS = (Fraction("0.2"), Fraction("0.5"), Fraction("0.8"))

MATCH_REQUIRED = {
    "case_key",
    "case_histology",
    "case_variant_uid",
    "case_gene",
    "case_D2_subclass",
    "case_allele_class",
    "case_boundary_bin",
    "case_coverage_bin",
    "control_key",
    "control_histology",
    "control_variant_uid",
    "control_gene",
    "control_D2_subclass",
    "control_allele_class",
    "control_boundary_bin",
    "control_coverage_bin",
    "control_proxy_status",
    "same_gene",
    "match_rank",
    "match_rule",
}
FROZEN_REQUIRED = {"histology", "variant_uid", "gene", "junction_id", "cluster_id", "selected_transcript"}
PREDICTOR_REQUIRED = {
    "variant_uid",
    "status",
    "D500_status",
}


class ComparisonError(RuntimeError):
    """Raised for an input-integrity or pre-specified-analysis violation."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_size(path: Path) -> int:
    return int(path.stat().st_size)


def stable_json(value: object) -> str:
    def default(item: object) -> object:
        if isinstance(item, Fraction):
            return fraction_text(item)
        if isinstance(item, Path):
            return str(item)
        raise TypeError(f"unsupported JSON value: {type(item)!r}")

    return json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False, default=default) + "\n"


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ComparisonError(f"missing TSV header: {path}")
        return [{str(k): ("" if v is None else str(v)) for k, v in row.items()} for row in reader]


def require_columns(path: Path, rows: Sequence[Mapping[str, str]], columns: Iterable[str]) -> None:
    if not rows:
        raise ComparisonError(f"empty input table: {path}")
    missing = sorted(set(columns).difference(rows[0]))
    if missing:
        raise ComparisonError(f"{path} missing columns: {missing}")


def as_int(value: object, label: str) -> int:
    text = str(value).strip()
    try:
        number = int(text)
    except ValueError as exc:
        raise ComparisonError(f"invalid integer {label}: {text!r}") from exc
    return number


def as_fraction(value: object, label: str) -> Fraction:
    text = str(value).strip()
    if not text:
        raise ComparisonError(f"blank score {label}")
    try:
        result = Fraction(text)
    except (ValueError, ZeroDivisionError) as exc:
        raise ComparisonError(f"invalid score {label}: {text!r}") from exc
    if not math.isfinite(float(result)):
        raise ComparisonError(f"non-finite score {label}: {text!r}")
    return result


def fraction_float(value: Fraction | None) -> float | None:
    if value is None:
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ComparisonError(f"non-finite result from Fraction: {value}")
    return result


def format_fraction(value: Fraction | None) -> str:
    """Stable decimal output for Fraction-derived summaries."""

    if value is None:
        return "NA"
    result = float(value)
    if not math.isfinite(result):
        raise ComparisonError(f"non-finite summary: {value}")
    return format(result, ".17g")


def fraction_text(value: Fraction) -> str:
    """Exact rational text used only in self-test/audit diagnostics."""

    return f"{value.numerator}/{value.denominator}"


def parse_bool(value: object, label: str) -> bool:
    text = str(value).strip().upper()
    if text == "TRUE":
        return True
    if text == "FALSE":
        return False
    raise ComparisonError(f"invalid boolean {label}: {value!r}")


def deterministic_gzip_write(path: Path, header: Sequence[str], rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            zipped.write(("\t".join(header) + "\n").encode("utf-8"))
            for row in rows:
                values: list[str] = []
                for column in header:
                    value = row.get(column, "")
                    if value is None:
                        value = ""
                    elif isinstance(value, bool):
                        value = "TRUE" if value else "FALSE"
                    elif isinstance(value, float):
                        if not math.isfinite(value):
                            raise ComparisonError(f"non-finite output {column}: {value}")
                        value = format(value, ".17g")
                    values.append(str(value))
                zipped.write(("\t".join(values) + "\n").encode("utf-8"))
    os.replace(temp, path)


def deterministic_tsv_write(path: Path, header: Sequence[str], rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            values: list[str] = []
            for column in header:
                value = row.get(column, "")
                if value is None:
                    value = ""
                elif isinstance(value, bool):
                    value = "TRUE" if value else "FALSE"
                elif isinstance(value, float):
                    if not math.isfinite(value):
                        raise ComparisonError(f"non-finite output {column}: {value}")
                    value = format(value, ".17g")
                values.append(str(value))
            handle.write("\t".join(values) + "\n")
    os.replace(temp, path)


def load_frozen_cases(path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(path)
    require_columns(path, rows, FROZEN_REQUIRED)
    cases: dict[str, dict[str, str]] = {}
    for row in rows:
        hist = row["histology"].strip()
        variant = row["variant_uid"].strip()
        key = f"{hist}|{variant}"
        if not hist or not variant or key in cases:
            raise ComparisonError(f"invalid/duplicate frozen case key: {key!r}")
        cases[key] = row
    if set(cases) != set(EXPECTED_CASE_K):
        raise ComparisonError(f"frozen9 key set conflict: {sorted(cases)}")
    return cases


def load_match_table(path: Path, cases: Mapping[str, Mapping[str, str]]) -> tuple[dict[str, list[dict[str, str]]], dict[str, object]]:
    observed_hash = sha256_file(path)
    if observed_hash != EXPECTED_MATCH_SHA256:
        raise ComparisonError(f"M6 match table SHA256 conflict: {observed_hash} != {EXPECTED_MATCH_SHA256}")
    rows = read_tsv(path)
    require_columns(path, rows, MATCH_REQUIRED)
    if len(rows) != EXPECTED_MATCH_ROWS:
        raise ComparisonError(f"M6 match row count conflict: {len(rows)} != {EXPECTED_MATCH_ROWS}")
    grouped: dict[str, list[dict[str, str]]] = {key: [] for key in cases}
    control_uids: list[str] = []
    for row in rows:
        case_key = row["case_key"].strip()
        if case_key not in grouped:
            raise ComparisonError(f"match table contains unknown case: {case_key}")
        if row["case_key"] != case_key or row["case_histology"] != cases[case_key]["histology"] or row["case_variant_uid"] != cases[case_key]["variant_uid"]:
            raise ComparisonError(f"case identity mismatch in match row: {case_key}")
        if row["case_gene"] != cases[case_key]["gene"]:
            raise ComparisonError(f"case gene mismatch in match row: {case_key}")
        if row["case_histology"] != row["control_histology"]:
            raise ComparisonError(f"histology matching violation: {case_key}")
        for pair in (
            ("case_D2_subclass", "control_D2_subclass"),
            ("case_allele_class", "control_allele_class"),
            ("case_boundary_bin", "control_boundary_bin"),
            ("case_coverage_bin", "control_coverage_bin"),
        ):
            if row[pair[0]] != row[pair[1]]:
                raise ComparisonError(f"mandatory matching violation {pair} in {case_key}")
        if row["control_proxy_status"] != "R0_PROXY_NO_SUPPORTED_PHASE0_EVENT":
            raise ComparisonError(f"control proxy status conflict in {case_key}")
        if row["match_rule"] != "M6_EXACT_MANDATORY_NO_RELAXATION":
            raise ComparisonError(f"match rule conflict in {case_key}")
        parse_bool(row["same_gene"], f"same_gene {case_key}")
        as_int(row["match_rank"], f"match_rank {case_key}")
        control_uids.append(row["control_variant_uid"])
        grouped[case_key].append(row)
    observed_k = {key: len(value) for key, value in grouped.items()}
    if observed_k != EXPECTED_CASE_K:
        raise ComparisonError(f"matched-set K conflict: {observed_k} != {EXPECTED_CASE_K}")
    if len(control_uids) != len(set(control_uids)):
        raise ComparisonError("global control variant_uid reuse detected")
    for key, value in grouped.items():
        ranks = sorted(as_int(row["match_rank"], f"match_rank {key}") for row in value)
        if ranks != list(range(1, len(value) + 1)):
            raise ComparisonError(f"non-contiguous match ranks for {key}: {ranks}")
        value.sort(key=lambda row: as_int(row["match_rank"], f"match_rank {key}"))
    audit = {
        "rows": len(rows),
        "case_count": len(grouped),
        "case_K": {key: observed_k[key] for key in sorted(observed_k)},
        "unique_control_variant_uids": len(set(control_uids)),
        "global_no_replacement": len(control_uids) == len(set(control_uids)),
    }
    return grouped, audit


def load_match_audit(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    for field in EXPECTED_MATCH_AUDIT_FIELDS:
        if field not in payload:
            raise ComparisonError(f"M6 match audit missing field: {field}")
    if payload["status"] != "PASS":
        raise ComparisonError(f"M6 match audit status is not PASS: {payload['status']!r}")
    if payload["outcome_fields_used"] != [] or payload["predictor_fields_used"] != []:
        raise ComparisonError("M6 matching audit is not outcome/predictor blind")
    design = payload["design"]
    if not design.get("global_control_no_replacement") or design.get("no_relaxation") is not True:
        raise ComparisonError("M6 matching audit design flags conflict")
    qc = payload["global_matching_qc"]
    if not qc.get("global_control_variant_uid_no_replacement") or not qc.get("mandatory_exact_matching_pass") or not qc.get("no_relaxation"):
        raise ComparisonError("M6 matching audit QC flags conflict")
    if qc.get("selected_rows") != EXPECTED_MATCH_ROWS or qc.get("selected_unique_control_variant_uids") != EXPECTED_MATCH_ROWS:
        raise ComparisonError("M6 matching audit row/QC count conflict")
    return payload


def load_predictors(path: Path, kind: str) -> dict[str, dict[str, object]]:
    rows = read_tsv(path)
    required = PREDICTOR_REQUIRED | ({"ds_max", "D500_ds_max"} if kind == "SpliceAI" else {"max_abs_delta", "D500_max_abs_delta"})
    require_columns(path, rows, required)
    by_uid: dict[str, dict[str, object]] = {}
    for row in rows:
        uid = row["variant_uid"].strip()
        if not uid or uid in by_uid:
            raise ComparisonError(f"{kind} duplicate/blank variant_uid: {uid!r}")
        if row["status"].strip() != "SCORED" or row["D500_status"].strip() != "SCORED":
            raise ComparisonError(f"{kind} unscored predictor row: {uid}")
        primary_column = "ds_max" if kind == "SpliceAI" else "max_abs_delta"
        sensitivity_column = "D500_ds_max" if kind == "SpliceAI" else "D500_max_abs_delta"
        primary_text = row[primary_column].strip()
        sensitivity_text = row[sensitivity_column].strip()
        by_uid[uid] = {
            "primary": as_fraction(primary_text, f"{kind} D50 {uid}"),
            "sensitivity": as_fraction(sensitivity_text, f"{kind} D500 {uid}"),
            "primary_text": primary_text,
            "sensitivity_text": sensitivity_text,
        }
    if len(by_uid) != 2374:
        raise ComparisonError(f"{kind} predictor row count conflict: {len(by_uid)} != 2374")
    return by_uid


def score_values_for_sets(
    grouped: Mapping[str, Sequence[Mapping[str, str]]],
    predictors: Mapping[str, Mapping[str, object]],
    score_label: str,
) -> tuple[dict[str, list[Fraction]], dict[str, str], dict[str, list[str]]]:
    values: dict[str, list[Fraction]] = {}
    score_text: dict[str, str] = {}
    control_text: dict[str, list[str]] = {}
    for case_key, rows in grouped.items():
        if not rows:
            raise ComparisonError(f"empty matched set: {case_key}")
        case_uid = rows[0]["case_variant_uid"]
        if case_uid not in predictors:
            raise ComparisonError(f"missing {score_label} score for case: {case_uid}")
        case_score = predictors[case_uid][score_label]
        if not isinstance(case_score, Fraction):
            raise ComparisonError(f"internal score type conflict: {score_label}")
        controls: list[Fraction] = []
        texts: list[str] = []
        for row in rows:
            uid = row["control_variant_uid"]
            if uid not in predictors:
                raise ComparisonError(f"missing {score_label} score for control: {uid}")
            score = predictors[uid][score_label]
            if not isinstance(score, Fraction):
                raise ComparisonError(f"internal control score type conflict: {score_label}")
            controls.append(score)
            texts.append(str(predictors[uid][f"{score_label}_text"]))
        values[case_key] = [case_score, *controls]
        score_text[case_key] = str(predictors[case_uid][f"{score_label}_text"])
        control_text[case_key] = texts
    return values, score_text, control_text


def median_fraction(values: Sequence[Fraction]) -> Fraction:
    if not values:
        raise ComparisonError("median of empty score list")
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def set_theta(case_score: Fraction, controls: Sequence[Fraction]) -> Fraction:
    if not controls:
        raise ComparisonError("matched set has zero controls")
    greater = sum(1 for value in controls if case_score > value)
    ties = sum(1 for value in controls if case_score == value)
    return (Fraction(greater) + Fraction(ties, 2)) / len(controls)


def observed_summary(values: Mapping[str, Sequence[Fraction]], selected_cases: Sequence[str]) -> dict[str, object]:
    if not selected_cases:
        raise ComparisonError("no matched sets remain for analysis")
    set_thetas: list[Fraction] = []
    mean_differences: list[Fraction] = []
    median_differences: list[Fraction] = []
    pair_count = 0
    set_details: dict[str, dict[str, object]] = {}
    for key in selected_cases:
        all_values = values[key]
        case_score = all_values[0]
        controls = all_values[1:]
        theta = set_theta(case_score, controls)
        mean_control = sum(controls, Fraction(0)) / len(controls)
        median_control = median_fraction(controls)
        mean_diff = case_score - mean_control
        median_diff = case_score - median_control
        set_thetas.append(theta)
        mean_differences.append(mean_diff)
        median_differences.append(median_diff)
        pair_count += len(controls)
        set_details[key] = {
            "K": len(controls),
            "theta": theta,
            "case_score": case_score,
            "mean_control": mean_control,
            "median_control": median_control,
            "case_minus_mean_controls": mean_diff,
            "case_minus_median_controls": median_diff,
        }
    overall_theta = sum(set_thetas, Fraction(0)) / len(set_thetas)
    mean_delta = sum(mean_differences, Fraction(0)) / len(mean_differences)
    median_delta = median_fraction(mean_differences)
    return {
        "n_sets": len(selected_cases),
        "pair_count": pair_count,
        "theta": overall_theta,
        "rank_biserial": 2 * overall_theta - 1,
        "mean_case_minus_mean_controls": mean_delta,
        "median_case_minus_median_controls": median_delta,
        "set_thetas": set_thetas,
        "set_details": set_details,
    }


def pseudo_case_distribution(scores: Sequence[Fraction]) -> Counter[Fraction]:
    if len(scores) < 2:
        raise ComparisonError("permutation block requires K+1 >= 2 scores")
    distribution: Counter[Fraction] = Counter()
    for pseudo_index, pseudo_score in enumerate(scores):
        controls = [score for index, score in enumerate(scores) if index != pseudo_index]
        distribution[set_theta(pseudo_score, controls)] += 1
    return distribution


def convolve_distributions(distributions: Sequence[Counter[Fraction]]) -> Counter[Fraction]:
    total: Counter[Fraction] = Counter({Fraction(0): 1})
    for distribution in distributions:
        updated: Counter[Fraction] = Counter()
        for left, left_weight in total.items():
            for right, right_weight in distribution.items():
                updated[left + right] += left_weight * right_weight
        total = updated
    return total


def exact_permutation(values: Mapping[str, Sequence[Fraction]], selected_cases: Sequence[str]) -> dict[str, object]:
    if not selected_cases:
        raise ComparisonError("cannot permute zero matched sets")
    distributions = [pseudo_case_distribution(values[key]) for key in selected_cases]
    total_distribution = convolve_distributions(distributions)
    observed = observed_summary(values, selected_cases)
    observed_sum = observed["theta"] * len(selected_cases)
    numerator = sum(weight for total, weight in total_distribution.items() if total >= observed_sum)
    denominator = math.prod(len(values[key]) for key in selected_cases)
    if sum(total_distribution.values()) != denominator:
        raise ComparisonError("DP permutation mass does not equal product(K+1)")
    p_value = Fraction(numerator, denominator)
    return {
        "exact_p": p_value,
        "exact_numerator": numerator,
        "exact_denominator": denominator,
        "dp_states_by_set": [len(distribution) for distribution in distributions],
        "dp_final_states": len(total_distribution),
        "total_combinations": denominator,
        "observed_sum_theta": observed_sum,
    }


def threshold_counts(values: Mapping[str, Sequence[Fraction]], selected_cases: Sequence[str]) -> dict[str, int | float]:
    cases = [values[key][0] for key in selected_cases]
    controls = [score for key in selected_cases for score in values[key][1:]]
    result: dict[str, int | float] = {
        "case_n": len(cases),
        "control_n": len(controls),
    }
    for threshold in THRESHOLDS:
        label = format(threshold.numerator / threshold.denominator, "g").replace(".", "_")
        case_n = sum(score >= threshold for score in cases)
        control_n = sum(score >= threshold for score in controls)
        result[f"case_ge_{label}_n"] = int(case_n)
        result[f"case_ge_{label}_prop"] = case_n / len(cases) if cases else math.nan
        result[f"control_ge_{label}_n"] = int(control_n)
        result[f"control_ge_{label}_prop"] = control_n / len(controls) if controls else math.nan
    return result


def holm_adjust(p_values: Mapping[str, Fraction]) -> dict[str, Fraction]:
    if not p_values:
        return {}
    ordered = sorted(p_values.items(), key=lambda item: (item[1], item[0]))
    m = len(ordered)
    adjusted: dict[str, Fraction] = {}
    running = Fraction(0)
    # Holm step-down adjustment: after sorting p-values increasingly,
    # adjusted p_(i) is the running maximum of (m-i) * p_(i).
    # Iterating backwards would incorrectly propagate the largest raw
    # p-value to every smaller p-value.
    for index in range(m):
        label, p_value = ordered[index]
        candidate = min(Fraction(1), p_value * (m - index))
        running = max(running, candidate)
        adjusted[label] = running
    return adjusted


def loo_summary(values: Mapping[str, Sequence[Fraction]], case_keys: Sequence[str], full_theta: Fraction) -> dict[str, object]:
    records: list[tuple[str, dict[str, object]]] = []
    for omitted in case_keys:
        remaining = [key for key in case_keys if key != omitted]
        if not remaining:
            raise ComparisonError("LOO would leave zero matched sets")
        observed = observed_summary(values, remaining)
        permutation = exact_permutation(values, remaining)
        records.append((omitted, {"observed": observed, "permutation": permutation}))
    theta_values = [entry["observed"]["theta"] for _, entry in records]
    p_values = [entry["permutation"]["exact_p"] for _, entry in records]
    influence = sorted(
        records,
        key=lambda item: (-abs(item[1]["observed"]["theta"] - full_theta), item[0]),
    )[0]
    return {
        "records": records,
        "theta_min": min(theta_values),
        "theta_max": max(theta_values),
        "p_min": min(p_values),
        "p_max": max(p_values),
        "most_influential_case": influence[0],
        "most_influential_abs_theta_change": abs(influence[1]["observed"]["theta"] - full_theta),
    }


def score_output_rows(
    grouped: Mapping[str, Sequence[Mapping[str, str]]],
    predictor_maps: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case_key in sorted(grouped):
        matches = grouped[case_key]
        case_uid = matches[0]["case_variant_uid"]
        for match in matches:
            control_uid = match["control_variant_uid"]
            row: dict[str, object] = {
                "case_key": case_key,
                "case_histology": match["case_histology"],
                "case_variant_uid": case_uid,
                "case_gene": match["case_gene"],
                "control_key": match["control_key"],
                "control_variant_uid": control_uid,
                "control_gene": match["control_gene"],
                "match_rank": match["match_rank"],
                "same_gene": match["same_gene"],
                "set_K": len(matches),
            }
            for predictor in ("SpliceAI", "Pangolin"):
                for scope, suffix in (("primary", "D50"), ("sensitivity", "D500")):
                    row[f"{predictor.lower()}_{suffix}_case_score"] = predictor_maps[predictor][case_uid][f"{scope}_text"]
                    row[f"{predictor.lower()}_{suffix}_control_score"] = predictor_maps[predictor][control_uid][f"{scope}_text"]
            rows.append(row)
    return rows


SCORE_OUTPUT_HEADER = (
    "case_key",
    "case_histology",
    "case_variant_uid",
    "case_gene",
    "control_key",
    "control_variant_uid",
    "control_gene",
    "match_rank",
    "same_gene",
    "set_K",
    "spliceai_D50_case_score",
    "spliceai_D50_control_score",
    "spliceai_D500_case_score",
    "spliceai_D500_control_score",
    "pangolin_D50_case_score",
    "pangolin_D50_control_score",
    "pangolin_D500_case_score",
    "pangolin_D500_control_score",
)


STATS_HEADER = (
    "predictor",
    "score_scope",
    "analysis_role",
    "status",
    "n_sets",
    "n_pairs",
    "theta",
    "rank_biserial",
    "mean_case_minus_mean_controls",
    "median_case_minus_median_controls",
    "exact_p",
    "exact_p_numerator",
    "exact_p_denominator",
    "holm_p",
    "holm_p_numerator",
    "holm_p_denominator",
    "same_gene_excluded_case_keys",
    "loo_theta_min",
    "loo_theta_max",
    "loo_exact_p_min",
    "loo_exact_p_min_numerator",
    "loo_exact_p_min_denominator",
    "loo_exact_p_max",
    "loo_exact_p_max_numerator",
    "loo_exact_p_max_denominator",
    "loo_most_influential_case",
    "loo_most_influential_abs_theta_change",
    "case_n",
    "control_n",
    "case_ge_0_2_n",
    "case_ge_0_2_prop",
    "control_ge_0_2_n",
    "control_ge_0_2_prop",
    "case_ge_0_5_n",
    "case_ge_0_5_prop",
    "control_ge_0_5_n",
    "control_ge_0_5_prop",
    "case_ge_0_8_n",
    "case_ge_0_8_prop",
    "control_ge_0_8_n",
    "control_ge_0_8_prop",
    "reason",
)


def na_row(**values: object) -> dict[str, object]:
    row = {column: "NA" for column in STATS_HEADER}
    row.update(values)
    return row


def stats_row(
    predictor: str,
    score_scope: str,
    analysis_role: str,
    values: Mapping[str, Sequence[Fraction]],
    selected_cases: Sequence[str],
    *,
    holm: Fraction | None = None,
    same_gene_excluded: Sequence[str] = (),
    include_loo: bool = False,
) -> tuple[dict[str, object], dict[str, object]]:
    observed = observed_summary(values, selected_cases)
    permutation = exact_permutation(values, selected_cases)
    thresholds = threshold_counts(values, selected_cases)
    row = na_row(
        predictor=predictor,
        score_scope=score_scope,
        analysis_role=analysis_role,
        status="PASS",
        n_sets=observed["n_sets"],
        n_pairs=observed["pair_count"],
        theta=format_fraction(observed["theta"]),
        rank_biserial=format_fraction(observed["rank_biserial"]),
        mean_case_minus_mean_controls=format_fraction(observed["mean_case_minus_mean_controls"]),
        median_case_minus_median_controls=format_fraction(observed["median_case_minus_median_controls"]),
        exact_p=format_fraction(permutation["exact_p"]),
        exact_p_numerator=permutation["exact_numerator"],
        exact_p_denominator=permutation["exact_denominator"],
        same_gene_excluded_case_keys=",".join(sorted(same_gene_excluded)) if same_gene_excluded else "",
        case_n=thresholds["case_n"],
        control_n=thresholds["control_n"],
        case_ge_0_2_n=thresholds["case_ge_0_2_n"],
        case_ge_0_2_prop=thresholds["case_ge_0_2_prop"],
        control_ge_0_2_n=thresholds["control_ge_0_2_n"],
        control_ge_0_2_prop=thresholds["control_ge_0_2_prop"],
        case_ge_0_5_n=thresholds["case_ge_0_5_n"],
        case_ge_0_5_prop=thresholds["case_ge_0_5_prop"],
        control_ge_0_5_n=thresholds["control_ge_0_5_n"],
        control_ge_0_5_prop=thresholds["control_ge_0_5_prop"],
        case_ge_0_8_n=thresholds["case_ge_0_8_n"],
        case_ge_0_8_prop=thresholds["case_ge_0_8_prop"],
        control_ge_0_8_n=thresholds["control_ge_0_8_n"],
        control_ge_0_8_prop=thresholds["control_ge_0_8_prop"],
        reason="Continuous matched-set comparison; threshold columns are descriptive only.",
    )
    if holm is not None:
        row["holm_p"] = format_fraction(holm)
        row["holm_p_numerator"] = holm.numerator
        row["holm_p_denominator"] = holm.denominator
    loo_payload: dict[str, object] = {}
    if include_loo:
        loo = loo_summary(values, selected_cases, observed["theta"])
        row.update(
            {
                "loo_theta_min": format_fraction(loo["theta_min"]),
                "loo_theta_max": format_fraction(loo["theta_max"]),
                "loo_exact_p_min": format_fraction(loo["p_min"]),
                "loo_exact_p_min_numerator": loo["p_min"].numerator,
                "loo_exact_p_min_denominator": loo["p_min"].denominator,
                "loo_exact_p_max": format_fraction(loo["p_max"]),
                "loo_exact_p_max_numerator": loo["p_max"].numerator,
                "loo_exact_p_max_denominator": loo["p_max"].denominator,
                "loo_most_influential_case": loo["most_influential_case"],
                "loo_most_influential_abs_theta_change": format_fraction(loo["most_influential_abs_theta_change"]),
            }
        )
        loo_payload = loo
    return row, {"observed": observed, "permutation": permutation, "thresholds": thresholds, "loo": loo_payload}


def exhaustive_permutation(values: Mapping[str, Sequence[Fraction]], selected_cases: Sequence[str]) -> Counter[Fraction]:
    total = Counter()
    choices = [range(len(values[key])) for key in selected_cases]
    for selected_indices in itertools.product(*choices):
        theta_sum = Fraction(0)
        for key, index in zip(selected_cases, selected_indices):
            scores = values[key]
            pseudo = scores[index]
            controls = [score for j, score in enumerate(scores) if j != index]
            theta_sum += set_theta(pseudo, controls)
        total[theta_sum] += 1
    return total


def run_self_test() -> dict[str, object]:
    """Test ties, Fraction DP, separation, Holm and LOO without project files."""

    tie_values = {
        "A": [Fraction(1), Fraction(1), Fraction(0)],
        "B": [Fraction(2), Fraction(2), Fraction(2)],
    }
    dp = convolve_distributions([pseudo_case_distribution(tie_values[key]) for key in ("A", "B")])
    brute = exhaustive_permutation(tie_values, ("A", "B"))
    if dp != brute:
        raise ComparisonError(f"Fraction DP differs from exhaustive enumeration: {dp} != {brute}")
    all_ties = {"A": [Fraction(1), Fraction(1), Fraction(1)], "B": [Fraction(2), Fraction(2)]}
    all_ties_p = exact_permutation(all_ties, ("A", "B"))
    if all_ties_p["exact_p"] != 1:
        raise ComparisonError(f"all-ties p is not 1: {all_ties_p}")
    separated = {"A": [Fraction(1), Fraction(0), Fraction(0), Fraction(0)]}
    separated_p = exact_permutation(separated, ("A",))
    if separated_p["exact_p"] != Fraction(1, 4):
        raise ComparisonError(f"strict separation p conflict: {separated_p}")
    holm = holm_adjust({"a": Fraction(1, 4), "b": Fraction(1, 2)})
    if holm != {"a": Fraction(1, 2), "b": Fraction(1, 2)}:
        raise ComparisonError(f"Holm self-test conflict: {holm}")
    holm_separated = holm_adjust({"a": Fraction(1, 100), "b": Fraction(1, 5)})
    if holm_separated != {"a": Fraction(1, 50), "b": Fraction(1, 5)}:
        raise ComparisonError(f"Holm separated-p self-test conflict: {holm_separated}")
    loo_values = {
        "A": [Fraction(3), Fraction(0), Fraction(0)],
        "B": [Fraction(2), Fraction(0), Fraction(0)],
        "C": [Fraction(0), Fraction(0), Fraction(0)],
    }
    loo = loo_summary(loo_values, ("A", "B", "C"), observed_summary(loo_values, ("A", "B", "C"))["theta"])
    if len(loo["records"]) != 3 or loo["theta_min"] > loo["theta_max"] or loo["p_min"] > loo["p_max"]:
        raise ComparisonError(f"LOO self-test conflict: {loo}")
    return {
        "status": "PASS",
        "fraction_dp_matches_exhaustive": True,
        "ties_midrank": True,
        "all_ties_p": "1",
        "strict_separation_p": "1/4",
        "holm": {key: fraction_text(value) for key, value in sorted(holm.items())},
        "loo_case_count": len(loo["records"]),
    }


def project_paths(root: Path) -> dict[str, Path]:
    return {
        "match": root / "results/tables/V23_M6_MATCHED_CONTROLS.tsv.gz",
        "match_audit": root / "results/audit/m4_m5_phase1/V23_M6_MATCH_FREEZE_AUDIT.json",
        "frozen9": root / "results/audit/V23_POST_M2_M3_FROZEN9_KEYS.tsv",
        "spliceai": root / "results/tables/V23_M4_SPLICEAI_ALL_VARIANTS.tsv.gz",
        "pangolin": root / "results/tables/V23_M5_PANGOLIN_ALL_VARIANTS.tsv.gz",
        "predictor_audit": root / "results/audit/m4_m5_phase1/V23_PREDICTOR_RAW_FREEZE_AUDIT.json",
        "scores": root / "results/stats/m6_predictor/V23_M6_MATCHED_SET_SCORES.tsv.gz",
        "stats": root / "results/stats/V23_M6_PREDICTOR_COMPARISON.tsv",
        "audit": root / "results/audit/m4_m5_phase1/V23_M6_PREDICTOR_COMPARISON_AUDIT.json",
        "manifest": root / "results/audit/m4_m5_phase1/V23_M6_PREDICTOR_COMPARISON_SHA256.tsv",
        "script": root / "scripts/32_analyze_m6_predictor_comparison.py",
    }


def verify_paths(paths: Mapping[str, Path]) -> None:
    required = ("match", "match_audit", "frozen9", "spliceai", "pangolin", "predictor_audit", "script")
    missing = [label for label in required if not paths[label].is_file()]
    if missing:
        raise ComparisonError(f"missing required input(s): {missing}")
    for relpath, expected in EXPECTED_STEP7.items():
        path = paths["script"].parents[1] / relpath
        observed = sha256_file(path)
        if observed != expected:
            raise ComparisonError(f"frozen Step7 hash conflict for {relpath}: {observed} != {expected}")
    if any(paths[label].exists() for label in ("scores", "stats", "audit", "manifest")):
        existing = [str(paths[label]) for label in ("scores", "stats", "audit", "manifest") if paths[label].exists()]
        raise ComparisonError(f"formal output already exists; refusing overwrite: {existing}")


def build_stats(
    grouped: Mapping[str, Sequence[Mapping[str, str]]],
    predictor_maps: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    case_keys = tuple(sorted(grouped))
    values_by_metric: dict[tuple[str, str], dict[str, list[Fraction]]] = {}
    for predictor in ("SpliceAI", "Pangolin"):
        for scope in ("primary", "sensitivity"):
            values, _, _ = score_values_for_sets(grouped, predictor_maps[predictor], scope)
            values_by_metric[(predictor, scope)] = values
    primary_p = {
        predictor: exact_permutation(values_by_metric[(predictor, "primary")], case_keys)["exact_p"]
        for predictor in ("SpliceAI", "Pangolin")
    }
    holm = holm_adjust(primary_p)
    rows: list[dict[str, object]] = []
    details: dict[str, object] = {"primary": {}, "D500": {}, "same_gene": {}, "LOO": {}, "DP": {}}
    for predictor in ("SpliceAI", "Pangolin"):
        row, detail = stats_row(
            predictor,
            "D50",
            "PRIMARY_D50",
            values_by_metric[(predictor, "primary")],
            case_keys,
            holm=holm[predictor],
            include_loo=True,
        )
        rows.append(row)
        details["primary"][predictor] = detail
        details["LOO"][predictor] = detail["loo"]
        details["DP"][f"{predictor}_D50"] = detail["permutation"]
    for predictor in ("SpliceAI", "Pangolin"):
        row, detail = stats_row(
            predictor,
            "D500",
            "SENSITIVITY_D500",
            values_by_metric[(predictor, "sensitivity")],
            case_keys,
        )
        rows.append(row)
        details["D500"][predictor] = detail
        details["DP"][f"{predictor}_D500"] = detail["permutation"]
    same_gene_grouped: dict[str, list[dict[str, str]]] = {}
    excluded: list[str] = []
    for key in case_keys:
        selected = [row for row in grouped[key] if parse_bool(row["same_gene"], f"same_gene {key}")]
        if selected:
            same_gene_grouped[key] = selected
        else:
            excluded.append(key)
    same_gene_values_by_metric: dict[tuple[str, str], dict[str, list[Fraction]]] = {}
    same_gene_keys = tuple(sorted(same_gene_grouped))
    if not same_gene_keys:
        raise ComparisonError("same-gene-only sensitivity has zero matched sets")
    for predictor in ("SpliceAI", "Pangolin"):
        values, _, _ = score_values_for_sets(same_gene_grouped, predictor_maps[predictor], "primary")
        same_gene_values_by_metric[(predictor, "primary")] = values
        row, detail = stats_row(
            predictor,
            "D50",
            "SENSITIVITY_SAME_GENE_D50",
            values,
            same_gene_keys,
            same_gene_excluded=excluded,
        )
        rows.append(row)
        details["same_gene"][predictor] = {"n_sets": len(same_gene_keys), "excluded_case_keys": excluded, **detail}
        details["DP"][f"{predictor}_same_gene_D50"] = detail["permutation"]
    rows.append(
        na_row(
            predictor="GEOMETRY_CONCORDANCE",
            score_scope="NA",
            analysis_role="PRIMARY_D50",
            status="NOT_ESTIMABLE_UNDER_PRIMARY_R0_CONTROL_DESIGN",
            reason="Controls are R0 proxy variants without observed altered RNA geometry; no geometry outcome exists for comparison, so effect sizes and p values are NA.",
        )
    )
    details["geometry"] = {
        "status": "NOT_ESTIMABLE_UNDER_PRIMARY_R0_CONTROL_DESIGN",
        "effect_and_p_values": "NA",
        "reason": "Controls are R0 proxy variants without observed altered RNA geometry; candidate junctions are not used and controls are not coded as zero.",
    }
    return rows, details


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=None, help="v2.3 project root")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args(argv)
    if args.self_test:
        print(stable_json(run_self_test()), end="")
        return 0
    root = (args.root or Path(__file__).resolve().parents[1]).resolve()
    paths = project_paths(root)
    verify_paths(paths)
    cases = load_frozen_cases(paths["frozen9"])
    grouped, match_counts = load_match_table(paths["match"], cases)
    match_audit = load_match_audit(paths["match_audit"])
    predictor_maps = {
        "SpliceAI": load_predictors(paths["spliceai"], "SpliceAI"),
        "Pangolin": load_predictors(paths["pangolin"], "Pangolin"),
    }
    relevant_uids = {row["case_variant_uid"] for rows in grouped.values() for row in rows} | {row["control_variant_uid"] for rows in grouped.values() for row in rows}
    for predictor, mapping in predictor_maps.items():
        missing = sorted(relevant_uids.difference(mapping))
        if missing:
            raise ComparisonError(f"{predictor} missing frozen case/control scores: {missing}")
    score_rows = score_output_rows(grouped, predictor_maps)
    stats_rows, details = build_stats(grouped, predictor_maps)
    self_test = run_self_test()
    predictor_fields = [
        "SpliceAI.ds_max",
        "SpliceAI.D500_ds_max",
        "Pangolin.max_abs_delta",
        "Pangolin.D500_max_abs_delta",
        "status",
        "D500_status",
    ]
    input_paths = {
        "match": paths["match"],
        "match_audit": paths["match_audit"],
        "frozen9": paths["frozen9"],
        "spliceai": paths["spliceai"],
        "pangolin": paths["pangolin"],
        "predictor_audit": paths["predictor_audit"],
        "script": paths["script"],
    }
    input_sha = {label: sha256_file(path) for label, path in input_paths.items()}
    audit = {
        "status": "PASS",
        "analysis": "M6_PREDICTOR_COMPARISON_PHASE2",
        "outcome_fields_used": [],
        "RNA_outcome_fields_used": [],
        "predictor_fields_used": predictor_fields,
        "input_sha256": input_sha,
        "design": {
            "analysis_unit": "9 non-overlapping histology x variant_uid matched sets; case equal-weighted",
            "case_count": EXPECTED_CASE_COUNT,
            "case_K": {key: EXPECTED_CASE_K[key] for key in sorted(EXPECTED_CASE_K)},
            "primary": {
                "SpliceAI": "ds_max D50",
                "Pangolin": "max_abs_delta D50",
            },
            "sensitivity": {
                "D500": {"SpliceAI": "D500_ds_max", "Pangolin": "D500_max_abs_delta"},
                "same_gene_only": "retain only same_gene TRUE; exclude sets with zero same-gene controls",
                "leave_one_case_out": "primary D50; report theta and exact p ranges; influence is largest absolute change in overall theta",
            },
            "theta": "per-set (#case>control + 0.5*#tie)/K; overall theta is unweighted mean across sets; rank_biserial=2*theta-1",
            "continuous_effects": "mean across sets of (case - mean controls), and median across sets of (case - median controls)",
            "permutation": "one pseudo-case uniformly chosen from K+1 units in each set; midrank ties; Counter/Fraction convolution; one-sided H1 case higher; no pseudocount",
            "primary_multiplicity": "Holm adjustment over exactly the two primary D50 p values; D500 and other sensitivity rows excluded",
            "descriptive_thresholds": ">=0.2, >=0.5, >=0.8 case/control proportions only; never replace continuous test",
        },
        "matched_sets": match_counts,
        "match_audit_status": match_audit["status"],
        "score_integrity": {
            "predictor_variant_rows": {predictor: len(mapping) for predictor, mapping in predictor_maps.items()},
            "all_rows_d50_d500_scored": True,
            "relevant_case_control_uids": len(relevant_uids),
            "relevant_case_control_uids_scored": len(relevant_uids),
            "global_unique_variant_uid_per_predictor": True,
        },
        "results": details,
        "self_test": self_test,
        "geometry": details["geometry"],
        "outputs": {},
    }
    deterministic_gzip_write(paths["scores"], SCORE_OUTPUT_HEADER, score_rows)
    deterministic_tsv_write(paths["stats"], STATS_HEADER, stats_rows)
    output_targets = [
        ("output", "V23_M6_MATCHED_SET_SCORES.tsv.gz", paths["scores"]),
        ("output", "V23_M6_PREDICTOR_COMPARISON.tsv", paths["stats"]),
    ]
    audit["outputs"] = {
        label: {"rows": len(score_rows) if label.endswith("SCORES.tsv.gz") else len(stats_rows), "sha256": sha256_file(path), "size_bytes": file_size(path)}
        for _, label, path in output_targets
    }
    paths["audit"].parent.mkdir(parents=True, exist_ok=True)
    paths["audit"].write_text(stable_json(audit), encoding="utf-8")
    output_targets.append(("output", "V23_M6_PREDICTOR_COMPARISON_AUDIT.json", paths["audit"]))
    manifest_rows: list[dict[str, object]] = []
    for label, path in input_paths.items():
        manifest_rows.append({"kind": "input", "label": label, "path": str(path), "sha256": sha256_file(path), "size_bytes": file_size(path)})
    for _, label, path in output_targets:
        manifest_rows.append({"kind": "output", "label": label, "path": str(path), "sha256": sha256_file(path), "size_bytes": file_size(path)})
    manifest_rows.sort(key=lambda row: (row["kind"], row["label"], row["path"]))
    deterministic_tsv_write(paths["manifest"], ("kind", "label", "path", "sha256", "size_bytes"), manifest_rows)
    summary = {
        "status": "PASS",
        "stats": str(paths["stats"]),
        "scores": str(paths["scores"]),
        "audit": str(paths["audit"]),
        "manifest": str(paths["manifest"]),
        "primary": {
            predictor: {
                "theta": details["primary"][predictor]["observed"]["theta"],
                "rank_biserial": details["primary"][predictor]["observed"]["rank_biserial"],
                "exact_p": details["primary"][predictor]["permutation"]["exact_p"],
                "holm_p": holm_adjust({p: details["primary"][p]["permutation"]["exact_p"] for p in ("SpliceAI", "Pangolin")})[predictor],
            }
            for predictor in ("SpliceAI", "Pangolin")
        },
        "sha256": {
            "stats": sha256_file(paths["stats"]),
            "scores": sha256_file(paths["scores"]),
            "audit": sha256_file(paths["audit"]),
            "manifest": sha256_file(paths["manifest"]),
        },
    }
    print(stable_json(summary), end="")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ComparisonError as exc:
        print(f"M6_PREDICTOR_COMPARISON_FAIL: {exc}", file=sys.stderr)
        raise SystemExit(2)
