#!/usr/bin/env python3
"""Outcome-blind Step 8 RNA--predictor geometry classification for frozen9.

This script joins the already frozen predictor tables to the nine post-M2/M3
events and classifies only splice-site geometry.  It intentionally does not
read patient outcomes, usage, M2, M3, or any other effect/result field.  The
raw junction endpoints are intron interval endpoints; they are converted to
predictor splice-site coordinates before any coordinate comparison:

    plus:  left donor = left - 1; right acceptor = right + 1
    minus: left acceptor = left - 1; right donor = right + 1

The 0.2 threshold is a descriptive working rule, not a validated clinical
cutoff.  Continuous scores are retained as the primary predictor result and
0.5/0.8 bins are descriptive only.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Iterable


RULE_VERSION = "V23_STEP8_GEOMETRY_RULES_v2.3.1"
LOW_THRESHOLD = 0.2
EXPECTED_EVENTS = 9
EXPECTED_PREDICTOR_VARIANTS = 2374
EXPECTED_STEP7 = {
    "scripts/29_freeze_predictor_outputs_qc.py": "5a02495034208e1a22ba280244bfe0a3702c8f7715ee085847f44b50a214ad68",
    "results/audit/V23_PREDICTOR_QC.tsv": "41f4c29dcae7fa8fec01bf4a272fd269963c08af9689ccab15d1a3e908286b07",
    "results/audit/m4_m5_phase1/V23_PREDICTOR_RAW_FREEZE_AUDIT.json": "cea4ff166808319d8b2998d629a904a37b6f6874a50adab72fb613d6922447d9",
    "results/audit/m4_m5_phase1/V23_PREDICTOR_RAW_FREEZE_SHA256.tsv": "c9e103ea0e9aae26f0973d40f90588f86ec8943afaa5b55a9254e036aaae4791",
    "results/tables/V23_M4_SPLICEAI_ALL_VARIANTS.tsv.gz": "525e896f358ee1d6c7a2e83ae54207c353c31a0e6861968af4dfb420358f68f6",
    "results/tables/V23_M5_PANGOLIN_ALL_VARIANTS.tsv.gz": "1c9f75d93c48e9d67f031ceeb5c571bf3ca2c916fd3e6fd8a3936ee31797dc2d",
}

GEOMETRY_CLASSES = {
    "EXACT_GEOMETRY_CONCORDANT",
    "NEAR_GEOMETRY_CONCORDANT",
    "CANONICAL_SITE_LOSS_CONCORDANT",
    "ALTERNATIVE_SITE_GAIN_COMPATIBLE",
}
VALID_CLASSES = GEOMETRY_CLASSES | {
    "SPLICE_EFFECT_ONLY",
    "PREDICTOR_LOW",
    "UNSCORED",
}

FROZEN_REQUIRED = ("histology", "variant_uid", "gene", "junction_id", "cluster_id", "selected_transcript")
EVIDENCE_ALLOWED = ("variant_uid", "gene", "junction_id", "v23_altered_splice_type")
COMPETITOR_ALLOWED = (
    "variant_uid",
    "gene",
    "junction_id",
    "primary_competitor_junction_id",
    "secondary_competitor_junction_id",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative_path(root: Path, path: Path) -> str:
    return path.resolve().relative_to(root.resolve()).as_posix()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"empty TSV header: {path}")
        return list(reader)


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"empty gzip TSV header: {path}")
        return list(reader)


def require_columns(rows: list[dict[str, str]], columns: Iterable[str], path: Path) -> None:
    # DictReader stores the header on the first row only; an empty table is
    # rejected by callers, while this check protects all non-empty inputs.
    if not rows:
        raise ValueError(f"empty input table: {path}")
    missing = sorted(set(columns) - set(rows[0]))
    if missing:
        raise ValueError(f"{path} missing columns: {missing}")


def finite_float(text: str, label: str, *, allow_blank: bool = False) -> float | None:
    if text == "" and allow_blank:
        return None
    try:
        value = float(text)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid {label}: {text!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"non-finite {label}: {text!r}")
    return value


def int_field(text: str, label: str, *, allow_blank: bool = False) -> int | None:
    if text == "" and allow_blank:
        return None
    try:
        return int(text)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid {label}: {text!r}") from exc


def fmt_float(value: float | None) -> str:
    return "" if value is None else format(value, ".17g")


def score_bin(value: float | None) -> str:
    if value is None:
        return "UNSCORED"
    if value < 0.2:
        return "LOW_LT_0.2"
    if value < 0.5:
        return "DESCRIPTIVE_0.2_TO_LT_0.5"
    if value < 0.8:
        return "DESCRIPTIVE_0.5_TO_LT_0.8"
    return "DESCRIPTIVE_GE_0.8"


JUNCTION_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<left>[0-9]+)-(?P<right>[0-9]+):(?P<strand>[+-])$")
VARIANT_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<pos>[0-9]+):(?P<ref>[^:]+):(?P<alt>[^:]+)$")


def parse_variant_uid(value: str) -> tuple[str, int]:
    match = VARIANT_RE.fullmatch(value)
    if not match:
        raise ValueError(f"invalid variant_uid: {value!r}")
    return match.group("chrom"), int(match.group("pos"))


def splice_site_coordinate(side: str, raw: int) -> int:
    if side not in {"left", "right"}:
        raise ValueError(f"invalid junction side: {side!r}")
    # Both left and right are intron interval boundaries in the frozen tables.
    # Their adjacent exonic splice-site nucleotide is one base outside the
    # interval, independent of strand.
    return raw - 1 if side == "left" else raw + 1


def parse_junction(value: str) -> dict[str, Any]:
    match = JUNCTION_RE.fullmatch(value)
    if not match:
        raise ValueError(f"invalid junction_id: {value!r}")
    left = int(match.group("left"))
    right = int(match.group("right"))
    if left >= right:
        raise ValueError(f"junction endpoints are not increasing: {value!r}")
    strand = match.group("strand")
    if strand == "+":
        roles = {"left": "donor", "right": "acceptor"}
    else:
        roles = {"left": "acceptor", "right": "donor"}
    endpoints: list[dict[str, Any]] = []
    for side, raw in (("left", left), ("right", right)):
        endpoints.append(
            {
                "side": side,
                "raw": raw,
                "role": roles[side],
                "site": splice_site_coordinate(side, raw),
            }
        )
    return {
        "chrom": match.group("chrom"),
        "left": left,
        "right": right,
        "strand": strand,
        "endpoints": endpoints,
    }


def endpoint_map(junction: dict[str, Any], label: str) -> dict[int, dict[str, Any]]:
    result: dict[int, dict[str, Any]] = {}
    for endpoint in junction["endpoints"]:
        raw = endpoint["raw"]
        if raw in result:
            raise ValueError(f"duplicate endpoint in {label}: {raw}")
        result[raw] = endpoint
    return result


def join_key_rows(rows: list[dict[str, str]], label: str) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        key = row.get("variant_uid", "")
        if not key:
            raise ValueError(f"{label} contains blank variant_uid")
        if key in result:
            raise ValueError(f"{label} has duplicate variant_uid: {key}")
        result[key] = row
    return result


def canonical_endpoint_map(junctions: list[tuple[str, dict[str, Any]]]) -> dict[int, dict[str, Any]]:
    result: dict[int, dict[str, Any]] = {}
    for label, junction in junctions:
        for endpoint in junction["endpoints"]:
            raw = endpoint["raw"]
            previous = result.get(raw)
            if previous is not None and (
                previous["role"] != endpoint["role"] or previous["site"] != endpoint["site"]
            ):
                raise ValueError(f"canonical endpoint role/site conflict at raw coordinate {raw}")
            result[raw] = endpoint
    return result


def prepare_context(key: dict[str, str], event: dict[str, str], competitor: dict[str, str]) -> dict[str, Any]:
    if key["variant_uid"] != event["variant_uid"] or key["variant_uid"] != competitor["variant_uid"]:
        raise ValueError("key join mismatch")
    if key["gene"] != event["gene"] or key["gene"] != competitor["gene"]:
        raise ValueError(f"gene mismatch for {key['variant_uid']}")
    if key["junction_id"] != event["junction_id"] or key["junction_id"] != competitor["junction_id"]:
        raise ValueError(f"junction mismatch for {key['variant_uid']}")
    event_type = event["v23_altered_splice_type"]
    if event_type not in {"EXON_SKIP_1", "CRYPTIC_ACCEPTOR_TRUNCATION_13NT"}:
        raise ValueError(f"unsupported event type for {key['variant_uid']}: {event_type}")
    observed = parse_junction(key["junction_id"])
    primary_text = competitor["primary_competitor_junction_id"]
    if not primary_text:
        raise ValueError(f"missing primary competitor for {key['variant_uid']}")
    primary = parse_junction(primary_text)
    secondary_text = competitor["secondary_competitor_junction_id"]
    secondary = parse_junction(secondary_text) if secondary_text else None
    for label, other in (("primary", primary), ("secondary", secondary)):
        if other is None:
            continue
        if other["chrom"] != observed["chrom"] or other["strand"] != observed["strand"]:
            raise ValueError(f"{label} competitor chrom/strand mismatch for {key['variant_uid']}")
    if observed["chrom"] != parse_variant_uid(key["variant_uid"])[0]:
        raise ValueError(f"junction and variant chromosome mismatch for {key['variant_uid']}")
    observed_by_raw = endpoint_map(observed, "observed junction")
    primary_by_raw = endpoint_map(primary, "primary competitor")
    canonical = canonical_endpoint_map(
        [("primary", primary)] + ([('secondary', secondary)] if secondary is not None else [])
    )
    context: dict[str, Any] = {
        "key": key,
        "event_type": event_type,
        "observed": observed,
        "primary": primary,
        "secondary": secondary,
        "observed_by_raw": observed_by_raw,
        "primary_by_raw": primary_by_raw,
        "canonical_by_raw": canonical,
        "primary_text": primary_text,
        "secondary_text": secondary_text,
    }
    if event_type == "EXON_SKIP_1":
        if secondary is None:
            raise ValueError(f"EXON_SKIP_1 requires secondary competitor: {key['variant_uid']}")
        lost_raw = sorted(set(canonical) - set(observed_by_raw))
        if len(lost_raw) != 2:
            raise ValueError(f"EXON_SKIP_1 must derive exactly two lost raw endpoints for {key['variant_uid']}: {lost_raw}")
        lost = [canonical[raw] for raw in lost_raw]
        roles = Counter(endpoint["role"] for endpoint in lost)
        if roles != Counter({"donor": 1, "acceptor": 1}):
            raise ValueError(f"EXON_SKIP_1 lost endpoints must contain one donor and one acceptor: {key['variant_uid']}")
        context["lost_by_role"] = {endpoint["role"]: endpoint for endpoint in lost}
    else:
        altered_raw = sorted(set(observed_by_raw) - set(primary_by_raw))
        replaced_raw = sorted(set(primary_by_raw) - set(observed_by_raw))
        if len(altered_raw) != 1 or len(replaced_raw) != 1:
            raise ValueError(
                f"cryptic event must derive one altered and one replaced endpoint for {key['variant_uid']}: "
                f"altered={altered_raw}, replaced={replaced_raw}"
            )
        altered = observed_by_raw[altered_raw[0]]
        replaced = primary_by_raw[replaced_raw[0]]
        if altered["role"] != "acceptor" or replaced["role"] != "acceptor":
            raise ValueError(f"cryptic acceptor event endpoint roles are not acceptors: {key['variant_uid']}")
        context["altered"] = altered
        context["replaced"] = replaced
    return context


def validate_predictor_target(row: dict[str, str], context: dict[str, Any], model: str) -> None:
    key = context["key"]
    expected_chrom, expected_pos = parse_variant_uid(key["variant_uid"])
    if row.get("variant_uid") != key["variant_uid"]:
        raise ValueError(f"{model} variant mismatch for {key['variant_uid']}")
    if row.get("gene") != key["gene"]:
        raise ValueError(f"{model} gene mismatch for {key['variant_uid']}")
    if row.get("chrom") != expected_chrom or int_field(row.get("pos_1based", ""), f"{model} pos") != expected_pos:
        raise ValueError(f"{model} variant coordinates mismatch for {key['variant_uid']}")
    if row.get("target_transcript") != key["selected_transcript"]:
        raise ValueError(f"{model} transcript mismatch for {key['variant_uid']}")
    if row.get("strand") != context["observed"]["strand"]:
        raise ValueError(f"{model} strand mismatch for {key['variant_uid']}")
    if row.get("outcome_blind_input", "").upper() != "TRUE":
        raise ValueError(f"{model} outcome_blind_input is not TRUE for {key['variant_uid']}")
    if row.get("status") not in {"SCORED", "UNSCORED", "FAILED", "MISSING"}:
        raise ValueError(f"unexpected {model} status for {key['variant_uid']}: {row.get('status')!r}")
    if row.get("D500_status") not in {"SCORED", "UNSCORED", "FAILED", "MISSING"}:
        raise ValueError(f"unexpected {model} D500_status for {key['variant_uid']}: {row.get('D500_status')!r}")


def component_prefix(window: str) -> str:
    return "" if window == "D50" else "D500_"


def parse_spliceai_metrics(row: dict[str, str], window: str) -> dict[str, Any]:
    prefix = component_prefix(window)
    status = row[f"{prefix}status"] if prefix else row["status"]
    if status != "SCORED":
        return {"status": status, "score": None, "components": {}, "coords": {}}
    score_key = f"{prefix}ds_max" if prefix else "ds_max"
    score = finite_float(row[score_key], f"SpliceAI {window} ds_max")
    components: dict[str, float] = {}
    coords: dict[str, int] = {}
    for name in ("AG", "AL", "DG", "DL"):
        score_name = f"{prefix}DS_{name}"
        coord_name = f"{prefix}COORD_{name}"
        value = finite_float(row[score_name], f"SpliceAI {window} DS_{name}")
        coord = int_field(row[coord_name], f"SpliceAI {window} COORD_{name}")
        assert value is not None and coord is not None
        if value < 0:
            raise ValueError(f"negative SpliceAI {window} DS_{name}")
        components[name] = value
        coords[name] = coord
    if score is None or score < 0:
        raise ValueError(f"negative SpliceAI {window} ds_max")
    if not math.isclose(score, max(components.values()), rel_tol=1e-9, abs_tol=1e-9):
        raise ValueError(f"SpliceAI {window} ds_max is inconsistent with component scores")
    type_key = f"{prefix}ds_max_type" if prefix else "ds_max_type"
    coord_key = f"{prefix}ds_max_coordinate" if prefix else "ds_max_coordinate"
    max_type = row[type_key]
    max_coordinate = int_field(row[coord_key], f"SpliceAI {window} ds_max_coordinate")
    if max_type not in {"DS_AG", "DS_AL", "DS_DG", "DS_DL"} or max_coordinate is None:
        raise ValueError(f"invalid SpliceAI {window} max metadata")
    return {
        "status": status,
        "score": score,
        "score_text": row[score_key],
        "components": components,
        "component_text": {name: row[f"{prefix}DS_{name}"] for name in ("AG", "AL", "DG", "DL")},
        "coords": coords,
        "max_type": max_type,
        "max_coordinate": max_coordinate,
        "max_coordinate_text": row[coord_key],
    }


def nearest_to(coords: list[int], targets: list[int]) -> tuple[int | None, int | None]:
    if not coords or not targets:
        return None, None
    distance, coordinate, _target = min(
        (abs(coordinate - target), coordinate, target) for coordinate in coords for target in targets
    )
    return distance, coordinate


def classify_spliceai(row: dict[str, str], context: dict[str, Any], window: str) -> dict[str, Any]:
    metrics = parse_spliceai_metrics(row, window)
    prefix = "spliceai_" + window.lower() + "_"
    result: dict[str, Any] = {
        prefix + "status": metrics["status"],
        prefix + "score": metrics.get("score_text", ""),
        prefix + "score_bin": score_bin(metrics.get("score")),
        prefix + "max_type": metrics.get("max_type", ""),
        prefix + "max_coordinate": str(metrics.get("max_coordinate", "")) if metrics.get("max_coordinate") is not None else "",
        prefix + "relevant_score": "",
        prefix + "relevant_score_bin": "UNSCORED" if metrics.get("score") is None else "",
        prefix + "raw_coordinate_match_class": "UNSCORED" if metrics.get("score") is None else "RAW_NO_COORDINATE_MATCH",
        prefix + "geometry_class": "UNSCORED" if metrics.get("score") is None else "",
        prefix + "geometry_match_component": "",
        prefix + "geometry_match_role": "",
        prefix + "geometry_match_site": "",
        prefix + "geometry_distance_nt": "",
        prefix + "altered_ag_distance_nt": "",
        prefix + "replaced_al_distance_nt": "",
        prefix + "loss_distance_nt": "",
        prefix + "gain_distance_nt": "",
    }
    if metrics.get("score") is None:
        return result
    score = metrics["score"]
    components: dict[str, float] = metrics["components"]
    coords: dict[str, int] = metrics["coords"]
    if context["event_type"] == "CRYPTIC_ACCEPTOR_TRUNCATION_13NT":
        altered_site = context["altered"]["site"]
        replaced_site = context["replaced"]["site"]
        altered_distance = abs(coords["AG"] - altered_site)
        replaced_distance = abs(coords["AL"] - replaced_site)
        result[prefix + "altered_ag_distance_nt"] = str(altered_distance)
        result[prefix + "replaced_al_distance_nt"] = str(replaced_distance)
        if altered_distance == 0:
            raw_class = "RAW_EXACT_COORDINATE_MATCH"
        elif altered_distance <= 2:
            raw_class = "RAW_NEAR_COORDINATE_MATCH_2NT"
        else:
            raw_class = "RAW_NO_COORDINATE_MATCH"
        result[prefix + "raw_coordinate_match_class"] = raw_class
        candidate = components["AG"] if altered_distance <= 2 else None
        if candidate is not None:
            result[prefix + "relevant_score"] = metrics["component_text"]["AG"]
            result[prefix + "relevant_score_bin"] = score_bin(candidate)
        if candidate is not None and candidate >= LOW_THRESHOLD:
            result[prefix + "geometry_class"] = (
                "EXACT_GEOMETRY_CONCORDANT" if altered_distance == 0 else "NEAR_GEOMETRY_CONCORDANT"
            )
            result[prefix + "geometry_match_component"] = "DS_AG"
            result[prefix + "geometry_match_role"] = "acceptor"
            result[prefix + "geometry_match_site"] = str(altered_site)
            result[prefix + "geometry_distance_nt"] = str(altered_distance)
        else:
            result[prefix + "geometry_class"] = "SPLICE_EFFECT_ONLY" if score >= LOW_THRESHOLD else "PREDICTOR_LOW"
        return result

    lost = context["lost_by_role"]
    observed = {endpoint["role"]: endpoint for endpoint in context["observed"]["endpoints"]}
    loss_specs = {"donor": "DL", "acceptor": "AL"}
    gain_specs = {"donor": "DG", "acceptor": "AG"}
    loss_distances: list[int] = []
    gain_distances: list[int] = []
    exact_loss: list[tuple[str, str, int, float]] = []
    exact_gain: list[tuple[str, str, int, float]] = []
    for role in ("donor", "acceptor"):
        component = loss_specs[role]
        target = lost[role]["site"]
        distance = abs(coords[component] - target)
        loss_distances.append(distance)
        if distance == 0:
            exact_loss.append((component, role, target, components[component]))
        component = gain_specs[role]
        target = observed[role]["site"]
        distance = abs(coords[component] - target)
        gain_distances.append(distance)
        if distance == 0:
            exact_gain.append((component, role, target, components[component]))
    result[prefix + "loss_distance_nt"] = str(min(loss_distances))
    result[prefix + "gain_distance_nt"] = str(min(gain_distances))
    raw_parts: list[str] = []
    if exact_loss:
        raw_parts.append("RAW_CANONICAL_LOSS_EXACT")
    if exact_gain:
        raw_parts.append("RAW_ALTERNATIVE_GAIN_EXACT")
    result[prefix + "raw_coordinate_match_class"] = "+".join(raw_parts) if raw_parts else "RAW_NO_COORDINATE_MATCH"
    meaningful_loss = [item for item in exact_loss if item[3] >= LOW_THRESHOLD]
    meaningful_gain = [item for item in exact_gain if item[3] >= LOW_THRESHOLD]
    if meaningful_loss:
        component, role, target, value = meaningful_loss[0]
        result[prefix + "relevant_score"] = metrics["component_text"][component]
        result[prefix + "relevant_score_bin"] = score_bin(value)
        result[prefix + "geometry_class"] = "CANONICAL_SITE_LOSS_CONCORDANT"
        result[prefix + "geometry_match_component"] = "DS_" + component
        result[prefix + "geometry_match_role"] = role
        result[prefix + "geometry_match_site"] = str(target)
        result[prefix + "geometry_distance_nt"] = "0"
    elif meaningful_gain:
        component, role, target, value = meaningful_gain[0]
        result[prefix + "relevant_score"] = metrics["component_text"][component]
        result[prefix + "relevant_score_bin"] = score_bin(value)
        result[prefix + "geometry_class"] = "ALTERNATIVE_SITE_GAIN_COMPATIBLE"
        result[prefix + "geometry_match_component"] = "DS_" + component
        result[prefix + "geometry_match_role"] = role
        result[prefix + "geometry_match_site"] = str(target)
        result[prefix + "geometry_distance_nt"] = "0"
    else:
        # Preserve the relevant score even for a low-score raw coordinate hit;
        # if neither relevant coordinate is exact, the overall score remains
        # the appropriate continuous fallback.
        candidate = exact_loss[0] if exact_loss else (exact_gain[0] if exact_gain else None)
        if candidate is not None:
            result[prefix + "relevant_score"] = metrics["component_text"][candidate[0]]
            result[prefix + "relevant_score_bin"] = score_bin(candidate[3])
        result[prefix + "geometry_class"] = "SPLICE_EFFECT_ONLY" if score >= LOW_THRESHOLD else "PREDICTOR_LOW"
        result[prefix + "geometry_distance_nt"] = str(min(loss_distances + gain_distances))
    return result


def parse_pangolin_metrics(row: dict[str, str], window: str) -> dict[str, Any]:
    prefix = component_prefix(window)
    status = row[f"{prefix}status"] if prefix else row["status"]
    if status != "SCORED":
        return {"status": status, "score": None}
    score_key = f"{prefix}max_abs_delta" if prefix else "max_abs_delta"
    signed_key = f"{prefix}signed_delta" if prefix else "signed_delta"
    coord_key = f"{prefix}predicted_coordinate" if prefix else "predicted_coordinate"
    direction_key = f"{prefix}direction" if prefix else "direction"
    site_type_key = f"{prefix}site_type" if prefix else "site_type"
    score = finite_float(row[score_key], f"Pangolin {window} max_abs_delta")
    signed = finite_float(row[signed_key], f"Pangolin {window} signed_delta")
    coord = int_field(row[coord_key], f"Pangolin {window} predicted_coordinate")
    if score is None or signed is None or coord is None or score < 0:
        raise ValueError(f"invalid Pangolin {window} metrics")
    if not math.isclose(abs(signed), score, rel_tol=1e-6, abs_tol=1e-6):
        raise ValueError(f"Pangolin {window} max_abs_delta inconsistent with signed_delta")
    if row[direction_key] not in {"GAIN", "LOSS"}:
        raise ValueError(f"invalid Pangolin {window} direction")
    return {
        "status": status,
        "score": score,
        "score_text": row[score_key],
        "signed_text": row[signed_key],
        "coordinate": coord,
        "coordinate_text": row[coord_key],
        "direction": row[direction_key],
        "site_type": row[site_type_key],
    }


def classify_pangolin(row: dict[str, str], context: dict[str, Any], window: str) -> dict[str, Any]:
    metrics = parse_pangolin_metrics(row, window)
    prefix = "pangolin_" + window.lower() + "_"
    result: dict[str, Any] = {
        prefix + "status": metrics["status"],
        prefix + "score": metrics.get("score_text", ""),
        prefix + "signed_delta": metrics.get("signed_text", ""),
        prefix + "predicted_coordinate": str(metrics.get("coordinate", "")) if metrics.get("coordinate") is not None else "",
        prefix + "direction": metrics.get("direction", ""),
        prefix + "site_type": metrics.get("site_type", ""),
        prefix + "score_bin": score_bin(metrics.get("score")),
        prefix + "raw_coordinate_match_class": "UNSCORED" if metrics.get("score") is None else "RAW_NO_COORDINATE_MATCH",
        prefix + "geometry_class": "UNSCORED" if metrics.get("score") is None else "",
        prefix + "geometry_match_role": "",
        prefix + "geometry_match_site": "",
        prefix + "geometry_distance_nt": "",
        prefix + "altered_distance_nt": "",
        prefix + "replaced_distance_nt": "",
        prefix + "loss_distance_nt": "",
        prefix + "gain_distance_nt": "",
    }
    if metrics.get("score") is None:
        return result
    score = metrics["score"]
    coordinate = metrics["coordinate"]
    direction = metrics["direction"]
    if context["event_type"] == "CRYPTIC_ACCEPTOR_TRUNCATION_13NT":
        altered_site = context["altered"]["site"]
        replaced_site = context["replaced"]["site"]
        altered_distance = abs(coordinate - altered_site)
        replaced_distance = abs(coordinate - replaced_site)
        result[prefix + "altered_distance_nt"] = str(altered_distance)
        result[prefix + "replaced_distance_nt"] = str(replaced_distance)
        if direction == "GAIN" and altered_distance == 0:
            raw_class = "RAW_EXACT_COORDINATE_MATCH"
        elif direction == "GAIN" and altered_distance <= 2:
            raw_class = "RAW_NEAR_COORDINATE_MATCH_2NT"
        elif direction == "LOSS" and replaced_distance == 0:
            raw_class = "RAW_REPLACED_LOSS_EXACT_SUPPORTING"
        else:
            raw_class = "RAW_NO_COORDINATE_MATCH"
        result[prefix + "raw_coordinate_match_class"] = raw_class
        if direction == "GAIN" and altered_distance <= 2 and score >= LOW_THRESHOLD:
            result[prefix + "geometry_class"] = "EXACT_GEOMETRY_CONCORDANT" if altered_distance == 0 else "NEAR_GEOMETRY_CONCORDANT"
            result[prefix + "geometry_match_role"] = "acceptor"
            result[prefix + "geometry_match_site"] = str(altered_site)
            result[prefix + "geometry_distance_nt"] = str(altered_distance)
        else:
            result[prefix + "geometry_class"] = "SPLICE_EFFECT_ONLY" if score >= LOW_THRESHOLD else "PREDICTOR_LOW"
        return result

    observed = {endpoint["role"]: endpoint for endpoint in context["observed"]["endpoints"]}
    lost = context["lost_by_role"]
    loss_sites = [lost["donor"]["site"], lost["acceptor"]["site"]]
    gain_sites = [observed["donor"]["site"], observed["acceptor"]["site"]]
    loss_distance = min(abs(coordinate - site) for site in loss_sites)
    gain_distance = min(abs(coordinate - site) for site in gain_sites)
    result[prefix + "loss_distance_nt"] = str(loss_distance)
    result[prefix + "gain_distance_nt"] = str(gain_distance)
    loss_match = direction == "LOSS" and coordinate in loss_sites
    gain_match = direction == "GAIN" and coordinate in gain_sites
    raw_parts: list[str] = []
    if loss_match:
        raw_parts.append("RAW_CANONICAL_LOSS_EXACT")
    if gain_match:
        raw_parts.append("RAW_ALTERNATIVE_GAIN_EXACT")
    result[prefix + "raw_coordinate_match_class"] = "+".join(raw_parts) if raw_parts else "RAW_NO_COORDINATE_MATCH"
    if loss_match and score >= LOW_THRESHOLD:
        role = "donor" if coordinate == lost["donor"]["site"] else "acceptor"
        result[prefix + "geometry_class"] = "CANONICAL_SITE_LOSS_CONCORDANT"
        result[prefix + "geometry_match_role"] = role
        result[prefix + "geometry_match_site"] = str(coordinate)
        result[prefix + "geometry_distance_nt"] = "0"
    elif gain_match and score >= LOW_THRESHOLD:
        role = "donor" if coordinate == observed["donor"]["site"] else "acceptor"
        result[prefix + "geometry_class"] = "ALTERNATIVE_SITE_GAIN_COMPATIBLE"
        result[prefix + "geometry_match_role"] = role
        result[prefix + "geometry_match_site"] = str(coordinate)
        result[prefix + "geometry_distance_nt"] = "0"
    else:
        result[prefix + "geometry_class"] = "SPLICE_EFFECT_ONLY" if score >= LOW_THRESHOLD else "PREDICTOR_LOW"
        result[prefix + "geometry_distance_nt"] = str(min(loss_distance, gain_distance))
    return result


def consensus_class(spliceai_class: str, pangolin_class: str) -> str:
    if spliceai_class == "UNSCORED" or pangolin_class == "UNSCORED":
        return "UNSCORED"
    if spliceai_class not in VALID_CLASSES or pangolin_class not in VALID_CLASSES:
        raise ValueError(f"unknown model class combination: {spliceai_class}, {pangolin_class}")
    spliceai_geometry = spliceai_class in GEOMETRY_CLASSES
    pangolin_geometry = pangolin_class in GEOMETRY_CLASSES
    if spliceai_geometry and pangolin_geometry:
        return "PREDICTOR_FULLY_CONCORDANT"
    if spliceai_geometry or pangolin_geometry:
        return "PREDICTOR_PARTIALLY_CONCORDANT"
    if spliceai_class == "SPLICE_EFFECT_ONLY" or pangolin_class == "SPLICE_EFFECT_ONLY":
        return "PREDICTOR_SPLICE_ONLY_WRONG_OR_INCOMPLETE_GEOMETRY"
    if spliceai_class == "PREDICTOR_LOW" and pangolin_class == "PREDICTOR_LOW":
        return "RNA_ADDS_GEOMETRY_INFORMATION"
    raise ValueError(f"unhandled consensus class combination: {spliceai_class}, {pangolin_class}")


BASE_COLUMNS = [
    "rule_version", "threshold_note", "variant_uid", "gene", "event_type", "junction_id", "strand", "selected_transcript",
    "observed_raw_left", "observed_raw_right", "observed_site_left", "observed_site_right",
    "observed_donor_raw", "observed_donor_site", "observed_acceptor_raw", "observed_acceptor_site",
    "primary_competitor_junction_id", "secondary_competitor_junction_id",
    "lost_donor_raw", "lost_donor_site", "lost_acceptor_raw", "lost_acceptor_site",
    "cryptic_altered_raw", "cryptic_altered_site", "cryptic_replaced_raw", "cryptic_replaced_site",
]


def model_columns(prefix: str, spliceai: bool) -> list[str]:
    if spliceai:
        return [
            prefix + "status", prefix + "score", prefix + "score_bin", prefix + "max_type", prefix + "max_coordinate",
            prefix + "DS_AG", prefix + "DS_AL", prefix + "DS_DG", prefix + "DS_DL",
            prefix + "COORD_AG", prefix + "COORD_AL", prefix + "COORD_DG", prefix + "COORD_DL",
            prefix + "relevant_score", prefix + "relevant_score_bin", prefix + "raw_coordinate_match_class",
            prefix + "geometry_class", prefix + "geometry_match_component", prefix + "geometry_match_role",
            prefix + "geometry_match_site", prefix + "geometry_distance_nt", prefix + "altered_ag_distance_nt",
            prefix + "replaced_al_distance_nt", prefix + "loss_distance_nt", prefix + "gain_distance_nt",
        ]
    return [
        prefix + "status", prefix + "score", prefix + "signed_delta", prefix + "predicted_coordinate", prefix + "direction",
        prefix + "site_type", prefix + "score_bin", prefix + "raw_coordinate_match_class", prefix + "geometry_class",
        prefix + "geometry_match_role", prefix + "geometry_match_site", prefix + "geometry_distance_nt",
        prefix + "altered_distance_nt", prefix + "replaced_distance_nt", prefix + "loss_distance_nt", prefix + "gain_distance_nt",
    ]


OUTPUT_COLUMNS = BASE_COLUMNS + model_columns("spliceai_d50_", True) + model_columns("spliceai_d500_", True) + model_columns("pangolin_d50_", False) + model_columns("pangolin_d500_", False) + [
    "consensus_d50_class", "consensus_d50_rule", "interpretation_boundary",
]


def context_output(context: dict[str, Any]) -> dict[str, str]:
    observed = context["observed"]
    observed_by_role = {endpoint["role"]: endpoint for endpoint in observed["endpoints"]}
    row: dict[str, str] = {
        "rule_version": RULE_VERSION,
        "threshold_note": "0.2 is a descriptive working threshold only; not a validated clinical cutoff; continuous scores are primary",
        "variant_uid": context["key"]["variant_uid"],
        "gene": context["key"]["gene"],
        "event_type": context["event_type"],
        "junction_id": context["key"]["junction_id"],
        "strand": observed["strand"],
        "selected_transcript": context["key"]["selected_transcript"],
        "observed_raw_left": str(observed["left"]),
        "observed_raw_right": str(observed["right"]),
        "observed_site_left": str(next(endpoint["site"] for endpoint in observed["endpoints"] if endpoint["side"] == "left")),
        "observed_site_right": str(next(endpoint["site"] for endpoint in observed["endpoints"] if endpoint["side"] == "right")),
        "observed_donor_raw": str(observed_by_role["donor"]["raw"]),
        "observed_donor_site": str(observed_by_role["donor"]["site"]),
        "observed_acceptor_raw": str(observed_by_role["acceptor"]["raw"]),
        "observed_acceptor_site": str(observed_by_role["acceptor"]["site"]),
        "primary_competitor_junction_id": context["primary_text"],
        "secondary_competitor_junction_id": context["secondary_text"],
        "lost_donor_raw": "", "lost_donor_site": "", "lost_acceptor_raw": "", "lost_acceptor_site": "",
        "cryptic_altered_raw": "", "cryptic_altered_site": "", "cryptic_replaced_raw": "", "cryptic_replaced_site": "",
    }
    if context["event_type"] == "EXON_SKIP_1":
        for role, raw_key, site_key in (("donor", "lost_donor_raw", "lost_donor_site"), ("acceptor", "lost_acceptor_raw", "lost_acceptor_site")):
            row[raw_key] = str(context["lost_by_role"][role]["raw"])
            row[site_key] = str(context["lost_by_role"][role]["site"])
    else:
        row["cryptic_altered_raw"] = str(context["altered"]["raw"])
        row["cryptic_altered_site"] = str(context["altered"]["site"])
        row["cryptic_replaced_raw"] = str(context["replaced"]["raw"])
        row["cryptic_replaced_site"] = str(context["replaced"]["site"])
    return row


def predictor_output_fields(row: dict[str, str], model: str, window: str) -> dict[str, str]:
    if model == "spliceai":
        prefix = "spliceai_" + window.lower() + "_"
        input_prefix = component_prefix(window)
        values = {
            prefix + "DS_AG": row[f"{input_prefix}DS_AG"], prefix + "DS_AL": row[f"{input_prefix}DS_AL"],
            prefix + "DS_DG": row[f"{input_prefix}DS_DG"], prefix + "DS_DL": row[f"{input_prefix}DS_DL"],
            prefix + "COORD_AG": row[f"{input_prefix}COORD_AG"], prefix + "COORD_AL": row[f"{input_prefix}COORD_AL"],
            prefix + "COORD_DG": row[f"{input_prefix}COORD_DG"], prefix + "COORD_DL": row[f"{input_prefix}COORD_DL"],
        }
    else:
        values = {}
    return values


def load_step7_validation(root: Path) -> dict[str, Any]:
    fingerprint: dict[str, Any] = {}
    for rel, expected in EXPECTED_STEP7.items():
        path = root / rel
        if not path.exists():
            raise ValueError(f"missing Step7 frozen input: {rel}")
        observed = sha256_file(path)
        if observed != expected:
            raise ValueError(f"Step7 hash mismatch for {rel}: {observed} != {expected}")
        fingerprint[rel] = {"sha256": observed, "size_bytes": path.stat().st_size}
    audit_path = root / "results/audit/m4_m5_phase1/V23_PREDICTOR_RAW_FREEZE_AUDIT.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit.get("status") != "PASS" or audit.get("outcome_blind") is not True or audit.get("outcome_fields_used") != []:
        raise ValueError("Step7 raw freeze audit is not PASS/outcome-blind")
    if audit.get("n_variants") != EXPECTED_PREDICTOR_VARIANTS:
        raise ValueError("Step7 predictor variant count mismatch")
    if audit.get("primary", {}).get("distance") != 50 or audit.get("sensitivity", {}).get("distance") != 500:
        raise ValueError("Step7 distance freeze mismatch")
    if audit.get("primary", {}).get("continuous_score_primary") is not True:
        raise ValueError("Step7 primary continuous score flag missing")
    fingerprint["step7_audit_values"] = {
        "status": audit["status"], "outcome_blind": audit["outcome_blind"], "outcome_fields_used": audit["outcome_fields_used"],
        "n_variants": audit["n_variants"], "primary_distance": audit["primary"]["distance"], "sensitivity_distance": audit["sensitivity"]["distance"],
    }
    return fingerprint


def read_inputs(root: Path) -> tuple[list[dict[str, str]], dict[str, dict[str, str]], dict[str, dict[str, str]], dict[str, dict[str, str]], dict[str, dict[str, str]], dict[str, Any]]:
    keys_path = root / "results/audit/V23_POST_M2_M3_FROZEN9_KEYS.tsv"
    evidence_path = root / "results/tables/V23_R2_EVIDENCE_MATRIX_POST_M2_M3.tsv"
    competitors_path = root / "results/tables/V23_CANONICAL_COMPETITORS_POST_M3.tsv"
    spliceai_path = root / "results/tables/V23_M4_SPLICEAI_ALL_VARIANTS.tsv.gz"
    pangolin_path = root / "results/tables/V23_M5_PANGOLIN_ALL_VARIANTS.tsv.gz"
    keys = read_tsv(keys_path)
    require_columns(keys, FROZEN_REQUIRED, keys_path)
    if len(keys) != EXPECTED_EVENTS:
        raise ValueError(f"frozen9 key table must contain {EXPECTED_EVENTS} rows")
    key_map = join_key_rows(keys, "frozen9 keys")
    evidence_rows = read_tsv(evidence_path)
    require_columns(evidence_rows, EVIDENCE_ALLOWED, evidence_path)
    evidence_subset = [row for row in evidence_rows if row.get("variant_uid") in key_map]
    evidence = join_key_rows(
        [{column: row[column] for column in EVIDENCE_ALLOWED} for row in evidence_subset],
        "frozen9 evidence subset",
    )
    competitors_rows = read_tsv(competitors_path)
    require_columns(competitors_rows, COMPETITOR_ALLOWED, competitors_path)
    competitors_subset = [row for row in competitors_rows if row.get("variant_uid") in key_map]
    competitors = join_key_rows(
        [{column: row[column] for column in COMPETITOR_ALLOWED} for row in competitors_subset],
        "frozen9 competitor subset",
    )
    if set(evidence) != set(key_map) or set(competitors) != set(key_map):
        raise ValueError("frozen9/evidence/competitor key sets are not exactly equal")
    spliceai_rows = read_gzip_tsv(spliceai_path)
    pangolin_rows = read_gzip_tsv(pangolin_path)
    if len(spliceai_rows) != EXPECTED_PREDICTOR_VARIANTS or len(pangolin_rows) != EXPECTED_PREDICTOR_VARIANTS:
        raise ValueError("Step7 predictor tables do not contain 2374 variants")
    spliceai = join_key_rows(spliceai_rows, "SpliceAI table")
    pangolin = join_key_rows(pangolin_rows, "Pangolin table")
    if set(spliceai) != set(pangolin):
        raise ValueError("SpliceAI/Pangolin complete variant universes differ")
    step7 = load_step7_validation(root)
    return keys, evidence, competitors, spliceai, pangolin, step7


def build_rows(root: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    keys, evidence, competitors, spliceai, pangolin, step7 = read_inputs(root)
    rows: list[dict[str, str]] = []
    spliceai_counts: dict[str, Counter[str]] = {"D50": Counter(), "D500": Counter()}
    pangolin_counts: dict[str, Counter[str]] = {"D50": Counter(), "D500": Counter()}
    consensus_counts: Counter[str] = Counter()
    for key in keys:
        context = prepare_context(key, evidence[key["variant_uid"]], competitors[key["variant_uid"]])
        spliceai_row = spliceai[key["variant_uid"]]
        pangolin_row = pangolin[key["variant_uid"]]
        validate_predictor_target(spliceai_row, context, "SpliceAI")
        validate_predictor_target(pangolin_row, context, "Pangolin")
        if spliceai_row.get("target_gene_id") != pangolin_row.get("target_gene_id"):
            raise ValueError(f"SpliceAI/Pangolin target_gene_id mismatch for {key['variant_uid']}")
        out = context_output(context)
        for window in ("D50", "D500"):
            out.update(predictor_output_fields(spliceai_row, "spliceai", window))
            out.update(predictor_output_fields(pangolin_row, "pangolin", window))
            out.update(classify_spliceai(spliceai_row, context, window))
            out.update(classify_pangolin(pangolin_row, context, window))
        s_class = out["spliceai_d50_geometry_class"]
        p_class = out["pangolin_d50_geometry_class"]
        out["consensus_d50_class"] = consensus_class(s_class, p_class)
        out["consensus_d50_rule"] = "UNSCORED if either model unscored; both geometry -> FULLY; one -> PARTIALLY; else any SPLICE_EFFECT_ONLY -> SPLICE_ONLY; both LOW -> RNA_ADDS_GEOMETRY_INFORMATION"
        out["interpretation_boundary"] = "Outcome-blind descriptive predictor geometry only; no causal, clinical, or protein-level inference"
        rows.append({column: out.get(column, "") for column in OUTPUT_COLUMNS})
        spliceai_counts["D50"][s_class] += 1
        spliceai_counts["D500"][out["spliceai_d500_geometry_class"]] += 1
        pangolin_counts["D50"][p_class] += 1
        pangolin_counts["D500"][out["pangolin_d500_geometry_class"]] += 1
        consensus_counts[out["consensus_d50_class"]] += 1
    if len(rows) != EXPECTED_EVENTS or {row["variant_uid"] for row in rows} != {key["variant_uid"] for key in keys}:
        raise ValueError("Step8 output row/key validation failed")
    audit = {
        "status": "PASS",
        "rule_version": RULE_VERSION,
        "outcome_blind": True,
        "outcome_fields_used": [],
        "allowed_evidence_fields": list(EVIDENCE_ALLOWED),
        "allowed_competitor_fields": list(COMPETITOR_ALLOWED),
        "coordinate_rule": "junction raw intron endpoints converted to predictor splice-site coordinates: left raw - 1, right raw + 1; plus left donor/right acceptor, minus left acceptor/right donor",
        "thresholds": {
            "working_low_threshold": LOW_THRESHOLD,
            "descriptive_bins": [0.2, 0.5, 0.8],
            "clinical_cutoff": False,
            "continuous_score_primary": True,
            "geometry_relevant_score_rule": "a geometry coordinate match is concordant only if its relevant component (SpliceAI AG/AL/DG/DL or Pangolin max_abs_delta) is >=0.2; raw coordinate matches remain separately reported",
        },
        "n_rows": len(rows),
        "categories": {
            "spliceai_d50": dict(sorted(spliceai_counts["D50"].items())),
            "spliceai_d500": dict(sorted(spliceai_counts["D500"].items())),
            "pangolin_d50": dict(sorted(pangolin_counts["D50"].items())),
            "pangolin_d500": dict(sorted(pangolin_counts["D500"].items())),
            "consensus_d50": dict(sorted(consensus_counts.items())),
        },
        "step7_fingerprint": step7,
    }
    return rows, audit


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def self_test() -> None:
    # Conversion checks cover both sides and both strands.
    plus = parse_junction("chr1:100-200:+")
    minus = parse_junction("chr1:100-200:-")
    assert {(e["side"], e["role"], e["site"]) for e in plus["endpoints"]} == {
        ("left", "donor", 99), ("right", "acceptor", 201)
    }
    assert {(e["side"], e["role"], e["site"]) for e in minus["endpoints"]} == {
        ("left", "acceptor", 99), ("right", "donor", 201)
    }
    def fake_context(event_type: str, observed: str, primary: str, secondary: str = "") -> dict[str, Any]:
        key = {"variant_uid": "chr1:150:A:G", "gene": "X", "junction_id": observed, "selected_transcript": "ENSTX"}
        event = {"variant_uid": key["variant_uid"], "gene": "X", "junction_id": observed, "v23_altered_splice_type": event_type}
        comp = {"variant_uid": key["variant_uid"], "gene": "X", "junction_id": observed, "primary_competitor_junction_id": primary, "secondary_competitor_junction_id": secondary}
        return prepare_context(key, event, comp)
    exon_plus = fake_context("EXON_SKIP_1", "chr1:100-200:+", "chr1:100-130:+", "chr1:170-200:+")
    assert exon_plus["lost_by_role"]["acceptor"]["site"] == 131
    assert exon_plus["lost_by_role"]["donor"]["site"] == 169
    exon_minus = fake_context("EXON_SKIP_1", "chr1:100-200:-", "chr1:100-130:-", "chr1:170-200:-")
    assert exon_minus["lost_by_role"]["donor"]["site"] == 131
    assert exon_minus["lost_by_role"]["acceptor"]["site"] == 169
    cryptic_exact = fake_context("CRYPTIC_ACCEPTOR_TRUNCATION_13NT", "chr1:100-210:+", "chr1:100-200:+")
    assert cryptic_exact["altered"]["site"] == 211 and cryptic_exact["replaced"]["site"] == 201
    cryptic_near = fake_context("CRYPTIC_ACCEPTOR_TRUNCATION_13NT", "chr1:110-200:-", "chr1:100-200:-")
    assert cryptic_near["altered"]["site"] == 109 and cryptic_near["replaced"]["site"] == 99
    # Synthetic SpliceAI cryptic exact/near/non-match classes at meaningful and
    # low relevant scores; the raw match class remains available for the low case.
    def fake_splice_row(coords: dict[str, int], values: dict[str, float]) -> dict[str, str]:
        row = {"status": "SCORED", "D500_status": "SCORED", "outcome_blind_input": "TRUE"}
        for prefix in ("", "D500_"):
            row[prefix + "ds_max_type"] = "DS_AG"
            row[prefix + "ds_max_coordinate"] = str(coords["AG"])
            for name in ("AG", "AL", "DG", "DL"):
                row[prefix + "DS_" + name] = str(values[name])
                row[prefix + "COORD_" + name] = str(coords[name])
            row[prefix + "ds_max"] = str(max(values.values()))
        return row
    exact_row = fake_splice_row({"AG": 211, "AL": 201, "DG": 211, "DL": 201}, {"AG": 0.4, "AL": 0.0, "DG": 0.0, "DL": 0.0})
    assert classify_spliceai(exact_row, cryptic_exact, "D50")["spliceai_d50_geometry_class"] == "EXACT_GEOMETRY_CONCORDANT"
    near_row = fake_splice_row({"AG": 111, "AL": 99, "DG": 111, "DL": 99}, {"AG": 0.4, "AL": 0.0, "DG": 0.0, "DL": 0.0})
    assert classify_spliceai(near_row, cryptic_near, "D50")["spliceai_d50_geometry_class"] == "NEAR_GEOMETRY_CONCORDANT"
    nonmatch_row = fake_splice_row({"AG": 300, "AL": 99, "DG": 300, "DL": 99}, {"AG": 0.4, "AL": 0.0, "DG": 0.0, "DL": 0.0})
    assert classify_spliceai(nonmatch_row, cryptic_near, "D50")["spliceai_d50_geometry_class"] == "SPLICE_EFFECT_ONLY"
    low_match_row = fake_splice_row({"AG": 109, "AL": 99, "DG": 300, "DL": 99}, {"AG": 0.1, "AL": 0.0, "DG": 0.0, "DL": 0.0})
    low_result = classify_spliceai(low_match_row, cryptic_near, "D50")
    assert low_result["spliceai_d50_geometry_class"] == "PREDICTOR_LOW"
    assert low_result["spliceai_d50_raw_coordinate_match_class"] == "RAW_EXACT_COORDINATE_MATCH"
    print("SELF_TEST_PASS")


def run(root: Path) -> dict[str, Any]:
    rows, audit = build_rows(root)
    table_path = root / "results/tables/V23_PREDICTOR_RNA_GEOMETRY_FROZEN9.tsv"
    audit_path = root / "results/audit/m4_m5_phase1/V23_PREDICTOR_RNA_GEOMETRY_FROZEN9_AUDIT.json"
    manifest_path = root / "results/audit/m4_m5_phase1/V23_PREDICTOR_RNA_GEOMETRY_FROZEN9_SHA256.tsv"
    write_tsv(table_path, rows, OUTPUT_COLUMNS)
    output_hash = sha256_file(table_path)
    audit["output"] = {"path": relative_path(root, table_path), "size_bytes": table_path.stat().st_size, "sha256": output_hash}
    audit["script"] = {"path": relative_path(root, Path(__file__)), "size_bytes": Path(__file__).stat().st_size, "sha256": sha256_file(Path(__file__))}
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    manifest_rows = []
    for label, path in (("step8_output", table_path), ("step8_audit", audit_path), ("step8_script", Path(__file__))):
        manifest_rows.append({"label": label, "path": relative_path(root, path), "size_bytes": str(path.stat().st_size), "sha256": sha256_file(path)})
    write_tsv(manifest_path, manifest_rows, ["label", "path", "size_bytes", "sha256"])
    # The audit is deliberately written before the manifest; its content and
    # hash are therefore stable.  The manifest records that exact audit hash.
    return audit


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, help="v2.3 project root")
    parser.add_argument("--self-test", action="store_true", help="run synthetic geometry tests and exit")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.self_test:
        self_test()
        return 0
    if args.root is None:
        raise SystemExit("--root is required unless --self-test is used")
    audit = run(args.root.resolve())
    print(json.dumps({"status": audit["status"], "n_rows": audit["n_rows"], "categories": audit["categories"]}, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
