#!/usr/bin/env python3
"""M3 quantitative reciprocal-splicing analysis for the frozen nine events.

The canonical competitor is read from the already frozen GENCODE50
``primary_competitor_junction_id`` column.  This script does not reconstruct,
rank, or otherwise select competitors from RNA measurements.  It evaluates
the altered junction and that one pre-defined competitor in the frozen
cluster-total coverage pool and performs the single pre-specified sign test
on the canonical carrier-minus-noncarrier deltas.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


KEY = ["histology", "variant_uid", "gene", "junction_id", "cluster_id", "selected_transcript"]
EVENT_KEY = KEY[:5]
PRIMARY = "primary_competitor_junction_id"
EXPECTED_EVENTS = 9
EXPECTED_M3_STATUS = "RECIPROCAL_DIRECTION_CONCORDANT_DESCRIPTIVE"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_histology(values: pd.Series) -> pd.Series:
    return values.astype(str).str.replace("TCGA-", "", regex=False)


def quantile_linear(values: np.ndarray, q: float) -> float:
    """NumPy linear quantile, with compatibility for older NumPy versions."""
    values = np.asarray(values, dtype=float)
    try:
        return float(np.quantile(values, q, method="linear"))
    except TypeError:  # pragma: no cover - compatibility path
        return float(np.quantile(values, q, interpolation="linear"))


def finite_nonnegative(values: pd.Series, label: str) -> None:
    numeric = pd.to_numeric(values, errors="coerce")
    if not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise ValueError(f"{label} contains non-finite values")
    if (numeric.to_numpy(dtype=float) < 0).any():
        raise ValueError(f"{label} contains negative values")


def exact_binomial_upper_tail(n: int, k: int, probability: float = 0.5) -> float:
    """Return P[X >= k] for X~Binomial(n, probability) without scipy."""
    if n <= 0:
        return 1.0
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    q = 1.0 - probability
    return float(sum(math.comb(n, i) * probability**i * q ** (n - i) for i in range(k, n + 1)))


def load_frozen_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    frozen = pd.read_csv(args.frozen_keys, sep="\t", dtype=str, low_memory=False).fillna("")
    if len(frozen) != EXPECTED_EVENTS:
        raise ValueError(f"frozen keys must contain {EXPECTED_EVENTS} events, found {len(frozen)}")
    missing = sorted(set(KEY) - set(frozen.columns))
    if missing:
        raise ValueError(f"frozen keys missing required columns: {missing}")
    if frozen[KEY].duplicated().any() or frozen[EVENT_KEY].duplicated().any():
        raise ValueError("frozen keys are not unique")
    if frozen[PRIMARY].eq("").any() if PRIMARY in frozen.columns else False:
        raise ValueError("frozen keys unexpectedly contain an empty primary competitor")

    competitors = pd.read_csv(args.competitors, sep="\t", dtype=str, low_memory=False).fillna("")
    missing = sorted(set(EVENT_KEY + [PRIMARY]) - set(competitors.columns))
    if missing:
        raise ValueError(f"canonical competitor table missing required columns: {missing}")
    if competitors[EVENT_KEY].duplicated().any():
        raise ValueError("canonical competitor table has duplicate event keys")
    if competitors[PRIMARY].str.strip().eq("").any():
        raise ValueError("one or more frozen events has an empty primary competitor")

    frozen_key_frame = frozen[KEY].sort_values(KEY, kind="mergesort").reset_index(drop=True)
    competitor_key_frame = competitors[KEY].sort_values(KEY, kind="mergesort").reset_index(drop=True)
    if not frozen_key_frame[KEY].equals(competitor_key_frame[KEY]):
        raise ValueError("frozen keys and canonical competitor keys are not exactly equal")

    existing = pd.read_csv(args.existing_m3_summary, sep="\t", dtype=str, low_memory=False).fillna("")
    if len(existing) != EXPECTED_EVENTS:
        raise ValueError(f"existing M3 summary must contain {EXPECTED_EVENTS} rows, found {len(existing)}")
    if not set(EVENT_KEY).issubset(existing.columns) or "m3_status" not in existing.columns:
        raise ValueError("existing M3 summary lacks frozen keys or m3_status")
    if existing[EVENT_KEY].duplicated().any():
        raise ValueError("existing M3 summary has duplicate event keys")
    if set(map(tuple, existing[EVENT_KEY].astype(str).to_numpy())) != set(
        map(tuple, frozen[EVENT_KEY].astype(str).to_numpy())
    ):
        raise ValueError("existing M3 summary keys do not exactly match frozen keys")
    if set(existing["m3_status"]) != {EXPECTED_M3_STATUS}:
        raise ValueError(
            "existing M3 summary is not descriptive reciprocal concordance for all nine events: "
            + json.dumps(existing["m3_status"].value_counts().to_dict(), ensure_ascii=False)
        )

    events = frozen[KEY].merge(
        competitors[EVENT_KEY + [PRIMARY, "secondary_competitor_junction_id", "all_reciprocal_reference_junctions", "m3_status"]],
        on=EVENT_KEY,
        how="left",
        validate="one_to_one",
    )
    if events[PRIMARY].isna().any() or events[PRIMARY].astype(str).str.strip().eq("").any():
        raise ValueError("canonical competitor merge produced an empty primary competitor")
    if set(events["m3_status"]) != {EXPECTED_M3_STATUS}:
        raise ValueError("canonical competitor table does not retain the required descriptive M3 status")
    return frozen, competitors, existing


def load_coverage_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, dict[tuple[str, str, str, str, str], float]]:
    totals = pd.read_csv(args.cluster_totals, sep="\t", low_memory=False)
    required_totals = {"histology", "patient", "gene", "cluster_id", "cluster_total"}
    missing = sorted(required_totals - set(totals.columns))
    if missing:
        raise ValueError(f"cluster totals missing required columns: {missing}")
    totals["histology"] = clean_histology(totals["histology"])
    totals["patient"] = totals["patient"].astype(str)
    totals["gene"] = totals["gene"].astype(str)
    totals["cluster_id"] = totals["cluster_id"].astype(str)
    finite_nonnegative(totals["cluster_total"], "cluster_total")
    total_key = ["histology", "patient", "gene", "cluster_id"]
    if totals.duplicated(total_key).any():
        raise ValueError("cluster totals are not unique at histology/patient/gene/cluster key")
    total_lookup = {
        (str(row.histology), str(row.patient), str(row.gene), str(row.cluster_id)): float(row.cluster_total)
        for row in totals.itertuples(index=False)
    }

    sparse_cols = ["histology", "patient", "gene", "junction_id", "cluster_id", "count", "cluster_total"]
    sparse = pd.read_csv(args.sparse_counts, sep="\t", usecols=sparse_cols, low_memory=False)
    sparse["histology"] = clean_histology(sparse["histology"])
    sparse["patient"] = sparse["patient"].astype(str)
    sparse["gene"] = sparse["gene"].astype(str)
    sparse["junction_id"] = sparse["junction_id"].astype(str)
    sparse["cluster_id"] = sparse["cluster_id"].astype(str)
    finite_nonnegative(sparse["count"], "sparse count")
    finite_nonnegative(sparse["cluster_total"], "sparse cluster_total")
    count_key = ["histology", "patient", "gene", "junction_id", "cluster_id"]
    if sparse.duplicated(count_key).any():
        raise ValueError("sparse counts are not unique at histology/patient/gene/junction/cluster key")
    if (sparse["count"].to_numpy(dtype=float) > sparse["cluster_total"].to_numpy(dtype=float) + 1e-12).any():
        raise ValueError("sparse count exceeds its recorded cluster_total")
    # The sparse table carries the same cluster total used for each count row;
    # verify it agrees with the frozen totals table wherever a sparse row exists.
    sparse_total = sparse.merge(
        totals[total_key + ["cluster_total"]].rename(columns={"cluster_total": "cluster_total_from_totals"}),
        on=total_key,
        how="left",
        validate="many_to_one",
    )
    if sparse_total["cluster_total_from_totals"].isna().any():
        raise ValueError("sparse count row lacks a matching frozen cluster total")
    if not np.allclose(
        sparse_total["cluster_total"].to_numpy(dtype=float),
        sparse_total["cluster_total_from_totals"].to_numpy(dtype=float),
        rtol=0,
        atol=1e-12,
    ):
        raise ValueError("sparse and cluster-total coverage values disagree")
    count_lookup = {
        (str(row.histology), str(row.patient), str(row.gene), str(row.junction_id), str(row.cluster_id)): float(row.count)
        for row in sparse.itertuples(index=False)
    }
    return totals, sparse, count_lookup


def run_analysis(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    started = time.perf_counter()
    frozen, competitors, existing = load_frozen_inputs(args)
    totals, sparse, count_lookup = load_coverage_inputs(args)
    sparse_presence = {
        (str(row.histology), str(row.gene), str(row.junction_id), str(row.cluster_id))
        for row in sparse[["histology", "gene", "junction_id", "cluster_id"]].itertuples(index=False)
    }

    carriers = pd.read_csv(args.carriers, sep="\t", compression="gzip", dtype=str, low_memory=False).fillna("")
    required_carriers = {"histology", "patient", "variant_uid", "gene"}
    missing = sorted(required_carriers - set(carriers.columns))
    if missing:
        raise ValueError(f"carrier table missing required columns: {missing}")
    carriers["histology"] = clean_histology(carriers["histology"])
    carriers["patient"] = carriers["patient"].astype(str)
    carriers["variant_uid"] = carriers["variant_uid"].astype(str)
    carriers["gene"] = carriers["gene"].astype(str)
    carrier_sets: dict[tuple[str, str, str], set[str]] = {}
    for row in carriers[["histology", "patient", "variant_uid", "gene"]].itertuples(index=False):
        carrier_sets.setdefault((str(row.histology), str(row.variant_uid), str(row.gene)), set()).add(str(row.patient))

    # Keep only lookup rows relevant to the frozen event keys.  Missing sparse
    # rows are explicitly interpreted as count=0 below.
    events = frozen.merge(
        competitors[EVENT_KEY + [PRIMARY, "secondary_competitor_junction_id", "m3_status"]],
        on=EVENT_KEY,
        how="left",
        validate="one_to_one",
    )
    total_groups = {
        (str(h), str(g), str(c)): group[["patient", "cluster_total"]].copy()
        for (h, g, c), group in totals.groupby(["histology", "gene", "cluster_id"], sort=False)
    }
    result_rows: list[dict[str, Any]] = []
    secondary_audit: list[dict[str, str]] = []
    delta_records: list[float] = []
    for event in events.to_dict("records"):
        hist = str(event["histology"])
        variant = str(event["variant_uid"])
        gene = str(event["gene"])
        altered = str(event["junction_id"])
        cluster = str(event["cluster_id"])
        canonical = str(event[PRIMARY])
        base = total_groups.get((hist, gene, cluster))
        if base is None:
            raise ValueError(f"missing cluster total group for {hist}/{gene}/{cluster}")
        pool = base[base["cluster_total"].astype(float) >= 10].copy()
        pool["patient"] = pool["patient"].astype(str)
        if pool.empty:
            raise ValueError(f"empty evaluable candidate pool for {hist}/{gene}/{cluster}")
        carrier_set = carrier_sets.get((hist, variant, gene), set())
        if len(carrier_set) not in {1, 2}:
            raise ValueError(f"frozen event carrier set must contain 1 or 2 patients, found {len(carrier_set)} for {hist}/{variant}/{gene}")
        carrier_patients = sorted(set(pool["patient"]) & carrier_set)
        noncarrier_patients = sorted(set(pool["patient"]) - carrier_set)
        if len(carrier_patients) < 1 or len(noncarrier_patients) < 20:
            raise ValueError(
                f"evaluable pool invariant failed for {hist}/{variant}/{gene}: "
                f"carriers={len(carrier_patients)}, noncarriers={len(noncarrier_patients)}"
            )

        totals_by_patient = dict(zip(pool["patient"], pool["cluster_total"].astype(float)))
        def usage(junction: str, patient: str) -> float:
            total = float(totals_by_patient[patient])
            count = float(count_lookup.get((hist, patient, gene, junction, cluster), 0.0))
            if not np.isfinite(count) or count < 0 or count > total + 1e-12:
                raise ValueError(f"invalid count/coverage for {hist}/{patient}/{gene}/{junction}/{cluster}")
            return count / total

        altered_carrier = np.asarray([usage(altered, patient) for patient in carrier_patients], dtype=float)
        altered_noncarrier = np.asarray([usage(altered, patient) for patient in noncarrier_patients], dtype=float)
        canonical_carrier = np.asarray([usage(canonical, patient) for patient in carrier_patients], dtype=float)
        canonical_noncarrier = np.asarray([usage(canonical, patient) for patient in noncarrier_patients], dtype=float)
        for label, values in {
            "altered_carrier": altered_carrier,
            "altered_noncarrier": altered_noncarrier,
            "canonical_carrier": canonical_carrier,
            "canonical_noncarrier": canonical_noncarrier,
        }.items():
            if not np.isfinite(values).all():
                raise ValueError(f"non-finite usage in {label} for {hist}/{variant}/{gene}")

        altered_carrier_usage = float(np.mean(altered_carrier))
        altered_noncarrier_median = float(np.median(altered_noncarrier))
        altered_noncarrier_q99 = quantile_linear(altered_noncarrier, 0.99)
        canonical_carrier_usage = float(np.mean(canonical_carrier))
        canonical_noncarrier_median = float(np.median(canonical_noncarrier))
        canonical_noncarrier_q10 = quantile_linear(canonical_noncarrier, 0.10)
        canonical_noncarrier_q99 = quantile_linear(canonical_noncarrier, 0.99)
        canonical_carrier_percentile = float(np.mean(canonical_noncarrier <= canonical_carrier_usage))
        canonical_delta = canonical_carrier_usage - canonical_noncarrier_median
        strong = bool(
            altered_carrier_usage > altered_noncarrier_q99
            and canonical_carrier_usage < canonical_noncarrier_q10
        )
        directional = bool(
            altered_carrier_usage > altered_noncarrier_q99
            and canonical_carrier_usage < canonical_noncarrier_median
        )
        classification = "RECIPROCAL_STRONG" if strong else ("RECIPROCAL_DIRECTIONAL" if directional else "NO_RECIPROCAL")
        delta_records.append(canonical_delta)
        secondary_audit.append({
            "histology": hist,
            "variant_uid": variant,
            "gene": gene,
            "junction_id": altered,
            "secondary_competitor_junction_id": str(event.get("secondary_competitor_junction_id", "")),
        })
        result_rows.append({
            "histology": hist,
            "variant_uid": variant,
            "gene": gene,
            "junction_id": altered,
            "cluster_id": cluster,
            "selected_transcript": str(event["selected_transcript"]),
            "primary_competitor_junction_id": canonical,
            "n_carriers_total": len(carrier_set),
            "n_carriers_evaluable": len(carrier_patients),
            "n_noncarriers_evaluable": len(noncarrier_patients),
            "aberrant_carrier_usage": altered_carrier_usage,
            "aberrant_noncarrier_median": altered_noncarrier_median,
            "aberrant_noncarrier_q99": altered_noncarrier_q99,
            "aberrant_delta": altered_carrier_usage - altered_noncarrier_median,
            "canonical_carrier_usage": canonical_carrier_usage,
            "canonical_noncarrier_median": canonical_noncarrier_median,
            "canonical_noncarrier_q10": canonical_noncarrier_q10,
            "canonical_noncarrier_q99": canonical_noncarrier_q99,
            "canonical_delta": canonical_delta,
            "canonical_carrier_percentile": canonical_carrier_percentile,
            "canonical_percentile_ties_included": True,
            "primary_competitor_present_in_cluster": bool((hist, gene, canonical, cluster) in sparse_presence),
            "classification": classification,
            "m3_status": str(event["m3_status"]),
        })

    quantitative = pd.DataFrame(result_rows).sort_values(EVENT_KEY, kind="mergesort").reset_index(drop=True)
    if len(quantitative) != EXPECTED_EVENTS or quantitative[EVENT_KEY].duplicated().any():
        raise ValueError("quantitative output does not contain exactly nine unique frozen events")
    deltas = np.asarray(delta_records, dtype=float)
    negative = int(np.sum(deltas < 0))
    positive = int(np.sum(deltas > 0))
    zero = int(np.sum(deltas == 0))
    n_nonzero = negative + positive
    p_value = exact_binomial_upper_tail(n_nonzero, negative, 0.5) if n_nonzero else 1.0
    sign_test = pd.DataFrame([{
        "n_events": int(len(deltas)),
        "n_negative": negative,
        "n_positive": positive,
        "n_zero": zero,
        "n_nonzero": n_nonzero,
        "p_value": p_value,
        "test": "exact_binomial_sign_test",
        "alternative": "canonical_delta < 0",
        "null": "P(negative)=P(positive)=0.5; zeros excluded",
    }])

    elapsed = time.perf_counter() - started
    audit: dict[str, Any] = {
        "status": "PASS",
        "scope": "frozen_9_R2",
        "events": EXPECTED_EVENTS,
        "reference_release": "GENCODE50_index_v2.2 / GRCh38 primary assembly",
        "competitor_selection": "frozen primary_competitor_junction_id from V23_CANONICAL_COMPETITORS_POST_M3.tsv; no data-driven re-selection",
        "secondary_competitor": "retained only in audit metadata; excluded from the single-event summary metrics",
        "coverage_rule": "histology/gene/cluster candidate pool with cluster_total >= 10",
        "usage_definition": "junction count / cluster_total; missing sparse count is zero",
        "carrier_summary_definition": "arithmetic mean across evaluable carrier patients",
        "noncarrier_quantiles": "numpy linear quantile; median is q=0.5, q10/q99 as named",
        "canonical_percentile": "mean(noncarrier_usage <= canonical_carrier_usage); ties included",
        "classification_rule": "RECIPROCAL_STRONG iff aberrant_carrier_usage > aberrant_noncarrier_q99 and canonical_carrier_usage < canonical_noncarrier_q10; RECIPROCAL_DIRECTIONAL iff aberrant_carrier_usage > aberrant_noncarrier_q99 and canonical_carrier_usage < canonical_noncarrier_median; else NO_RECIPROCAL",
        "sign_test": sign_test.iloc[0].to_dict(),
        "class_counts": quantitative["classification"].value_counts().to_dict(),
        "m3_status_counts": quantitative["m3_status"].value_counts().to_dict(),
        "secondary_competitors": secondary_audit,
        "input_sha256": {
            "frozen_keys": sha256_file(args.frozen_keys),
            "canonical_competitors": sha256_file(args.competitors),
            "existing_m3_summary": sha256_file(args.existing_m3_summary),
            "carriers": sha256_file(args.carriers),
            "sparse_counts": sha256_file(args.sparse_counts),
            "cluster_totals": sha256_file(args.cluster_totals),
            "script": sha256_file(Path(__file__).resolve()),
        },
        "counts_total_validation": "PASS: finite/non-negative counts and totals; count <= cluster_total for every sparse row; sparse totals agree with cluster-total table",
        "event_pool_validation": "PASS: each event has >=1 evaluable carrier and >=20 evaluable noncarriers",
        "existing_m3_status_validation": "PASS: all nine existing statuses are RECIPROCAL_DIRECTION_CONCORDANT_DESCRIPTIVE",
        "elapsed_seconds": elapsed,
    }
    return quantitative, sign_test, audit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frozen-keys", required=True, type=Path)
    parser.add_argument("--competitors", required=True, type=Path)
    parser.add_argument("--existing-m3-summary", required=True, type=Path)
    parser.add_argument("--carriers", required=True, type=Path)
    parser.add_argument("--sparse-counts", required=True, type=Path)
    parser.add_argument("--cluster-totals", required=True, type=Path)
    parser.add_argument("--project-root", required=True, type=Path)
    args = parser.parse_args()
    for path in [args.frozen_keys, args.competitors, args.existing_m3_summary, args.carriers, args.sparse_counts, args.cluster_totals]:
        if not path.is_file():
            raise FileNotFoundError(path)
    table_path = args.project_root / "results" / "tables" / "V23_M3_RECIPROCAL_QUANTITATIVE.tsv"
    sign_path = args.project_root / "results" / "stats" / "V23_M3_RECIPROCAL_SIGNTEST.tsv"
    audit_path = args.project_root / "results" / "audit" / "V23_M3_QUANTITATIVE_AUDIT.json"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    sign_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    quantitative, sign_test, audit = run_analysis(args)
    quantitative.to_csv(table_path, sep="\t", index=False)
    sign_test.to_csv(sign_path, sep="\t", index=False)
    audit["outputs"] = {
        "quantitative": str(table_path),
        "sign_test": str(sign_path),
        "audit": str(audit_path),
    }
    audit["output_sha256"] = {
        "quantitative": sha256_file(table_path),
        "sign_test": sha256_file(sign_path),
    }
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": audit["status"],
        "events": len(quantitative),
        "sign_test": sign_test.iloc[0].to_dict(),
        "class_counts": audit["class_counts"],
        "elapsed_seconds": audit["elapsed_seconds"],
        "outputs": audit["outputs"],
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
