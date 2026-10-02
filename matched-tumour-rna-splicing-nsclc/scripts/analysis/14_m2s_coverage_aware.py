#!/usr/bin/env python3
"""Coverage-aware M2 sensitivity using the frozen v2.2 rare-outlier call.

Pseudo-carriers are sampled independently for every candidate pair.  Sampling
is matched to the observed carrier log1p(cluster_total) decile within the
pair's histology/gene/cluster pool, with deterministic adjacent-decile
fallback.  Outcome fields are not used for sampling.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


V22 = None
STATE: dict[str, Any] = {}
PAIR_KEY = ["histology", "variant_uid", "gene", "junction_id", "cluster_id"]
CHECKPOINT_SCHEMA = "V23_M2S_CHECKPOINT_V2"


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("v22_rare_impl", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import v2.2 rare implementation: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def choose_pseudo_carriers(task_index: int, iteration: int) -> tuple[set[str], np.ndarray]:
    task = STATE["tasks"][task_index]
    pool = STATE["pools"][task["pool_key"]]
    patients = pool["patients"]
    labels = pool["deciles"]
    desired = task["carrier_deciles"]
    rng = np.random.default_rng(
        np.random.SeedSequence([STATE["seed"], int(iteration), int(task_index)])
    )
    selected: list[int] = []
    excluded = np.asarray(task["real_carrier_indices"], dtype=np.int64)
    fallback = np.zeros(10, dtype=np.int64)
    for target in desired:
        chosen = None
        for distance in range(10):
            allowed = {int(target) - distance, int(target) + distance}
            allowed = {value for value in allowed if 0 <= value <= 9}
            candidates = np.flatnonzero(np.isin(labels, list(allowed)))
            candidates = candidates[~np.isin(candidates, excluded)]
            if selected:
                candidates = candidates[~np.isin(candidates, np.asarray(selected, dtype=np.int64))]
            if candidates.size:
                chosen = int(rng.choice(candidates))
                fallback[distance] += 1
                break
        if chosen is None:
            raise RuntimeError(f"no pseudo-carrier available for task {task_index}")
        selected.append(chosen)
    pseudo = {str(patients[index]) for index in selected}
    if len(pseudo) != len(desired):
        raise RuntimeError(f"pseudo-carrier cardinality failure for task {task_index}")
    return pseudo, fallback


def run_iteration(iteration: int) -> dict[str, Any]:
    n_tasks = len(STATE["tasks"])
    rare_vector = np.zeros(n_tasks, dtype=np.uint8)
    fallback = np.zeros(10, dtype=np.int64)
    coverage_failures = 0
    target_failures = 0
    carrier_failures = 0
    rare_tests = 0

    for index, task in enumerate(STATE["tasks"]):
        pseudo, used = choose_pseudo_carriers(index, iteration)
        fallback += used
        if len(pseudo) != task["n_carriers"]:
            carrier_failures += 1
        key = (task["histology"], task["link"]["variant_uid"], task["link"]["gene"])
        V22.CARRIER_GROUPS = {key: pseudo}
        recurrent, rare_rows = V22.analyze_task((task["link"], task["histology"]))
        observed_evaluable = int(recurrent.get("n_carriers_evaluable", 0)) + int(
            recurrent.get("n_noncarriers_evaluable", 0)
        )
        if observed_evaluable != task["expected_evaluable"]:
            coverage_failures += 1
        if recurrent.get("model_status") == "FAIL_TARGET_EXCEEDS_CLUSTER_TOTAL":
            target_failures += 1
        rare_tests += len(rare_rows)
        if any(as_bool(row.get("rare_call", False)) for row in rare_rows):
            rare_vector[index] = 1

    packed = np.packbits(rare_vector, bitorder="little").tobytes()
    r2_calls = int(rare_vector[STATE["r2_flags"]].sum())
    summary = {
        "iteration": int(iteration),
        "n_rare_calls": int(rare_vector.sum()),
        "n_r2_rare_calls": r2_calls,
        "rare_carrier_tests": int(rare_tests),
        "carrier_count_pass": bool(carrier_failures == 0),
        "carrier_count_failures": int(carrier_failures),
        "coverage_invariance_pass": bool(coverage_failures == 0),
        "coverage_failures": int(coverage_failures),
        "target_leq_cluster_total_pass": bool(target_failures == 0),
        "target_exceeds_cluster_total": int(target_failures),
        "rare_call_vector_sha256": hashlib.sha256(packed).hexdigest(),
    }
    for distance, count in enumerate(fallback.tolist()):
        summary[f"fallback_distance_{distance}"] = int(count)
    return {"packed": packed, "summary": summary}


def save_checkpoint(path: Path, summaries: list[dict[str, Any]], fingerprint: dict[str, Any],
                    reproducibility_pass: bool) -> None:
    atomic_json(path, {
        "schema": CHECKPOINT_SCHEMA,
        "run_fingerprint": fingerprint,
        "reproducibility_pass": bool(reproducibility_pass),
        "summaries": sorted(summaries, key=lambda row: int(row["iteration"])),
    })


def main() -> int:
    global V22, STATE
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-functional", required=True, type=Path)
    parser.add_argument("--carriers", required=True, type=Path)
    parser.add_argument("--sparse-counts", required=True, type=Path)
    parser.add_argument("--cluster-totals", required=True, type=Path)
    parser.add_argument("--v22-rare-implementation", required=True, type=Path)
    parser.add_argument("--expected-v22-sha256", required=True)
    parser.add_argument("--frozen-r2-keys", required=True, type=Path)
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--stats-output", required=True, type=Path)
    parser.add_argument("--audit-output", required=True, type=Path)
    parser.add_argument("--iterations", required=True, type=int)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--processes", required=True, type=int)
    parser.add_argument("--checkpoint-every", type=int, default=25)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--repro-check", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.processes <= 8:
        raise ValueError("processes must be in 1..8")
    if args.iterations < 1 or args.checkpoint_every < 1:
        raise ValueError("iterations and checkpoint-every must be positive")
    if not args.repro_check:
        raise ValueError("--repro-check is mandatory for every M2-S run")
    if "fork" not in mp.get_all_start_methods():
        raise RuntimeError("M2-S requires Linux fork to share the frozen read-only count tables")
    started = time.time()
    args.workdir.mkdir(parents=True, exist_ok=True)
    args.stats_output.parent.mkdir(parents=True, exist_ok=True)
    args.audit_output.parent.mkdir(parents=True, exist_ok=True)

    V22 = load_module(args.v22_rare_implementation)
    actual_v22_sha256 = sha256_file(args.v22_rare_implementation)
    if actual_v22_sha256.lower() != args.expected_v22_sha256.strip().lower():
        raise RuntimeError(
            f"frozen v2.2 implementation SHA256 mismatch: {actual_v22_sha256} "
            f"!= {args.expected_v22_sha256}"
        )
    candidates = pd.read_csv(
        args.candidate_functional, sep="\t", compression="gzip", low_memory=False
    ).fillna("")
    required = set(PAIR_KEY) | {
        "DNA_class", "primary_consequence", "selected_transcript", "junction_class",
        "link_rule", "n_carriers_total", "R2_FUNCTIONALLY_ELIGIBLE",
    }
    missing = sorted(required - set(candidates.columns))
    if missing:
        raise ValueError(f"candidate table missing columns: {missing}")
    candidates["histology"] = candidates["histology"].astype(str).str.replace("TCGA-", "", regex=False)
    candidates = candidates.sort_values(PAIR_KEY, kind="mergesort").reset_index(drop=True)
    if len(candidates) < 1 or candidates[PAIR_KEY].duplicated().any():
        raise ValueError("candidate universe must contain at least one unique pair")
    for column in ["n_carriers_total"]:
        numeric = pd.to_numeric(candidates[column], errors="coerce")
        if not np.isfinite(numeric.to_numpy()).all():
            raise ValueError(f"non-finite candidate field: {column}")
        if not np.isin(numeric.to_numpy(dtype=float), [1.0, 2.0]).all():
            raise ValueError("n_carriers_total must be an integer in {1,2} for every candidate pair")

    frozen = pd.read_csv(args.frozen_r2_keys, sep="\t", dtype=str).fillna("")
    if sorted(set(PAIR_KEY) - set(frozen.columns)) or len(frozen) < 1 or frozen[PAIR_KEY].duplicated().any():
        raise ValueError("frozen R2 key file must contain at least one unique PAIR_KEY row")
    frozen["histology"] = frozen["histology"].astype(str).str.replace("TCGA-", "", regex=False)
    if frozen[PAIR_KEY].duplicated().any():
        raise ValueError("frozen R2 keys are duplicated after histology normalization")
    frozen_key_set = {tuple(row[column] for column in PAIR_KEY) for row in frozen.to_dict("records")}
    candidate_key_set = {tuple(str(row[column]) for column in PAIR_KEY) for row in candidates.to_dict("records")}
    if not frozen_key_set.issubset(candidate_key_set):
        raise RuntimeError("one or more frozen R2 keys is absent from the candidate-pair universe")
    candidate_flagged = {
        tuple(str(row[column]) for column in PAIR_KEY)
        for row in candidates.to_dict("records") if as_bool(row["R2_FUNCTIONALLY_ELIGIBLE"])
    }
    if not frozen_key_set.issubset(candidate_flagged):
        raise RuntimeError("one or more frozen R2 keys is absent from the prespecified eligible null stratum")
    observed_r2 = len(frozen_key_set)

    carriers = pd.read_csv(args.carriers, sep="\t", compression="gzip", dtype=str).fillna("")
    carrier_required = {"histology", "variant_uid", "gene", "patient"}
    if carrier_required - set(carriers.columns):
        raise ValueError(f"carrier table missing columns: {sorted(carrier_required - set(carriers.columns))}")
    carrier_sets: dict[tuple[str, str, str], set[str]] = {}
    for row in carriers.itertuples(index=False):
        hist = str(row.histology).replace("TCGA-", "")
        carrier_sets.setdefault((hist, str(row.variant_uid), str(row.gene)), set()).add(str(row.patient))

    sparse = pd.read_csv(args.sparse_counts, sep="\t", compression="gzip", low_memory=False)
    totals = pd.read_csv(args.cluster_totals, sep="\t", compression="gzip", low_memory=False)
    sparse_required = {"histology", "junction_id", "patient", "count"}
    totals_required = {"histology", "gene", "cluster_id", "patient", "cluster_total"}
    if sparse_required - set(sparse.columns):
        raise ValueError(f"sparse table missing columns: {sorted(sparse_required - set(sparse.columns))}")
    if totals_required - set(totals.columns):
        raise ValueError(f"totals table missing columns: {sorted(totals_required - set(totals.columns))}")
    sparse["histology"] = sparse["histology"].astype(str).str.replace("TCGA-", "", regex=False)
    totals["histology"] = totals["histology"].astype(str).str.replace("TCGA-", "", regex=False)
    for frame, columns in [(sparse, ["count"]), (totals, ["cluster_total"])]:
        for column in columns:
            values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
            if not np.isfinite(values).all():
                raise ValueError(f"non-finite input values in {column}")
    total_unique_key = ["histology", "patient", "gene", "cluster_id"]
    if totals[total_unique_key].duplicated().any():
        raise ValueError("cluster totals are not unique by histology/patient/gene/cluster")
    sparse_lookup_key = ["histology", "junction_id", "patient"]
    sparse_duplicate_rows = sparse[sparse.duplicated(sparse_lookup_key, keep=False)].copy()
    sparse_duplicate_keys = int(sparse_duplicate_rows[sparse_lookup_key].drop_duplicates().shape[0])
    if sparse_duplicate_keys:
        conflicting = sparse_duplicate_rows.groupby(sparse_lookup_key, sort=False)["count"].nunique()
        if (conflicting > 1).any():
            raise ValueError("sparse target counts conflict within histology/junction/patient lookup key")
    if (sparse["count"].astype(float) < 0).any() or (totals["cluster_total"].astype(float) < 0).any():
        raise ValueError("negative count or coverage detected")
    V22.TARGET_GROUPS = {
        (hist, junction): dict(zip(group["patient"].astype(str), group["count"].astype(float)))
        for (hist, junction), group in sparse.groupby(["histology", "junction_id"], sort=False)
    }
    V22.TOTAL_GROUPS = {
        (hist, gene, cluster): group[["patient", "cluster_total"]].copy()
        for (hist, gene, cluster), group in totals.groupby(["histology", "gene", "cluster_id"], sort=False)
    }
    samples_by_histology = {
        str(hist): np.asarray(sorted(group["patient"].astype(str).unique()), dtype=object)
        for hist, group in totals.groupby("histology", sort=False)
    }

    pools: dict[tuple[str, str, str], dict[str, Any]] = {}
    tasks: list[dict[str, Any]] = []
    r2_flags: list[bool] = []
    real_carrier_zero_imputed_assignments = 0
    pairs_with_real_carrier_zero_imputed = 0
    model_columns = [
        "variant_uid", "gene", "junction_id", "cluster_id", "DNA_class",
        "primary_consequence", "selected_transcript", "junction_class", "link_rule",
    ]
    for row in candidates.to_dict("records"):
        hist = str(row["histology"]).replace("TCGA-", "")
        pool_key = (hist, str(row["gene"]), str(row["cluster_id"]))
        base = V22.TOTAL_GROUPS.get(pool_key)
        if base is None or base.empty:
            raise RuntimeError(f"missing cluster pool: {pool_key}")
        if pool_key not in pools:
            patients = samples_by_histology.get(hist)
            if patients is None or patients.size == 0:
                raise RuntimeError(f"missing histology-wide sample universe: {hist}")
            observed = base.assign(patient=base["patient"].astype(str)).set_index("patient")["cluster_total"]
            coverage_series = pd.Series(patients, index=patients).map(observed)
            zero_imputed = coverage_series.isna().to_numpy()
            coverage = coverage_series.fillna(0.0).to_numpy(dtype=float)
            transformed = np.log1p(coverage)
            edges = np.quantile(transformed, np.linspace(0.0, 1.0, 11), method="linear")
            labels = np.searchsorted(edges[1:-1], transformed, side="right").astype(np.int8)
            pools[pool_key] = {
                "patients": patients, "coverage": coverage, "edges": edges, "deciles": labels,
                "zero_imputed": zero_imputed,
                "patient_index": {patient: index for index, patient in enumerate(patients)},
            }
        pool = pools[pool_key]
        carrier_key = (hist, str(row["variant_uid"]), str(row["gene"]))
        real = sorted(carrier_sets.get(carrier_key, set()))
        expected_n = int(float(row["n_carriers_total"]))
        if len(real) != expected_n or expected_n not in {1, 2}:
            raise RuntimeError(f"carrier count mismatch for {carrier_key}: {len(real)} vs {expected_n}")
        try:
            real_indices = [int(pool["patient_index"][patient]) for patient in real]
            carrier_deciles = [int(pool["deciles"][index]) for index in real_indices]
        except KeyError as exc:
            raise RuntimeError(f"real carrier absent from candidate pair pool: {carrier_key} {exc}") from exc
        zero_imputed_real = int(np.sum(pool["zero_imputed"][real_indices]))
        real_carrier_zero_imputed_assignments += zero_imputed_real
        pairs_with_real_carrier_zero_imputed += int(zero_imputed_real > 0)
        tasks.append({
            "histology": hist,
            "pool_key": pool_key,
            "link": {column: row[column] for column in model_columns},
            "n_carriers": expected_n,
            "carrier_deciles": carrier_deciles,
            "real_carrier_indices": real_indices,
            "real_carrier_zero_imputed": zero_imputed_real,
            "expected_evaluable": int((pool["coverage"] >= 10).sum()),
        })
        r2_flags.append(as_bool(row["R2_FUNCTIONALLY_ELIGIBLE"]))

    STATE = {"seed": int(args.seed), "tasks": tasks, "pools": pools, "r2_flags": np.asarray(r2_flags, dtype=bool)}
    input_hashes = {
        str(path): sha256_file(path) for path in [
            args.candidate_functional, args.carriers, args.sparse_counts,
            args.cluster_totals, args.v22_rare_implementation, args.frozen_r2_keys,
        ]
    }
    fingerprint = {
        "schema": CHECKPOINT_SCHEMA,
        "iterations": int(args.iterations),
        "seed": int(args.seed),
        "candidate_pairs": int(len(tasks)),
        "observed_r2": int(observed_r2),
        "process_script_sha256": sha256_file(Path(__file__)),
        "input_sha256": input_hashes,
        "expected_v22_sha256": args.expected_v22_sha256.strip().lower(),
    }
    checkpoint = args.workdir / "V23_M2_COVERAGE_SENSITIVITY.checkpoint.json"
    summaries: list[dict[str, Any]] = []
    stored_reproducibility = False
    if args.resume and checkpoint.is_file():
        saved = json.loads(checkpoint.read_text(encoding="utf-8"))
        if saved.get("schema") != CHECKPOINT_SCHEMA or saved.get("run_fingerprint") != fingerprint:
            raise RuntimeError("checkpoint fingerprint mismatch; refusing unsafe resume")
        summaries = saved.get("summaries", [])
        stored_reproducibility = bool(saved.get("reproducibility_pass", False))
        ids = [int(row["iteration"]) for row in summaries]
        if len(ids) != len(set(ids)) or any(value < 0 or value >= args.iterations for value in ids):
            raise RuntimeError("checkpoint has duplicate or out-of-range iteration IDs")
    completed = {int(row["iteration"]) for row in summaries}
    pending = [iteration for iteration in range(args.iterations) if iteration not in completed]

    context = mp.get_context("fork")
    with context.Pool(processes=args.processes) as pool:
        if args.repro_check:
            duplicate = list(pool.imap(run_iteration, [0, 0], chunksize=1))
            reproducibility_pass = (
                duplicate[0]["packed"] == duplicate[1]["packed"]
                and duplicate[0]["summary"] == duplicate[1]["summary"]
            )
        else:
            reproducibility_pass = stored_reproducibility if completed else True
        for result in pool.imap_unordered(run_iteration, pending, chunksize=1):
            summary = result["summary"]
            summaries.append(summary)
            summaries.sort(key=lambda row: int(row["iteration"]))
            completed.add(int(summary["iteration"]))
            print(
                f"M2-S iteration={int(summary['iteration']) + 1}/{args.iterations} "
                f"r2={summary['n_r2_rare_calls']} exact={summary['fallback_distance_0']}",
                flush=True,
            )
            if len(completed) % args.checkpoint_every == 0 or len(completed) == args.iterations:
                save_checkpoint(checkpoint, summaries, fingerprint, reproducibility_pass)

    if len(completed) != args.iterations:
        raise RuntimeError(f"incomplete iterations: {len(completed)}/{args.iterations}")
    summary_df = pd.DataFrame(summaries).sort_values("iteration")
    required_pass = ["carrier_count_pass", "coverage_invariance_pass", "target_leq_cluster_total_pass"]
    if not reproducibility_pass or not all(bool(summary_df[column].all()) for column in required_pass):
        raise RuntimeError("M2-S reproducibility or iteration QC failed")
    summary_df.to_csv(args.workdir / "V23_M2_COVERAGE_SENSITIVITY_ITERATIONS.tsv", sep="\t", index=False)

    null = summary_df["n_r2_rare_calls"].astype(int).to_numpy()
    exceedance = int(np.sum(null >= observed_r2))
    empirical_p = float((exceedance + 1) / (args.iterations + 1))
    null_mean = float(np.mean(null))
    ratio = float(observed_r2 / null_mean) if null_mean > 0 else None
    result_row = {
        "analysis": "M2-S_COVERAGE_AWARE",
        "B": int(args.iterations),
        "MASTER_SEED": int(args.seed),
        "candidate_pairs": int(len(tasks)),
        "observed_R2": int(observed_r2),
        "null_mean": null_mean,
        "null_median": float(np.median(null)),
        "null_q95": float(np.quantile(null, 0.95, method="linear")),
        "null_q99": float(np.quantile(null, 0.99, method="linear")),
        "exceedance_ge9": exceedance,
        "empirical_p": empirical_p,
        "observed_expected_ratio": ratio,
        "ROBUST": "YES" if empirical_p <= 0.05 else "NO",
        "quantile_definition": "numpy linear interpolation at (n-1)*q",
    }
    pd.DataFrame([result_row]).to_csv(args.stats_output, sep="\t", index=False)
    fallback_totals = {
        str(distance): int(summary_df[f"fallback_distance_{distance}"].sum()) for distance in range(10)
    }
    audit = {
        "status": "PASS",
        "analysis": "M2-S_COVERAGE_AWARE",
        "result": result_row,
        "inputs": input_hashes,
        "run_fingerprint": fingerprint,
        "sampling_contract": {
            "scope": "independent within histology and within same candidate pair",
            "carrier_cardinality": "preserve observed N_carriers (1 or 2)",
            "coverage_transform": "log1p(cluster_total)",
            "coverage_pool": "all patients in the histology-wide frozen totals universe; absent pair-specific cluster rows are cluster_total=0",
            "zero_coverage_semantics": "matches frozen v2.2 CARRIER_CLUSTER_NOT_OBSERVED behavior and preserves non-evaluable zero-coverage labels",
            "decile_edges": "per unique (histology,gene,cluster_id), numpy linear quantiles q=0..1 by 0.1",
            "boundary_and_tie_rule": "searchsorted on internal edges with side=right; boundary ties enter upper bin; repeated edges may yield empty bins",
            "fallback": "symmetric adjacent-decile expansion by integer distance; uniform random choice within first non-empty distance",
            "without_replacement_within_pair": True,
            "original_real_carriers_eligible": False,
            "pseudo_carrier_pool": "observed noncarriers in the same candidate pair coverage pool",
            "rng": "numpy Generator(SeedSequence([MASTER_SEED, iteration, stable_pair_index]))",
        },
        "v22_function_reused_unchanged": "analyze_task",
        "v22_implementation_sha256": actual_v22_sha256,
        "frozen_r2_keys_sha256": sha256_file(args.frozen_r2_keys),
        "sampling_fields": ["histology", "gene", "cluster_id", "patient", "cluster_total", "n_carriers_total"],
        "outcome_fields_excluded_from_sampling": [
            "robust_z", "carrier_junction_usage", "junction_usage", "functional_consensus",
            "R2_FUNCTIONALLY_ELIGIBLE", "high_confidence", "Phase-0", "M2_event_level_result",
        ],
        "r2_flag_use": "post-call aggregation only; never used for pseudo-carrier selection or v2.2 task payload",
        "null_r2_eligible_pairs": int(np.sum(STATE["r2_flags"])),
        "frozen_nine_use": "observed count and membership/eligibility validation only; never defines the null aggregation mask",
        "candidate_sort": PAIR_KEY,
        "unique_coverage_pools": int(len(pools)),
        "zero_imputed_histology_pair_pool_entries": int(sum(np.sum(pool["zero_imputed"]) for pool in pools.values())),
        "real_carrier_pair_assignments_zero_imputed": int(real_carrier_zero_imputed_assignments),
        "candidate_pairs_with_real_carrier_zero_imputed": int(pairs_with_real_carrier_zero_imputed),
        "sparse_identical_duplicate_lookup_keys": sparse_duplicate_keys,
        "sparse_conflicting_duplicate_lookup_keys": 0,
        "sparse_duplicate_handling": "identical counts across gene/cluster aliases retained; v2.2 dict lookup is unchanged",
        "matching_selection_counts_by_fallback_distance": fallback_totals,
        "iteration_qc": {
            "rows": int(len(summary_df)),
            "ids_exact_0_to_B_minus_1": summary_df["iteration"].astype(int).tolist() == list(range(args.iterations)),
            "all_carrier_count_pass": bool(summary_df["carrier_count_pass"].all()),
            "all_coverage_invariance_pass": bool(summary_df["coverage_invariance_pass"].all()),
            "all_target_validity_pass": bool(summary_df["target_leq_cluster_total_pass"].all()),
        },
        "reproducibility_seed_check": bool(reproducibility_pass),
        "processes": int(args.processes),
        "elapsed_seconds": float(time.time() - started),
        "interpretation_rule": "p <= 0.05 => robust to coverage-aware sensitivity; otherwise downgrade wording and do not add another null",
        "outputs": [str(args.stats_output), str(args.audit_output)],
    }
    atomic_json(args.audit_output, audit)
    checkpoint.unlink(missing_ok=True)
    print(json.dumps(result_row, indent=2, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
