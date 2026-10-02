#!/usr/bin/env python3
"""v2.3 M2 carrier-label permutation using the frozen v2.2 rare-call code.

The null randomizes carrier identities within histology while preserving the
observed carrier count for every (histology, variant, gene).  The v2.2
``analyze_task`` implementation is imported unchanged; functional annotations
are retained only for post-hoc stratification and are never passed as model
inputs.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import multiprocessing as mp
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


V22 = None
STATE: dict[str, Any] = {}


def load_v22_module(path: Path):
    spec = importlib.util.spec_from_file_location("v22_rare_impl", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import v2.2 rare implementation: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def run_one_iteration(iteration: int) -> dict[str, Any]:
    """Run all candidate tasks for one deterministic permutation.

    This function is executed in forked workers.  All reference counts and
    candidate tasks are copy-on-write globals; only carrier labels change per
    iteration and are assigned to the imported v2.2 module before evaluation.
    """
    rng = np.random.default_rng(np.random.SeedSequence([STATE["seed"], int(iteration)]))
    permuted = {}
    carrier_count_ok = True
    for key, n_carriers in STATE["observed_carrier_counts"].items():
        histology = key[0]
        sample_pool = STATE["samples_by_histology"][histology]
        if int(n_carriers) > len(sample_pool):
            carrier_count_ok = False
            permuted[key] = set()
            continue
        selected = rng.choice(sample_pool, size=int(n_carriers), replace=False)
        permuted[key] = {str(value) for value in selected.tolist()}
        carrier_count_ok = carrier_count_ok and len(permuted[key]) == int(n_carriers)

    V22.CARRIER_GROUPS = permuted
    n_tasks = len(STATE["tasks"])
    rare_calls = np.zeros(n_tasks, dtype=np.uint8)
    functional_rare_calls = np.zeros(n_tasks, dtype=np.uint8)
    r2_rare_calls = np.zeros(n_tasks, dtype=np.uint8)
    coverage_failures = 0
    target_exceeds_total = 0
    rare_carrier_tests = 0

    for index, (link, histology) in enumerate(STATE["tasks"]):
        recurrent, rare_rows = V22.analyze_task((link, histology))
        expected_evaluable = int(STATE["expected_evaluable"][index])
        observed_evaluable = int(recurrent.get("n_carriers_evaluable", 0)) + int(
            recurrent.get("n_noncarriers_evaluable", 0)
        )
        if observed_evaluable != expected_evaluable:
            coverage_failures += 1
        if recurrent.get("model_status") == "FAIL_TARGET_EXCEEDS_CLUSTER_TOTAL":
            target_exceeds_total += 1
        rare_carrier_tests += len(rare_rows)
        call = any(as_bool(row.get("rare_call", False)) for row in rare_rows)
        if call:
            rare_calls[index] = 1
            if STATE["functional_flags"][index]:
                functional_rare_calls[index] = 1
            if STATE["r2_flags"][index]:
                r2_rare_calls[index] = 1

    rare_packed = np.packbits(rare_calls, bitorder="little").tobytes()
    functional_packed = np.packbits(functional_rare_calls, bitorder="little").tobytes()
    r2_packed = np.packbits(r2_rare_calls, bitorder="little").tobytes()
    summary = {
        "iteration": int(iteration),
        "n_rare_calls": int(rare_calls.sum()),
        "n_functional_rare_calls": int(functional_rare_calls.sum()),
        "n_r2_rare_calls": int(r2_rare_calls.sum()),
        "rare_carrier_tests": int(rare_carrier_tests),
        "carrier_count_pass": bool(carrier_count_ok),
        "coverage_invariance_pass": bool(coverage_failures == 0),
        "coverage_failures": int(coverage_failures),
        "target_leq_cluster_total_pass": bool(target_exceeds_total == 0),
        "target_exceeds_cluster_total": int(target_exceeds_total),
        "rare_call_vector_sha256": sha256_bytes(rare_packed),
    }
    return {
        "iteration": int(iteration),
        "rare_packed": rare_packed,
        "functional_packed": functional_packed,
        "r2_packed": r2_packed,
        "summary": summary,
    }


def unpacked(value: bytes, n_tasks: int) -> np.ndarray:
    return np.unpackbits(np.frombuffer(value, dtype=np.uint8), bitorder="little")[:n_tasks]


def checkpoint_write(path: Path, rare_counts: np.ndarray, functional_counts: np.ndarray,
                    r2_counts: np.ndarray, completed: set[int], summaries: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    # Passing an open handle prevents numpy from appending a second ``.npz``
    # suffix to the temporary checkpoint path.
    with tmp.open("wb") as handle:
        np.savez_compressed(
            handle,
            rare_counts=rare_counts,
            functional_counts=functional_counts,
            r2_counts=r2_counts,
            completed=np.array(sorted(completed), dtype=np.int32),
        )
    os.replace(tmp, path)
    path.with_suffix(".summary.json").write_text(
        json.dumps(summaries, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def load_checkpoint(path: Path, n_tasks: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, set[int], list[dict[str, Any]]]:
    with np.load(path, allow_pickle=False) as loaded:
        rare = loaded["rare_counts"].astype(np.uint32, copy=True)
        functional = loaded["functional_counts"].astype(np.uint32, copy=True)
        r2 = loaded["r2_counts"].astype(np.uint32, copy=True)
        completed = {int(x) for x in loaded["completed"].tolist()}
    if any(array.shape != (n_tasks,) for array in (rare, functional, r2)):
        raise ValueError("checkpoint task dimension does not match current candidate universe")
    summary_path = path.with_suffix(".summary.json")
    summaries = json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.is_file() else []
    return rare, functional, r2, completed, summaries


def main() -> int:
    global V22, STATE
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidate-functional", required=True, type=Path)
    ap.add_argument("--carriers", required=True, type=Path)
    ap.add_argument("--sparse-counts", required=True, type=Path)
    ap.add_argument("--cluster-totals", required=True, type=Path)
    ap.add_argument("--v22-rare-implementation", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--iterations", required=True, type=int)
    ap.add_argument("--seed", required=True, type=int)
    ap.add_argument("--processes", required=True, type=int)
    ap.add_argument("--label", required=True)
    ap.add_argument("--checkpoint-every", type=int, default=10)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--repro-check", action="store_true")
    args = ap.parse_args()
    if not 1 <= args.processes <= 8:
        raise ValueError("processes must be within the single-task server limit 1..8")
    if args.iterations < 1:
        raise ValueError("iterations must be positive")
    if args.checkpoint_every < 1:
        raise ValueError("checkpoint-every must be positive")
    if args.repro_check and args.iterations < 1:
        raise ValueError("reproducibility check requires at least one iteration")
    args.outdir.mkdir(parents=True, exist_ok=True)

    V22 = load_v22_module(args.v22_rare_implementation)
    candidates = pd.read_csv(args.candidate_functional, sep="\t", compression="gzip", low_memory=False).fillna("")
    required = {
        "histology", "variant_uid", "gene", "junction_id", "cluster_id", "DNA_class",
        "primary_consequence", "selected_transcript", "junction_class", "link_rule",
        "n_carriers_total", "functional_consensus", "R2_FUNCTIONALLY_ELIGIBLE",
    }
    missing = sorted(required - set(candidates.columns))
    if missing:
        raise ValueError(f"candidate functional table missing columns: {missing}")
    candidate_key = ["histology", "variant_uid", "gene", "junction_id", "cluster_id"]
    if candidates[candidate_key].duplicated().any():
        raise ValueError("candidate functional table is not unique at histology/pair/cluster level")

    # This is the only task payload passed to v2.2.  Functional columns are
    # explicitly kept outside it and are used only after rare-call evaluation.
    model_columns = [
        "variant_uid", "gene", "junction_id", "cluster_id", "DNA_class",
        "primary_consequence", "selected_transcript", "junction_class", "link_rule",
    ]
    forbidden_model_columns = {
        "carrier_junction_usage", "junction_usage", "robust_z", "R2_FUNCTIONALLY_ELIGIBLE",
        "functional_consensus", "high_confidence", "Phase-0", "phase0",
    }
    if forbidden_model_columns.intersection(model_columns):
        raise AssertionError("functional/RNA evidence leaked into permutation model payload")
    tasks = []
    functional_flags = []
    r2_flags = []
    for row in candidates.to_dict("records"):
        histology = str(row["histology"])
        link = {column: row[column] for column in model_columns}
        tasks.append((link, histology))
        functional_flags.append(as_bool(row["functional_consensus"]))
        r2_flags.append(as_bool(row["R2_FUNCTIONALLY_ELIGIBLE"]))

    carriers = pd.read_csv(args.carriers, sep="\t", compression="gzip", dtype=str).fillna("")
    carrier_sets = {}
    for row in carriers.itertuples(index=False):
        histology = str(row.histology).replace("TCGA-", "")
        key = (histology, str(row.variant_uid), str(row.gene))
        carrier_sets.setdefault(key, set()).add(str(row.patient))
    observed_carrier_counts = {}
    for row in candidates[["histology", "variant_uid", "gene", "n_carriers_total"]].drop_duplicates().itertuples(index=False):
        key = (str(row.histology), str(row.variant_uid), str(row.gene))
        n_carriers = int(float(row.n_carriers_total))
        if key in observed_carrier_counts and observed_carrier_counts[key] != n_carriers:
            raise ValueError(f"inconsistent candidate carrier count for {key}")
        observed_carrier_counts[key] = n_carriers
    carrier_count_mismatches = {
        str(key): {"candidate": int(observed_carrier_counts[key]), "v22": len(carrier_sets.get(key, set()))}
        for key in observed_carrier_counts
        if len(carrier_sets.get(key, set())) != int(observed_carrier_counts[key])
    }
    if carrier_count_mismatches:
        raise ValueError(f"candidate/v2.2 carrier counts disagree: {list(carrier_count_mismatches.items())[:3]}")

    sparse = pd.read_csv(args.sparse_counts, sep="\t", compression="gzip", low_memory=False)
    totals = pd.read_csv(args.cluster_totals, sep="\t", compression="gzip", low_memory=False)
    V22.TARGET_GROUPS = {
        (hist, junction): dict(zip(group["patient"].astype(str), group["count"].astype(float)))
        for (hist, junction), group in sparse.groupby(["histology", "junction_id"], sort=False)
    }
    V22.TOTAL_GROUPS = {
        (hist, gene, cluster): group[["patient", "cluster_total"]].copy()
        for (hist, gene, cluster), group in totals.groupby(["histology", "gene", "cluster_id"], sort=False)
    }
    samples_by_histology = {
        str(hist): np.array(sorted(group["patient"].astype(str).unique()), dtype=object)
        for hist, group in totals.groupby("histology", sort=False)
    }
    expected_evaluable = []
    for link, histology in tasks:
        base = V22.TOTAL_GROUPS.get((histology, link["gene"], link["cluster_id"]))
        expected_evaluable.append(0 if base is None else int((base["cluster_total"].astype(float) >= 10).sum()))

    STATE = {
        "seed": int(args.seed),
        "tasks": tasks,
        "observed_carrier_counts": observed_carrier_counts,
        "samples_by_histology": samples_by_histology,
        "expected_evaluable": np.array(expected_evaluable, dtype=np.int32),
        "functional_flags": np.array(functional_flags, dtype=bool),
        "r2_flags": np.array(r2_flags, dtype=bool),
    }
    n_tasks = len(tasks)
    if n_tasks < 1:
        raise ValueError("candidate functional universe is empty")
    if any(key[0] not in samples_by_histology for key in observed_carrier_counts):
        raise ValueError("carrier permutation sample universe is missing a histology")
    if (totals["cluster_total"].astype(float) < 0).any() or (sparse["count"].astype(float) < 0).any():
        raise ValueError("negative coverage/count input detected")
    # The model payload is structural evidence for the leakage audit.
    model_input_columns = sorted(model_columns)
    functional_columns_excluded = sorted(forbidden_model_columns | {"functional_consensus", "R2_FUNCTIONALLY_ELIGIBLE"})

    outdir = args.outdir
    checkpoint = outdir / f"{args.label}.checkpoint.npz"
    rare_counts = np.zeros(n_tasks, dtype=np.uint32)
    functional_counts = np.zeros(n_tasks, dtype=np.uint32)
    r2_counts = np.zeros(n_tasks, dtype=np.uint32)
    completed: set[int] = set()
    summaries: list[dict[str, Any]] = []
    if args.resume and checkpoint.is_file():
        rare_counts, functional_counts, r2_counts, completed, summaries = load_checkpoint(checkpoint, n_tasks)
    pending = [idx for idx in range(args.iterations) if idx not in completed]

    ctx = mp.get_context("fork")
    with ctx.Pool(processes=args.processes) as pool:
        if args.repro_check and not completed:
            reproducibility = list(pool.imap(run_one_iteration, [0, 0], chunksize=1))
            reproducibility_pass = all(
                left["rare_packed"] == right["rare_packed"]
                and left["functional_packed"] == right["functional_packed"]
                and left["r2_packed"] == right["r2_packed"]
                and left["summary"]["rare_call_vector_sha256"] == right["summary"]["rare_call_vector_sha256"]
                for left, right in [(reproducibility[0], reproducibility[1])]
            )
        else:
            reproducibility_pass = True
        for result in pool.imap_unordered(run_one_iteration, pending, chunksize=1):
            iteration = int(result["iteration"])
            rare_counts += unpacked(result["rare_packed"], n_tasks).astype(np.uint32)
            functional_counts += unpacked(result["functional_packed"], n_tasks).astype(np.uint32)
            r2_counts += unpacked(result["r2_packed"], n_tasks).astype(np.uint32)
            completed.add(iteration)
            summaries.append(result["summary"])
            summaries.sort(key=lambda row: int(row["iteration"]))
            print(
                f"{args.label} iteration={iteration + 1}/{args.iterations} "
                f"rare={result['summary']['n_rare_calls']} "
                f"coverage_pass={result['summary']['coverage_invariance_pass']}",
                flush=True,
            )
            if len(completed) % args.checkpoint_every == 0 or len(completed) == args.iterations:
                checkpoint_write(checkpoint, rare_counts, functional_counts, r2_counts, completed, summaries)

    if len(completed) != args.iterations:
        raise RuntimeError(f"permutation incomplete: {len(completed)}/{args.iterations}")
    if any(not bool(row["carrier_count_pass"]) for row in summaries):
        raise RuntimeError("carrier-count invariance failed")
    if any(not bool(row["coverage_invariance_pass"]) for row in summaries):
        raise RuntimeError("coverage invariance failed")
    if any(not bool(row["target_leq_cluster_total_pass"]) for row in summaries):
        raise RuntimeError("target/cluster-total validity failed")

    prefix = args.label
    summary_df = pd.DataFrame(summaries).sort_values("iteration")
    summary_df.to_csv(outdir / f"{prefix}_ITERATION_SUMMARY.tsv", sep="\t", index=False)
    pair_out = candidates[candidate_key + ["DNA_class", "primary_consequence", "selected_transcript", "junction_class", "link_rule", "functional_consensus", "R2_FUNCTIONALLY_ELIGIBLE"]].copy()
    pair_out["null_rare_call_count"] = rare_counts
    pair_out["null_functional_rare_call_count"] = functional_counts
    pair_out["null_r2_rare_call_count"] = r2_counts
    pair_out["n_permutations"] = int(args.iterations)
    pair_out["null_rare_call_rate"] = pair_out["null_rare_call_count"] / float(args.iterations)
    pair_out["null_functional_rare_call_rate"] = pair_out["null_functional_rare_call_count"] / float(args.iterations)
    pair_out["null_r2_rare_call_rate"] = pair_out["null_r2_rare_call_count"] / float(args.iterations)
    pair_out.to_csv(outdir / f"{prefix}_PAIR_COUNTS.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})

    audit = {
        "status": "PASS",
        "label": prefix,
        "iterations": int(args.iterations),
        "seed": int(args.seed),
        "candidate_pairs": int(n_tasks),
        "candidate_variant_gene_histology_groups": int(len(observed_carrier_counts)),
        "carrier_count_distribution": {str(k): int(v) for k, v in pd.Series(list(observed_carrier_counts.values())).value_counts().sort_index().items()},
        "v22_rare_implementation": str(args.v22_rare_implementation),
        "v22_rare_implementation_sha256": hashlib.sha256(args.v22_rare_implementation.read_bytes()).hexdigest(),
        "permutation_scheme": "within-histology patient carrier-label permutation; observed carrier count preserved per (histology, variant, gene); labels generated independently with SeedSequence(seed, iteration)",
        "model_input_columns": model_input_columns,
        "functional_columns_excluded_from_model": functional_columns_excluded,
        "leakage_test": "PASS: functional/RNA/Phase-0 fields absent from v2.2 model task payload and only used for post-hoc strata",
        "coverage_test": "PASS: cluster-total evaluability invariant and target_count <= cluster_total for every iteration/task",
        "carrier_count_test": "PASS: every permutation preserved observed carrier cardinality",
        "reproducibility_test": "PASS" if reproducibility_pass else "FAIL",
        "reproducibility_seed_check": bool(reproducibility_pass),
        "all_iteration_coverage_pass": bool(summary_df["coverage_invariance_pass"].all()),
        "all_iteration_carrier_count_pass": bool(summary_df["carrier_count_pass"].all()),
        "all_iteration_target_validity_pass": bool(summary_df["target_leq_cluster_total_pass"].all()),
        "processes": int(args.processes),
        "m2_formal_allowed": bool(prefix.endswith("FORMAL")),
        "m3_not_released_by_this_script": True,
    }
    (outdir / f"{prefix}_AUDIT.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    checkpoint.unlink(missing_ok=True)
    checkpoint.with_suffix(".summary.json").unlink(missing_ok=True)
    print(json.dumps(audit, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
