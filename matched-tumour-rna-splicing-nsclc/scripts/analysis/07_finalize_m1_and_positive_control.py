#!/usr/bin/env python3
"""Finalize M1 after independent function reconstruction and 9/9 control."""

from __future__ import annotations

import csv
import gzip
import json
import os
from pathlib import Path

import pandas as pd


V23_ROOT = Path(__file__).resolve().parents[1]
V22_ROOT = Path(os.environ.get("V22_ROOT", str(V23_ROOT.parent / "v2.2")))
TABLES = V23_ROOT / "results" / "tables"
AUDIT = V23_ROOT / "results" / "audit"


def main() -> int:
    raw_path = TABLES / "V23_RARE_D2_NULL_UNIVERSE.tsv.gz"
    function_path = TABLES / "V23_FULL_D2_FUNCTIONAL_POTENTIAL.tsv.gz"
    # V22_ROOT is the project root; v2.2 final_freeze is kept read-only at root level.
    freeze_candidates = [
        V22_ROOT / "final_freeze" / "v2.2" / "A_PHASE0_R2_UNMASKED_V22.tsv",
        V22_ROOT.parent / "final_freeze" / "v2.2" / "A_PHASE0_R2_UNMASKED_V22.tsv",
    ]
    freeze_r2_path = next((path for path in freeze_candidates if path.is_file()), freeze_candidates[0])
    if not raw_path.is_file() or not function_path.is_file():
        raise FileNotFoundError("M1 raw universe or independent functional table is missing")

    v22 = pd.read_csv(freeze_r2_path, sep="\t", dtype=str, low_memory=False).fillna("")
    function = pd.read_csv(function_path, sep="\t", compression="gzip", dtype=str, low_memory=False).fillna("")
    key = ["variant_uid", "gene", "junction_id"]
    compare = ["event_class", "altered_splice_type", "frame_consequence", "nmd_status", "functional_consensus"]
    signatures = function.groupby(key)[compare].nunique(dropna=False)
    if (signatures > 1).any().any():
        raise RuntimeError("functional signature differs across histology for the same candidate pair")
    function_unique = function.sort_values(["gene", "variant_uid", "junction_id", "histology"]).drop_duplicates(key)
    control = v22[key + compare].merge(function_unique[key + compare], on=key, how="left", suffixes=("_v22", "_v23"), validate="one_to_one")
    control_rows = []
    for row in control.to_dict("records"):
        per_field = {field: str(row[f"{field}_v22"]) == str(row[f"{field}_v23"]) for field in compare}
        control_rows.append({
            **{field: row[field] for field in key},
            "v22_r2_present": bool(row["event_class_v22"]),
            "v23_function_present": bool(row["event_class_v23"]),
            "all_function_fields_match": all(per_field.values()),
            **{f"match_{field}": value for field, value in per_field.items()},
        })
    control_df = pd.DataFrame(control_rows)
    control_pass = len(control_df) == 9 and bool(control_df["v23_function_present"].all()) and bool(control_df["all_function_fields_match"].all())
    control_df["control_status"] = "PASS" if control_pass else "FAIL"
    control_df.to_csv(AUDIT / "V23_R2_POSITIVE_CONTROL.tsv", sep="\t", index=False)

    merge_key = ["histology", "variant_uid", "gene", "junction_id", "cluster_id"]
    function_cols = merge_key + [
        "event_class", "altered_splice_type", "frame_consequence", "nmd_status", "functional_consensus",
        "R2_FUNCTIONALLY_ELIGIBLE", "v23_functional_potential_status", "v23_domain_status",
    ]
    function_merge = function[function_cols].drop_duplicates(merge_key)
    function_index = function_merge.set_index(merge_key)
    temp_path = raw_path.with_suffix(raw_path.suffix + ".functional_tmp")
    expected_rows = 0
    missing_function_rows = 0
    first = True
    with gzip.open(temp_path, "wt", encoding="utf-8", newline="") as handle:
        for chunk in pd.read_csv(raw_path, sep="\t", compression="gzip", chunksize=250_000, low_memory=False):
            expected_rows += len(chunk)
            merged = chunk.merge(function_merge, on=merge_key, how="left", suffixes=("", "_functional"), validate="many_to_one")
            # Count missing independent annotations before overwriting the raw placeholder
            # event_class column.  The raw M1 table intentionally contains blank/HOLD
            # functional fields, so checking merged["event_class"] here would falsely
            # reject every row.
            missing_function_rows += int(merged["event_class_functional"].isna().sum())
            for field in ["event_class", "altered_splice_type", "frame_consequence", "nmd_status", "functional_consensus"]:
                merged[field] = merged[f"{field}_functional"]
                merged.drop(columns=[f"{field}_functional"], inplace=True)
            merged["R2_FUNCTIONALLY_ELIGIBLE"] = merged["R2_FUNCTIONALLY_ELIGIBLE_functional"]
            merged["functional_potential_status"] = merged["v23_functional_potential_status"]
            merged["functional_reconstruction_independent"] = True
            merged.drop(columns=["R2_FUNCTIONALLY_ELIGIBLE_functional", "v23_functional_potential_status", "v23_domain_status"], inplace=True)
            merged.to_csv(handle, sep="\t", index=False, header=first)
            first = False
    if not control_pass or missing_function_rows != 0:
        temp_path.unlink(missing_ok=True)
        raise RuntimeError(f"M1 gate failed: positive_control={control_pass}, missing_function_rows={missing_function_rows}")
    os.replace(temp_path, raw_path)

    gate = {
        "m1_status": "PASS",
        "m2_allowed": True,
        "raw_universe_rows": expected_rows,
        "functional_candidate_rows": int(len(function)),
        "functional_unique_pair_keys": int(function_unique[key].drop_duplicates().shape[0]),
        "r2_positive_control": "PASS_9_OF_9",
        "positive_control_fields": compare,
        "missing_function_rows_after_merge": missing_function_rows,
        "v22_final_freeze_read_only": True,
        "m2_prerequisite": "independent functional potential PASS; 9/9 frozen R2 exact reconstruction PASS",
    }
    (AUDIT / "V23_M1_FUNCTIONAL_GATE.json").write_text(json.dumps(gate, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(gate, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
