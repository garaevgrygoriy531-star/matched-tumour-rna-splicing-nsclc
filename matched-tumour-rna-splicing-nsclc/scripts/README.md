# Analysis scripts

The scripts in `analysis/` correspond to the computational components described in the manuscript and Supplementary Methods.

| Analysis component | Main scripts |
|---|---|
| TCGA DNA/RNA sample mapping and cohort construction | `08_map_tcga_dna_rna.py`, `30_finalize_phase0_strict68.py` |
| Candidate DNA–RNA relationships and event eligibility | `02_build_rare_d2_null_universe.py`, `03_validate_m1_universe.py`, `06_reconstruct_functional_potential.py`, `07_finalize_m1_and_positive_control.py`, `14_select_primary_dna_consequence.py` |
| Primary conditional null | `08_m2_permutation.py` |
| Coverage-matched conditional null | `14_m2s_coverage_aware.py` |
| Reciprocal splice-route redistribution | `09_m3_reciprocal_splice.py`, `15_m3_quantitative_reciprocal.py` |
| SpliceAI/Pangolin annotation and RNA-geometry comparison | `16_build_gencode50_predictor_annotation.py`, `17_validate_m4_m5_variant_input.py`, `25_prepare_m4_m5_transport_vcf.py`, `27_merge_validate_predictor_vcfs.py`, `29_freeze_predictor_outputs_qc.py`, `30_classify_predictor_rna_geometry_frozen9.py` |
| Matched-control predictor comparison | `31_freeze_m6_matched_controls.py`, `32_analyze_m6_predictor_comparison.py` |
| Transcript consequence and normal-lung context | `33_refine_r2_protein_consequences.py`, `34_freeze_gtex_route_context.py`, `35_build_final_r2_evidence_matrix.py`, `36_select_representative_r2_events.py` |
| Molecular reference-state modelling | `04_run_corec_corrected_qc.R` |
| CPTAC STK11 reference-state analysis | `05_run_cptac_stk11_validation.R` |

Several scripts retain their original chronological filenames so that code provenance is preserved. Scientific definitions and thresholds should be taken from the manuscript, Supplementary Methods and `metadata/analysis_parameters.csv`.

The scripts require project-specific input tables derived from the public source datasets. These data are not redistributed here. Input requirements are exposed through command-line arguments where available and are summarised in `docs/code_scope.md`.
