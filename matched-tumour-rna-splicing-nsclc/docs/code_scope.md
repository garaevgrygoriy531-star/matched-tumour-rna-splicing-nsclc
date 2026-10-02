# Repository scope and data requirements

## Included code

The repository contains the available source code for the principal computational analyses reported in the manuscript: cohort mapping, candidate DNA–RNA relationships, conditional null models, splice-route redistribution, splice-prediction analysis, matched-control inference, transcript-context analysis, normal-lung route context, molecular reference-state modelling and the CPTAC STK11 reference-state analysis.

## Data not distributed

No participant-level source data are included. The following inputs must be obtained from the original public resources or reconstructed from them according to the Methods:

- TCGA/GDC sample metadata, MC3 somatic calls and ASCAT3 copy-number segments;
- recount3/Snaptron tumour junction counts and sample metadata;
- GENCODE/MANE transcript annotations and GRCh38 reference sequence;
- SpliceAI and Pangolin model assets;
- GTEx lung junction data queried through Snaptron;
- CPTAC/LinkedOmicsKB mutation, copy-number, RNA and proteomic data.

The repository does not contain patient-level crosswalks, expression matrices, junction-count matrices, event-level result tables or final manuscript figures.

## Components not represented by a complete standalone source module

The final analysis archive did not contain a standalone source module for the ASCAT3 allelic-context assignment used in Figure 5a, nor a verified final figure-assembly workflow. The transcript-context code included here covers the available reconstruction and refinement modules, but the original project-specific intermediate tables are not redistributed. The CPTAC R script contains the STK11/LUAD analysis used for the orthogonal reference-state results.

These limits concern code packaging rather than changes to the reported analyses. Missing components have not been recreated or replaced with new implementations.

## Re-running the analysis

The source scripts were written for structured intermediate tables produced from the public datasets. A complete rerun therefore requires reconstruction of those inputs using the versions and definitions specified in the manuscript and Supplementary Methods. `metadata/data_sources.csv`, `metadata/software_versions.csv` and `metadata/analysis_parameters.csv` provide the public-resource and parameter information available in the paper.
