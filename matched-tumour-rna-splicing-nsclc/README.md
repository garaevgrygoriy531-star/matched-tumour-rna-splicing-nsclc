# Matched tumour RNA and somatic splice outcomes in NSCLC

Analysis code accompanying the manuscript **“Matched tumour RNA resolves somatic splice outcomes beyond sequence prediction in non-small-cell lung cancer.”**

## Study overview

The study integrates matched somatic DNA and tumour RNA from 1,016 TCGA primary tumours (516 lung adenocarcinomas and 500 lung squamous cell carcinomas) across a clinically anchored 68-gene analysis set. The reported analysis contains 2,373 RNA-evaluable variant–histology instances, 22,469 candidate DNA–RNA pairs, 78 supported DNA–RNA transitions and nine RNA-resolved events.

The repository contains code for the principal computational analyses, including cohort construction, RNA-supported splice-event analysis, conditional null models, local splice-route redistribution, sequence-prediction comparison, matched-control inference, transcript-context analysis, normal-lung route context, molecular reference-state modelling and the CPTAC STK11 reference-state analysis.

## Repository structure

```text
data/
  analysis_gene_set.csv       Clinically anchored 68-gene analysis set
metadata/
  analysis_parameters.csv     Main thresholds and analysis settings
  data_sources.csv            Public data resources and versions
  software_versions.csv       Software and reference versions reported in the Methods
scripts/analysis/             Analysis scripts
scripts/README.md              Script-to-analysis map
docs/
  workflow.md                 Analysis sequence
  code_scope.md               Repository scope and data requirements
  code_availability.md        Suggested manuscript code-availability wording
requirements.txt              Python package requirements
R_packages.txt                R package requirements
CITATION.cff                  Citation metadata
LICENSE                       MIT software licence
```

## Data

No participant-level or sample-level source data are distributed in this repository. The analyses use public resources from TCGA/GDC, MC3, recount3/Snaptron, GTEx, CIViC and CPTAC/LinkedOmicsKB together with public reference annotations and splice-prediction tools. Source information is listed in `metadata/data_sources.csv`.

Users should retrieve source data from the original providers and follow the corresponding access and use conditions. Derived patient-level tables, expression matrices, copy-number segments and junction-count matrices are intentionally excluded from version control.

## Software

The manuscript specifies the versions of the principal annotation and prediction tools used in the analysis. These are recorded in `metadata/software_versions.csv`. Python package requirements are listed in `requirements.txt`; R package requirements are listed in `R_packages.txt`. Package versions not reported in the Methods are not inferred here.

## Running the code

The scripts are organised by analysis component rather than as a single command-line workflow. Most Python modules expose their required input paths through command-line arguments. Where supported, use:

```bash
python scripts/analysis/<script>.py --help
```

R scripts take positional input paths as described in `scripts/README.md` and in the source files.

The principal numerical settings used in the manuscript—including event-support thresholds, 10,000 conditional-null iterations, splice-prediction windows and reference-state modelling criteria—are summarised in `metadata/analysis_parameters.csv`.

## Reproducibility scope

The repository provides the available source code for the main computational analyses but does not redistribute the underlying public datasets or project-specific intermediate tables. A small number of project-specific components are not included because the corresponding source code was not available in the final analysis archive; these are documented in `docs/code_scope.md`. No replacement implementation has been introduced for missing components.

## Citation

Please cite the associated manuscript and this software archive. Citation metadata are provided in `CITATION.cff`.

## Licence

This repository is released under the MIT License. External datasets and third-party software remain subject to their own licences and terms of use.
