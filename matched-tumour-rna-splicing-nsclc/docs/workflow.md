# Analysis workflow

The analysis proceeds through the following components.

1. **Cohort and gene-set definition.** Select TCGA LUAD and LUSC primary tumours with matched DNA and RNA data and apply the clinically anchored 68-gene set.
2. **DNA–RNA relationship construction.** Map somatic variants to local splice-junction observations, define RNA-evaluable variant–histology instances and construct candidate DNA–RNA pairs.
3. **RNA-supported event identification.** Apply recurrent and rare-event RNA-support criteria and reconstruct interpretable splice consequences in the selected transcript.
4. **Conditional enrichment.** Compare the observed number of RNA-resolved events with the primary carrier-structure null and the coverage-matched conditional null, each using 10,000 iterations.
5. **Local splice-route redistribution.** Compare tumour-observed abnormal-junction usage with the corresponding canonical competitor.
6. **Sequence prediction and matched controls.** Evaluate SpliceAI and Pangolin, compare predicted splice-site geometry with tumour-observed RNA, and estimate within-set probability of superiority against matched controls.
7. **Transcript and normal-lung context.** Reconstruct selected and major alternative transcripts and quantify tumour-observed splice-route prevalence in GTEx lung.
8. **Allelic and molecular context.** Position event-bearing tumours relative to gene- and histology-specific canonical loss-of-function and clean-wild-type reference states.
9. **Orthogonal CPTAC reference state.** Characterise the canonical STK11-loss state in CPTAC LUAD using RNA-axis and protein measurements.

Key thresholds and model settings are listed in `../metadata/analysis_parameters.csv`.
