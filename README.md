# Supplementary Code and Source Data

This repository contains the analysis scripts and processed source tables supporting the intracranial aneurysm scRNA-seq manuscript.

## External public data
- GSE54083: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54083

## Contents
- `09_Additional_Bioinformatics/`: VSMC state sensitivity, pseudotime dynamics, and GSE54083 co-expression module source tables.
- `10_Formal_Bioinformatics/`: formal VSMC substate, Milo differential-abundance, CellChat, DoRothEA/VIPER, SCENIC, pseudobulk, leave-one-out, pseudo-fate, and AMPK-SIRT1 compensation source tables.
- `01_QC_Integration/`, `02_Clustering_Annotation/`, `04_Trajectory_Analysis/`, `05_Virtual_Perturbation/`, `06_Network_Perturbation/`, `07_VSMC_DE_Enrichment/`, `08_GSE54083_Validation/`, and `GSE54083验证/`: processed source tables used for manuscript figures and analyses.

## Scope and data-sharing note
Large intermediate `.rds` objects, locally installed R package folders, and raw human scRNA-seq matrices are not included in this supplementary archive. Public release of raw matrices and full processed single-cell objects should follow the final institutional ethics/data-sharing approval and journal requirements.

## Reuse
The scripts are provided for peer review and reproducibility. Some scripts may require local path adjustment before rerunning in a new computing environment.

## Large table handling
Large text source tables over 20 MB were gzip-compressed to keep individual repository files within GitHub upload limits.

Compressed files:
- `10_Formal_Bioinformatics/SCENIC_VSMC_RegulonAUC_Long.csv.gz` (5.27 MB)
