# External cohort validation summary

Generated: 2026-06-06 22:21:24

## Included cohorts
   dataset      group n_samples
 GSE122897    Control        16
 GSE122897   Ruptured        21
 GSE122897 Unruptured        21
  GSE13353   Ruptured        11
  GSE13353 Unruptured         8
  GSE15629    Control         5
  GSE15629   Ruptured         8
  GSE15629 Unruptured         6
  GSE75436         IA        15
  GSE75436        STA        15

## Main significant signature-level findings (global BH-FDR < 0.05)
- GSE122897 IA_all_vs_Control: ECM_Remodeling higher (delta=0.88, FDR=0.0112)
- GSE122897 IA_all_vs_Control: Inflammatory higher (delta=0.57, FDR=0.014)
- GSE122897 RIA_vs_Control: ECM_Remodeling higher (delta=0.92, FDR=0.0112)
- GSE122897 RIA_vs_Control: Inflammatory higher (delta=0.51, FDR=0.0429)
- GSE122897 UIA_vs_Control: ECM_Remodeling higher (delta=0.85, FDR=0.014)
- GSE122897 UIA_vs_Control: Inflammatory higher (delta=0.63, FDR=0.0172)
- GSE13353 RIA_vs_UIA: Inflammatory higher (delta=0.76, FDR=0.0206)
- GSE15629 IA_all_vs_Control: Contractile lower (delta=-1.20, FDR=0.0335)
- GSE15629 IA_all_vs_Control: Inflammatory higher (delta=0.98, FDR=0.0112)
- GSE15629 RIA_vs_Control: Inflammatory higher (delta=0.81, FDR=0.0172)
- GSE15629 UIA_vs_Control: Contractile lower (delta=-1.60, FDR=0.0392)
- GSE15629 UIA_vs_Control: Inflammatory higher (delta=1.20, FDR=0.027)
- GSE75436 IA_vs_paired_STA: AMPK_SIRT1 lower (delta=-0.40, FDR=0.014)
- GSE75436 IA_vs_paired_STA: Contractile lower (delta=-1.33, FDR=0.0112)
- GSE75436 IA_vs_paired_STA: ECM_Remodeling higher (delta=1.13, FDR=0.0112)

## Interpretation guardrails
- These datasets are bulk aneurysm-wall or arterial-wall transcriptomic cohorts, so they validate tissue-level reproducibility rather than VSMC-specific single-cell states.
- Directional AMPK-SIRT1 changes should be interpreted together with inflammatory and ECM-remodeling scores, because stress compensation can raise individual metabolic-defense genes in inflamed tissue.
- Ruptured versus unruptured comparisons are not uniformly powered across cohorts; effect-size consistency is more informative than any single nominal p-value.

## Output files
- processed/signature_scores_all_external_cohorts.csv
- processed/signature_statistics_external_cohorts.csv
- processed/signature_trend_statistics_external_cohorts.csv
- processed/AMPK_signature_correlation_external_cohorts.csv
- processed/key_gene_statistics_external_cohorts.csv
- figures/Figure_External_Validation_Composite.pdf/png/tiff
- figures/Figure_External_Validation_A_signature_scores.pdf/png/tiff
- figures/Figure_External_Validation_B_effect_heatmap.pdf/png/tiff
- figures/Figure_External_Validation_C_AMPK_correlations.pdf/png/tiff
