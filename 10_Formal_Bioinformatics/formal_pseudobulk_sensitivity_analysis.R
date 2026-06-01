set.seed(20260524)
base_dir <- '/Users/wangxiaolong/Desktop/VSMC单细胞方法与结果'
lib_dir <- file.path(base_dir, 'R_formal_packages')
.libPaths(c(lib_dir, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

out_dir <- file.path(base_dir, '10_Formal_Bioinformatics')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
obj_path <- file.path(out_dir, 'VSMC_Formal_Substate_Annotated.rds')
obj <- readRDS(obj_path)
meta <- obj@meta.data %>% tibble::rownames_to_column('cell')
meta$status <- dplyr::case_when(
  meta$stage == 'STA' ~ 'STA',
  meta$stage %in% c('Ruptured_de_novo', 'Ruptured_recurrent') ~ 'Ruptured',
  TRUE ~ 'Unruptured'
)
meta_use <- meta %>% filter(status %in% c('STA', 'Ruptured'))

get_layer <- function(obj, assay='RNA', layer='counts') {
  tryCatch(GetAssayData(obj, assay=assay, layer=layer), error=function(e) GetAssayData(obj, assay=assay, slot=layer))
}
counts <- get_layer(obj, 'RNA', 'counts')[, meta_use$cell, drop=FALSE]

sample_info <- meta_use %>% distinct(sample, stage, status) %>% arrange(status, sample)
samples <- sample_info$sample
pb_counts <- sapply(samples, function(s) {
  cells <- meta_use$cell[meta_use$sample == s]
  Matrix::rowSums(counts[, cells, drop=FALSE])
})
colnames(pb_counts) <- samples
rownames(pb_counts) <- rownames(counts)
write.csv(sample_info, file.path(out_dir, 'VSMC_Pseudobulk_Sample_Info.csv'), row.names=FALSE)
write.csv(as.data.frame(pb_counts), file.path(out_dir, 'VSMC_Pseudobulk_Counts.csv'))

keep_samples <- sample_info %>% filter(status %in% c('STA','Ruptured'))
y <- DGEList(counts=pb_counts[, keep_samples$sample, drop=FALSE], group=keep_samples$status)
keep_genes <- filterByExpr(y, group=keep_samples$status, min.count=10)
y <- y[keep_genes,, keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~0 + factor(keep_samples$status, levels=c('STA','Ruptured')))
colnames(design) <- c('STA','Ruptured')
# limma-voom handles small sample sizes better than cell-level DE here; use trend moderation.
v <- voom(y, design, plot=FALSE)
fit <- lmFit(v, design)
contrast <- makeContrasts(Ruptured_vs_STA = Ruptured - STA, levels=design)
fit2 <- eBayes(contrasts.fit(fit, contrast), trend=TRUE, robust=TRUE)
de <- topTable(fit2, coef='Ruptured_vs_STA', number=Inf, sort.by='P') %>% tibble::rownames_to_column('gene')
write.csv(de, file.path(out_dir, 'VSMC_Pseudobulk_DE_Ruptured_vs_STA_limma_voom.csv'), row.names=FALSE)

sig_sets <- list(
  Contractile_VSMC = c('ACTA2','MYH11','TAGLN','CNN1','MYLK','CALD1','TPM1','TPM2','DES','LPP','CSRP1','MYL9'),
  Synthetic_VSMC = c('VCAN','FN1','COL1A1','COL1A2','COL3A1','THBS1','LGALS3','SPP1','VIM','MGP'),
  Inflammatory_VSMC = c('CCL2','CXCL8','IL6','NFKBIA','TNFAIP3','IRF1','STAT1','STAT3','ICAM1','CXCL2','CXCL3','HLA-DRA'),
  ECM_Remodeling = c('COL1A1','COL1A2','COL3A1','COL5A1','COL6A1','MMP2','MMP9','MMP14','TIMP1','FN1','POSTN','LOX'),
  AMPK_SIRT1_Protective = c('PRKAA1','PRKAA2','SIRT1','PPARGC1A','FOXO1','FOXO3','CAT','SOD2','NFE2L2','ADIPOQ','NAMPT'),
  Oxidative_Stress = c('HMOX1','DDIT3','TXNIP','SOD2','NQO1','ATF3','JUN','FOS','DNAJB1','HSPA1A'),
  NFkB_AP1 = c('NFKB1','RELA','NFKBIA','TNFAIP3','JUN','FOS','JUNB','FOSB','IER3','DUSP1'),
  CellChat_ECM_Axes = c('FN1','CD44','ITGAV','ITGB1','ITGA1','COL1A1','COL1A2','COL3A1','LAMC1','LAMB1')
)
ranked <- de$t
names(ranked) <- de$gene
ranked <- sort(ranked[is.finite(ranked)], decreasing=TRUE)

gsea_one <- function(genes, stats, nperm=10000) {
  genes <- intersect(unique(genes), names(stats))
  N <- length(stats); Nh <- length(genes)
  if (Nh < 3) return(c(NES=NA_real_, pvalue=NA_real_, size=Nh, direction=NA_real_))
  hits <- names(stats) %in% genes
  weights <- abs(stats)
  Phit <- cumsum(ifelse(hits, weights / sum(weights[hits]), 0))
  Pmiss <- cumsum(ifelse(!hits, 1/(N-Nh), 0))
  running <- Phit - Pmiss
  ES <- running[which.max(abs(running))]
  perm_es <- replicate(nperm, {
    ph <- rep(FALSE, N); ph[sample.int(N, Nh)] <- TRUE
    phit <- cumsum(ifelse(ph, weights / sum(weights[ph]), 0))
    pmiss <- cumsum(ifelse(!ph, 1/(N-Nh), 0))
    rr <- phit - pmiss
    rr[which.max(abs(rr))]
  })
  if (ES >= 0) {
    denom <- mean(perm_es[perm_es >= 0], na.rm=TRUE)
    p <- (sum(perm_es >= ES) + 1) / (sum(perm_es >= 0) + 1)
  } else {
    denom <- abs(mean(perm_es[perm_es < 0], na.rm=TRUE))
    p <- (sum(perm_es <= ES) + 1) / (sum(perm_es < 0) + 1)
  }
  NES <- ES / denom
  c(NES=NES, pvalue=p, size=Nh, direction=sign(ES))
}

gsea <- bind_rows(lapply(names(sig_sets), function(nm) {
  res <- gsea_one(sig_sets[[nm]], ranked, nperm=5000)
  data.frame(pathway=nm, NES=as.numeric(res['NES']), pvalue=as.numeric(res['pvalue']), size=as.integer(res['size']), direction=as.numeric(res['direction']))
})) %>% mutate(padj=p.adjust(pvalue, method='BH'), interpretation=ifelse(NES > 0, 'Higher in ruptured pseudobulk', 'Higher in STA pseudobulk')) %>% arrange(pvalue)
write.csv(gsea, file.path(out_dir, 'VSMC_Pseudobulk_Preranked_GSEA_Custom_Signatures.csv'), row.names=FALSE)

top_genes <- unique(c(
  head(de$gene[de$logFC > 0], 20),
  head(de$gene[de$logFC < 0], 20),
  unique(unlist(sig_sets))
))
top_genes <- intersect(top_genes, rownames(v$E))
expr_df <- as.data.frame(v$E[top_genes, , drop=FALSE]) %>% tibble::rownames_to_column('gene') %>% pivot_longer(-gene, names_to='sample', values_to='expr') %>% left_join(sample_info, by='sample')
# Show a compact subset: top 15 positive/negative DE plus key mechanistic genes available.
plot_genes <- unique(c(head(de$gene[de$logFC > 0], 15), head(de$gene[de$logFC < 0], 15), intersect(c('SIRT1','PRKAA1','PRKAA2','PPARGC1A','ACTA2','MYH11','TAGLN','CCL2','IL6','FN1','COL1A1','MMP2','MMP9','RELA','JUN','FOS'), top_genes)))
plot_genes <- intersect(plot_genes, top_genes)
plot_df <- expr_df %>% filter(gene %in% plot_genes) %>% group_by(gene) %>% mutate(z=as.numeric(scale(expr))) %>% ungroup()
plot_df$gene <- factor(plot_df$gene, levels=rev(plot_genes))
plot_df$sample <- factor(plot_df$sample, levels=sample_info$sample)
p_heat <- ggplot(plot_df, aes(sample, gene, fill=z)) + geom_tile(color='white', linewidth=0.25) +
  scale_fill_gradient2(low='#2A6FBB', mid='white', high='#C65A46', midpoint=0, name='z-score') +
  labs(title='VSMC pseudobulk expression: ruptured vs STA', x=NULL, y=NULL) +
  theme_classic(base_size=9) + theme(plot.title=element_text(face='bold'), axis.text.x=element_text(angle=35, hjust=1))
ggsave('/tmp/FigureS11_Pseudobulk_DE_GSEA.png', p_heat, width=10.5, height=9, dpi=300, bg='white')
file.copy('/tmp/FigureS11_Pseudobulk_DE_GSEA.png', file.path(out_dir, 'FigureS11_Pseudobulk_DE_GSEA.png'), overwrite=TRUE)

p_gsea <- ggplot(gsea, aes(reorder(pathway, NES), NES, fill=NES>0)) + geom_col(width=0.72) + coord_flip() +
  scale_fill_manual(values=c('TRUE'='#C65A46','FALSE'='#2A6FBB'), guide='none') +
  geom_hline(yintercept=0, linewidth=0.25) + labs(title='Custom preranked enrichment from VSMC pseudobulk DE', x=NULL, y='NES (Ruptured vs STA)') +
  theme_classic(base_size=9) + theme(plot.title=element_text(face='bold'))
ggsave('/tmp/FigureS11b_Pseudobulk_GSEA_Bar.png', p_gsea, width=7.5, height=4.8, dpi=300, bg='white')
file.copy('/tmp/FigureS11b_Pseudobulk_GSEA_Bar.png', file.path(out_dir, 'FigureS11b_Pseudobulk_GSEA_Bar.png'), overwrite=TRUE)

# Leave-one-ruptured-sample-out sensitivity for sample-level state proportion and VSMC module scores.
ruptured_samples <- sample_info$sample[sample_info$status == 'Ruptured']
score_cols <- c('AMPK_SIRT1_Score','Inflammation_Score','vsmc_Contractile','vsmc_Synthetic','vsmc_Inflammatory','vsmc_ECM_Remodeling','vsmc_Metabolic_Stress')
score_cols <- intersect(score_cols, colnames(meta_use))
base_state <- meta_use %>% count(sample, status, vsmc_state, name='cells') %>% group_by(sample, status) %>% mutate(total_vsmc=sum(cells), proportion=cells/total_vsmc) %>% ungroup()
base_scores <- meta_use %>% group_by(sample, status) %>% summarise(across(all_of(score_cols), ~mean(.x, na.rm=TRUE)), n_cells=n(), .groups='drop')

sensitivity <- bind_rows(lapply(c('none', ruptured_samples), function(drop_sample) {
  dat <- meta_use
  if (drop_sample != 'none') dat <- dat %>% filter(sample != drop_sample)
  state_summary <- dat %>% count(status, vsmc_state, name='cells') %>% group_by(status) %>% mutate(total_vsmc=sum(cells), proportion=cells/total_vsmc) %>% ungroup()
  inf_r <- state_summary %>% filter(status=='Ruptured', vsmc_state=='inflammatory VSMC') %>% pull(proportion)
  inf_s <- state_summary %>% filter(status=='STA', vsmc_state=='inflammatory VSMC') %>% pull(proportion)
  score_summary <- dat %>% group_by(status) %>% summarise(across(all_of(score_cols), ~mean(.x, na.rm=TRUE)), .groups='drop')
  get_score_diff <- function(col) {
    r <- score_summary %>% filter(status=='Ruptured') %>% pull(all_of(col))
    s <- score_summary %>% filter(status=='STA') %>% pull(all_of(col))
    ifelse(length(r)==1 & length(s)==1, r-s, NA_real_)
  }
  data.frame(
    analysis=ifelse(drop_sample=='none','all_ruptured_samples', paste0('drop_', drop_sample)),
    dropped_sample=ifelse(drop_sample=='none', NA, drop_sample),
    inflammatory_prop_ruptured=ifelse(length(inf_r)==1, inf_r, NA_real_),
    inflammatory_prop_STA=ifelse(length(inf_s)==1, inf_s, NA_real_),
    inflammatory_prop_delta=ifelse(length(inf_r)==1 & length(inf_s)==1, inf_r-inf_s, NA_real_),
    AMPK_SIRT1_delta=get_score_diff('AMPK_SIRT1_Score'),
    Inflammation_delta=get_score_diff('Inflammation_Score'),
    Contractile_delta=get_score_diff('vsmc_Contractile'),
    Synthetic_delta=get_score_diff('vsmc_Synthetic'),
    ECM_Remodeling_delta=get_score_diff('vsmc_ECM_Remodeling'),
    Metabolic_Stress_delta=get_score_diff('vsmc_Metabolic_Stress')
  )
}))
write.csv(base_state, file.path(out_dir, 'VSMC_Sensitivity_State_Proportions_BySample.csv'), row.names=FALSE)
write.csv(base_scores, file.path(out_dir, 'VSMC_Sensitivity_ModuleScores_BySample.csv'), row.names=FALSE)
write.csv(sensitivity, file.path(out_dir, 'VSMC_LeaveOneRupturedOut_Sensitivity.csv'), row.names=FALSE)

sens_long <- sensitivity %>% select(analysis, inflammatory_prop_delta, AMPK_SIRT1_delta, Inflammation_delta, Contractile_delta, ECM_Remodeling_delta) %>% pivot_longer(-analysis, names_to='metric', values_to='delta')
sens_long$analysis <- factor(sens_long$analysis, levels=sensitivity$analysis)
p_sens <- ggplot(sens_long, aes(analysis, delta, fill=metric)) + geom_col(position=position_dodge(width=0.78), width=0.68) +
  geom_hline(yintercept=0, linewidth=0.25) + coord_flip() +
  scale_fill_manual(values=c('inflammatory_prop_delta'='#A33E8A','AMPK_SIRT1_delta'='#2A6FBB','Inflammation_delta'='#C65A46','Contractile_delta'='#3A8F62','ECM_Remodeling_delta'='#D19A2E')) +
  labs(title='Leave-one-ruptured-sample-out sensitivity', x=NULL, y='Ruptured minus STA delta', fill=NULL) + theme_classic(base_size=9) + theme(plot.title=element_text(face='bold'))
ggsave('/tmp/FigureS12_LeaveOneOut_Sensitivity.png', p_sens, width=10, height=5.8, dpi=300, bg='white')
file.copy('/tmp/FigureS12_LeaveOneOut_Sensitivity.png', file.path(out_dir, 'FigureS12_LeaveOneOut_Sensitivity.png'), overwrite=TRUE)

summary_lines <- c(
  'VSMC pseudobulk and sensitivity analysis complete.',
  paste('Pseudobulk samples:', paste(keep_samples$sample, collapse=', ')),
  paste('Genes tested:', nrow(de)),
  paste('FDR<0.10 genes:', sum(de$adj.P.Val < 0.10, na.rm=TRUE)),
  paste('Top ruptured-up genes:', paste(head(de$gene[order(de$logFC, decreasing=TRUE)], 10), collapse=', ')),
  paste('Top STA-up genes:', paste(head(de$gene[order(de$logFC, decreasing=FALSE)], 10), collapse=', ')),
  paste('Custom enrichment top pathways:', paste(head(gsea$pathway, 5), collapse=', ')),
  paste('Leave-one-out inflammatory prop delta range:', paste(range(sensitivity$inflammatory_prop_delta, na.rm=TRUE), collapse=' to ')),
  paste('Leave-one-out AMPK-SIRT1 delta range:', paste(range(sensitivity$AMPK_SIRT1_delta, na.rm=TRUE), collapse=' to ')),
  paste('Leave-one-out inflammation score delta range:', paste(range(sensitivity$Inflammation_delta, na.rm=TRUE), collapse=' to '))
)
writeLines(summary_lines, file.path(out_dir, 'VSMC_Pseudobulk_Sensitivity_Run_Summary.txt'))
cat(paste(summary_lines, collapse='\n'), '\n')
