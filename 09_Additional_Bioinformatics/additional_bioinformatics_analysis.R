Sys.setenv(LANG = "C")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

base_dir <- "/Users/wangxiaolong/Desktop/VSMC单细胞方法与结果"
out_dir <- file.path(base_dir, "09_Additional_Bioinformatics")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stage_order <- c("STA", "Unruptured", "Ruptured_de_novo", "Ruptured_recurrent")
stage_code_map <- c(STA = 1, Unruptured = 2, Ruptured_de_novo = 3, Ruptured_recurrent = 4)

sig_sets <- list(
  Contractile = c("ACTA2", "TAGLN", "MYH11", "CNN1", "MYLK", "SMTN", "LMOD1", "TPM2", "ACTG2"),
  AMPK_SIRT1 = c("SIRT1", "PRKAA1", "PRKAA2", "PPARGC1A", "TFAM", "FOXO3", "NAMPT", "SOD2", "CAT"),
  Inflammatory = c("IL1B", "CCL2", "CXCL8", "CXCL2", "NFKBIA", "TNFAIP3", "STAT1", "JUN", "FOS"),
  ECM_Remodeling = c("COL1A1", "COL1A2", "COL3A1", "FN1", "VCAN", "LUM", "DCN", "MMP2", "MMP9")
)

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 4 || length(unique(x[keep])) < 2 || length(unique(y[keep])) < 2) {
    return(c(rho = NA_real_, p = NA_real_))
  }
  res <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  c(rho = unname(res$estimate), p = res$p.value)
}

score_from_fetch <- function(obj, genes) {
  present <- intersect(genes, rownames(obj))
  if (length(present) == 0) return(rep(NA_real_, ncol(obj)))
  expr <- FetchData(obj, vars = present)
  rowMeans(expr, na.rm = TRUE)
}

theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#EEF3F7", color = "#A8B4BF"),
      plot.title = element_text(face = "bold", size = base_size + 1),
      axis.text.x = element_text(angle = 35, hjust = 1)
    )
}

message("Loading VSMC Seurat object...")
vsmc <- readRDS(file.path(base_dir, "03_VSMC_Macrophage_Analysis", "VSMC_Subset_Harmony.rds"))
emb <- as.data.frame(Embeddings(vsmc, "umap"))
colnames(emb) <- c("UMAP_1", "UMAP_2")
meta <- vsmc@meta.data
meta$cell <- rownames(meta)
meta$stage <- factor(meta$stage, levels = stage_order)
meta$stage_code <- as.numeric(stage_code_map[as.character(meta$stage)])
meta$seurat_cluster <- as.character(meta$seurat_clusters)

for (nm in names(sig_sets)) {
  meta[[paste0(nm, "_Score2")]] <- score_from_fetch(vsmc, sig_sets[[nm]])
}

# The VSMC subset RDS contains the core module scores, while the trajectory
# workflow stores pseudotime and, in some runs, the ECM score in a separate CSV.
pt_meta <- read.csv(file.path(base_dir, "04_Trajectory_Analysis", "VSMC_Pseudotime_Metadata.csv"),
                    check.names = FALSE)
pt_keep <- intersect(c("cell", "pseudotime", "ECM_Remodeling_Score"), colnames(pt_meta))
meta <- meta %>% left_join(pt_meta[, pt_keep, drop = FALSE], by = "cell", suffix = c("", "_trajectory"))
if (!"pseudotime" %in% colnames(meta) && "pseudotime_trajectory" %in% colnames(meta)) {
  meta$pseudotime <- meta$pseudotime_trajectory
}
if (!"ECM_Remodeling_Score" %in% colnames(meta) && "ECM_Remodeling_Score_trajectory" %in% colnames(meta)) {
  meta$ECM_Remodeling_Score <- meta$ECM_Remodeling_Score_trajectory
}
if (!"ECM_Remodeling_Score" %in% colnames(meta)) {
  meta$ECM_Remodeling_Score <- meta$ECM_Remodeling_Score2
}
if (!"pseudotime" %in% colnames(meta)) {
  meta$pseudotime <- NA_real_
}
state_scores <- cbind(meta, emb[meta$cell, , drop = FALSE])
write.csv(state_scores, file.path(out_dir, "VSMC_State_Signature_Scores.csv"), row.names = FALSE)

cluster_summary <- state_scores %>%
  group_by(seurat_cluster) %>%
  summarise(
    n_cells = n(),
    Contractile = mean(Contractile_Score2, na.rm = TRUE),
    AMPK_SIRT1 = mean(AMPK_SIRT1_Score2, na.rm = TRUE),
    Inflammatory = mean(Inflammatory_Score2, na.rm = TRUE),
    ECM_Remodeling = mean(ECM_Remodeling_Score2, na.rm = TRUE),
    AMPK_SIRT1_existing = mean(AMPK_SIRT1_Score, na.rm = TRUE),
    Inflammation_existing = mean(Inflammation_Score, na.rm = TRUE),
    ECM_existing = mean(ECM_Remodeling_Score, na.rm = TRUE),
    pseudotime = mean(pseudotime, na.rm = TRUE),
    .groups = "drop"
  )

cluster_long <- cluster_summary %>%
  select(seurat_cluster, n_cells, Contractile, AMPK_SIRT1, Inflammatory, ECM_Remodeling) %>%
  tidyr::pivot_longer(c(Contractile, AMPK_SIRT1, Inflammatory, ECM_Remodeling),
                      names_to = "signature", values_to = "mean_score") %>%
  group_by(signature) %>%
  mutate(z_score = as.numeric(scale(mean_score))) %>%
  ungroup()

write.csv(cluster_summary, file.path(out_dir, "VSMC_State_Cluster_Summary.csv"), row.names = FALSE)
write.csv(cluster_long, file.path(out_dir, "VSMC_State_Cluster_Signature_Long.csv"), row.names = FALSE)

sample_summary <- state_scores %>%
  group_by(sample, stage) %>%
  summarise(
    n_cells = n(),
    Contractile = mean(Contractile_Score2, na.rm = TRUE),
    AMPK_SIRT1_gene = mean(AMPK_SIRT1_Score2, na.rm = TRUE),
    Inflammatory_gene = mean(Inflammatory_Score2, na.rm = TRUE),
    ECM_gene = mean(ECM_Remodeling_Score2, na.rm = TRUE),
    AMPK_SIRT1_module = mean(AMPK_SIRT1_Score, na.rm = TRUE),
    Inflammation_module = mean(Inflammation_Score, na.rm = TRUE),
    ECM_module = mean(ECM_Remodeling_Score, na.rm = TRUE),
    pseudotime = mean(pseudotime, na.rm = TRUE),
    .groups = "drop"
  )
sample_summary$stage_code <- as.numeric(stage_code_map[as.character(sample_summary$stage)])
write.csv(sample_summary, file.path(out_dir, "VSMC_Sample_Level_Sensitivity.csv"), row.names = FALSE)

sample_stats <- lapply(c("Contractile", "AMPK_SIRT1_gene", "Inflammatory_gene", "ECM_gene",
                         "AMPK_SIRT1_module", "Inflammation_module", "ECM_module", "pseudotime"), function(v) {
  tr <- safe_spearman(sample_summary$stage_code, sample_summary[[v]])
  ru <- sample_summary %>% filter(stage %in% c("STA", "Ruptured_de_novo", "Ruptured_recurrent"))
  ru$binary <- ifelse(ru$stage == "STA", "STA", "Rupture_like")
  p_w <- NA_real_
  if (length(unique(ru$binary)) == 2) {
    p_w <- suppressWarnings(wilcox.test(ru[[v]] ~ ru$binary, exact = FALSE)$p.value)
  }
  data.frame(metric = v, sample_stage_spearman_rho = tr["rho"], sample_stage_spearman_p = tr["p"],
             STA_vs_rupture_like_wilcox_p = p_w)
})
sample_stats <- bind_rows(sample_stats)
write.csv(sample_stats, file.path(out_dir, "VSMC_Sample_Level_Sensitivity_Stats.csv"), row.names = FALSE)

pseudotime_stats <- lapply(c("AMPK_SIRT1_Score", "Inflammation_Score", "ECM_Remodeling_Score",
                             "Contractile_Score2", "AMPK_SIRT1_Score2",
                             "Inflammatory_Score2", "ECM_Remodeling_Score2"), function(v) {
  tr <- safe_spearman(state_scores$pseudotime, state_scores[[v]])
  data.frame(metric = v, pseudotime_spearman_rho = tr["rho"], pseudotime_spearman_p = tr["p"])
})
pseudotime_stats <- bind_rows(pseudotime_stats)
write.csv(pseudotime_stats, file.path(out_dir, "VSMC_Pseudotime_Module_Dynamics_Stats.csv"), row.names = FALSE)

bin_df <- state_scores %>%
  mutate(pseudotime_bin = ntile(pseudotime, 20)) %>%
  group_by(pseudotime_bin) %>%
  summarise(
    pseudotime = mean(pseudotime, na.rm = TRUE),
    AMPK_SIRT1 = mean(AMPK_SIRT1_Score2, na.rm = TRUE),
    Inflammation = mean(Inflammatory_Score2, na.rm = TRUE),
    ECM = mean(ECM_Remodeling_Score2, na.rm = TRUE),
    Contractile = mean(Contractile_Score2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(c(AMPK_SIRT1, Inflammation, ECM, Contractile),
                      names_to = "program", values_to = "mean_score")
bin_df <- bin_df %>% group_by(program) %>%
  mutate(mean_score_z = as.numeric(scale(mean_score))) %>%
  ungroup()
write.csv(bin_df, file.path(out_dir, "VSMC_Pseudotime_Binned_Module_Trajectories.csv"), row.names = FALSE)

message("Creating VSMC supplementary figure...")
score_plot_df <- sample_summary %>%
  tidyr::pivot_longer(c(Contractile, AMPK_SIRT1_gene, Inflammatory_gene, ECM_gene, pseudotime),
                      names_to = "metric", values_to = "value")
score_plot_df$metric <- factor(
  score_plot_df$metric,
  levels = c("Contractile", "AMPK_SIRT1_gene", "Inflammatory_gene", "ECM_gene", "pseudotime"),
  labels = c("Contractile", "AMPK-SIRT1 genes", "Inflammatory genes", "ECM genes", "Pseudotime")
)

p1 <- ggplot(cluster_long, aes(signature, seurat_cluster, fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", limits = c(-2.5, 2.5)) +
  labs(title = "A. Cluster-level VSMC state-score heatmap", x = NULL, y = "cluster", fill = "z") +
  theme_pub(9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

p2 <- ggplot(score_plot_df, aes(stage, value, color = stage)) +
  geom_point(aes(size = n_cells), alpha = 0.85, position = position_jitter(width = 0.08, height = 0)) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c(STA = "#4C78A8", Unruptured = "#59A14F",
                                Ruptured_de_novo = "#F28E2B", Ruptured_recurrent = "#E15759")) +
  labs(title = "B. Sample-level sensitivity summary", x = NULL, y = "sample mean", size = "VSMC n") +
  theme_pub(8) +
  theme(legend.position = "bottom")

p3 <- ggplot(bin_df, aes(pseudotime, mean_score_z, color = program)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  scale_color_manual(values = c(AMPK_SIRT1 = "#4C78A8", Inflammation = "#E15759",
                                ECM = "#59A14F", Contractile = "#8E6C8A")) +
  labs(title = "C. Binned gene-set dynamics along pseudotime", x = "pseudotime", y = "z-scored mean") +
  theme_pub(9) +
  theme(legend.position = "bottom")

stat_plot_df <- sample_stats %>%
  filter(metric %in% c("Contractile", "AMPK_SIRT1_gene", "Inflammatory_gene", "ECM_gene", "pseudotime"))
stat_plot_df$metric <- factor(
  stat_plot_df$metric,
  levels = c("Contractile", "AMPK_SIRT1_gene", "Inflammatory_gene", "ECM_gene", "pseudotime"),
  labels = c("Contractile", "AMPK-SIRT1 genes", "Inflammatory genes", "ECM genes", "Pseudotime")
)
p4 <- ggplot(stat_plot_df, aes(metric, sample_stage_spearman_rho, fill = sample_stage_spearman_rho)) +
  geom_col(width = 0.62) +
  geom_hline(yintercept = 0, color = "#555555", linewidth = 0.4) +
  geom_text(aes(label = sprintf("rho=%.2f\nP=%.3g", sample_stage_spearman_rho, sample_stage_spearman_p)),
            vjust = ifelse(stat_plot_df$sample_stage_spearman_rho >= 0, -0.15, 1.15), size = 2.8) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", limits = c(-1, 1)) +
  labs(title = "D. Sample-level monotonic trends across stage", x = NULL, y = "Spearman rho") +
  theme_pub(9) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))

png(file.path(out_dir, "FigureS5_VSMC_State_Sensitivity_Pseudotime.png"), width = 3600, height = 2600, res = 300)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

# Use a grid layout by explicitly composing through grid because patchwork/gridExtra may not be installed.
png(file.path(out_dir, "FigureS5_VSMC_State_Sensitivity_Pseudotime_panel.png"), width = 3600, height = 2600, res = 300)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p3, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p4, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

message("Running external GSE54083 co-expression module analysis...")
expr <- read.csv(file.path(base_dir, "GSE54083验证", "GSE54083_gene_expression_matrix.csv"), check.names = FALSE)
meta_bulk <- read.csv(file.path(base_dir, "GSE54083验证", "GSE54083_Module_Scores_Extended.csv"), check.names = FALSE)
rownames(expr) <- make.unique(expr$gene)
expr$gene <- NULL
common_samples <- intersect(colnames(expr), meta_bulk$sample)
expr <- expr[, common_samples]
meta_bulk <- meta_bulk[match(common_samples, meta_bulk$sample), ]

vars <- apply(expr, 1, var, na.rm = TRUE)
vars <- vars[is.finite(vars)]
top_n <- min(2000, length(vars))
top_genes <- names(sort(vars, decreasing = TRUE))[seq_len(top_n)]
expr_top <- as.matrix(expr[top_genes, ])
expr_z <- t(scale(t(expr_top)))
expr_z[!is.finite(expr_z)] <- 0

cor_gene <- cor(t(expr_z), method = "pearson")
cor_gene[!is.finite(cor_gene)] <- 0
dist_gene <- as.dist(1 - cor_gene)
hc <- hclust(dist_gene, method = "average")
k <- 6
modules <- cutree(hc, k = k)
module_ids <- paste0("M", modules)

gene_modules <- data.frame(gene = rownames(expr_z), module = module_ids, variance = vars[rownames(expr_z)])
module_scores <- sapply(sort(unique(module_ids)), function(m) {
  genes <- gene_modules$gene[gene_modules$module == m]
  colMeans(expr_z[genes, , drop = FALSE], na.rm = TRUE)
})
module_scores <- as.data.frame(module_scores)
module_scores$sample <- rownames(module_scores)
module_scores <- module_scores[, c("sample", sort(setdiff(colnames(module_scores), "sample")))]
write.csv(gene_modules, file.path(out_dir, "GSE54083_Coexpression_Gene_Modules.csv"), row.names = FALSE)
write.csv(module_scores, file.path(out_dir, "GSE54083_Coexpression_Module_Scores.csv"), row.names = FALSE)

traits <- meta_bulk[, c("sample", "stage_code", "AMPK_SIRT1_score", "Inflammation_score", "ECM_score")]
module_trait_rows <- list()
for (m in sort(unique(module_ids))) {
  ms <- module_scores[[m]]
  for (trait in setdiff(colnames(traits), "sample")) {
    st <- safe_spearman(ms, traits[[trait]])
    module_trait_rows[[length(module_trait_rows) + 1]] <- data.frame(
      module = m, trait = trait, spearman_rho = st["rho"], spearman_p = st["p"],
      n_genes = sum(gene_modules$module == m)
    )
  }
}
module_trait <- bind_rows(module_trait_rows)
write.csv(module_trait, file.path(out_dir, "GSE54083_Coexpression_Module_Trait_Correlations.csv"), row.names = FALSE)

signature_universe <- rownames(expr_z)
enrich_rows <- list()
for (m in sort(unique(module_ids))) {
  genes_m <- gene_modules$gene[gene_modules$module == m]
  for (sig in names(sig_sets)) {
    sig_genes <- intersect(sig_sets[[sig]], signature_universe)
    overlap <- intersect(genes_m, sig_genes)
    p <- if (length(sig_genes) > 0) {
      phyper(length(overlap) - 1, length(sig_genes), length(signature_universe) - length(sig_genes),
             length(genes_m), lower.tail = FALSE)
    } else {
      NA_real_
    }
    enrich_rows[[length(enrich_rows) + 1]] <- data.frame(
      module = m, signature = sig, overlap_n = length(overlap),
      module_n = length(genes_m), signature_n = length(sig_genes),
      p_value = p, overlap_genes = paste(overlap, collapse = ";")
    )
  }
}
enrich <- bind_rows(enrich_rows)
enrich$p_adj <- p.adjust(enrich$p_value, method = "BH")
write.csv(enrich, file.path(out_dir, "GSE54083_Coexpression_Signature_Enrichment.csv"), row.names = FALSE)

top_gene_rows <- list()
for (m in sort(unique(module_ids))) {
  genes <- gene_modules$gene[gene_modules$module == m]
  ms <- module_scores[[m]]
  cors <- sapply(genes, function(g) cor(as.numeric(expr_z[g, ]), ms, method = "spearman"))
  top <- head(names(sort(abs(cors), decreasing = TRUE)), 20)
  top_gene_rows[[length(top_gene_rows) + 1]] <- data.frame(module = m, gene = top, module_correlation = cors[top])
}
top_gene_df <- bind_rows(top_gene_rows)
write.csv(top_gene_df, file.path(out_dir, "GSE54083_Coexpression_Top_Module_Genes.csv"), row.names = FALSE)

mod_stage <- module_trait %>% filter(trait == "stage_code") %>% arrange(desc(abs(spearman_rho)))
top_pos <- module_trait %>% filter(trait %in% c("Inflammation_score", "ECM_score")) %>%
  group_by(module) %>% summarise(mean_abs = mean(abs(spearman_rho)), mean_rho = mean(spearman_rho), .groups = "drop") %>%
  arrange(desc(mean_abs)) %>% slice(1) %>% pull(module)
top_neg <- module_trait %>% filter(trait == "AMPK_SIRT1_score") %>%
  arrange(desc(abs(spearman_rho))) %>% slice(1) %>% pull(module)

mt_plot <- module_trait
mt_plot$trait <- factor(mt_plot$trait, levels = c("stage_code", "AMPK_SIRT1_score", "Inflammation_score", "ECM_score"))
p5 <- ggplot(mt_plot, aes(trait, module, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 2.7) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", limits = c(-1, 1)) +
  labs(title = "A. GSE54083 co-expression module-trait correlations", x = NULL, y = "module", fill = "rho") +
  theme_pub(8)

plot_module_df <- module_scores %>%
  select(sample, all_of(unique(c(top_neg, top_pos)))) %>%
  left_join(meta_bulk[, c("sample", "group", "stage_code")], by = "sample") %>%
  tidyr::pivot_longer(all_of(unique(c(top_neg, top_pos))), names_to = "module", values_to = "module_score")
plot_module_df$group <- factor(plot_module_df$group, levels = c("STA", "UIA", "RIA"))
p6 <- ggplot(plot_module_df, aes(group, module_score, color = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.15) +
  geom_point(size = 1.8, alpha = 0.9, position = position_jitter(width = 0.08, height = 0)) +
  facet_wrap(~module, scales = "free_y") +
  scale_color_manual(values = c(STA = "#4C78A8", UIA = "#59A14F", RIA = "#E15759")) +
  labs(title = "B. Representative module scores by disease group", x = NULL, y = "mean z-expression") +
  theme_pub(8) +
  theme(legend.position = "bottom")

enrich_plot <- enrich
enrich_plot$neg_log10_p <- -log10(pmax(enrich_plot$p_value, 1e-12))
p7 <- ggplot(enrich_plot, aes(signature, module, fill = neg_log10_p)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = overlap_n), size = 2.7) +
  scale_fill_gradient(low = "white", high = "#6A3D9A") +
  labs(title = "C. Signature overlap enrichment across modules", x = NULL, y = "module",
       fill = "-log10(P)") +
  theme_pub(8)

tg <- top_gene_df %>% group_by(module) %>% slice_head(n = 8) %>% ungroup()
tg$gene <- factor(tg$gene, levels = rev(unique(tg$gene)))
p8 <- ggplot(tg, aes(abs(module_correlation), gene, fill = module)) +
  geom_col() +
  facet_wrap(~module, scales = "free_y") +
  labs(title = "D. Top hub-like genes by module correlation", x = "|Spearman rho with module score|", y = NULL) +
  theme_pub(8) +
  theme(legend.position = "none")

png(file.path(out_dir, "FigureS6_GSE54083_Coexpression_Modules.png"), width = 3600, height = 2600, res = 300)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
print(p5, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p6, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p7, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p8, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

module_summary <- module_trait %>%
  select(module, trait, spearman_rho, spearman_p, n_genes) %>%
  tidyr::pivot_wider(names_from = trait, values_from = c(spearman_rho, spearman_p)) %>%
  left_join(top_gene_df %>% group_by(module) %>%
              summarise(top_genes = paste(head(gene, 10), collapse = ";"), .groups = "drop"),
            by = "module") %>%
  arrange(desc(abs(spearman_rho_stage_code)))
write.csv(module_summary, file.path(out_dir, "GSE54083_Coexpression_Module_Summary.csv"), row.names = FALSE)

message("Completed additional bioinformatics analysis.")
