Sys.setenv(LC_ALL = "C")
set.seed(20260522)

base_dir <- "/Users/wangxiaolong/Desktop/VSMC单细胞方法与结果"
lib_dir <- file.path(base_dir, "R_formal_packages")
.libPaths(c(lib_dir, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(patchwork)
  library(SingleCellExperiment)
  library(miloR)
  library(CellChat)
  library(dorothea)
  library(viper)
})

out_dir <- file.path(base_dir, "10_Formal_Bioinformatics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")
  flush.console()
}

save_png <- function(plot, filename, width = 11, height = 7.5, dpi = 300) {
  target <- file.path(out_dir, filename)
  tmp <- tempfile(pattern = "formal_plot_", fileext = ".png")
  ggsave(tmp, plot, width = width, height = height, dpi = dpi, bg = "white")
  ok <- file.copy(tmp, target, overwrite = TRUE)
  unlink(tmp)
  if (!ok) stop("Could not copy plot to ", target)
}

get_layer <- function(obj, assay = "RNA", layer = "data") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
}

present_genes <- function(obj, genes) {
  intersect(unique(genes), rownames(obj))
}

score_signatures <- function(obj, signatures, prefix = "sig") {
  DefaultAssay(obj) <- "RNA"
  for (nm in names(signatures)) {
    genes <- present_genes(obj, signatures[[nm]])
    if (length(genes) < 2) {
      obj[[paste0(prefix, "_", nm)]] <- NA_real_
      next
    }
    obj <- AddModuleScore(
      object = obj,
      features = list(genes),
      name = paste0(prefix, "_", nm),
      assay = "RNA",
      search = FALSE,
      seed = 20260522
    )
    generated <- tail(grep(paste0("^", prefix, "_", nm), colnames(obj@meta.data), value = TRUE), 1)
    obj[[paste0(prefix, "_", nm)]] <- obj@meta.data[[generated]]
  }
  obj
}

plot_empty <- function(label) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = label, size = 4) +
    theme_void()
}

signature_sets <- list(
  Contractile = c("ACTA2", "MYH11", "TAGLN", "CNN1", "MYLK", "SMTN", "TPM2", "DES", "LMOD1", "MYL9"),
  Synthetic = c("FN1", "VCAN", "THBS1", "SPARC", "SPP1", "LGALS3", "VIM", "CD44", "COL8A1", "ITGA5"),
  Inflammatory = c("CCL2", "CXCL8", "IL6", "NFKBIA", "TNFAIP3", "ICAM1", "VCAM1", "CXCL2", "CXCL3", "SOCS3"),
  ECM_Remodeling = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "MMP2", "MMP9", "MMP14", "TIMP1", "LOX", "POSTN", "FN1"),
  Metabolic_Stress = c("HMOX1", "DDIT3", "ATF3", "SOD2", "TXNIP", "PDK4", "LDHA", "SLC2A1", "HIF1A", "PPARGC1A")
)

global_marker_sets <- list(
  VSMC = c("ACTA2", "MYH11", "TAGLN", "CNN1", "MYL9", "TPM2", "COL1A1", "FN1"),
  Macrophage = c("LYZ", "C1QA", "C1QB", "C1QC", "CD68", "AIF1", "LST1", "FCGR3A", "MSR1"),
  EC = c("PECAM1", "VWF", "KDR", "CLDN5", "RAMP2", "ESAM", "CDH5"),
  Fibroblast = c("DCN", "LUM", "COL1A1", "COL3A1", "PDGFRA", "FBLN1", "COL6A1"),
  T_NK = c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1", "IL7R"),
  B_Plasma = c("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHG1"),
  Mast = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2")
)

state_colors <- c(
  "contractile VSMC" = "#2A6FBB",
  "synthetic VSMC" = "#C65A46",
  "inflammatory VSMC" = "#A33E8A",
  "ECM-remodeling VSMC" = "#3A8F62",
  "metabolic-stress VSMC" = "#D19A2E",
  "mixed/other VSMC" = "#7F7F7F"
)

celltype_colors <- c(
  VSMC = "#2A6FBB",
  Macrophage = "#C65A46",
  EC = "#3A8F62",
  Fibroblast = "#D19A2E",
  T_NK = "#7A5AA6",
  B_Plasma = "#4E9A9A",
  Mast = "#A66B3F",
  Other = "#8A8A8A"
)

run_vsmc_substates <- function() {
  log_msg("Loading VSMC object")
  cache_path <- file.path(out_dir, "VSMC_Formal_Substate_Annotated.rds")
  meta_path <- file.path(out_dir, "VSMC_Substate_Cell_Metadata.csv")
  prop_path <- file.path(out_dir, "VSMC_State_Proportions_by_Sample.csv")
  if (file.exists(cache_path) && file.exists(meta_path) && file.exists(prop_path)) {
    log_msg("Using cached VSMC substate object:", cache_path)
    vsmc <- readRDS(cache_path)
    emb <- read.csv(meta_path, check.names = FALSE)
    props <- read.csv(prop_path, check.names = FALSE)
    emb$vsmc_state <- factor(emb$vsmc_state, levels = names(state_colors))
    props$vsmc_state <- factor(props$vsmc_state, levels = names(state_colors))

    p_umap <- ggplot(emb, aes(vUMAP_1, vUMAP_2, color = vsmc_state)) +
      geom_point(size = 0.18, alpha = 0.65) +
      scale_color_manual(values = state_colors, drop = FALSE) +
      labs(x = "VSMC UMAP 1", y = "VSMC UMAP 2", color = "VSMC state", title = "VSMC substates") +
      theme_classic(base_size = 10) +
      theme(plot.title = element_text(face = "bold"), legend.position = "right")

    p_prop <- props %>%
      dplyr::mutate(sample = factor(sample, levels = unique(sample[order(stage, sample)]))) %>%
      ggplot(aes(sample, proportion, fill = vsmc_state)) +
      geom_col(width = 0.78) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
      scale_fill_manual(values = state_colors, drop = FALSE) +
      labs(x = NULL, y = "VSMC fraction", fill = "VSMC state", title = "Sample-level VSMC composition") +
      theme_classic(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"))
    return(list(vsmc = vsmc, p_umap = p_umap, p_prop = p_prop, props = props))
  }

  vsmc_path <- file.path(base_dir, "03_VSMC_Macrophage_Analysis", "VSMC_Subset_Harmony.rds")
  vsmc <- readRDS(vsmc_path)
  DefaultAssay(vsmc) <- "RNA"
  if (ncol(get_layer(vsmc, "RNA", "data")) == 0) {
    vsmc <- NormalizeData(vsmc, verbose = FALSE)
  }

  dims_use <- seq_len(min(20, ncol(Embeddings(vsmc, "harmony"))))
  vsmc <- FindNeighbors(
    vsmc,
    reduction = "harmony",
    dims = dims_use,
    graph.name = c("vsmc_nn", "vsmc_snn"),
    verbose = FALSE
  )

  chosen <- NULL
  for (res in c(0.25, 0.4, 0.55, 0.7, 0.9, 1.1)) {
    tmp <- FindClusters(vsmc, graph.name = "vsmc_snn", resolution = res, algorithm = 1, verbose = FALSE)
    ncl <- length(unique(as.character(Idents(tmp))))
    chosen <- list(obj = tmp, res = res, ncl = ncl)
    if (ncl >= 5 && ncl <= 12) break
  }
  vsmc <- chosen$obj
  vsmc$vsmc_subcluster <- as.character(Idents(vsmc))
  log_msg("VSMC reclustering resolution", chosen$res, "yielded", chosen$ncl, "clusters")

  vsmc <- RunUMAP(
    vsmc,
    reduction = "harmony",
    dims = dims_use,
    reduction.name = "umap_vsmc_sub",
    reduction.key = "vUMAP_",
    verbose = FALSE
  )

  vsmc <- score_signatures(vsmc, signature_sets, prefix = "vsmc")
  sig_cols <- paste0("vsmc_", names(signature_sets))
  score_df <- vsmc@meta.data[, c("vsmc_subcluster", sig_cols), drop = FALSE]
  z_scores <- as.data.frame(scale(score_df[, sig_cols, drop = FALSE]))
  z_scores$vsmc_subcluster <- score_df$vsmc_subcluster
  cluster_scores <- z_scores %>%
    dplyr::group_by(vsmc_subcluster) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(sig_cols), ~ mean(.x, na.rm = TRUE)), n_cells = dplyr::n(), .groups = "drop")
  max_sig <- names(signature_sets)[max.col(as.matrix(cluster_scores[, sig_cols, drop = FALSE]), ties.method = "first")]
  state_map <- setNames(
    dplyr::recode(
      max_sig,
      Contractile = "contractile VSMC",
      Synthetic = "synthetic VSMC",
      Inflammatory = "inflammatory VSMC",
      ECM_Remodeling = "ECM-remodeling VSMC",
      Metabolic_Stress = "metabolic-stress VSMC"
    ),
    cluster_scores$vsmc_subcluster
  )
  vsmc$vsmc_state <- unname(state_map[vsmc$vsmc_subcluster])
  vsmc$vsmc_state[is.na(vsmc$vsmc_state)] <- "mixed/other VSMC"
  vsmc$vsmc_state <- factor(vsmc$vsmc_state, levels = names(state_colors))

  cluster_scores$assigned_state <- unname(state_map[cluster_scores$vsmc_subcluster])
  write.csv(cluster_scores, file.path(out_dir, "VSMC_Subcluster_Signature_Scores.csv"), row.names = FALSE)

  props <- vsmc@meta.data %>%
    as_tibble(rownames = "cell") %>%
    dplyr::count(sample, stage, vsmc_state, name = "cells") %>%
    dplyr::group_by(sample, stage) %>%
    dplyr::mutate(total_vsmc = sum(cells), proportion = cells / total_vsmc) %>%
    dplyr::ungroup()
  write.csv(props, file.path(out_dir, "VSMC_State_Proportions_by_Sample.csv"), row.names = FALSE)

  emb <- as.data.frame(Embeddings(vsmc, "umap_vsmc_sub"))
  emb$cell <- rownames(emb)
  emb <- cbind(emb, vsmc@meta.data[emb$cell, c("sample", "stage", "vsmc_subcluster", "vsmc_state")])
  write.csv(emb, file.path(out_dir, "VSMC_Substate_Cell_Metadata.csv"), row.names = FALSE)
  saveRDS(vsmc, file.path(out_dir, "VSMC_Formal_Substate_Annotated.rds"))

  marker_path <- file.path(out_dir, "VSMC_State_Markers_FindAllMarkers.csv")
  if (!file.exists(marker_path)) {
    Idents(vsmc) <- "vsmc_state"
    markers <- FindAllMarkers(
      vsmc,
      only.pos = TRUE,
      min.pct = 0.1,
      logfc.threshold = 0.25,
      test.use = "wilcox",
      max.cells.per.ident = 1000,
      random.seed = 20260522,
      verbose = FALSE
    )
    write.csv(markers, marker_path, row.names = FALSE)
  } else {
    log_msg("Using existing VSMC marker table:", marker_path)
  }

  p_umap <- ggplot(emb, aes(vUMAP_1, vUMAP_2, color = vsmc_state)) +
    geom_point(size = 0.18, alpha = 0.65) +
    scale_color_manual(values = state_colors, drop = FALSE) +
    labs(x = "VSMC UMAP 1", y = "VSMC UMAP 2", color = "VSMC state", title = "VSMC substates") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")

  p_prop <- props %>%
    dplyr::mutate(sample = factor(sample, levels = unique(sample[order(stage, sample)]))) %>%
    ggplot(aes(sample, proportion, fill = vsmc_state)) +
    geom_col(width = 0.78) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
    scale_fill_manual(values = state_colors, drop = FALSE) +
    labs(x = NULL, y = "VSMC fraction", fill = "VSMC state", title = "Sample-level VSMC composition") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"))

  list(vsmc = vsmc, p_umap = p_umap, p_prop = p_prop, props = props)
}

run_milo_da <- function(vsmc, p_umap, p_prop) {
  log_msg("Running Milo differential abundance")
  da_path <- file.path(out_dir, "Milo_VSMC_Neighborhood_DA_Ruptured_vs_STA.csv")
  milo_path <- file.path(out_dir, "Milo_VSMC_Object.rds")
  if (file.exists(da_path) && file.exists(milo_path)) {
    log_msg("Using cached Milo results:", da_path)
    da <- read.csv(da_path, check.names = FALSE)
    milo <- readRDS(milo_path)
  } else {
    vsmc_da <- subset(vsmc, subset = stage %in% c("STA", "Ruptured_de_novo", "Ruptured_recurrent"))
    vsmc_da$status <- ifelse(vsmc_da$stage == "STA", "STA", "Ruptured")
    vsmc_da$status <- factor(vsmc_da$status, levels = c("STA", "Ruptured"))

    dims_use <- seq_len(min(20, ncol(Embeddings(vsmc_da, "harmony"))))
    sce <- as.SingleCellExperiment(vsmc_da, assay = "RNA")
    SingleCellExperiment::reducedDim(sce, "HARMONY") <- Embeddings(vsmc_da, "harmony")[colnames(vsmc_da), dims_use, drop = FALSE]
    SingleCellExperiment::reducedDim(sce, "UMAP") <- Embeddings(vsmc_da, "umap_vsmc_sub")[colnames(vsmc_da), , drop = FALSE]

    milo <- Milo(sce)
    milo <- buildGraph(milo, k = 30, d = length(dims_use), reduced.dim = "HARMONY")
    milo <- makeNhoods(milo, prop = 0.1, k = 30, d = length(dims_use), refined = TRUE, reduced_dims = "HARMONY")
    milo <- countCells(milo, samples = "sample", meta.data = as.data.frame(colData(milo)))

    design_df <- as.data.frame(colData(milo)) %>%
      as_tibble() %>%
      dplyr::select(sample, status) %>%
      dplyr::distinct() %>%
      dplyr::mutate(status = factor(status, levels = c("STA", "Ruptured"))) %>%
      as.data.frame()
    rownames(design_df) <- design_df$sample
    design_df <- design_df[colnames(nhoodCounts(milo)), , drop = FALSE]

    da <- testNhoods(
      milo,
      design = ~ status,
      design.df = design_df,
      model.contrasts = "statusRuptured",
      fdr.weighting = "none",
      reduced.dim = "HARMONY"
    )
    da <- annotateNhoods(milo, da, coldata_col = "vsmc_state")
    write.csv(da, da_path, row.names = FALSE)
    saveRDS(milo, milo_path)
  }

  da$DA_FDR <- ifelse(is.na(da$SpatialFDR), da$FDR, da$SpatialFDR)
  da$significant <- da$DA_FDR < 0.1
  write.csv(da, da_path, row.names = FALSE)

  p_milo <- ggplot(da, aes(x = vsmc_state, y = logFC, color = significant)) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.35) +
    geom_jitter(width = 0.25, height = 0, alpha = 0.55, size = 0.9) +
    scale_color_manual(values = c(`TRUE` = "#C65A46", `FALSE` = "#7F7F7F"), na.value = "#7F7F7F") +
    labs(title = "Milo DA: ruptured aneurysm vs STA", x = NULL, y = "Neighborhood logFC", color = "FDR < 0.1") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.title = element_text(face = "bold"))

  panel <- (p_umap | p_prop) / p_milo + plot_layout(heights = c(1, 0.8))
  save_png(panel, "FigureS7_VSMC_Substates_Milo_DA.png", width = 13.5, height = 10)
  list(milo = milo, da = da, p_milo = p_milo)
}

annotate_global_celltypes <- function(obj) {
  log_msg("Annotating global cell types")
  DefaultAssay(obj) <- "RNA"
  if (ncol(get_layer(obj, "RNA", "data")) == 0) {
    obj <- NormalizeData(obj, verbose = FALSE)
  }
  dims_use <- seq_len(min(30, ncol(Embeddings(obj, "harmony"))))
  obj <- FindNeighbors(obj, reduction = "harmony", dims = dims_use, graph.name = c("global_nn", "global_snn"), verbose = FALSE)
  obj <- FindClusters(obj, graph.name = "global_snn", resolution = 0.45, algorithm = 1, verbose = FALSE)
  obj$global_subcluster <- as.character(Idents(obj))
  obj <- RunUMAP(obj, reduction = "harmony", dims = dims_use, reduction.name = "umap_global_formal", reduction.key = "gUMAP_", verbose = FALSE)

  obj <- score_signatures(obj, global_marker_sets, prefix = "ct")
  ct_cols <- paste0("ct_", names(global_marker_sets))
  z <- as.data.frame(scale(obj@meta.data[, ct_cols, drop = FALSE]))
  z$global_subcluster <- obj$global_subcluster
  cluster_scores <- z %>%
    dplyr::group_by(global_subcluster) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(ct_cols), ~ mean(.x, na.rm = TRUE)), n_cells = dplyr::n(), .groups = "drop")
  best <- names(global_marker_sets)[max.col(as.matrix(cluster_scores[, ct_cols, drop = FALSE]), ties.method = "first")]
  max_score <- apply(as.matrix(cluster_scores[, ct_cols, drop = FALSE]), 1, max, na.rm = TRUE)
  assigned <- ifelse(max_score < 0.05, "Other", best)
  cluster_scores$assigned_celltype <- assigned
  write.csv(cluster_scores, file.path(out_dir, "Global_Cluster_Celltype_Signature_Scores.csv"), row.names = FALSE)

  ct_map <- setNames(cluster_scores$assigned_celltype, cluster_scores$global_subcluster)
  obj$celltype_formal <- unname(ct_map[obj$global_subcluster])
  obj$celltype_formal[is.na(obj$celltype_formal)] <- "Other"
  obj$celltype_formal <- factor(obj$celltype_formal, levels = names(celltype_colors))
  obj
}

run_cellchat <- function() {
  log_msg("Running CellChat")
  int_path <- file.path(base_dir, "01_QC_Integration", "Integrated_Seurat_Harmony.rds")
  annot_path <- file.path(out_dir, "Integrated_Formal_Celltype_Annotated.rds")
  cellchat_path <- file.path(out_dir, "CellChat_STA_Ruptured_List.rds")
  if (
    file.exists(cellchat_path) &&
      file.exists(file.path(out_dir, "CellChat_All_Communications.csv")) &&
      file.exists(file.path(out_dir, "CellChat_VSMC_Targeted_Axes.csv")) &&
      file.exists(file.path(out_dir, "FigureS8_CellChat_VSMC_Interactions.png"))
  ) {
    log_msg("Using cached CellChat results:", cellchat_path)
    return(list(
      cellchat_path = cellchat_path,
      p_focus = NULL,
      p_net = NULL,
      annotated_path = annot_path
    ))
  }
  if (file.exists(annot_path)) {
    log_msg("Using cached global celltype annotation:", annot_path)
    obj <- readRDS(annot_path)
  } else {
    obj <- readRDS(int_path)
    obj <- annotate_global_celltypes(obj)
    saveRDS(obj, annot_path)
  }
  obj$status <- dplyr::case_when(
    obj$stage == "STA" ~ "STA",
    obj$stage %in% c("Ruptured_de_novo", "Ruptured_recurrent") ~ "Ruptured",
    TRUE ~ "Unruptured"
  )

  meta_all <- obj@meta.data %>%
    as_tibble(rownames = "cell") %>%
    dplyr::filter(status %in% c("STA", "Ruptured"), celltype_formal %in% c("VSMC", "Macrophage", "EC", "Fibroblast", "T_NK"))

  meta_all <- meta_all %>%
    dplyr::group_split(status, celltype_formal) %>%
    lapply(function(df) {
      df[sample(seq_len(nrow(df)), size = min(1000, nrow(df)), replace = FALSE), , drop = FALSE]
    }) %>%
    dplyr::bind_rows()

  data_mat <- get_layer(obj, "RNA", "data")[, meta_all$cell, drop = FALSE]
  meta_use <- meta_all %>%
    dplyr::select(cell, sample, stage, status, celltype = celltype_formal) %>%
    as.data.frame()
  rownames(meta_use) <- meta_use$cell
  meta_use$samples <- meta_use$sample
  write.csv(meta_use, file.path(out_dir, "CellChat_Input_Cell_Metadata.csv"), row.names = FALSE)

  run_one <- function(status_name) {
    cells <- rownames(meta_use)[meta_use$status == status_name]
    meta <- meta_use[cells, , drop = FALSE]
    meta$celltype <- droplevels(factor(meta$celltype))
    mat <- data_mat[, cells, drop = FALSE]
    cc <- createCellChat(object = mat, meta = meta, group.by = "celltype", datatype = "RNA")
    cc@DB <- CellChatDB.human
    cc <- subsetData(cc)
    future::plan("sequential")
    cc <- identifyOverExpressedGenes(cc, do.fast = FALSE, min.cells = 10)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc, type = "triMean", population.size = TRUE, nboot = 50, seed.use = 20260522)
    cc <- filterCommunication(cc, min.cells = 10)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    cc
  }

  cellchat_list <- list(STA = run_one("STA"), Ruptured = run_one("Ruptured"))
  saveRDS(cellchat_list, cellchat_path)

  comm <- bind_rows(lapply(names(cellchat_list), function(nm) {
    df <- subsetCommunication(cellchat_list[[nm]], thresh = 0.1)
    if (is.null(df) || nrow(df) == 0) return(data.frame())
    df$status <- nm
    df
  }))
  write.csv(comm, file.path(out_dir, "CellChat_All_Communications.csv"), row.names = FALSE)

  focus_sources <- c("Macrophage", "EC", "Fibroblast")
  focus_pathways <- c("TNF", "IL1", "CCL", "TGFb", "PDGF", "COLLAGEN", "LAMININ", "FN1", "MIF", "SPP1", "MMP")
  focus <- comm %>%
    dplyr::filter(source %in% focus_sources, target == "VSMC") %>%
    dplyr::mutate(
      axis = paste(source, "->", target),
      lr_axis = paste0(ligand, "-", receptor),
      focus_hit = pathway_name %in% focus_pathways | grepl("TNF|IL1|CCL|TGF|PDGF|COLLAGEN|MMP|FN1|SPP1|MIF", pathway_name, ignore.case = TRUE)
    ) %>%
    dplyr::arrange(dplyr::desc(focus_hit), dplyr::desc(prob)) %>%
    dplyr::group_by(status) %>%
    dplyr::slice_head(n = 35) %>%
    dplyr::ungroup()
  write.csv(focus, file.path(out_dir, "CellChat_VSMC_Targeted_Axes.csv"), row.names = FALSE)

  if (nrow(focus) == 0) {
    p_focus <- plot_empty("No focused CellChat interactions passed the threshold")
  } else {
    top_lr <- focus %>%
      dplyr::group_by(lr_axis) %>%
      dplyr::summarise(max_prob = max(prob, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(max_prob)) %>%
      dplyr::slice_head(n = 25) %>%
      pull(lr_axis)
    p_focus <- focus %>%
      dplyr::filter(lr_axis %in% top_lr) %>%
      dplyr::mutate(lr_axis = factor(lr_axis, levels = rev(top_lr))) %>%
      ggplot(aes(status, lr_axis, size = prob, color = pathway_name)) +
      geom_point(alpha = 0.85) +
      facet_wrap(~ axis, scales = "free_y") +
      scale_size_continuous(range = c(1.5, 6)) +
      labs(x = NULL, y = "Ligand-receptor pair", size = "CellChat probability", color = "Pathway", title = "CellChat axes targeting VSMC") +
      theme_classic(base_size = 9) +
      theme(plot.title = element_text(face = "bold"), legend.position = "right")
  }

  net_summary <- bind_rows(lapply(names(cellchat_list), function(nm) {
    net <- cellchat_list[[nm]]@net$count
    as.data.frame(as.table(net)) %>%
      setNames(c("source", "target", "n_interactions")) %>%
      dplyr::mutate(status = nm)
  }))
  write.csv(net_summary, file.path(out_dir, "CellChat_Net_Count_Summary.csv"), row.names = FALSE)
  p_net <- net_summary %>%
    dplyr::filter(source %in% c("Macrophage", "EC", "Fibroblast", "VSMC"), target %in% c("VSMC", "Macrophage", "EC", "Fibroblast")) %>%
    ggplot(aes(source, target, fill = n_interactions)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = n_interactions), size = 2.6) +
    facet_wrap(~ status) +
    scale_fill_gradient(low = "#F4F4F4", high = "#C65A46") +
    labs(x = "Source", y = "Target", fill = "Count", title = "CellChat interaction counts") +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.title = element_text(face = "bold"))

  save_png(p_focus / p_net + plot_layout(heights = c(1.2, 0.9)), "FigureS8_CellChat_VSMC_Interactions.png", width = 13, height = 10)
  saveRDS(obj, annot_path)
  list(cellchat = cellchat_list, p_focus = p_focus, p_net = p_net, annotated = obj)
}

run_tf_activity <- function(vsmc) {
  log_msg("Running DoRothEA/VIPER TF activity")
  meta <- vsmc@meta.data %>%
    as_tibble(rownames = "cell") %>%
    dplyr::mutate(status = dplyr::case_when(
      stage == "STA" ~ "STA",
      stage %in% c("Ruptured_de_novo", "Ruptured_recurrent") ~ "Ruptured",
      TRUE ~ "Unruptured"
    ))

  group_counts <- meta %>%
    dplyr::count(sample, status, vsmc_state, name = "cells") %>%
    dplyr::filter(cells >= 25) %>%
    dplyr::mutate(tf_group = paste(gsub("_", "-", sample), status, vsmc_state, sep = "--"))
  keep_cells <- meta %>%
    dplyr::inner_join(group_counts, by = c("sample", "status", "vsmc_state")) %>%
    dplyr::select(cell, tf_group)

  vsmc$tf_group <- NA_character_
  vsmc$tf_group[keep_cells$cell] <- keep_cells$tf_group
  sub <- subset(vsmc, cells = keep_cells$cell)
  avg <- AggregateExpression(sub, assays = "RNA", group.by = "tf_group", slot = "data", return.seurat = FALSE, verbose = FALSE)$RNA

  regulons <- dorothea::dorothea_hs %>%
    dplyr::filter(confidence %in% c("A", "B", "C")) %>%
    dorothea::df2regulon()
  tf_mat <- viper::viper(as.matrix(avg), regulons, minsize = 5, verbose = FALSE)
  tf_df <- as.data.frame(tf_mat) %>%
    tibble::rownames_to_column("tf") %>%
    pivot_longer(-tf, names_to = "tf_group", values_to = "NES") %>%
    dplyr::left_join(group_counts, by = "tf_group")
  write.csv(tf_df, file.path(out_dir, "DoRothEA_VIPER_VSMC_State_TF_Activity.csv"), row.names = FALSE)

  tf_interest <- c("KLF4", "STAT1", "STAT3", "RELA", "NFKB1", "JUN", "FOS", "PPARGC1A", "FOXO1", "FOXO3", "SIRT1", "ATF3", "KLF2")
  tf_plot_df <- tf_df %>%
    dplyr::filter(tf %in% tf_interest, status %in% c("STA", "Ruptured")) %>%
    dplyr::group_by(tf, status, vsmc_state) %>%
    dplyr::summarise(mean_NES = mean(NES, na.rm = TRUE), .groups = "drop")

  if (nrow(tf_plot_df) == 0) {
    p_tf <- plot_empty("No DoRothEA/VIPER TF activities available for selected TFs")
  } else {
    p_tf <- tf_plot_df %>%
      dplyr::mutate(tf = factor(tf, levels = rev(tf_interest))) %>%
      ggplot(aes(status, tf, fill = mean_NES)) +
      geom_tile(color = "white", linewidth = 0.35) +
      facet_wrap(~ vsmc_state) +
      scale_fill_gradient2(low = "#2A6FBB", mid = "white", high = "#C65A46", midpoint = 0) +
      labs(x = NULL, y = "TF", fill = "Mean NES", title = "DoRothEA/VIPER TF activity in VSMC states") +
      theme_classic(base_size = 9) +
      theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 35, hjust = 1))
  }
  save_png(p_tf, "FigureS9_DoRothEA_VIPER_TF_Activity.png", width = 12, height = 8.5)
  list(tf_activity = tf_df, p_tf = p_tf)
}

main <- function() {
  vsmc_res <- run_vsmc_substates()
  milo_res <- run_milo_da(vsmc_res$vsmc, vsmc_res$p_umap, vsmc_res$p_prop)
  tf_res <- run_tf_activity(vsmc_res$vsmc)
  cellchat_res <- run_cellchat()

  summary_lines <- c(
    "Formal bioinformatics analysis complete.",
    paste("Output directory:", out_dir),
    "Main output files:",
    "FigureS7_VSMC_Substates_Milo_DA.png",
    "FigureS8_CellChat_VSMC_Interactions.png",
    "FigureS9_DoRothEA_VIPER_TF_Activity.png",
    "VSMC_Formal_Substate_Annotated.rds",
    "Milo_VSMC_Neighborhood_DA_Ruptured_vs_STA.csv",
    "CellChat_VSMC_Targeted_Axes.csv",
    "DoRothEA_VIPER_VSMC_State_TF_Activity.csv"
  )
  writeLines(summary_lines, file.path(out_dir, "Formal_Bioinformatics_Run_Summary.txt"))
  log_msg("Formal bioinformatics analysis complete")
}

main()
