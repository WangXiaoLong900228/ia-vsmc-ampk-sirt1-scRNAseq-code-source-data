stable_tmp_dir <- "/Users/wangxiaolong/Desktop/VSMC单细胞方法与结果/10_Formal_Bioinformatics/R_tmp"
dir.create(stable_tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(LC_ALL = "C", TMPDIR = stable_tmp_dir, TMP = stable_tmp_dir, TEMP = stable_tmp_dir)
set.seed(20260522)

base_dir <- "/Users/wangxiaolong/Desktop/VSMC单细胞方法与结果"
lib_dir <- file.path(base_dir, "R_formal_packages")
.libPaths(c(lib_dir, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(data.table)
})

out_dir <- file.path(base_dir, "10_Formal_Bioinformatics")
scenic_dir <- "/Users/wangxiaolong/Desktop/IA_SCENIC_Run_targeted"
db_dir <- "/Users/wangxiaolong/Desktop/IA_SCENIC_databases_mc9nr"
dir.create(scenic_dir, recursive = TRUE, showWarnings = FALSE)
data("motifAnnotations_hgnc", package = "RcisTarget")

patch_scenic_motif_annotations <- function() {
  ns <- asNamespace("SCENIC")
  patched <- function(scenicOptions) {
    dbAnnotFiles <- scenicOptions@settings$db_annotFiles
    if (!is.null(dbAnnotFiles)) {
      motifAnnotations <- NULL
      for (annotPath in dbAnnotFiles) {
        motifAnnot <- data.table::fread(annotPath)
        motifAnnot$annotationSource <- factor(motifAnnot$annotationSource)
        colnames(motifAnnot)[1] <- "motif"
        levels(motifAnnot$annotationSource) <- c(levels(motifAnnot$annotationSource),
          c("directAnnotation", "inferredBy_Orthology",
            "inferredBy_MotifSimilarity", "inferredBy_MotifSimilarity_n_Orthology"))
        motifAnnotations <- rbind(motifAnnotations, motifAnnot)
      }
    } else {
      org <- SCENIC::getDatasetInfo(scenicOptions, "org")
      if (is.na(org)) stop("Please provide an organism (scenicOptions@inputDatasetInfo$org).")
      if (!org %in% c("hgnc", "mgi", "dmel")) stop("Organism not recognized (scenicOptions@inputDatasetInfo$org).")
      motifAnnotName <- switch(org,
        hgnc = "motifAnnotations_hgnc",
        mgi = "motifAnnotations_mgi",
        dmel = "motifAnnotations_dmel"
      )
      if (!is.null(scenicOptions@settings$db_mcVersion)) {
        if (scenicOptions@settings$db_mcVersion == "v8") motifAnnotName <- paste0(motifAnnotName, "_v8")
      }
      data(list = motifAnnotName, package = "RcisTarget", verbose = FALSE)
      if (exists(motifAnnotName, inherits = FALSE)) {
        motifAnnotations <- get(motifAnnotName)
      } else if (exists("motifAnnotations", inherits = FALSE)) {
        motifAnnotations <- get("motifAnnotations")
      } else {
        v9_name <- paste0(motifAnnotName, "_v9")
        data(list = v9_name, package = "RcisTarget", verbose = FALSE)
        motifAnnotations <- get(v9_name)
      }
    }
    motifAnnotations
  }
  unlockBinding("getDbAnnotations", ns)
  assign("getDbAnnotations", patched, envir = ns)
  lockBinding("getDbAnnotations", ns)
}

patch_scenic_motif_annotations()

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")
  flush.console()
}

ensure_tempdir <- function() {
  dir.create(stable_tmp_dir, recursive = TRUE, showWarnings = FALSE)
  td <- tempdir()
  if (!dir.exists(td)) dir.create(td, recursive = TRUE, showWarnings = FALSE)
  invisible(td)
}

save_png <- function(plot, filename, width = 11, height = 7.5, dpi = 300) {
  ensure_tempdir()
  target <- file.path(out_dir, filename)
  tmp <- tempfile(pattern = "scenic_plot_", tmpdir = "/tmp", fileext = ".png")
  ggsave(tmp, plot, width = width, height = height, dpi = dpi, bg = "white")
  ok <- file.copy(tmp, target, overwrite = TRUE)
  unlink(tmp)
  if (!ok) stop("Could not copy plot to ", target)
}

get_feather_genes <- function(db_file) {
  rf <- arrow::ReadableFile$create(db_file)
  fr <- arrow::FeatherReader$create(rf)
  setdiff(names(fr), "features")
}

resume_scenic_2_from_cached_motif_enrichment <- function(scenic_options, minGenes = 20) {
  log_msg("Resuming SCENIC step 2 from cached motif enrichment with per-database gene filtering")
  nCores <- getSettings(scenic_options, "nCores")
  db_names <- getDatabases(scenic_options)
  tf_modules <- loadInt(scenic_options, "tfModules_forEnrichment")
  motif_enrichment <- loadInt(scenic_options, "motifEnrichment_full")

  motif_enrichment_self <- motif_enrichment[which(motif_enrichment$TFinDB != ""), , drop = FALSE]
  if (nrow(motif_enrichment_self) == 0) {
    stop("None of the co-expression modules present enrichment of the TF motif: There are no regulons.")
  }

  met_by_db <- split(motif_enrichment_self, motif_enrichment_self$motifDb)
  for (db in names(met_by_db)) {
    met <- met_by_db[[db]]
    met <- split(met, factor(met$highlightedTFs))
    met <- lapply(met, function(x) {
      data.table::rbindlist(lapply(split(x, x$motif), function(y) y[which.max(y$NES), ]))
    })
    met_by_db[[db]] <- data.table::rbindlist(met)
  }
  motif_enrichment_self <- data.table::rbindlist(met_by_db)
  log_msg("Pruning cached self-motif enrichments:", nrow(motif_enrichment_self), "rows")

  motif_db_keys <- unique(as.character(motif_enrichment_self$motifDb))
  motif_enrichment_self_w_genes <- lapply(names(db_names), function(motif_db_name) {
    db_path <- db_names[[motif_db_name]]
    motif_db_key <- motif_db_keys[endsWith(basename(db_path), motif_db_keys)]
    if (length(motif_db_key) != 1) {
      motif_db_key <- motif_db_name
    }
    db_genes <- get_feather_genes(db_path)
    tf_modules_db <- lapply(tf_modules, function(module_genes) {
      intersect(as.character(module_genes), db_genes)
    })
    tf_modules_db <- tf_modules_db[lengths(tf_modules_db) >= minGenes]
    results_db <- motif_enrichment_self[motif_enrichment_self$motifDb == motif_db_key, , drop = FALSE]
    results_db <- results_db[results_db$geneSet %in% names(tf_modules_db), , drop = FALSE]
    if (nrow(results_db) == 0 || length(tf_modules_db) == 0) return(NULL)

    all_genes_db <- sort(unique(unlist(tf_modules_db, use.names = FALSE)))
    log_msg("Pruning database", motif_db_key, "with", length(all_genes_db), "genes and", nrow(results_db), "motif rows")
    ranking <- RcisTarget::importRankings(db_path, columns = all_genes_db)
    RcisTarget::addSignificantGenes(
      resultsTable = results_db,
      geneSets = tf_modules_db,
      rankings = ranking,
      plotCurve = FALSE,
      maxRank = 5000,
      method = "aprox",
      nMean = 100,
      nCores = nCores
    )
  })
  motif_enrichment_self_w_genes <- data.table::rbindlist(Filter(Negate(is.null), motif_enrichment_self_w_genes))
  if (nrow(motif_enrichment_self_w_genes) == 0) {
    stop("No motifs remained after per-database pruning. Check motifDb/database name matching.")
  }
  saveRDS(motif_enrichment_self_w_genes, file = getIntName(scenic_options, "motifEnrichment_selfMotifs_wGenes"))
  log_msg("Motifs supporting regulons after pruning:", nrow(motif_enrichment_self_w_genes))

  if (!file.exists("output")) dir.create("output")
  write.table(
    motif_enrichment_self_w_genes,
    file = getOutName(scenic_options, "s2_motifEnrichment"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  motif_enrichment_as_incid_list <- apply(motif_enrichment_self_w_genes, 1, function(one_motif_row) {
    genes <- strsplit(one_motif_row["enrichedGenes"], ";")[[1]]
    one_motif_row <- data.frame(rbind(one_motif_row), stringsAsFactors = FALSE)
    data.frame(
      one_motif_row[rep(1, length(genes)), c("NES", "motif", "highlightedTFs", "TFinDB", "geneSet", "motifDb")],
      genes,
      stringsAsFactors = FALSE
    )
  })
  motif_enrichment_as_incid_list <- data.table::rbindlist(motif_enrichment_as_incid_list)
  colnames(motif_enrichment_as_incid_list)[which(colnames(motif_enrichment_as_incid_list) == "highlightedTFs")] <- "TF"
  colnames(motif_enrichment_as_incid_list)[which(colnames(motif_enrichment_as_incid_list) == "TFinDB")] <- "annot"
  colnames(motif_enrichment_as_incid_list)[which(colnames(motif_enrichment_as_incid_list) == "genes")] <- "gene"
  motif_enrichment_as_incid_list <- data.frame(motif_enrichment_as_incid_list, stringsAsFactors = FALSE)

  regulon_targets_info <- lapply(split(motif_enrichment_as_incid_list, motif_enrichment_as_incid_list$TF), function(tf_targets) {
    tf_table <- as.data.frame(do.call(rbind, lapply(split(tf_targets, tf_targets$gene), function(enr_one_gene) {
      high_conf_annot <- "**" %in% enr_one_gene$annot
      enr_one_gene_by_annot <- enr_one_gene
      if (high_conf_annot) enr_one_gene_by_annot <- enr_one_gene_by_annot[which(enr_one_gene$annot == "**"), ]
      best_motif <- which.max(enr_one_gene_by_annot$NES)
      tf <- unique(enr_one_gene$TF)
      cbind(
        TF = tf,
        gene = unique(enr_one_gene$gene),
        highConfAnnot = high_conf_annot,
        nMotifs = nrow(enr_one_gene),
        bestMotif = as.character(enr_one_gene_by_annot[best_motif, "motif"]),
        NES = as.numeric(enr_one_gene_by_annot[best_motif, "NES"]),
        motifDb = as.character(enr_one_gene_by_annot[best_motif, "motifDb"]),
        coexModule = gsub(paste0(tf, "_"), "", as.character(enr_one_gene_by_annot[best_motif, "geneSet"]), fixed = TRUE)
      )
    })), stringsAsFactors = FALSE)
    tf_table[order(tf_table$NES, decreasing = TRUE), ]
  })
  regulon_targets_info <- data.table::rbindlist(regulon_targets_info)

  corr_mat <- loadInt(scenic_options, "corrMat", ifNotExists = "null")
  if (!is.null(corr_mat)) {
    regulon_targets_info$spearCor <- NA_real_
    for (tf in unique(regulon_targets_info$TF)) {
      rows <- which(regulon_targets_info$TF == tf)
      genes <- unlist(regulon_targets_info[rows, "gene"])
      present <- genes %in% colnames(corr_mat)
      if (tf %in% rownames(corr_mat) && any(present)) {
        regulon_targets_info[rows[present], "spearCor"] <- corr_mat[tf, genes[present]]
      }
    }
  } else {
    warning("It was not possible to add the correlation to the regulonTargetsInfo table.")
  }

  link_list <- loadInt(scenic_options, "genie3ll", ifNotExists = "null")
  if (!is.null(link_list) && ("weight" %in% colnames(link_list))) {
    if (data.table::is.data.table(link_list)) link_list <- as.data.frame(link_list)
    unique_pairs <- nrow(unique(link_list[, c("TF", "Target")]))
    if (unique_pairs == nrow(link_list)) {
      link_list <- link_list[which(link_list$weight >= getSettings(scenic_options, "modules/weightThreshold")), ]
      rownames(link_list) <- paste(link_list$TF, link_list$Target, sep = "__")
      regulon_targets_info <- cbind(
        regulon_targets_info,
        CoexWeight = link_list[paste(regulon_targets_info$TF, regulon_targets_info$gene, sep = "__"), "weight"]
      )
    } else {
      warning("There are duplicated regulator-target pairs in the co-expression link list.")
    }
  } else {
    warning("It was not possible to add the weight to the regulonTargetsInfo table.")
  }

  saveRDS(regulon_targets_info, file = getIntName(scenic_options, "regulonTargetsInfo"))
  write.table(
    regulon_targets_info,
    file = getOutName(scenic_options, "s2_regulonTargetsInfo"),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )

  regulon_targets_by_annot <- split(regulon_targets_info, regulon_targets_info$highConfAnnot)
  regulons <- list()
  if (!is.null(regulon_targets_by_annot[["TRUE"]])) {
    regulons <- lapply(split(regulon_targets_by_annot[["TRUE"]], regulon_targets_by_annot[["TRUE"]][, "TF"]), function(x) {
      sort(as.character(unlist(x[, "gene"])))
    })
  }
  regulons_extended <- list()
  if (!is.null(regulon_targets_by_annot[["FALSE"]])) {
    regulons_extended <- lapply(split(regulon_targets_by_annot[["FALSE"]], regulon_targets_by_annot[["FALSE"]][, "TF"]), function(x) {
      unname(unlist(x[, "gene"]))
    })
    regulons_extended <- setNames(lapply(names(regulons_extended), function(tf) {
      sort(unique(c(regulons[[tf]], unlist(regulons_extended[[tf]]))))
    }), names(regulons_extended))
    names(regulons_extended) <- paste(names(regulons_extended), "_extended", sep = "")
  }
  regulons <- c(regulons, regulons_extended)
  regulons <- regulons[lengths(regulons) >= minGenes]
  if (length(regulons) == 0) stop("No regulons passed the minimum gene threshold after pruning.")

  saveRDS(regulons, file = getIntName(scenic_options, "regulons"))
  incid_list <- reshape2::melt(regulons)
  incid_mat <- table(incid_list[, 2], incid_list[, 1])
  saveRDS(incid_mat, file = getIntName(scenic_options, "regulons_incidMat"))
  scenic_options@status$current <- 2
  invisible(scenic_options)
}

get_layer <- function(obj, assay = "RNA", layer = "counts") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
}

state_colors <- c(
  "contractile VSMC" = "#2A6FBB",
  "synthetic VSMC" = "#C65A46",
  "inflammatory VSMC" = "#A33E8A",
  "ECM-remodeling VSMC" = "#3A8F62",
  "metabolic-stress VSMC" = "#D19A2E",
  "mixed/other VSMC" = "#7F7F7F"
)

main <- function() {
  summary_path <- file.path(out_dir, "SCENIC_VSMC_Run_Summary.txt")
  auc_path <- file.path(out_dir, "SCENIC_VSMC_RegulonAUC_ByCell.csv")
  if (file.exists(auc_path) && file.exists(file.path(out_dir, "FigureS10_SCENIC_Regulons.png"))) {
    log_msg("Using cached SCENIC outputs")
    writeLines(c("SCENIC outputs already exist.", auc_path), summary_path)
    return(invisible(NULL))
  }

  dbs <- c(
    "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",
    "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
  )
  missing_dbs <- dbs[!file.exists(file.path(db_dir, dbs))]
  if (length(missing_dbs) > 0) stop("Missing SCENIC database files: ", paste(missing_dbs, collapse = ", "))

  log_msg("Loading VSMC substate object")
  vsmc <- readRDS(file.path(out_dir, "VSMC_Formal_Substate_Annotated.rds"))
  vsmc$status <- dplyr::case_when(
    vsmc$stage == "STA" ~ "STA",
    vsmc$stage %in% c("Ruptured_de_novo", "Ruptured_recurrent") ~ "Ruptured",
    TRUE ~ "Unruptured"
  )
  meta <- vsmc@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::filter(status %in% c("STA", "Ruptured")) %>%
    dplyr::group_split(status, vsmc_state) %>%
    lapply(function(df) {
      df[sample(seq_len(nrow(df)), size = min(300, nrow(df)), replace = FALSE), , drop = FALSE]
    }) %>%
    dplyr::bind_rows()
  write.csv(meta, file.path(out_dir, "SCENIC_Input_Cell_Metadata.csv"), row.names = FALSE)

  log_msg("Preparing expression matrix for", nrow(meta), "cells")
  expr_mat <- get_layer(vsmc, "RNA", "counts")[, meta$cell, drop = FALSE]
  expr_mat <- as.matrix(expr_mat)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(scenic_dir)

  log_msg("Initializing SCENIC")
  scenic_options <- initializeScenic(
    org = "hgnc",
    dbDir = db_dir,
    dbs = dbs,
    datasetTitle = "IA_VSMC_SCENIC",
    nCores = 4,
    dbIndexCol = "features"
  )
  saveRDS(scenic_options, file.path(out_dir, "SCENIC_VSMC_Options_initial.rds"))

  log_msg("Filtering genes")
  genes_kept <- geneFiltering(
    expr_mat,
    scenicOptions = scenic_options,
    minCountsPerGene = 3 * 0.01 * ncol(expr_mat),
    minSamples = 0.01 * ncol(expr_mat)
  )
  expr_filt <- expr_mat[genes_kept, , drop = FALSE]

  marker_path <- file.path(out_dir, "VSMC_State_Markers_FindAllMarkers.csv")
  marker_genes <- character()
  if (file.exists(marker_path)) {
    marker_tbl <- read.csv(marker_path)
    fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(marker_tbl))[1]
    marker_tbl <- marker_tbl %>%
      dplyr::filter(!is.na(gene)) %>%
      dplyr::arrange(cluster, p_val_adj, dplyr::desc(.data[[fc_col]])) %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_head(n = 150) %>%
      dplyr::ungroup()
    marker_genes <- unique(marker_tbl$gene)
  }
  mechanism_genes <- unique(c(
    "KLF4", "STAT1", "STAT3", "RELA", "NFKB1", "JUN", "FOS", "PPARGC1A",
    "FOXO1", "FOXO3", "SIRT1", "ATF3", "KLF2", "PRKAA1", "PRKAA2",
    "ACTA2", "MYH11", "TAGLN", "CNN1", "FN1", "VCAN", "SPP1", "CCL2",
    "IL6", "NFKBIA", "TNFAIP3", "COL1A1", "COL1A2", "COL3A1", "MMP2",
    "MMP9", "MMP14", "HMOX1", "DDIT3", "SOD2", "TXNIP"
  ))
  db_tfs <- tryCatch(SCENIC::getDbTfs(scenic_options), error = function(e) character())
  target_genes <- unique(c(db_tfs, marker_genes, mechanism_genes))
  target_genes <- intersect(target_genes, rownames(expr_filt))
  db_gene_intersection <- Reduce(intersect, lapply(file.path(db_dir, dbs), get_feather_genes))
  target_genes <- intersect(target_genes, db_gene_intersection)
  expr_filt <- expr_filt[target_genes, , drop = FALSE]
  write.csv(data.frame(gene = rownames(expr_filt)), file.path(out_dir, "SCENIC_Targeted_Genes.csv"), row.names = FALSE)
  saveRDS(genes_kept, file.path(out_dir, "SCENIC_VSMC_GenesKept.rds"))
  log_msg("Genes retained for targeted SCENIC/GENIE3:", nrow(expr_filt))

  corr_path <- file.path(scenic_dir, "int", "1.2_corrMat.Rds")
  if (file.exists(corr_path)) {
    log_msg("Using cached SCENIC correlation matrix:", corr_path)
  } else {
    log_msg("Running SCENIC correlation")
    runCorrelation(expr_filt, scenic_options)
  }

  genie3_link_path <- file.path(scenic_dir, "int", "1.4_GENIE3_linkList.Rds")
  if (file.exists(genie3_link_path)) {
    log_msg("Using cached GENIE3 link list:", genie3_link_path)
  } else {
    log_msg("Running GENIE3")
    runGenie3(expr_filt, scenic_options, nParts = 20, resumePreviousRun = TRUE)
  }

  log_msg("Building co-expression modules")
  ensure_tempdir()
  scenic_options <- runSCENIC_1_coexNetwork2modules(scenic_options)

  log_msg("Creating regulons with RcisTarget")
  ensure_tempdir()
  regulons_path <- file.path(scenic_dir, getIntName(scenic_options, "regulons"))
  motif_enrichment_path <- file.path(scenic_dir, getIntName(scenic_options, "motifEnrichment_full"))
  if (file.exists(regulons_path)) {
    log_msg("Using cached SCENIC regulons:", regulons_path)
    scenic_options@status$current <- 2
  } else if (file.exists(motif_enrichment_path)) {
    scenic_options <- resume_scenic_2_from_cached_motif_enrichment(scenic_options, minGenes = 20)
  } else {
    scenic_options <- runSCENIC_2_createRegulons(scenic_options, minGenes = 20, dbIndexCol = "features")
  }

  log_msg("Scoring cells with AUCell")
  ensure_tempdir()
  scenic_options@settings$nCores <- 1
  scenic_options <- runSCENIC_3_scoreCells(
    scenic_options,
    expr_filt,
    skipBinaryThresholds = TRUE,
    skipHeatmap = TRUE,
    skipTsne = TRUE
  )
  saveRDS(scenic_options, file.path(out_dir, "SCENIC_VSMC_Options_final.rds"))

  regulon_auc <- loadInt(scenic_options, "aucell_regulonAUC")
  auc <- as.data.frame(t(AUCell::getAUC(regulon_auc)))
  auc <- tibble::rownames_to_column(auc, "cell")
  write.csv(auc, auc_path, row.names = FALSE)

  auc_long <- auc %>%
    tidyr::pivot_longer(-cell, names_to = "regulon", values_to = "AUC") %>%
    dplyr::left_join(meta %>% dplyr::select(cell, sample, stage, status, vsmc_state), by = "cell")
  write.csv(auc_long, file.path(out_dir, "SCENIC_VSMC_RegulonAUC_Long.csv"), row.names = FALSE)

  state_summary <- auc_long %>%
    dplyr::group_by(regulon, status, vsmc_state) %>%
    dplyr::summarise(mean_AUC = mean(AUC, na.rm = TRUE), .groups = "drop")
  write.csv(state_summary, file.path(out_dir, "SCENIC_VSMC_RegulonAUC_ByState.csv"), row.names = FALSE)

  tf_interest <- c("KLF4", "STAT1", "STAT3", "RELA", "NFKB1", "JUN", "FOS", "PPARGC1A", "FOXO1", "FOXO3", "SIRT1", "ATF3", "KLF2")
  plot_df <- state_summary %>%
    dplyr::mutate(tf = sub("\\(.*", "", regulon)) %>%
    dplyr::filter(tf %in% tf_interest)
  write.csv(plot_df, file.path(out_dir, "SCENIC_VSMC_TF_Interest_Regulons.csv"), row.names = FALSE)

  if (nrow(plot_df) == 0) {
    top_regs <- state_summary %>%
      dplyr::group_by(regulon) %>%
      dplyr::summarise(dynamic_range = max(mean_AUC, na.rm = TRUE) - min(mean_AUC, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(dynamic_range)) %>%
      dplyr::slice_head(n = 20) %>%
      dplyr::pull(regulon)
    plot_df <- state_summary %>% dplyr::filter(regulon %in% top_regs)
  }

  p <- plot_df %>%
    dplyr::mutate(
      regulon = factor(regulon, levels = rev(unique(regulon))),
      vsmc_state = factor(vsmc_state, levels = names(state_colors))
    ) %>%
    ggplot(aes(status, regulon, fill = mean_AUC)) +
    geom_tile(color = "white", linewidth = 0.35) +
    facet_wrap(~ vsmc_state, scales = "free_y") +
    scale_fill_gradient(low = "#F4F4F4", high = "#A33E8A") +
    labs(x = NULL, y = "SCENIC regulon", fill = "Mean AUC", title = "SCENIC regulon activity in VSMC states") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 35, hjust = 1))
  save_png(p, "FigureS10_SCENIC_Regulons.png", width = 13, height = 9)

  writeLines(c(
    "SCENIC VSMC regulon analysis complete.",
    paste("Cells used:", nrow(meta)),
    paste("Genes retained in targeted SCENIC/GENIE3:", nrow(expr_filt)),
    paste("Regulons scored:", ncol(auc) - 1),
    "Target set: RcisTarget human TFs plus VSMC state marker genes and prespecified AMPK/SIRT1-inflammatory-remodeling genes.",
    paste("Databases:", paste(dbs, collapse = "; "))
  ), summary_path)
  log_msg("SCENIC analysis complete")
}

main()
