#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

project_dir <- "/Users/wangxiaolong/Desktop/VSMC单细胞方法与结果/11_External_Cohort_Validation"
raw_dir <- file.path(project_dir, "raw")
processed_dir <- file.path(project_dir, "processed")
figure_dir <- file.path(project_dir, "figures")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

signature_sets <- list(
  AMPK_SIRT1 = c("SIRT1", "PRKAA1", "PRKAA2", "PPARGC1A", "TFAM", "FOXO1", "FOXO3", "NAMPT", "SOD2", "CAT", "NFE2L2"),
  Contractile = c("ACTA2", "TAGLN", "MYH11", "CNN1", "MYLK", "SMTN", "LMOD1", "TPM2", "ACTG2", "CALD1", "DES", "MYL9"),
  Inflammatory = c("IL1B", "CCL2", "CXCL8", "CXCL2", "CXCL3", "IL6", "NFKBIA", "TNFAIP3", "STAT1", "STAT3", "JUN", "FOS", "ICAM1", "VCAM1", "HLA-DRA"),
  ECM_Remodeling = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "FN1", "VCAN", "LUM", "DCN", "MMP2", "MMP9", "MMP14", "TIMP1", "POSTN", "LOX")
)

key_genes <- unique(c(
  signature_sets$AMPK_SIRT1,
  "ACTA2", "MYH11", "TAGLN", "KLF4", "SPP1",
  "IL1B", "CCL2", "CXCL8", "NFKBIA", "STAT1", "JUN", "FOS",
  "COL1A1", "COL3A1", "FN1", "MMP2", "MMP9"
))

message_time <- function(...) {
  message(format(Sys.time(), "%H:%M:%S"), " | ", ...)
}

clean_geo_values <- function(x) {
  x <- gsub('^"|"$', "", x)
  x <- gsub('\\"', '"', x, fixed = TRUE)
  x
}

parse_geo_line <- function(lines, tag) {
  hit <- grep(paste0("^", tag, "\t"), lines, value = TRUE)
  if (!length(hit)) return(character())
  clean_geo_values(strsplit(hit[1], "\t", fixed = TRUE)[[1]][-1])
}

parse_geo_metadata <- function(lines) {
  geo <- parse_geo_line(lines, "!Sample_geo_accession")
  title <- parse_geo_line(lines, "!Sample_title")
  source <- parse_geo_line(lines, "!Sample_source_name_ch1")

  meta <- data.frame(
    geo_accession = geo,
    title = title,
    source_name = source,
    stringsAsFactors = FALSE
  )

  characteristic_lines <- grep("^!Sample_characteristics_ch1\t", lines, value = TRUE)
  for (line in characteristic_lines) {
    values <- clean_geo_values(strsplit(line, "\t", fixed = TRUE)[[1]][-1])
    keys <- str_trim(sub(":.*$", "", values))
    vals <- str_trim(sub("^[^:]+:\\s*", "", values))
    key <- keys[which.max(table(keys))]
    if (!length(key) || is.na(key) || key == values[1]) {
      key <- paste0("characteristic_", ncol(meta))
      vals <- values
    }
    key <- make.names(tolower(key))
    if (key %in% names(meta)) key <- make.unique(c(names(meta), key))[length(names(meta)) + 1]
    meta[[key]] <- vals
  }

  meta
}

read_geo_series_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  begin <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  if (!length(begin) || !length(end) || end <= begin + 1) {
    stop("No expression table found in ", basename(path))
  }
  matrix_text <- paste(lines[(begin + 1):(end - 1)], collapse = "\n")
  expr_dt <- fread(text = matrix_text, data.table = FALSE, check.names = FALSE)
  stopifnot("ID_REF" %in% names(expr_dt))
  ids <- as.character(expr_dt$ID_REF)
  expr_mat <- as.matrix(expr_dt[, setdiff(names(expr_dt), "ID_REF"), drop = FALSE])
  suppressWarnings(storage.mode(expr_mat) <- "numeric")
  rownames(expr_mat) <- ids
  list(expr = expr_mat, meta = parse_geo_metadata(lines))
}

clean_symbol <- function(x) {
  x <- as.character(x)
  x <- str_split_fixed(x, "///", 2)[, 1]
  x <- str_trim(x)
  x <- toupper(x)
  x[x %in% c("", "---", "NA", "N/A", "NULL")] <- NA_character_
  x
}

read_platform_annotation <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  begin <- grep("^!platform_table_begin", lines)
  end <- grep("^!platform_table_end", lines)
  if (!length(begin) || !length(end) || end <= begin + 1) {
    stop("No platform table found in ", basename(path))
  }
  platform_text <- paste(lines[(begin + 1):(end - 1)], collapse = "\n")
  annot <- fread(text = platform_text, data.table = FALSE, check.names = FALSE, quote = "")
  symbol_col <- grep("^Gene symbol$", names(annot), value = TRUE)[1]
  if (is.na(symbol_col)) stop("Cannot find Gene symbol column in ", basename(path))
  out <- data.frame(
    ID_REF = as.character(annot$ID),
    symbol = clean_symbol(annot[[symbol_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$symbol), ]
  out <- out[!duplicated(out$ID_REF), ]
  out
}

collapse_matrix_to_symbols <- function(expr, annot, id_col = "ID_REF") {
  row_ids <- rownames(expr)
  map <- annot %>%
    transmute(ID_REF = as.character(.data[[id_col]]), symbol = toupper(symbol)) %>%
    filter(!is.na(symbol), symbol != "")
  idx <- match(row_ids, map$ID_REF)
  keep <- !is.na(idx)
  expr <- expr[keep, , drop = FALSE]
  symbols <- map$symbol[idx[keep]]
  collapsed_sum <- rowsum(expr, group = symbols, reorder = FALSE, na.rm = TRUE)
  counts <- as.numeric(table(factor(symbols, levels = rownames(collapsed_sum))))
  sweep(collapsed_sum, 1, counts, "/")
}

maybe_log2_microarray <- function(expr) {
  q99 <- suppressWarnings(as.numeric(quantile(expr, 0.99, na.rm = TRUE)))
  mx <- suppressWarnings(max(expr, na.rm = TRUE))
  if (is.finite(q99) && (q99 > 100 || mx > 1000)) {
    expr <- log2(expr + 1)
  }
  expr
}

make_human_ensembl_symbol_map <- function() {
  out_path <- file.path(processed_dir, "human_ensembl_symbol_map.tsv")
  if (file.exists(out_path) && file.info(out_path)$size > 1000) {
    return(fread(out_path, data.table = FALSE))
  }

  gene2ensembl <- file.path(raw_dir, "gene2ensembl.gz")
  gene_info <- file.path(raw_dir, "Homo_sapiens.gene_info.gz")
  human_gene2ensembl <- file.path(processed_dir, "human_gene2ensembl.tsv")
  human_gene_info <- file.path(processed_dir, "human_gene_info.tsv")

  if (!file.exists(human_gene2ensembl)) {
    cmd <- sprintf(
      "gzip -dc %s | awk 'BEGIN{FS=OFS=\"\\t\"} $1==9606 && $3!=\"-\" {print $2,$3}' | sort -u > %s",
      shQuote(gene2ensembl), shQuote(human_gene2ensembl)
    )
    system(cmd)
  }
  if (!file.exists(human_gene_info)) {
    cmd <- sprintf(
      "gzip -dc %s | awk 'BEGIN{FS=OFS=\"\\t\"} $1==9606 {print $2,$3}' | sort -u > %s",
      shQuote(gene_info), shQuote(human_gene_info)
    )
    system(cmd)
  }

  e2g <- fread(human_gene2ensembl, col.names = c("GeneID", "ensembl_id"), data.table = FALSE)
  ginfo <- fread(human_gene_info, col.names = c("GeneID", "symbol"), data.table = FALSE)
  map <- e2g %>%
    inner_join(ginfo, by = "GeneID") %>%
    mutate(
      ensembl_id = sub("\\..*$", "", ensembl_id),
      symbol = clean_symbol(symbol)
    ) %>%
    filter(!is.na(symbol), symbol != "") %>%
    distinct(ensembl_id, symbol)
  fwrite(map, out_path, sep = "\t")
  map
}

read_gse122897_counts <- function() {
  message_time("Reading GSE122897 RNA-seq counts")
  counts_path <- file.path(raw_dir, "GSE122897_readCounts_raw.txt.gz")
  counts <- fread(counts_path, data.table = FALSE, check.names = FALSE)
  names(counts)[1] <- "ensembl_id"
  counts$ensembl_id <- sub("\\..*$", "", counts$ensembl_id)
  mat <- as.matrix(counts[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- counts$ensembl_id

  map <- make_human_ensembl_symbol_map()
  idx <- match(rownames(mat), map$ensembl_id)
  keep <- !is.na(idx)
  mat <- mat[keep, , drop = FALSE]
  symbols <- map$symbol[idx[keep]]
  mat_sum <- rowsum(mat, group = symbols, reorder = FALSE, na.rm = TRUE)
  gene_counts <- as.numeric(table(factor(symbols, levels = rownames(mat_sum))))
  mat_mean <- sweep(mat_sum, 1, gene_counts, "/")

  series <- read_geo_series_matrix(file.path(raw_dir, "GSE122897_series_matrix.txt.gz"))
  meta <- series$meta
  meta$sample_id <- meta$title
  meta$group <- case_when(
    str_detect(str_to_lower(meta$tissue.status), "unknown") ~ "Unknown",
    str_detect(str_to_lower(meta$tissue.status), "intracranial cortical artery") ~ "Control",
    str_detect(str_to_lower(meta$tissue.status), "unruptured") ~ "Unruptured",
    str_detect(str_to_lower(meta$tissue.status), "ruptured") ~ "Ruptured",
    TRUE ~ "Unknown"
  )
  meta$dataset <- "GSE122897"
  meta$pair_id <- NA_character_
  meta <- meta %>% filter(sample_id %in% colnames(mat_mean), group != "Unknown")
  mat_mean <- mat_mean[, meta$sample_id, drop = FALSE]

  dge <- DGEList(counts = round(mat_mean), group = meta$group)
  keep_genes <- filterByExpr(dge, group = meta$group, min.count = 5)
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  log_cpm <- cpm(dge, log = TRUE, prior.count = 2)

  log_cpm <- log_cpm[, meta$sample_id, drop = FALSE]
  list(expr = log_cpm, meta = meta)
}

read_microarray_dataset <- function(dataset, matrix_file, annot_file, platform, group_fun) {
  message_time("Reading ", dataset, " microarray matrix")
  series <- read_geo_series_matrix(file.path(raw_dir, matrix_file))
  annot <- read_platform_annotation(file.path(raw_dir, annot_file))
  expr <- maybe_log2_microarray(series$expr)
  expr <- collapse_matrix_to_symbols(expr, annot)
  expr <- normalizeBetweenArrays(expr, method = "quantile")

  meta <- series$meta
  meta$sample_id <- meta$geo_accession
  meta$group <- group_fun(meta)
  meta$dataset <- dataset
  meta$platform <- platform
  if (!"pair_id" %in% names(meta)) meta$pair_id <- NA_character_
  meta <- meta %>% filter(sample_id %in% colnames(expr), !is.na(group), group != "Unknown")
  expr <- expr[, meta$sample_id, drop = FALSE]
  list(expr = expr, meta = meta)
}

zscore_rows <- function(expr) {
  expr <- expr[rowSums(is.na(expr)) < ncol(expr), , drop = FALSE]
  sds <- apply(expr, 1, sd, na.rm = TRUE)
  expr <- expr[is.finite(sds) & sds > 0, , drop = FALSE]
  z <- t(scale(t(expr)))
  z[!is.finite(z)] <- NA_real_
  z
}

score_signatures <- function(expr, meta, dataset) {
  rownames(expr) <- toupper(rownames(expr))
  z <- zscore_rows(expr)
  scores <- lapply(names(signature_sets), function(sig) {
    genes <- intersect(signature_sets[[sig]], rownames(z))
    if (length(genes) < 2) {
      return(data.frame())
    }
    data.frame(
      dataset = dataset,
      sample_id = colnames(z),
      signature = sig,
      score = colMeans(z[genes, , drop = FALSE], na.rm = TRUE),
      present_genes = length(genes),
      total_genes = length(signature_sets[[sig]]),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  gene_z <- z[intersect(key_genes, rownames(z)), , drop = FALSE]
  gene_long <- as.data.frame(t(gene_z)) %>%
    tibble::rownames_to_column("sample_id") %>%
    pivot_longer(-sample_id, names_to = "gene", values_to = "z_expression") %>%
    mutate(dataset = dataset)

  meta_for_join <- meta %>% select(-any_of("dataset"))
  scores <- scores %>% left_join(meta_for_join, by = "sample_id")
  gene_long <- gene_long %>% left_join(meta_for_join, by = "sample_id")
  coverage <- data.frame(
    dataset = dataset,
    signature = names(signature_sets),
    present_genes = sapply(signature_sets, function(g) length(intersect(g, rownames(z)))),
    total_genes = sapply(signature_sets, length),
    present_gene_list = sapply(signature_sets, function(g) paste(intersect(g, rownames(z)), collapse = ";")),
    stringsAsFactors = FALSE
  )

  list(scores = scores, genes = gene_long, coverage = coverage)
}

wilcox_two_group <- function(df, ref, test, comparison) {
  out <- lapply(unique(df$signature), function(sig) {
    x <- df %>% filter(signature == sig, group %in% c(ref, test))
    a <- x$score[x$group == ref]
    b <- x$score[x$group == test]
    if (length(a) < 2 || length(b) < 2) return(NULL)
    p <- suppressWarnings(wilcox.test(b, a, exact = FALSE)$p.value)
    data.frame(
      dataset = unique(df$dataset),
      comparison = comparison,
      reference_group = ref,
      test_group = test,
      paired = FALSE,
      signature = sig,
      n_reference = length(a),
      n_test = length(b),
      delta_mean = mean(b, na.rm = TRUE) - mean(a, na.rm = TRUE),
      median_reference = median(a, na.rm = TRUE),
      median_test = median(b, na.rm = TRUE),
      p_value = p,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(out)
}

wilcox_combined_group <- function(df, ref, test_groups, test_label, comparison) {
  tmp <- df %>%
    filter(group == ref | group %in% test_groups) %>%
    mutate(group2 = if_else(group == ref, ref, test_label))
  out <- lapply(unique(tmp$signature), function(sig) {
    x <- tmp %>% filter(signature == sig)
    a <- x$score[x$group2 == ref]
    b <- x$score[x$group2 == test_label]
    if (length(a) < 2 || length(b) < 2) return(NULL)
    p <- suppressWarnings(wilcox.test(b, a, exact = FALSE)$p.value)
    data.frame(
      dataset = unique(df$dataset),
      comparison = comparison,
      reference_group = ref,
      test_group = test_label,
      paired = FALSE,
      signature = sig,
      n_reference = length(a),
      n_test = length(b),
      delta_mean = mean(b, na.rm = TRUE) - mean(a, na.rm = TRUE),
      median_reference = median(a, na.rm = TRUE),
      median_test = median(b, na.rm = TRUE),
      p_value = p,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(out)
}

wilcox_paired <- function(df, ref, test, comparison) {
  out <- lapply(unique(df$signature), function(sig) {
    wide <- df %>%
      filter(signature == sig, group %in% c(ref, test), !is.na(pair_id)) %>%
      select(pair_id, group, score) %>%
      distinct() %>%
      pivot_wider(names_from = group, values_from = score) %>%
      filter(!is.na(.data[[ref]]), !is.na(.data[[test]]))
    if (nrow(wide) < 3) return(NULL)
    p <- suppressWarnings(wilcox.test(wide[[test]], wide[[ref]], paired = TRUE, exact = FALSE)$p.value)
    data.frame(
      dataset = unique(df$dataset),
      comparison = comparison,
      reference_group = ref,
      test_group = test,
      paired = TRUE,
      signature = sig,
      n_reference = nrow(wide),
      n_test = nrow(wide),
      delta_mean = mean(wide[[test]] - wide[[ref]], na.rm = TRUE),
      median_reference = median(wide[[ref]], na.rm = TRUE),
      median_test = median(wide[[test]], na.rm = TRUE),
      p_value = p,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(out)
}

run_signature_statistics <- function(scores) {
  out <- list()
  for (ds in unique(scores$dataset)) {
    df <- scores %>% filter(dataset == ds)
    if (ds %in% c("GSE122897", "GSE15629")) {
      out[[paste0(ds, "_UIA_vs_Control")]] <- wilcox_two_group(df, "Control", "Unruptured", "UIA_vs_Control")
      out[[paste0(ds, "_RIA_vs_UIA")]] <- wilcox_two_group(df, "Unruptured", "Ruptured", "RIA_vs_UIA")
      out[[paste0(ds, "_RIA_vs_Control")]] <- wilcox_two_group(df, "Control", "Ruptured", "RIA_vs_Control")
      out[[paste0(ds, "_IA_all_vs_Control")]] <- wilcox_combined_group(df, "Control", c("Unruptured", "Ruptured"), "IA_all", "IA_all_vs_Control")
    } else if (ds == "GSE13353") {
      out[[paste0(ds, "_RIA_vs_UIA")]] <- wilcox_two_group(df, "Unruptured", "Ruptured", "RIA_vs_UIA")
    } else if (ds == "GSE75436") {
      out[[paste0(ds, "_IA_vs_STA")]] <- wilcox_paired(df, "STA", "IA", "IA_vs_paired_STA")
    }
  }
  stats <- bind_rows(out)
  stats %>%
    group_by(comparison) %>%
    mutate(p_adj_within_comparison = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(p_adj_global = p.adjust(p_value, method = "BH"))
}

run_trend_statistics <- function(scores) {
  trend_datasets <- scores %>%
    filter(dataset %in% c("GSE122897", "GSE15629"), group %in% c("Control", "Unruptured", "Ruptured")) %>%
    mutate(stage_code = case_when(
      group == "Control" ~ 0,
      group == "Unruptured" ~ 1,
      group == "Ruptured" ~ 2,
      TRUE ~ NA_real_
    ))
  trend_datasets %>%
    group_by(dataset, signature) %>%
    summarize(
      n = n(),
      rho = suppressWarnings(cor(stage_code, score, method = "spearman", use = "pairwise.complete.obs")),
      p_value = suppressWarnings(cor.test(stage_code, score, method = "spearman", exact = FALSE)$p.value),
      .groups = "drop"
    ) %>%
    group_by(dataset) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup()
}

run_correlation_statistics <- function(scores) {
  scores %>%
    select(dataset, sample_id, signature, score) %>%
    pivot_wider(names_from = signature, values_from = score) %>%
    group_by(dataset) %>%
    group_modify(function(.x, .y) {
      targets <- setdiff(names(signature_sets), "AMPK_SIRT1")
      bind_rows(lapply(targets, function(target) {
        ct <- suppressWarnings(cor.test(.x$AMPK_SIRT1, .x[[target]], method = "spearman", exact = FALSE))
        data.frame(
          target_signature = target,
          n = sum(complete.cases(.x$AMPK_SIRT1, .x[[target]])),
          rho = unname(ct$estimate),
          p_value = ct$p.value,
          stringsAsFactors = FALSE
        )
      }))
    }) %>%
    ungroup() %>%
    group_by(dataset) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup()
}

run_gene_statistics <- function(gene_long) {
  gene_long <- gene_long %>%
    rename(score = z_expression, signature = gene)
  run_signature_statistics(gene_long) %>%
    rename(gene = signature)
}

save_plot_all_formats <- function(plot, name, width, height) {
  stage_dir <- "/tmp/ia_external_validation_publication_figures"
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  targets <- file.path(figure_dir, paste0(name, c(".pdf", ".png", ".tiff")))
  staged <- file.path(stage_dir, paste0(name, c(".pdf", ".png", ".tiff")))
  ggsave(staged[1], plot, width = width, height = height, units = "in", device = pdf)
  ggsave(staged[2], plot, width = width, height = height, units = "in", dpi = 600)
  ggsave(staged[3], plot, width = width, height = height, units = "in", dpi = 600, compression = "lzw")
  ok <- file.copy(staged, targets, overwrite = TRUE)
  if (!all(ok)) {
    for (i in seq_along(staged)) {
      system2("cp", c(staged[i], targets[i]))
    }
  }
}

p_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

make_figures <- function(scores, stats, trends, cors) {
  group_levels <- c("Control", "STA", "Unruptured", "IA", "Ruptured")
  group_labels <- c(Control = "Control", STA = "STA", Unruptured = "UIA", IA = "IA", Ruptured = "RIA")
  group_colors <- c(Control = "#7A7A7A", STA = "#7A7A7A", Unruptured = "#2A9D8F", IA = "#277DA1", Ruptured = "#D1495B")
  sig_levels <- c("AMPK_SIRT1", "Contractile", "Inflammatory", "ECM_Remodeling")
  sig_labels <- c(
    AMPK_SIRT1 = "AMPK-SIRT1",
    Contractile = "Contractile VSMC",
    Inflammatory = "Inflammatory",
    ECM_Remodeling = "ECM remodeling"
  )

  plot_scores <- scores %>%
    mutate(
      group = factor(group, levels = group_levels),
      signature = factor(signature, levels = sig_levels, labels = sig_labels[sig_levels]),
      dataset = factor(dataset, levels = c("GSE122897", "GSE13353", "GSE15629", "GSE75436"))
    )

  p_scores <- ggplot(plot_scores, aes(group, score, color = group, fill = group)) +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.25) +
    geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.16, linewidth = 0.38) +
    geom_jitter(width = 0.14, height = 0, size = 1.35, alpha = 0.78, stroke = 0) +
    facet_grid(signature ~ dataset, scales = "free_x", space = "free_x") +
    scale_x_discrete(labels = group_labels, drop = TRUE) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(x = NULL, y = "Module score (mean gene-level z-score)") +
    theme_bw(base_size = 9, base_family = "Helvetica") +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey95", color = "grey72", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = 8.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 7.5),
      axis.text.y = element_text(size = 7.5),
      axis.title.y = element_text(size = 9)
    )
  save_plot_all_formats(p_scores, "Figure_External_Validation_A_signature_scores", 12.2, 8.2)

  comparison_order <- c(
    "GSE122897 | IA_all_vs_Control",
    "GSE122897 | UIA_vs_Control",
    "GSE122897 | RIA_vs_UIA",
    "GSE122897 | RIA_vs_Control",
    "GSE13353 | RIA_vs_UIA",
    "GSE15629 | IA_all_vs_Control",
    "GSE15629 | UIA_vs_Control",
    "GSE15629 | RIA_vs_UIA",
    "GSE15629 | RIA_vs_Control",
    "GSE75436 | IA_vs_paired_STA"
  )
  heat <- stats %>%
    mutate(
      row_label = paste(dataset, comparison, sep = " | "),
      row_label = factor(row_label, levels = rev(comparison_order)),
      signature = factor(signature, levels = sig_levels, labels = sig_labels[sig_levels]),
      label = paste0(sprintf("%.2f", delta_mean), p_stars(p_adj_global))
    )

  p_heat <- ggplot(heat, aes(signature, row_label, fill = delta_mean)) +
    geom_tile(color = "white", linewidth = 0.55) +
    geom_text(aes(label = label), size = 2.65, family = "Helvetica") +
    scale_fill_gradient2(
      low = "#3B6EA8", mid = "white", high = "#C44536", midpoint = 0,
      name = "Mean\ndifference"
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 9, base_family = "Helvetica") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", size = 8),
      axis.text.y = element_text(size = 7.7),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7.5)
    )
  save_plot_all_formats(p_heat, "Figure_External_Validation_B_effect_heatmap", 8.6, 4.8)

  cor_plot <- cors %>%
    mutate(
      target_signature = factor(target_signature, levels = c("Contractile", "Inflammatory", "ECM_Remodeling"), labels = sig_labels[c("Contractile", "Inflammatory", "ECM_Remodeling")]),
      dataset = factor(dataset, levels = c("GSE122897", "GSE13353", "GSE15629", "GSE75436")),
      label = paste0(sprintf("%.2f", rho), p_stars(p_adj))
    )
  p_cor <- ggplot(cor_plot, aes(dataset, target_signature, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.55) +
    geom_text(aes(label = label), size = 3.0, family = "Helvetica") +
    scale_fill_gradient2(low = "#3B6EA8", mid = "white", high = "#C44536", midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho") +
    labs(x = NULL, y = "Correlation with AMPK-SIRT1 score") +
    theme_minimal(base_size = 9, base_family = "Helvetica") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 9)
    )
  save_plot_all_formats(p_cor, "Figure_External_Validation_C_AMPK_correlations", 6.4, 2.8)

  p_composite <- (p_scores / (p_heat | p_cor)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 12, family = "Helvetica"))
  save_plot_all_formats(p_composite, "Figure_External_Validation_Composite", 13.2, 12.0)
}

make_report <- function(sample_table, coverage, stats, trends, cors) {
  key_stats <- stats %>%
    filter(signature %in% c("AMPK_SIRT1", "Inflammatory", "ECM_Remodeling", "Contractile")) %>%
    mutate(direction = if_else(delta_mean > 0, "higher", "lower")) %>%
    arrange(dataset, comparison, signature)

  sig_hits <- key_stats %>%
    filter(p_adj_global < 0.05) %>%
    mutate(text = paste0(dataset, " ", comparison, ": ", signature, " ", direction, " (delta=", sprintf("%.2f", delta_mean), ", FDR=", signif(p_adj_global, 3), ")")) %>%
    pull(text)
  if (!length(sig_hits)) sig_hits <- "No signature-level comparisons passed global BH-FDR < 0.05; inspect effect sizes and within-comparison FDR."

  lines <- c(
    "# External cohort validation summary",
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## Included cohorts",
    paste(capture.output(print(sample_table, row.names = FALSE)), collapse = "\n"),
    "",
    "## Main significant signature-level findings (global BH-FDR < 0.05)",
    paste0("- ", sig_hits),
    "",
    "## Interpretation guardrails",
    "- These datasets are bulk aneurysm-wall or arterial-wall transcriptomic cohorts, so they validate tissue-level reproducibility rather than VSMC-specific single-cell states.",
    "- Directional AMPK-SIRT1 changes should be interpreted together with inflammatory and ECM-remodeling scores, because stress compensation can raise individual metabolic-defense genes in inflamed tissue.",
    "- Ruptured versus unruptured comparisons are not uniformly powered across cohorts; effect-size consistency is more informative than any single nominal p-value.",
    "",
    "## Output files",
    "- processed/signature_scores_all_external_cohorts.csv",
    "- processed/signature_statistics_external_cohorts.csv",
    "- processed/signature_trend_statistics_external_cohorts.csv",
    "- processed/AMPK_signature_correlation_external_cohorts.csv",
    "- processed/key_gene_statistics_external_cohorts.csv",
    "- figures/Figure_External_Validation_Composite.pdf/png/tiff",
    "- figures/Figure_External_Validation_A_signature_scores.pdf/png/tiff",
    "- figures/Figure_External_Validation_B_effect_heatmap.pdf/png/tiff",
    "- figures/Figure_External_Validation_C_AMPK_correlations.pdf/png/tiff"
  )
  writeLines(lines, file.path(project_dir, "External_Cohort_Validation_Summary.md"))
}

message_time("Starting external validation")

gse122897 <- read_gse122897_counts()
gse13353 <- read_microarray_dataset(
  "GSE13353",
  "GSE13353_series_matrix.txt.gz",
  "GPL570.annot.gz",
  "GPL570",
  function(meta) {
    case_when(
      str_detect(str_to_lower(meta$source_name), "unruptured") ~ "Unruptured",
      str_detect(str_to_lower(meta$source_name), "ruptured") ~ "Ruptured",
      TRUE ~ "Unknown"
    )
  }
)
gse15629 <- read_microarray_dataset(
  "GSE15629",
  "GSE15629_series_matrix.txt.gz",
  "GPL6244.annot.gz",
  "GPL6244",
  function(meta) {
    case_when(
      str_detect(str_to_lower(meta$type), "control") ~ "Control",
      str_detect(str_to_lower(meta$type), "unruptured") ~ "Unruptured",
      str_detect(str_to_lower(meta$type), "ruptured") ~ "Ruptured",
      TRUE ~ "Unknown"
    )
  }
)
gse75436 <- read_microarray_dataset(
  "GSE75436",
  "GSE75436_series_matrix.txt.gz",
  "GPL570.annot.gz",
  "GPL570",
  function(meta) {
    case_when(
      str_detect(str_to_lower(meta$tissue), "superficial temporal artery") ~ "STA",
      str_detect(str_to_lower(meta$tissue), "intracranial aneurysm") ~ "IA",
      TRUE ~ "Unknown"
    )
  }
)
gse75436$meta$pair_id <- gse75436$meta$individual

datasets <- list(GSE122897 = gse122897, GSE13353 = gse13353, GSE15629 = gse15629, GSE75436 = gse75436)

sample_table <- bind_rows(lapply(datasets, function(x) x$meta)) %>%
  count(dataset, group, name = "n_samples") %>%
  arrange(dataset, group)
fwrite(sample_table, file.path(processed_dir, "sample_counts_external_cohorts.csv"))

score_objects <- lapply(names(datasets), function(ds) {
  message_time("Scoring signatures for ", ds)
  score_signatures(datasets[[ds]]$expr, datasets[[ds]]$meta, ds)
})
names(score_objects) <- names(datasets)

scores <- bind_rows(lapply(score_objects, `[[`, "scores"))
gene_long <- bind_rows(lapply(score_objects, `[[`, "genes"))
coverage <- bind_rows(lapply(score_objects, `[[`, "coverage"))

stats <- run_signature_statistics(scores)
trends <- run_trend_statistics(scores)
cors <- run_correlation_statistics(scores)
gene_stats <- run_gene_statistics(gene_long)

fwrite(scores, file.path(processed_dir, "signature_scores_all_external_cohorts.csv"))
fwrite(gene_long, file.path(processed_dir, "key_gene_zscores_all_external_cohorts.csv"))
fwrite(coverage, file.path(processed_dir, "signature_gene_coverage_external_cohorts.csv"))
fwrite(stats, file.path(processed_dir, "signature_statistics_external_cohorts.csv"))
fwrite(trends, file.path(processed_dir, "signature_trend_statistics_external_cohorts.csv"))
fwrite(cors, file.path(processed_dir, "AMPK_signature_correlation_external_cohorts.csv"))
fwrite(gene_stats, file.path(processed_dir, "key_gene_statistics_external_cohorts.csv"))

message_time("Generating publication figures")
make_figures(scores, stats, trends, cors)
make_report(sample_table, coverage, stats, trends, cors)

message_time("External validation complete")
