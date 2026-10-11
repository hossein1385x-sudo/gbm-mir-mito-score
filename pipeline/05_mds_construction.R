#!/usr/bin/env Rscript
# =============================================================================
# 05_mds_construction_train.R - PCA-based version
# -----------------------------------------------------------------------------
# MDS = PC1 of the z-scored matrix of mitochondrial genes
# The sign of PC1 is set automatically so that it points in the same direction as mean mito expression
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(matrixStats)
})

log_msg <- function(...) cat(sprintf("[%s] %s\n",
                                     format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 ||
                             (length(a) == 1 && is.na(a))) b else a

normalize_tcga_id <- function(ids) {
  ids <- toupper(as.character(ids))
  ids <- gsub("\\.", "-", ids)
  ids <- trimws(ids)
  m <- regmatches(ids, regexpr("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", ids))
  out <- rep(NA_character_, length(ids))
  matched <- nchar(m) > 0
  out[matched] <- m[matched]
  out
}

run_mds_construction <- function(input_root,
                                 sample_col  = NULL,
                                 gene_col    = "ensembl_id",
                                 force       = FALSE,
                                 min_samples = 10,
                                 method      = c("pca", "ssgsea", "mean")) {
  method <- match.arg(method)

  log_msg("==================== Step 05: MDS Construction (Train Only) ====================")
  log_msg(sprintf("MDS construction method: %s", method))

  ## ---- 1) train ----
  train_path <- file.path(input_root, "train_set.csv")
  if (!file.exists(train_path))
    train_path <- file.path(input_root, "train_set_ORIGINAL_BACKUP.csv")
  if (!file.exists(train_path)) stop("train_set.csv not found.")

  train_meta <- fread(train_path)
  log_msg(sprintf("reading %s: %d records", basename(train_path), nrow(train_meta)))

  if (is.null(sample_col)) {
    candidates <- c("patient_barcode", "barcode", "sample_id", "patient", "case_id")
    sample_col <- intersect(candidates, names(train_meta))[1]
    if (is.na(sample_col)) stop("sample column not found.")
  }

  train_keys <- unique(normalize_tcga_id(train_meta[[sample_col]]))
  train_keys <- train_keys[!is.na(train_keys)]
  log_msg(sprintf("  unique train patients: %d", length(train_keys)))

  ## ---- 2) RNA-seq ----
  rna_path <- file.path(input_root, "raw_rnaseq_counts_matrix.rds")
  if (!file.exists(rna_path))
    rna_path <- file.path(input_root, "raw_rnaseq.rds")
  if (!file.exists(rna_path)) stop("RNA-seq file not found.")

  raw_counts <- readRDS(rna_path)
  if (!is.matrix(raw_counts)) raw_counts <- as.matrix(raw_counts)
  storage.mode(raw_counts) <- "numeric"
  log_msg(sprintf("  RNA-seq matrix: %d genes × %d samples",
                  nrow(raw_counts), ncol(raw_counts)))

  col_keys <- normalize_tcga_id(colnames(raw_counts))
  common_keys <- intersect(train_keys, col_keys[!is.na(col_keys)])
  log_msg(sprintf("  shared patients: %d of %d", length(common_keys), length(train_keys)))

  if (length(common_keys) < min_samples) stop("too few shared patients.")

  col_idx <- match(common_keys, col_keys)
  counts_train_all <- raw_counts[, col_idx, drop = FALSE]
  colnames(counts_train_all) <- common_keys

  ## ---- 3) Mito filter ----
  mito_genes <- fread(file.path(input_root, "mito_gene_list.csv"))
  if (!gene_col %in% names(mito_genes))
    stop(sprintf("column '%s' is not in mito_gene_list.csv.", gene_col))

  mito_ids <- gsub("\\.[0-9]+$", "", unique(na.omit(as.character(mito_genes[[gene_col]]))))
  rn_norm  <- gsub("\\.[0-9]+$", "", rownames(counts_train_all))
  mito_idx <- which(rn_norm %in% mito_ids)

  if (length(mito_idx) == 0) stop("no mitochondrial gene found.")

  counts_train_mito <- counts_train_all[mito_idx, , drop = FALSE]
  log_msg(sprintf("  %d mitochondrial genes will be used.", nrow(counts_train_mito)))

  ## ---- 4) DESeq2 normalisation ----
  log_msg("normalising with DESeq2 ...")
  counts_int <- round(counts_train_mito)
  counts_int[counts_int < 0] <- 0

  zero_genes <- rowSums(counts_int) == 0
  if (sum(zero_genes) > 0) {
    log_msg(sprintf("  removed %d genes with zero total.", sum(zero_genes)))
    counts_int <- counts_int[!zero_genes, , drop = FALSE]
  }

  col_data <- data.frame(row.names = colnames(counts_int),
                         dummy = factor(rep("A", ncol(counts_int))))
  dds <- DESeqDataSetFromMatrix(counts_int, col_data, ~ 1)
  dds <- estimateSizeFactors(dds)
  size_factors <- sizeFactors(dds)
  norm_counts  <- counts(dds, normalized = TRUE)
  log_norm     <- log2(norm_counts + 1)

  ## ---- 5) z-score ----
  log_msg("computing z-scores ...")
  gene_means <- rowMeans(log_norm)
  gene_sds   <- rowSds(log_norm)
  gene_sds[gene_sds == 0 | is.na(gene_sds)] <- 1e-8
  z <- (log_norm - gene_means) / gene_sds

  ## ---- 6) Build the MDS with the selected method ----
  if (method == "pca") {
    log_msg("building the MDS with PCA (PC1) ...")
    # PCA on samples (rows = samples, columns = genes)
    pca <- prcomp(t(z), center = FALSE, scale. = FALSE)
    mds_raw <- pca$x[, 1]

    # orient PC1 so that its direction agrees with the mean mito expression
    ref <- colMeans(log_norm)
    if (cor(mds_raw, ref) < 0) {
      log_msg("  PC1 sign was flipped to agree with the mean mito expression.")
      mds_raw <- -mds_raw
    }

    # explained variance
    var_explained <- (pca$sdev^2 / sum(pca$sdev^2))[1]
    log_msg(sprintf("  variance explained by PC1: %.1f%%", 100 * var_explained))

    mds_scores <- data.table(sample_id = colnames(z), MDS = as.numeric(mds_raw))
    params_extra <- list(
      method            = "pca",
      pca_rotation      = pca$rotation,
      pca_sdev          = pca$sdev,
      pca_center        = pca$center,
      pca_scale         = pca$scale,
      sign_corrected    = TRUE,
      var_explained_pc1 = var_explained
    )

  } else if (method == "ssgsea") {
    log_msg("building the MDS with ssGSEA ...")
    if (!requireNamespace("GSVA", quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("GSVA", update = FALSE, ask = FALSE)
    }
    library(GSVA)
    gene_set <- list(MitoDysfunction = rownames(z))
    ss <- gsva(as.matrix(z), gene_set, method = "ssgsea",
               verbose = FALSE, kcdf = "Gaussian")
    mds_raw <- as.numeric(ss[1, ])

    mds_scores <- data.table(sample_id = colnames(z), MDS = mds_raw)
    params_extra <- list(method = "ssgsea")

  } else {
    log_msg("building the MDS with the mean z-score ...")
    mds_raw <- colMeans(z)
    mds_scores <- data.table(sample_id = names(mds_raw), MDS = as.numeric(mds_raw))
    params_extra <- list(method = "mean")
  }

  log_msg(sprintf("MDS computed for %d patients.", nrow(mds_scores)))
  log_msg(sprintf("   mean: %.4f | SD: %.4f | range: [%.4f, %.4f]",
                  mean(mds_scores$MDS), sd(mds_scores$MDS),
                  min(mds_scores$MDS), max(mds_scores$MDS)))

  ## ---- 7) Save ----
  out_scores <- file.path(input_root, "mds_scores_train.csv")
  out_params <- file.path(input_root, "mds_transform_params.rds")

  if (!force && file.exists(out_scores) && file.exists(out_params)) {
    log_msg("files already exist. Pass force=TRUE to overwrite.")
    return(invisible(NULL))
  }

  fwrite(mds_scores, out_scores)

  params <- c(
    list(
      mito_genes_used    = rownames(counts_int),
      gene_means_train   = gene_means,
      gene_sds_train     = gene_sds,
      size_factors_train = size_factors,
      train_samples      = colnames(counts_int),
      gene_col_source    = gene_col,
      n_genes            = nrow(counts_int),
      n_samples          = ncol(counts_int),
      fit_timestamp      = as.character(Sys.time()),
      r_version          = R.version.string
    ),
    params_extra
  )
  saveRDS(params, out_params)

  log_msg(sprintf("saved: %s", out_scores))
  log_msg(sprintf("saved: %s", out_params))
  log_msg("Step 05 complete.")

  invisible(list(scores = mds_scores, params = params))
}

## ---- CLI ----
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  parse_args <- function(argv) {
    opts <- list(); i <- 1
    while (i <= length(argv)) {
      if (substr(argv[i], 1, 2) == "--") {
        key <- sub("^--", "", argv[i])
        if (i < length(argv) && substr(argv[i + 1], 1, 2) != "--") {
          opts[[key]] <- argv[i + 1]; i <- i + 2
        } else { opts[[key]] <- TRUE; i <- i + 1 }
      } else i <- i + 1
    }
    opts
  }
  opts <- parse_args(args)
  default_root <- WORK_DIR

  run_mds_construction(
    input_root  = opts[["input-root"]] %||% default_root,
    sample_col  = opts[["sample-col"]] %||% NULL,
    gene_col    = opts[["gene-col"]]   %||% "ensembl_id",
    force = TRUE,
    min_samples = as.integer(opts[["min-samples"]] %||% 10),
    method      = opts[["method"]]     %||% "pca"
  )
}
