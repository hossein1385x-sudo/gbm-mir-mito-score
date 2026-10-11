#!/usr/bin/env Rscript
# =============================================================================
# 06_DE_analysis_train_only.R — FIXED v2
# -----------------------------------------------------------------------------
# DE analysis between MDS-high and MDS-low patients (train set only)
#
# Inputs: mirna_isoform_counts_matrix.csv  (preferred)
#        mirna_isoform_rpm_matrix.csv     (fallback, limma-voom)
#        mds_scores_train.csv             (output of step 5)
#
# Outputs: DE_results_train.csv
#        candidate_miRNAs_isoform.csv
#        de_manifest.json
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(jsonlite)
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
  ok <- nchar(m) > 0
  out[ok] <- m[ok]
  out
}

looks_like_mirna <- function(x) {
  mean(grepl("^hsa-", x, ignore.case = TRUE)) > 0.5
}

## ---- Read the miRNA matrix + detect counts/RPM -----------------------------
read_mirna_matrix <- function(input_root) {
  counts_candidates <- c(
    "mirna_isoform_counts_matrix.csv",
    "mirna_isoform_counts_matrix.rds"
  )
  rpm_candidates <- c(
    "mirna_isoform_rpm_matrix.csv",
    "mirna_isoform_rpm_matrix.rds"
  )

  found <- NULL
  is_rpm <- FALSE
  for (f in counts_candidates) {
    p <- file.path(input_root, f)
    if (file.exists(p)) { found <- p; is_rpm <- FALSE; break }
  }
  if (is.null(found)) {
    for (f in rpm_candidates) {
      p <- file.path(input_root, f)
      if (file.exists(p)) { found <- p; is_rpm <- TRUE; break }
    }
  }
  if (is.null(found)) stop("no miRNA matrix found.")

  log_msg(sprintf("reading %s (type: %s) ...", basename(found),
                  if (is_rpm) "RPM" else "counts"))

  if (grepl("\\.rds$", found, ignore.case = TRUE)) {
    raw <- readRDS(found)
    if (!is.matrix(raw)) raw <- as.matrix(raw)
  } else {
    raw <- as.matrix(read.csv(found, row.names = 1, check.names = FALSE))
  }
  storage.mode(raw) <- "numeric"

  ## ---- detect the matrix orientation ----
  rn_mir <- looks_like_mirna(rownames(raw))
  cn_mir <- looks_like_mirna(colnames(raw))

  if (rn_mir && !cn_mir) {
    mir_mat <- raw
  } else if (cn_mir && !rn_mir) {
    mir_mat <- t(raw)
  } else if (rn_mir && cn_mir) {
    stop("both rows and columns are miRNAs — invalid matrix.")
  } else {
    stop("neither rows nor columns are miRNAs.")
  }

  list(mat = mir_mat, is_rpm = is_rpm, file = found)
}

## ---- Statistical tests: DESeq2 on counts, limma on RPM ---------------------
run_deseq2 <- function(counts_mat, groups, padj_thr, log2fc_thr) {
  log_msg("running DESeq2 (on counts) ...")
  col_data <- data.frame(row.names = colnames(counts_mat),
                         mds_group = groups)
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_mat),
    colData   = col_data,
    design    = ~ mds_group
  )
  # filter: at least 3 samples with counts ≥ 10
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  log_msg(sprintf("  after low-count filter: %d isoforms", nrow(dds)))

  if (nrow(dds) < 50) stop("the number of isoforms after filtering is very low.")

  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds,
                 contrast = c("mds_group", "high", "low"),
                 alpha = padj_thr)
  res_df <- as.data.frame(res)
  res_df$isoform_id <- rownames(res_df)
  res_df
}

run_limma_voom <- function(counts_mat, groups, padj_thr) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("limma", update = FALSE, ask = FALSE)
  }
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    BiocManager::install("edgeR", update = FALSE, ask = FALSE)
  }
  library(limma); library(edgeR)

  log_msg("running limma-voom (on RPM) ...")
  # filter
  keep <- rowSums(counts_mat > 1) >= 3
  counts_mat <- counts_mat[keep, , drop = FALSE]
  log_msg(sprintf("  after filter: %d isoforms", nrow(counts_mat)))

  dge <- DGEList(counts = counts_mat)
  # RPM values are already library-size normalised, so TMM normalisation is skipped
  dge <- calcNormFactors(dge, method = "none")

  design <- model.matrix(~ groups)
  colnames(design) <- c("Intercept", "mds_grouphigh")

  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  res_df <- topTable(fit, coef = "mds_grouphigh", number = Inf,
                     adjust.method = "BH", sort.by = "P")
  res_df$isoform_id <- rownames(res_df)
  res_df
}

## =============================================================================
run_de_analysis <- function(input_root,
                            padj_threshold   = 0.05,
                            log2fc_threshold = 1.0,
                            split_method     = c("median", "tertile"),
                            force            = TRUE) {

  split_method <- match.arg(split_method)

  mds_file <- file.path(input_root, "mds_scores_train.csv")
  if (!file.exists(mds_file)) stop("mds_scores_train.csv not found.")

  log_msg("==================== Step 06: DE Analysis (Train Only) ====================")

  ## ---- 1) Read the miRNA matrix ----
  mir_info <- read_mirna_matrix(input_root)
  mir_mat  <- mir_info$mat
  is_rpm   <- mir_info$is_rpm

  log_msg(sprintf("  dimensions: %d isoforms × %d samples", nrow(mir_mat), ncol(mir_mat)))

  ## ---- 2) Collapse barcodes to patient level ----
  sample_ids_norm <- normalize_tcga_id(colnames(mir_mat))
  if (anyNA(sample_ids_norm)) {
    stop("some barcodes could not be normalised.")
  }
  if (anyDuplicated(sample_ids_norm)) {
    log_msg("  collapse aliquot → patient ...")
    mir_dt <- as.data.table(t(mir_mat))
    mir_dt[, patient_id := sample_ids_norm]
    # for counts: sum. for RPM: mean.
    agg_fun <- if (is_rpm) mean else sum
    mir_dt <- mir_dt[, lapply(.SD, function(x) agg_fun(x, na.rm = TRUE)),
                     by = patient_id,
                     .SDcols = setdiff(names(mir_dt), "patient_id")]
    mir_mat <- t(as.matrix(mir_dt[, -1, with = FALSE]))
    colnames(mir_mat) <- mir_dt$patient_id
    log_msg(sprintf("  after collapse: %d samples", ncol(mir_mat)))
  } else {
    colnames(mir_mat) <- sample_ids_norm
  }

  ## ---- 3) MDS ----
  log_msg(sprintf("reading %s ...", basename(mds_file)))
  mds_scores <- fread(mds_file)
  if (!all(c("sample_id", "MDS") %in% names(mds_scores))) {
    stop("columns sample_id and MDS are required in mds_scores_train.csv.")
  }
  mds_scores[, sample_id := normalize_tcga_id(sample_id)]

  matched <- intersect(mds_scores$sample_id, colnames(mir_mat))
  log_msg(sprintf("  overlap with MDS: %d samples of %d",
                  length(matched), nrow(mds_scores)))

  if (length(matched) < 20) {
    stop("fewer than 20 shared samples — analysis is not possible.")
  }

  missing <- setdiff(mds_scores$sample_id, matched)
  if (length(missing) > 0) {
    log_msg(sprintf("  %d MDS samples are not in the miRNA matrix.", length(missing)))
  }

  mir_train <- mir_mat[, matched, drop = FALSE]
  mds_scores <- mds_scores[match(matched, mds_scores$sample_id), ]
  stopifnot(all(mds_scores$sample_id == colnames(mir_train)))

  ## ---- 4) Grouping ----
  if (split_method == "median") {
    cutoff <- median(mds_scores$MDS)
    mds_group <- factor(ifelse(mds_scores$MDS >= cutoff, "high", "low"),
                        levels = c("low", "high"))
  } else {
    q <- quantile(mds_scores$MDS, probs = c(1/3, 2/3))
    mds_group <- cut(mds_scores$MDS,
                     breaks = c(-Inf, q[1], q[2], Inf),
                     labels = c("low", "mid", "high"))
    keep <- mds_group != "mid"
    mir_train <- mir_train[, keep, drop = FALSE]
    mds_group <- droplevels(mds_group[keep])
    mds_scores <- mds_scores[keep, ]
  }

  log_msg(sprintf("  groups: low=%d, high=%d",
                  sum(mds_group == "low"), sum(mds_group == "high")))

  if (min(table(mds_group)) < 10) {
    stop("fewer than 10 samples per group — analysis is not meaningful.")
  }

  ## ---- 5) Test ----
  if (!is_rpm) {
    res_df <- run_deseq2(mir_train, mds_group, padj_threshold, log2fc_threshold)
    method_used <- "DESeq2 (counts)"
    padj_col    <- "padj"
    lfc_col     <- "log2FoldChange"
  } else {
    res_df <- run_limma_voom(mir_train, mds_group, padj_threshold)
    method_used <- "limma-voom (RPM)"
    padj_col    <- "adj.P.Val"
    lfc_col     <- "logFC"
  }

  ## ---- 6) Filter candidates ----
  res_df <- res_df[order(res_df[[padj_col]]), ]

  keep_idx <- !is.na(res_df[[padj_col]]) &
    res_df[[padj_col]] < padj_threshold &
    abs(res_df[[lfc_col]]) >= log2fc_threshold
  candidates <- res_df[keep_idx, , drop = FALSE]

  log_msg(sprintf("  %d significant isoforms (padj<%.2f, |log2FC|>=%.1f)",
                  nrow(candidates), padj_threshold, log2fc_threshold))

  # report the families that have both arms
  if (nrow(candidates) > 0) {
    base <- gsub("-(5p|3p)$", "", candidates$isoform_id)
    dup  <- base[duplicated(base) | duplicated(base, fromLast = TRUE)]
    if (length(unique(dup)) > 0) {
      log_msg(sprintf("  %d families with both arms: %s",
                      length(unique(dup)),
                      paste(head(unique(dup), 5), collapse = ", ")))
    }
  }

  ## ---- 7) Save ----
  out_de    <- file.path(input_root, "DE_results_train.csv")
  out_cand  <- file.path(input_root, "candidate_miRNAs_isoform.csv")
  out_man   <- file.path(input_root, "de_manifest.json")

  if (!force && file.exists(out_de)) {
    log_msg("files already exist. Use force=TRUE to overwrite.")
    return(invisible(NULL))
  }

  for (f in c(out_de, out_cand)) {
    if (file.exists(f)) {
      tryCatch(Sys.chmod(f, mode = "0666"), error = function(e) NULL)
      file.remove(f)
    }
  }

  fwrite(as.data.table(res_df), out_de)
  fwrite(as.data.table(candidates), out_cand)

  manifest <- list(
    input_mirna_file  = mir_info$file,
    data_type         = if (is_rpm) "RPM" else "counts",
    method_used       = method_used,
    split_method      = split_method,
    padj_threshold    = padj_threshold,
    log2fc_threshold  = log2fc_threshold,
    n_train_samples   = length(matched),
    n_low             = sum(mds_group == "low"),
    n_high            = sum(mds_group == "high"),
    n_isoforms_tested = nrow(res_df),
    n_candidates      = nrow(candidates),
    timestamp         = as.character(Sys.time()),
    r_version         = R.version.string
  )
  write_json(manifest, out_man, pretty = TRUE, auto_unbox = TRUE)

  log_msg(sprintf("saved: %s", out_de))
  log_msg(sprintf("saved: %s", out_cand))
  log_msg(sprintf("saved: %s", out_man))
  log_msg("Step 06 complete.")

  invisible(list(de = res_df, candidates = candidates, manifest = manifest))
}

## ---- Entry point ----
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

  run_de_analysis(
    input_root       = opts[["input-root"]] %||% default_root,
    padj_threshold   = as.numeric(opts[["padj"]]   %||% 0.05),
    log2fc_threshold = as.numeric(opts[["log2fc"]] %||% 1.0),
    split_method     = opts[["split"]] %||% "median",
    force = TRUE
  )
}
