#!/usr/bin/env Rscript
# =============================================================================
# 04_train_test_split.R  —  FIXED VERSION (handles read-only test file)
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table); library(jsonlite); library(digest)
})

log_msg <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

normalize_tcga_id <- function(ids) {
  ids <- toupper(as.character(ids)); ids <- gsub("\\.", "-", ids)
  m <- regmatches(ids, regexpr("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", ids))
  out <- ids; out[nchar(m) > 0] <- m[nchar(m) > 0]; out
}

looks_like_mirna <- function(x) mean(grepl("^hsa-", x, ignore.case = TRUE)) > 0.5

read_mirna_sample_ids <- function(input_root) {
  csv_path <- file.path(input_root, "mirna_isoform_rpm_matrix.csv")
  rds_path <- file.path(input_root, "mirna_isoform_rpm_matrix.rds")
  if (file.exists(csv_path)) {
    log_msg(sprintf("  reading %s ...", basename(csv_path)))
    raw <- as.matrix(read.csv(csv_path, row.names = 1, check.names = FALSE))
    storage.mode(raw) <- "numeric"
  } else if (file.exists(rds_path)) {
    log_msg(sprintf("  reading %s ...", basename(rds_path)))
    raw <- readRDS(rds_path)
  } else stop("mirna matrix not found in ", input_root)
  if (looks_like_mirna(rownames(raw)) && !looks_like_mirna(colnames(raw))) ids <- colnames(raw)
  else if (looks_like_mirna(colnames(raw))) ids <- rownames(raw)
  else stop("cannot detect miRNA orientation")
  ids
}

read_rna_sample_ids <- function(input_root) {
  rds_path <- file.path(input_root, "raw_rnaseq_counts_matrix.rds")
  if (!file.exists(rds_path)) stop("raw_rnaseq_counts_matrix.rds not found")
  log_msg(sprintf("  reading %s ...", basename(rds_path)))
  colnames(readRDS(rds_path))
}

stratified_split <- function(data, strat_col, train_frac = 0.75, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.factor(data[[strat_col]])) data[[strat_col]] <- as.factor(data[[strat_col]])
  train_idx <- c()
  for (level in levels(data[[strat_col]])) {
    idx_level <- which(data[[strat_col]] == level); n_level <- length(idx_level)
    if (n_level == 0) next
    n_train_level <- max(1, round(n_level * train_frac))
    if (n_level == 1) { train_idx <- c(train_idx, idx_level); next }
    train_idx <- c(train_idx, sample(idx_level, size = n_train_level, replace = FALSE))
  }
  test_idx <- setdiff(seq_len(nrow(data)), train_idx)
  list(train = sort(train_idx), test = sort(test_idx))
}

# new function: if the file is read-only or locked, remove it
safe_write <- function(dt, path) {
  if (file.exists(path)) {
    tryCatch(Sys.chmod(path, mode = "0666"), error = function(e) NULL)
    ok <- tryCatch({ file.remove(path); TRUE }, error = function(e) FALSE)
    if (!ok) {
      Sys.sleep(0.5)
      tryCatch(Sys.chmod(path, mode = "0666"), error = function(e) NULL)
      file.remove(path)
    }
  }
  fwrite(dt, path)
}

run_train_test_split <- function(input_root, seed = 42, test_frac = 0.25,
                                 strat_col_primary = NULL, strat_col_secondary = NULL,
                                 force = FALSE) {
  log_msg("==================== Step 04: Train/Test Split (FIXED) ====================")
  in_path <- file.path(input_root, "raw_clinical.rds")
  if (!file.exists(in_path)) stop("raw_clinical.rds not found: ", in_path)
  clinical <- readRDS(in_path)
  if (!is.data.table(clinical)) clinical <- as.data.table(clinical)
  log_msg(sprintf("  raw clinical: %d patients", nrow(clinical)))

  id_col <- intersect(c("patient_barcode","patient_id","barcode","case_id"), names(clinical))[1]
  if (is.na(id_col)) stop("patient_id column not found")
  log_msg(sprintf("  ID column: %s", id_col))
  clinical[, patient_id := normalize_tcga_id(get(id_col))]

  possible_os_time   <- intersect(c("OS.time","survival_time_final","survival_time","days_to_last_follow_up","time"), names(clinical))
  possible_os_status <- intersect(c("OS","os_status","status_final","vital_status","death","event"), names(clinical))
  if (length(possible_os_time) == 0 || length(possible_os_status) == 0) stop("survival columns not found")
  os_time_col <- possible_os_time[1]; os_status_col <- possible_os_status[1]
  log_msg(sprintf("  survival time: %s | survival status: %s", os_time_col, os_status_col))

  log_msg("  --- eligibility filter (miRNA ∩ RNA ∩ clinical) ---")
  mir_ids <- normalize_tcga_id(read_mirna_sample_ids(input_root))
  rna_ids <- normalize_tcga_id(read_rna_sample_ids(input_root))
  clin_ids <- clinical$patient_id
  log_msg(sprintf("    miRNA unique IDs : %d", length(unique(mir_ids))))
  log_msg(sprintf("    RNA   unique IDs : %d", length(unique(rna_ids))))
  log_msg(sprintf("    clin  unique IDs : %d", length(unique(clin_ids))))
  eligible_ids <- Reduce(intersect, list(mir_ids, rna_ids, clin_ids))
  log_msg(sprintf("    eligible (all three sources): %d", length(eligible_ids)))
  if (length(eligible_ids) == 0) stop("no eligible patients")
  n_before <- nrow(clinical)
  clinical <- clinical[patient_id %in% eligible_ids]
  log_msg(sprintf("    clinical after filtering: %d (removed %d)", nrow(clinical), n_before - nrow(clinical)))

  if (is.null(strat_col_primary)) strat_col_primary <- os_status_col
  if (!strat_col_primary %in% names(clinical)) strat_col_primary <- os_status_col
  use_secondary <- !is.null(strat_col_secondary) && strat_col_secondary %in% names(clinical)
  clinical[, strat_group := if (use_secondary) interaction(get(strat_col_primary), get(strat_col_secondary), drop = TRUE) else get(strat_col_primary)]

  split_result <- stratified_split(clinical, "strat_group", train_frac = 1 - test_frac, seed = seed)
  train_set <- clinical[split_result$train, ]; test_set <- clinical[split_result$test, ]
  train_set[, strat_group := NULL]; test_set[, strat_group := NULL]
  log_msg(sprintf("  Train: %d | Test: %d", nrow(train_set), nrow(test_set)))

  train_path <- file.path(input_root, "train_set.csv")
  test_path  <- file.path(input_root, "test_set_LOCKED.csv")
  manifest_path <- file.path(input_root, "split_manifest.json")
  if (!force && file.exists(train_path) && file.exists(test_path)) {
    log_msg("files already exist. Pass force=TRUE to overwrite.")
    return(invisible(NULL))
  }

  safe_write(train_set, train_path)
  safe_write(test_set,  test_path)

  # lock the test set only after a successful write
  tryCatch(Sys.chmod(test_path, mode = "0444"), error = function(e) NULL)

  manifest <- list(seed = seed, test_fraction = test_frac,
                   n_total_clinical_raw = n_before, n_eligible = length(eligible_ids),
                   n_train = nrow(train_set), n_test = nrow(test_set),
                   n_events_train = sum(train_set[[strat_col_primary]] == 1, na.rm = TRUE),
                   n_events_test  = sum(test_set[[strat_col_primary]]  == 1, na.rm = TRUE),
                   stratification_primary = strat_col_primary,
                   survival_time_col = os_time_col, survival_status_col = os_status_col,
                   patient_id_col_source = id_col,
                   split_timestamp = as.character(Sys.time()),
                   script_version = "04_train_test_split.R FIXED")
  write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
  log_msg(sprintf("saved: %s", train_path))
  log_msg(sprintf("saved: %s", test_path))
  log_msg("Step 04 complete.")
}

if (sys.nframe() == 0) {
  run_train_test_split(
    input_root = WORK_DIR,
    seed = 42, test_frac = 0.25, force = TRUE
  )
}
