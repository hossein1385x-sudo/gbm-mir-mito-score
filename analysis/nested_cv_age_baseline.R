#!/usr/bin/env Rscript
# 01_age_baseline_deltaC.R
#
# Question answered:
#   Do the miRNAs add anything to AGE alone, when the WHOLE pipeline
#   (filter -> MDS -> DESeq2 -> univariate screen -> Elastic Net ->
#   frequency consolidation -> refit) is re-run inside every outer fold?

# =============================================================================
# FIX 0: set the temporary directory before anything else (fixes the temp-dir error)
# =============================================================================

source(here::here("config.R"))  # defines WORK_DIR (see config.R)

output_dir <- paste0(WORK_DIR, "/")
tmp_dir    <- file.path(output_dir, "r_temp")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir)
Sys.setenv(TEMP   = tmp_dir)
Sys.setenv(TMP    = tmp_dir)
options(tmpdir = tmp_dir)

# check that temp works
.test_tmp <- tempfile()
if (!file.create(.test_tmp)) {
  stop("cannot create a temporary file. Check the temp path manually: ", tmp_dir)
}
file.remove(.test_tmp)
cat(sprintf("Temp dir OK: %s\n", tmp_dir))

# =============================================================================
suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(DESeq2)
  library(parallel)
})

# ============================== 0. CONFIG ====================================
PATH_MIRNA_COUNTS <- file.path(output_dir, "mirna_isoform_counts_matrix.csv")
PATH_RNA_COUNTS   <- file.path(output_dir, "raw_rnaseq_counts_matrix.rds")
PATH_CLINICAL     <- file.path(output_dir, "raw_clinical_cdr_integrated.rds")
PATH_MITO_GENES   <- file.path(output_dir, "mito_gene_list.csv")
PATH_TRAIN_IDS    <- file.path(output_dir, "train_set.csv")
OUT_DIR           <- file.path(output_dir, "01_age_baseline_deltaC")

K_OUTER   <- 5
R_REPEATS <- 10
SEED      <- 2026

MIN_COUNT     <- 10
MIN_SAMPLES   <- 3
DE_FDR        <- 0.05
DE_LFC        <- 1
CANDIDATE_MODE <- "DE"
UNI_P         <- 0.05
UNI_MIN_KEEP  <- 5
UNI_MAX_KEEP  <- 200
ALPHA         <- 0.5
INNER_K       <- 5
LAMBDA_CHOICE <- "lambda.min"
AGE_PENALIZED <- TRUE

USE_CONSOLIDATION <- TRUE
R_INNER      <- 20
SUBSAMPLE    <- 0.8
FREQ_THR     <- 0.70
MAX_FEATURES <- 10
MIN_FALLBACK <- 3

N_CORES         <- 4
DEBUG_N_TASKS   <- NULL
PERMUTE_OUTCOME <- FALSE

dir.create(file.path(OUT_DIR, "folds"), recursive = TRUE, showWarnings = FALSE)

# ============================== 1. LOAD DATA =================================
read_mat <- function(p) {
  if (!file.exists(p)) stop(sprintf("File not found: %s", p))
  m <- if (grepl("\\.rds$", p, ignore.case = TRUE)) readRDS(p) else as.matrix(read.csv(p, row.names = 1, check.names = FALSE))
  if (is.data.frame(m)) m <- as.matrix(m)
  storage.mode(m) <- "numeric"
  m
}

read_clin <- function(p) {
  if (!file.exists(p)) stop(sprintf("File not found: %s", p))
  df <- if (grepl("\\.rds$", p, ignore.case = TRUE)) readRDS(p) else read.csv(p, stringsAsFactors = FALSE)
  df <- as.data.frame(df)

  if (!"patient_id" %in% colnames(df)) {
    pc <- intersect(c("patient_barcode", "patient_id", "barcode", "case_id"), colnames(df))[1]
    if (!is.na(pc)) df$patient_id <- df[[pc]]
  }
  if (!"time" %in% colnames(df)) {
    tc <- intersect(c("survival_time_final", "survival_time", "OS.time", "time", "days_to_death"), colnames(df))[1]
    if (!is.na(tc)) df$time <- df[[tc]]
  }
  if (!"event" %in% colnames(df)) {
    ec <- intersect(c("status_final", "status", "OS", "event", "vital_status"), colnames(df))[1]
    if (!is.na(ec)) {
      ev <- df[[ec]]
      df$event <- if (is.character(ev)) ifelse(tolower(ev) %in% c("dead", "deceased", "1"), 1L, 0L) else as.integer(ev)
    }
  }
  if (!"age" %in% colnames(df)) {
    ac <- intersect(c("age_at_diagnosis", "age"), colnames(df))[1]
    if (!is.na(ac)) df$age <- df[[ac]]
  }
  df
}

normalize_tcga_id <- function(ids) {
  ids <- toupper(as.character(ids))
  ids <- gsub("\\.", "-", ids)
  m <- regmatches(ids, regexpr("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", ids))
  out <- ids
  out[nchar(m) > 0] <- m[nchar(m) > 0]
  out
}

# ---- Load ----
mir  <- read_mat(PATH_MIRNA_COUNTS)
rna  <- read_mat(PATH_RNA_COUNTS)
clin <- read_clin(PATH_CLINICAL)

cat(sprintf("mir (raw): %d × %d\n", nrow(mir), ncol(mir)))
cat(sprintf("rna: %d × %d\n", nrow(rna), ncol(rna)))
cat(sprintf("clin: %d rows\n", nrow(clin)))

# ============================== 1a. TRANSPOSE MIR ============================
looks_like_mirna <- function(x) {
  mean(grepl("^hsa-", x, ignore.case = TRUE)) > 0.5
}

if (looks_like_mirna(colnames(mir)) && !looks_like_mirna(rownames(mir))) {
  cat("mir is transposed: transposing.\n")
  mir <- t(mir)
} else if (looks_like_mirna(rownames(mir))) {
  cat("mir is already in the correct orientation.\n")
} else {
  stop("cannot determine which dimension of mir holds the miRNA names.")
}

cat(sprintf("mir (after transpose): %d × %d\n", nrow(mir), ncol(mir)))

# ============================== 1b. NORMALIZE IDs ============================
mir_ids_norm  <- normalize_tcga_id(colnames(mir))
rna_ids_norm  <- normalize_tcga_id(colnames(rna))
clin_ids_norm <- normalize_tcga_id(clin$patient_id)

common_ids <- Reduce(intersect, list(mir_ids_norm, rna_ids_norm, clin_ids_norm))
cat(sprintf("Common IDs: %d\n", length(common_ids)))

if (length(common_ids) == 0) stop("no shared ID found.")

mir  <- mir[,  match(common_ids, mir_ids_norm),  drop = FALSE]; colnames(mir) <- common_ids
rna  <- rna[,  match(common_ids, rna_ids_norm),  drop = FALSE]; colnames(rna) <- common_ids
clin <- clin[match(common_ids, clin_ids_norm), , drop = FALSE]; clin$patient_id <- common_ids

# ---- Filter by train_ids ----
if (!is.null(PATH_TRAIN_IDS) && file.exists(PATH_TRAIN_IDS)) {
  train_df <- read.csv(PATH_TRAIN_IDS, stringsAsFactors = FALSE)
  id_col <- intersect(c("patient_barcode", "patient_id", "barcode", "case_id"), colnames(train_df))[1]
  if (!is.na(id_col)) {
    train_ids <- normalize_tcga_id(train_df[[id_col]])
    ids_keep <- intersect(common_ids, train_ids)
    cat(sprintf("Restricted to train set: %d patients\n", length(ids_keep)))

    mir  <- mir[,  ids_keep, drop = FALSE]
    rna  <- rna[,  ids_keep, drop = FALSE]
    clin <- clin[match(ids_keep, clin$patient_id), , drop = FALSE]
    clin$patient_id <- ids_keep
  }
}

# ---- Filter incomplete clinical rows ----
complete_rows <- complete.cases(clin[, c("patient_id", "time", "event", "age")])
clin <- clin[complete_rows, , drop = FALSE]
ids_final <- clin$patient_id

mir <- mir[, ids_final, drop = FALSE]
rna <- rna[, ids_final, drop = FALSE]

stopifnot(all(clin$time > 0), all(clin$event %in% 0:1))
cat(sprintf("Patients: %d | events: %d\n", nrow(clin), sum(clin$event)))

# ============================== 1c. MITO GENES ===============================
cat("matching the mitochondrial genes...\n")

rna_rownames_full <- rownames(rna)
rna_rownames_noversion <- gsub("\\.[0-9]+$", "", rna_rownames_full)

cat(sprintf("   RNA sample (with version):    %s\n", paste(head(rna_rownames_full, 3), collapse = ", ")))
cat(sprintf("   RNA sample (no version): %s\n", paste(head(rna_rownames_noversion, 3), collapse = ", ")))

mito_df <- read.csv(PATH_MITO_GENES, stringsAsFactors = FALSE)
cat(sprintf("   mito columns: %s\n", paste(colnames(mito_df), collapse = ", ")))

mito_ensembl <- unique(na.omit(mito_df$ensembl_id))
common_ensembl <- intersect(rna_rownames_noversion, mito_ensembl)

mito_symbols <- unique(na.omit(mito_df$gene_symbol))
common_symbols <- intersect(rna_rownames_noversion, mito_symbols)

cat(sprintf("   match on ensembl_id:  %d\n", length(common_ensembl)))
cat(sprintf("   match on gene_symbol: %d\n", length(common_symbols)))

if (length(common_ensembl) >= length(common_symbols) && length(common_ensembl) >= 5) {
  cat("   using ensembl_id\n")
  mito_in_rna <- common_ensembl
  rna_match_idx <- match(mito_in_rna, rna_rownames_noversion)
} else if (length(common_symbols) >= 5) {
  cat("   using gene_symbol\n")
  mito_in_rna <- common_symbols
  rna_match_idx <- match(mito_in_rna, rna_rownames_noversion)
} else {
  cat("   no sufficient match found.\n")
  stop("the number of matched mitochondrial genes is below 5.")
}

rna <- rna[rna_match_idx, , drop = FALSE]
rownames(rna) <- mito_in_rna

cat(sprintf("mito genes matched: %d\n", length(mito_in_rna)))
cat(sprintf("Final: %d patients | %d miRNAs | %d mito genes\n",
            nrow(clin), nrow(mir), nrow(rna)))

if (PERMUTE_OUTCOME) {
  set.seed(SEED + 999)
  p <- sample(nrow(clin))
  clin$time  <- clin$time[p]
  clin$event <- clin$event[p]
  cat(">>> PERMUTED-OUTCOME NULL RUN <<<\n")
}

# ============================== 2. HELPERS ===================================
cindex <- function(time, event, lp) {
  as.numeric(concordance(Surv(time, event) ~ lp,
                         data = data.frame(time = time, event = event, lp = lp),
                         reverse = TRUE)$concordance)
}

geo_nz <- function(x) if (all(x == 0)) NA_real_ else exp(sum(log(x[x > 0])) / length(x))
ref_geomeans <- function(m) apply(m, 1, geo_nz)
size_factors <- function(m, g) {
  ok0 <- !is.na(g)
  sf <- apply(m, 2, function(x) { ok <- ok0 & x > 0; median(x[ok] / g[ok]) })
  sf[!is.finite(sf) | sf <= 0] <- 1
  sf
}

transform_mirna <- function(counts, geo) {
  sf <- size_factors(counts, geo)
  log2(sweep(counts, 2, sf, "/") + 1)
}

scale_train <- function(X) {
  sdv <- apply(X, 2, sd); sdv[sdv == 0 | is.na(sdv)] <- 1
  list(mu = colMeans(X), sd = sdv)
}
apply_scale <- function(X, s) {
  sweep(sweep(X, 2, s$mu[colnames(X)], "-"), 2, s$sd[colnames(X)], "/")
}

make_folds <- function(event, K, seed) {
  set.seed(seed)
  fold <- integer(length(event))
  for (ev in unique(event)) {
    idx <- sample(which(event == ev))
    fold[idx] <- rep(seq_len(K), length.out = length(idx))
  }
  fold
}

uni_cox_p <- function(X, time, event) {
  apply(X, 2, function(x) {
    if (sd(x) == 0) return(1)
    s <- tryCatch(summary(coxph(Surv(time, event) ~ x)), error = function(e) NULL)
    if (is.null(s)) 1 else s$coefficients[1, "Pr(>|z|)"]
  })
}

en_select <- function(Xc, time, event, age_z = NULL) {
  x <- if (is.null(age_z)) Xc else cbind(Xc, age = age_z)
  if (ncol(x) == 1) {
    b <- tryCatch(coef(coxph(Surv(time, event) ~ x[, 1])), error = function(e) 0)
    return(setNames(as.numeric(b), colnames(x)))
  }
  pf <- rep(1, ncol(x))
  if (!AGE_PENALIZED && !is.null(age_z)) pf[ncol(x)] <- 0
  fit <- tryCatch(
    cv.glmnet(x, cbind(time = time, status = event), family = "cox",
              alpha = ALPHA, nfolds = INNER_K, penalty.factor = pf),
    error = function(e) NULL)
  if (is.null(fit)) return(setNames(rep(0, ncol(x)), colnames(x)))
  cf <- as.matrix(coef(fit, s = LAMBDA_CHOICE))
  setNames(cf[, 1], rownames(cf))
}

get_candidates <- function(mir_tr, rna_tr) {
  keep <- rowSums(mir_tr >= MIN_COUNT) >= MIN_SAMPLES
  mir_f <- mir_tr[keep, , drop = FALSE]
  if (CANDIDATE_MODE == "all") return(list(cand = rownames(mir_f), keep = rownames(mir_f)))

  sfr <- tryCatch(estimateSizeFactorsForMatrix(rna_tr),
                  error = function(e) size_factors(rna_tr, ref_geomeans(rna_tr)))
  lr  <- log2(sweep(rna_tr, 2, sfr, "/") + 1)
  z   <- t(scale(t(lr)))
  z   <- z[complete.cases(z), , drop = FALSE]
  mds <- colMeans(z)
  grp <- factor(ifelse(mds > median(mds), "high", "low"), levels = c("low", "high"))

  dds <- DESeqDataSetFromMatrix(round(mir_f), data.frame(group = grp, row.names = colnames(mir_f)), ~ group)
  dds <- tryCatch(DESeq(dds, sfType = "poscounts", quiet = TRUE),
                  error = function(e) DESeq(dds, sfType = "poscounts", fitType = "mean", quiet = TRUE))
  res <- as.data.frame(results(dds, contrast = c("group", "high", "low")))
  ok  <- !is.na(res$padj) & res$padj < DE_FDR & abs(res$log2FoldChange) >= DE_LFC
  list(cand = rownames(res)[ok], keep = rownames(mir_f))
}

# ============================== 3. PIPELINE ==================================
fit_pipeline <- function(idx_tr, seed) {
  set.seed(seed)
  mir_tr <- mir[, idx_tr, drop = FALSE]
  rna_tr <- rna[, idx_tr, drop = FALSE]
  cl     <- clin[idx_tr, ]
  n      <- nrow(cl)

  cc <- get_candidates(mir_tr, rna_tr)
  if (length(cc$cand) < 2) stop("fewer than 2 candidates in this fold")

  geo    <- ref_geomeans(mir_tr[cc$keep, , drop = FALSE])
  Xall   <- t(transform_mirna(mir_tr[cc$keep, , drop = FALSE], geo))
  Xtr    <- Xall[, cc$cand, drop = FALSE]
  sc     <- scale_train(Xtr)
  Xtr_z  <- apply_scale(Xtr, sc)
  age_mu <- mean(cl$age); age_sd <- sd(cl$age)
  age_z  <- (cl$age - age_mu) / age_sd

  select_once <- function(rows) {
    p <- uni_cox_p(Xtr_z[rows, , drop = FALSE], cl$time[rows], cl$event[rows])
    keepf <- names(p)[p < UNI_P]
    if (length(keepf) < UNI_MIN_KEEP) keepf <- names(sort(p))[seq_len(min(UNI_MIN_KEEP, length(p)))]
    if (length(keepf) > UNI_MAX_KEEP) keepf <- names(sort(p[keepf]))[seq_len(UNI_MAX_KEEP)]
    b <- en_select(Xtr_z[rows, keepf, drop = FALSE], cl$time[rows], cl$event[rows], age_z[rows])
    names(b)[b != 0 & names(b) != "age"]
  }

  if (USE_CONSOLIDATION) {
    picks <- unlist(lapply(seq_len(R_INNER), function(i)
      select_once(sample(n, round(SUBSAMPLE * n)))))
    fr  <- sort(table(factor(picks, levels = cc$cand)) / R_INNER, decreasing = TRUE)
    sel <- names(fr)[fr >= FREQ_THR]
    sel <- head(sel, MAX_FEATURES)
    if (length(sel) == 0) sel <- names(fr)[seq_len(min(MIN_FALLBACK, length(fr)))]
  } else {
    sel <- select_once(seq_len(n))
    if (length(sel) == 0) stop("Elastic Net selected nothing in this fold")
  }

  b_full <- en_select(Xtr_z[, sel, drop = FALSE], cl$time, cl$event, age_z)
  b_mir  <- en_select(Xtr_z[, sel, drop = FALSE], cl$time, cl$event, age_z = NULL)
  b_age  <- coef(coxph(Surv(cl$time, cl$event) ~ age_z))

  list(sel = sel, keep = cc$keep, geo = geo, sc = sc,
       age_mu = age_mu, age_sd = age_sd,
       b_full = b_full, b_mir = b_mir, b_age = as.numeric(b_age))
}

predict_pipeline <- function(m, idx_va) {
  Xva   <- t(transform_mirna(mir[m$keep, idx_va, drop = FALSE], m$geo))[, m$sel, drop = FALSE]
  Xva_z <- apply_scale(Xva, m$sc)
  a_z   <- (clin$age[idx_va] - m$age_mu) / m$age_sd
  list(
    lp_full = drop(Xva_z %*% m$b_full[m$sel]) + a_z * m$b_full[["age"]],
    lp_mir  = drop(Xva_z %*% m$b_mir[m$sel]),
    lp_age  = a_z * m$b_age
  )
}

run_task <- function(task) {
  f <- file.path(OUT_DIR, "folds", sprintf("rep%02d_fold%d.rds", task$rep, task$fold))
  if (file.exists(f)) return(readRDS(f))
  out <- tryCatch({
    m <- fit_pipeline(task$idx_tr, task$seed)
    p <- predict_pipeline(m, task$idx_va)
    tv <- clin$time[task$idx_va]; ev <- clin$event[task$idx_va]
    data.frame(rep = task$rep, fold = task$fold, ok = TRUE, msg = "",
               n_train = length(task$idx_tr), n_valid = length(task$idx_va),
               events_valid = sum(ev), n_sel = length(m$sel),
               C_age   = cindex(tv, ev, p$lp_age),
               C_mirna = cindex(tv, ev, p$lp_mir),
               C_full  = cindex(tv, ev, p$lp_full),
               sel = paste(m$sel, collapse = ";"), stringsAsFactors = FALSE)
  }, error = function(e)
    data.frame(rep = task$rep, fold = task$fold, ok = FALSE, msg = conditionMessage(e),
               n_train = NA, n_valid = NA, events_valid = NA, n_sel = NA,
               C_age = NA, C_mirna = NA, C_full = NA, sel = "", stringsAsFactors = FALSE))
  saveRDS(out, f)
  out
}

# ============================== 4. RUN =======================================
tasks <- list()
for (r in seq_len(R_REPEATS)) {
  fold <- make_folds(clin$event, K_OUTER, SEED + r)
  for (k in seq_len(K_OUTER)) {
    tasks[[length(tasks) + 1]] <- list(rep = r, fold = k,
                                       idx_tr = which(fold != k), idx_va = which(fold == k),
                                       seed = SEED * 100 + r * 10 + k)
  }
}
if (!is.null(DEBUG_N_TASKS)) tasks <- tasks[seq_len(DEBUG_N_TASKS)]
cat(sprintf("\nRunning %d outer folds on %d core(s)...\n", length(tasks), N_CORES))

t0 <- Sys.time()
if (N_CORES > 1) {
  cl_ <- makeCluster(N_CORES, outfile = "")
  clusterExport(cl_, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
  clusterEvalQ(cl_, {
    suppressPackageStartupMessages({
      library(survival); library(glmnet); library(DESeq2)
    })
    Sys.setenv(TMPDIR = file.path(WORK_DIR, "r_temp"))
    Sys.setenv(TEMP   = file.path(WORK_DIR, "r_temp"))
    Sys.setenv(TMP    = file.path(WORK_DIR, "r_temp"))
    options(tmpdir = file.path(WORK_DIR, "r_temp"))
    NULL
  })
  results <- parLapplyLB(cl_, tasks, run_task)
  stopCluster(cl_)
} else {
  results <- lapply(seq_along(tasks), function(i) {
    cat(sprintf("  fold %d/%d\n", i, length(tasks))); run_task(tasks[[i]]) })
}
cat("Elapsed:", format(Sys.time() - t0), "\n")

res <- do.call(rbind, results)
write.csv(res, file.path(OUT_DIR, "fold_level_results.csv"), row.names = FALSE)
if (any(!res$ok)) { cat("FAILED folds:\n"); print(res[!res$ok, c("rep", "fold", "msg")]) }
res_ok <- res[res$ok, ]

# ============================== 5. SUMMARY ===================================
res_ok$dC_full_vs_age  <- res_ok$C_full  - res_ok$C_age
res_ok$dC_mirna_vs_age <- res_ok$C_mirna - res_ok$C_age

summ_nb <- function(d, K = K_OUTER) {
  J <- length(d); m <- mean(d); s <- sd(d)
  se_naive <- s / sqrt(J)
  se_nb    <- sqrt((1 / J + 1 / (K - 1)) * s^2)
  tq <- qt(0.975, J - 1)
  data.frame(mean = m, sd = s, se_naive = se_naive, se_NB = se_nb,
             ci_lo = m - tq * se_nb, ci_hi = m + tq * se_nb,
             t = m / se_nb, p = 2 * pt(-abs(m / se_nb), J - 1),
             prop_pos = mean(d > 0), n_folds = J)
}

summary_tbl <- do.call(rbind, lapply(
  c("C_age", "C_mirna", "C_full", "dC_mirna_vs_age", "dC_full_vs_age"),
  function(v) cbind(metric = v, summ_nb(res_ok[[v]]))))
print(summary_tbl, digits = 3, row.names = FALSE)
write.csv(summary_tbl, file.path(OUT_DIR, "summary_NB.csv"), row.names = FALSE)

rep_means <- aggregate(cbind(C_age, C_mirna, C_full, dC_full_vs_age) ~ rep, res_ok, mean)
write.csv(rep_means, file.path(OUT_DIR, "per_repeat_means.csv"), row.names = FALSE)

sel_freq <- sort(table(unlist(strsplit(res_ok$sel, ";"))), decreasing = TRUE) / nrow(res_ok)
write.csv(data.frame(miRNA = names(sel_freq), freq = as.numeric(sel_freq)),
          file.path(OUT_DIR, "selection_frequency_fully_nested.csv"), row.names = FALSE)

png(file.path(OUT_DIR, "deltaC.png"), width = 1600, height = 800, res = 150)
par(mfrow = c(1, 2))
boxplot(res_ok[, c("C_age", "C_mirna", "C_full")], names = c("Age only", "miRNA only", "Age + miRNA"),
        ylab = "Outer-fold C-index", main = "Fully nested CV"); abline(h = 0.5, lty = 3)
hist(res_ok$dC_full_vs_age, breaks = 20, col = "grey85",
     xlab = "Delta C (Age+miRNA minus Age-only)", main = "Incremental value of miRNAs")
abline(v = 0, lty = 2, col = "red"); abline(v = mean(res_ok$dC_full_vs_age), col = "blue", lwd = 2)
dev.off()

cat("\nDone. Key file: summary_NB.csv (row 'dC_full_vs_age').\n")
