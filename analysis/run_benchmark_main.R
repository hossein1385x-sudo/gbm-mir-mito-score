#!/usr/bin/env Rscript
# =============================================================================
# run_benchmark_main.R — Benchmark of the miR-Mito model against published signatures
# CGGA expression matrix contains only the real miRNAs
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(readr)
})

# ---------------------------------------------------------------- 0) paths
source(here::here("config.R"))   # defines WORK_DIR (see README)
OUT <- WORK_DIR
if (!dir.exists(OUT)) stop("OUT does not exist: ", OUT)
cat("Working in:", OUT, "\n")

src_file <- here::here("analysis", "benchmark", "benchmark_signatures.R")
if (!file.exists(src_file)) stop("benchmark_signatures.R not found: ", src_file)
source(src_file)
cat("Loaded functions from:", src_file, "\n\n")

BENCH_OUT <- file.path(OUT, "benchmark_out")
dir.create(BENCH_OUT, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- override check_scale
check_scale <- function(cohorts) {
  cat("--- expression scale check ---\n")
  for (nm in names(cohorts)) {
    x <- cohorts[[nm]]$expr
    if (is.null(x) || length(x) == 0) {
      cat(sprintf("%-8s EMPTY — skipped\n", nm)); next
    }
    if (all(is.na(x))) {
      cat(sprintf("%-8s all NA — skipped (no expression; signature eval impossible)\n", nm))
      next
    }
    q <- quantile(x, c(0, .5, .99, 1), na.rm = TRUE)
    cat(sprintf("%-8s n=%d p=%d | min=%.2f median=%.2f p99=%.2f max=%.2f | NA=%.1f%%\n",
                nm, nrow(x), ncol(x), q[1], q[2], q[3], q[4], 100 * mean(is.na(x))))
    if (!is.na(q[4]) && q[4] > 50)
      warning(nm, ": max expression > 50 — probably NOT log2-scale")
  }
  invisible(NULL)
}

# ---------------------------------------------------------------- 1) Train
cat("=== 1) Loading Train ===\n")
tr  <- readRDS(file.path(OUT, "train_full_preprocessed.rds"))
sig <- readRDS(file.path(OUT, "signature_locked_latest.rds"))
beta <- setNames(sig$coefficients$Coefficient, sig$coefficients$Feature)
X <- as.matrix(tr$X_full); y <- tr$y_full
cat("  Xtr dim:", nrow(X), "x", ncol(X), "\n")

ts <- read_csv(file.path(OUT, "train_set.csv"), show_col_types = FALSE)
ts$Age_years <- ts$age_at_diagnosis / 365.25
age_tr <- setNames(ts$Age_years, ts$patient_barcode)

z1 <- (X[,"hsa-miR-196a-5p"] - sig$scaling$mirna$center["hsa-miR-196a-5p"]) /
  sig$scaling$mirna$scale ["hsa-miR-196a-5p"]
z2 <- (X[,"hsa-miR-675-5p"]  - sig$scaling$mirna$center["hsa-miR-675-5p"])  /
  sig$scaling$mirna$scale ["hsa-miR-675-5p"]
zA <- (age_tr[rownames(X)] - sig$scaling$clinical$center["Age"]) /
  sig$scaling$clinical$scale ["Age"]
train_lp <- setNames(
  as.numeric(beta["hsa-miR-196a-5p"]*z1 + beta["hsa-miR-675-5p"]*z2 + beta["Age"]*zA),
  rownames(X))

ok <- is.finite(y[,"time"]) & y[,"time"] > 0 & !is.na(age_tr[rownames(X)])
cohorts <- list(Train = list(
  expr  = X[ok, , drop = FALSE],
  time  = as.numeric(y[ok, "time"]),
  event = as.integer(y[ok, "status"]),
  age   = as.numeric(age_tr[rownames(X)[ok]]) ))
own_scores <- list(Train = train_lp[rownames(cohorts$Train$expr)])
cat("  Train n =", nrow(cohorts$Train$expr),
    "| events =", sum(cohorts$Train$event), "\n\n")

# ---------------------------------------------------------------- 2) TCGA matrix
cat("=== 2) Loading TCGA expression ===\n")
raw_path <- file.path(OUT, "GBM_mirna_matrix.csv")
if (!file.exists(raw_path)) stop("GBM_mirna_matrix.csv not found in ", OUT)

Xraw <- read.csv(raw_path, check.names = FALSE, row.names = 1)
cat("  raw dim:", nrow(Xraw), "x", ncol(Xraw), "\n")

looks_mirna_rows   <- mean(grepl("^hsa-", rownames(Xraw))) > 0.9
looks_patient_rows <- mean(grepl("^TCGA", rownames(Xraw))) > 0.9
if (looks_mirna_rows && !looks_patient_rows) {
  cat("  → transposing (rows are miRNAs)\n")
  Xtcga <- t(as.matrix(Xraw))
} else {
  Xtcga <- as.matrix(Xraw)
}
cat("  oriented: ", nrow(Xtcga), " samples x", ncol(Xtcga), " features\n")

if (grepl("\\.", rownames(Xtcga)[1]))
  rownames(Xtcga) <- gsub("\\.", "-", rownames(Xtcga))
Xtcga <- Xtcga[!duplicated(rownames(Xtcga)), , drop = FALSE]

common <- intersect(colnames(Xtcga), colnames(X))
cat("  overlap with Xtr:", length(common), "miRNA\n")
cat("  n mature (-5p/-3p):", sum(grepl("-[35]p$", colnames(Xtcga))), "\n")
cat("  sample id example:", rownames(Xtcga)[1], "\n\n")

# ---------------------------------------------------------------- 3) Test
cat("=== 3) Loading Test (TCGA test set) ===\n")
tp <- read_csv(file.path(OUT, "test_predictions.csv"), show_col_types = FALSE)
cat("  test_predictions: n =", nrow(tp), "| events =", sum(tp$event), "\n")

te_ids_in_expr <- intersect(tp$patient_id, rownames(Xtcga))
cat("  test IDs found in Xtcga:", length(te_ids_in_expr), "/", nrow(tp), "\n")

use_test_expr <- length(common) >= 100 && length(te_ids_in_expr) >= 20
if (use_test_expr) {
  Xte <- Xtcga[te_ids_in_expr, common, drop = FALSE]
  Xte_full <- matrix(NA_real_, nrow = nrow(Xte), ncol = ncol(X),
                     dimnames = list(rownames(Xte), colnames(X)))
  Xte_full[, common] <- Xte
  test_time  <- tp$survival_time[match(te_ids_in_expr, tp$patient_id)]
  test_event <- tp$event        [match(te_ids_in_expr, tp$patient_id)]
  cohorts$Test <- list(
    expr  = Xte_full,
    time  = as.numeric(test_time),
    event = as.integer(test_event),
    age   = NULL )
  cat("  → using real expression for Test\n\n")
} else {
  cohorts$Test <- list(
    expr  = matrix(NA_real_, nrow = nrow(tp), ncol = ncol(X),
                   dimnames = list(tp$patient_id, colnames(X))),
    time  = as.numeric(tp$survival_time),
    event = as.integer(tp$event),
    age   = NULL )
  cat("  → WARNING: expression overlap too low; Test will be OWN-only\n\n")
}
own_scores$Test <- setNames(tp$risk_score, tp$patient_id)

# ---------------------------------------------------------------- 4) CGGA
cat("=== 4) Loading CGGA ===\n")
cgga_path <- file.path(OUT, "cgga_model_input/cgga_flat.csv")
if (file.exists(cgga_path)) {
  cg <- read_csv(cgga_path, show_col_types = FALSE)
  cat("  cgga_flat: n =", nrow(cg), "| events =", sum(cg$status, na.rm = TRUE), "\n")
  mirna_cols_cg <- grep("^hsa-", names(cg), value = TRUE)
  cat("  available miRNAs:", paste(mirna_cols_cg, collapse = ", "), "\n")

  # OWN 2-miRNA score in CGGA
  z1cg <- (cg$`hsa-miR-196a-5p` - sig$scaling$mirna$center["hsa-miR-196a-5p"]) /
    sig$scaling$mirna$scale ["hsa-miR-196a-5p"]
  z2cg <- (cg$`hsa-miR-675-5p`  - sig$scaling$mirna$center["hsa-miR-675-5p"])  /
    sig$scaling$mirna$scale ["hsa-miR-675-5p"]
  lp_cg <- beta["hsa-miR-196a-5p"]*z1cg + beta["hsa-miR-675-5p"]*z2cg

  # expr matrix with the real miRNAs only (not 2196 NA-filled columns)
  ok_cg <- is.finite(cg$time) & cg$time > 0
  Xcg_real <- as.matrix(cg[ok_cg, mirna_cols_cg, drop = FALSE])
  rownames(Xcg_real) <- cg$id[ok_cg]

  cohorts$CGGA <- list(
    expr  = Xcg_real,
    time  = as.numeric(cg$time[ok_cg]),
    event = as.integer(cg$status[ok_cg]),
    age   = as.numeric(cg$age[ok_cg]) )
  own_scores$CGGA <- setNames(lp_cg[ok_cg], cg$id[ok_cg])
  cat("  CGGA ready: n =", nrow(cohorts$CGGA$expr),
      "| p =", ncol(cohorts$CGGA$expr), "miRNAs\n\n")
} else {
  cat("  cgga_flat not found; skipping CGGA\n\n")
}

# ---------------------------------------------------------------- 5) signatures
cat("=== 5) Reading published signatures ===\n")
sig_path <- file.path(OUT, "signatures_published.csv")
if (!file.exists(sig_path)) sig_path <- here::here("analysis", "benchmark", "signatures_published.csv")
if (!file.exists(sig_path))
  stop("signatures_published.csv not found (looked in WORK_DIR and analysis/benchmark/). Columns: signature, feature, coef, alt_names")
signatures <- read_signatures(sig_path)
cat("  signatures:", paste(unique(signatures$signature), collapse = ", "), "\n")
cat("  total rows:", nrow(signatures), "\n\n")

# ---------------------------------------------------------------- 6) benchmark
cat("=== 6) Running benchmark ===\n")
valid_cohorts <- names(cohorts)[vapply(cohorts, function(co) {
  !is.null(co$expr) && !all(is.na(co$expr))
}, logical(1))]
cat("Eligible cohorts for signature benchmark:",
    paste(valid_cohorts, collapse = ", "), "\n")

eval_names <- intersect(valid_cohorts, c("Train", "Test", "CGGA"))

res <- NULL
if (length(eval_names) > 0) {
  res <- run_benchmark(
    cohorts      = cohorts[eval_names],
    own_scores   = own_scores[eval_names],
    signatures   = signatures,
    train_cohort = if ("Train" %in% eval_names) "Train" else NULL,
    eval_cohorts = eval_names,
    min_coverage = 0.8,
    B            = 1000,
    out_dir      = BENCH_OUT,
    seed         = 2026 )
} else {
  cat("  ⚠ no cohort with real expression; skipping run_benchmark\n")
}

# ---------------------------------------------------------------- 7) OWN out-of-sample
cat("\n=== 7) OWN model C-index (all cohorts) ===\n")
cidx <- function(time, event, score) {
  ok <- is.finite(time) & is.finite(event) & is.finite(score) & time > 0
  if (sum(ok) < 10 || sum(event[ok]) < 5) return(NA_real_)
  survival::concordance(survival::Surv(time[ok], event[ok]) ~ score[ok],
                        reverse = TRUE)$concordance
}

own_rows <- list()
for (cn in c("Train", "CGGA")) {
  if (!is.null(cohorts[[cn]]) && !is.null(own_scores[[cn]])) {
    sc <- own_scores[[cn]]; co <- cohorts[[cn]]
    own_rows[[length(own_rows)+1]] <- data.frame(
      cohort = cn,
      n      = length(sc),
      events = sum(co$event, na.rm = TRUE),
      C_index = round(cidx(co$time, co$event, sc), 3),
      source  = if (cn == "Train") "locked score (in-sample)" else "cgga_flat.csv (out-of-sample)",
      stringsAsFactors = FALSE)
  }
}
if (!is.null(tp)) {
  own_rows[[length(own_rows)+1]] <- data.frame(
    cohort = "Test",
    n      = nrow(tp),
    events = sum(tp$event, na.rm = TRUE),
    C_index = round(cidx(tp$survival_time, tp$event, tp$risk_score), 3),
    source  = "test_predictions.csv (out-of-sample)",
    stringsAsFactors = FALSE)
}
own_summary <- do.call(rbind, own_rows)
if (!is.null(own_summary)) {
  print(own_summary, row.names = FALSE)
  write.csv(own_summary,
            file.path(BENCH_OUT, "OWN_out_of_sample_summary.csv"),
            row.names = FALSE)
}

# ---------------------------------------------------------------- 8) report
cat("\n=========================================================\n")
cat("DONE — outputs in:", BENCH_OUT, "\n")
cat("=========================================================\n")
for (f in list.files(BENCH_OUT)) cat("  -", f, "\n")
cat("\nInterpretation notes:\n")
cat("  - benchmark_coverage.csv         : miRNA coverage per cohort × signature\n")
cat("  - benchmark_performance.csv      : C-index + bootstrap CI\n")
cat("  - benchmark_delta.csv            : paired difference OWN vs each competitor\n")
cat("  - OWN_out_of_sample_summary.csv  : OWN C-index all cohorts\n")
cat("  - benchmark_forest.pdf           : forest plot\n")
cat("\n  Limitations (state in paper):\n")
cat("    * Test : miRNA expression not saved in pipeline → only OWN evaluated\n")
cat("    * CGGA : only 5 miRNAs available → published signatures coverage < 0.8\n")
cat("    * Train: full benchmark possible, but in-sample (optimistic)\n")
