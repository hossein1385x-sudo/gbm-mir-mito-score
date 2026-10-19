#!/usr/bin/env Rscript
# =============================================================================
# benchmark_li2014.R - benchmark of the 5 subtype-specific signatures of Li et al. 2014
# Cox coefficients were not published -> unit weights (+1 risky / -1 protective)
# Evaluated in all three cohorts: Train, Test, CGGA
# Output: WORK_DIR/benchmark_out_li2014/  (read by combine_benchmarks.R)
# =============================================================================

suppressPackageStartupMessages({ library(survival); library(readr) })

source(here::here("config.R"))   # defines WORK_DIR (see README)
OUT <- WORK_DIR
BENCH_OUT <- file.path(OUT, "benchmark_out_li2014")
dir.create(BENCH_OUT, showWarnings = FALSE, recursive = TRUE)

source(here::here("analysis", "benchmark", "benchmark_signatures.R"))

rd <- function(p) utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
norm_id <- function(x) toupper(substr(gsub("[._]", "-", x), 1, 12))
norm_cgga <- function(x) {
  x <- toupper(gsub("[._]", "-", x))
  ifelse(grepl("^CGGA", x), x, paste0("CGGA-", x))
}

# =============================================================================
# 1) Train
# =============================================================================
cat("=== 1) Train ===\n")
tr  <- readRDS(file.path(OUT, "train_full_preprocessed.rds"))
sig <- readRDS(file.path(OUT, "signature_locked_latest.rds"))
beta <- setNames(sig$coefficients$Coefficient, sig$coefficients$Feature)
M1 <- "hsa-miR-196a-5p"; M2 <- "hsa-miR-675-5p"
X <- as.matrix(tr$X_full); y <- tr$y_full

ts <- rd(file.path(OUT, "train_set.csv"))
age_tr_raw <- ts$age_at_diagnosis
if (median(age_tr_raw, na.rm = TRUE) > 120) age_tr_raw <- age_tr_raw / 365.25
age_tr <- setNames(age_tr_raw, ts$patient_barcode)

zmi  <- function(v, f) (v - sig$scaling$mirna$center[f]) / sig$scaling$mirna$scale[f]
zage <- function(a) (a - sig$scaling$clinical$center["Age"]) / sig$scaling$clinical$scale["Age"]
lp_mirna <- function(m1, m2) as.numeric(beta[M1] * zmi(m1, M1) + beta[M2] * zmi(m2, M2))

ok <- is.finite(y[, "time"]) & y[, "time"] > 0 & !is.na(age_tr[rownames(X)])
cohorts <- list(Train = list(
  expr  = X[ok, , drop = FALSE],
  time  = as.numeric(y[ok, "time"]),
  event = as.integer(y[ok, "status"]),
  age   = as.numeric(age_tr[rownames(X)[ok]])))
rn <- rownames(cohorts$Train$expr)
lp_tr_mi <- setNames(lp_mirna(X[rn, M1], X[rn, M2]), rn)
lp_tr    <- lp_tr_mi + as.numeric(beta["Age"] * zage(age_tr[rn]))
own_scores <- list(Train = lp_tr)
extra      <- list(Train = list(OWN_miRNA_only = lp_tr_mi))
cat(sprintf("  n=%d events=%d\n", nrow(cohorts$Train$expr), sum(cohorts$Train$event)))

# =============================================================================
# 2) Test - with the mature-miRNA matrix
# =============================================================================
cat("\n=== 2) Test ===\n")
tp <- rd(file.path(OUT, "test_predictions.csv"))
tp$pid <- norm_id(tp$patient_id)

Xmat_path <- file.path(OUT, "mirna_isoform_rpm_matrix.csv")
if (file.exists(Xmat_path)) {
  Xr <- as.matrix(read.csv(Xmat_path, check.names = FALSE, row.names = 1))
  ids_raw <- gsub("[._]", "-", rownames(Xr))
  is_tumor <- nchar(ids_raw) <= 12 | substr(ids_raw, 14, 15) == "01"
  Xr <- Xr[is_tumor, , drop = FALSE]
  rownames(Xr) <- norm_id(rownames(Xr))
  Xr <- Xr[!duplicated(rownames(Xr)), , drop = FALSE]

  te_ids <- intersect(tp$pid, rownames(Xr))
  cat(sprintf("  test patients matched: %d / %d\n", length(te_ids), nrow(tp)))

  if (length(te_ids) >= 20) {
    tt <- tp[match(te_ids, tp$pid), ]
    Xte <- Xr[te_ids, , drop = FALSE]

    # Age
    age_te <- NULL
    for (f in c("test_set_LOCKED.csv", "test_set.csv")) {
      pth <- file.path(OUT, f)
      if (!file.exists(pth)) next
      ct <- rd(pth)
      idc <- intersect(c("patient_barcode", "patient_id"), names(ct))[1]
      if (is.na(idc) || !"age_at_diagnosis" %in% names(ct)) next
      a <- ct$age_at_diagnosis
      if (median(a, na.rm = TRUE) > 120) a <- a / 365.25
      age_te <- a[match(te_ids, norm_id(ct[[idc]]))]
      break
    }

    cohorts$Test <- list(
      expr  = Xte,
      time  = as.numeric(tt$survival_time),
      event = as.integer(tt$event),
      age   = age_te)
    own_scores$Test <- setNames(tt$risk_score, te_ids)

    if (all(c(M1, M2) %in% colnames(Xte))) {
      lp_te_mi <- setNames(lp_mirna(Xte[, M1], Xte[, M2]), te_ids)
      extra$Test <- list(OWN_miRNA_only = lp_te_mi)
    }
    cat(sprintf("  Test ready: n=%d events=%d p=%d\n", nrow(Xte), sum(tt$event), ncol(Xte)))
  }
} else cat("  [warn] mature matrix not found; Test skipped\n")

# =============================================================================
# 3) CGGA - with the full matrix (cgga_full_matrix.csv)
# =============================================================================
cat("\n=== 3) CGGA ===\n")
cg_full_path <- file.path(OUT, "cgga_full_matrix.csv")
cg_flat_path <- file.path(OUT, "cgga_model_input", "cgga_flat.csv")

if (file.exists(cg_full_path) && file.exists(cg_flat_path)) {
  cg_full <- as.matrix(read.csv(cg_full_path, check.names = FALSE, row.names = 1))
  cg_flat <- rd(cg_flat_path)
  flat_ids <- norm_cgga(cg_flat$id)
  full_ids <- norm_cgga(rownames(cg_full))
  common <- intersect(flat_ids, full_ids)
  cat(sprintf("  matched %d / %d | p=%d miRNAs\n", length(common), length(flat_ids), ncol(cg_full)))

  if (length(common) >= 100) {
    fi <- match(common, flat_ids)
    mi <- match(common, full_ids)

    cg_age <- cg_flat$age[fi]
    if (median(cg_age, na.rm = TRUE) > 120) cg_age <- cg_age / 365.25

    Xcg <- cg_full[mi, , drop = FALSE]
    rownames(Xcg) <- common

    t_cg <- cg_flat$time[fi]; e_cg <- cg_flat$status[fi]
    okc <- is.finite(t_cg) & t_cg > 0
    Xcg <- Xcg[okc, , drop = FALSE]
    ids_cg <- rownames(Xcg)

    flat_row <- match(ids_cg, flat_ids)
    lp_cg_mi <- setNames(
      lp_mirna(cg_flat[[M1]][flat_row], cg_flat[[M2]][flat_row]),
      ids_cg)
    lp_cg <- lp_cg_mi + as.numeric(beta["Age"] * zage(cg_age[okc]))

    cohorts$CGGA <- list(
      expr  = Xcg,
      time  = as.numeric(t_cg[okc]),
      event = as.integer(e_cg[okc]),
      age   = cg_age[okc])
    own_scores$CGGA <- lp_cg[ids_cg]
    extra$CGGA <- list(OWN_miRNA_only = lp_cg_mi[ids_cg])
    cat(sprintf("  CGGA ready: n=%d events=%d p=%d\n",
                length(ids_cg), sum(cohorts$CGGA$event), ncol(Xcg)))
  }
}

# =============================================================================
# 4) Li2014 signatures - unit weights
# =============================================================================
cat("\n=== 4) Building Li2014 signatures ===\n")

li <- data.frame(
  signature = c(
    rep("Li2014_Classical", 7),
    rep("Li2014_Neural", 8),
    rep("Li2014_Mesenchymal", 5),
    rep("Li2014_ProneuralGCIMP", 3),
    rep("Li2014_ProneuralNonGCIMP", 10)),
  feature = c(
    # Classical
    "hsa-miR-26a-5p","hsa-miR-767-3p","hsa-miR-153-3p","hsa-miR-31-5p","hsa-miR-222-3p",
    "hsa-miR-654-3p","hsa-miR-378c",
    # Neural
    "hsa-miR-222-3p","hsa-miR-422a","hsa-miR-662","hsa-miR-566","hsa-miR-24-3p",
    "hsa-miR-370-3p","hsa-miR-492","hsa-miR-629-5p",
    # Mesenchymal
    "hsa-miR-373-3p","hsa-miR-296-5p","hsa-miR-191-5p","hsa-miR-602","hsa-miR-223-3p",
    # Proneural-G-CIMP
    "hsa-miR-582-5p","hsa-miR-130a-3p","hsa-miR-195-5p",
    # Proneural-nonG-CIMP
    "hsa-miR-335-5p","hsa-miR-34a-5p","hsa-miR-581","hsa-miR-21-5p",
    "hsa-miR-361-5p","hsa-miR-145-5p","hsa-miR-143-5p","hsa-miR-378a-3p",
    "hsa-miR-182-5p","hsa-miR-183-5p"),
  coef = c(
    # Classical: 5 risky + 2 protective
    1,1,1,1,1, -1,-1,
    # Neural: 1 risky + 7 protective
    1, -1,-1,-1,-1,-1,-1,-1,
    # Mesenchymal: 4 risky + 1 protective
    1,1,1,1, -1,
    # Proneural-G-CIMP: 1 risky + 2 protective
    1, -1,-1,
    # Proneural-nonG-CIMP: 4 risky + 6 protective
    1,1,1,1, -1,-1,-1,-1,-1,-1),
  alt_names = c(
    # Classical
    "hsa-miR-26a","hsa-miR-767-3p","hsa-miR-153","hsa-miR-31","hsa-miR-222",
    "hsa-miR-654","hsa-miR-422b",
    # Neural
    "hsa-miR-222","hsa-miR-422a","hsa-miR-662","hsa-miR-566","hsa-miR-24",
    "hsa-miR-370","hsa-miR-492","hsa-miR-629",
    # Mesenchymal
    "hsa-miR-373","hsa-miR-296","hsa-miR-191","hsa-miR-602","hsa-miR-223",
    # Proneural-G-CIMP
    "hsa-miR-582","hsa-miR-130a","hsa-miR-195",
    # Proneural-nonG-CIMP
    "hsa-miR-335","hsa-miR-34a","hsa-miR-581","hsa-miR-21",
    "hsa-miR-361","hsa-miR-145","hsa-miR-143","hsa-miR-378",
    "hsa-miR-182","hsa-miR-183"),
  stringsAsFactors = FALSE)

# Coverage before running
cat("\n--- Coverage of Li2014 signatures ---\n")
for (cn in names(cohorts)) {
  cat(sprintf("\n[%s]\n", cn))
  avail <- colnames(cohorts[[cn]]$expr)
  for (s in unique(li$signature)) {
    sub <- li[li$signature == s, ]
    res <- resolve_features(sub, avail)
    n_ok <- sum(!is.na(res$matched))
    cat(sprintf("  %-28s %d / %d\n", s, n_ok, nrow(sub)))
  }
}

# =============================================================================
# 5) Benchmark
# =============================================================================
cat("\n=== 5) Running benchmark ===\n")
valid <- names(cohorts)[vapply(cohorts, function(co)
  !is.null(co$expr) && ncol(co$expr) > 0 && !all(is.na(co$expr)), logical(1))]

res <- run_benchmark(
  cohorts      = cohorts[valid],
  own_scores   = own_scores[valid],
  signatures   = li,
  train_cohort = "Train",
  eval_cohorts = valid,
  extra_scores = extra[valid],
  min_coverage = 0.8,
  B            = 1000,
  out_dir      = BENCH_OUT,
  seed         = 2026)

# =============================================================================
# 6) Summary table
# =============================================================================
cat("\n=== Summary ===\n")
perf <- res$performance
perf_sub <- perf[grepl("Li2014", perf$model) |
                 perf$model %in% c("OWN","OWN_miRNA_only","AGE_ONLY"), ]
print(perf_sub[, c("cohort","model","n","events","C","C_lo","C_hi")], row.names = FALSE)

cat("\nDone:", BENCH_OUT, "\n")
