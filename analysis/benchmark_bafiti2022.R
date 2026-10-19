#!/usr/bin/env Rscript
# benchmark_bafiti2022.R — 3-miRNA signature (miR-20a, miR-21, miR-10a)
# All risky (+). Weights: relative RF importance (the paper did not publish betas).

suppressPackageStartupMessages({ library(survival); library(readr) })
source(here::here("config.R"))   # defines WORK_DIR (see README)
OUT <- WORK_DIR
BENCH_OUT <- file.path(OUT, "benchmark_out_bafiti2022")
dir.create(BENCH_OUT, showWarnings = FALSE, recursive = TRUE)

source(here::here("analysis", "benchmark", "benchmark_signatures.R"))

# ---- Train ----
tr  <- readRDS(file.path(OUT, "train_full_preprocessed.rds"))
sig <- readRDS(file.path(OUT, "signature_locked_latest.rds"))
beta <- setNames(sig$coefficients$Coefficient, sig$coefficients$Feature)
X <- as.matrix(tr$X_full); y <- tr$y_full
ts <- read_csv(file.path(OUT, "train_set.csv"), show_col_types = FALSE)
age_tr <- setNames(ts$age_at_diagnosis / 365.25, ts$patient_barcode)

z1 <- (X[,"hsa-miR-196a-5p"] - sig$scaling$mirna$center["hsa-miR-196a-5p"]) / sig$scaling$mirna$scale["hsa-miR-196a-5p"]
z2 <- (X[,"hsa-miR-675-5p"]  - sig$scaling$mirna$center["hsa-miR-675-5p"])  / sig$scaling$mirna$scale["hsa-miR-675-5p"]
zA <- (age_tr[rownames(X)] - sig$scaling$clinical$center["Age"]) / sig$scaling$clinical$scale["Age"]
train_lp <- setNames(as.numeric(beta["hsa-miR-196a-5p"]*z1 + beta["hsa-miR-675-5p"]*z2 + beta["Age"]*zA), rownames(X))

ok <- is.finite(y[,"time"]) & y[,"time"] > 0 & !is.na(age_tr[rownames(X)])
cohorts <- list(Train = list(
  expr  = X[ok, , drop = FALSE],
  time  = as.numeric(y[ok, "time"]),
  event = as.integer(y[ok, "status"]),
  age   = as.numeric(age_tr[rownames(X)[ok]])))
own_scores <- list(Train = train_lp[rownames(cohorts$Train$expr)])

# ---- Bafiti2022 signature: relative RF importance weights ----
imp <- c("hsa-miR-20a-5p"=9.41, "hsa-miR-21-5p"=7.02, "hsa-miR-10a-5p"=13.02)
bafiti <- data.frame(
  signature = "Bafiti2022",
  feature   = names(imp),
  coef      = as.numeric(imp / sum(imp)),   # normalized to sum=1
  alt_names = c("hsa-miR-20a", "hsa-miR-21", "hsa-miR-10a"),
  stringsAsFactors = FALSE)

cat("=== Bafiti2022 weights ===\n")
print(bafiti[, c("feature","coef")], row.names = FALSE)

# coverage check
found <- bafiti$feature %in% colnames(X)
cat(sprintf("\nCoverage in Train: %d / 3\n", sum(found)))
if (!all(found)) cat("Missing:", paste(bafiti$feature[!found], collapse=", "), "\n")

# ---- Run ----
res <- run_benchmark(
  cohorts      = cohorts,
  own_scores   = own_scores,
  signatures   = bafiti,
  train_cohort = "Train",
  eval_cohorts = c("Train"),
  min_coverage = 0.66,   # 2/3 = 0.667 (still evaluated if one miRNA is missing)
  B            = 1000,
  out_dir      = BENCH_OUT,
  seed         = 2026)

cat("\nDone:", BENCH_OUT, "\n")
