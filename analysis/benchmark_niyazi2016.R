#!/usr/bin/env Rscript
# benchmark_niyazi2016.R — 4-miRNA signature (let-7a, let-7b, miR-125a, miR-615)
# Coefficients taken verbatim from Table 2 of the paper.

suppressPackageStartupMessages({ library(survival); library(readr) })
source(here::here("config.R"))   # defines WORK_DIR (see README)
OUT <- WORK_DIR
BENCH_OUT <- file.path(OUT, "benchmark_out_niyazi2016")
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

# ---- Niyazi2016: coefficients exactly from Table 2 ----
niyazi <- data.frame(
  signature = "Niyazi2016",
  feature   = c("hsa-let-7a-5p","hsa-let-7b-5p","hsa-miR-125a-5p","hsa-miR-615-5p"),
  coef      = c(0.5059587, -0.9669152, -0.2821517, 0.3254795),
  alt_names = c("", "", "", ""),
  stringsAsFactors = FALSE)

cat("=== Niyazi2016 coefficients (from Table 2) ===\n")
print(niyazi[, c("feature","coef")], row.names = FALSE)

# coverage
found <- niyazi$feature %in% colnames(X)
cat(sprintf("\nCoverage in Train: %d / 4\n", sum(found)))
if (!all(found)) cat("  missing:", paste(niyazi$feature[!found], collapse=", "), "\n")

# ---- Run ----
res <- run_benchmark(
  cohorts      = cohorts,
  own_scores   = own_scores,
  signatures   = niyazi,
  train_cohort = "Train",
  eval_cohorts = c("Train"),
  min_coverage = 0.75,   # 3/4 = 0.75
  B            = 1000,
  out_dir      = BENCH_OUT,
  seed         = 2026)

cat("\nDone:", BENCH_OUT, "\n")
