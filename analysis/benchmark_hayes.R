#!/usr/bin/env Rscript
# benchmark_hayes.R — Hayes 2015 (9-miRNA) benchmark vs miR-Mito
# Exact coefficients from the paper (Table 2 + formula in section 3.2)

suppressPackageStartupMessages({ library(survival); library(readr) })

source(here::here("config.R"))   # defines WORK_DIR (see README)
OUT <- WORK_DIR
BENCH_OUT <- file.path(OUT, "benchmark_out_hayes")
dir.create(BENCH_OUT, showWarnings = FALSE, recursive = TRUE)

source(here::here("analysis", "benchmark", "benchmark_signatures.R"))

# ---------------------------------------------------------------- 1) Train
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

# ---------------------------------------------------------------- 2) Test (OWN only)
tp <- read_csv(file.path(OUT, "test_predictions.csv"), show_col_types = FALSE)

# ---------------------------------------------------------------- 3) CGGA (5 miRNAs)
cg <- read_csv(file.path(OUT, "cgga_model_input/cgga_flat.csv"), show_col_types = FALSE)
mirna_cg <- grep("^hsa-", names(cg), value = TRUE)
Xcg <- as.matrix(cg[, mirna_cg, drop = FALSE]); rownames(Xcg) <- cg$id
ok_cg <- is.finite(cg$time) & cg$time > 0
z1cg <- (cg$`hsa-miR-196a-5p` - sig$scaling$mirna$center["hsa-miR-196a-5p"]) / sig$scaling$mirna$scale["hsa-miR-196a-5p"]
z2cg <- (cg$`hsa-miR-675-5p`  - sig$scaling$mirna$center["hsa-miR-675-5p"])  / sig$scaling$mirna$scale["hsa-miR-675-5p"]
lp_cg <- as.numeric(beta["hsa-miR-196a-5p"]*z1cg + beta["hsa-miR-675-5p"]*z2cg)
cohorts$CGGA <- list(
  expr  = Xcg[ok_cg, , drop = FALSE],
  time  = cg$time[ok_cg],
  event = cg$status[ok_cg],
  age   = cg$age[ok_cg])
own_scores$CGGA <- setNames(lp_cg[ok_cg], cg$id[ok_cg])

# ---------------------------------------------------------------- 4) Hayes2015
hayes_map <- data.frame(
  article = c("miR-370","miR-124a","miR-145","miR-34a","miR-10b",
              "miR-148a","miR-222","miR-9","miR-182"),
  mature  = c("hsa-miR-370-3p","hsa-miR-124-3p","hsa-miR-145-5p","hsa-miR-34a-5p","hsa-miR-10b-5p",
              "hsa-miR-148a-3p","hsa-miR-222-3p","hsa-miR-9-5p","hsa-miR-182-5p"),
  coef    = c(-0.044, 0.062, -0.066, 0.005, 0.015, 0.092, 0.162, -0.032, -0.021),
  stringsAsFactors = FALSE)

miss <- !hayes_map$mature %in% colnames(X)
if (any(miss)) {
  cat("Some Hayes mature names not found exactly. Searching alternatives:\n")
  for (i in which(miss)) {
    base <- sub("([0-9]+)[a-z]$", "\\1", hayes_map$article[i])
    hits <- grep(paste0("^hsa-", base, "(-[35]p)?$"), colnames(X), value = TRUE, ignore.case = TRUE)
    cat(sprintf("  %-10s (want %s) -> %s\n", hayes_map$article[i], hayes_map$mature[i],
                if (length(hits) == 0) "NOT FOUND" else paste(hits, collapse = ", ")))
    if (length(hits) > 0) hayes_map$mature[i] <- hits[1]
  }
}

hayes <- data.frame(
  signature = "Hayes2015",
  feature   = hayes_map$mature,
  coef      = hayes_map$coef,
  alt_names = hayes_map$article,
  stringsAsFactors = FALSE)
hayes <- hayes[hayes$feature %in% colnames(X), ]

cat(sprintf("\n=== Hayes2015: %d / 9 miRNAs available in Train ===\n", nrow(hayes)))
print(hayes[, c("feature","coef")], row.names = FALSE)
cat("\n")

# ---------------------------------------------------------------- 5) Benchmark (Train + CGGA)
res <- run_benchmark(
  cohorts      = cohorts[c("Train","CGGA")],
  own_scores   = own_scores[c("Train","CGGA")],
  signatures   = hayes,
  train_cohort = "Train",
  eval_cohorts = c("Train","CGGA"),
  min_coverage = 0.8, B = 1000,
  out_dir      = BENCH_OUT, seed = 2026)

# ---------------------------------------------------------------- 6) OWN summary (all 3)
cidx <- function(t,e,s) {
  ok <- is.finite(t) & is.finite(e) & is.finite(s) & t > 0
  if (sum(ok) < 10 || sum(e[ok]) < 5) return(NA_real_)
  survival::concordance(survival::Surv(t[ok], e[ok]) ~ s[ok], reverse = TRUE)$concordance
}
own_tbl <- rbind(
  data.frame(cohort = "Train", n = nrow(cohorts$Train$expr), events = sum(cohorts$Train$event),
             C_index = round(cidx(cohorts$Train$time, cohorts$Train$event, own_scores$Train), 3),
             source = "locked (in-sample)"),
  data.frame(cohort = "CGGA", n = nrow(cohorts$CGGA$expr), events = sum(cohorts$CGGA$event),
             C_index = round(cidx(cohorts$CGGA$time, cohorts$CGGA$event, own_scores$CGGA), 3),
             source = "cgga_flat (out-of-sample)"),
  data.frame(cohort = "Test", n = nrow(tp), events = sum(tp$event),
             C_index = round(cidx(tp$survival_time, tp$event, tp$risk_score), 3),
             source = "test_predictions (out-of-sample)"))

cat("\n=== OWN model C-index ===\n"); print(own_tbl, row.names = FALSE)
write.csv(own_tbl, file.path(BENCH_OUT, "OWN_summary.csv"), row.names = FALSE)

cat("\nDone. Output:", BENCH_OUT, "\n")
