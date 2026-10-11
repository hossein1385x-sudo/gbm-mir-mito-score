#!/usr/bin/env Rscript
# =============================================================================
# Cross-platform sensitivity analysis of the signature miRNAs in CGGA (corrected C-index)
# =============================================================================
rm(list = ls())

source(here::here("config.R"))  # defines WORK_DIR (see config.R)

ROOT <- WORK_DIR
OUT  <- file.path(ROOT, "CROSSPLATFORM_SENS")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(survival); library(data.table) })

d <- as.data.frame(fread(file.path(ROOT, "cgga_model_input/cgga_flat.csv")))
miRNAs <- c("hsa-miR-196a-5p","hsa-miR-675-5p","hsa-miR-221-3p",
            "hsa-miR-31-5p","hsa-miR-181a-2-3p")

d <- d[complete.cases(d[, c("time","status", miRNAs)]), ]
cat("n =", nrow(d), " events =", sum(d$status), "\n\n")

# TCGA-locked coefficients
beta <- c(
  "hsa-miR-196a-5p"    =  0.2731,
  "hsa-miR-675-5p"     =  0.1878,
  "hsa-miR-221-3p"     =  0.07998,
  "hsa-miR-31-5p"      =  0.1731,
  "hsa-miR-181a-2-3p"  = -0.1546
)

mu_tcga <- c(4.249, 1.684, 11.234, 8.085, 10.206)
sd_tcga <- c(2.627, 1.703, 1.397, 1.914, 0.925)

# ---------- Method A: TCGA-locked ----------
z_locked <- sapply(seq_along(miRNAs), function(i) {
  (log2(d[[miRNAs[i]]] + 1) - mu_tcga[i]) / sd_tcga[i]
})
colnames(z_locked) <- miRNAs
risk_locked <- as.numeric(z_locked %*% beta[miRNAs])

# ---------- Method B: within-CGGA z ----------
z_inner <- scale(log2(as.matrix(d[, miRNAs]) + 1))
risk_inner <- as.numeric(z_inner %*% beta[miRNAs])

# ---------- Method C: rank-based ----------
rank_mat <- apply(as.matrix(d[, miRNAs]), 2, function(x) rank(x)/length(x))
risk_rank <- as.numeric(rank_mat %*% beta[miRNAs])

# ---------- C-index: use reverse = TRUE ----------
# Because risk score: higher = worse survival (shorter time)
cind <- function(r, y) {
  concordance(Surv(y[,1], y[,2]) ~ r, reverse = TRUE)$concordance
}

Y <- cbind(d$time, d$status)
res <- data.frame(
  method   = c("TCGA-locked scaling", "Within-CGGA z", "Rank-based"),
  C_index  = round(c(
    cind(risk_locked, Y),
    cind(risk_inner,  Y),
    cind(risk_rank,   Y)
  ), 4)
)
print(res)

# ---------- Verification against pipeline risk score ----------
rs_pipe <- read.csv(file.path(ROOT, "supplementary_final/CGGA_sample_scores.csv"),
                    stringsAsFactors = FALSE)
merged <- merge(data.frame(id = d$id, our_risk = risk_locked),
                rs_pipe, by.x = "id", by.y = "CGGA_ID")
cat("\nCorrelation with pipeline risk_score:",
    round(cor(merged$our_risk, merged$risk_score, use = "complete.obs"), 4), "\n")
cat("C-index of pipeline risk_score (reverse=TRUE):",
    round(concordance(Surv(OS, status) ~ risk_score,
                      data = merged, reverse = TRUE)$concordance, 4), "\n")

write.csv(res, file.path(OUT, "crossplatform_C_index_corrected.csv"),
          row.names = FALSE)
