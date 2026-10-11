#!/usr/bin/env Rscript
# =============================================================================
# Supplementary — Figure S5 (KM for CGGA) — Diagnostic + Robust
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table); library(survival); library(survminer); library(ggplot2)
})

input_dir <- paste0(WORK_DIR, "/")
sig_dir   <- file.path(input_dir, "signature_v2")
out_dir   <- file.path(input_dir, "supplementary_final")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(sprintf("[%s] %s\n",
                                     format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

# --- Load ---
expr_file <- file.path(input_dir,
                       "cgga_raw/exp cgga/CGGA.microRNA_array_198_gene_level.20200506.txt")
clin_file <- file.path(input_dir,
                       "cgga_raw/clin cgga/CGGA.microRNA_array_198_clinical.20200506.txt")
sc_file   <- file.path(sig_dir, "signature_v2_scaling.rds")
if (!file.exists(sc_file)) sc_file <- file.path(input_dir, "signature_scaling_params.rds")
feat_file <- file.path(sig_dir, "signature_features.csv")

expr <- read.delim(expr_file, header = TRUE, sep = "\t",
                   check.names = FALSE, stringsAsFactors = FALSE)
clin <- read.delim(clin_file, header = TRUE, sep = "\t",
                   check.names = FALSE, stringsAsFactors = FALSE)
sc   <- readRDS(sc_file)
feat <- fread(feat_file)

log_msg(sprintf("   expr: %d × %d | clin: %d × %d",
                nrow(expr), ncol(expr), nrow(clin), ncol(clin)))

# --- mapping CGGA ---
mirna_mapping <- list(
  "hsa-miR-196a-5p"   = "hsa-miR-196a",
  "hsa-miR-675-5p"    = "hsa-miR-675",
  "hsa-miR-221-3p"    = "hsa-miR-221#",
  "hsa-miR-31-5p"     = "hsa-miR-31",
  "hsa-miR-181a-2-3p" = "hsa-miR-181a-2#"
)

# --- shared samples ---
sample_ids <- as.character(colnames(expr)[-1])
clin$CGGA_ID <- as.character(clin$CGGA_ID)
common_ids <- intersect(sample_ids, clin$CGGA_ID)
log_msg(sprintf("   common IDs: %d", length(common_ids)))

# --- debug mapping ---
log_msg("")
log_msg("check of the mapping and coefficients:")
for (m in names(mirna_mapping)) {
  probe <- mirna_mapping[[m]]
  has_probe   <- probe %in% expr$microRNA_ID
  has_center  <- m %in% names(sc$mirna$center)
  has_scale   <- m %in% names(sc$mirna$scale)
  has_beta    <- m %in% feat$Feature
  cat(sprintf("  %-22s probe=%-20s | probe=%s | center=%s | scale=%s | beta=%s\n",
              m, probe,
              ifelse(has_probe,  "OK", "MISSING"),
              ifelse(has_center, "OK", "MISSING"),
              ifelse(has_scale,  "OK", "MISSING"),
              ifelse(has_beta,   "OK", "MISSING")))
}

# --- z-score ---
beta <- setNames(feat$Coefficient, feat$Feature)
mirnas_used <- intersect(names(mirna_mapping), names(beta))
log_msg(sprintf("   miRNAs used (%d): %s",
                length(mirnas_used), paste(mirnas_used, collapse = ", ")))

if (length(mirnas_used) == 0) {
  stop("no miRNA is shared between the mapping and the coefficients. Check names(feat$Feature).")
}

z_mat <- matrix(NA_real_, nrow = length(common_ids),
                ncol = length(mirnas_used),
                dimnames = list(common_ids, mirnas_used))

for (m in mirnas_used) {
  probe <- mirna_mapping[[m]]
  i <- which(expr$microRNA_ID == probe)
  if (length(i) != 1) {
    log_msg(sprintf("   probe '%s' is not unique (n=%d)", probe, length(i)))
    next
  }
  x_all <- as.numeric(expr[i, -1]); names(x_all) <- sample_ids
  x <- x_all[common_ids]
  x_log <- log2(x + 1)
  z_mat[, m] <- (x_log - sc$mirna$center[m]) / sc$mirna$scale[m]
}

log_msg(sprintf("   z_mat: NA count = %d of %d",
                sum(is.na(z_mat)), length(z_mat)))

# --- risk score ---
z_sub <- z_mat[, mirnas_used, drop = FALSE]
valid_rows <- complete.cases(z_sub)
log_msg(sprintf("   complete rows: %d of %d", sum(valid_rows), nrow(z_sub)))

if (sum(valid_rows) < 20) {
  stop(sprintf("only %d rows with a complete z — fewer than 20.", sum(valid_rows)))
}

risk <- rep(NA_real_, nrow(z_sub))
risk[valid_rows] <- as.numeric(z_sub[valid_rows, , drop = FALSE] %*%
                                 beta[mirnas_used])

# --- clinical ---
clin$OS     <- suppressWarnings(as.numeric(clin$OS))
clin$status <- suppressWarnings(as.integer(clin[["Censor (alive=0; dead=1)"]]))

# --- df_cgga ---
df_cgga <- data.frame(
  CGGA_ID    = common_ids,
  risk_score = risk,
  OS         = clin$OS[match(common_ids, clin$CGGA_ID)],
  status     = clin$status[match(common_ids, clin$CGGA_ID)],
  stringsAsFactors = FALSE
)

log_msg("")
log_msg("debug of df_cgga:")
log_msg(sprintf("   n = %d", nrow(df_cgga)))
log_msg(sprintf("   OS NA     : %d", sum(is.na(df_cgga$OS))))
log_msg(sprintf("   status NA : %d", sum(is.na(df_cgga$status))))
log_msg(sprintf("   risk NA   : %d", sum(is.na(df_cgga$risk_score))))
log_msg(sprintf("   OS > 0    : %d", sum(df_cgga$OS > 0, na.rm = TRUE)))

df_cgga <- df_cgga[!is.na(df_cgga$OS) & !is.na(df_cgga$status) &
                     df_cgga$OS > 0 & !is.na(df_cgga$risk_score), ]

log_msg(sprintf("   final n: %d | events = %d",
                nrow(df_cgga), sum(df_cgga$status == 1)))

if (nrow(df_cgga) < 20) {
  stop("not enough samples for KM. Look at the debug output.")
}

# --- KM ---
med <- median(df_cgga$risk_score, na.rm = TRUE)
df_cgga$risk_group <- factor(
  ifelse(df_cgga$risk_score >= med, "High risk", "Low risk"),
  levels = c("Low risk", "High risk"))

surv_obj <- Surv(df_cgga$OS, df_cgga$status)
fit <- survfit(surv_obj ~ risk_group, data = df_cgga)
lr  <- survdiff(surv_obj ~ risk_group, data = df_cgga)
lr_p <- 1 - pchisq(lr$chisq, df = 1)
log_msg(sprintf("   log-rank P = %.4g", lr_p))

p_label <- if (lr_p < 0.0001) "p < 0.0001" else sprintf("p = %.4g", lr_p)

km <- ggsurvplot(
  fit, data = df_cgga,
  pval = p_label, pval.method = TRUE,
  risk.table = TRUE, risk.table.col = "strata",
  palette = c("#2C7FB8", "#D7301F"),
  title = "CGGA External Cohort — Kaplan-Meier (median split)",
  legend.title = "Risk group",
  legend.labs = c("Low risk", "High risk"),
  xlab = "Time (days)", ylab = "Overall Survival Probability",
  ggtheme = theme_bw(base_size = 12),
  tables.theme = theme_cleantable(),
  break.time.by = 365, xlim = c(0, 3650))

ggsave(file.path(out_dir, "FigureS5_KM_CGGA.png"),
       plot = km$plot, width = 8, height = 6, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "FigureS5_KM_CGGA.pdf"),
       plot = km$plot, width = 8, height = 6)
write.csv(df_cgga,
          file.path(out_dir, "CGGA_sample_scores.csv"),
          row.names = FALSE)

log_msg("   FigureS5_KM_CGGA.png / .pdf")
