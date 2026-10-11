#!/usr/bin/env Rscript
# =============================================================================
# Supplementary Batch 3 — Table S15, S19, S20 + Figure S8 + Data S1-S3
# -----------------------------------------------------------------------------
# Table S15  — MDS gene list + PC1 loadings
# Table S19  — Age-only benchmark vs signature
# Table S20  — CGGA probe mapping
# Figure S8  — z-score distribution TCGA vs CGGA
# Data S1    — barcodes train/test
# Data S2    — risk scores test + CGGA
# Data S3    — session info
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(jsonlite)
})

input_dir <- paste0(WORK_DIR, "/")
sig_dir   <- file.path(input_dir, "signature_v2")
out_dir   <- file.path(input_dir, "supplementary_final")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(sprintf("[%s] %s\n",
                                     format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

norm_id <- function(x) substr(toupper(gsub("\\.", "-", as.character(x))), 1, 12)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 ||
                             (length(a) == 1 && is.na(a))) b else a

log_msg("============================================================")
log_msg("Supplementary Batch 3 — S15, S19, S20, S8, Data S1-S3")
log_msg("============================================================")


# =============================================================================
# TABLE S15 — MDS gene list + PC1 loadings
# =============================================================================
log_msg("")
log_msg("Table S15 — MDS gene list")

mds_param_path <- file.path(input_dir, "mds_transform_params.rds")
mito_path      <- file.path(input_dir, "mito_gene_list.csv")

if (!file.exists(mds_param_path)) {
  log_msg("   mds_transform_params.rds not found — reporting only mito_gene_list")
  mito <- fread(mito_path)
  write.csv(mito, file.path(out_dir, "TableS15_mito_genes_only.csv"),
            row.names = FALSE)
} else {
  mds_param <- readRDS(mds_param_path)

  # --- list of genes used ---
  gene_ids <- mds_param$mito_genes_used %||%
    mds_param$gene_ids %||%
    names(mds_param$gene_means_train) %||%
    NULL

  if (is.null(gene_ids)) {
    log_msg("   gene list not found in mds_transform_params.rds.")
  } else {
    # --- PC1 loadings (if saved) ---
    loadings <- mds_param$pca_rotation[, 1] %||%
      mds_param$pc1_loadings %||%
      NULL

    if (!is.null(loadings)) {
      tbl_s15 <- data.frame(
        gene_ensembl = gene_ids,
        PC1_loading  = round(as.numeric(loadings), 5),
        train_mean   = round(as.numeric(mds_param$gene_means_train), 4),
        train_sd     = round(as.numeric(mds_param$gene_sds_train), 4),
        stringsAsFactors = FALSE)
    } else {
      tbl_s15 <- data.frame(
        gene_ensembl = gene_ids,
        train_mean   = round(as.numeric(mds_param$gene_means_train), 4),
        train_sd     = round(as.numeric(mds_param$gene_sds_train), 4),
        stringsAsFactors = FALSE)
    }

    tbl_s15 <- tbl_s15[order(-abs(tbl_s15$PC1_loading %||%
                                    tbl_s15$train_mean)), ]

    write.csv(tbl_s15,
              file.path(out_dir, "TableS15_MDS_genes_loadings.csv"),
              row.names = FALSE)

    log_msg(sprintf("   TableS15 saved (%d genes)", nrow(tbl_s15)))
  }
}


# =============================================================================
# TABLE S19 — Age-only benchmark
# =============================================================================
log_msg("")
log_msg("Table S19 — Age-only vs signature")

pred_path <- file.path(sig_dir, "signature_v2_test_predictions.csv")
pred <- fread(pred_path)

y_test <- Surv(pred$survival_time, pred$event)

# --- Model 1: age only ---
cox_age <- coxph(y_test ~ scale(pred$age_years))

# --- Model 2: risk score only ---
cox_risk <- coxph(y_test ~ scale(pred$risk_score_full))

# --- Model 3: age + risk score ---
cox_both <- coxph(y_test ~ scale(pred$age_years) + scale(pred$risk_score_full))

# --- C-index ---
c_age  <- summary(cox_age)$concordance[1]
c_risk <- summary(cox_risk)$concordance[1]
c_both <- summary(cox_both)$concordance[1]

# --- LRT: age vs age + risk ---
lrt <- anova(cox_age, cox_both, test = "LRT")
lrt_p     <- lrt[2, "Pr(>|Chi|)"]
lrt_chisq <- lrt$Chisq[2]
lrt_df    <- lrt$Df[2]

# --- HR risk (adjusted for age) ---
hr_risk_adj <- summary(cox_both)$coefficients["scale(pred$risk_score_full)", "exp(coef)"]
hr_lo       <- exp(summary(cox_both)$coefficients["scale(pred$risk_score_full)", "coef"] -
                     1.96 * summary(cox_both)$coefficients["scale(pred$risk_score_full)", "se(coef)"])
hr_hi       <- exp(summary(cox_both)$coefficients["scale(pred$risk_score_full)", "coef"] +
                     1.96 * summary(cox_both)$coefficients["scale(pred$risk_score_full)", "se(coef)"])
p_risk_adj  <- summary(cox_both)$coefficients["scale(pred$risk_score_full)", "Pr(>|z|)"]

tbl_s19 <- data.frame(
  Model          = c("Age only", "Risk score only", "Age + Risk score"),
  C_index        = round(c(c_age, c_risk, c_both), 4),
  AIC            = round(c(AIC(cox_age), AIC(cox_risk), AIC(cox_both)), 2),
  stringsAsFactors = FALSE)

tbl_s19_extra <- data.frame(
  Metric          = c("Risk-score HR (adjusted for age)",
                      "95% CI (adjusted)",
                      "P-value (adjusted)",
                      "LRT (age-only vs age + risk)",
                      "ΔAIC (age-only → age + risk)"),
  Value           = c(round(hr_risk_adj, 3),
                      sprintf("%.3f\u2013%.3f", hr_lo, hr_hi),
                      sprintf("%.4f", p_risk_adj),
                      sprintf("chi-sq=%.3f, df=%d, P=%.4f", lrt_chisq, lrt_df, lrt_p),
                      round(AIC(cox_both) - AIC(cox_age), 2)),
  stringsAsFactors = FALSE)

write.csv(tbl_s19,
          file.path(out_dir, "TableS19_age_benchmark.csv"),
          row.names = FALSE)
write.csv(tbl_s19_extra,
          file.path(out_dir, "TableS19_extra.csv"),
          row.names = FALSE)

log_msg("   TableS19 saved")
print(tbl_s19, row.names = FALSE)


# =============================================================================
# TABLE S20 — CGGA probe mapping
# =============================================================================
log_msg("")
log_msg("Table S20 — CGGA probe mapping")

# --- check the probes available in CGGA ---
expr_path <- file.path(input_dir,
                       "cgga_raw/exp cgga/CGGA.microRNA_array_198_gene_level.20200506.txt")
expr_cgga <- read.delim(expr_path, header = TRUE, sep = "\t",
                        check.names = FALSE, stringsAsFactors = FALSE)

probe_map <- data.frame(
  TCGA_miRNA = c("hsa-miR-196a-5p", "hsa-miR-675-5p",
                 "hsa-miR-221-3p", "hsa-miR-31-5p",
                 "hsa-miR-181a-2-3p", "hsa-miR-4660-5p"),
  CGGA_probe = c("hsa-miR-196a", "hsa-miR-675",
                 "hsa-miR-221#", "hsa-miR-31",
                 "hsa-miR-181a-2#", NA),
  stringsAsFactors = FALSE)

probe_map$In_CGGA <- sapply(probe_map$CGGA_probe, function(p) {
  if (is.na(p)) return("No")
  if (p %in% expr_cgga$microRNA_ID) "Yes" else "No"
})

probe_map$Note <- c(
  "Direct match (5p arm)",
  "Direct match (5p arm)",
  "# indicates 3p arm in CGGA annotation",
  "Direct match (5p arm)",
  "# indicates 3p arm in CGGA annotation",
  "Not represented on microRNA_array_198 platform"
)

write.csv(probe_map,
          file.path(out_dir, "TableS20_CGGA_probe_mapping.csv"),
          row.names = FALSE)

log_msg("   TableS20 saved")
print(probe_map, row.names = FALSE)


# =============================================================================
# FIGURE S8 — z-score distribution TCGA vs CGGA
# =============================================================================
log_msg("")
log_msg("Figure S8 — z-score TCGA vs CGGA")

# --- CGGA z-scores ---
sc_file <- file.path(sig_dir, "signature_v2_scaling.rds")
if (!file.exists(sc_file)) sc_file <- file.path(input_dir, "signature_scaling_params.rds")
sc <- readRDS(sc_file)

mirna_mapping <- c(
  "hsa-miR-196a-5p"   = "hsa-miR-196a",
  "hsa-miR-675-5p"    = "hsa-miR-675",
  "hsa-miR-221-3p"    = "hsa-miR-221#",
  "hsa-miR-31-5p"     = "hsa-miR-31",
  "hsa-miR-181a-2-3p" = "hsa-miR-181a-2#"
)

sample_ids <- as.character(colnames(expr_cgga)[-1])

z_cgga_list <- list()
for (m in names(mirna_mapping)) {
  probe <- mirna_mapping[[m]]
  i <- which(expr_cgga$microRNA_ID == probe)
  if (length(i) != 1) next
  x <- as.numeric(expr_cgga[i, -1])
  x_log <- log2(x + 1)
  z <- (x_log - sc$mirna$center[m]) / sc$mirna$scale[m]
  z_cgga_list[[m]] <- data.frame(miRNA = m, z = z, cohort = "CGGA")
}
z_cgga <- rbindlist(z_cgga_list)

# --- TCGA z-scores ---
mirna_train <- read.csv(file.path(input_dir, "mirna_isoform_rpm_matrix.csv"),
                        row.names = 1, check.names = FALSE)
if (mean(grepl("^TCGA-", rownames(mirna_train))) > 0.5) {
  mirna_train <- t(as.matrix(mirna_train))
} else {
  mirna_train <- as.matrix(mirna_train)
}
# restrict to the shared miRNAs
common_mirnas <- intersect(names(mirna_mapping), rownames(mirna_train))
z_tcga_list <- list()
for (m in common_mirnas) {
  x <- as.numeric(mirna_train[m, ])
  x_log <- log2(x + 1)
  z <- (x_log - sc$mirna$center[m]) / sc$mirna$scale[m]
  z_tcga_list[[m]] <- data.frame(miRNA = m, z = z, cohort = "TCGA")
}
z_tcga <- rbindlist(z_tcga_list)

z_all <- rbindlist(list(z_tcga, z_cgga))

# --- plot ---
p_s8 <- ggplot(z_all, aes(x = z, fill = cohort, color = cohort)) +
  geom_density(alpha = 0.30, linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ miRNA, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("TCGA" = "#2166AC", "CGGA" = "#B2182B")) +
  scale_color_manual(values = c("TCGA" = "#2166AC", "CGGA" = "#B2182B")) +
  labs(
    title = "Z-score distribution of signature miRNAs: TCGA vs CGGA",
    subtitle = "TCGA-derived centering/scaling applied identically to both cohorts",
    x = "Standardized expression (z-score)",
    y = "Density") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom")

ggsave(file.path(out_dir, "FigureS8_zscore_comparison.png"),
       p_s8, width = 10, height = 5, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "FigureS8_zscore_comparison.pdf"),
       p_s8, width = 10, height = 5)

writeLines(c(
  "Figure S8. Density distributions of standardized expression (z-scores) ",
  "for the five signature miRNAs with available CGGA probes, comparing TCGA ",
  "and CGGA cohorts. Both cohorts were transformed using the identical ",
  "log2(x + 1) operation followed by direct application of the TCGA-locked ",
  "centering and scaling parameters. Differences in the shape and location ",
  "of the z-score distributions between cohorts reflect platform-level ",
  "differences in absolute expression range and should be considered when ",
  "interpreting cross-cohort performance."
), file.path(out_dir, "FigureS8_caption.txt"))

log_msg("   FigureS8 saved")


# =============================================================================
# DATA S1 — Barcodes
# =============================================================================
log_msg("")
log_msg("Data S1 — Barcodes")

train_df <- fread(file.path(input_dir, "train_set.csv"))
test_df  <- fread(file.path(input_dir, "test_set_LOCKED.csv"))

barcodes <- rbind(
  data.frame(patient_barcode = train_df$patient_barcode, split = "train",
             seed = 42, stringsAsFactors = FALSE),
  data.frame(patient_barcode = test_df$patient_barcode, split = "test",
             seed = 42, stringsAsFactors = FALSE))

write.csv(barcodes,
          file.path(out_dir, "DataS1_barcodes.csv"),
          row.names = FALSE)

log_msg(sprintf("   DataS1 saved (train=%d, test=%d)",
                nrow(train_df), nrow(test_df)))


# =============================================================================
# DATA S2 — Risk scores
# =============================================================================
log_msg("")
log_msg("Data S2 — Risk scores")

pred <- fread(file.path(sig_dir, "signature_v2_test_predictions.csv"))

data_s2 <- data.frame(
  patient_id    = pred$patient_id,
  cohort        = "TCGA_Test",
  risk_score    = pred$risk_score_full,
  risk_score_miRNA_only = pred$risk_score_mirna,
  survival_time = pred$survival_time,
  event         = pred$event,
  stringsAsFactors = FALSE)

# if CGGA sample scores exist, add them
cgga_path <- file.path(out_dir, "CGGA_sample_scores.csv")
if (file.exists(cgga_path)) {
  cgga <- fread(cgga_path)
  cgga_data <- data.frame(
    patient_id    = cgga$CGGA_ID,
    cohort        = "CGGA",
    risk_score    = cgga$risk_score,
    risk_score_miRNA_only = NA_real_,
    survival_time = cgga$OS,
    event         = cgga$status,
    stringsAsFactors = FALSE)
  data_s2 <- rbind(data_s2, cgga_data)
}

write.csv(data_s2,
          file.path(out_dir, "DataS2_risk_scores.csv"),
          row.names = FALSE)

log_msg(sprintf("   DataS2 saved (%d rows)", nrow(data_s2)))


# =============================================================================
# DATA S3 — Session info
# =============================================================================
log_msg("")
log_msg("Data S3 — Session info")

si <- capture.output(sessionInfo())

writeLines(c(
  "============================================================",
  "Session Information — miR-Mito Score (Signature v2)",
  "============================================================",
  "",
  sprintf("Date: %s", as.character(Sys.time())),
  sprintf("R version: %s", R.version.string),
  "",
  "------------------- Package versions -------------------",
  "",
  si
), file.path(out_dir, "DataS3_sessionInfo.txt"))

# --- versions of the key packages in a CSV ---
key_pkgs <- c("survival", "glmnet", "DESeq2", "timeROC",
              "survminer", "ggplot2", "dplyr", "data.table",
              "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi",
              "jsonlite", "digest")
pkg_versions <- data.frame(
  Package = key_pkgs,
  Version = sapply(key_pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE))
      as.character(packageVersion(p))
    else "not_installed"
  }),
  stringsAsFactors = FALSE)

write.csv(pkg_versions,
          file.path(out_dir, "DataS3_package_versions.csv"),
          row.names = FALSE)

log_msg("   DataS3 saved")
print(pkg_versions, row.names = FALSE)


# =============================================================================
# summary
# =============================================================================
log_msg("")
log_msg("============================================================")
log_msg("all outputs in: ", out_dir)
log_msg("   TableS15_MDS_genes_loadings.csv  (or _mito_genes_only)")
log_msg("   TableS19_age_benchmark.csv + TableS19_extra.csv")
log_msg("   TableS20_CGGA_probe_mapping.csv")
log_msg("   FigureS8_zscore_comparison.png / .pdf")
log_msg("   DataS1_barcodes.csv")
log_msg("   DataS2_risk_scores.csv")
log_msg("   DataS3_sessionInfo.txt + DataS3_package_versions.csv")
log_msg("============================================================")
