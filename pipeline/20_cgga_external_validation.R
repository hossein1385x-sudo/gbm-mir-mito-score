#!/usr/bin/env Rscript
# =============================================================================
# CGGA External Validation — New (non-leaky) Signatures
# -----------------------------------------------------------------------------
# Validation of the new (non-leaky) signatures
#
# Two signatures:
#   Main:  196a-5p + 675-5p + Age
#   v2:    196a-5p + 675-5p + 221-3p + 31-5p + 181a-2-3p + Age
#          (4660-5p is not present in CGGA -> excluded from validation)
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(survival)
  library(jsonlite)
})

output_dir <- paste0(WORK_DIR, "/")

expr_file <- file.path(output_dir,
                       "cgga_raw/exp cgga/CGGA.microRNA_array_198_gene_level.20200506.txt")
clin_file <- file.path(output_dir,
                       "cgga_raw/clin cgga/CGGA.microRNA_array_198_clinical.20200506.txt")

sig_main_path <- file.path(output_dir, "signature_locked_latest.rds")
sig_v2_path   <- file.path(output_dir, "signature_v2/signature_v2_locked_latest.rds")

log_msg <- function(...) cat(sprintf("[%s] %s\n",
                                     format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

log_msg("============================================================")
log_msg("CGGA External Validation — New Signatures")
log_msg("============================================================")

# ======================================================================
# 1) Load CGGA
# ======================================================================
expr <- read.delim(expr_file, header = TRUE, sep = "\t",
                   check.names = FALSE, stringsAsFactors = FALSE)
clin <- read.delim(clin_file, header = TRUE, sep = "\t",
                   check.names = FALSE, stringsAsFactors = FALSE)

sample_ids <- colnames(expr)[-1]
all_mirnas <- expr$microRNA_ID

log_msg(sprintf("   CGGA: %d miRNAs × %d samples", length(all_mirnas), length(sample_ids)))

# ======================================================================
# 2) Explicit TCGA → CGGA mapping (manual decision based on the array structure)
# ======================================================================
# Rule:
#   - if the suffix is -5p → in CGGA usually without #
#   - if the suffix is -3p → in CGGA usually with #
#   - if both exist, prefer the one that matches the rule
mirna_mapping <- list(
  "hsa-miR-196a-5p"    = "hsa-miR-196a",
  "hsa-miR-675-5p"     = "hsa-miR-675",
  "hsa-miR-221-3p"     = "hsa-miR-221#",
  "hsa-miR-31-5p"      = "hsa-miR-31",
  "hsa-miR-181a-2-3p"  = "hsa-miR-181a-2#"
)

# verify the mapping
for (m in names(mirna_mapping)) {
  probe <- mirna_mapping[[m]]
  ok <- probe %in% all_mirnas
  log_msg(sprintf("   %-22s → %-22s %s", m, probe,
                  ifelse(ok, "OK", "NOT FOUND")))
  if (!ok) stop("probe not found: ", probe)
}

# ======================================================================
# 3) Extraction and z-score
# ======================================================================
extract_z <- function(probe, tcga_name, sig) {
  i <- which(all_mirnas == probe)
  if (length(i) != 1) stop("multiple or zero probes: ", probe)

  x_raw <- as.numeric(expr[i, -1])
  names(x_raw) <- sample_ids

  x_log <- log2(x_raw + 1)

  center_m <- sig$scaling$mirna$center[tcga_name]
  scale_m  <- sig$scaling$mirna$scale[tcga_name]

  if (is.na(center_m) || is.na(scale_m) || scale_m == 0)
    stop("center/scale missing for: ", tcga_name)

  (x_log - center_m) / scale_m
}

# ======================================================================
# 4) Prepare the clinical data
# ======================================================================
clin$CGGA_ID <- as.character(clin$CGGA_ID)
clin$OS      <- as.numeric(clin$OS)
clin$status  <- as.numeric(clin[["Censor (alive=0; dead=1)"]])

# Age
if ("age_at_diagnosis" %in% names(clin)) {
  age_raw <- as.numeric(clin$age_at_diagnosis)
} else if ("Age" %in% names(clin)) {
  age_raw <- as.numeric(clin$Age)
} else {
  age_raw <- NA_real_
}

# merge
idx <- match(sample_ids, clin$CGGA_ID)
clin_df <- data.frame(
  CGGA_ID = sample_ids,
  OS      = clin$OS[idx],
  status  = clin$status[idx],
  age_raw = age_raw[idx],
  stringsAsFactors = FALSE
)

# ======================================================================
# 5) Validation function
# ======================================================================
validate_signature <- function(sig, label, available_mirnas) {
  log_msg("")
  log_msg(sprintf("[%s] validating with %d miRNAs",
                  label, length(available_mirnas)))

  # build the z matrix
  z_mat <- matrix(NA_real_, nrow = length(sample_ids),
                  ncol = length(available_mirnas))
  colnames(z_mat) <- available_mirnas
  rownames(z_mat) <- sample_ids

  for (m in available_mirnas) {
    z_mat[, m] <- extract_z(mirna_mapping[[m]], m, sig)
  }

  # coefficients
  coefs <- setNames(sig$coefficients$Coefficient, sig$coefficients$Feature)
  beta_m <- coefs[available_mirnas]

  # risk = Σ β·z  (miRNA)
  risk_mirna <- as.numeric(z_mat %*% beta_m)

  # Age, if it is in the signature
  has_age <- "Age" %in% sig$clinical_features
  if (has_age) {
    age_years <- if (median(clin_df$age_raw, na.rm = TRUE) > 200)
      clin_df$age_raw / 365.25 else clin_df$age_raw
    age_years[is.na(age_years)] <- median(age_years, na.rm = TRUE)

    age_z <- (age_years - sig$scaling$clinical$center["Age"]) /
      sig$scaling$clinical$scale["Age"]
    risk <- risk_mirna + coefs["Age"] * age_z
    risk_mirna_only <- risk_mirna
  } else {
    risk <- risk_mirna
    risk_mirna_only <- risk_mirna
  }

  # final dataset
  df <- data.frame(
    CGGA_ID = clin_df$CGGA_ID,
    risk_full = risk,
    risk_mirna_only = risk_mirna_only,
    OS = clin_df$OS,
    status = clin_df$status,
    stringsAsFactors = FALSE
  )
  df <- df[complete.cases(df[, c("OS", "status", "risk_full")]), ]
  df <- df[df$OS > 0, ]

  if (nrow(df) < 20) {
    log_msg(sprintf("   n=%d is too small", nrow(df)))
    return(NULL)
  }

  y <- Surv(df$OS, df$status)

  # --- Full score ---
  c_full <- survival::survConcordance(y ~ df$risk_full)
  cox_full <- coxph(y ~ scale(df$risk_full))
  s_full <- summary(cox_full)

  # --- miRNA-only ---
  c_mir <- survival::survConcordance(y ~ df$risk_mirna_only)
  cox_mir <- coxph(y ~ scale(df$risk_mirna_only))
  s_mir <- summary(cox_mir)

  # --- KM / logrank median split ---
  med <- median(df$risk_full)
  grp <- factor(ifelse(df$risk_full >= med, "High", "Low"),
                levels = c("Low", "High"))
  lr <- survdiff(y ~ grp)
  lr_p <- 1 - pchisq(lr$chisq, df = 1)

  data.frame(
    signature            = label,
    n                    = nrow(df),
    events               = sum(df$status == 1),
    n_mirnas_available   = length(available_mirnas),
    mirnas_used          = paste(available_mirnas, collapse = "; "),
    mirnas_missing       = paste(setdiff(sig$mirna_features, available_mirnas),
                                 collapse = "; "),

    c_full               = as.numeric(c_full$concordance),
    c_full_SE            = as.numeric(c_full$std.err),
    c_full_CI_low        = max(0, c_full$concordance - 1.96 * c_full$std.err),
    c_full_CI_high       = min(1, c_full$concordance + 1.96 * c_full$std.err),

    HR_full_per_SD       = exp(s_full$coefficients[1, "coef"]),
    HR_full_CI_low       = exp(s_full$coefficients[1, "coef"] -
                                 1.96 * s_full$coefficients[1, "se(coef)"]),
    HR_full_CI_high      = exp(s_full$coefficients[1, "coef"] +
                                 1.96 * s_full$coefficients[1, "se(coef)"]),
    p_full               = s_full$coefficients[1, "Pr(>|z|)"],

    c_mirna_only         = as.numeric(c_mir$concordance),
    HR_mirna_only_per_SD = exp(s_mir$coefficients[1, "coef"]),
    HR_mirna_only_CI_low = exp(s_mir$coefficients[1, "coef"] -
                                 1.96 * s_mir$coefficients[1, "se(coef)"]),
    HR_mirna_only_CI_high= exp(s_mir$coefficients[1, "coef"] +
                                 1.96 * s_mir$coefficients[1, "se(coef)"]),
    p_mirna_only         = s_mir$coefficients[1, "Pr(>|z|)"],

    logrank_median_p     = lr_p,

    stringsAsFactors = FALSE
  )
}

# ======================================================================
# 6) Load the signatures and run the validation
# ======================================================================
log_msg("")
log_msg("loading the signatures ...")

sig_main <- if (file.exists(sig_main_path)) readRDS(sig_main_path) else NULL
sig_v2   <- if (file.exists(sig_v2_path))   readRDS(sig_v2_path)   else NULL

if (!is.null(sig_main)) {
  log_msg(sprintf("   Main: %d miRNA — %s",
                  length(sig_main$mirna_features),
                  paste(sig_main$mirna_features, collapse = ", ")))
}
if (!is.null(sig_v2)) {
  log_msg(sprintf("   v2:   %d miRNA — %s",
                  length(sig_v2$mirna_features),
                  paste(sig_v2$mirna_features, collapse = ", ")))
}

# miRNAs available in CGGA for each signature
cgga_available <- names(mirna_mapping)

res_main <- if (!is.null(sig_main)) {
  avail_main <- intersect(sig_main$mirna_features, cgga_available)
  validate_signature(sig_main, "Main", avail_main)
} else NULL

res_v2 <- if (!is.null(sig_v2)) {
  avail_v2 <- intersect(sig_v2$mirna_features, cgga_available)
  validate_signature(sig_v2, "Sensitivity_v2", avail_v2)
} else NULL

all_res <- rbind(res_main, res_v2)
print(all_res)

# ======================================================================
# 7) Save
# ======================================================================
write.csv(all_res,
          file.path(output_dir, "CGGA_new_signature_validation.csv"),
          row.names = FALSE)

manifest <- list(
  validation_type     = "external_CGGA",
  main_signature      = if (!is.null(sig_main)) sig_main$mirna_features else NULL,
  v2_signature        = if (!is.null(sig_v2))   sig_v2$mirna_features   else NULL,
  cgga_available      = cgga_available,
  mirna_mapping       = mirna_mapping,
  transformation      = "log2(x+1) + TCGA locked center/scale",
  n_cgga_total        = length(sample_ids),
  timestamp           = as.character(Sys.time())
)
write_json(manifest,
           file.path(output_dir, "CGGA_new_signature_manifest.json"),
           pretty = TRUE, auto_unbox = TRUE)

log_msg("")
log_msg("============================================================")
log_msg("summary")
log_msg("============================================================")
for (i in seq_len(nrow(all_res))) {
  r <- all_res[i, ]
  log_msg(sprintf("[%s] n=%d events=%d | C_full=%.3f (%.3f–%.3f) | HR_full=%.3f (p=%.4g)",
                  r$signature, r$n, r$events,
                  r$c_full, r$c_full_CI_low, r$c_full_CI_high,
                  r$HR_full_per_SD, r$p_full))
  log_msg(sprintf("        mirnas_missing: %s",
                  ifelse(r$mirnas_missing == "", "none", r$mirnas_missing)))
}
log_msg("")
log_msg("outputs:")
log_msg("   CGGA_new_signature_validation.csv")
log_msg("   CGGA_new_signature_manifest.json")
log_msg("DONE")
