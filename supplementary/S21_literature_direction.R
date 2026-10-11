#!/usr/bin/env Rscript
# =============================================================================
# Supplementary — Table S21 (Literature direction-of-effect concordance)
# -----------------------------------------------------------------------------
# Build the pre-filled table for the 6 miRNAs of the v2 signature
# Outputs in supplementary_final/
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table)
})

input_dir <- paste0(WORK_DIR, "/")
sig_dir   <- file.path(input_dir, "signature_v2")
out_dir   <- file.path(input_dir, "supplementary_final")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(sprintf("[%s] %s\n",
                                     format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

log_msg("============================================================")
log_msg("Table S21 — Literature direction concordance")
log_msg("============================================================")

# ======================================================================
# 1) Coefficients of the signature miRNAs
# ======================================================================
feat_path <- file.path(sig_dir, "signature_features.csv")
if (!file.exists(feat_path)) stop("signature_features.csv not found.")

feat <- fread(feat_path)
feat_mirna <- feat[grepl("^hsa-", Feature)]
feat_mirna[, model_direction := ifelse(Coefficient > 0, "Risk", "Protective")]
feat_mirna[, model_beta      := round(Coefficient, 4)]

# ======================================================================
# 2) Literature information (manual, from PubMed)
# ======================================================================
literature <- data.table(
  Feature = c(
    "hsa-miR-196a-5p",
    "hsa-miR-675-5p",
    "hsa-miR-221-3p",
    "hsa-miR-31-5p",
    "hsa-miR-181a-2-3p",
    "hsa-miR-4660-5p"
  ),
  Representative_references = c(
    "PMID 33469841; PMID 33289788",
    "PMID 28970140",
    "PMID 42695418",
    "PMID 32906592; PMID 31933922",
    "PMID 41721414; PMID 40093871; PMID 36232448",
    "—"
  ),
  Cancer_type = c(
    "Glioma / GBM",
    "Glioma",
    "Glioma",
    "Glioma",
    "GBM / glioma (5p arm only)",
    "—"
  ),
  Reported_direction = c(
    "Upregulated, oncogenic",
    "Upregulated, oncogenic",
    "Upregulated, oncogenic",
    "Tumor-suppressive (downregulated in tumors)",
    "Studies address miR-181a-5p arm; 2-3p–specific evidence absent",
    "No GBM-specific literature"
  ),
  Concordance = c(
    "Concordant",
    "Concordant",
    "Concordant",
    "Discordant",
    "Insufficient (arm mismatch)",
    "Insufficient (4xxx family, recently annotated)"
  ),
  Notes = c(
    "Independent reports link miR-196a upregulation to malignant glioma phenotype and hypoxia response.",
    "miR-675 promotes glioma proliferation/motility via RB1 repression.",
    "miR-221/222 promote glioma progression via Akt pathway activation.",
    "Published evidence supports a tumor-suppressive role for miR-31 in glioma, opposite to the positive model coefficient.",
    "Literature on miR-181a in glioma pertains to the 5p arm; the seed sequence and targetome of 2-3p differ substantially.",
    "No GBM-specific functional or prognostic literature identified within the scope of the present search."
  )
)

# ======================================================================
# 3) Merge with the model coefficients
# ======================================================================
tbl_s21 <- merge(
  feat_mirna[, .(miRNA = Feature,
                 Model_direction = model_direction,
                 Model_beta = model_beta)],
  literature,
  by.x = "miRNA", by.y = "Feature",
  all.x = TRUE)

# order by |beta|
tbl_s21 <- tbl_s21[order(-abs(Model_beta))]

# order the columns
setcolorder(tbl_s21, c(
  "miRNA",
  "Model_direction",
  "Model_beta",
  "Representative_references",
  "Cancer_type",
  "Reported_direction",
  "Concordance",
  "Notes"
))

# ======================================================================
# 4) Save
# ======================================================================
out_path <- file.path(out_dir, "TableS21_literature_direction.csv")
write.csv(tbl_s21, out_path, row.names = FALSE)

# short version (without Notes) for inclusion in the manuscript
tbl_s21_short <- tbl_s21[, .(
  miRNA,
  `Model direction` = Model_direction,
  `β (model)` = Model_beta,
  `Representative reference(s)` = Representative_references,
  `Reported direction in glioma` = Reported_direction,
  Concordance
)]
write.csv(tbl_s21_short,
          file.path(out_dir, "TableS21_literature_direction_short.csv"),
          row.names = FALSE)

# ======================================================================
# 5) Count summary
# ======================================================================
counts <- tbl_s21[, .N, by = Concordance][order(-N)]
write.csv(counts,
          file.path(out_dir, "TableS21_concordance_summary.csv"),
          row.names = FALSE)

# ======================================================================
# 6) Console output
# ======================================================================
log_msg("")
log_msg("Table S21 created:")
print(tbl_s21[, .(miRNA, Model_direction, Model_beta, Concordance)],
      row.names = FALSE)

log_msg("")
log_msg("Concordance counts:")
print(counts, row.names = FALSE)

log_msg("")
log_msg("outputs:")
log_msg("   TableS21_literature_direction.csv          (full, with Notes)")
log_msg("   TableS21_literature_direction_short.csv    (short version)")
log_msg("   TableS21_concordance_summary.csv           (counts)")
log_msg("")
log_msg("DONE")
