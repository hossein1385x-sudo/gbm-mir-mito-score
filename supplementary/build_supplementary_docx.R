#!/usr/bin/env Rscript
# =============================================================================
# Build Supplementary_Materials.docx - version compatible with officer 0.7.6
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(officer)
  library(flextable)
  library(readr)
  library(dplyr)
})

# --- page-break helper (compatible with different officer versions) ---
add_page_break <- function(doc) {
  if (exists("body_add_pagebreak", mode = "function")) {
    doc <- body_add_pagebreak(doc)
  } else if (exists("body_add_break", mode = "function")) {
    doc <- body_add_break(doc, pos = "after")
  }
  doc
}

supp_dir <- file.path(WORK_DIR, "supplementary_final")
out_file <- file.path(supp_dir, "Supplementary_Materials.docx")

log_msg <- function(...) cat(sprintf("[%s] %s\n",
                                     format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

log_msg("============================================================")
log_msg("building Supplementary_Materials.docx")
log_msg("============================================================")

# --- styles ---
tryCatch({
  set_flextable_defaults(
    font.family  = "Calibri",
    font.size    = 9,
    border.color = "grey70",
    padding      = 4
  )
}, error = function(e) NULL)

# --- helper functions ---
csv_to_flextable <- function(csv_path, max_rows = NULL, font_size = 9) {
  df <- read_csv(csv_path, show_col_types = FALSE)
  note <- NULL
  if (!is.null(max_rows) && nrow(df) > max_rows) {
    note <- sprintf("(Showing first %d of %d rows)", max_rows, nrow(df))
    df <- head(df, max_rows)
  }
  ft <- flextable(df) %>%
    fontsize(size = font_size, part = "all") %>%
    font(fontname = "Calibri", part = "all") %>%
    bold(part = "header") %>%
    bg(bg = "#EDF3F7", part = "header") %>%
    align(align = "left", part = "all") %>%
    valign(valign = "top", part = "all") %>%
    autofit()
  list(ft = ft, note = note, n_rows = nrow(df))
}

add_table_with_caption <- function(doc, title, csv_path, max_rows = NULL) {
  if (!file.exists(csv_path)) {
    return(doc %>%
             body_add_par(title, style = "heading 2") %>%
             body_add_par("[File not found — skipped]", style = "Normal"))
  }
  info <- csv_to_flextable(csv_path, max_rows = max_rows)
  doc <- doc %>%
    body_add_par(title, style = "heading 2") %>%
    body_add_flextable(info$ft)
  if (!is.null(info$note)) {
    doc <- doc %>% body_add_par(info$note, style = "Normal")
  }
  doc %>% body_add_par("", style = "Normal")
}

add_figure <- function(doc, title, png_path, caption = NULL, width_in = 6.5) {
  if (!file.exists(png_path)) {
    return(doc %>%
             body_add_par(title, style = "heading 2") %>%
             body_add_par("[Image not found — skipped]", style = "Normal"))
  }
  doc <- doc %>%
    body_add_par(title, style = "heading 2") %>%
    body_add_img(src = png_path, width = width_in, height = width_in * 0.7)
  if (!is.null(caption)) {
    doc <- doc %>% body_add_par(caption, style = "Normal")
  }
  doc %>% body_add_par("", style = "Normal")
}

# =============================================================================
# Build the DOCX
# =============================================================================
doc <- read_docx()

# ---- title page ----
doc <- doc %>%
  body_add_par("Supplementary Materials", style = "heading 1") %>%
  body_add_par("A Mitochondrial Dysfunction-Associated miRNA Signature for Prognostic Stratification of Glioblastoma",
               style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_par(sprintf("Generated: %s", as.character(Sys.time())), style = "Normal") %>%
  body_add_par("", style = "Normal")

doc <- add_page_break(doc)

# ---- table of contents ----
doc <- doc %>%
  body_add_par("Contents", style = "heading 1") %>%
  body_add_par("Tables", style = "heading 2") %>%
  body_add_par("Table S1   — Participant flow", style = "Normal") %>%
  body_add_par("Table S12  — miRNA–target concordance", style = "Normal") %>%
  body_add_par("Table S13  — GO/KEGG enrichment", style = "Normal") %>%
  body_add_par("Table S14  — Locked model specification", style = "Normal") %>%
  body_add_par("Table S15  — MDS gene list and PC1 loadings", style = "Normal") %>%
  body_add_par("Table S18  — Baseline characteristics", style = "Normal") %>%
  body_add_par("Table S19  — Age-only benchmark", style = "Normal") %>%
  body_add_par("Table S20  — CGGA probe mapping", style = "Normal") %>%
  body_add_par("Table S21  — Literature direction concordance", style = "Normal") %>%
  body_add_par("Table S22  — MDS tertile sensitivity", style = "Normal") %>%
  body_add_par("Figures", style = "heading 2") %>%
  body_add_par("Figure S4 — miRNA–target network (Cytoscape)", style = "Normal") %>%
  body_add_par("Figure S5 — Kaplan-Meier (CGGA)", style = "Normal") %>%
  body_add_par("Figure S6 — Schoenfeld residuals", style = "Normal") %>%
  body_add_par("Figure S7 — Spline hazard ratio", style = "Normal") %>%
  body_add_par("Figure S8 — Z-score TCGA vs CGGA", style = "Normal") %>%
  body_add_par("Figure S9 — Enrichment dot plot", style = "Normal") %>%
  body_add_par("Data and Methods", style = "heading 2") %>%
  body_add_par("Data S1–S3, Supplementary Methods, TRIPOD+AI Checklist",
               style = "Normal")

doc <- add_page_break(doc)

# ---- tables ----
doc <- doc %>% body_add_par("Tables", style = "heading 1")

doc <- add_table_with_caption(doc, "Table S1. Participant flow",
                              file.path(supp_dir, "TableS1_participant_flow.csv"))

doc <- add_table_with_caption(doc, "Table S12. miRNA–target concordance (per-miRNA)",
                              file.path(supp_dir, "TableS12_concordance_per_mirna.csv"))

doc <- add_table_with_caption(doc, "Table S13. GO/KEGG enrichment",
                              file.path(supp_dir, "TableS13_enrichment.csv"),
                              max_rows = 40)

doc <- add_table_with_caption(doc, "Table S14a. Locked signature coefficients",
                              file.path(supp_dir, "TableS14a_coefficients.csv"))
doc <- add_table_with_caption(doc, "Table S14b. Training mean and SD (log2)",
                              file.path(supp_dir, "TableS14b_mean_sd_log2.csv"))
doc <- add_table_with_caption(doc, "Table S14d. Model manifest",
                              file.path(supp_dir, "TableS14d_manifest.csv"))

doc <- add_table_with_caption(doc, "Table S15. MDS gene list and PC1 loadings",
                              file.path(supp_dir, "TableS15_MDS_genes_loadings.csv"),
                              max_rows = 100)

doc <- add_table_with_caption(doc, "Table S18. Baseline characteristics",
                              file.path(supp_dir, "TableS18_baseline_characteristics.csv"))

doc <- add_table_with_caption(doc, "Table S19. Age-only benchmark",
                              file.path(supp_dir, "TableS19_age_benchmark.csv"))
doc <- add_table_with_caption(doc, "Table S19 (extra). Adjusted HR and LRT",
                              file.path(supp_dir, "TableS19_extra.csv"))

doc <- add_table_with_caption(doc, "Table S20. CGGA probe mapping",
                              file.path(supp_dir, "TableS20_CGGA_probe_mapping.csv"))

doc <- add_table_with_caption(doc, "Table S21. Literature direction concordance",
                              file.path(supp_dir, "TableS21_literature_direction.csv"))

doc <- add_table_with_caption(doc, "Table S22. MDS tertile sensitivity",
                              file.path(supp_dir, "TableS22_MDS_tertile_sensitivity.csv"))

doc <- add_page_break(doc)

# ---- figures ----
doc <- doc %>% body_add_par("Figures", style = "heading 1")

fig_s4 <- file.path(supp_dir, "miR-Mito Score Network.png")
doc <- add_figure(doc, "Figure S4. miRNA–target interaction network",
                  fig_s4,
                  caption = paste0("Network of the six signature miRNAs and their ",
                                   "predicted and validated target genes, constructed ",
                                   "in Cytoscape."))

fig_s5 <- file.path(supp_dir, "FigureS5_KM_CGGA.png")
cap_s5 <- file.path(supp_dir, "FigureS5_KM_CGGA_caption.txt")
capt <- if (file.exists(cap_s5)) paste(readLines(cap_s5, warn = FALSE), collapse = " ") else NULL
doc <- add_figure(doc, "Figure S5. Kaplan-Meier (CGGA external cohort)",
                  fig_s5, caption = capt, width_in = 6.5)

fig_s6 <- file.path(supp_dir, "FigureS6_Schoenfeld.png")
cap_s6 <- file.path(supp_dir, "FigureS6_Schoenfeld_caption.txt")
capt <- if (file.exists(cap_s6)) paste(readLines(cap_s6, warn = FALSE), collapse = " ") else NULL
doc <- add_figure(doc, "Figure S6. Schoenfeld residuals (test set)",
                  fig_s6, caption = capt, width_in = 5.5)

fig_s7 <- file.path(supp_dir, "FigureS7_Spline_HR.png")
cap_s7 <- file.path(supp_dir, "FigureS7_Spline_HR_caption.txt")
capt <- if (file.exists(cap_s7)) paste(readLines(cap_s7, warn = FALSE), collapse = " ") else NULL
doc <- add_figure(doc, "Figure S7. Spline-based hazard ratio",
                  fig_s7, caption = capt, width_in = 6)

fig_s8 <- file.path(supp_dir, "FigureS8_zscore_comparison.png")
cap_s8 <- file.path(supp_dir, "FigureS8_caption.txt")
capt <- if (file.exists(cap_s8)) paste(readLines(cap_s8, warn = FALSE), collapse = " ") else NULL
doc <- add_figure(doc, "Figure S8. Z-score distribution TCGA vs CGGA",
                  fig_s8, caption = capt, width_in = 6.5)

fig_s9 <- file.path(supp_dir, "FigureS9_enrichment_dotplot.png")
cap_s9 <- file.path(supp_dir, "FigureS9_caption.txt")
capt <- if (file.exists(cap_s9)) paste(readLines(cap_s9, warn = FALSE), collapse = " ") else NULL
doc <- add_figure(doc, "Figure S9. GO/KEGG enrichment dot plot",
                  fig_s9, caption = capt, width_in = 6.5)

doc <- add_page_break(doc)

# ---- Data ----
doc <- doc %>% body_add_par("Data Files", style = "heading 1")

doc <- add_table_with_caption(doc, "Data S1. Barcodes by split",
                              file.path(supp_dir, "DataS1_barcodes.csv"),
                              max_rows = 30)
doc <- add_table_with_caption(doc, "Data S2. Risk scores",
                              file.path(supp_dir, "DataS2_risk_scores.csv"),
                              max_rows = 30)
doc <- add_table_with_caption(doc, "Data S3. Package versions",
                              file.path(supp_dir, "DataS3_package_versions.csv"))

doc <- add_page_break(doc)

# ---- Supplementary Methods ----
doc <- doc %>% body_add_par("Supplementary Methods", style = "heading 1")

methods_path <- file.path(supp_dir, "SupplementaryMethods_draft.txt")
if (file.exists(methods_path)) {
  lines <- readLines(methods_path, warn = FALSE)
  for (line in lines) {
    if (grepl("^S[0-9]+\\.", line)) {
      doc <- doc %>% body_add_par(line, style = "heading 2")
    } else if (grepl("^[-=]+$", line) || !nzchar(trimws(line))) {
      next
    } else {
      doc <- doc %>% body_add_par(line, style = "Normal")
    }
  }
}

doc <- add_page_break(doc)

# ---- TRIPOD+AI ----
doc <- doc %>% body_add_par("TRIPOD+AI Checklist", style = "heading 1")
doc <- doc %>% body_add_par(
  "The TRIPOD+AI checklist (BMJ 2024) was completed for the development + evaluation (D;E) pathway.",
  style = "Normal")
doc <- doc %>% body_add_par("", style = "Normal")

tripod_path <- file.path(supp_dir, "TRIPOD_AI_checklist_prefilled.csv")
if (file.exists(tripod_path)) {
  info <- csv_to_flextable(tripod_path, font_size = 8)
  doc <- doc %>% body_add_flextable(info$ft)
}

# =============================================================================
# save
# =============================================================================
print(doc, target = out_file)

log_msg("")
log_msg("============================================================")
log_msg("DOCX built: ", out_file)
log_msg("   size: ", round(file.info(out_file)$size / 1024, 1), " KB")
log_msg("============================================================")
