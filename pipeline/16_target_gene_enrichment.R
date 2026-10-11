#!/usr/bin/env Rscript
# =============================================================================
# STEP 15 (Signature v2) — Target Gene GO/KEGG Enrichment
# -----------------------------------------------------------------------------
# GO (BP/MF/CC) and KEGG enrichment of the target genes of the v2 signature miRNAs
#
# Inputs:
#   signature_v2/signature_v2_locked_latest.rds
#   raw_rnaseq_counts_matrix.rds
#   reference/mirtarbase_MTI.csv
#   reference/miRDB_v6.0_prediction_result.txt
#
# Outputs in signature_v2/15_target_genes/
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

# ======================================================================
# 0) Settings
# ======================================================================
output_dir <- paste0(WORK_DIR, "/")
sig_dir    <- file.path(output_dir, "signature_v2")
out_dir    <- file.path(sig_dir, "15_target_genes")
ref_dir    <- file.path(output_dir, "reference")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

path_signature   <- file.path(sig_dir, "signature_v2_locked_latest.rds")
path_rnaseq_expr <- file.path(output_dir, "raw_rnaseq_counts_matrix.rds")
path_mirtarbase  <- file.path(ref_dir, "mirtarbase_MTI.csv")
path_mirdb       <- file.path(ref_dir, "miRDB_v6.0_prediction_result.txt")

MIRDB_SCORE_CUTOFF <- 80
TOP_N_PER_CATEGORY <- 10

log_msg <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

log_msg("============================================================")
log_msg("STEP 15 (Signature v2) — Target Gene Enrichment")
log_msg("============================================================")

for (f in c(path_signature, path_rnaseq_expr, path_mirtarbase, path_mirdb)) {
  if (!file.exists(f)) stop("not found: ", f)
}

# ======================================================================
# 1) Load the v2 signature miRNAs
# ======================================================================
log_msg("Loading the v2 signature ...")

signature_locked <- readRDS(path_signature)

signature_mirnas <- if (!is.null(signature_locked$mirna_features)) {
  signature_locked$mirna_features
} else if (!is.null(signature_locked$selected_features)) {
  signature_locked$selected_features
} else {
  stop("field mirna_features is missing from signature_v2_locked_latest.rds.")
}

log_msg(sprintf("   miRNAs in the signature (%d): %s",
                length(signature_mirnas),
                paste(signature_mirnas, collapse = ", ")))

# light normalisation
normalize_mirna <- function(x) {
  x <- str_replace(x, "^hsa-mir-", "hsa-miR-")
  x
}
signature_mirnas_mature <- normalize_mirna(signature_mirnas)

# build alternative forms (with/without locus, with/without arm)
build_mirna_candidates <- function(x) {
  locus_stripped <- str_replace(
    x, "^(hsa-(?:miR|let)-[A-Za-z0-9]+)-[0-9]+-(3p|5p)$", "\\1-\\2")
  arm_locus_stripped <- str_replace(
    x, "^(hsa-(?:miR|let)-[A-Za-z0-9]+)(?:-[0-9]+)?-(3p|5p)$", "\\1")
  unique(c(x, locus_stripped, arm_locus_stripped))
}

mirna_candidate_map <- purrr::map_dfr(signature_mirnas_mature, function(m) {
  data.frame(canonical_mirna = m,
             candidate       = build_mirna_candidates(m),
             stringsAsFactors = FALSE)
})

# ======================================================================
# 2) miRTarBase
# ======================================================================
log_msg("")
log_msg("miRTarBase ...")

mirtarbase_raw <- read_csv(path_mirtarbase, show_col_types = FALSE)

# species filter — both "Homo sapiens" and "hsa"
species_col <- intersect(c("Species (miRNA)", "Species"), names(mirtarbase_raw))[1]
if (!is.na(species_col)) {
  mirtarbase_raw <- mirtarbase_raw %>%
    filter(.data[[species_col]] %in% c("Homo sapiens", "hsa"))
}
log_msg(sprintf("   after species filter: %d rows", nrow(mirtarbase_raw)))

if (nrow(mirtarbase_raw) == 0)
  warning("miRTarBase is empty after filtering.")

mirtarbase_filtered <- mirtarbase_raw %>%
  filter(`Support Type` == "Functional MTI") %>%
  inner_join(mirna_candidate_map, by = c("miRNA" = "candidate")) %>%
  transmute(mirna = canonical_mirna,
            target_gene_symbol = `Target Gene`,
            evidence_source = "miRTarBase") %>%
  distinct()

log_msg(sprintf("   final miRTarBase: %d miRNA-target pairs", nrow(mirtarbase_filtered)))

# ======================================================================
# 3) miRDB v6.0
# ======================================================================
log_msg("")
log_msg("miRDB v6.0 ...")

mirdb_raw <- read_tsv(
  path_mirdb,
  col_names = c("mirna", "refseq_transcript", "score"),
  show_col_types = FALSE
)
log_msg(sprintf("   raw miRDB: %d rows", nrow(mirdb_raw)))

mirdb_matched <- mirdb_raw %>%
  filter(score >= MIRDB_SCORE_CUTOFF) %>%
  inner_join(mirna_candidate_map, by = c("mirna" = "candidate")) %>%
  mutate(refseq_clean = str_remove(refseq_transcript, "\\.[0-9]+$"))

log_msg(sprintf("   after score>=%d filter: %d rows (%d unique RefSeq)",
                MIRDB_SCORE_CUTOFF, nrow(mirdb_matched),
                length(unique(mirdb_matched$refseq_clean))))

# RefSeq → SYMBOL
refseq_to_symbol <- tryCatch(
  AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = unique(mirdb_matched$refseq_clean),
    column = "SYMBOL",
    keytype = "REFSEQ",
    multiVals = "first"
  ),
  error = function(e) { warning(conditionMessage(e)); NA }
)

mirdb_filtered <- mirdb_matched %>%
  mutate(target_gene_symbol = refseq_to_symbol[refseq_clean]) %>%
  filter(!is.na(target_gene_symbol)) %>%
  transmute(mirna = canonical_mirna,
            target_gene_symbol,
            evidence_source = "miRDB") %>%
  distinct()

log_msg(sprintf("   final miRDB: %d pairs", nrow(mirdb_filtered)))

# ======================================================================
# 4) Combine and report coverage
# ======================================================================
target_gene_table <- bind_rows(mirtarbase_filtered, mirdb_filtered) %>%
  distinct(mirna, target_gene_symbol, evidence_source)

unique_target_genes <- target_gene_table %>%
  distinct(target_gene_symbol) %>%
  pull(target_gene_symbol)

log_msg("")
log_msg(sprintf("total pairs: %d | unique target genes: %d",
                nrow(target_gene_table), length(unique_target_genes)))

if (length(unique_target_genes) < 10)
  warning("the number of unique target genes is very low.")

# per-miRNA coverage
per_mirna_coverage <- target_gene_table %>%
  count(mirna, evidence_source) %>%
  pivot_wider(names_from = evidence_source, values_from = n, values_fill = 0)

missing_mirnas <- setdiff(signature_mirnas_mature, per_mirna_coverage$mirna)
if (length(missing_mirnas) > 0) {
  log_msg(sprintf("miRNAs without targets (%d): %s",
                  length(missing_mirnas), paste(missing_mirnas, collapse = ", ")))
  per_mirna_coverage <- bind_rows(
    per_mirna_coverage,
    data.frame(mirna = missing_mirnas)
  )
}

write_csv(per_mirna_coverage,
          file.path(out_dir, "mirna_target_coverage.csv"))
write_csv(target_gene_table,
          file.path(out_dir, "signature_target_genes.csv"))

# ======================================================================
# 5) Background universe from RNA-seq
# ======================================================================
log_msg("")
log_msg("building the background ...")

rnaseq_obj <- readRDS(path_rnaseq_expr)

raw_ids <- if (is.matrix(rnaseq_obj) || is.data.frame(rnaseq_obj)) {
  # gene × sample matrix
  if (nrow(rnaseq_obj) >= ncol(rnaseq_obj)) rownames(rnaseq_obj) else colnames(rnaseq_obj)
} else {
  stop("format of raw_rnaseq_counts_matrix.rds not recognised.")
}
raw_ids <- raw_ids[!is.na(raw_ids) & raw_ids != ""]

# remove the version suffix
raw_ids <- gsub("\\.[0-9]+$", "", raw_ids)
raw_ids <- unique(raw_ids)

log_msg(sprintf("   RNA-seq IDs: %d", length(raw_ids)))

# ======================================================================
# 6) Convert to Entrez
# ======================================================================
convert_to_entrez <- function(ids, keytype) {
  ids <- trimws(as.character(ids))
  ids <- ids[!is.na(ids) & ids != ""]
  if (length(ids) == 0) return(data.frame(ORIGINAL = character(), ENTREZID = character()))

  mapped <- tryCatch(
    AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids,
                          column = "ENTREZID", keytype = keytype,
                          multiVals = "first"),
    error = function(e) { warning(conditionMessage(e)); NA }
  )
  valid <- !is.na(mapped)
  if (sum(valid) == 0) return(data.frame(ORIGINAL = character(), ENTREZID = character()))
  data.frame(ORIGINAL = names(mapped)[valid],
             ENTREZID = as.character(mapped[valid]),
             stringsAsFactors = FALSE)
}

target_entrez     <- convert_to_entrez(unique_target_genes, "SYMBOL")
background_entrez <- convert_to_entrez(raw_ids,             "ENSEMBL")

log_msg(sprintf("   targets mapped to Entrez: %d of %d",
                nrow(target_entrez), length(unique_target_genes)))
log_msg(sprintf("   background mapped to Entrez: %d", nrow(background_entrez)))

if (nrow(target_entrez) == 0)
  stop("no target could be converted to Entrez.")
if (nrow(background_entrez) == 0) {
  warning("Background is empty — using the whole human genome.")
  background_entrez <- convert_to_entrez(
    keys(org.Hs.eg.db, keytype = "SYMBOL"), "SYMBOL")
}

# ======================================================================
# 7) GO enrichment
# ======================================================================
log_msg("")
log_msg("GO enrichment ...")

run_go <- function(ont_code, ont_label) {
  res <- enrichGO(
    gene          = target_entrez$ENTREZID,
    universe      = background_entrez$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = ont_code,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  as.data.frame(res) %>% mutate(ontology_source = ont_label)
}

go_bp <- run_go("BP", "GO_BP")
go_mf <- run_go("MF", "GO_MF")
go_cc <- run_go("CC", "GO_CC")

# ======================================================================
# 8) KEGG enrichment
# ======================================================================
log_msg("KEGG enrichment ...")

kegg_res <- tryCatch(
  enrichKEGG(
    gene          = target_entrez$ENTREZID,
    universe      = background_entrez$ENTREZID,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  ),
  error = function(e) { warning(conditionMessage(e)); NULL }
)

kegg_df <- if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
  as.data.frame(kegg_res) %>% mutate(ontology_source = "KEGG")
} else NULL

# ======================================================================
# 9) Combine and save
# ======================================================================
all_enrichment <- bind_rows(go_bp, go_mf, go_cc, kegg_df)

write_csv(all_enrichment,
          file.path(out_dir, "go_kegg_enrichment.csv"))

# summary
n_sig <- all_enrichment %>%
  filter(p.adjust < 0.05) %>%
  count(ontology_source)

log_msg("")
log_msg("Significant terms (FDR<0.05) per ontology:")
print(as.data.frame(n_sig))

# ======================================================================
# 10) Dotplot
# ======================================================================
top_terms <- all_enrichment %>%
  filter(p.adjust < 0.05) %>%
  group_by(ontology_source) %>%
  slice_min(order_by = p.adjust, n = TOP_N_PER_CATEGORY) %>%
  ungroup()

if (nrow(top_terms) > 0) {
  p_dot <- ggplot(top_terms,
                  aes(x = -log10(p.adjust),
                      y = reorder(Description, -log10(p.adjust)),
                      size = Count,
                      color = ontology_source)) +
    geom_point() +
    facet_wrap(~ ontology_source, scales = "free_y", ncol = 1) +
    labs(
      title    = "GO/KEGG Enrichment — Signature v2 Target Genes",
      subtitle = sprintf("Sources: miRTarBase (validated) + miRDB v6.0 (score>=%d)",
                         MIRDB_SCORE_CUTOFF),
      x = "-log10(adjusted p-value)",
      y = NULL, size = "Gene count", color = "Ontology"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))

  ggsave(file.path(out_dir, "go_kegg_dotplot.pdf"),
         p_dot, width = 9, height = max(6, 3 * length(unique(top_terms$ontology_source))))
  ggsave(file.path(out_dir, "go_kegg_dotplot_600dpi.png"),
         p_dot, width = 9, height = max(6, 3 * length(unique(top_terms$ontology_source))),
         dpi = 600)
  log_msg("   go_kegg_dotplot.pdf + .png")
} else {
  log_msg("No significant term (FDR<0.05) found.")
}

# ======================================================================
# 11) Summary
# ======================================================================
log_msg("")
log_msg("============================================================")
log_msg("Summary of STEP 15 (Signature v2)")
log_msg("============================================================")
log_msg(sprintf("miRNAs in the signature: %d", length(signature_mirnas_mature)))
log_msg(sprintf("miRNA-target pairs: %d", nrow(target_gene_table)))
log_msg(sprintf("unique target genes: %d", length(unique_target_genes)))
log_msg(sprintf("miRNAs without targets: %d", length(missing_mirnas)))
log_msg(sprintf("Significant terms FDR<0.05: %d", sum(all_enrichment$p.adjust < 0.05, na.rm = TRUE)))
log_msg("")
log_msg("outputs in: ", out_dir)
log_msg("STEP 15 (Signature v2) COMPLETED")
log_msg("This analysis is hypothesis-generating.")
