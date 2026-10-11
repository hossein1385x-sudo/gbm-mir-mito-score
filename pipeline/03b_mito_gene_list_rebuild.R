#!/usr/bin/env Rscript
# =============================================================================
# Build the mitochondrial gene list from a local MitoCarta 3.0 file and clean the RNA-seq matrix
# (supersedes the earlier step 03, which downloaded MitoCarta and had fallback sources)
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table)
})

log_msg <- function(...) cat(sprintf("[%s] %s\n",
                                     format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

OUT_DIR <- WORK_DIR

## =============================================================================
## 1) Locate the MitoCarta file
## =============================================================================
log_msg("===== Locating the MitoCarta file =====")

candidates <- list.files(OUT_DIR, pattern = "MitoCarta",
                         full.names = TRUE, ignore.case = TRUE)
if (length(candidates) == 0) stop("MitoCarta file not found.")

xls_files <- candidates[grepl("\\.(xls|xlsx)$", candidates, ignore.case = TRUE)]
mito_file <- if (length(xls_files) > 0) xls_files[1] else candidates[1]
log_msg(sprintf("file: %s (%.1f KB)", basename(mito_file),
                file.info(mito_file)$size / 1024))

## =============================================================================
## 2) Read the correct sheet
## =============================================================================
log_msg("")
log_msg("===== Reading the MitoCarta sheet =====")

if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}

sheets <- readxl::excel_sheets(mito_file)
log_msg(sprintf("available sheets: %s", paste(sheets, collapse = " | ")))

# select the sheet: A Human MitoCarta3.0
target <- sheets[grepl("^A\\s", sheets)][1]
if (is.na(target)) target <- sheets[grepl("MitoCarta3", sheets)][1]
if (is.na(target)) target <- sheets[2]   # default: second sheet
log_msg(sprintf("reading sheet: %s", target))

mito_raw <- as.data.table(readxl::read_excel(mito_file, sheet = target))
log_msg(sprintf("number of rows: %d", nrow(mito_raw)))
log_msg(sprintf("columns: %s", paste(head(names(mito_raw), 10), collapse = ", ")))

# Symbol column
symbol_col <- intersect(c("Symbol", "Gene Symbol", "gene_symbol"), names(mito_raw))[1]
if (is.na(symbol_col)) stop("Symbol column not found.")

# Entrez column
entrez_col <- intersect(c("EntrezID", "Entrez ID", "entrez_id"), names(mito_raw))[1]

symbols <- unique(na.omit(as.character(mito_raw[[symbol_col]])))
symbols <- symbols[nzchar(symbols)]
log_msg(sprintf("number of unique symbols: %d", length(symbols)))

## =============================================================================
## 3) Install/load org.Hs.eg.db — forced install
## =============================================================================
log_msg("")
log_msg("===== Preparing org.Hs.eg.db =====")

# check that the install is real: try to load it
orgdb_ok <- tryCatch({
  suppressPackageStartupMessages(library(org.Hs.eg.db, character.only = TRUE))
  TRUE
}, error = function(e) FALSE)

if (!orgdb_ok) {
  log_msg("org.Hs.eg.db is not available — reinstalling ...")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  # remove the previous incomplete install (if any)
  tryCatch(
    remove.packages("org.Hs.eg.db"),
    error = function(e) NULL
  )

  BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE, force = TRUE)

  orgdb_ok <- tryCatch({
    suppressPackageStartupMessages(library(org.Hs.eg.db, character.only = TRUE))
    TRUE
  }, error = function(e) FALSE)
}

if (!orgdb_ok) {
  log_msg("org.Hs.eg.db could not be installed. Trying biomaRt ...")
}

## =============================================================================
## 4) Mapping
## =============================================================================
log_msg("")
log_msg("===== Mapping to ensembl_id =====")

mapping <- NULL
method_used <- NA_character_

## Route 1: org.Hs.eg.db
if (orgdb_ok) {
  log_msg("using org.Hs.eg.db ...")
  res <- tryCatch(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = symbols,
      columns = c("SYMBOL", "ENSEMBL", "ENTREZID"),
      keytype = "SYMBOL"
    ),
    error = function(e) { log_msg("  error: ", conditionMessage(e)); NULL }
  )
  if (!is.null(res) && nrow(res) > 0) {
    mapping <- as.data.table(res)
    mapping <- mapping[!is.na(ENSEMBL) & ENSEMBL != ""]
    setnames(mapping, c("SYMBOL", "ENSEMBL", "ENTREZID"),
             c("gene_symbol", "ensembl_id", "entrez_id"))
    mapping[, ensembl_id := gsub("\\.[0-9]+$", "", ensembl_id)]
    mapping <- unique(mapping, by = "gene_symbol")
    method_used <- "org.Hs.eg.db"
    log_msg(sprintf("  %d mappings", nrow(mapping)))
  }
}

## Route 2: biomaRt (if org.Hs.eg.db is unavailable)
if (is.null(mapping) && requireNamespace("biomaRt", quietly = TRUE)) {
  log_msg("using biomaRt ...")
  for (h in c("https://e116.ensembl.org",
              "https://www.ensembl.org",
              "https://e115.ensembl.org",
              "https://e114.ensembl.org",
              "https://useast.ensembl.org")) {
    log_msg(sprintf("  trying: %s", h))
    ens <- tryCatch(
      biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                       dataset = "hsapiens_gene_ensembl", host = h),
      error = function(e) NULL)
    if (is.null(ens)) { log_msg("    connection failed"); next }

    res <- tryCatch(
      biomaRt::getBM(
        attributes = c("external_gene_name", "ensembl_gene_id", "entrezgene_id"),
        filters    = "external_gene_name",
        values     = symbols, mart = ens),
      error = function(e) NULL)
    if (!is.null(res) && nrow(res) > 0) {
      mapping <- as.data.table(res)
      setnames(mapping, c("external_gene_name", "ensembl_gene_id", "entrezgene_id"),
               c("gene_symbol", "ensembl_id", "entrez_id"))
      mapping[, ensembl_id := gsub("\\.[0-9]+$", "", ensembl_id)]
      mapping <- unique(mapping, by = "gene_symbol")
      method_used <- paste0("biomaRt ", h)
      log_msg(sprintf("  %d mappings", nrow(mapping)))
      break
    }
  }
}

## Route 3: convertid package
if (is.null(mapping)) {
  log_msg("using convertid ...")
  if (!requireNamespace("convertid", quietly = TRUE)) {
    install.packages("convertid")
  }
  if (requireNamespace("convertid", quietly = TRUE)) {
    res <- tryCatch(
      convertid::convertId2(symbols, from = "SYMBOL",
                            to = "ENSEMBL", species = "Human"),
      error = function(e) NULL)
    if (!is.null(res) && length(res) > 0) {
      res <- as.data.table(res)
      if (all(c("from", "to") %in% names(res))) {
        setnames(res, c("from", "to"), c("gene_symbol", "ensembl_id"))
        res[, ensembl_id := gsub("\\.[0-9]+$", "", ensembl_id)]
        mapping <- res[!is.na(ensembl_id)]
        method_used <- "convertid"
        log_msg(sprintf("  %d mappings", nrow(mapping)))
      }
    }
  }
}

if (is.null(mapping) || nrow(mapping) == 0) {
  stop("No mapping method worked. ",
       "Manual fix: put a gene_symbol→ensembl_id mapping file in ",
       file.path(OUT_DIR, "gene_symbol_to_ensembl.csv"),
       " with two columns, gene_symbol and ensembl_id.")
}

log_msg(sprintf("mapping with %s: %d genes", method_used, nrow(mapping)))

## =============================================================================
## 5) Build the final mito_gene_list
## =============================================================================
mito_final <- merge(
  data.table(gene_symbol = symbols),
  mapping[, .(gene_symbol, ensembl_id, entrez_id)],
  by = "gene_symbol", all.x = TRUE
)

mito_final <- mito_final[!is.na(ensembl_id) & ensembl_id != ""]
mito_final[, source := "MitoCarta3.0"]
mito_final[, evidence := NA_character_]

setcolorder(mito_final, c("gene_symbol", "entrez_id",
                          "source", "evidence", "ensembl_id"))

log_msg(sprintf("final mito_gene_list: %d genes (from %d initial symbols)",
                nrow(mito_final), length(symbols)))

mito_out <- file.path(OUT_DIR, "mito_gene_list.csv")
fwrite(mito_final, mito_out)
log_msg(sprintf("saved: %s", mito_out))

## =============================================================================
## 6) Clean the RNA-seq matrix
## =============================================================================
log_msg("")
log_msg("===== Cleaning the RNA-seq matrix =====")

rna_path_in  <- file.path(OUT_DIR, "raw_rnaseq_counts_matrix.rds")
rna_path_out <- file.path(OUT_DIR, "raw_rnaseq_counts_matrix_clean.rds")

if (!file.exists(rna_path_in)) stop("raw_rnaseq_counts_matrix.rds not found.")

rna <- readRDS(rna_path_in)
if (!is.matrix(rna)) rna <- as.matrix(rna)
storage.mode(rna) <- "numeric"

log_msg(sprintf("before: %d genes × %d samples", nrow(rna), ncol(rna)))

star_stats <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
keep <- !(rownames(rna) %in% star_stats)
if (sum(!keep) > 0) {
  log_msg(sprintf("removed %d STAR rows: %s", sum(!keep),
                  paste(rownames(rna)[!keep], collapse = ", ")))
}

rna_clean <- rna[keep, , drop = FALSE]
log_msg(sprintf("after: %d genes × %d samples", nrow(rna_clean), ncol(rna_clean)))

saveRDS(rna_clean, rna_path_out)
log_msg(sprintf("saved: %s", rna_path_out))

## =============================================================================
## 7) Final check
## =============================================================================
log_msg("")
log_msg("===== Final check =====")

rn_norm  <- gsub("\\.[0-9]+$", "", rownames(rna_clean))
mito_ids <- gsub("\\.[0-9]+$", "", mito_final$ensembl_id)
common   <- intersect(mito_ids, rn_norm)

log_msg(sprintf("  mito genes: %d", length(mito_ids)))
log_msg(sprintf("  RNA-seq genes: %d", length(rn_norm)))
log_msg(sprintf("  overlap: %d genes", length(common)))

if (length(common) >= 800) {
  log_msg("Everything is ready.")
} else {
  log_msg("Overlap is lower than expected.")
  log_msg("  mito sample: ", paste(head(mito_ids, 3), collapse = ", "))
  log_msg("  RNA sample:  ", paste(head(rn_norm, 3), collapse = ", "))
}

log_msg("")
log_msg("Next step: put the clean file in place of the original and run Step 05:")
log_msg('  file.copy("raw_rnaseq_counts_matrix_clean.rds",')
log_msg('            "raw_rnaseq_counts_matrix.rds", overwrite = TRUE)')
log_msg('  run_mds_construction(input_root = "...", force = TRUE)')
