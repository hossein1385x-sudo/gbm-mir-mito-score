#!/usr/bin/env Rscript
# =============================================================================
# STEP 15b (Signature v2) — Cytoscape Input Generator
# -----------------------------------------------------------------------------
# Build input files for Cytoscape from the miRNA-target network
#
# Inputs:
#   signature_v2/15_target_genes/signature_target_genes.csv
#   signature_v2/signature_features.csv
#
# Outputs in signature_v2/15_target_genes/cytoscape/:
#   nodes_full.csv           (all miRNAs and genes)
#   edges_full.csv           (all pairs, with evidence)
#   nodes_validated.csv      (miRTarBase-validated targets only)
#   edges_validated.csv
#   nodes_compact.csv        (top-N targets per miRNA, for a readable network)
#   edges_compact.csv
#   README_cytoscape.txt     (loading instructions)
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

input_dir <- paste0(WORK_DIR, "/")
sig_dir   <- file.path(input_dir, "signature_v2")
out_dir   <- file.path(sig_dir, "15_target_genes/cytoscape")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

TOP_N_PER_MIRNA <- 15   # for the compact network

log_msg("============================================================")
log_msg("STEP 15b (Signature v2) — Cytoscape Input Generator")
log_msg("============================================================")

# ======================================================================
# 1) Load the miRNA-target pairs
# ======================================================================
targets_path <- file.path(sig_dir, "15_target_genes/signature_target_genes.csv")
if (!file.exists(targets_path)) stop("not found: ", targets_path)

targets <- fread(targets_path)
log_msg(sprintf("   pairs: %d | columns: %s",
                nrow(targets), paste(names(targets), collapse = ", ")))

# ======================================================================
# 2) Load the coefficients of the signature miRNAs
# ======================================================================
coef_path <- file.path(sig_dir, "signature_features.csv")
if (!file.exists(coef_path)) stop("not found: ", coef_path)

coefs <- fread(coef_path)
log_msg(sprintf("   coefficients: %d features", nrow(coefs)))

# separate miRNA from clinical
coefs_mirna <- coefs[grepl("^hsa-", Feature)]
coefs_mirna[, direction := ifelse(Coefficient > 0, "Risk", "Protective")]

# ======================================================================
# 3) Build the edges (Full)
# ======================================================================
edges <- as.data.table(targets)
setnames(edges,
         c("mirna", "target_gene_symbol", "evidence_source"),
         c("source", "target", "evidence"),
         skip_absent = TRUE)

edges[, interaction := "miRNA-target"]
edges[, weight := ifelse(evidence == "miRTarBase", 2, 1)]  # weight: validated > predicted

# de-duplicate (if a miRNA-target pair appears in both miRTarBase and miRDB, keep both rows)
edges <- unique(edges, by = c("source", "target", "evidence"))
edges <- edges[order(source, target)]

write.csv(edges, file.path(out_dir, "edges_full.csv"), row.names = FALSE)
log_msg(sprintf("   edges_full.csv: %d rows", nrow(edges)))

# ======================================================================
# 4) Build the nodes (Full)
# ======================================================================
mirna_nodes <- data.table(
  id            = coefs_mirna$Feature,
  node_type     = "miRNA",
  coefficient   = round(coefs_mirna$Coefficient, 4),
  direction     = coefs_mirna$direction,
  n_targets_all = as.integer(table(factor(edges$source,
                                          levels = coefs_mirna$Feature))),
  stringsAsFactors = FALSE)

gene_nodes <- data.table(
  id            = unique(edges$target),
  node_type     = "Gene",
  coefficient   = NA_real_,
  direction     = NA_character_,
  n_targets_all = NA_integer_)

nodes <- rbindlist(list(mirna_nodes, gene_nodes), fill = TRUE)
nodes <- unique(nodes, by = "id")

# add the number of miRTarBase hits per gene (useful for filtering)
gene_evidence_count <- edges[, .(
  n_validated = sum(evidence == "miRTarBase"),
  n_predicted = sum(evidence == "miRDB"),
  n_mirnas    = uniqueN(source)
), by = .(id = target)]

nodes <- merge(nodes, gene_evidence_count, by = "id", all.x = TRUE)
nodes[is.na(n_validated), n_validated := 0]
nodes[is.na(n_predicted), n_predicted := 0]
nodes[is.na(n_mirnas),    n_mirnas    := 0]

write.csv(nodes, file.path(out_dir, "nodes_full.csv"), row.names = FALSE)
log_msg(sprintf("   nodes_full.csv: %d rows (%d miRNAs + %d genes)",
                nrow(nodes),
                sum(nodes$node_type == "miRNA"),
                sum(nodes$node_type == "Gene")))

# ======================================================================
# 5) Validated network (miRTarBase only)
# ======================================================================
edges_validated <- edges[evidence == "miRTarBase"]
if (nrow(edges_validated) > 0) {
  genes_validated <- unique(edges_validated$target)
  nodes_validated <- nodes[id %in% c(coefs_mirna$Feature, genes_validated)]

  write.csv(edges_validated,
            file.path(out_dir, "edges_validated.csv"),
            row.names = FALSE)
  write.csv(nodes_validated,
            file.path(out_dir, "nodes_validated.csv"),
            row.names = FALSE)
  log_msg(sprintf("   validated: %d edge | %d node",
                  nrow(edges_validated), nrow(nodes_validated)))
}

# ======================================================================
# 6) Compact network (top-N targets per miRNA, combining miRTarBase + miRDB)
# ======================================================================
# priority: miRTarBase first, then miRDB
edges_compact <- edges[order(source, evidence, target)]
edges_compact <- edges_compact[, .SD[seq_len(min(.N, TOP_N_PER_MIRNA))],
                                by = source]

# if a miRNA has both miRTarBase and miRDB targets, miRTarBase takes priority
edges_compact <- unique(edges_compact, by = c("source", "target"))

genes_compact <- unique(edges_compact$target)
nodes_compact <- nodes[id %in% c(coefs_mirna$Feature, genes_compact)]

write.csv(edges_compact,
          file.path(out_dir, "edges_compact.csv"),
          row.names = FALSE)
write.csv(nodes_compact,
          file.path(out_dir, "nodes_compact.csv"),
          row.names = FALSE)
log_msg(sprintf("   compact (top %d per miRNA): %d edge | %d node",
                TOP_N_PER_MIRNA,
                nrow(edges_compact), nrow(nodes_compact)))

# ======================================================================
# 7) Statistical summary
# ======================================================================
summary_dt <- data.table(
  network        = c("Full", "Validated", "Compact"),
  n_nodes        = c(nrow(nodes),
                     if (exists("nodes_validated")) nrow(nodes_validated) else 0,
                     nrow(nodes_compact)),
  n_edges        = c(nrow(edges),
                     nrow(edges_validated),
                     nrow(edges_compact)),
  n_miRNAs       = c(sum(nodes$node_type == "miRNA"),
                     sum(nodes_validated$node_type == "miRNA"),
                     sum(nodes_compact$node_type == "miRNA")),
  n_genes        = c(sum(nodes$node_type == "Gene"),
                     sum(nodes_validated$node_type == "Gene"),
                     sum(nodes_compact$node_type == "Gene"))
)
write.csv(summary_dt,
          file.path(out_dir, "network_summary.csv"),
          row.names = FALSE)

log_msg("")
log_msg("============================================================")
log_msg("network summary")
log_msg("============================================================")
print(summary_dt)

# ======================================================================
# 8) README for Cytoscape
# ======================================================================
readme <- c(
  "============================================================",
  "Cytoscape Import Guide — miR-Mito Score Target Network",
  "============================================================",
  "",
  "FILES:",
  "  nodes_full.csv / edges_full.csv            — all pairs",
  "  nodes_validated.csv / edges_validated.csv  — miRTarBase only",
  "  nodes_compact.csv / edges_compact.csv      — top-N targets per miRNA",
  "",
  "IMPORT (Cytoscape):",
  "  1. File → Import → Network → File...",
  "     - Select: edges_full.csv",
  "     - Source node column: source",
  "     - Target node column: target",
  "     - Edge attribute columns: evidence, interaction, weight",
  "  2. File → Import → Table → File...",
  "     - Select: nodes_full.csv",
  "     - Key column for network: id",
  "     - Import all attributes",
  "",
  "NODE ATTRIBUTES (for styling):",
  "  node_type    : 'miRNA' or 'Gene'",
  "  coefficient  : Cox coefficient (miRNAs only)",
  "  direction    : 'Risk' (β>0) or 'Protective' (β<0)",
  "  n_targets_all: number of targets (miRNAs only)",
  "  n_validated  : number of miRTarBase-validated edges (genes)",
  "  n_predicted  : number of miRDB-predicted edges (genes)",
  "  n_mirnas     : number of connected signature miRNAs (genes)",
  "",
  "EDGE ATTRIBUTES (for styling):",
  "  evidence     : 'miRTarBase' (validated) or 'miRDB' (predicted)",
  "  interaction  : 'miRNA-target'",
  "  weight       : 2 (validated) or 1 (predicted)",
  "",
  "SUGGESTED STYLING:",
  "  - Node shape: ellipse for miRNA, round-rectangle for gene",
  "  - Node fill: red = Risk miRNA; blue = Protective miRNA; grey = gene",
  "  - Node size: map to n_targets_all (miRNAs) or n_mirnas (genes)",
  "  - Edge color: solid black for miRTarBase, dashed grey for miRDB",
  "  - Edge width: map to weight",
  "",
  "EXPORT:",
  "  - File → Export → Network to Image (PDF, 600 dpi)",
  "  - For publication, use: PostScript or SVG export",
  "",
  "============================================================"
)
writeLines(readme, file.path(out_dir, "README_cytoscape.txt"))

log_msg("")
log_msg("outputs in: ", out_dir)
log_msg("DONE")
