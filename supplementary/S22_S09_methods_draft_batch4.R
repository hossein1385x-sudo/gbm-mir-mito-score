#!/usr/bin/env Rscript
# =============================================================================
# Supplementary Batch 4 — Table S22 + Figure S9 + Supplementary Methods draft
# -----------------------------------------------------------------------------
# Table S22  — MDS tertile sensitivity analysis
# Figure S9  — GO/KEGG enrichment dot plot (v2)
# Supplementary Methods — draft text (auto-generated from current pipeline)
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(DESeq2)
})

input_dir <- paste0(WORK_DIR, "/")
sig_dir   <- file.path(input_dir, "signature_v2")
out_dir   <- file.path(input_dir, "supplementary_final")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) cat(sprintf("[%s] %s\n",
                                     format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

norm_id <- function(x) substr(toupper(gsub("\\.", "-", as.character(x))), 1, 12)


# =============================================================================
# TABLE S22 — MDS tertile sensitivity
# =============================================================================
log_msg("")
log_msg("Table S22 — MDS tertile sensitivity")

mds_path  <- file.path(input_dir, "mds_scores_train.csv")
mir_path  <- file.path(input_dir, "mirna_isoform_counts_matrix.csv")

if (!file.exists(mds_path) || !file.exists(mir_path)) {
  log_msg("   required files are missing — Table S22 skipped.")
} else {

  mds <- fread(mds_path)
  mds[, key := norm_id(sample_id)]
  mds <- mds[!is.na(MDS)]

  # --- tertile split ---
  qs <- quantile(mds$MDS, probs = c(1/3, 2/3), na.rm = TRUE)
  mds[, tertile := ifelse(MDS <= qs[1], "low",
                          ifelse(MDS >= qs[2], "high", "mid"))]
  mds_keep <- mds[tertile != "mid"]

  log_msg(sprintf("   tertile: low=%d | mid=%d | high=%d",
                  sum(mds$tertile == "low"),
                  sum(mds$tertile == "mid"),
                  sum(mds$tertile == "high")))

  # --- miRNA matrix ---
  mir <- as.matrix(read.csv(mir_path, row.names = 1, check.names = FALSE))
  storage.mode(mir) <- "numeric"
  if (mean(grepl("^hsa-", rownames(mir))) < 0.5) mir <- t(mir)

  colnames(mir) <- norm_id(colnames(mir))
  if (anyDuplicated(colnames(mir))) {
    dt <- as.data.table(t(mir)); dt[, key := colnames(mir)]
    dt <- dt[, lapply(.SD, function(x) sum(x, na.rm = TRUE)), by = key,
             .SDcols = setdiff(names(dt), "key")]
    mir <- t(as.matrix(dt[, -1, with = FALSE]))
    colnames(mir) <- dt$key
  }

  common <- intersect(mds_keep$key, colnames(mir))
  log_msg(sprintf("   shared samples: %d", length(common)))

  if (length(common) >= 20) {
    mir_sub  <- round(mir[, common, drop = FALSE])
    mds_sub  <- mds_keep[match(common, mds_keep$key)]
    grp      <- factor(mds_sub$tertile, levels = c("low", "high"))

    col_data <- data.frame(row.names = common, tertile = grp)

    dds <- DESeqDataSetFromMatrix(
      countData = mir_sub,
      colData   = col_data,
      design    = ~ tertile)

    keep <- rowSums(counts(dds) >= 10) >= 3
    dds  <- dds[keep, ]
    dds  <- DESeq(dds, quiet = TRUE)

    res <- results(dds, contrast = c("tertile", "high", "low"), alpha = 0.05)
    res_df <- as.data.frame(res)
    res_df$isoform_id <- rownames(res_df)
    res_df <- res_df[order(res_df$padj), ]

    n_sig <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 &
                   abs(res_df$log2FoldChange) >= 1)

    # --- summary of the comparison with the median split ---
    summary_tbl <- data.frame(
      Metric = c("Split method",
                 "n low", "n high",
                 "Isoforms tested",
                 "Candidates (FDR<0.05 & |LFC|>=1)",
                 "Top up isoform",
                 "Top up log2FC",
                 "Top down isoform",
                 "Top down log2FC"),
      Value = c("Tertile (33rd/67th percentile)",
                sum(grp == "low"), sum(grp == "high"),
                nrow(res_df), n_sig,
                res_df$isoform_id[which.max(res_df$log2FoldChange)],
                round(max(res_df$log2FoldChange, na.rm = TRUE), 3),
                res_df$isoform_id[which.min(res_df$log2FoldChange)],
                round(min(res_df$log2FoldChange, na.rm = TRUE), 3)),
      stringsAsFactors = FALSE)

    write.csv(summary_tbl,
              file.path(out_dir, "TableS22_MDS_tertile_sensitivity.csv"),
              row.names = FALSE)

    log_msg("   TableS22 saved")
    print(summary_tbl, row.names = FALSE)
  }
}


# =============================================================================
# FIGURE S9 — GO/KEGG dot plot (v2)
# =============================================================================
log_msg("")
log_msg("Figure S9 — Enrichment dot plot")

go_path <- file.path(sig_dir, "15_target_genes/go_kegg_enrichment.csv")
if (!file.exists(go_path)) {
  log_msg("   go_kegg_enrichment.csv is missing — Figure S9 skipped.")
} else {
  go <- fread(go_path)

  # filter and sort
  go_sig <- go[!is.na(p.adjust) & p.adjust < 0.05]
  go_sig <- go_sig[order(p.adjust)]
  go_sig <- go_sig[, .SD[1:min(.N, 10)], by = ontology_source]

  # labels
  go_sig[, ontology_label := fcase(
    ontology_source == "GO_BP", "GO Biological Process",
    ontology_source == "GO_MF", "GO Molecular Function",
    ontology_source == "GO_CC", "GO Cellular Component",
    ontology_source == "KEGG",  "KEGG Pathway",
    default = ontology_source)]

  # GeneRatio to numeric
  go_sig[, GeneRatio_num := as.numeric(sub("/.*", "", GeneRatio)) /
           as.numeric(sub(".*/", "", GeneRatio))]

  p_s9 <- ggplot(go_sig,
                 aes(x = -log10(p.adjust),
                     y = reorder(Description, -log10(p.adjust)),
                     size = Count,
                     color = ontology_label)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = c(2, 7),
                          name = "Gene count") +
    scale_color_manual(
      values = c("GO Biological Process" = "#B2182B",
                 "GO Molecular Function" = "#2166AC",
                 "GO Cellular Component" = "#4DAF4A",
                 "KEGG Pathway" = "#FF7F00"),
      name = NULL) +
    facet_wrap(~ ontology_label, scales = "free_y", ncol = 1) +
    labs(
      title = "GO/KEGG enrichment of signature miRNA target genes",
      subtitle = "Top 10 terms per ontology (FDR < 0.05)",
      x = expression(-log[10]~"(FDR-adjusted P)"),
      y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "grey30"))

  height_calc <- max(6, 3 * length(unique(go_sig$ontology_label)))

  ggsave(file.path(out_dir, "FigureS9_enrichment_dotplot.png"),
         p_s9, width = 10, height = height_calc, dpi = 600, bg = "white")
  ggsave(file.path(out_dir, "FigureS9_enrichment_dotplot.pdf"),
         p_s9, width = 10, height = height_calc)

  writeLines(c(
    "Figure S9. Dot plot of the top 10 significantly enriched terms ",
    "(FDR-adjusted P < 0.05) per ontology category for the combined ",
    "predicted and validated target genes of the six signature miRNAs. ",
    "Point size indicates the number of genes in each term; color indicates ",
    "the ontology category (GO Biological Process, Molecular Function, ",
    "Cellular Component, or KEGG Pathway). Enrichment was performed using ",
    "clusterProfiler with the full set of genes present in the RNA-seq ",
    "dataset as the background universe. Terms were ranked by FDR-adjusted ",
    "P-value within each ontology category."
  ), file.path(out_dir, "FigureS9_caption.txt"))

  log_msg("   FigureS9 saved")
}


# =============================================================================
# SUPPLEMENTARY METHODS — draft text
# =============================================================================
log_msg("")
log_msg("Supplementary Methods — draft")

methods_text <- c(
  "============================================================",
  "Supplementary Methods",
  "============================================================",
  "",
  "S1. Mitochondrial Dysfunction Score (MDS)",
  "------------------------------------------------------------",
  "A curated list of nuclear- and mitochondrial-encoded mitochondrial",
  "genes was compiled from MitoCarta 3.0 and matched to the RNA-seq",
  "expression matrix. After removal of genes with zero counts and",
  "application of DESeq2 median-of-ratios size-factor normalization",
  "followed by log2(count + 1) transformation, gene-level z-score",
  "standardization was performed using training-set parameters only.",
  "The MDS for each patient was defined as the first principal component",
  "(PC1) of the standardized mitochondrial gene expression matrix,",
  "capturing 20.6% of total variance. PC1 was sign-oriented so that",
  "higher MDS values corresponded to higher overall mitochondrial gene",
  "expression. Patients were dichotomized at the training-set median",
  "MDS into MDS-high and MDS-low groups for downstream differential",
  "expression analysis. A Cox univariate test relating MDS to overall",
  "survival in the training set was not significant (HR = 0.92 per SD,",
  "P = 0.36), confirming that the MDS dichotomization operated as an",
  "outcome-independent enrichment filter.",
  "",
  "S2. Differential Expression Analysis",
  "------------------------------------------------------------",
  "miRNA isoforms in the training set were compared between MDS-high",
  "and MDS-low groups using DESeq2. Features with fewer than 10 counts",
  "in at least three samples were removed prior to testing. Candidate",
  "miRNAs were defined by Benjamini-Hochberg FDR-adjusted P < 0.05 and",
  "|log2FC| >= 1. This step was performed exclusively within the",
  "training set.",
  "",
  "S3. Elastic Net–Cox Model and Nested Cross-Validation",
  "------------------------------------------------------------",
  "An Elastic Net-penalized Cox proportional-hazards model (alpha = 0.5)",
  "was fit by minimizing the penalized negative partial log-likelihood:",
  "",
  "  beta_hat = argmin { ell(beta) + lambda * [ (1-alpha)/2 * ||beta||_2^2",
  "                                            + alpha * ||beta||_1 ] }",
  "",
  "where ell(beta) denotes the Cox partial log-likelihood. The",
  "regularization parameter lambda was selected via cv.glmnet using",
  "5-fold internal cross-validation. Model validation was performed",
  "using a repeated nested cross-validation scheme (5 outer folds x",
  "5 inner folds) repeated over 10 iterations, yielding up to 50 outer-",
  "fold performance evaluations. In each outer fold, the inner folds",
  "were used to select lambda, the model was refit on the full outer-",
  "training partition, and discrimination was assessed on the held-out",
  "outer-test partition using Harrell's concordance index (C-index).",
  "",
  "S4. Feature Selection Stability Analysis",
  "------------------------------------------------------------",
  "To assess selection stability, 50 independent 80/20 splits of the",
  "training set were performed. In each iteration, univariate Cox",
  "screening (p < 0.01, with a fallback to p < 0.05 if fewer than 20",
  "candidates survived) was restricted to the training subset, scaling",
  "parameters were derived exclusively from that subset, and an Elastic",
  "Net-Cox model (alpha = 0.5, 5-fold CV) was fit. Up to 10 miRNAs were",
  "retained per split (ranked by absolute coefficient magnitude), and",
  "the selected set was compared against the locked signature. Stability",
  "was quantified by the overlap count, the Jaccard similarity index,",
  "and the training/validation C-index.",
  "",
  "S5. Bootstrap Optimism Correction",
  "------------------------------------------------------------",
  "To estimate and correct for performance optimism, 1,000 bootstrap",
  "replicates were generated from the training set. In each replicate,",
  "resampling with replacement was performed and the entire modeling",
  "pipeline (screening, scaling, lambda selection, and Elastic Net-Cox",
  "fitting) was re-executed. The C-index was computed both on the",
  "bootstrap sample (C_boot) and on the original training set using",
  "the bootstrap-derived model (C_orig), and the difference was taken",
  "as the optimism. The mean optimism across 1,000 replicates was",
  "subtracted from the apparent (in-sample) C-index to obtain the",
  "optimism-corrected estimate. The test set was not involved at any",
  "stage of this procedure.",
  "",
  "S6. Locked Test Set Evaluation",
  "------------------------------------------------------------",
  "The final model, comprising six miRNAs and age, was locked with a",
  "SHA-256 fingerprint immediately after finalization. The signature",
  "was applied to the test set with no refitting, no re-selection of",
  "features, and no re-estimation of scaling parameters. For each test-",
  "set patient, a linear predictor (risk score) was computed using the",
  "fixed model coefficients. Discrimination was assessed with Harrell's",
  "C-index; hazard ratios were computed per standard deviation of the",
  "risk score using univariate and age-adjusted Cox regression.",
  "Time-dependent discrimination was evaluated using time-dependent",
  "AUC(t) at 6-month, 1-year, and 2-year horizons.",
  "",
  "S7. Calibration and Decision Curve Analysis",
  "------------------------------------------------------------",
  "Predicted survival probabilities were computed using the baseline",
  "survival function estimated from the training set. Calibration was",
  "assessed by grouping test-set patients into quartiles of predicted",
  "probability and comparing Kaplan-Meier-observed survival against",
  "model-predicted probability within each quartile. The IPCW-corrected",
  "time-dependent Brier score was computed at each horizon and compared",
  "against a Kaplan-Meier reference. Decision curve analysis (DCA) was",
  "performed across threshold probabilities of 0.01-0.50, comparing",
  "the signature against Treat-All and Treat-None reference strategies.",
  "",
  "S8. External Validation in CGGA",
  "------------------------------------------------------------",
  "The CGGA cohort (microRNA_array_198) was used as an independent",
  "external cohort. Because the array platform did not include all six",
  "signature miRNAs, only five components were evaluable (hsa-miR-196a-",
  "5p, hsa-miR-675-5p, hsa-miR-221-3p, hsa-miR-31-5p, and hsa-miR-181a-",
  "2-3p); hsa-miR-4660-5p was not represented. CGGA expression values",
  "were transformed using log2(x + 1), and the identical TCGA-locked",
  "centering/scaling parameters were applied without re-estimation.",
  "Risk scores were computed using the TCGA-locked coefficients, and",
  "discrimination was assessed using the C-index, univariate and age-",
  "adjusted Cox regression, and Kaplan-Meier analysis.",
  "",
  "S9. Functional Enrichment Analysis",
  "------------------------------------------------------------",
  "Target genes of the six signature miRNAs were identified by",
  "integrating experimentally validated interactions (miRTarBase) and",
  "computationally predicted interactions (miRDB v6.0, score >= 80).",
  "The background gene universe was defined as all genes present in the",
  "RNA-seq dataset. Gene Ontology (BP, MF, CC) and KEGG pathway",
  "enrichment were performed using clusterProfiler with the Benjamini-",
  "Hochberg correction. These analyses were interpreted as exploratory",
  "and hypothesis-generating.",
  "",
  "S10. miRNA–Target Expression Concordance",
  "------------------------------------------------------------",
  "To evaluate directional consistency between signature miRNAs and",
  "their predicted targets, Spearman correlations were computed between",
  "each miRNA and each target gene in the training set. Multiple-testing",
  "correction was applied using the Benjamini-Hochberg method. A pair",
  "was classified as concordant if the correlation was negative",
  "(rho < 0) and significant after FDR correction (FDR-adjusted",
  "P < 0.05).",
  "",
  "S11. Software and Reproducibility",
  "------------------------------------------------------------",
  "All analyses were performed in R (see Data S3 for version and",
  "package specifications). The full pipeline was organized as a",
  "sequence of numbered scripts (01-18) deposited in the project",
  "repository, and random seeds were fixed at 42 for the primary",
  "analyses and 1001-1050 for the 50-seed stability analysis. The",
  "SHA-256 fingerprint of the locked signature is recorded in Table S14.",
  "",
  "============================================================"
)

writeLines(methods_text,
           file.path(out_dir, "SupplementaryMethods_draft.txt"))

log_msg("   SupplementaryMethods_draft.txt saved")


# =============================================================================
# summary
# =============================================================================
log_msg("")
log_msg("============================================================")
log_msg("Supplementary Batch 4 complete")
log_msg("   TableS22_MDS_tertile_sensitivity.csv")
log_msg("   FigureS9_enrichment_dotplot.png / .pdf")
log_msg("   FigureS9_caption.txt")
log_msg("   SupplementaryMethods_draft.txt")
log_msg("============================================================")
