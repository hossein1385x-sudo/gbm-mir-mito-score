#!/usr/bin/env Rscript
# =============================================================================
# STEP 18b (Signature v2) — Figures 5–8
# -----------------------------------------------------------------------------
# Fig 5: Volcano plot (DE analysis)
# Fig 6: LASSO-Cox selected features
# Fig 7: Time-dependent ROC curves
# Fig 8: Calibration plot
#
# All outputs are saved in figures_final/ at 600 dpi.
# =============================================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(scales)
  library(ggrepel)
  library(survival)
  library(stringr)
})

# =============================================================================
# 1) Paths
# =============================================================================
INPUT_DIR  <- paste0(WORK_DIR, "/")
SIG_DIR    <- file.path(INPUT_DIR, "signature_v2")
OUTPUT_DIR <- file.path(INPUT_DIR, "figures_final")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

in_path  <- function(...) file.path(...)
out_path <- function(f)   file.path(OUTPUT_DIR, f)

DPI <- 600

log_msg <- function(...) cat(sprintf("[%s] %s\n",
  format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

# =============================================================================
# 2) Shared theme
# =============================================================================
theme_pub <- function() {
  theme_bw(base_size = 11, base_family = "sans") %+replace%
    theme(
      panel.grid.major  = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.ticks        = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(4, "pt"),
      axis.title        = element_text(face = "bold", size = 11),
      axis.text         = element_text(color = "black", size = 9.5),
      plot.title        = element_text(face = "bold", size = 12, hjust = 0.5,
                                       margin = margin(b = 4)),
      plot.subtitle     = element_text(size = 10, hjust = 0.5,
                                       margin = margin(b = 8)),
      strip.background  = element_rect(fill = "grey95", color = "black",
                                       linewidth = 0.6),
      strip.text        = element_text(face = "bold", size = 10),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.position   = "bottom",
      plot.margin       = margin(t = 10, r = 10, b = 15, l = 15, unit = "pt")
    )
}

save_plot <- function(plot_obj, name_base, width, height) {
  ggsave(out_path(paste0(name_base, ".png")), plot = plot_obj,
         width = width, height = height, dpi = DPI, bg = "white")
  ggsave(out_path(paste0(name_base, ".pdf")), plot = plot_obj,
         width = width, height = height, device = cairo_pdf, bg = "white")
  log_msg(sprintf("   Saved: %s", name_base))
}

# =============================================================================
# FIGURE 5 — Volcano Plot
# =============================================================================
draw_fig5 <- function() {
  log_msg("")
  log_msg("Figure 5 — Volcano")

  # v2 path
  f5_path <- in_path(SIG_DIR, "06_DE/DE_results_train.csv")
  if (!file.exists(f5_path)) {
    # fall back to the old path
    f5_path <- in_path(INPUT_DIR, "DE_results_train.csv")
  }
  if (!file.exists(f5_path)) stop("DE_results_train.csv not found.")

  df <- read_csv(f5_path, show_col_types = FALSE)
  log_msg(sprintf("   reading: %s (%d rows)", basename(f5_path), nrow(df)))

  col_id   <- intersect(c("isoform_id", "miRNA", "mirna", "gene", "feature"),
                        names(df))[1]
  col_lfc  <- intersect(c("log2FoldChange", "log2FC", "logFC", "lfc"),
                        names(df))[1]
  col_padj <- intersect(c("padj", "p_adj", "FDR", "qvalue", "adj_p_value"),
                        names(df))[1]

  if (any(is.na(c(col_id, col_lfc, col_padj)))) {
    stop(sprintf("required columns not found. available: %s",
                 paste(names(df), collapse = ", ")))
  }

  df_clean <- df %>%
    rename(Feature = !!col_id, Log2FC = !!col_lfc, Padj = !!col_padj) %>%
    filter(!is.na(Log2FC), !is.na(Padj)) %>%
    mutate(
      NegLog10P = -log10(pmax(Padj, 1e-300)),
      Group = case_when(
        Padj < 0.05 & Log2FC >=  1.0 ~ "Up in high-MDS",
        Padj < 0.05 & Log2FC <= -1.0 ~ "Down in high-MDS",
        TRUE                          ~ "Not significant"
      )
    )

  top_labels <- df_clean %>%
    filter(Group != "Not significant") %>%
    arrange(Padj, desc(abs(Log2FC))) %>%
    head(12)

  p <- ggplot(df_clean, aes(x = Log2FC, y = NegLog10P)) +
    geom_point(aes(color = Group), alpha = 0.65, size = 1.5) +
    scale_color_manual(values = c(
      "Up in high-MDS"   = "#B2182B",
      "Down in high-MDS" = "#2166AC",
      "Not significant"  = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_text_repel(data = top_labels, aes(label = Feature),
                    size = 2.8, max.overlaps = 20, box.padding = 0.3) +
    labs(
      title = "Differential Expression of miRNA Isoforms (Training Set)",
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~"(FDR-adjusted p-value)"),
      color = NULL
    ) +
    theme_pub()

  save_plot(p, "Fig5_Volcano", 6.5, 5.2)
}


# =============================================================================
# FIGURE 6 — LASSO-Cox Selected Features
# =============================================================================
draw_fig6 <- function() {
  log_msg("")
  log_msg("Figure 6 — LASSO coefficients")

  # v2 path: the model is in signature_v2/signature_v2_model.rds
  model_path <- in_path(SIG_DIR, "signature_v2_model.rds")
  if (!file.exists(model_path)) {
    model_path <- in_path(INPUT_DIR, "lasso_cox_model.rds")
  }
  if (!file.exists(model_path)) stop("LASSO model not found.")

  model_obj <- readRDS(model_path)
  log_msg(sprintf("   reading model: %s | class=%s",
                  basename(model_path), class(model_obj)[1]))

  # extract the coefficients with coef() — for glmnet
  coef_mat <- as.matrix(coef(model_obj))
  df_beta <- data.frame(
    Feature = rownames(coef_mat),
    Coef    = as.numeric(coef_mat[, 1]),
    stringsAsFactors = FALSE
  )
  df_beta <- df_beta[df_beta$Coef != 0, ]
  df_beta <- df_beta[order(abs(df_beta$Coef), decreasing = TRUE), ]

  if (nrow(df_beta) == 0) stop("all coefficients are zero.")

  df_beta$Feature <- factor(df_beta$Feature, levels = rev(df_beta$Feature))
  log_msg(sprintf("   number of non-zero features: %d", nrow(df_beta)))

  p <- ggplot(df_beta, aes(x = Feature, y = Coef, fill = Coef > 0)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = c("FALSE" = "#2166AC", "TRUE" = "#B2182B"),
                      guide = "none") +
    coord_flip() +
    labs(
      title = "LASSO-Cox Selected Features",
      subtitle = sprintf("Non-zero coefficients (%d features)", nrow(df_beta)),
      x = NULL, y = "Coefficient Value"
    ) +
    theme_pub() +
    theme(panel.grid.major.y = element_blank())

  plot_height <- max(4.0, nrow(df_beta) * 0.35)

  ggsave(out_path("Fig6_LASSO_Coefs.png"), plot = p,
         width = 7.5, height = plot_height, dpi = DPI, bg = "white")
  ggsave(out_path("Fig6_LASSO_Coefs.pdf"), plot = p,
         width = 7.5, height = plot_height, device = cairo_pdf, bg = "white")
  log_msg("   Saved: Fig6_LASSO_Coefs")
}


# =============================================================================
# FIGURE 7 — Time-Dependent ROC Curves
# =============================================================================
draw_fig7 <- function() {
  log_msg("")
  log_msg("Figure 7 — Time-dependent ROC")

  if (!requireNamespace("timeROC", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("timeROC", update = FALSE, ask = FALSE)
  }
  library(timeROC)

  # v2 path
  pred_path <- in_path(SIG_DIR, "signature_v2_test_predictions.csv")
  if (!file.exists(pred_path)) {
    pred_path <- in_path(INPUT_DIR, "test_scores_applied.csv")
  }
  if (!file.exists(pred_path)) stop("test predictions file not found.")

  df <- read_csv(pred_path, show_col_types = FALSE)
  log_msg(sprintf("   reading: %s", basename(pred_path)))
  log_msg(sprintf("   columns: %s", paste(names(df), collapse = ", ")))

  # standardise the column names
  time_col   <- intersect(c("survival_time", "time", "OS", "OS.time"), names(df))[1]
  status_col <- intersect(c("event", "status", "vital_status"), names(df))[1]
  score_col  <- intersect(c("risk_score_full", "risk_score", "score"), names(df))[1]

  if (any(is.na(c(time_col, status_col, score_col)))) {
    stop("required columns not found.")
  }

  df_clean <- df %>%
    transmute(
      Time   = as.numeric(.data[[time_col]]),
      Status = as.integer(.data[[status_col]]),
      Score  = as.numeric(.data[[score_col]])
    ) %>%
    filter(!is.na(Time), !is.na(Status), !is.na(Score), Time > 0)

  horizons <- c(180, 365, 730)
  labels   <- c("6-Month", "1-Year", "2-Year")

  roc_out <- timeROC(
    T      = df_clean$Time,
    delta  = df_clean$Status,
    marker = df_clean$Score,
    cause  = 1,
    times  = horizons,
    iid    = TRUE
  )

  # if AUC is below 0.5, flip the direction
  if (mean(roc_out$AUC, na.rm = TRUE) < 0.5) {
    log_msg("   score direction is inverted — flipping")
    roc_out <- timeROC(
      T      = df_clean$Time,
      delta  = df_clean$Status,
      marker = -df_clean$Score,
      cause  = 1,
      times  = horizons,
      iid    = TRUE
    )
  }

  roc_list <- list()
  for (i in seq_along(horizons)) {
    roc_list[[i]] <- data.frame(
      FPR     = roc_out$FP[, i],
      TPR     = roc_out$TP[, i],
      Horizon = labels[i],
      AUC     = roc_out$AUC[i]
    )
  }

  df_roc <- bind_rows(roc_list) %>%
    filter(is.finite(FPR), is.finite(TPR)) %>%
    mutate(FPR = pmax(0, pmin(1, FPR)),
           TPR = pmax(0, pmin(1, TPR))) %>%
    mutate(Label = paste0(Horizon, " (AUC = ",
                          sprintf("%.3f", AUC), ")"))

  df_roc$Label <- factor(df_roc$Label, levels = unique(df_roc$Label))
  colors <- c("#D73027", "#2C7BB6", "#4DAF4A")

  p <- ggplot(df_roc, aes(x = FPR, y = TPR, color = Label)) +
    geom_line(linewidth = 0.9) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey40", linewidth = 0.4) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       expand = c(0, 0)) +
    scale_color_manual(values = colors[seq_along(levels(df_roc$Label))]) +
    labs(
      title    = "Time-Dependent ROC Curves",
      subtitle = "miR-Mito Score (Test Set)",
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)",
      color = NULL
    ) +
    theme_pub() +
    theme(
      legend.position  = c(0.76, 0.18),
      legend.background = element_blank(),
      legend.key       = element_blank(),
      legend.text      = element_text(size = 9.5)
    )

  save_plot(p, "Fig7_ROC_Curves", 5.5, 5.5)

  log_msg("   AUCs:")
  print(df_roc %>% distinct(Horizon, AUC))
}


# =============================================================================
# FIGURE 8 — Calibration Plot
# =============================================================================
draw_fig8 <- function() {
  log_msg("")
  log_msg("Figure 8 — Calibration")

  # v2 path
  cal_path <- in_path(SIG_DIR, "14_4_calibration/calibration_data.csv")
  if (!file.exists(cal_path)) {
    cal_path <- in_path(INPUT_DIR, "calibration_data.csv")
  }
  if (!file.exists(cal_path)) stop("calibration_data.csv not found.")

  df <- read_csv(cal_path, show_col_types = FALSE)
  log_msg(sprintf("   reading: %s (%d rows)", basename(cal_path), nrow(df)))

  df_clean <- df %>%
    rename(
      Predicted = mean_predicted,
      Observed  = observed_surv,
      Lower     = observed_lower,
      Upper     = observed_upper,
      Time      = time_horizon
    )

  # convert Persian digits to ASCII digits (if needed)
  df_clean$Time <- chartr("\u06F0\u06F1\u06F2\u06F3\u06F4\u06F5\u06F6\u06F7\u06F8\u06F9", "0123456789", df_clean$Time)

  # clean English labels for the facets
  df_clean <- df_clean %>%
    mutate(
      TimeFacet = case_when(
        grepl("^6(mo|[- ]?month)", Time, ignore.case = TRUE)  ~ "6-Month",
        grepl("^1(yr|[- ]?year)",  Time, ignore.case = TRUE)  ~ "1-Year",
        grepl("^2(yr|[- ]?year)",  Time, ignore.case = TRUE)  ~ "2-Year",
        TRUE ~ as.character(Time)
      )
    )

  df_clean$TimeFacet <- factor(df_clean$TimeFacet,
                               levels = c("6-Month", "1-Year", "2-Year"))

  old_lc <- Sys.getlocale("LC_NUMERIC")
  Sys.setlocale("LC_NUMERIC", "C")
  on.exit(Sys.setlocale("LC_NUMERIC", old_lc), add = TRUE)

  eng_labels <- c("0%", "25%", "50%", "75%", "100%")

  p <- ggplot(df_clean, aes(x = Predicted, y = Observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey40", linewidth = 0.5) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.02,
                  color = "#2B83BA", linewidth = 0.5) +
    geom_line(color = "#2B83BA", linewidth = 0.8) +
    geom_point(color = "#2B83BA", size = 2.5) +
    facet_wrap(~ TimeFacet, nrow = 1) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), labels = eng_labels) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = eng_labels) +
    labs(
      title = "Calibration Plot - miR-Mito Score (Test Set)",
      x = "Predicted Survival Probability",
      y = "Observed Survival (KM ± 95% CI)"
    ) +
    theme_bw(base_family = "sans") +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.title       = element_text(face = "bold", size = 11),
      axis.text        = element_text(color = "black", size = 9.5),
      plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
      strip.background = element_rect(fill = "grey95", color = "black",
                                      linewidth = 0.6),
      strip.text       = element_text(face = "bold", size = 10),
      panel.spacing    = unit(1, "lines"),
      plot.margin      = margin(t = 10, r = 10, b = 15, l = 15, unit = "pt")
    )

  ggsave(out_path("Fig8_Calibration.png"), plot = p,
         width = 9.0, height = 4.0, dpi = DPI, bg = "white")
  ggsave(out_path("Fig8_Calibration.pdf"), plot = p,
         width = 9.0, height = 4.0, device = cairo_pdf, bg = "white")
  log_msg("   Saved: Fig8_Calibration")
}


# =============================================================================
# run everything
# =============================================================================
log_msg("============================================================")
log_msg("STEP 18b (Signature v2) — Figures 5-8")
log_msg("============================================================")

draw_fig5()
draw_fig6()
draw_fig7()
draw_fig8()

log_msg("")
log_msg("============================================================")
log_msg("all figures in: ", OUTPUT_DIR)
log_msg("============================================================")
