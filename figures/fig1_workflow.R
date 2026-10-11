#!/usr/bin/env Rscript
# ============================================================
#  Figure 1 workflow — 600 DPI, single panel
# ============================================================
rm(list = ls())

source(here::here("config.R"))  # defines WORK_DIR (see config.R)

OUT <- file.path(WORK_DIR, "figures_final")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

# 4-step workflow boxes
make_box <- function(title, subtitle, x, y, w = 0.22, h = 0.55, fill = "#E8F0FE") {
  data.frame(x = x, y = y, w = w, h = h, title = title, subtitle = subtitle, fill = fill)
}

boxes <- rbind(
  make_box("TCGA-GBM Patients",
           "MitoCarta 3.0\n1,057 genes\nMDS = PC1",
           0.12, 0.5),
  make_box("6-miRNA + Age\nSignature",
           "Elastic Net-Cox\n(alpha = 0.5)\nLocked before test",
           0.37, 0.5),
  make_box("Overall Survival",
           "HR = 1.64 / SD\n(TCGA test)\nHR = 2.22 / SD\n(CGGA external)",
           0.62, 0.5),
  make_box("Target Enrichment",
           "Hypoxia-response\nGliogenesis\nNeurogenesis\nCell-matrix adhesion",
           0.87, 0.5, fill = "#FCE8E6")
)

p <- ggplot() +
  # boxes
  geom_rect(data = boxes,
            aes(xmin = x - w/2, xmax = x + w/2,
                ymin = y - h/2, ymax = y + h/2,
                fill = fill),
            colour = "black", linewidth = 0.6) +
  scale_fill_identity() +
  # titles
  geom_text(data = boxes,
            aes(x = x, y = y + h/2 - 0.09, label = title),
            fontface = "bold", size = 3.8, lineheight = 0.95) +
  # subtitles
  geom_text(data = boxes,
            aes(x = x, y = y - 0.06, label = subtitle),
            size = 3.0, colour = "grey20", lineheight = 1.1) +
  # arrows
  annotate("segment", x = 0.235, xend = 0.255, y = 0.5, yend = 0.5,
           arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.6) +
  annotate("segment", x = 0.485, xend = 0.505, y = 0.5, yend = 0.5,
           arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.6) +
  annotate("segment", x = 0.735, xend = 0.755, y = 0.5, yend = 0.5,
           arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.6) +
  # panel
  coord_fixed(ratio = 0.55, xlim = c(0, 1), ylim = c(0.1, 0.9), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(15, 15, 15, 15))

ggsave(file.path(OUT, "Fig1_Workflow_600dpi.png"), p,
       width = 10, height = 3.2, dpi = 600, bg = "white")
ggsave(file.path(OUT, "Fig1_Workflow_600dpi.pdf"), p,
       width = 10, height = 3.2, bg = "white")

cat("Saved to:", OUT, "\n")
