#!/usr/bin/env Rscript
# ============================================================
# Combined Robustness Forest Plot for miR-Mito score in GBM
# Layout: HORIZONTAL  (Panel A | Panel B)
# Output: 600 DPI TIFF + PNG
# ============================================================


source(here::here("config.R"))  # defines WORK_DIR (see config.R)

library(ggplot2)
library(patchwork)
library(dplyr)
library(scales)

# ------------------------------------------------------------
# 1) Panel A data  →  Hazard Ratio (HR) with 95% CI
# ------------------------------------------------------------
df_hr <- data.frame(
  Model = c(
    "TCGA-GBM: Age + Sex",
    "TCGA-GBM: + Stage",
    "TCGA-GBM: + MGMT",
    "TCGA-GBM: + IDH",
    "TCGA-GBM: Fully adjusted",
    "CGGA (external validation)"
  ),
  HR    = c(1.92, 1.78, 1.65, 1.58, 1.52, 1.47),
  Lower = c(1.45, 1.32, 1.22, 1.17, 1.12, 1.08),
  Upper = c(2.54, 2.40, 2.23, 2.13, 2.06, 2.00),
  Group = c(rep("TCGA-GBM (adjusted)", 5), "CGGA (external)")
)

df_hr$Model <- factor(df_hr$Model, levels = rev(df_hr$Model))

# ------------------------------------------------------------
# 2) Panel B data  →  C-index with 95% CI
# ------------------------------------------------------------
df_c <- data.frame(
  Metric = c(
    "Apparent (TCGA)",
    "Bootstrap-corrected (TCGA)",
    "CGGA external validation"
  ),
  Cindex = c(0.74, 0.71, 0.68),
  Lower  = c(0.70, 0.67, 0.63),
  Upper  = c(0.78, 0.75, 0.73),
  Group  = c("TCGA-GBM", "TCGA-GBM", "CGGA")
)

df_c$Metric <- factor(df_c$Metric, levels = rev(df_c$Metric))

# ------------------------------------------------------------
# 3) Panel A: hazard-ratio forest plot
# ------------------------------------------------------------
panel_A <- ggplot(df_hr, aes(x = HR, y = Model, color = Group)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.18, linewidth = 0.7) +
  geom_point(size = 3.2, shape = 18) +
  geom_text(aes(label = sprintf("%.2f (%.2f–%.2f)", HR, Lower, Upper)),
            x = max(df_hr$Upper) * 1.25,
            hjust = 0, size = 3.0, color = "black") +
  scale_color_manual(values = c("TCGA-GBM (adjusted)" = "#1f77b4",
                                "CGGA (external)"      = "#d62728")) +
  scale_x_continuous(
    trans  = "log10",
    breaks = c(0.5, 1, 1.5, 2, 3),
    labels = c("0.5", "1.0", "1.5", "2.0", "3.0"),
    limits = c(0.8, max(df_hr$Upper) * 2.2)
  ) +
  labs(
    title = "A.  Hazard-ratio robustness",
    x = "Hazard Ratio (95% CI, log scale)",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position   = "bottom",
    legend.title      = element_blank(),
    plot.title        = element_text(face = "bold", size = 12),
    axis.text.y       = element_text(size = 9.5)
  )

# ------------------------------------------------------------
# 4) Panel B: C-index forest plot
# ------------------------------------------------------------
panel_B <- ggplot(df_c, aes(x = Cindex, y = Metric, color = Group)) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.15, linewidth = 0.7) +
  geom_point(size = 3.2, shape = 18) +
  geom_text(aes(label = sprintf("%.2f (%.2f–%.2f)", Cindex, Lower, Upper)),
            x = max(df_c$Upper) * 1.05,
            hjust = 0, size = 3.0, color = "black") +
  scale_color_manual(values = c("TCGA-GBM" = "#2ca02c",
                                "CGGA"     = "#d62728")) +
  scale_x_continuous(
    breaks = seq(0.5, 0.9, 0.05),
    limits = c(0.5, max(df_c$Upper) * 1.9)
  ) +
  labs(
    title = "B.  Discrimination robustness (C-index)",
    x = "C-index (95% CI)",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position   = "bottom",
    legend.title      = element_blank(),
    plot.title        = element_text(face = "bold", size = 12),
    axis.text.y       = element_text(size = 9.5)
  )

# ------------------------------------------------------------
# 5) Combine the panels horizontally  →  A | B  (side by side)
#    note: the | operator in patchwork means a horizontal layout
# ------------------------------------------------------------
combined <- (panel_A | panel_B) +
  plot_annotation(
    title    = "Figure 9. Combined robustness forest plot for the miR-Mito score in GBM",
    subtitle = "Left: adjusted TCGA-GBM models and CGGA external validation   |   Right: apparent, bootstrap-corrected, and CGGA validation",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0)
    )
  ) &
  theme(plot.margin = margin(6, 12, 6, 6))

# ------------------------------------------------------------
# 6) Save at DPI = 600  (wider, shorter)
# ------------------------------------------------------------
out_dir <- file.path(WORK_DIR, "signature_v2")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(
  filename = file.path(out_dir, "Figure9_combined_forest_600dpi.tiff"),
  plot     = combined,
  width    = 14, height = 6, units = "in",   # ← horizontal layout: larger width, smaller height
  dpi      = 600, compression = "lzw"
)

ggsave(
  filename = file.path(out_dir, "Figure9_combined_forest_600dpi.png"),
  plot     = combined,
  width    = 14, height = 6, units = "in",
  dpi      = 600, bg = "white"
)

message("Figure 9 (horizontal layout) saved at 600 DPI ")
