#!/usr/bin/env Rscript
# =============================================================================
# Combine the median-split and quartile-split Kaplan-Meier panels into one 600-dpi image
# =============================================================================

source(here::here("config.R"))  # defines WORK_DIR (see config.R)

library(magick)

# image paths
img1_path <- file.path(WORK_DIR, "signature_v2/signature_v2_km_median_600dpi.png")
img2_path <- file.path(WORK_DIR, "signature_v2/signature_v2_km_quartile_600dpi.png")
out_path  <- file.path(WORK_DIR, "signature_v2/signature_v2_combined_600dpi.png")

# read the images
img1 <- image_read(img1_path)
img2 <- image_read(img2_path)

# paste side by side (horizontally)
combined <- image_append(c(img1, img2), stack = FALSE)

# save at 600 DPI
image_write(
  combined,
  path = out_path,
  format = "png",
  density = "600x600"   # ← set the DPI when saving
)

# check the final DPI
info <- image_info(combined)
print(info$density)
