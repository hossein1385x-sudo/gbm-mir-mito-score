#!/usr/bin/env Rscript
# ============================================================
#  Sensitivity: consensus threshold 50%, 60%, 70%
# ============================================================

source(here::here("config.R"))  # defines WORK_DIR (see config.R)

freq <- read.csv(file.path(WORK_DIR, "signature_v2/signature_v2_selection_frequency.csv"))

for (thr in c(0.50, 0.60, 0.70)) {
  sel <- freq$miRNA[freq$frequency >= thr]
  cat(sprintf("Threshold = %.0f%%  ->  %d miRNAs selected: %s\n\n",
              thr*100, length(sel), paste(sel, collapse = ", ")))
}
out <- do.call(rbind, lapply(c(0.50, 0.60, 0.70), function(thr) {
  sel <- freq$miRNA[freq$frequency >= thr]
  data.frame(
    threshold       = sprintf("%.0f%%", thr*100),
    n_miRNAs        = length(sel),
    miRNAs_selected = paste(sel, collapse = "; "),
    includes_locked = all(c("hsa-miR-196a-5p","hsa-miR-675-5p","hsa-miR-4660-5p",
                            "hsa-miR-221-3p","hsa-miR-31-5p","hsa-miR-181a-2-3p") %in% sel)
  )
}))
print(out)
write.csv(out, file.path(WORK_DIR, "supplementary_final/TableS25_threshold_sensitivity.csv"),
          row.names = FALSE)
