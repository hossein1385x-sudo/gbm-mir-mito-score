#!/usr/bin/env Rscript
# =============================================================================
# PubMed hit counts for each signature miRNA (supports Table S21)
# =============================================================================

source(here::here("config.R"))  # defines WORK_DIR (see config.R)

if (!requireNamespace("rentrez", quietly = TRUE))
  install.packages("rentrez")
library(rentrez)

mirnas <- c("hsa-miR-196a-5p" = "miR-196a",
            "hsa-miR-675-5p"  = "miR-675",
            "hsa-miR-4660-5p" = "miR-4660",
            "hsa-miR-221-3p"  = "miR-221",
            "hsa-miR-31-5p"   = "miR-31",
            "hsa-miR-181a-2-3p" = "miR-181a")

results <- list()
for (name in names(mirnas)) {
  q <- paste0('(', mirnas[name], '[Title/Abstract]) AND ',
              '(glioblastoma[Title/Abstract] OR glioma[Title/Abstract]) AND ',
              '(prognosis[Title/Abstract] OR survival[Title/Abstract] OR expression[Title/Abstract])')
  hits <- entrez_search(db = "pubmed", term = q, retmax = 5)

  if (length(hits$ids) > 0) {
    summ <- entrez_summary(db = "pubmed", id = hits$ids)
    results[[name]] <- data.frame(
      miRNA = name,
      PMID  = hits$ids,
      Title = sapply(summ, function(s) s$title),
      Year  = sapply(summ, function(s) s$pubdate),
      stringsAsFactors = FALSE)
  }
}

df_s21 <- do.call(rbind, results)
write.csv(df_s21,
          file.path(WORK_DIR, "supplementary_final/TableS21_pubmed_hits.csv"),
          row.names = FALSE)
print(df_s21, row.names = FALSE)
