# Installs every package used by the scripts (not run automatically).
# Optional: step 03b has a fallback that uses the 'convertid' package; install it only if you need that path.
cran <- c("data.table",
         "jsonlite",
         "digest",
         "httr",
         "stringr",
         "dplyr",
         "tidyr",
         "purrr",
         "readr",
         "readxl",
         "forcats",
         "scales",
         "ggplot2",
         "ggrepel",
         "patchwork",
         "gridExtra",
         "magick",
         "survival",
         "survminer",
         "glmnet",
         "timeROC",
         "dcurves",
         "matrixStats",
         "broom",
         "officer",
         "flextable",
         "rentrez",
         "here")
bioc <- c("DESeq2",
         "edgeR",
         "limma",
         "GSVA",
         "biomaRt",
         "org.Hs.eg.db",
         "AnnotationDbi",
         "clusterProfiler",
         "TCGAbiolinks")

install.packages(setdiff(cran, rownames(installed.packages())))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE, ask = FALSE)
