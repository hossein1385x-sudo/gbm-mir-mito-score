# miR-Mito Score - miRNA-based prognostic signature for glioblastoma (TCGA-GBM)

Analysis code for the manuscript **"A miRNA Signature for Prognostic Stratification of Glioblastoma"**.
The pipeline builds the miR-Mito Score in TCGA-GBM and validates it externally in CGGA.

## Design principle

The patient-level **test set is split off and locked before any feature selection**; every selection,
scaling and tuning step uses the training set only (step 04 -> 08). The locked signature is applied to the
test set without refitting (step 09).

## Repository layout

```
config.R                  # the only place where paths are defined
install_packages.R        # installs all required packages
pipeline/                 # steps 00-20, run in numeric order
analysis/                 # robustness and sensitivity analyses
figures/                  # manuscript figures
supplementary/            # supplementary tables, data files and the supplementary .docx
```

## Setup

1. Install the packages: `Rscript install_packages.R`
2. Choose a data/results folder **outside** the repository and point to it (copy `.Renviron.example` to `.Renviron`):
   `GBM_WORK_DIR="D:/gbm_work"`. All inputs and outputs are read from / written to this folder.
3. Run scripts from anywhere inside the repository folder, e.g. `Rscript pipeline/07_signature_lasso_cox.R`.

### External inputs (not included; place them in `GBM_WORK_DIR`)

| Resource | Expected location |
|---|---|
| TCGA-GBM RNA-seq and miRNA-seq (GDC, via `pipeline/00_download_gdc.ps1`) | `rnaseq/`, `mirna/` |
| miRBase `mature.fa` | `mature.fa.txt` |
| MitoCarta 3.0 (Excel file whose name contains `MitoCarta`) | `GBM_WORK_DIR/` (read by step 03b) |
| miRTarBase, miRDB v6.0 | `reference/mirtarbase_MTI.csv`, `reference/miRDB_v6.0_prediction_result.txt` |
| CGGA microRNA array (198 samples): expression and clinical | `cgga_raw/exp cgga/...`, `cgga_raw/clin cgga/...` |

## Pipeline

| Script | Purpose |
|---|---|
| `00_download_gdc.ps1` | Resumable download of GDC manifest parts (RNA-seq / miRNA-seq) with gdc-client |
| `01_parse_gdc_downloads.R` | Parse raw downloads, attach barcodes and clinical data, assert every file is TCGA-GBM |
| `02_isoform_quantification.R` | miRNA isoform-level quantification (MIMAT -> name mapping from miRBase) |
| `03b_mito_gene_list_rebuild.R` | Build the mitochondrial gene list from MitoCarta 3.0 (symbol -> Ensembl ID) and clean the RNA-seq matrix |
| `04_train_test_split.R` | Patient-level train/test split; test set locked |
| `05_mds_construction.R` | Mitochondrial score (MDS): PC1 of z-scored mito genes, train only |
| `06_differential_expression.R` | DE of miRNA isoforms, MDS-high vs MDS-low (train only) |
| `07_signature_lasso_cox.R` | LASSO-Cox with nested CV and final model (6 miRNA + Age, signature v2) |
| `08_signature_lock.R` | Freeze the signature (SHA-256 hash and manifest) |
| `09_apply_to_test.R` | Apply the locked signature to the test set, without refitting |
| `10_performance_evaluation.R` | Test-set performance: C-index, HR, time-dependent AUC |
| `11_kaplan_meier.R` | Kaplan-Meier stratification and spline check |
| `12_bootstrap_optimism.R` | Bootstrap optimism correction (full pipeline rebuilt per replicate) |
| `13_verhaak_subtype.R` | Association with Verhaak mRNA subtype; incremental value; interaction |
| `14_calibration.R` | Calibration at 6 months / 1 year / 2 years |
| `15_brier_score.R` | Brier score (corrected, Graf formula) |
| `16_target_gene_enrichment.R` | GO / KEGG enrichment of miRNA target genes |
| `17_cytoscape_inputs.R` | Node / edge tables for Cytoscape |
| `18_decision_curve_analysis.R` | Decision curve analysis |
| `19_multiseed_stability.R` | Stability across 50 repeated 80/20 splits |
| `20_cgga_external_validation.R` | External validation in CGGA |

The scripts in `analysis/`, `figures/` and `supplementary/` depend on the outputs of the pipeline; their
header comments list inputs and outputs.

## Reproducibility

- Every script reads and writes through `config.R`, so no path needs to be edited inside the code.
- Random seeds are fixed with `set.seed()` inside each script that uses randomness (train/test split, nested CV, bootstrap, multi-seed stability).
- The signature is frozen with a SHA-256 hash and manifest (step 08) before it touches the test set.

## Citation

If you use this code, please cite the associated manuscript:
*A miRNA Signature for Prognostic Stratification of Glioblastoma* (citation details will be added on publication).

## License

Released under the [MIT License](LICENSE).
