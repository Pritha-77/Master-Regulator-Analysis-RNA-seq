# Master Regulator Analysis of Adipogenic Differentiation (Day0–Day8)

This repository contains an end-to-end **Master Regulator Analysis (MRA)** pipeline for bulk RNA-seq data across adipogenic differentiation time points (Day0, Day4, Day8). The workflow integrates **DESeq2**, **RegEnrich**, and **RTN** to robustly identify transcription factors (TFs) that drive late-stage adipogenesis, with a focus on high-confidence regulators supported by multiple methods.

---

## Overview of the Workflow

The analysis follows these major steps:

1. **Input preprocessing**

   * Import raw gene count matrix
   * Clean gene symbols (remove duplicates, NAs)
   * Construct sample metadata

2. **Filtering & normalization**

   * Filter lowly expressed genes (≥80% samples)
   * Variance stabilizing transformation (DESeq2)
   * TMM normalization + voom (edgeR/limma)

3. **Differential expression analysis**

   * Day4 vs Day0 (early differentiation)
   * Day8 vs Day0 (late/terminal differentiation)

4. **Master Regulator Analysis (three approaches)**

   * **RegEnrich (COEN + GSEA)**
   * **RTN + Fisher’s Exact Test (FET)**
   * **RTN + GSEA-based enrichment**

5. **Integration & visualization**

   * Overlap of master regulators across methods
   * Venn diagram of shared TFs
   * Identification of high-confidence master regulators

---

## Input Data

### Required files

* `raw.genecount.xlsx`

  * Gene-level raw counts
  * First column: **Gene symbols**
  * Remaining columns: samples

### Sample design

The script assumes **9 samples**:

| Condition | Replicates |
| --------- | ---------- |
| Day0      | 3          |
| Day4      | 3          |
| Day8      | 3          |

If your experimental design differs, update the `conditions` vector in the script accordingly.

---

## Software Requirements

### R version

* R ≥ 4.2 recommended

### Required R packages

```r
RegEnrich
RTN
DESeq2
limma
edgeR
org.Hs.eg.db
RedeR
igraph
clusterProfiler
pheatmap
ggplot2
dplyr
tidyr
data.table
VennDiagram
openxlsx
grid
```

Install missing packages before running the script.

---


## Output Files

| File                                             | Description                       |
| ------------------------------------------------ | --------------------------------- |
| `sample_info.txt`                                | Sample metadata                   |
| `res_mra_regenrich.csv`                          | RegEnrich master regulator scores |
| `RTN_MRA_Day8_vs_Day0_FET.csv`                   | RTN MRA (Fisher’s Exact Test)     |
| `RTN_MRA_GSEA_Day8_vs_Day0.csv`                  | RTN MRA (GSEA-based)              |
| `Venn_Master_Regulators_Day8_vs_Day0.png`        | Overlap of MRs across methods     |
| `Common_Master_Regulators_All3_Day8_vs_Day0.txt` | High-confidence shared MRs        |

---

## Interpretation

* **RegEnrich** identifies regulators based on co-expression networks and enrichment of differential signals.
* **RTN-FET** highlights TFs whose regulons are overrepresented among significantly DE genes.
* **RTN-GSEA** detects TFs whose regulons show coordinated expression shifts.

Transcription factors shared across all three methods represent **high-confidence master regulators** of late-stage adipogenesis.

---

## Notes & Customization

* Adjust filtering thresholds (`padj`, `log2FoldChange`, regulon size) as needed.
* Increase RTN permutations (`nPermutations`) for higher confidence at the cost of runtime.
* TF list can be replaced with a custom, tissue-specific TF annotation.


