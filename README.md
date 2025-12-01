# ***rnaseqDegy*** <img src="https://img.shields.io/badge/status-active-success" align="right"/>

### 🚀 Unified Bulk RNA‑seq QC, Differential Expression, and Enrichment Pipeline
**Author:** Dr. Eric Bareke (Majewski Lab, McGill University)  
**Repository:** https://github.com/ebareke/rnaseqDegy

---

`rnaseqDegy` is an R package providing a **single, fully automated pipeline** for bulk RNA‑seq analysis — covering **QC**, **differential expression**, **GO ORA**, and **GSEA** using **user‑supplied GMT files**, for both **human and mouse** datasets.  

The package is designed for **research labs**, **HPC clusters**, and **automated workflows**, offering clean, reproducible, publication‑ready outputs.

---

## ✨ Key Features

### 🔍 QC & Normalization
- Count filtering (≥10 counts in ≥50% samples)
- PCA (with optional batch correction via `limma::removeBatchEffect`)
- Sample correlation heatmap (Pearson/Spearman/Kendall)
- Top 500 most variable gene heatmap
- Auto‑dimensioned plots with **publication‑ready styling**

### ⚔ Differential Expression
- Choose **DESeq2** or **edgeR**
- Batch included in the design (`~ Batch + Condition`)
- Export full results + normalized counts

### 🌋 Volcano Plots
- Highlights top 10 up/down DEGs
- Shows counts of Up/Down genes in subtitle
- Custom colour scheme and bold theme

### 🧬 Functional Enrichment
#### ORA (GO BP/MF/CC)
- Uses **clusterProfiler**
- Simplification by semantic similarity
- Dotplot, barplot, cnetplot, emapplot

#### GSEA (offline, user GMT)
- Accepts any custom GMT file
- Dotplot, barplot, cnetplot, emapplot
- Auto‑filtered term subset via `gmt_terms`

---

## 📦 Installation

This package requires **Bioconductor dependencies**. Install everything cleanly with:

```r
# Install devtools if not present
install.packages("devtools")

# Install rnaseqDegy from GitHub
devtools::install_github("ebareke/rnaseqDegy")
```

---

## 🧠 Quick Start

````r
library(rnaseqDegy)

run_rnaseqDegy(
  counts = "counts.tsv",
  design = "design.tsv",
  species = "human",
  de_method = "DESeq2",
  condition_reference = "Control",
  condition_test = "Mutant",
  batch_correction = TRUE,
  run_ora = TRUE,
  run_gsea = TRUE,
  gmt_dir = "gmt/",
  gmt_files = c("HALLMARK.gmt", "GO_BP.gmt"),
  output_dir = "results/Control_vs_Mutant",
  prefix = "C_vs_M"
)
````

This will generate:
- QC plots
- DE tables + volcano
- ORA (BP/MF/CC)
- GSEA (for each GMT)

Outputs are auto‑named and saved as **PDF + PNG**.

---

## 📁 Input File Format

### 1. **Count Matrix** (TSV)
```
GeneID    S1   S2   S3   S4
ENSG0001  50   23   11   9
ENSG0002  140  180  90   76
...
```
- First column: **GeneID** (Ensembl or Symbol)
- Following columns: raw integer counts

### 2. **Design File** (TSV)
```
Sample   Batch   Condition
S1       B1      Control
S2       B1      Control
S3       B2      Mutant
S4       B2      Mutant
```

---

## 🧬 Output Structure
```
results/
 ├── QC_PCA_*.pdf/png
 ├── QC_sampleCorrelation.pdf/png
 ├── QC_topVar500_heatmap.pdf/png
 ├── DE_results_all.tsv
 ├── normalized_counts.tsv
 ├── DE_volcano.pdf/png
 ├── ORA_GO_BP_*.tsv + plots
 ├── ORA_GO_MF_*.tsv + plots
 ├── ORA_GO_CC_*.tsv + plots
 ├── GSEA_<gmt>_*.tsv + plots
```

---

## 🧪 Example Folder Layout
```
project/
 ├── counts.tsv
 ├── design.tsv
 ├── gmt/
 │    ├── HALLMARK.gmt
 │    └── GO_BP.gmt
 └── rnaseqDegy_results/
```

---

## ⚙ Advanced Options
- `filter_min_count`, `filter_min_prop`
- `batch_correction = TRUE/FALSE`
- `gmt_terms = c("GOBP_APOPTOSIS","MYC_TARGETS_V2")`
- `corr_method = "spearman"`
- `output_dir`, `prefix`

All parameters mirror the original CLI pipeline.

---

## 📘 Citation
If you use **rnaseqDegy** in your research, please cite:
```
Dr. Eric Bareke, Majewski Lab, McGill University.
```
(A full CITATION.cff file is included in this repository.)

---

## 🤝 Contributing
See **CONTRIBUTING.md** for guidelines. Pull requests are welcome.

---

## 📝 License
MIT License — see `LICENSE` file.
