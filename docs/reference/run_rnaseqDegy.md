# run_rnaseqDegy

Unified bulk RNA-seq pipeline for QC, DE (DESeq2/edgeR), ORA and GSEA.
Accepts file paths and parameters as R function arguments (wrapper
around the original CLI pipeline). All outputs are written into the
chosen directory.

## Usage

``` r
run_rnaseqDegy(
  counts,
  design,
  species = "human",
  de_method = "DESeq2",
  condition_reference,
  condition_test,
  filter_min_count = 10,
  filter_min_prop = 0.5,
  batch_correction = TRUE,
  corr_method = "pearson",
  alpha = 0.05,
  lfc_cutoff = 1,
  run_ora = TRUE,
  run_gsea = TRUE,
  gmt_dir = NULL,
  gmt_files = NULL,
  gmt_terms = NULL,
  output_dir = "rnaseqDegy_output",
  prefix = NULL
)
```

## Arguments

- counts:

  Path to count matrix (TSV)

- design:

  Path to design matrix (TSV)

- species:

  "human" or "mouse"

- de_method:

  "DESeq2" or "edgeR"

- condition_reference:

  Reference condition

- condition_test:

  Test/contrast condition

- filter_min_count:

  Minimum count per sample (default 10)

- filter_min_prop:

  Minimum proportion of samples (default 0.5)

- batch_correction:

  Logical, PCA batch correction (TRUE/FALSE)

- corr_method:

  Correlation method: "pearson", "spearman", "kendall"

- alpha:

  FDR cutoff (default 0.05)

- lfc_cutoff:

  Log2FC threshold (default 1)

- run_ora:

  Logical, run ORA

- run_gsea:

  Logical, run GSEA

- gmt_dir:

  Directory containing GMT files

- gmt_files:

  Character vector of GMT filenames

- gmt_terms:

  Optional character vector of GO/GMT terms

- output_dir:

  Output directory

- prefix:

  Optional prefix for output files

## Value

Invisibly returns the output directory path.
