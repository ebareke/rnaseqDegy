# Advanced rnaseqDegy Workflows

## Overview

This vignette shows **advanced usage patterns** of `rnaseqDegy`,
including:

- Custom filtering thresholds
- Running the pipeline with **edgeR** instead of DESeq2
- Enabling **GO ORA** and **GSEA** with user-supplied GMTs
- Running multiple contrasts
- Tips for HPC and scripting

We assume you already read the **Introduction** vignette.

## Setup

``` r
library(rnaseqDegy)

counts <- system.file("extdata", "example_counts.tsv", package = "rnaseqDegy")
design <- system.file("extdata", "example_design.tsv", package = "rnaseqDegy")
```

## 1. Custom Filtering and Correlation Method

You can relax or tighten the expression filtering criteria and change
the correlation method:

``` r
run_rnaseqDegy(
  counts = counts,
  design = design,
  species = "human",
  de_method = "DESeq2",
  condition_reference = "Control",
  condition_test = "Mutant",
  filter_min_count = 5,     # instead of default 10
  filter_min_prop  = 0.4,   # instead of default 0.5
  corr_method      = "spearman",
  batch_correction = TRUE,
  run_ora = FALSE,
  run_gsea = FALSE,
  output_dir = "rnaseqDegy_advanced_filter"
)
```

## 2. Using edgeR as the DE Engine

To use **edgeR** for differential expression:

``` r
run_rnaseqDegy(
  counts = counts,
  design = design,
  species = "human",
  de_method = "edgeR",
  condition_reference = "Control",
  condition_test = "Mutant",
  batch_correction = TRUE,
  alpha = 0.01,      # stricter FDR cutoff
  lfc_cutoff = 1.5,  # higher fold-change threshold
  run_ora = FALSE,
  run_gsea = FALSE,
  output_dir = "rnaseqDegy_edgeR"
)
```

`rnaseqDegy` will automatically:

- Fit a GLM with `~ batch + group`
- Compute QL F-tests
- Harmonize output columns to `log2FoldChange` and `padj`

## 3. Running GO ORA (BP/MF/CC)

To enable GO ORA, simply set `run_ora = TRUE`:

``` r
run_rnaseqDegy(
  counts = counts,
  design = design,
  species = "human",
  de_method = "DESeq2",
  condition_reference = "Control",
  condition_test = "Mutant",
  run_ora = TRUE,
  run_gsea = FALSE,
  output_dir = "rnaseqDegy_ORA"
)
```

Internally, the pipeline:

- Uses all expressed genes as **background** (universe)
- Selects significant DEGs using `alpha` and `lfc_cutoff`
- Maps symbols to **ENTREZ IDs**
- Runs
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html)
  for **BP**, **MF**, and **CC**
- Simplifies terms by semantic similarity (`simplify()`)

The following files are created:

- `ORA_GO_BP.tsv`, `ORA_GO_MF.tsv`, `ORA_GO_CC.tsv`
- `ORA_GO_<ONT>_dotplot.pdf/png`
- `ORA_GO_<ONT>_barplot.pdf/png`
- `ORA_GO_<ONT>_cnetplot.pdf/png`
- `ORA_GO_<ONT>_emapplot.pdf/png`

## 4. Running GSEA with Offline GMT Files

To use GSEA you must provide:

- `gmt_dir`: directory containing GMT files
- `gmt_files`: character vector of GMT filenames

Example GMT layout:

``` text
gmt/
 ├── HALLMARK.v2023.1.Hs.symbols.gmt
 └── C5.GO.BP.v2023.1.Hs.symbols.gmt
```

Example call:

``` r
run_rnaseqDegy(
  counts = counts,
  design = design,
  species = "human",
  de_method = "DESeq2",
  condition_reference = "Control",
  condition_test = "Mutant",
  run_ora = FALSE,
  run_gsea = TRUE,
  gmt_dir = "gmt",
  gmt_files = c(
    "HALLMARK.v2023.1.Hs.symbols.gmt",
    "C5.GO.BP.v2023.1.Hs.symbols.gmt"
  ),
  output_dir = "rnaseqDegy_GSEA"
)
```

For each GMT file, the pipeline:

- Ranks all genes by log2FC
- Runs
  [`clusterProfiler::GSEA()`](https://rdrr.io/pkg/clusterProfiler/man/GSEA.html)
- Writes `GSEA_<gmt_name>.tsv`
- Creates dotplot, barplot, cnetplot, emapplot

### Focusing on specific terms

Optionally, you can restrict plots to a subset of terms using
`gmt_terms`:

``` r
run_rnaseqDegy(
  counts = counts,
  design = design,
  species = "human",
  de_method = "DESeq2",
  condition_reference = "Control",
  condition_test = "Mutant",
  run_ora = TRUE,
  run_gsea = TRUE,
  gmt_dir = "gmt",
  gmt_files = c("HALLMARK.v2023.1.Hs.symbols.gmt"),
  gmt_terms = c("HALLMARK_APOPTOSIS", "HALLMARK_INFLAMMATORY_RESPONSE"),
  output_dir = "rnaseqDegy_focus_terms"
)
```

## 5. Multiple Contrasts in a Script

`rnaseqDegy` focuses on **one contrast per call**, but you can script
multiple contrasts:

``` r
contrasts <- list(
  list(ref = "Control", test = "Mutant"),
  list(ref = "Control", test = "Treatment2")
)

for (c in contrasts) {
  outdir <- sprintf("rnaseqDegy_%s_vs_%s", c$ref, c$test)

  run_rnaseqDegy(
    counts = counts,
    design = design,
    species = "human",
    de_method = "DESeq2",
    condition_reference = c$ref,
    condition_test = c$test,
    batch_correction = TRUE,
    run_ora = TRUE,
    run_gsea = FALSE,
    output_dir = outdir,
    prefix = sprintf("%s_vs_%s", c$ref, c$test)
  )
}
```

This pattern is particularly useful on **HPC clusters**, where each
contrast can be launched as a separate job.

## 6. Tips for HPC and Automation

- Use absolute paths for `counts`, `design`, `gmt_dir`, and
  `output_dir`.
- Keep `output_dir` unique for each job/contrast.
- You can wrap
  [`run_rnaseqDegy()`](https://ebareke.github.io/rnaseqDegy/reference/run_rnaseqDegy.md)
  inside a small R script and call it with `Rscript` from a scheduler
  (SLURM, PBS, etc.).

Example SLURM command (shell-level):

``` bash
Rscript -e "rnaseqDegy::run_rnaseqDegy(\
  counts = 'counts.tsv', \
  design = 'design.tsv', \
  species = 'human', \
  de_method = 'DESeq2', \
  condition_reference = 'Control', \
  condition_test = 'Mutant', \
  batch_correction = TRUE, \
  run_ora = TRUE, \
  run_gsea = FALSE, \
  output_dir = 'rnaseqDegy_job1')"
```

## 7. Troubleshooting

- **Installation issues**: Ensure that all Bioconductor dependencies
  (`DESeq2`, `edgeR`, `clusterProfiler`, etc.) are installed with
  [`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html).
- **No significant DEGs**: Try relaxing `alpha` or `lfc_cutoff`, or
  verify your design and sample size.
- **ORA/GSEA empty**: Ensure enough DEGs for ORA (e.g., ≥ 10) and that
  GMT identifiers match your gene identifiers (symbols).

## Conclusion

This vignette illustrated advanced control over filtering, DE engines,
enrichment options, and multi-contrast workflows. Combined with the
introductory vignette, you should now be able to integrate `rnaseqDegy`
into both interactive analyses and production pipelines.

For questions and feature requests, please open an issue at:
<https://github.com/ebareke/rnaseqDegy/issues>
