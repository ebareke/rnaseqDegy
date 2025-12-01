# Changelog

All notable changes to **rnaseqDegy** will be documented in this file.

The format follows **Keep a Changelog**
(<https://keepachangelog.com/en/1.1.0/>)  
and the project adheres to **Semantic Versioning**
(<https://semver.org/>).

------------------------------------------------------------------------

## \[0.1.0\] - 2025-02-01

### Added

- Initial release of **rnaseqDegy**, a unified bulk RNA-seq pipeline R
  package.
- End-to-end workflow covering:
  - Count filtering (≥10 counts in ≥50% samples).
  - QC: PCA (with/without batch correction), sample correlation heatmap,
    top 500 variable-gene heatmap.
  - Differential expression (DESeq2 or edgeR) with batch correction
    integrated in the design.
  - Enhanced volcano plots with gene counts and top 10 labels.
  - ORA (GO BP/MF/CC) via clusterProfiler (dotplot, barplot, cnetplot,
    emapplot).
  - GSEA using user-supplied offline GMT files (dotplot, barplot,
    cnetplot, emapplot).
- Added
  [`run_rnaseqDegy()`](https://ebareke.github.io/rnaseqDegy/reference/run_rnaseqDegy.md)
  high-level wrapper function for R workflow integration.
- Internal modular engine located under
  `inst/scripts/engine_internal.R`.
- Added comprehensive documentation and roxygen tags.
- Added metadata files:
  - `README.md`
  - `LICENSE` (MIT)
  - `CITATION.cff`
  - `CONTRIBUTING.md`
  - `NAMESPACE` and `DESCRIPTION`

### Notes

- This is the first stable functional version of the package.
- Future versions will expand enrichment options (KEGG, Reactome),
  introduce multi-contrast support, and add automated HTML reporting.

------------------------------------------------------------------------

## \[Unreleased\]

### Planned

- Support for KEGG/Reactome ORA + GSEA.
- Optional batch-aware normalization choices (e.g., RUVseq, ComBat‑Seq).
- Multi-contrast DE support.
- Automated QC/DE report (HTML via rmarkdown).
- Extending species support (rat, drosophila).
- Shiny app wrapper for interactive use.
