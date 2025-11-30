#' run_rnaseqDegy
#'
#' Unified bulk RNA-seq pipeline for QC, DE (DESeq2/edgeR), ORA and GSEA.
#' Accepts file paths and parameters as R function arguments (wrapper around
#' the original CLI pipeline). All outputs are written into the chosen directory.
#'
#' @param counts Path to count matrix (TSV)
#' @param design Path to design matrix (TSV)
#' @param species "human" or "mouse"
#' @param de_method "DESeq2" or "edgeR"
#' @param condition_reference Reference condition
#' @param condition_test Test/contrast condition
#' @param filter_min_count Minimum count per sample (default 10)
#' @param filter_min_prop Minimum proportion of samples (default 0.5)
#' @param batch_correction Logical, PCA batch correction (TRUE/FALSE)
#' @param corr_method Correlation method: "pearson", "spearman", "kendall"
#' @param alpha FDR cutoff (default 0.05)
#' @param lfc_cutoff Log2FC threshold (default 1)
#' @param run_ora Logical, run ORA
#' @param run_gsea Logical, run GSEA
#' @param gmt_dir Directory containing GMT files
#' @param gmt_files Character vector of GMT filenames
#' @param gmt_terms Optional character vector of GO/GMT terms
#' @param output_dir Output directory
#' @param prefix Optional prefix for output files
#'
#' @return Invisibly returns the output directory path.
#' @export
run_rnaseqDegy <- function(
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
  lfc_cutoff = 1.0,
  run_ora = TRUE,
  run_gsea = TRUE,
  gmt_dir = NULL,
  gmt_files = NULL,
  gmt_terms = NULL,
  output_dir = "rnaseqDegy_output",
  prefix = NULL
) {
  message("===== rnaseqDegy :: Unified Bulk RNA-seq Pipeline =====")

  # Validate file paths
  if (!file.exists(counts)) stop("Count matrix not found: ", counts)
  if (!file.exists(design)) stop("Design file not found: ", design)

  if (!dir.exists(output_dir)) {
    message("Creating output directory: ", output_dir)
    dir.create(output_dir, recursive = TRUE)
  }

  # Convert gmt_files from vector to comma-separated string for compatibility
  gmt_files_str <- if (!is.null(gmt_files)) paste(gmt_files, collapse = ",") else ""
  gmt_terms_str <- if (!is.null(gmt_terms)) paste(gmt_terms, collapse = ",") else ""

  # Construct command-line-like argument list
  args <- list(
    counts = counts,
    design = design,
    species = species,
    de_method = de_method,
    condition_reference = condition_reference,
    condition_test = condition_test,
    filter_min_count = filter_min_count,
    filter_min_prop = filter_min_prop,
    batch_correction = batch_correction,
    corr_method = corr_method,
    alpha = alpha,
    lfc_cutoff = lfc_cutoff,
    run_ora = run_ora,
    run_gsea = run_gsea,
    gmt_dir = ifelse(is.null(gmt_dir), "", gmt_dir),
    gmt_files = gmt_files_str,
    gmt_terms = gmt_terms_str,
    output_dir = output_dir,
    prefix = ifelse(is.null(prefix), "", prefix)
  )

  message("Running rnaseqDegy pipeline...")

  # Source the internal engine
  # The engine will be created as inst/scripts/engine_internal.R
  internal_engine <- system.file("scripts", "engine_internal.R", package = "rnaseqDegy")

  if (!file.exists(internal_engine)) {
    stop("Internal engine 'engine_internal.R' not found in package inst/scripts/")
  }

  # Load the engine
  source(internal_engine, local = TRUE)

  # Call internal function .run_engine (defined inside engine_internal.R)
  .run_engine(args)

  message("rnaseqDegy completed successfully.")
  invisible(output_dir)
}
