# How to run it ?
# Rscript GSEA.r -f /path/to/your/input/files -c VEH,T_DM1 -r VEH -t 2 -p 0.05 -m 10
set.seed(54321)
# List of required packages
libraries <- c(
  "DESeq2",
  "org.Hs.eg.db",
  "tidyverse",
  "clusterProfiler",
#  "GseaVis",
  "RColorBrewer",
  "patchwork",
 # "msigdbr",
  "ggplot2",
  "ggrepel",
  "reshape2",
  "argparse",
  "here",
  "data.table",
  "grDevices"
)

# Function to check, install, and load packages
check_install_load <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    # Install the package if not installed
#    if (package %in% BiocManager::installed()) {
#      BiocManager::install(package)
#   } else {
#      install.packages(package)
#    }
  }

  # Load the package
  suppressPackageStartupMessages(library(package, character.only = TRUE))
}

# Apply the function to each package in the list
lapply(libraries, check_install_load)


# Function to parse command line arguments
parse_args <- function() {
  parser <- ArgumentParser()
  parser$add_argument("-f", "--filePath", type = "character", help = "Path to the input files")
  parser$add_argument("-c", "--conditions", type = "character", help = "Treatments to compare (comma-separated)")
  parser$add_argument("-r", "--refStr", type = "character", help = "Reference condition")
  parser$add_argument("-t", "--fc_cutoff", type = "numeric", help = "Fold change cutoff", default = 2)
  parser$add_argument("-p", "--adjusted_p_value_cutoff", type = "numeric", help = "Adjusted p-value cutoff", default = 0.05)
  parser$add_argument("-m", "--min_read_cutoff", type = "numeric", help = "Minimum read count cutoff", default = 10)

  args <- parser$parse_args()
  return(args)
}

# Parse command line arguments
opt <- parse_args()


# Function Definitions

# Function for saving plots in multiple formats
save_plot <- function(plot, filename) {
  svg_filename <- paste0(filename, ".svg")
  pdf_filename <- paste0(filename, ".pdf")

  # Save as SVG
  ggsave(svg_filename, plot = plot, path = outputPath, width = 8, height = 6, dpi = 600, device = 'svg')

  # Save as PDF
  ggsave(pdf_filename, plot = plot, path = outputPath, width = 8, height = 6, dpi = 600, device = 'pdf')
}

# Write results to Excel-compatible files with raw counts and Gene Symbol
write_results <- function(data, filename) {
  write.table(data, file = here::here(outputPath, filename), sep = "\t", row.names = FALSE, quote = FALSE)
}

save_results <- function(data, filename) {
  fwrite(data, file = here::here(outputPath, filename), sep = "\t")
}

# Assign values to variables
filePath <- opt$filePath
conditions <- unlist(strsplit(opt$conditions, ","))
refStr <- opt$refStr
FC_cutoff <- as.numeric(opt$fc_cutoff)
adjusted_p_value_cutoff <- as.numeric(opt$adjusted_p_value_cutoff)
min_read_cutoff <- as.numeric(opt$min_read_cutoff)

# Set input and output paths
inputPath <- file.path(filePath, "inputFiles")
outputPath <- file.path(filePath, "gseaFiles")

# Create inputFiles folder if it doesn't exist
if (!dir.exists(inputPath)) {
  dir.create(inputPath, recursive = TRUE)
} else {
  cat("Input directory already exists. Skipping creation.\n")
}

# Create outputFiles folder if it doesn't exist
if (!dir.exists(outputPath)) {
  dir.create(outputPath, recursive = TRUE)
} else {
  cat("Output directory already exists. Skipping creation.\n")
}


# Copy count matrix and sample information files to inputFiles folder
countMatrixFiles <- list.files(pattern = "*Counts.t.*")
infoFiles <- list.files(pattern = "*Info.t.*")

if (length(countMatrixFiles) > 0) {
  # Specify the destination folder for count matrix files
  destination_count_matrix <- file.path(inputPath, countMatrixFiles)

  # Copy count matrix files with overwrite
  file.copy(countMatrixFiles, destination_count_matrix, overwrite = TRUE)
}

if (length(infoFiles) > 0) {
  # Specify the destination folder for sample information files
  destination_info_files <- file.path(inputPath, infoFiles)

  # Copy sample information files with overwrite
  file.copy(infoFiles, destination_info_files, overwrite = TRUE)
}

countMatrixFile <- list.files(path = inputPath, pattern = ".*Counts.t.*")
sampleInfoFile <- list.files(path = inputPath, pattern = ".*Info.t.*")

countMatrix <- read.table(file.path(inputPath, countMatrixFile), header = TRUE, row.names = 1)
sampleInfo <- read.table(file.path(inputPath, sampleInfoFile), header = TRUE, row.names = 1)

# Filter rows based on minimum read count cutoff
countMatrix <- countMatrix[rowSums(countMatrix) >= min_read_cutoff, ]

# Subset countMatrix to include only samples present in sampleInfo
countMatrix <- countMatrix[colnames(countMatrix) %in% row.names(sampleInfo), ]

# DEseq - RNA seq analysis (Note: may need to change design element depending on comparison)
dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = sampleInfo, design = ~ Treatment)
dds$Treatment <- relevel(dds$Treatment, ref = refStr)

# Perform DESeq analysis
dds <- DESeq(dds)

# Extract DESeq results
res <- results(dds, tidy = TRUE)
#head(res)
gene_list <- res$log2FoldChange
names(gene_list) = res$row
head(gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
#entrez_gene_list <- mapIds(org.Hs.eg.db, keys = names(gene_list), column = "ENTREZID", keytype = "SYMBOL")
#entrez_gene_list <- entrez_gene_list[!is.na(entrez_gene_list)]
#entrez_gene_list = sort(entrez_gene_list, decreasing = TRUE)
#head(entrez_gene_list)

# List of GMT files
gmt_files <- c(
 # "C2.KEGG_Legacy.v2023.2.Hs.symbols.gmt",
 # "C3.TFT.v2023.2.Hs.symbols.gmt",
 # "C5.GO.BP.v2023.2.Hs.symbols.gmt",
 # "C5.GO.MF.v2023.2.Hs.symbols.gmt",
  "C6.ONCOGENIC.v2023.2.Hs.symbols.gmt",
  "C7.IMMUNOLOGIC.v2023.2.Hs.symbols.gmt",
  "HALLMARK.v2023.2.Hs.symbols.gmt",
  "PDAC.v2023.2.Hs.symbols.gmt"
)

# Loop over GMT files
for (gmt_file in gmt_files) {
  pathway_path <- paste("/Users/ebareke/Desktop/GOZANI/RNAseq/Signatures/GSEA/", gmt_file, sep = "")
  pathway_set <- read.gmt(pathway_path)
  #enrich_result <- enricher(gene = gene_list, TERM2GENE = pathway_set, universe = names(org.Hs.eg.db), pvalueCutoff = adjusted_p_value_cutoff, pAdjustMethod = "BH")
  ## Perform fgsea analysis
  gsea_result = GSEA(gene_list, TERM2GENE = pathway_set, pvalueCutoff = adjusted_p_value_cutoff, pAdjustMethod = "BH")
  ## generate GSEA plot
  #enrich_dot_plot <- dotplot(enrich_result, showCategory = 15) + ggtitle("Over-representation analysis : Day 9")
  #gsea_dot_plot <- dotplotGsea(data = gsea_result, topn = 15)
  gsea_dot_plot <- dotplot(gsea_result, showCategory = 10, split = ".sign") + facet_grid(.~.sign) + ggtitle("GeneSet Enrichment Analysis : Day 9") + theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
  ## Save the results in EXCEL-compatible format
  #enrich_result_df <- as.data.frame(enrich_result)
  #save_results(enrich_result_df, paste(gsub("\\.v2023\\.2\\.Hs\\.symbols\\.gmt$", "", gmt_file), "_oraResult.txt", sep = ""))
  gsea_result_df <- as.data.frame(gsea_result)
  save_results(gsea_result_df, paste(gsub("\\.v2023\\.2\\.Hs\\.symbols\\.gmt$", "", gmt_file), "_gseaResult.txt", sep = ""))
  ## Save enrichment plots
  #save_plot(enrich_dot_plot, paste(gsub("\\.v2023\\.2\\.Hs\\.symbols\\.gmt$", "", gmt_file), "_oraPlot", sep = ""))
  save_plot(gsea_dot_plot, paste(gsub("\\.v2023\\.2\\.Hs\\.symbols\\.gmt$", "", gmt_file), "_gseaPlot", sep = ""))
}
