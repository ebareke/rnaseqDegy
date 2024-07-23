# How to run it ?
# Rscript GSEA.r -f /path/to/your/input/files -c VEH,T_DM1 -r VEH -t 2 -p 0.05 -m 10

# List of required packages
libraries <- c(
  "DESeq2",
  "org.Hs.eg.db",
  "tidyverse",
  "fgsea",
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
  parser$add_argument("-t", "--fc_cutoff", type = "numeric", help = "Fold change cutoff", default = 1.5)
  parser$add_argument("-p", "--adjusted_p_value_cutoff", type = "numeric", help = "Adjusted p-value cutoff", default = 0.05)
  parser$add_argument("-m", "--min_read_cutoff", type = "numeric", help = "Minimum read count cutoff", default = 20)

  args <- parser$parse_args()
  return(args)
}

# Parse command line arguments
opt <- parse_args()


# Function Definitions

# Function for saving plots in multiple formats
save_plot <- function(plot, filename) {
  png_filename <- paste0(filename, ".png")
  #pdf_filename <- paste0(filename, ".pdf")

  # Save as PNG
  ggsave(png_filename, plot = plot, path = outputPath, width = 16, height = 12, dpi = 600, device = 'png')

  # Save as PDF
  #ggsave(pdf_filename, plot = plot, path = outputPath,  width = 16, height = 12, dpi = 600, device = 'pdf')
}

# Write results to Excel-compatible files with raw counts and Gene Symbol
write_results <- function(data, filename) {
  write.table(data, file = here::here(outputPath, filename), sep = "\t", row.names = FALSE, quote = FALSE)
}

save_results <- function(data, filename) {
  fwrite(data, file = here::here(outputPath, filename), sep = "\t")
}

#### Perform GSEA
GSEA = function(gene_list, GMT_file, pval) {
  set.seed(54321)
  library(dplyr)
  library(fgsea)

  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGMT = fgsea::gmtPathways(GMT_file)

  fgRes <- fgsea::fgseaMultilevel(pathways = myGMT,
                           stats = gene_list,
                           minSize=20, ## minimum gene set size
                           maxSize=2000) %>%
                  as.data.frame() %>%
                  dplyr::filter(padj < !!pval) %>%
                  arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))

  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = myGMT,
                                      stats = gene_list)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))

  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 20),
                  tail(fgRes, n = 20 ))

  total_up = sum(fgRes$Enrichment == "Up-regulated")
  total_down = sum(fgRes$Enrichment == "Down-regulated")
  header = paste0("Top 40 Significant Enrichments (FDR = 0.1) - [Total: Up = ", total_up, ", Down = ", total_down, "]")

  colos = setNames(c("firebrick2", "dodgerblue2"),
                 c("Up-regulated", "Down-regulated"))

 g1= ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
  geom_point( aes(fill = Enrichment, size = size), shape=21) +
  scale_fill_manual(values = colos ) +
  scale_size_continuous(range = c(2,10)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x="Enrichments", y="Normalized Enrichment Score",
       title=header) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))



  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
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

gene_list <- res$log2FoldChange
names(gene_list) = res$row
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

# List of GMT files
gmt_files <- c(
  "M10431_GRUETZMANN_PANCREATIC_CANCER_DN.v2023.2.Hs.gmt",
  "M15193_GRUETZMANN_PANCREATIC_CANCER_UP.v2023.2.Hs.gmt",
  "M19938_REACTOME_REGULATION_OF_GENE_EXPRESSION_IN_EARLY_PANCREATIC_PRECURSOR_CELLS.v2023.2.Hs.gmt",
  "M27123_SU_PANCREAS.v2023.2.Hs.gmt",
  "M9726_KEGG_PANCREATIC_CANCER.v2023.2.Hs.gmt"
)

# Loop over GMT files
for (gmt_file in gmt_files) {
  # Construct the pathway object
  pathways <- paste("/Users/ebareke/Desktop/Or_PDAC/Signatures/", gmt_file, sep = "")

  # Perform fgsea analysis
  fgsea_res = GSEA(gene_list, pathways, pval = adjusted_p_value_cutoff)

  # Extract GSEA result
  res <- fgsea_res$Results

  # Save the results in EXCEL-compatible format
   save_results(res, paste(gsub("\\.gmt$", "", gmt_file), "_GSEA.txt", sep = ""))

  # Create a plot
  plot <- fgsea_res$Plot

  # Save the plot in both formats
  save_plot(plot, paste(gsub("\\.gmt$", "", gmt_file), "_GSEA", sep = ""))
}

# Delete the Rplots.pdf file
# file.remove("Rplots.pdf")
