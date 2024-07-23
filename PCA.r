# MIT License
# How to run it ?
# Rscript PCA.r -f /path/to/your/input/files -c condition1,condition2 -r reference_condition -t 4 -p 0.05 -m 30

# List of required packages
libraries <- c(
  "ggplot2",
  "DESeq2",
  "grDevices",
  "RColorBrewer",
  "pheatmap",
  "ggrepel",
  "dplyr",
  "argparse",
  "here"
)

# Function to check, install, and load packages
check_install_load <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    # Install the package if not installed
    if (package %in% BiocManager::installed()) {
      BiocManager::install(package)
    } else {
      install.packages(package)
    }
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

############# Function Definitions #############

# Function for saving plots in multiple formats
save_plot <- function(plot, filename) {
  png_filename <- paste0(filename, ".svg")
  pdf_filename <- paste0(filename, ".pdf")

  # Save as PNG
  ggsave(png_filename, plot = plot, path = outputPath, dpi = 600, device = 'svg')

  # Save as PDF
  ggsave(pdf_filename, plot = plot, path = outputPath, dpi = 600, device = 'pdf')
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
outputPath <- file.path(filePath, "outputFiles")

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
countMatrixFiles <- list.files(pattern = "*Counts.t*")
infoFiles <- list.files(pattern = "*Info.t*")

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

# Access the normalized counts from DESeq2
#counts <- counts(dds, normalized = FALSE)
normalized_counts <- counts(dds, normalized = TRUE)

## Calculate sample correlation matrix
#sample_cor_matrix <- cor(normalized_counts)

# Step 2: Compute the correlation matrix
cor_matrix <- cor(normalized_counts)

# Step 3: Perform hierarchical clustering
hc <- hclust(as.dist(1 - cor_matrix), method = "complete")

# Step 4: Extract sample correlation matrix based on hierarchical clustering order
ordered_data <- normalized_counts[, hc$order]  # reorder columns based on clustering
sample_cor <- cor(ordered_data)   # compute correlation matrix of reordered data

# Create a heatmap using pheatmap
Samples <- pheatmap(sample_cor,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    #clustering_method = "complete",
                    color = colorRampPalette(c("blue", "white", "red"))(100),
                    main = "Correlation Heatmap with Hierarchical Clustering")

# Save the plot in both formats
save_plot(Samples, "correlationPlot")


# PCA
rld <- rlogTransformation(dds)
rv <- rowVars(assay(rld))
vsd <- varianceStabilizingTransformation(dds)
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pc <- prcomp(t(assay(vsd)[select, ]))
condition <- sampleInfo$Treatment
scores <- data.frame(pc$x, condition)
samples <- rownames(sampleInfo)
explained_variance <- summary(pc)$importance[2, ] * 100

AllSamples <- ggplot(scores, aes(x = PC1, y = PC2, col = factor(condition))) +
  xlab(paste0("PC1 Variance Explained: ", explained_variance[1], "%")) +
  ylab(paste0("PC2 Variance Explained: ", explained_variance[2], "%")) +
  geom_point(size = 5) +
  geom_label_repel(aes(label = samples), show.legend = FALSE, max.overlaps = 50) +
  ggtitle("Principal Component Analysis") +
  scale_colour_brewer(name = " ", palette = "Set1") +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold'),
    legend.position = "right",
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.background = element_rect(color = 'black', fill = NA),
    axis.line.y = element_line(linewidth = 0.5, color = "black"),
    axis.line.x = element_line(linewidth = 0.5, color = "black"),
    axis.text = element_text(color = "black"),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),
    panel.grid.major.x = element_line(color = "grey", linetype = "dashed")
  )

save_plot(AllSamples, "pcaPlot")

# EnrichedHeatMap

# Create custom_colors dynamically
custom_colors <- list(
  Treatment = setNames(
    c("dodgerblue2", "firebrick2"),
    c(refStr, setdiff(conditions, refStr))
  )
)

countMatrix <- data.matrix(countMatrix)
top_varied_genes <- rownames(countMatrix)[order(rowVars(countMatrix, useNames = TRUE), decreasing = TRUE)][1:100]
MostVariedGenes <- pheatmap(
  countMatrix[top_varied_genes, ],
  main = "Top 100 most variable genes across samples",
  annotation_col = sampleInfo,
  annotation_colors = custom_colors,
  scale = "row",
  show_rownames = FALSE,
  legend_title = "Treatment Time",
  color = brewer.pal(8, "Dark2")
)

save_plot(MostVariedGenes, "heatmapPlot")
