# MIT License
#
# Copyright (c) 2023 Eric BAREKE (eric.bareke@mcgill.ca), Emma Carlson (emma.carlson@mail.mcgill.ca), and Majewski Lab (McGill University)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# How to run it ?
# Rscript DGEA.r -f /path/to/your/input/files -c condition1,condition2 -r reference_condition -t 4 -p 0.05 -m 30
#awk -F"\t" -v OFS="\t" '{print $NF,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25}' Desktop/Gozani-DeSeq/Day5/inputFiles/annotated_combined.counts| awk -F"\t" -v OFS="\t" '!seen[$1]++' > Desktop/counts.tsv

# List of required packages
libraries <- c(
  "clusterProfiler",
  "enrichplot",
  "ggplot2",
  "DESeq2",
  "EnhancedVolcano",
  "grDevices",
  "RColorBrewer",
  "pheatmap",
  "ggrepel",
  "dplyr",
  "argparse",
  "here",
  "org.Hs.eg.db",
  "pathfindR"
)

# Function to check, install, and load packages
check_install_load <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
#    # Install the package if not installed
#    if (package %in% BiocManager::installed()) {
#      BiocManager::install(package)
#    } else {
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


############# Function Definitions #############

# Function for saving plots in multiple formats
save_plot <- function(plot, filename) {
  png_filename <- paste0(filename, ".png")
#  pdf_filename <- paste0(filename, ".pdf")

  # Save as PNG
  ggsave(png_filename, plot = plot, path = outputPath, dpi = 600, device = 'png')

  # Save as PDF
#  ggsave(pdf_filename, plot = plot, path = outputPath, dpi = 600, device = 'pdf')
}

# Perform inner joins
merge_and_rename <- function(data, subset_data, suffix) {
  merged_data <- merge(counts(dds, normalized = FALSE), subset_data, by = "row.names", all = FALSE)
  colnames(merged_data)[1] <- "GeneSymbol"
  colnames(merged_data)[-1] <- paste0(colnames(merged_data)[-1], suffix)
  return(merged_data)
}

# Write results to Excel-compatible files with raw counts and Gene Symbol
write_results <- function(data, filename) {
  write.table(data, file = here::here(outputPath, filename), sep = "\t", row.names = FALSE, quote = FALSE)
}

# Function to perform GO enrichment analysis
perform_GO_enrichment <- function(gene_list, direction, output_path, adjusted_p_value_cutoff = 0.1) {
  GO_results <- enrichGO(
    gene = gene_list,
    universe = background,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "fdr",
    minGSSize = 10,
    maxGSSize = 2000,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  #Change so GO_results is run once
  if (!is.null(GO_results)) {
    if (dim(GO_results)[1] > 0) {
      # Plot and save charts...
      write.csv(GO_results, file = paste0(output_path, "/GO_enrichment_results_", direction, ".csv"), row.names = FALSE)
      bpResults <- filter(GO_results, ONTOLOGY == "BP", p.adjust < adjusted_p_value_cutoff)
      if (dim(bpResults)[1] > 0) {
        save_enrichment_plots(bpResults, direction, "GO", "BP", output_path)
      }
      mfResults <- filter(GO_results, ONTOLOGY == "MF", p.adjust < adjusted_p_value_cutoff)
      if (dim(mfResults)[1] > 0) {
        save_enrichment_plots(mfResults, direction, "GO", "MF", output_path)
      }
      ccResults <- filter(GO_results, ONTOLOGY == "CC", p.adjust < adjusted_p_value_cutoff)
      if (dim(ccResults)[1] > 0) {
        save_enrichment_plots(ccResults, direction, "GO", "CC", output_path)
      }
    }
  }

}

perform_KEGG_enrichment <- function(gene_list, direction, output_path, adjusted_p_value_cutoff = 0.1) {
  KEGG_results <- enrichKEGG(
    gene = gene_list,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = adjusted_p_value_cutoff,
    pAdjustMethod = "BH",
    universe = background,
    minGSSize = 10,
    maxGSSize = 2000,
    qvalueCutoff = adjusted_p_value_cutoff,
    use_internal_data = FALSE
  )

  # Plot and save charts...
  if(!is.null(KEGG_results)) {
	if (dim(KEGG_results)[1] > 0) {
    save_enrichment_plots(KEGG_results, direction, "KEGG", "KEGG", output_path)
  		}
	}
}

# Save enrichment plots
#save_enrichment_plots <- function(results, direction, analysis_type, ont_type, output_path) {
#  bar_plot <- barplot(results, showCategory = 15)
#  dot_plot <- dotplot(results, showCategory = 15)
#  write.csv(results, file = paste0(output_path, analysis_type, ifelse(!is.null(direction), paste0("_", direction), ""), ifelse(!is.null(ont_type), paste0("_", ont_type), "")), "_results.csv"), row.names = FALSE)
#    save_plot(bar_plot, paste0("barPlot_", analysis_type, ifelse(!is.null(direction), paste0("_", direction), ""), ifelse(!is.null(ont_type), paste0("_", ont_type), "")))
#    save_plot(dot_plot, paste0("dotPlot_", analysis_type, ifelse(!is.null(direction), paste0("_", direction), ""), ifelse(!is.null(ont_type), paste0("_", ont_type), "")))
#
#}

save_enrichment_plots <- function(results, direction, analysis_type, ont_type, output_path) {
  if (ont_type != "KEGG") {
    SimResults <- simplify(results, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
    bar_plot <- barplot(SimResults, showCategory = 15)
    dot_plot <- dotplot(SimResults, showCategory = 15)
    write.csv(SimResults@result, file = paste0(output_path, analysis_type, ifelse(!is.null(direction), paste0("_", direction), ""), ifelse(!is.null(ont_type), paste0("_", ont_type), ""), "_results.csv"), row.names = FALSE)
  } else {
    bar_plot <- barplot(results, showCategory = 15)
    dot_plot <- dotplot(results, showCategory = 15)
    write.csv(results, file = paste0(output_path, analysis_type, ifelse(!is.null(direction), paste0("_", direction), ""), ifelse(!is.null(ont_type), paste0("_", ont_type), ""), "_results.csv"), row.names = FALSE)
  }

  save_plot(bar_plot, paste0("barPlot_", analysis_type, ifelse(!is.null(direction), paste0("_", direction), ""), ifelse(!is.null(ont_type), paste0("_", ont_type), "")))
  save_plot(dot_plot, paste0("dotPlot_", analysis_type, ifelse(!is.null(direction), paste0("_", direction), ""), ifelse(!is.null(ont_type), paste0("_", ont_type), "")))
}



# Function to perform pathfindR analysis
perform_pathfindR_analysis <- function(input_data, direction, output_path) {
  # Format data
  input_df <- data.frame(
    Gene.symbol = input_data$GeneSymbol,
    logFC = input_data$log2FoldChange,
    adj.P.Val = input_data$padj
  )

  # Define gene sets
  gene_sets <- c("KEGG", "GO-BP", "GO-MF")

  for (gene_set in gene_sets) {
    # Run pathfindR and save results
    output <- run_pathfindR(input_df, gene_sets = gene_set)
    clustered_output <- cluster_enriched_terms(output)

    # Create enrichment chart
    enrichment_chart(
      result_df = output,
      top_terms = 10
    )

    # Create clustered enrichment chart
    clustered_chart <- enrichment_chart(clustered_output, plot_by_cluster = TRUE)

    # Save plots in different formats with unique names
    save_plot(clustered_chart, paste0("clusterPlot_", direction, gene_set))
  }
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
# countMatrix <- countMatrix[rowSums(countMatrix) >= min_read_cutoff, ]

# Subset countMatrix to include only samples present in sampleInfo
countMatrix <- countMatrix[colnames(countMatrix) %in% row.names(sampleInfo), ]

# DEseq - RNA seq analysis (Note: may need to change design element depending on comparison)
dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = sampleInfo, design = ~ Treatment)

# Determine the number of samples
num_samples <- ncol(counts(dds))

# Filter out lowly expressed genes
keep <- rowSums(counts(dds) >= min_read_cutoff) >= 0.5 * num_samples
dds <- dds[keep,]

dds$Treatment <- relevel(dds$Treatment, ref = refStr)

# Perform DESeq analysis
dds <- DESeq(dds)

# Access the normalized counts from DESeq2
normalized_counts <- counts(dds, normalized = TRUE)

# Calculate sample correlation matrix
sample_cor_matrix <- cor(normalized_counts)

# Create a heatmap using pheatmap
Samples <- pheatmap(sample_cor_matrix,
                    clustering_distance_rows = "correlation",
                    clustering_distance_cols = "correlation",
                    color = colorRampPalette(c("blue", "white", "red"))(50),
                    main = "Sample Correlation Heatmap")

# Save the plot in both formats
save_plot(Samples, "correlationPlot")

# Extract DESeq results
res0 <- results(dds)

bpl0 <- sapply( levels(dds$Treatment), function(levels) rowMeans( counts(dds,normalized=TRUE)[,dds$Treatment==levels]))
mrg0 <- merge(bpl0, res0, by = "row.names", all = TRUE)
res <- mrg0[,-1]
row.names(res) <- mrg0[,1]
res <- as.data.frame(res)

#  # Filter significant DEGs based on adjusted p-value and fold change
resSig <- res[!is.na(res$padj) & (res$padj < adjusted_p_value_cutoff) & (abs(res$log2FoldChange) > log2(FC_cutoff)), ]

# Separate upregulated and downregulated genes
upregulatedResSig <- resSig[resSig$log2FoldChange > log2(FC_cutoff), ]
downregulatedResSig <- resSig[resSig$log2FoldChange < -log2(FC_cutoff), ]

# Merge and rename
merged_upregulatedResSig <- merge_and_rename(counts(dds, normalized = FALSE), upregulatedResSig, "_up")
merged_downregulatedResSig <- merge_and_rename(counts(dds, normalized = FALSE), downregulatedResSig, "_down")
merged_resSig <- merge_and_rename(counts(dds, normalized = FALSE), resSig, "_significant")
merged_res <- merge_and_rename(counts(dds, normalized = FALSE), res, "_raw")

# Write the merged dataframes to files
write_results(merged_upregulatedResSig, "upregulated_results.txt")
write_results(merged_downregulatedResSig, "downregulated_results.txt")
write_results(merged_resSig, "significant_results.txt")
write_results(merged_res, "unfiltered_results.txt")


# PCA
rld <- rlogTransformation(dds)
rv <- rowVars(assay(rld), useNames = TRUE)
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

custom_colors <- list(
  Treatment = setNames(
    c("darkgreen", "magenta"),
    levels(sampleInfo$Treatment)
  )
)


# EnrichedHeatMap

# Create custom_colors dynamically
# custom_colors <- list(
#  Treatment = setNames(
#    c("darkgreen", "magenta"),
#    c(refStr, setdiff(conditions, refStr))
#  )
#)

countMatrix <- data.matrix(countMatrix)
top_varied_genes <- rownames(countMatrix)[order(rowVars(countMatrix, useNames = TRUE), decreasing = TRUE)][1:500]
MostVariedGenes <- pheatmap(
  countMatrix[top_varied_genes, ],
  main = "Top 500 most variable genes across samples",
  annotation_col = sampleInfo,
  #annotation_colors = custom_colors,
  scale = "row",
  show_rownames = FALSE,
  legend_title = "Treatment Time",
  color = colorRampPalette(c("blue", "white", "red"))(50)
)

save_plot(MostVariedGenes, "heatmapPlot")


# Volcano plot of up and down regulated genes
keyvals <- ifelse(
  (res$log2FoldChange < -log2(FC_cutoff) & res$padj < adjusted_p_value_cutoff), 'blue',
  ifelse((res$log2FoldChange > log2(FC_cutoff) & res$padj < adjusted_p_value_cutoff), 'red',
         'white')
)
keyvals[is.na(keyvals)] <- 'white'
names(keyvals)[keyvals == 'red'] <- 'Up-regulated'
names(keyvals)[keyvals == 'white'] <- ''
names(keyvals)[keyvals == 'blue'] <- 'Down-regulated'

# Set xlim and ylim dynamically
DEG <- EnhancedVolcano(res,
                                lab = "",
                                x = 'log2FoldChange',
                                y = 'padj',
                                title = 'Differentially Expressed Genes',
                                subtitle = paste('Down-regulated Genes:', sum(keyvals == 'blue'),
                                                 ' -- Up-regulated Genes:', sum(keyvals == 'red')),
                                caption = paste("FC cutoff =", FC_cutoff, "; p-adj cutoff =", adjusted_p_value_cutoff),
                                xlab = bquote(~Log[2]~ 'Fold Change'),
                                xlim = c(-10,10),
                                ylim = c(0,20),
                                pCutoff = 0.1,
                                FCcutoff = 0.0,
                                pointSize = 2.0,
                                colCustom = keyvals,
                                legendPosition = 'bottom',
                                legendLabSize = 14,
                                legendIconSize = 4.0,
                                gridlines.major = TRUE,
                                gridlines.minor = FALSE,
                                border = 'partial',
                                borderWidth = 1.2,
                                borderColour = 'black', max.overlaps = Inf
)

DEG <- DEG + theme(plot.subtitle = element_text(hjust = 0, size = 12))

save_plot(DEG, "volcanoPlot")


# Define input parameters as needed...

universe <- merged_res$GeneSymbol
background <- mapIds(org.Hs.eg.db, keys = universe, column = "ENTREZID", keytype = "SYMBOL")

background <- background[!is.na(background)]

listDN <- merged_downregulatedResSig$GeneSymbol
listDN <- mapIds(org.Hs.eg.db, keys = listDN, column = "ENTREZID", keytype = "SYMBOL")
listDN <- listDN[!is.na(names(listDN))]

# Main analysis for downregulated genes
perform_GO_enrichment(listDN, "Down", outputPath, adjusted_p_value_cutoff)
#perform_GO_enrichment(listDN, "Down", "MF", outputPath, adjusted_p_value_cutoff)
#perform_GO_enrichment(listDN, "Down", "CC", outputPath, adjusted_p_value_cutoff)
perform_KEGG_enrichment(listDN, "Down", outputPath, adjusted_p_value_cutoff)
#perform_pathfindR_analysis(merged_downregulatedResSig, "Down", outputPath)

# Define other input parameters as needed

listUP <- merged_upregulatedResSig$GeneSymbol
listUP <- mapIds(org.Hs.eg.db, keys = listUP, column = "ENTREZID", keytype = "SYMBOL")
listUP <- listUP[!is.na(names(listUP))]

# Main analysis for upregulated genes
perform_GO_enrichment(listUP, "Up", outputPath, adjusted_p_value_cutoff)
#perform_GO_enrichment(listUP, "Up", "MF", outputPath, adjusted_p_value_cutoff)
#perform_GO_enrichment(listUP, "Up", "CC", outputPath, adjusted_p_value_cutoff)
perform_KEGG_enrichment(listUP, "Up", outputPath,adjusted_p_value_cutoff)
#perform_pathfindR_analysis(merged_upregulatedResSig, "Up", outputPath)
