.run_engine <- function(args) {
  ###########################################################################
  # Internal engine for rnaseqDegy
  #
  # This is a function version of the original CLI pipeline. It expects a
  # named list `args` constructed by run_rnaseqDegy().
  # All heavy lifting (QC, DE, ORA, GSEA) happens here.
  ###########################################################################

  # Convenience local variables ------------------------------------------------
  run_batch_correction <- isTRUE(args$batch_correction)
  run_ora              <- isTRUE(args$run_ora)
  run_gsea             <- isTRUE(args$run_gsea)

  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive = TRUE)
  }

  prefix_tag <- ifelse(is.null(args$prefix) || args$prefix == "", "", paste0(args$prefix, "_"))

  # Helper theme + saving utilities -------------------------------------------

  base_theme <- function() {
    ggplot2::theme_bw(base_size = 14) +
      ggplot2::theme(
        text = ggplot2::element_text(size = 14, face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.title = ggplot2::element_text(face = "bold"),
        legend.text = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(face = "bold")
      )
  }

  save_gg <- function(p, filename, width = 7, height = 5) {
    ggplot2::ggsave(
      filename = file.path(args$output_dir, paste0(prefix_tag, filename, ".pdf")),
      plot = p, width = width, height = height
    )
    ggplot2::ggsave(
      filename = file.path(args$output_dir, paste0(prefix_tag, filename, ".png")),
      plot = p, width = width, height = height, dpi = 300
    )
  }

  message("[1/8] Reading input files...")

  # Data import ---------------------------------------------------------------

  counts_df <- data.table::fread(args$counts)
  if (ncol(counts_df) < 3) stop("Count matrix must have >= 2 samples.")

  colnames(counts_df)[1] <- "GeneID"

  gene_ids  <- counts_df$GeneID
  count_mat <- as.matrix(counts_df[, -1, with = FALSE])
  rownames(count_mat) <- gene_ids
  mode(count_mat) <- "numeric"

  design_df <- data.table::fread(args$design)
  required_cols <- c("Sample", "Batch", "Condition")
  if (!all(required_cols %in% colnames(design_df))) {
    stop("Design file must contain columns: Sample, Batch, Condition")
  }

  # Match samples -------------------------------------------------------------

  samples_in_both <- intersect(colnames(count_mat), design_df$Sample)
  if (length(samples_in_both) < 2) stop("Need at least 2 overlapping samples between counts and design.")

  count_mat <- count_mat[, samples_in_both, drop = FALSE]
  design_df <- design_df[match(samples_in_both, design_df$Sample), ]
  rownames(design_df) <- design_df$Sample

  # Subset to conditions of interest -----------------------------------------

  keep_cond <- design_df$Condition %in% c(args$condition_reference, args$condition_test)
  count_mat <- count_mat[, keep_cond, drop = FALSE]
  design_df <- design_df[keep_cond, ]

  if (nrow(design_df) < 2) stop("Not enough samples for selected conditions.")

  design_df$Condition <- factor(design_df$Condition, levels = c(args$condition_reference, args$condition_test))
  design_df$Batch     <- factor(design_df$Batch)

  # Species / annotation helpers ---------------------------------------------

  if (args$species == "human") {
    OrgDb <- org.Hs.eg.db::org.Hs.eg.db
    ens_prefix <- "ENSG"
  } else {
    OrgDb <- org.Mm.eg.db::org.Mm.eg.db
    ens_prefix <- "ENSMUSG"
  }

  infer_id_type <- function(ids, ens_prefix) {
    ens_frac <- mean(grepl(paste0("^", ens_prefix), ids, ignore.case = TRUE))
    if (ens_frac > 0.5) return("ENSEMBL")
    return("SYMBOL")
  }

  id_type <- infer_id_type(gene_ids, ens_prefix)
  message("Detected ID type: ", id_type)

  map_to_symbols <- function(ids) {
    if (id_type == "SYMBOL") return(ids)
    mapped <- AnnotationDbi::mapIds(
      OrgDb,
      keys = ids,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    as.character(mapped)
  }

  map_to_entrez <- function(symbols) {
    AnnotationDbi::mapIds(
      OrgDb,
      keys = symbols,
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
  }

  symbols <- map_to_symbols(rownames(count_mat))

  # Filtering -----------------------------------------------------------------

  message("[2/8] Filtering low-count genes...")

  n_samples <- ncol(count_mat)
  keep_gene <- rowSums(count_mat >= args$filter_min_count) >= (args$filter_min_prop * n_samples)
  filtered_counts  <- count_mat[keep_gene, , drop = FALSE]
  filtered_symbols <- symbols[keep_gene]

  message(sum(keep_gene), " genes retained after filtering.")

  # Normalisation + QC matrices ----------------------------------------------

  message("[3/8] Normalising counts and preparing QC matrices...")

  coldata <- design_df

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData   = coldata,
    design    = ~ Batch + Condition
  )

  dds$Condition <- stats::relevel(dds$Condition, ref = args$condition_reference)

  vst_obj <- DESeq2::vst(dds, blind = TRUE)
  vst_mat <- SummarizedExperiment::assay(vst_obj)

  # Batch correction ---------------------------------------------------------

  if (run_batch_correction) {
    message("[4/8] Applying batch correction for QC...")
    batch <- coldata$Batch
    vst_batch_corrected <- limma::removeBatchEffect(vst_mat, batch = batch)
  } else {
    vst_batch_corrected <- NULL
  }

  # PCA helper ---------------------------------------------------------------

  pca_plot <- function(mat, title) {
    rv <- matrixStats::rowVars(mat)
    sel <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
    mat_sel <- t(mat[sel, , drop = FALSE])
    pc <- stats::prcomp(mat_sel, scale. = TRUE)
    var_expl <- (pc$sdev^2) / sum(pc$sdev^2) * 100

    df <- data.frame(pc$x[, 1:2, drop = FALSE])
    df$Sample    <- rownames(df)
    df$Condition <- coldata[df$Sample, "Condition"]
    df$Batch     <- coldata[df$Sample, "Batch"]

    ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, colour = Condition, shape = Batch)) +
      ggplot2::geom_point(size = 4) +
      ggrepel::geom_text_repel(ggplot2::aes(label = Sample), size = 3.5,
                               show.legend = FALSE, max.overlaps = Inf) +
      ggplot2::labs(
        title = title,
        x = sprintf("PC1 (%.1f%%)", var_expl[1]),
        y = sprintf("PC2 (%.1f%%)", var_expl[2])
      ) +
      ggplot2::scale_colour_brewer(palette = "Set1") +
      base_theme()
  }

  # PCA plots -----------------------------------------------------------------

  if (!run_batch_correction) {
    message("[5/8] Generating PCA (no batch correction, as requested)...")
    p_pca <- pca_plot(vst_mat, "PCA (no batch correction)")
    save_gg(p_pca, "QC_PCA_noBatch")
  } else {
    message("[5/8] Generating PCA after batch correction...")
    p_pca_bc <- pca_plot(vst_batch_corrected, "PCA (after batch correction)")
    save_gg(p_pca_bc, "QC_PCA_batchCorrected")
  }

  # Sample correlation heatmap -----------------------------------------------

  message("[6/8] Sample correlation heatmap...")

  cor_mat <- stats::cor(vst_mat, method = args$corr_method)

  p_cor <- pheatmap::pheatmap(
    cor_mat,
    main = sprintf("Sample correlation (%s)", args$corr_method),
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    color = grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(50),
    fontsize = 11,
    fontsize_row = 11,
    fontsize_col = 11,
    border_color = NA
  )

  pdf(file.path(args$output_dir, paste0(prefix_tag, "QC_sampleCorrelation.pdf")),
      width = max(6, 0.35 * ncol(cor_mat)),
      height = max(6, 0.35 * ncol(cor_mat)))
  print(p_cor)
  grDevices::dev.off()

  png(file.path(args$output_dir, paste0(prefix_tag, "QC_sampleCorrelation.png")),
      width = max(6, 0.35 * ncol(cor_mat)),
      height = max(6, 0.35 * ncol(cor_mat)),
      units = "in", res = 300)
  print(p_cor)
  grDevices::dev.off()

  # Top 500 most variable genes heatmap --------------------------------------

  message("[7/8] Heatmap of top 500 most variable genes...")

  v_var <- matrixStats::rowVars(vst_mat)
  ord <- order(v_var, decreasing = TRUE)
  ng <- min(500, length(ord))
  sel_genes <- ord[seq_len(ng)]

  ann <- data.frame(
    Batch = coldata$Batch,
    Condition = coldata$Condition
  )
  rownames(ann) <- rownames(coldata)

  batch_levels <- unique(ann$Batch)
  cond_levels  <- unique(ann$Condition)

  ha_colors <- list(
    Batch = stats::setNames(
      RColorBrewer::brewer.pal(max(3, length(batch_levels)), "Set2")[seq_along(batch_levels)],
      batch_levels
    ),
    Condition = stats::setNames(
      RColorBrewer::brewer.pal(max(3, length(cond_levels)), "Set1")[seq_along(cond_levels)],
      cond_levels
    )
  )

  p_hm <- pheatmap::pheatmap(
    vst_mat[sel_genes, , drop = FALSE],
    scale = "row",
    main = sprintf("Top %d most variable genes", ng),
    annotation_col = ann,
    annotation_colors = ha_colors,
    show_rownames = FALSE,
    color = grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(50),
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 10,
    border_color = NA
  )

  pdf(file.path(args$output_dir, paste0(prefix_tag, "QC_topVar500_heatmap.pdf")),
      width = max(7, 0.3 * ncol(vst_mat)),
      height = 8)
  print(p_hm)
  grDevices::dev.off()

  png(file.path(args$output_dir, paste0(prefix_tag, "QC_topVar500_heatmap.png")),
      width = max(7, 0.3 * ncol(vst_mat)),
      height = 8,
      units = "in", res = 300)
  print(p_hm)
  grDevices::dev.off()

  # Differential expression ---------------------------------------------------

  message("[8/8] Running differential expression (", args$de_method, ")...")

  res_table   <- NULL
  norm_counts <- NULL

  if (identical(args$de_method, "DESeq2")) {

    dds_de <- dds
    dds_de <- DESeq2::DESeq(dds_de)

    res <- DESeq2::results(dds_de, contrast = c("Condition", args$condition_test, args$condition_reference))
    res <- as.data.frame(res)
    res$GeneID <- rownames(res)
    res$Symbol <- filtered_symbols[match(res$GeneID, rownames(filtered_counts))]

    norm_counts <- DESeq2::counts(dds_de, normalized = TRUE)
    norm_counts <- as.data.frame(norm_counts)
    norm_counts$GeneID <- rownames(norm_counts)
    norm_counts$Symbol <- filtered_symbols[match(norm_counts$GeneID, rownames(filtered_counts))]

    res_table <- res

  } else if (identical(args$de_method, "edgeR")) {

    group <- design_df$Condition
    batch <- design_df$Batch

    y <- edgeR::DGEList(counts = filtered_counts, group = group)
    y <- edgeR::calcNormFactors(y)

    design_mat <- stats::model.matrix(~ batch + group)
    y <- edgeR::estimateDisp(y, design_mat)
    fit <- edgeR::glmQLFit(y, design_mat)

    col_idx <- which(colnames(design_mat) == paste0("group", args$condition_test))
    if (length(col_idx) == 0) stop("Could not find appropriate contrast column for edgeR.")
    contrast_vec <- rep(0, ncol(design_mat))
    contrast_vec[col_idx] <- 1

    qlf <- edgeR::glmQLFTest(fit, contrast = contrast_vec)
    res <- edgeR::topTags(qlf, n = Inf)$table

    res$GeneID <- rownames(res)
    res$Symbol <- filtered_symbols[match(res$GeneID, rownames(filtered_counts))]

    res$log2FoldChange <- res$logFC
    res$padj <- stats::p.adjust(res$PValue, method = "BH")

    norm_counts <- edgeR::cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
    norm_counts <- as.data.frame(norm_counts)
    norm_counts$GeneID <- rownames(norm_counts)
    norm_counts$Symbol <- filtered_symbols[match(norm_counts$GeneID, rownames(filtered_counts))]

    res_table <- res
  } else {
    stop("Unsupported de_method: ", args$de_method)
  }

  # Save DE tables ------------------------------------------------------------

  data.table::fwrite(res_table, file.path(args$output_dir, paste0(prefix_tag, "DE_results_all.tsv")), sep = "\t")
  data.table::fwrite(norm_counts, file.path(args$output_dir, paste0(prefix_tag, "normalized_counts.tsv")), sep = "\t")

  # DEG subsets + volcano -----------------------------------------------------

  sig_idx <- which(!is.na(res_table$padj) & res_table$padj <= args$alpha &
                     abs(res_table$log2FoldChange) >= args$lfc_cutoff)
  res_sig <- res_table[sig_idx, , drop = FALSE]

  up_idx   <- which(res_table$padj <= args$alpha & res_table$log2FoldChange >= args$lfc_cutoff)
  down_idx <- which(res_table$padj <= args$alpha & res_table$log2FoldChange <= -args$lfc_cutoff)

  up_genes   <- res_table$Symbol[up_idx]
  down_genes <- res_table$Symbol[down_idx]

  sel_up <- head(up_idx[order(res_table$log2FoldChange[up_idx], decreasing = TRUE)], 10)
  sel_dn <- head(down_idx[order(res_table$log2FoldChange[down_idx], decreasing = FALSE)], 10)

  select_lab <- unique(stats::na.omit(c(res_table$Symbol[sel_up], res_table$Symbol[sel_dn])))

  keyvals <- ifelse(res_table$padj > args$alpha | is.na(res_table$padj), "grey80",
                    ifelse(res_table$log2FoldChange >= args$lfc_cutoff, "#d73027",
                           ifelse(res_table$log2FoldChange <= -args$lfc_cutoff, "#4575b4", "grey80")))

  names(keyvals)[keyvals == "#d73027"] <- "Up"
  names(keyvals)[keyvals == "#4575b4"] <- "Down"
  names(keyvals)[keyvals == "grey80"]  <- "NS"

  volcano_subtitle <- sprintf("Up: %d  |  Down: %d  (|log2FC| >= %.2f, FDR <= %.3f)",
                              length(up_idx), length(down_idx), args$lfc_cutoff, args$alpha)

  volcano <- EnhancedVolcano::EnhancedVolcano(
    res_table,
    lab = res_table$Symbol,
    x = "log2FoldChange",
    y = "padj",
    pCutoff = args$alpha,
    FCcutoff = args$lfc_cutoff,
    xlab = bquote(~Log[2]~" fold change"),
    title = paste0("DE: ", args$condition_test, " vs ", args$condition_reference),
    subtitle = volcano_subtitle,
    caption = paste("DE method:", args$de_method),
    selectLab = select_lab,
    colCustom = keyvals,
    legendPosition = "bottom",
    legendLabSize = 14,
    legendIconSize = 4,
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    border = "partial",
    borderWidth = 1,
    borderColour = "black",
    max.overlaps = Inf
  ) + base_theme()

  save_gg(volcano, "DE_volcano")

  # Enrichment: ORA (GO BP/MF/CC) --------------------------------------------

  if (run_ora) {
    message("[ORA] Running GO ORA (BP/MF/CC)...")

    deg_symbols <- unique(stats::na.omit(res_sig$Symbol))
    bg_symbols  <- unique(stats::na.omit(res_table$Symbol))

    bg_entrez <- map_to_entrez(bg_symbols)
    bg_entrez <- unique(bg_entrez[!is.na(bg_entrez)])

    deg_entrez <- map_to_entrez(deg_symbols)
    deg_entrez <- unique(deg_entrez[!is.na(deg_entrez)])

    if (length(deg_entrez) >= 10) {
      for (ont in c("BP", "MF", "CC")) {
        go_res <- clusterProfiler::enrichGO(
          gene          = deg_entrez,
          universe      = bg_entrez,
          OrgDb         = OrgDb,
          keyType       = "ENTREZID",
          ont           = ont,
          pAdjustMethod = "BH",
          pvalueCutoff  = 1,
          qvalueCutoff  = 1,
          readable      = TRUE
        )

        if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
          go_res <- clusterProfiler::simplify(go_res, cutoff = 0.7, by = "p.adjust", select_fun = min)

          go_df <- as.data.frame(go_res)
          data.table::fwrite(go_df, file.path(args$output_dir,
                                              paste0(prefix_tag, "ORA_GO_", ont, ".tsv")), sep = "\t")

          if (!is.null(args$gmt_terms) && args$gmt_terms != "") {
            terms_vec <- strsplit(args$gmt_terms, ",")[[1]]
            go_res <- go_res[go_res$ID %in% terms_vec | go_res$Description %in% terms_vec]
            if (nrow(as.data.frame(go_res)) == 0) next
          }

          p_dot <- enrichplot::dotplot(go_res, showCategory = 20) + base_theme() +
            ggplot2::ggtitle(paste("GO", ont, "ORA - dotplot"))
          save_gg(p_dot, paste0("ORA_GO_", ont, "_dotplot"))

          p_bar <- enrichplot::barplot(go_res, showCategory = 20) + base_theme() +
            ggplot2::ggtitle(paste("GO", ont, "ORA - barplot"))
          save_gg(p_bar, paste0("ORA_GO_", ont, "_barplot"))

          p_cnet <- enrichplot::cnetplot(go_res, showCategory = 10, foldChange = NULL) + base_theme() +
            ggplot2::ggtitle(paste("GO", ont, "ORA - cnetplot"))
          save_gg(p_cnet, paste0("ORA_GO_", ont, "_cnetplot"))

          go_res2 <- enrichplot::pairwise_termsim(go_res)
          p_ema <- enrichplot::emapplot(go_res2, showCategory = 20, cex_label_category = 0.8) + base_theme() +
            ggplot2::ggtitle(paste("GO", ont, "ORA - emapplot"))
          save_gg(p_ema, paste0("ORA_GO_", ont, "_emapplot"))
        }
      }
    } else {
      message("[ORA] Not enough DEG for ORA (need >= 10 entrez IDs). Skipping.")
    }
  }

  # Enrichment: GSEA ----------------------------------------------------------

  if (run_gsea) {
    if (is.null(args$gmt_dir) || args$gmt_dir == "" || is.null(args$gmt_files) || args$gmt_files == "") {
      stop("run_gsea=TRUE but gmt_dir or gmt_files not provided.")
    }

    message("[GSEA] Running GSEA on GMT files...")

    gmt_files <- strsplit(args$gmt_files, ",")[[1]]
    gmt_files <- trimws(gmt_files)

    ranked_df <- res_table
    ranked_df$Symbol_final <- ifelse(is.na(ranked_df$Symbol) | ranked_df$Symbol == "",
                                     ranked_df$GeneID, ranked_df$Symbol)

    gene_list <- ranked_df$log2FoldChange
    names(gene_list) <- ranked_df$Symbol_final
    gene_list <- sort(gene_list, decreasing = TRUE)

    gene_list <- gene_list[!is.na(gene_list)]
    gene_list <- gene_list[!duplicated(names(gene_list))]

    for (gmt in gmt_files) {
      gmt_path <- file.path(args$gmt_dir, gmt)
      if (!file.exists(gmt_path)) {
        warning("GMT file not found: ", gmt_path, " — skipping.")
        next
      }

      message("[GSEA] ", gmt)
      gmt_df <- clusterProfiler::read.gmt(gmt_path)

      gsea_res <- tryCatch({
        clusterProfiler::GSEA(
          geneList      = gene_list,
          TERM2GENE     = gmt_df,
          pAdjustMethod = "BH",
          pvalueCutoff  = args$alpha,
          minGSSize     = 10,
          maxGSSize     = 5000,
          verbose       = FALSE
        )
      }, error = function(e) {
        warning("GSEA failed for ", gmt, ": ", conditionMessage(e))
        NULL
      })

      if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
        message("[GSEA] No significant terms for ", gmt)
        next
      }

      gsea_df <- as.data.frame(gsea_res)
      data.table::fwrite(gsea_df, file.path(args$output_dir,
                                            paste0(prefix_tag, "GSEA_", gsub("\\\\.gmt$", "", gmt), ".tsv")),
                         sep = "\t")

      if (!is.null(args$gmt_terms) && args$gmt_terms != "") {
        terms_vec <- strsplit(args$gmt_terms, ",")[[1]]
        gsea_res <- gsea_res[gsea_res$ID %in% terms_vec | gsea_res$Description %in% terms_vec]
        if (nrow(as.data.frame(gsea_res)) == 0) next
      }

      base_name <- gsub("\\\\.gmt$", "", gmt)

      p_dot <- enrichplot::dotplot(gsea_res, showCategory = 20, split = ".sign") +
        ggplot2::facet_grid(. ~ .sign) + base_theme() +
        ggplot2::ggtitle(paste("GSEA -", base_name))
      save_gg(p_dot, paste0("GSEA_", base_name, "_dotplot"))

      p_bar <- enrichplot::barplot(gsea_res, showCategory = 20) + base_theme() +
        ggplot2::ggtitle(paste("GSEA -", base_name))
      save_gg(p_bar, paste0("GSEA_", base_name, "_barplot"))

      p_cnet <- enrichplot::cnetplot(gsea_res, showCategory = 10, foldChange = gene_list) + base_theme() +
        ggplot2::ggtitle(paste("GSEA -", base_name, " cnetplot"))
      save_gg(p_cnet, paste0("GSEA_", base_name, "_cnetplot"))

      gsea_res2 <- enrichplot::pairwise_termsim(gsea_res)
      p_ema <- enrichplot::emapplot(gsea_res2, showCategory = 20, cex_label_category = 0.8) + base_theme() +
        ggplot2::ggtitle(paste("GSEA -", base_name, " emapplot"))
      save_gg(p_ema, paste0("GSEA_", base_name, "_emapplot"))
    }
  }

  message("Pipeline completed successfully. Outputs written to: ", args$output_dir)

  invisible(NULL)
}
