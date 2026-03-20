library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  counts_df <- read_tsv(counts_csv, show_col_types = FALSE)

  meta <- read_csv(metafile_csv, show_col_types = FALSE)

  # put vP0 first so it becomes the reference level in DESeq2
  time_levels <- c("vP0", "vAd")
  time_levels <- time_levels[time_levels %in% unique(meta$timepoint)]

  meta_sub <- meta %>%
    filter(timepoint %in% selected_times) %>%
    mutate(timepoint = factor(timepoint, levels = time_levels)) %>%
    arrange(timepoint, samplename)

  sample_ids <- meta_sub$samplename
  # first column is the gene id
  gene_ids <- counts_df[[1]]
  counts_only <- counts_df %>%
    dplyr::select(all_of(sample_ids)) # biomaRt masks plain select()

  counts_mat <- as.matrix(counts_only)
  rownames(counts_mat) <- gene_ids
  storage.mode(counts_mat) <- "integer"

  # data.frame so rownames are ok for SummarizedExperiment
  coldata <- meta_sub %>%
    dplyr::select(samplename, timepoint) %>%
    as.data.frame()
  rownames(coldata) <- coldata$samplename

  se <- SummarizedExperiment(
    assays = list(counts = counts_mat),
    colData = coldata
  )
  se
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  dds <- DESeqDataSet(se, design = design)
  dds <- DESeq(dds)
  res_df <- results(dds)
  res_df <- as.data.frame(res_df)
  list(res_df, dds)
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  labeled <- deseq2_res %>%
    as_tibble(rownames = "genes")

  labeled <- labeled %>%
    mutate(
      volc_plot_status = case_when(
        is.na(padj) ~ "NS",
        is.na(log2FoldChange) ~ "NS",
        padj < padj_threshold & log2FoldChange > 0 ~ "UP",
        padj < padj_threshold & log2FoldChange < 0 ~ "DOWN",
        TRUE ~ "NS"
      )
    )
  labeled
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  plot_dat <- labeled_results %>%
    filter(!is.na(pvalue))

  ggplot(plot_dat, aes(x = pvalue)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white", boundary = 0) +
    theme_bw() +
    labs(
      x = "raw p-value",
      y = "number of genes",
      title = "Histogram of raw p-values (all genes)"
    )
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  sig <- labeled_results %>%
    filter(!is.na(padj), padj < padj_threshold)

  ggplot(sig, aes(x = log2FoldChange)) +
    geom_histogram(bins = 40, fill = "coral", color = "white") +
    theme_bw() +
    labs(
      x = "log2 fold change (adult vs P0)",
      y = "number of genes",
      title = "log2FC for genes significant at chosen padj cutoff"
    )
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes) {
  ranked <- labeled_results %>%
    mutate(padj_sort = ifelse(is.na(padj), Inf, padj)) %>%
    arrange(padj_sort, genes)

  top <- ranked %>%
    slice_head(n = num_genes)

  top_ids <- top$genes

  norm_mat <- counts(dds_obj, normalized = TRUE)
  sub <- norm_mat[top_ids, , drop = FALSE]

  long_df <- as.data.frame(sub)
  long_df$genes <- rownames(long_df)
  long_df <- long_df %>%
    pivot_longer(-genes, names_to = "samplename", values_to = "norm_count")

  ggplot(long_df, aes(
    x = samplename,
    y = log10(norm_count + 1),
    color = genes
  )) +
    geom_point(position = position_jitter(width = 0.15, height = 0), size = 2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "sample",
      y = "log10(normalized count + 1)",
      title = paste("Top", num_genes, "genes by lowest adjusted p-value"),
      color = "gene"
    )
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  vp <- labeled_results %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      padj_plot = pmax(padj, .Machine$double.xmin),
      neg_log10_padj = -log10(padj_plot)
    )

  ggplot(vp, aes(
    x = log2FoldChange,
    y = neg_log10_padj,
    color = volc_plot_status
  )) +
    geom_point(alpha = 0.35, size = 0.6) +
    scale_color_manual(
      values = c(UP = "#d73027", DOWN = "#313695", NS = "gray55")
    ) +
    theme_bw() +
    labs(
      x = "log2 fold change",
      y = "-log10(adjusted p-value)",
      color = "status",
      title = "Volcano plot (adult vs P0)"
    )
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  idmap <- read_tsv(
    id2gene_path,
    col_names = c("ensembl_id", "symbol"),
    show_col_types = FALSE
  )

  with_sym <- labeled_results %>%
    left_join(idmap, by = c("genes" = "ensembl_id")) %>%
    filter(!is.na(symbol))

  # average lfc if the same symbol shows up more than once
  collapsed <- with_sym %>%
    group_by(symbol) %>%
    summarise(log2FoldChange = mean(log2FoldChange, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(log2FoldChange)) %>%
    arrange(desc(log2FoldChange))

  rnk <- collapsed$log2FoldChange
  names(rnk) <- collapsed$symbol
  rnk
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  pathways <- gmtPathways(gmt_file_path)

  fg_out <- fgsea(
    pathways = pathways,
    stats = rnk_list,
    minSize = min_size,
    maxSize = max_size
  )

  fg_out <- as.data.frame(fg_out) %>%
    as_tibble()
  fg_out
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths) {
  clean <- fgsea_results %>%
    filter(!is.na(NES))

  hi <- clean %>%
    slice_max(order_by = NES, n = num_paths, with_ties = FALSE) %>%
    mutate(direction = "positive NES")

  lo <- clean %>%
    slice_min(order_by = NES, n = num_paths, with_ties = FALSE) %>%
    mutate(direction = "negative NES")

  both <- bind_rows(hi, lo) %>%
    mutate(
      lab = str_wrap(pathway, width = 45),
      nes_dir = ifelse(NES >= 0, "positive", "negative")
    )

  ggplot(both, aes(x = reorder(lab, NES), y = NES, fill = nes_dir)) +
    geom_col(width = 0.85) +
    coord_flip() +
    scale_fill_manual(values = c(positive = "#e34a33", negative = "#2c7fb8")) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(
      x = NULL,
      y = "NES",
      fill = NULL,
      title = "Top enriched / depleted pathways (fgsea)"
    )
}
