
# BINF6210 Software Tools
# Assignment 5
# Script File 1 - Preprocessing, Gene Expression and Pathway Enrichment for Biomarker Identification in R
# Isha Baxi
# Optimizing Gene Expression Analysis for Alzheimer's Disease with Machine Learning for Biomarker Identification and Validation

# Load Libraries
library(ggplot2)
library(dplyr)
library(pheatmap)
library(reshape2)
library(gridExtra)
library(DESeq2)
library(clusterProfiler)
library(edgeR)
library(org.Hs.eg.db)

# Code Section 1 – Data Acquisition, Filtering, and Quality Control -----

## Load the Data ----
# The RNA Seq was taken from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153873, from a study that looks at samples of human hippocampus from control subjects and those with Alzheimer's disease (AD).
file_path <- "../Data/GSE153873_summary_count.star.txt"  
data <- read.delim(file_path, row.names = 1, check.names = FALSE)

# Initial Dataset Dimensions
cat("Initial dataset dimensions:", dim(data), "\n")

# Generate Metadata (Condition Assignment)
condition <- sub(".*-(AD|Old|Young)", "\\1", colnames(data))
condition <- ifelse(condition %in% c("Old", "Young"), "Control", condition)
metadata <- data.frame(Sample = colnames(data), Condition = condition)
summary(metadata)

# Check: Ensure all samples are valid
valid_samples <- metadata$Condition %in% c("AD", "Control")
data <- data[, valid_samples]
metadata <- metadata[valid_samples, ]
cat("Filtered dataset dimensions (valid samples only):", dim(data), "\n")

# Define a color palette for AD vs Control to keep consistency throughout analysis
color_palette <- c("AD" = "#1f77b4", "Control" = "#ff7f0e")  # Blue for AD, Orange for Control

## Define Functions for modularity and efficiency in pre-processing ----
filter_low_counts <- function(data, threshold = 10, min_samples = 2) {
  data[rowSums(data >= threshold) >= min_samples, ]
}

filter_high_zero_counts <- function(data, zero_threshold = 80) {
  zero_counts <- rowSums(data == 0) / ncol(data) * 100
  data[zero_counts <= zero_threshold, ]
}

# Function to generate library size plots and display them side by side
generate_library_size_plots <- function(metadata) {
  
  colors <- color_palette
  
  # Internal function to create individual plots
  visualize_library_sizes <- function(metadata, colors) {
    # Bar plot of library sizes by sample
    p1 <- ggplot(metadata, aes(x = reorder(Sample, LibrarySize), y = LibrarySize, fill = Condition)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      labs(title = "Library Sizes by Sample", x = "Samples", y = "Library Size") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    # Boxplot of library sizes by condition
    p2 <- ggplot(metadata, aes(x = Condition, y = LibrarySize, fill = Condition)) +
      geom_boxplot() +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      labs(title = "Library Sizes by Condition", x = "Condition", y = "Library Size")
    
    list(sample_plot = p1, condition_plot = p2)
  }
  
  # Generate the plots
  library_size_plots <- visualize_library_sizes(metadata, colors)
  
  # Arrange and display the plots side by side
  grid.arrange(
    library_size_plots$sample_plot,
    library_size_plots$condition_plot,
    ncol = 2
  )
}

## Quality Control and Filtering: Library Size and Handling Data ----
metadata$LibrarySize <- colSums(data)
generate_library_size_plots(metadata)
# This figure provides an overview of the library sizes for each sample and their distribution across conditions (AD vs. Control). It ensures that library sizes are relatively consistent, identifying any significant discrepancies or batch effects that may impact downstream analyses and impact computational speed. 

# Handling Lowly Expressed Genes:
filtered_data <- filter_low_counts(data, threshold = 10, min_samples = 2)
cat("Filtered dataset dimensions (lowly expressed genes removed):", dim(filtered_data), "\n")

# In RNA-Seq analysis, lowly expressed genes (i.e., genes with low counts across most samples) are often excluded
# to reduce noise and improve statistical power. Genes that are expressed in fewer than 2 samples (threshold of 10 counts per sample) will be removed to focus on genes that show variation across the samples, reducing the impact of noise.

# Filtering: High Zero Counts
filtered_data <- filter_high_zero_counts(filtered_data, zero_threshold = 80)
cat("Filtered dataset dimensions (high zero count genes removed):", dim(filtered_data), "\n")

# Handling Missing Data:
missing_threshold <- 0.2  # Set threshold for missing values (e.g., 20%)
gene_missing_percentage <- rowSums(is.na(filtered_data)) / ncol(filtered_data)
filtered_data <- filtered_data[gene_missing_percentage <= missing_threshold, ]
cat("Filtered dataset dimensions (after handling missing data):", dim(filtered_data), "\n")

# In RNA-Seq data, missing values can occur due to technical issues or biological reasons. Genes with more than 20% missing values across samples will be filtered and excluded from the analysis.
# This threshold is chosen to ensure that genes with sufficient data is used. 

## Exploratory Data Analysis ----

# Log2 Transform Counts for Data Normalization
log_counts <- log2(filtered_data + 1)

# Density plot of gene expression
# The density plot for gene expression provides a high-level overview of the distribution of expression values across a dataset, helping to assess data quality, detect outliers, and identify technical biases. It can be a good indicator of evaluating the effectiveness of normalization techniques and ensuring consistency between experimental groups. 
melted_data <- melt(as.data.frame(log_counts))
ggplot(melted_data, aes(value)) +
  geom_density(fill = "lightblue", alpha = 0.6) +
  theme_minimal() +
  labs(title = "Density Plot of Gene Expression (Log2 Transformed)", x = "Log2 Counts", y = "Density")

# PCA: Check for clustering by condition
pca <- prcomp(t(log_counts))
pca_data <- data.frame(Sample = metadata$Sample,
                       PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                       Condition = metadata$Condition)

ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  theme_minimal() +
  labs(title = "PCA of Samples", x = "Principal Component 1", y = "Principal Component 2")

# Outlier Detection
mad_scores <- apply(log_counts, 2, mad)
boxplot(mad_scores, main = "Sample-wise MAD Scores", ylab = "MAD",
        col = "gray", border = "black", horizontal = TRUE)

mad_threshold <- mean(mad_scores) + 2 * sd(mad_scores)
# MAD scores detect outliers by measuring the spread of data around the median, making them more robust to extreme values than methods like standard deviation. High MAD scores indicate data points that deviate significantly from the central tendency, helping to identify potential outliers, especially in skewed datasets.
outliers <- names(which(mad_scores > mad_threshold))
cat("Potential outliers detected:", outliers, "\n")

summary(filtered_data)
summary(filtered_data$`4-13A-Young`)
summary(filtered_data$`8-15A-Young`)

# Highlight outliers in PCA plot
pca_data$Outlier <- ifelse(pca_data$Sample %in% outliers, "Outlier", "Normal")
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, shape = Outlier)) +
  geom_point(size = 4) +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = c("Normal" = 16, "Outlier" = 17)) +
  theme_minimal() +
  labs(title = "PCA of Samples with Outliers Highlighted",
       x = "Principal Component 1", y = "Principal Component 2")

# This PCA plot captures variability in the dataset, effectively differentiating AD and Control samples while highlighting outliers. This is useful for identifying patterns in the data and detecting samples that may need to be excluded to improve data quality. 
# For this analysis, both samples were kept as the whole datset is relatively small, but it is to be known that the dataset may be slightly skewed due to this. This acts as a limitation in the analysis. 

# Heatmap of Top 50 Most Variable Genes
variances <- apply(log_counts, 1, var)
top_genes <- names(sort(variances, decreasing = TRUE))[1:50]

annotation_col <- metadata
rownames(annotation_col) <- metadata$Sample
log_counts <- log_counts[, rownames(annotation_col)]  # Ensure match

# Define colors for the annotation
annotation_colors <- list(
  Condition = c("AD" = "#1f77b4", "Control" = "#ff7f0e")  # Blue for AD, Orange for Control
)

# Heatmap of Top 50 Most Variable Genes
pheatmap(
  log_counts[top_genes, ],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col[, "Condition", drop = FALSE],
  annotation_colors = annotation_colors,  # Apply the custom colors
  main = "Heatmap of Top 50 Most Variable Genes",
  color = colorRampPalette(c("#edf8fb", "#2b8cbe"))(50)  # Blue gradient for expression
)
# This heatmap provides an initial look at the genes driving variability in the dataset, which may include biomarkers or genes associated with disease pathology. The clustering of samples by condition (AD vs. Control) suggests the presence of biologically relevant differences. Moreover, hierarchical clustering of genes can identify groups of co-expressed genes, which may reveal pathways or gene networks of interest for further analysis.

summary(metadata)

# Save Filtered Data and Metadata
write.csv(filtered_data, "filtered_gene_expression.csv")
write.csv(metadata, "metadata.csv")

# Final Dataset Dimensions
cat("Final dataset dimensions (after filtering):", dim(filtered_data), "\n")

# Code 2 - Main Analysis ------ 

# Initialization and Metadata
metadata_norm <- data.frame(row.names = colnames(filtered_data), condition = condition)

## Functions ---- 
# Normalizations Functions using DSeq2 and EdgeR 
# This function normalizes RNA-seq count data using DESeq2's default median ratio method for size factor estimation, which adjusts for sequencing depth across samples. It also uses the design formula "condition" by default to account for the experimental condition during normalization.
normalize_deseq2 <- function(count_data, metadata, design_formula = "~ condition") {
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = as.formula(design_formula))
  
  time_deseq2 <- system.time({
    dds <- estimateSizeFactors(dds)  # Median ratio method for size factor estimation
    normalized_counts <- counts(dds, normalized = TRUE)  # Get normalized counts
  })
  
  return(list(normalized_counts = normalized_counts, time = time_deseq2["elapsed"], dds = dds))
}

# The function normalizes count data using edgeR's default TMM (Trimmed Mean of M-values) method, which accounts for sequencing depth and composition biases. The normalized data is output as Counts Per Million (CPM), a standard metric for RNA-seq data normalization.
normalize_edger <- function(count_data, metadata, method = "TMM") {
  dge <- DGEList(counts = count_data, group = metadata$condition)
  
  time_edger <- system.time({
    dge <- calcNormFactors(dge, method = method)  # Apply TMM normalization
    normalized_counts <- cpm(dge)  # Get counts per million (CPM)
  })
  
  return(list(normalized_counts = normalized_counts, time = time_edger["elapsed"], dge = dge))
}

# Visualization Functions - Normalization Comparison Plots and PCA of Normalizations
visualize_normalizations <- function(raw_data, deseq2_data, edger_data) {
  par(mfrow = c(1, 3))
  
  boxplot(log2(raw_data + 1), main = "Raw Counts", col = "skyblue", las = 2)
  boxplot(log2(deseq2_data + 1), main = "DESeq2 Normalized", col = "lightgreen", las = 2)
  boxplot(log2(edger_data + 1), main = "edgeR Normalized", col = "orange", las = 2)
  
  par(mfrow = c(1, 1)) # Reset layout
}

visualize_pca <- function(normalized_data, metadata, title = "PCA Plot") {
  # PCA calculation
  pca <- prcomp(t(log2(normalized_data + 1)))
  pca_data <- as.data.frame(pca$x)
  pca_data$Condition <- metadata$condition
  
  # Plot
  ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 4) +
    labs(title = title, x = "PC1", y = "PC2") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}

# Differential Expression Functions
# This function performs differential expression analysis using DESeq2, applying a specified contrast between conditions (e.g., AD vs. Control) and returns the significant genes based on an adjusted p-value threshold (default 0.05).
run_deseq2 <- function(dds, contrast = c("condition", "AD", "Control"), pvalue_threshold = 0.05) {
  dds <- DESeq(dds)
  res <- results(dds, contrast = contrast)
  sig_genes <- subset(as.data.frame(res), padj < pvalue_threshold & !is.na(padj))
  return(sig_genes)
}

# This function performs differential expression analysis using edgeR, fitting a negative binomial model to the count data, and returns significant genes based on a false discovery rate (FDR) threshold (default 0.05).
run_edger <- function(dge, metadata, coef = 2, fdr_threshold = 0.05) {
  design <- model.matrix(~ condition, data = metadata)
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, coef = coef)
  
  res <- topTags(lrt, n = Inf)$table
  sig_genes <- subset(res, FDR < fdr_threshold)
  return(sig_genes)
}

# This function performs KEGG pathway enrichment analysis by mapping significant gene symbols to Entrez IDs, then testing for overrepresented pathways with a p-value cutoff of 0.05.
pathway_enrichment <- function(sig_genes) {
  # Ensure sig_genes is non-empty
  if (nrow(sig_genes) == 0) {
    message("No significant genes provided.")
    return(NULL)
  }
  
  # Extract gene symbols and map to Entrez IDs
  gene_list <- rownames(sig_genes)
  gene_ids <- mapIds(
    org.Hs.eg.db, 
    keys = gene_list, 
    column = "ENTREZID", 
    keytype = "SYMBOL", 
    multiVals = "first"
  )
  
  # Remove NA and duplicate Entrez IDs
  gene_ids <- na.omit(unique(gene_ids))
  
  # Perform KEGG pathway enrichment analysis
  enrich_result <- enrichKEGG(
    gene = gene_ids, 
    organism = 'hsa', 
    pvalueCutoff = 0.05
  )
  
  return(enrich_result)
}

# This function generates a bar plot of significant pathways from the KEGG enrichment results, displaying the top pathways based on the specified number of categories (default 10).
visualize_pathway_enrichment <- function(enrich_result, title, showCategory = 10) {
  # Ensure enrich_result is not null and has rows
  if (!is.null(enrich_result) && nrow(as.data.frame(enrich_result)) > 0) {
    # Generate the barplot directly
    barplot(
      enrich_result,
      showCategory = showCategory
    ) +
      ggtitle(title) +
      theme_minimal()
  } else {
    message("No significant pathways found.")
  }
}

# Pathway enrichment analysis is significant because it helps identify biological pathways that are disproportionately associated with a set of genes of interest, often revealing underlying molecular mechanisms or disease processes. For this project, it can pathways relevant to AD, offering insights into potential biomarkers or therapeutic targets for detection. 

# Accuracy Metric for Normalization
# This function computes the coefficient of variation (CV) for each gene across samples in the `normalized_data`, which is defined as the ratio of the standard deviation to the mean. It returns the mean of these CV values across all genes, providing an overall measure of variability in the dataset. A higher CV = greater variability among the samples for that gene, which could be useful for identifying genes with high or low expression variability in the dataset.
calculate_accuracy <- function(normalized_data) {
  cv <- apply(normalized_data, 1, function(x) sd(x) / mean(x))
  return(mean(cv, na.rm = TRUE))
}

compare_methods <- function(deseq2_time, edger_time, deseq2_data, edger_data) {
  deseq2_accuracy <- calculate_accuracy(deseq2_data)
  edger_accuracy <- calculate_accuracy(edger_data)
  
  comparison <- data.frame(
    Method = c("DESeq2", "edgeR"),
    Computation_Time_Seconds = c(deseq2_time, edger_time),
    Mean_Coefficient_of_Variation = c(deseq2_accuracy, edger_accuracy)
  )
  print(comparison)
}

## Main Workflow and Implementation of Functions ----

# Normalize data using DESeq2 and EdgeR
deseq2_result <- normalize_deseq2(filtered_data, metadata_norm)
edger_result <- normalize_edger(filtered_data, metadata_norm)

# Visualizations for Normalizations
visualize_normalizations(filtered_data, deseq2_result$normalized_counts, edger_result$normalized_counts)
# The figure demonstrates how normalization reduces variability and adjusts for library size differences, improving the reliability of downstream analyses. In the Raw Counts panel, the variability is higher, with uneven medians, indicating potential biases due to differing sequencing depths. The DESeq2 Normalized and edgeR Normalized panels show more consistent medians across samples, demonstrating successful normalization. This consistency is critical for reducing technical noise and ensuring the biological signal drives subsequent analyses. 

visualize_pca(deseq2_result$normalized_counts, metadata_norm, title = "PCA Plot (DESeq2 Normalized)")
visualize_pca(edger_result$normalized_counts, metadata_norm, title = "PCA Plot (EdgeR Normalized)")

# Conduct Differential Expression
sig_genes_deseq2 <- run_deseq2(deseq2_result$dds)
sig_genes_edger <- run_edger(edger_result$dge, metadata_norm)

# Extract gene names from both results, and find the common genes to use for biomarker identification and validation
genes_deseq2 <- rownames(sig_genes_deseq2)
genes_edger <- rownames(sig_genes_edger)
common_genes <- intersect(genes_deseq2, genes_edger)
cat("Number of common genes: ", length(common_genes), "\n")

## Pathway Enrichment ----
enrich_deseq2 <- pathway_enrichment(sig_genes_deseq2)
enrich_edger <- pathway_enrichment(sig_genes_edger)

visualize_pathway_enrichment(enrich_deseq2, title = "Top Pathways (DESeq2)")
# This bar plot illustrates the top pathways enriched in the Alzheimer's dataset based on DESeq2-normalized gene expression data. Pathways like "Alzheimer disease," "Parkinson disease," and "Pathways of neurodegeneration - multiple diseases" show significant enrichment, highlighting their strong association with neurodegenerative processes. The color gradient corresponds to adjusted p-values, with darker red representing higher significance, emphasizing the importance of these pathways in the context of Alzheimer's disease research.

visualize_pathway_enrichment(enrich_edger, title = "Top Pathways (edgeR)")
# Similar to the previous plot, this one displays the top pathways using the edgeR package. Both display the same diseases and processes, but the order and p-values differ, indicating and displaying the different normalization algorithms (median ratio method and TMM) effects. 

## Compare DESeq2 and EdgeR -----
# DESeq2 is computationally faster (0.384 seconds) than edgeR (0.888 seconds), while edgeR achieves a slightly lower mean coefficient of variation (0.453 vs. 0.456), indicating greater consistency in normalized values. This trade-off suggests that DESeq2 is preferable for speed, whereas edgeR might be slightly better for stability in gene expression normalization.
compare_methods(deseq2_result$time, edger_result$time, deseq2_result$normalized_counts, edger_result$normalized_counts)

# Prepare output file for Python ML Implementation
common_gene_data <- sig_genes_deseq2[common_genes, ]
ml_input <- t(common_gene_data) # Transpose the df to have samples as rows and genes as columns
ml_input <- data.frame(ml_input, Condition = metadata_norm$condition) # Add condition labels to the df

# Save Results as CSV to be used. All files can be found in the Data folder. ----
#write.csv(sig_genes_deseq2, "DEG_results_DESeq2.csv")
#write.csv(sig_genes_edger, "DEG_results_edgeR.csv")
#write.csv(ml_input, "ml_input_common_genes.csv", row.names = TRUE)
