# RNAseq_DESeq2.R
library(DESeq2)
library(tidyverse)

# Set working directory
setwd("/projectnb/bf528/students/addisony/ATACseq-and-Differential-Chromatin-Accessibility-Analysis-Pipeline/")

# Read your RNA-seq raw counts
# Assuming you have files: cDC1_Raw_counts.tsv, cDC2_Raw_counts.tsv
cDC1_counts <- read.table("rnaseq_counts/cDC1_Raw_counts.tsv", 
                          header = TRUE, sep = "\t", row.names = 1)
cDC2_counts <- read.table("rnaseq_counts/cDC2_Raw_counts.tsv", 
                          header = TRUE, sep = "\t", row.names = 1)

# Create sample information
# Assuming your columns are: WT_1, WT_2, KO_1, KO_2
create_deseq_dataset <- function(counts_df, cell_type) {
  # Extract sample names from column names
  sample_names <- colnames(counts_df)
  
  # Create sample info
  sample_info <- data.frame(
    sample = sample_names,
    condition = ifelse(grepl("WT", sample_names), "WT", "KO"),
    cell_type = cell_type
  )
  rownames(sample_info) <- sample_names
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = counts_df,
    colData = sample_info,
    design = ~ condition
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results
  res <- results(dds, contrast = c("condition", "KO", "WT"))
  
  # Convert to dataframe
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Save results
  write.csv(res_df, paste0("rnaseq_results/", cell_type, "_DESeq2_results.csv"))
  
  # Also save normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  write.csv(norm_counts, paste0("rnaseq_results/", cell_type, "_normalized_counts.csv"))
  
  return(res_df)
}

# Create output directory
dir.create("rnaseq_results", showWarnings = FALSE)

# Run for both cell types
cDC1_res <- create_deseq_dataset(cDC1_counts, "cDC1")
cDC2_res <- create_deseq_dataset(cDC2_counts, "cDC2")

# Summary
cat("cDC1 results:\n")
cat("  Total genes:", nrow(cDC1_res), "\n")
cat("  Significant (p < 0.05):", sum(cDC1_res$padj < 0.05, na.rm = TRUE), "\n")
cat("  Upregulated (log2FC > 0):", sum(cDC1_res$log2FoldChange > 0 & cDC1_res$padj < 0.05, na.rm = TRUE), "\n")
cat("  Downregulated (log2FC < 0):", sum(cDC1_res$log2FoldChange < 0 & cDC1_res$padj < 0.05, na.rm = TRUE), "\n\n")

cat("cDC2 results:\n")
cat("  Total genes:", nrow(cDC2_res), "\n")
cat("  Significant (p < 0.05):", sum(cDC2_res$padj < 0.05, na.rm = TRUE), "\n")
cat("  Upregulated (log2FC > 0):", sum(cDC2_res$log2FoldChange > 0 & cDC2_res$padj < 0.05, na.rm = TRUE), "\n")
cat("  Downregulated (log2FC < 0):", sum(cDC2_res$log2FoldChange < 0 & cDC2_res$padj < 0.05, na.rm = TRUE), "\n")