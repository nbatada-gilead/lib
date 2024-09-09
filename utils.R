library(ggplot2)
library(dplyr)    # alternatively, this also loads %>%



rename_genes_ensembl_to_hgnc <- function(sobj, species = "human", mirror = "useast") {
   # usage: 
   # sobj <- rename_genes_ensembl_to_hgnc(sobj, species = "human", mirror = "useast")
   
   library(biomaRt)
   
   # Choose the correct dataset and use a mirror
   if (species == "human") {
      ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = mirror)
   } else if (species == "mouse") {
      ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", mirror = mirror)
   } else {
      stop("Unsupported species. Use 'human' or 'mouse'.")
   }
   
   # Extract Ensembl IDs
   ensembl_ids <- rownames(sobj)
   
   # Query biomaRt to get HGNC symbols for the Ensembl IDs
   gene_info <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = "ensembl_gene_id",
      values = ensembl_ids,
      mart = ensembl
   )
   
   # Create a named vector for mapping
   ensembl_to_hgnc <- setNames(gene_info$hgnc_symbol, gene_info$ensembl_gene_id)
   
   # Replace Ensembl IDs with HGNC symbols, keeping the original Ensembl ID if no symbol is found
   gene_symbols <- ensembl_to_hgnc[ensembl_ids]
   gene_symbols[is.na(gene_symbols)] <- ensembl_ids[is.na(gene_symbols)]
   
   # Ensure uniqueness of rownames after replacement
   unique_gene_symbols <- make.unique(gene_symbols)
   
   # Update rownames in counts and data slots
   rownames(sobj@assays$RNA@counts) <- unique_gene_symbols
   rownames(sobj@assays$RNA@data) <- unique_gene_symbols
   
   # Ensure metadata consistency if present
   if ("meta.features" %in% slotNames(sobj@assays$RNA)) {
      sobj@assays$RNA@meta.features <- sobj@assays$RNA@meta.features[unique_gene_symbols, , drop = FALSE]
   }
   
   # Return the updated Seurat object
   return(sobj)
}


stats_check_if_gene_is_de <- function(df_expr, df_meta, gene, condition) {
   
   # Check if the gene exists in df_expr
   if (!gene %in% rownames(df_expr)) {
      stop(sprintf("Gene '%s' not found in the expression matrix.", gene))
   }
   
   # Check if the condition exists in df_meta
   if (!condition %in% colnames(df_meta)) {
      stop(sprintf("Condition '%s' not found in the metadata.", condition))
   }
   
   # Ensure there are matching samples in both df_expr and df_meta
   matching_samples <- intersect(colnames(df_expr), rownames(df_meta))
   gene_counts <- df_expr[gene, matching_samples]
   
   # Extract the condition values
   condition_values <- df_meta[matching_samples, condition]
   
   # Check if condition has at least two unique levels
   if (length(unique(condition_values)) < 2) {
      stop(sprintf("Condition '%s' does not have at least two unique levels.", condition))
   }
   
   # Perform the Wilcoxon test
   group1 <- gene_counts[condition_values == unique(condition_values)[1]]
   group2 <- gene_counts[condition_values == unique(condition_values)[2]]
   
   wilcox_result <- wilcox.test(group1, group2)
   sprintf("%s (condition: %s) p.value = %s", gene,condition,wilcox_result$p.value)
   
   # Return the p-value
   return(wilcox_result)
}

plot_expr_of_gene_by_condition <- function(df_expr, df_meta, gene, condition) {
   
   # Check if the gene exists in df_expr
   if (!gene %in% rownames(df_expr)) {
      stop(sprintf("--\n<gene> '%s' not found in the expression matrix.\n--", gene))
   }
   
   # Check if the condition exists in df_meta
   if (!condition %in% colnames(df_meta)) {
      stop(sprintf("--\n<condition> '%s' not found in the metadata.\n--", condition))
   }
   
  
  # Extract gene expression values
  gene_counts <- df_expr[gene, ]

  # Ensure there are matching samples in both df_expr and df_meta
  matching_samples <- intersect(colnames(df_expr), rownames(df_meta))
  gene_counts <- gene_counts[matching_samples]
  
  # Create a data frame for plotting
  plot_data <- data.frame(
    Sample = matching_samples,
    Expression = gene_counts,
    Condition = df_meta[matching_samples, condition]
  )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.7) +
    labs(title = sprintf("%s Expression by %s", gene, condition),
         x = condition,
         y = "Normalized Expression") +
    theme_minimal()

  wilcox_result <- stats_check_if_gene_is_de(df_expr, df_meta, gene, condition)  
    
  # Add p-value to the plot
  p <- p + 
    annotate("text", x = 1.5, y = max(plot_data$Expression, na.rm = TRUE), 
             label = sprintf("p = %.4f", wilcox_result$p.value), 
             vjust = -0.5)
  
  print(p)
  
  # Return the p-value as well
  return(wilcox_result$p.value)
}

library(ggplot2)
