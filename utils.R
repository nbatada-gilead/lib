library(ggplot2)
library(dplyr)    # alternatively, this also loads %>%


seurat_view_metadata <-function(sobj){
    lapply(sobj@meta.data[, sapply(sobj@meta.data, is.character)], unique) 
}


whos <- function() {
    # like matlab: list variables and their sizes
  object_names <- ls(envir = .GlobalEnv)
  
  object_sizes <- sapply(object_names, function(x) {
    size_in_bytes <- object.size(get(x, envir = .GlobalEnv))
    size_in_MB <- size_in_bytes / (1024^2)  # Convert bytes to megabytes
    round(size_in_MB, 0)  # Round to three decimal places for cleaner output
  }, simplify = "vector")
  
  named_sizes <- setNames(object_sizes, object_names)
  
  sorted_sizes <- sort(named_sizes, decreasing = TRUE)
  
  print(sorted_sizes)
}

rename_genes_ensembl_to_hgnc <- function(sobj){
    filename = "~/data/geneinfo/genemaps/genesymbol_to_ensgeneid/hgnc/genesymbol2geneid_hgnc.map"
   # Load the gene symbol to Ensembl ID mapping file
   gene_map <- read.table(filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

   # Extract Ensembl IDs and Gene Symbols from the file
   ensembl_ids_in_file <- gene_map$ensembl_gene_id
   gene_symbols_in_file <- setNames(gene_map$symbol, gene_map$ensembl_gene_id)

   # Extract Ensembl IDs from the Seurat object
   ensembl_ids_in_sobj <- rownames(sobj)

   # Filter Ensembl IDs to keep only those present in the local file
   valid_ensembl_ids <- ensembl_ids_in_sobj[ensembl_ids_in_sobj %in% ensembl_ids_in_file]

   # If no valid Ensembl IDs are found, stop with a message
   if (length(valid_ensembl_ids) == 0) {
      stop("No matching Ensembl IDs found in the local file.")
   }

   # Replace Ensembl IDs with Gene Symbols using the file
   gene_symbols <- gene_symbols_in_file[valid_ensembl_ids]

   # Ensure uniqueness of rownames after replacement
   unique_gene_symbols <- make.unique(gene_symbols)

   # Check if fewer gene symbols than Ensembl IDs
   if (length(unique_gene_symbols) < length(valid_ensembl_ids)) {
      message("Creating a new Seurat object with reduced gene symbols...")

      # Subset original Seurat object
      sobj <- sobj[valid_ensembl_ids, ]

      # Create a new Seurat object based on the available counts or data
      if (!is.null(SeuratObject::GetAssayData(sobj, slot = "counts"))) {
         new_counts <- SeuratObject::GetAssayData(sobj, slot = "counts")
         rownames(new_counts) <- unique_gene_symbols
      } else {
         new_counts <- NULL
         warning("Layer 'counts' isn't present, skipping 'counts'.")
      }

      if (!is.null(SeuratObject::GetAssayData(sobj, slot = "data"))) {
         new_data <- SeuratObject::GetAssayData(sobj, slot = "data")
         rownames(new_data) <- unique_gene_symbols
      } else {
         new_data <- NULL
         warning("Layer 'data' isn't present, skipping 'data'.")
      }

      # Initialize a new Seurat object
      new_sobj <- CreateSeuratObject(counts = new_counts, data = new_data)

      # Copy over other information like metadata and dimensional reductions
      new_sobj@meta.data <- sobj@meta.data

      # Copy over dimensional reductions if they exist
      for (reduction in names(sobj@reductions)) {
         new_sobj@reductions[[reduction]] <- sobj@reductions[[reduction]]
      }

      # Copy over meta.features if available
      if ("meta.features" %in% slotNames(sobj@assays$RNA)) {
         new_sobj@assays$RNA@meta.features <- sobj@assays$RNA@meta.features[valid_ensembl_ids, ]
      }

      # Return the new Seurat object
      return(new_sobj)
   } else {
      # Update rownames in counts and data layers if no new object is created
      if (!is.null(SeuratObject::GetAssayData(sobj, slot = "counts"))) {
         counts_layer <- SeuratObject::GetAssayData(sobj, slot = "counts")
         if (nrow(counts_layer) == length(unique_gene_symbols)) {
            rownames(counts_layer) <- unique_gene_symbols
            sobj <- SeuratObject::SetAssayData(sobj, slot = "counts", new.data = counts_layer)
         } else {
            stop("Length of unique gene symbols does not match RNA counts matrix.")
         }
      }

      if (!is.null(SeuratObject::GetAssayData(sobj, slot = "data"))) {
         data_layer <- SeuratObject::GetAssayData(sobj, slot = "data")
         if (nrow(data_layer) == length(unique_gene_symbols)) {
            rownames(data_layer) <- unique_gene_symbols
            sobj <- SeuratObject::SetAssayData(sobj, slot = "data", new.data = data_layer)
         } else {
            stop("Length of unique gene symbols does not match RNA data matrix.")
         }
      }

      return(sobj)
   }
}


rename_genes_ensembl_to_hgnc_remote <- function(sobj, species = "human", mirror = "useast") {
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

   # Ensure uniqueness by taking the first occurrence of each HGNC symbol
   unique_indices <- !duplicated(gene_symbols)
   unique_gene_symbols <- gene_symbols[unique_indices]

   # Subset the Seurat object to include only the rows that correspond to the unique gene symbols
   sobj <- sobj[unique_indices, ]

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


rename_genes_ensembl_to_hgnc_old <- function(sobj, species = "human", mirror = "useast") {
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

