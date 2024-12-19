# scrds_utils.R
#
# Wrapper to work with scRDS ingestion specifically for Abio annotated data
# Nizar Batada
# Version 1
# 12-12-2024

library(scRDS)
library(Seurat)
library(dplyr)
library(stringr)

validate_meta_data <- function(sobj) {
  # Define required columns
  required_columns <- c("nCount_RNA", "nFeature_RNA", "percent.mito")
  flexible_columns <- c("condition", "fraction", "clusters")
  
  # Generate column mapping
  column_mapping <- setNames(paste0("column_", seq_along(colnames(sobj@meta.data))), colnames(sobj@meta.data))
  
  # Log user decisions for reproducibility
  reproducibility_log <- list()
  
  # Check for fixed required columns
  for (required_col in required_columns) {
    if (!(required_col %in% colnames(sobj@meta.data))) {
      # Suggest alternatives for percent.mito
      if (required_col == "percent.mito") {
        alternatives <- grep("percent.mt|percent.mito", colnames(sobj@meta.data), value = TRUE, ignore.case = TRUE)
        cat(paste("\nThe required column", required_col, "is missing.\n"))
        
        # Display equivalent columns with their unique values
        cat("Available columns in meta.data and their top 5 most frequent or unique values:\n")
        for (col in colnames(sobj@meta.data)) {
          unique_vals <- unique(sobj@meta.data[[col]])
          top_5_vals <- unique_vals[1:min(5, length(unique_vals))]
          cat(paste0(column_mapping[col], " (", col, "): ", paste(top_5_vals, collapse = ", "), " (top 5)\n"))
        }
        
        # Show detected alternatives
        if (length(alternatives) > 0) {
          cat("Alternative column(s) detected:\n")
          for (i in seq_along(alternatives)) {
            cat(paste0(i, ": ", alternatives[i], "\n"))
          }
        } else {
          cat("No alternative column(s) detected.\n")
        }
        
        # Prompt user input
        user_input <- readline(prompt = "Enter the number of the alternative column to use, its name, or 'none' to create and populate with a default value: ")
        
        # Handle user input
        if (tolower(user_input) == "none") {
          default_value <- as.numeric(readline(prompt = paste("Enter the default value to initialize", required_col, ": ")))
          sobj@meta.data[[required_col]] <- default_value
          reproducibility_log[[required_col]] <- paste0("sobj@meta.data[['", required_col, "']] <- ", default_value)
          cat(paste("Initialized", required_col, "with default value:", default_value, "\n"))
        } else if (suppressWarnings(as.numeric(user_input)) %in% seq_along(alternatives)) {
          selected_index <- as.numeric(user_input)
          selected_col <- alternatives[selected_index]
          sobj@meta.data[[required_col]] <- sobj@meta.data[[selected_col]]
          reproducibility_log[[required_col]] <- paste0("sobj@meta.data[['", required_col, "']] <- sobj@meta.data[['", selected_col, "']]")
          cat(paste("\nUsing column", selected_col, "for", required_col, ".\n"))
        } else if (user_input %in% colnames(sobj@meta.data)) {
          sobj@meta.data[[required_col]] <- sobj@meta.data[[user_input]]
          reproducibility_log[[required_col]] <- paste0("sobj@meta.data[['", required_col, "']] <- sobj@meta.data[['", user_input, "']]")
          cat(paste("\nUsing column", user_input, "for", required_col, ".\n"))
        } else if (user_input %in% column_mapping) {
          original_col <- names(column_mapping[column_mapping == user_input])
          sobj@meta.data[[required_col]] <- sobj@meta.data[[original_col]]
          reproducibility_log[[required_col]] <- paste0("sobj@meta.data[['", required_col, "']] <- sobj@meta.data[['", original_col, "']]")
          cat(paste("\nUsing column", original_col, "for", required_col, ".\n"))
        } else {
          stop("Invalid selection. Please ensure you enter a valid number, column name, or 'none'.")
        }
      } else {
        stop(paste("The required column", required_col, "is missing in meta.data."))
      }
    } else {
      cat(paste("Column", required_col, "is present in meta.data.\n"))
    }
  }
  
  # Flexible column detection and user interaction
  for (flex_col in flexible_columns) {
    detected <- NULL
    
    # Display description and meta.data columns
    cat("\n----\n")
    cat(paste("Specify the column for", flex_col, ".\n"))
    cat("----\n")
    cat("Available columns in meta.data and their top 5 most frequent or unique values:\n")
    for (col in colnames(sobj@meta.data)) {
      unique_vals <- unique(sobj@meta.data[[col]])
      top_5_vals <- unique_vals[1:min(5, length(unique_vals))]
      cat(paste0(column_mapping[col], " (", col, "): ", paste(top_5_vals, collapse = ", "), " (top 5)\n"))
    }
    
    # Handle user input
    column_name <- readline(prompt = paste("Specify the column name, mapped name, or enter 'none' to create and populate a default column for", flex_col, ": "))
    if (tolower(column_name) == "none") {
      sobj@meta.data[[flex_col]] <- paste("unknown", flex_col, sep = "_")
      reproducibility_log[[flex_col]] <- paste0("sobj@meta.data[['", flex_col, "']] <- 'unknown_", flex_col, "'")
      cat(paste("Initialized", flex_col, "with default value: unknown_", flex_col, "\n", sep = ""))
    } else if (column_name %in% colnames(sobj@meta.data)) {
      sobj@meta.data[[flex_col]] <- sobj@meta.data[[column_name]]
      reproducibility_log[[flex_col]] <- paste0("sobj@meta.data[['", flex_col, "']] <- sobj@meta.data[['", column_name, "']]")
      cat(paste("\nUsing column", column_name, "for", flex_col, ".\n"))
    } else if (column_name %in% column_mapping) {
      original_col <- names(column_mapping[column_mapping == column_name])
      sobj@meta.data[[flex_col]] <- sobj@meta.data[[original_col]]
      reproducibility_log[[flex_col]] <- paste0("sobj@meta.data[['", flex_col, "']] <- sobj@meta.data[['", original_col, "']]")
      cat(paste("\nUsing column", original_col, "for", flex_col, ".\n"))
    } else {
      stop("Invalid selection. Please ensure you enter a valid column name or 'none'.")
    }
  }
  
  # Print reproducibility log
  cat("\n----\nReproducibility Code:\n----\n")
  for (log in reproducibility_log) {
    cat(log, "\n")
  }
  
  cat("\nAll required columns are now present and validated in meta.data.\n")
  return(sobj)
}

validate_misc_experiment <- function(sobj) {
  # Ensure sobj@misc exists
  if (is.null(sobj@misc)) {
    sobj@misc <- list()
  }
  
  # Ensure sobj@misc$experiment exists
  if (is.null(sobj@misc$experiment)) {
    sobj@misc$experiment <- list()
  }
  
  # Define required fields in experiment
  required_fields <- c(
    "experiment", "study", "type", "class", "order", "species", 
    "tissue", "disease", "platform", "cell.selection", "info", 
    "reference", "ExpressionMatrixType", "link", "contact", 
    "description", "notes"
  )
  
  # Log for reproducibility
  reproducibility_log <- list()
  
  # Iterate through required fields
  for (field in required_fields) {
    if (is.null(sobj@misc$experiment[[field]])) {
      # Field is missing; handle based on type
      if (field %in% c("experiment", "study", "type", "class", "species", "tissue", "disease", 
                       "platform", "cell.selection", "ExpressionMatrixType", "reference", 
                       "notes", "link", "contact", "description")) {
        # Ask user to manually specify a string
        user_input <- readline(prompt = paste("Enter a value for", field, "(string): "))
        sobj@misc$experiment[[field]] <- user_input
        reproducibility_log[[field]] <- paste0("sobj@misc$experiment[['", field, "']] <- '", user_input, "'")
        cat(paste("\nSet", field, "to:", user_input, "\n"))
      
      } else if (field == "info") {
        # Suggest and calculate default value for 'info'
        default_info <- paste(
          nrow(sobj), "cells from",
          paste(unique(sobj@meta.data$condition), collapse = ", ")
        )
        user_input <- readline(prompt = paste("Enter a value for", field, "(or press Enter to use default: '", default_info, "'): "))
        if (user_input == "") {
          user_input <- default_info
        }
        sobj@misc$experiment[[field]] <- user_input
        reproducibility_log[[field]] <- paste0("sobj@misc$experiment[['", field, "']] <- '", user_input, "'")
        cat(paste("\nSet", field, "to:", user_input, "\n"))
      
      } else if (field == "order") {
        # Ask user to specify a numeric value
        user_input <- as.numeric(readline(prompt = paste("Enter a numeric value for", field, ": ")))
        sobj@misc$experiment[[field]] <- user_input
        reproducibility_log[[field]] <- paste0("sobj@misc$experiment[['", field, "']] <- ", user_input)
        cat(paste("\nSet", field, "to:", user_input, "\n"))
      }
    } else {
      cat(paste("Field", field, "is already present in misc$experiment.\n"))
    }
  }
  
  # Print reproducibility log
  cat("\n----\nReproducibility Code:\n----\n")
  for (log in reproducibility_log) {
    cat(log, "\n")
  }
  
  cat("\nAll required fields are now present in misc$experiment.\n")
  return(sobj)
}

validate_misc_celltypes <- function(sobj) {
  # Ensure sobj@misc exists
  if (is.null(sobj@misc)) {
    sobj@misc <- list()
  }
  
  # Ensure sobj@misc$celltypes exists
  if (is.null(sobj@misc$celltypes)) {
    sobj@misc$celltypes <- list()
  }
  
  # Define required fields in celltypes
  required_fields <- c("category", "celltype", "description", "order")
  optional_fields <- c("marker")
  
  # Log for reproducibility
  reproducibility_log <- list()
  
  # Helper: Display meta.data columns and their top 5 unique values
  display_meta_data <- function() {
    column_mapping <- setNames(paste0("column_", seq_along(colnames(sobj@meta.data))), colnames(sobj@meta.data))
    cat("Available columns in meta.data and their top 5 unique values:\n")
    for (col in colnames(sobj@meta.data)) {
      unique_vals <- unique(sobj@meta.data[[col]])
      top_5_vals <- unique_vals[1:min(5, length(unique_vals))]
      cat(paste0(column_mapping[col], " (", col, "): ", paste(top_5_vals, collapse = ", "), " (top 5)\n"))
    }
  }
  
  # Iterate through required fields
  for (field in required_fields) {
    if (is.null(sobj@misc$celltypes[[field]])) {
      if (field == "category") {
        # category must match a column in meta.data
        cat("\nThe required field 'category' is missing in misc$celltypes.\n")
        cat("Annotation type (e.g., 'Cell Type (Internal)'). Must match a column in meta.data.\n")
        display_meta_data()
        
        user_input <- readline(prompt = "Specify the column name or mapped name for 'category': ")
        if (user_input %in% colnames(sobj@meta.data)) {
          sobj@misc$celltypes[[field]] <- user_input
          reproducibility_log[[field]] <- paste0("sobj@misc$celltypes[['", field, "']] <- '", user_input, "'")
          cat(paste("\nSet 'category' to:", user_input, "\n"))
        } else if (user_input %in% names(column_mapping)) {
          original_col <- names(column_mapping[column_mapping == user_input])
          sobj@misc$celltypes[[field]] <- original_col
          reproducibility_log[[field]] <- paste0("sobj@misc$celltypes[['", field, "']] <- '", original_col, "'")
          cat(paste("\nSet 'category' to:", original_col, "\n"))
        } else {
          stop("Invalid selection for 'category'. Please ensure it matches a column in meta.data.")
        }
      } else if (field == "celltype") {
        # celltype should be derived from the selected category in meta.data
        category_column <- sobj@misc$celltypes$category
        if (is.null(category_column) || !(category_column %in% colnames(sobj@meta.data))) {
          stop("'category' must be defined and match a column in meta.data before setting 'celltype'.")
        }
        sobj@misc$celltypes[[field]] <- unique(sobj@meta.data[[category_column]])
        reproducibility_log[[field]] <- paste0(
          "sobj@misc$celltypes[['", field, "']] <- unique(sobj@meta.data[['", category_column, "']])"
        )
        cat(paste("\nSet 'celltype' to unique values from:", category_column, "\n"))
      } else if (field == "description") {
        # description is manually specified
        cat("\nDescription of the annotation (e.g., 'Cell type annotations derived from meta.data').\n")
        user_input <- readline(prompt = "Enter a value for 'description': ")
        sobj@misc$celltypes[[field]] <- user_input
        reproducibility_log[[field]] <- paste0("sobj@misc$celltypes[['", field, "']] <- '", user_input, "'")
        cat(paste("\nSet 'description' to:", user_input, "\n"))
      } else if (field == "order") {
        # order is based on the unique values in celltype
        celltype_values <- sobj@misc$celltypes$celltype
        if (is.null(celltype_values)) {
          stop("'celltype' must be defined before setting 'order'.")
        }
        sobj@misc$celltypes[[field]] <- seq_along(celltype_values)
        reproducibility_log[[field]] <- paste0(
          "sobj@misc$celltypes[['", field, "']] <- seq_along(sobj@misc$celltypes[['celltype']])"
        )
        cat("\nSet 'order' to the sequence of unique cell types.\n")
      }
    } else {
      cat(paste("Field", field, "is already present in misc$celltypes.\n"))
    }
  }
  
  # Iterate through optional fields
  for (field in optional_fields) {
    if (is.null(sobj@misc$celltypes[[field]])) {
      # marker is manually specified
      cat("\nMarker genes associated with the cell type or cluster (e.g., 'CD3E, CD19, CD14').\n")
      user_input <- readline(prompt = paste("Enter a value for '", field, "' (optional, press Enter to skip): ", sep = ""))
      if (user_input != "") {
        sobj@misc$celltypes[[field]] <- user_input
        reproducibility_log[[field]] <- paste0("sobj@misc$celltypes[['", field, "']] <- '", user_input, "'")
        cat(paste("\nSet '", field, "' to:", user_input, "\n", sep = ""))
      }
    } else {
      cat(paste("Field", field, "is already present in misc$celltypes.\n"))
    }
  }
  
  # Print reproducibility log
  cat("\n----\nReproducibility Code:\n----\n")
  for (log in reproducibility_log) {
    cat(log, "\n")
  }
  
  cat("\nAll required fields are now present in misc$celltypes.\n")
  return(sobj)
}


##----------------
##----------------

ingestAbioData <- function(rds_in, rds_out, DESCRIPTION = "") {
    message("[INFO]: Reading Seurat object: ", rds_in)
    sobj <- readRDS(rds_in)

    # Validate Seurat object
    message("[INFO]: Seurat object version: ", paste(Version(sobj), collapse = "."))
    if (DefaultAssay(sobj) != "RNA") {
        message("[ERROR]: Default assay must be 'RNA'")
        return(NULL)
    }
    if (!all(c('data', 'counts') %in% Layers(sobj))) {
        message("[ERROR]: Missing 'counts' or 'data' layers")
        return(NULL)
    }

    required_meta <- c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'condition', 'fraction')
    if (!all(required_meta %in% colnames(sobj@meta.data))) {
        message("[ERROR]: Missing required metadata columns")
        return(NULL)
    }

    if (!any(c('umap', 'tsne') %in% names(sobj@reductions))) {
        message("[ERROR]: UMAP or t-SNE reduction required")
        return(NULL)
    }

    if (!all(c('experiment', 'celltypes') %in% names(Misc(sobj)))) {
        message("[ERROR]: Missing 'experiment' or 'celltypes' in misc")
        return(NULL)
    }

    # Fix column names
    sobj <- correct_colnames_in_misc_celltypes(sobj, 'celltypes', 'celltype')
    sobj <- correct_colnames_in_misc_celltypes(sobj, 'catagory', 'category')

    # Update experiment metadata
    experiment_name <- sub("\\.rds$", "", gsub("_", " ", rds_out))
    sobj@misc$experiment$experiment <- experiment_name
    sobj@misc$experiment$study <- experiment_name
    sobj@misc$experiment$order[1] <- 0
    sobj@misc$experiment$class <- 'oncology'
    sobj@misc$experiment$type <- 'singlecell'
    sobj@misc$experiment$species <- 'HS'
    PUBMEDID <- strsplit(rds_out, '_')[[1]][3]
    sobj@misc$experiment$reference <- PUBMEDID
    sobj@misc$experiment$notes <- ifelse(DESCRIPTION == "", "None", DESCRIPTION)

    if (!all(unique(sobj@misc$celltypes$category) %in% colnames(sobj@meta.data))) {
        message("[ERROR]: Categories in misc celltypes must match metadata column names")
        return(NULL)
    }

    # Validate Seurat object
    if (!isValid(sobj, unique(sobj@misc$celltypes$category))) {
        sobj[['RNA']]@counts <- sobj[['RNA']]@data
        sobj[['RNA']]@counts@x <- (exp(sobj[['RNA']]@data@x) - 1) / 10000
        sobj@meta.data$nCount_RNA <- 1

        if (!isValid(sobj, unique(sobj@misc$celltypes$category))) {
            message("[ERROR]: Validation failed (after fixing counts slot). Investigate manually.")
            return(NULL)
        }
    }

    # Write updated Seurat object
    message("[INFO]: Writing updated Seurat object: ", rds_out)
    saveRDS(sobj, file = rds_out)

    # Extract subset_name
    subset_name <- strsplit(sub("\\.rds$", "", rds_in), "_")[[1]][2]

    # Ingest and return job ID
    unique_categories <- paste(unique(sobj@misc$celltypes$category), collapse = ",")
    jobid <- ingest(rds_out, subset_name, unique_categories)

    message("[INFO]: Ingest successful. Job ID: ", jobid)
    return(jobid)
}

# Ingest Abio Data
ingestAbioData_old <- function(rds_in, rds_out, DESCRIPTION = "") {
    tryCatch({
        message("[INFO]: Reading Seurat object: ", rds_in)
        sobj <- readRDS(rds_in)


        message("[INFO]: Seurat object version: ", paste(Version(sobj), collapse = "."))
        if (DefaultAssay(sobj) != "RNA") stop("[ERROR]: Default assay must be 'RNA'")
        if (!all(c('data', 'counts') %in% Layers(sobj))) stop("[ERROR]: Missing 'counts' or 'data' layers")
        required_meta <- c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'condition', 'fraction')
        if (!all(required_meta %in% colnames(sobj@meta.data))) {
            warning("[ERROR]: Missing required metadata columns")
        }
        if (!any(c('umap', 'tsne') %in% names(sobj@reductions))) stop("[ERROR]: UMAP or t-SNE reduction required")
        if (!all(c('experiment', 'celltypes') %in% names(Misc(sobj)))) stop("[ERROR]: Missing 'experiment' or 'celltypes' in misc")

        sobj <- correct_colnames_in_misc_celltypes(sobj, 'celltypes', 'celltype')
        sobj <- correct_colnames_in_misc_celltypes(sobj, 'catagory', 'category')

        experiment_name <- sub("\\.rds$", "", gsub("_", " ", rds_out))
        sobj@misc$experiment$experiment <- experiment_name
        sobj@misc$experiment$study <- experiment_name
        sobj@misc$experiment$order[1] <- 0
        sobj@misc$experiment$class <- 'oncology'
        sobj@misc$experiment$type <- 'singlecell'
        sobj@misc$experiment$species <- 'HS'
        PUBMEDID <- strsplit(rds_out, '_')[[1]][3]
        sobj@misc$experiment$reference <- PUBMEDID
        sobj@misc$experiment$notes <- ifelse(DESCRIPTION == "", "None", DESCRIPTION)

        if (!all(unique(sobj@misc$celltypes$category) %in% colnames(sobj@meta.data))) {
            warning("[ERROR]: Categories in misc celltypes must match metadata column names")
        }

        if (!isValid(sobj, unique(sobj@misc$celltypes$category))) {
            sobj[['RNA']]@counts <- sobj[['RNA']]@data
            sobj[['RNA']]@counts@x <- (exp(sobj[['RNA']]@data@x) - 1) / 10000
            sobj@meta.data$nCount_RNA <- 1

            if (!isValid(sobj, unique(sobj@misc$celltypes$category))) {
                warning("[ERROR]: Validation failed (after fixing counts slot). Investigate manually.")
            }
        }

        message("[INFO]: Writing updated Seurat object: ", rds_out)
        saveRDS(sobj, file = rds_out)

        
        jobid <- tryCatch({
            unique_categories <- paste(unique(sobj@misc$celltypes$category), collapse = ",")
            
            subset_name <- strsplit(sub("\\.rds$", "", rds_in), "_")[[1]][5]

            ingest(rds_out, subset_name, unique_categories)
            
            message("[INFO]: Ingest successful. Job ID: ", jobid)
            return(jobid)
        }, error = function(e) {
            warning("[ERROR]: Ingestion failed for ", rds_out, " - ", e$message)
        })


    }, error = function(e) {
        message("[ERROR]: ", e$message)
        return(NULL)
    })
}

# Helper: Correct column names in misc celltypes
correct_colnames_in_misc_celltypes <- function(sobj, oldname, newname) {
    if (oldname %in% names(sobj@misc$celltypes)) {
        sobj@misc$celltypes[[newname]] <- sobj@misc$celltypes[[oldname]]
        sobj@misc$celltypes[[oldname]] <- NULL
    }
    return(sobj)
}

# Retrieve log for a job ID
jobidToLog <- function(jid) {
    log_output <- capture.output(printLog(jobId = jid))
    if (length(log_output) == 0) stop("[ERROR]: Log output is empty.")
    return(log_output)
}

# Retrieve experiment ID from job ID
jobidToExperimentid <- function(jobid) {
    log_output <- jobidToLog(jobid)
    experiment_id_line <- grep("experiment_id:", log_output, value = TRUE)
    if (length(experiment_id_line) == 0) stop("[ERROR]: 'experiment_id' not found in log output.")
    experiment_id <- sub(".*experiment_id:\\s*([^ ]+).*", "\\1", experiment_id_line)
    if (nchar(experiment_id) == 0) stop("[ERROR]: Failed to extract 'experiment_id'")
    return(experiment_id)
}

# Search experiment IDs by query
searchExperimentids <- function(query) {
    experiments <- getTestExperiments() %>%
        filter(if_any(everything(), ~ str_detect(as.character(.), regex(query, ignore_case = TRUE))))
    if (nrow(experiments) == 0) stop("No matching experiments found")
    list(experiment_ids = experiments$experiment_id, details = experiments)
}

# Retrieve experiment information
experimentidToInfo <- function(eid) {
    experiments <- getTestExperiments() %>% filter(experiment_id == eid)
    if (nrow(experiments) == 0) stop("[ERROR]: No experiment found for the given ID.")
    t(data.table(experiments))
}

# Cleanup all succeeded jobs
suceededJobsCleanupAll <- function(createdAt_prefix = NULL) {
    experiment_ids <- suceededJobsGetExperimentIds(createdAt_prefix)
    lapply(experiment_ids, function(expid) {
        message(sprintf("Executing cleanup for experiment ID: %s", expid))
        cleanup(experiment_id = expid)
    })
}

# Get experiment IDs for succeeded jobs
suceededJobsGetExperimentIds <- function(createdAt_prefix = NULL) {
    succeeded_jobs <- listJobs(hours=96) %>%
        filter(
            status == "SUCCEEDED",
            if (!is.null(createdAt_prefix)) startsWith(as.character(createdAt), createdAt_prefix) else TRUE
        )
    experiment_ids <- lapply(succeeded_jobs$jobId, jobidToExperimentid)
    if (length(experiment_ids) == 0) stop("[ERROR]: No experiment IDs could be retrieved.")
    return(unlist(experiment_ids))
}

