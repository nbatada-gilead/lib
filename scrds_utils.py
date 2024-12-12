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

# Ingest Abio Data
ingestAbioData <- function(rds_in, rds_out, DESCRIPTION = "") {
    tryCatch({
        message("[INFO]: Reading Seurat object: ", rds_in)
        sobj <- readRDS(rds_in)

        message("[INFO]: Seurat object version: ", paste(Version(sobj), collapse = "."))
        if (DefaultAssay(sobj) != "RNA") stop("[ERROR]: Default assay must be 'RNA'")
        if (!all(c('data', 'counts') %in% Layers(sobj))) stop("[ERROR]: Missing 'counts' or 'data' layers")
        required_meta <- c('nCount_RNA', 'nFeature_RNA', 'percent.mito', 'condition', 'fraction')
        if (!all(required_meta %in% colnames(sobj@meta.data))) {
            stop("[ERROR]: Missing required metadata columns")
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
                stop("[ERROR]: Validation failed. Investigate manually.")
            }
        }

        message("[INFO]: Writing updated Seurat object: ", rds_out)
        saveRDS(sobj, file = rds_out)

        jobid <- tryCatch({
            unique_categories <- paste(unique(sobj@misc$celltypes$category), collapse = ",")
            ingest(rds_out, subset_name, unique_categories)
        }, error = function(e) {
            stop("[ERROR]: Ingestion failed for ", rds_out, " - ", e$message)
        })

        message("[INFO]: Ingest successful. Job ID: ", jobid)
        return(jobid)
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
    succeeded_jobs <- listJobs() %>%
        filter(
            status == "SUCCEEDED",
            if (!is.null(createdAt_prefix)) startsWith(as.character(createdAt), createdAt_prefix) else TRUE
        )
    experiment_ids <- lapply(succeeded_jobs$jobId, jobidToExperimentid)
    if (length(experiment_ids) == 0) stop("[ERROR]: No experiment IDs could be retrieved.")
    return(unlist(experiment_ids))
}

