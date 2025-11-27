library("ggplot2")

library(anndata)
library(Matrix)
library(data.table)
library(sva)
library(Seurat)
library(pheatmap)
library(here)
library(jsonlite)
library(fs)
library(glue)

library(RColorBrewer) # Define a Brewer palette 
palette_patient <- brewer.pal(n = 8, name = "Accent")
palette_label <- c( "invasive cancer" = "red", "cancer in situ" = "orange", "immune infiltrate" = "yellow", "breast glands" = "green", "connective tissue" = "blue", "adipose tissue" = "cyan", "undetermined" = "lightgrey" )


rds_file <- here("results", "TNBC", "molecular", "filtering", "TNBC_filtered_seurat_object.rds")
source(here("digitalhistopathology", "molecular_helpers.R"))
setwd("./")
Sys.setenv(R_USER = "/./")


# Read the Seurat object from the RDS file 
seurat_object <- readRDS(rds_file) 

# Rename cell names: remove batch_...parquet. prefix and replace X with spot
old_names <- colnames(seurat_object)
new_names <- gsub("^batch_\\d+\\.parquet\\.", "", old_names)  # Remove batch_XX.parquet.
new_names <- gsub("_X", "_spot", new_names)  # Replace _X with _spot
new_names <- gsub("_", "-", new_names)  # Replace remaining underscores with hyphens

# Update cell names in all assays
for (assay_name in names(seurat_object@assays)) {
  colnames(seurat_object@assays[[assay_name]]) <- new_names
}

# Update cell names in metadata
rownames(seurat_object@meta.data) <- new_names

# Update cell names in reductions if they exist
if (length(seurat_object@reductions) > 0) {
  for (reduction_name in names(seurat_object@reductions)) {
    rownames(seurat_object@reductions[[reduction_name]]@cell.embeddings) <- new_names
  }
}

cat("Cell names updated. Example:", head(colnames(seurat_object), 3), "\n")


n <- 5
# Directory to search for label files
file <- here("results", "transferability_exp", "TNBC_predicted_labels_from_HER2_space.csv")


# Load true labels DGE results
res_true_labels <- load_dge_pathways_analysis_per_clusters(directory_name = here("results", "TNBC", "molecular", "tissue_DGE"), add_name = "_upregulated_true_labels_all")

# Get pathways heatmaps for upregulated pathways
results_upregulated <- get_pathways_heatmaps(labels_clusters_uni_file = file, seurat_object = seurat_object, res_true_labels = res_true_labels, upregulated = TRUE, add_name = paste0("_", n, "_clusters_all"))
