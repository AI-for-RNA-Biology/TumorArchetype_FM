
# if (!require("BiocManager", quietly = TRUE))
#   BiocManager::install("BiocManager")
# if (!require("DESeq2", quietly = TRUE))
#   BiocManager::install("DESeq2")
# library("DESeq2")
# if (!require("ape", quietly = TRUE))
#   BiocManager::install("ape")
# library("ape")
# if (!require("limma", quietly = TRUE))
#   BiocManager::install("limma")
# library("limma")
# if (!require("geneplotter", quietly = TRUE))
#   BiocManager::install("geneplotter")
# library("geneplotter")
# if (!require("gprofiler2", quietly = TRUE))
#   BiocManager::install("gprofiler2")
# library("gprofiler2")
# if (!require("edgeR", quietly = TRUE))
#   BiocManager::install("edgeR")
# library("edgeR")
# if (!require("AMR", quietly = TRUE))
#   BiocManager::install("AMR")
# library("AMR")
# if (!require("Seurat", quietly = TRUE))
#   BiocManager::install("Seurat")
# library("Seurat")
# if (!require("dplyr", quietly = TRUE))
#   install.packages("dplyr")
# library("dplyr")
# if (!require("patchwork", quietly = TRUE))
#   BiocManager::install("patchwork")
# library("patchwork")
# if (!require("Matrix", quietly = TRUE))
#   BiocManager::install("h")
# library("Matrix")
# if (!require("ggplot2", quietly = TRUE))
#   BiocManager::install("ggplot2")
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


# Get model name from command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: Model name must be provided as an argument.")
}

model <- args[1]
cat("Processing model:", model, "\n")

# Load config file
config <- fromJSON(here("config", "config_notebooks_TNBC.json"))


if (model %in% config$retrained_model_list1) {
  # Remove all ../ from the beginning of the path
  relative_path <- gsub("^(\\.\\./)+", "", config$retrained_benchmark_folder1)
  folder <- here(relative_path)
} else {
  # Remove all ../ from the beginning of the path
  relative_path <- gsub("^(\\.\\./)+", "", config$retrained_benchmark_folder2)
  folder <- here(relative_path)
}


# Build path to opti_clusters.csv
path_opti_clusters <- fs::path(folder, "invasive_cancer_clustering", 
                                config$invasive_cancer_clustering_algo, "optimal_clusters.csv")

# Read the CSV
opti_clusters <- read.csv(path_opti_clusters)

# Extract number of clusters for the current model
# Assuming the first column contains the model names
n <- opti_clusters$n_clusters[opti_clusters[[1]] == model]

# Directory to search for label files
label_dir <- fs::path(folder, "invasive_cancer_clustering", 
                      config$invasive_cancer_clustering_algo, model)

# Build the pattern with the current value of `n`
pattern <- glue("invasive_labels_{n}_clusters_umap.*\\.csv$")

# Find matching files
path_labels <- dir(label_dir, pattern = pattern, full.names = TRUE)

if (length(path_labels) > 1) {
  warning(glue("More than one file found for model {model}"))
}

# Take the first matching file
file <- path_labels[1]

# Load true labels DGE results
res_true_labels <- load_dge_pathways_analysis_per_clusters(directory_name = here("results", "TNBC", "molecular", "tissue_DGE"), add_name = "_upregulated_true_labels")

# Get number of clusters from filename
n <- gsub(".*labels_(\\d+)_clusters.*", "\\1", file)

# Get pathways heatmaps for upregulated pathways
results_upregulated <- get_pathways_heatmaps(labels_clusters_uni_file = file, seurat_object = seurat_object, res_true_labels = res_true_labels, upregulated = TRUE, add_name = paste0("_", n, "_clusters"))
