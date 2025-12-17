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


rds_file <- here("results", "HER2", "molecular",  "gene_embeddings_HER2.rds")
output_dir <- here("results", "HER2", "molecular", "tissue_DGE")
figures_dir <- file.path(output_dir, "Figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

upregulated <- TRUE

direction <- ifelse(upregulated, "upregulated", "downregulated")

source(here("digitalhistopathology", "molecular_helpers.R"))
setwd("./")
Sys.setenv(R_USER = "/./")


# Read the Seurat object from the RDS file 
seurat_object <- readRDS(rds_file) 

seurat_object_labeled <- seurat_object[ , rownames(seurat_object@meta.data)[seurat_object@meta.data$label != ""]]
seurat_object_labeled <- seurat_object_labeled[ , rownames(seurat_object_labeled@meta.data)[seurat_object_labeled@meta.data$label != "undetermined"]]
Idents(seurat_object_labeled) <- seurat_object_labeled@meta.data$label
cat("Seurat object with labeled cells loaded from", rds_file, "\n")

# Get model name from command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: Model name must be provided as an argument.")
}

model <- args[1]
cat("Processing model:", model, "\n")

# Load config file
config <- fromJSON(here("config", "config_notebooks_HER2.json"))


if (model %in% config$retrained_model_list1) {
  # Remove all ../ from the beginning of the path
  relative_path <- gsub("^(\\.\\./)+", "", config$retrained_benchmark_folder1)
  folder <- here(relative_path)
} else {
  # Remove all ../ from the beginning of the path
  relative_path <- gsub("^(\\.\\./)+", "", config$retrained_benchmark_folder2)
  folder <- here(relative_path)
}

n <- 7
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
# res_true_labels <- load_dge_pathways_analysis_per_clusters(directory_name = here("results", "HER2", "molecular", "tissue_DGE"), add_name = "_upregulated_true_labels_all")
res_true_labels <- get_clusters_DGE_BPs(seurat_object_labeled)
save_dge_pathways_analysis_per_clusters(res = res_true_labels, directory_name = output_dir, add_name = paste0("_", direction, "_true_labels_all"))
result_matrix <- get_pathway_scores_across_all_clusters(res = res_true_labels)
heatmap_pathways(result_matrix = result_matrix, display_numbers = TRUE, directory_name = figures_dir, name = paste0("pathway_heatmap_", direction, "_known_tissues_all"))