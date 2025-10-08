#!/usr/bin/env Rscript

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
palette_label <- c( "invasive cancer" = "red", "cancer in situ" = "orange", "immune infiltrate" = "yellow", "breast glands" = "green", "connective tissue" = "blue", "adipose tissue" = "cyan", "undetermined" = "lightgrey" )

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  # Default values
  rds_file <- "/storage/research/dbmr_luisierlab/temp/lfournier/repositories/TumorArchetype-FM/results/molecular/TNBC_processed/TNBC_filtered_seurat_object.rds"
  output_dir <- here("Figures", "tissue_DGE")
  upregulated <- TRUE
} else if (length(args) == 1) {
  rds_file <- args[1]
  output_dir <- here("Figures", "tissue_DGE")
  upregulated <- TRUE
} else if (length(args) == 2) {
  rds_file <- args[1]
  output_dir <- args[2]
  upregulated <- TRUE
} else if (length(args) == 3) {
  rds_file <- args[1]
  output_dir <- args[2]
  upregulated <- as.logical(args[3])
} else {
  cat("Usage: Rscript 4_tissue_DGE.R [rds_file] [output_dir] [upregulated]\n")
  cat("  rds_file: Path to the Seurat object RDS file\n")
  cat("  output_dir: Directory to save figures (will create subdirectories)\n")
  cat("  upregulated: TRUE for upregulated genes, FALSE for downregulated genes\n")
  quit(status = 1)
}

# Create output directory for figures
output_dir <- file.path(output_dir, "tissue_DGE")
figures_dir <- file.path(output_dir, "Figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

direction <- ifelse(upregulated, "upregulated", "downregulated")

cat("=== TNBC Tissue DGE Analysis ===\n")
cat("Analysis type:", direction, "genes\n")
cat("RDS file:", rds_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Figures will be saved to:", figures_dir, "\n")

# Check if files exist
if (!file.exists(rds_file)) {
  stop("RDS file does not exist: ", rds_file)
}

source("/storage/research/dbmr_luisierlab/temp/lfournier/repositories/TumorArchetype-FM/digitalhistopathology/molecular_helpers.R")

create_pseudobulk <- function(obj) {
  pseudobulk <- Matrix::rowSums(obj@assays$RNA$counts)
  return(pseudobulk)
}

seurat_object <- readRDS(rds_file) 
cat("Seurat object loaded from", rds_file, "\n") 
# Extract patient ID from spot names (format: filename.parquet.TNBC1_X...)
# Split by "." and take the last part before "_", then split by "_" and take first part
seurat_object$patient <- sapply(colnames(seurat_object), function(x) {
  # Remove everything up to the last "parquet." 
  after_parquet <- sub(".*parquet\\.", "", x)
  # Split by "_" and take the first part (TNBC1)
  strsplit(after_parquet, "_")[[1]][1]
})

# Debug: Check patient extraction
cat("Sample spot names:", paste(head(colnames(seurat_object), 3), collapse = ", "), "\n")
cat("Extracted patients:", paste(head(unique(seurat_object$patient), 5), collapse = ", "), "\n")
# Apply the function to each patient and store the results in a list
# Alternative to SplitObject for compatibility
unique_patients <- unique(seurat_object$patient)
cat("Found", length(unique_patients), "unique patients:", paste(unique_patients, collapse = ", "), "\n")

seurat_list <- list()
for (patient in unique_patients) {
  patient_cells <- colnames(seurat_object)[seurat_object$patient == patient]
  seurat_list[[patient]] <- subset(seurat_object, cells = patient_cells)
  cat("Patient", patient, ":", length(patient_cells), "spots\n")
}

pseudobulk_list <- lapply(seurat_list, create_pseudobulk)
pseudobulk_list_log <- lapply(pseudobulk_list, function(x) log2(x + 1))

# Plot number of genes per UMI for each cell to metadata
seurat_object$log10GenesPerUMI <- log10(seurat_object$nFeature_RNA) / log10(seurat_object$nCount_RNA)

# Calculate mitochondrial gene percentage manually
mito_genes <- grep("^MT-", rownames(seurat_object), value = TRUE)
if (length(mito_genes) > 0) {
  mito_counts <- Matrix::colSums(seurat_object@assays$RNA$counts[mito_genes, , drop = FALSE])
  seurat_object$mitoRatio <- mito_counts / seurat_object$nCount_RNA
  cat("Found", length(mito_genes), "mitochondrial genes\n")
} else {
  seurat_object$mitoRatio <- 0
  cat("No mitochondrial genes found with pattern '^MT-'\n")
}

# Save smoothScatter plot
png(file.path(figures_dir, "genes_vs_UMI_scatter.png"), width = 800, height = 600)
smoothScatter(seurat_object$nFeature_RNA, seurat_object$nCount_RNA, las=1, main="", xlab="# genes", ylab="# UMI")
dev.off()
cat("Saved genes vs UMI scatter plot\n")


seurat_object_labeled <- seurat_object[ , rownames(seurat_object@meta.data)[seurat_object@meta.data$label != ""]]
seurat_object_labeled <- seurat_object_labeled[ , rownames(seurat_object_labeled@meta.data)[seurat_object_labeled@meta.data$label != "undetermined"]]
# Set the identities in the Seurat object to your custom labels
Idents(seurat_object_labeled) <- seurat_object_labeled@meta.data$label

cat("Starting DGE analysis with MAST on", ncol(seurat_object_labeled), "cells...\n")
cat("Analysis type:", direction, "genes\n")
cat("This may take several hours with MAST and large dataset.\n")

# Show tissue distribution
tissue_counts <- table(seurat_object_labeled$label)
cat("Tissue distribution:\n")
print(tissue_counts)

# Set options for better MAST performance
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB limit for future objects

start_time <- Sys.time()

res_dge <- get_clusters_DGE_BPs(seurat_object_labeled, upregulated = upregulated)
save_dge_pathways_analysis_per_clusters(res = res_dge, directory_name = output_dir, add_name = paste0("_", direction, "_true_labels"))

end_time <- Sys.time()
cat("DGE analysis completed in", round(as.numeric(end_time - start_time, units = "hours"), 2), "hours\n")

# Save the DGE results
cat("Saving DGE analysis results...\n")

for (cluster in names(res_dge$gprofiler_results)) {
    ggsave(filename = file.path(figures_dir, paste0("pathway_cluster_", direction, "_", cluster, ".png")), 
           plot = plot_pathways(cluster, res_dge$gprofiler_results[[cluster]]), 
           width = 10, height = 8)
}
cat("Saved pathway plots for", length(names(res_dge$gprofiler_results)), "clusters\n")

result_matrix <- get_pathway_scores_across_all_clusters(res = res_dge)
heatmap_pathways(result_matrix = result_matrix, display_numbers = TRUE, name = paste0("_", direction), directory_name = figures_dir, name = paste0("pathway_heatmap_", direction, "_known_tissues"))