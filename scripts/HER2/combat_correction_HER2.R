# Title: combat_correction
# Date: 2025-01-17

# Load required libraries
if (!require("BiocManager", quietly = TRUE))
  BiocManager::install("BiocManager")
if (!require("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")
library("DESeq2")
if (!require("ape", quietly = TRUE))
  BiocManager::install("ape")
library("ape")
if (!require("limma", quietly = TRUE))
  BiocManager::install("limma")
library("limma")
if (!require("geneplotter", quietly = TRUE))
  BiocManager::install("geneplotter")
library("geneplotter")
if (!require("gprofiler2", quietly = TRUE))
  BiocManager::install("gprofiler2")
library("gprofiler2")
if (!require("edgeR", quietly = TRUE))
  BiocManager::install("edgeR")
library("edgeR")
if (!require("AMR", quietly = TRUE))
  BiocManager::install("AMR")
library("AMR")
if (!require("Seurat", quietly = TRUE))
  BiocManager::install("Seurat")
library("Seurat")
if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
library("dplyr")
if (!require("patchwork", quietly = TRUE))
  BiocManager::install("patchwork")
library("patchwork")
if (!require("Matrix", quietly = TRUE))
  BiocManager::install("h")
library("Matrix")
if (!require("ggplot2", quietly = TRUE))
  BiocManager::install("ggplot2")
library("ggplot2")

library(anndata)
library(Matrix)
library(data.table)
library(sva)
library(harmony)
library(crescendo)
library(pheatmap)
library(RColorBrewer)

# Define palettes
palette_patient <- brewer.pal(n = 8, name = "Accent")
palette_label <- c(
  "invasive cancer" = "red", "cancer in situ" = "orange", 
  "immune infiltrate" = "yellow", "breast glands" = "green", 
  "connective tissue" = "blue", "adipose tissue" = "cyan", 
  "undetermined" = "lightgrey"
)

# Load helper functions
source("../digitalhistopathology/molecular_helpers.R")

# Set working directory and environment
setwd("./")
Sys.setenv(R_USER = "/./")

# Function to create a pseudobulk profile by summing counts
create_pseudobulk <- function(obj) {
  pseudobulk <- Matrix::rowSums(obj@assays$RNA$counts)
  return(pseudobulk)
}

# Load or create Seurat object
rds_file <- "../results/molecular/gene_embeddings_HER2.rds"
if (file.exists(rds_file)) {
  # Load Seurat object
  seurat_object <- readRDS(rds_file)
  cat("Seurat object loaded from", rds_file, "\n")
  
  # Split by patient and create pseudobulk profiles
  seurat_list <- SplitObject(seurat_object, split.by = "patient")
  pseudobulk_list <- lapply(seurat_list, create_pseudobulk)
  pseudobulk_list_log <- lapply(pseudobulk_list, function(x) log2(x + 1))
  
  # Plot pseudobulk densities
  par(mfrow = c(2, 4))
  for (i in seq_along(pseudobulk_list_log)) {
    pseudobulk_log <- pseudobulk_list_log[[i]]
    pseudobulk_log <- pseudobulk_log[pseudobulk_log > 0]
    bimdens <- mclust::densityMclust(data = pseudobulk_log, G = 2, plot = FALSE)
    lim <- qnorm(0.99, mean = bimdens$parameters$mean[1], sd = sqrt(bimdens$parameters$variance$sigmasq[1]))
    if (names(pseudobulk_list_log)[i] == 'H') { lim <- 2 }
    plot(density(pseudobulk_log), main = names(pseudobulk_list_log)[i], xlab = "log2(# genes per spots + 1)", ylab = "Density", col = palette_patient[i])
    abline(v = lim, col = "black", lwd = 2.5)
  }
} else {
  # Create Seurat object from CSV files
  csv_file <- "../results/molecular/gene_embedding_HER2.csv"
  spots_metadata_file <- "../results/molecular/spots_metadata.csv"
  m <- read.csv(csv_file, row.names = "index")
  spots_metadata <- read.csv(spots_metadata_file, row.names = "index")
  rownames(m) <- gsub("_", "-", rownames(m))
  rownames(spots_metadata) <- gsub("_", "-", rownames(spots_metadata))
  seurat_object <- CreateSeuratObject(counts = t(m), meta.data = spots_metadata)
  seurat_object$patient <- substr(seurat_object$name_origin, 1, 1)
  
  # Filter and clean data
  total_counts_per_cell <- colSums(seurat_object@assays$RNA$counts)
  cells_to_keep <- total_counts_per_cell > 0
  seurat_object <- seurat_object[, cells_to_keep]
  seurat_object <- seurat_object[, rownames(seurat_object@meta.data)[seurat_object@meta.data$patient != "A"]]
  
  # Split by patient and create pseudobulk profiles
  seurat_list <- SplitObject(seurat_object, split.by = "patient")
  pseudobulk_list <- lapply(seurat_list, create_pseudobulk)
  pseudobulk_list_log <- lapply(pseudobulk_list, function(x) log2(x + 1))
  
  # Plot pseudobulk densities and filter genes
  par(mfrow = c(2, 4))
  limits <- list()
  genes_above_limit <- list()
  for (i in seq_along(pseudobulk_list_log)) {
    pseudobulk_log <- pseudobulk_list_log[[i]]
    pseudobulk_log <- pseudobulk_log[pseudobulk_log > 0]
    bimdens <- mclust::densityMclust(data = pseudobulk_log, G = 2, plot = FALSE)
    lim <- qnorm(0.99, mean = bimdens$parameters$mean[1], sd = sqrt(bimdens$parameters$variance$sigmasq[1]))
    if (names(pseudobulk_list_log)[i] == 'H') { lim <- 2 }
    plot(density(pseudobulk_log), main = names(pseudobulk_list_log)[i], xlab = "log2(# genes per spots + 1)", ylab = "Density", col = palette_patient[i])
    abline(v = lim, col = "black", lwd = 2.5)
    genes_above_limit[[names(pseudobulk_list_log)[i]]] <- names(pseudobulk_log)[pseudobulk_log > lim]
  }
  common_genes <- Reduce(intersect, genes_above_limit)
  seurat_object <- subset(seurat_object, features = common_genes)
  
  # Normalize and scale data
  seurat_object <- NormalizeData(seurat_object, assay = "RNA")
  seurat_object <- FindVariableFeatures(seurat_object, assay = "RNA", selection.method = "vst", nfeatures = 19000)
  seurat_object <- ScaleData(seurat_object, assay = "RNA")
  seurat_object <- RunPCA(seurat_object, assay = "RNA")
  seurat_object <- RunUMAP(seurat_object, assay = "RNA", dims = 1:30)
  
  # Save Seurat object
  saveRDS(seurat_object, file = "../results/molecular/gene_embeddings_HER2.rds")
  cat("Seurat object created and saved to", rds_file, "\n")
}

# Save filtered gene expression
filtered_gene_expression <- as.matrix(seurat_object@assays$RNA$counts)
write.csv(filtered_gene_expression, file = "../results/molecular/filtered_gene_expression.csv", row.names = TRUE)

# Save normalized data
normalized_data <- as.matrix(seurat_object@assays$RNA$data)
write.csv(normalized_data, file = "filtered_normalized_gene_expression.csv")

# Apply ComBat batch correction
scaled_data <- as.matrix(seurat_object@assays$RNA$scale.data)
batch <- seurat_object@meta.data$patient
combat_corrected_data <- ComBat(dat = scaled_data, batch = batch)
write.csv(combat_corrected_data, file = "../results/molecular/combat_corrected_filtered_counts.csv")

# Create Seurat object with batch-corrected data
seurat_object_combat_corrected <- CreateSeuratObject(counts = seurat_object@assays$RNA$counts)
seurat_object_combat_corrected@assays$RNA$scale.data <- combat_corrected_data
seurat_object_combat_corrected@meta.data <- seurat_object@meta.data
seurat_object_combat_corrected <- RunPCA(seurat_object_combat_corrected, features = VariableFeatures(seurat_object))
seurat_object_combat_corrected <- RunUMAP(seurat_object_combat_corrected, dims = 1:10)

# Plot UMAPs
DimPlot(seurat_object, reduction = "umap", group.by = "patient", cols = palette_patient)
DimPlot(seurat_object_combat_corrected, reduction = "umap", group.by = "patient", cols = palette_patient)
seurat_object_combat_corrected_labeled <- seurat_object_combat_corrected[, rownames(seurat_object_combat_corrected@meta.data)[seurat_object_combat_corrected@meta.data$label != ""]]
DimPlot(seurat_object_combat_corrected_labeled, reduction = "umap", group.by = "label", cols = palette_label)
