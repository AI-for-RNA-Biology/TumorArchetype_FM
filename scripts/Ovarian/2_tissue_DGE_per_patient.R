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



source(here("digitalhistopathology", "molecular_helpers.R"))

upregulated <- TRUE

palette_patient <- brewer.pal(n = 8, name = "Accent")
palette_label <- c("invasive cancer" = "red", "immune infiltrate" = "yellow", "epithelial" = "green", "connective tissue" = "blue", "undetermined" = "lightgrey" )

output_dir <- here("results", "Ovarian", "molecular")

# Create output directory for figures
output_dir <- file.path(output_dir, "tissue_DGE")
figures_dir <- file.path(output_dir, "Figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

direction <- ifelse(upregulated, "upregulated", "downregulated")

rds_file <- here("results", "Ovarian", "molecular", "gene_embeddings.rds")


cat("=== TNBC Tissue DGE Analysis ===\n")
cat("Analysis type:", direction, "genes\n")
cat("RDS file:", rds_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Figures will be saved to:", figures_dir, "\n")

# Check if the RDS file exists 

  # Read the Seurat object from the RDS file 
seurat_object <- readRDS(rds_file) 
cat("Seurat object loaded from", rds_file, "\n") 

# Get unique patients
unique_patients <- unique(seurat_object$patient)
cat("Found", length(unique_patients), "unique patients:", paste(unique_patients, collapse = ", "), "\n")

# Iterate over each patient
for (patient_id in unique_patients) {
  cat("\n========================================\n")
  cat("Processing patient:", patient_id, "\n")
  cat("========================================\n")
  
  # Create patient-specific output directories
  patient_output_dir <- file.path(output_dir, "per_patient", paste0("patient_", patient_id))
  patient_figures_dir <- file.path(patient_output_dir, "Figures")
  dir.create(patient_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(patient_figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Subset Seurat object for current patient
  seurat_patient <- seurat_object[, seurat_object$patient == patient_id]
  cat("Patient", patient_id, "has", ncol(seurat_patient), "cells/spots\n")
  
  # QC plots for this patient
  pdf(file.path(patient_figures_dir, paste0("QC_patient_", patient_id, ".pdf")), width = 12, height = 6)
  layout(matrix(ncol=2, nrow=1, c(1:2), byrow = TRUE))
  
  # Add number of genes per UMI for each cell to metadata
  seurat_patient$log10GenesPerUMI <- log10(seurat_patient$nFeature_RNA) / log10(seurat_patient$nCount_RNA)
  
  # Compute percent mito ratio
  seurat_patient$mitoRatio <- PercentageFeatureSet(object = seurat_patient, pattern = "^MT-") / 100
  
  # Compare #genes versus UMI
  smoothScatter(seurat_patient$nFeature_RNA, seurat_patient$nCount_RNA, 
                las=1, main=paste("Patient", patient_id), xlab="# genes", ylab="# UMI")
  dev.off()
  
  # Filter for labeled samples (exclude empty and undetermined)
  seurat_patient_labeled <- seurat_patient[, rownames(seurat_patient@meta.data)[seurat_patient@meta.data$label != ""]]
  seurat_patient_labeled <- seurat_patient_labeled[, rownames(seurat_patient_labeled@meta.data)[seurat_patient_labeled@meta.data$label != "undetermined"]]
  
  cat("After filtering, patient", patient_id, "has", ncol(seurat_patient_labeled), "labeled cells/spots\n")
  
  # Check if there are enough labeled samples
  if (ncol(seurat_patient_labeled) < 10) {
    cat("WARNING: Patient", patient_id, "has too few labeled samples (<10). Skipping DGE analysis.\n")
    next
  }
  
  # Check if there are multiple tissue labels
  label_counts <- table(seurat_patient_labeled$label)
  cat("Label distribution for patient", patient_id, ":\n")
  print(label_counts)
  
  if (length(label_counts) < 2) {
    cat("WARNING: Patient", patient_id, "has only one tissue label. Skipping DGE analysis.\n")
    next
  }
  
  # Set the identities in the Seurat object to your custom labels
  Idents(seurat_patient_labeled) <- seurat_patient_labeled@meta.data$label
  
  # Run DGE analysis
  start_time <- Sys.time()


  res_dge <- get_clusters_DGE_BPs(seurat_patient_labeled, upregulated = upregulated)
  save_dge_pathways_analysis_per_clusters(res = res_dge, 
                                          directory_name = patient_output_dir, 
                                          add_name = paste0("_", direction, "_patient_", patient_id))
  
  end_time <- Sys.time()
  cat("DGE analysis for patient", patient_id, "completed in", 
      round(as.numeric(end_time - start_time, units = "mins"), 2), "minutes\n")
  
  # Save pathway plots for this patient
  cat("Saving pathway plots for patient", patient_id, "...\n")
  for (cluster in names(res_dge$gprofiler_results)) {
    ggsave(filename = file.path(patient_figures_dir, 
                                paste0("pathway_cluster_", direction, "_", cluster, "_patient_", patient_id, ".png")), 
            plot = plot_pathways(cluster, res_dge$gprofiler_results[[cluster]]), 
            width = 10, height = 8)
  }
  cat("Saved pathway plots for", length(names(res_dge$gprofiler_results)), "clusters\n")
  
  # Create pathway heatmap for this patient
  result_matrix <- get_pathway_scores_across_all_clusters(res = res_dge)
  heatmap_pathways(result_matrix = result_matrix, 
                  display_numbers = TRUE, 
                  directory_name = patient_figures_dir, 
                  name = paste0("pathway_heatmap_", direction, "_patient_", patient_id))
  
  cat("Patient", patient_id, "analysis completed successfully!\n")
    
}

cat("\n========================================\n")
cat("All patients processed!\n")
cat("========================================\n")
