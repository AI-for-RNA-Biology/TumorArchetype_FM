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
if (file.exists(rds_file)) { 
  # Read the Seurat object from the RDS file 
  seurat_object <- readRDS(rds_file) 
  cat("Seurat object loaded from", rds_file, "\n") 
    # Apply the function to each patient and store the results in a list
  
} else { 
  # Create a new Seurat object 
  # Read the .csv file
    csv_file <- here("results", "Ovarian", "molecular", "df_spots_filtered.csv")
    m <- read.csv(csv_file, row.names="name")
    # Read the metadata
    spots_metadata_file <- here("results", "Ovarian", "molecular", "spots_metadata.csv") 
    spots_metadata <- read.csv(spots_metadata_file, row.names="name.1")
    
    rownames(m) <- gsub("_", "-", rownames(m))
    rownames(spots_metadata) <- gsub("_", "-", rownames(spots_metadata))
    
    # Create the Seurat object
    seurat_object <- CreateSeuratObject(counts = t(m), meta.data = spots_metadata)

    # Primary cleaning
    # Sum the counts for each cell 
    total_counts_per_cell <- colSums(seurat_object@assays$RNA$counts) 
    # Identify cells with zero counts 
    cells_to_keep <- total_counts_per_cell > 0 
    # Filter out cells with zero counts 
    seurat_object <- seurat_object[, cells_to_keep]
    
    seurat_object <- NormalizeData(seurat_object, assay = "RNA")
    seurat_object <- FindVariableFeatures(seurat_object, assay = "RNA", selection.method = "vst", nfeatures = 15000)
    seurat_object <- ScaleData(seurat_object, assay = "RNA")
    
    seurat_object <- RunPCA(seurat_object, assay = "RNA")
    seurat_object <- RunUMAP(seurat_object, assay = "RNA", dims = 1:30)
    
    
    # Save the Seurat object to a file 
    saveRDS(seurat_object, file = here("results", "Ovarian", "molecular", "gene_embeddings.rds"))
  
    cat("Seurat object created and saved to", here("results", "Ovarian", "molecular", "gene_embeddings.rds"), "\n") 

    # Assuming seurat_object is your Seurat object
    rna_data <- seurat_object@assays$RNA$data

    # Write the RNA data to a CSV file
    write.csv(as.matrix(rna_data), file = here("results", "Ovarian", "molecular", "filtered_normalized_gene_expression.csv"))
    cat("Filtered normalized gene expression data saved to", here("results", "Ovarian", "molecular", "filtered_normalized_gene_expression.csv"), "\n")
}


layout(matrix(ncol=2,nrow=1,c(1:2),byrow = TRUE))
# Add number of genes per UMI for each cell to metadata
seurat_object$log10GenesPerUMI <- log10(seurat_object$nFeature_RNA) / log10(seurat_object$nCount_RNA)
#another notation: GE[["log10GenesPerUMI"]]        <- log10(GE$nFeature_RNA) / log10(GE$nCount_RNA)
# Compute percent mito ratio
seurat_object$mitoRatio        <- PercentageFeatureSet(object = seurat_object, pattern = "^MT-")/100
#another notation: GE[["mitoRatio"]]        <- PercentageFeatureSet(object = GE, pattern = "^MT-")/100

#Compare #genes versus UMI
smoothScatter(seurat_object$nFeature_RNA,seurat_object$nCount_RNA,las=1,main="",xlab="# genes",ylab="# UMI")

seurat_object_labeled <- seurat_object[ , rownames(seurat_object@meta.data)[seurat_object@meta.data$label != ""]]
seurat_object_labeled <- seurat_object_labeled[ , rownames(seurat_object_labeled@meta.data)[seurat_object_labeled@meta.data$label != "undetermined"]]
# Set the identities in the Seurat object to your custom labels
Idents(seurat_object_labeled) <- seurat_object_labeled@meta.data$label


start_time <- Sys.time()

res_dge <- get_clusters_DGE_BPs(seurat_object_labeled, upregulated = upregulated)
save_dge_pathways_analysis_per_clusters(res = res_dge, directory_name = output_dir, add_name = paste0("_", direction, "_true_labels_all"))

end_time <- Sys.time()
cat("DGE analysis completed in", round(as.numeric(end_time - start_time, units = "hours"), 2), "hours\n")

# Save the DGE results
cat("Saving DGE analysis results...\n")

for (cluster in names(res_dge$gprofiler_results)) {
    ggsave(filename = file.path(figures_dir, paste0("pathway_cluster_", direction, "_", cluster, "_all.png")), 
           plot = plot_pathways(cluster, res_dge$gprofiler_results[[cluster]]), 
           width = 10, height = 8)
}
cat("Saved pathway plots for", length(names(res_dge$gprofiler_results)), "clusters\n")

result_matrix <- get_pathway_scores_across_all_clusters(res = res_dge)
heatmap_pathways(result_matrix = result_matrix, display_numbers = TRUE, directory_name = figures_dir, name = paste0("pathway_heatmap_", direction, "_known_tissues_all"))
