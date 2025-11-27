library(Seurat)
library(sva)
library(arrow)
library(Matrix)
library(here)


rds_file <- here("results", "Ovarian", "molecular", "gene_embeddings.rds")

# Check if the RDS file exists 

  # Read the Seurat object from the RDS file 
seurat_object <- readRDS(rds_file) 

# Save filtered gene expression
filtered_gene_expression <- as.matrix(seurat_object@assays$RNA$counts)

# Save normalized data
normalized_data <- as.matrix(seurat_object@assays$RNA$data)

# Apply ComBat batch correction
scaled_data <- as.matrix(seurat_object@assays$RNA$scale.data)
batch <- seurat_object@meta.data$patient
combat_corrected_data <- ComBat(dat = scaled_data, batch = batch)

combat_corrected_file <-  here("results", "Ovarian", "molecular", "combat", "combat_corrected_filtered_counts.csv")
write.csv(combat_corrected_data, file = combat_corrected_file)
