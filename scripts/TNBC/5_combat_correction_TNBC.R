


#!/usr/bin/env Rscript

library(Seurat)
library(sva)
library(arrow)
library(Matrix)
library(here)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript 5_combat_correction_TNBC.R <input_rds> <output_parquet> [batch_column] [use_scale_data]\n")
  cat("  input_rds: Path to the Seurat object RDS file\n")
  cat("  output_parquet: Path to output parquet file\n")
  cat("  batch_column: Column name for batch correction (default: 'patient')\n")
  cat("  use_scale_data: Use scaled data (TRUE) or raw counts (FALSE, default: FALSE)\n")
  quit(status = 1)
}

# Parse arguments
input_rds <- args[1]
output_parquet <- args[2]
batch_column <- ifelse(length(args) >= 3, args[3], "patient")
use_scale_data <- ifelse(length(args) >= 4, as.logical(args[4]), FALSE)

cat("=== ComBat Batch Correction ===\n")
cat("Input RDS file:", input_rds, "\n")
cat("Output parquet file:", output_parquet, "\n")
cat("Batch column:", batch_column, "\n")
cat("Use scaled data:", use_scale_data, "\n")

# Check if input file exists
if (!file.exists(input_rds)) {
  stop("Input RDS file does not exist: ", input_rds)
}

# Create output directory if needed
output_dir <- dirname(output_parquet)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Load Seurat object
cat("Loading Seurat object...\n")
seurat_object <- readRDS(input_rds)
cat("Loaded Seurat object with", ncol(seurat_object), "cells and", nrow(seurat_object), "genes\n")

# Ensure the patient columns is properly set 
seurat_object$patient <- sapply(colnames(seurat_object), function(x) {
  # Remove everything up to the last "parquet." 
  after_parquet <- sub(".*parquet\\.", "", x)
  # Split by "_" and take the first part (TNBC1)
  strsplit(after_parquet, "_")[[1]][1]
})

# Check if batch column exists
if (!batch_column %in% colnames(seurat_object@meta.data)) {
  stop("Batch column '", batch_column, "' not found in metadata. Available columns: ", 
       paste(colnames(seurat_object@meta.data), collapse = ", "))
}

# Extract batch information
batch <- seurat_object@meta.data[[batch_column]]
cat("Found", length(unique(batch)), "unique batches:", paste(unique(batch), collapse = ", "), "\n")

# Check if we have enough batches for ComBat
if (length(unique(batch)) < 2) {
  stop("ComBat requires at least 2 batches, but found only ", length(unique(batch)))
}

# Extract data matrix
if (use_scale_data) {
  # Check for Seurat v5 vs v4 compatibility
  if (inherits(seurat_object@assays$RNA, "Assay5")) {
    # Seurat v5 - use LayerData function (this may take several minutes for large datasets)
    cat("Extracting scaled data from Seurat v5 object (this may take a while)...\n")
    scaled_data <- LayerData(seurat_object, assay = "RNA", layer = "scale.data")
    if (is.null(scaled_data) || nrow(scaled_data) == 0) {
      cat("Scaled data not available, scaling data first...\n")
      seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
      cat("Re-extracting scaled data after scaling...\n")
      scaled_data <- LayerData(seurat_object, assay = "RNA", layer = "scale.data")
    }
    cat("Converting to matrix format...\n")
    data_matrix <- as.matrix(scaled_data)
  } else {
    # Seurat v4 - use traditional slot access
    if (is.null(seurat_object@assays$RNA@scale.data) || nrow(seurat_object@assays$RNA@scale.data) == 0) {
      cat("Scaled data not available, scaling data first...\n")
      seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
    }
    data_matrix <- as.matrix(seurat_object@assays$RNA@scale.data)
  }
  cat("Using scaled data matrix with", nrow(data_matrix), "genes\n")
} else {
  # Use normalized counts or raw counts
  if (inherits(seurat_object@assays$RNA, "Assay5")) {
    # Seurat v5 - use LayerData function (this may take several minutes for large datasets)
    cat("Extracting normalized data from Seurat v5 object (this may take a while)...\n")
    norm_data <- LayerData(seurat_object, assay = "RNA", layer = "data")
    if (!is.null(norm_data) && sum(norm_data) > 0) {
      cat("Converting normalized data to matrix format...\n")
      data_matrix <- as.matrix(norm_data)
      cat("Using normalized count data\n")
    } else {
      cat("Extracting raw count data from Seurat v5 object...\n")
      count_data <- LayerData(seurat_object, assay = "RNA", layer = "counts")
      cat("Converting raw count data to matrix format...\n")
      data_matrix <- as.matrix(count_data)
      cat("Using raw count data\n")
    }
  } else {
    # Seurat v4 - use traditional slot access
    if (!is.null(seurat_object@assays$RNA@data) && sum(seurat_object@assays$RNA@data) > 0) {
      data_matrix <- as.matrix(seurat_object@assays$RNA@data)
      cat("Using normalized count data\n")
    } else {
      data_matrix <- as.matrix(seurat_object@assays$RNA@counts)
      cat("Using raw count data\n")
    }
  }
}

cat("Data matrix dimensions:", nrow(data_matrix), "genes x", ncol(data_matrix), "cells\n")

# Apply ComBat batch correction
cat("Applying ComBat batch correction...\n")
cat("This may take several minutes for large datasets...\n")

start_time <- Sys.time()

# ComBat correction
combat_corrected_data <- ComBat(
  dat = data_matrix, 
  batch = batch,
  mod = NULL,  # No covariates to preserve
  par.prior = TRUE,  # Use parametric adjustments
  prior.plots = FALSE  # Don't create plots
)

end_time <- Sys.time()
cat("ComBat correction completed in", round(as.numeric(end_time - start_time, units = "mins"), 2), "minutes\n")

# Convert to data.frame for parquet export
cat("Converting to data.frame for parquet export...\n")
combat_df <- as.data.frame(t(combat_corrected_data))  # Transpose so cells are rows
combat_df$spot_id <- rownames(combat_df)  # Add spot IDs as a column

# Reorder columns to have spot_id first
combat_df <- combat_df[, c("spot_id", setdiff(colnames(combat_df), "spot_id"))]

cat("Final data.frame dimensions:", nrow(combat_df), "cells x", ncol(combat_df)-1, "genes\n")

# Save as parquet
cat("Saving ComBat-corrected data to parquet...\n")
write_parquet(combat_df, output_parquet)

# Verify the saved file
if (file.exists(output_parquet)) {
  file_size <- file.size(output_parquet) / (1024^2)  # Size in MB
  cat("✅ Successfully saved ComBat-corrected data to:", output_parquet, "\n")
  cat("File size:", round(file_size, 2), "MB\n")
  
  # Quick verification by reading a few rows
  cat("Verification: Reading first few rows...\n")
  test_read <- read_parquet(output_parquet, n_max = 3)
  cat("Sample data shape:", nrow(test_read), "x", ncol(test_read), "\n")
  cat("Sample spot IDs:", paste(head(test_read$spot_id, 3), collapse = ", "), "\n")
} else {
  stop("❌ Failed to create output file: ", output_parquet)
}

cat("ComBat batch correction completed successfully!\n")