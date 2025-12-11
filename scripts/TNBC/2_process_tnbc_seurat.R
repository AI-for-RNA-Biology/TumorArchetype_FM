#!/usr/bin/env Rscript

# Memory-efficient processing of TNBC gene expression data
# Processes patient-wise filtering, bimodal fitting, and Seurat object creation

# Set CRAN mirror first
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install missing packages
if (!require("arrow", quietly = TRUE)) {
  install.packages("arrow")
}
if (!require("here", quietly = TRUE)) {
  install.packages("here")
}
if (!require("stringi", quietly = TRUE)) {
  install.packages("stringi")
}

# Check R library paths and fix permissions
cat("Checking R library paths...\n")
lib_paths <- .libPaths()
cat("Current .libPaths():", paste(lib_paths, collapse = ", "), "\n")

# Ensure user library directory exists and is writable
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") {
  user_lib <- file.path(Sys.getenv("HOME"), "R", paste0("R-", getRversion()[1,1], ".", getRversion()[1,2]), "library")
}
cat("User library path:", user_lib, "\n")

if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  cat("Created user library directory\n")
}

# Add user library to the beginning of library paths
if (!user_lib %in% .libPaths()) {
  .libPaths(c(user_lib, .libPaths()))
  cat("Added user library to .libPaths()\n")
}

# Try to load Seurat and reinstall if needed
seurat_loaded <- FALSE
tryCatch({
  library(Seurat)
  seurat_loaded <- TRUE
  cat("Seurat loaded successfully!\n")
}, error = function(e) {
  cat("Error loading Seurat:", e$message, "\n")
  cat("Attempting comprehensive Seurat reinstallation...\n")
  
  # Clean up any corrupted installations
  try({
    remove.packages(c("Seurat", "SeuratObject"))
    cat("Removed existing Seurat packages\n")
  }, silent = TRUE)
  
  # Update package repositories
  options(repos = c(CRAN = "https://cloud.r-project.org", 
                   BioCsoft = "https://bioconductor.org/packages/release/bioc",
                   BioCann = "https://bioconductor.org/packages/release/data/annotation"))
  
  # Install system dependencies first
  essential_deps <- c("Rcpp", "Matrix", "stringi", "stringr")
  cat("Installing essential dependencies...\n")
  for (pkg in essential_deps) {
    tryCatch({
      install.packages(pkg, dependencies = TRUE, type = "source")
      cat("Installed", pkg, "\n")
    }, error = function(e2) {
      cat("Failed to install", pkg, "- trying binary version\n")
      try(install.packages(pkg, dependencies = TRUE), silent = TRUE)
    })
  }
  
  # Install Seurat dependencies in order
  seurat_deps <- c("sp", "future", "future.apply", "progressr", "ggplot2", 
                  "cowplot", "patchwork", "rlang", "scales", "tibble", "dplyr",
                  "RcppArmadillo", "RcppEigen", "RcppProgress", "irlba", "uwot")
  
  cat("Installing Seurat dependencies...\n")
  for (pkg in seurat_deps) {
    tryCatch({
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        cat("Installed", pkg, "\n")
      }
    }, error = function(e2) {
      cat("Warning: Could not install", pkg, "\n")
    })
  }
  
  # Install SeuratObject first
  cat("Installing SeuratObject...\n")
  tryCatch({
    install.packages("SeuratObject", dependencies = TRUE)
    library(SeuratObject)
    cat("SeuratObject installed and loaded successfully\n")
  }, error = function(e3) {
    cat("Failed to install SeuratObject:", e3$message, "\n")
    stop("Cannot proceed without SeuratObject")
  })
  
  # Install Seurat
  cat("Installing Seurat...\n")
  tryCatch({
    install.packages("Seurat", dependencies = TRUE)
    library(Seurat)
    seurat_loaded <<- TRUE
    cat("Seurat installed and loaded successfully!\n")
  }, error = function(e4) {
    cat("Failed to install Seurat:", e4$message, "\n")
    
    # Last resort: try installing from GitHub
    cat("Trying installation from GitHub as last resort...\n")
    if (!require("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    try({
      remotes::install_github("satijalab/seurat", dependencies = TRUE)
      library(Seurat)
      seurat_loaded <<- TRUE
      cat("Seurat installed from GitHub successfully!\n")
    })
  })
})

if (!seurat_loaded) {
  stop("Failed to install and load Seurat after multiple attempts. Please check your R environment and conda setup.")
}
library(data.table)
library(arrow)
library(mclust)
library(Matrix)
library(here)
library(dplyr)

# Configuration
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  # Default paths
  input_method <- "parquet"  # or "csv"
  input_path <- here("results", "molecular", "TNBC", "TNBC_all_samples_coding_genes.parquet")
  intermediate_dir <- here("results", "molecular", "TNBC", "intermediate")
  output_dir <- here("results", "molecular", "TNBC")
  min_patient_fraction <- 0.10  # 10% of patients
  force_recompute <- FALSE  # Set to TRUE to force recomputation of gene filtering
  spots_filter_file <- here("results", "compute_patches", "TNBC", "spots_labels.csv")
} else {
  input_method <- args[1]
  input_path <- args[2] 
  intermediate_dir <- args[3]
  output_dir <- args[4]
  min_patient_fraction <- as.numeric(args[5])
  force_recompute <- ifelse(length(args) >= 6, as.logical(args[6]), FALSE)
  spots_filter_file <- ifelse(length(args) >= 7, args[7], here("results", "compute_patches", "TNBC", "spots_labels.csv"))
}

cat("=== TNBC Gene Expression Processing ===\n")
cat("Input method:", input_method, "\n")
cat("Input path:", input_path, "\n")
cat("Output directory:", output_dir, "\n")
cat("Minimum patient fraction:", min_patient_fraction, "\n")
cat("Force recompute gene filtering:", force_recompute, "\n")
cat("Spots filter file:", spots_filter_file, "\n")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Check if gene selection results already exist
gene_results_file <- file.path(output_dir, "gene_selection_results.csv")
skip_gene_filtering <- file.exists(gene_results_file) && !force_recompute

if (skip_gene_filtering) {
  cat("\n=== Gene filtering results already exist ===\n")
  cat("Loading existing results from:", gene_results_file, "\n")
  cat("To recompute gene filtering, set force_recompute=TRUE or delete the results file\n")
}

# Load spots filter if file exists
valid_spots <- NULL
cat("\nDEBUG: Checking spots filter file:", spots_filter_file, "\n")
cat("DEBUG: File exists:", file.exists(spots_filter_file), "\n")

if (file.exists(spots_filter_file)) {
  cat("\n=== Loading spots filter ===\n")
  cat("Spots filter file:", spots_filter_file, "\n")
  spots_filter_df <- read.csv(spots_filter_file, stringsAsFactors = FALSE)
  cat("DEBUG: spots_filter_df dimensions:", nrow(spots_filter_df), "x", ncol(spots_filter_df), "\n")
  cat("DEBUG: Column names:", paste(colnames(spots_filter_df), collapse = ", "), "\n")
  valid_spots <- spots_filter_df[, 1]  # First column contains spot names
  cat("Found", length(valid_spots), "valid spots in filter file\n")
  cat("Sample spot names:", paste(head(valid_spots, 5), collapse = ", "), "...\n")
  cat("DEBUG: valid_spots is null:", is.null(valid_spots), "\n")
  cat("DEBUG: valid_spots length:", length(valid_spots), "\n")
} else {
  cat("\n=== No spots filter file found ===\n")
  cat("Processing all spots (no filtering applied)\n")
}

# Function to extract patient ID from spot names
extract_patient_id <- function(spot_names) {
  # Assuming format like "TNBC64_spot_id" -> patient is "TNBC64"
  sapply(strsplit(spot_names, "_"), function(x) x[1])
}

# Function to normalize spot names for consistent matching
normalize_spot_names <- function(spot_names) {
  # Convert from "TNBC1_X2x14" to "TNBC1_spot2x14" format
  # Or from "TNBC1_spot2x14" to "TNBC1_X2x14" format
  normalized <- gsub("_X([0-9]+)x([0-9]+)", "_spot\\1x\\2", spot_names)
  return(normalized)
}

# Function to filter spots based on valid spots list
filter_valid_spots <- function(spot_names, valid_spots_list = NULL) {
  if (is.null(valid_spots_list)) {
    return(rep(TRUE, length(spot_names)))  # Keep all spots if no filter
  }
  
  # Normalize both spot names to the same format for comparison
  normalized_spot_names <- normalize_spot_names(spot_names)
  normalized_valid_spots <- normalize_spot_names(valid_spots_list)
  
  # Check which spots are in the valid list after normalization
  is_valid <- normalized_spot_names %in% normalized_valid_spots
  
  if (sum(is_valid) == 0) {
    cat("WARNING: No spots match the filter list!\n")
    cat("Sample current spots:", paste(head(spot_names, 3), collapse = ", "), "\n")
    cat("Sample normalized current:", paste(head(normalized_spot_names, 3), collapse = ", "), "\n")
    cat("Sample filter spots:", paste(head(valid_spots_list, 3), collapse = ", "), "\n")
    cat("Sample normalized filter:", paste(head(normalized_valid_spots, 3), collapse = ", "), "\n")
  } else {
    cat("DEBUG: Found", sum(is_valid), "matching spots after normalization\n")
  }
  
  return(is_valid)
}

# Function to create pseudobulk profile by summing counts
create_pseudobulk <- function(count_matrix) {
  if (is(count_matrix, "Matrix")) {
    pseudobulk <- Matrix::rowSums(count_matrix)
  } else {
    pseudobulk <- rowSums(count_matrix)
  }
  return(pseudobulk)
}

# Function to find bimodal cutoff
find_bimodal_cutoff <- function(pseudobulk_log, patient_id) {
  # Keep original vector with gene names for final filtering
  original_pseudobulk <- pseudobulk_log
  
  # Remove zeros AND NAs for bimodal fitting
  pseudobulk_filtered <- pseudobulk_log[pseudobulk_log > 0 & !is.na(pseudobulk_log)]
  
  if (length(pseudobulk_filtered) < 50) {
    cat("Warning: Patient", patient_id, "has too few non-zero values\n")
    return(list(cutoff = 1, genes_above = character(0)))
  }
  
  tryCatch({
    # Fit bimodal distribution on filtered data
    bimdens <- densityMclust(data = pseudobulk_filtered, G = 2, plot = FALSE)
    
    # Identify cutoff (90th percentile of first component)
    cutoff <- qnorm(0.99, 
                   mean = bimdens$parameters$mean[1], 
                   sd = sqrt(bimdens$parameters$variance$sigmasq[1]))
    
    # Debug: Check the data before applying cutoff
    cat("    DEBUG - Original vector length:", length(original_pseudobulk), 
        ", range:", round(range(original_pseudobulk, na.rm = TRUE), 3), 
        ", cutoff:", round(cutoff, 3), "\n")
    cat("    DEBUG - Names present:", !is.null(names(original_pseudobulk)), 
        ", names length:", length(names(original_pseudobulk)), "\n")
    
    # Get genes above cutoff from ORIGINAL vector (with gene names intact)
    above_cutoff_indices <- which(original_pseudobulk > cutoff & !is.na(original_pseudobulk))
    genes_above <- names(original_pseudobulk)[above_cutoff_indices]
    # Remove any NAs from gene names
    genes_above <- genes_above[!is.na(genes_above)]
    
    cat("    DEBUG - Indices above cutoff:", length(above_cutoff_indices), 
        ", values above cutoff:", sum(original_pseudobulk > cutoff, na.rm = TRUE), "\n")
    
    cat("Patient", patient_id, ": cutoff =", round(cutoff, 3), 
        ", genes above =", length(genes_above), "\n")
    
    return(list(cutoff = cutoff, genes_above = genes_above))
    
  }, error = function(e) {
    cat("Error fitting bimodal distribution for patient", patient_id, ":", e$message, "\n")
    cat("Using median as cutoff\n")
    cutoff <- median(pseudobulk_filtered)
    
    # Debug info for error case too
    cat("    DEBUG (ERROR) - Original vector length:", length(original_pseudobulk), 
        ", cutoff:", round(cutoff, 3), "\n")
    
    # Use original vector for gene selection
    above_cutoff_indices <- which(original_pseudobulk > cutoff & !is.na(original_pseudobulk))
    genes_above <- names(original_pseudobulk)[above_cutoff_indices]
    genes_above <- genes_above[!is.na(genes_above)]
    
    cat("    DEBUG (ERROR) - Values above cutoff:", length(above_cutoff_indices), "\n")
    
    return(list(cutoff = cutoff, genes_above = genes_above))
  })
}

# Method 1: Process from intermediate parquet files
process_from_parquet_files <- function(intermediate_dir) {
  cat("\n=== Processing from intermediate parquet files ===\n")
  
  # Find all batch files
  batch_files <- list.files(intermediate_dir, pattern = "batch_.*\\.parquet$", full.names = TRUE)
  cat("Found", length(batch_files), "batch files\n")
  
  if (length(batch_files) == 0) {
    stop("No batch files found in ", intermediate_dir)
  }
  
  # Process each batch and collect patient-wise results
  all_patients <- list()
  
  for (batch_file in batch_files) {
    cat("Processing", basename(batch_file), "...\n")
    
    # Read batch and convert to data.frame (parquet often returns tibbles)
    batch_df <- as.data.frame(read_parquet(batch_file))
    
    cat("  Batch shape:", nrow(batch_df), "x", ncol(batch_df), "\n")
    cat("  Column names (first 5):", paste(head(colnames(batch_df), 5), collapse = ", "), "...\n")
    
    # Handle index column as row names if present
    index_cols <- grep("^__index_level_0__|^index$|^Unnamed", colnames(batch_df), value = TRUE)
    if (length(index_cols) > 0) {
      cat("  Found index column(s):", paste(index_cols, collapse = ", "), "\n")
      # Use first index column as row names
      rownames(batch_df) <- as.character(batch_df[[index_cols[1]]])
      # Remove index columns from data
      batch_df <- batch_df[, !colnames(batch_df) %in% index_cols, drop = FALSE]
      cat("  After removing index column(s):", nrow(batch_df), "x", ncol(batch_df), "\n")
    }
    
    # Ensure all remaining columns are numeric
    numeric_cols <- sapply(batch_df, is.numeric)
    if (!all(numeric_cols)) {
      cat("  Warning: Non-numeric data columns detected, converting...\n")
      non_numeric_cols <- names(batch_df)[!numeric_cols]
      cat("  Non-numeric columns:", paste(head(non_numeric_cols, 5), collapse = ", "), "\n")
      
      # Convert non-numeric columns to numeric (coerce to 0 if can't convert)
      for (col in non_numeric_cols) {
        batch_df[[col]] <- as.numeric(as.character(batch_df[[col]]))
        batch_df[[col]][is.na(batch_df[[col]])] <- 0
      }
    }
    
    # Extract patient IDs from row names
    spot_names <- rownames(batch_df)
    patient_ids <- extract_patient_id(spot_names)
    
    # Process each patient in this batch
    unique_patients <- unique(patient_ids)
    cat("  Processing", length(unique_patients), "patients:", paste(unique_patients, collapse = ", "), "\n")
    
    for (i in seq_along(unique_patients)) {
      patient_id <- unique_patients[i]
      cat("    Patient", i, "of", length(unique_patients), ":", patient_id, "... ")
      
      patient_spots <- spot_names[patient_ids == patient_id]
      cat(length(patient_spots), "spots... ")
      
      # Use more efficient subsetting - avoid as.matrix which can be slow
      patient_data <- batch_df[patient_spots, , drop = FALSE]
      
      # Quick check - if already numeric data.frame, skip conversion
      if (is.data.frame(patient_data) && all(sapply(patient_data, is.numeric))) {
        # Convert to matrix more efficiently
        patient_matrix <- as.matrix(patient_data)
      } else {
        cat("converting...")
        patient_matrix <- as.matrix(patient_data)
        if (!is.numeric(patient_matrix)) {
          patient_matrix <- apply(patient_matrix, 2, as.numeric)
          patient_matrix[is.na(patient_matrix)] <- 0
        }
      }
      
      # Remove spots with zero counts - use faster method with NA handling
      spot_sums <- .rowSums(patient_matrix, nrow(patient_matrix), ncol(patient_matrix), na.rm = TRUE)
      keep_spots <- !is.na(spot_sums) & spot_sums > 0
      
      if (sum(keep_spots, na.rm = TRUE) == 0) {
        cat("no valid spots\n")
        next
      }
      
      patient_matrix <- patient_matrix[keep_spots, , drop = FALSE]
      cat("filtered to", nrow(patient_matrix), "spots... ")
      
      # Create pseudobulk more efficiently and preserve gene names
      pseudobulk <- .colSums(patient_matrix, nrow(patient_matrix), ncol(patient_matrix))
      names(pseudobulk) <- colnames(patient_matrix)  # Preserve gene names!
      
      # Debug: check for issues before log transformation
      if (any(is.na(pseudobulk))) {
        cat("WARNING: NAs in pseudobulk for", patient_id, "- converting to 0\n")
        pseudobulk[is.na(pseudobulk)] <- 0
      }
      if (any(pseudobulk < 0)) {
        cat("WARNING: Negative values in pseudobulk for", patient_id, "- converting to 0\n")
        pseudobulk[pseudobulk < 0] <- 0
      }
      
      pseudobulk_log <- log2(pseudobulk + 1)
      
      # Store for this patient
      if (patient_id %in% names(all_patients)) {
        # Adding to existing patient data
        combined_data <- all_patients[[patient_id]] + pseudobulk_log
        if (any(is.na(combined_data))) {
          cat("WARNING: NAs after combining data for", patient_id, "\n")
          combined_data[is.na(combined_data)] <- 0
        }
        all_patients[[patient_id]] <- combined_data
      } else {
        all_patients[[patient_id]] <- pseudobulk_log
      }
      
      cat("done\n")
      
      # Force garbage collection every 10 patients to prevent memory buildup
      if (i %% 10 == 0) {
        gc()
      }
    }
    
    # Clear memory
    rm(batch_df)
    gc()
  }
  
  return(all_patients)
}

# Method 2: Process from big CSV/parquet file (more memory intensive)
process_from_big_file <- function(input_path) {
  cat("\n=== Processing from big file ===\n")
  cat("Reading file:", input_path, "\n")
  
  # Read the big matrix
  if (grepl("\\.parquet$", input_path)) {
    big_df <- as.data.frame(read_parquet(input_path))
  } else {
    big_df <- fread(input_path, data.table = FALSE)
    rownames(big_df) <- big_df[, 1]
    big_df <- big_df[, -1]
  }
  
  cat("Matrix dimensions:", nrow(big_df), "x", ncol(big_df), "\n")
  
  # Handle index column as row names if present (for parquet files)
  if (grepl("\\.parquet$", input_path)) {
    index_cols <- grep("^__index_level_0__|^index$|^Unnamed", colnames(big_df), value = TRUE)
    if (length(index_cols) > 0) {
      cat("Found index column(s):", paste(index_cols, collapse = ", "), "\n")
      # Use first index column as row names
      rownames(big_df) <- as.character(big_df[[index_cols[1]]])
      # Remove index columns from data
      big_df <- big_df[, !colnames(big_df) %in% index_cols, drop = FALSE]
      cat("After removing index column(s):", nrow(big_df), "x", ncol(big_df), "\n")
    }
  }
  
  # Ensure all remaining columns are numeric
  numeric_cols <- sapply(big_df, is.numeric)
  if (!all(numeric_cols)) {
    cat("Warning: Non-numeric data columns detected, converting...\n")
    non_numeric_cols <- names(big_df)[!numeric_cols]
    cat("Non-numeric columns:", paste(head(non_numeric_cols, 5), collapse = ", "), "\n")
    
    # Convert non-numeric columns to numeric
    for (col in non_numeric_cols) {
      big_df[[col]] <- as.numeric(as.character(big_df[[col]]))
      big_df[[col]][is.na(big_df[[col]])] <- 0
    }
  }
  
  # Extract patient IDs
  spot_names <- rownames(big_df)
  patient_ids <- extract_patient_id(spot_names)
  unique_patients <- unique(patient_ids)
  
  cat("Found", length(unique_patients), "patients:", paste(unique_patients, collapse = ", "), "\n")
  
  # Process each patient
  all_patients <- list()
  
  for (patient_id in unique_patients) {
    cat("Processing patient", patient_id, "...\n")
    
    patient_spots <- spot_names[patient_ids == patient_id]
    patient_data <- big_df[patient_spots, , drop = FALSE]
    
    # Ensure patient_data is numeric matrix
    patient_data <- as.matrix(patient_data)
    if (!is.numeric(patient_data)) {
      cat("  Warning: Converting patient data to numeric for", patient_id, "\n")
      patient_data <- apply(patient_data, 2, as.numeric)
      patient_data[is.na(patient_data)] <- 0
    }
    
    # Remove spots with zero counts
    spot_sums <- rowSums(patient_data, na.rm = TRUE)
    patient_data <- patient_data[spot_sums > 0, , drop = FALSE]
    
    if (nrow(patient_data) == 0) {
      cat("Warning: No valid spots for patient", patient_id, "\n")
      next
    }
    
    cat("  Patient", patient_id, ":", nrow(patient_data), "spots\n")
    
    # Create pseudobulk
    pseudobulk <- create_pseudobulk(t(patient_data))
    pseudobulk_log <- log2(pseudobulk + 1)
    
    all_patients[[patient_id]] <- pseudobulk_log
    
    # Clear memory
    rm(patient_data)
    gc()
  }
  
  return(all_patients)
}

################## Main processing ##################
cat("\n=== Starting processing ===\n")

if (skip_gene_filtering) {
  # Load existing gene selection results
  cat("Loading existing gene selection results...\n")
  gene_selection_results <- read.csv(gene_results_file, stringsAsFactors = FALSE)
  selected_genes <- gene_selection_results$gene[gene_selection_results$selected]
  cat("Loaded", length(selected_genes), "selected genes from existing results\n")
  
  # Also load patient cutoffs if available
  cutoff_file <- file.path(output_dir, "patient_cutoffs.csv")
  if (file.exists(cutoff_file)) {
    cutoff_results <- read.csv(cutoff_file, stringsAsFactors = FALSE)
    cat("Loaded patient cutoffs for", nrow(cutoff_results), "patients\n")
  }
  
} else {
  # Perform gene filtering computation
  cat("Computing gene filtering from scratch...\n")

  # Get patient-wise pseudobulk data
  if (input_method == "parquet_batch") {
    all_patients <- process_from_parquet_files(intermediate_dir)
  } else {
    all_patients <- process_from_big_file(input_path)
  }

  # Find cutoffs and genes above limit for each patient
  cat("\n=== Finding bimodal cutoffs ===\n")
  patient_results <- list()
  genes_above_limit <- list()

  # Create PDF for filtering plots
  pdf(file.path(output_dir, "TNBC_filtering_plots.pdf"), width = 16, height = 12)
  n_patients <- length(all_patients)

  # Plot maximum 20 patients per page (4x5 grid)
  plots_per_page <- 20
  n_pages <- ceiling(n_patients / plots_per_page)
  patient_names <- names(all_patients)

  for (page in 1:n_pages) {
    # Set up grid for this page (4 rows x 5 columns)
    par(mfrow = c(4, 5), mar = c(4, 4, 2, 1))
    
    # Calculate patient indices for this page
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, n_patients)
    
    cat("Creating plots for page", page, "of", n_pages, "(patients", start_idx, "to", end_idx, ")\n")
    
    for (i in start_idx:end_idx) {
      patient_id <- patient_names[i]
    result <- find_bimodal_cutoff(all_patients[[patient_id]], patient_id)
    patient_results[[patient_id]] <- result
    genes_above_limit[[patient_id]] <- result$genes_above
    
    # Create the density plot for this patient
    pseudobulk_log <- all_patients[[patient_id]]
    # Remove zeros AND NAs for plotting
    pseudobulk_log <- pseudobulk_log[pseudobulk_log > 0 & !is.na(pseudobulk_log)]
    
    if (length(pseudobulk_log) > 0 && !any(is.na(pseudobulk_log))) {
      # Plot the density
      plot(density(pseudobulk_log), 
           main = patient_id,
           xlab = "log2(# genes per spots + 1)", 
           ylab = "Density",
           col = "darkblue",
           lwd = 2)
      
      # Add the vertical line for the cutoff
      abline(v = result$cutoff, col = "red", lwd = 2.5, lty = 2)
      
      # Add text with cutoff value
      text(x = result$cutoff, y = max(density(pseudobulk_log)$y) * 0.8, 
           labels = paste("Cutoff:", round(result$cutoff, 3)), 
           pos = 4, col = "red", cex = 0.8)
    } else {
      # Create empty plot if no valid data
      plot(1, 1, type = "n", main = paste(patient_id, "(no valid data)"), 
           xlab = "log2(# genes per spots + 1)", ylab = "Density")
      text(1, 1, "No valid data for plotting", cex = 1.2, col = "red")
    }
    }
    
    # Fill remaining spaces on last page with empty plots if needed
    if (page == n_pages && end_idx < page * plots_per_page) {
      remaining_slots <- plots_per_page - (end_idx - start_idx + 1)
      for (j in 1:remaining_slots) {
        plot.new()
      }
    }
  }

  dev.off()
  cat("Filtering plots saved to:", file.path(output_dir, "TNBC_filtering_plots.pdf"), "\n")

  # Find genes above limit in at least min_patient_fraction of patients
  cat("\n=== Selecting common genes ===\n")
  n_patients <- length(all_patients)
  min_patients <- ceiling(n_patients * min_patient_fraction)

  cat("Total patients:", n_patients, "\n")
  cat("Minimum patients for gene inclusion:", min_patients, "\n")

  # Count how many patients each gene appears in
  all_genes <- unique(unlist(genes_above_limit))
  gene_counts <- sapply(all_genes, function(gene) {
    sum(sapply(genes_above_limit, function(patient_genes) gene %in% patient_genes))
  })

  # Select genes that appear in at least min_patients
  selected_genes <- names(gene_counts)[gene_counts >= min_patients]
  cat("Selected", length(selected_genes), "genes (out of", length(all_genes), "total)\n")

  # Save gene selection results
  gene_selection_results <- data.frame(
    gene = names(gene_counts),
    n_patients = gene_counts,
    selected = gene_counts >= min_patients
  )
  write.csv(gene_selection_results, file.path(output_dir, "gene_selection_results.csv"))

  # Save patient cutoffs
  cutoff_results <- data.frame(
    patient = names(patient_results),
    cutoff = sapply(patient_results, function(x) x$cutoff),
    n_genes_above = sapply(patient_results, function(x) length(x$genes_above))
  )
  write.csv(cutoff_results, file.path(output_dir, "patient_cutoffs.csv"))

}

# Now create the filtered matrix with selected genes
cat("\n=== Creating filtered matrix ===\n")

# Re-read the data and filter for selected genes only
if (input_method == "parquet_batch") {
  # Process batch files again, but only keep selected genes
  cat("Re-processing batch files for selected genes...\n")
  
  filtered_data_list <- list()
  
  for (batch_file in list.files(intermediate_dir, pattern = "batch_.*\\.parquet$", full.names = TRUE)) {
    cat("Re-processing", basename(batch_file), "...\n")
    
    batch_df <- as.data.frame(read_parquet(batch_file))
    
    # Handle index column as row names if present
    index_cols <- grep("^__index_level_0__|^index$|^Unnamed", colnames(batch_df), value = TRUE)
    if (length(index_cols) > 0) {
      cat("  Found index column(s):", paste(index_cols, collapse = ", "), "\n")
      # Use first index column as row names
      rownames(batch_df) <- as.character(batch_df[[index_cols[1]]])
      # Remove index columns from data
      batch_df <- batch_df[, !colnames(batch_df) %in% index_cols, drop = FALSE]
    }
    
    # Ensure all columns are numeric
    numeric_cols <- sapply(batch_df, is.numeric)
    if (!all(numeric_cols)) {
      non_numeric_cols <- names(batch_df)[!numeric_cols]
      for (col in non_numeric_cols) {
        batch_df[[col]] <- as.numeric(as.character(batch_df[[col]]))
        batch_df[[col]][is.na(batch_df[[col]])] <- 0
      }
    }
    
    # Keep all spots for now - filtering will happen later
    
    # Keep only selected genes
    available_genes <- intersect(selected_genes, colnames(batch_df))
    batch_df_filtered <- batch_df[, available_genes, drop = FALSE]
    
    # Add missing genes as zeros
    missing_genes <- setdiff(selected_genes, colnames(batch_df_filtered))
    if (length(missing_genes) > 0) {
      missing_matrix <- matrix(0, nrow = nrow(batch_df_filtered), ncol = length(missing_genes))
      colnames(missing_matrix) <- missing_genes
      rownames(missing_matrix) <- rownames(batch_df_filtered)
      batch_df_filtered <- cbind(batch_df_filtered, missing_matrix)
    }
    
    # Ensure column order
    batch_df_filtered <- batch_df_filtered[, selected_genes, drop = FALSE]
    
    # Apply final spot filtering here (after gene filtering)
    cat("  DEBUG: valid_spots is null:", is.null(valid_spots), "\n")
    cat("  DEBUG: About to check spot filtering condition\n")
    if (!is.null(valid_spots)) {
      cat("  DEBUG: Entering spot filtering loop\n")
      batch_spot_names <- rownames(batch_df_filtered)
      cat("  DEBUG: batch has", length(batch_spot_names), "spots\n")
      valid_mask <- filter_valid_spots(batch_spot_names, valid_spots)
      if (sum(valid_mask) > 0) {
        batch_df_filtered <- batch_df_filtered[valid_mask, , drop = FALSE]
        cat("  Final spot filtering: kept", nrow(batch_df_filtered), "valid spots (from", length(valid_mask), "total)\n")
      } else {
        cat("  WARNING: No valid spots found in this batch after filtering\n")
      }
    } else {
      cat("  DEBUG: valid_spots is NULL, skipping spot filtering\n")
    }
    
    filtered_data_list[[basename(batch_file)]] <- batch_df_filtered
    
    rm(batch_df)
    gc()
  }
  
  # Combine all batches
  cat("Combining filtered batches...\n")
  filtered_matrix <- do.call(rbind, filtered_data_list)
  
} else {
  # Filter the big matrix
  cat("Filtering big matrix...\n")
  
  if (grepl("\\.parquet$", input_path)) {
    big_df <- read_parquet(input_path)
  } else {
    big_df <- fread(input_path, data.table = FALSE)
    rownames(big_df) <- big_df[, 1]
    big_df <- big_df[, -1]
  }
  
  # Keep all spots for now - filtering will happen later
  
  # Keep only selected genes
  available_genes <- intersect(selected_genes, colnames(big_df))
  filtered_matrix <- big_df[, available_genes, drop = FALSE]
  
  # Add missing genes as zeros
  missing_genes <- setdiff(selected_genes, colnames(filtered_matrix))
  if (length(missing_genes) > 0) {
    missing_matrix <- matrix(0, nrow = nrow(filtered_matrix), ncol = length(missing_genes))
    colnames(missing_matrix) <- missing_genes
    rownames(missing_matrix) <- rownames(filtered_matrix)
    filtered_matrix <- cbind(filtered_matrix, missing_matrix)
  }
  
  # Ensure column order
  filtered_matrix <- filtered_matrix[, selected_genes, drop = FALSE]
  
  # Apply final spot filtering here (after gene filtering)
  cat("DEBUG: valid_spots is null:", is.null(valid_spots), "\n")
  if (!is.null(valid_spots)) {
    cat("DEBUG: Entering spot filtering for big file method\n")
    matrix_spot_names <- rownames(filtered_matrix)
    cat("DEBUG: matrix has", length(matrix_spot_names), "spots\n")
    valid_mask <- filter_valid_spots(matrix_spot_names, valid_spots)
    if (sum(valid_mask) > 0) {
      filtered_matrix <- filtered_matrix[valid_mask, , drop = FALSE]
      cat("Final spot filtering: kept", nrow(filtered_matrix), "valid spots (from", length(valid_mask), "total)\n")
    } else {
      cat("WARNING: No valid spots found after filtering\n")
    }
  } else {
    cat("DEBUG: valid_spots is NULL, skipping spot filtering in big file method\n")
  }
  
  rm(big_df)
  gc()
}

cat("Filtered matrix dimensions:", nrow(filtered_matrix), "x", ncol(filtered_matrix), "\n")

# Clean up any NaN values that may have been introduced during concatenation
cat("Checking for NaN/NA values in filtered matrix...\n")
na_count <- sum(is.na(filtered_matrix))
nan_count <- sum(is.nan(as.matrix(filtered_matrix)))
if (na_count > 0 || nan_count > 0) {
  cat("Found", na_count, "NA values and", nan_count, "NaN values - replacing with zeros\n")
  filtered_matrix[is.na(filtered_matrix)] <- 0
  filtered_matrix[is.nan(as.matrix(filtered_matrix))] <- 0
} else {
  cat("No NA/NaN values found\n")
}

# Ensure all values are numeric
if (!is.numeric(as.matrix(filtered_matrix))) {
  cat("Converting filtered matrix to fully numeric format...\n")
  filtered_matrix <- apply(filtered_matrix, 2, as.numeric)
  filtered_matrix[is.na(filtered_matrix)] <- 0
}

matrix_nrow <- as.numeric(nrow(filtered_matrix))
matrix_ncol <- as.numeric(ncol(filtered_matrix))
matrix_size_gb <- matrix_nrow * matrix_ncol * 8 / 1e9
cat("Estimated memory usage:", round(matrix_size_gb, 2), "GB\n")

# Check if matrix is too large for available memory
if (!is.na(matrix_size_gb) && matrix_size_gb > 25) {  # Conservative limit to leave room for Seurat operations
  cat("WARNING: Matrix is very large (", round(matrix_size_gb, 1), "GB)\n")
  cat("Consider increasing memory or reducing the number of genes/spots\n")
  # Sample spots to reduce memory if needed
}

# Remove spots with zero counts after filtering (memory-efficient)
cat("Removing zero-count spots...\n")
# Use more memory-efficient row sum calculation
spot_sums <- .rowSums(as.matrix(filtered_matrix), nrow(filtered_matrix), ncol(filtered_matrix), na.rm = TRUE)
# Handle potential NAs and filter in one step
keep_spots <- !is.na(spot_sums) & spot_sums > 0
spots_removed <- sum(keep_spots) < nrow(filtered_matrix)
if (spots_removed) {
  cat("Removing", sum(!keep_spots), "spots with zero counts\n")
  filtered_matrix <- filtered_matrix[keep_spots, , drop = FALSE]
  gc()  # Force garbage collection after filtering
} else {
  cat("No zero-count spots to remove\n")
}
cat("After removing zero-count spots:", nrow(filtered_matrix), "x", ncol(filtered_matrix), "\n")

# Remove genes with zero variance to prevent Seurat issues (BEFORE creating Seurat object)
cat("Checking gene variance and expression...\n")
gene_sums <- .colSums(as.matrix(filtered_matrix), nrow(filtered_matrix), ncol(filtered_matrix), na.rm = TRUE)
gene_vars <- apply(as.matrix(filtered_matrix), 2, var, na.rm = TRUE)
keep_genes <- !is.na(gene_vars) & gene_vars > 0 & gene_sums > 0
genes_removed <- sum(!keep_genes)
if (genes_removed > 0) {
  cat("Removing", genes_removed, "genes with zero variance or zero expression\n")
  filtered_matrix <- filtered_matrix[, keep_genes, drop = FALSE]
  cat("After removing problematic genes:", nrow(filtered_matrix), "x", ncol(filtered_matrix), "\n")
} else {
  cat("All genes have variance > 0\n")
}

# Final check - ensure we have enough genes for analysis
if (ncol(filtered_matrix) < 100) {
  cat("WARNING: Very few genes remaining (", ncol(filtered_matrix), "). This may cause issues with downstream analysis.\n")
}

# Force garbage collection before Seurat operations
gc()

# Create metadata
spot_names <- rownames(filtered_matrix)
patient_ids <- extract_patient_id(spot_names)

metadata <- data.frame(
  spot_id = spot_names,
  patient = patient_ids,
  row.names = spot_names
)

# Create Seurat object
cat("\n=== Creating Seurat object ===\n")

# Convert to sparse matrix to save memory
library(Matrix)

# Save filtered count matrix WITH SPOT NAMES PRESERVED
cat("DEBUG: Saving filtered counts with spot names...\n")
cat("DEBUG: filtered_matrix has rownames:", !is.null(rownames(filtered_matrix)), "\n")
cat("DEBUG: Sample rownames:", paste(head(rownames(filtered_matrix), 3), collapse = ", "), "\n")

filtered_counts_df <- as.data.frame(filtered_matrix)
filtered_counts_df$spot_id <- rownames(filtered_matrix)  # Add spot names as explicit column

# Additional safety check
if (is.null(filtered_counts_df$spot_id) || any(is.na(filtered_counts_df$spot_id))) {
  cat("WARNING: spot_id column has issues, using row indices instead\n")
  filtered_counts_df$spot_id <- paste0("spot_", 1:nrow(filtered_counts_df))
}

filtered_file <- file.path(output_dir, "TNBC_filtered_counts.parquet")
write_parquet(filtered_counts_df, filtered_file)
cat("Filtered counts saved to:", filtered_file, "(with", nrow(filtered_counts_df), "spots)\n")
cat("DEBUG: spot_id column sample:", paste(head(filtered_counts_df$spot_id, 3), collapse = ", "), "\n")



# For Seurat object creation, we need to exclude the spot_id column and transpose
# Remove spot_id column for matrix operations
filtered_matrix_for_seurat <- filtered_counts_df[, !colnames(filtered_counts_df) %in% "spot_id", drop = FALSE]

# Convert to matrix and transpose for Seurat (genes x spots)
filtered_matrix_seurat <- as.matrix(filtered_matrix_for_seurat)
filtered_sparse <- as(t(filtered_matrix_seurat), "dgCMatrix")

# Debug: Check data integrity
cat("DEBUG: Matrix for Seurat dimensions:", dim(filtered_sparse), "(genes x spots)\n")
cat("DEBUG: Any NAs in sparse matrix:", any(is.na(filtered_sparse@x)), "\n")
cat("DEBUG: Sample gene names:", paste(head(rownames(filtered_sparse), 3), collapse = ", "), "\n")
cat("DEBUG: Sample spot names:", paste(head(colnames(filtered_sparse), 3), collapse = ", "), "\n")

seurat_object <- CreateSeuratObject(
  counts = filtered_sparse,
  meta.data = metadata,
  project = "TNBC"
)

# Clean up intermediate objects
rm(filtered_matrix_for_seurat, filtered_matrix_seurat, filtered_sparse)
gc()

cat("Seurat object created with", ncol(seurat_object), "spots and", nrow(seurat_object), "genes\n")

# Clear intermediate objects
rm(metadata)
gc()

# Normalize data (log normalization, not scaled)
cat("Normalizing data...\n")
seurat_object <- NormalizeData(seurat_object, assay = "RNA")


# Save normalized data (not scaled) WITH SPOT NAMES PRESERVED
cat("DEBUG: Saving normalized expression with spot names...\n")
normalized_data <- as.matrix(seurat_object@assays$RNA$data)
cat("DEBUG: normalized_data dimensions:", dim(normalized_data), "(genes x spots)\n")
cat("DEBUG: colnames sample:", paste(head(colnames(normalized_data), 3), collapse = ", "), "\n")

# Note: normalized data is transposed (genes x spots), so we need column names as spot_id
normalized_df <- as.data.frame(t(normalized_data))  # Transpose back to spots x genes
normalized_df$spot_id <- rownames(normalized_df)    # Add spot names as explicit column

# Additional safety check
if (is.null(normalized_df$spot_id) || any(is.na(normalized_df$spot_id))) {
  cat("WARNING: spot_id column has issues in normalized data\n")
  cat("DEBUG: Using colnames from original normalized_data\n")
  normalized_df$spot_id <- colnames(normalized_data)  # Use original column names
}

normalized_file <- file.path(output_dir, "TNBC_filtered_normalized_expression.parquet")
write_parquet(normalized_df, normalized_file)
cat("Normalized expression saved to:", normalized_file, "(with", nrow(normalized_df), "spots)\n")
cat("DEBUG: spot_id column sample:", paste(head(normalized_df$spot_id, 3), collapse = ", "), "\n")



# Find variable features
cat("Finding variable features...\n")
tryCatch({
  seurat_object <- FindVariableFeatures(seurat_object, assay = "RNA", selection.method = "vst", nfeatures = min(19000, nrow(seurat_object)))
}, error = function(e) {
  cat("Error with VST method, trying dispersion method instead...\n")
  cat("Original error:", e$message, "\n")
  seurat_object <<- FindVariableFeatures(seurat_object, assay = "RNA", selection.method = "dispersion", nfeatures = min(19000, nrow(seurat_object)))
})

# Scale data for PCA
cat("Scaling data...\n")
# Check if we have variable features
var_features <- VariableFeatures(seurat_object)
if (length(var_features) == 0) {
  cat("WARNING: No variable features found. Using top 2000 most variable genes instead.\n")
  # Calculate variance manually and select top genes
  all_genes <- rownames(seurat_object)
  gene_vars <- apply(as.matrix(seurat_object@assays$RNA$data), 1, var, na.rm = TRUE)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(2000, length(gene_vars))]
  VariableFeatures(seurat_object) <- top_genes
  cat("Set", length(top_genes), "most variable genes as variable features\n")
}

tryCatch({
  seurat_object <- ScaleData(seurat_object, assay = "RNA")
}, error = function(e) {
  cat("Error in ScaleData:", e$message, "\n")
  cat("Trying to scale only variable features...\n")
  seurat_object <<- ScaleData(seurat_object, assay = "RNA", features = VariableFeatures(seurat_object))
})

# Add tissue labels to Seurat object metadata
cat("\n=== Adding tissue labels to metadata ===\n")

if (file.exists(spots_filter_file)) {
  # Read the spots labels file
  spots_labels_df <- read.csv(spots_filter_file, stringsAsFactors = FALSE)
  cat("Loaded spots labels from:", spots_filter_file, "\n")
  cat("Spots labels dimensions:", nrow(spots_labels_df), "x", ncol(spots_labels_df), "\n")
  cat("Spots labels columns:", paste(colnames(spots_labels_df), collapse = ", "), "\n")
  
  # Debug: Check what's in each column to find the tissue names
  cat("DEBUG: First few values in each column:\n")
  for (i in 1:ncol(spots_labels_df)) {
    col_name <- colnames(spots_labels_df)[i]
    sample_values <- paste(head(as.character(spots_labels_df[, i]), 3), collapse = ", ")
    cat("  Column", i, "(", col_name, "):", sample_values, "\n")
  }
  
  # Look for a column that might contain tissue names instead of numbers
  # Check if there's a column with text values that look like tissue names
  tissue_col_idx <- NULL
  for (i in 1:ncol(spots_labels_df)) {
    col_values <- as.character(spots_labels_df[, i])
    # Check if this column contains non-numeric tissue-like values
    unique_vals <- unique(col_values[!is.na(col_values)])
    if (length(unique_vals) < 20 && any(grepl("[a-zA-Z]", unique_vals))) {  # Contains letters and not too many unique values
      cat("DEBUG: Column", i, "(", colnames(spots_labels_df)[i], ") might contain tissue names:", paste(head(unique_vals, 5), collapse = ", "), "\n")
      tissue_col_idx <- i
      break
    }
  }
  
  # If we didn't find a tissue name column, create a mapping from numbers to tissue names
  if (is.null(tissue_col_idx)) {
    cat("DEBUG: No tissue name column found, creating mapping from label numbers\n")
    # Create a simple mapping - you may need to adjust this based on your data
    label_to_tissue <- c(
      "2" = "tissue_type_2", "3" = "tissue_type_3", "4" = "tissue_type_4", "5" = "tissue_type_5",
      "6" = "tissue_type_6", "7" = "tissue_type_7", "8" = "tissue_type_8", "9" = "tissue_type_9",
      "10" = "tissue_type_10", "11" = "tissue_type_11", "12" = "tissue_type_12", "13" = "tissue_type_13"
    )
    
    # Use the first column as spot names and the 'label' column (column 9) as labels
    spot_names_labels <- as.character(spots_labels_df[, 1])
    numeric_labels <- as.character(spots_labels_df$label)  # Use the 'label' column
    
    # Convert numeric labels to tissue names using mapping
    tissue_labels <- label_to_tissue[numeric_labels]
    tissue_labels[is.na(tissue_labels)] <- paste0("tissue_type_", numeric_labels[is.na(tissue_labels)])  # fallback for unmapped values
    names(tissue_labels) <- spot_names_labels
  } else {
    # Use the found tissue column
    spot_names_labels <- as.character(spots_labels_df[, 1])
    tissue_labels <- as.character(spots_labels_df[, tissue_col_idx])
    names(tissue_labels) <- spot_names_labels
  }
  
  cat("Sample labels:", paste(head(names(tissue_labels), 3), collapse = ", "), "\n")
  cat("Sample tissue types:", paste(head(tissue_labels, 3), collapse = ", "), "\n")
  
  # Convert Seurat spot names to match spots_labels format
  seurat_spot_names <- colnames(seurat_object)
  cat("Sample Seurat spot names:", paste(head(seurat_spot_names, 3), collapse = ", "), "\n")
  
  # Transform spot names: split by 'parquet.' and take the part after, then replace 'X' with 'spot'
  converted_spot_names <- sapply(seurat_spot_names, function(x) {
    if (grepl("parquet\\.", x)) {
      after_parquet <- sub(".*parquet\\.", "", x)
      converted <- gsub("X([0-9]+)x([0-9]+)", "spot\\1x\\2", after_parquet)
      return(converted)
    } else {
      return(x)  # Return original if no 'parquet.' found
    }
  })
  
  cat("Sample converted spot names:", paste(head(converted_spot_names, 3), collapse = ", "), "\n")
  
  # Map labels to converted spot names
  spot_labels <- tissue_labels[converted_spot_names]
  
  # Debug: Check what we got from the mapping
  cat("DEBUG: Class of spot_labels before conversion:", class(spot_labels), "\n")
  cat("DEBUG: First few spot_labels before conversion:", paste(head(spot_labels, 3), collapse = ", "), "\n")
  
  # Force conversion to character to avoid factor issues
  spot_labels <- as.character(spot_labels)
  
  # Set missing labels to "undetermined"
  spot_labels[is.na(spot_labels)] <- "undetermined"
  names(spot_labels) <- seurat_spot_names  # Use original Seurat spot names as names
  
  # Debug: Check final result
  cat("DEBUG: Class of spot_labels after conversion:", class(spot_labels), "\n")
  cat("DEBUG: First few spot_labels after conversion:", paste(head(spot_labels, 3), collapse = ", "), "\n")
  
  # Add labels to Seurat metadata
  seurat_object$label <- spot_labels
  
  # Summary of label assignment
  label_counts <- table(spot_labels)
  cat("Label assignment summary:\n")
  for (label in names(label_counts)) {
    cat("  ", label, ":", label_counts[label], "spots\n")
  }
  
  # Check how many spots were successfully matched
  matched_spots <- sum(spot_labels != "undetermined")
  total_spots <- length(spot_labels)
  cat("Successfully matched", matched_spots, "out of", total_spots, "spots (",
      round(100 * matched_spots / total_spots, 1), "%)\n")
  
} else {
  cat("No spots labels file found, all spots will be labeled as 'undetermined'\n")
  seurat_object$label <- "undetermined"
}

# Add patient IDs to metadata (extracted from spot names)
seurat_object$patient <- sapply(colnames(seurat_object), function(x) {
  # Remove everything up to the last "parquet." 
  after_parquet <- sub(".*parquet\\.", "", x)
  # Split by "_" and take the first part (TNBC1)
  strsplit(after_parquet, "_")[[1]][1]
})

# Save Seurat object
seurat_file <- file.path(output_dir, "TNBC_filtered_seurat_object.rds")
saveRDS(seurat_object, file = seurat_file)
cat("Seurat object saved before PCA to:", seurat_file, "\n")

# Run PCA
cat("Running PCA...\n")
tryCatch({
  seurat_object <- RunPCA(seurat_object, assay = "RNA")
}, error = function(e) {
  cat("Error in RunPCA:", e$message, "\n")
  cat("Trying PCA with fewer features...\n")
  # Use top 1000 variable features for PCA
  top_features <- head(VariableFeatures(seurat_object), 1000)
  seurat_object <<- RunPCA(seurat_object, assay = "RNA", features = top_features)
})

# Run UMAP
cat("Running UMAP...\n")
# Check how many PCs we actually have
n_pcs <- min(30, ncol(seurat_object@reductions$pca@cell.embeddings))
cat("Using", n_pcs, "PCs for UMAP\n")

tryCatch({
  seurat_object <- RunUMAP(seurat_object, assay = "RNA", dims = 1:n_pcs)
}, error = function(e) {
  cat("Error in RunUMAP:", e$message, "\n")
  cat("Trying UMAP with fewer dimensions...\n")
  n_pcs_reduced <- min(10, n_pcs)
  seurat_object <<- RunUMAP(seurat_object, assay = "RNA", dims = 1:n_pcs_reduced)
})



# Save results
cat("\n=== Saving results ===\n")

# Save Seurat object
seurat_file <- file.path(output_dir, "TNBC_filtered_seurat_object.rds")
saveRDS(seurat_object, file = seurat_file)
cat("Seurat object saved after PCA to:", seurat_file, "\n")

cat("\n=== Processing complete ===\n")
cat("Summary:\n")
if (!is.null(valid_spots)) {
  cat("- Spots in filter file:", length(valid_spots), "\n")
}
cat("- Final spots in Seurat object:", ncol(seurat_object), "\n")
cat("- Final genes in Seurat object:", nrow(seurat_object), "\n")
cat("- Patients:", length(unique(patient_ids)), "\n")
cat("- Min patient fraction:", min_patient_fraction, "\n")
cat("- Output directory:", output_dir, "\n")
if (!is.null(valid_spots) && ncol(seurat_object) < length(valid_spots)) {
  cat("\nNote: Final spot count is less than filter file because:\n")
  cat("  1. Some spots may have zero expression for all selected genes\n")
  cat("  2. Zero-count spots are removed during QC filtering\n")
}