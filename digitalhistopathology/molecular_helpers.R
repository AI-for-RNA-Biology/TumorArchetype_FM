#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#

library(tibble)
library(gprofiler2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)



retrieve_genes_from_GO_term <- function(go_term){
  result <- gost(query = go_term, organism = "hsapiens", sources = "GO:BP")
  gene_ensgs <- result$meta$genes_metadata$query$query_1$ensgs
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_ensgs, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  gene_names <- as.list(mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
  return(paste(gene_names, collapse = ";"))
}

plot_pathways <- function(cluster, results, gprofiler_results_true_labels = NULL) {
  ggplot(results, aes(x = reorder(term_name, -p_value), y = -log10(p_value))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = paste("Cluster", cluster), x = "Pathway", y = "-log10(p-value)") +
    theme_minimal()
}


get_pathway_scores_across_all_clusters <- function(res, res_true_labels = NULL) {
  # Get the union of all biological pathways
  
  if (!is.null(res_true_labels)) {
    gprofiler_results <- c(res$gprofiler_results, res_true_labels$gprofiler_results)
    markers <- bind_rows(res$markers, res_true_labels$markers)
  } else {
    gprofiler_results <- res$gprofiler_results  # Assuming a shallow copy is sufficient
    markers <- res$markers
  }
  
  
  all_pathways <- unique(unlist(lapply(gprofiler_results, function(df) df$term_name)))
  
  # Initialize an empty dataframe to store the -log10 p-values
  result_df <- data.frame(term_name = all_pathways)
  
  # Iterate over each cluster's dataframe in gprofiler_results
  for (cluster in names(gprofiler_results)) {
    # Get the current cluster's dataframe
    df <- gprofiler_results[[cluster]]
    
    # Calculate -log10 p-values for existing pathways
    pvalues_df <- df %>%
      mutate(log_p_value = -log10(p_value)) %>%
      select(term_name, log_p_value)
    
    # Identify missing pathways
    missing_pathways <- setdiff(all_pathways, df$term_name)
    
    # Recompute significance for missing pathways using g:Profiler
    if (length(missing_pathways) > 0) {
      missing_genes <- markers %>%
        dplyr::filter(cluster == !!cluster) %>%
        dplyr::pull(gene)
      
      gp_result <- gost(query = missing_genes, organism = "hsapiens", correction_method = "g_SCS", sources='GO:BP')
      
      if (length(gp_result) > 0){
        
        missing_pvalues_df <- gp_result$result %>%
          filter(term_name %in% missing_pathways) %>%
          mutate(log_p_value = -log10(p_value)) %>%
          select(term_name, log_p_value)
        
        # Combine existing and missing p-values
        pvalues_df <- bind_rows(pvalues_df, missing_pvalues_df)
      }
      
      # Merge with the result dataframe
      result_df <- result_df %>%
        left_join(pvalues_df, by = "term_name") %>%
        rename(!!cluster := log_p_value)
    }
  }
  
  # Replace NA values with 0
  result_df[is.na(result_df)] <- 0
  
  # Replace infinite values with 0
  # Replace infinite values with 0
  result_df[] <- lapply(result_df, function(x) {
    if (is.numeric(x)) {
      x[is.infinite(x)] <- 0
    }
    return(x)
  })  
  
  # Convert the dataframe to a matrix
  result_matrix <- as.matrix(result_df[,-1])  # Exclude the first column if it contains pathway names
  rownames(result_matrix) <- result_df$term_name
  
  return(result_matrix)
  
}

heatmap_pathways <- function(result_matrix, display_numbers=FALSE, directory_name="./", name="", method = 'pearson') {
  
  # Check for rows or columns with zero variance and remove them because correlation for distance matrix cannot be computed
  constant_cols <- apply(result_matrix, 2, sd) == 0
  constant_rows <- apply(result_matrix, 1, sd) == 0
  
  print(which(constant_cols))
  print(which(constant_rows))
  
  result_matrix <- result_matrix[, !constant_cols]
  result_matrix <- result_matrix[!constant_rows, ]
  
  # Calculate the distance matrix for columns using 1 - Spearman correlation
  
  if (method == 'pearson'){
    dist_matrix_cols <- as.dist(1 - cor(result_matrix, method = "pearson"))
    
    # Calculate the distance matrix for rows using 1 - pearson correlation
    dist_matrix_rows <- as.dist(1 - cor(t(result_matrix), method = "pearson"))
  } else if (method == 'spearman') {
    dist_matrix_cols <- as.dist(1 - cor(result_matrix, method = "pearson"))
    
    # Calculate the distance matrix for rows using 1 - Spearman correlation
    dist_matrix_rows <- as.dist(1 - cor(t(result_matrix), method = "spearman"))
  } else {
    dist_matrix_cols <- dist(t(result_matrix), method = "euclidean")
    dist_matrix_rows <- dist(result_matrix, method = "euclidean")
  }

  
  # Create a custom color palette
  custom_palette <- colorRampPalette(c("white", "black"))(100)
  custom_palette <- topo.colors(100)
  custom_palette <- colorRampPalette(c("gray", "deepskyblue", "magenta"))(100)
  custom_palette <- colorRampPalette(c('#FFFFFF','#69C4E0',  '#242B7A',  '#CF48DB', '#ffe800'))(100)
  # pdf(paste0(directory_name, "/", "heatmap_BPs", name, ".pdf"), width = 7, height = 10)
  
  # Set the background to transparent
  # par(bg = NA)
  # Plot the heatmap with the specified clustering method
  pheatmap(result_matrix, 
           cluster_rows = hclust(dist_matrix_rows, method = "complete"), 
           cluster_cols = hclust(dist_matrix_cols, method = "complete"), 
           display_numbers = display_numbers, 
           fontsize_number = 6, 
           main = "-log10(p-values) Heatmap",
           color = custom_palette,
           number_color = "black", 
           filename = paste0(directory_name, "/", "heatmap_BPs", name, ".pdf"),
           height=10,
           width = 7,
           fontsize_row=9)
  
  # dev.off()
  

}


format_seurat_with_predicted_csv <- function(seurat_object, path_to_predicted_clusters){
  predicted_clusters_uni_metadata <- read.csv(path_to_predicted_clusters, row.names = 1)
  predicted_clusters_uni_metadata <- predicted_clusters_uni_metadata %>% select(-label)
  rownames(predicted_clusters_uni_metadata) <- gsub("_", "-", rownames(predicted_clusters_uni_metadata))
  
  predicted_clusters_uni_metadata <- predicted_clusters_uni_metadata %>% rownames_to_column(var = "rowname")
  
  
  meta_data <- seurat_object@meta.data %>%
    rownames_to_column(var = "rowname")
  
  merged_meta_data <- merge(predicted_clusters_uni_metadata, meta_data, by = "rowname")
  rownames(merged_meta_data) <- merged_meta_data$rowname
  merged_meta_data <- merged_meta_data %>% select(-rowname)
  
  seurat_object_predicted <- seurat_object[ , rownames(merged_meta_data)]
  seurat_object_predicted@meta.data <- merged_meta_data
  # Set the identities in the Seurat object to your custom labels
  Idents(seurat_object_predicted) <- seurat_object_predicted@meta.data$predicted_label
  
  return(seurat_object_predicted)
  
}

get_clusters_DGE_BPs <- function(seurat_object, upregulated = TRUE){

  
  # Use FindAllMarkers to find markers for each cluster
  if (upregulated == TRUE){
    markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, latent.vars = "patient", test.use='MAST')
  } else {
    markers <- FindAllMarkers(seurat_object, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, latent.vars = "patient", test.use='MAST')
  }
  
  
  # Print the top markers
  print(head(markers))

  # Initialize a list to store the results
  gprofiler_results <- list()
  
  # Loop through each cluster and run g:Profiler
  for (cluster in unique(markers$cluster)) {
    # Get the markers for the current cluster
    print(cluster)
    
    if (upregulated == TRUE){
      cluster_genes <- markers %>%
        dplyr::filter(cluster == !!cluster) %>% dplyr::filter(avg_log2FC > 0) %>% dplyr::pull(gene)
    } else {
      cluster_genes <- markers %>%
        dplyr::filter(cluster == !!cluster) %>% dplyr::filter(avg_log2FC < 0) %>% dplyr::pull(gene)
    }
    
    # Run g:Profiler for the current cluster
    gp_result <- gost(query = cluster_genes, organism = "hsapiens", correction_method = "g_SCS", sources='GO:BP')
    
    # Select top 10 pathways based on p-value 
    if (length(gp_result) > 0) {
      specific_terms <- gp_result$result %>% dplyr::filter(term_size < 300) %>% dplyr::filter(intersection_size > 5)  %>% arrange(p_value) %>% slice_head(n = 10) 
      
      gprofiler_results[[cluster]] <- specific_terms
      
      # Store the results
    }
  }
  
  return(list(markers = markers, gprofiler_results = gprofiler_results))
}


plot_pathways <- function(cluster, results) {
  ggplot(results, aes(x = reorder(term_name, -p_value), y = -log10(p_value))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = paste("Cluster", cluster), x = "Pathway", y = "-log10(p-value)") +
    theme_minimal()
}


save_dge_pathways_analysis_per_clusters <- function(res, directory_name, add_name = ""){
  write.csv(res$markers, file = paste0(directory_name, "/", "marker_genes.csv"))
  
  res$gprofiler_results <- lapply(res$gprofiler_results, function(df) {
    if (length(df) > 0) {
      df$GO_genes <- sapply(df$term_id, retrieve_genes_from_GO_term)
    }
    return(df)
  })
  
  res$gprofiler_results <- lapply(res$gprofiler_results, function(df) {
    df$parents <- sapply(df$parents, function(x) paste(x, collapse = ";"))
    return(df)
  })
  
  
  lapply(names(res$gprofiler_results), function(name) {
    write.csv(res$gprofiler_results[[name]], file = paste0(directory_name, "/", "pathways_results_cluster_", name, add_name, ".csv"), row.names = FALSE)
  })
}
