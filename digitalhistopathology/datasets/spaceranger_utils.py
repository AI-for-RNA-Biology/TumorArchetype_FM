#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
"""
Utility functions for processing Space Ranger Visium HD data.

This module provides functions for loading and aggregating spatial transcriptomics data
from Space Ranger outputs, including filtering for protein-coding genes and aggregating
spots to different resolutions.
"""

import os
import scanpy as sc
import pandas as pd
import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix
from genomemanager.gtf_utils import read_gtf
import scipy.sparse as sp
    


def load_spaceranger_coding_genes(visium_output_dir, 
                                  gtf_path,
                                  resolution="002um", 
                                  filtered=True):
    """
    Load Space Ranger binned outputs and filter for protein-coding genes.
    
    Parameters:
    -----------
    visium_output_dir : str
        Path to the Space Ranger output directory (e.g., .../outs)
    resolution : str
        Resolution bin size (e.g., "002um", "008um", "016um")
    gtf_path : str
        Path to the GTF annotation file
        
    Returns:
    --------
    adata_coding : AnnData
        AnnData object filtered for protein-coding genes
    """
    # Construct path to binned outputs
    square_dir = os.path.join(visium_output_dir, "binned_outputs", f"square_{resolution}")
    positions = pd.read_parquet(os.path.join(square_dir, "spatial", "tissue_positions.parquet"))
    
    if not os.path.exists(square_dir):
        raise FileNotFoundError(f"Directory not found: {square_dir}")
    
    # Load the matrix with scanpy
    print(f"Loading {resolution} resolution matrix...")
    if filtered:
        adata = sc.read_10x_h5(os.path.join(square_dir, "filtered_feature_bc_matrix.h5"))
    else:
        adata = sc.read_10x_h5(os.path.join(square_dir, "raw_feature_bc_matrix.h5"))
    print(f"Loaded {resolution} matrix:")
    print(f"  Shape: {adata.shape} (spots × genes)")
    print(f"  Total UMIs: {adata.X.sum():,.0f}")
    

    
    # Load GTF annotation to get protein-coding genes
    if gtf_path is not None:
        print(f"\nLoading GTF annotation from: {gtf_path}")
        hg38 = read_gtf(gtf_path)
        hg38_coding = hg38[hg38["gene_type"] == "protein_coding"]
        hg38_coding_genes = list(hg38_coding['gene_name'].unique())
        print(f"Found {len(hg38_coding_genes)} protein-coding genes in GTF")
        
        # Fix duplicate gene names in adata
        print(f"\nChecking for duplicate gene names...")
        print(f"  Total genes: {len(adata.var_names)}")
        print(f"  Unique genes: {len(adata.var_names.unique())}")
        
        if len(adata.var_names) != len(adata.var_names.unique()):
            print("  Making gene names unique...")
            adata.var_names_make_unique()
            print("  Done!")
        
        # Filter to retain only protein-coding genes
        genes_to_keep = [g for g in adata.var_names if g in hg38_coding_genes]
        print(f"\nFound {len(genes_to_keep)} coding genes in adata (of {len(hg38_coding_genes)} total coding genes)")
        
        if len(genes_to_keep) == 0:
            raise RuntimeError("No coding genes were found in adata.var_names")
        
        # Subset to coding genes
        adata_coding = adata[:, genes_to_keep].copy()
        
        print(f"\nFinal filtered data:")
        print(f"  Shape: {adata_coding.shape}")
        print(f"  Total UMIs: {adata_coding.X.sum():,.0f}")
        
        adata_coding.obs = adata_coding.obs.merge(positions, left_index=True, right_on="barcode")
        print(f"Added spatial positions to adata.obs")
        
        return adata_coding

    else:
        print("No GTF path provided, returning unfiltered adata")
        adata.obs = adata.obs.merge(positions, left_index=True, right_on="barcode")

        return adata


def aggregate_spots_to_resolution(adata, current_resolution_um, target_resolution_um, verbose=True):
    """
    Aggregate spatial transcriptomics spots to a larger resolution.
    
    This function takes spots at a fine resolution and aggregates them into larger spots
    by summing expression counts. It uses grid-based binning with proper order preservation.
    Only large spots containing at least 95% of the expected number of small spots are retained.
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object with spatial coordinates in .obs
        Must contain: 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres'
    current_resolution_um : float
        Current resolution in micrometers (e.g., 2 for 2um spots)
    target_resolution_um : float
        Target resolution in micrometers (e.g., 100 for 100um spots)
    verbose : bool
        Whether to print progress information
        
    Returns:
    --------
    adata_aggregated : AnnData
        New AnnData object with aggregated spots. Only includes large spots that contain
        at least 95% of the expected number of small spots.
        
    Example:
    --------
    >>> # Aggregate 2um spots to 100um spots
    >>> large_adata = aggregate_spots_to_resolution(adata_002um, 
    ...                                              current_resolution_um=2, 
    ...                                              target_resolution_um=100)
    """
    # Calculate bin size (how many spots to aggregate)
    bin_size = int(target_resolution_um / current_resolution_um)
    
    if verbose:
        print(f"Aggregating from {current_resolution_um}um to {target_resolution_um}um resolution")
        print(f"Bin size: {bin_size}x{bin_size} spots")
    
    # Get array coordinates
    coords = adata.obs[['array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']].copy()
    
    # Create large spot grid by binning array coordinates
    coords['large_grid_row'] = (coords['array_row'] // bin_size).astype(int)
    coords['large_grid_col'] = (coords['array_col'] // bin_size).astype(int)
    
    # Create unique identifier for each large spot
    coords['large_spot_id'] = coords['large_grid_row'].astype(str) + '_' + coords['large_grid_col'].astype(str)

    spots_per_large = coords.groupby('large_spot_id').size()
    if verbose:
        print(f"\nOriginal number of small spots: {len(coords)}")
        print(f"Number of large spots: {coords['large_spot_id'].nunique()}")
        print(f"Average small spots per large spot: {len(coords) / coords['large_spot_id'].nunique():.2f}")
        
        # Show distribution of spots per large spot
        
        print(f"\nDistribution of small spots per large spot:")
        print(spots_per_large.value_counts().sort_index())

    ## Keep only the large spots that contains 95% of expected number of small spots
    expected_small_spots = bin_size**2

    if verbose:
        print(f"\nExpected small spots per large spot: {expected_small_spots}")
        print(f"Threshold (95%): {expected_small_spots*0.95:.1f}")
        print(f"Filtering to keep only large spots with > {expected_small_spots*0.95:.1f} small spots...")

    spots_per_large_filtered = spots_per_large[spots_per_large > expected_small_spots*0.95]
    
    if verbose:
        print(f"Large spots before filtering: {len(spots_per_large)}")
        print(f"Large spots after filtering: {len(spots_per_large_filtered)}")
        print(f"Removed {len(spots_per_large) - len(spots_per_large_filtered)} large spots")

    # Filter coords to only include spots that pass the threshold
    coords_filtered = coords[coords['large_spot_id'].isin(spots_per_large_filtered.index)].copy()
    
    if verbose:
        print(f"Small spots after filtering: {len(coords_filtered)}")

    expr_matrix = adata.X
    
    # Create mapping from original spots to large spots (using filtered coords)
    spot_to_large = coords_filtered['large_spot_id'].values
    
    # Create a mapping matrix: rows = large spots, columns = small spots
    unique_large_spots = coords_filtered['large_spot_id'].unique()
    large_spot_to_idx = {spot: idx for idx, spot in enumerate(unique_large_spots)}
    
    # Create index arrays for sparse matrix
    # We need to map the filtered coords back to their original indices in adata
    original_indices = coords_filtered.index.values
    large_spot_indices = np.array([large_spot_to_idx[spot] for spot in spot_to_large])
    
    data = np.ones(len(coords_filtered))
    aggregation_matrix = csr_matrix((data, (large_spot_indices, original_indices)), 
                                     shape=(len(unique_large_spots), len(adata)))
    
    # Aggregate expression: sum counts from all small spots in each large spot
    if verbose:
        print(f"\nShapes: aggregation matrix shape: {aggregation_matrix.shape}")
        print(f"Shapes: expression matrix shape: {expr_matrix.shape}")
    large_spot_expr = aggregation_matrix @ expr_matrix
    
    if verbose:
        print(f"\nOriginal expression matrix shape: {expr_matrix.shape}")
        print(f"Aggregated expression matrix shape: {large_spot_expr.shape}")
    
    # Create new AnnData object for large spots
    adata_aggregated = ad.AnnData(X=large_spot_expr, var=adata.var.copy())
    
    # Calculate center coordinates for large spots
    large_spot_coords = coords_filtered.groupby('large_spot_id').agg({
        'pxl_row_in_fullres': 'mean',
        'pxl_col_in_fullres': 'mean',
        'large_grid_row': 'first',
        'large_grid_col': 'first',
        'array_row': 'min',  # Keep track of original array coords
        'array_col': 'min'
    }).reset_index()
    
    # CRITICAL FIX: Reorder large_spot_coords to match the order of unique_large_spots used in aggregation
    large_spot_coords_ordered = large_spot_coords.set_index('large_spot_id').loc[unique_large_spots].reset_index()
    
    if verbose:
        print("\nOrder matches after reordering:", 
              (list(unique_large_spots) == large_spot_coords_ordered['large_spot_id'].tolist()))
    
    # Update adata_aggregated.obs with the correctly ordered coordinates
    adata_aggregated.obs = large_spot_coords_ordered.set_index('large_spot_id')
    adata_aggregated.obs['n_genes'] = (adata_aggregated.X > 0).sum(axis=1).A1
    adata_aggregated.obs['total_counts'] = adata_aggregated.X.sum(axis=1).A1
    
    if verbose:
        print(f"\nAggregated AnnData created:")
        print(f"  Shape: {adata_aggregated.shape}")
        print(f"  Total UMIs: {adata_aggregated.X.sum():,.0f}")
        print(f"  Mean UMIs per spot: {adata_aggregated.obs['total_counts'].mean():.1f}")
        print(f"  Mean genes per spot: {adata_aggregated.obs['n_genes'].mean():.1f}")
    
    return adata_aggregated



def assign_barcodes_to_patches_fast(df_patches, adata_small_spots, margin):
    """
    Fast assignment of barcodes to patches.
    Strategy:
    - sort spots by row once (O(n log n))
    - for each patch binary-search the row range (O(log n))
    - then filter only the small slice by columns (fast)
    Returns a dict mapping patch name -> numpy array of barcodes.
    """
    spot_rows = adata_small_spots.obs['pxl_row_in_fullres'].to_numpy()
    spot_cols = adata_small_spots.obs['pxl_col_in_fullres'].to_numpy()
    barcodes = adata_small_spots.obs['barcode'].to_numpy()

    # sort by row coordinate once
    order = np.argsort(spot_rows)
    rows_sorted = spot_rows[order]
    cols_sorted = spot_cols[order]
    bar_sorted = barcodes[order]

    result = {}

    for idx, patch in df_patches.iterrows():
        row_min = patch['start_height_origin'] + margin
        row_max = patch['start_height_origin'] + patch['shape_pixel'] - margin
        col_min = patch['start_width_origin'] + margin
        col_max = patch['start_width_origin'] + patch['shape_pixel'] - margin

        # locate row slice with binary search
        i0 = np.searchsorted(rows_sorted, row_min, side='left')
        i1 = np.searchsorted(rows_sorted, row_max, side='right')

        if i0 >= i1:
            result[patch['name']] = np.array([], dtype=bar_sorted.dtype)
            continue

        # filter the small slice by columns
        cols_slice = cols_sorted[i0:i1]
        mask = (cols_slice >= col_min) & (cols_slice <= col_max)
        result[patch['name']] = bar_sorted[i0:i1][mask]

    return result


def create_pseudobulk_for_patches_from_smaller_bins(adata_source, df_patches, save_path=None):
    """
    Create pseudobulk AnnData from spatial transcriptomics data by aggregating spots within patches.
    
    Parameters
    ----------
    adata_source : AnnData
        Source AnnData object with spatial transcriptomics data (spots x genes).
        Must have either 'barcode' column in obs or barcodes as obs_names.
    df_patches : pd.DataFrame
        DataFrame with patch information. Must contain:
        - 'name': patch names
        - 'barcodes': iterable of barcode lists for each patch
        Other columns will be preserved in the output obs.
    save_path : str, optional
        If provided, save the resulting AnnData to this path.
        
    Returns
    -------
    adata_pseudobulk : AnnData
        Pseudobulk AnnData object (patches x genes) with:
        - X: sparse matrix with summed expression per patch
        - obs: patch metadata from df_patches
        - var: gene metadata from adata_source
        - uns: metadata about pseudobulk creation
    """

    # 1) Build mapping barcode -> obs index in adata
    if 'barcode' in adata_source.obs.columns:
        source_barcodes = adata_source.obs['barcode'].astype(str).to_numpy()
    else:
        source_barcodes = adata_source.obs_names.astype(str)
    
    barcode_to_idx = {b: i for i, b in enumerate(source_barcodes)}
    
    n_genes = adata_source.n_vars
    dtype = getattr(adata_source.X, 'dtype', None) or np.float32
    
    rows = []          # will hold 1xG sparse rows (csr)
    valid_patch_names = []  # keep track of patch names aligned to rows
    
    for patch_name, barcodes in zip(df_patches['name'], df_patches['barcodes']):
        # barcodes may be numpy array / list / empty; ensure strings
        if barcodes is None or len(barcodes) == 0:
            # empty row
            rows.append(sp.csr_matrix((1, n_genes), dtype=dtype))
            valid_patch_names.append(patch_name)
            continue
        
        # map to indices (skip barcodes not found)
        idxs = [barcode_to_idx[b] for b in map(str, barcodes) if b in barcode_to_idx]
        if len(idxs) == 0:
            rows.append(sp.csr_matrix((1, n_genes), dtype=dtype))
            valid_patch_names.append(patch_name)
            continue
        
        # subset the AnnData X for these obs indices
        subX = adata_source.X[idxs]
        # compute pseudobulk sum over axis=0
        # handle sparse/dense
        if sp.issparse(subX):
            summed = subX.sum(axis=0)            # returns 1 x n_genes (sparse matrix or numpy.matrix)
            summed = sp.csr_matrix(summed)       # ensure csr 1xG
        else:
            # dense numpy array
            summed = np.asarray(subX).sum(axis=0, keepdims=True)
            summed = sp.csr_matrix(summed)
        # ensure 1 x n_genes shape
        if summed.shape[1] != n_genes:
            summed = summed.reshape(1, n_genes)
        rows.append(summed)
        valid_patch_names.append(patch_name)
    
    # 2) Stack rows to form (n_patches x n_genes) sparse matrix
    if len(rows) == 0:
        raise ValueError("No patches found in df_patches['barcodes'].")
    pseudobulk_X = sp.vstack(rows).tocsr()
    
    # 3) Build AnnData: obs <- df_patches (keep metadata), var <- adata_source.var
    obs = df_patches.copy()
    # set index to patch name for clarity
    obs = obs.set_index('name', drop=False).loc[valid_patch_names]
    
    # create AnnData
    adata_pseudobulk = ad.AnnData(X=pseudobulk_X, obs=obs, var=adata_source.var.copy())
    
    # 4) Store source info
    adata_pseudobulk.uns['pseudobulk_from'] = 'spatial_transcriptomics'
    adata_pseudobulk.uns['n_original_spots'] = int(adata_source.n_obs)
    adata_pseudobulk.uns['n_patches'] = len(valid_patch_names)
    
    # 5) Save to disk if desired
    if save_path is not None:
        adata_pseudobulk.write_h5ad(save_path)
        print(f"Saved pseudobulk AnnData to {save_path}")
    
    return adata_pseudobulk