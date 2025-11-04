#!/usr/bin/env python3
"""
Visualize gene expression UMAP for random patients.

This script:
1. Loads gene expression data from a parquet file
2. Loads labels from a CSV file
3. Computes UMAP on gene expression data
4. Randomly selects 10 patients
5. Plots UMAP colored by patient and by label
"""

import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

PALETTE = {
    "invasive cancer": "red",
    "cancer in situ": "orange",
    "immune infiltrate": "yellow",
    "breast glands": "green",
    "connective tissue": "blue",
    "adipose tissue": "cyan",
    "undetermined": "lightgrey",
}

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Visualize gene expression UMAP for random patients'
    )
    parser.add_argument(
        '--gene_expression',
        type=str,
        required=True,
        help='Path to parquet file containing gene expression data'
    )
    parser.add_argument(
        '--labels',
        type=str,
        required=True,
        help='Path to CSV file containing labels'
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        default='.',
        help='Directory to save output plots (default: current directory)'
    )
    parser.add_argument(
        '--n_patients',
        type=int,
        default=10,
        help='Number of random patients to select (default: 10)'
    )
    parser.add_argument(
        '--name',
        type=str,
        default='',
        help='Name for the output files (default: empty string)'
    )


    
    return parser.parse_args()


def format_spot_names(spot_names):
    """
    Format spot names from batch_01.parquet.TNBC1_X2x12 to TNBC1_spot2x12.
    
    Args:
        spot_names: Series or list of spot names
        
    Returns:
        Formatted spot names
    """
    formatted = pd.Series(spot_names).astype(str)
    # Remove batch_XX.parquet. prefix
    formatted = formatted.str.replace(r'^batch_\d+\.parquet\.', '', regex=True)
    # Replace _X with _spot
    formatted = formatted.str.replace('_X', '_spot')
    
    return formatted


def extract_patient_id(spot_name):
    """
    Extract patient ID from spot name.
    
    Args:
        spot_name: Spot name (e.g., 'TNBC1-spot2x12')
        
    Returns:
        Patient ID (e.g., 'TNBC1')
    """
    # Extract prefix before the first hyphen or underscore
    import re
    match = re.match(r'^([A-Za-z0-9]+)', spot_name)
    if match:
        return match.group(1)
    return spot_name


def load_gene_expression(filepath):
    """
    Load gene expression data from parquet file.
    
    Args:
        filepath: Path to parquet file
        
    Returns:
        DataFrame with gene expression data
    """
    print(f"Loading gene expression data from {filepath}...")
    df = pd.read_parquet(filepath)
    
    # Set spot_id as index if it exists
    if 'spot_id' in df.columns:
        df = df.set_index('spot_id')
        print(f"Set 'spot_id' as index")
    
    # Format spot names
    print("Formatting spot names...")
    df.index = format_spot_names(df.index)
    
    print(f"Loaded gene expression data: {df.shape[0]} spots x {df.shape[1]} genes")
    
    return df


def load_labels(filepath):
    """
    Load labels from CSV file.
    
    Args:
        filepath: Path to CSV file
        
    Returns:
        DataFrame with labels
    """
    print(f"\nLoading labels from {filepath}...")
    df = pd.read_csv(filepath, index_col=0)
    
    # Format spot names
    print("Formatting label index...")
    df.index = format_spot_names(df.index)
    
    print(f"Loaded labels for {df.shape[0]} spots")
    print(f"Label columns: {df.columns.tolist()}")
    
    return df


def compute_umap(gene_expr_df):
    """
    Compute UMAP on gene expression data.
    
    Args:
        gene_expr_df: DataFrame with gene expression data
        
    Returns:
        AnnData object with UMAP coordinates
    """
    print("\nComputing UMAP...")
    
    # Create AnnData object
    adata = sc.AnnData(gene_expr_df)
    
    print(f"Computing highly variable genes (top 19000)...")
    sc.pp.highly_variable_genes(adata, n_top_genes=19000)
    
    print(f"Computing PCA...")
    sc.tl.pca(adata, use_highly_variable=True)
    
    print(f"Computing neighbors...")
    sc.pp.neighbors(adata, n_pcs=10)
    
    print(f"Computing UMAP...")
    sc.tl.umap(adata)
    
    print("UMAP computation complete!")
    
    return adata


def select_random_patients(patient_ids, n_patients=10, seed=42):
    """
    Randomly select n patients.
    
    Args:
        patient_ids: Series of patient IDs
        n_patients: Number of patients to select
        seed: Random seed
        
    Returns:
        List of selected patient IDs
    """
    unique_patients = patient_ids.unique()
    n_available = len(unique_patients)
    
    print(f"\nFound {n_available} unique patients")
    
    if n_available < n_patients:
        print(f"Warning: Only {n_available} patients available, using all")
        return unique_patients.tolist()
    
    np.random.seed(seed)
    selected = np.random.choice(unique_patients, size=n_patients, replace=False)
    
    print(f"Randomly selected {n_patients} patients: {', '.join(selected)}")
    
    return selected.tolist()


def plot_umap(adata, labels_df, selected_patients, output_dir, name=""):
    """
    Plot UMAP colored by patient and by label.
    
    Args:
        adata: AnnData object with UMAP coordinates (already filtered to spots with labels)
        labels_df: DataFrame with labels
        selected_patients: List of selected patient IDs
        output_dir: Directory to save plots
    """
    # Extract patient IDs
    adata.obs['patient'] = [extract_patient_id(spot) for spot in adata.obs.index]
    
    # Add labels to adata
    label_column = 'label' if 'label' in labels_df.columns else labels_df.columns[0]
    adata.obs[label_column] = labels_df.loc[adata.obs.index, label_column]
    
    # Filter for selected patients
    mask = adata.obs['patient'].isin(selected_patients)
    adata_subset = adata[mask].copy()
    
    print(f"Subset contains {adata_subset.n_obs} labeled spots from {len(selected_patients)} patients")
    
    # Check for any remaining NaN labels (shouldn't happen now)
    has_label = ~adata_subset.obs[label_column].isna()
    n_with_labels = has_label.sum()
    if n_with_labels < adata_subset.n_obs:
        print(f"Warning: {adata_subset.n_obs - n_with_labels} spots still have missing labels")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Plot UMAP
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Plot 1: Colored by patient
    print("\nPlotting UMAP colored by patient...")
    umap_coords = adata_subset.obsm['X_umap']
    
    # Create color palette for patients
    patient_palette = dict(zip(
        selected_patients,
        sns.color_palette('tab10', len(selected_patients))
    ))
    
    for patient in selected_patients:
        patient_mask = adata_subset.obs['patient'] == patient
        axes[0].scatter(
            umap_coords[patient_mask, 0],
            umap_coords[patient_mask, 1],
            c=[patient_palette[patient]],
            label=patient,
            s=5,
            alpha=0.7
        )
    
    axes[0].set_xlabel('UMAP 1')
    axes[0].set_ylabel('UMAP 2')
    axes[0].set_title('UMAP colored by patient', fontweight='bold')
    axes[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left', markerscale=2)
    
    # Plot 2: Colored by label
    print("Plotting UMAP colored by label...")
    
    unique_labels = adata_subset.obs[label_column].dropna().unique()

    
    # Plot spots with labels
    for label in unique_labels:
        label_mask = adata_subset.obs[label_column] == label
        axes[1].scatter(
            umap_coords[label_mask, 0],
            umap_coords[label_mask, 1],
            c=[PALETTE[label]],
            label=label,
            s=5,
            alpha=0.7
        )
    
    axes[1].set_xlabel('UMAP 1')
    axes[1].set_ylabel('UMAP 2')
    axes[1].set_title('UMAP colored by label', fontweight='bold')
    axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', markerscale=2)
    
    plt.tight_layout()
    
    # Save figure
    output_file = output_path / f'gene_expression_umap_{name}.pdf'
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    print(f"\nSaved plot to {output_file}")
    
    plt.show()


def plot_umap_kde(adata, labels_df, output_dir, name=""):
    """
    Plot UMAP using KDE (density) plots colored by label for all patients.
    
    Args:
        adata: AnnData object with UMAP coordinates (already filtered to spots with labels)
        labels_df: DataFrame with labels
        output_dir: Directory to save plots
        name: Name suffix for output file
    """
    # Extract patient IDs
    adata.obs['patient'] = [extract_patient_id(spot) for spot in adata.obs.index]
    
    # Add labels to adata
    label_column = 'label' if 'label' in labels_df.columns else labels_df.columns[0]
    adata.obs[label_column] = labels_df.loc[adata.obs.index, label_column]
    
    print(f"\nPlotting KDE UMAP for all {adata.n_obs} spots")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Get UMAP coordinates
    umap_coords = adata.obsm['X_umap']
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    print("Plotting KDE UMAP colored by label...")
    
    # Get unique labels
    unique_labels = adata.obs[label_column].dropna().unique()
    
    # Plot KDE for each label
    for label in unique_labels:
        label_mask = adata.obs[label_column] == label
        if label_mask.sum() > 1:  # Need at least 2 points for KDE
            sns.kdeplot(
                x=umap_coords[label_mask, 0],
                y=umap_coords[label_mask, 1],
                color=PALETTE[label],
                label=label,
                fill=False,
                ax=ax
            )
    
    ax.set_xlabel('UMAP 1', fontsize=12)
    ax.set_ylabel('UMAP 2', fontsize=12)
    ax.set_title('UMAP KDE colored by label (all patients)', fontweight='bold', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    
    # Save figure
    output_file = output_path / f'gene_expression_umap_kde_{name}.pdf'
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    print(f"\nSaved KDE plot to {output_file}")
    
    plt.show()


def main():
    """Main function."""
    args = parse_args()
    
    # Load data
    gene_expr_df = load_gene_expression(args.gene_expression)
    labels_df = load_labels(args.labels)
    
    # Filter gene expression to keep only spots with labels
    common_spots = gene_expr_df.index.intersection(labels_df.index)

    print(f"\nFiltering to {len(common_spots)} spots with labels out of {gene_expr_df.shape[0]} total spots")
    gene_expr_df = gene_expr_df.loc[common_spots]
    
    # Compute UMAP on filtered data
    adata = compute_umap(gene_expr_df)
    
    # Extract patient IDs
    patient_ids = pd.Series([extract_patient_id(spot) for spot in adata.obs.index])
    
    # Plot KDE for all patients
    print("\n" + "="*50)
    print("Creating KDE plot for all patients")
    print("="*50)
    plot_umap_kde(adata, labels_df, args.output_dir, name=args.name)
    
    # Select random patients and create scatter plots
    for seed in range(10):
        print("\n" + "="*50)
        print(f"Creating scatter plot for seed {seed}")
        print("="*50)
        selected_patients = select_random_patients(patient_ids, args.n_patients, seed)

        # Plot UMAP
        plot_umap(adata, labels_df, selected_patients, args.output_dir, name=args.name + f"_seed{seed}")

    # Plot scatter for all patients
    print("Plotting UMAP for all patients")
    plot_umap_kde(adata, labels_df, args.output_dir, name=args.name + "_all_patients")
        
    print("\nDone!")


if __name__ == '__main__':
    main()
