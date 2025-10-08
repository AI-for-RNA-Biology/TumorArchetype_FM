# Warning: this script will need a lot of memory to run

# Import required libraries
import pandas as pd
import os
import numpy as np
import sys
sys.path.append("../")
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import glob
import scanpy as sc
from digitalhistopathology.datasets.real_datasets import HER2Dataset
import matplotlib
import json

# Configure matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

# Load spots metadata
spots_metadata = pd.read_csv("../results/molecular/spots_metadata.csv", index_col=0)

# Load gene expression data
gene_raw = pd.read_csv('../results/molecular/gene_embedding_HER2.csv', index_col=0)
emb_raw = ad.AnnData(gene_raw)

gene_filtered = pd.read_csv('../results/molecular/filtered_gene_expression.csv', index_col=0)
emb_filtered = ad.AnnData(gene_filtered.T)

gene_filtered_normalized = pd.read_csv('../results/molecular/filtered_normalized_gene_expression.csv', index_col=0)
emb_filtered_normalized = ad.AnnData(gene_filtered_normalized.T)

combat_filtered = pd.read_csv('../results/molecular/combat_corrected_filtered_counts.csv', index_col=0)
emb_combat = ad.AnnData(combat_filtered.T)

figures_folder = "../results/Figures"

# Process and visualize embeddings
for emb, name in zip([emb_raw, emb_filtered, emb_filtered_normalized, emb_combat], 
                     ["gene_embedding_HER2", "filtered_gene_expression", "filtered_normalized_gene_expression", "combat_corrected_filtered_counts"]):

    emb.obs.index = [idx.replace('-', '_') for idx in emb.obs.index]
    emb.obs['tumor'] = [x.split("_")[0] for x in emb.obs.index]
    emb.obs = emb.obs.merge(spots_metadata, left_index=True, right_index=True, how='left')
    emb.X = emb.X.astype(float)
    emb.layers['counts'] = emb.X.copy()

    # Perform dimensionality reduction
    sc.pp.highly_variable_genes(emb, n_top_genes=19000)
    sc.tl.pca(emb, use_highly_variable=True)
    sc.pp.neighbors(emb, n_pcs=10)
    sc.tl.umap(emb)

    # Plot UMAPs
    plt.figure(figsize=(20, 5))
    plt.subplot(1, 3, 1)
    sns.scatterplot(x=emb.obsm['X_umap'][:, 0], 
                    y=emb.obsm['X_umap'][:, 1], 
                    hue=emb.obs['tumor'], 
                    palette='Accent',
                    s=5,
                    alpha=0.7)
    sns.despine()
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.title("Labeled by tumor (all spots)", weight='bold')

    plt.subplot(1, 3, 2)
    sub_emb = emb[~emb.obs['label'].isna()]
    sns.scatterplot(x=sub_emb.obsm['X_umap'][:, 0],
                    y=sub_emb.obsm['X_umap'][:, 1],
                    hue=sub_emb.obs['tumor'],
                    palette="Accent",
                    s=10,
                    alpha=0.7)
    plt.title("Labeled by tumor (labeled spots)", weight='bold')

    plt.subplot(1, 3, 3)
    sns.scatterplot(x=emb.obsm['X_umap'][:, 0], 
                    y=emb.obsm['X_umap'][:, 1], 
                    hue=emb.obs['label'], 
                    palette=HER2Dataset.PALETTE,
                    s=10,
                    alpha=0.7)
    sns.despine()
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.title("Labeled by tissue type", weight='bold')
    plt.suptitle("UMAP of corrected gene expression", weight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_folder, "Figure4", "umap_combat_corrected_gene_expression.pdf"), bbox_inches='tight')

    # Plot KDE UMAPs
    plt.figure(figsize=(20, 5))
    plt.subplot(1, 3, 1)
    sns.kdeplot(x=emb.obsm['X_umap'][:, 0], 
                y=emb.obsm['X_umap'][:, 1], 
                hue=emb.obs['tumor'], 
                palette='Accent',
                alpha=0.7)
    sns.despine()
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.title("Labeled by tumor (all spots)", weight='bold')

    plt.subplot(1, 3, 2)
    sns.kdeplot(x=sub_emb.obsm['X_umap'][:, 0],
                y=sub_emb.obsm['X_umap'][:, 1],
                hue=sub_emb.obs['tumor'],
                palette="Accent",
                alpha=0.7)
    plt.title("Labeled by tumor (labeled spots)", weight='bold')

    plt.subplot(1, 3, 3)
    sns.kdeplot(x=emb.obsm['X_umap'][:, 0], 
                y=emb.obsm['X_umap'][:, 1], 
                hue=emb.obs['label'], 
                palette=HER2Dataset.PALETTE,
                alpha=0.7)
    sns.despine()
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.title("Labeled by tissue type", weight='bold')
    plt.suptitle("UMAP of corrected gene expression", weight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_folder, "Figure4", "umap_kde_combat_corrected_gene_expression.pdf"), bbox_inches='tight')

    # Save processed data
    emb.write_h5ad(f"../results/molecular/{name}.h5ad")
