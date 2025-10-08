import pandas as pd
import numpy as np
import multiprocessing as mp
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import os


def process_one_spot_fast(spot, gene_list, n_bootstrap=1000):
    # Filter valid genes
    genes_in_data = [gene for gene in gene_list if gene in spot.index]
    gene_count = len(genes_in_data)

    # If no genes matched, return NaN or another placeholder
    if gene_count == 0:
        return np.nan

    # Convert to NumPy for speed
    spot_values = spot.values
    spot_genes = np.array(spot.index)
    
    # Get indices for target genes
    gene_idx = np.isin(spot_genes, genes_in_data)
    observed_mean = spot_values[gene_idx].mean()

    # Generate bootstrap samples: (n_bootstrap × gene_count) indices
    sampled_indices = np.random.choice(len(spot), size=(n_bootstrap, gene_count), replace=True)

    # Gather values and compute means (vectorized)
    sampled_values = spot_values[sampled_indices]
    bootstrapped_means = sampled_values.mean(axis=1)

    # Compute two-sided p-value
    boot_mean = bootstrapped_means.mean()
    obs_diff = np.abs(observed_mean - boot_mean)
    extreme_count = np.sum(np.abs(bootstrapped_means - boot_mean) >= obs_diff)
    p_value = (extreme_count + 1) / (n_bootstrap + 1)

    # Return signed log10 p-value
    return -np.log10(p_value) if observed_mean > boot_mean else np.log10(p_value)

def compute_pathway_safe(pathway_name, gene_list, data_T, n_bootstrap=1000):
    np.random.seed()  # reseed RNG in each subprocess
    return pathway_name, data_T.apply(
        lambda spot: process_one_spot_fast(spot, gene_list, n_bootstrap), axis=1
    )

if __name__ == '__main__':
    pathways = pd.read_csv("/idiap/group/genomics/annotation/KEGG/human/KEGG_human_agg.csv", index_col=0)
    pathways['Genes'] = pathways['Gene_Symbols'].apply(lambda x: x.split(', '))
    pathways.drop(columns=['Gene_Symbols', 'Symbol'], inplace=True)

    data = pd.read_csv("/idiap/group/genomics/lfournier/digitalhistopathology/results/molecular/filtered_normalized_gene_expression.csv", index_col=0)

    # Transpose to spots x genes
    data_T = data.T  # shape: (n_spots, n_genes)
    spot_names = data_T.index
    pathway_names = pathways.index
    
    n_jobs = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
    print(f"Using {n_jobs} processes")


    pathway_args = [
        (name, pathways.loc[name, 'Genes'], data_T)
        for name in pathways.index
    ]

    results = []
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        futures = [executor.submit(compute_pathway_safe, *args) for args in pathway_args]
        for future in tqdm(as_completed(futures), total=len(futures)):
            results.append(future.result())

    # Reconstruct result DataFrame
    results_df = pd.concat([r[1].rename(r[0]) for r in results], axis=1)

    results_df.to_csv("/idiap/group/genomics/lfournier/digitalhistopathology/results/molecular/KEGG_pathways_HER2.csv")