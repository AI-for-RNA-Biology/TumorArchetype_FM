#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import sys
sys.path.append("../")  # Adjust the path to include the parent directory
import os 
from digitalhistopathology.datasets.real_datasets import HER2Dataset, TNBCDataset
import argparse




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Load molecular data and save gene embeddings.")
    parser.add_argument("--gene_embedding_saving_path", type=str, default="../results/molecular", help="Path to save gene embeddings.")
    parser.add_argument("--patches_folder", type=str, default="../results/compute_patches/her2_final_without_A", help="Path to the folder containing patches.")
    parser.add_argument("--dataset", type=str, default="HER2", help="Dataset name (default: HER2).")
    parser.add_argument("--gene_counts_path", type=str, default="../results/TNBC/count-matrices", help="Path to the gene counts file.")
    parser.add_argument("--dataDir", type=str, default="../data", help="Path to the data directory.")
    args = parser.parse_args()

    gene_embedding_saving_path = args.gene_embedding_saving_path
    patches_folder = args.patches_folder

    if not os.path.exists(os.path.join(gene_embedding_saving_path, f"gene_embedding_{args.dataset}.csv")):

        print("Start loading gene embeddings")
        
        if args.dataset == "HER2":
            d = HER2Dataset(patches_folder=patches_folder, saving_emb_folder=gene_embedding_saving_path)
        else:
            d = TNBCDataset(patches_folder=patches_folder, 
                            saving_emb_folder=gene_embedding_saving_path, 
                            gene_counts_path=args.gene_counts_path,
                            dataDir=args.dataDir)
            
        ge = d.get_gene_embeddings(compute_emb=True)
        ge.emb.to_df().to_csv(os.path.join(d.molecular_results_folder, f"gene_embedding_{args.dataset}.csv"))

        # spots metadata
        spots_metadata = ge.emb.obs[['label']]
        spots_metadata['name_origin'] = [f"{idx[0]}{idx[5]}" for idx in spots_metadata.index]
        spots_metadata.to_csv(os.path.join(d.molecular_results_folder, "spots_metadata.csv"))