#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import sys
import os 
sys.path.append(os.getcwd())
from digitalhistopathology.datasets.real_datasets import HER2Dataset 
import argparse




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Load molecular data and save gene embeddings.")
    parser.add_argument("--gene_embedding_saving_path", type=str, default="results/molecular", help="Path to save gene embeddings.")
    parser.add_argument("--patches_folder", type=str, default="results/HER2/compute_patches", help="Path to the folder containing patches.")
    args = parser.parse_args()

    gene_embedding_saving_path = args.gene_embedding_saving_path
    patches_folder = args.patches_folder

    if not os.path.exists(os.path.join(gene_embedding_saving_path, f"gene_embedding_HER2.csv")):

        print("Start loading gene embeddings")
        
        d = HER2Dataset(patches_folder=patches_folder, saving_emb_folder=gene_embedding_saving_path)
        ge = d.get_gene_embeddings(compute_emb=True)
        os.makedirs(gene_embedding_saving_path, exist_ok=True)
        ge.emb.to_df().to_csv(os.path.join(gene_embedding_saving_path, f"gene_embedding_HER2.csv"))

        # spots metadata
        spots_metadata = ge.emb.obs[['label']]
        spots_metadata['name_origin'] = [f"{idx[0]}{idx[5]}" for idx in spots_metadata.index]
        spots_metadata.to_csv(os.path.join(gene_embedding_saving_path, "spots_metadata.csv"))