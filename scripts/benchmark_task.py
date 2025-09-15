#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import argparse
import pandas as pd
import sys
import seaborn as sns
import json
import numpy as np
import os
sys.path.append("../")

import glob

from digitalhistopathology.benchmark.benchmark_shannon import BenchmarkShannon
from digitalhistopathology.benchmark.benchmark_clustering import BenchmarkClustering
from digitalhistopathology.benchmark.benchmark_regression import BenchmarkRegression
from digitalhistopathology.benchmark.benchmark_invasive_cancer import BenchmarkInvasive


def benchmark(task,
              path_to_pipeline=None, 
              pipelines_list=None, 
              dataset_name=None,
              results_folder="../results",
              emb_df_csv_path=None,
              saving_folder=None,
              image_embedding_name="image_embedding.h5ad",
              engineered_features_saving_folder=None,
              engineered_features_type='scMTOP',
              extension='png',
              group='tumor',
              min_comp=512,
              algo='kmeans',
              pct_variance=0.9,
              molecular_emb_path=None,
              molecular_name='gene',
              ref_model_emb=None,
              regression_type="linear",
              n_splits=5,
              alpha_reg=1,
              on_invasive=False,
              min_cluster=4,
              max_cluster=10,
              cluster_step=1,
              dataset="HER2"):
        
    if task == 'shannon_entropy':

        benchmark_obj = BenchmarkShannon(path_to_pipeline=path_to_pipeline, 
                         pipelines_list=pipelines_list, 
                         dataset_name=dataset_name,
                         results_folder=results_folder,
                         emb_df_csv_path=emb_df_csv_path,
                         saving_folder=saving_folder,
                         image_embedding_name=image_embedding_name,
                         engineered_features_saving_folder=engineered_features_saving_folder,
                         engineered_features_type=engineered_features_type,
                         extension=extension,
                         group=group,
                         min_comp=min_comp,
                         pct_variance=pct_variance,
                         dataset=dataset
                         )
        
    elif task == 'unsupervised_clustering_ARI':

        benchmark_obj = BenchmarkClustering(path_to_pipeline=path_to_pipeline, 
                         pipelines_list=pipelines_list, 
                         dataset_name=dataset_name,
                         results_folder=results_folder,
                         emb_df_csv_path=emb_df_csv_path,
                         saving_folder=saving_folder,
                         image_embedding_name=image_embedding_name,
                         engineered_features_saving_folder=engineered_features_saving_folder,
                         engineered_features_type=engineered_features_type,
                         extension=extension,
                         group=group,
                         algo=algo,
                         dataset=dataset)
        
    elif task == 'regression':

        benchmark_obj = BenchmarkRegression(path_to_pipeline=path_to_pipeline,
                         pipelines_list=pipelines_list, 
                         dataset_name=dataset_name,
                         results_folder=results_folder,
                         emb_df_csv_path=emb_df_csv_path,
                         saving_folder=saving_folder,
                         image_embedding_name=image_embedding_name,
                         engineered_features_saving_folder=engineered_features_saving_folder,
                         engineered_features_type=engineered_features_type,
                         extension=extension,
                         group=group,
                         regression_type=regression_type,
                         n_splits=n_splits,
                         alpha_reg=alpha_reg,
                         on_invasive=on_invasive,
                         dataset=dataset)
        
    elif task == 'invasive_cancer_clustering':

        benchmark_obj = BenchmarkInvasive(path_to_pipeline=path_to_pipeline,
                         pipelines_list=pipelines_list, 
                         dataset_name=dataset_name,
                         results_folder=results_folder,
                         emb_df_csv_path=emb_df_csv_path,
                         saving_folder=saving_folder,
                         image_embedding_name=image_embedding_name,
                         engineered_features_saving_folder=engineered_features_saving_folder,
                         engineered_features_type=engineered_features_type,
                         extension=extension,
                         group=group,
                         molecular_emb_path=molecular_emb_path,
                         molecular_name=molecular_name,
                         ref_model_emb=ref_model_emb,
                         algo=algo,
                         min_cluster=min_cluster,
                         max_cluster=max_cluster,
                         cluster_step=cluster_step,
                         dataset=dataset)
    else:
        raise ValueError("Task not recognized. Please choose between 'shannon_entropy', 'unsupervised_clustering_ARI', 'regression', and 'invasive_cancer_clustering'.")

    benchmark_obj.saving_folder = os.path.join(benchmark_obj.saving_folder, task)
    benchmark_obj.compute_image_embeddings()
    benchmark_obj.load_engineered_features()

    benchmark_obj.execute_pipeline()



def main():
    parser = argparse.ArgumentParser(description="Benchmarking script for digital histopathology pipelines")
    parser.add_argument('--config', type=str, default=None, help='Path to the config file')
    parser.add_argument('--benchmark_task', type=str, default='regression', help='Type of benchmark task to perform')
    parser.add_argument('--path_to_pipeline', type=str, nargs='+', required=True, help='Path to the pipeline directory')
    parser.add_argument('--pipelines_list', type=str, nargs='+', required=True, help='List of pipeline models to include')
    parser.add_argument('--saving_folder', type=str, default="../results/benchmark", help='Folder to save benchmark results')
    parser.add_argument('--emb_df_csv_path', type=str, required=False, default=None, help='Path to the csv file containing the engineered features')
    parser.add_argument('--dataset_name', type=str, default='her2_final_without_A', help='Dataset name')
    parser.add_argument('--engineered_features_saving_folder', type=str, default=None, help='Folder where engineered features are stored')
    parser.add_argument('--extension', type=str, default='png', help='Extension of the image to save')
    parser.add_argument('--image_embedding_name', type=str, default='image_embedding.h5ad', help='Name of the image embedding file')
    parser.add_argument('--results_folder', type=str, default="../results", help='Base of the results folder')
    parser.add_argument('--dataset', type=str, default="HER2", help='Dataset to use for the benchmark')
    
    # For the 2 clustering tasks
    parser.add_argument('--clustering_algo', type=str, default='kmeans', help='Clustering algorithm to use')

    # Specific to Shannon entropy
    parser.add_argument('--pct_variance', type=float, default=0.9, help='Percentage of variance to explain')
    parser.add_argument('--min_comp', type=int, default=512, help='Minimum number of components to use')

    # Specific to regression of engineered features
    parser.add_argument('--regression_type', type=str, default='linear', help='Type of regression to perform')
    parser.add_argument('--n_splits', type=int, default=5, help='Number of splits for cross validation')
    parser.add_argument('--alpha_reg', type=float, default=1, help='Alpha parameter (regularization) for regression')
    parser.add_argument('--regress_only_invasive', action="store_true", help='Regress only on invasive cancer')

    ## Specific to invasive cancer clustering
    parser.add_argument('--molecular_emb_path', type=str, default='../results/molecular/combat_corrected_embeddings.h5ad', help='Molecular embeddings to use')
    parser.add_argument('--molecular_name', type=str, default='gene', help='Name of the molecular feature')
    parser.add_argument('--ref_model_emb', type=str, default='../pipeline/uni/image_embeddings.h5ad', help='Reference image embeddings to use')
    parser.add_argument('--min_cluster', type=int, default=4, help='Minimum number of clusters to use')
    parser.add_argument('--max_cluster', type=int, default=10, help='Maximum number of clusters to use')
    parser.add_argument('--cluster_step', type=int, default=1, help='Step size for clustering')
    

    ## Plot colors
    parser.add_argument('--shades_palette', action="store_true", help='Use shades palette')
    parser.add_argument('--color', type=str, default='Blues', help='Verbose mode')

    args = parser.parse_args()

    # Set your colorpalette
    if args.shades_palette:
        a = sns.color_palette(args.color, n_colors=10)
    else:
        a = sns.color_palette("Set3")[0:1] + ['powderblue'] + ['lightpink'] + sns.color_palette("Set3")[1:4] + [
            'lightgray', 'mediumspringgreen', 'cadetblue', 'teal', 'aquamarine', 'aqua', 'orangered', 
            'firebrick', 'maroon', 'sienna', 'deeppink'
        ]

    # If provided, load the config file
    if args.config is not None:
        with open(args.config, 'r') as f:
            config = json.load(f)
            args.benchmark_task = config.get('benchmark_task', 'shannon_entropy')
            args.path_to_pipeline = config.get('path_to_pipeline', [])
            args.pipelines_list = config.get('pipelines_list', [])
            args.dataset_name = config.get('dataset_name', 'her2_final_without_A')
            args.results_folder = config.get('results_folder', "../results")
            args.emb_df_csv_path = config.get('emb_df_csv_path', None)
            args.saving_folder = config.get('saving_folder', "../results/benchmark")
            args.image_embedding_name = config.get('image_embedding_name', 'image_embedding.h5ad')
            args.engineered_features_saving_folder = config.get('engineered_features_saving_folder', None)
            args.extension = config.get('extension', 'png')
            args.clustering_algo = config.get('clustering_algo', 'kmeans')
            args.pct_variance = config.get('pct_variance', 0.9)
            args.min_comp = config.get('min_comp', 512)
            args.molecular_emb_path = config.get('molecular_emb_path', '../results/molecular/combat_corrected_embeddings.h5ad')
            args.molecular_name = config.get('molecular_name', 'gene')
            args.ref_model_emb = config.get('ref_model_emb', '../pipeline/uni/image_embeddings.h5ad')
            args.min_cluster = config.get('min_cluster', 4)
            args.max_cluster = config.get('max_cluster', 10)
            args.cluster_step = config.get('cluster_step', 1)
            args.dataset = config.get('dataset', 'HER2')
            
    sns.set_palette(a)

    if args.benchmark_task == "all":
        # Run all benchmark tasks
        for task in ['shannon_entropy', 'unsupervised_clustering_ARI', 'regression', 'invasive_cancer_clustering']:
            benchmark(
                task=task,
                path_to_pipeline=args.path_to_pipeline,
                pipelines_list=args.pipelines_list,
                dataset_name=args.dataset_name,
                results_folder=args.results_folder,
                emb_df_csv_path=args.emb_df_csv_path,
                saving_folder=args.saving_folder,
                image_embedding_name=args.image_embedding_name,
                engineered_features_saving_folder=args.engineered_features_saving_folder,
                engineered_features_type='scMTOP',
                extension=args.extension,
                group='tumor',
                min_comp=args.min_comp,
                algo=args.clustering_algo,
                pct_variance=args.pct_variance,
                molecular_emb_path=args.molecular_emb_path,
                molecular_name=args.molecular_name,
                ref_model_emb=args.ref_model_emb,
                regression_type=args.regression_type,
                n_splits=args.n_splits,
                alpha_reg=args.alpha_reg,
                on_invasive=args.regress_only_invasive,
                min_cluster=args.min_cluster,
                max_cluster=args.max_cluster,
                cluster_step=args.cluster_step,
                dataset=args.dataset,
            )
    else:
        # Run the specified benchmark task
        benchmark(
            task=args.benchmark_task,
            path_to_pipeline=args.path_to_pipeline,
            pipelines_list=args.pipelines_list,
            dataset_name=args.dataset_name,
            results_folder=args.results_folder,
            emb_df_csv_path=args.emb_df_csv_path,
            saving_folder=args.saving_folder,
            image_embedding_name=args.image_embedding_name,
            engineered_features_saving_folder=args.engineered_features_saving_folder,
            engineered_features_type='scMTOP',
            extension=args.extension,
            group='tumor',
            min_comp=args.min_comp,
            algo=args.clustering_algo,
            pct_variance=args.pct_variance,
            molecular_emb_path=args.molecular_emb_path,
            molecular_name=args.molecular_name,
            ref_model_emb=args.ref_model_emb,
            regression_type=args.regression_type,
            n_splits=args.n_splits,
            alpha_reg=args.alpha_reg,
            on_invasive=args.regress_only_invasive,
            min_cluster=args.min_cluster,
            max_cluster=args.max_cluster,
            cluster_step=args.cluster_step,
            dataset=args.dataset,
        )

if __name__ == "__main__":
    main()

