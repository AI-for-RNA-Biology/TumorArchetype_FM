#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#

import pandas as pd
import seaborn as sns
import numpy as np

import json
import matplotlib.pyplot as plt

from sklearn.metrics import adjusted_rand_score
import glob 
import os


def plot_ari_scores_all_patients(clustering_dict):

    ari_scores = {}
    for model in clustering_dict.keys():
        ari_scores[model] = {}
        for patient in clustering_dict[model].keys():
            if (patient != 'all') and (patient != 'mean') and (patient != f'ARI_tumor'):
                ari_scores[model][patient] = clustering_dict[model][patient]['ari']
    df_aris = pd.DataFrame.from_dict(ari_scores)
    df_aris_melted = pd.melt(df_aris, var_name='model', value_name='ari')
    df_aris_melted['patient'] = list(df_aris.index)*len(df_aris.columns)

    sns.boxplot(data=df_aris_melted, x='model', y='ari', color='white', linewidth=2)
    sns.stripplot(data=df_aris_melted, x='model', y='ari', jitter=True, dodge=True, linewidth=1, hue='patient', palette='Accent')
    plt.xticks(rotation=90)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    sns.despine()
    plt.title('ARI scores for unsupervised clustering', weight='bold')


def barplot_correlation_with_eng_features(correlation_dict, title):
    plt.figure(figsize=(15, 15))
    i=1
    len(correlation_dict['stat'])/5
    for pc, pc_dict_stat in correlation_dict['stat'].items():
        absolute_dict = {k: abs(v) for k, v in pc_dict_stat.items()}
        plt.subplot(int(np.ceil(len(correlation_dict['stat'])/3)), 3, i)
        sns.barplot(x=list(absolute_dict.keys()), y=list(absolute_dict.values()), hue=list(absolute_dict.keys()))
        i+=1
        plt.xticks(rotation=90)
        plt.title(pc, weight='bold')
        plt.ylim(0,0.8)
        for j, model in enumerate(correlation_dict['p_value'][pc].keys()):
            plt.text(s=f"P={correlation_dict['p_value'][pc][model]:.2e}", y=0.75, x=j-0.35)
        sns.despine()

    plt.suptitle(title, weight='bold')
    plt.tight_layout()

# def get_optimal_cluster_number_one_model(files):

#     # Read the first file to get the min distances

#     with open(files[0]) as f:
#         model_clusters = json.load(f)
#         min_dists = list(model_clusters.keys())
#         min_dists = [float(x) for x in min_dists if x != 'samples' and x != 'labels']
#         min_dists.sort()

#     df_all_clusters = []

#     for min_dist in min_dists:

#         dfs_clusters = []

#         for file in files:
#             k = int(file.split("_clusters")[0].split("_")[-1])
#             with open(file) as f:
#                 model_clusters = json.load(f)
                
#                 model_clusters['patient_label'] = [x.split('_')[0] for x in model_clusters['samples']]
#                 # Convert labels to integers
#                 # unique_labels = sorted(set(model_clusters['patient_label']))
#                 # label_to_int = {label: idx for idx, label in enumerate(unique_labels)}
#                 # model_clusters['patient_label_int'] = [label_to_int[label] for label in model_clusters['patient_label']]
                
#                 df_model_clusters = pd.DataFrame(model_clusters[str(min_dist)]).T
#                 df_model_clusters = df_model_clusters.reset_index().rename(columns={"index": "n_neighbors"})
                
#                 df_model_clusters["n_clusters"] = k
#                 # df_model_clusters["ARI_patient"] = df_model_clusters['labels'].apply(lambda x: adjusted_rand_score(model_clusters['patient_label_int'], x))
#                 dfs_clusters.append(df_model_clusters)
                        
#         dfs_clusters = pd.concat(dfs_clusters) 
#         dfs_clusters["min_dist"] = min_dist
#         df_all_clusters.append(dfs_clusters)

#     df_all = pd.concat(df_all_clusters)
#     df_all['silhouette_score'] = pd.to_numeric(df_all["silhouette_score"], errors='coerce')

#     df_all = df_all.dropna(subset=["silhouette_score"]).reset_index(drop=True)

#     idx = df_all.groupby(['n_clusters'])["silhouette_score"].idxmax()

#     df_sil_ARI = df_all.loc[idx]

#     df_sil_ARI['euclidian_dist_to_optimal'] = df_sil_ARI.apply(lambda row: np.sqrt((row['silhouette_score'] - 1)**2 + (row['ARI_patient'] - 0)**2), axis=1)


#     opti_clusters = df_sil_ARI.loc[df_sil_ARI['euclidian_dist_to_optimal'].idxmin()]

#     return opti_clusters['n_clusters']

def select_best_umap_parameters_based_on_silhouette(saving_folder=None,
                                         files=None):
    
    """
    Select the best UMAP parameters based on silhouette score. Then, 
    select the best cluster number based on both silhouette score and ARI for patient batches.
    Algorithm first searches for the best UMAP parameters combination (min_dist, n_neighbors) for each 
    number of cluster that maximize the silhouette score. Then, for across best combinations, it searches 
    for the number of clusters that minimize the euclidean distance to the optimal point 
    (silhouette score = 1, ARI_patient = 0).
    Args:
        saving_folder (str): Folder where the JSON files are saved.
        files (list): List of JSON files to load. If None, will load from saving_folder.
    Returns:
        dict: Best UMAP parameters and corresponding scores."""

    if files is None:
        if saving_folder is not None:
            files = glob.glob(os.path.join(saving_folder, "scores_umap_across_parameters_*_clusters.json"))
        else: 
            raise Exception("No files to load")
    

    
    best_UMAP_params_per_file = {}
    
    for file in files:
        best_silhouette_score = 0
        best_n_neighbors = None
        best_min_dist = None

        n = int(file.split("_clusters")[0].split("_")[-1])
        with open(file) as f:
            results = json.load(f)
            
        
        for min_dist, n_neighbors_dict in results.items():
            if min_dist == 'samples':
                continue
            else:
                for n_neighbors, scores in n_neighbors_dict.items():
                    # euclidean_dist_to_optimal = np.sqrt((scores['silhouette_score'] - 1)**2 + (scores['ARI_patient'] - 0)**2)
                    if scores['silhouette_score'] > best_silhouette_score:
                        best_silhouette_score = scores['silhouette_score']
                        best_n_neighbors = int(n_neighbors)
                        best_min_dist = float(min_dist)
                        ari_patient = scores['ARI_patient']
                        # n_clusters = int(file.split("_clusters")[0].split("_")[-1])
                        labels = scores['labels']
        best_UMAP_params_per_file[n] = {'n_neighbors': best_n_neighbors, 
                                        'min_dist': best_min_dist,
                                        'silhouette_score': best_silhouette_score,
                                        'ari_patient': ari_patient,
                                        'labels': labels}
        
    df_best_params = pd.DataFrame.from_dict(best_UMAP_params_per_file, orient='index')
    df_best_params['euclidian_dist_to_optimal'] = df_best_params.apply(lambda row: np.sqrt((row['silhouette_score'] - 1)**2 + (row['ari_patient'] - 0)**2), axis=1)
    
    optimal_clusters = df_best_params.loc[df_best_params['euclidian_dist_to_optimal'].idxmin()]
    
    best_n_neighbors = optimal_clusters['n_neighbors']
    best_min_dist = optimal_clusters['min_dist']
    n_clusters = optimal_clusters.name
    best_score = optimal_clusters['euclidian_dist_to_optimal']
    silhouette_score = optimal_clusters['silhouette_score']
    ari_patient = optimal_clusters['ari_patient']
    labels = optimal_clusters['labels']
                    
    best_params = {'n_neighbors': best_n_neighbors, 
                    'min_dist': best_min_dist, 
                    'n_clusters': n_clusters, 
                    'euclidean_dist_to_optimal': best_score,
                    'silhouette_score': silhouette_score,
                    'ARI_patient': ari_patient,
                    'labels': labels,
                    'samples': results['samples']}
    
    return best_params

def select_best_umap_parameters_based_on_silhouette_and_batch_effect(saving_folder=None,
                                                                     files=None):
    
    """
    Select the best UMAP and kmeans parameters based on silhouette score and ARI for patient batches.
    Algorithm searches for the parameters combination (min_dist, n_neighbors, n_clusters) that minimize 
    the euclidean distance to the optimal point (silhouette score = 1, ARI_patient = 0).
    Args:
        saving_folder (str): Folder where the JSON files are saved.
        files (list): List of JSON files to load. If None, will load from saving_folder.
    Returns:
        dict: Best UMAP parameters and corresponding scores.
    """
    
    
    if files is None:
        if saving_folder is not None:
            files = glob.glob(os.path.join(saving_folder, "scores_umap_across_parameters_*_clusters.json"))
        else: 
            raise Exception("No files to load")
    
    best_score = 1000
    best_n_neighbors = None
    best_min_dist = None
    n_clusters = None
    
    for file in files:
        with open(file) as f:
            results = json.load(f)
            
        
        for min_dist, n_neighbors_dict in results.items():
            if min_dist == 'samples':
                continue
            else:
                for n_neighbors, scores in n_neighbors_dict.items():
                    euclidean_dist_to_optimal = np.sqrt((scores['silhouette_score'] - 1)**2 + (scores['ARI_patient'] - 0)**2)
                    if euclidean_dist_to_optimal < best_score:
                        best_score = euclidean_dist_to_optimal
                        best_n_neighbors = int(n_neighbors)
                        best_min_dist = float(min_dist)
                        n_clusters = int(file.split("_clusters")[0].split("_")[-1])
                        silhouette_score = scores['silhouette_score']
                        ari_patient = scores['ARI_patient']
                        labels = scores['labels']
                    
    best_params = {'n_neighbors': best_n_neighbors, 
                    'min_dist': best_min_dist, 
                    'n_clusters': n_clusters, 
                    'euclidean_dist_to_optimal': best_score,
                    'silhouette_score': silhouette_score,
                    'ARI_patient': ari_patient,
                    'labels': labels,
                    'samples': results['samples']}
    
    return best_params


def select_best_UMAP_and_kmeans_parameters(saving_folder=None,
                                           files=None, 
                                           strategy_silhouette_then_batch_effect=True):
    
    if strategy_silhouette_then_batch_effect:
        print("Selecting best UMAP parameters based on silhouette score first, then batch effect for clusters selection...")
        return select_best_umap_parameters_based_on_silhouette(saving_folder=saving_folder,
                                                               files=files)
    else:
        print("Selecting best UMAP parameters and clusters based on both silhouette score and batch effect...")
        return select_best_umap_parameters_based_on_silhouette_and_batch_effect(saving_folder=saving_folder,
                                                                                files=files)