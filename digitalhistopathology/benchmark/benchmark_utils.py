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

def get_optimal_cluster_number_one_model(files):

    # Read the first file to get the min distances

    with open(files[0]) as f:
        model_clusters = json.load(f)
        min_dists = list(model_clusters.keys())
        min_dists = [float(x) for x in min_dists if x != 'samples' and x != 'labels']
        min_dists.sort()

    df_all_clusters = []

    for min_dist in min_dists:

        dfs_clusters = []

        for file in files:
            k = int(file.split("_clusters")[0].split("_")[-1])
            with open(file) as f:
                model_clusters = json.load(f)
                
                model_clusters['patient_label'] = [x.split('_')[0] for x in model_clusters['samples']]
                # Convert labels to integers
                # unique_labels = sorted(set(model_clusters['patient_label']))
                # label_to_int = {label: idx for idx, label in enumerate(unique_labels)}
                # model_clusters['patient_label_int'] = [label_to_int[label] for label in model_clusters['patient_label']]
                
                df_model_clusters = pd.DataFrame(model_clusters[str(min_dist)]).T
                df_model_clusters = df_model_clusters.reset_index().rename(columns={"index": "n_neighbors"})
                
                df_model_clusters["n_clusters"] = k
                # df_model_clusters["ARI_patient"] = df_model_clusters['labels'].apply(lambda x: adjusted_rand_score(model_clusters['patient_label_int'], x))
                dfs_clusters.append(df_model_clusters)
                        
        dfs_clusters = pd.concat(dfs_clusters) 
        dfs_clusters["min_dist"] = min_dist
        df_all_clusters.append(dfs_clusters)

    df_all = pd.concat(df_all_clusters)
    df_all['silhouette_score'] = pd.to_numeric(df_all["silhouette_score"], errors='coerce')

    df_all = df_all.dropna(subset=["silhouette_score"]).reset_index(drop=True)

    idx = df_all.groupby(['n_clusters'])["silhouette_score"].idxmax()

    df_sil_ARI = df_all.loc[idx]

    df_sil_ARI['euclidian_dist_to_optimal'] = df_sil_ARI.apply(lambda row: np.sqrt((row['silhouette_score'] - 1)**2 + (row['ARI_patient'] - 0)**2), axis=1)

    # TODO: REMOVE AFTER THIS TEST
    # df_sil_ARI = df_sil_ARI[df_sil_ARI['n_clusters'] == 5]

    opti_clusters = df_sil_ARI.loc[df_sil_ARI['euclidian_dist_to_optimal'].idxmin()]

    return opti_clusters['n_clusters']

