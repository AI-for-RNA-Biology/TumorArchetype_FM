#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#

import numpy as np
from sklearn.cluster import KMeans, MeanShift, SpectralClustering, AgglomerativeClustering, Birch
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
from hdbscan import HDBSCAN
import scanpy as sc
import anndata as ad
from GraphST import utils
import pandas as pd
from scipy.cluster.hierarchy import dendrogram
import os
import matplotlib.pyplot as plt
import seaborn as sns

def clustering_scores(matrix, 
                        n_clusters, 
                        name="raw", 
                        algo='kmeans', 
                        min_cluster_size_hdbscan=None,
                        metric_hdbscan='euclidean',
                        assign_labels_spectral="cluster_qr",):
    
    
    """Run kmeans clustering and compute metrics from it.

    Args:
        matrix (np.array): Matrix on which kmeans is fit.
        n_clusters (int): Number of clusters.
        name (str, optional): Name of the data. Defaults to "raw".

    Returns:
        dict: Dictionnary with metrics and some other informations.
    """

    if algo == "kmeans":
        model = KMeans(n_clusters=n_clusters, random_state=42, n_init="auto")
    elif algo == "gaussian":
        model = GaussianMixture(n_components=n_clusters, random_state=42)
    elif algo == "meanshift":
        model = MeanShift(bandwidth=None)
    elif algo == "hdbscan":
        if min_cluster_size_hdbscan is None:
            min_cluster_size_hdbscan = int(0.05 * matrix.shape[0])
            print("Min cluster size = {}".format(min_cluster_size_hdbscan))
        model = HDBSCAN(min_cluster_size=min_cluster_size_hdbscan, metric=metric_hdbscan)
    elif algo == "spectral":
        model = SpectralClustering(n_clusters=n_clusters, assign_labels=assign_labels_spectral, random_state=42)
    elif algo == "birch":
        model = Birch(threshold=0.5, n_clusters=n_clusters)
    elif algo == "ward":
        model = AgglomerativeClustering(n_clusters=n_clusters, linkage="ward")
    elif algo == "leiden" or algo == "louvain":
        model = None
        an = ad.AnnData(matrix)
        an.obsm["cluster"] = matrix
        res = utils.search_res(an, n_clusters, use_rep="cluster", method=algo, start=0.1, end=3.0, increment=0.01)
        if algo == "leiden":
            sc.tl.leiden(an, random_state=0, resolution=res)
            labels = an.obs["leiden"]
        else:
            sc.tl.louvain(an, random_state=0, resolution=res)
            labels = an.obs["louvain"]
    else:
        raise Exception("This clustering algorithm is not implemented.")

    if model is not None:
        labels = model.fit_predict(matrix)
    

    # model = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    # model.fit(matrix)
    silhouette = silhouette_score(matrix, labels)
    # inertia = model.inertia_

    # Compute inertia
    centroids = np.array([matrix[labels == i].mean(axis=0) for i in range(labels.max() + 1)])
    inertia = 0
    for i in range(len(matrix)):
        cluster = labels[i]
        inertia += np.sum((matrix[i] - centroids[cluster]) ** 2)

    db = davies_bouldin_score(matrix, labels)
    cal = calinski_harabasz_score(matrix, labels)
    return {
        "data": name,
        "labels": labels.tolist(),
        "cluster_number": int(n_clusters),
        "silhouette_score": float(silhouette),
        "inertia": float(inertia),
        "davies_bouldin_score": float(db),
        "calinski_harabasz_score": float(cal),
    }

def plot_dendrogram(model, **kwargs):
    """Create linkage matrix and then plot the dendrogram
    From: https://scikit-learn.org/stable/auto_examples/cluster/plot_agglomerative_dendrogram.html

    Args:
        model: AgglomerativeClustering scikit-learn model.
    """
    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_, counts]).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)


def clustering_boxplot_per_patient(data, name, folder, ylabel="ari", saving_plot=True):
    """Boxplot of the per patient clustering ari score per algorithm.

    Args:
        data (Dictionary): Data dictionary containing the per patient ylabel metric per algorithm
        name (string): Title of the boxplot also used for saving the figure
        folder (string): Path to the saving folder for the boxplot.
        ylabel (string, optional): Metric used to evaluate the clustering. Default to ari
        saving_plot (boolean, optional): If the plot is saved to folder. Default to True.
    """
    data = {
        algo: {
            patient: {
                ylabel: v[ylabel],
                "number_of_categories": len(v["contingency_matrix"]),
            }
            for patient, v in d.items()
            if patient != "mean"
        }
        for algo, d in data.items()
    }
    transformed_data = {}
    for algo, patients in data.items():
        for patient, metrics in patients.items():
            num_categories = metrics["number_of_categories"]
            if len(data.items()) == 1:
                key = f"{num_categories}"
            else:
                key = f"{algo} {num_categories}"
            if key not in transformed_data:
                transformed_data[key] = {}
            transformed_data[key][patient] = metrics[ylabel]

    df = pd.DataFrame.from_records(transformed_data)
    df.boxplot(
        grid=False,
        ylabel=ylabel,
    )
    plt.xlabel("Number of categories")
    plt.title("Boxplot per patient of the {} score for the {} per algorithm".format(ylabel, name))

    colors = plt.get_cmap("tab20")(np.linspace(0, 1, 100))
    for idx, index in enumerate(df.index):
        plt.scatter(
            np.arange(1, len(list(df.columns)) + 1),
            df.loc[index, :],
            color=colors[idx],
            label=index,
            alpha=0.8,
        )
    plt.legend(title="Patient", bbox_to_anchor=(1.05, 1), loc="upper left")
    if saving_plot:
        plt.savefig(
            os.path.join(
                folder,
                "boxplot_per_patient_{}_score_{}_per_algo.pdf".format(ylabel, name.replace(" ", "_")),
            ),
            bbox_inches="tight",
        )
    plt.close()

def plot_ari_scores_all_patients(clustering_dict, model_list=None, stripplot=True, color_boxes=False, palette=None):

    if model_list is None:
        model_list = clustering_dict.keys()

    ari_scores = {}
    for model in model_list:
        ari_scores[model] = {}
        for patient in clustering_dict[model].keys():
            if (patient != 'all') and (patient != 'mean'):
                ari_scores[model][patient] = clustering_dict[model][patient]['ari']
    df_aris = pd.DataFrame.from_dict(ari_scores)
    df_aris_melted = pd.melt(df_aris, var_name='model', value_name='ari')
    df_aris_melted['patient'] = df_aris.index.to_list()*len(df_aris.columns)

    if color_boxes:
        sns.boxplot(data=df_aris_melted, x='model', y='ari', hue='model', linewidth=2, showfliers=False, palette=palette)
    else:
        sns.boxplot(data=df_aris_melted, x='model', y='ari', color='white', linewidth=2, showfliers=False)
    if stripplot:
        # sns.stripplot(data=df_aris_melted, x='model', y='ari', jitter=True, dodge=True, linewidth=1, hue='patient', palette='Accent')
        sns.stripplot(data=df_aris_melted, x='model', y='ari', jitter=True, dodge=True, linewidth=1, hue='patient', palette='Accent')
    plt.xticks(rotation=90)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    sns.despine()
    plt.title('ARI scores for unsupervised clustering', weight='bold')