#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
from digitalhistopathology.dimred.dimensionality_reduction import DimRed
from GraphST import utils
from digitalhistopathology.clustering.clustering_utils import plot_dendrogram
import numpy as np
from sklearn.cluster import KMeans, MeanShift, SpectralClustering, Birch, AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from hdbscan import HDBSCAN
import scanpy as sc
import os
import json
import glob
import warnings
from matplotlib import cm
from scipy.cluster.hierarchy import dendrogram
from concurrent.futures import ThreadPoolExecutor, as_completed
from kneed import KneeLocator
import pandas as pd
import anndata as ad

from sklearn.metrics import (
    adjusted_rand_score,
    adjusted_mutual_info_score,
    fowlkes_mallows_score,
)

from sklearn.metrics.cluster import contingency_matrix

from digitalhistopathology.clustering.clustering_utils import clustering_scores

import matplotlib.pyplot as plt

class FindClusters(DimRed):
    """Class to compute clustering on the embedding matrix.
    It inherits from EmbDimensionalityReduction class.
    """

    def __init__(self, emb=None, saving_plots=False, result_saving_folder=None):
        DimRed.__init__(self, emb=emb, saving_plots=saving_plots, result_saving_folder=result_saving_folder)
        self._unsupervised_clustering_score_files = []
        self.info_optimal_number_of_clusters = dict()

    def set_unsupervised_clustering_score_files(self, files_list):
        self._unsupervised_clustering_score_files = files_list


    def unsupervised_clustering(
        self,
        n_clusters=6,
        true_label_col="label",
        algo="kmeans",
        layer=None,
        obsm=None,
        var_ratio_threshold_for_svd=0.9,
        u_comp_list=None,
        denoised_for_svd=False,
        metric_hdbscan="euclidean",
        min_cluster_size_hdbscan=None,
        assign_labels_spectral="cluster_qr",
        center_before_svd=False,
        scale_before_svd=False,
        multiply_by_variance=False,
    ):
        """Compute ARI between the ground truth labels and the labels predicted from the chosen clustering algorithm.

        Args:
            n_clusters (int, optional): Number of clusters to consider for the clustering algorithm. Defaults to 6.
            true_label_col (str, optional): Column name in emb.obs  where the ground truth labels are. Defaults to "label".
            algo (str, optional): Clustering algorithm to use. Defaults to "kmeans".
            layer (str, optional): Layer name in emb.layers on which to run the clustering algorithm. If None and obsm is None, the raw emb is taken. If "svd", u matrix from svd analysis is used. Defaults to None.
            obsm (str, optional): Layer name in emb.obsm on which to run the clustering algorithm. If None and layer is None, the raw emb is taken. Defaults to None.
            var_ratio_threshold_for_svd (float, optional): If you use layer = "svd", it is the threshold of cumulative variance explained ratio until which we keep the ui. Between 0 and 1. Defaults to 0.9.
            u_comp_list (list, optional): If layer = "svd", it is the list of the ui you want to keep for the matrix clustering. If None, the var_ratio_threshold_for_svd is used. Defaults None.
            denoised_for_svd (bool, optional): If you use layer = "svd", it chooses between the svd computed on the raw or the denoised matrix. Defaults to False.
            metric_hdbscan (str, optional): HDBSCAN metric parameter. Defaults to "euclidean".
            min_cluster_size_hdbscan (int, optional): HDBSCAN min_cluster_size parameter. If None, the 5% of the number of samples is taken as min_cluster_size. Default to None.
            assign_labels_spectral (str, optional): SpectralClustering assign_labels parameter. "kmeans", "discretize", "cluster_qr". Defaults to "discretize".
        Raises:
            Exception: The chosen algorithm is not available in this method.

        Returns:
            dict: Dictionary of metrics: ari, ami, fmi
        """

        if layer is None and obsm is None:
            matrix = self.emb.X.copy()
        elif layer == "svd":
            if self.svd["denoised"] is None or self.svd["denoised"] != denoised_for_svd:
                self.compute_svd(denoised=denoised_for_svd, center=center_before_svd, scale=scale_before_svd)

            # Equivalent to using PCs if we multiply by the variance
            if multiply_by_variance:
                U_df = self.svd["U_df"] * self.svd["S"]
            else:
                U_df = self.svd["U_df"].copy()

            if u_comp_list is None:
                if var_ratio_threshold_for_svd is not None:
                    cumulative_variance = np.cumsum(self.get_explained_variance_ratio_list())
                    threshold_u_number = np.argwhere(cumulative_variance > var_ratio_threshold_for_svd)[0][0]
                    matrix = np.array(U_df.iloc[:, 0:threshold_u_number].copy())
                    print(
                        "Cut u matrix until {}% explained variance, u shape = {}".format(
                            round(var_ratio_threshold_for_svd * 100, 2), matrix.shape
                        )
                    )
                else:
                    matrix = np.array(U_df)
                    print("Using the full SVD matrix", flush=True)
            else:
                matrix = np.array(U_df.loc[:, u_comp_list].copy())
        elif layer is not None and obsm is None:
            if (layer == 'umap') and ('umap' not in self.emb.obsm.keys()):
                print(f"Computing UMAP for unsupervised clustering", flush=True)
                self.compute_umap()
            if layer in self.emb.layers.keys():
                matrix = self.emb.layers[layer].copy()
            elif layer in self.emb.obsm.keys():
                matrix = self.emb.obsm[layer].copy()
            else:
                raise Exception("Layer {} not found in emb.layers or emb.obsm".format(layer))
        elif obsm is not None and layer is None:
            if (obsm == 'umap') and ('umap' not in self.emb.obsm.keys()):
                print(f"Computing UMAP for unsupervised clustering", flush=True)
                self.compute_umap()
            if obsm in self.emb.obsm.keys():
                matrix = self.emb.obsm[obsm].copy()
            elif obsm in self.emb.layers.keys():
                matrix = self.emb.layers[obsm].copy()
            else:
                raise Exception("Obsm {} not found in emb.obsm or emb.layers".format(obsm))
        else:
            raise Exception("Cannot have both obsm and layer, or layer or obsm requested do not exist")

        print("Matrix shape = {}".format(matrix.shape), flush=True)

        # If there are nan values in the matrix, we replace them by the median of the column

        if np.isnan(matrix).any():
            for i in range(matrix.shape[1]):
                col = matrix[:,i]
                median = np.nanmedian(col)
                col[np.isnan(col)] = median  
                print("Nan values in column {} replaced by median".format(i), flush=True)

        if "predicted_label" in self.emb.obs.columns:
            self.emb.obs.drop(columns=["predicted_label"], inplace=True)

        algo = algo.lower()
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
        elif algo == "mclust":
            model = None
            if matrix.shape[0] <= matrix.shape[1]:
                raise Exception(
                    "Cannot do mclust with a matrix with more rows than columns. Matrix shape = {}".format(matrix.shape)
                )
            self.emb.obsm["cluster"] = np.array(matrix)
            self.emb = utils.mclust_R(self.emb, used_obsm="cluster", num_cluster=n_clusters)
            self.emb.obs = self.emb.obs.rename(columns={"mclust": "predicted_label"})
        elif algo == "leiden":
            model = None
            self.emb.obsm["cluster"] = matrix
            res = utils.search_res(self.emb, n_clusters, use_rep="cluster", method=algo, start=0.1, end=3.0, increment=0.01)
            sc.tl.leiden(self.emb, random_state=0, resolution=res)
            self.emb.obs = self.emb.obs.rename(columns={"leiden": "predicted_label"})
        elif algo == "louvain":
            model = None
            self.emb.obsm["cluster"] = matrix
            res = utils.search_res(self.emb, n_clusters, use_rep="cluster", method=algo, start=0.1, end=3.0, increment=0.01)
            sc.tl.louvain(self.emb, random_state=0, resolution=res)
            self.emb.obs = self.emb.obs.rename(columns={"louvain": "predicted_label"})
        else:
            raise Exception("This clustering algorithm is not implemented.")

        self.emb.uns.update({"predicted_label": algo})
        if model is not None:
            self.emb.obs["predicted_label"] = model.fit_predict(matrix)
        else:
            del self.emb.obsm["cluster"]

        results = None
        if true_label_col is not None:
            # Filter labeled data that is not undetermined
            subset = self.emb.obs[~self.emb.obs["label"].isna()]
            subset = subset[subset["label"] != "nan"]
            subset = subset[subset["label"] != "undetermined"]
            results = {}
            results["ari"] = adjusted_rand_score(subset["label"].values, subset["predicted_label"].values)
            results["ami"] = adjusted_mutual_info_score(subset["label"].values, subset["predicted_label"].values)
            results["fmi"] = fowlkes_mallows_score(subset["label"].values, subset["predicted_label"].values)
            results["contingency_matrix"] = contingency_matrix(
                subset["label"].values, subset["predicted_label"].values
            ).tolist()

        return results

    def unsupervised_clustering_per_sample(
        self,
        sample="tumor",
        layer=None,
        obsm=None,
        algo="kmeans",
        n_clusters=6,
        var_ratio_threshold_for_svd=0.9,
        denoised_for_svd=False,
        u_comp_list=None,
        filter_label=True,
        assign_labels_spectral="cluster_qr",
        on_all_data=True,
        tsne_perplexity=30,
        tsne_n_components=2,
        umap_n_components=2,
        umap_n_neighbors=15,
        umap_min_dist=0.1,
        recompute_dim_red_on_subset=False,
        center_before_svd=False,
        scale_before_svd=False,
        multiply_by_variance=False,
    ):
        """Compute the ARI score for each group defined by the categorical column given by sample and also compute the ARI on the whole data.

        Args:
            sample (str, optional): Categorical column from emb.obs used to group the data. Defaults to "tumor".
            layer (str, optional): Layer name in emb.layers on which to run the clustering algorithm. If None and obsm is None, the raw emb is taken. Defaults to None.
            obsm (str, optional): Layer name in emb.obsm on which to run the clustering algorithm. If None and layer is None, the raw emb is taken. Defaults to None.
            algo (str, optional): Clustering algorithm to use. Defaults to "kmeans".
            n_clusters (int, optional): Number of clusters to consider for the clustering algorithm. If None, it is adapted with the number of different labels a sample has. Defaults to 6.
            var_ratio_threshold_for_svd (float, optional): If you use layer = "svd", it is the threshold of cumulative variance explained ratio until which we keep the ui. Between 0 and 1. Defaults to 0.9.
            denoised_for_svd (bool, optional): If you use layer = "svd", it chooses between the svd computed on the raw or the denoised matrix. Defaults to False.
            u_comp_list (list, optional): If layer = "svd", it is the list of the ui you want to keep for the matrix clustering. If None, the var_ratio_threshold_for_svd is used. Defaults None.
            filter_label (bool, optional): If we keep only embeddings with label (in emb.obs) that is not na. Defaults to True.
            assign_labels_spectral (str, optional): SpectralClustering assign_labels parameter. "kmeans", "discretize", "cluster_qr". Defaults to "discretize".
            on_all_data (bool, optional): If the clustering is also performed on all the data. Defaults to True.
        Returns:
            dict: Dictionary of all the metrics results with keys corresponding to groups and values to disctionary of metrics.
        """
        if filter_label:
            print("Filtering data with label different from nan...", flush=True)
            emb_without_nan = self.emb[self.emb.obs["label"] != "nan"].copy()
            print(f"Removed {len(self.emb) - len(emb_without_nan)} elements that were 'nan'", flush=True)
            emb_without_nan = emb_without_nan[~emb_without_nan.obs["label"].isna()]
            print(f"Removed {len(self.emb) - len(emb_without_nan)} elements that were np.NaN", flush=True)

            print(emb_without_nan.obs.value_counts("label"))
            for idx in emb_without_nan.obs.value_counts("label").index:
                print(f"{idx}: oftype {type(idx)}")
            groups = emb_without_nan[~emb_without_nan.obs["label"].isna()].obs.groupby(by=sample)
        else:
            groups = self.emb.obs.groupby(by=sample)
        ari_list, ami_list, fmi_list = [], [], []
        results = {}
        for name, group in groups:
            try:
                emb_subset = FindClusters()
                emb_subset.emb = self.emb[group.index, :]
                print("Processing sample {} with {} elements".format(name, len(emb_subset.emb)))

                if recompute_dim_red_on_subset:
                    # Recompute the SVD, TSNE or UMAP on the subset!!
                    if layer == 'svd':
                        print(f"Computing SVD on subset {name}", flush=True)
                        emb_subset.compute_svd(denoised=denoised_for_svd, center=center_before_svd, scale=scale_before_svd)
                    elif obsm == 'tsne':
                        print(f"Computing TSNE on subset {name}", flush=True)
                        emb_subset.compute_tsne(n_components=tsne_n_components, perplexity=tsne_perplexity)
                    elif obsm == 'umap':
                        print(f"Computing UMAP on subset {name}", flush=True)
                        emb_subset.compute_umap(n_components=umap_n_components, n_neighbors=umap_n_neighbors, min_dist=umap_min_dist)
                    else:
                        print(f"Using raw data on subset {name}", flush=True)

                if n_clusters is None:
                    nc = list(emb_subset.emb.obs["label"].unique())
                    nc = [e for e in nc if isinstance(e, str) and e != "undetermined" and e != "nan" and e != np.nan]
                    print(nc)
                    nc = len(nc)
                    print("Number of clusters = {}".format(nc))
                else:
                    nc = n_clusters
                metrics = emb_subset.unsupervised_clustering(
                    true_label_col="label",
                    algo=algo,
                    n_clusters=nc,
                    layer=layer,
                    obsm=obsm,
                    var_ratio_threshold_for_svd=var_ratio_threshold_for_svd,
                    u_comp_list=u_comp_list,
                    denoised_for_svd=denoised_for_svd,
                    center_before_svd=center_before_svd,
                    scale_before_svd=scale_before_svd,
                    multiply_by_variance=multiply_by_variance,
                )
                print("Metrics for {} = {}".format(name, metrics))
                results[name] = metrics
                ari_list.append(metrics["ari"])
                ami_list.append(metrics["ami"])
                fmi_list.append(metrics["fmi"])
                # TODO: Add the emb_subset.emb.obs["predicted_label"] to self.emb ?
            except Exception as e:
                print("Cannot compute clustering on {} sample: {}".format(name, e))

        if len(ari_list) > 0:
            results["mean"] = {"ari": np.mean(ari_list), "ami": np.mean(ami_list), "fmi": np.mean(fmi_list)}
        if on_all_data:
            try:
                results["all"] = self.unsupervised_clustering(
                    true_label_col="label",
                    algo=algo,
                    layer=layer,
                    obsm=obsm,
                    var_ratio_threshold_for_svd=var_ratio_threshold_for_svd,
                    u_comp_list=u_comp_list,
                    denoised_for_svd=denoised_for_svd,
                    assign_labels_spectral=assign_labels_spectral,
                )
            except Exception as e:
                print("Cannot compute clustering on all samples: {}".format(e))
        return results


    def hierarchical_clustering_dendrogram(self, ax=None, colored_by_cluster=True, palette=plt.get_cmap("tab10")):
        """Plot a dendrogram with 1 - spearman correlation as distance measure of all rows of the embedding. It is mainly used for embedding with low
        number of rows such as pseudo-bulk matrix.

        Args:
            ax (matplotlib.ax, optional): Matplotlib axe. Defaults to None.
            colored_by_cluster (bool, optional): If the y labels are colored according to predicted_label obs column. If False it is colored by patient (tumor emb.obs column). Defaults to True.
            palette (plt.cmap/dict, optional): Matplotlib palette, can also be a dictionnary with cluster/patient as key and color as values. Defaults to plt.get_cmap("tab10").

        Returns:
            matplotlib.ax: axe
        """
        show = False
        if ax is None:
            show = True
            ax = plt.gca()

        distance_matrix = 1 - self.emb.to_df().T.corr("spearman")

        model = AgglomerativeClustering(distance_threshold=0, n_clusters=None, metric="precomputed", linkage="average")
        model = model.fit(distance_matrix)
        ax.set_title("Hierarchical clustering dendrogram")

        # Extract labels
        labels = list(self.emb.obs.index)

        # Plot the dendrogram
        dendrogram_data = plot_dendrogram(
            model,
            truncate_mode="level",
            p=100,
            labels=labels,
            color_threshold=0,
            above_threshold_color="gray",
            orientation="right",
            ax=ax,
        )

        ax.set_xlabel("1 - Spearman correlation")
        # Define colors for clusters or patients
        try:
            if isinstance(palette, dict):
                colors_dict = palette
            else:
                if colored_by_cluster:
                    clusters = list(self.emb.obs["predicted_label"].unique())
                    colors_dict = {cluster: palette(i % 10) for i, cluster in enumerate(clusters)}
                else:
                    patients = list(self.emb.obs["tumor"].unique())
                    colors_dict = {patient: palette(i % 10) for i, patient in enumerate(patients)}

            y_labels = ax.get_ymajorticklabels()

            for label in y_labels:
                text = label.get_text()
                if colored_by_cluster:
                    cluster = text.split("_patient")[0]
                    label.set_color(colors_dict.get(cluster, "black"))
                    add_figname = "_colored_by_cluster"
                else:
                    patient = text.split("patient_")[-1]
                    label.set_color(colors_dict.get(patient, "black"))
                    add_figname = "_colored_by_patient"
        except Exception as e:
            add_figname = ""
            print("Y labels colored could not be set: {}".format(e))

        if show:
            if self.saving_plots:
                plt.savefig(
                    os.path.join(self.result_saving_folder, "hierarchical_clustering_dendrogram{}.pdf".format(add_figname)),
                    bbox_inches="tight",
                )
                plt.close()
            else:
                plt.show()
        else:
            return ax

    def hierarchical_clustering_dendrogram_per_group(self, group, ax=None, palette=plt.get_cmap("tab10")):
        """Plot a dendrogram with 1 - spearman correlation as distance measure of all rows of the embedding. It is mainly used for embedding with low
        number of rows such as pseudo-bulk matrix.

        Args:
            ax (matplotlib.ax, optional): Matplotlib axe. Defaults to None.
            colored_by_cluster (bool, optional): If the y labels are colored according to predicted_label obs column. If False it is colored by patient (tumor emb.obs column). Defaults to True.
            palette (plt.cmap/dict, optional): Matplotlib palette, can also be a dictionnary with cluster/patient as key and color as values. Defaults to plt.get_cmap("tab10").

        Returns:
            matplotlib.ax: axe
        """

        if ax is None:
            show = True
            ax = plt.gca()

        for g in self.emb.obs[group].unique():
            df = self.emb.obs[self.emb.obs[group] == g]

            distance_matrix = 1 - df.T.corr("spearman")

            model = AgglomerativeClustering(distance_threshold=0, n_clusters=None, metric="precomputed", linkage="average")
            model = model.fit(distance_matrix)
            ax.set_title("Hierarchical clustering dendrogram")

            # Extract labels
            labels = list(self.emb.obs.index)

            # Plot the dendrogram
            dendrogram_data = plot_dendrogram(
                model,
                truncate_mode="level",
                p=100,
                labels=labels,
                color_threshold=0,
                above_threshold_color="gray",
                orientation="right",
                ax=ax,
            )

            ax.set_xlabel("1 - Spearman correlation")
            # Define colors for clusters or patients
            try:
                if isinstance(palette, dict):
                    colors_dict = palette
                else:
                    clusters = list(self.emb.obs["predicted_label"].unique())
                    colors_dict = {cluster: palette(i % 10) for i, cluster in enumerate(clusters)}

                y_labels = ax.get_ymajorticklabels()

                for label in y_labels:
                    text = label.get_text()
                    cluster = text.split("_patient")[0]
                    label.set_color(colors_dict.get(cluster, "black"))
                    add_figname = "_colored_by_cluster"

            except Exception as e:
                add_figname = ""
                print("Y labels colored could not be set: {}".format(e))

            if show:
                if self.saving_plots:
                    plt.savefig(
                        os.path.join(self.result_saving_folder, "hierarchical_clustering_dendrogram{}.pdf".format(add_figname)),
                        bbox_inches="tight",
                    )
                    plt.close()
                else:
                    plt.show()
            else:
                return ax
            
    def clustering_across_umap_parameters(
        self,
        n_neighbors_list=[10, 30, 50, 100, 150, 200, 250, 300, 350, 400],
        min_dist_list=[0.001, 0.1],
        n_components=2,
        algo="kmeans",
        n_clusters=6,
        saving_folder=None,
    ):
        """Parallelized version of clustering_across_umap_parameters."""
        results = {}

        def process_combination(n_neighbors, min_dist):
            """Helper function to process a single combination of parameters."""
            print(f"n_neighbors: {n_neighbors}, min_dist: {min_dist}, n_components: {n_components}", flush=True)
            self.compute_umap(n_components=n_components, n_neighbors=n_neighbors, min_dist=min_dist)
            matrix = self.emb.obsm["umap"].copy()
            scores = clustering_scores(matrix=matrix, n_clusters=n_clusters, algo=algo, name="umap")
            return min_dist, n_neighbors, scores

        # Use ThreadPoolExecutor for parallel processing
        with ThreadPoolExecutor(max_workers=6) as executor:
            future_to_params = {
                executor.submit(process_combination, n_neighbors, min_dist): (n_neighbors, min_dist)
                for n_neighbors in n_neighbors_list
                for min_dist in min_dist_list
            }

            for future in as_completed(future_to_params):
                n_neighbors, min_dist = future_to_params[future]
                try:
                    min_dist, n_neighbors, scores = future.result()
                    if min_dist not in results:
                        results[min_dist] = {}
                    results[min_dist][n_neighbors] = scores
                except Exception as e:
                    print(f"Error processing n_neighbors={n_neighbors}, min_dist={min_dist}: {e}")

        results['samples'] = list(self.emb.obs.index)

        if saving_folder is None:
            return results
        else:
            filename = f"scores_umap_across_parameters_{n_clusters}_clusters.json"
            saving_path = os.path.join(saving_folder, filename)
            with open(saving_path, "w") as fp:
                json.dump(results, fp)
            print("Save results to {}".format(saving_path), flush=True)
                
    def select_best_umap_parameters(self, 
                                    saving_folder=None,
                                    files=None,
                                    score='silhouette_score'):
        
        
        if files is None:
            if saving_folder is not None:
                files = glob.glob(os.path.join(saving_folder, "scores_umap_across_parameters_*_clusters.json"))
            else: 
                raise Exception("No files to load")
        
        best_score = 0
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
                        if scores[score] > best_score:
                            best_score = scores[score]
                            best_n_neighbors = int(n_neighbors)
                            best_min_dist = float(min_dist)
                            n_clusters = int(file.split("_clusters")[0].split("_")[-1])
                            labels = scores['labels']
                        
        best_params = {'n_neighbors': best_n_neighbors, 
                        'min_dist': best_min_dist, 
                        'n_clusters': n_clusters, 
                        score: best_score,
                        'labels': labels,
                        'samples': results['samples']}
        
        return best_params
            

    def get_all_unsupervised_clustering_score_df(self, all_svd=True):
        """Compute a dataframe from the unsupervised_clustering_score_files with one line for each experiment.

        Args:
            all_svd (bool, optional): Add a data column to identify the different svd files. Defaults to True.

        Returns:
            pd.DataFrame: Dataframe with the following columns: "data", "cluster_number", "silhouette_score", "inertia", "davies_bouldin_score", "calinski_harabasz_score".
        """
        all_scores_df = None
        for file in self._unsupervised_clustering_score_files:
            with open(file) as f:
                scores_df = pd.DataFrame(json.load(f))
            if all_svd:
                scores_df["data"] = "svd" + file.split("_svd")[-1].split(".json")[0].split("_")[-1]
            if all_scores_df is None:
                all_scores_df = scores_df
            else:
                all_scores_df = pd.concat((all_scores_df, scores_df), axis=0)
        # print(all_scores_df.shape)
        return all_scores_df

    def plot_unsupervised_clustering_score(
        self,
        all_svd=True,
        algo="kmeans",
        suptitle=True,
        return_ax=False,
        add_best_cluster_vis=False,
        c=None,
    ):
        """One plot per metric (silhouette, inertia, DBI and CHI) and in each plot all the different experiments are plotted ("data" column of get_all_unsupervised_clustering_score_df).

        Args:
            all_svd (bool, optional): If we have several experiment. Defaults to True.
            algo (str, optional): Clustering algorithm used. Defaults to "kmeans".
            suptitle (bool, optional): Plot suptitle. Defaults to True.
            return_ax (bool, optional): If we want the plot ax to be returned. Defaults to False.
            add_best_cluster_vis (bool, optional): Add crosses on each best score of each experiment. Defaults to False.
            c (str/int, optional): Color of the lines and points. Defaults to None.

        Returns:
            matplotlib.ax: axe
        """
        fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(12, 12), sharex=True, layout="constrained")
        all_scores_df = self.get_all_unsupervised_clustering_score_df(all_svd=all_svd)
        print(f"All scores df: {all_scores_df}")
        for label, df in all_scores_df.groupby(by="data"):
            print(f"Plotting unsupervised clustering scores for {label}", flush=True)
            max_silh = df["silhouette_score"].max()
            cluster_number_max_silh = df.iloc[df["silhouette_score"].idxmax()]["cluster_number"]
            print(
                "Label {}: max silhouette score = {} for cluster number = {}".format(label, max_silh, cluster_number_max_silh)
            )
            min_db = df["davies_bouldin_score"].min()
            cluster_number_min_db = df.iloc[df["davies_bouldin_score"].idxmin()]["cluster_number"]
            print(
                "Label {}: min davies bouldin score = {} for cluster number = {}".format(label, min_db, cluster_number_min_db)
            )
            max_ch = df["calinski_harabasz_score"].max()
            cluster_number_max_ch = df.iloc[df["calinski_harabasz_score"].idxmax()]["cluster_number"]
            print(
                "Label {}: max calinski harabasz score = {} for cluster number = {}".format(
                    label, max_ch, cluster_number_max_ch
                )
            )
            ax[0, 0] = df.plot(x="cluster_number", y="silhouette_score", marker="o", ax=ax[0, 0], label=label, c=c)
            if "inertia" in df.columns:
                kn = KneeLocator(df["cluster_number"].values, y=df["inertia"].values, curve="convex", direction="decreasing")
                print("Label {}: elbow for {} clusters".format(label, kn.knee))
                if add_best_cluster_vis:
                    ax[0, 1] = df[df["cluster_number"] == kn.knee].plot(
                        x="cluster_number",
                        y="inertia",
                        marker="x",
                        ax=ax[0, 1],
                        c=c,
                        kind="scatter",
                        s=100,
                        label="",
                        linewidths=4,
                    )
                df.plot(x="cluster_number", y="inertia", marker="o", ax=ax[0, 1], label=label, c=c)
                ax[0, 1].get_legend().remove()
            ax[1, 0] = df.plot(x="cluster_number", y="davies_bouldin_score", marker="o", ax=ax[1, 0], label=label, c=c)
            ax[1, 1] = df.plot(x="cluster_number", y="calinski_harabasz_score", marker="o", ax=ax[1, 1], label=label, c=c)
            print("\n")
            if add_best_cluster_vis:
                df[df["cluster_number"] == cluster_number_max_silh].plot(
                    x="cluster_number",
                    y="silhouette_score",
                    marker="x",
                    ax=ax[0, 0],
                    c=c,
                    kind="scatter",
                    s=100,
                    label="",
                    linewidths=4,
                )
                df[df["cluster_number"] == cluster_number_min_db].plot(
                    x="cluster_number",
                    y="davies_bouldin_score",
                    marker="x",
                    ax=ax[1, 0],
                    c=c,
                    kind="scatter",
                    s=100,
                    label="",
                    linewidths=4,
                )
                df[df["cluster_number"] == cluster_number_max_ch].plot(
                    x="cluster_number",
                    y="calinski_harabasz_score",
                    marker="x",
                    ax=ax[1, 1],
                    c=c,
                    kind="scatter",
                    s=100,
                    label="",
                    linewidths=4,
                )
        ax[0, 0].set_ylabel("Silhouette score")
        ax[0, 1].set_ylabel("Inertia")
        ax[1, 0].set_ylabel("Davies Bouldin score")
        ax[1, 1].set_ylabel("Calinski Harabasz score")
        ax[1, 1].set_xlabel("")
        ax[1, 0].set_xlabel("")
        handles, labels = ax[0, 0].get_legend_handles_labels()
        labels = ["raw" if "raw" in l else l for l in labels]
        fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(1.1, 1), fancybox=True, shadow=True)
        ax[0, 0].get_legend().remove()
        ax[1, 0].get_legend().remove()
        # ax[1, 1].get_legend().remove()
        fig.supxlabel("Number of clusters")
        if suptitle:
            fig.suptitle("{} scores on the images embeddings to choose the ideal number of clusters".format(algo.capitalize()))
        if self.saving_plots:

            plt.savefig(
                os.path.join(
                    self.result_saving_folder,
                    "{}_benchmark_all_scores.pdf".format(algo.lower()),
                ),
                bbox_inches="tight",
            )
            plt.close()
        else:
            if return_ax:
                return ax
            else:
                plt.show()

    def compute_optimal_number_of_clusters(self, all_svd=False):
        """From get_all_unsupervised_clustering_score_df, the function goes through all experiments to find the best one. It first extract the optimal number of clutsers for all the metrics and all the
        experiments. It chooses the best number of svd component based on the lower standard deviation between the optimal number of clusters of the metrics. Finally, it calculates the final optimal number
        of clusters by computing the median of the optimal number of clusters from each metric and it rounds it according to mean of the optimal number of clusters from each metric. The results are saved
        in the info_optimal_number_of_clusters attribute.
        """
        all_scores_df = self.get_all_unsupervised_clustering_score_df(all_svd=all_svd)
        all_predicted_cluster = []
        for label, df in all_scores_df.groupby(by="data"):
            current_dict = dict()
            current_dict["data"] = label
            max_silh = df["silhouette_score"].max()
            current_dict["silhouette"] = df.iloc[df["silhouette_score"].idxmax()]["cluster_number"]

            min_db = df["davies_bouldin_score"].min()
            current_dict["DBI"] = df.iloc[df["davies_bouldin_score"].idxmin()]["cluster_number"]

            max_ch = df["calinski_harabasz_score"].max()
            current_dict["CHI"] = df.iloc[df["calinski_harabasz_score"].idxmax()]["cluster_number"]

            if "inertia" in df.columns:
                inertia = KneeLocator(
                    df["cluster_number"].values, y=df["inertia"].values, curve="convex", direction="decreasing"
                ).knee
                if inertia is not None:
                    current_dict["inertia"] = inertia
                    # metric_list = ["silhouette", "DBI", "CHI", "inertia"]
                else:
                    warnings.warn("Elbow score was impossible to compute.")
                    # metric_list = ["silhouette", "DBI", "CHI"]

            all_predicted_cluster.append(current_dict)

        all_predicted_cluster = pd.DataFrame(all_predicted_cluster).set_index("data")
        metric_list = all_predicted_cluster.columns

        all_predicted_cluster["std"] = all_predicted_cluster.apply(lambda r: np.std(r[metric_list].values[~np.isnan(r[metric_list].values)]), axis=1)
        all_predicted_cluster = all_predicted_cluster.sort_values(by="std")
        all_predicted_cluster["mean"] = all_predicted_cluster.apply(lambda r: np.mean(r[metric_list].values[~np.isnan(r[metric_list].values)]), axis=1)
        all_predicted_cluster["median"] = all_predicted_cluster.apply(lambda r: np.median(r[metric_list].values[~np.isnan(r[metric_list].values)]), axis=1)

        if all_svd:
            best_number_svd_comp = int(all_predicted_cluster.index[0].split("svd")[-1])
            print(
                "The number of svd components with the least variance in the number of optimal clusters is: {}".format(
                    best_number_svd_comp
                )
            )

            optimal_cluster_number = all_predicted_cluster.iloc[0, :]["median"]
            print(
                "\nThe median of the number of optimal clusters for svd {} is {}".format(
                    best_number_svd_comp, optimal_cluster_number
                )
            )
        else:
            best_exp = all_predicted_cluster.index[0]
            print("The experiment with the least variance in the number of optimal clusters is: {}".format(best_exp))
            optimal_cluster_number = all_predicted_cluster.iloc[0, :]["median"]
            print(
                "\nThe median of the number of optimal clusters for exp {} is {}".format(
                    best_exp, optimal_cluster_number
                )
            )

        if not optimal_cluster_number.is_integer():
            mean = all_predicted_cluster.iloc[0, :]["mean"]
            if optimal_cluster_number > mean:
                optimal_cluster_number = np.floor(optimal_cluster_number)
            elif optimal_cluster_number < mean:
                optimal_cluster_number = np.ceil(optimal_cluster_number)
            else:
                optimal_cluster_number = round(optimal_cluster_number)
            print("As the median was not a integer, the number of optimal cluster was rounded depending on the mean")

        print("The optimal number of cluster is {}".format(optimal_cluster_number))

        self.info_optimal_number_of_clusters["number"] = int(optimal_cluster_number)

        if all_svd:
            self.info_optimal_number_of_clusters["svd_comp"] = best_number_svd_comp
        else:
            self.info_optimal_number_of_clusters["exp"] = best_exp

        self.info_optimal_number_of_clusters_all_exp = all_predicted_cluster
    

    def clustering_across_different_n_clusters(
        self, 
        clusters_list=np.arange(2, 21, 1), 
        algo='kmeans',
        layer="svd", 
        denoised_svd=False, 
        svd_component_number=5, 
        saving_folder=None, 
        multiply_by_variance=True
    ):
        """Save json results file or return a dictionnary containing different metrics (Silhouette, CHI, DBI, inertia) to evaluate the kmeans clustering on different number of clusters k. The kmeans clustering
        can be applied on the U svd matrix with a given number of components.

        Args:
            clusters_list (list, optional): List containing the different number of clusters to test. Defaults to np.arange(2, 21, 1).
            layer (str, optional): Layer on which the kmeans clustering is applied. It can be a anndata.layers or svd. If None the raw matrix is used. Defaults to "svd".
            denoised_svd (bool, optional): If the denoised svd is used if layer svd is chosen. Defaults to False.
            svd_component_number (int, optional): The number of svd components if layer svd is chosen. Defaults to 5.
            saving_folder (str, optional): Saving folder where to store the json results file. If None, the results are returned as a dictionnary. Defaults to None.

        Raises:
            NotImplementedError: If a layer is is implemented.

        Returns:
            dict: Dictionnary of results if saving_folder is None.
        """
        if layer is None:
            print(f"Start {algo} clustering on the raw embedding matrix")
            matrix = self.emb.X
            add_to_filename = "_raw"
        elif layer == "svd":
            if svd_component_number is None:
                svd_component_number = U_df.shape[1]
            print(f"Start {algo} clustering for {'denoised' if denoised_svd else ''} svd with {svd_component_number} components")

            self.compute_svd(denoised=denoised_svd)

            if multiply_by_variance:
                U_df = self.svd["U_df"] * self.svd["S"]
            else:
                U_df = self.svd["U_df"].copy()

            matrix = U_df.iloc[:, 0:svd_component_number].to_numpy()
            add_to_filename = "_svd_{}".format(svd_component_number)
            if denoised_svd:
                add_to_filename = "_denoised" + add_to_filename
        elif layer in list(self.emb.layers.keys()):
            matrix = self.emb.layers[layer].to_numpy()
            add_to_filename = "_{}".format(layer)
        elif layer in list(self.emb.obsm.keys()):
            matrix = self.emb.obsm[layer].to_numpy()
            add_to_filename = "_{}".format(layer)
        else:
            raise NotImplementedError

        results = []
        for c in clusters_list:
            print("c: {}/{}".format(c, clusters_list[-1]), flush=True)
            current_results = clustering_scores(
                matrix=matrix, n_clusters=c, algo=algo, name="raw" if layer is None else add_to_filename[1:]
            )
            current_results["samples"] = list(self.emb.obs.index)
            results.append(current_results)
        if saving_folder is None:
            return results
        else:
            filename = "scores_{}_{}{}.json".format(np.array(clusters_list).min(), np.array(clusters_list).max(), add_to_filename)
            saving_path = os.path.join(saving_folder, filename)
            with open(saving_path, "w") as fp:
                json.dump(results, fp)
            print("Save results to {}".format(saving_path), flush=True)

    def UMAP_validation_unsupervised_clustering(self, 
                                            algo='kmeans', 
                                            n_neighbors_list=[10, 30, 50, 100, 150, 200, 250, 300, 350, 400],
                                            min_dist_list=[0.001, 0.1], 
                                            n_components_list=[2], 
                                            n_clusters=None):
        """
        Perform UMAP validation using unsupervised clustering to find the best parameters.
        This method iterates over combinations of UMAP parameters (n_neighbors, min_dist, n_components)
        and performs unsupervised clustering using the specified algorithm. It returns the best
        Adjusted Rand Index (ARI) score and the corresponding parameters.
        Parameters:
        algo (str): The clustering algorithm to use.
        n_neighbors_list (list): List of n_neighbors values to try.
        min_dist_list (list): List of min_dist values to try.
        n_components_list (list): List of n_components values to try.
        n_clusters (int, optional): The number of clusters to form. Default is None.
        Returns:
        tuple: A tuple containing the best ARI score (float) and a dictionary of the best parameters.
        """

        all_runs = []
        best_ari = -1
        best_params = {'n_neighbors': None, 'min_dist': None, 'n_components': None}
        for n_neighbors in n_neighbors_list:
            for min_dist in min_dist_list:
                print(f"n_neighbors: {n_neighbors}, min_dist: {min_dist}")
                for n_components in n_components_list:
                    self.compute_umap(n_components=n_components, n_neighbors=n_neighbors, min_dist=min_dist)
                    ari = self.unsupervised_clustering(n_clusters=n_clusters, algo=algo, obsm='umap')['ari']

                    if ari > best_ari:
                        best_ari = ari
                        best_params = {'n_neighbors': n_neighbors, 'min_dist': min_dist, 'n_components': n_components}
                        
                    all_runs.append([n_neighbors, min_dist, n_components, ari])
                    
        df_all_runs = pd.DataFrame(all_runs, columns=['n_neighbors', 'min_dist', 'n_components', 'ari'])
        return best_ari, best_params, df_all_runs
    
    def unsupervised_clustering_no_labels(self,
                                            clusters_list=[3, 4, 5, 6, 7, 8, 9, 10], 
                                            layer="svd", 
                                            algo = "kmeans",
                                            denoised_svd=False, 
                                            svd_component_number=5, 
                                            saving_folder=None, 
                                            multiply_by_variance=True,
                                            extension='pdf'):

        if saving_folder is not None:
            if not os.path.exists(saving_folder):
                os.makedirs(saving_folder)

        if layer is None:
            add_to_filename = ""
        elif layer == "svd":
            add_to_filename = "_svd_{}".format(svd_component_number)
            if denoised_svd:
                add_to_filename = "_denoised" + add_to_filename
        elif layer in list(self.emb.layers.keys()):
            add_to_filename = "_{}".format(layer)
        elif layer in list(self.emb.obsm.keys()):
            add_to_filename = "_{}".format(layer)
        else:
            raise NotImplementedError

        if 'start_width_origin' in self.emb.obs.columns:
            self.emb.obs['start_width_origin'] = self.emb.obs['start_width_origin'].astype(float)

        if 'start_height_origin' in self.emb.obs.columns:
            self.emb.obs['start_height_origin'] = self.emb.obs['start_height_origin'].astype(float)

        filename = "scores_{}_{}{}.json".format(np.array(clusters_list).min(), np.array(clusters_list).max(), add_to_filename)

        if os.path.exists(os.path.join(saving_folder, filename)):
            print("k-means for n clusters already computed already exists")
        else:

            # Unsupervised clustering
            self.clustering_across_different_n_clusters(clusters_list=clusters_list,
                                                        layer=layer,
                                                        algo=algo,
                                                        denoised_svd=denoised_svd,
                                                        svd_component_number=svd_component_number,
                                                        saving_folder=saving_folder,
                                                        multiply_by_variance=multiply_by_variance)

        self._unsupervised_clustering_score_files = [os.path.join(saving_folder, filename)]

        # Plot scores
        self.plot_unsupervised_clustering_score(return_ax=True, all_svd=False, algo=algo)
        plt.savefig(os.path.join(saving_folder,  f"unsupervised_scores_{np.array(clusters_list).min()}_{np.array(clusters_list).max()}{add_to_filename}.{extension}"), bbox_inches='tight')

        # Compute optimal number of clusters
        self.compute_optimal_number_of_clusters(all_svd=False)

        self.unsupervised_clustering(n_clusters=self.info_optimal_number_of_clusters['number'], 
                                     layer=layer,
                                     algo=algo,
                                     multiply_by_variance=multiply_by_variance, 
                                     u_comp_list=[f"u{i+1}" for i in range(svd_component_number)],
                                     denoised_for_svd=denoised_svd,
                                     true_label_col=None)

