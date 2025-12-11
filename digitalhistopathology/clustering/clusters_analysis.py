#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances
from digitalhistopathology.utils import quantized_wasserstein
import os
import seaborn as sns
import matplotlib.pyplot as plt

class ClustersAnalysis:

    def __init__(self, emb=None, saving_plots=False, result_saving_folder=None):
        """
        Initialize the ClustersAnalysis class.

        Parameters:
        -----------
        emb : EmbDimensionalityReduction
            An instance of the EmbDimensionalityReduction class.
        saving_plots : bool, optional
            Whether to save plots or not. Default is False.
        result_saving_folder : str, optional
            The folder where results will be saved. Default is None.
        """
        super().__init__(emb=emb, saving_plots=saving_plots, result_saving_folder=result_saving_folder)

    def extract_representative_samples(self, 
                                       layer='umap', 
                                       n_samples=10, 
                                       label_col='predicted_label'):
        
        if layer is None:
            matrix = self.emb.X
        elif layer in list(self.emb.layers.keys()):
            matrix = self.emb.layers[layer]
        elif layer in list(self.emb.obsm.keys()):
            matrix = self.emb.obsm[layer]
        else:
            raise NotImplementedError
        
        
        labels = self.emb.obs[label_col].values
        
        df_matrix = pd.DataFrame(matrix, index=self.emb.obs.index)
        df_matrix['label'] = labels 
        
        centroids = df_matrix.groupby('label').mean()
        
        distances = pairwise_distances(df_matrix[[col for col in df_matrix if col != 'label']], centroids) 
        
        df_matrix['distance_to_centroid'] = distances.min(axis=1)
        
        closest_samples_dict = {}
        for label in df_matrix['label'].unique():
            cluster_samples = df_matrix[df_matrix['label'] == label]
            closest_samples = cluster_samples.nsmallest(n_samples, 'distance_to_centroid')
            closest_samples_dict[int(label)] = closest_samples.index.tolist()
        
        return closest_samples_dict
    
    
    def compute_quantized_wasserstein_distance_between_clusters(self, 
                                                                cluster_col, 
                                                                layer=None, 
                                                                ref_space=None, 
                                                                k=5000):
        
        """
        Compute the Sinkhorn distance between clusters in the embedding matrix.
        Parameters:
        -----------
        cluster_col : str
            The column name in `self.emb.obs` that contains cluster labels.
        layer : str, optional
            The layer of the embedding matrix to use. If None, use `self.emb.X`.
            If the layer is in `self.emb.layers`, use that layer.
            If the layer is in `self.emb.obsm`, use that layer.
        normalize_by_max_diameter : bool, optional
            Whether to normalize the Sinkhorn distance by the maximum Euclidean distance
            between points in the embedding matrix. Default is True.
        Returns:
        --------
        df_sinkhorn : pd.DataFrame
            A DataFrame containing the Sinkhorn distances between clusters.
        df_sinkhorn_std : pd.DataFrame
            A DataFrame containing the standard deviations of the Sinkhorn distances
            between clusters.
        """


        if layer is None:
            matrix = self.emb.X
        elif layer in list(self.emb.layers.keys()):
            matrix = self.emb.layers[layer].to_numpy()
        elif layer in list(self.emb.obsm.keys()):
            matrix = self.emb.obsm[layer].to_numpy()
        else:
            raise NotImplementedError

        all_clusters = self.emb.obs[cluster_col].unique()
        
        min_sample_size = self.emb.obs[cluster_col].value_counts().min()
        
        if k > min_sample_size:
            k = min_sample_size
               
        df_w = pd.DataFrame(index=all_clusters, columns=all_clusters)
        
        for i in range(len(all_clusters)):
            for j in range(i+1, len(all_clusters)):
                if i != j:
                    cluster1 = all_clusters[i]
                    cluster2 = all_clusters[j]

                    print(f"Compute quantized wasserstein distance between {cluster1} and {cluster2}", flush=True)

                    samples_cluster1 = self.emb.obs[self.emb.obs[cluster_col] == cluster1].index
                    idx_samples_cluster1 = [self.emb.obs.index.get_loc(x) for x in samples_cluster1]
                    samples_cluster2 = self.emb.obs[self.emb.obs[cluster_col] == cluster2].index
                    idx_samples_cluster2 = [self.emb.obs.index.get_loc(x) for x in samples_cluster2]
                    
                    print(f"Number of samples in {cluster1}: {len(idx_samples_cluster1)}")
                    print(f"Number of samples in {cluster2}: {len(idx_samples_cluster2)}")

                    w = quantized_wasserstein(matrix=matrix, 
                                              idx_samples_cluster1=idx_samples_cluster1,
                                              idx_samples_cluster2=idx_samples_cluster2,
                                              ref_space=ref_space,
                                              k=k)
                        
                    df_w.loc[cluster1, cluster2] = w
                    df_w.loc[cluster2, cluster1] = w

        df_w.replace(np.nan, 0, inplace=True)

        return df_w
    
    def get_predicted_clusters_infos(
        self, groupby="predicted_label", fraction=True, predicted_clusters_to_consider=None, label_column="predicted_label"
    ):
        """Get information (size, labeled size, spot size, fraction of patients/predicted_cluster, fraction of labels) about each predicted clutser for groupby predicted_clutser or each patient
        for groupby tumor from unsupervised clustering.

        Args:
            groupby (str, optional): If the information are grouped per predicted cluster ("predicted_label") or patient ("tumor"). Other emb.obs column can also be used. Defaults to "predicted_label".
            fraction (bool, optional): If the amount of each patient or label is computed in term of fraction. If false, counts are returned. Defaults to True.
            predicted_clusters_to_consider (list, optional): List of predicted clusters to consider, it filters out all predicted clusters not present in the list. Defaults to None.
            label_column (str, optional): Column in emb.obs where the labels are. Defaults to "predicted_label".

        Returns:
            pd.DataFrame: DataFrame with all the information as columns and each predicted cluster or patient (groupby) as rows.
        """
        if groupby not in list(self.emb.obs.columns):
            raise Exception(
                "You must have the groupby column in the emb.obs column generated from clustering, see unsupervised_clutsering function from Embedding class"
            )

        groups = self.emb.obs.groupby(by=groupby)
        first_fraction_groupby = "tumor" if "predicted_label" in groupby else label_column
        results = []
        for name, group in groups:
            if predicted_clusters_to_consider is not None:
                group = group[group[label_column].apply(lambda l: l in predicted_clusters_to_consider)]
            first_fraction_dict = (
                (group.groupby(by=first_fraction_groupby)["name_origin"].count() / len(group)).to_dict()
                if fraction
                else (group.groupby(by=first_fraction_groupby)["name_origin"].count()).to_dict()
            )

            labeled_group = group[(~group["label"].isna()) & (group["label"] != "undetermined")]
            label_fraction_dict = (
                (labeled_group.groupby(by="label")["name_origin"].count() / len(labeled_group)).to_dict()
                if fraction
                else (labeled_group.groupby(by="label")["name_origin"].count()).to_dict()
            )
            current_results = {
                "name": name,
                "size": len(group),
                "spots_size": len(group),
                "labeled_size": len(labeled_group),
            }
            current_results.update(label_fraction_dict)
            current_results.update(first_fraction_dict)
            results.append(current_results)

        cluster_info_df = pd.DataFrame(results)

        if "predicted_label" in groupby:
            count = (
                self.emb.obs[(~self.emb.obs["label"].isna()) & (self.emb.obs["label"] != "undetermined")]
                .groupby(by="label")
                .count()["name_origin"]
            )
            # if "invasive cancer" not in count.index:
            #     count["invasive cancer"] = 0
            random_clustering_threshold_invasive_cancer = count["invasive cancer"] / count.sum()
            print(
                "Fraction of invasive cancer in cluster if random cluster assignement = {}".format(
                    random_clustering_threshold_invasive_cancer
                )
            )
            self._random_clustering_threshold_invasive_cancer = random_clustering_threshold_invasive_cancer
            try:
                if isinstance(cluster_info_df.name[0], int):
                    cluster_info_df.index = cluster_info_df.name
                elif isinstance(cluster_info_df.name[0], str):
                    cluster_info_df.index = cluster_info_df.name.apply(lambda s: int(s.split(" ")[-1]))
            except Exception as e:
                print("Impossible to rename the clusters: {}".format(e))

        return cluster_info_df

    def pie_charts_each_patient_across_selected_clusters(
        self, selected_clusters_list, label_column="predicted_label", palette=None
    ):
        """Pie chart for each patient according to their selected predicted cluster content.

        Args:
            selected_clusters_list (list): List of the selected cluster names.
            label_column (str, optional): Column in emb.obs where the labels are. Defaults to "predicted_label".
            palette (dict, optional): Palette dictionary with label as key and color as item. Defaults to None.
        """
        nrows = int(np.ceil((len(self.emb.obs.tumor.unique()) / 4)))
        fig, axes = plt.subplots(
            nrows=nrows, ncols=4, figsize=(20, 5 * nrows), layout="constrained"
        )
        axes = axes.reshape(-1)
        # To get the same colors for each patient
        if selected_clusters_list is not None and palette is None:
            palette = sns.color_palette()
            palette = {c: palette[i] for i, c in enumerate(selected_clusters_list)}
        cluster_df_patient = self.get_predicted_clusters_infos(
            groupby="tumor",
            fraction=True,
            predicted_clusters_to_consider=selected_clusters_list,
            label_column=label_column,
        )
        labels_list = self.emb.obs.label.dropna().unique().tolist()
        labels_list = [l for l in labels_list if l != "undetermined"]
        for index, row in cluster_df_patient.iterrows():
            # radius = row["size"] / cluster_df_patient["size"].sum()
            cancer_row = row[selected_clusters_list].dropna()
            axes[index].pie(
                cancer_row.values,
                labels=cancer_row.index,
                autopct="%.0f%%",
                textprops={"size": "smaller"},
                colors=[palette[d] for d in list(cancer_row.index)] if palette is not None else None,
            )
            axes[index].set_title("Patient {}, n={}".format(row["name"], int(row["size"])))

        for i in range(index + 1, axes.shape[0]):
            axes[i].remove()

        if self.saving_plots:
            plt.savefig(
                os.path.join(
                    self.result_saving_folder,
                    "unsupervised_clustering_pie_chart_each_patient_across_selected_clusters.pdf",
                ),
                bbox_inches="tight",
            )
            plt.close()
        else:
            plt.show()

    def barplots_predicted_clusters(self, groupby="predicted_cluster", palette_label=None, threshold_invasive_dashed_line=True):
        """Two barplots. The top one shows the distibution of patients (tumor column) across the predicted clusters and the bottom plot
        shows the distibution of labels (label column) across the predicted clusters.

        Args:
            groupby (str, optional): self.emb.obs categorical column to group the data. Defaults to predicted_cluster.
            palette_label (dict, optional): Palette dictionary with label as key and color as item. Defaults to None.
            threshold_invasive_dashed_line (bool, optional): If a horizontal dashed line is plotted at the height of the fraction of invasive cancer labeled patches of the whole datatset. Defaults to True.
        """
        self.emb.obs.replace("nan", np.NaN, inplace=True)

        cluster_info_df = self.get_predicted_clusters_infos(groupby=groupby, fraction=True, predicted_clusters_to_consider=None)
        whole_dataset_serie_patient = self.emb.obs.groupby(by="tumor")["name_origin"].count() / self.emb.shape[0]
        cluster_info_df.loc["Whole", whole_dataset_serie_patient.index] = whole_dataset_serie_patient.values

        labeled_emb = self.emb[~self.emb.obs.label.isna(), :]
        labeled_emb = labeled_emb[labeled_emb.obs.label != "undetermined", :]
        whole_dataset_serie_label = self.emb.obs.groupby(by="label")["name_origin"].count() / labeled_emb.shape[0]
        cluster_info_df.loc["Whole", whole_dataset_serie_label.index] = whole_dataset_serie_label.values
        cluster_info_df.loc["Whole", "size"] = self.emb.shape[0]
        cluster_info_df.loc["Whole", "labeled_size"] = labeled_emb.shape[0]

        cluster_info_df.index = ["{} (n={})".format(index, int(row["size"])) for index, row in cluster_info_df.iterrows()]
        patient_list = self.emb.obs.tumor.dropna().unique().tolist()
        cluster_number = len(cluster_info_df)
        high_patient_number_adaptation = len(patient_list) > 16
        n_rows = cluster_number + 1 if high_patient_number_adaptation else 2
        fig, axes = plt.subplots(
            nrows=n_rows, figsize=(12, n_rows * 6), layout="constrained", sharey=True
        )
        if not high_patient_number_adaptation:
            cluster_info_df[patient_list].plot.bar(ax=axes[0])
            axes[0].set_title("Fraction of patient in each cluster")
            axes[0].set_xlabel("Clusters")
            axes[0].set_ylabel("Fraction")
            axes[0].set_ylim(0, cluster_info_df[patient_list].max().max() * 1.1)

            axes[0].legend(
                loc="upper center",
                bbox_to_anchor=(1.05, 1),
                fancybox=True,
                shadow=True,
                ncol=1,
                title="Patient",
            )

        else:
            # Dynamically generate colors based on the number of patients
            num_patients = len(patient_list)
            colors = ["black"] * num_patients
            for i, (cluster_index, row) in enumerate(cluster_info_df.iterrows()):
                row = row[patient_list]
                row.plot.bar(ax=axes[i], color=colors)
                axes[i].set_title("Cluster {} ".format(cluster_index))
                axes[i].set_ylim(0, row.max() * 1.1)
        label_list = self.emb.obs.label.dropna().unique().tolist()
        label_list = [l for l in label_list if l != "undetermined"]
        cluster_info_df.index = [
            "{} (n={})".format(index.split(" ")[0], int(row["labeled_size"]))
            for index, row in cluster_info_df.iterrows()
        ]
        cluster_info_df[label_list].plot.bar(ax=axes[-1], color=palette_label)

        if threshold_invasive_dashed_line:
            axes[-1].hlines(
                y=self._random_clustering_threshold_invasive_cancer,
                xmin=-10,
                xmax=30,
                color="black",
                linestyles="dashed",
            )
        axes[-1].set_title("Fraction of label in each cluster")
        axes[-1].set_xlabel("Clusters")
        axes[-1].set_ylabel("Fraction")
        axes[-1].legend(
            loc="upper center",
            bbox_to_anchor=(1.1, 1),
            fancybox=True,
            shadow=True,
            ncol=1,
            title="Label",
        )

        if self.saving_plots:
            plt.savefig(
                os.path.join(
                    self.result_saving_folder,
                    "unsupervised_clustering_fraction_in_each_cluster_{}".format(groupby) + ".pdf",
                ),
                bbox_inches="tight",
            )
            plt.close()
        else:
            plt.show()

 
