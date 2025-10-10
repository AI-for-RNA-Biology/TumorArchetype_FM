#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#

import os
import glob
import json
import numpy as np
import pandas as pd
import seaborn as sns
import anndata
from digitalhistopathology.embeddings.image_embedding import ImageEmbedding
from digitalhistopathology.embeddings.gene_embedding import GeneEmbedding
from digitalhistopathology.datasets.real_datasets import HER2Dataset, TNBCDataset
from digitalhistopathology.benchmark.benchmark_base import BenchmarkBase
from digitalhistopathology.benchmark.benchmark_utils import get_optimal_cluster_number_one_model

import matplotlib.pyplot as plt

from sklearn.metrics import adjusted_rand_score

class BenchmarkInvasive(BenchmarkBase):
    def __init__(self, 
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
                 molecular_emb_path=None,
                 molecular_name='gene',
                 ref_model_emb=None,
                 algo='kmeans',
                 min_cluster=4,
                 max_cluster=10,
                 cluster_step=1,
                 dataset="HER2"):
        
        super().__init__(path_to_pipeline=path_to_pipeline, 
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
                         dataset=dataset)
        
        self.molecular_name = molecular_name
        self.algo = algo

        if dataset == "HER2":
            self.dataset = HER2Dataset()
        elif dataset == "TNBC":
            self.dataset = TNBCDataset()

        if molecular_emb_path is not None:
            if molecular_emb_path.endswith('.h5ad'):
                molecular_emb = anndata.read_h5ad(molecular_emb_path)
                self.molecular_emb = GeneEmbedding(emb=molecular_emb)
            elif molecular_emb_path.endswith('.parquet'):
                self.molecular_emb = self.dataset.get_gene_embeddings(molecular_emb_path=molecular_emb_path)
            else:
                print(f"Loading molecular embeddings from {molecular_emb_path} not implemented")

            print(f"Loaded molecular embeddings with shape {self.molecular_emb.emb.X.shape}")
            print(f"Molecular embeddings have labels: {self.molecular_emb.emb.obs.columns.tolist()}")

        else:
            self.molecular_emb = None

        if ref_model_emb is not None:
            self.ref_model_emb = None
        else:
            self.ref_model_emb = anndata.read_h5ad(ref_model_emb).X
        
        self.min_cluster = min_cluster
        self.max_cluster = max_cluster
        self.cluster_step = cluster_step

    
    def compare_invasive_patches_selection(self):
        """
        Compare the invasive patches selection across different models.
        This method calculates and visualizes the percentage and number of common invasive patches 
        between different models' embeddings. It also plots the number of invasive patches per model.
        The method performs the following steps:
        1. Loads invasive image embeddings if not already loaded.
        2. Initializes matrices to store the percentage and number of common patches between models.
        3. Iterates over pairs of models to compute the common and not common patches.
        4. Fills the matrices with the calculated values.
        5. Computes the number of invasive patches for each model.
        6. Plots and saves heatmaps for the percentage and number of common patches.
        7. Plots and saves a bar plot for the number of invasive patches per model.
        The generated plots are saved in the specified saving folder with the given file extension.
        Raises:
            ValueError: If `self.invasive_image_embeddings` is None and cannot be loaded.
        """

        print("Comparing invasive patches selection across models...", flush=True)

        if self.invasive_image_embeddings is None:
            self.load_invasive_image_embeddings()

        percentage_common_patches_matrix = np.zeros((len(self.pipelines_list), len(self.pipelines_list)))
        number_common_patches_matrix = np.zeros((len(self.pipelines_list), len(self.pipelines_list)))   
        for idx_model1 in range(len(self.pipelines_list)):
            for idx_model2 in range(len(self.pipelines_list)):
                invasive_patches1 = self.invasive_image_embeddings[self.pipelines_list[idx_model1]].emb.obs_names
                invasive_patches2 = self.invasive_image_embeddings[self.pipelines_list[idx_model2]].emb.obs_names
                common_patches = set(invasive_patches1).intersection(set(invasive_patches2))
                not_common_patches = set(invasive_patches1).union(set(invasive_patches2)) - common_patches
                percentage_common_patches_matrix[idx_model1, idx_model2] = len(common_patches)*100 / (len(not_common_patches) + len(common_patches))
                number_common_patches_matrix[idx_model1, idx_model2] = len(common_patches)

        number_invasive_patches = {}
        for model in self.pipelines_list:
            number_invasive_patches[model] = len(self.invasive_image_embeddings[model].emb.obs_names)
            print(f"Model {model} has {len(self.invasive_image_embeddings[model].emb.obs_names)} patches")

        plt.figure()
        sns.heatmap(percentage_common_patches_matrix, annot=True, xticklabels=self.pipelines_list, yticklabels=self.pipelines_list, cmap='Reds')
        plt.title("Percentage of common patches between invasive cancer embeddings")
        plt.savefig(os.path.join(self.saving_folder, f"percentage_common_patches_between_invasive_embeddings.{self.extension}"), bbox_inches='tight')

        plt.figure()
        sns.heatmap(number_common_patches_matrix, annot=True, xticklabels=self.pipelines_list, yticklabels=self.pipelines_list, cmap='Reds')
        plt.title("Number of common patches between invasive cancer embeddings")
        plt.savefig(os.path.join(self.saving_folder, f"number_common_patches_between_invasive_embeddings.{self.extension}"), bbox_inches='tight')

        plt.figure()
        sns.barplot(x=number_invasive_patches.keys(), y=number_invasive_patches.values())
        plt.title("Number of invasive patches per model")
        plt.xticks(rotation=90)
        plt.savefig(os.path.join(self.saving_folder, f"number_invasive_patches_per_model.{self.extension}"), bbox_inches='tight')


    
    def sub_invasive_cancer_clustering_pipeline(self,
                                                file, 
                                                model):
        
        ## Retrieve the number of clusters
        n_cluster = int(file.split("_clusters")[0].split("_")[-1])
        
        ## Retrieve the best UMAP parameters
        best_params = self.invasive_image_embeddings[model].select_best_umap_parameters(files=[file])
        
        ## Get the labels
        self.invasive_image_embeddings[model].emb.obs['predicted_label'] = best_params['labels']
        
        ## Save the labels as a DataFrame
        labels = pd.DataFrame(best_params["labels"], best_params["samples"])
        labels.columns = ["predicted_label"]
        labels['label'] = ""
        labels.to_csv(os.path.join(self.saving_folder, model, f"invasive_labels_{n_cluster}_clusters_umap_min_dist_{best_params['min_dist']}_n_neighbors_{best_params['n_neighbors']}.csv"))
        
        
        ## Load the molecular embeddings: the full one and the one with the labels
        molecular_emb = GeneEmbedding()
        molecular_emb.emb = self.molecular_emb.emb.copy()
        molecular_emb.emb.obs['predicted_label'] = [labels.loc[idx, 'predicted_label'] if idx in labels.index else "not invasive" for idx in molecular_emb.emb.obs.index]
        
        molecular_emb_labeled = GeneEmbedding()
        molecular_emb_labeled.emb = molecular_emb.emb[~molecular_emb.emb.obs['predicted_label'].isna()]
        molecular_emb_labeled.emb = molecular_emb_labeled.emb[~(molecular_emb_labeled.emb.obs['predicted_label'] == 'not invasive')]
        
        if not 'umap' in molecular_emb_labeled.emb.obsm.keys():
            print("Computing UMAP for the molecular embedding...", flush=True)
            molecular_emb_labeled.compute_umap(palette=self.dataset.PALETTE)
        
        if "tumor" not in molecular_emb_labeled.emb.obs.columns:
            molecular_emb_labeled.emb.obs['tumor'] = [idx.split("_")[0] for idx in molecular_emb_labeled.emb.obs.index]
        
        ## Save and plot the representative patches for each cluster
        representative_patches = self.invasive_image_embeddings[model].extract_representative_samples(layer='umap', n_samples=10, label_col='predicted_label')
        with open(os.path.join(self.saving_folder, model, f"representative_patches_{n_cluster}_clusters_umap.json"), 'w') as f:
            json.dump(representative_patches, f)
            
        # plot representative patches
        self.invasive_image_embeddings[model].plot_representative_samples(representative_samples=representative_patches, 
                                                                        saving_filename=os.path.join(self.saving_folder, model, f"representative_patches_{n_cluster}_clusters_umap.{self.extension}"))
            

        ## Plot UMAPs
        for label in ['predicted_label', self.group]:
            if label == 'predicted_label':
                
                if self.invasive_image_embeddings[model].emb.obs["predicted_label"].nunique() > 5:
                    palette = sns.color_palette("Accent", n_colors=self.invasive_image_embeddings[model].emb.obs["predicted_label"].nunique())
                else:
                    palette = ['#F9A11B', '#31C4F3', '#ACB5B6', '#EC2A90', '#66BB46', 'white']
                palette = {c: palette[i] for i, c in enumerate(sorted(
                    self.invasive_image_embeddings[model].emb.obs["predicted_label"][
                        ~self.invasive_image_embeddings[model].emb.obs["predicted_label"].isna()
                    ]
                    .unique()
                    .tolist()
                ))}
            else:
                palette = 'Accent'
                
            
                
            print(f"Plotting UMAP kmeans for model {model} - n_cluster {n_cluster} - colored by {label}...", flush=True)
            plt.figure()
            sns.scatterplot(x=self.invasive_image_embeddings[model].emb.obsm["umap"][:,0],
                            y=self.invasive_image_embeddings[model].emb.obsm["umap"][:,1],
                            hue=self.invasive_image_embeddings[model].emb.obs[label],
                            palette=palette,
                            s=10)
            plt.title(f"UMAP kmeans - {model} - colored by {label}", weight='bold')
            sns.despine()
            plt.legend(loc='upper right')
            plt.savefig(os.path.join(self.saving_folder, model, f"umap_kmeans_{label}_{model}_{n_cluster}_clusters.{self.extension}"), bbox_inches='tight')


            plt.figure()
            sns.kdeplot(x=self.invasive_image_embeddings[model].emb.obsm["umap"][:,0],
                            y=self.invasive_image_embeddings[model].emb.obsm["umap"][:,1],
                            hue=self.invasive_image_embeddings[model].emb.obs[label],
                            palette=palette)
            plt.title(f"UMAP kmeans - {model} - colored by {label}", weight='bold')
            sns.despine()
            plt.legend(loc='upper right')
            plt.savefig(os.path.join(self.saving_folder, model, f"umap_kde_kmeans_{label}_{model}_{n_cluster}_clusters.{self.extension}"), bbox_inches='tight')
        
        
            plt.figure()
            sns.scatterplot(x=molecular_emb_labeled.emb.obsm["umap"][:,0],
                            y=molecular_emb_labeled.emb.obsm["umap"][:,1],
                            hue=molecular_emb_labeled.emb.obs[label],
                            palette=palette,
                            s=10)
            plt.title(f"UMAP kmeans molecular - {model} - colored by {label}", weight='bold')
            sns.despine()
            # plt.legend(loc='upper right')
            plt.savefig(os.path.join(self.saving_folder, model, f"umap_{self.molecular_name}_kmeans_{label}_{model}_{n_cluster}_clusters.{self.extension}"), bbox_inches='tight')


            plt.figure()
            sns.kdeplot(x=molecular_emb_labeled.emb.obsm["umap"][:,0],
                            y=molecular_emb_labeled.emb.obsm["umap"][:,1],
                            hue=molecular_emb_labeled.emb.obs[label],
                            palette=palette)
            
            plt.title(f"UMAP kmeans molecular - {model} - colored by {label}", weight='bold')
            sns.despine()
            # plt.legend(loc='upper right')
            plt.savefig(os.path.join(self.saving_folder, model, f"umap_{self.molecular_name}_kde_kmeans_{label}_{model}_{n_cluster}_clusters.{self.extension}"), bbox_inches='tight')
        ## barplots and pie charts 
        
            self.invasive_image_embeddings[model].emb.obs.replace("nan", np.NaN, inplace=True)
            self.invasive_image_embeddings[model].barplots_predicted_clusters(
                groupby="predicted_label",
                threshold_invasive_dashed_line=False,
                palette_label=sns.color_palette("Accent"),
            )
            plt.savefig(os.path.join(self.saving_folder, model, f"barplot_predicted_clusters_umap_model_{model}_{n_cluster}_clusters.{self.extension}"), bbox_inches='tight')

            self.invasive_image_embeddings[model].pie_charts_each_patient_across_selected_clusters(
                    selected_clusters_list=sorted(
                        self.invasive_image_embeddings[model].emb.obs["predicted_label"][
                            ~self.invasive_image_embeddings[model].emb.obs["predicted_label"].isna()
                        ]
                        .unique()
                        .tolist()
                    ),
                    label_column="predicted_label",
                )
            plt.savefig(os.path.join(self.saving_folder, model, f"piechart_predicted_clusters_umap_model_{model}_{n_cluster}_clusters.{self.extension}"), bbox_inches='tight')
            
        # Compute quantized wasserstein on the UMAP results

            
        # if not os.path.exists(os.path.join(self.saving_folder, model, f"quantized_wasserstein_distance_{model}_{n_cluster}_clusters.csv")):
        df_w = self.invasive_image_embeddings[model].compute_quantized_wasserstein_distance_between_clusters(cluster_col='predicted_label', 
                                                                                                                    layer=None, 
                                                                                                                    ref_space=self.ref_model_emb)
            
        df_w.to_csv(os.path.join(self.saving_folder, model, f"quantized_wasserstein_distance_{model}_{n_cluster}_clusters.csv"))

        for patient in molecular_emb_labeled.emb.obs[self.group].unique():
            
            ## Computing wasserstein per patient in the image space
            
            print(f"Computing quantized wasserstein distance in the image space for patient {patient}...", flush=True)
            
            image_emb_patient = ImageEmbedding()
            image_emb_patient.emb = self.invasive_image_embeddings[model].emb[self.invasive_image_embeddings[model].emb.obs[self.group] == patient]
            
            if not os.path.exists(os.path.join(self.saving_folder, model, f"quantized_wasserstein_distance_image_{model}_{n_cluster}_clusters_patient_{patient}.csv")):
                df_w = image_emb_patient.compute_quantized_wasserstein_distance_between_clusters(cluster_col="predicted_label", 
                                                                                                    layer=None, 
                                                                                                    ref_space=None, 
                                                                                                    k=5000)
                
                df_w.to_csv(os.path.join(self.saving_folder, model, f"quantized_wasserstein_distance_image_{model}_{n_cluster}_clusters_patient_{patient}.csv"))
            
            
            ## Computing wasserstein per patient in the molecular space
            
            print(f"Computing quantized wasserstein distance for patient {patient}...", flush=True)
            
            molecular_emb_labeled_patient = GeneEmbedding()
            molecular_emb_labeled_patient.emb = molecular_emb_labeled.emb[molecular_emb_labeled.emb.obs[self.group] == patient]

            # Compute wasserstein distance in the molecular embedding
            if not os.path.exists(os.path.join(self.saving_folder, model, f"quantized_wasserstein_distance_molecular_{self.molecular_name}_{model}_{n_cluster}_clusters_patient_{patient}.csv")):
                df_w = molecular_emb_labeled_patient.compute_quantized_wasserstein_distance_between_clusters(cluster_col="predicted_label", 
                                                                                                        layer=None, 
                                                                                                        ref_space=None, 
                                                                                                        k=5000)
                
                df_w.to_csv(os.path.join(self.saving_folder, model, f"quantized_wasserstein_distance_molecular_{self.molecular_name}_{model}_{n_cluster}_clusters_patient_{patient}.csv"))
        
        
    
        

    def execute_pipeline(self):
        
        range_clusters = np.arange(self.min_cluster, self.max_cluster, self.cluster_step)

        self.saving_folder = os.path.join(self.saving_folder, self.algo)
        
        self.load_invasive_image_embeddings()
        self.initialize_model_saving_folders()


        for model in self.pipelines_list:
            print(f"Computing invasive cancer clustering for model {model}...", flush=True)
            print(self.invasive_image_embeddings[model].emb.X)
            

            self.invasive_image_embeddings[model].palette = self.dataset.PALETTE    
            self.invasive_image_embedding_raw_svd_umap(image_embedding=self.invasive_image_embeddings[model], 
                                                  saving_folder=self.saving_folder,
                                                  algo=self.algo,
                                                  model_name=model,
                                                  min_cluster=self.min_cluster,
                                                  max_cluster=self.max_cluster)
            
            files = glob.glob(os.path.join(self.saving_folder, model, "scores_umap_across_parameters_*_clusters.json"))
            
            files = [file for file in files if int(file.split("_clusters")[0].split("_")[-1]) in range_clusters]
            
            opti_cluster_number = get_optimal_cluster_number_one_model(files)
            opti_file = os.path.join(self.saving_folder, model, f"scores_umap_across_parameters_{opti_cluster_number}_clusters.json")

            self.sub_invasive_cancer_clustering_pipeline(opti_file, model)

    # @staticmethod
    # def sub_invasive_cancer_clustering_pipeline_multiprocess(args):
    #     benchmark_obj, file, model, ref_model = args
    #     benchmark_obj.sub_invasive_cancer_clustering_pipeline(file, model, ref_model=ref_model)



    @staticmethod
    def invasive_image_embedding_raw_svd_umap(image_embedding, 
                                              saving_folder, 
                                              model_name, 
                                              multiply_by_variance=True,
                                              algo='kmeans',
                                              min_cluster=3,
                                              max_cluster=10,
                                              svd_comp=5,
                                              cluster_step=1):

        clustering_files = sorted(glob.glob(os.path.join(saving_folder, model_name, "*.json")))

        print(f"Files found in {os.path.join(saving_folder, model_name)}: {clustering_files}", flush=True)



        basenames = [os.path.basename(f) for f in clustering_files]
        clusters = np.arange(min_cluster, max_cluster, cluster_step)

        print(f"Files found in {os.path.join(saving_folder, model_name)}: {basenames}", flush=True)


        if not f"scores_{clusters.min()}_{clusters.max()}_raw.json" in basenames:
            
            # raw
            print(f"scores_{clusters.min()}_{clusters.max()}_raw.json not found in {os.path.join(saving_folder, model_name)}")
            print(f"Computing {algo} clustering for raw data", flush=True)
            image_embedding.clustering_across_different_n_clusters(clusters_list=np.arange(min_cluster, max_cluster, cluster_step),
                                                        algo=algo,
                                                        layer=None,
                                                        saving_folder=os.path.join(saving_folder, model_name),)

        if not f"scores_{clusters.min()}_{clusters.max()}_svd_{svd_comp}.json" in basenames:
            # svd
            print(f"Computing {algo} clustering for svd data", flush=True)

            image_embedding.clustering_across_different_n_clusters(clusters_list=np.arange(min_cluster, max_cluster, cluster_step),
                                                        algo=algo,
                                                        layer='svd', 
                                                        denoised_svd=False, 
                                                        svd_component_number=svd_comp, 
                                                        saving_folder=os.path.join(saving_folder, model_name),
                                                        multiply_by_variance=multiply_by_variance)

        # umap
        for n_clust in clusters:
            try: 
                with open(os.path.join(saving_folder, model_name, f"scores_umap_across_parameters_{n_clust}_clusters.json"), 'r') as f:
                    results = json.load(f)
            except:
            # if not f"scores_umap_across_parameters_{n_clust}_clusters.json" in basenames:
                print(f"Computing {algo} clustering for umap with n={n_clust} clusters", flush=True)
                image_embedding.clustering_across_umap_parameters(n_neighbors_list=[10, 30, 50, 100, 150, 200, 250, 300, 350, 400],
                                                       min_dist_list=[0.001, 0.1],
                                                       n_components=2,
                                                       algo=algo, 
                                                       n_clusters=n_clust, 
                                                       saving_folder=os.path.join(saving_folder, model_name))
                
        best_umap_params = image_embedding.select_best_umap_parameters(saving_folder = os.path.join(saving_folder, model_name),
                                                            score='silhouette_score')
        



        print("Computing umap with best parameters", flush=True)
        print(f"optimal number of clusters {best_umap_params['n_clusters']}")
        image_embedding.compute_umap(n_neighbors=best_umap_params['n_neighbors'],
                          min_dist=best_umap_params['min_dist'])


        clustering_files = sorted(glob.glob(os.path.join(saving_folder, model_name, "*.json")))
        clustering_files = [file for file in clustering_files if f"{clusters.min()}_{clusters.max()}" in file]
        print(f"Clustering files: {clustering_files}")

        # Plot the scores individually (only for raw and svd; as plots are different for UMAPs)
        for clustering_file in clustering_files:
            if "svd" in clustering_file:
                add_filename = f"svd_{svd_comp}"
            elif "umap" in clustering_file:
                continue
            else:
                add_filename = "raw"
            image_embedding.set_unsupervised_clustering_score_files([clustering_file])

            print(image_embedding._unsupervised_clustering_score_files)
            image_embedding.plot_unsupervised_clustering_score(
                all_svd=False,
                algo=algo,
                suptitle=True,
                return_ax=True,
                add_best_cluster_vis=False,
                c=None,
            )
            plt.savefig(os.path.join(saving_folder, f"unsupervised_clustering_score_{model_name}_{add_filename}.png"))

        image_embedding.set_unsupervised_clustering_score_files(clustering_files)


        image_embedding.compute_optimal_number_of_clusters()

        for experiment in ["raw", f"svd_{svd_comp}", "umap"]:
        
            print(f"################# {experiment} #################", flush=True)
            

            nb_svd_component = 0
            if experiment == "raw":
                layer = None
            elif "svd" in experiment:
                layer = "svd"
                nb_svd_component = int(experiment.split("svd_")[-1])
            elif experiment == "umap":
                layer = "umap"
                experiment = f"umap_min_dist_{best_umap_params['min_dist']}_n_neighbors_{best_umap_params['n_neighbors']}"
            else:
                layer = None
                
            if layer != "umap":

                print(image_embedding.info_optimal_number_of_clusters_all_exp.columns)
                opt_k = int(image_embedding.info_optimal_number_of_clusters_all_exp.loc[experiment, "silhouette"])
                
            else:
                opt_k = best_umap_params['n_clusters']
            

            print("Optimal number of clusters: ", opt_k, flush=True)
            image_embedding.unsupervised_clustering(
            algo=algo,
            n_clusters=opt_k,
            u_comp_list=[f"u{i+1}" for i in range(nb_svd_component)],
            layer=layer
            )

            image_embedding.emb.obs[["label", "predicted_label"]].to_csv(os.path.join(saving_folder, model_name, f"invasive_labels_{opt_k}_clusters_among_{clusters.min()}_{clusters.max()}_{experiment}.csv"))


            image_embedding.emb.obs.replace("nan", np.NaN, inplace=True)
            image_embedding.barplots_predicted_clusters(
                groupby="predicted_label",
                threshold_invasive_dashed_line=False,
                palette_label=image_embedding.palette,
            )
            plt.savefig(os.path.join(saving_folder, f"barplot_predicted_clusters_{experiment}_model_{model_name}.png"))

            image_embedding.pie_charts_each_patient_across_selected_clusters(
                    selected_clusters_list=sorted(
                        image_embedding.emb.obs["predicted_label"][
                            ~image_embedding.emb.obs["predicted_label"].isna()
                        ]
                        .unique()
                        .tolist()
                    ),
                    label_column="predicted_label",
                )
            plt.savefig(os.path.join(saving_folder, f"piechart_predicted_clusters_{experiment}_model_{model_name}.png"))
            