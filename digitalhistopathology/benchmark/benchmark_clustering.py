#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#

import glob
import matplotlib.pyplot as plt
from digitalhistopathology.benchmark.benchmark_base import BenchmarkBase
import os
import json
import numpy as np
import seaborn as sns
from sklearn.metrics import adjusted_rand_score
from digitalhistopathology.embeddings.embedding import Embedding
from digitalhistopathology.datasets.real_datasets import HER2Dataset
from digitalhistopathology.benchmark.benchmark_utils import plot_ari_scores_all_patients



class BenchmarkClustering(BenchmarkBase):
    """
    """


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
                 algo='kmeans',
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

        self.algo = algo


    def unsupervised_clustering_benchmark(self,
                                          n_clusters=None, 
                                          obsm=None, 
                                          layer=None, 
                                          recompute_dim_red_on_subset=False, 
                                          var_ratio_threshold_for_svd=None,
                                          u_comp_list=None,
                                          center_before_svd=False,
                                          scale_before_svd=False,
                                          multiply_by_variance=False,
                                          clustering_name='raw',
                                          **kwargs):
        """ 
        Run unsupervised clustering on the image embeddings and handcrafted features.
        Parameters:
        -----------
        n_clusters : int, optional
            Number of clusters to use for the clustering algorithm. Default is None.
        Returns:
        --------
        unsupervised_clustering_results : dict
        
        """

        unsupervised_clustering_results = {}


        for model in self.pipelines_list:
            
            if os.path.exists(os.path.join(self.saving_folder, f'{clustering_name}_clustering_results_{model}.json')):
                print(f"Loading unsupervised clustering results for model {model}...", flush=True)
                with open(os.path.join(self.saving_folder, f'{clustering_name}_clustering_results_{model}.json'), 'r') as f:
                    unsupervised_clustering_results[model] = json.load(f)[model]
                                    
            else:
            

                print(f"Computing {clustering_name} clustering for model {model}, obsm {obsm}", flush=True)
                if "label" not in self.image_embeddings[model].emb.obs:
                    print(f"Adding label to the embeddings of {model}...", flush=True)
                    self.image_embeddings[model].add_label()
                # self.image_embeddings[model].emb.obs.replace("nan", np.nan, inplace=True)
                unsupervised_clustering_results[model] = self.image_embeddings[model].unsupervised_clustering_per_sample(n_clusters=n_clusters, 
                                                                                                                        sample='name_origin',
                                                                                                                        obsm=obsm, 
                                                                                                                        layer=layer, 
                                                                                                                        recompute_dim_red_on_subset=recompute_dim_red_on_subset,
                                                                                                                        var_ratio_threshold_for_svd=var_ratio_threshold_for_svd,
                                                                                                                        u_comp_list=u_comp_list,
                                                                                                                        center_before_svd=center_before_svd,
                                                                                                                        scale_before_svd=scale_before_svd,
                                                                                                                        algo=self.algo,
                                                                                                                        multiply_by_variance=multiply_by_variance,
                                                                                                                        filter_label=False,
                                                                                                                        **kwargs)
                
                
                unsupervised_clustering_results[model][f'ARI_{self.group}'] = adjusted_rand_score(self.image_embeddings[model].emb.obs['tumor'], 
                                                                                                self.image_embeddings[model].emb.obs['predicted_label'])
            
                with open(os.path.join(self.saving_folder, f'{clustering_name}_clustering_results_{model}.json'), 'w') as f:
                    json.dump({model: unsupervised_clustering_results[model]}, f, indent=4)
        # For handcrafted features
        if self.ef is not None:
            
            if os.path.exists(os.path.join(self.saving_folder, f'{clustering_name}_clustering_results_handcrafted_features.json')):
                print(f"Loading {clustering_name} clustering results for handcrafted features...", flush=True)
                with open(os.path.join(self.saving_folder, f'{clustering_name}_clustering_results_handcrafted_features.json'), 'r') as f:
                    unsupervised_clustering_results['handcrafted_features'] = json.load(f)['handcrafted_features']
            else:
                if "label" not in self.ef.emb.obs:
                    self.ef.add_label()


                unsupervised_clustering_results['handcrafted_features'] = self.ef.unsupervised_clustering_per_sample(n_clusters=n_clusters, 
                                                                                                                    sample='name_origin',
                                                                                                                    obsm=obsm, 
                                                                                                                    layer=layer, 
                                                                                                                    recompute_dim_red_on_subset=recompute_dim_red_on_subset,
                                                                                                                    var_ratio_threshold_for_svd=var_ratio_threshold_for_svd,
                                                                                                                    u_comp_list=u_comp_list,
                                                                                                                    center_before_svd=center_before_svd,
                                                                                                                    scale_before_svd=scale_before_svd,
                                                                                                                    algo=self.algo,
                                                                                                                    multiply_by_variance=multiply_by_variance,
                                                                                                                    filter_label=False,
                                                                                                                    **kwargs)
                

                
                unsupervised_clustering_results['handcrafted_features'][f'ARI_{self.group}'] = adjusted_rand_score(self.ef.emb.obs['tumor'],
                                                                                                                    self.ef.emb.obs['predicted_label'])
                
                with open(os.path.join(self.saving_folder, f'{clustering_name}_clustering_results_handcrafted_features.json'), 'w') as f:
                    json.dump({'handcrafted_features': unsupervised_clustering_results['handcrafted_features']}, f, indent=4)
            
        return unsupervised_clustering_results


    def set_best_UMAP_overall(self, annotated_only=False):

        if annotated_only:
            
            for model in self.annotated_embeddings.keys():
                print(f"Setting best UMAP parameters for all annotated images for model {model}", flush=True)
                if os.path.exists(os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_all_annotated_only.json')):
                    with open(os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_all_annotated_only.json'), 'r') as f:
                        best_ari_umap = json.load(f)
                else:
                    print(f"File {os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_all_annotated_only.json')} doesn't exist")
                    raise NotImplementedError(f"First compute the best UMAP parameters on annotated only for model {model}.")
            
                self.annotated_embeddings[model].compute_umap(n_neighbors=best_ari_umap['params']['n_neighbors'],
                                                            min_dist=best_ari_umap['params']['min_dist'])

        else:
            
            for model in self.image_embeddings.keys():
                print(f"Setting best UMAP parameters for all images for model {model}", flush=True)

                if os.path.exists(os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_all.json')):
                    with open(os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_all.json'), 'r') as f:
                        best_ari_umap = json.load(f)
                else:
                    raise NotImplementedError(f"First compute the best UMAP parameters for model {model}.")
                
                self.image_embeddings[model].compute_umap(n_neighbors=best_ari_umap['params']['n_neighbors'],
                                                            min_dist=best_ari_umap['params']['min_dist'])
            
            if self.ef is not None:

                if os.path.exists(os.path.join(self.saving_folder, f'best_umap_ari_model_handcrafted_features_all.json')):
                    with open(os.path.join(self.saving_folder, f'best_umap_ari_model_handcrafted_features_all.json'), 'r') as f:
                        best_ari_umap = json.load(f)
                else:
                    raise NotImplementedError(f"First compute the best UMAP parameters for handcrafted features.")
                
                self.ef.compute_umap(n_neighbors=best_ari_umap['params']['n_neighbors'],
                                    min_dist=best_ari_umap['params']['min_dist'])
            

            

    
    def set_best_UMAP_per_slide(self):
            
        embeddings_per_slide_umap = {}
        for model in self.embeddings_per_slide.keys():
            embeddings_per_slide_umap[model] = {}
            for slide, subset_emb in self.embeddings_per_slide[model].items():
                

                n_labeled_patches = subset_emb.emb.obs[(subset_emb.emb.obs['label'] != 'nan') & ~(subset_emb.emb.obs['label'].isna())].shape[0]
                if n_labeled_patches > 0:
                    embeddings_per_slide_umap[model][slide] = {}
                    subset_emb_umap = Embedding()
                    subset_emb_umap.emb = subset_emb.emb.copy()

                    if os.path.exists(os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_patient_{slide}.json')):
                        with open(os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_patient_{slide}.json'), 'r') as f:
                            best_ari_umap = json.load(f)
                    else:
                        # best_ari_umap = self.get_best_UMAP_ari_subset(subset_emb, 
                        #                                                 model, 
                        #                                                 patient, 
                        #                                                 algo='kmeans', 
                        #                                                 n_neighbors_list=[5, 10, 20, 30, 40, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300], 
                        #                                                 min_dist_list=[0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1], 
                        #                                                 n_components_list=[2])
                        raise NotImplementedError(f"First compute the best UMAP parameters per patient for slide {slide} and model {model}.")
                    
                    subset_emb_umap.compute_umap(n_neighbors=best_ari_umap['params']['n_neighbors'], 
                                            min_dist=best_ari_umap['params']['min_dist'])
                    
                    embeddings_per_slide_umap[model][slide] = subset_emb_umap
                else:
                    embeddings_per_slide_umap[model][slide] = subset_emb

        self.embeddings_per_slide = embeddings_per_slide_umap




    def UMAP_best_params_clustering_visualization_overall(self, 
                                                          hue='predicted_label', 
                                                          annotated_only=False, 
                                                          palette=None):

        add_title = ""

        if annotated_only:
            embeddings = list(self.annotated_embeddings.values())
            models = list(self.annotated_embeddings.keys())
            add_filename = "_annotated_only"

        else:
            embeddings = list(self.image_embeddings.values())
            models = list(self.image_embeddings.keys())
            add_filename = ""
            
            if self.ef is not None:
                embeddings.append(self.ef)
                models.append('handcrafted_features')


        for model, embedding in zip(models, embeddings):
            print(f"Plotting UMAP for model {model}...", flush=True)
            if hue == 'label':
                embedding.emb.obs['label'] = embedding.emb.obs['label'].replace(np.nan, 'nan')

            if hue == 'predicted_label':
                n_tissues = len([label for label in embedding.emb.obs['label'].unique() if label not in ['nan', 'undetermined', np.nan]])
                
                # Note: here, for HER2, there was an error and n_clusters was at 6 instead of 7.
                r = embedding.unsupervised_clustering(n_clusters=n_tissues, obsm='umap', algo=self.algo)
                ARI = f"{r['ari']:.2f}"
                add_title = f" - ARI: {ARI}"
                
                
            if hue == 'tumor':
                # number of patients
                n_patients = embedding.emb.obs['tumor'].nunique()

                if n_patients <= 10:
                    r = embedding.unsupervised_clustering(n_clusters=n_patients, obsm='umap', algo=self.algo)
                    ARI_patient = adjusted_rand_score(embedding.emb.obs['tumor'], embedding.emb.obs['predicted_label'])
                else:
                    aris_bootstrapped = []
                    if n_patients > 10:
                        r = embedding.unsupervised_clustering(n_clusters=10, obsm='umap', algo=self.algo)
                        for i in range(20):
                            random_patients = np.random.choice(embedding.emb.obs['tumor'].unique(), size=10, replace=False)
                            sub_emb = embedding.emb[embedding.emb.obs['tumor'].isin(random_patients)]
                            ari = adjusted_rand_score(sub_emb.obs['tumor'], sub_emb.obs['predicted_label'])
                            aris_bootstrapped.append(ari)
                        ARI_patient = np.mean(aris_bootstrapped)
                        
                add_title = f" - patient ARI (batch): {ARI_patient:.2f}"
                

            
            plt.figure(figsize=(10,10))
            sns.scatterplot(x=embedding.emb.obsm['umap'][:,0],
                            y=embedding.emb.obsm['umap'][:,1],
                            hue=embedding.emb.obs[hue],
                            palette=palette)
            sns.despine()
            plt.title(f"UMAP - {model} - colored by {hue} -{add_filename.replace('_', ' ')}{add_title}", weight='bold')
            plt.savefig(os.path.join(self.saving_folder, f"UMAP_best_params_{model}{add_filename}_colored_by_{hue}.{self.extension}"), bbox_inches='tight')
            
            plt.figure(figsize=(10,10))
            sns.kdeplot(x=embedding.emb.obsm['umap'][:,0],
                            y=embedding.emb.obsm['umap'][:,1],
                            hue=embedding.emb.obs[hue],
                            palette=palette)
            sns.despine()
            plt.title(f"UMAP - {model} - colored by {hue} -{add_filename.replace('_', ' ')}{add_title}", weight='bold')
            plt.savefig(os.path.join(self.saving_folder, f"UMAP_kde_best_params_{model}{add_filename}_colored_by_{hue}.{self.extension}"), bbox_inches='tight')


    def UMAP_best_params_clustering_many_patients_batch_old(self, 
                                                          annotated_only=False, 
                                                          palette=None,
                                                          n_patients=10):

        add_title = ""

        if annotated_only:
            embeddings = list(self.annotated_embeddings.values())
            models = list(self.annotated_embeddings.keys())
            add_filename = "_annotated_only"

        else:
            embeddings = list(self.image_embeddings.values())
            models = list(self.image_embeddings.keys())
            add_filename = ""
            
            if self.ef is not None:
                embeddings.append(self.ef)
                models.append('handcrafted_features')

        boostraped_ARI_patients = {}
        for model, embedding in zip(models, embeddings):
            
            if os.path.exists(os.path.join(self.saving_folder, f'boostraped_ARI_patient_{model}_{n_patients}_patients{add_filename}.json')):
                print(f"Loading boostraped ARI for model {model}...", flush=True)
                with open(os.path.join(self.saving_folder, f'boostraped_ARI_patient_{model}_{n_patients}_patients{add_filename}.json'), 'r') as f:
                    boostraped_ARI_patients[model] = json.load(f)[model]
            else:
                all_patients = embedding.emb.obs['tumor'].unique()

                boostraped_ARI_patients[model] = []
                for i in range(20):
                    patients = np.random.choice(all_patients, size=n_patients, replace=False)
                    subset_emb = Embedding()
                    subset_emb.emb = embedding.emb[embedding.emb.obs['tumor'].isin(patients)].copy()
                    r = subset_emb.unsupervised_clustering(n_clusters=n_patients, obsm='umap', algo=self.algo)
                    ARI_patient = adjusted_rand_score(subset_emb.emb.obs['tumor'], subset_emb.emb.obs['predicted_label'])
                    add_title = f" - patient ARI (batch): {ARI_patient:.2f} - random {len(patients)} patients - round {i+1}"
                    boostraped_ARI_patients[model].append(ARI_patient)
                
                    plt.figure(figsize=(10,10))
                    sns.scatterplot(x=subset_emb.emb.obsm['umap'][:,0],
                                    y=subset_emb.emb.obsm['umap'][:,1],
                                    hue=subset_emb.emb.obs['tumor'],
                                    palette=palette)
                    sns.despine()
                    plt.title(f"UMAP - {model} - colored by tumor -{add_filename.replace('_', ' ')}{add_title}", weight='bold')
                    plt.savefig(os.path.join(self.saving_folder, f"UMAP_best_params_{model}{add_filename}_random_{n_patients}_{i}_colored_by_tumor.{self.extension}"), bbox_inches='tight')
                    
                    plt.figure(figsize=(10,10))
                    sns.kdeplot(x=subset_emb.emb.obsm['umap'][:,0],
                                    y=subset_emb.emb.obsm['umap'][:,1],
                                    hue=subset_emb.emb.obs['tumor'],
                                    palette=palette)
                    sns.despine()
                    plt.title(f"UMAP - {model} - colored by tumor -{add_filename.replace('_', ' ')}{add_title}", weight='bold')
                    plt.savefig(os.path.join(self.saving_folder, f"UMAP_kde_best_params_{model}{add_filename}_random_{n_patients}_{i}_colored_by_tumor.{self.extension}"), bbox_inches='tight')
                
                    plt.figure(figsize=(10,10))
                    subset_emb.emb = subset_emb.emb[~subset_emb.emb.obs['label'].isna()]
                    subset_emb.emb = subset_emb.emb[subset_emb.emb.obs['label'] != 'nan']
                    sns.scatterplot(x=subset_emb.emb.obsm['umap'][:,0],
                                    y=subset_emb.emb.obsm['umap'][:,1],
                                    hue=subset_emb.emb.obs['label'],
                                    palette=self.dataset.PALETTE)
                    sns.despine()
                    plt.title(f"UMAP - {model} - colored by label -{add_filename.replace('_', ' ')}{add_title}", weight='bold')
                    plt.savefig(os.path.join(self.saving_folder, f"UMAP_best_params_{model}{add_filename}_random_{n_patients}_{i}_colored_by_label.{self.extension}"), bbox_inches='tight')

                    plt.figure(figsize=(10,10))
                    sns.kdeplot(x=subset_emb.emb.obsm['umap'][:,0],
                                    y=subset_emb.emb.obsm['umap'][:,1],
                                    hue=subset_emb.emb.obs['label'],
                                    palette=self.dataset.PALETTE)
                    sns.despine()
                    plt.title(f"UMAP - {model} - colored by label -{add_filename.replace('_', ' ')}{add_title}", weight='bold')
                    plt.savefig(os.path.join(self.saving_folder, f"UMAP_kde_best_params_{model}{add_filename}_random_{n_patients}_{i}_colored_by_label.{self.extension}"), bbox_inches='tight')
                
                with open(os.path.join(self.saving_folder, f'boostraped_ARI_patient_{model}_{n_patients}_patients{add_filename}.json'), 'w') as f:
                    json.dump({model: boostraped_ARI_patients[model]}, f, indent=4)
                
        with open(os.path.join(self.saving_folder, f'boostraped_ARI_patient_{n_patients}_patients{add_filename}.json'), 'w') as f:
            json.dump(boostraped_ARI_patients, f, indent=4)

    def UMAP_best_params_clustering_many_patients_batch(self, 
                                                          annotated_only=False, 
                                                          palette=None,
                                                          n_patients=10):

        add_title = ""

        if annotated_only:
            embeddings = list(self.annotated_embeddings.values())
            models = list(self.annotated_embeddings.keys())
            add_filename = "_annotated_only"

        else:
            embeddings = list(self.image_embeddings.values())
            models = list(self.image_embeddings.keys())
            add_filename = ""
            
            if self.ef is not None:
                embeddings.append(self.ef)
                models.append('handcrafted_features')

        boostraped_ARI_patients = {}
        for model, embedding in zip(models, embeddings):
            
            # if os.path.exists(os.path.join(self.saving_folder, f'boostraped_ARI_patient_{model}_{n_patients}_patients{add_filename}.json')):
            #     print(f"Loading boostraped ARI for model {model}...", flush=True)
            #     with open(os.path.join(self.saving_folder, f'boostraped_ARI_patient_{model}_{n_patients}_patients{add_filename}.json'), 'r') as f:
            #         boostraped_ARI_patients[model] = json.load(f)[model]
            # else:
            all_patients = embedding.emb.obs['tumor'].unique()
            embedding.unsupervised_clustering(n_clusters=n_patients, obsm='umap', algo=self.algo)
            boostraped_ARI_patients[model] = []
            for i in range(20):
                patients = np.random.choice(all_patients, size=n_patients, replace=False)
                subset_emb = Embedding()
                subset_emb.emb = embedding.emb[embedding.emb.obs['tumor'].isin(patients)].copy()
                # r = subset_emb.unsupervised_clustering(n_clusters=n_patients, obsm='umap', algo=self.algo)
                ARI_patient = adjusted_rand_score(subset_emb.emb.obs['tumor'], subset_emb.emb.obs['predicted_label'])
                add_title = f" - patient ARI (batch): {ARI_patient:.2f} - random {len(patients)} patients - round {i+1}"
                boostraped_ARI_patients[model].append(ARI_patient)
            
                plt.figure(figsize=(10,10))
                sns.scatterplot(x=subset_emb.emb.obsm['umap'][:,0],
                                y=subset_emb.emb.obsm['umap'][:,1],
                                hue=subset_emb.emb.obs['tumor'],
                                palette=palette)
                sns.despine()
                plt.title(f"UMAP - {model} - colored by tumor -{add_filename.replace('_', ' ')}{add_title}", weight='bold')
                plt.savefig(os.path.join(self.saving_folder, f"UMAP_best_params_{model}{add_filename}_random_{n_patients}_{i}_colored_by_tumor.{self.extension}"), bbox_inches='tight')
                
                plt.figure(figsize=(10,10))
                sns.kdeplot(x=subset_emb.emb.obsm['umap'][:,0],
                                y=subset_emb.emb.obsm['umap'][:,1],
                                hue=subset_emb.emb.obs['tumor'],
                                palette=palette)
                sns.despine()
                plt.title(f"UMAP - {model} - colored by tumor -{add_filename.replace('_', ' ')}{add_title}", weight='bold')
                plt.savefig(os.path.join(self.saving_folder, f"UMAP_kde_best_params_{model}{add_filename}_random_{n_patients}_{i}_colored_by_tumor.{self.extension}"), bbox_inches='tight')
            
                plt.figure(figsize=(10,10))
                subset_emb.emb = subset_emb.emb[~subset_emb.emb.obs['label'].isna()]
                subset_emb.emb = subset_emb.emb[subset_emb.emb.obs['label'] != 'nan']
                sns.scatterplot(x=subset_emb.emb.obsm['umap'][:,0],
                                y=subset_emb.emb.obsm['umap'][:,1],
                                hue=subset_emb.emb.obs['label'],
                                palette=self.dataset.PALETTE)
                sns.despine()
                plt.title(f"UMAP - {model} - colored by label -{add_filename.replace('_', ' ')}{add_title}", weight='bold')
                plt.savefig(os.path.join(self.saving_folder, f"UMAP_best_params_{model}{add_filename}_random_{n_patients}_{i}_colored_by_label.{self.extension}"), bbox_inches='tight')

                plt.figure(figsize=(10,10))
                sns.kdeplot(x=subset_emb.emb.obsm['umap'][:,0],
                                y=subset_emb.emb.obsm['umap'][:,1],
                                hue=subset_emb.emb.obs['label'],
                                palette=self.dataset.PALETTE)
                sns.despine()
                plt.title(f"UMAP - {model} - colored by label -{add_filename.replace('_', ' ')}{add_title}", weight='bold')
                plt.savefig(os.path.join(self.saving_folder, f"UMAP_kde_best_params_{model}{add_filename}_random_{n_patients}_{i}_colored_by_label.{self.extension}"), bbox_inches='tight')
            
            with open(os.path.join(self.saving_folder, f'boostraped_ARI_patient_{model}_{n_patients}_patients{add_filename}.json'), 'w') as f:
                json.dump({model: boostraped_ARI_patients[model]}, f, indent=4)
                
        with open(os.path.join(self.saving_folder, f'boostraped_ARI_patient_{n_patients}_patients{add_filename}.json'), 'w') as f:
            json.dump(boostraped_ARI_patients, f, indent=4)

    def UMAP_best_params_clustering_visualization_per_slide(self, slide_id_col='name_origin'):
        if self.embeddings_per_slide is None:
            self.get_embeddings_per_slide()

        
        if self.ef is not None:
            for slide in self.ef.emb.obs[slide_id_col].unique():
                plt.figure(figsize=(6, 2*len(self.pipelines_list)))
                for i, model in enumerate(self.pipelines_list):

                    subset_emb = Embedding()
                    subset_emb.emb = self.embeddings_per_slide[model][slide].emb.copy()

                    nc = list(subset_emb.emb.obs["label"].unique())
                    nc = [e for e in nc if isinstance(e, str) and e != "undetermined" and e != "nan" and e != np.nan]
                    print(nc)
                    nc = len(nc)
                    print("Number of clusters = {}".format(nc))

                    if nc > 1:

                        r = subset_emb.unsupervised_clustering(n_clusters=nc, obsm='umap', algo=self.algo)

                        # Remove nan labels for the plot
                        subset_emb.emb = subset_emb.emb[~subset_emb.emb.obs['label'].isna()]
                        subset_emb.emb = subset_emb.emb[subset_emb.emb.obs['label'] != 'nan']

                        plt.subplot(len(self.pipelines_list), 2, (list(self.pipelines_list).index(model))*2+1)
                        sns.scatterplot(x=subset_emb.emb.obsm['umap'][:,0],
                                        y=subset_emb.emb.obsm['umap'][:,1],
                                        hue=subset_emb.emb.obs['label'],
                                        palette=HER2Dataset.PALETTE, 
                                        s=10)   
                        plt.title(f"{model}", weight='bold')
                        sns.despine()
                        if i != len(self.pipelines_list)-1:
                            plt.legend().remove()
                        else:
                            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

                        plt.subplot(len(self.pipelines_list), 2, (list(self.pipelines_list).index(model))*2+2)


                        sns.scatterplot(x=subset_emb.emb.obsm['umap'][:,0],
                                        y=subset_emb.emb.obsm['umap'][:,1],
                                        hue=subset_emb.emb.obs['predicted_label'],
                                        palette='Accent',
                                        s=10)
                        plt.title(f"ARI: {r['ari']:.2f}", weight='bold')
                        sns.despine()
                        if i != len(self.pipelines_list)-1:
                            plt.legend().remove()
                        else:
                            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

                    plt.tight_layout()
                    plt.suptitle(f"patient {slide}", weight='bold')
                    #plt.suptitle(f"UMAP kmeans - patient {patient}", weight='bold')
                    plt.savefig(os.path.join(self.saving_folder, f"UMAP_kmeans_best_params_patient_{slide}.{self.extension}"), bbox_inches='tight')



    def get_best_UMAP_ari_per_slide(self, 
                                    n_neighbors_list=[10, 30, 50, 100, 150, 200, 250, 300, 350, 400],
                                    min_dist_list=[0.001, 0.1], 
                                    n_components_list=[2],
                                    force_loading=False):

        best_ari_umap = {}
        print(f"Embeddings per patient: {self.embeddings_per_slide}", flush=True)
        # For each model
        # for model in self.embeddings_per_slide.keys():
        for model in list(self.embeddings_per_slide.keys()):
            print(f"Model: {model}", flush=True)
            best_ari_umap[model] = {}
            if force_loading:
                print(f"Force loading best UMAP parameters for model {model} for existing patients...", flush=True)
                for file in glob.glob(os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_patient_*.json')):
                    patient_name = os.path.basename(file).split("_patient_")[1].split('.json')[0]
                    with open(file, 'r') as f:
                        best_ari_umap[model][patient_name] = json.load(f)

            # elif self.embeddings_per_slide is not None:
            #     if len(glob.glob(os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_patient_*.json'))) == len(self.embeddings_per_slide[model]):
            #         print(f"Loading best UMAP parameters for model {model}...", flush=True)
            #         for file in glob.glob(os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_patient_*.json')):
            #             with open(file, 'r') as f:
            #                 best_ari_umap[model][os.path.basename(file).split('_')[6]] = json.load(f)
                
            else:
                for slide, subset_emb in self.embeddings_per_slide[model].items():

                    n_labeled_patches = subset_emb.emb.obs[(subset_emb.emb.obs['label'] != 'nan') & ~(subset_emb.emb.obs['label'].isna())].shape[0]
                    if n_labeled_patches > 0:

                        print(f"Computing best UMAP parameters for patient {slide} ({n_labeled_patches} labeled patches) with model {model}...", flush=True)

                        if os.path.exists(os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_patient_{slide}.json')):
                            print(f"Loading best UMAP parameters for model {model} and patient {slide}...", flush=True)
                            with open(os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_patient_{slide}.json'), 'r') as f:
                                best_ari_umap[model][slide] = json.load(f)

                        else:
                            print(f"File {os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_patient_{slide}.json')} doesn't exist", flush=True)
                            best_ari_umap[model][slide] = self.get_best_UMAP_ari_subset(subset_emb, 
                                                                                        model, 
                                                                                        slide, 
                                                                                        n_neighbors_list=n_neighbors_list, 
                                                                                        min_dist_list=min_dist_list, 
                                                                                        n_components_list=n_components_list)
                            
                            with open(os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_patient_{slide}.json'), 'w') as f:
                                json.dump(best_ari_umap[model][slide], f, indent=4)

        return best_ari_umap
    
    
    

    def get_best_UMAP_ari_subset(self, 
                                 subset_emb,
                                 model_name, 
                                 patient_name,
                                n_neighbors_list=[10, 30, 50, 100, 150, 200, 250, 300, 350, 400],
                                min_dist_list=[0.001, 0.1], 
                                 n_components_list=[2]):


        nc = subset_emb.emb.obs["label"].unique().tolist()
        nc = [e for e in nc if isinstance(e, str) and e != "undetermined" and e != "nan" and e != np.nan]
        print(nc)
        nc = len(nc)
        print(f"Number of clusters = {nc}")

        best_ari, best_params, df_all_runs = subset_emb.UMAP_validation_unsupervised_clustering(n_clusters=nc,
                                                                                    algo=self.algo,
                                                                                    n_neighbors_list=n_neighbors_list,
                                                                                    min_dist_list=min_dist_list,
                                                                                    n_components_list=n_components_list)

        best_ari_umap = {'ari': best_ari, 'params': best_params}
        
        df_all_runs.to_csv(os.path.join(self.saving_folder, f'UMAP_validation_{model_name}_{patient_name}.csv'))


        with open(os.path.join(self.saving_folder, f'best_umap_ari_model_{model_name}_patient_{patient_name}.json'), 'w') as f:
            json.dump(best_ari_umap, f)


        return best_ari_umap
    

    def compute_umap_ARI_patient_for_batch_control(self):
        if self.image_embeddings is not None and len(self.image_embeddings) > 0:
            n_patients = self.image_embeddings[self.pipelines_list[0]].emb.obs[self.group].nunique()
        else:
            n_patients = 0
        aris_patient = {}
        for model in self.pipelines_list:
            
            if os.path.exists(os.path.join(self.saving_folder, f'ARI_patient_{model}.json')):
                print(f"Loading ARI for model {model}...", flush=True)
                with open(os.path.join(self.saving_folder, f'ARI_patient_{model}.json'), 'r') as f:
                    aris_patient[model] = json.load(f)[model]
            else:
                # recompute unsupervised_clustering overall
                self.image_embeddings[model].unsupervised_clustering(n_clusters=n_patients, 
                                                                    algo=self.algo,
                                                                    obsm='umap')
                
                # compute ARI with patient
                aris_patient[model] = adjusted_rand_score(self.image_embeddings[model].emb.obs[self.group], self.image_embeddings[model].emb.obs['predicted_label'])
                
                with open(os.path.join(self.saving_folder, f'ARI_patient_{model}.json'), 'w') as f:
                    json.dump({model: aris_patient[model]}, f)   
                 
        if self.ef is not None:
            
            if os.path.exists(os.path.join(self.saving_folder, f'ARI_patient_handcrafted_features.json')):
                print(f"Loading ARI for handcrafted features...", flush=True)
                with open(os.path.join(self.saving_folder, f'ARI_patient_handcrafted_features.json'), 'r') as f:
                    aris_patient['handcrafted_features'] = json.load(f)['handcrafted_features']
            else:
            # handcrafed features
                
                self.ef.unsupervised_clustering(n_clusters=6,
                                                algo=self.algo,
                                                obsm='umap')
                

                
                aris_patient['handcrafted_features'] = adjusted_rand_score(self.ef.emb.obs[self.group], self.ef.emb.obs['predicted_label'])
                
                with open(os.path.join(self.saving_folder, f'ARI_patient_{model}.json'), 'w') as f:
                    json.dump({'handcrafted_features': aris_patient['handcrafted_features']}, f)   

        with open(os.path.join(self.saving_folder, 'ARI_patient.json'), 'w') as f:
            json.dump(aris_patient, f)
                
        return aris_patient

    
    def get_best_UMAP_ari_overall(self, n_clusters=6, 
                                      n_neighbors_list=[10, 30, 50, 100, 150, 200, 250, 300, 350, 400],
                                      min_dist_list=[0.001, 0.1], 
                                 n_components_list=[2],
                                 annotated_only=False,
                                 slides_id_col='name_origin'):
        
        models = self.pipelines_list.copy()
        embeddings = self.pipelines_list.copy()
        if annotated_only:
            # embeddings = list(self.annotated_embeddings.values())
            # models = list(self.annotated_embeddings.keys())
            if self.annotated_embeddings is not None and len(self.annotated_embeddings) > 0:
                embeddings = [self.annotated_embeddings[model] for model in models if self.annotated_embeddings is not None]
            add_filename = "_annotated_only"
        else:
            # embeddings = list(self.image_embeddings.values())
            # models = list(self.image_embeddings.keys())
            if self.image_embeddings is not None and len(self.image_embeddings) > 0:
                embeddings = [self.image_embeddings[model] for model in models if self.image_embeddings[model] is not None]
            add_filename = ""

            if self.ef is not None:
                embeddings.append(self.ef)
                models.append('handcrafted_features')


        best_ari_umap = {}
        for model, embedding in zip(models, embeddings):
            best_ari_umap[model] = {}

            print(f"Computing best UMAP parameters for model {model}...", flush=True)

            if os.path.exists(os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_all{add_filename}.json')):

                print(f"Loading best UMAP parameters for model {model}...", flush=True)
                with open(os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_all{add_filename}.json'), 'r') as f:
                    best_ari_umap[model] = json.load(f)
            else:

                print(f"File {os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_all{add_filename}.json')} doesn't exist", flush=True)

                best_ari, best_params, df_all_runs = embedding.UMAP_validation_unsupervised_clustering_non_parallel(n_clusters=n_clusters,
                                                                                          n_neighbors_list=n_neighbors_list,
                                                                                          min_dist_list=min_dist_list,
                                                                                          n_components_list=n_components_list,
                                                                                          algo=self.algo)
                
                df_all_runs.to_csv(os.path.join(self.saving_folder, f'UMAP_validation_{model}_all{add_filename}.csv'))

                best_ari_umap[model]['ari'] = best_ari
                best_ari_umap[model]['params'] = best_params

                with open(os.path.join(self.saving_folder, f'best_umap_ari_model_{model}_all{add_filename}.json'), 'w') as f:
                    json.dump(best_ari_umap[model], f)

        return best_ari_umap
    


    def execute_pipeline(self):
        # Compute unsupervised clustering on raw data
        self.saving_folder = os.path.join(self.saving_folder, self.algo)

        if not os.path.exists(self.saving_folder):
            os.makedirs(self.saving_folder, exist_ok=True)

        print(f"Saving folder: {self.saving_folder}", flush=True)

        if self.ef is not None:
            self.ef.add_label(dataset=self.dataset_name)


            
        unsupervised_clustering_results_optk = self.unsupervised_clustering_benchmark(n_clusters=None, clustering_name='raw')

        
        plt.figure()
        plot_ari_scores_all_patients(unsupervised_clustering_results_optk)
        plt.savefig(os.path.join(self.saving_folder, f'ari_scores_optk.{self.extension}'), bbox_inches='tight')




        svd5_multiplied_by_S = self.unsupervised_clustering_benchmark(n_clusters=None, 
                                                                             layer='svd', 
                                                                             var_ratio_threshold_for_svd=None,
                                                                               u_comp_list=[f"u{i+1}" for i in range(5)],
                                                                               recompute_dim_red_on_subset=True,
                                                                               multiply_by_variance=True,
                                                                               center_before_svd=False,
                                                                               scale_before_svd=False,
                                                                               clustering_name='svd5')
            

    
        plt.figure()
        plot_ari_scores_all_patients(svd5_multiplied_by_S)
        plt.savefig(os.path.join(self.saving_folder, f'ari_scores_svd5.{self.extension}'), bbox_inches='tight')



        self.get_embeddings_per_slide()
        print(f"Embeddings per slides loaded for models: {self.embeddings_per_slide.keys()}")
        best_umap_ari_per_slide = self.get_best_UMAP_ari_per_slide()
        with open(os.path.join(self.saving_folder, 'best_umap_ari_per_slide.json'), 'w') as f:
            json.dump(best_umap_ari_per_slide, f)

        best_umap_ari_overall = self.get_best_UMAP_ari_overall()
        with open(os.path.join(self.saving_folder, 'best_umap_ari_overall.json'), 'w') as f:
            json.dump(best_umap_ari_overall, f)

        self.get_annotated_embeddings()
        best_umap_ari_overall_annotated_only = self.get_best_UMAP_ari_overall(annotated_only=True)
        with open(os.path.join(self.saving_folder, 'best_umap_ari_overall_annotated_only.json'), 'w') as f:
            json.dump(best_umap_ari_overall_annotated_only, f)

        # Visualization UMAP-kmeans per slide
        self.set_best_UMAP_per_slide()
        # self.UMAP_best_params_clustering_visualization_per_slide()
        self.plot_slides_with_clusters(obsm='umap')


        
        self.set_best_UMAP_overall()

        self.UMAP_best_params_clustering_visualization_overall(hue='predicted_label', palette='Set3')
        self.UMAP_best_params_clustering_visualization_overall(hue='tumor', palette='Accent')
        self.UMAP_best_params_clustering_visualization_overall(hue='name_origin')

        palette_with_nans = self.dataset.PALETTE.copy()
        palette_with_nans['nan'] = 'black'
        self.UMAP_best_params_clustering_visualization_overall(hue='label', palette=palette_with_nans)

        # Visiualization UMAP k-means on annotated slides only
        self.set_best_UMAP_overall(annotated_only=True)
        self.UMAP_best_params_clustering_visualization_overall(hue='predicted_label', annotated_only=True, palette='Set3')
        self.UMAP_best_params_clustering_visualization_overall(hue='tumor', annotated_only=True, palette='Accent')
        self.UMAP_best_params_clustering_visualization_overall(hue='name_origin', annotated_only=True, palette='Accent')
        self.UMAP_best_params_clustering_visualization_overall(hue='label', annotated_only=True, palette=HER2Dataset.PALETTE)
        
        # Visualization UMAP k-means overall all slides
        n_patients = self.image_embeddings[self.pipelines_list[0]].emb.obs[self.group].nunique()
        if n_patients <= 10:
            ari_patients = self.compute_umap_ARI_patient_for_batch_control()
        else:
            self.UMAP_best_params_clustering_many_patients_batch(annotated_only=False, palette='Set3', n_patients=10)

        # Plot ARI scores
        plt.figure()
        plot_ari_scores_all_patients(best_umap_ari_per_slide)
        plt.savefig(os.path.join(self.saving_folder, f'ari_scores_best_umap.{self.extension}'), bbox_inches='tight')


    def plot_slides_with_clusters(self, obsm=None, layer=None):

        if obsm is not None:
            layer_name = obsm
        elif layer is not None:
            layer_name = layer
        else:
            layer_name = 'raw'
                
        if self.embeddings_per_slide is None:
            self.get_embeddings_per_slide()

        # Count number of annotated slides
        n_slides = len(self.embeddings_per_slide[self.pipelines_list[0]])

        if n_slides <= 10:
            emb_slides = self.embeddings_per_slide.copy()

        else:
            emb_slides = {}
            for model in self.pipelines_list:
                emb_slides[model] = {}
                slides = list(self.embeddings_per_slide[model].keys())
                selected_slides = list(self.embeddings_per_slide[self.pipelines_list[0]].keys())[:10]
                for slide in selected_slides:
                    emb_slides[model][slide] = self.embeddings_per_slide[model][slide]

            print(f"Number of slides > 10, randomly selecting 10 slides: {list(emb_slides[self.pipelines_list[0]].keys())}", flush=True)
        
        for model in self.pipelines_list:

            for slide, ie_subset_img in emb_slides[model].items():
                if ie_subset_img.emb.obs[(ie_subset_img.emb.obs['label'] != 'nan') & ~(ie_subset_img.emb.obs['label'].isna())].shape[0] > 0:
                    # ie_subset_img.emb.obs['path_origin'] = ie_subset_img.emb.obs['path_origin'].apply(lambda x: x.replace('ghaefliger', 'lfournier/repositories'))
                    ie_subset_img.emb.obs['start_width_origin'] = ie_subset_img.emb.obs['start_width_origin'].astype(int)
                    ie_subset_img.emb.obs['start_height_origin'] = ie_subset_img.emb.obs['start_height_origin'].astype(int)



                    nc = list(ie_subset_img.emb.obs["label"].unique())
                    nc = [e for e in nc if isinstance(e, str) and e != "undetermined" and e != "nan" and e != np.nan]
                    print(nc)
                    nc = len(nc)
                    print("Number of clusters = {}".format(nc))

                    if nc > 1:

                        r = ie_subset_img.unsupervised_clustering(n_clusters=nc, obsm=obsm, algo=self.algo, layer=layer)
                        
                        # annotated_patient = ie_subset_img.emb.obs['name_origin'].unique()[0] # find the name of the annotated sample
                        subset_without_nans = Embedding()
                        subset_without_nans.emb = ie_subset_img.emb[~ie_subset_img.emb.obs['label'].isna()]
                        subset_without_nans.emb = subset_without_nans.emb[subset_without_nans.emb.obs['label'] != 'nan']
                        
                        subset_without_nans.compare_two_labels_plot(sample_name=slide,
                                                            label_obs_name1="label", 
                                                            label_obs_name2="predicted_label",
                                                            palette_1=HER2Dataset.PALETTE,
                                                            )
                        


                        plt.suptitle(f"Comparison of labels for {slide} - {model} - ARI: {r['ari']:.2f}", weight='bold')
                        plt.savefig(os.path.join(self.saving_folder, f"comparison_labels_{layer_name}_{self.algo}_{slide}_{model}.{self.extension}"), bbox_inches='tight')

