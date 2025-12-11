#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#

import os
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
from scipy.cluster.hierarchy import dendrogram, fcluster
import json

import glob
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from digitalhistopathology.benchmark.benchmark_base import BenchmarkBase

from matplotlib.colors import LinearSegmentedColormap

colors_heatmaps = [(1, 0, 1), (0, 0, 0), (0, 1, 1)]  # RGB for magenta, black, cyan
custom_diverging_map = LinearSegmentedColormap.from_list('custom_diverging', colors_heatmaps, N=256)




class BenchmarkCorrelation(BenchmarkBase):
    """
    A class to benchmark image embeddings and engineered features.
    Attributes:
    -----------
    path_to_pipeline : str or list of str
        Path to the pipeline directory.
    pipelines_list : list
        List of pipeline models.
    dataset_name : str
        Folder containing the dataset.
    patches_filenames : list
        List of patch filenames.
    patches_info_filename : str
        Filename containing patch information.
    results_folder : str
        Folder to save results.
    name_engineered_features : str
        Name of the engineered features.
    saving_folder : str
        Folder to save intermediate results.
    image_embedding_name : str
        Name of the image embedding file.
    image_embeddings : dict
        Dictionary to store image embeddings.
    ef : EngineeredFeatures
        Instance of EngineeredFeatures class.
    Methods:
    --------
    __init__(self, path_to_pipeline=None, pipelines_list=None, dataset_name=None, patches_filenames=None, patches_info_filename=None, results_folder="../results", name_engineered_features=None, saving_folder=None, image_embedding_name=None):
        Initializes the Benchmark class with the given parameters.
    compute_image_embeddings(self, denoised=False):
        Computes image embeddings for each model in the pipeline list.
    compute_correlation_btw_svd_components_and_one_engineered_feature(self, pcs=[0], denoised=False, engineered_feature_name='nuclei_number'):
        Computes the correlation between SVD components and a single engineered feature.
    compute_correlation_btw_svd_components_and_multiple_engineered_features(self, pcs=[0], denoised=False, engineered_feature_names=None):
        Computes the correlation between SVD components and multiple engineered features.
    plot_heatmaps_correlations(self, nb_pcs=10, denoised=False):
        Plots heatmaps of the correlations between SVD components and engineered features.
    plot_scree_plots(self):
        Plots scree plots for each model in the pipeline list.
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
                 cluster_threshold=0.5,
                 corr_type='spearman',
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
        

        self.corr_type = corr_type
        self.cluster_threshold = cluster_threshold

 

    def compute_correlation_btw_svd_components_and_one_feature(self,
                                                                 pcs=None, 
                                                                 denoised=False,
                                                                 feature_values=None,
                                                                 feature_name=None,
                                                                 pct=None):
        
        """
            Computes the correlation between specified SVD components and multiple engineered features.
            Parameters:
            -----------
            pcs : list of int, optional
                List of principal components to compute correlations for. Default is 10.
            denoised : bool, optional
                If True, use denoised image embeddings. Default is False.
            engineered_feature_names : list of str, optional
                List of engineered feature names to compute correlations with. If None, all features in `self.ef.emb_df` are used. Default is None.
            Returns:
            --------
            correlations : dict
                Dictionary containing correlation statistics and p-values for each model and principal component.
                Structure:
                {
                'stat': {
                    model_name: {
                    'PC{pc+1}': {
                        engineered_feature_name: spearman_statistic
                    }
                    }
                },
                'p_value': {
                    model_name: {
                    'PC{pc+1}': {
                        engineered_feature_name: p_value
                    }
                    }
                }
                }
        """
        
        if denoised:
            matrix_name = "denoised"
        else:
            matrix_name = "raw"

    
        correlations = {}
        correlations['stat'] = {}
        correlations['p_value'] = {}

        ordered_models = self.pipelines_list.copy()

        if self.image_embeddings[ordered_models[0]].svd['denoised'] != denoised:
            self.compute_svd_image_embeddings(denoised)
            
        for model in ordered_models:
            correlations['stat'][model] = {}
            correlations['p_value'][model] = {}
            if len(self.pipelines_list) != 0:
                if model not in self.pipelines_list:
                    continue

            svd = self.image_embeddings[model].svd

            if pct is None:
                if pcs is None or pcs > svd['U_df'].shape[1]:
                    pcs = svd['U_df'].shape[1]
            else:
                pcs = np.argmax(np.cumsum(self.image_embeddings[model].get_explained_variance_ratio_list()) > pct)

            for pc in range(pcs):

                pc_values = svd['U_df'][f"u{pc+1}"].values
                if self.corr_type == 'spearman':
                    s, p = spearmanr(pc_values, feature_values, nan_policy='omit')
                else:
                    nas = np.logical_or(np.isnan(pc_values), np.isnan(feature_values))
                    s, p = pearsonr(pc_values[~nas], feature_values[~nas])
                correlations['stat'][model][f"PC{pc+1}"] = s
                correlations['p_value'][model][f"PC{pc+1}"]= p
        
        return correlations
    

    def compute_correlation_btw_svd_components_and_multiple_engineered_features(self,
                                                                                pcs=None, 
                                                                                denoised=False,
                                                                                engineered_feature_names=None, 
                                                                                pct=None):
        
        """
            Computes the correlation between specified SVD components and multiple engineered features.
            Parameters:
            -----------
            pcs : list of int, optional
                List of principal components to compute correlations for. Default is 10.
            denoised : bool, optional
                If True, use denoised image embeddings. Default is False.
            engineered_feature_names : list of str, optional
                List of engineered feature names to compute correlations with. If None, all features in `self.ef.emb_df` are used. Default is None.
            Returns:
            --------
            correlations : dict
                Dictionary containing correlation statistics and p-values for each model and principal component.
                Structure:
                {
                'stat': {
                    model_name: {
                    'PC{pc+1}': {
                        engineered_feature_name: spearman_statistic
                    }
                    }
                },
                'p_value': {
                    model_name: {
                    'PC{pc+1}': {
                        engineered_feature_name: p_value
                    }
                    }
                }
                }
        """
        
        if denoised:
            matrix_name = "denoised"
        else:
            matrix_name = "raw"

        if engineered_feature_names is None:
            engineered_feature_names = list(self.ef.emb_df.columns)
    
    
        ordered_models = self.pipelines_list.copy()

        if self.image_embeddings[ordered_models[0]].svd['denoised'] != denoised:
            self.compute_svd_image_embeddings(denoised)
            
        for model in ordered_models:
            if model not in self.correlations_image_pcs_and_handcrafted_features[matrix_name]['stat'].keys():
                self.correlations_image_pcs_and_handcrafted_features[matrix_name]['stat'][model] = {}
                self.correlations_image_pcs_and_handcrafted_features[matrix_name]['p_value'][model] = {}
                # if len(self.pipelines_list) != 0:
                #     if model not in self.pipelines_list:
                #         continue

            svd = self.image_embeddings[model].svd

            if pct is None:
                if pcs is None or pcs > svd['U_df'].shape[1]:
                    pcs = svd['U_df'].shape[1]
            else:
                pcs = np.argmax(np.cumsum(self.image_embeddings[model].get_explained_variance_ratio_list()) > pct)

            for pc in range(pcs):
                if f"PC{pc+1}" not in self.correlations_image_pcs_and_handcrafted_features[matrix_name]['stat'][model].keys():
                    self.correlations_image_pcs_and_handcrafted_features[matrix_name]['stat'][model][f"PC{pc+1}"] = {}
                    self.correlations_image_pcs_and_handcrafted_features[matrix_name]['p_value'][model][f"PC{pc+1}"] = {}
                
                eng_features_to_compute = [feature for feature in engineered_feature_names if feature not in self.correlations_image_pcs_and_handcrafted_features[matrix_name]['stat'][model][f"PC{pc+1}"].keys()]

                pc_values = svd['U_df'].loc[list(self.ef.emb_df.index)][f"u{pc+1}"].values
                print(f"Computing correlations for PC{pc+1} and {model}...", flush=True)
                print(f"Engineered features to compute: {eng_features_to_compute}: ", flush=True)
                for engineered_feature_name in eng_features_to_compute:
                    if self.corr_type == 'spearman':
                        s, p = spearmanr(pc_values, self.ef.emb_df[engineered_feature_name], nan_policy='omit')
                    else:
                        nas = np.logical_or(np.isnan(pc_values), np.isnan(self.ef.emb_df[engineered_feature_name]))
                        s, p = pearsonr(pc_values[~nas], self.ef.emb_df[engineered_feature_name][~nas])
                    self.correlations_image_pcs_and_handcrafted_features[matrix_name]['stat'][model][f"PC{pc+1}"][engineered_feature_name] = s
                    self.correlations_image_pcs_and_handcrafted_features[matrix_name]['p_value'][model][f"PC{pc+1}"][engineered_feature_name] = p
        
    
    def compute_explanation_scores(self, denoised=False, pct_variance=0.9):
        """
        Compute explanation scores for image embeddings.
        This method calculates the explanation scores for image embeddings by 
        computing the correlation between SVD components and multiple engineered 
        features. The scores are weighted by the explained variance ratio of the 
        components.
        Args:
            denoised (bool): If True, use denoised image embeddings; otherwise, use raw image embeddings.
            pct_variance (float): The percentage of variance to be explained by the selected components.
        Returns:
            None: The method updates the `explanation_scores` attribute of the class with the computed scores.
        """


        if denoised:
            matrix_name = "denoised"
        else:
            matrix_name = "raw"
        
        ordered_models = self.pipelines_list.copy()

        for model in ordered_models:
            # if model not in self.explanation_scores.keys():
            # s = self.image_embeddings[matrix_name][model].svd['S']
            explained_variance_ratio = self.image_embeddings[model].get_explained_variance_ratio_list()
            n_comp = np.argmax(np.cumsum(explained_variance_ratio) > pct_variance)

            # if f"PC{n_comp+1}" not in self.correlations_image_pcs_and_handcrafted_features[matrix_name]['stat'][model].keys():
            #     self.compute_correlation_btw_svd_components_and_multiple_engineered_features(pcs=n_comp, denoised=denoised)

            # Get the correlation scores for these components
            corr_scores = abs(pd.DataFrame.from_dict(self.correlations_image_pcs_and_handcrafted_features[matrix_name]['stat'][model]))[[f"PC{i+1}" for i in range(n_comp)]]
            vars = explained_variance_ratio[:n_comp]

            scores = corr_scores.mul(vars, axis=1).sum(axis=1).to_dict()
            self.explanation_scores[model] = scores
        


    @staticmethod
    def SVD_correlation_two_pipelines(emb1, emb2):
        """
        Compute the Pearson correlation coefficients and p-values between the singular vectors
        of two embeddings obtained from Singular Value Decomposition (SVD).
        Parameters:
        emb1 (object): The first embedding object containing SVD results with a 'U_df' DataFrame.
        emb2 (object): The second embedding object containing SVD results with a 'U_df' DataFrame.
        Returns:
        tuple: A tuple containing two lists:
            - correlation_vector (list): A list of absolute Pearson correlation coefficients between 
              the singular vectors of emb1 and emb2.
            - p_value_vector (list): A list of p-values corresponding to the Pearson correlation coefficients.
        """

        correlation_vector = []
        p_value_vector = []
        for i in range(50):
            for j in range(50):
                s, p = pearsonr(emb1.svd['U_df'][f"u{i+1}"], emb2.svd['U_df'][f"u{j+1}"])
                correlation_vector.append(abs(s))
                p_value_vector.append(p)
        return correlation_vector, p_value_vector       
    
    def compute_SVD_correlations_between_pipelines(self):
        """
        Compute Singular Value Decomposition (SVD) correlations between different pipelines and visualize the results.
        This method calculates the maximum and median SVD correlations between pairs of pipelines and stores the results 
        in dataframes. It then generates and saves heatmaps and boxplots to visualize these correlations.
        The method performs the following steps:
        1. Initializes dataframes to store maximum and median SVD correlations.
        2. Iterates over pairs of pipelines to compute SVD correlations.
        3. Stores the maximum and median correlations in the respective dataframes.
        4. Plots and saves heatmaps of the maximum and median correlations.
        5. Plots and saves boxplots of the SVD correlations, both with and without outliers, and with and without intra-pipeline correlations.
        Parameters:
        None
        Returns:
        None
        """

        df_SVD_max_corrs = pd.DataFrame(columns=self.pipelines_list, index=self.pipelines_list)
        df_SVD_median_corrs = pd.DataFrame(columns=self.pipelines_list, index=self.pipelines_list)
        df_corrs = pd.DataFrame()
        for i in range(len(self.pipelines_list)):
            for j in range(i, len(self.pipelines_list)):
                corrs, _ = self.SVD_correlation_two_pipelines(self.image_embeddings[self.pipelines_list[i]],
                                                            self.image_embeddings[self.pipelines_list[j]])

                df_SVD_max_corrs.loc[self.pipelines_list[i], self.pipelines_list[j]] = max(corrs)
                df_SVD_max_corrs.loc[self.pipelines_list[j], self.pipelines_list[i]] = max(corrs)
                df_SVD_median_corrs.loc[self.pipelines_list[i], self.pipelines_list[j]] = np.median(corrs)
                df_SVD_median_corrs.loc[self.pipelines_list[j], self.pipelines_list[i]] = np.median(corrs)
                df_corrs[f"{self.pipelines_list[i]}_{self.pipelines_list[j]}"] = corrs
        
        # Plot the results
        df_SVD_max_corrs = df_SVD_max_corrs.apply(pd.to_numeric, errors='coerce')
        # Plot the heatmap of max correlations
        plt.figure(figsize=(10, 8))
        sns.heatmap(df_SVD_max_corrs, 
                    annot=True, 
                    cmap=custom_diverging_map, 
                    vmin=df_SVD_max_corrs.min().min(), 
                    vmax=df_SVD_max_corrs.max().max(), 
                    center=df_SVD_max_corrs.min().min(), 
                    cbar=True)

        plt.title('Heatmap of SVD maximal Correlations', weight='bold', fontsize=12)
        plt.savefig(os.path.join(self.saving_folder, f"heatmap_SVD_max_correlations.{self.extension}"), bbox_inches='tight')

        # Plot the heatmap of median correlations
        df_SVD_median_corrs = df_SVD_median_corrs.apply(pd.to_numeric, errors='coerce')
        plt.figure(figsize=(10, 8))
        sns.heatmap(df_SVD_median_corrs, 
                    annot=True, 
                    cmap=custom_diverging_map, 
                    vmin=df_SVD_median_corrs.min().min(), 
                    vmax=df_SVD_median_corrs.max().max(), 
                    center=df_SVD_median_corrs.min().min(), 
                    cbar=True)
        plt.title('Heatmap of SVD median Correlations', weight='bold', fontsize=12)
        plt.savefig(os.path.join(self.saving_folder, f"heatmap_SVD_median_correlations.{self.extension}"), bbox_inches='tight')

        # Plot the distribution of correlations
        plt.figure()
        sns.boxplot(df_corrs, 
            showfliers=True, 
            linewidth=2.5,
            color='salmon')
        plt.xticks(rotation=90)
        sns.despine()
        plt.title('Boxplot of SVD correlations between pipelines', weight='bold', fontsize=12)
        plt.savefig(os.path.join(self.saving_folder, f"boxplot_SVD_correlations_with_intra_with_outliers.{self.extension}"), bbox_inches='tight')

        plt.figure()
        sns.boxplot(df_corrs.drop(columns=[f"{model}_{model}" for model in self.pipelines_list]), 
            showfliers=False, 
            linewidth=2.5, color='salmon')
        
        sns.despine()
        plt.xticks(rotation=90)
        plt.title('Boxplot of SVD correlations between pipelines', weight='bold', fontsize=12)
        plt.savefig(os.path.join(self.saving_folder, f"boxplot_SVD_correlations_without_intra_without_outliers.{self.extension}"), bbox_inches='tight')
    
    def plot_heatmaps_correlations(self, 
                                   nb_pcs=10, 
                                   denoised=False, 
                                   engineered_feature_names=None, 
                                   group_features=False, 
                                   stat=True):
        
        """
        Plots heatmaps of the correlations between SVD components and engineered features.
        This method computes the image embeddings and the correlations between the SVD components 
        and multiple engineered features. It then generates and saves heatmaps of these correlations.
        Parameters:
        -----------
        nb_pcs : int, optional
            The number of principal components to consider for the SVD. Default is 10.
        denoised : bool, optional
            If True, use denoised image embeddings. Otherwise, use raw image embeddings. Default is False.
        Returns:
        --------
        None
        """
        self.initialize_model_saving_folders()

        if denoised:
            matrix_name = "denoised"
        else:
            matrix_name = "raw"
        
        # if matrix_name not in self.correlations_image_pcs_and_handcrafted_features.keys():
        #     self.compute_correlation_btw_svd_components_and_multiple_engineered_features(pcs=list(range(nb_pcs)), denoised=denoised, engineered_feature_names=engineered_feature_names)

        if engineered_feature_names is None:
            engineered_feature_names = list(self.ef.emb_df.columns)

        features_per_type = {}
        for col in engineered_feature_names:
            if col.split("_")[0] not in features_per_type:
                features_per_type[col.split("_")[0]] = []
            features_per_type[col.split("_")[0]].append(col)

        if group_features:
            print("Plotting heatmaps per feature group...", flush=True)
            for feature_type, feature_list in features_per_type.items():
                if stat:
                    for model, model_dict in self.correlations_image_pcs_and_handcrafted_features[matrix_name]['stat'].items():
                        plt.figure(figsize=(10, 10))
                        df = pd.DataFrame.from_dict(model_dict)[[f"PC{i+1}" for i in range(nb_pcs)]]
                        df.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
                        df = df.loc[feature_list]
                        sns.clustermap(df, cmap=custom_diverging_map, vmin=-1, vmax=1, yticklabels=True, xticklabels=True, row_cluster=True, col_cluster=True)
                        plt.title(f"Correlation between {matrix_name} SVD components and {feature_type} features \n {model}", weight='bold') 
                        plt.savefig(os.path.join(self.saving_folder, model, f"correlation_{matrix_name}_svd_{feature_type}_features_{model}_.{self.extension}"), bbox_inches='tight')
                else:
                    for model, model_dict in self.correlations_image_pcs_and_handcrafted_features[matrix_name]['p_value'].items():
                        plt.figure(figsize=(10, 10))
                        df = -np.log10(pd.DataFrame.from_dict(model_dict)[[f"PC{i+1}" for i in range(nb_pcs)]])
                        df = df.loc[feature_list]
                        df.replace([np.inf, -np.inf], 300, inplace=True)
                        df.replace(np.nan, 0, inplace=True)
                        sns.clustermap(df, cmap=custom_diverging_map, center=0, vmin=df.min().min(), vmax=df.max().max(), yticklabels=True, xticklabels=True, row_cluster=True, col_cluster=True)
                        plt.title(f"Correlation p-value between {matrix_name} SVD components and {feature_type} features \n {model}", weight='bold') 
                        plt.savefig(os.path.join(self.saving_folder, model, f"p_value_correlation_{matrix_name}_svd_{feature_type}_features_{model}_.{self.extension}"), bbox_inches='tight')
        else:   
                ## Heatmap of the correlations
            print("Plotting whole heatmap...", flush=True)
            if stat:
                for model, model_dict in self.correlations_image_pcs_and_handcrafted_features[matrix_name]['stat'].items():
                    plt.figure(figsize=(10, 10))
                    df = pd.DataFrame.from_dict(model_dict)[[f"PC{i+1}" for i in range(nb_pcs)]]
                    df.replace([np.inf, -np.inf, np.nan], 0, inplace=True)


                    if self.ef.feature_type_color_dict is None:

                        sns.clustermap(df, cmap=custom_diverging_map, vmin=-1, vmax=1, row_cluster=True, col_cluster=False, yticklabels=False, xticklabels=True)
                        plt.title(f"Correlation between raw SVD components and handcrafted features \n {model}", weight='bold') 
                        plt.savefig(os.path.join(self.saving_folder, model, f"correlation_{matrix_name}_svd_handcrafted_features_{model}_.{self.extension}"), bbox_inches='tight')

                    else:
                        print(df.index)
                        colors_ = [self.ef.feature_type_color_dict[col.split("_")[0]] for col in df.index]

                        sns.clustermap(df, cmap=custom_diverging_map, vmin=-1, vmax=1, row_cluster=True, col_cluster=False, yticklabels=False, xticklabels=True, row_colors=colors_)
                        plt.title(f"Correlation between raw SVD components and handcrafted features \n {model}", weight='bold') 
                        plt.savefig(os.path.join(self.saving_folder, model, f"correlation_{matrix_name}_svd_handcrafted_features_{model}_.{self.extension}"), bbox_inches='tight')
            else:
                for model, model_dict in self.correlations_image_pcs_and_handcrafted_features[matrix_name]['p_value'].items():
                    plt.figure(figsize=(10, 10))
                    df = -np.log10(pd.DataFrame.from_dict(model_dict)[[f"PC{i+1}" for i in range(nb_pcs)]])
                    df.replace([np.inf, -np.inf], 300, inplace=True)
                    df.replace(np.nan, 0, inplace=True)


                    if self.ef.feature_type_color_dict is None:
                        sns.clustermap(df, cmap=custom_diverging_map, center=0, vmin=df.min().min(), vmax=df.max().max(), row_cluster=True, col_cluster=False, yticklabels=False, xticklabels=True)
                        plt.title(f"Correlation p-value between raw SVD components and handcrafted features \n {model}", weight='bold') 
                        plt.savefig(os.path.join(self.saving_folder, model, f"p_value_correlation_{matrix_name}_svd_handcrafted_features_{model}_.{self.extension}"), bbox_inches='tight')

                    else:
                        colors = [self.ef.feature_type_color_dict[col.split("_")[0]] for col in df.index]

                        sns.clustermap(df, cmap=custom_diverging_map, center=0, vmin=df.min().min(), vmax=df.max().max(), row_cluster=True, col_cluster=False, yticklabels=False, xticklabels=True, row_colors=colors)
                        plt.title(f"Correlation p-value between raw SVD components and handcrafted features \n {model}", weight='bold') 
                        plt.savefig(os.path.join(self.saving_folder, model, f"p_value_correlation_{matrix_name}_svd_handcrafted_features_{model}_.{self.extension}"), bbox_inches='tight')



    def boxplot_correlations_per_feature_type(self, nb_pcs=5, denoised=False, engineered_feature_names=None):

        if denoised:
            matrix_name = "denoised"
        else:
            matrix_name = "raw"

        # if matrix_name not in self.correlations_image_pcs_and_handcrafted_features.keys():
        #     self.compute_correlation_btw_svd_components_and_multiple_engineered_features(pcs=list(range(nb_pcs)), denoised=denoised, engineered_feature_names=engineered_feature_names)

        if engineered_feature_names is None:
            engineered_feature_names = list(self.ef.emb_df.columns)

        features_per_type = {}
        for col in engineered_feature_names:
            if col.split("_")[0] not in features_per_type:
                features_per_type[col.split("_")[0]] = []
            features_per_type[col.split("_")[0]].append(col)

        print("Plotting boxplots per feature group...", flush=True)
        for pc in range(nb_pcs):
            df_models = pd.DataFrame()
            for model_name, model_dict in self.correlations_image_pcs_and_handcrafted_features[matrix_name]['stat'].items():
                df_models[model_name] = model_dict[f'PC{pc+1}']

            for feature_type, feature_list in features_per_type.items():
                print(f"Plotting boxplot for {feature_type} features...", flush=True)
                print(f"length of the dataframe: {len(df_models.loc[feature_list])}")

                print(f"{feature_type} features: {feature_list}")
                print("Dataframe length: ", len(abs(df_models.loc[feature_list])))
                print("Dataframe columns: ", abs(df_models.loc[feature_list]).columns)
                print("Dataframe head: ", abs(df_models.loc[feature_list]).head())

                plt.figure(figsize=(2*len(self.correlations_image_pcs_and_handcrafted_features[matrix_name]['stat']), 6))
                sns.boxplot(data=abs(df_models.loc[feature_list]))
                
                plt.title(f"Absolute correlation between PC{pc+1} and {feature_type} features", weight='bold') 
                plt.xticks(rotation=90)
                plt.ylabel("Absolute correlation")
                plt.ylim(0,1)
                plt.savefig(os.path.join(self.saving_folder, f"boxplot_correlation_{feature_type}_features_with_PC{pc+1}_per_model.{self.extension}"), bbox_inches='tight')
    

    def heatmap_correlations_per_feature_type_per_pc(self, nb_pcs=5, denoised=False, engineered_feature_names=None):

        if denoised:
            matrix_name = "denoised"
        else:
            matrix_name = "raw"

        if matrix_name not in self.correlations_image_pcs_and_handcrafted_features.keys():
            self.compute_correlation_btw_svd_components_and_multiple_engineered_features(pcs=list(range(nb_pcs)), denoised=denoised, engineered_feature_names=engineered_feature_names)

        if engineered_feature_names is None:
            engineered_feature_names = list(self.ef.emb_df.columns)

        features_per_type = {}
        for col in engineered_feature_names:
            if col.split("_")[0] not in features_per_type:
                features_per_type[col.split("_")[0]] = []
            features_per_type[col.split("_")[0]].append(col)

        print("Plotting boxplots per feature group...", flush=True)
        for pc in range(nb_pcs):
            df_models = pd.DataFrame()
            for model_name, model_dict in self.correlations_image_pcs_and_handcrafted_features[matrix_name]['stat'].items():
                df_models[model_name] = model_dict[f'PC{pc+1}']

            for feature_type, feature_list in features_per_type.items():
                plt.figure(figsize=(10, int(np.ceil(0.2*len(feature_list)))))
                sns.heatmap(data=df_models.loc[feature_list], cmap=custom_diverging_map, vmin=-1, vmax=1, yticklabels=True)
                plt.title(f"Correlation between PC{pc+1} and {feature_type} features", weight='bold') 
                plt.savefig(os.path.join(self.saving_folder, f"heatmap_correlation_{feature_type}_features_with_PC{pc+1}_per_model.{self.extension}"), bbox_inches='tight')


    def heatmap_explanation_scores(self, group_features=True, threshold=0.15, features_subset=None):
        df_explanation_scores = pd.DataFrame.from_dict(self.explanation_scores)
        feature_types = list(set([idx.split("_")[0] for idx in df_explanation_scores.index]))
        print(feature_types)
        if features_subset is not None:
            df_explanation_scores_to_keep = [idx for idx in df_explanation_scores.index if idx.split("_")[0] in features_subset]
            df_explanation_scores = df_explanation_scores.loc[df_explanation_scores_to_keep]
            features_names = "_".join(features_subset)
        else:
            features_subset = feature_types
            features_names = "all"

        if group_features:
            for feature_type in feature_types:
                plt.figure()
                sns.clustermap(df_explanation_scores.loc[[idx for idx in df_explanation_scores.index if idx.split("_")[0] == feature_type]], 
                            cmap=sns.color_palette("Blues", as_cmap=True),
                            vmin=df_explanation_scores.min().min(), 
                            vmax=df_explanation_scores.max().max(), 
                            row_cluster=True,
                                col_cluster=True,
                            )
                plt.title(f"Explanation scores for {feature_type} features", weight='bold')
                plt.savefig(os.path.join(self.saving_folder, f"heatmap_explanation_scores_{feature_type}.{self.extension}"), bbox_inches='tight')
        else:
            if self.ef.feature_type_color_dict is None:
                plt.figure()
                sns.clustermap(df_explanation_scores,
                            cmap=custom_diverging_map,
                                vmin=df_explanation_scores.min().min(),
                                vmax=df_explanation_scores.max().max(),
                                row_cluster=True,
                                col_cluster=True,
                                xticklabels=True, 
                                center=0)
                plt.title(f"Explanation scores for all features", weight='bold')
                plt.savefig(os.path.join(self.saving_folder, f"heatmap_explanation_scores_{features_names}_features.{self.extension}"), bbox_inches='tight')
            
            else:
                plt.figure()
                colors = [self.ef.feature_type_color_dict[col.split("_")[0]] for col in df_explanation_scores.index]
                g = sns.clustermap(df_explanation_scores,
                            cmap=custom_diverging_map,
                                vmin=df_explanation_scores.min().min(),
                                vmax=df_explanation_scores.max().max(),
                                row_cluster=True,
                                col_cluster=True,
                                row_colors=colors,
                                xticklabels=True, 
                                center=0,
                                figsize=(6,10),
                                yticklabels=False)
                plt.title(f"Explanation scores for all features", weight='bold')
                plt.savefig(os.path.join(self.saving_folder, f"heatmap_explanation_scores_{features_names}_features.{self.extension}"), bbox_inches='tight')

                row_linkage = g.dendrogram_row.linkage
                plt.figure(figsize=(10, 7))
                dendrogram(row_linkage, labels=df_explanation_scores.index)
                plt.title('Dendrogram')
                plt.xlabel('Features')
                plt.ylabel('Distance')
                plt.savefig(os.path.join(self.saving_folder, f"dendrogram_explanation_scores_{features_names}_features.{self.extension}"), bbox_inches='tight')

                row_clusters = fcluster(row_linkage, threshold, criterion='distance')

                # Create a color mapping for the clusters
                unique_clusters = np.unique(row_clusters)
                cluster_colors = sns.color_palette('tab20', len(unique_clusters))
                row_colors = [cluster_colors[cluster_id - 1] for cluster_id in row_clusters]

                # Map clusters to row names
                clustered_features = {}
                for cluster_id in unique_clusters:
                    clustered_features[int(cluster_id)] = df_explanation_scores.index[row_clusters == cluster_id].tolist()

                # Print the clustered features with their associated colors
                for cluster_id, features in clustered_features.items():
                    print(f"Cluster {cluster_id}: {features}")

                # Create custom legend

                # Replot the ClusterMap with row colors
                plt.figure()
                g = sns.clustermap(df_explanation_scores, 
                                   cmap=custom_diverging_map, 
                                   row_cluster=True, 
                                   col_cluster=True, 
                                   row_colors=row_colors, 
                                   center=0, 
                                   yticklabels=False, 
                                   figsize=(6,10))
                handles = [mpatches.Patch(color=cluster_colors[i], label=f'Cluster {i+1}') for i in range(len(unique_clusters))]
                plt.legend(handles=handles, title='Clusters', bbox_to_anchor=(20, 1.4), loc='upper left', borderaxespad=0.)
                plt.savefig(os.path.join(self.saving_folder, f"heatmap_explanation_scores_{features_names}_features_clustered.{self.extension}"), bbox_inches='tight')


                # Quantify the proportion of features starting with each key for each cluster
                proportions = []
                for cluster_id, features in clustered_features.items():
                    cluster_proportions = {}
                    for key in features_subset:
                        count_key_features = sum(1 for feature in features if feature.startswith(key))
                        proportion_key_features = count_key_features / len(features)
                        cluster_proportions[key] = proportion_key_features
                    proportions.append((cluster_id, cluster_proportions))

                # Prepare data for bar plots
                proportion_data = []
                for cluster_id, cluster_proportions in proportions:
                    for key, proportion in cluster_proportions.items():
                        proportion_data.append((cluster_id, key, proportion))

                proportion_df = pd.DataFrame(proportion_data, columns=['Cluster', 'Feature', 'Proportion'])

                # Plot bar plots
                plt.figure()
                sns.set_style(style="whitegrid")
                barplot = sns.barplot(x='Cluster', y='Proportion', hue='Feature', data=proportion_df, palette=self.ef.feature_type_color_dict)
                barplot.set_title('Proportion of Features per type in each cluster', weight='bold')
                plt.legend(title='Feature', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
                plt.savefig(os.path.join(self.saving_folder, f"barplot_proportion_features_per_cluster_from_heatmap_{features_names}.{self.extension}"), bbox_inches='tight')


                total_features = {key: sum(1 for feature in df_explanation_scores.index if feature.startswith(key)) for key in features_subset}

                # Quantify the proportion of features of each type in each cluster
                proportions = []
                for cluster_id, features in clustered_features.items():
                    cluster_proportions = {}
                    for key in features_subset:
                        count_key_features = sum(1 for feature in features if feature.startswith(key))
                        proportion_key_features = count_key_features / total_features[key]
                        cluster_proportions[key] = proportion_key_features
                    proportions.append((cluster_id, cluster_proportions))

                # Prepare data for bar plots
                proportion_data = []
                for cluster_id, cluster_proportions in proportions:
                    for key, proportion in cluster_proportions.items():
                        proportion_data.append((cluster_id, key, proportion))

                proportion_df = pd.DataFrame(proportion_data, columns=['Cluster', 'Feature', 'Proportion'])

                # Plot bar plots
                plt.figure()
                sns.set_style(style="whitegrid")
                barplot = sns.barplot(x='Feature', y='Proportion', hue='Cluster', data=proportion_df, palette=cluster_colors)
                barplot.set_title('Features representation', weight='bold')
                plt.xticks(rotation=90)
                plt.legend(title='Cluster', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
                plt.savefig(os.path.join(self.saving_folder, f"barplot_proportion_features_type_in_clusters_from_heatmap_{features_names}.{self.extension}"), bbox_inches='tight')

                # Save the clustered features
                with open(os.path.join(self.saving_folder, f'clustered_features_from_heatmap_{features_names}.json'), 'w') as f:
                    json.dump(clustered_features, f)


    def boxplot_explanation_scores(self):
        df_explanation_scores = pd.DataFrame.from_dict(self.explanation_scores)
        feature_types = list(set([idx.split("_")[0] for idx in df_explanation_scores.index]))

        for feature_type in feature_types:
            plt.figure(figsize=(2*len(self.explanation_scores), 6))
            sns.boxplot(data=df_explanation_scores.loc[[idx for idx in df_explanation_scores.index if idx.split("_")[0] == feature_type]])


            plt.title(f"Explanation scores for {feature_type} features", weight='bold')
            plt.xticks(rotation=90)
            plt.ylabel("Explanation score")
            plt.ylim(df_explanation_scores.min().min(),df_explanation_scores.max().max())
            plt.savefig(os.path.join(self.saving_folder, f"boxplot_explanation_scores_{feature_type}.{self.extension}"), bbox_inches='tight')


    def get_pc_correlation_image_embeddings(self, pc=1):
        df_pcs_correlations = pd.DataFrame()

        for model1 in self.pipelines_list:
            for model2 in self.pipelines_list:
                if model1 != model2:
                    s, _ = spearmanr(self.image_embeddings[model1].svd['U_df'][f'u{pc}'], self.image_embeddings[model2].svd['U_df'][f'u{pc}'])
                    df_pcs_correlations.loc[model1, model2] = abs(s)
                else:
                    df_pcs_correlations.loc[model1, model2] = 1

        plt.figure()
        sns.heatmap(df_pcs_correlations, annot=True, cmap=custom_diverging_map, center=0, vmin=0, vmax=1)
        plt.title(f"Correlation between SVD components (image embeddings): u{pc}", weight='bold')
        plt.savefig(os.path.join(self.saving_folder, f"heatmap_correlation_u{pc}.{self.extension}"), bbox_inches='tight')


    def correlations_between_svd_components_and_engineered_features_pipeline(self, pct_variance=0.9, n_pcs_plot=10):

        # Compute correlations with engineered features 
        if os.path.exists(os.path.join(self.saving_folder, f'correlations_{self.corr_type}.json')):
            with open(os.path.join(self.saving_folder, f'correlations_{self.corr_type}.json'), 'r') as f:
                self.correlations_image_pcs_and_handcrafted_features = json.load(f)
        
        # Will compute the correlations if the json file does not exist or if models or features are missing
        self.compute_correlation_btw_svd_components_and_multiple_engineered_features(pct=pct_variance, denoised=False)



                # ## Avoid recomputing these correlations and taking the json files previously computed into account
                # ## The json file needs to be put in the right folder beforehand
                # if len([model for model in self.pipelines_list if model not in models_in_corr_json]) > 0:
                #     print("Models in the benchmark that are not in the correlations json file. Computing correlations for these models...")


        with open(os.path.join(self.saving_folder, f'correlations_{self.corr_type}.json'), 'w') as f:
            json.dump(self.correlations_image_pcs_and_handcrafted_features, f)

        models_in_corr_json = list(self.correlations_image_pcs_and_handcrafted_features['raw']['stat'].keys())


        if set(models_in_corr_json) != set(self.pipelines_list):
            models_to_remove = [model for model in models_in_corr_json if model not in self.pipelines_list]
            for stat_type in self.correlations_image_pcs_and_handcrafted_features["raw"].keys():
                for model in models_to_remove:
                    del self.correlations_image_pcs_and_handcrafted_features["raw"][stat_type][model]


        # Plots associated with correlation
        nb_pcs_min = min([len(self.correlations_image_pcs_and_handcrafted_features['raw']['stat'][model].keys()) for model in self.pipelines_list])
        nb_pcs_plot = min(n_pcs_plot, nb_pcs_min)
        self.plot_heatmaps_correlations(nb_pcs=nb_pcs_plot, denoised=False, group_features=False, stat=True)
        self.plot_heatmaps_correlations(nb_pcs=nb_pcs_plot, denoised=False, group_features=True, stat=False)
        self.plot_heatmaps_correlations(nb_pcs=nb_pcs_plot, denoised=False, group_features=True, stat=True)
        self.plot_heatmaps_correlations(nb_pcs=nb_pcs_plot, denoised=False, group_features=False, stat=False)


        self.boxplot_correlations_per_feature_type(nb_pcs=nb_pcs_plot, denoised=False)
        self.heatmap_correlations_per_feature_type_per_pc(nb_pcs=nb_pcs_plot, denoised=False)


    def explanation_scores_per_model_per_engineered_feature_pipeline(self, pct_variance=0.9, enforce_recompute=False):


        if enforce_recompute:
            explanation_scores_handcrafted = self.ef.compute_explanation_scores()
            explanation_scores_handcrafted.to_csv(os.path.join(self.saving_folder, f"explanation_scores_handcrafted.csv"))

            self.compute_explanation_scores(denoised=False, pct_variance=pct_variance)

            with open(os.path.join(self.saving_folder, f'explanation_scores_{self.corr_type}.json'), 'w') as f:
                json.dump(self.explanation_scores, f)

        else:

            # Compute the correlations scores from the handcrafted features matrix
            if not os.path.exists(os.path.join(self.saving_folder, f"explanation_scores_handcrafted.csv")):
                explanation_scores_handcrafted = self.ef.compute_explanation_scores()
                explanation_scores_handcrafted.to_csv(os.path.join(self.saving_folder, f"explanation_scores_handcrafted.csv"))


            if os.path.exists(os.path.join(self.saving_folder, f'explanation_scores_{self.corr_type}.json')):
                with open(os.path.join(self.saving_folder, f'explanation_scores_{self.corr_type}.json'), 'r') as f:
                    self.explanation_scores = json.load(f)
                    
                    # Remove models not used in the benchmark
                    if set(self.explanation_scores.keys()) != set(self.pipelines_list):
                        models_to_remove = [model for model in self.explanation_scores.keys() if model not in self.pipelines_list]
                        for model in models_to_remove:
                            del self.explanation_scores[model]

                        # Compute explanation score if not already done
                        models_in_exp_json = list(self.explanation_scores.keys())
                        if len([model for model in self.pipelines_list if model not in models_in_exp_json]) > 0:
                            print("Models in the benchmark that are not in the explanation scores json file. Computing explanation scores for these models...")

                            with open(os.path.join(self.saving_folder, f'explanation_scores_{self.corr_type}.json'), 'w') as f:
                                json.dump(self.explanation_scores, f)

            else:
                self.compute_explanation_scores(denoised=False, pct_variance=pct_variance)

                with open(os.path.join(self.saving_folder, f'explanation_scores_{self.corr_type}.json'), 'w') as f:
                    json.dump(self.explanation_scores, f)


        # Plot explanation scores
        self.boxplot_explanation_scores()
        self.heatmap_explanation_scores()
        self.heatmap_explanation_scores(group_features=True)
        self.heatmap_explanation_scores(group_features=False, threshold=self.cluster_threshold, features_subset=None)

    