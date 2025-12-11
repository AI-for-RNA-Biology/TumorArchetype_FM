#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#

import os
import glob
import numpy as np
import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from digitalhistopathology.benchmark.benchmark_base import BenchmarkBase


class BenchmarkShannon(BenchmarkBase):

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
                 min_comp=512,
                 pct_variance=0.9,
                 dataset="HER2",

                 ):
        
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
        
        self.min_comp = min_comp  # Minimum number of components to keep for the SVD
        self.pct_variance = pct_variance  # Percentage of variance to keep for the SVD



    def nb_components_explaining_variance(self, pct=0.9):

        nb_components = {}
        for model in self.pipelines_list:
            explained_variance_ratio = self.image_embeddings[model].get_explained_variance_ratio_list()
            n_comp = np.argmax(np.cumsum(explained_variance_ratio) > pct)
            nb_components[model] = int(n_comp)
        return nb_components


    def shannon_entropies_and_kl_divergence_all_models(self, pct=None, denoised=False, n_comp=None, rescale=False, group=None, retrieve_from_file=False, group_name_in_file='per_patient'):

        if pct is not None:
            name_comp = f'comp_explaining_{int(np.round(pct*100))}%_variance'
        elif n_comp is None:
            name_comp = 'all_comps'
        else:
            name_comp = f'{n_comp}_comps'

        if denoised:
            name_denoised = 'denoised'
        else:
            name_denoised = 'raw'

        if rescale:
            name_rescale = '_rescaled'
        else:
            name_rescale = ''
            
        if group is not None:

            entropy_models = {}
            kl_divergence_models = {}

            for path_to_pipeline, model in zip(self.path_to_pipelines, self.pipelines_list):
                if retrieve_from_file:
                    print(f"Warning... Retrieving from file with method n_comp=512 and rescale=False")
                    with open(os.path.join(path_to_pipeline, "shannon_entropy", "shannon_entropy.json"), 'r') as f:
                        entropy_dict = json.load(f)[group_name_in_file][name_denoised]
                else:
                    entropy_dict = self.image_embeddings[model].get_shannon_entropy_per_group(group=group, denoised=denoised, n_comp=n_comp, rescale=rescale, pct=pct)
                
                entropy_models[model] = entropy_dict
                kl_divergence_models[model] = self.image_embeddings[model].get_kl_divergence_per_group(group=group, denoised=denoised, n_comp=n_comp, pct=pct)

            # Add handcrafted features
            if self.ef is not None:
                entropy_models['handcrafted'] = self.ef.get_shannon_entropy_per_group(group=group, denoised=denoised, n_comp=n_comp, rescale=rescale, pct=pct)
                kl_divergence_models['handcrafted'] = self.ef.get_kl_divergence_per_group(group=group, denoised=denoised, n_comp=n_comp, pct=pct)
            
            print(f"Entropy models: {entropy_models}")
            ## Boxlots
            for dict_, name in zip([entropy_models, kl_divergence_models], ['shannon_entropy', 'kl_divergence']):
                plt.figure(figsize=(2*len(dict_), 5))
                print(f"Dict: {dict_}")
                df = pd.DataFrame.from_dict(dict_)

                df_ = pd.melt(df, var_name='model', value_name=name)
                df_[group] = list(df.index) * len(df.columns)

                plt.figure(figsize=(10, 6))
                sns.boxplot(data=df_, x='model', y=name, hue='model')
                sns.stripplot(data=df_, x='model', y=name, hue=group, 
                            dodge=True, palette=sns.color_palette("colorblind"), edgecolor="black", linewidth=1)
                plt.title(f"Boxplot of {name.replace('_', ' ')} per {group} \n {name_comp.replace('_', ' ')} - {name_denoised} SVD -{name_rescale.replace('_', ' ')}", weight='bold')
                plt.ylabel(name.replace('_', ' ').capitalize())
                plt.xlabel("Model")
                plt.xticks(rotation=90)
                sns.despine()

                # print(f"Saving shannon entropy plot to {os.path.join(self.saving_folder, f"boxplot_{name}_per_{group}_{name_comp}_{name_denoised}{name_rescale}.{self.extension}")}")
                plt.savefig(os.path.join(self.saving_folder, f"boxplot_{name}_per_{group}_{name_comp}_{name_denoised}{name_rescale}.{self.extension}"), bbox_inches='tight')


                ## Barplot
                plt.figure(figsize=(2*len(entropy_dict), 6))
                sns.barplot(data=df_, x=group, y=name, hue='model')
                plt.title(f"Barplot of {name.replace('_', ' ')} per {group} \n {name_comp.replace('_', ' ')} - {name_denoised} SVD -{name_rescale.replace('_', ' ')}", weight='bold')
                plt.ylabel(name.replace('_', ' ').capitalize())
                plt.xlabel(f"{group.capitalize()}")
                sns.despine()
                plt.legend(loc='upper right', bbox_to_anchor=(1.8, 1))
                plt.savefig(os.path.join(self.saving_folder, f"barplot_{name}_per_{group}_{name_comp}_{name_denoised}{name_rescale}.{self.extension}"), bbox_inches='tight')

                with open(os.path.join(self.saving_folder, f"{name}_per_{group}_{name_comp}_{name_denoised}{name_rescale}.json"), 'w') as f:
                    json.dump(dict_, f)
        else:
            entropies = {}
            kl_divergences = {}

            for path_to_pipeline, model in zip(self.path_to_pipelines, self.pipelines_list):
                if retrieve_from_file:
                    print(f"Warning... Retrieving from file with method n_comp=512 and rescale=False")
                    with open(os.path.join(path_to_pipeline, "shannon_entropy", "shannon_entropy.json"), 'r') as f:
                        entropy_dict = json.load(f)['whole_matrix'][name_denoised]
                else:
                    entropies[model] = self.image_embeddings[model].get_shannon_entropy(denoised=denoised, n_comp=n_comp, rescale=rescale, pct=pct)
                self.image_embeddings[model].scree_plot(matrix_name=name_denoised, ylim=(0, 0.8))
                kl_divergences[model] = self.image_embeddings[model].get_kl_divergence(denoised=denoised, n_comp=n_comp, pct=pct)

            # Add handcrafted features
            if self.ef is not None:
                entropies['handcrafted'] = self.ef.get_shannon_entropy(denoised=denoised, n_comp=n_comp, rescale=rescale, pct=pct)
                kl_divergences['handcrafted'] = self.ef.get_kl_divergence(denoised=denoised, n_comp=n_comp, pct=pct)
            
            for dict_, name in zip([entropies, kl_divergences], ['shannon_entropy', 'kl_divergence']):
                plt.figure(figsize=(2*len(self.pipelines_list), 5))
                sns.barplot(x=dict_.keys(), y=dict_.values(), hue=dict_.keys())
                plt.title(f"{name.replace('_', ' ').capitalize()}  \n {name_comp.replace('_', ' ')} - {name_denoised} SVD -{name_rescale.replace('_', ' ')}", weight='bold')
                plt.ylabel(f"{name.replace('_', ' ').capitalize()}")
                sns.despine()
                plt.xticks(rotation=90)
                plt.savefig(os.path.join(self.saving_folder, f"{name}_whole_matrix_{name_comp}_{name_denoised}{name_rescale}.{self.extension}"), bbox_inches='tight')

                with open(os.path.join(self.saving_folder, f"{name}_{name_comp}_{name_denoised}{name_rescale}.json"), 'w') as f:
                    json.dump(dict_, f)

        


    def get_all_shannon_entropy_and_kl_divergence_plots_whole_matrix(self, denoised=False):


        print(f"Plotting Shannon entropy... denoised: {denoised}", flush=True)
        # Shannon entropy on all components
        self.shannon_entropies_and_kl_divergence_all_models(denoised=denoised, group=None, n_comp=None, pct=None)

        # Shannon entropy on the first 512 components
        # self.shannon_entropies_all_models(denoised=denoised, group=None, n_comp=self.min_comp, pct=None, retrieve_from_file=True)

        # Shannon entropy on the first 90% of variance
        # self.shannon_entropies_all_models(denoised=denoised, group=None, n_comp=None, pct=0.9)

        # Shannon entropy on the first 512 components rescaled
        self.shannon_entropies_and_kl_divergence_all_models(denoised=denoised, group=None, n_comp=self.min_comp, pct=None, rescale=True)
    
    def get_all_shannon_entropy_and_kl_divergence_plots_per_group(self, group, denoised=False):

        print(f"Plotting Shannon entropy per group... denoised: {denoised}", flush=True)

        # Shannon entropy on all components
        self.shannon_entropies_and_kl_divergence_all_models(denoised=denoised, group=group, n_comp=None, pct=None)

        # Shannon entropy on the first 512 components
        # self.shannon_entropies_all_models(denoised=denoised, group=group, n_comp=self.min_comp, pct=None, retrieve_from_file=True)

        # Shannon entropy on the first 90% of variance
        # self.shannon_entropies_all_models(denoised=denoised, group=group, n_comp=None, pct=0.9)

        # Shannon entropy on the first 512 components rescaled
       
        self.shannon_entropies_and_kl_divergence_all_models(denoised=denoised, group=group, n_comp=self.min_comp, pct=None, rescale=True)

    def execute_pipeline(self):
        self.get_all_shannon_entropy_and_kl_divergence_plots_whole_matrix(denoised=False)
        self.get_all_shannon_entropy_and_kl_divergence_plots_per_group(denoised=False, group=self.group)
        # benchmark_obj.get_all_shannon_entropy_and_kl_divergence_plots_whole_matrix(denoised=True)
        # benchmark_obj.get_all_shannon_entropy_and_kl_divergence_plots_per_group(denoised=True, group=benchmark_obj.group)
        nb_components = self.nb_components_explaining_variance(pct=self.pct_variance)

        with open(os.path.join(self.saving_folder, f'nb_comp_explaining_{int(np.round(self.pct_variance*100))}%_variance.json'), 'w') as f:
            json.dump(nb_components, f)
    