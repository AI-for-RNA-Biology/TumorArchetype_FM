#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#

import os
import numpy as np
import anndata
from digitalhistopathology.embeddings.image_embedding import ImageEmbedding
from digitalhistopathology.engineered_features.engineered_features import EngineeredFeatures, scMTOP_EngineeredFeatures
from digitalhistopathology.embeddings.embedding import Embedding
from digitalhistopathology.datasets.real_datasets import HER2Dataset, TNBCDataset
import glob

class BenchmarkBase:
    """
    Base class for all benchmarks.
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
                 dataset="HER2",
                 ):
        

        if len(path_to_pipeline) == 1 and len(pipelines_list) > 1:
            self.pipelines_list = pipelines_list
            self.path_to_pipelines = [os.path.join(path_to_pipeline, model) for model in pipelines_list]
        else:
            self.path_to_pipelines = path_to_pipeline
            print(pipelines_list)
            print(type(pipelines_list))
            self.pipelines_list = pipelines_list
        

        self.dataset_name = dataset_name
        self.image_embedding_name = image_embedding_name
        self.results_folder = results_folder
        self.image_embeddings = {}


        self.saving_folder = saving_folder
        self.engineered_features_type = engineered_features_type
        self.emb_df_csv_path = emb_df_csv_path
        self.engineered_features_saving_folder = engineered_features_saving_folder
        self.label_files = None
        self.extension = extension

        if not os.path.exists(self.saving_folder):
            print(f"Creating folder {self.saving_folder}...", flush=True)
            os.makedirs(self.saving_folder, exist_ok=True)


        self.patches_filenames = sorted(glob.glob(os.path.join(self.results_folder, "{}/compute_patches/all/*.tiff".format(self.dataset_name))))
        self.patches_info_filename = os.path.join(self.results_folder, "{}/compute_patches/all/patches_info.pkl.gz".format(self.dataset_name))
        
        if dataset == "HER2":
            self.dataset = HER2Dataset()
        elif dataset == "TNBC":
            self.dataset = TNBCDataset()
        else:
            raise ValueError(f"Dataset {dataset} not recognized. Please use 'HER2' or 'TNBC'.")
        
        self.label_files = self.dataset.label_filenames
        print(f"Label files: {self.label_files}")
        self.group = group
        self.embeddings_per_slide = None
        self.annotated_embeddings = None
        self.ef = None
        
    def initialize_model_saving_folders(self):
        """
        Creates directories for saving model outputs if they do not already exist.
        This method iterates over the list of models specified in `self.pipelines_list`
        and checks if a directory for each model exists within the `self.saving_folder`.
        If a directory does not exist, it creates one.
        Returns:
            None
        """

        for model in self.pipelines_list:
            if not os.path.exists(os.path.join(self.saving_folder, model)):
                os.makedirs(os.path.join(self.saving_folder, model))
    

    def load_engineered_features(self):
        """
        Load engineered features based on the specified type.
        This method loads engineered features according to the type specified in 
        `self.engineered_features_type`. It supports 'scMTOP'. The loaded features are then processed to ensure they are 
        ordered with the same indexes as the image embeddings and any NaN values 
        are handled.
        Attributes:
            self.engineered_features_type (str): Type of engineered features to load.
            self.patches_info_filename (str): Filename containing patches information.
            self.emb_df_csv_path (str): Path to the CSV file containing embedding data.
            self.dataset_name (str): Name of the dataset.
            self.engineered_features_saving_folder (str): Folder to save the engineered features.
            self.label_files (list): List of label files.
            self.name_engineered_feature (str): Name of the engineered feature.
            self.image_embeddings (dict): Dictionary containing image embeddings.
            self.pipelines_list (list): List of pipelines.
        Raises:
            ValueError: If `self.engineered_features_type` is not recognized.
        Returns:
            None
        """

        if self.engineered_features_saving_folder is None:
            self.ef = None
        else:
            if self.engineered_features_type == 'scMTOP':
                # We load the engineered features
                print("Loading scMTOP engineered features...", flush=True)
                ef = scMTOP_EngineeredFeatures(patches_info_filename=self.patches_info_filename,
                                    saving_plots=True,
                                    emb_df_csv=self.emb_df_csv_path,
                                    dataset_name=self.dataset_name,
                                    result_saving_folder=self.engineered_features_saving_folder,
                                    label_files=self.label_files)
                
            else:
                ef = EngineeredFeatures(patches_info_filename=self.patches_info_filename,
                                        result_saving_folder=self.engineered_features_saving_folder,
                                        saving_plots=True,
                                        emb_df_csv=self.emb_df_csv_path,
                                        dataset_name=self.dataset_name, 
                                        label_files=self.label_files)
            ef.load_emb_df()
            ef.fill_emb()


            self.ef = ef

            # Check that the engineered features are ordered with the same indexes
            indexes = self.image_embeddings[self.pipelines_list[0]].emb.obs_names
            self.ef.emb = self.ef.emb[indexes,:]

            # Remove nans
            df = self.ef.emb.to_df()
            df.replace([np.inf, -np.inf], np.nan, inplace=True)
            df.fillna(df.mean(), inplace=True)
            self.ef.emb.X = df

    
    def compute_image_embeddings(self):
        """
        Computes image embeddings for each model in the pipelines list and stores them in the instance variable `image_embeddings`.
        Parameters:
        denoised (bool): If True, applies denoising during the SVD computation of the image embeddings. Default is False.
        Returns:
        None
        This method performs the following steps:
        1. Iterates over each model in the `pipelines_list`.
        2. Reads the image embedding data from an .h5ad file.
        3. Initializes an `ImageEmbedding` object with the read data and other relevant parameters.
        4. Computes the SVD of the image embeddings, optionally applying denoising.
        5. Stores the computed `ImageEmbedding` object in the `image_embeddings` dictionary with the model name as the key.
        6. Prints a message indicating that SVD computation is complete for each model.
        7. Plots the scree plot for each model.
        """

        
        image_embeddings = {}
        for path_to_pipeline, model in zip(self.path_to_pipelines, self.pipelines_list):
            print(f"Loading image embeddings for model {model}... at {os.path.join(path_to_pipeline, self.image_embedding_name)}", flush=True)
            emb = anndata.read_h5ad(os.path.join(path_to_pipeline, self.image_embedding_name))
            image_emb = ImageEmbedding(patches_filenames=self.patches_filenames, 
                                       patches_info_filename=self.patches_info_filename, 
                                       emb=emb, 
                                       # result_saving_folder=os.path.join(self.saving_folder, model),
                                       result_saving_folder=os.path.join(self.saving_folder),
                                       saving_plots=True,
                                       label_files=self.label_files)            
            image_embeddings[model] = image_emb

        self.image_embeddings = image_embeddings

        # Check that the image embeddings are all ordered with the same indexes
        indexes = self.image_embeddings[self.pipelines_list[0]].emb.obs_names


        for model in self.pipelines_list[1:]:
            self.image_embeddings[model].emb = self.image_embeddings[model].emb[indexes,:]


    def compute_svd_image_embeddings(self, denoised=False, scale=False, center=False):
        """
        Compute Singular Value Decomposition (SVD) for image embeddings.
        This method computes the SVD for image embeddings stored in the `image_embeddings` attribute
        for each model in the `pipelines_list`. The SVD can be computed on either the raw or denoised
        data, and can optionally be centered and/or scaled.
        Parameters:
        denoised (bool): If True, use denoised data for SVD computation. Default is False.
        scale (bool): If True, scale the data before SVD computation. Default is False.
        center (bool): If True, center the data before SVD computation. Default is False.
        Returns:
        None
        """

        if denoised:
            matrix_name = "denoised"
        else:
            matrix_name = "raw"
    
        for model in self.pipelines_list:
            self.image_embeddings[model].compute_svd(denoised=denoised, center=center, scale=scale)
            print(f"SVD computed for image embeddings: {model}", flush=True)
            self.image_embeddings[model].scree_plot(matrix_name=matrix_name)

    def get_annotated_embeddings(self):
        annotated_embeddings = {}
        for model, whole_emb in self.image_embeddings.items():
            subset_emb = ImageEmbedding()
            subset_emb.emb = whole_emb.emb[~whole_emb.emb.obs['label'].isna()]
            subset_emb.emb = subset_emb.emb[subset_emb.emb.obs['label'] != 'nan']
            annotated_embeddings[model] = subset_emb
        
        self.annotated_embeddings = annotated_embeddings

        if self.ef is not None:
            subset_emb = Embedding()
            subset_emb.emb = self.ef.emb[~self.ef.emb.obs['label'].isna()]
            subset_emb.emb = subset_emb.emb[subset_emb.emb.obs['label'] != 'nan']
            self.annotated_embeddings['handcrafted_features'] = subset_emb


    def get_embeddings_per_slide(self, slide_id_col='name_origin', only_labeled=False):
        embeddings_per_slide = {}
        for model in self.pipelines_list:
            embeddings_per_slide[model] = {}
            for slide in self.image_embeddings[model].emb.obs[slide_id_col].unique():
                print(f"Getting embeddings for patient {slide} for model {model}...", flush=True)

                if "label" not in self.image_embeddings[model].emb.obs:
                    self.image_embeddings[model].add_label()

                subset_emb = ImageEmbedding()
                subset_emb.emb = self.image_embeddings[model].emb[self.image_embeddings[model].emb.obs[slide_id_col] == slide]
                if only_labeled:
                    subset_emb.emb = subset_emb.emb[~subset_emb.emb.obs['label'].isna()]
                    subset_emb.emb = subset_emb.emb[subset_emb.emb.obs['label'] != 'nan']

                embeddings_per_slide[model][slide] = subset_emb
        
        # Do it also for handcrafted features
        if self.ef is not None:
            embeddings_per_slide['handcrafted_features'] = {}
            for slide in self.ef.emb.obs[slide_id_col].unique():

                subset_emb = ImageEmbedding()
                subset_emb.emb = self.ef.emb[self.ef.emb.obs[slide_id_col] == slide]

                if only_labeled:
                    subset_emb.emb = subset_emb.emb[~subset_emb.emb.obs['label'].isna()]
                    subset_emb.emb = subset_emb.emb[subset_emb.emb.obs['label'] != 'nan']

                embeddings_per_slide['handcrafted_features'][slide] = subset_emb

        self.embeddings_per_slide = embeddings_per_slide

    def load_invasive_image_embeddings(self):
        image_embeddings = {}

        for path_to_pipeline, model in zip(self.path_to_pipelines, self.pipelines_list):
            emb = anndata.read_h5ad(os.path.join(path_to_pipeline, "invasive_image_embedding.h5ad"))
            image_emb = ImageEmbedding(emb=emb, 
                                       result_saving_folder=os.path.join(self.saving_folder, model),
                                       label_files=self.label_files,
                                       patches_info_filename=self.patches_info_filename)
            # image_emb.compute_svd(denoised=denoised)
            
            image_embeddings[model] = image_emb
            # print(f"SVD computed for image embeddings: {model}", flush=True)
            # image_emb.scree_plot(matrix_name=matrix_name)

        self.invasive_image_embeddings = image_embeddings

    def execute_pipeline(self):
        pass