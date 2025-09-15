#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import gzip
import pickle
import random
import json
import glob
import argparse

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import torch
import sys

sys.path.append("../")
import pyvips 
from pathlib import Path
import warnings
import openslide
import mahotas
import cv2



from skimage.color import rgb2gray
from skimage.feature import graycomatrix, graycoprops

from skimage.measure import shannon_entropy
from scipy.stats import spearmanr
from scipy.stats import skew, kurtosis
import mahotas as mh

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

from cellpose import models
from cellpose.io import imread
from PIL import Image, ImageDraw
from skimage.measure import regionprops_table
import skimage.io as io
from skimage.measure import label, regionprops
from skimage.color import rgb2gray
import datetime
from tqdm import tqdm


from multiprocessing import Pool


Image.MAX_IMAGE_PIXELS = None

from digitalhistopathology.embeddings.image_embedding import ImageEmbedding
from digitalhistopathology.engineered_features.engineered_utils import (
    apply_mask,
    get_whole_image_texture_features,
    get_color_features,
    get_patch_zernike_moments,
    get_zernike_moments_per_nuclei
)

import os

CELL_TYPE_DICT = {1: "T", 2: "I", 3: "S", 5: "N"}


class EngineeredFeatures(ImageEmbedding):
    def __init__(self,
                 saving_plots=False,
                 label_files=None,
                 emb=None,
                 result_saving_folder="../results/engineered_features",
                 name="",
                 emb_df_csv=None,
                 patches_info_filename=None,
                 dataset_name="dataset", 
                 feature_type_color_dict=None):
        """EngineeredFeatures class contains dimensionally reduction, clustering and visualization techniques to analyze engineered features embeddings. It inherits from Embeddings class.

        Args:
            emb (anndata, optional): Embeddings anndata object. Defaults to None.
            result_saving_folder (str, optional): Result folder in which the results are saved. Defaults to "../results".
            name (str, optional): Name of the embeddings. Defaults to "".
            saving_plots (bool, optional): If the plots are saved to the result folder or not. Defaults to False.
            label_files (list, optional): List of files containing label of each spot, "x" column corresponds to the first spot coordinate, "y" column corresponds to the second. Can be csv, tsv with gzip compression or not. One per sample or the name of the file contain the sample name at the beginning with a "_" just after. Defaults to None.
            emb_df (pd.DataFrame, optional): DataFrame containing embeddings. Defaults to None.
        """
        ImageEmbedding.__init__(self,
                                emb=emb, 
                                result_saving_folder=result_saving_folder, 
                                name=name, saving_plots=saving_plots, 
                                label_files=label_files,
                                patches_info_filename=patches_info_filename,)
        self.dataset_name = dataset_name

        self.emb_df_csv = emb_df_csv
        self.feature_type_color_dict = feature_type_color_dict

        if self.result_saving_folder is not None:
            if not os.path.exists(self.result_saving_folder):
                os.makedirs(self.result_saving_folder)
            
            self.save_path = os.path.join(self.result_saving_folder, self.dataset_name)

            if not os.path.exists(self.save_path):
                os.makedirs(self.save_path)
                
        else:
            self.save_path = None

        if patches_info_filename is not None:
            with gzip.open(patches_info_filename) as file:
                self.patches_info = pickle.load(file)
        else:
            self.patches_info = None


    def fill_emb(self):
        """Fill emb anndata with emb_df attribute.

        Raises:
            Exception: Attribute emb_df is empty
        """
        print("Filling emb...")
        try:
            if len(self.emb_df) > 0:
                print("Emb_df is not empty.")
                self.emb = ad.AnnData(self.emb_df)
                if self.patches_info_filename is not None:
                    print("Loading patches info...")
                    patches_info_df = self.load_patches_infos()
                    self.emb.obs = self.emb.obs.merge(patches_info_df, right_index=True, left_index=True)
            else:
                self.emb_df = None
        except:
            self.emb_df = None

        
    
    def load_emb_df(self):
        """Load emb_df from a file.

        Args:
            filename (str): Filename to load the dataframe.
        """
        try:
            self.emb_df = pd.read_csv(self.emb_df_csv, index_col=0)
        except:
            print("File storing the emf_df doesn't exist - You need to compute the engineered features first.")
            self.emb_df = None
        # self.emb_df['detection_mask'] = self.emb_df['detection_mask'].apply(lambda x: pickle.loads(eval(x)))

    def save_emb(self):
        """Save emb to a .h5ad file.
        """
        self.save_embeddings(saving_path=os.path.join(self.save_path, self.dataset_name + "_scMTOP.h5ad"))

    def compute_explanation_scores(self, pct=0.9, pcs=None):
        df = self.emb.to_df()
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        df.fillna(df.mean(), inplace=True)
        self.emb.X = df

        print("Computing SVD...")
        self.compute_svd(denoised=False, center=True, scale=True)

        correlations = {}
        correlations['stat'] = {}
        correlations['p_value'] = {}

        engineered_feature_names = list(self.emb_df.columns)

        if pct is None:
            if pcs is None or pcs > self.svd['U_df'].shape[1]:
                pcs = self.svd['U_df'].shape[1]
        else:
            pcs = np.argmax(np.cumsum(self.get_explained_variance_ratio_list()) > pct)

        for pc in range(pcs):
            correlations['stat'][f"PC{pc+1}"] = {}
            correlations['p_value'][f"PC{pc+1}"] = {}
            pc_values = self.svd['U_df'].loc[list(self.emb_df.index)][f"u{pc+1}"].values
            for engineered_feature_name in engineered_feature_names:
                s, p = spearmanr(pc_values, self.emb_df[engineered_feature_name], nan_policy='omit')
                correlations['stat'][f"PC{pc+1}"][engineered_feature_name] = s
                correlations['p_value'][f"PC{pc+1}"][engineered_feature_name] = p

        explained_variance_ratio = self.get_explained_variance_ratio_list()
        vars = explained_variance_ratio[:pcs]

        explanation_scores_handcrafted = pd.DataFrame(abs(pd.DataFrame(correlations['stat'])).mul(vars, axis=1).sum(axis=1), columns=['explanation_score'])
        explanation_scores_handcrafted['feature_type'] = [idx.split('_')[0] for idx in explanation_scores_handcrafted.index]
        return explanation_scores_handcrafted
    


class scMTOP_EngineeredFeatures(EngineeredFeatures):
    """
    scMTOP_EngineeredFeatures is a class that extends the Embedding class to handle engineered features for scMTOP.
    Attributes:
        path_to_wsis (str): Path to whole slide images.
        list_wsi_filenames (list): List of WSI filenames.
        patches_info_filename (str): Filename for patches information.
        path_to_cellvit_folder (str): Path to the CellVit folder.
        result_saving_folder (str): Folder to save results.
        emb_df (pd.DataFrame): DataFrame containing embeddings.
        dataset_name (str): Name of the dataset.
        temporary_folder (str): Temporary folder for intermediate files.
        wsi_names (list): List of WSI names.
        save_path (str): Path to save the results.
    Methods:
        __init__(self, emb=None, result_saving_folder="../results/engineered_features/scMTOP", name="", saving_plots=False, label_files=None, patches_info_filename=None, path_to_cellvit_folder=None, list_wsi_filenames=None, path_to_wsis=None, emb_df=None, dataset_name="dataset", temporary_folder="../results/tmp"):
            Initializes the scMTOP_EngineeredFeatures class with the given parameters.
        convertjson_from_cellvit_to_scMTOP(self, cellvit_json_filename):
            Converts JSON files from CellVit format to scMTOP format.
        check_and_convert_wsifile(self, path_to_wsi):
            Checks and converts WSI files to the appropriate format.
        compute_scMTOP_features(self):
            Computes scMTOP features for the given WSI files.
        map_cells_and_patches(self, wsi_name):
            Maps cells and patches for a given WSI name.
        get_patch_features(self, wsi_name, patch_to_cell):
            Gets patch features for a given WSI name and patch-to-cell mapping.
        get_whole_scMTOP_embedding(self, save_individual_wsi=False, save_whole=True):
            Computes and saves the whole scMTOP embedding.
        load_emb_df_from_individual_wsi(self):
            Loads embedding DataFrame from individual WSI files.
        load_emb_df_from_whole(self):
            Loads embedding DataFrame from the whole dataset.
        fill_emb(self):
            Fills the embedding attribute with the DataFrame.
    """
    
    def __init__(self, 
                 emb=None, 
                 result_saving_folder="../results/engineered_features/scMTOP", 
                 name="", 
                 saving_plots=False, 
                 label_files=None,
                 patches_info_filename=None,
                 path_to_cellvit_folder=None,
                 list_wsi_filenames=None,
                 path_to_wsis=None,
                 emb_df_csv=None,
                 dataset_name="dataset",
                 temporary_folder="../results/tmp",
                 feature_type_color_dict=None):
        """Initializes the scMTOP_EngineeredFeatures class with the given parameters.

        Args:
            emb (anndata, optional): Embeddings anndata object. Defaults to None.
            result_saving_folder (str, optional): Folder to save results. Defaults to "../results/engineered_features/scMTOP".
            name (str, optional): Name of the embeddings. Defaults to "".
            saving_plots (bool, optional): If the plots are saved to the result folder or not. Defaults to False.
            label_files (list, optional): List of files containing label of each spot, "x" column corresponds to the first spot coordinate, "y" column corresponds to the second. Can be csv, tsv with gzip compression or not. One per sample or the name of the file contain the sample name at the beginning with a "_" just after. Defaults to None.
            patches_info_filename (str, optional): Filename for patches information. Defaults to None.
            path_to_cellvit_folder (str, optional): Path to the CellVit folder. Defaults to None.
            list_wsi_filenames (list, optional): List of WSI filenames. Defaults to None.
            path_to_wsis (str, optional): Path to whole slide images. Defaults to None.
            emb_df (pd.DataFrame, optional): DataFrame containing embeddings. Defaults to None.
            dataset_name (str, optional): Name of the dataset. Defaults to "dataset".
            temporary_folder (str, optional): Temporary folder for intermediate files. Defaults to "../results/tmp".
        """
        
        
        super().__init__(
            emb=emb, 
            result_saving_folder=result_saving_folder, 
            name=name, 
            saving_plots=saving_plots, 
            label_files=label_files,
            dataset_name=dataset_name,
            patches_info_filename=patches_info_filename,
            feature_type_color_dict=feature_type_color_dict
        )
        

        self.path_to_wsis = path_to_wsis
        self.list_wsi_filenames = list_wsi_filenames
        self.patches_info_filename = patches_info_filename

        if self.path_to_wsis is not None:
            self.list_wsi_filenames = glob.glob(os.path.join(path_to_wsis, '*'))
        elif self.list_wsi_filenames is not None:
            self.list_wsi_filenames = list_wsi_filenames
        elif self.patches_info is not None:
            self.list_wsi_filenames = list(set([patch['path_origin'] for patch in self.patches_info]))
            
        else:
            warnings.warn("path_to_wsis, list_wsi_filenames and patches_info_filename cannot be all None. If you want \
                          to compute the scMTOP features, you need to provide at least one of these parameters.")

        self.temporary_folder = temporary_folder
        self.path_to_cellvit_folder = path_to_cellvit_folder
        self.path_to_wsis = path_to_wsis
        self.result_saving_folder =result_saving_folder
        self.emb_df_csv = emb_df_csv

        if self.list_wsi_filenames is not None:
            self.list_wsi_filenames = [self.check_and_convert_wsifile(path_to_wsi) for path_to_wsi in self.list_wsi_filenames]  
            self.wsi_names = [Path(path_to_wsi).stem for path_to_wsi in self.list_wsi_filenames]
            print(self.wsi_names, flush=True)
        else:
            self.wsi_names = None


        self.load_emb_df()


        if self.path_to_cellvit_folder is not None:
            cellvit_paths = glob.glob(os.path.join(self.path_to_cellvit_folder, '*'))


            wsi_cell_vit_paths = []
            for wsi_name in self.wsi_names:
                for cellvit_path in cellvit_paths:
                    if wsi_name in os.listdir(cellvit_path):
                        wsi_cell_vit_paths.append(cellvit_path)
                        break
            self.wsi_cellvit_paths = wsi_cell_vit_paths

        self.feature_type_color_dict = {'Nuclei-Composition': 'mediumpurple', 
                                        'Nuclei-Morph': 'tomato', 
                                        'Nuclei-Texture': 'lightgreen', 
                                        'ExtraCell-Morph': 'darkorange', 
                                        'ExtraCell-Color': 'lightpink', 
                                        'WholePatch-Texture': 'brown',
                                        'WholePatch-Morph': 'lightblue'}
        
        if self.save_path is not None:
            if not os.path.exists(os.path.join(self.save_path, "zernike_cell")):
                os.makedirs(os.path.join(self.save_path, "zernike_cells"), exist_ok=True)
            
            if not os.path.exists(os.path.join(self.save_path, "patch_to_cell")):
                os.makedirs(os.path.join(self.save_path, "patch_to_cell"), exist_ok=True)
                              

    def convertjson_from_cellvit_to_scMTOP(self, cellvit_json_filename):
        """
        Converts a JSON file from the CellViT format to the scMTOP format.
        This method reads a JSON file containing cell data in the CellViT format,
        extracts the relevant information, and writes it to a new JSON file in the
        scMTOP format.
        As a reminder the CellViT format is as follows:
        {   "wsi_metadata": {
                "magnification": 40
            },
            "cells": [
                {
                    "cell_id": 0,
                    "centroid": [100, 200],
                    "bbox": [90, 190, 110, 210],
                    "cell_type": 1
                },
                ...
            ]
        }
        The scMTOP format is as follows:
        {  "mag": 40,
            "nuc": {
                0: {
                    "centroid": [100, 200],
                    "bbox": [90, 190, 110, 210],
                    "cell_type": 1
                },
                ...
            }
        }
        Args:
            cellvit_json_filename (str): The path to the input JSON file in CellViT format.
        Returns:
            None
        """

        dirname = os.path.dirname(cellvit_json_filename)
        with open(cellvit_json_filename, 'r') as f:
            cellvit_json = json.load(f)

        mag = cellvit_json['wsi_metadata']['magnification']
        nuc = {}
        for i, item in enumerate(cellvit_json['cells']):
            nuc[i] = item
        
        scMTOP_json = {'mag': mag, 'nuc': nuc}

        with open(os.path.join(dirname, 'cells_for_scMTOP.json'), 'w') as f:
            json.dump(scMTOP_json, f)
    
    def check_and_convert_wsifile(self, path_to_wsi):
        """
        Checks the WSI file extension and converts it to the appropriate format if necessary.
        Args:
            path_to_wsi (str): The path to the WSI file.
        Returns:
            str: The path to the WSI file in the appropriate format.
        Raises:
            ValueError: If the WSI file extension is not supported.
        """

        
        extension = os.path.basename(path_to_wsi).split('.')[-1]
        filename = os.path.basename(path_to_wsi).split(f'.{extension}')[0]

        if extension in ["svs",
                         "tiff",
                         "tif",
                         "bif",
                         "scn",
                         "ndpi",
                         "vms",
                         "vmu"]:
            return path_to_wsi
        elif extension == 'jpg':
            image = pyvips.Image.new_from_file(path_to_wsi, access='sequential')
            path_new_tif = os.path.join(self.temporary_folder, self.dataset_name, "wsi", f"{filename}.tif")

            if not os.path.exists(os.path.dirname(path_new_tif)):
                os.makedirs(os.path.dirname(path_new_tif))

            image.tiffsave(path_new_tif,
                           compression="jpeg",
                           Q=75,
                           tile=True,
                           tile_width=256,
                           tile_height=256,
                           pyramid=True)
            return path_new_tif
        else:
            raise ValueError("The WSI file extension is not supported.")
    
    def compute_scMTOP_features(self):
        """
        Computes scMTOP features for the given WSI files.
        This method reads the JSON files containing cell data in the CellViT format,
        converts them to the scMTOP format, and computes the scMTOP features.
        The scMTOP features are saved in the specified directory.
        Args:
            None
        Returns:
            None
        """

        from sc_MTOP.F3_FeatureExtract import fun3

        ## Convert the json files
        for wsi_name, wsi_path, wsi_cellvit_path in zip(self.wsi_names, self.list_wsi_filenames, self.wsi_cellvit_paths):
            print(f"Computing scMTOP features for {wsi_name}")
            cellvit_json_filename = os.path.join(wsi_cellvit_path, wsi_name, 'cell_detection', 'cells.json')
            self.convertjson_from_cellvit_to_scMTOP(cellvit_json_filename)
            fun3(os.path.join(wsi_cellvit_path, wsi_name, 'cell_detection', 'cells_for_scMTOP.json'), wsi_path, self.save_path)
    
    def map_cells_and_patches(self, wsi_name, wsi_cellvit_path):
        """
        Maps cells and patches for a given WSI name.
        This method reads the JSON file containing cell data in the scMTOP format,
        extracts the relevant information, and maps the cells to the patches.
        Args:
            wsi_name (str): The name of the WSI.
        Returns:
            dict: A dictionary mapping patches to cells.
        """

        path_to_scMTOP_json = os.path.join(wsi_cellvit_path, wsi_name, 'cell_detection', 'cells_for_scMTOP.json')
        with open(path_to_scMTOP_json, 'r') as f:
            scMTOP_json = json.load(f)
        
        nuc = scMTOP_json['nuc']
        patch_to_cell = {}

        print(f"Length of the scMTOP json: {len(nuc)}")

        # patches_wsi = [patch for patch in patches_info if os.path.basename(patch['path_origin'].split('.')[0]) == wsi_name]
        patches_wsi = []
        filenames = []
        for patch in self.patches_info:
            extension = os.path.basename(patch['path_origin']).split('.')[-1]
            filename = os.path.basename(patch['path_origin']).split(f'.{extension}')[0]
            filenames.append(filename)
            if filename == wsi_name:
                patches_wsi.append(patch)
        print(f"Number of patches in {wsi_name}: {len(patches_wsi)}")
        print(list(set(filenames)))

        df_nuc = pd.DataFrame.from_dict(nuc, orient='index')

        df_nuc = df_nuc[['centroid']]
        df_nuc['centroid_x'] = df_nuc['centroid'].apply(lambda x: x[0])
        df_nuc['centroid_y'] = df_nuc['centroid'].apply(lambda x: x[1])


        for patch in patches_wsi:
            patch_shape = patch['shape_pixel']
            start_width_origin = patch['start_width_origin']
            start_height_origin = patch['start_height_origin']
            end_width_origin = start_width_origin + patch_shape
            end_height_origin = start_height_origin + patch_shape

            cells = df_nuc[(df_nuc['centroid_x'] >= start_width_origin) & \
                    (df_nuc['centroid_x'] <= end_width_origin) & \
                    (df_nuc['centroid_y'] >= start_height_origin) & \
                    (df_nuc['centroid_y'] <= end_height_origin)]
            
            patch_to_cell[patch['name']] = cells.index.tolist()

        return patch_to_cell
    
    def process_patch_features(patch, 
                                patch_to_cell, 
                                all_cells, 
                                self, 
                                wsi_path,
                                wsi_name,
                                wsi_cellvit_path, 
                                zernike_cells=None):
        emb_df = pd.DataFrame()

        print(f"Computing patch features for {patch}...", flush=True)

        # Mean and std of the features of the cells in the patch
        cells_ids = [int(cell_id) for cell_id in patch_to_cell[patch]]

        # Check for dead cells
        cells_dead = [i for i in cells_ids if i not in all_cells.index]
        if len(cells_dead) > 0:
            print(f"Warning: there is {len(cells_dead)} dead cells in patch {patch}. These cells are ignored.", flush=True)

        # Keep only cells not dead
        cells_ids = [i for i in cells_ids if i in all_cells.index]

        cells_features = all_cells.loc[cells_ids]
        numeric_feats = [col for col in cells_features.columns if col not in ['Bbox', 'Centroid', 'CellType']]
        mean_columns = [f"{col}_mean" for col in numeric_feats]
        emb_df.loc[patch, mean_columns] = cells_features[numeric_feats].mean().values
        std_columns = [f"{col}_std" for col in numeric_feats]
        emb_df.loc[patch, std_columns] = cells_features[numeric_feats].std().values

        # Number of cells in the patch
        number_of_cells = len(cells_ids)
        emb_df.loc[patch, 'Nuclei-Composition_number_of_cells'] = number_of_cells

        # Density of cells in the patch
        shape_patch = [p for p in self.patches_info if p['name'] == patch][0]['shape_pixel']
        cell_density = number_of_cells / (shape_patch ** 2)
        emb_df.loc[patch, 'Nuclei-Composition_total_cell_density'] = cell_density

        # Number, proportion and density of cells per type in the patch
        number_of_cells_per_type = cells_features['CellType'].value_counts()
        for id, cell_type in CELL_TYPE_DICT.items():
            if id in number_of_cells_per_type.index:
                emb_df.loc[patch, f"Nuclei-Composition_prop_of_{cell_type}_cells"] = number_of_cells_per_type[id] / number_of_cells
                emb_df.loc[patch, f"Nuclei-Composition_density_of_{cell_type}_cells"] = number_of_cells_per_type[id] / (shape_patch ** 2)
            else:
                emb_df.loc[patch, f"Nuclei-Composition_prop_of_{cell_type}_cells"] = 0
                emb_df.loc[patch, f"Nuclei-Composition_density_of_{cell_type}_cells"] = 0

        # Extracellular matrix
        non_cell_pixels, patch_mask, patch_img = self.get_extracellular_pixels(wsi_path=wsi_path,
                                                                            wsi_name=wsi_name,
                                                                            wsi_cellvit_path=wsi_cellvit_path,
                                                                            patch_name=patch)
        
        patch_pixels = np.array(patch_img)[:, :, :3].reshape(-1, 3)
        whole_patch_color_features = get_color_features(patch_pixels)
        
        for feature_name, feature_value in whole_patch_color_features.items():
            feature_name = f"WholePatch-Color_{feature_name}"
            emb_df.loc[patch, feature_name] = feature_value
        
        

        # Extracellular matrix features -- color and texture
        extracellular_color_features = get_color_features(non_cell_pixels)
        extracellular_texture_features = self.get_extracellular_texture_features(patch_img, patch_mask)

        for feature_name, feature_value in extracellular_color_features.items():
            feature_name = f"ExtraCell-Color_{feature_name}"
            emb_df.loc[patch, feature_name] = feature_value

        for feature_name, feature_value in extracellular_texture_features.items():
            feature_name = f"ExtraCell-Texture_{feature_name}"
            emb_df.loc[patch, feature_name] = feature_value

        # Whole patch features
        whole_patch_texture_features = get_whole_image_texture_features(patch_img)

        for feature_name, feature_value in whole_patch_texture_features.items():
            feature_name = f"WholePatch-Texture_{feature_name}"
            emb_df.loc[patch, feature_name] = feature_value

        if zernike_cells is not None:
            nuclei_patch_zernike_features = get_patch_zernike_moments(patch, patch_to_cell, zernike_cells)

            for feature_name, feature_value in nuclei_patch_zernike_features.items():
                feature_name = f"Nuclei-Morph_{feature_name}"
                emb_df.loc[patch, feature_name] = feature_value

        return emb_df

    @staticmethod
    def process_patch_features_wrapper(args):
        return scMTOP_EngineeredFeatures.process_patch_features(*args)

    def get_patch_features(self, wsi_path, wsi_name, wsi_cellvit_path, patch_to_cell, zernike_cells=None, num_cores=None):
        print(f"Computing patch features for {wsi_name}...", flush=True)
        T_cells = pd.read_csv(os.path.join(self.save_path, wsi_name, f"{wsi_name}_Feats_T.csv"), index_col=0)
        I_cells = pd.read_csv(os.path.join(self.save_path, wsi_name, f"{wsi_name}_Feats_I.csv"), index_col=0)
        S_cells = pd.read_csv(os.path.join(self.save_path, wsi_name, f"{wsi_name}_Feats_S.csv"), index_col=0)
        N_cells = pd.read_csv(os.path.join(self.save_path, wsi_name, f"{wsi_name}_Feats_N.csv"), index_col=0)

        all_cells = pd.concat([T_cells, I_cells, S_cells, N_cells])

        # Remove graph features
        all_cells = all_cells[[col for col in all_cells.columns if 'Graph' not in col]]

        # Remove stromaBlocker features
        all_cells = all_cells[[col for col in all_cells.columns if 'stromaBlocker' not in col]]

        patches = list(patch_to_cell.keys())
        args = [(patch, patch_to_cell, all_cells, self, wsi_path, wsi_name, wsi_cellvit_path, zernike_cells) for patch in patches]

        with Pool(processes=num_cores) as pool:
            results = list(tqdm(pool.imap(scMTOP_EngineeredFeatures.process_patch_features_wrapper, args), total=len(patches)))

        emb_df = pd.concat(results, axis=0)

        #color_cols = [col for col in emb_df if "Color" in col or "Transparency" in col]
        
        
        emb_df.rename(columns=lambda x: x.replace("Morph", "Nuclei-Morph") if x.startswith("Morph") else x, inplace=True)
        emb_df.rename(columns=lambda x: x.replace("Texture", "Nuclei-Texture") if x.startswith("Texture") else x, inplace=True)
        
        color_cols = [col for col in emb_df if "Color" in col or "Transparency" in col]
        emb_df.rename(columns=lambda x: x.replace("Nuclei-Texture", "Nuclei-Color") if x in color_cols else x, inplace=True)

        return emb_df
        
    

    def get_extracellular_pixels(self, wsi_path, wsi_cellvit_path, wsi_name, patch_name):

        # Load cellvit data
        path_to_scMTOP_json = os.path.join(wsi_cellvit_path, wsi_name, 'cell_detection', 'cells_for_scMTOP.json')

        with open(path_to_scMTOP_json, 'r') as f:
            scMTOP_json = json.load(f)

        # Load wsi
        wsi = openslide.OpenSlide(wsi_path)

        # Load patches info

        patch_info = [patch_info for patch_info in self.patches_info if patch_info['name'] == patch_name][0]
        patch_shape = patch_info['shape_pixel']
        start_width_origin = patch_info['start_width_origin']
        start_height_origin = patch_info['start_height_origin']

        # Get the mask
        mask = Image.new('L', wsi.level_dimensions[0], 255)
        draw = ImageDraw.Draw(mask)
        for cell in scMTOP_json['nuc']:

            draw.polygon([tuple(point) for point in scMTOP_json['nuc'][cell]['contour']], outline=1, fill=1)

        mask = np.array(mask)

        patch_mask = mask[start_height_origin:start_height_origin+patch_shape, start_width_origin:start_width_origin+patch_shape]

        patch_img = wsi.read_region((start_width_origin, start_height_origin), 0, (patch_shape, patch_shape))

        non_cell_pixels = np.array(patch_img)[patch_mask == 255]

        return non_cell_pixels, patch_mask, patch_img
    
    def get_extracellular_texture_features(self, patch_img, patch_mask):
        
        gray_image = (rgb2gray(np.array(patch_img)[:, :, :3]) * 255).astype(np.uint8)
        gray_image[patch_mask == 1] = gray_image.min()
        
        grayco_features = {}
        angles = [0, np.pi/4, np.pi/2, 3*np.pi/4]

        # Compute the gray-level co-occurrence matrix (GLCM)
        glcm = graycomatrix(gray_image, [1], angles, 256, symmetric=True, normed=True)

        properties = ['ASM', 'contrast', 'dissimilarity', 'homogeneity', 'energy', 'correlation']
        grayco_results = {}

        for prop in properties:
            grayco_results[prop] = graycoprops(glcm, prop)[0]

            grayco_features[f'{prop}'] = np.mean(grayco_results[prop])

        return grayco_features        



    def get_whole_scMTOP_embedding(self, 
                                   save_individual_wsi=False, 
                                   save_whole=True, 
                                   zernike_per_nuclei=False,
                                   num_cores=None):
        """
        Computes and saves the whole scMTOP embedding.
        This method computes the scMTOP embedding for each WSI file, concatenates
        the embeddings, and saves the whole embedding to a CSV file.
        Args:
            save_individual_wsi (bool): Whether to save the scMTOP embeddings for each WSI. Defaults to False.
            save_whole (bool): Whether to save the whole scMTOP embeddings. Defaults to True.
        Returns:
            None
        """

        emb_dfs = []

        for wsi_path, wsi_name, wsi_cellvit_path in zip(self.list_wsi_filenames, self.wsi_names, self.wsi_cellvit_paths):
            print(f"Computing the scMTOP embedding for {wsi_name}")
            print(f"Path to the WSI: {wsi_path}")
            print(f"Path to the CellViT folder: {wsi_cellvit_path}")
           # if not os.path.exists(os.path.join(self.save_path, wsi_name)):
            patch_to_cell = self.map_cells_and_patches(wsi_name, wsi_cellvit_path)

            if zernike_per_nuclei:
                if len(patch_to_cell) > 0:
                    if os.path.exists(os.path.join(self.save_path, "zernike_cells", f"{wsi_name}_zernike_cells.json")):
                        with open(os.path.join(self.save_path, "zernike_cells", f"{wsi_name}_zernike_cells.json"), 'r') as f:
                            zernike_cells = json.load(f)
                    else:
                        zernike_cells = get_zernike_moments_per_nuclei(wsi_path, wsi_cellvit_path, wsi_name, num_cores=num_cores)
                        zernike_cells_serializable = {cell: moments.tolist() for cell, moments in zernike_cells.items()}
                        
                        with open(os.path.join(self.save_path, "zernike_cells", f"{wsi_name}_zernike_cells.json"), 'w') as f:
                            json.dump(zernike_cells_serializable, f)


                    if os.path.exists(os.path.join(self.save_path, "patch_to_cell", f"{wsi_name}_patch_to_cell.json")):
                        with open(os.path.join(self.save_path, "patch_to_cell", f"{wsi_name}_patch_to_cell.json"), 'r') as f:
                            patch_to_cell = json.load(f)
                    else:
                        
                        patch_to_cell = self.map_cells_and_patches(wsi_name, wsi_cellvit_path)
                        
                        with open(os.path.join(self.save_path, "patch_to_cell", f"{wsi_name}_patch_to_cell.json"), 'w') as f:
                            json.dump(patch_to_cell, f)
                else: 
                    zernike_cells = None
            else:
                zernike_cells = None

            print(f"Length of patch to cell: {len(patch_to_cell)}", flush=True)
            emb_df = self.get_patch_features(wsi_path=wsi_path, 
                                             wsi_name=wsi_name, 
                                             wsi_cellvit_path=wsi_cellvit_path, 
                                             patch_to_cell=patch_to_cell, 
                                             zernike_cells=zernike_cells,
                                             num_cores=num_cores)
         
            emb_dfs.append(emb_df)
            print(len(emb_dfs), flush=True)
            print(emb_df.shape, flush=True)

            if len(self.wsi_names) == 1 or save_individual_wsi:
                print("Saving the scMTOP embedding at: ", os.path.join(self.save_path, f"{wsi_name}_scMTOP.csv"))
                emb_df.to_csv(os.path.join(self.save_path, f"{wsi_name}_scMTOP.csv"))
        
        self.emb_df = pd.concat(emb_dfs)
        
        if save_whole:
            print("Saving the whole scMTOP embedding at: ", os.path.join(self.save_path, f"{self.dataset_name}_scMTOP.csv"))
            self.emb_df.to_csv(os.path.join(self.save_path, f"{self.dataset_name}_scMTOP.csv"))

    def get_extracellular_morphological_features(self, patch_img, patch_mask):

        masked_image = apply_mask(image=patch_img, mask=patch_mask)

        gray_image = (rgb2gray(np.array(masked_image)[:, :, :3]) * 255).astype(np.uint8)
        # Label the connected regions in the mask
        label_image = label(patch_mask > 0)
        regions = regionprops(label_image, intensity_image=gray_image)

        features = []

        for region in regions:
            features.append({
                'Area': region.area,
                'Perimeter': region.perimeter,
                'Major Axis Length': region.major_axis_length,
                'Minor Axis Length': region.minor_axis_length,
                'Eccentricity': region.eccentricity,
                'Solidity': region.solidity
            })

        all_features = {"nb_regions": len(features),
                    "Area_mean": np.mean([ft['Area'] for ft in features]),
                    "Perimeter_mean": np.mean([ft['Perimeter'] for ft in features]),
                    "Major_Axis_Length_mean": np.mean([ft['Major Axis Length'] for ft in features]),
                    "Minor_Axis_Length_mean": np.mean([ft['Minor Axis Length'] for ft in features]),
                    "Eccentricity_mean": np.mean([ft['Eccentricity'] for ft in features]),
                    "Solidity_mean": np.mean([ft['Solidity'] for ft in features])
        }

        return all_features

    def load_emb_df_from_individual_wsi(self, savewhole=True):
        """
        Loads embedding dataframes from individual Whole Slide Images (WSI) and concatenates them into a single dataframe.
        This method reads CSV files corresponding to each WSI name stored in `self.wsi_names`. Each CSV file is expected 
        to be located in the directory specified by `self.save_path` and named in the format "{wsi_name}_scMTOP.csv". 
        The dataframes are then concatenated into a single dataframe and stored in `self.emb_df`.
        Returns:
            None
        """

        emb_dfs = []

        csv_files = glob.glob(os.path.join(self.save_path, "*scMTOP.csv"))
        for file in csv_files:
            sample_name = os.path.basename(file).split('_scMTOP.csv')[0]
            if sample_name != self.dataset_name:
                print(f"Loading the scMTOP embedding for {sample_name}")
                emb_df = pd.read_csv(file, index_col=0)
                emb_dfs.append(emb_df)
                    
        self.emb_df = pd.concat(emb_dfs)
        print(f"Shape of the whole embedding: {self.emb_df.shape}", flush=True)

        if savewhole:
            self.emb_df.to_csv(os.path.join(self.save_path, f"{self.dataset_name}_scMTOP.csv"))

    def load_emb_df_from_whole(self):
        """
        Loads embedding dataframe from the whole dataset.
        This method reads the CSV file corresponding to the whole dataset stored in `self.save_path` and named in the format
        "{self.dataset_name}_scMTOP.csv". The dataframe is then stored in `self.emb_df`.
        Returns:
            None
        """

        self.emb_df = pd.read_csv(os.path.join(self.save_path, f"{self.dataset_name}_scMTOP.csv"), index_col=0)

    def load_emb_df(self):
        try:
            if os.path.exists(os.path.join(self.save_path, f"{self.dataset_name}_scMTOP.csv")):
                print("Loading the whole scMTOP embedding from ", os.path.join(self.save_path, f"{self.dataset_name}_scMTOP.csv"))
                self.load_emb_df_from_whole()
            else:
                try:
                    self.load_emb_df_from_individual_wsi()
                except:
                    warnings.warn("No scMTOP embedding found. Embedding will be empty.")
                    self.emb_df = None
        except:
            self.emb_df = None



def compute_features_from_scMTOP(patches_info_filename=None,
                                path_to_cellvit_folder=None,
                                list_wsi_filenames=None,
                                path_to_wsis=None,
                                result_saving_folder="../results/engineered_features/scMTOP",
                                dataset_name="dataset",
                                temporary_folder="../results/tmp", 
                                save_individual_wsi=False,
                                save_whole=True, 
                                zernike_per_nuclei=False,
                                num_cores=None):
    """
    Computes scMTOP features for the given WSI files.
    Args:
        patches_info_filename (str): Filename for patches information.
        path_to_cellvit_folder (str): Path to the CellVit folder.
        list_wsi_filenames (list): List of WSI filenames.
        path_to_wsis (str): Path to whole slide images.
        result_saving_folder (str): Folder to save results.
        emb_df (pd.DataFrame): DataFrame containing embeddings.
        dataset_name (str): Name of the dataset.
        temporary_folder (str): Temporary folder for intermediate files.
        save_individual_wsi (bool): Whether to save the scMTOP embeddings for each WSI. Defaults to False.
        save_whole (bool): Whether to save the whole scMTOP embeddings. Defaults to True.
    Returns:
        None
    """
    
    scmtop = scMTOP_EngineeredFeatures(patches_info_filename=patches_info_filename,
                                       path_to_cellvit_folder=path_to_cellvit_folder,
                                       list_wsi_filenames=list_wsi_filenames,
                                       path_to_wsis=path_to_wsis,
                                       result_saving_folder=result_saving_folder,
                                       dataset_name=dataset_name,
                                       temporary_folder=temporary_folder)
    scmtop.compute_scMTOP_features()
    scmtop.get_whole_scMTOP_embedding(save_individual_wsi=save_individual_wsi, 
                                      save_whole=save_whole, 
                                      zernike_per_nuclei=zernike_per_nuclei,
                                      num_cores=num_cores)




def compute_engineered_features(method="scMTOP",
                                patches_info_filename=None,
                                patches_filename=None,
                                path_to_cellvit_folder=None,
                                list_wsi_filenames=None,
                                path_to_wsis=None,
                                result_saving_folder="../results/engineered_features/scMTOP",
                                dataset_name="dataset",
                                temporary_folder="../results/tmp",
                                save_individual_wsi=False,
                                save_whole=True, 
                                zernike_per_nuclei=False,
                                num_cores=None):
    """
    Computes engineered features for digital histopathology.
    Args:
        method (str): Method to compute engineered features.
        patches_info_filename (str): Filename that contains all patches information.
        path_to_cellvit_folder (str): Path to the CellViT folder.
        list_wsi_filenames (list): List of WSI filenames.
        path_to_wsis (str): Path to the WSIs.
        result_saving_folder (str): Output path for the results.
        dataset_name (str): Name of the dataset.
        temporary_folder (str): Temporary folder for intermediate files.
        save_individual_wsi (bool): Save the scMTOP embeddings for each WSI.
        save_whole (bool): Save the whole scMTOP embeddings.
    Returns:
        None
    """
    if method == "scMTOP":
        compute_features_from_scMTOP(patches_info_filename=patches_info_filename,
                                path_to_cellvit_folder=path_to_cellvit_folder,
                                list_wsi_filenames=list_wsi_filenames,
                                path_to_wsis=path_to_wsis,
                                result_saving_folder=result_saving_folder,
                                dataset_name=dataset_name,
                                temporary_folder=temporary_folder,
                                save_individual_wsi=save_individual_wsi,
                                save_whole=save_whole, 
                                zernike_per_nuclei=zernike_per_nuclei,
                                num_cores=num_cores)
        
    else:
        raise ValueError("Method not supported.")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compute engineered features for digital histopathology.")
    parser.add_argument("--method", type=str, default="scMTOP", help="Method to compute engineered features.")
    parser.add_argument("--patches_filename", type=str, help="Path to all patches.")
    parser.add_argument("--patches_info_filename", type=str, help="Filename that contains all patches information.")
    parser.add_argument("--path_to_cellvit_folder", type=str, help="Path to the CellViT folder.")
    parser.add_argument("--list_wsi_filenames", type=str, nargs='+', help="List of WSI filenames.")
    parser.add_argument("--path_to_wsis", type=str, help="Path to the WSIs.")
    parser.add_argument("--result_saving_folder", type=str, default="../results/engineered_features/scMTOP", help="Output path for the results.")
    parser.add_argument("--dataset_name", type=str, default="dataset", help="Name of the dataset.")
    parser.add_argument("--temporary_folder", type=str, default="../results/tmp", help="Temporary folder for intermediate files.")
    parser.add_argument("--save_individual_wsi", action="store_true", help="Save the scMTOP embeddings for each WSI.")
    parser.add_argument("--save_whole", action="store_true", help="Save the whole scMTOP embeddings.")
    parser.add_argument("--zernike_per_nuclei", action="store_true", help="Compute Zernike moments per nuclei.")
    parser.add_argument("--num_cores", type=int, default=None, help="Number of cores to use for parallel processing.")
    args = parser.parse_args()


    compute_engineered_features(
        method=args.method,
        patches_info_filename=args.patches_info_filename,
        path_to_cellvit_folder=args.path_to_cellvit_folder,
        list_wsi_filenames=args.list_wsi_filenames,
        path_to_wsis=args.path_to_wsis,
        result_saving_folder=args.result_saving_folder,
        dataset_name=args.dataset_name,
        temporary_folder=args.temporary_folder,
        save_individual_wsi=args.save_individual_wsi,
        save_whole=args.save_whole,
        zernike_per_nuclei=args.zernike_per_nuclei,
        num_cores=args.num_cores
    )