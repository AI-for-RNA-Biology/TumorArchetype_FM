#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import ast
import os

from PIL import Image

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

from digitalhistopathology.classification.classification import Classification
from digitalhistopathology.clustering.clustering import Clustering
from digitalhistopathology.visualization.spatial_viz import SpatialViz

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

Image.MAX_IMAGE_PIXELS = None


class Embedding(Classification, Clustering, SpatialViz):

    def __init__(
        self,
        emb=None,
        result_saving_folder=None,
        name="",
        saving_plots=False,
        label_files=None,
        palette=None,
    ):
        """Embedding class contains dimensionally reduction, clustering and visualization techniques to analyze the different types of embeddings that a spatial
        transcriptomics dataset can offer (image, genetic, engineered).

        Args:
            emb (anndata, optional): Embeddings anndata object. Defaults to None.
            result_saving_folder (str, optional): Result folder in which the results are saved. Defaults to None.
            name (str, optional): Name of the embeddings. Defaults to "".
            saving_plots (bool, optional): If the plots are saved to the result folder or not. Defaults to False.
            label_files (list, optional): List of files containing label of each spot, "x" column corresponds to the first spot coordinate, "y" column corresponds to the second. Can be csv, tsv with gzip compression or not. One per sample or the name of the file contain the sample name at the beginning with a "_" just after. Defaults to None.
        """
        Classification.__init__(self, emb, saving_plots=saving_plots, result_saving_folder=result_saving_folder)
        Clustering.__init__(self, emb, saving_plots=saving_plots, result_saving_folder=result_saving_folder)
        # SpatialViz.__init__(self, emb, saving_plots=saving_plots, result_saving_folder=result_saving_folder)

        self.name = name
        self.label_files = label_files

        self.init_result_saving_folder()
        self.palette = palette


    def init_result_saving_folder(self):
        """Create result_saving_folder if given."""
        if self.result_saving_folder is not None and not os.path.exists(self.result_saving_folder):
            try:
                os.mkdir(self.result_saving_folder)
            except Exception as e:
                print("Cannot create the results saving folder {}: {}".format(self.result_saving_folder, e))

    def compute_embeddings(self):
        """Abstract method; it is implemented in the children classes."""
        pass

    def save_embeddings(
        self,
        saving_path="../results/embeddings.h5ad",
    ):
        """Save embeddings anndata into a .h5ad format.

        Args:
            saving_path (str, optional): Path where you want to save the embeddings, it has to ends with .h5ad. Defaults to "../results/embeddings.h5ad".

        Raises:
            Exception: If the saving_path does not ends with the correct extension (.h5ad).
        """
        if not saving_path.endswith(".h5ad"):
            raise Exception("The saving path must finish with .h5ad extension")
        emb_copy = self.emb.copy()
        emb_copy.obs = emb_copy.obs.dropna(axis=1, how="all")
        emb_copy.obs = emb_copy.obs.applymap(str)
        emb_copy.write_h5ad(saving_path)


    def load_embeddings(
        self,
        data_path,
        columns_numeric=[
            "shape_pixel",
            "start_height_origin",
            "start_width_origin",
            "mean_intensity",
            "median_intensity",
            "std_intensity",
            "entropy_intensity",
            "mpp_height",
            "mpp_width",
            "shape_micron",
            "overlap_pixel",
            "radius_pixel",
            "x",
            "y",
            "predicted_label",
        ],
    ):
        """Load the embeddings into an anndata format from .h5ad file. Default values are specific to image embeddings.

        Args:
            data_path (_type_): Path the .h5ad file.
            columns_numeric (list, optional): List of the columns of the embeddings that need to be changed to the numeric format. Defaults to [ "shape_pixel", "start_height_origin", "start_width_origin", "mean_intensity", "median_intensity", "std_intensity", "entropy_intensity", "mpp_height", "mpp_width", "shape_micron", "overlap_pixel", "radius_pixel", "x", "y", ].

        Raises:
            Exception: If the data_path does not ends with the correct extension (.h5ad).
        """
        if not data_path.endswith(".h5ad"):
            raise Exception("The saving path must finish with .h5ad extension")
        self.emb = ad.read_h5ad(data_path)
        self.emb.obs = self.emb.obs.applymap(str)

        for col in columns_numeric:
            if col in self.emb.obs.columns:
                self.emb.obs[col] = pd.to_numeric(self.emb.obs[col], errors="coerce")
        self.emb.obs.replace("nan", np.nan, inplace=True)


    def add_label(self, dataset='HER2'):
        """Add the label to emb.obs from the label_files. 
        Designed specifically for the format of label files from the HER2-positive breast cancer dataset.

        Raises:
            Exception: If there is no label_files.
            Exception: If there is no spots info in emb.obs columns.
            Exception: If there is no name_origin in emb.obs columns.
        """
        if dataset == 'HER2':
            cols = self.emb.obs.columns
            if self.label_files is None:
                raise Exception("No label files")
            if "spots_info" not in cols and not ("x" in cols and "y" in cols):
                raise Exception("No spots info in emb.obs columns")
            if "name_origin" not in cols:
                raise Exception("No name_origin in emb.obs columns")

            print("Start adding labels to patches with {} files".format(len(self.label_files)))

            all_labels_df = pd.DataFrame()
            for file in self.label_files:
                compression = "gzip" if file.endswith("gz") else None
                sep = "\t" if ".tsv" in file.split("/")[-1] else ","
                current_df = pd.read_csv(file, sep=sep, compression=compression)

                if self.label_files[0].split("/")[-3].split("_")[0] == "HER2":
                    current_df["name_origin"] = file.split("/")[-1].split("_")[0]
                else:
                    current_df["name_origin"] = file.split("/")[-1].split(".")[0]

                if len(all_labels_df) == 0:
                    all_labels_df = current_df.copy()
                else:
                    print("Concatenating {} with {} rows".format(file, len(current_df)))
                    all_labels_df = pd.concat((all_labels_df, current_df), axis=0)

            all_labels_df = all_labels_df.dropna(axis=0)
            all_labels_df["x"] = all_labels_df["x"].apply(lambda x: round(x))
            all_labels_df["y"] = all_labels_df["y"].apply(lambda y: round(y))
            # one problem with one file of her2 dataset
            all_labels_df.loc[all_labels_df["label"] == "immune infiltrate⁄", "label"] = "immune infiltrate"

            if not ("x" in cols and "y" in cols):
                print("Adding spots info to emb.obs")
                spot_df = self.emb[~self.emb.obs["spots_info"].isna()].obs
                spot_df = spot_df[~(spot_df["spots_info"] == 'nan')]

                # Convert the 'spots_info' column to string
                spot_df['spots_info'] = spot_df['spots_info'].astype(str)

                # Parse the string representation into dictionaries
                spot_df['spots_info'] = spot_df['spots_info'].apply(ast.literal_eval)

                self.emb.obs = self.emb.obs.merge(
                    pd.json_normalize(spot_df["spots_info"]).set_index(spot_df.index),
                    how="left",
                    left_index=True,
                    right_index=True,
                )

            if "label" in self.emb.obs.columns:
                print("Here was the error") 
                self.emb.obs.drop(columns=["label"], inplace=True)
                print(self.emb.obs.columns)

            self.emb.obs = (
                self.emb.obs.reset_index()
                .merge(
                    all_labels_df[["x", "y", "name_origin", "label"]],
                    how="left",
                    on=["x", "y", "name_origin"],
                )
                .set_index("index")
            )

            print(
                "Added labels to {} / {} patches".format(
                    len(self.emb.obs) - self.emb.obs["label"].isna().sum(),
                    len(self.emb.obs),
                )
            )
            
        elif dataset == 'TNBC':
            
            all_labels_df = []
            for file in self.label_files:
                
                l = pd.read_csv(file, index_col=0)
                all_labels_df.append(l)
                
            label_file = pd.concat(all_labels_df, axis=0)
                
            # ids_common = set(label_file.index).intersection(set(self.emb.obs.index))
            
            # Filter the image embedding accordingly
            # self.emb = self.emb[self.emb.obs.index.isin(ids_common)]
            
            # Merge the label file with the image embedding
            # self.emb.obs = self.emb.obs.merge(label_file, left_index=True, right_index=True)
            self.emb.obs = self.emb.obs.merge(label_file, left_index=True, right_index=True, how="left")

            
            self.emb.obs.replace("undetermined", np.nan, inplace=True)
            
            print(f"Added TNBC labels to {len(self.emb.obs)} patches. Patches with no label were set to NaNs.")
        else:
            raise Exception("Dataset not supported for label addition.")   


