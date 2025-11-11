#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import glob
import math
from multiprocessing.resource_tracker import getfd
import os
import pandas as pd
import numpy as np
import seaborn as sns
from PIL import Image
from json import load
from digitalhistopathology.datasets.spaceranger_utils import aggregate_spots_to_resolution, load_spaceranger_coding_genes
from natsort import natsorted
import openslide

Image.MAX_IMAGE_PIXELS = None

from digitalhistopathology.embeddings.embedding import Embedding
from digitalhistopathology.embeddings.gene_embedding import GeneEmbedding
from digitalhistopathology.embeddings.image_embedding import ImageEmbedding
from digitalhistopathology.engineered_features.engineered_features import EngineeredFeatures
from digitalhistopathology.datasets.spatial_dataset import SpatialDataset


# R imports for TNBC  
import rpy2.robjects as ro
import rpy2.robjects.vectors as rvec
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()
EBImage = importr('EBImage')



class HER2Dataset(SpatialDataset):
    """HER2-positive breast cancer dataset from https://www.nature.com/articles/s41467-021-26271-2."""

    PALETTE = {
        "invasive cancer": "red",
        "cancer in situ": "orange",
        "immune infiltrate": "yellow",
        "breast glands": "green",
        "connective tissue": "blue",
        "adipose tissue": "cyan",
        "undetermined": "lightgrey",
    }
    ORDER_LABEL = [
        "invasive cancer",
        "cancer in situ",
        "immune infiltrate",
        "breast glands",
        "connective tissue",
        "adipose tissue",
        "undetermined",
    ]

    def __init__(self, patches_folder=None, saving_emb_folder=None):
        """
        Args:
            patches_folder (str, optional): Folder that contains patches. Defaults to None.
            saving_emb_folder (str, optional): Path to folder where to save the embeddings. Defaults to None.
        """

        super().__init__(patches_folder)
        self.saving_emb_folder = saving_emb_folder
        self.name = "HER2"
        super().__init__(patches_folder)
        # Image embeddings
        if patches_folder is not None:
            self.init_patches_filenames()
        self.label_filenames = sorted(
            glob.glob("data/HER2/meta/*.tsv")
        )

        # Gene embeddings
        AREA_SPOT_HER2_PIXEL2 = 8000  # from QPath
        self.spot_diameter = 2 * int(np.sqrt(AREA_SPOT_HER2_PIXEL2 / math.pi))
        self.genes_count_filenames = sorted(
            glob.glob("data/HER2/count-matrices/*")
        )
        self.genes_spots_filenames = sorted(
            glob.glob("data/HER2/spot-selections/*")
        )
        self.images_filenames = sorted(
            glob.glob("data/HER2/images/HE/*")
        )
        self.samples_names = [
            f.split("/")[-1].split(".")[0][0:2] for f in self.images_filenames
        ]

    def get_image_embeddings(self, model, filename="ie", emb_path=None):
        """Compute the image embedding or load it if it already exists.

        Args:
            model (models.PretrainedModel): Pretrained deep learning vision encoder.
            filename (str, optional): Filename of the .h5ad image embedding. Defaults to "ie".

        Returns:
            image_embedding.ImageEmbedding: image embedding
        """
        self.init_patches_filenames()
        ie = ImageEmbedding(
            patches_filenames=self.patches_filenames,
            patches_info_filename=self.patches_info_filename,
            label_files=self.label_filenames,
            pretrained_model=model,
            name=self.name + "_" + model.name,
        )
        try:
            if os.path.exists(self.saving_emb_folder):
                loading_file = os.path.join(
                    self.saving_emb_folder, "{}.h5ad".format(filename)
                )
            elif emb_path is not None:
                loading_file = os.path.join(emb_path, "{}.h5ad".format(filename))
            else:
                loading_file = (
                    "results/HER2/embeddings/images_embeddings/save/{}/{}.h5ad".format(
                        model.name.lower(), self.saving_emb_folder
                    )
                )
            ie.load_embeddings(loading_file)
            if "label" not in ie.emb.obs.columns:
                ie.add_label()
            print("Fill emb with: {}".format(loading_file))
            print(ie.emb)
        except Exception as e:
            print("Cannot load images embeddings: {}".format(e))
        return ie

    def get_gene_embeddings(
        self,
        compute_emb=False,
        preprocessing=False,
        spot_normalization=False,
        load=False,
        filename="ge",
    ):
        """Compute the gene embedding or load it if it already exists.

        Args:
            compute_emb (bool, optional): If we need compute the gene embedding. Defaults to False.
            preprocessing (bool, optional): If preprocessing is done. Defaults to False.
            spot_normalization (bool, optional): If spots normalization is done if preprocessing. Defaults to False.
            load (bool, optional): If we need to load gene embedding from file. Defaults to False.
            filename (str, optional): Filename of the .h5ad gene embedding. Defaults to "ge".

        Returns:
            gene_embedding.GeneEmebdding: gene embedding
        """
        # Filter only spot patches
        patches_filenames = [
            f for f in self.patches_filenames if "spot" in f.split("/")[-1]
        ]
        ge = GeneEmbedding(
            spot_diameter_fullres=self.spot_diameter,
            samples_names=self.samples_names,
            genes_count_filenames=self.genes_count_filenames,
            spots_filenames=self.genes_spots_filenames,
            image_filenames=self.images_filenames,
            label_files=self.label_filenames,
            patches_filenames=patches_filenames,
            name=self.name,
            st_method="old_st",
        )
        if compute_emb:
            ge.compute_embeddings()
            if preprocessing:
                ge.emb = ge.preprocessing(spot_norm=spot_normalization)
            print(ge.emb)
        elif load:
            try:
                loading_file = os.path.join(
                    self.saving_emb_folder, "{}.h5ad".format(filename)
                )
                ge.load_embeddings(
                    loading_file,
                    columns_ast=[],
                    columns_numeric=[
                        "imagecol",
                        "imagerow",
                        "x",
                        "y",
                        "total_genes_count",
                        "predicted_label",
                    ],
                )
                if "label" not in ge.emb.obs.columns:
                    ge.add_label()
                print("Fill emb with: {}".format(loading_file))
                print(ge.emb)
            except Exception as e:
                print("Cannot load genes embeddings: {}".format(e))
        return ge

    def get_engineered_features(self, remove_nan=True, filename="ef", emb_path=None):
        """Compute the engineered features or load it if it already exists.

        Args:
            remove_nan (bool, optional): If rows with NaN. Defaults to True.
            filename (str, optional): Filename of the .h5ad engineered features. Defaults to "ef".

        Returns:
            engineered_features.EngineeredFeatures: engineered features
        """
        ef = EngineeredFeatures(
            patches_info_filename=self.patches_info_filename,
            name=self.name,
            label_files=self.label_filenames,
        )

        if (self.saving_emb_folder is not None) and (
            os.path.exists(self.saving_emb_folder)
        ):
            loading_file = os.path.join(
                self.saving_emb_folder, "{}.h5ad".format(filename)
            )
        elif emb_path is not None:
            loading_file = os.path.join(emb_path, "{}.h5ad".format(filename))
        else:
            loading_file = (
                "results/HER2/embeddings/engineered_features/save/{}.h5ad".format(
                    self.saving_emb_folder
                )
            )

        try:
            if (self.saving_emb_folder is not None) and (
                os.path.exists(self.saving_emb_folder)
            ):
                loading_file = os.path.join(
                    self.saving_emb_folder, "{}.h5ad".format(filename)
                )
            elif emb_path is not None:
                loading_file = os.path.join(emb_path, "{}.h5ad".format(filename))
            else:
                loading_file = (
                    "results/HER2/embeddings/engineered_features/save/{}.h5ad".format(
                        self.saving_emb_folder
                    )
                )

            print(f"loading_file: {loading_file}", flush=True)
            ef.load_embeddings(loading_file)

            if remove_nan:
                # remove rows with nan
                print("Remove Nan...")
                print(ef.emb.shape)
                ef.emb = ef.emb[~np.isnan(ef.emb.X).any(axis=1), :]
                print(ef.emb.shape)
            if "label" not in ef.emb.obs.columns:
                ef.add_label()
            print("Fill emb with: {}".format(loading_file))
            print(ef.emb)
        except Exception as e:
            print("Cannot load engineered features: {}".format(e))
        return ef

    def get_palette_2():
        palette = sns.color_palette(palette="bright")
        palette_2 = dict()
        palette_2[0] = palette[0]
        palette_2[1] = palette[5]
        palette_2[2] = palette[1]
        palette_2[3] = palette[2]
        palette_2[4] = palette[6]
        palette_2[5] = palette[3]
        palette_2[6] = palette[4]
        palette_2[7] = palette[7]
        palette_2[8] = palette[8]
        palette_2[9] = palette[9]
        palette_2[10] = sns.color_palette(palette="colorblind")[5]
        return palette_2
    
    def add_label_from_emb(self, emb):
        
        cols = emb.emb.obs.columns
        if emb.emb.label_files is None:
            raise Exception("No label files")
        if "spots_info" not in cols and not ("x" in cols and "y" in cols):
            raise Exception("No spots info in emb.obs columns")
        if "name_origin" not in cols:
            raise Exception("No name_origin in emb.obs columns")

        print("Start adding labels to patches with {} files".format(len(emb.label_files)))

        all_labels_df = pd.DataFrame()
        for file in emb.label_files:
            compression = "gzip" if file.endswith("gz") else None
            sep = "\t" if ".tsv" in file.split("/")[-1] else ","
            current_df = pd.read_csv(file, sep=sep, compression=compression)

            if emb.label_files[0].split("/")[-3].split("_")[0] == "HER2":
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
            spot_df = emb.emb[~emb.emb.obs["spots_info"].isna()].obs
            spot_df = spot_df[~(spot_df["spots_info"] == 'nan')]

            # Convert the 'spots_info' column to string
            spot_df['spots_info'] = spot_df['spots_info'].astype(str)

            # Parse the string representation into dictionaries
            spot_df['spots_info'] = spot_df['spots_info'].apply(ast.literal_eval)

            emb.emb.obs = emb.emb.obs.merge(
                pd.json_normalize(spot_df["spots_info"]).set_index(spot_df.index),
                how="left",
                left_index=True,
                right_index=True,
            )

        if "label" in emb.emb.obs.columns:
            print("Here was the error") 
            emb.emb.obs.drop(columns=["label"], inplace=True)
            print(emb.emb.obs.columns)

        emb.emb.obs = (
            emb.emb.obs.reset_index()
            .merge(
                all_labels_df[["x", "y", "name_origin", "label"]],
                how="left",
                on=["x", "y", "name_origin"],
            )
            .set_index("index")
        )

        print(
            "Added labels to {} / {} patches".format(
                len(emb.emb.obs) - emb.emb.obs["label"].isna().sum(),
                len(emb.emb.obs),
            )
        )



class TNBCDataset(SpatialDataset):
    # PALETTE = {
    #     "invasive cancer": "#017801",
    #     "Necrosis": "#000000",
    #     "Fat tissue": "#000080",
    #     "Vessels": "#dc0000",
    #     "Lactiferous duct": "#9980e6",
    #     "in situ": "#ccffcc",
    #     "Lymphoid nodule": "#80801a",
    #     "Lymphocyte": "#c4417f",
    #     "Stroma": "#ff9980",
    #     "Nerve": "#4d8080",
    #     "Heterologous elements": "#808080",
    # }

    # ORDER_LABEL = [
    #     "invasive cancer",
    #     "Necrosis",
    #     "Fat tissue",
    #     "Vessels",
    #     "Lactiferous duct",
    #     "in situ",
    #     "Lymphoid nodule",
    #     "Lymphocyte",
    #     "Stroma",
    #     "Nerve",
    #     "Heterologous elements",
    # ]
    
    PALETTE = {
        "invasive cancer": "red",
        "cancer in situ": "orange",
        "immune infiltrate": "yellow",
        "breast glands": "green",
        "connective tissue": "blue",
        "adipose tissue": "cyan",
        "undetermined": "lightgrey",
    }
    ORDER_LABEL = [
        "invasive cancer",
        "cancer in situ",
        "immune infiltrate",
        "breast glands",
        "connective tissue",
        "adipose tissue",
        "undetermined",
    ]

    def __init__(self, 
                 patches_folder=None, 
                 saving_emb_folder=None, 
                 dataDir=None,
                 gene_counts_path=None):
        """
        Args:
            patches_folder (str, optional): Folder that contains patches. Defaults to None.
            saving_emb_folder (str, optional): Path to folder where to save the embeddings. Defaults to None.
            dataDir (str, optional): Path to the folder that contains the images. Defaults to None.
        """
        self.dataDir = dataDir
        super().__init__(patches_folder)
        self.saving_emb_folder = saving_emb_folder
        self.name = "TNBC"
        # Image embeddings
        if patches_folder is not None:
            self.init_patches_filenames()
        self.label_filenames = sorted(glob.glob("results/TNBC/compute_patches/all/*labels.csv"))

        self.spot_diameter = 253
        
        if gene_counts_path is not None:
            self.genes_count_filenames = sorted(
                glob.glob(os.path.join(gene_counts_path, "*.csv"))
            )
        else:
            self.genes_count_filenames = sorted(
                glob.glob("results/TNBC/count-matrices/all/*.csv")
            )
        self.genes_spots_filenames = self.label_filenames
        self.images_filenames = [
            os.path.join(
                self.dataDir, "Images", "imagesHD", os.path.basename(name)[:-4] + ".jpg"
            )
            for name in self.genes_count_filenames
        ]
        self.samples_names = [
            f.split("/")[-1].split(".")[0] for f in self.images_filenames
        ]


    def preprocess(self, annotated_only=True):

        dataDir = self.dataDir
    
        # Set the image names
        if annotated_only:
            images_name_list = natsorted([os.path.basename(f).replace('.png', '.jpg') for f in glob.glob(f"{dataDir}/Images/imageAnnotations/*.png")])
            images_filenames_spot_patches = [f"{dataDir}/Images/imagesHD/{image}" for image in images_name_list]
        else:
            images_filenames_spot_patches = natsorted(glob.glob(f"{dataDir}/Images/imagesHD/*.jpg"))
            images_name_list = [os.path.basename(f) for f in images_filenames_spot_patches]

        def convert_to_df(data):
            # Check if data is an R Matrix
            if isinstance(data, rvec.Matrix):
                return pandas2ri.rpy2py(ro.r['as.data.frame'](data))
            elif isinstance(data, rvec.Vector):
                return pandas2ri.rpy2py(data)
            # Check if data is an R DataFrame
            elif isinstance(data, rvec.DataFrame):
                return pandas2ri.rpy2py(data)

            # Check if data is an EBImage object (matrix or image)
            elif isinstance(data, ro.rinterface.Sexp):
                # Check if it's an EBImage object (usually matrix)
                try:
                    ro.r('library(EBImage)')
                    image_array = ro.r['imageData'](data)  # Convert EBImage object to matrix
                    if isinstance(image_array, rvec.Matrix):
                        # Convert the image matrix to a DataFrame (pixel values)
                        return pd.DataFrame(pandas2ri.rpy2py(image_array))
                except Exception as e:
                    print(f"Error converting EBImage object: {e}")
                    return None

            # Return None for unsupported data types
            else:
                return None

        # Set the spots paths
        spots_path = [f"{dataDir}/byArray/{image.split('_')[1]}/{image.split('_')[2][:2]}/allSpots.RDS" for image in images_name_list]
        
        # Initialize an empty list to store DataFrames
        spots_dfs = []

        # Process each RDS file
        idxs_to_remove = []
        for i, spots in enumerate(spots_path):

            if not os.path.exists(spots):
                print(f"File not found: {spots}")
                idxs_to_remove.append(i)
                continue

            try:
                # Read the RDS file
                spots_df = ro.r['readRDS'](spots)
                spots_df = convert_to_df(spots_df)
                # Add a column to track the source file
                spots_df['pixel_x'] = spots_df['pixel_x'] * 10.0 / 3.0 # This is to get the coordinates in the HD image (31*744*31744 px). Reminder: the large image is 9523*9523px.
                spots_df['pixel_y'] = spots_df['pixel_y'] * 10.0 / 3.0 # This is to get the coordinates in the HD image (31*744*31744 px). Reminder: the large image is 9523*9523px.
                # Append to the list of DataFrames
                spots_dfs.append(spots_df)
            
            except Exception as e:
                print(f"Error processing {spots}: {e}")
        
        # Remove the missing RDS files 
        for idx in sorted(idxs_to_remove, reverse=True):
            del spots_path[idx]
            del images_filenames_spot_patches[idx]
            del images_name_list[idx]

        # Compute spot diameter
        inter_distances = []

        for i in range(2, 65):
            inter_distances.append(spots_df.groupby("x").mean().loc[i+1]['pixel_x'] - spots_df.groupby("x").mean().loc[i]['pixel_x'])
            
        for i in range(2, 62):
            inter_distances.append(spots_df.groupby("y").mean().loc[i+1]['pixel_y'] - spots_df.groupby("y").mean().loc[i]['pixel_y'])
            
        center_to_center_distance = np.mean(inter_distances)

        # Center to center distance: 150um (according to the paper)
        resolution = 150/center_to_center_distance
        print(f"Resolution: {resolution:.2f}um/px")

        # Spot diameter: 100um (according to the paper)
        spot_diameter = round(100/resolution)
        print(f"Spot diameter: {spot_diameter} px")
        
        # Generate all RDS file paths
        all_dataframes = []
        annotsBySpot = natsorted(glob.glob(f"{dataDir}/Robjects/annotsBySpot/TNBC*.RDS"))

        if not annotated_only:
            annotated_names = natsorted([os.path.basename(f).replace('.png', '') for f in glob.glob(f"{dataDir}/Images/imageAnnotations/TNBC*.png")])

            # Process each RDS file
            for idx, img_name in enumerate(images_name_list):
                spots_df = spots_dfs[idx]
                
                if img_name.replace('.jpg', '') not in annotated_names:
                    spots_df['source_file'] = img_name.replace('.jpg', '')
                    all_dataframes.append(spots_df)
                else:
                    annots = f"{dataDir}/Robjects/annotsBySpot/{img_name.replace('.jpg', '').split('_')[0]}.RDS"

                    try:
                        # Read the RDS file
                        loaded_data = ro.r['readRDS'](annots)

                        # Extract 'annots'
                        names = loaded_data.names
                        for i, name in enumerate(names):
                            if name == 'annots':
                                annots_rds = loaded_data[i]
                        annots_df = convert_to_df(annots_rds)
                        # annots_df[
                        #     'source_file'] = os.path.basename(annots).replace('.RDS', '')
                        
                        merged_dfs = pd.merge(annots_df, spots_df, left_index=True, right_index=True, suffixes=('_spot', '_annot'), how='outer')
                        merged_dfs['source_file'] = img_name.replace('.jpg', '')

                        # Add a column to track the source file
                        
                        # Append to the list of DataFrames
                        all_dataframes.append(merged_dfs)

                    except Exception as e:
                        print(f"Error processing {spots}: {e}")# Combine all DataFrames into one giant DataFrame
            
        else:
            
            # Process each RDS file
            for spots_df, annots in zip(spots_dfs, annotsBySpot):

                try:
                    # Read the RDS file
                    loaded_data = ro.r['readRDS'](annots)

                    # Extract 'annots'
                    names = loaded_data.names
                    for i, name in enumerate(names):
                        if name == 'annots':
                            annots_rds = loaded_data[i]
                    annots_df = convert_to_df(annots_rds)

                    spots_df = spots_df[spots_df['selected'] == 1]  # Filter spots to only include selected ones
                    merged_dfs = pd.merge(annots_df, spots_df, left_index=True, right_index=True, suffixes=('_spot', '_annot'), how='outer')
                    merged_dfs['source_file'] = os.path.basename(annots).replace('.RDS', '')
                    
                    # Append to the list of DataFrames
                    all_dataframes.append(merged_dfs)

                except Exception as e:
                    print(f"Error processing {spots}: {e}")
        
        # Save the results
        self.spot_diameter = spot_diameter
        self.images_filenames = images_filenames_spot_patches
        self.all_dataframes = all_dataframes
        

    def get_image_embeddings(self, model, filename="ie", emb_path=None):
        """Compute the image embedding or load it if it already exists.

        Args:
            model (models.PretrainedModel): Pretrained deep learning vision encoder.
            filename (str, optional): Filename of the .h5ad image embedding. Defaults to "ie".

        Returns:
            image_embedding.ImageEmbedding: image embedding
        """
        self.init_patches_filenames()
        ie = ImageEmbedding(
            patches_filenames=self.patches_filenames,
            patches_info_filename=self.patches_info_filename,
            label_files=self.label_filenames,
            pretrained_model=model,
            name=self.name + "_" + model.name,
        )
        try:
            if os.path.exists(self.saving_emb_folder):
                loading_file = os.path.join(
                    self.saving_emb_folder, "{}.h5ad".format(filename)
                )
            elif emb_path is not None:
                loading_file = os.path.join(emb_path, "{}.h5ad".format(filename))
            else:
                loading_file = (
                    "results/TNBC/embeddings/images_embeddings/save/{}/{}.h5ad".format(
                        model.name.lower(), self.saving_emb_folder
                    )
                )
            ie.load_embeddings(loading_file)
            if "label" not in ie.emb.obs.columns:
                ie.add_label()
            print("Fill emb with: {}".format(loading_file))
            print(ie.emb)
        except Exception as e:
            print("Cannot load images embeddings: {}".format(e))
        return ie

    def get_gene_embeddings(
        self,
        compute_emb=False,
        preprocessing=False,
        spot_normalization=False,
        load=False,
        filename="ge",
        molecular_emb_path=None,
    ):
        """Compute the gene embedding or load it if it already exists.

        Args:
            compute_emb (bool, optional): If we need compute the gene embedding. Defaults to False.
            preprocessing (bool, optional): If preprocessing is done. Defaults to False.
            spot_normalization (bool, optional): If spots normalization is done if preprocessing. Defaults to False.
            load (bool, optional): If we need to load gene embedding from file. Defaults to False.
            filename (str, optional): Filename of the .h5ad gene embedding. Defaults to "ge".
            molecular_emb_path (str, optional): Path to pre-processed parquet file. Defaults to None.

        Returns:
            gene_embedding.GeneEmebdding: gene embedding
        """
        ge = GeneEmbedding(
            spot_diameter_fullres=self.spot_diameter,
            samples_names=self.samples_names,
            genes_count_filenames=self.genes_count_filenames,
            spots_filenames=self.genes_spots_filenames,
            image_filenames=self.images_filenames,
            label_files=self.label_filenames,
            patches_filenames=getattr(self, 'patches_filenames', None),
            name=self.name,
            st_method="old_st",
        )
        
        if molecular_emb_path is not None:
            # Load from pre-processed parquet file (raw, normalized, or combat-corrected)
            print(f"Loading gene embeddings from parquet: {molecular_emb_path}")
            self._load_from_parquet(ge, molecular_emb_path)
            
        elif compute_emb:
            ge.compute_embeddings()
            if preprocessing:
                ge.emb = ge.preprocessing(spot_norm=spot_normalization)
            print(ge.emb)
            
        elif load:
            try:
                loading_file = os.path.join(
                    self.saving_emb_folder, "{}.h5ad".format(filename)
                )
                ge.load_embeddings(
                    loading_file,
                    columns_ast=[],
                    columns_numeric=[
                        "imagecol",
                        "imagerow",
                        "x",
                        "y",
                        "total_genes_count",
                        "predicted_label",
                    ],
                )
                if "label" not in ge.emb.obs.columns:
                    ge.add_label()
                print("Fill emb with: {}".format(loading_file))
                print(ge.emb)
            except Exception as e:
                print("Cannot load genes embeddings: {}".format(e))
        return ge

    def get_engineered_features(self, remove_nan=True, filename="ef", emb_path=None):
        """Compute the engineered features or load it if it already exists.

        Args:
            remove_nan (bool, optional): If rows with NaN. Defaults to True.
            filename (str, optional): Filename of the .h5ad engineered features. Defaults to "ef".

        Returns:
            engineered_features.EngineeredFeatures: engineered features
        """
        ef = EngineeredFeatures(
            patches_info_filename=self.patches_info_filename,
            name=self.name,
            label_files=self.label_filenames,
        )

        if (self.saving_emb_folder is not None) and (
            os.path.exists(self.saving_emb_folder)
        ):
            loading_file = os.path.join(
                self.saving_emb_folder, "{}.h5ad".format(filename)
            )
        elif emb_path is not None:
            loading_file = os.path.join(emb_path, "{}.h5ad".format(filename))
        else:
            loading_file = (
                "results/TNBC/embeddings/engineered_features/save/{}.h5ad".format(
                    self.saving_emb_folder
                )
            )

        try:
            if (self.saving_emb_folder is not None) and (
                os.path.exists(self.saving_emb_folder)
            ):
                loading_file = os.path.join(
                    self.saving_emb_folder, "{}.h5ad".format(filename)
                )
            elif emb_path is not None:
                loading_file = os.path.join(emb_path, "{}.h5ad".format(filename))
            else:
                loading_file = (
                    "results/TNBC/embeddings/engineered_features/save/{}.h5ad".format(
                        self.saving_emb_folder
                    )
                )

            print(f"loading_file: {loading_file}", flush=True)
            ef.load_embeddings(loading_file)

            if remove_nan:
                # remove rows with nan
                print("Remove Nan...")
                print(ef.emb.shape)
                ef.emb = ef.emb[~np.isnan(ef.emb.X).any(axis=1), :]
                print(ef.emb.shape)
            if "label" not in ef.emb.obs.columns:
                ef.add_label()
            print("Fill emb with: {}".format(loading_file))
            print(ef.emb)
        except Exception as e:
            print("Cannot load engineered features: {}".format(e))
        return ef

    def _load_from_parquet(self, ge, molecular_emb_path):
        """Load gene embeddings from processed parquet file"""
        import anndata
        
        try:
            # Load parquet file
            molecular_emb_df = pd.read_parquet(molecular_emb_path)
            print(f"Loaded parquet with shape: {molecular_emb_df.shape}")
            
            # Handle different parquet formats
            if 'spot_id' in molecular_emb_df.columns:
                # spot_id is a column - set as index
                molecular_emb_df.set_index('spot_id', inplace=True)
                print("Using 'spot_id' column as index")
            elif molecular_emb_df.columns[0] not in molecular_emb_df.index.names:
                # First column might be spot names - set as index
                molecular_emb_df.set_index(molecular_emb_df.columns[0], inplace=True)
                print(f"Using first column as index")

            # Change index to match same convention than before
            if "parquet" in molecular_emb_df.index[0]:
                molecular_emb_df.index = [idx.split('.parquet.')[1].replace('_X', '_spot') for idx in molecular_emb_df.index]
                print("Modified index to match our convention")
            
            # Create AnnData object (spots x genes)
            molecular_emb = anndata.AnnData(
                X=molecular_emb_df.values, 
                obs=pd.DataFrame(index=molecular_emb_df.index), 
                var=pd.DataFrame(index=molecular_emb_df.columns)
            )
            
            # Add sample information to obs
            # molecular_emb.obs['sample'] = [idx.split('_')[0] for idx in molecular_emb.obs.index]
            molecular_emb.obs['patient'] = [idx.split('_')[0] for idx in molecular_emb.obs.index]
            molecular_emb.obs['tumor'] = [idx.split('_')[0] for idx in molecular_emb.obs.index]
            molecular_emb.obs['name_origin'] = [idx.split('_spot')[0] for idx in molecular_emb.obs.index]

            # Set the embedding in GeneEmbedding object
            ge.emb = molecular_emb
            
            print(f"Created AnnData object with shape: {ge.emb.shape}")
            print(f"Sample spots: {list(ge.emb.obs.index[:3])}")
            print(f"Sample genes: {list(ge.emb.var.index[:5])}")
            
            # Add labels if label files are available
            if hasattr(self, 'label_filenames') and self.label_filenames:
                try:
                    ge.add_label(dataset=self.name)
                    print("Added tissue labels to gene embedding")
                except Exception as e:
                    print(f"Could not add labels: {e}")
                    
        except Exception as e:
            print(f"Error loading from parquet: {e}")
            raise

    def get_palette_2():
        palette = sns.color_palette(palette="bright")
        palette_2 = dict()
        palette_2[0] = palette[0]
        palette_2[1] = palette[5]
        palette_2[2] = palette[1]
        palette_2[3] = palette[2]
        palette_2[4] = palette[6]
        palette_2[5] = palette[3]
        palette_2[6] = palette[4]
        palette_2[7] = palette[7]
        palette_2[8] = palette[8]
        palette_2[9] = palette[9]
        palette_2[10] = sns.color_palette(palette="colorblind")[5]
        return palette_2
    

class VisiumHDdataset(SpatialDataset):
    def __init__(self, 
                 patches_folder=None, 
                 saving_emb_folder=None, 
                 raw_images_dir=None,
                 spaceranger_dir=None, 
                 spot_diameter_micron=100,
                 name="VisiumHD"):
        
        # For ovarian
        # name = "Ovarian"
        # raw_images_dir = "/storage/research/dbmr_luisierlab/database/Ovarian_Visium_GTOP/Visium_HD"
        # spaceranger_dir = "/storage/research/dbmr_luisierlab/temp/lfournier/Data/Ovarian_Visium_GTOP/hg38/spaceranger"
        
        """
        Args:
            patches_folder (str, optional): Folder that contains patches. Defaults to None.
            saving_emb_folder (str, optional): Path to folder where to save the embeddings. Defaults to None.
        """
        super().__init__(patches_folder)
        self.saving_emb_folder = saving_emb_folder
        self.name = name
        self.spaceranger_dir = spaceranger_dir
        self.spot_diameter_micron = spot_diameter_micron
        self.raw_images_dir = raw_images_dir
        
        # Image embeddings
        if patches_folder is not None:
            self.init_patches_filenames()

        self.images_filenames = glob.glob(os.path.join(self.raw_images_dir, "**", "*.qptiff"), recursive=True)

        self.spot_diameter = self.get_spot_diameter()
        print(f"Spot diameter (in pixels): {self.spot_diameter}")
        
        
    def get_spot_diameter(self):
        slide = openslide.OpenSlide(self.images_filenames[0])
        x_resolution = float(slide.properties.get('tiff.XResolution', 0))
        microns_per_pixel_x = 10000 / x_resolution if x_resolution > 0 else None
        spot_diameter = round(self.spot_diameter_micron / microns_per_pixel_x) if microns_per_pixel_x else None
        return spot_diameter
    
    def preprocess(self, gtf_path=None):
        sample_dirs = [d for d in os.listdir(self.spaceranger_dir) if os.path.isdir(os.path.join(self.spaceranger_dir, d))] 
        
        # Ensure to have the same order than self.images_filenames
        filenames = [os.path.basename(f).replace('.qptiff', '').replace('-', '_') for f in self.images_filenames]
        all_dataframes = []
        
        
        for filename in filenames:
            sample_dirs_simplified = ['_'.join(d.split('_')[1:]) for d in sample_dirs]
            # find the sample dir that matches the filename
            matched_sample_simplified_dir = [d for d in sample_dirs_simplified if d in filename]
            if len(matched_sample_simplified_dir) == 0:
                raise Exception(f"No matching sample dir found for image filename: {filename}") 
            elif len(matched_sample_simplified_dir) > 1:
                raise Exception(f"Multiple matching sample dirs found for image filename: {filename}: {matched_sample_simplified_dir}")
            else:
                sample_dir = sample_dirs[sample_dirs_simplified.index(matched_sample_simplified_dir[0])]
            
                adata_002um_coding = load_spaceranger_coding_genes(visium_output_dir=os.path.join(self.spaceranger_dir, sample_dir, "outs"),
                                                                gtf_path=gtf_path, 
                                                                resolution="002um", 
                                                                filtered=False)
                
                adata_agg = aggregate_spots_to_resolution(adata_002um_coding, 
                                                            current_resolution_um=2, 
                                                            target_resolution_um=self.spot_diameter_micron)

                adata_agg.obs["source_file"] = sample_dir
                adata_agg.obs["pixel_x"] = adata_agg.obs["pxl_row_in_fullres"]
                adata_agg.obs["pixel_y"] = adata_agg.obs["pxl_col_in_fullres"]
                adata_agg.obs["x"] = adata_agg.obs["large_grid_row"]
                adata_agg.obs["y"] = adata_agg.obs["large_grid_col"]
            
            all_dataframes.append(adata_agg.obs)
        
        self.all_dataframes = all_dataframes
        
    

