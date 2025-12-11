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

from digitalhistopathology.embeddings.image_embedding import ImageEmbedding
from digitalhistopathology.engineered_features.engineered_features import EngineeredFeatures


class SpatialDataset:
    PALETTE = None
    ORDER_LABEL = None

    def __init__(self, 
                 patches_folder=None, 
                 saving_emb_folder=None, 
                 name=None, 
                 patches_filenames=None, 
                 patches_info_filename=None, 
                 label_filenames=None):
        
        self.patches_folder = patches_folder
        self.saving_emb_folder = saving_emb_folder
        self.name = name
        self.patches_filenames = patches_filenames
        self.patches_info_filename = patches_info_filename
        self.label_filenames = label_filenames

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
                    f"results/{self.name}/embeddings/images_embeddings/save/{model.name.lower()}/{self.saving_emb_folder}.h5ad"
                )
            ie.load_embeddings(loading_file)
            if "label" not in ie.emb.obs.columns:
                ie.add_label(dataset=self.name)
            print("Fill emb with: {}".format(loading_file))
            print(ie.emb)
        except Exception as e:
            print("Cannot load images embeddings: {}".format(e))
        return ie

    def get_gene_embeddings(self):
        pass

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
                f"results/{self.name}/embeddings/engineered_features/save/{self.saving_emb_folder}.h5ad"
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
                ef.add_label(dataset=self.name)
            print("Fill emb with: {}".format(loading_file))
            print(ef.emb)
        except Exception as e:
            print("Cannot load engineered features: {}".format(e))
        return ef

    def init_patches_filenames(self):

        if glob.glob(self.patches_folder + "/*.tiff"):
            self.patches_filenames = sorted(glob.glob(self.patches_folder + "/*.tiff"))
        else:
            self.patches_filenames = glob.glob(
                os.path.join(self.patches_folder, "*.hdf5")
            )
        self.patches_info_filename = os.path.join(
            self.patches_folder, "patches_info.pkl.gz"
        )


class MixedImageDataset(SpatialDataset):
    """Dataset with patches from different dataset."""

    def __init__(self, folder, saving_emb_folder, label_filenames=None, name="Mixed image dataset"):
        """
        Args:
            folder (str): Path to folder with the patches.
            saving_emb_folder (str): Path to folder where to save the image embedding.
            label_filenames (list, optional): List of filenames with labels. Defaults to None.
            name (str, optional): Name of the dataset. Defaults to "Mixed image dataset".
        """
        self.folder = folder
        self.saving_emb_folder = saving_emb_folder
        self.patches_filenames = sorted(glob.glob("{}/*.tiff".format(folder)))
        self.patches_info_filename = os.path.join(folder, "patches_info.pkl.gz")
        self.name = name
        self.label_filenames = label_filenames

    def get_image_embeddings(self, model, filename="ie"):
        """Compute the image embedding or load it if it already exists.

        Args:
            model (models.PretrainedModel): Pretrained deep learning vision encoder.
            filename (str, optional): Filename of the .h5ad image embedding. Defaults to "ie".

        Returns:
            image_embedding.ImageEmbedding: image embedding
        """
        ie = ImageEmbedding(
            patches_filenames=self.patches_filenames,
            patches_info_filename=self.patches_info_filename,
            pretrained_model=model,
            label_files=self.label_filenames,
            name=self.name + "_" + model.name,
        )
        try:
            if os.path.exists(self.folder):
                loading_file = os.path.join(self.saving_emb_folder, "{}.h5ad".format(filename))
            else:
                loading_file = "../results/embeddings/images_embeddings/save/{}/{}.h5ad".format(
                    model.name.lower(), self.saving_emb_folder
                )
            ie.load_embeddings(loading_file)
        except Exception as e:
            print("Cannot load images embeddings: {}".format(e))
        return ie