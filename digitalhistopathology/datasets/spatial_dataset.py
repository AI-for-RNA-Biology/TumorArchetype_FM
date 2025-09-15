#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import os
import glob

from digitalhistopathology.embeddings.image_embedding import ImageEmbedding


class SpatialDataset:
    PALETTE = None
    ORDER_LABEL = None

    def __init__(self, patches_folder=None):
        self.patches_folder = patches_folder
        pass

    def get_image_embeddings(self, model):
        pass

    def get_gene_embeddings(self):
        pass

    def get_engineered_features(self):
        pass

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