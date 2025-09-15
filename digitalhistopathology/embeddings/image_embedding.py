#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import gzip
import pickle

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import torch

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

from digitalhistopathology.embeddings.embedding import Embedding
from digitalhistopathology.datasets.image_embedding_dataset import ImageEmbeddingDataset
from digitalhistopathology.visualization.spatial_viz import SpatialViz
from tqdm import tqdm


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")



class ImageEmbedding(Embedding, SpatialViz):

    def __init__(
        self,
        patches_filenames=None,
        patches_info_filename=None,
        emb_path=None,
        pretrained_model=None,
        emb=None,
        result_saving_folder="../results",
        name="",
        saving_plots=False,
        label_files=None,
        palette=None,
    ):
        """ImageEmbedding class contains dimensionally reduction, clustering and visualization techniques to analyze patches embeddings. It inherits from Embeddings class.

        Args:
            patches_filenames (_type_, optional): List of all patches paths. Defaults to None.
            patches_info_filename (_type_, optional): Filename that finished with .plk.gz that contains all patches information. Defaults to None.
            emb_path (str, optional): Path to an embedding csv. Defaults to None.
            pretrained_model (PretrainedModel, optional): Pretrained model to compute the embeddings, see models.py. Defaults to None.
            emb (anndata, optional): Embeddings anndata object. Defaults to None.
            result_saving_folder (str, optional): Result folder in which the results are saved. Defaults to "../results".
            name (str, optional): Name of the embeddings. Defaults to "".
            saving_plots (bool, optional): If the plots are saved to the result folder or not. Defaults to False.
            label_files (list, optional): List of files containing label of each spot, "x" column corresponds to the first spot coordinate, "y" column corresponds to the second. Can be csv, tsv with gzip compression or not. One per sample or the name of the file contain the sample name at the beginning with a "_" just after. Defaults to None.
            spot coordinate and "label" column the labels of each spot. Defaults to None.
        """
        Embedding.__init__(self,
            emb=emb,
            result_saving_folder=result_saving_folder,
            name=name,
            saving_plots=saving_plots,
            label_files=label_files,
            palette=palette,
        )
        SpatialViz.__init__(self, emb, saving_plots=saving_plots, result_saving_folder=result_saving_folder)

        self.patches_filenames = patches_filenames
        self.patches_info_filename = patches_info_filename

        self.emb_path = emb_path
        self.pretrained_model = pretrained_model
        self._random_clustering_threshold_invasive_cancer = None
        self._unsupervised_clustering_score_files = []
        self.info_optimal_number_of_clusters = dict()

    def compute_embeddings(self):
        """Compute the image embeddings. Fill emb anndata.

        Raises:
            Exception: If there is not enough information to compute the embeddings
        """


        if (self.pretrained_model is not None and self.patches_filenames is not None) or self.emb_path is not None:
            if self.emb_path is not None:
                embeddings_df = pd.read_csv(self.emb_path, index_col=0)
            else:
                embeddings_df = self.compute_embeddings_df()
            self.emb = ad.AnnData(embeddings_df)
            if self.patches_info_filename is not None:
                patches_info_df = self.load_patches_infos()
                self.emb.obs = self.emb.obs.merge(patches_info_df, right_index=True, left_index=True)
        else:
            raise Exception("Not enough information to compute the image embeddings")

    def compute_embeddings_df(self, batch_size=256, num_workers=0):
        """Compute the image embeddings with the chosen pretrained model.
        Arg:
            batch_size (int, optional): batch size of the data loader

        Returns:
            pd.DataFrame: Embedding dataframe with patches names as index.
        """
        print("Load data...")
        if self.patches_filenames[0].endswith(".hdf5"):
            print("Loading patches from hdf5 file", flush=True)
            images_dataset = ImageEmbeddingDataset(hdf5_file=self.patches_filenames[0], model_name=self.pretrained_model.name)
        else:
            print("Loading patches from folder", flush=True)
            images = pd.DataFrame(sorted(self.patches_filenames), columns=["filename"])
            images_dataset = ImageEmbeddingDataset(images_filenames_pd=images, model_name=self.pretrained_model.name)

        print("Number of images = {}".format(len(images_dataset)))

        database_loader = torch.utils.data.DataLoader(
            images_dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers
        )

        print("Compute features...")
        loaded_pretrained_model = self.pretrained_model.model
        loaded_pretrained_model.eval()
        with torch.no_grad():
            features_list = []
            features_names = []
            for i, batch in enumerate(tqdm(database_loader, desc="Computing features")):

                batch = batch.to(device)
                features = loaded_pretrained_model(batch)
                features = features.cpu().detach().numpy()
                # we can do this as shuffle = false in the dataloader
                for x, y in enumerate(features):
                    features_names.append(images_dataset.name(i * batch_size + x))
                    features_list.append(y.tolist())

        embeddings_df = pd.DataFrame(features_list, index=features_names)
        print("Embedding df shape: {}".format(embeddings_df.shape))
        return embeddings_df
    
    def load_patches_infos(self):
        """Load the patches information from the patches_info_filename attribute.

        Returns:
            pd.DataFrame: Dataframe with the patches information.
        """
        if self.patches_info_filename.endswith(".pkl.gz"):
            with gzip.open(self.patches_info_filename) as file:
                patches_info = pickle.load(file)
            patches_info = (
                pd.DataFrame.from_records(patches_info)
                .sort_values(by="name")
                .reset_index(drop=True)
                .rename(columns={"name": "index"})
                .set_index("index")
            )
        elif self.patches_info_filename.endswith(".csv"):
            patches_info = pd.read_csv(self.patches_info_filename).set_index("index")
        else:
            raise Exception("The file format for self.patches_info_filename is not supported. It needs a .pkl.gz or .csv!")
        return patches_info   


    def plot_representative_samples(self, representative_samples, saving_filename=None):
        """
        Plots representative samples from different clusters.
        Parameters:
        -----------
        representative_samples : dict
            A dictionary where keys are cluster identifiers and values are lists of sample names 
            representing the samples in each cluster.
        Returns:
        --------
        None
        """
        with gzip.open(self.patches_info_filename) as file:
            patches_info = pickle.load(file)
        # Create a dictionary to map patch names to their paths
        patches_dict = {patch['name']: patch['path'] for patch in patches_info}

        # Determine the maximum number of samples in any cluster
        max_samples = max(len(samples) for samples in representative_samples.values())

        # Create a single figure with multiple rows (one for each cluster)
        plt.figure(figsize=(20, 5 * len(representative_samples)))
        # Plot the images
        for row, (key, samples) in enumerate(representative_samples.items()):
            for col in range(max_samples):
                plt.subplot(len(representative_samples), max_samples, row * max_samples + col + 1)
                if col < len(samples):
                    sample = samples[col]
                    img_path = patches_dict.get(sample)
                    if img_path:
                        img = mpimg.imread(img_path)
                        plt.imshow(img)
                        plt.title(sample)
                plt.axis('off')
                
                if col == int(np.floor(max_samples/2) - 1):
                    plt.text(s=f'Cluster {key}', x=0, y=-20, ha='left', va='bottom', fontsize=20, weight='bold')

        plt.tight_layout()
        
        if saving_filename is not None:
            plt.savefig(saving_filename, bbox_inches='tight')
        else:
            plt.show()

