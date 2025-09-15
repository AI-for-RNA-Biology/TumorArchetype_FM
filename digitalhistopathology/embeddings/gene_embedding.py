#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import gzip
import itertools
import math
import os
import glob
import gc
import tempfile
import warnings

import anndata as ad
import mantel
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import scvi
import seaborn as sns
import stlearn as st
import stringdb
import torch
import tqdm

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
from collections import Counter

import matplotlib.pyplot as plt
import plotly.express as px
import rpy2.robjects as robjects
from GraphST import GraphST
from matplotlib_venn import venn2
from PIL import Image
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.packages import importr
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.mixture import BayesianGaussianMixture, GaussianMixture

Image.MAX_IMAGE_PIXELS = None

from digitalhistopathology.embeddings.embedding import Embedding
from digitalhistopathology.helpers import (
    bubbleHeatmap,
    intersection,
    parsing_oldST,
    read_visium,
)

from digitalhistopathology.utils import calculate_mahalanobis_distance, calculate_z_score

import multiprocessing as mp


class GeneEmbedding(Embedding):
    def __init__(
        self,
        samples_names=None,
        genes_count_filenames=None,
        spots_filenames=None,
        spot_diameter_fullres=None,
        image_filenames=None,
        label_files=None,
        emb=None,
        result_saving_folder="../results",
        name="",
        saving_plots=False,
        patches_filenames=None,
        folders_visium=None,
        st_method="old_st",
    ):
        """GeneEmbedding class contains dimensionally reduction, clustering and visualization techniques to analyze genes embeddings. It inherits from Embeddings class.
        Be careful that samples_names, genes_count_filenames, spots_filenames, image_filenames are of the same length and order (each index of each list corresponding to
        a single sample).

        Args:
            samples_names (list): List of all samples names. Defaults to None.
            genes_count_filenames (list): List of all genes counts filenames. One for each sample. For old_st tsv.gz or tsv format. Defaults to None.
            spots_filenames (list, optional): List of the spots locations on the original image filenames. One for each sample. For old_st tsv.gz or tsv format. Defaults to None.
            spot_diameter_fullres (float, optional): Spot diameter on the original images. Defaults to None.
            image_filenames (list, optional): List of the images of origin filenames. One for each sample. Defaults to None.
            label_files (list, optional): List of files containing label of each spot, "x" column corresponds to the first spot coordinate, "y" column corresponds to the second. Can be csv, tsv with gzip compression or not. One per sample or the name of the file contain the sample name at the beginning with a "_" just after. Defaults to None.
            emb (anndata, optional): Embeddings anndata object. Defaults to None.
            result_saving_folder (str, optional): Result folder in which the results are saved. Defaults to "../results".
            name (str, optional): Name of the embeddings. Defaults to "".
            saving_plots (bool, optional): If the plots are saved to the result folder or not. Defaults to False.
            patches_filenames (list, optional): List of the patches filenames corresponding to the spot. Defaults to None.
            folders_visium (str, optional): Path to folder that contains the images for the visium data. Defaults to None.
            st_method (str, optional): Spatial transcriptomics method used, "old_st" or "visium_10x". Defaults to "old_st".
        """
        super().__init__(
            emb=emb,
            result_saving_folder=result_saving_folder,
            name=name,
            saving_plots=saving_plots,
            label_files=label_files,
        )
        self.samples_names = samples_names
        self.genes_count_filenames = genes_count_filenames
        self.spots_filenames = spots_filenames
        self.image_filenames = image_filenames
        self.spot_diameter_fullres = spot_diameter_fullres
        self.patches_filenames = patches_filenames
        self.folders_visium = folders_visium
        self.st_method = st_method
        self.__adjacency_matrix_params = {
            "name": "adjacency_cosine_sim",
            "n_neighbors": None,
            "continuous_adj": False,
            "patches_embeddings": None,
        }
        self.marker_genes_df = None
        self.marker_genes_infos = {"p_val_thresh": None, "logfoldchanges_thresh": None}
        self.enriched_pathways = {"down": None, "up": None}


    def add_label_to_adata(self, data_st, file_label):
        """Add the label to data_st.obs from the file_label.

        Args:
            data_st (anndata): Gene embeddings data for only one sample.
            file_label (str): Label file corresponding to this sample.

        Returns:
            anndata: Gene embeddings data from the sample with label in .obs.
        """
        print("Add label to data")

        if "her2" in self.name.lower():
            compression = "gzip" if file_label.endswith("gz") else None
            sep = "\t" if ".tsv" in file_label.split("/")[-1] else ","
            label_df = pd.read_csv(file_label, sep=sep, compression=compression)
            label_df = label_df.dropna(axis=0)

            label_df["x"] = label_df["x"].apply(lambda x: round(x))
            label_df["y"] = label_df["y"].apply(lambda y: round(y))

            # Problem in some immune infiltrate labels
            label_df.loc[label_df["label"] == "immune infiltrate⁄", "label"] = "immune infiltrate"

            data_st.obs["x"] = [int(ind.split("x")[0]) for ind in list(data_st.obs_names)]
            data_st.obs["y"] = [int(ind.split("x")[-1]) for ind in list(data_st.obs_names)]

            data_st.obs = (
                data_st.obs.reset_index()
                .merge(label_df[["x", "y", "label"]], how="left", on=["x", "y"])
                .drop(columns=["x", "y"])
                .set_index("index")
            )
            assert data_st.obs["label"].isna().sum() == 0, "Problem: some spots do not have labels"

        elif "dlpfc" in self.name.lower():
            df_meta = pd.read_csv(file_label, sep="\t")
            data_st.obs["label"] = df_meta["layer_guess"].values
        elif "tnbc" in self.name.lower():
            df_meta = pd.read_csv(file_label)
            data_st.obs["label"] = df_meta["label"].values
        else:
            warnings.warn("No implementation to add label to this dataset")

        return data_st

    def add_patches_path_to_adata(self, gene_emb):
        """Add patches path to gene_emb.obs.

        Args:
            gene_emb (anndata): Genes embeddings anndata.
        """
        print("Add patch path to data")
        df_filenames = pd.DataFrame()
        df_filenames["path"] = self.patches_filenames
        df_filenames["name_origin"] = df_filenames["path"].apply(lambda f: f.split("/")[-1].split("_spot")[0])
        # In patches names that I added rep_ from origin name
        if self.name is not None and "her2" in self.name.lower():
            df_filenames["name_origin"] = df_filenames["name_origin"].apply(lambda o: o.split("_rep")[0] + o.split("_rep")[1])
        print(f"Path: {df_filenames['path'].values[0].split('/')[-1]}", flush=True)
        print(df_filenames.columns, flush=True)
        
        if df_filenames['path'].values[0].endswith('hdf5'):
            
            df_filenames["x"] = "none"
            df_filenames["y"] = "none"
            
            gene_emb.obs['path'] = df_filenames['path'].values[0]
            
        else:
            
            df_filenames["x"] = df_filenames["path"].apply(lambda f: int(f.split("/")[-1].split("_spot")[-1].split("x")[0]))
            
            df_filenames["y"] = df_filenames["path"].apply(
                lambda f: int(f.split("/")[-1].split("_spot")[-1].split("x")[1].split(".")[0])
            )
            gene_emb.obs["path"] = gene_emb.obs.merge(df_filenames, how="left", on=["name_origin", "x", "y"])["path"].values

    def read_old_st_with_gzip_compression(self, index):
        """Read old ST sample (Ex: her2-positive breast cancer data).

        Args:
            index (int): Index of the sample from the dataset you want to read, >=0.

        Returns:
            anndata: Anndata of the sample with gene count in X and additional information in obs.
        """
        genes_count_file = self.genes_count_filenames[index]
        spots_file = self.spots_filenames[index]
        image_file = self.image_filenames[index]
        # Maybe not a label file for all samples
        if self.label_files is not None:
            try:
                label_file = [f for f in self.label_files if self.samples_names[index] in f.split("/")[-1]][0]
            except Exception as e:
                print("No label file was found for sample {}: {}".format(self.samples_names[index], e))
                label_file = None

        # Automatically remove temp_dir outside of the with
        with tempfile.TemporaryDirectory() as temp_dir:
            if genes_count_file.endswith(".gz"):
                temp_file_path_count = os.path.join(temp_dir, "tempfile_count")
                # Decompress the Gzip file and save it to the temporary file
                with gzip.open(genes_count_file, "rb") as gz_file, open(temp_file_path_count, "wb") as temp_file:
                    temp_file.write(gz_file.read())
            else:
                temp_file_path_count = genes_count_file

            if spots_file.endswith(".gz"):
                temp_file_path_spatial = os.path.join(temp_dir, "tempfile")
                with gzip.open(spots_file, "rb") as gz_file, open(temp_file_path_spatial, "wb") as temp_file:
                    temp_file.write(gz_file.read())
            else:
                temp_file_path_spatial = spots_file

            # st.ReadOldST filters all genes whose total gene count is 0, so I adapted the parsing function
            data_st = sc.read_text(temp_file_path_count)
            if genes_count_file.endswith(".csv"):
                data_st = sc.read_csv(genes_count_file)
            else:
                data_st = sc.read_text(temp_file_path_count)
                
            data_st = parsing_oldST(data_st, coordinates_file=temp_file_path_spatial)
            st.add.image(
                data_st,
                library_id="OldST",
                quality="hires",
                imgpath=image_file,
                scale=1.0,
                spot_diameter_fullres=self.spot_diameter_fullres,
            )

        if label_file is not None:
            data_st = self.add_label_to_adata(data_st, label_file)
        else:
            data_st.obs["label"] = np.nan

        # Add some other useful informations
        data_st.obs["x"] = [int(ind.split("x")[0]) for ind in list(data_st.obs_names)]
        data_st.obs["y"] = [int(ind.split("x")[-1]) for ind in list(data_st.obs_names)]
        data_st.obs["name_origin"] = self.samples_names[index]
        data_st.obs["path_origin"] = self.image_filenames[index]
        data_st.obs["total_genes_count"] = data_st.X.sum(axis=1)
        data_st.obs_names = data_st.obs.apply(lambda r: "spot{}".format(r.name), axis=1)

        if self.patches_filenames is not None:
            self.add_patches_path_to_adata(data_st)

        return data_st

    def read_visium_10x(self, index, image_info_folder=""):
        """Read 10x visium data sample.

        Args:
            index (int): Index of the sample from the dataset you want to read, >=0.
            image_info_folder (str, optional): Folder in which there is the visium images. Defaults to "".

        Returns:
            anndata: Anndata of the sample with gene count in X and additional information in obs.
        """
        data = read_visium(
            self.folders_visium[index],
            count_file=self.genes_count_filenames[index].split("/")[-1],
            image_info_folder=image_info_folder,
        )
        data.var_names_make_unique()
        data.X = data.X.toarray()
        data.obs["imagecol"] = data.obsm["spatial"][:, 0]
        data.obs["imagerow"] = data.obsm["spatial"][:, 1]
        data.obs = data.obs.rename(columns={"array_col": "x", "array_row": "y"})
        data.obs["index"] = data.obs.apply(lambda r: "spot{}x{}".format(int(r["x"]), int(r["y"])), axis=1)
        data.obs = data.obs.reset_index(drop=True).set_index("index")
        data.obs["name_origin"] = self.samples_names[index]
        data.obs["path_origin"] = self.image_filenames[index]
        data.obs["total_genes_count"] = data.X.sum(axis=1)

        if self.patches_filenames is not None:
            self.add_patches_path_to_adata(data)

        # Maybe not a label file for all samples
        if self.label_files is not None:
            try:
                label_file = [f for f in self.label_files if self.samples_names[index] in f.split("/")[-2]][0]
            except Exception as e:
                print("No label file was found for sample {}: {}".format(self.samples_names[index], e))
                label_file = None

        if label_file is not None:
            data = self.add_label_to_adata(data, label_file)
        else:
            data.obs["label"] = np.nan

        return data

    def get_anndata(self, index, whole_data=True):
        """Get the anndata sample from index irrespective of the method.

        Args:
            index (int): Index of the sample from the dataset you want to read, >=0.
            whole_data (bool, optional): If you want the whole anndata (with obs, ...) or a simple anndata with only count in X. Defaults to True.

        Raises:
            Exception: Not a correct spatial transcriptomics method set: old_st or visium_10x

        Returns:
            anndata: Anndata of the sample with gene count in X and additional information in obs if whole_data is set to True.
        """
        print("Read data for sample {}".format(self.samples_names[index]))
        if whole_data:
            if self.st_method == "old_st":
                adata = self.read_old_st_with_gzip_compression(index=index)
            elif self.st_method == "visium_10x":
                adata = self.read_visium_10x(index=index)
            else:
                raise Exception("Not a correct spatial transcriptomics method: old_st or visium_10x")

        else:
            if self.st_method == "old_st":
                delimiter = "\t" if ".tsv" in self.genes_count_filenames[index].split("/")[-1] else ","
                adata = sc.read(self.genes_count_filenames[index], delimiter=delimiter)
            elif self.st_method == "visium_10x":
                adata = sc.read_10x_h5(self.genes_count_filenames[index])
                adata.X = adata.X.toarray()
            else:
                raise Exception("Not a correct spatial transcriptomics method: old_st or visium_10x")
        return adata


    def sample_threshold(
        self,
        index=None,
        adata=None,
        quantile_threshold=0.1,
        intersection=True,
        density=False,
        bins=50,
        baysian_gmm=True,
        scaling_to_max=True,
        show_plot=True,
        matrix_name="",
        ax=None,
    ):
        """Compute the log2(bulk genes counts + 1) threshold for a sample by fitting a gaussian mixture model to identify a bimodal behavior in the distribution.

        Args:
            index (int): Index of the sample from the dataset you want to read, >=0. Defaults to None.
            adata (anndata): Anndata from which we want to plot bulk (samples x genes). Defaults to None.
            quantile_threshold (float, optional): Quantile threshold to apply on the foreground gaussian. Defaults to 0.1.
            intersection (bool, optional): To display the intersection between the two gaussian distributions. Defaults to True.
            density (bool, optional): If the histogram y axis is in density or count. Defaults to False.
            bins (int, optional): Number of bins for the histogram. Defaults to 50.
            baysian_gmm (bool, optional): To use BayesianGaussianMixture instead of GaussianMixture. Defaults to True.
            scaling_to_max (bool, optional): To scale the gaussians to their max between mean-std and mean+std. If False it is scaled to their median between mean-std and mean+std. Defaults to True.
            show_plot (bool, optional): If you want to show the plot or not. Defaults to True.
            matrix_name (str, optional): Name of the data that appear in the title of the plot. Defaults to "".
            ax (matplotlib.ax, optional): Possibility to plot the figure on a axe. Defaults to None.

        Returns:
            float: Log2(bulk genes counts + 1) threshold
        """
        print(
            "Start computing threshold for {}".format(
                "sample " + self.samples_names[index] if index is not None else matrix_name
            )
        )

        if adata is None:
            adata = self.get_anndata(index=index, whole_data=False)

        if baysian_gmm:
            gmm = BayesianGaussianMixture(n_components=2, random_state=0, max_iter=1000)
        else:
            gmm = GaussianMixture(n_components=2, random_state=0, max_iter=500)

        adata.var["bulk"] = np.reshape(adata.X.sum(axis=0), (adata.shape[1], 1))
        adata.var["bulk_log"] = np.log2(adata.var["bulk"].values + 1)

        gmm.fit(np.reshape(adata.var["bulk_log"].values, (adata.var["bulk_log"].values.shape[0], 1)))

        ax_none = ax is None
        if ax_none:
            fig, ax = plt.subplots(figsize=(6, 5))

        y_hist, x_hist, _ = ax.hist(
            adata.var["bulk_log"].values,
            density=density,
            alpha=0.5,
            color="gray",
            bins=bins,
        )
        x_hist = (x_hist[1:] + x_hist[:-1]) / 2
        hist_df = pd.DataFrame(np.array([y_hist, x_hist]).T, columns=["y", "x"])

        x = np.linspace(adata.var["bulk_log"].min(), adata.var["bulk_log"].max(), 1000).reshape(-1, 1)

        if gmm.means_[0][0] > gmm.means_[1][0]:
            ind_list = [1, 0]
        else:
            ind_list = [0, 1]

        pdfs = {}
        for i, ind in enumerate(ind_list):
            mean = gmm.means_[ind][0]
            covariance = gmm.covariances_[ind][0, 0]
            dist = scipy.stats.norm(loc=mean, scale=np.sqrt(covariance))
            pdf = dist.pdf(x)

            # quantile_threshold on the foreground
            if i == 1:
                threshold = dist.ppf(quantile_threshold)
                color = "green"
            else:
                color = "orange"

            gauss_max = np.max(pdf)
            hist_max = hist_df[
                np.logical_and(
                    hist_df["x"] > mean - np.sqrt(covariance),
                    hist_df["x"] < mean + np.sqrt(covariance),
                )
            ]["y"]
            if scaling_to_max:
                hist_max = hist_max.max()
            else:
                hist_max = np.nanmean(hist_max.apply(lambda y: float("nan") if y == 0 else y))
            scaling_factor = hist_max / gauss_max
            ax.plot(x, pdf * scaling_factor, label=f"Component {i + 1}", color=color)
            pdfs[f"Component {i + 1}"] = np.reshape(pdf * scaling_factor, -1)

        ax.axvline(
            threshold,
            color="red",
            label=f"{quantile_threshold*100:.0f}th: {threshold:.3f}",
        )
        print("Threshold based on the {}th quantile: {}".format(quantile_threshold * 100, np.round(threshold, 3)))

        if intersection:
            intersection_ind = np.argwhere(np.diff(np.sign(pdfs["Component 1"] - pdfs["Component 2"]))).flatten()
            intersection_point = x[intersection_ind].flatten()
            ax.scatter(
                x[intersection_ind],
                pdfs["Component 1"][intersection_ind],
                c="red",
                label="Intersection: {}".format(np.round(intersection_point.max(), 3)),
                marker="X",
            )
            print("Intersection x point: {}".format(intersection_point[0]))

        ax.legend()
        ax.set_xlabel("Log2(total number of reads + 1)")
        if density:
            ax.set_ylabel("Genes density")
        else:
            ax.set_ylabel("Genes count")

        ax.set_title(
            "Gaussian mixture modeling of the total number of reads per gene of {}".format(
                "sample " + self.samples_names[index] if index is not None else matrix_name
            )
        )

        if self.saving_plots:
            plt.savefig(
                os.path.join(
                    self.result_saving_folder,
                    "gene_count_threshold_bimodal_fit_{}.pdf".format(
                        "sample" + self.samples_names[index] if index is not None else matrix_name.replace(" ", "_")
                    ),
                ),
                bbox_inches="tight",
            )
        if not ax_none:
            return ax
        else:
            if show_plot:
                plt.show()
            else:
                plt.close()
            return threshold

    def filter_genes(
        self, adata, quantile_threshold=0.1, show_plot=False, matrix_name="", baysian_gmm=True, scaling_to_max=True
    ):
        """Filter genes according to the pseudo-bulk bimodal fitting from an andata genes counts object.

        Args:
            adata (anndata): Anndata genes counts object.
            quantile_threshold (float, optional): Quantile threshold to apply on the foreground gaussian. Defaults to 0.1.
            show_plot (bool, optional): If you want to show the bimodal fitting plot or not. Defaults to False.
            baysian_gmm (bool, optional): To use BayesianGaussianMixture instead of GaussianMixture. Defaults to True.
            scaling_to_max (bool, optional): To scale the gaussians to their max between mean-std and mean+std. If False it is scaled to their median between mean-std and mean+std. Defaults to True.
            matrix_name (str, optional): Name of the data that appear in the title of the plot. Defaults to "".

        Returns:
            anndata: Anndata genes counts object without filtered genes.
        """
        print("Filering genes")
        threshold = self.sample_threshold(
            adata=adata,
            quantile_threshold=quantile_threshold,
            show_plot=show_plot,
            matrix_name=matrix_name,
            baysian_gmm=baysian_gmm,
            scaling_to_max=scaling_to_max,
        )
        print("Shape before filtering genes = {}".format(adata.X.shape))
        adata.var["bulk"] = np.reshape(adata.X.sum(axis=0), (adata.shape[1], 1))
        adata.var["bulk_log"] = np.log2(adata.var["bulk"].values + 1)
        discard_genes = adata.var["bulk_log"].apply(lambda e: e < threshold)
        percentage_discarded_genes = (discard_genes).sum() / len(discard_genes) * 100
        print("Percentage of discarded genes: {}%".format(round(percentage_discarded_genes, 3)))
        adata = adata[:, ~discard_genes]
        print("Shape after filtering genes = {}".format(adata.X.shape))
        adata.uns["filter_genes_log2_threshold"] = threshold
        return adata

    # def compute_raw_embeddings(self):
    #     """Fill emb with raws genes counts that are aligned across the different samples."""
    #     genes_present_df = None
    #     for i in range(len(self.genes_count_filenames)):
    #         gene_emb = self.get_anndata(i, whole_data=True)
    #         genes_list = list(gene_emb.var_names)
    #         del gene_emb.uns["spatial"]

    #         if "her2" in self.name.lower():
    #             gene_emb.obs["index"] = gene_emb.obs.apply(
    #                 lambda r: r["name_origin"][0] + "_rep{}_".format(r["name_origin"][1]) + str(r.name), axis=1
    #             )
    #             gene_emb.obs["tumor"] = gene_emb.obs["name_origin"].apply(lambda n: n[0])
    #         else:
    #             gene_emb.obs["index"] = gene_emb.obs.apply(lambda r: r["name_origin"] + "_" + str(r.name), axis=1)

    #         gene_emb.obs = gene_emb.obs.reset_index(drop=True).set_index("index")

    #         if i == 0:
    #             self.emb = gene_emb.copy()
    #             genes_present_df = pd.DataFrame(index=genes_list)
    #             genes_present_df[self.samples_names[i]] = np.ones(len(genes_list))
    #         else:
    #             # concat X, obs, uns, obsm (var deleted)
    #             self.emb = ad.concat([self.emb, gene_emb], axis=0, join="outer", uns_merge="first")
    #             # replace nan with 0
    #             self.emb.X = np.nan_to_num(self.emb.X, nan=0.0)
    #             current_genes_present_df = pd.DataFrame(index=genes_list)
    #             current_genes_present_df[self.samples_names[i]] = np.ones(len(genes_list))
    #             genes_present_df = genes_present_df.merge(
    #                 current_genes_present_df, left_index=True, right_index=True, how="outer"
    #             )

    #     self.emb.var = pd.concat((self.emb.var, genes_present_df.fillna(0).astype(bool)), axis=1)
    
    def compute_raw_embeddings(self, chunk_size=10, temp_dir="./tmp_gene_emb_chunks"):
        
        from gtfparse import read_gtf

        # gencode = read_gtf("/idiap/group/genomics/annotation/hg38/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf")

        os.makedirs(temp_dir, exist_ok=True)
        genes_present_df = None
        chunk_embs = []
        for i in range(len(self.genes_count_filenames)):
            gene_emb = self.get_anndata(i, whole_data=True)
            genes_list = list(gene_emb.var_names)
            
            # genes_list = list(set(gencode['gene_id'].unique()) & set(gene_emb.var_names))
            # gene_emb = gene_emb[:, genes_list]
            
            del gene_emb.uns["spatial"]

            if "her2" in self.name.lower():
                gene_emb.obs["index"] = gene_emb.obs.apply(
                    lambda r: r["name_origin"][0] + "_rep{}_".format(r["name_origin"][1]) + str(r.name), axis=1
                )
                gene_emb.obs["tumor"] = gene_emb.obs["name_origin"].apply(lambda n: n[0])
            else:
                gene_emb.obs["index"] = gene_emb.obs.apply(lambda r: r["name_origin"] + "_" + str(r.name), axis=1)

            gene_emb.obs = gene_emb.obs.reset_index(drop=True).set_index("index")

            chunk_embs.append(gene_emb)
            if genes_present_df is None:
                genes_present_df = pd.DataFrame(index=genes_list)
            else:
                # Add any new genes to the index
                new_genes = set(genes_list) - set(genes_present_df.index)
                if new_genes:
                    genes_present_df = genes_present_df.reindex(
                        list(genes_present_df.index) + list(new_genes),
                        fill_value=0
                    )
            genes_present_df[self.samples_names[i]] = pd.Series(1, index=genes_list)

            # Save chunk and clear memory
            if (i + 1) % chunk_size == 0 or (i + 1) == len(self.genes_count_filenames):
                chunk_file = os.path.join(temp_dir, f"chunk_{i // chunk_size}.h5ad")
                ad.concat(chunk_embs, axis=0, join="outer", uns_merge="first").write(chunk_file)
                chunk_embs = []
                gc.collect()

        # Merge all chunks at the end
        chunk_files = sorted(glob.glob(os.path.join(temp_dir, "chunk_*.h5ad")))
        print(f"Merging {len(chunk_files)} chunks from {temp_dir}...")
        all_embs = [ad.read_h5ad(f) for f in chunk_files]
        print(f"Total number of chunks: {len(all_embs)}")
        self.emb = ad.concat(all_embs, axis=0, join="outer", uns_merge="first")
        print(f"Shape of concatenated embeddings: {self.emb.shape}")
        self.emb.X = np.nan_to_num(self.emb.X, nan=0.0)
        print("Replacing NaN values with 0 in embeddings matrix.")
        self.emb.var = pd.concat((self.emb.var, genes_present_df.fillna(0).astype(bool)), axis=1)
        print("Gene embeddings computed and stored in self.emb.")
        
    def preprocessing(
        self,
        filter_genes_groupby="tumor",
        filter_genes_fraction_of_presence_across_group=0.5,
        filter_genes_intersection_groupby=None,
        filter_empty_spots=True,
        spot_norm=False,
        spot_norm_by="name_origin",
        log=True,
        scaling=False,
        adata=None,
        inplace=False,
    ):
        """Preprocessing of raw genes counts data.

        Args:
            filter_genes_groupby (str, optional): Obs column with which we group data to filter genes (bimodal filtering). Defaults to "tumor".
            filter_genes_fraction_of_presence_across_group (float, optional): Fraction of groupby value (ex: patient) that have this genes expressed after pseudo-bulk bimodal filtering. Defaults to 0.5.
            filter_genes_intersection_groupby (str, optional): In each group, do an intersection of the reliably expressed genes of each subgroup defined by intersection_groupby which is a self.emb.obs column. For example, if groupby is tumor (patient) you can decide to do the intersection in each group based on the replicates (name_origin). If None, it is not used.
            filter_empty_spots (bool, optional): If the filtering of empty spots (rows) is done. Defaults to True.
            spot_norm (bool, optional): If the normalization with target sum of 1e6 of spots (rows) is done. Defaults to False.
            spot_norm_by (str, optional): If spot_norm is True, if the spot normalization is individually performed on groups of the data (emb.obs column). If None, spot normalization is performed on the whole data. Defaults to name_origin.
            log (bool, optional): If the log2(x + 1) transformed is applied. Defaults to True.
            scaling (bool, optional): If the scaling is applied with zero_centered=False, omit zero-centering variables, which allows to handle sparse input efficiently. Defaults to False.
            inplace (bool, optional): If the self.emb is replaced by the preprocessed emb in this function. Defaults to False.
        """
        if adata is None:
            adata = self.emb.copy()

        reliably_expressed_genes = self.get_reliably_expressed_genes(
            fraction_of_presence_across_patient_threshold=filter_genes_fraction_of_presence_across_group,
            groupby=filter_genes_groupby,
            show_plot=False,
            intersection_groupby=filter_genes_intersection_groupby,
        )
        adata = adata[:, reliably_expressed_genes]
        print("Filtering genes done, shape = {}".format(adata.shape))

        if filter_empty_spots:
            adata = adata[(adata.X.sum(axis=1) != 0), :]
            print("\nShape after empty spots filtering: {}".format(adata.shape))

        if spot_norm:
            if spot_norm_by is not None:
                groups = adata.obs.groupby(by=spot_norm_by)
                for name, group in groups:
                    sub_adata = adata[group.index, :]
                    sc.pp.normalize_total(sub_adata, target_sum=1e6)
                    adata[group.index, :].X = np.array(sub_adata.X)
            else:
                sc.pp.normalize_total(adata, target_sum=1e6)
            print("\nSpots normalization done.")

        if log:
            st.pp.log1p(adata, base=2)
            print("\nLog transformed done.")

        if scaling:
            # If zero_centered=False, omit zero-centering variables, which allows to handle sparse input efficiently.
            sc.pp.scale(adata, zero_center=False)
            print("\nScaling done.")

        adata.uns.update({"preprocessing": True})
        adata.uns.update(
            {
                "filter_genes_groupby": filter_genes_groupby,
                "filter_genes_fraction_of_presence_across_group": filter_genes_fraction_of_presence_across_group,
                "filter_empty_spots": filter_empty_spots,
                "spot_norm": spot_norm,
                "spot_norm_by": spot_norm_by,
                "log": log,
                "scaling": scaling,
            }
        )

        if inplace:
            self.emb = adata
        else:
            return adata

    def compute_embeddings(self):
        """Compute gene embeddings for each spot of each sample of the dataset. Fill emb anndata."""
        print("Start computing genes embeddings")
        self.compute_raw_embeddings()


    def get_reliably_expressed_genes(
        self,
        fraction_of_presence_across_patient_threshold=0.5,
        quantile_threshold=0.1,
        groupby="tumor",
        show_plot=False,
        baysian_gmm=True,
        scaling_to_max=True,
        intersection_groupby=None,
        threshold_subgroup_size=0,
        OR_groupby_gene=True,
    ):
        """Compute the reliably expressed genes across the samples to filter genes.

        Args:
            fraction_of_presence_across_patient_threshold (float, optional): Fraction of patient that have this genes expressed after pseudo-bulk bimodal filtering. Defaults to 0.75.
            quantile_threshold (float, optional): Quantile threshold to apply on the foreground gaussian. Defaults to 0.1.
            groupby (str, optional): How to group the data to filter it. Defaults to "tumor".
            show_plot (bool, optional): If you want to show the bimodal fitting plot or not. Defaults to False.
            baysian_gmm (bool, optional): To use BayesianGaussianMixture instead of GaussianMixture. Defaults to True.
            scaling_to_max (bool, optional): To scale the gaussians to their max between mean-std and mean+std. If False it is scaled to their median between mean-std and mean+std. Defaults to True.
            intersection_groupby (str, optional): In each group, do an intersection of the reliably expressed genes of each subgroup defined by intersection_groupby which is a self.emb.obs column. For example, if groupby is tumor (patient) you can decide to do the intersection in each group based on the replicates (name_origin). If None, it is not used.
            threshold_subgroup_size (int, optional): The number of samples a subgroup need to have to be counted. Default to 0.
            OR_groupby_gene (bool, optional): If the reliably expressed genes of each groupby group are fused with OR. If False, AND is used. Defaults to True.

        Returns:
            list: List with all the reliably expressed genes
        """
        log2_threshold_per_group = dict()
        reliably_expressed_gene_per_group = dict()
        groups = self.emb.obs.groupby(by=groupby)
        genes_list = []
        for name, group in groups:
            if intersection_groupby is not None:
                reliably_expressed_gene_per_group[name] = dict()
                log2_threshold_per_group[name] = dict()
                # Intersection between same patient slides or same cluster patients
                subgroups = group.groupby(by=intersection_groupby)
                sub_genes_list = []
                for subname, subgroup in subgroups:
                    if subgroup.shape[0] < threshold_subgroup_size:
                        print(
                            "Do not count for group {} x {} because it has less than {} samples".format(
                                name, subname, threshold_subgroup_size
                            )
                        )
                        continue
                    # Get only the genes present in the specific sample file
                    if intersection_groupby == "name_origin":
                        subdata = self.emb[subgroup.index, self.emb.var[subname]]
                    else:
                        subdata = self.emb[subgroup.index, :]
                    subdata = self.filter_genes(
                        subdata,
                        quantile_threshold=quantile_threshold,
                        show_plot=show_plot,
                        matrix_name="{} {} x {} {}".format(
                            groupby if groupby != "predicted_label" else "cluster",
                            name,
                            intersection_groupby if intersection_groupby != "tumor" else "patient",
                            subname,
                        ),
                        baysian_gmm=baysian_gmm,
                        scaling_to_max=scaling_to_max,
                    )
                    reliably_expressed_gene_per_group[name][subname] = len(list(subdata.var_names))
                    log2_threshold_per_group[name][subname] = subdata.uns["filter_genes_log2_threshold"]
                    if len(sub_genes_list) == 0:
                        sub_genes_list = list(subdata.var_names)
                    else:
                        sub_genes_list = list(set(sub_genes_list) & set(list(subdata.var_names)))

                if OR_groupby_gene or len(genes_list) == 0:
                    genes_list.extend(sub_genes_list)
                else:
                    genes_list = list(set(genes_list) & set(list(sub_genes_list)))

                reliably_expressed_gene_per_group[name]["all"] = len(sub_genes_list)
            else:
                if name != "undetermined":
                    adata = self.filter_genes(
                        self.emb[group.index, :],
                        quantile_threshold=quantile_threshold,
                        show_plot=show_plot,
                        matrix_name=name,
                        baysian_gmm=baysian_gmm,
                        scaling_to_max=scaling_to_max,
                    )
                    if OR_groupby_gene or len(genes_list) == 0:
                        genes_list.extend(list(adata.var_names))
                    else:
                        genes_list = list(set(genes_list) & set(list(adata.var_names)))
                    reliably_expressed_gene_per_group[name] = len(list(adata.var_names))
                    log2_threshold_per_group[name] = adata.uns["filter_genes_log2_threshold"]

        reliably_expressed_gene_per_group["all_before_patient_thresholding"] = len(set(genes_list))
        counter_genes_df = pd.DataFrame.from_dict(dict(Counter(genes_list)), orient="index", columns=["count"])
        print("Length before number of patient filtering = {}".format(len(set(genes_list))))
        threshold = fraction_of_presence_across_patient_threshold * len(groups)
        print("Threshold = {}".format(threshold))
        reliably_expressed_genes = sorted(list(counter_genes_df[counter_genes_df["count"] >= threshold].index))
        print("Length after number of patient filtering = {}".format(len(reliably_expressed_genes)))

        reliably_expressed_gene_per_group["all"] = len(reliably_expressed_genes)
        self.emb.uns["reliably_expressed_gene_per_group"] = reliably_expressed_gene_per_group
        print("Final number of identified reliable genes = {}".format(len(reliably_expressed_genes)))

        self.emb.uns["filter_genes_log2_threshold_per_group"] = log2_threshold_per_group

        return reliably_expressed_genes


    @staticmethod
    def get_z_score_one_pathway(df_spots, GO_genes, df_background=None, mahalanobis=False):
        
        # print(f"Calculating z-score for pathway with {len(GO_genes)} genes", flush=True)
        
        if df_background is None:
            df_background = df_spots[[g for g in GO_genes if g in df_spots.columns]]
        else:
            df_background = df_background[[g for g in GO_genes if g in df_background.columns and g in df_spots.columns]]
        

        spots_pathway = df_spots[[g for g in GO_genes if g in df_spots.columns]]        
        # print(f"Spots pathway shape: {spots_pathway.shape}", flush=True)
        ## Positive values are up-regulated, negative values are down-regulated
        
        
        if mahalanobis:
            print("Calculating mahalanobis distance", flush=True)
            results = spots_pathway.apply(lambda col: calculate_mahalanobis_distance(point=col, data=df_background, both_sided=True), axis=1)
        else:
            results = spots_pathway.apply(lambda col: calculate_z_score(point=col, data=df_background), axis=1)
            
        return results
    
    @staticmethod
    def compute_z_score(row, df_spots, genes_list_column='GO_genes', term_name_column='term_name', df_background=None, mahalanobis=False):
        GO_genes = row[genes_list_column]
        term = row[term_name_column]
        print(f"Calculating z-score for pathway {term}", flush=True)
        z_scores = GeneEmbedding.get_z_score_one_pathway(df_spots, GO_genes, df_background, mahalanobis)
        print(f"Finished calculating z-score for pathway {term}", flush=True)
        return term, z_scores
         
    
    def get_spots_z_scores_multiple_pathways(self, 
                                                 df_pathways, 
                                                 term_name_column='term_name',
                                                 genes_list_column='GO_genes',
                                                 num_processes=1,
                                                 df_background=None, 
                                                 mahalanobis=False):

        df_spots = self.emb.to_df()

        with mp.Pool(processes=num_processes) as pool:
            results = [pool.apply_async(GeneEmbedding.compute_z_score, args=(row, df_spots, genes_list_column, term_name_column, df_background, mahalanobis)) for row in df_pathways.to_dict(orient='records')]
            pool.close()
            pool.join()

        results = [result.get() for result in results]
        distances = pd.DataFrame.from_dict({term_name: distance for term_name, distance in results}, orient='index')
        return distances

    
    

    @staticmethod
    def get_spots_log10_pvalue_one_pathway(df_spots, GO_genes, df_background=None):
        
        if df_background is None:
            mean_spots = df_spots.median(axis=0)
        else:
            mean_spots = df_background.median(axis=0)

        spots_pat1 = df_spots[[g for g in GO_genes if g in df_spots.columns]]
        
        ## Positive values are up-regulated, negative values are down-regulated
        
        

        results = spots_pat1.apply(lambda col: -np.log10(scipy.stats.ttest_ind(col, mean_spots, equal_var=False)[1]), axis=1)

        return results
    

    @staticmethod
    def compute_log10_pvalue(row, df_spots, genes_list_column='GO_genes', term_name_column='term_name', df_background=None):
        GO_genes = row[genes_list_column]
        return row[term_name_column], GeneEmbedding.get_spots_log10_pvalue_one_pathway(df_spots, GO_genes, df_background)


    def get_spots_log10_pvalue_multiple_pathways(self, 
                                                 df_pathways, 
                                                 term_name_column='term_name',
                                                 genes_list_column='GO_genes',
                                                 num_processes=1,
                                                 df_background=None):

        df_spots = self.emb.to_df()

        def update_progress(result):
            nonlocal progress
            progress.update()

        with mp.Pool(processes=num_processes) as pool:
            progress = tqdm.tqdm(total=len(df_pathways), desc="Processing pathways")
            results = []
            for row in df_pathways.to_dict(orient='records'):
                result = pool.apply_async(GeneEmbedding.compute_log10_pvalue, args=(row, df_spots, genes_list_column, term_name_column, df_background), callback=update_progress)
                results.append(result)
            pool.close()
            pool.join()

        results = [result.get() for result in results]
        pvalues = {term_name: log10_pvalue for term_name, log10_pvalue in results}  
        log10_pvalues = pd.DataFrame(pvalues)
        log10_pvalues = pd.DataFrame.from_dict(pvalues, orient='index')

        return log10_pvalues
