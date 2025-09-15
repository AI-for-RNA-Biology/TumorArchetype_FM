#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import anndata as ad
import matplotlib
import pandas as pd

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

from PIL import Image
from sklearn.preprocessing import scale

Image.MAX_IMAGE_PIXELS = None

from digitalhistopathology.embeddings.embedding import Embedding


class MergedEmbedding(Embedding):
    def __init__(
        self,
        gene_emb,
        img_emb,
        ef_emb=None,
        emb=None,
        result_saving_folder="../results/embeddings/merged_embeddings",
        name="",
        saving_plots=False,
        method="concat",
    ):
        """MergedEmbedding class contains dimensionally reduction, clustering and visualization techniques to analyze merged embeddings. It inherits from Embeddings class. It merged embeddings from genes and images originated from the same spot (location).

        Args:
            gene_emb (anndata): Genes embeddings.
            img_emb (anndata): Patches/images embeddings.
            emb (anndata, optional): Embeddings anndata object. Defaults to None.
            result_saving_folder (str, optional): Result folder in which the results are saved. Defaults to "../results/merged_embeddings".
            name (str, optional): Name of the embeddings. Defaults to "".
            saving_plots (bool, optional): If the plots are saved to the result folder or not. Defaults to False.
        """
        super().__init__(
            emb=emb,
            result_saving_folder=result_saving_folder,
            name=name,
            saving_plots=saving_plots,
        )
        self.gene_emb = gene_emb
        self.img_emb = img_emb
        self.ef_emb = ef_emb
        self.method = method
        if emb is not None:
            self.emb = emb

    def prepare_data(self, only_labeled_spots=False, only_cancer=False):
        """Align and prepare the data. All embeddings have the same number of samples in the same order.

        Args:
            only_labeled_spots (bool, optional): If you want only data with a label entry in .obs. Defaults to False.
            only_cancer (bool, optional): If you want only data with a label entry in .obs that is "invasive cancer" or "cancer in situ". Defaults to False.
        """
        # Alignement and prepare the data (standardization)
        df, _ = self.img_emb.obs.align(self.gene_emb.obs, join="inner", axis=0)
        if self.ef_emb is not None:
            df, _ = df.align(self.ef_emb.obs, join="inner", axis=0)
        if only_labeled_spots:
            df = df.loc[~df["label"].isna(), :]
        if only_cancer:
            df = df.loc[df["label"].apply(lambda l: l in ["invasive cancer", "cancer in situ"]), :]
        index = df.index
        self.img_emb = self.img_emb[index, :]
        self.img_emb.X = scale(self.img_emb.X, axis=1)
        self.gene_emb = self.gene_emb[index, :]
        self.gene_emb.X = scale(self.gene_emb.X, axis=1)
        if self.ef_emb is not None:
            self.ef_emb = self.ef_emb[index, :]
            self.ef_emb.X = scale(self.ef_emb.X, axis=1)
            assert self.img_emb.shape[0] == self.ef_emb.shape[0], "Problem of shape"

        assert self.img_emb.shape[0] == self.gene_emb.shape[0], "Problem of shape"

        # Prepare image emb
        self.img_emb = self.img_emb[~self.img_emb.obs["spots_info"].isna()]
        if "x" not in self.img_emb.obs.columns or "y" not in self.img_emb.obs.columns:
            self.img_emb.obs = pd.concat(
                (
                    self.img_emb.obs.reset_index(),
                    pd.json_normalize(self.img_emb.obs["spots_info"]).reset_index(drop=True),
                ),
                axis=1,
            ).set_index("index")

    def concat_all_emb(self):
        """Concatenate all the embeddings along the columns."""
        self.emb = ad.concat([self.img_emb, self.gene_emb], axis=1, join="outer", uns_merge="first")
        emb_names = ["image_emb_{}".format(i) for i in self.img_emb.var_names]
        emb_names.extend(["gene_emb_{}".format(i) for i in self.gene_emb.var_names])
        if self.ef_emb is not None:
            self.emb = ad.concat([self.emb, self.ef_emb], axis=1, join="outer", uns_merge="first")
            emb_names.extend(["ef_emb_{}".format(i) for i in self.ef_emb.var_names])

        different_cols = self.gene_emb.obs.columns.difference(self.img_emb.obs.columns)
        self.emb.obs = self.img_emb.obs.join(self.gene_emb.obs[different_cols], how="outer")

        assert len(emb_names) == self.emb.shape[1], "Problem with var names shape"
        self.emb.var["index"] = emb_names
        self.emb.var = self.emb.var.set_index("index")

    def addition_all_emb(self):
        """Add all the embeddings."""
        assert self.img_emb.shape[1] == self.gene_emb.shape[1], "Image emb and gene emb must have the same number of columns"
        self.emb = ad.AnnData(self.img_emb.X + self.gene_emb.X)
        if self.ef_emb is not None:
            assert (
                self.img_emb.shape[1] == self.ef_emb.shape[1]
            ), "Image emb, gene emb and ef emb must have the same number of columns"
            self.emb.X = self.emb.X + self.ef_emb.X

        different_cols = self.gene_emb.obs.columns.difference(self.img_emb.obs.columns)
        self.emb.obs = self.img_emb.obs.join(self.gene_emb.obs[different_cols], how="outer")

    def compute_embeddings(self, only_labeled_spots=False, only_cancer=False):
        """Compute merged embeddings. It fill emb anndata.

        Args:
            only_labeled_spots (bool, optional): If you want only data with a label entry in .obs. Defaults to False.
            only_cancer (bool, optional): If you want only data with a label entry in .obs that is "invasive cancer" or "cancer in situ". Defaults to False.
        """
        self.prepare_data(only_labeled_spots=only_labeled_spots, only_cancer=only_cancer)

        if self.method == "concat":
            self.concat_all_emb()
        elif self.method == "cca":
            pass
        elif self.method == "addition":
            self.addition_all_emb()
        else:
            raise Exception("Method not implemented")
