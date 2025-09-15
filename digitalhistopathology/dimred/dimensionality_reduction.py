#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import warnings
import os

import matplotlib.pyplot as plt
from digitalhistopathology.helpers import corr2_coeff, plot_bokeh
from digitalhistopathology.utils import kl_divergence
from sklearn.manifold import TSNE
from umap import UMAP
import plotly.express as px

class DimRed:

    def __init__(self, emb=None, saving_plots=False, result_saving_folder=None):
        self.emb = emb
        self.svd = {"denoised": None, "U_df": None, "S": None, "V_df": None}
        self.saving_plots = saving_plots
        self.result_saving_folder = result_saving_folder

    def compute_svd(self, denoised=False, center=False, scale=False):
        """Compute svd from the embeddings anndata and store U, S, Vt into the svd variable.

        Args:
            denoised (bool, optional): If True it will use the denoised_emb layer of emb, otherwise it will use emb.X. Defaults to False.
        """
        if center or scale:
            scaler = StandardScaler(with_mean=center, with_std=scale)
        if denoised:
            if "denoised_emb" not in list(self.emb.layers.keys()):
                self.compute_denoised_embeddings(center=center, scale=scale)
            U, S, Vt = np.linalg.svd(self.emb.layers["denoised_emb"], full_matrices=False)
        else:
            if center or scale:
                U, S, Vt = np.linalg.svd(scaler.fit_transform(self.emb.X), full_matrices=False)
            else:
                U, S, Vt = np.linalg.svd(self.emb.X, full_matrices=False)
        self.svd["denoised"] = denoised
        self.svd["U_df"] = pd.DataFrame(
            U,
            columns=["u{}".format(i) for i in range(1, U.shape[1] + 1)],
            index=self.emb.obs_names,
        )
        self.svd["S"] = S
        self.svd["V_df"] = pd.DataFrame(Vt.T, columns=["v{}".format(i) for i in range(1, Vt.T.shape[1] + 1)])
        if self.emb.shape[0] < self.emb.shape[1]:
            warnings.warn(
                "Be careful your embedding matrix has more columns than rows. So, U matrix does not have the same size as your embedding matrix. The other functions that use svd could have some problems..."
            )

    def compute_tsne(self, n_components=2, perplexity=30):
        self.emb.obsm["tsne"] = TSNE(
            n_components=n_components,
            learning_rate="auto",
            init="random",
            perplexity=perplexity,
        ).fit_transform(self.emb.X)[:, 0:2]

    def compute_umap(self, n_components=2, n_neighbors=15, min_dist=0.1, metric='euclidean'):

        umap_reducer = UMAP(n_components=n_components, n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, random_state=42)
        self.emb.obsm["umap"] = umap_reducer.fit_transform(self.emb.X)[:, 0:2]

    def get_svd_components(self, components_nb, denoised=False, center=False, scale=False):
        """Return the svd matrix corresponding to components_nb. si * Ui @ Vi.

        Args:
            components_nb (int): The number to choose the singular vectors and the singular values indexes. >=1.
            denoised (bool, optional): If True it will use the denoised svd matrices. Defaults to False.

        Returns:
            np.array: Matrix corresponding to components_nb of the same size as emb.X.
        """
        if center or scale:
            scaler = StandardScaler(with_mean=center, with_std=scale)

        if denoised:
            U, S, Vt = np.linalg.svd(self.emb.layers["denoised_emb"], full_matrices=False)
        else:
            if center or scale:
                U, S, Vt = np.linalg.svd(scaler.fit_transform(self.emb.X), full_matrices=False)
            else:
                U, S, Vt = np.linalg.svd(self.emb.X, full_matrices=False)
        u1 = U[:, components_nb - 1]
        u1 = np.reshape(u1, (u1.shape[0], 1))
        v1 = Vt[components_nb - 1, :]
        v1 = np.reshape(v1, (1, v1.shape[0]))
        s1 = S[components_nb - 1]
        return s1 * (u1 @ v1)
    
    def compute_denoised_embeddings(self, center=False, scale=False):
        """Compute svd on the embeddings from which component 1 is subtracted. The matrix is stored into the emb.layers["denoised_emb"]."""
        svd1 = self.get_svd_components(1, denoised=False, center=center, scale=scale)   
        assert svd1.shape == self.emb.X.shape
        self.emb.layers["denoised_emb"] = self.emb.X - svd1

    def get_explained_variance_ratio_list(self, denoised=False):
        """Get the explained variance ratio list from the S matrix that was computed during svd.

        Returns:
            list: List of the explained variance ratio for each svd component.
        """
        if self.svd["S"] is None:
            print(f"SVD was not computed. Computing it now with denoised={denoised}...")
            self.compute_svd(denoised=denoised)

        eigenvalues = self.svd["S"] ** 2
        explained_variance_ratio = eigenvalues / np.sum(eigenvalues)
        return explained_variance_ratio
    
    def get_shannon_entropy(self, n_comp=None, denoised=None, rescale=False, pct=None):
        """Get the Shannon entropy from the svd. The higher it is, the more spread is the variance across components.

        Args:
            n_comp (int, optional): Number of components on which to calculate the shannon entropy. If None, it takes the max number of components. Default to None.
            denoised (bool, optional): If True, the denoised svd is used. Defaults to None.
            rescale (bool, optional): If True, the explained variance ratio is rescaled to sum up to 1. Useful when computing Shannon entropy not on all components. Defaults to False.
            pct (float, optional): Percentage of variance to keep. If None, it takes the max number of components. Default to None.
        Returns:
            float: Shannon entropy
        """

        if denoised is not None:
            if self.svd["denoised"] != denoised:
                print("Computing svd in shannon entropy function with denoised={}".format(denoised))
                self.compute_svd(denoised=denoised)
        if (n_comp is None) and (pct is None):
            n_comp = min(self.emb.shape[0], self.emb.shape[1])

        if pct is not None:
            explained_variance_ratio = self.get_explained_variance_ratio_list()
            n_comp = np.argmax(np.cumsum(explained_variance_ratio) > pct)

        explained_variance_ratio = self.get_explained_variance_ratio_list()[0:n_comp]

        if rescale:
            variance_total = np.sum(explained_variance_ratio)
            explained_variance_ratio = explained_variance_ratio / variance_total

        return (
            -1 * np.sum(explained_variance_ratio * np.log10(explained_variance_ratio)) / np.log10(len(explained_variance_ratio))
        )

    def get_kl_divergence(self, n_comp=None, denoised=None, pct=None):
        """Get the Kullback-Leibler divergence from the svd. The higher it is, the less spread is the variance across components.

        Args:
            n_comp (int, optional): Number of components on which to calculate the shannon entropy. If None, it takes the max number of components. Default to None.
            denoised (bool, optional): If True, the denoised svd is used. Defaults to None.
            pct (float, optional): Percentage of variance to keep. If None, it takes the max number of components. Default to None.

        Returns:
            float: Kullback-Leibler divergence
        """

        if denoised is not None:
            if self.svd["denoised"] != denoised:
                print("Computing svd in shannon entropy function with denoised={}".format(denoised))
                self.compute_svd(denoised=denoised)

        if (n_comp is None) and (pct is None):
            n_comp = min(self.emb.shape[0], self.emb.shape[1])

        if pct is not None:
            explained_variance_ratio = self.get_explained_variance_ratio_list()
            n_comp = np.argmax(np.cumsum(explained_variance_ratio) > pct)

        explained_variance_ratio = self.get_explained_variance_ratio_list()[0:n_comp]

        # Uniform distribution
        q = np.ones_like(explained_variance_ratio) / len(explained_variance_ratio)

        # Calculate KL divergence
        kl_div = kl_divergence(explained_variance_ratio, q)

        return kl_div

    def get_shannon_entropy_per_group(
        self,
        group="batch",
        denoised=False,
        n_comp=None,
        pct=None,
        rescale=False,
        cancer_patches=False,
    ):
        """Get the Shannon entropy for each group of data.

        Args:
            group (str, optional): The columns of emb.obs with which the data are grouped. Defaults to "batch".
            denoised (bool, optional): If True, the denoised svd is used. Defaults to False.
            n_comp (int, optional): Number of components on which to calculate the shannon entropy. If None, it takes the max number of components. Default to None.

        Returns:
            dict: Keys are the name of the groups and values are the associated Shannon entropy.
        """
        if cancer_patches:
            embeddings = self.emb.obs[
                self.emb.obs["knn_predicted_label"] == "invasive cancer"
            ]
        else:
            embeddings = self.emb.obs
        groups = embeddings.groupby(by=group).groups
        batch_list = list(groups.keys())

        entropy_dict = {}
        for batch in batch_list:
            subset = DimRed()
            subset.emb = self.emb[groups[batch], :]
            subset.compute_svd(denoised=denoised)
            entropy_dict[batch] = subset.get_shannon_entropy(n_comp=n_comp, denoised=denoised, pct=pct, rescale=rescale)
        return entropy_dict

    def get_kl_divergence_per_group(
        self,
        group="batch",
        denoised=False,
        n_comp=None,
        pct=None,
        rescale=False,
        cancer_patches=False,
    ):
        """Get the Kullback-Leibler divergence for each group of data.

        Args:
            group (str, optional): The columns of emb.obs with which the data are grouped. Defaults to "batch".
            denoised (bool, optional): If True, the denoised svd is used. Defaults to False.
            n_comp (int, optional): Number of components on which to calculate the shannon entropy. If None, it takes the max number of components. Default to None.

        Returns:
            dict: Keys are the name of the groups and values are the associated Kullback-Leibler divergence.
        """
        if cancer_patches:
            embeddings = self.emb.obs[
                self.emb.obs["knn_predicted_label"] == "invasive cancer"
            ]
        else:
            embeddings = self.emb.obs
        groups = embeddings.groupby(by=group).groups
        batch_list = list(groups.keys())

        kl_dict = {}
        for batch in batch_list:
            subset = DimRed()
            subset.emb = self.emb[groups[batch], :]
            subset.compute_svd(denoised=denoised)
            kl_dict[batch] = subset.get_kl_divergence(n_comp=n_comp, denoised=denoised, pct=pct)
        return kl_dict
    
    def boxplot_value_per_group(
        self,
        dict_value_per_group,
        value="Shannon entropy",
        group_name="Patient",
        xlabel="SVD",
        invert_x_axis=True,
    ):
        df = pd.DataFrame.from_records(dict_value_per_group)
        df.loc["mean"] = df.mean()
        df.boxplot(grid=False)
        plt.xlabel(xlabel)
        plt.ylabel(value)
        plt.title("Boxplot of the {} per {} for each {}".format(value, group_name.lower(), xlabel.lower()))

        colors = plt.cm.get_cmap("tab10", len(df.index))

        for idx, index in enumerate(df.index):
            plt.scatter(
                np.arange(1, len(list(df.columns)) + 1),
                df.loc[index, :],
                color=colors(idx),
                label=index,
                alpha=0.6,
            )

        if invert_x_axis:
            plt.gca().invert_xaxis()  # Invert the x-axis

        # Rotate x-axis labels for better visibility
        plt.xticks(rotation=45, ha="right")

        plt.legend(title=group_name, bbox_to_anchor=(1.05, 1), loc="upper left")

        if self.saving_plots:
            plt.savefig(
                os.path.join(
                    self.result_saving_folder,
                    "boxplot_{}_per_{}_for_each_{}.pdf".format(value, group_name, xlabel).lower().replace(" ", "_"),
                ),
                bbox_inches="tight",
            )
            plt.close()
        else:
            plt.show()

    def scree_plot(self, components_number_to_display=50, matrix_name="raw", ylim=None):
        """Plot the scree plot from the svd.

        Args:
            components_number_to_display (int, optional): Number of svd components to display. Defaults to 50.
            matrix_name (str, optional): Name of the matrix to be displayed in the plot title. Defaults to "raw".
        """
        entropy = self.get_shannon_entropy(denoised=self.svd["denoised"])
        print("Shannon entropy: {}".format(entropy))

        explained_variance_ratio = self.get_explained_variance_ratio_list()

        plt.subplots(figsize=(10, 4))
        plt.bar(
            np.arange(1, components_number_to_display + 1),
            explained_variance_ratio[0:components_number_to_display],
        )
        plt.xticks(
            np.arange(1, components_number_to_display + 1, 2),
            np.arange(1, components_number_to_display + 1, 2),
        )
        plt.xlabel("SVD components")
        plt.ylabel("Fraction of explained variance")
        plt.title(
            "Fraction of explained variance captured by {} first SVD components from the {} matrix".format(
                components_number_to_display, matrix_name
            ), weight='bold'
        )
        if ylim is not None:
            plt.ylim(ylim[0], ylim[1])

        if self.saving_plots:
            plt.savefig(
                os.path.join(
                    self.result_saving_folder,
                    "svd_scree_plot_" + matrix_name.replace(" ", "_") + ".pdf",
                ),
                bbox_inches="tight",
            )
            plt.close()
        else:
            plt.show()

    def cumulative_variance_plot(self, matrix_name="raw"):
        """Plot the cumulative variance from the svd.

        Args:
            matrix_name (str, optional): Name of the matrix to be displayed in the plot title. Defaults to "raw".
        """
        explained_variance_ratio = self.get_explained_variance_ratio_list()

        cumulative_explained_variance_ratio = np.cumsum(explained_variance_ratio)
        print(
            "90 % of the variance is captured by {} components".format(
                np.argwhere(cumulative_explained_variance_ratio > 0.9)[0][0] + 1
            )
        )

        plt.figure(figsize=(8, 5))
        plt.title("Cumulative variance explained from the SVD on the {} matrix".format(matrix_name))
        plt.plot(
            np.arange(1, explained_variance_ratio.shape[0] + 1),
            cumulative_explained_variance_ratio,
        )
        plt.xlabel("SVD Components")
        plt.ylabel("Variance explained")
        plt.grid(True)
        if self.saving_plots:
            plt.savefig(
                os.path.join(
                    self.result_saving_folder,
                    "svd_cumulative_variance_plot_" + matrix_name.replace(" ", "_") + ".pdf",
                ),
                bbox_inches="tight",
            )
            plt.close()
        else:
            plt.show()

    def compute_corr_matrix_u(self):
        """Compute the correlation matrix between the left singular vectors from the svd (ui) and the embedding columns (features).

        Returns:
            pd.DataFrame: Correlation matrix the left singular vectors and the embedding columns.
        """
        correlation_matrix = corr2_coeff(self.emb.X.T, self.svd["U_df"].values.T)
        return pd.DataFrame(correlation_matrix, columns=self.svd["U_df"].columns)

    def compute_corr_ui(self, i):
        """Compute the correlation matrix between the left singular vectors from the svd (ui) and the embedding rows (spots).

        Args:
            i (int): Svd component number to select the left singular vector ui, i>=1

        Returns:
            pd.DataFrame: Correlation between the ith left singular vector and the embedding rows.
        """
        ui = self.svd["U_df"].iloc[:, i - 1]
        correlation_matrix = corr2_coeff(self.emb.X.T, ui.values.reshape((-1, 1)).T)
        return pd.DataFrame(correlation_matrix, columns=[ui.name])

    def compute_corr_matrix_v(self):
        """Compute the correlation matrix between the right singular vectors from the svd (vi) and the embedding rows (spots).

        Returns:
            pd.DataFrame: Correlation matrix between the right singular vectors and the embedding rows.
        """
        correlation_matrix = corr2_coeff(self.emb.X, self.svd["V_df"].values.T)
        return pd.DataFrame(correlation_matrix, columns=self.svd["V_df"].columns, index=self.emb.obs_names)

    def compute_corr_vi(self, i):
        """Compute the correlation matrix between the right singular vector from the svd (vi) and the embedding rows (spots).

        Args:
            i (int): Svd component number to select the right singular vector vi, i>=1

        Returns:
            pd.DataFrame: Correlation between the ith right singular vector and the embedding rows.
        """
        vi = self.svd["V_df"].iloc[:, i - 1]
        correlation_matrix = corr2_coeff(self.emb.X, vi.values.reshape((-1, 1)).T)
        return pd.DataFrame(correlation_matrix, columns=[vi.name], index=self.emb.obs_names)

    @staticmethod
    def extract_top_features(df, col_ind=0, top=20):
        """Return the row index of the top positive and negative contributor to the col_ind feature.

        Args:
            df (pd.DataFrame): Dataframe with features as columns and samples as rows.
            col_ind (int, optional): Column index for the feature. Defaults to 0.
            top (int, optional): Number of top features to return. Defaults to 20.

        Returns:
            list: List of the top positive indexes.
            list: List of the top negative indexes.
        """
        vector = df.iloc[:, col_ind]
        sorted_vector = vector.sort_values()
        sorted_index = list(sorted_vector.index)
        # [::-1] to have the top positive in descending order
        return sorted_index[0:top], sorted_index[-top:][::-1]
    
    def get_ui_p_val_dict(self, cluster_name, linear_mixed=False, max_number_comp=None, obs_col_name="predicted_label"):
        """
        Compute many linear model tests and gather all resulting pvalues into one dictionnary.
        Args:
            cluster_name (str): Cluster name from predicted_label column to test.
            linear_mixed (bool, optional): If linear mixed model is used. If False, linear model is used. Defaults to False.
            max_number_comp (int, optional): Number of first svd components to test. If None, all the available componenents are tested. Defaults to None.
            obs_col_name: (str, optional): self.emb.obs column name that contain cluster label. Defaults to "predicted_label".

        Returns:
            dict: Dictionnary that map each svd component to the p-value corresponding to linear model testing associated with cluster_name.
        """
        p_vals_dict = dict()
        if max_number_comp is None:
            max_number_comp = self.svd["U_df"].shape[1]
        for i in range(max_number_comp):
            if linear_mixed:
                p_vals_dict["u{}".format(i + 1)] = self.lmer_per_ui(
                    comp=i + 1, cluster_name=cluster_name, obs_col_name=obs_col_name, verbose=False
                )
            else:
                p_vals_dict["u{}".format(i + 1)] = self.lm_per_ui(
                    comp=i + 1, cluster_name=cluster_name, obs_col_name=obs_col_name, verbose=False
                )
        return p_vals_dict

    def compute_umap_plot(
        self,
        color_obs_column="name_origin",
        matrix_name="embeddings",
        title_font_size=20,
        recompute_umap=True,
        layer=None,
        palette=None,
        bokeh_plot=False,
    ):
        """Compute UMAP from emb or another emb layer and plot the two components. Save UMAP results in emb.obsm["umap"].

        Args:
            color_obs_column (str, optional): Column name in emb.obs for color. Defaults to "name_origin".
            matrix_name (str, optional): Name of the matrix to be displayed in the plot title. Defaults to "embeddings".
            title_font_size (int, optional): Font size of the title of the plot. Defaults to 20.
            recompute_umap (bool, optional): If the UMAP must be computed. If False, it will take the UMAP results in emb.obsm["umap"]. Defaults to True.
            layer (str, optional): Layer name in emb.layers to use to compute the UMAP. If None, the raw emb is taken. Defaults to None.
            palette (dict, optional): Dictionnary that maps each color_obs_name to a color. Defaults None.
            bokeh_plot (bool, optional): If we plot bokeh scatter plot with patches on hover data instead of regular plotly plots. Defaults to False.
        """
        if recompute_umap:
            reducer = UMAP()
            if layer is None:
                matrix = self.emb.X.copy()
            elif layer in list(self.emb.layers.keys()):
                matrix = self.emb.layers[layer].copy()
            elif layer in list(self.emb.obsm.keys()):
                matrix = self.emb.obsm[layer].copy()
            self.emb.obsm["umap"] = reducer.fit_transform(matrix)

        hover_data = ["name"]
        umap_df = pd.DataFrame(self.emb.obsm["umap"], columns=["comp1", "comp2"])
        umap_df["name"] = list(self.emb.obs_names)
        if "label" in self.emb.obs.columns:
            hover_data.append("label")
            umap_df["label"] = list(self.emb.obs.label)

        if color_obs_column and color_obs_column not in umap_df.columns:
            umap_df[color_obs_column] = self.emb.obs[color_obs_column].values

        if bokeh_plot:
            hover_data.append(color_obs_column)
            umap_df["path"] = list(self.emb.obs["path"])
            saving_path = (
                None
                if not self.saving_plots
                else os.path.join(
                    self.result_saving_folder,
                    "compute_umap_bokeh_plot_" + matrix_name.replace(" ", "_") + ".html",
                )
            )
            plot_bokeh(
                umap_df,
                x="comp1",
                y="comp2",
                img_path="path",
                other_hover_data=hover_data,
                title="UMAP on the {}".format(matrix_name),
                color_col=color_obs_column,
                palette=palette,
                saving_path=saving_path,
            )
            return

        if color_obs_column:
            umap_df[color_obs_column] = self.emb.obs[color_obs_column].values
            fig = px.scatter(
                umap_df,
                x="comp1",
                y="comp2",
                hover_data=hover_data,
                color=color_obs_column,
                width=900,
                height=700,
                color_discrete_map=palette,
            )
        else:
            fig = px.scatter(
                umap_df,
                x="comp1",
                y="comp2",
                hover_data=hover_data,
                width=900,
                height=700,
                color_discrete_map=palette,
            )
            fig.update_traces(marker_size=7)

        fig.update_layout(
            title="UMAP on the {}".format(matrix_name),
            xaxis=dict(showgrid=False, ticks="inside", linecolor="#BCCCDC"),
            yaxis=dict(showgrid=False, ticks="inside", linecolor="#BCCCDC"),
            plot_bgcolor="#FFF",
            title_font_size=title_font_size,
        )

        if self.saving_plots:
            fig.write_image(
                os.path.join(
                    self.result_saving_folder,
                    "compute_umap_plot_" + matrix_name.replace(" ", "_") + "_colored_by_" + color_obs_column + ".pdf",
                )
            )
            plt.close()
        else:
            fig.show()
