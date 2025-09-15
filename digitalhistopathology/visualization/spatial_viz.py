#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import os
import numpy as np
import pandas as pd
import seaborn as sns
from PIL import Image

import openslide


import matplotlib.pyplot as plt

class SpatialViz():

    def __init__(self, emb, saving_plots=False, result_saving_folder=None):
        """Class to visualize the spatial embedding of the spots.

        Args:
            emb (Emb): Embedding object containing the data to visualize.
            saving_plots (bool, optional): If True, save the plots. Defaults to False.
            result_saving_folder (str, optional): Folder to save the plots. Defaults to None.
        """
        self.emb = emb
        self.saving_plots = saving_plots
        self.result_saving_folder = result_saving_folder

    def compare_two_labels_plot(
        self,
        sample_name="A1",
        label_obs_name1="label",
        label_obs_name2="predicted_label",
        path_origin_col="path_origin",
        x_col="start_width_origin",
        y_col="start_height_origin",
        palette_1=None,
        add_raw_image=True,
        ax_list=None,
        legend_pos_next=False,
    ):
        """Plot two different labels of spots on the image to be visually compared.

        Args:
            sample_name (str, optional): Name of the sample to plot. One of the string in image_embedding.emb.name_origin. Defaults to A1.
            label_obs_name1 (str, optional): First emb.obs column to compare. Defaults to "label".
            label_obs_name2 (str, optional): Second emb.obs column to compare. Defaults to "predicted_label".
            path_origin_col (str, optional): Column in emb.obs that contains the path to origin image. Defaults to "path_origin".
            x_col (str, optional): Column in emb.obs that contains x spot location on the origin image. Defaults to "start_width_origin".
            y_col (str, optional): Column in emb.obs that contains y spot location on the origin image. Defaults to "start_height_origin".
            add_raw_image (bool, optional): Add the raw image in the first ax. Defaults to True.
            ax_list (list, optional): List of axes to plot on it. Defaults to None.
            legend_pos_next (bool, optional): Legend next to the plot upper right instead that under the plot. Defaults to False.
        """
        data = self.emb[~self.emb.obs[label_obs_name1].isna()]
        data = data[~data.obs[label_obs_name2].isna()]
        data = data[data.obs.name_origin == sample_name]

        # nc_predicted = len([e for e in data.obs[label_obs_name2].unique() if isinstance(e, str) and e != "undetermined"])

        print(data)
        if palette_1 is None:
            palette_1 = sns.color_palette("tab10")
        # if palette_2 is None:
        #     palette_2 = sns.color_palette("Accent")

        true_labels_list = data.obs[label_obs_name1].unique()
        true_labels_list = [label for label in true_labels_list if label != "undetermined"]

        palette_1 = {label: palette_1[label] for label in true_labels_list}

        print(true_labels_list)
        ## Create the matching palette for the predicted labels
        df_palette = pd.DataFrame(index = true_labels_list)
        print(data.obs[label_obs_name2].unique())

        for i, lab in enumerate(data.obs[label_obs_name2].unique()):
            subset_label = data.obs[data.obs[label_obs_name2] == lab]
            subset_label = subset_label[subset_label[label_obs_name1] != "undetermined"]
            # df = (subset_label[label_obs_name1].value_counts()/len(subset_label)).loc[true_labels_list]
            df_count = (subset_label[label_obs_name1].value_counts())
            for true_label in true_labels_list:
                if true_label not in df_count.index:
                    df_count[true_label] = 0
            df = df_count.loc[true_labels_list]
            df_palette[lab] = df.values
        
        palette_2 = {}
        assigned_true_label = []
        assigned_predicted_label = []
        assigned_colors = set()

        for predicted_label in df_palette.columns:
            highest_label = df_palette[predicted_label].idxmax()
            highest_label_proportion = df_palette.loc[highest_label, predicted_label]
            if highest_label_proportion != df_palette.loc[highest_label].max():
                pass
            else:
                if palette_1[highest_label] not in assigned_colors:
                    palette_2[predicted_label] = palette_1[highest_label]
                    assigned_true_label.append(highest_label)
                    assigned_predicted_label.append(predicted_label)
                    assigned_colors.add(palette_1[highest_label])

        for predicted_label in df_palette.columns:
            if predicted_label not in assigned_predicted_label:
                remaining_labels = [true_label for true_label in df_palette.index if true_label not in assigned_true_label]
                max_label = df_palette.loc[remaining_labels][predicted_label].idxmax()
                if palette_1[max_label] not in assigned_colors:
                    palette_2[predicted_label] = palette_1[max_label]
                    assigned_true_label.append(max_label)
                    assigned_predicted_label.append(predicted_label)
                    assigned_colors.add(palette_1[max_label])
                else:
                    # Find the next available color
                    for label in remaining_labels:
                        if palette_1[label] not in assigned_colors:
                            palette_2[predicted_label] = palette_1[label]
                            assigned_true_label.append(label)
                            assigned_predicted_label.append(predicted_label)
                            assigned_colors.add(palette_1[label])
                            break
        print(palette_2)

        ncols = 3 if add_raw_image else 2
        if ax_list is None:
            fig, ax = plt.subplots(ncols=ncols, figsize=(5 * ncols, 5))
        else:
            ax = ax_list
        for i, label in enumerate([label_obs_name1, label_obs_name2]):
            palette = palette_1 if i == 0 else palette_2
            if add_raw_image:
                i = i + 1
            sns.scatterplot(
                data=data.obs, x=x_col, y=y_col, hue=label, s=20, alpha=0.9, ax=ax[i], palette=palette, hue_order=palette.keys()
            )
            bbox_to_anchor_1 = 1.45
            if label == "predicted_label":
                markers = [plt.Line2D([0, 0], [0, 0], color=color, marker="o", linestyle="") for color in palette_2.values()]
                ax[i].legend(markers, palette_2.keys(), numpoints=1)
                bbox_to_anchor_1 = 1.35
            if legend_pos_next:
                sns.move_legend(
                    ax[i],
                    "upper right",
                    bbox_to_anchor=(bbox_to_anchor_1, 1),
                    ncols=1,
                    title=label.replace("_", " ").capitalize(),
                )
            else:
                ncols == 3 if label == "predicted_label" else 2
                sns.move_legend(
                    ax[i], "lower center", bbox_to_anchor=(0.5, -0.4), ncols=ncols, title=label.replace("_", " ").capitalize()
                )
            plt.gca().invert_yaxis()
            image = np.array(Image.open(data.obs[path_origin_col].unique()[0]).convert("RGB"))
            ax[i].imshow(image)
            ax[i].set_title(sample_name + " " + label.replace("_", " "))
            ax[i].axis("off")
        if add_raw_image:
            ax[0].imshow(image)
            ax[0].set_title("{} raw image".format(sample_name))
            ax[0].axis("off")
        fig.show()
        if self.saving_plots:
            print("Saving plot in {}".format(self.result_saving_folder))
            plt.savefig(
                os.path.join(
                    self.result_saving_folder,
                    "compare_two_labels_plot_{}.pdf".format(sample_name),
                ),
                bbox_inches="tight",
            )
            plt.close()
        else:
            if ax_list is None:
                if "name_origin" in data.obs.columns:
                    fig.suptitle("Comparaison between two labelings for sample {}".format(data.obs["name_origin"].unique()[0]))
                fig.show()
            else:
                return ax_list

    def plot_spot_location_with_color_on_origin_image(
        self,
        path_origin_col="path_origin",
        x_col="start_width_origin",
        y_col="start_height_origin",
        color="label",
        list_name_origin=None,
        palette=None,
        alpha=0.8,
        s=20,
    ):
        """Plot image of origin with a scatter plot with specified positions and with a specified color.

        Args:
            path_origin_col (str, optional): Column name of emb.obs where the paths to origin image are. Defaults to "path_origin".
            x_col (str, optional): Column name of emb.obs where the position on the x axis are. Defaults to "start_width_origin".
            y_col (str, optional): Column name of emb.obs where the position on the y axis are. Defaults to "start_height_origin".
            color (str, optional): Column name of emb.obs for color. Defaults to "label".
            list_name_origin (list, optional): List of the name_origin of the images to plot. If None, all the images are plotted. Defaults to None.
            palette (dict, optional): Dictionnary that maps each color_obs_name to a color. Defaults to None.
            alpha (float, optional): Intensity of the scatter points on the image, between 0 and 1. Defaults to 0.8.
        """
        groups = self.emb.obs.groupby(by=path_origin_col)
        labels = list(self.emb.obs[color].unique())
        if list_name_origin is None:
            list_name_origin = list(self.emb.obs["name_origin"].unique())
        if np.nan in labels:
            labels.remove(np.nan)
        for name, group in groups:
            group = group[~group[color].isna()]
            if len(group) > 0 and group["name_origin"].values[0] in list_name_origin:
                sns.relplot(
                    data=group,
                    x=x_col,
                    y=y_col,
                    hue=color,
                    hue_order=labels,
                    s=s,
                    alpha=alpha,
                    palette=palette,
                    edgecolor=None
                )
                plt.gca().invert_yaxis()
                plt.imshow(np.array(Image.open(name).convert("RGB")))
                if "name_origin" in group.columns:
                    plt.title(group["name_origin"].unique()[0])
                # plt.show()

    def plot_unsupervised_clustering_no_labels(self,
                                            slide_name,
                                            slide_col_id='tumor',
                                            subsample_level=2,
                                            saving_folder=None,
                                            add_to_filename="",
                                            extension='pdf',
                                            palette=None,
                                            s=4, 
                                            path_origin_col="path_origin",
                                            x_col="start_width_origin",
                                            y_col="start_height_origin",):

        if saving_folder is not None:
            if not os.path.exists(saving_folder):
                os.makedirs(saving_folder)

        slide_emb = SpatialViz()
        slide_emb.emb = self.emb[self.emb.obs[slide_col_id] == slide_name]

        path_to_slide = slide_emb.emb.obs[path_origin_col].unique()[0]

        ## Load the slide
        slide = openslide.OpenSlide(path_to_slide)

        level = subsample_level
        downsampled_width, downsampled_height = slide.level_dimensions[level]

        # Read the region at the downsampled level
        downsampled_region = slide.read_region((0, 0), level, (downsampled_width, downsampled_height))

        # Convert to numpy array
        image = np.array(downsampled_region)

        slide_emb.emb.obs[f'start_width_origin_to_plot'] = slide_emb.emb.obs[x_col].apply(lambda x: x/ float(slide.properties[f'openslide.level[{level}].downsample'])) 
        slide_emb.emb.obs[f'start_height_origin_to_plot'] = slide_emb.emb.obs[y_col].apply(lambda x: x/ float(slide.properties[f'openslide.level[{level}].downsample']))

        # Display the slide with the predicted labels

        proportion_size = image.shape[0] / image.shape[1]

        plt.figure(figsize=(20, 10 * proportion_size))

        plt.subplot(1, 2, 1)
        plt.imshow(image)
        plt.title("Slide image", weight='bold')

        plt.subplot(1, 2, 2)
        plt.imshow(image)

        sns.scatterplot(data=slide_emb.emb.obs, 
                        x=f'start_width_origin_to_plot', 
                        y=f'start_height_origin_to_plot', 
                        s=s,
                        hue=slide_emb.emb.obs['predicted_label'],
                        palette=palette,
                        edgecolor=None)
        plt.ylabel('')
        plt.xlabel('')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.title(f"Slide image with predicted labels \n (20'000 random patches)", weight='bold')

        plt.tight_layout()

        if saving_folder is not None:
            plt.savefig(os.path.join(saving_folder,  f"unsupervised_clustering_no_labels_{slide_name}{add_to_filename}.{extension}"), bbox_inches='tight', dpi=300)
