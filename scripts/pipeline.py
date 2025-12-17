
#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import os
import argparse
import copy
import glob
import json
import time
import warnings
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
from huggingface_hub import logout

import sys
sys.path.append(os.getcwd())

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

from digitalhistopathology.clustering.clustering_utils import clustering_boxplot_per_patient
from digitalhistopathology.datasets.real_datasets import HER2Dataset, TNBCDataset
# from digitalhistopathology.datasets.real_datasets import VisiumHDdataset (preliminary analysis for Ovarian)
from digitalhistopathology.datasets.spatial_dataset import MixedImageDataset, SpatialDataset
from digitalhistopathology.embeddings.image_embedding import ImageEmbedding
from digitalhistopathology.helpers import NumpyEncoder
from digitalhistopathology.models import Model, load_model
from digitalhistopathology.monitor import Monitor


class Pipeline:

    def __init__(
        self,
        name: str,
        model: Model,
        results_folder="results",
        patches_folder=None,
        selected_invasive_cancer_file=None,
        dataset=None,
        svd_components_number_shannon_entropy=512,  # to be fixed to the smallest number of deep features across all the tested models
        gene_embedding_spot_normalization=False,
        stats_dge_within_patient="wilcox",
    ):
        """A pipeline was developed to test different deep learning encoders for digital histopathology images 
        to find the ones that reflect the most the molecular data

        Args:
            name (str): Name of the pipeline
            model (PretrainedModel): Pretrained model from which computes the image embedding.
            results_folder (str, optional): Where to save the pipeline results. Defaults to "../results/HER2/pipeline".
            patches_folder (str, optional): Path to patches folder. If None or path does not exists, patches are computed to "../results/HER2/pipeline/patches". Defaults to None.
            selected_invasive_cancer_file (str, optional): Csv file for selected invasive cancer. Defaults to None.
            dataset (SpatialDataset, optional): Dataset. Defaults to HER2Dataset().
            svd_components_number_shannon_entropy (int, optional): Number of svd component to take for the Shannon entropy computation. To be fixed to the lower embedding dimension of all the models you want to benchmark. Defaults to 512.
            n_clustes_to_test_kmeans (list/1D-array, optional): List of number of the number of clusters k to test during kmeans clustering. Defaults to np.arange(2, 21, 1).
            gene_embedding_spot_normalization (bool, optional): If CPM spot normalization is done in the gene embedding or not. Defaults to False.
            stats_dge_within_patient (str, optional): Statistic for FindMarkers or FindAllMarkers from Seurat R package. Defaults to "wilcox".
        """
        self.name = name
        self.model = model

        self.stats_dge_within_patient = stats_dge_within_patient
        self.results_folder = results_folder
        self.__run_compute_patches = False
        self.patches_folder = patches_folder
        self.__run_select_invasive_cancer_patches = False
        self.selected_invasive_cancer_file = selected_invasive_cancer_file
        self.dataset: SpatialDataset = dataset
        self.saving_folders_and_files_init()

        self.dataset.patches_folder = self.patches_folder
        self.dataset.saving_emb_folder = self.model_results_folder

        self.svd_components_number_shannon_entropy = svd_components_number_shannon_entropy

        self.gene_embedding_spot_normalization = gene_embedding_spot_normalization

        self.image_embedding = None
        self.invasive_image_embedding = None
        self.invasive_gene_embedding = None

    def __create_folder(self, name, parent_folder):
        folder = os.path.join(parent_folder, name)
        if not os.path.exists(folder):
            os.mkdir(folder)
        return folder

    def __set_figure_folder(self, folder):
        if self.image_embedding is not None:
            self.image_embedding.result_saving_folder = folder
        if self.invasive_image_embedding is not None:
            self.invasive_image_embedding.result_saving_folder = folder
        if self.invasive_gene_embedding is not None:
            self.invasive_gene_embedding.result_saving_folder = folder

    def saving_folders_and_files_init(self):
        # Ensure the parent directory exists
        parent_dir = os.path.dirname(self.results_folder)
        print(self.results_folder)
        if not os.path.exists(parent_dir):
            os.makedirs(parent_dir)

        # Create the results folder
        if not os.path.exists(self.results_folder):
            print("Saving results to:", self.results_folder)
            os.makedirs(self.results_folder)

        if self.patches_folder is None:
            self.__run_compute_patches = True
            self.patches_folder = self.__create_folder(name="patches", parent_folder=self.results_folder)
            if len(glob.glob(os.path.join(self.patches_folder, "*.tiff"))) > 0:
                print("The folder {} is already filled with .tiff images.".format(self.patches_folder))
                raise Exception(
                    "Problem with patches folder {}: it is already filled, specify the patches_folder of the pipeline or rename/delete the folder to recompute patches".format(
                        self.patches_folder
                    )
                )

        self.model_results_folder = self.results_folder

        # TODO: Transformed the print into logs (logging) in pipeline.py ? For now, copying the sbatch logs to the model folder
        self.logs_folder = self.__create_folder(name="logs", parent_folder=self.model_results_folder)

        self.image_embedding_saving_path = os.path.join(self.model_results_folder, "image_embedding.h5ad")

        self.shannon_entropy_result_parent_folder = self.__create_folder(
            name="shannon_entropy", parent_folder=self.model_results_folder
        )

        self.shannon_entropy_results_folder = self.__create_folder(
            name=self.dataset.name, parent_folder=self.shannon_entropy_result_parent_folder
        )
        self.shannon_entropy_results_path = os.path.join(self.shannon_entropy_results_folder, "shannon_entropy.json")

        self.select_invasive_cancer_results_folder = self.__create_folder(
            name="select_invasive_cancer", parent_folder=self.model_results_folder
        )
        if self.selected_invasive_cancer_file is None:
            self.selected_invasive_cancer_file = os.path.join(
                self.select_invasive_cancer_results_folder, "selected_invasive_cancer.csv"
            )
            self.__run_select_invasive_cancer_patches = True
        self.mlp_classification_folder = self.__create_folder(
            name="mlp_classification", parent_folder=self.model_results_folder
        )
        self.mlp_classification_file = os.path.join(
            self.mlp_classification_folder, "confusion_matrix_avg.pdf"
        )

        self.enrichement_score_from_knn_results_path = os.path.join(
            self.select_invasive_cancer_results_folder, "enrichement_score_from_knn.json"
        )

        self.invasive_image_embedding_saving_path = os.path.join(self.model_results_folder, "invasive_image_embedding.h5ad")

        self.gene_embedding_saving_path = os.path.join(self.model_results_folder, "gene_embedding.h5ad")
        
        self.invasive_gene_embedding_saving_path = os.path.join(self.model_results_folder, "invasive_gene_embedding.h5ad")

        self.pipeline_settings_path = os.path.join(self.model_results_folder, "pipeline_settings.json")

    def run(self, knn_relabeling):

        print("START running pipeline\n")

        self.load_and_save_whole_images_embeddings()

        self.select_invasive_cancer_patches_pipeline(relabeling=knn_relabeling)

        self.shannon_entropy_pipeline()

        self.__set_figure_folder(folder=self.model_results_folder)

        logout()

        print("\nEND pipeline")

    def load_and_save_whole_images_embeddings(self, monitor_delay=None, benchmarking=False):
        """Load and and save whole image emebedding as image_embedding.h5ad. Save cpu.json, ram.json and gpu.json to compute_image_embedding_monitoring folder.

        Args:            monitor_delay (int, optional): Monitor delay to store cpu, ram and gpu. Defaults to 5.
            benchmarking (boolean, optional): True if computing the image embeddings from the benchmarking dataset. False if computing the image embeddings from self.dataset. Defaults to False.
        """
        if benchmarking:
            self.benchmark_dataset = MixedImageDataset(
                folder=self.benchmark_patches_folder, saving_emb_folder=self.model_results_folder, name="benchmark"
            )

        def select_folder():
            return self.benchmark_patches_folder if benchmarking else self.patches_folder

        def select_image_embedding_saving_path():
            return self.benchmark_image_embedding_saving_path if benchmarking else self.image_embedding_saving_path

        def select_dataset():
            print(f"Dataset: {self.dataset}", flush=True)
            return self.benchmark_dataset if benchmarking else self.dataset

        if not os.path.exists(select_image_embedding_saving_path()):
            print("\nStart saving {} embeddings for {}".format(self.model.name, select_folder()))

            patches_filenames = (
                self.dataset.patches_filenames
                if not benchmarking
                else self.benchmark_dataset.patches_filenames
            )
            patches_info_filename = os.path.join(select_folder(), "patches_info.pkl.gz")

            self.image_embedding = ImageEmbedding(
                patches_filenames=patches_filenames,
                patches_info_filename=patches_info_filename,
                pretrained_model=self.model,
                label_files=select_dataset().label_filenames,
                name=select_dataset().name + "_" + self.model.name,
                result_saving_folder=self.model_results_folder,
                saving_plots=True,
                palette=select_dataset().PALETTE,
            )
            
            print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)
            print(f"Label files: {self.image_embedding.label_files}", flush=True)

            if monitor_delay is not None:
                monitor = Monitor(delay=monitor_delay)
                monitor.start()
            
            print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), flush=True)

            start = time.time()
            self.image_embedding.compute_embeddings()
            end = time.time()

            print("Time to compute embeddings for {}: {}".format(self.model.name, end - start))

            if monitor_delay is not None:
                monitor.stop()
                monitor_saving_path = os.path.join(self.model_results_folder, "compute_image_embedding_monitoring")
                if not os.path.exists(monitor_saving_path):
                    os.mkdir(monitor_saving_path)
                monitor.save_results(monitor_saving_path)

            if not benchmarking:
                print("\nAdd label to image embedding", flush=True)
                self.image_embedding.add_label(dataset=self.dataset.name)
                
            if self.dataset.name == "TNBC":
                self.image_embedding.emb.obs['tumor'] = self.image_embedding.emb.obs['tumor'].apply(lambda x: x.split('_')[0])

            if self.dataset.name == "Ovarian":
                self.image_embedding.emb.obs['tumor'] = self.image_embedding.emb.obs['tumor'].apply(lambda x: x.replace('_', '-'))
                self.image_embedding.emb.obs['tumor'] = self.image_embedding.emb.obs['tumor'].apply(lambda x: x.split('HUG-')[1])
                self.image_embedding.emb.obs['tumor'] = self.image_embedding.emb.obs['tumor'].apply(lambda x: x.split('-')[0] if '-' in x else x)
            self.image_embedding.save_embeddings(saving_path=select_image_embedding_saving_path())
            print("Saving to {} ok".format(select_image_embedding_saving_path()))
        
        else:
            print("\nLoad image embedding from: {}".format(select_image_embedding_saving_path()))
            self.image_embedding = select_dataset().get_image_embeddings(
                self.model, filename=select_image_embedding_saving_path().split("/")[-1].split(".")[0]
            )
            self.image_embedding.result_saving_folder = self.model_results_folder
            self.image_embedding.saving_plots = True
            
            if "label" not in self.image_embedding.emb.obs.columns:
                print("Add label to image embedding", flush=True)
                self.image_embedding.add_label(dataset=self.dataset.name)
                

            if self.dataset.name == "Ovarian":
                self.image_embedding.emb.obs['tumor'] = self.image_embedding.emb.obs['tumor'].apply(lambda x: x.replace('_', '-'))
                self.image_embedding.emb.obs['tumor'] = self.image_embedding.emb.obs['tumor'].apply(lambda x: x.split('HUG-')[1])
                self.image_embedding.emb.obs['tumor'] = self.image_embedding.emb.obs['tumor'].apply(lambda x: x.split('-')[0] if '-' in x else x)
            # # Just for this round

            self.image_embedding.save_embeddings(saving_path=select_image_embedding_saving_path())
            
    def shannon_entropy_pipeline(self):
        """Compute the shannon entropy from denoised and non-denoised SVD on the whole image_embedding and separated by patient (shannon_entropy.json). Scree plots from the SVD on
        the whole matrix and boxplot of the shannon entropy across patient are saved. All results can be found in the shannon_entropy/HER2 folder.
        """
        self.__set_figure_folder(folder=self.shannon_entropy_results_folder)
        print("\nStart Shannon entropy pipeline")
        if not os.path.exists(self.shannon_entropy_results_path):
            print("Start Shannon entropy analysis")
            shannon_entropies_dict = dict()
            shannon_entropies_dict["svd_components_number"] = self.svd_components_number_shannon_entropy
            shannon_entropies_dict["per_patient"] = dict()
            shannon_entropies_dict["whole_matrix"] = dict()

            shannon_entropies_dict["whole_matrix"]["raw_all_patches"] = self.image_embedding.get_shannon_entropy(
                n_comp=None, denoised=False
            )
            self.image_embedding.scree_plot(components_number_to_display=15, matrix_name="{} raw_all_patche".format(self.model.name))

            shannon_entropies_dict["whole_matrix"]["raw_cancer_patches"] = self.invasive_image_embedding.get_shannon_entropy(
                n_comp=None, denoised=False
            )
            self.image_embedding.scree_plot(
                components_number_to_display=15, matrix_name="{} raw_cancer_patches".format(self.model.name)
            )
            shannon_entropies_dict["per_patient"]["raw_cancer_patches"] = (
                self.image_embedding.get_shannon_entropy_per_group(
                    group="tumor",
                    denoised=False,
                    n_comp=None,
                    cancer_patches=True,
                )
            )
            shannon_entropies_dict["per_patient"]["raw_all_patches"] = (
                self.image_embedding.get_shannon_entropy_per_group(
                    group="tumor",
                    denoised=False,
                    n_comp=None,
                )
            )
            self.image_embedding.boxplot_value_per_group(dict_value_per_group=shannon_entropies_dict["per_patient"])

            self.image_embedding.emb.uns["shannon_entropy_results"] = shannon_entropies_dict
            self.image_embedding.save_embeddings(saving_path=self.image_embedding_saving_path)

            with open(self.shannon_entropy_results_path, "w") as fp:
                json.dump(shannon_entropies_dict, fp)
            print("Save Shannon entropy results to {}".format(self.shannon_entropy_results_path))
        elif "shannon_entropy_results" not in list(self.image_embedding.emb.uns.keys()):
            print("Load Shannon entropy results")
            with open(self.shannon_entropy_results_path) as f:
                self.image_embedding.emb.uns["shannon_entropy_results"] = json.load(f)
        else:
            print("Shannon entropy results already stored in the image embedding object.")

    def enrichement_score_from_knn(self):
        """Compute enrichement scores by compting the fraction of already labeled data that was relabeled in the same label for all label and save results
        in enrichement_score_from_knn.json.
        """
        results = self.image_embedding.get_enrichement_score(groupby="knn_predicted_label")

        with open(self.enrichement_score_from_knn_results_path, "w") as fp:
            json.dump(results, fp)
            print("Save enrichement score results to {}".format(self.enrichement_score_from_knn_results_path))

        self.image_embedding.emb.uns["enrichement_score_from_knn_results"] = results

    def select_invasive_cancer_patches_pipeline(self, relabeling=False):
        """Select invasive cancer patches by computing KNN on image_embeddding by using label column as training data. At the end, we selected only invasive_cancer label to
        go more in depth into tumor heterogeneity. We also computed enrichement scores by compting the fraction of already labeled data that was relabeled in
        the same label for all label (enrichement_score_from_knn.json) which is computed only if relabeling=True. The selected invasive cancer patches are stored in a csv
        (selected_invasive_cancer.csv). Boxplot of the F1 and accuracy across patient to choose the optimal k of knn and barplots with patient and label fraction in each knn
        cluster are saved. All results can be found in the select_invasive_cancer folder.

        Args:
            relabeling (bool, optional): If the already labeled patches are relabeled by KNN. If False, no enrichement score is computed. Defaults to False.
        """
        self.__set_figure_folder(folder=self.select_invasive_cancer_results_folder)
        print("\nStart select invasive cancer patches pipeline")
        if not os.path.exists(self.selected_invasive_cancer_file):
            print("Start select invasive cancer patches analysis")
            optimal_k = self.image_embedding.get_optimal_k_for_knn_per_patient(
                svd_comp=None,
                weights="uniform",
                metric="minkowski",
                k_list=[3, 5, 7, 9, 11, 13, 15, 17, 19, 21],
            )
            if self.image_embedding.emb.obs["label"].isnull().sum() > 0:
                self.image_embedding.predict_labels_with_knn(
                    n_neighbors=optimal_k, reclassify_all_data=relabeling, svd_comp=None, weights="uniform", metric="minkowski"
                )
                self.image_embedding.barplots_predicted_clusters(
                    groupby="knn_predicted_label", palette_label=self.dataset.PALETTE, threshold_invasive_dashed_line=True
                )
            else:
                self.image_embedding.emb.obs["knn_predicted_label"] = self.image_embedding.emb.obs["label"]
            if relabeling:
                self.enrichement_score_from_knn()
            self.image_embedding.save_embeddings(saving_path=self.image_embedding_saving_path)
            if "knn_predicted_label" in self.image_embedding.emb.obs.columns:
                self.image_embedding.emb[self.image_embedding.emb.obs["knn_predicted_label"] == self.dataset.ORDER_LABEL[0], :].obs.to_csv(
                    self.selected_invasive_cancer_file
                )

        elif "knn_predicted_label" not in self.image_embedding.emb.obs.columns:
            print("Load select invasive cancer patches results")
            selected_invasive_cancer_df = pd.read_csv(self.selected_invasive_cancer_file, index_col="index")
            # Better to align like this instead with index if we want to give an external selected_invasive_cancer_file whose indexes do not correspond
            self.image_embedding.emb.obs["knn_predicted_label"] = self.image_embedding.emb.obs.merge(
                selected_invasive_cancer_df[
                    ["start_height_origin", "start_width_origin", "name_origin", "knn_predicted_label"]
                ],
                on=["start_height_origin", "start_width_origin", "name_origin"],
                how="left",
            )["knn_predicted_label"]
            self.enrichement_score_from_knn()
        else:
            print("Select invasive cancer patches results already stored in the image embedding object.")

        self.invasive_image_embedding = copy.deepcopy(self.image_embedding)
        self.invasive_image_embedding.emb = self.invasive_image_embedding.emb[
            self.invasive_image_embedding.emb.obs["knn_predicted_label"] == self.dataset.ORDER_LABEL[0], :
        ]

        self.invasive_image_embedding.emb.layers = {}
        self.invasive_image_embedding.emb.uns = {}
        self.invasive_image_embedding.svd = {"denoised": None, "U_df": None, "S": None, "V_df": None}

        print(f"Saving invasive_image_embedding to {self.invasive_image_embedding_saving_path} -- in select invasive cancer pipeline")
        self.invasive_image_embedding.save_embeddings(saving_path=self.invasive_image_embedding_saving_path)

    def mlp_classification_pipeline(self):
        """Run MLP classification on the whole labeled dataset.
        The classification results are saved in the classification folder.
        """
        print("\nStart MLP classification pipeline")
        if not os.path.exists(self.mlp_classification_file):
            print("Start MLP classification analysis")
            self.image_embedding.mlp_classification(
                self.mlp_classification_folder,
            )
        else:
            print("MLP classification results already stored.")

    def load_and_save_rnaseq_data(self):
        """Load and save RNAseq data: gene_embedding.h5ad, invasive_gene_embedding.h5ad. It aligns RNAseq data with their associated patches by assigning the predicted cluster from kmeans to the RNAseq data.
        Then, invasive_gene_embedding is computed by keeping only the RNAseq data related to invasive cancer knn label.
        """
        if not os.path.exists(self.invasive_gene_embedding_saving_path):
            print("\nStart saving invasive RNAseq data")
            gene_embedding = self.dataset.get_gene_embeddings(compute_emb=True, preprocessing=False)
            gene_embedding.emb.obs = gene_embedding.emb.obs.merge(
                self.image_embedding.emb.obs[["predicted_label"]], how="left", left_index=True, right_index=True
            )
            gene_embedding.save_embeddings(self.gene_embedding_saving_path)
            self.invasive_gene_embedding = copy.deepcopy(gene_embedding)
            self.invasive_gene_embedding.emb = self.invasive_gene_embedding.emb[
                ~self.invasive_gene_embedding.emb.obs["predicted_label"].isna(), :
            ]
            self.invasive_gene_embedding.emb.var = pd.DataFrame(index=self.invasive_gene_embedding.emb.var_names)
            self.invasive_gene_embedding.emb.obs["predicted_label_str"] = self.invasive_gene_embedding.emb.obs[
                "predicted_label"
            ].apply(lambda l: "cluster_{}".format(int(l)))
            self.invasive_gene_embedding.save_embeddings(saving_path=self.invasive_gene_embedding_saving_path)
        else:
            print("\nLoading invasive RNAseq data from file: {}".format(self.invasive_gene_embedding_saving_path))
            self.invasive_gene_embedding = self.dataset.get_gene_embeddings(
                compute_emb=False,
                preprocessing=False,
                load=True,
                filename=self.invasive_gene_embedding_saving_path.split("/")[-1].split(".")[0],
            )

        self.invasive_gene_embedding.emb.obs["predicted_label"] = self.invasive_gene_embedding.emb.obs["predicted_label"].apply(
            lambda l: int(l)
        )
        self.invasive_gene_embedding.result_saving_folder = self.model_results_folder
        self.invasive_gene_embedding.saving_plots = True

        print("\nShape before empty spots filtering: {}".format(self.invasive_gene_embedding.emb.shape))
        self.invasive_gene_embedding.emb = self.invasive_gene_embedding.emb[
            (self.invasive_gene_embedding.emb.X.sum(axis=1) != 0), :
        ]
        print("Shape after empty spots filtering: {}".format(self.invasive_gene_embedding.emb.shape))

        if self.gene_embedding_spot_normalization:
            sc.pp.normalize_total(self.invasive_gene_embedding.emb, target_sum=1e6)
            print("\nSpots normalization done.")

    def save_pipeline_settings(self):
        """Save the pipeline attributes that are not of the digitalhistopathology classes and not dataframe to pipeline_setting.json"""
        dict_results = {
            key: value
            for key, value in vars(self).items()
            if (not key.startswith("_"))
            and (not str(value.__class__).split("'")[1].startswith("digitalhistopathology"))
            and (not isinstance(value, pd.DataFrame))
        }
        dict_results["model"] = {"name": self.model.name, "pretrained_model_path": self.model.pretrained_model_path}
        dict_results["image_embedding"] = {
            "name": self.image_embedding.name,
            "shape": self.image_embedding.emb.shape,
        }
        dict_results["invasive_image_embedding"] = {
            "shape": self.invasive_image_embedding.emb.shape,
            "info_optimal_number_of_clusters": self.invasive_image_embedding.info_optimal_number_of_clusters,
        }
        dict_results["invasive_gene_embedding"] = {
            "shape": self.invasive_gene_embedding.emb.shape,
        }

        with open(self.pipeline_settings_path, "w") as fp:
            json.dump(dict_results, fp, cls=NumpyEncoder)
        print("\nSave pipeline settings to {}".format(self.pipeline_settings_path))

    def load_pipeline_settings(self):
        pass


def main():
    parser = argparse.ArgumentParser(description="Pipeline")
    parser.add_argument(
        "--model_name",
        "-m",
        default="uni",
        help="Model name",
        type=str,
        )
    parser.add_argument(
        "--dataset",
        "-d",
        default="HER2",
        help="Dataset name",
        type=str,
        )
    parser.add_argument(
        "--retrained_model_path",
        "-rp",
        default="",
        help="Retrained model path",
        type=str,
        )
    parser.add_argument(
        "--patches_folder",
        "-pf",
        default="results/HER2/compute_patches/all",
        help="Patches folder",
        type=str,
        )
    parser.add_argument(
        "--pipeline_name",
        "-n",
        help="Pipeline name",
        type=str,
        )
    parser.add_argument(
        "--results_folder",
        "-rf",
        default="results",
        help="Pipeline results folder",
        type=str,
        )
    parser.add_argument(
        "--knn_relabeling",
        "-kr",
        default=False,
        help="KNN relabeling",
        type=bool,
        )


    args = parser.parse_args()
    if not args.pipeline_name:
        args.pipeline_name = input("Please enter the pipeline name: ")
    if args.retrained_model_path == "" and args.model_name == "vit":
        args.retrained_model_path = input("Please enter the retrained model path: ")

    model_name = args.model_name.lower()
    model = load_model(model_name, args.retrained_model_path)

    if args.dataset == "TNBC":
        dataset = TNBCDataset(patches_folder=args.patches_folder)
    elif args.dataset == "HER2":
        dataset = HER2Dataset(
            patches_folder=args.patches_folder,
        )

    else:
        dataset = VisiumHDdataset(
            patches_folder=args.patches_folder,
            name=args.dataset
        )



    # Get rid of patient A (outlier)
    if dataset is HER2Dataset:
        dataset.images_filenames = dataset.images_filenames[6:]
        dataset.genes_spots_filenames = dataset.genes_spots_filenames[6:]
        dataset.genes_count_filenames = dataset.genes_count_filenames[6:]
        dataset.samples_names = dataset.samples_names[6:]
        dataset.label_filenames = dataset.label_filenames[1:]

    if not os.path.exists(args.patches_folder):
        warnings.warn("Patches folder does not exist. Patches will be computed.")
    
    results_folder = os.path.join(args.results_folder, args.dataset, "pipeline", args.pipeline_name)
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    
    pipeline = Pipeline(
        name=args.pipeline_name,
        model=model,
        results_folder=results_folder,
        patches_folder=(
            args.patches_folder if os.path.exists(args.patches_folder) else None
        ),
        selected_invasive_cancer_file=None,
        dataset=dataset,
        svd_components_number_shannon_entropy=SVD_COMPONENT_NUMBER,  # to be fixed to the smallest number of deep features across all the tested models
    )

    pipeline.run(knn_relabeling=args.knn_relabeling)

SVD_COMPONENT_NUMBER = 512
if __name__ == "__main__":
    main()
