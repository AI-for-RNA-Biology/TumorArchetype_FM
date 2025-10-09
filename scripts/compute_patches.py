#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#

import os
import glob
import argparse

import sys
sys.path.append(os.getcwd())
from digitalhistopathology.datasets.real_datasets import HER2Dataset, TNBCDataset
from digitalhistopathology.patch_generator import PatchGenerator



def create_patch_generator(dataset_name):

    print(f"Creating patch generator for {dataset_name}...")

    if dataset_name == "HER2":

        # Use the dataset class to extract spot diameter
        saving_folder = "results/HER2/compute_patches/all/"
        dataset = HER2Dataset(saving_folder)
        spot_diameter = dataset.spot_diameter

        # Load the images and spots dataframes without patient A
        images_filenames = sorted(glob.glob("data/HER2/images/HE/*.jpg"))[6:]
        print(images_filenames)
        genes_spots_files = sorted(glob.glob("data/HER2/spot-selections/*"))[6:]
        print(genes_spots_files)

        # Set the other parameters
        patches_number = None
        spots_df = None

    elif dataset_name == "TNBC":

        saving_folder = "results/TNBC/compute_patches/all/"
        dataset = TNBCDataset(patches_folder=saving_folder, dataDir="data/TNBC")
        dataset.preprocess(annotated_only=True)

        images_filenames = dataset.images_filenames
        all_dataframes = dataset.all_dataframes
        spot_diameter = dataset.spot_diameter
        print(images_filenames)

        genes_spots_files = None

        patches_number = 900
        spots_df = all_dataframes

    else:
        raise Exception("Dataset not found")

    p = PatchGenerator(
        images_filenames_random_patches=images_filenames,
        patch_size_pixels=spot_diameter,
        patch_size_micron=None,
        patches_number=patches_number,
        overlap_pixels=0,
        extension="tiff",
        filter_background=True,
        saving_folder=saving_folder,
        name_patches_info_file="patches_info.pkl.gz",
        images_filenames_spot_patches=images_filenames,
        spots_filenames=genes_spots_files,
        spots_df=spots_df,
        spot_mask=False,
        spot_diameter=spot_diameter,
        filter_with_neighbours=False,
        log_file="logs/compute_patches.log",
    )

    return p


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Computing patches")
    parser.add_argument("--dataset", "-d", help="Dataset name", type=str, default="HER2")
    args = parser.parse_args()

    p = create_patch_generator(dataset_name=args.dataset)
    p.compute_all_patches()
