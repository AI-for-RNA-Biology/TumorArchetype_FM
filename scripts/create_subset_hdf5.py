#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import os
import argparse

import sys
sys.path.append(os.getcwd())

from digitalhistopathology.patch_generator import PatchGenerator


def main():
    parser = argparse.ArgumentParser(description="Create a subset HDF5 file with selected patches.")
    parser.add_argument(
        "--original_hdf5",
        "-o",
        required=True,
        help="Path to the original HDF5 file.",
        default="results/HER2/compute_patches/all/patches.hdf5",
        type=str,
    )
    parser.add_argument(
        "--csv_path",
        "-c",
        required=True,
        help="Path to the CSV file containing the names of the patches to include.",
        default="results/HER2/pipeline/uni/select_invasive_cancer/selected_invasive_cancer.csv",
        type=str,
    )
    parser.add_argument(
        "--new_hdf5",
        "-n",
        required=True,
        help="Path to save the new HDF5 file.",
        default="results/HER2/compute_patches/invasive/patches.hdf5",
        type=str,
    )

    args = parser.parse_args()

    patch_generator = PatchGenerator(saving_folder="")
    patch_generator.create_subset_hdf5(
        original_hdf5_path=args.original_hdf5,
        csv_path=args.csv_path,
        new_hdf5_path=args.new_hdf5,
    )

if __name__ == "__main__":
    main()
