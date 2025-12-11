#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import sys

sys.path.append("../")
import glob

from digitalhistopathology.datasets.real_datasets import HER2Dataset
from digitalhistopathology.patch_generator import PatchGenerator


images_filenames = sorted(glob.glob("../data/HER2_breast_cancer/images/HE/*.jpg"))[6:]
print(images_filenames)

genes_spots_files = sorted(glob.glob("../data/HER2_breast_cancer/spot-selections/*"))[6:]
print(genes_spots_files)
her2_dataset = HER2Dataset("../results/compute_patches/her2_final_without_A/")
p = PatchGenerator(
    images_filenames_random_patches=images_filenames,
    patch_size_pixels=her2_dataset.spot_diameter,
    patch_size_micron=None,
    patches_number=None,
    overlap_pixels=0,
    extension="tiff",
    filter_background=True,
    saving_folder="../results/compute_patches/her2_final_without_A/",
    name_patches_info_file="patches_info.pkl.gz",
    images_filenames_spot_patches=images_filenames,
    spots_filenames=genes_spots_files,
    spot_mask=False,
    spot_diameter=her2_dataset.spot_diameter,
    filter_with_neighbours=False,
    log_file="../logs/compute_patches.log",
)

p.compute_all_patches()
