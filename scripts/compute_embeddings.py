#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import os
import glob
import warnings
import numpy as np
import argparse
import time
from PIL import Image

import sys
sys.path.append(os.getcwd())

from digitalhistopathology.embeddings.image_embedding import ImageEmbedding
from digitalhistopathology.models import load_model



def resample_image(image_path, scale_factor, output_path):
    with Image.open(image_path) as img:
        new_dimensions = (int(img.width / scale_factor), int(img.height / scale_factor))
        downsampled_image = img.resize(new_dimensions, Image.LANCZOS)
        downsampled_image.save(output_path, 'TIFF')


def check_mpp_and_resample(directory, patches_info, target_mpp=1.0, output_directory=None):
    for filename in os.listdir(directory):
        if filename.endswith('.tiff'):
            file_path = os.path.join(directory, filename)
            mpp = [p for p in patches_info if p['name'] == filename.split('.tiff')[0]][0]['mpp_width']
            scale_factor = target_mpp / mpp
            output_path = os.path.join(output_directory, filename)
            resample_image(file_path, scale_factor, output_path)
            print(f"Resampled {filename} to {target_mpp} micron per pixel and saved as {output_path}")


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
        help="dataset name",
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
        default="../results/HER2/compute_patches/all",
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
        default="../results/",
        help="Pipeline results folder",
        type=str,
        )

    args = parser.parse_args()
    if not args.pipeline_name:
        args.pipeline_name = input("Please enter the pipeline name: ")
    if args.retrained_model_path == "" and args.model_name == "vit":
        args.retrained_model_path = input("Please enter the retrained model path: ")


    model_name = args.model_name.lower()
    model = load_model(model_name, args.retrained_model_path)

    if glob.glob(args.patches_folder + "/*.tiff"):
        patches_filenames = sorted(glob.glob(args.patches_folder + "/*.tiff"))
    else:
        patches_filenames = glob.glob(os.path.join(args.patches_folder, "*.hdf5"))
        
    print(f"Patches filenames: {patches_filenames}", flush=True)

    patches_info_filename = os.path.join(args.patches_folder, "patches_info.pkl.gz")
    
    results_folder = os.path.join(args.results_folder, args.dataset, "pipeline", args.pipeline_name)
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

    image_embedding = ImageEmbedding(patches_filenames=patches_filenames,
                                     patches_info_filename=patches_info_filename,
                                     pretrained_model=model,
                                     name=args.pipeline_name,
                                     result_saving_folder=results_folder,
                                     saving_plots=True)
    
    start = time.time()
    image_embedding.compute_embeddings()
    end = time.time()

    image_embedding.save_embeddings(saving_path=os.path.join(results_folder, "embeddings.h5ad"))

    print("Time to compute embeddings for {}: {}".format(model.name, end - start))



if __name__ == "__main__":
    main()