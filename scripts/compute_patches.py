#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#

import os
import glob
import array
import argparse
from PIL import Image
from json import load
import pandas as pd
import numpy as np
from natsort import natsorted

import sys
sys.path.append(os.getcwd())
from digitalhistopathology.datasets.real_datasets import HER2Dataset, TNBCDataset
from digitalhistopathology.patch_generator import PatchGenerator

Image.MAX_IMAGE_PIXELS = None

# --- R / EBImage imports for TNBC ---
import rpy2.robjects as ro
import rpy2.robjects.vectors as rvec
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

EBImage = importr('EBImage')
pandas2ri.activate()


### TNBC PREPROCESSING FUNCTIONS

def convert_to_df(data):
    # Check if data is an R Matrix
    if isinstance(data, rvec.Matrix):
        return pandas2ri.rpy2py(ro.r['as.data.frame'](data))
    elif isinstance(data, rvec.Vector):
        return pandas2ri.rpy2py(data)
    # Check if data is an R DataFrame
    elif isinstance(data, rvec.DataFrame):
        return pandas2ri.rpy2py(data)

    # Check if data is an EBImage object (matrix or image)
    elif isinstance(data, ro.rinterface.Sexp):
        # Check if it's an EBImage object (usually matrix)
        try:
            ro.r('library(EBImage)')
            image_array = ro.r['imageData'](data)  # Convert EBImage object to matrix
            if isinstance(image_array, rvec.Matrix):
                # Convert the image matrix to a DataFrame (pixel values)
                return pd.DataFrame(pandas2ri.rpy2py(image_array))
        except Exception as e:
            print(f"Error converting EBImage object: {e}")
            return None

    # Return None for unsupported data types
    else:
        return None


def preprocess_tnbc_images(dataDir, annotated_only=True, cancer_percentage_thresholed=0.5):
    
    # Set the image names
    if annotated_only:
        images_name_list = natsorted([os.path.basename(f).replace('.png', '.jpg') for f in glob.glob(f"{dataDir}/Images/imageAnnotations/*.png")])
        images_filenames_spot_patches = [f"{dataDir}/Images/imagesHD/{image}" for image in images_name_list]
    else:
        images_filenames_spot_patches = natsorted(glob.glob(f"{dataDir}/Images/imagesHD/*.jpg"))
        images_name_list = [os.path.basename(f) for f in images_filenames_spot_patches]

    def convert_to_df(data):
        # Check if data is an R Matrix
        if isinstance(data, rvec.Matrix):
            return pandas2ri.rpy2py(ro.r['as.data.frame'](data))
        elif isinstance(data, rvec.Vector):
            return pandas2ri.rpy2py(data)
        # Check if data is an R DataFrame
        elif isinstance(data, rvec.DataFrame):
            return pandas2ri.rpy2py(data)

        # Check if data is an EBImage object (matrix or image)
        elif isinstance(data, ro.rinterface.Sexp):
            # Check if it's an EBImage object (usually matrix)
            try:
                ro.r('library(EBImage)')
                image_array = ro.r['imageData'](data)  # Convert EBImage object to matrix
                if isinstance(image_array, rvec.Matrix):
                    # Convert the image matrix to a DataFrame (pixel values)
                    return pd.DataFrame(pandas2ri.rpy2py(image_array))
            except Exception as e:
                print(f"Error converting EBImage object: {e}")
                return None

        # Return None for unsupported data types
        else:
            return None

    # Set the spots paths
    spots_path = [f"{dataDir}/byArray/{image.split('_')[1]}/{image.split('_')[2][:2]}/allSpots.RDS" for image in images_name_list]
    
    # Initialize an empty list to store DataFrames
    spots_dfs = []

    # Process each RDS file
    idxs_to_remove = []
    for i, spots in enumerate(spots_path):

        if not os.path.exists(spots):
            print(f"File not found: {spots}")
            idxs_to_remove.append(i)
            continue

        try:
            # Read the RDS file
            spots_df = ro.r['readRDS'](spots)
            spots_df = convert_to_df(spots_df)
            # Add a column to track the source file
            spots_df['pixel_x'] = spots_df['pixel_x'] * 10.0 / 3.0 # This is to get the coordinates in the HD image (31*744*31744 px). Reminder: the large image is 9523*9523px.
            spots_df['pixel_y'] = spots_df['pixel_y'] * 10.0 / 3.0 # This is to get the coordinates in the HD image (31*744*31744 px). Reminder: the large image is 9523*9523px.
            # Append to the list of DataFrames
            spots_dfs.append(spots_df)
        
        except Exception as e:
            print(f"Error processing {spots}: {e}")
    
    # Remove the missing RDS files 
    for idx in sorted(idxs_to_remove, reverse=True):
        del spots_path[idx]
        del images_filenames_spot_patches[idx]
        del images_name_list[idx]

    # Compute spot diameter
    inter_distances = []
    spots_df = spots_dfs[0]

    for i in range(2, 65):
        inter_distances.append(spots_df.groupby("x").mean().loc[i+1]['pixel_x'] - spots_df.groupby("x").mean().loc[i]['pixel_x'])
        
    for i in range(2, 62):
        inter_distances.append(spots_df.groupby("y").mean().loc[i+1]['pixel_y'] - spots_df.groupby("y").mean().loc[i]['pixel_y'])
        
    center_to_center_distance = np.mean(inter_distances)

    # Center to center distance: 150um (according to the paper)
    resolution = 150/center_to_center_distance
    print(f"Resolution: {resolution:.2f}um/px")

    # Spot diameter: 100um (according to the paper)
    spot_diameter = round(100/resolution)
    print(f"Spot diameter: {spot_diameter} px")
    
    # Generate all RDS file paths
    all_dataframes = []
    annotsBySpot = natsorted(glob.glob(f"{dataDir}/Robjects/annotsBySpot/TNBC*.RDS"))

    if not annotated_only:
        annotated_names = natsorted([os.path.basename(f).replace('.png', '') for f in glob.glob(f"{dataDir}/Images/imageAnnotations/TNBC*.png")])

        # Process each RDS file
        for idx, img_name in enumerate(images_name_list):
            spots_df = spots_dfs[idx]
            
            if img_name.replace('.jpg', '') not in annotated_names:
                spots_df['source_file'] = img_name.replace('.jpg', '')
                all_dataframes.append(spots_df)
            else:
                annots = f"{dataDir}/Robjects/annotsBySpot/{img_name.replace('.jpg', '').split('_')[0]}.RDS"

                try:
                    # Read the RDS file
                    loaded_data = ro.r['readRDS'](annots)

                    # Extract 'annots'
                    names = loaded_data.names
                    for i, name in enumerate(names):
                        if name == 'annots':
                            annots_rds = loaded_data[i]
                    annots_df = convert_to_df(annots_rds)
                    # annots_df[
                    #     'source_file'] = os.path.basename(annots).replace('.RDS', '')
                    
                    merged_dfs = pd.merge(annots_df, spots_df, left_index=True, right_index=True, suffixes=('_spot', '_annot'), how='outer')
                    merged_dfs['source_file'] = img_name.replace('.jpg', '')

                    # Add a column to track the source file
                    
                    # Append to the list of DataFrames
                    all_dataframes.append(merged_dfs)

                except Exception as e:
                    print(f"Error processing {spots}: {e}")# Combine all DataFrames into one giant DataFrame
        
    else:
        
        # Process each RDS file
        for spots_df, annots in zip(spots_dfs, annotsBySpot):

            try:
                # Read the RDS file
                loaded_data = ro.r['readRDS'](annots)

                # Extract 'annots'
                names = loaded_data.names
                for i, name in enumerate(names):
                    if name == 'annots':
                        annots_rds = loaded_data[i]
                annots_df = convert_to_df(annots_rds)

                spots_df = spots_df[spots_df['selected'] == 1]  # Filter spots to only include selected ones
                merged_dfs = pd.merge(annots_df, spots_df, left_index=True, right_index=True, suffixes=('_spot', '_annot'), how='outer')
                merged_dfs['source_file'] = os.path.basename(annots).replace('.RDS', '')
                
                # Append to the list of DataFrames
                all_dataframes.append(merged_dfs)

            except Exception as e:
                print(f"Error processing {spots}: {e}")
        
    return images_filenames_spot_patches, all_dataframes, spot_diameter


### MAIN FUNCTION

def create_patch_generator(dataset_name):

    print(f"Creating patch generator for {dataset_name}...")

    if dataset_name == "HER2":

        # Use the dataset class to extract spot diameter
        saving_folder = "results/compute_patches/HER2/"
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

        saving_folder = "results/compute_patches/TNBC/"
        dataset = TNBCDataset(patches_folder=saving_folder, dataDir="data/TNBC")

        images_filenames = dataset.images_filenames
        all_dataframes = dataset.all_dataframes
        spot_diameter = dataset.spot_diameter


        #images_filenames, all_dataframes, spot_diameter = preprocess_tnbc_images(
         #   "data/TNBC", 
          #  annotated_only=True, 
           # cancer_percentage_thresholed=0.5)
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
    parser.add_argument("--dataset_name", "-m", help="Dataset name", type=str, default="HER2")
    args = parser.parse_args()

    p = create_patch_generator(dataset_name=args.dataset_name)
    p.compute_all_patches()
