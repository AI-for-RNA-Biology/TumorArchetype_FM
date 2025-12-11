#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import array
from PIL import Image
from json import load
import pandas as pd
import rpy2.robjects as ro
import rpy2.robjects.vectors as rvec
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import numpy as np
import matplotlib.pyplot as plt
Image.MAX_IMAGE_PIXELS = None
import os
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
from natsort import natsorted

import sys
sys.path.append("../")
from digitalhistopathology.patch_generator import PatchGenerator




# Load EBImage package from R
EBImage = importr('EBImage')

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
# Activate the automatic conversion of R objects to pandas objects
pandas2ri.activate()

annotated_only = True

#paths
dataDir = "/idiap/resource/database/ST_TNBC_v3"
cancer_percentage_thresholed = 0.5

if annotated_only:
    images_name_list = natsorted([os.path.basename(f).replace('.png', '.jpg') for f in glob.glob(f"{dataDir}/Images/imageAnnotations/*.png")])
    images_filenames_spot_patches = [f"{dataDir}/Images/imagesHD/{image}" for image in images_name_list]
else:
    images_filenames_spot_patches = natsorted(glob.glob(f"{dataDir}/Images/imagesHD/*.jpg"))
    images_name_list = [os.path.basename(f) for f in images_filenames_spot_patches]


spots_path = [f"{dataDir}/byArray/{image.split('_')[1]}/{image.split('_')[2][:2]}/allSpots.RDS" for image in images_name_list]
# Initialize an empty list to store DataFrames
spots_dfs = []

# Process each RDS file
idxs_to_remove = []
for spots in spots_path:
    if not os.path.exists(spots):
        print(f"File not found: {spots}")
        idx = spots_path.index(spots)
        idxs_to_remove.append(idx)

        continue

    try:
        # Read the RDS file
        spots_df = robjects.r['readRDS'](spots)

        spots_df = convert_to_df(spots_df)
        # Add a column to track the source file
        spots_df['pixel_x'] = spots_df['pixel_x'] * 10.0 / 3.0 # This is to get the coordinates in the HD image (31*744*31744 px). Reminder: the large image is 9523*9523px.
        spots_df['pixel_y'] = spots_df['pixel_y'] * 10.0 / 3.0 # This is to get the coordinates in the HD image (31*744*31744 px). Reminder: the large image is 9523*9523px.
        # Append to the list of DataFrames
        spots_dfs.append(spots_df)

    except Exception as e:
        print(f"Error processing {spots}: {e}")
spots_columns = spots_df.columns


for idx in idxs_to_remove:
    del spots_path[idx]
    del images_filenames_spot_patches[idx]
    del images_name_list[idx]
    
print(images_filenames_spot_patches )
print(spots_path)

inter_distances = []
for i in range(2, 65):
    inter_distances.append(spots_df.groupby("x").mean().loc[i+1]['pixel_x'] - spots_df.groupby("x").mean().loc[i]['pixel_x'])
    
for i in range(2, 62):
    inter_distances.append(spots_df.groupby("y").mean().loc[i+1]['pixel_y'] - spots_df.groupby("y").mean().loc[i]['pixel_y'])
    
center_to_center_distance = np.mean(inter_distances)

## Center to center distance: 150um (according to the paper)
resolution = 150/center_to_center_distance
print(f"Resolution: {resolution:.2f}um/px")

## Spot diameter: 100um (according to the paper)
spot_diameter = round(100/resolution)
print(f"Spot diameter: {round(spot_diameter)} px")

# Generate all RDS file paths
if not annotated_only:
    annotated_names = natsorted([os.path.basename(f).replace('.png', '') for f in glob.glob(f"{dataDir}/Images/imageAnnotations/TNBC*.png")])

    annotsBySpot = natsorted(glob.glob(f"{dataDir}/Robjects/annotsBySpot/TNBC*.RDS"))
    print(annotsBySpot)
    # Initialize an empty list to store DataFrames
    all_dataframes = []

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
                loaded_data = robjects.r['readRDS'](annots)

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
    # annots_collumn = annots_df.columns
    # # Ensure the relevant columns are present
    # cancer_columns = ['Tumor', 'in situ', 'Tumor region']
    # cancer_spots = []
    # for merged_dfs in all_dataframes:
    #     merged_dfs['total_pixels_excluding_nothing'] = merged_dfs.drop(columns=[*spots_columns,'Nothing','source_file']).sum(axis=1)

    #     # Sum the values in the cancer-related columns
    #     merged_dfs['cancer_pixels'] = merged_dfs[cancer_columns].sum(axis=1)

    #     cancer_spots.append(merged_dfs[ merged_dfs['cancer_pixels'] / merged_dfs['total_pixels_excluding_nothing'] > cancer_percentage_thresholed])

else:
    # Generate all RDS file paths
    annotated_names = natsorted([os.path.basename(f).replace('.png', '') for f in glob.glob(f"{dataDir}/Images/imageAnnotations/TNBC*.png")])

    annotsBySpot = natsorted(glob.glob(f"{dataDir}/Robjects/annotsBySpot/TNBC*.RDS"))
    print(annotsBySpot)
    # Initialize an empty list to store DataFrames
    all_dataframes = []

    # Process each RDS file
    for spots_df, annots in zip(spots_dfs,annotsBySpot):

        try:
            # Read the RDS file
            loaded_data = robjects.r['readRDS'](annots)

            # Extract 'annots'
            names = loaded_data.names
            for i, name in enumerate(names):
                if name == 'annots':
                    annots_rds = loaded_data[i]
            annots_df = convert_to_df(annots_rds)

            spots_df = spots_df[spots_df['selected'] == 1]  # Filter spots to only include selected ones
            merged_dfs = pd.merge(annots_df, spots_df, left_index=True, right_index=True, suffixes=('_spot', '_annot'), how='outer')
            merged_dfs['source_file'] = os.path.basename(annots).replace('.RDS', '')

            # Add a column to track the source file
            
            # Append to the list of DataFrames
            all_dataframes.append(merged_dfs)

        except Exception as e:
            print(f"Error processing {spots}: {e}")# Combine all DataFrames into one giant DataFrame
    annots_collumn = annots_df.columns
    # Ensure the relevant columns are present
    cancer_columns = ['Tumor', 'in situ', 'Tumor region']
    cancer_spots = []
    for merged_dfs in all_dataframes:
        merged_dfs['total_pixels_excluding_nothing'] = merged_dfs.drop(columns=[*spots_columns,'Nothing','source_file']).sum(axis=1)

        # Sum the values in the cancer-related columns
        merged_dfs['cancer_pixels'] = merged_dfs[cancer_columns].sum(axis=1)

        cancer_spots.append(merged_dfs[ merged_dfs['cancer_pixels'] / merged_dfs['total_pixels_excluding_nothing'] > cancer_percentage_thresholed])

#Here we generate the patches corresponding to the spots ONLY. We do not add random patches.
# Initialize the PatchGenerator
patch_generator = PatchGenerator(
    images_filenames_random_patches=images_filenames_spot_patches,
    patch_size_pixels=spot_diameter,
    patch_size_micron=None,
    patches_number=900,
    overlap_pixels=0,
    extension="tiff",
    filter_background=True,
    #saving_folder="../results/compute_patches/TNBC_fixed/",
    saving_folder="../results/compute_patches/TNBC/",
    name_patches_info_file="patches_info.pkl.gz",
    images_filenames_spot_patches=images_filenames_spot_patches,
    spots_df=all_dataframes,
    spot_mask=False,
    spot_diameter=spot_diameter,
    filter_with_neighbours=False,
    log_file="../logs/compute_patches.log",
)
patch_generator.compute_all_patches()