#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import numpy as np
from skimage.color import rgb2gray
from skimage.feature import graycomatrix, graycoprops
import os
import json
import datetime
from multiprocessing import Pool
from tqdm import tqdm
import pandas as pd
import mahotas
from PIL import Image, ImageDraw
import openslide
from scipy.stats import skew, kurtosis
from skimage.measure import shannon_entropy



def get_whole_image_texture_features(patch_img):
    grayco_features = {}
    angles = [0, np.pi/4, np.pi/2, 3*np.pi/4]

    gray_image = (rgb2gray(np.array(patch_img)[:, :, :3]) * 255).astype(np.uint8)
    # Compute the gray-level co-occurrence matrix (GLCM)
    glcm = graycomatrix(gray_image, [1], angles, 256, symmetric=True, normed=True)

    properties = ['ASM', 'contrast', 'dissimilarity', 'homogeneity', 'correlation']
    grayco_results = {}

    for prop in properties:
        grayco_results[prop] = graycoprops(glcm, prop)[0]

        grayco_features[f'{prop}'] = np.mean(grayco_results[prop])

    return grayco_features


def apply_mask(image, mask):
    # Ensure the mask is binary
    mask = mask > 0
    # Apply the mask
    masked_image = image * np.expand_dims(mask, axis=-1)
    return masked_image


def process_zernike_one_cell(args):
    cell, scMTOP_json, mask = args
    y1, x1 = scMTOP_json['nuc'][cell]['bbox'][0]
    y2, x2 = scMTOP_json['nuc'][cell]['bbox'][1]
    roi = np.array(mask)[y1:y2, x1:x2]
    zernike_moments = mahotas.features.zernike_moments(roi, radius=10, degree=5)
    return cell, zernike_moments


def get_zernike_moments_per_nuclei(wsi_path, wsi_cellvit_path, wsi_name, num_cores=None):


    # Load cellvit data
    path_to_scMTOP_json = os.path.join(wsi_cellvit_path, wsi_name, 'cell_detection', 'cells_for_scMTOP.json')

    with open(path_to_scMTOP_json, 'r') as f:
        scMTOP_json = json.load(f)

    # Load wsi
    wsi = openslide.OpenSlide(wsi_path)

    # Get the mask
    mask = Image.new('L', wsi.level_dimensions[0], 255)
    draw = ImageDraw.Draw(mask)

    zernike_cells = {}
    for cell in scMTOP_json['nuc']:
        draw.polygon([tuple(point) for point in scMTOP_json['nuc'][cell]['contour']], outline=1, fill=1)
    
    print(f"Finished to draw the mask for {wsi_name}: {datetime.datetime.now()}", flush=True)

    cells = list(scMTOP_json['nuc'].keys())
    args = [(cell, scMTOP_json, mask) for cell in cells]
    
    with Pool(processes=num_cores) as pool:
        results = list(tqdm(pool.imap(process_zernike_one_cell, args), total=len(cells)))

    zernike_cells = {cell: moments for cell, moments in results}
    
    return zernike_cells

def get_patch_zernike_moments(patch_name, patch_to_cell, zernike_cells):

    patch_zernike_moments = {}
    for cell in zernike_cells:
        if cell in patch_to_cell[patch_name]:
            patch_zernike_moments[cell] = zernike_cells[cell]
    
    df_zernike = pd.DataFrame(patch_zernike_moments)

    zernike_features = {}

    for i in range(df_zernike.shape[0]):
        zernike_features[f'zernike_moment_{i+1}_mean'] = df_zernike.iloc[i].mean()
        zernike_features[f'zernike_moment_{i+1}_std'] = df_zernike.iloc[i].std()

    return zernike_features

def get_color_features(pixels_array):
    
    """
    Extracts various color and intensity features from the given extracellular pixels.
    Parameters:
    -----------
    non_cell_pixels : np.ndarray
        A 2D array where each row represents a pixel with its color channels (R, G, B).
    Returns:
    --------
    dict
        A dictionary containing the following features:
        - mean_color_R: Mean value of the red color channel.
        - mean_color_G: Mean value of the green color channel.
        - mean_color_B: Mean value of the blue color channel.
        - mean_intensity: Mean intensity of the color channels.
        - std_color_R: Standard deviation of the red color channel.
        - std_color_G: Standard deviation of the green color channel.
        - std_color_B: Standard deviation of the blue color channel.
        - std_intensity: Standard deviation of the intensity.
        - skew_color_R: Skewness of the red color channel.
        - skew_color_G: Skewness of the green color channel.
        - skew_color_B: Skewness of the blue color channel.
        - kurtosis_color_R: Kurtosis of the red color channel.
        - kurtosis_color_G: Kurtosis of the green color channel.
        - kurtosis_color_B: Kurtosis of the blue color channel.
        - entropy_color_R: Entropy of the red color channel.
        - entropy_color_G: Entropy of the green color channel.
        - entropy_color_B: Entropy of the blue color channel.
    """

    
    # Mean color
    mean_color_R, mean_color_G, mean_color_B = np.mean(pixels_array, axis=0)

    # Mean intensity
    mean_intensity = np.mean(np.mean(pixels_array[:, :3], axis=1))
    
    # Standard deviation of color
    std_color_R, std_color_G, std_color_B = np.std(pixels_array, axis=0)
    
    # Variance of intensity
    std_intensity = np.std(np.mean(pixels_array[:, :3], axis=1))
    
    # Skewness of color channels
    skew_color_R, skew_color_G, skew_color_B = skew(pixels_array[:, :3], axis=0)
    
    # Kurtosis of color channels
    kurtosis_color_R, kurtosis_color_G, kurtosis_color_B = kurtosis(pixels_array[:, :3], axis=0)
    
    # Entropy
    entropy_R = shannon_entropy(pixels_array[:, 0])
    entropy_G = shannon_entropy(pixels_array[:, 1])
    entropy_B = shannon_entropy(pixels_array[:, 2])
        
    return {
        'mean_color_R': mean_color_R,
        'mean_color_G': mean_color_G,
        'mean_color_B': mean_color_B,
        'mean_intensity': mean_intensity,
        'std_color_R': std_color_R,
        'std_color_G': std_color_G,
        'std_color_B': std_color_B,
        'std_intensity': std_intensity,
        'skew_color_R': skew_color_R,
        'skew_color_G': skew_color_G,
        'skew_color_B': skew_color_B,
        'kurtosis_color_R': kurtosis_color_R,
        'kurtosis_color_G': kurtosis_color_G,
        'kurtosis_color_B': kurtosis_color_B,
        'entropy_color_R': entropy_R,
        'entropy_color_G': entropy_G,
        'entropy_color_B': entropy_B,
        }

