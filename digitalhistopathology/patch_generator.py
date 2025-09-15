#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import glob
import gzip
import itertools
import logging
from math import e
import os
import pickle
import random

import matplotlib
import numpy as np
import openslide
import pandas as pd
import scipy.stats as sps
import h5py

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from openslide.deepzoom import DeepZoomGenerator
from PIL import ExifTags, Image, ImageFile
from skimage.measure import shannon_entropy
from sklearn.mixture import BayesianGaussianMixture, GaussianMixture

Image.MAX_IMAGE_PIXELS = None
random.seed(42)
np.random.seed(42)

ImageFile.LOAD_TRUNCATED_IMAGES = True

class PatchGenerator:
    """Class to tile whole-slide images into patches."""

    SHANNON_ENTROPY_THRESHOLD = 6.451
    # SHANNON_ENTROPY_THRESHOLD = 5.5

    def __init__(
        self,
        images_filenames_random_patches=None,
        patch_size_pixels=224,
        patch_size_micron=None,
        patches_number=None,
        overlap_pixels=0,
        extension="tiff",
        filter_background=True,
        saving_folder="../results/patches_test/",
        name_patches_info_file="patches_info.pkl.gz",
        images_filenames_spot_patches=None,
        spots_filenames=None,
        spot_mask=True,
        spot_diameter=None,
        filter_with_neighbours=True,
        log_file="../logs/compute_patches.log",
        spots_df=None,
    ):
        """
        Args:
            images_filenames_random_patches (list): List of images files from which compute random patches
            patch_size_pixels (int, optional): Length in pixels of the square patch side, is taken into account only if patch_size_micron is None. Defaults to 224.
            patch_size_micron (float, optional): Length in microns of the square patch side. The patch_size_pixels is adapted from it. Defaults to None.
            patches_number (int, optional):  Number of patches you want. The patches number can be inferior to the asked number, because of lack os possibility with the overlap parameter. Defaults to None, which means all possible patches.
            overlap_pixels (int, optional): Number of overlapping pixels you allow for the patches computation. Better to give an even number because the pixels are shared between columns and rows. Defaults to 0.
            extension (str, optional): Extension of your patches files. Defaults to "tiff".
            filter_background (bool, optional): If you want to have only foreground patches. Defaults to True.
            saving_folder (str, optional): Where to save the patches and the information file. Defaults to "../results/patches_test/".
            name_patches_info_file (str, optional): Name of the information file. Need to end with .pkl.gz. Defaults to "patches_info.pkl.gz".
            images_filenames_spot_patches (list, optional): List of image files from which spot patches are computed. Defaults to None.
            spots_filenames (list, optional): List of spots location filename corresponding to the images_filenames_spot_patches. Defaults to None.
            spot_mask (bool, optional): Whether or not using the circular mask for the spot patches. Defaults to True.
            spot_diameter (float, optional): Diameter of the spot size on the image in pixels. Defaults to None.
            filter_with_neighbours (bool, optional): Whether to post-filter patches with neighbours threshold. Defaults to True.
            log_file (str, optional): Logging file. Defaults to "../logs/compute_patches.log".
            spots
        """

        self.spots_df = spots_df
        self.log_file = log_file
        if os.path.exists(self.log_file):
            os.remove(self.log_file)

        if not os.path.exists(os.path.dirname(self.log_file)):
            os.mkdir(os.path.dirname(self.log_file))
            
        logging.basicConfig(
            filename=self.log_file,
            encoding="utf-8",
            level=logging.INFO,
            format="%(asctime)s %(levelname)-8s %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

        self.images_filenames_random_patches = images_filenames_random_patches
        self.patch_size_pixels = patch_size_pixels
        self.patch_size_micron = patch_size_micron
        self.patches_number = patches_number
        self.overlap_pixels = overlap_pixels
        self.extension = extension
        self.filter_background = filter_background
        self.saving_folder = saving_folder
        self.name_patches_info_file = name_patches_info_file
        self.images_filenames_spot_patches = images_filenames_spot_patches
        self.spots_filenames = spots_filenames
        self.spot_mask = spot_mask
        self.spot_radius = int(spot_diameter / 2) if spot_diameter else None
        self.filter_with_neighbours = filter_with_neighbours
        self.all_patches_info = []

        if not os.path.exists(self.saving_folder):
            print("Creating saving folder: {}".format(self.saving_folder), flush=True)
            os.makedirs(self.saving_folder, exist_ok=True)

    @staticmethod
    def mpp_from_a_non_pyramidal_image(img, file):
        """Compute the microns per pixel resolutions from an pillow image

        Args:
            img (pillow): Pillow rgb image
            file (str): File of the image to assign correct resolution of there is None in function of the batch

        Returns:
            mmp_width (float): Microns per pixels for the width resolution of the image
            mmp_height (float): Microns per pixels for the height resolution of the image
        """
        try:
            exif = {ExifTags.TAGS[k]: v for k, v in img._getexif().items() if k in ExifTags.TAGS}
            resolution_unit = exif["ResolutionUnit"]
            resolution_width = exif["XResolution"]
            resolution_height = exif["YResolution"]
        except Exception as e:
            logging.warning("No exif resolution for this image, try with in another way: {}".format(e))
            resolution_unit = 2
            resolution_width = img.info.get("dpi", (72, 72))[0]
            resolution_height = img.info.get("dpi", (72, 72))[1]

        if resolution_width == 72 and resolution_height == 72:
            logging.warning("Default value for resolution only!!")
            if "human_melanoma" in file:
                logging.info("Set estimated resolution for human melanoma data from QPath")
                mpp_width = 0.303
                mpp_height = 0.303
            elif "HER2" in file:
                logging.info("Set estimated resolution for HER2 breast cancer data from QPath")
                mpp_width = 1
                mpp_height = 1
            elif "TNBC" in file:
                mpp_height = 0.40
                mpp_width = 0.40
        else:
            if resolution_unit == 2:
                # pixels per inch
                mpp_width = 1 / resolution_width * 25400
                mpp_height = 1 / resolution_height * 25400
            elif resolution_unit == 3:
                # pixels per cm
                mpp_width = 1 / resolution_width * 1e4
                mpp_height = 1 / resolution_height * 1e4
            else:
                logging.error("No resolution unit was found, impossible to compute patches with patche_size_micron not None.")
                return
        return mpp_width, mpp_height

    @staticmethod
    def mpp_from_a_pyramidal_image(slide):
        """Compute the microns per pixel resolutions from an openslide object

        Args:
            slide (openslide): Pyramidal image that can be processed with openslide library

        Returns:
            mmp_width (float): Microns per pixels for the width resolution of the image
            mmp_height (float): Microns per pixels for the height resolution of the image
        """
        slide_props = slide.properties
        mpp_width = float(slide_props["openslide.mpp-x"])
        mpp_height = float(slide_props["openslide.mpp-y"])
        return mpp_width, mpp_height

    @staticmethod
    def get_the_higher_mpp(images_filenames):
        """Return the lower resolution from a list of images files, so in term of microns per pixel, it means the higher one.

        Args:
            images_filenames (list): List of images filenames

        Returns:
            float: The higher microns per pixel
        """
        # TODO: Parallelize
        mpp_list = []
        for file in images_filenames:
            if file.endswith(".svs") or file.endswith(".tiff") or file.endswith(".tif"):
                slide = openslide.open_slide(file)
                mpp_width, mpp_height = mpp_from_a_pyramidal_image(slide)
                slide = None
            else:
                img = Image.open(file).convert("RGB")
                mpp_width, mpp_height = mpp_from_a_non_pyramidal_image(
                    img, human_melanoma_data=file.split("/")[-1].startswith("ST_mel")
                )
                img = None
            logging.debug("Mpp height = {}, Mpp width = {} for {}".format(mpp_height, mpp_width, file))
            assert mpp_width == mpp_height, "Mpp width = {} is different from mpp height = {}".format(mpp_width, mpp_height)
            mpp_list.append(mpp_width)
        return max(mpp_list)

    def random_patches(
        self,
        img_filename,
        debug_jupyter=False,
        benchmark_count=None,
        hdf5_file=None,
    ):
        """Compute random patches from an image (pyramidal or not)

        Args:
            img_filename (str): Image path
            debug_jupyter (bool, optional): To plot the discarded patches in jupyter notebook for debug. Defaults to False.
            benchmark_count: (int, optional): If the benchmark dataset is not computed, set to None. Defaults to None.

        Returns:
            patches_info_list (list): List of the patches additional informations
        """

        logging.info("Start computing patches for {}".format(os.path.abspath(img_filename)))

        self.overlap_pixels = int(self.overlap_pixels)

        pyramidal_image = False

        # Load image and resolution
        if img_filename.endswith(".svs") or img_filename.endswith(".tiff") or img_filename.endswith(".tif"):
            pyramidal_image = True
            slide = openslide.open_slide(img_filename)
            # shape of the whole slide
            shape_origin = (slide.dimensions[1], slide.dimensions[0], 3)
            # resolution
            mpp_width, mpp_height = self.mpp_from_a_pyramidal_image(slide)
        else:
            img = Image.open(img_filename).convert("RGB")
            # Resolution
            mpp_width, mpp_height = self.mpp_from_a_non_pyramidal_image(img, file=img_filename)
            img = np.array(img)
            # shape of the whole slide
            shape_origin = img.shape

        if mpp_height != mpp_width:
            logging.warning("The resolution is not the same for height and width!")
        logging.debug("Mpp height = {}, Mpp width = {}".format(mpp_height, mpp_width))

        # Adapt the pixel size if micron size is given
        if self.patch_size_micron is not None:
            assert mpp_height == mpp_width, "Not the same resolution in each axis"
            self.patch_size_pixels = int(self.patch_size_micron / mpp_height)
            logging.info(
                "Patch size in pixels for {} microns = {} pixels".format(self.patch_size_micron, self.patch_size_pixels)
            )

        # Generate the patches locations
        if pyramidal_image:
            # overlap = the number of extra pixels to add to each interior edge of a tile
            # tile_size = the width and height of a single tile = tile_size + 2 * overlap
            tiles = DeepZoomGenerator(
                slide,
                tile_size=self.patch_size_pixels - 2 * int(self.overlap_pixels / 2),
                overlap=int(self.overlap_pixels / 2),
                limit_bounds=False,
            )
            # tiles at the max resolution
            cols, rows = tiles.level_tiles[tiles.level_count - 1]

            # -1 because the border left and bottom tiles are usually smaller (because the slide is not perfectly divided by our pacth size)
            # we discard it, usually it will be only background, so it does not matter
            all_cols = list(np.arange(0, cols - 1))
            all_rows = list(np.arange(0, rows - 1))
        else:
            step = self.patch_size_pixels - int(self.overlap_pixels / 2)
            assert (
                step > 0
            ), "Step cannot be lower than 1, you fix an overlap pixels number bigger than the patch size: {} > {}".format(
                self.overlap_pixels, self.patch_size_pixels
            )

            # possible heights list
            all_rows = list(np.arange(0, img.shape[0] - self.patch_size_pixels, step))
            # possible widths list
            all_cols = list(np.arange(0, img.shape[1] - self.patch_size_pixels, step))

        # all possible patches location,
        all_tuples_row_col = list(itertools.product(all_rows, all_cols))
        logging.debug("Maximum number of patches without background filtering = {}".format(len(all_tuples_row_col)))
        patches_info_list = []

        if self.patches_number is None:
            patches_number = len(all_tuples_row_col)
        else:
            patches_number = self.patches_number
        
        print(f"Will compute {patches_number} random patches", flush=True)

        while len(patches_info_list) < patches_number and len(all_tuples_row_col) > 0:
            logging.debug(
                "len(patches_info_list) = {}, len(all_tuples_row_col) = {}".format(
                    len(patches_info_list), len(all_tuples_row_col)
                )
            )
            rand_ind_tuple_row_col = np.random.randint(0, len(all_tuples_row_col))
            rand_tuple_row_col = all_tuples_row_col.pop(rand_ind_tuple_row_col)
            logging.debug("(random row, random col) = {}".format(rand_tuple_row_col))

            if pyramidal_image:
                rand_row, rand_col = rand_tuple_row_col
                patch = tiles.get_tile(tiles.level_count - 1, (rand_col, rand_row)).convert("RGB")
                patch = np.array(patch)
            else:
                rand_start_h, rand_start_w = rand_tuple_row_col
                patch = img[
                    rand_start_h : rand_start_h + self.patch_size_pixels,
                    rand_start_w : rand_start_w + self.patch_size_pixels,
                    :,
                ]

            mean = patch.mean()
            std = patch.std()
            median = np.median(patch)
            entropy = shannon_entropy(patch)
            logging.debug(
                "patch mean = {}, patch std = {}, patch median = {}, patch entropy = {}".format(mean, std, median, entropy)
            )

            if self.filter_background:
                patch_ok = entropy > self.SHANNON_ENTROPY_THRESHOLD
            else:
                patch_ok = True

            if patch_ok:
                
                origin_name = img_filename.split("/")[-1].split("." + img_filename.split(".")[-1])[0]
                print(f"Origin name = {origin_name}", flush=True)
                # special case for HER2 breast cancer data
                if len(origin_name) == 2:
                    name = (
                        origin_name[0]
                        + "_rep"
                        + origin_name[1]
                        + "_patch{}".format(str(len(patches_info_list)).zfill(len(str(patches_number - 1))))
                    )
                else:
                    name = origin_name + "_patch{}".format(str(len(patches_info_list)).zfill(len(str(patches_number - 1))))

                if benchmark_count is not None:
                    benchmark_count += 1
                    name = origin_name + "_benchmark{}".format(str(benchmark_count).zfill(len(str(10000))))

                saving_path = os.path.join(self.saving_folder, name + "." + self.extension)

                # Convert to pillow to ensure the correct RGB format of our patches
                if hdf5_file is None:
                    Image.fromarray(patch).save(saving_path)
                    logging.debug("Patch saved to {}".format(os.path.abspath(saving_path)))
                else:
                    print(f"Saving random patch: {name}", flush=True)
                    hdf5_file.create_dataset(name, data=patch)
                    logging.debug("Patch saved to hdf5 file")

                if pyramidal_image:
                    rand_start_w, rand_start_h = tiles.get_tile_coordinates(tiles.level_count - 1, (rand_col, rand_row))[0]

                patch_data = {
                    "path": os.path.abspath(saving_path),
                    "name": name,
                    "mpp_height": mpp_height,
                    "mpp_width": mpp_width,
                    "shape_micron": self.patch_size_micron,
                    "shape_pixel": self.patch_size_pixels,
                    "overlap_pixel": self.overlap_pixels,
                    "path_origin": img_filename,
                    "name_origin": origin_name,
                    "shape_origin": shape_origin,
                    "start_height_origin": rand_start_h,
                    "start_width_origin": rand_start_w,
                    "extension_origin": img_filename.split(".")[-1],
                    "batch": img_filename.split("/data/")[-1].split("/")[0].replace("_", " ").capitalize(),
                    "tumor": name.split("_patch")[0].split(".")[0].split("_rep")[0],
                    "mean_intensity": mean,
                    "median_intensity": median,
                    "std_intensity": std,
                    "entropy_intensity": entropy,
                }
                patches_info_list.append(patch_data)
                logging.debug("patch info: {}".format(patch_data))
            else:
                if debug_jupyter:
                    print("patch mean = {}, patch std = {}, patch median = {}".format(mean, std, median))
                    print(patch.std())
                    plt.imshow(patch)
                    plt.show()

            patch = None

        logging.info("Finish computing {} patches for {}".format(len(patches_info_list), os.path.abspath(img_filename)))

        if benchmark_count is None:
            return patches_info_list
        else:
            return patches_info_list, benchmark_count

    def check_neighbours(self, t, all_values, neighbours_level=3):
        """Compute a score for a patch based on the number of direct neighbors

        Args:
            t (tuple): Top left patch coordinates in pixels (height, width)
            all_values (_type_): List of the tuples coordinates from all patches from the same image
            neighbours_level (int, optional): The max distance level between the current patch and its neighbours. Defaults to 3.

        Returns:
            int: the number of direct neighbors the patch has
        """
        step = self.patch_size_pixels - int(self.overlap_pixels / 2)
        x_list = [t[0] + i * step for i in range(-neighbours_level, neighbours_level + 1)]
        y_list = [t[1] + i * step for i in range(-neighbours_level, neighbours_level + 1)]
        neighbours = list(itertools.product(x_list, y_list))
        # remove the current patch position from neighbours
        neighbours.remove(t)
        score = 0
        for neighbour in neighbours:
            score += 1 * (neighbour in all_values)
        return score

    def patches_filtering_with_neighbours(self, plot=False, neighbours_level=3):
        """Filter patches based on the number of neighbors it has. Useful for not clean batch like human melanoma. Only applied for human melanoma dataset for the moment.

        Args:
            plot (bool, optional): Whether or not plotting the patches location colored in a different color if the patch is discard. Defaults to False.
            neighbours_level (int, optional): The max distance level between the current patch and its neighbours. Defaults to 3.
        """
        logging.info("Start filtering patches...")

        patches_info_df = pd.DataFrame.from_records(self.all_patches_info)
        shape_before = patches_info_df.shape

        # Do not take data from tcga database for this filtering step
        data = patches_info_df[~(patches_info_df["batch"] == "Tcga")]

        data["center_height_origin"] = data.apply(lambda r: r["start_height_origin"] + int(r["shape_pixel"] / 2), axis=1)
        data["center_width_origin"] = data.apply(lambda r: r["start_width_origin"] + int(r["shape_pixel"] / 2), axis=1)

        data["top_left_h_w"] = data.apply(lambda row: (row["start_height_origin"], row["start_width_origin"]), axis=1)

        groups = data.groupby(by="name_origin")
        patches_name_to_discard = []
        for name, group in groups:
            logging.info("\nImage {}".format(name))
            group["neighbour_score"] = group["top_left_h_w"].apply(
                lambda t: self.check_neighbours(
                    t=t,
                    all_values=list(group["top_left_h_w"].values),
                    neighbours_level=neighbours_level,
                )
            )
            group["outlier"] = 1 * (group["neighbour_score"] < 4 * neighbours_level)
            logging.info("Number of outliers = {}".format(group["outlier"].sum()))
            logging.info("Percentage of outliers = {} %".format(group["outlier"].sum() / len(group) * 100))

            patches_name_to_discard.extend(list(group[group["outlier"] == True]["name"].values))

            if plot:
                group.plot.scatter(x="center_width_origin", y="center_height_origin", c="outlier", cmap="viridis")
                plt.title(name)
                plt.gca().invert_yaxis()
                plt.show()

            for path in group[group["outlier"] == True]["path"].values:
                logging.info("Remove patch {}".format(path))
                os.remove(path)

            patches_info_df = patches_info_df[patches_info_df["name"].apply(lambda n: n not in patches_name_to_discard)]

        logging.info("Shape patches_info_df before = {} and after {}".format(shape_before, patches_info_df.shape))
        self.all_patches_info = patches_info_df.to_dict(orient="records")

    @staticmethod
    def create_circular_mask(h, w, center=None, radius=None):
        """Create a circular mask around center=center of radius=radius. Outside the circle = 0 and inside = 1.

        Args:
            h (int): height of the image
            w (int): width of the image
            center (float, optional): Center of the mask. Defaults to None.
            radius (float, optional): Radius of the mask. Defaults to None.

        Returns:
            np.array: mask with outside the circle = 0 and inside = 1.
        """
        if center is None:
            # the middle of the image
            center = (int(w / 2), int(h / 2))
        if radius is None:
            # the smallest distance between the center and image walls
            radius = min(center[0], center[1], w - center[0], h - center[1])
        Y, X = np.ogrid[:h, :w]
        dist_from_center = np.sqrt((X - center[0]) ** 2 + (Y - center[1]) ** 2)
        mask = dist_from_center <= radius
        return mask

    def spots_patches(self, img_filename, spots_df,hdf5_file=None):
        """Compute patches on the spatial transcriptomics spots.

        Args:
            img_filename (str): Image file
            spots_df (pd.DataFrame): Dataframe with the spots location
            hdf5_file (h5py.File, optional): If you want to save the patches in a hdf5 file. Defaults to None.
        """
        logging.info("Start computing spots patches for image: {}".format(img_filename))

        if self.spot_radius is None:
            raise Exception("spot_radius must be given to compute the spot patches!")

        img = np.array(Image.open(img_filename).convert("RGB"))
        patches_info_list = []

        for index, row in spots_df.iterrows():
            start_x = int(row["pixel_x"] - self.spot_radius)
            start_y = int(row["pixel_y"] - self.spot_radius)
            patch = img[start_y : start_y + 2 * self.spot_radius, start_x : start_x + 2 * self.spot_radius, :]

            # Create a spot mask on the patch
            if self.spot_mask:
                mask = self.create_circular_mask(patch.shape[0], patch.shape[1], radius=self.spot_radius)
                patch = np.expand_dims(mask, axis=2) * patch

            origin_name = img_filename.split("/")[-1].split(
                "." + img_filename.split(".")[-1]
            )[0]
            if spots_df.columns.str.contains("source_file").any():
                name = (
                    spots_df.loc[index, "source_file"]
                    + "_spot{}x{}".format(int(row["x"]), int(row["y"]))
                )
            else:
                name = (
                    origin_name[0]
                    + "_rep"
                    + origin_name[1]
                    + "_spot{}x{}".format(int(row["x"]), int(row["y"]))
                )
                
            print(f"Name = {name}", flush=True)
            
            saving_path = os.path.join(self.saving_folder, name + "." + self.extension)
            if hdf5_file is None:
                Image.fromarray(patch).save(saving_path)
                logging.debug("Patch saved to {}".format(os.path.abspath(saving_path)))
            else:
                hdf5_file.create_dataset(name, data=patch)
                logging.debug("Patch saved to hdf5 file")

            patch_data = {
                "path": os.path.abspath(saving_path),
                "name": name,
                "mpp_height": None,
                "mpp_width": None,
                "shape_micron": None,
                "shape_pixel": patch.shape[0],
                "overlap_pixel": None,
                "path_origin": img_filename,
                "name_origin": origin_name,
                "shape_origin": img.shape,
                "start_height_origin": start_y,
                "start_width_origin": start_x,
                "extension_origin": img_filename.split(".")[-1],
                "batch": img_filename.split("/data/")[-1].split("/")[0].replace("_", " ").capitalize(),
                "tumor": name.split("_spot")[0].split(".")[0].split("_rep")[0],
                "mean_intensity": np.mean(patch),
                "median_intensity": np.median(patch),
                "std_intensity": np.std(patch),
                "entropy_intensity": shannon_entropy(patch),
                "spots_info": {
                    "y": int(row["y"]),
                    "x": int(row["x"]),
                    "radius_pixel": self.spot_radius,
                },
            }
            patches_info_list.append(patch_data)

            patch = None
        logging.info("Finish saving {} spots patches".format(len(patches_info_list)))

        self.all_patches_info.extend(patches_info_list)

    def compute_all_patches(self):
        """Compute all patches from image files (random and spot patches)."""
        hdf5_file = h5py.File(os.path.join(self.saving_folder, "patches.hdf5"), "w")
        if self.images_filenames_random_patches is not None:

            for file in self.images_filenames_random_patches:
                print(f"Computing random patches for file {file}", flush=True)
                try:
                    patches_info = self.random_patches(
                        file,
                        debug_jupyter=False,
                        benchmark_count=None,
                        hdf5_file=hdf5_file,
                    )
                    self.all_patches_info.extend(patches_info)
                    print(
                        "Number of patches computed: {}".format(
                            len(self.all_patches_info)
                        )
                    )
                except Exception as e:
                    logging.error(
                        "Problem with image {}, no patch can be computed: {}".format(
                            file, e
                        )
                    )

            

            # Discard some patches with low neigbours count
            if self.patches_number is None and self.filter_with_neighbours:
                self.patches_filtering_with_neighbours()

        # Add the spots patches
        if self.images_filenames_spot_patches is not None and self.spots_filenames is not None:
            assert len(self.images_filenames_spot_patches) == len(
                self.spots_filenames
            ), "images_filenames_spot_patches and spots_filenames does not have the same size"
            for img_file, spots_file in zip(sorted(self.images_filenames_spot_patches), sorted(self.spots_filenames)):
                spots_df = pd.read_csv(spots_file, sep="\t", compression="gzip")
                self.spots_patches(img_file, spots_df, hdf5_file=hdf5_file)
                
        elif (
            self.images_filenames_spot_patches is not None and self.spots_df is not None
        ):
            assert len(self.images_filenames_spot_patches) == len(
                self.spots_df
            ), "images_filenames_spot_patches and spots_filenames does not have the same size"
            for img_file, spots_df in zip(
                self.images_filenames_spot_patches, self.spots_df
            ):
                self.spots_patches(img_file, spots_df, hdf5_file=hdf5_file)

        try:
            saving_path = os.path.join(self.saving_folder, self.name_patches_info_file)
            with gzip.open(saving_path, "wb") as file:
                pickle.dump(self.all_patches_info, file)
            print("Saving ok")
            logging.info("Saved patches info pkl file to {}".format(os.path.abspath(saving_path)))
        except Exception as e:
            logging.error("Problem with patches information pkl file saving: {}".format(e))
        hdf5_file.close()

    @staticmethod
    def get_random_balanced_benchmark_files():
        """Get random WSI with a balanced number from the human melanoma, HER2-positive breast cancer and TCGA dataset to create the benchmark dataset."""
        hm = sorted(glob.glob("../data/human_melanoma/raw_images/*.jpg"))
        tcga = sorted(glob.glob("../data/tcga/*/*.svs"))
        # filtering patient A
        her2 = sorted(glob.glob("../data/HER2_breast_cancer/images/HE/*.jpg"))[6:]
        files_number = int(min([len(hm), len(tcga), len(her2)]))
        benchmark_list = random.sample(hm, files_number)
        benchmark_list.extend(random.sample(tcga, files_number))
        benchmark_list.extend(random.sample(her2, files_number))
        print("Benchmark files list: {}".format(benchmark_list))
        return benchmark_list

    def benchmark_dataset(self):
        """Compute the benchmark patches dataset."""
        benchmark_count = 0
        patch_size_micron_list = [50, 75, 100, 125, 150, 175, 200, 225, 250]
        self.patches_number = 100
        self.overlap_pixel = 0
        self.extension = "tiff"
        self.filter_background = True
        all_patches_info = []
        for patch_size_micron in patch_size_micron_list:
            images_filenames = self.get_random_balanced_benchmark_files()
            for file in images_filenames:
                try:
                    self.patch_size_micron = patch_size_micron
                    patches_info, benchmark_count = self.random_patches(
                        file,
                        benchmark_count=benchmark_count,
                    )
                    all_patches_info.extend(patches_info)
                except Exception as e:
                    logging.error("Problem with image {}, no patch can be computed: {}".format(file, e))

        try:
            saving_path = os.path.join(self.saving_folder, self.name_patches_info_file)
            with gzip.open(saving_path, "wb") as file:
                pickle.dump(all_patches_info, file)
            print("Saving ok")
            logging.info("Saved patches info pkl file to {}".format(os.path.abspath(saving_path)))
        except Exception as e:
            logging.error("Problem with patches information pkl file saving: {}".format(e))

    def entropy_intensity_patches_threshold(
        self,
        patches_info_df=None,
        quantile_threshold=0.99,
        intersection=True,
        density=False,
        bins=50,
        bayesian_gmm=True,
        scaling_to_max=True,
        save_fig_name=None,
    ):
        """It finds the threshold of Shannon entropy that corresponds to the quantile threshold to distinguish foreground from background patches.
        It models two gaussian distributions from data coming from both foreground and background patches and applies the quantile threshold on
        the gaussian modeling background patches. It plots the histograms with the two gaussian distributions, the threshold and the intersection.

        Args:
            patches_info_df (pd.DataFrame): Dataframe with all the patches information from both foreground and background that is computed from the .pkl that outputs the compute_all_patches function. Defaults to None.
            quantile_threshold (float, optional): Quantile threshold for the gaussian that model the background patches distribution. Defaults to 0.995.
            intersection (bool, optional): To plot and print the intersection between the two gaussian distribution. Defaults to True.
            density (bool, optional): To have the histogram in density instead of count. Defaults to False.
            scaling_to_max (bool, optional): True to scale the gaussian to the max value between [mean-std;mean+std]. False to scale the gaussian to the mean value ignoring 0 between [mean-std;mean+std]. Defaults False.
            save_fig_name (str, optional): Name for the figure that will be saved, pdf format preferred

        Returns:
            threshold (float): Shannon entropy threshold that corresponds to the quantile_threshold of the background patches gaussian
            labels (list): List of labels for each patch of the patches_info_df (1=background, 2=foreground)
        """
        if patches_info_df is None:
            if self.all_patches_info is None:
                try:
                    path = os.path.join(self.saving_folder, self.name_patches_info_file)
                    with gzip.open(path) as file:
                        self.all_patches_info = pickle.load(file)
                except Exception as e:
                    logging.error("Problem to load the {}: {}".format(path, e))
            patches_info_df = pd.DataFrame.from_records(self.all_patches_info).sort_values(by="name").reset_index(drop=True)

        if bayesian_gmm:
            gmm = BayesianGaussianMixture(n_components=2, random_state=0, max_iter=500)
        else:
            gmm = GaussianMixture(n_components=2, random_state=0, max_iter=500)
        gmm.fit(
            np.reshape(patches_info_df["entropy_intensity"].values, (patches_info_df["entropy_intensity"].values.shape[0], 1))
        )

        fig, ax = plt.subplots(figsize=(6, 5))

        y_hist, x_hist, _ = plt.hist(
            patches_info_df["entropy_intensity"].values, density=density, bins=bins, alpha=0.5, color="gray"
        )
        x_hist = (x_hist[1:] + x_hist[:-1]) / 2
        hist_df = pd.DataFrame(np.array([y_hist, x_hist]).T, columns=["y", "x"])

        x = np.linspace(patches_info_df["entropy_intensity"].min(), patches_info_df["entropy_intensity"].max(), 1000).reshape(
            -1, 1
        )

        if gmm.means_[0][0] > gmm.means_[1][0]:
            ind_list = [1, 0]
        else:
            ind_list = [0, 1]

        pdfs = {}
        for i, ind in enumerate(ind_list):
            mean = gmm.means_[ind][0]
            covariance = gmm.covariances_[ind][0, 0]
            dist = sps.norm(loc=mean, scale=np.sqrt(covariance))
            pdf = dist.pdf(x)

            if i == 0:
                threshold = dist.ppf(quantile_threshold)
                color = "orange"
            else:
                color = "green"

            gauss_max = np.max(pdf)
            hist_max = hist_df[
                np.logical_and(hist_df["x"] > mean - np.sqrt(covariance), hist_df["x"] < mean + np.sqrt(covariance))
            ]["y"]
            if scaling_to_max:
                hist_max = hist_max.max()
            else:
                hist_max = np.nanmean(hist_max.apply(lambda y: float("nan") if y == 0 else y))
            scaling_factor = hist_max / gauss_max
            plt.plot(x, pdf * scaling_factor, label=f"Component {i + 1}", color=color)
            pdfs[f"Component {i + 1}"] = np.reshape(pdf * scaling_factor, -1)

        plt.axvline(threshold, color="red", label=f"{quantile_threshold*100:.0f}th: {threshold:.3f}")
        print("Threshold based on the {}th quantile: {}".format(quantile_threshold * 100, np.round(threshold, 3)))

        if intersection:
            intersection_ind = np.argwhere(np.diff(np.sign(pdfs["Component 1"] - pdfs["Component 2"]))).flatten()
            plt.scatter(
                x[intersection_ind],
                pdfs["Component 1"][intersection_ind],
                c="red",
                label="Intersection: {}".format(np.round(x[intersection_ind].flatten().max(), 3)),
                marker="X",
            )
            print("Intersection x point: {}".format(np.round(x[intersection_ind][0, 0], 3)))

        labels = patches_info_df["entropy_intensity"].apply(lambda e: 2 if e > threshold else 1)
        print("Percentage of patches in component 2: {}%".format(np.round((labels == 2).sum() / len(labels) * 100), 3))

        plt.xlabel("Shannon entropy")
        if density:
            plt.ylabel("Patches density")
        else:
            plt.ylabel("Patches count")
        plt.title("Bimodal fitting of the Shannon entropy of patches")

        handles, previous_labels = ax.get_legend_handles_labels()
        new_labels = ["Background", "Foreground"]
        new_labels.extend(previous_labels[2:])
        plt.legend(handles=handles, labels=new_labels)

        if save_fig_name is not None:
            plt.savefig(os.path.join(SAVE_FIG_FOLDER, save_fig_name), bbox_inches="tight")

        fig.show()

        return threshold, labels.to_list()

    def plot_patches_location(self, patches_info_df=None):
        """Subplot of a scatter of the patches location and the full image.

        Args:
            patches_info_df (pd.DataFrame): Dataframe containing the patches information that is computed from the .pkl that outputs the compute_all_patches function. Defaults to None.
        """
        if patches_info_df is None:
            if self.all_patches_info is None:
                try:
                    path = os.path.join(self.saving_folder, self.name_patches_info_file)
                    with gzip.open(path) as file:
                        self.all_patches_info = pickle.load(file)
                except Exception as e:
                    logging.error("Problem to load the {}: {}".format(path, e))
            patches_info_df = pd.DataFrame.from_records(self.all_patches_info).sort_values(by="name").reset_index(drop=True)

        patches_info_df["center_height_origin"] = patches_info_df.apply(
            lambda r: r["start_height_origin"] + int(r["shape_pixel"] / 2), axis=1
        )
        patches_info_df["center_width_origin"] = patches_info_df.apply(
            lambda r: r["start_width_origin"] + int(r["shape_pixel"] / 2), axis=1
        )

        groups = patches_info_df.groupby(by="name_origin")

        for name, group in groups:
            fig, ax = plt.subplots(1, 2, figsize=(11, 5), sharex=False, sharey=False, layout="constrained")

            group.plot.scatter(
                x="center_width_origin",
                y="center_height_origin",
                c="entropy_intensity",
                alpha=0.5,
                title="Patches position colored by entropy for image {}".format(group["name_origin"].iloc[0].split(".")[0]),
                ax=ax[0],
            )
            ax[0].set_xlim(0, group["shape_origin"].iloc[0][1])
            ax[0].set_ylim(0, group["shape_origin"].iloc[0][0])
            ax[0].invert_yaxis()

            if "path_origin" in group.columns:
                path_origin = group["path_origin"].iloc[0]
            else:
                path_origin = glob.glob("../data/HER2_breast_cancer/images/HE/{}.jpg".format(group["name_origin"].iloc[0]))
                path_origin.extend(glob.glob("../data/human_melanoma/raw_images/{}.jpg".format(group["name_origin"].iloc[0])))
                path_origin.extend(glob.glob("../data/tcga/*/{}*.svs".format(group["name_origin"].iloc[0])))
                assert len(path_origin) == 1
                path_origin = path_origin[0]

            slide = openslide.open_slide(path_origin)
            ax[1].imshow(slide.get_thumbnail((600, 600)), aspect="auto")
            ax[1].set_title("H&E image {}".format(group["name_origin"].iloc[0].split(".")[0]))
            ax[1].axis("off")

            plt.show()

    def create_subset_hdf5(self, original_hdf5_path, csv_path, new_hdf5_path):
        """
        Create a new HDF5 file containing only the patches listed in a CSV file.

        Args:
            original_hdf5_path (str): Path to the original HDF5 file.
            csv_path (str): Path to the CSV file containing the names of the patches to include.
            new_hdf5_path (str): Path to save the new HDF5 file.
        """
        import h5py
        import pandas as pd

        # Load the patch names from the CSV file
        patch_names = pd.read_csv(csv_path, index_col=0).index.tolist()

        if not os.path.exists(os.path.dirname(new_hdf5_path)):
            os.makedirs(os.path.dirname(new_hdf5_path))

        # Open the original HDF5 file and create a new HDF5 file
        with h5py.File(original_hdf5_path, "r") as original_hdf5, h5py.File(new_hdf5_path, "w") as new_hdf5:
            for patch_name in patch_names:
                if patch_name in original_hdf5:
                    # Copy the dataset from the original HDF5 to the new HDF5
                    original_hdf5.copy(patch_name, new_hdf5)
                else:
                    logging.warning(f"Patch {patch_name} not found in the original HDF5 file.")
        logging.info(f"New HDF5 file created at {new_hdf5_path} with selected patches.")
