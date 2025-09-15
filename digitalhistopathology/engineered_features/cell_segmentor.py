#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import pickle
import gzip
import pyvips
import os

import sys
import subprocess
import glob
from openslide import OpenSlide
sys.path.append("../")
sys.path.append("../CellViT")

from pathlib import Path
import matplotlib.pyplot as plt
import torch
from multiprocessing import Manager

print("No problem with the imports")



class CellSegmentor:
    def __init__(self,
                 path_to_wsis=None,
                 list_wsi_filenames=None,
                 patches_info_filename=None,
                 results_saving_folder="../results",
                 temporary_folder="../results/tmp",
                 dataset_name="dataset",
                 mpp=None,
                 magnification=None,
                 ):

        self.path_to_wsis = path_to_wsis
        self.patches_info = None

        if patches_info_filename is not None:
            with gzip.open(patches_info_filename) as file:
                self.patches_info = pickle.load(file)
            self.list_wsi_filenames = list(set([patch['path_origin'] for patch in self.patches_info]))
        elif path_to_wsis is not None:
            self.list_wsi_filenames = glob.glob(os.path.join(path_to_wsis, '*'))
            self.list_wsi_filenames = [path_to_wsi for path_to_wsi in self.list_wsi_filenames if os.path.isfile(path_to_wsi)]
        elif list_wsi_filenames is not None:
            self.list_wsi_filenames = list_wsi_filenames

        else:
            raise ValueError("path_to_wsis, list_wsi_filenames and patches_info_filename cannot be all None")

        self.results_saving_folder = results_saving_folder
        self.temporary_folder = temporary_folder
        self.dataset_name = dataset_name
        self.mpp = mpp
        self.magnification = magnification

        

        if not os.path.exists(self.results_saving_folder):
            os.makedirs(self.results_saving_folder)

        if not os.path.exists(self.temporary_folder):
            os.makedirs(self.temporary_folder)
    

            
    def segment(self, image):
        raise NotImplementedError("This method should be overridden by subclasses")
    

    def check_and_convert_wsifile(self, path_to_wsi):
        extension = os.path.basename(path_to_wsi).split('.')[-1]
        filename = os.path.basename(path_to_wsi).split(f'.{extension}')[0]

        print(f"Checking and converting WSI: {filename} with extension: {extension}")
        if extension in ["svs",
                         "tiff",
                         "tif",
                         "bif",
                         "scn",
                         "ndpi",
                         "vms",
                         "vmu"]:
            return path_to_wsi
        elif extension == 'jpg':
            image = pyvips.Image.new_from_file(path_to_wsi, access='sequential')
            path_new_tif = os.path.join(self.temporary_folder, "wsi", self.dataset_name, f"{filename}.tif")

            if not os.path.exists(os.path.dirname(path_new_tif)):
                os.makedirs(os.path.dirname(path_new_tif))

            image.tiffsave(path_new_tif,
                           compression="jpeg",
                           Q=75,
                           tile=True,
                           tile_width=256,
                           tile_height=256,
                           pyramid=True)
            return path_new_tif
        else:
            raise ValueError("The WSI file extension is not supported.")
    
    def get_all_magnifications(self):
        magnifications = []
        for path_to_wsi in self.list_wsi_filenames:
            slide = OpenSlide(str(path_to_wsi))
            if "openslide.objective-power" in slide.properties:
                slide_mag = float(slide.properties.get("openslide.objective-power"))
            elif self.magnification is not None:
                slide_mag = self.magnification
            else:
                raise NotImplementedError(
                    "MPP must be defined either by metadata or provided as an argument!")
            magnifications.append(slide_mag)
        self.magnifications = magnifications
    
    def get_all_mpps(self):


        mpps = []
        for path_to_wsi in self.list_wsi_filenames:
            print(f"Getting MPP for slide: {path_to_wsi}", flush=True)
            slide = OpenSlide(str(path_to_wsi))
            if "openslide.mpp-x" in slide.properties:
                slide_mpp = float(slide.properties.get("openslide.mpp-x"))
            elif self.mpp is not None:
                slide_mpp = self.mpp
            else:
                raise NotImplementedError(
                    "MPP must be defined either by metadata or provided as an argument!")
            mpps.append(slide_mpp)
        self.mpps = mpps

    def get_all_extensions(self):
        extensions = []
        for path_to_wsi in self.list_wsi_filenames:
            extension = os.path.basename(path_to_wsi).split('.')[-1]
            extensions.append(extension)
        self.extensions = extensions
    

class CellVIT_Segmentor(CellSegmentor):
    def __init__(self,
                 mpp=None, 
                 magnification=None,
                 path_to_wsis=None,                  
                 list_wsi_filenames=None,
                 patches_info_filename=None,
                 results_saving_folder="../results",
                 temporary_folder="../results/tmp",
                 model_path="../CellViT/models/pretrained",
                 gpu=0,
                 enforce_mixed_precision=False,
                 batch_size=8,
                 dataset_name="dataset"
                 ):
        super().__init__(path_to_wsis=path_to_wsis,
                         list_wsi_filenames=list_wsi_filenames,
                         patches_info_filename=patches_info_filename,
                         results_saving_folder=results_saving_folder,
                         temporary_folder=temporary_folder,
                         dataset_name=dataset_name,
                         mpp=mpp,
                         magnification=magnification)
        
        self.list_wsi_filenames = [self.check_and_convert_wsifile(path_to_wsi) for path_to_wsi in self.list_wsi_filenames]

        self.get_all_mpps()
        self.get_all_magnifications()
        self.get_all_extensions()

        self.patch_size = 1024
        self.patch_overlap = 6.25

        self.model_path = model_path
        self.gpu = gpu
        self.enforce_mixed_precision = enforce_mixed_precision
        self.batch_size = batch_size


    
    
    
    def cellvit_preprocess(self):
        from CellViT.preprocessing.patch_extraction.src.cli import PreProcessingConfig
        from CellViT.preprocessing.patch_extraction.src.patch_extraction import PreProcessor


        for path_to_wsi, mpp, magnification, extension in zip(self.list_wsi_filenames, self.mpps, self.magnifications, self.extensions):
            configuration = PreProcessingConfig(wsi_paths=path_to_wsi,
                                                output_path=os.path.join(self.results_saving_folder, 
                                                                             "segmentation/CellViT", self.dataset_name),
                                                patch_size=self.patch_size,
                                                patch_overlap=self.patch_overlap,
                                                wsi_properties={'slide_mpp': mpp,
                                                                'magnification': magnification},
                                                wsi_extension=extension)
            
            slide_processor = PreProcessor(slide_processor_config=configuration)
            slide_processor.sample_patches_dataset()


    def segment(self):
        # Implement the segmentation logic here
        from CellViT.cell_segmentation.inference.cell_detection import CellSegmentationInference
        from CellViT.cell_segmentation.inference.cell_detection import check_wsi
        from CellViT.datamodel.wsi_datamodel import WSI


        
        for path_to_wsi, magnification in zip(self.list_wsi_filenames, self.magnifications):
            if magnification not in [20, 40]:
                print(f"Skipping WSI {path_to_wsi} with magnification {magnification}. Only magnifications 20x and 40x are supported.")
                continue
            
            print(f"Segmenting WSI {path_to_wsi} with magnification {magnification} using model CellViT-SAM-H-x{magnification}...")
            cell_segmentation = CellSegmentationInference(model_path=os.path.join(self.model_path, f"CellViT-SAM-H-x{int(magnification)}.pth"),
                                                          gpu=self.gpu,
                                                          enforce_mixed_precision=self.enforce_mixed_precision)
            
            cell_segmentation.logger.info("Processing single WSI file")
            wsi_path = Path(path_to_wsi)
            wsi_name = wsi_path.stem
            wsi_file = WSI(
                name=wsi_name,
                patient=wsi_name,
                slide_path=wsi_path,
                patched_slide_path=os.path.join(self.results_saving_folder, "segmentation/CellViT", self.dataset_name, wsi_name),
            )
            check_wsi(wsi=wsi_file, magnification=magnification)
            cell_segmentation.process_wsi(
                wsi=wsi_file,
                subdir_name=None,
                geojson=True,
                batch_size=self.batch_size,
            )    


def cellVIT_segmentation(mpp=None, 
                         magnification=None,
                         path_to_wsis=None,      
                         list_wsi_filenames=None,
                         patches_info_filename=None,
                         results_saving_folder="../results",
                         temporary_folder="../results/tmp",
                         model_path="../CellViT/models/pretrained",
                         gpu=0,
                         enforce_mixed_precision=False,
                         dataset_name="dataset"):
    
    seg = CellVIT_Segmentor(mpp=mpp,
                            magnification=magnification,
                            path_to_wsis=path_to_wsis,
                            list_wsi_filenames=list_wsi_filenames,
                            patches_info_filename=patches_info_filename,
                            results_saving_folder=results_saving_folder,
                            temporary_folder=temporary_folder,
                            model_path=model_path,
                            gpu=gpu,
                            enforce_mixed_precision=enforce_mixed_precision,
                            dataset_name=dataset_name)
    
    seg.cellvit_preprocess()
    seg.segment()



def main_segmentation():
    import argparse

    parser = argparse.ArgumentParser(description='Segmentation')
    parser.add_argument('--segmentation_mode', type=str, default='cellvit', help='Which model to use for segmentation')
    parser.add_argument('--mpp', type=float, default=None, help='Microns per pixel')
    parser.add_argument('--magnification', type=float, default=None, help='Magnification of the WSI')
    parser.add_argument('--path_to_wsis', type=str, default=None, help='Path to the WSIs')
    parser.add_argument('--list_wsi_filenames', type=str, nargs='+', default=None, help='List of WSI filenames')
    parser.add_argument('--patches_info_filename', type=str, default=None, help='Patches info filename')
    parser.add_argument('--results_saving_folder', type=str, default="../results", help='Results saving folder')
    parser.add_argument('--temporary_folder', type=str, default="../results/tmp", help='Temporary folder')
    parser.add_argument('--model_path', type=str, default="../CellViT/models/pretrained/", help='Model path')
    parser.add_argument('--gpu', type=int, default=0, help='GPU')
    parser.add_argument('--enforce_mixed_precision', type=bool, default=False, help='Enforce mixed precision')
    parser.add_argument('--dataset_name', type=str, default="dataset", help='Dataset name')
    
    args = parser.parse_args()
    
    if args.segmentation_mode == 'cellvit':
        cellVIT_segmentation(mpp=args.mpp, 
                            magnification=args.magnification,
                            path_to_wsis=args.path_to_wsis,      
                            list_wsi_filenames=args.list_wsi_filenames,
                            patches_info_filename=args.patches_info_filename,
                            results_saving_folder=args.results_saving_folder,
                            temporary_folder=args.temporary_folder,
                            model_path=args.model_path,
                            gpu=args.gpu,
                            enforce_mixed_precision=args.enforce_mixed_precision,
                            dataset_name=args.dataset_name)

    else:
        raise ValueError("Segmentation mode not supported. Please choose among 'cellvit' and 'semantic'")

if __name__ == "__main__":

    main_segmentation()






    

