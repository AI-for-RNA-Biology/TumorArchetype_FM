#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import torch
import torchvision
import pandas as pd
from PIL import Image
import h5py



class ImageEmbeddingDataset(torch.utils.data.Dataset):
    def __init__(self,  model_name, images_filenames_pd=None, hdf5_file=None):
        """ImageEmbeddingDataset is a class used to load patches to be then evaluated 
        by a pretrained model.

        Args:
            images_filenames_pd (pd.DataFrame): Dataframe with one column named filename that contains all the patches filenames we want to evaluate.
        """
        super().__init__()
        # Transform from https://github.com/Xiyue-Wang/TransPath/blob/main/get_features_CTransPath.py
        # TODO: Test Vahadane’s color normalization ?
        if model_name == "CTransPath":
            self.transform = torchvision.transforms.Compose(
                [
                    torchvision.transforms.Resize(224, antialias=True),
                    torchvision.transforms.ToTensor(),
                    torchvision.transforms.Normalize(mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225)),
                ]
            )
        elif model_name == "SimCLR" or "simclr" in model_name:  # for retrained models
            # TODO: Normalization also for SimCLR ??? See https://github.com/ozanciga/self-supervised-histopathology/issues/2
            self.transform = torchvision.transforms.Compose(
                [
                    torchvision.transforms.Resize(224, antialias=True),
                    torchvision.transforms.ToTensor(),
                ]
            )
        elif model_name == "UNI":
            self.transform = torchvision.transforms.Compose(
                [
                    torchvision.transforms.Resize(224),
                    torchvision.transforms.ToTensor(),
                    torchvision.transforms.Normalize(mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225)),
                ]
            )
        elif model_name == "Virchow":
            # TODO: Normalization also for SimCLR ??? See https://github.com/ozanciga/self-supervised-histopathology/issues/2
            self.transform = torchvision.transforms.Compose(
                [
                    torchvision.transforms.Resize(224, antialias=True),
                    torchvision.transforms.ToTensor(),
                ]
            )
        elif model_name == "ProvGigaPath":
            self.transform = torchvision.transforms.Compose(
                [
                    torchvision.transforms.Resize(224, interpolation=torchvision.transforms.InterpolationMode.BICUBIC),
                    ##torchvision.transforms.CenterCrop(224) done if first resizing to 256
                    torchvision.transforms.ToTensor(),
                    torchvision.transforms.Normalize(mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225)),
                ]
            )
        elif model_name == "vit":
            self.transform = torchvision.transforms.Compose(
                [
                    torchvision.transforms.Resize(224),
                    torchvision.transforms.ToTensor(),
                    torchvision.transforms.Normalize(
                        mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225)
                    ),
                ]
            )
        elif model_name == "UNI2":
            self.transform = torchvision.transforms.Compose(
                [
                    torchvision.transforms.Resize(224),
                    torchvision.transforms.ToTensor(),
                    torchvision.transforms.Normalize(
                        mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225)
                    ),
                ]
            )
        else:
            raise Exception("Not a correct model name, see model.py: {}".format(model_name))

        self.images_filenames_pd = images_filenames_pd
        self.hdf5_file = hdf5_file

        if images_filenames_pd is not None:
            self.format = self.images_filenames_pd.filename[0].split(".")[-1]
        else:
            self.format = "hdf5"
            self.names = self.get_names_from_hdf5()            

    def __len__(self):
        if self.images_filenames_pd is not None:
            return len(self.images_filenames_pd)
        else:
            return len(self.names)

    def __getitem__(self, idx):
        if self.images_filenames_pd is not None:
            path = self.images_filenames_pd.filename[idx]
            image = Image.open(path).convert("RGB")
        else:
            with h5py.File(self.hdf5_file, "r") as f:
                image = f[self.names[idx]][:]
                image = image[:]
                image = Image.fromarray(image).convert("RGB")
                
        image = self.transform(image)
        
        return image

    def name(self, idx):
        if self.images_filenames_pd is not None:
            return self.images_filenames_pd.filename[idx].split("/")[-1].split("." + self.format)[0]    
        else:
            return self.names[idx]
            
    def get_names_from_hdf5(self):
        print(f"Loading names from hdf5 file: {self.hdf5_file}", flush=True)
        with h5py.File(self.hdf5_file, "r") as f:
            dataset_names = [] 
            def collect_datasets(name, obj): 
                if isinstance(obj, h5py.Dataset): 
                    dataset_names.append(name) 
            # Visit all items and collect dataset names
            f.visititems(collect_datasets) 
            # Count the number of datasets 
            return dataset_names
