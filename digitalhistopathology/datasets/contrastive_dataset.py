#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import numpy as np
from PIL import Image
from torch.utils.data import Dataset

import matplotlib.pyplot as plt
import torchvision.transforms as transforms

class ContrastiveDataset(Dataset):
    """Dataset class for retraining."""

    def __init__(
        self,
        image_files_list,
        debug=False,
        random_resized_crop=False,
        random_resized_scale=[0.01, 1.0],
        random_flip=True,
        random_rotate=True,
        random_grayscale=0.4,
        color_jitter=[0.4, 0.4, 0.4, 0.2],
        gaussian_blur=False,
        ratio_of_dataset=1,
    ):

        self.random_resized_crop = random_resized_crop
        self.random_resized_scale = random_resized_scale
        self.random_flip = random_flip
        self.random_rotate = random_rotate
        self.random_grayscale = random_grayscale
        self.color_jitter = color_jitter
        self.gaussian_blur = gaussian_blur
        self.image_files_list = np.random.choice(
            image_files_list,
            int(len(image_files_list) * ratio_of_dataset),
            replace=False,
        ).tolist()
        self.debug = debug
        self.transform = self.get_data_transformation(image_files_list)

    def __len__(self):
        return len(self.image_files_list)

    def __getitem__(self, idx):
        path = self.image_files_list[idx]
        image = Image.open(path).convert("RGB")
        sample = dict()
        sample["aug_copy_1"] = self.transform(image)
        sample["aug_copy_2"] = self.transform(image)

        if self.debug:
            fig, ax = plt.subplots(ncols=3, figsize=(15, 5))
            ax[0].imshow(image)
            ax[1].imshow(sample["aug_copy_1"].permute(1, 2, 0).numpy())
            ax[2].imshow(sample["aug_copy_2"].permute(1, 2, 0).numpy())
            fig.show()
        return sample

    def name(self, idx):
        file_format = self.image_files_list[idx].split(".")[-1]
        return self.image_files_list[idx].split("/")[-1].split("." + file_format)[0]

    def get_data_transformation(self, size=224):
        """Return a set of data augmentation transformations as described in the SimCLR paper. Inspired from:
        https://github.com/sthalles/SimCLR/blob/master/data_aug/contrastive_learning_dataset.py

        Args:
            size (int, optional): Target image size. Defaults to 224.

        Returns:
            torch.transforms: Tranformation for the dataset.
        """
        # TODO: normalization ?? No because torchvision.transforms.Normalize(mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225)) need to got tensor
        transform_list = []
        if self.random_resized_crop:
            transform_list.append(
                transforms.RandomResizedCrop(
                    size=size, scale=tuple(self.random_resized_scale)
                )
            )
        else:
            transform_list.append(transforms.Resize(224, antialias=True))
        if self.random_flip:
            transform_list.append(transforms.RandomHorizontalFlip())
            transform_list.append(transforms.RandomVerticalFlip())
        color_jitter = transforms.ColorJitter(
            self.color_jitter[0],
            self.color_jitter[1],
            self.color_jitter[2],
            self.color_jitter[3],
        )
        transform_list.extend(
            [
                transforms.RandomApply([color_jitter], p=0.8),
                transforms.RandomGrayscale(p=self.random_grayscale),
            ]
        )
        if self.random_rotate:
            transform_list.append(
                transforms.RandomApply([transforms.RandomRotation((90, 90))], p=0.5)
            )
        if self.gaussian_blur:
            # round a float up to next odd integer
            kernel_size = np.ceil(0.1 * size) // 2 * 2 + 1
            gaussian_blur = transforms.GaussianBlur(kernel_size=kernel_size)
            transform_list.append(transforms.RandomApply([gaussian_blur], p=0.5))
        transform_list.append(transforms.ToTensor())
        data_transforms = transforms.Compose(transform_list)
        return data_transforms

