#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#

from ast import mod
import math
import os
import pickle
from pyexpat import model
import re
import stat
from cycler import V
import yaml
import torch
from torch.utils.data import DataLoader
from tkinter import TRUE
from tracemalloc import stop
from arrow import get
from scipy import optimize
from sympy import Q
import timm.utils
import torch.cuda.amp as amp
import matplotlib.pyplot as plt
import numpy as np
import timm
import torch
import torch.nn as nn
import torchvision
from huggingface_hub import login
import warnings
import torch.nn.utils as utils
import logging
import sys
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
dinov2_path = os.path.join(project_root, "dinov2")
sys.path.append(dinov2_path)

from digitalhistopathology.access_token import READ_ONLY_HF_TOKEN


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Constants
MAX_NUM_OF_MEM_EVENTS_PER_SNAPSHOT = 100000
SNAPSHOT_FILE_PREFIX = "memory_snapshot"
warnings.filterwarnings("ignore", message="Applied workaround for CuDNN issue")
if (
    timm.__version__ == "0.5.4"
):  # when using CTransPath as it uses a modified version of the timm library
    from timm.models.layers.helpers import to_2tuple
from timm.optim.lamb import Lamb
from timm.optim.lars import Lars


from torch.optim import AdamW
import urllib.request

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

import tqdm
from enum import Enum


class Model(nn.Module):
    """Class for loading model from a filename."""

    def __init__(
        self,
        name,
        model,
        result_folder="../results",
    ):
        """
        Args:
            pretrained_model_path (str): Path to the model weights.
            name (str): Name of the pretrained model.
        """
        super().__init__()

        self.name = name
        self.model = model.to(device)

        self.result_folder = result_folder
        if not os.path.exists(self.result_folder):
            os.mkdir(self.result_folder)

    def model_head(self, model, head, loose_head=False):
        """Load the pretrained model and set it to the correct device to be ready to do inference.

        Args:
            head (str, optional): The head of the loaded model, can be "mlp" or nothing. Defaults to "".

        Returns:
            model: loaded model with added head
        """

        if head.lower() == "mlp":

            if hasattr(model, "fc"):
                in_features = model.fc.in_features
                if loose_head:
                    model.fc = nn.Identity()
                else:
                    model.fc = nn.Sequential(
                        nn.Linear(in_features, 512), nn.ReLU(), nn.Linear(512, 128)
                    )
            elif hasattr(model, "head"):
                if hasattr(model.head, "fc"):
                    in_features = model.head.fc.in_features
                    if loose_head:
                        model.head.fc = nn.Identity()
                    else:
                        model.head.fc = nn.Sequential(
                            nn.Linear(in_features, 512), nn.ReLU(), nn.Linear(512, 128)
                        )
                else:
                    if hasattr(model.head, "in_features"):
                        in_features = model.head.in_features
                    else:
                        # Create a dummy input (e.g., batch size of 1, input size 784)
                        dummy_input = torch.randn(1, 3, 224, 224).to(device)

                        # Get the output size
                        output = model(dummy_input)
                        in_features = output.size()[1]
                    if loose_head:
                        model.head = nn.Identity()
                    else:
                        model.head = nn.Sequential(
                            nn.Linear(in_features, 512), nn.ReLU(), nn.Linear(512, 128)
                        )
            else:
                raise AttributeError("Model does not have 'fc' or 'head' attribute.")
        if torch.cuda.is_available():
            model = model.cuda()
        return model

    def contrastive_loss(self, features):
        """

        Args:
            features (torch.tensor): Hidden feature representation of shape [b, 2, dim].

        Returns:
            float: loss
        """

        b, n, dim = features.size()
        assert n == 2
        mask = torch.eye(b, dtype=torch.float32).to(device)
        contrast_features = torch.cat(torch.unbind(features, dim=1), dim=0)
        anchor = features[:, 0]
        # Dot product
        dot_product = torch.matmul(anchor, contrast_features.T) / self.temperature
        # Log-sum trick for numerical stability
        logits_max, _ = torch.max(dot_product, dim=1, keepdim=True)
        logits = dot_product - logits_max.detach()
        mask = mask.repeat(1, 2)
        logits_mask = torch.scatter(
            torch.ones_like(mask), 1, torch.arange(b).view(-1, 1).to(device), 0
        )
        mask = mask * logits_mask
        # Log-softmax
        exp_logits = torch.exp(logits) * logits_mask
        log_prob = logits - torch.log(exp_logits.sum(1, keepdim=True))
        # Mean log-likelihood for positive
        loss = -((mask * log_prob).sum(1) / mask.sum(1)).mean()
        return loss

    def freeze_layers(self, model, first_freezed_layer_number, freeze_percentage):
        total_params = sum(p.numel() for p in model.parameters())
        params_to_freeze = int(total_params * freeze_percentage)
        frozen_params = 0

        for i, layer in enumerate(list(model.children())):
            if frozen_params >= params_to_freeze or i < first_freezed_layer_number:
                break
            print(f"Freezed layer {layer.__class__.__name__}")
            for param in layer.parameters():
                if frozen_params >= params_to_freeze or i < first_freezed_layer_number:

                    break
                param.requires_grad = False
                frozen_params += param.numel()
                # Check if the layer is a Transformer block (assuming your model uses Vision Transformers)
                if isinstance(layer, nn.TransformerEncoderLayer) or isinstance(
                    layer, nn.TransformerEncoder
                ):
                    # If LayerNorm is used, switch it to eval mode (similar to BatchNorm handling in ResNet)
                    if hasattr(layer, "norm1") and isinstance(
                        layer.norm1, nn.LayerNorm
                    ):
                        layer.norm1.eval()
                    if hasattr(layer, "norm2") and isinstance(
                        layer.norm2, nn.LayerNorm
                    ):
                        layer.norm2.eval()


    def train_one_epoch(self, model, dataloader, optimizer, epoch):

        scaler = amp.GradScaler()
        losses = []

        for i, batch in enumerate(
            tqdm.tqdm(
                dataloader,
                desc=f"Epoch {epoch+1}/{self.num_epochs}",
                unit="batch",
                leave=False,
            )
        ):
            aug_copy_1 = batch["aug_copy_1"]
            aug_copy_2 = batch["aug_copy_2"]
            b, c, h, w = aug_copy_1.size()
            input_ = torch.cat([aug_copy_1.unsqueeze(1), aug_copy_2.unsqueeze(1)], dim=1)
            input_ = input_.view(-1, c, h, w)

            if torch.cuda.is_available():
                input_ = input_.cuda(non_blocking=True)

            optimizer.zero_grad()

            with amp.autocast():

                output = model(input_).view(b, 2, -1)

                loss = self.contrastive_loss(output)

            if torch.isnan(loss) or torch.isinf(loss):
                print(f"NaN or Inf loss detected at epoch {epoch}, batch {i}")
                print(f"Loss: {loss}")
                print(f"Output: {output}")
                print(f"Input: {input_}")
                break

            scaler.scale(loss).backward()
            utils.clip_grad_norm_(model.parameters(), max_norm=1)
            scaler.step(optimizer)
            scaler.update()

            losses.append(loss.item())

            if torch.cuda.is_available():
                torch.cuda.empty_cache()

        return np.mean(losses)


# class CTransPath(Model):
#     """CTransPath pretrained model on digital histopathology datasets, from https://github.com/Xiyue-Wang/TransPath.
#     Get the backbone embeddings without the linear projection head."""
#     def __init__(
#         self,
#         pretrained_model_path="../pretrained_models/ctranspath.pth",
#         head="",
#         pretrained=False,
#         path_to_model=None,
#         lr=0.3,
#         beta1=0.9,
#         beta2=0.999,
#         weight_decay=1e-6,
#     ):
#         print("Load pretrained ctranspath model...")
#         model = torch.load("pytorch_model.bin", map_location=device)
#         if torch.cuda.is_available():
#             model = model.cuda()
#         super().__init__( name="CTransPath",model=model)

class PretrainedModel:
    """Class for loading pretrained model from a filename."""

    def __init__(self, pretrained_model_path, name):
        """
        Args:
            pretrained_model_path (str): Path to the model weights.
            name (str): Name of the pretrained model.
        """
        super().__init__()
        self.pretrained_model_path = pretrained_model_path
        self.name = name

    def load(self):
        """Load the pretrained model and set it to the correct device to be ready to do inference."""
        pass


class CTransPath(PretrainedModel):
    """CTransPath pretrained model on digital histopathology datasets, from https://github.com/Xiyue-Wang/TransPath.
    Get the backbone embeddings without the linear projection head."""

    def __init__(self):
        default_path = os.path.join(project_root, "pretrained_models", "ctranspath.pth")
        super().__init__(default_path, name="CTransPath")
        
        self.model = self.load()

    def load(self):
        print("Load pretrained ctranspath model...")
        model = timm.create_model("swin_tiny_patch4_window7_224", embed_layer=ConvStem, pretrained=False)
        model.head = nn.Identity()
        td = torch.load(self.pretrained_model_path, map_location=torch.device("cuda" if torch.cuda.is_available() else "cpu"), weights_only=False)
        model.load_state_dict(td["model"], strict=True)
        if torch.cuda.is_available():
            model = model.cuda()
        return model
    
def download_model(self, url, save_path):
    print(f"Downloading model from {url}...")
    urllib.request.urlretrieve(url, save_path)
    print(f"Model downloaded and saved to {save_path}")


class SimCLR(Model):
    """SimCLR pretrained model on digital histopathology datasets, from https://github.com/ozanciga/self-supervised-histopathology.
    Get the pre-activation embeddings (just after the CNN without the MLP)
    As we also retrained SimCLR, the class is more developed than for the other pretrained models.
    """

    def __init__(
        self,
        pretrained_model_path=None,
        name="SimCLR",
        batch_size=512,
        num_epochs=100,
        head="",
        pretrained=False,
        path_to_model=None,
        lr=0.3,
        weight_decay=1e-6,  # 1e-5
        beta1=0.9,
        beta2=0.999,
        temperature=0.1,
        epsilon=1e-08,
        loose_head=False,
    ):
        """
        Args:
            pretrained_model_path (str, optional): Path the pretrained or retrained model weigths. Defaults to None (uses default pretrained model).
            name (str, optional): Name of the model. Defaults to "SimCLR".
            training_dataset_files_list (list, optional): List with patches path to be used for the retraining. Defaults to [].

        """

        print("Load pretrained simclr model...")
        model = torchvision.models.__dict__["resnet18"](pretrained=True)
        super().model_head(model, head, loose_head=loose_head)
        if pretrained_model_path:
            state = torch.load(
                pretrained_model_path,
                map_location=device,
                weights_only=False
            )
        else:
            default_path = os.path.join(project_root, "pretrained_models", "tenpercent_resnet18.ckpt")
            state = torch.load(
                default_path,
                map_location=torch.device(
                    "cuda" if torch.cuda.is_available() else "cpu"
                ),
            )

        if "state_dict" in list(state.keys()):
            state_dict = state["state_dict"]
        else:
            state_dict = state["model_state"]
        if "hyper_parameters" in list(state.keys()):
            self.pretrained_hyper_params = state["hyper_parameters"]
        for key in list(state_dict.keys()):
            state_dict[key.replace("model.", "").replace("resnet.", "")] = (
                state_dict.pop(key)
            )
        model_dict = model.state_dict()
        weights = {k: v for k, v in state_dict.items() if k in model_dict}
        if weights == {}:
            print("No weight could be loaded..")
        model_dict.update(weights)
        model.load_state_dict(model_dict)
        super().__init__(name, model)


class UNI(Model):
    """Access must be granted to use this model, you must agree to the outlined terms of use: https://huggingface.co/MahmoodLab/UNI"""

    def __init__(
        self,
        head="",
        pretrained_model_path=None,
    ):
        login(token=READ_ONLY_HF_TOKEN)  # logout done at then end of pipeline.run
        if pretrained_model_path:
            print("Load pretrained uni model...")
            if device == torch.device("cpu"):
                model = torch.load(pretrained_model_path, map_location=device, weights_only=False)
            else:
                model = torch.load(pretrained_model_path, weights_only=False)
        else:
            print("Load MahmoodLab/UNI model from Hugging Face Hub...")
            model = timm.create_model(
                "hf-hub:MahmoodLab/uni",
                pretrained=True,
                init_values=1e-5,
                dynamic_img_size=True,
            )

        super().model_head(model, head)

        super().__init__(name="UNI", model=model)
        
class Uni2(Model):
    """Access must be granted to use this model, you must agree to the outlined terms of use: https://huggingface.co/MahmoodLab/UNI2-h"""

    def __init__(
        self,
        head="",
        pretrained_model_path=None,
    ):
        login(token=READ_ONLY_HF_TOKEN)  # logout done at then end of pipeline.run
        if pretrained_model_path:
            print("Load pretrained uni model...")
            if device == torch.device("cpu"):
                model = torch.load(pretrained_model_path, map_location=device, weights_only=False)
            else:
                model = torch.load(pretrained_model_path, weights_only=False)
        else:
            print("Load MahmoodLab/UNI model from Hugging Face Hub...")
            timm_kwargs = {
            'img_size': 224, 
            'patch_size': 14, 
            'depth': 24,
            'num_heads': 24,
            'init_values': 1e-5, 
            'embed_dim': 1536,
            'mlp_ratio': 2.66667*2,
            'num_classes': 0, 
            'no_embed_class': True,
            'mlp_layer': timm.layers.SwiGLUPacked, 
            'act_layer': torch.nn.SiLU, 
            'reg_tokens': 8, 
            'dynamic_img_size': True
        }
            model = timm.create_model("hf-hub:MahmoodLab/UNI2-h", pretrained=True, **timm_kwargs)

        super().model_head(model, head)

        super().__init__(name="UNI2", model=model)
        
class Virchow(Model):
    """Access must be granted to use this model, you must agree to the outlined terms of use: https://huggingface.co/MahmoodLab/Virchow"""

    def __init__(
        self,
        head="",
        pretrained_model_path=None,
    ):
        login(token=READ_ONLY_HF_TOKEN)  # logout done at the end of pipeline.run
        if pretrained_model_path:
            print("Load pretrained Virchow model...")
            if device == torch.device("cpu"):
                model = torch.load(pretrained_model_path, map_location=device, weights_only=False)
            else:
                model = torch.load(pretrained_model_path, weights_only=False)
        else:
            print("Load MahmoodLab/Virchow model from Hugging Face Hub...")
            # need to specify MLP layer and activation function for proper init
            model = timm.create_model("hf-hub:paige-ai/Virchow2", pretrained=True, mlp_layer=timm.layers.SwiGLUPacked, act_layer=torch.nn.SiLU)
        
        # Modify the forward method to return only the cls_token
        original_forward = model.forward

        def forward_with_cls_token(x):
            output = original_forward(x)
            return output[:, 0]  # Assuming the cls_token is the first token

        model.forward = forward_with_cls_token

        super().model_head(model, head)

        super().__init__(name="Virchow", model=model)
class ProvGigaPath(Model):
    """Access must be granted to use this model, you must agree to the outlined terms of use: https://huggingface.co/prov-gigapath/prov-gigapath"""

    def __init__(self):
        login(token=READ_ONLY_HF_TOKEN)  # logout done at the end of pipeline.run
        model = timm.create_model("hf_hub:prov-gigapath/prov-gigapath", pretrained=True)

        super().__init__(

            name="ProvGigaPath",
            model=model
        )


class ConvStem(nn.Module):
    """Embedding layer for CTransPath. Copied from https://github.com/Xiyue-Wang/TransPath/blob/main/ctran.py."""

    def __init__(
        self,
        img_size=224,
        patch_size=4,
        in_chans=3,
        embed_dim=768,
        norm_layer=None,
        flatten=True,
    ):
        super().__init__()

        assert patch_size == 4
        assert embed_dim % 8 == 0

        img_size = to_2tuple(img_size)
        patch_size = to_2tuple(patch_size)
        self.img_size = img_size
        self.patch_size = patch_size
        self.grid_size = (
            img_size[0] // patch_size[0],
            img_size[1] // patch_size[1],
        )
        self.num_patches = self.grid_size[0] * self.grid_size[1]
        self.flatten = flatten

        stem = []
        input_dim, output_dim = 3, embed_dim // 8
        for l in range(2):
            stem.append(
                nn.Conv2d(
                    input_dim,
                    output_dim,
                    kernel_size=3,
                    stride=2,
                    padding=1,
                    bias=False,
                )
            )
            stem.append(nn.BatchNorm2d(output_dim))
            stem.append(nn.ReLU(inplace=True))
            input_dim = output_dim
            output_dim *= 2
        stem.append(nn.Conv2d(input_dim, embed_dim, kernel_size=1))
        self.proj = nn.Sequential(*stem)

        self.norm = norm_layer(embed_dim) if norm_layer else nn.Identity()

    def forward(self, x):
        B, C, H, W = x.shape
        assert (
            H == self.img_size[0] and W == self.img_size[1]
        ), f"Input image size ({H}*{W}) doesn't match model ({self.img_size[0]}*{self.img_size[1]})."
        x = self.proj(x)
        if self.flatten:
            x = x.flatten(2).transpose(1, 2)  # BCHW -> BNC
        x = self.norm(x)
        return x
import collections


def remove_prefix(state_dict, prefix="student."):
    new_state_dict = collections.OrderedDict()
    for k, v in state_dict.items():
        if k.startswith(prefix):
            new_key = k[len(prefix) :]
            new_state_dict[new_key] = v
        else:
            new_state_dict[k] = v
    return new_state_dict
from peft import LoraModel

class Vit(Model):
    def __init__(self, retrained_model_path):
        model = torch.load(
            retrained_model_path,
            map_location=torch.device("cuda" if torch.cuda.is_available() else "cpu"),
        )

        model.to(device)
        if isinstance(model, LoraModel):
            model = model.merge_and_unload()
            torch.save(model, retrained_model_path)

        super().__init__(name="vit", model=model)


class ModelName(Enum):
    SIMCLR = "simclr"
    CTRANSPATH = "ctranspath"
    UNI = "uni"
    PROVGIGAPATH = "provgigapath"
    VIRCHOW = "virchow"
    UNI2 = "uni2"

def load_model(model_name, retrained_model_path):
    model_name = model_name.lower()
    if retrained_model_path != "" and not os.path.exists(retrained_model_path):
        raise Exception("Retrained model path does not exist: {}".format(retrained_model_path))

    if model_name == ModelName.SIMCLR.value:
        if retrained_model_path != "" and os.path.exists(retrained_model_path):
            return SimCLR(head="mlp", pretrained_model_path=retrained_model_path)
        else:
            return SimCLR()
    elif model_name == ModelName.CTRANSPATH.value:
        if retrained_model_path != "" and os.path.exists(retrained_model_path):
            return CTransPath(pretrained_model_path=retrained_model_path)
        else:
            return CTransPath()
    elif model_name == ModelName.UNI.value:
        if retrained_model_path != "" and os.path.exists(retrained_model_path):
            return UNI(
                pretrained_model_path=retrained_model_path,
                head="",
            )
        else:
            return UNI()
    elif model_name == ModelName.PROVGIGAPATH.value:
        if retrained_model_path != "" and os.path.exists(retrained_model_path):
            return ProvGigaPath(pretrained_model_path=retrained_model_path)
        else:
            return ProvGigaPath()
    elif model_name == ModelName.VIRCHOW.value:
        if retrained_model_path != "" and os.path.exists(retrained_model_path):
            return Virchow(pretrained_model_path=retrained_model_path)
        else:
            return Virchow()
    elif model_name == ModelName.UNI2.value:
        if retrained_model_path != "" and os.path.exists(retrained_model_path):
            return Uni2(pretrained_model_path=retrained_model_path)
        else:
            return Uni2()
    elif model_name == "vit":
        return Vit(retrained_model_path)
    else:
        raise Exception("Invalid model name")
