#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#

import torch.optim as optim
import torch.nn as nn

import torch
import numpy as np
from torch.utils.data import DataLoader
from sklearn.metrics import confusion_matrix
from itertools import cycle
from copy import deepcopy
from torch.optim.lr_scheduler import CosineAnnealingLR
import os
from sklearn.model_selection import KFold
from sklearn.preprocessing import LabelEncoder
from torch.utils.data import TensorDataset
import seaborn as sns
from matplotlib.colors import LogNorm
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve

import matplotlib.pyplot as plt
from digitalhistopathology.classification.mlp_classifier import Classification_mlp
from digitalhistopathology.classification.knn_classifier import Classification_knn


class Classification(Classification_mlp, Classification_knn):

    """
        EmbClassification is a utility class that provides functionality for performing 
        classification tasks using embeddings saved into AnnData format. It includes methods for multi-layer perceptron (MLP) 
        classification with 5-fold cross-validation and visualization of classification results.
        Attributes:
            emb (AnnData): AnnData object containing embeddings data and labels for classification in the "obs".
        Methods:
            __init__(emb=None):
                Initializes the EmbClassification with optional embeddings.
            mlp_classification(save_file, label_column="label"):
                Performs MLP classification on the embeddings using 5-fold cross-validation. 
                Saves the trained models, confusion matrix, and bar plots of prediction distributions.
                    save_file (str): Path to save the classification results, including models 
                                    and visualizations.
                    label_column (str): Column name in the AnnData object containing the labels for classification.
                    Default is "label".
        Example:
            emb_classification = EmbClassification(emb=ad_adata)
            emb_classification.mlp_classification(save_file="results", label_column="cell_type")
    """


    def __init__(self, emb=None, saving_plots=False, result_saving_folder=None):

        Classification_mlp.__init__(self, emb=emb)
        Classification_knn.__init__(self, emb=emb, 
                                       saving_plots=saving_plots, 
                                       result_saving_folder=result_saving_folder)
        
        

