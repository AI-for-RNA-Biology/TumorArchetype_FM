#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import torch
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import os 
from sklearn.metrics import adjusted_rand_score
from sklearn.model_selection import KFold, ParameterGrid
from sklearn.linear_model import LinearRegression
# from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.svm import SVR
from sklearn.kernel_ridge import KernelRidge
from sklearn.utils import resample
from sklearn.preprocessing import StandardScaler  
import anndata as ad
import argparse

class MLPRegressorPyTorch(nn.Module):
    def __init__(self, input_dim, hidden_dim=100):
        super(MLPRegressorPyTorch, self).__init__()
        self.hidden = nn.Linear(input_dim, hidden_dim)
        self.output = nn.Linear(hidden_dim, 1)
        self.relu = nn.ReLU()

    def forward(self, x):
        x = self.relu(self.hidden(x))
        x = self.output(x)
        return x


def train_mlp_model(X_train, y_train, input_dim, hidden_dim, learning_rate, epochs, batch_size, device, weight_decay=0.0):
    """
    Trains a Multi-Layer Perceptron (MLP) regression model using PyTorch.

    Args:
        X_train (torch.Tensor): The input features for training.
        y_train (torch.Tensor): The target values for training.
        input_dim (int): The number of input features.
        hidden_dim (int): The number of hidden units in the MLP.
        learning_rate (float): The learning rate for the optimizer.
        epochs (int): The number of epochs to train the model.
        batch_size (int): The batch size for training.

    Returns:
        MLPRegressorPyTorch: The trained MLP regression model.
    """
    model = MLPRegressorPyTorch(input_dim, hidden_dim).to(device)
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)

    train_dataset = torch.utils.data.TensorDataset(X_train, y_train)
    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=4)

    model.train()
    for epoch in range(epochs):
        for batch_X, batch_y in train_loader:
            batch_X, batch_y = batch_X.to(device), batch_y.to(device)
            optimizer.zero_grad()
            outputs = model(batch_X)
            loss = criterion(outputs, batch_y)
            loss.backward()
            optimizer.step()

    return model


def evaluate_model(model, X_test, device):
    """
    Evaluates a trained regression model on a test set.
    
    Args:
        model (MLPRegressorPyTorch): The trained regression model.
        X_test (torch.Tensor): The input features for testing.
        
    Returns:
        numpy.ndarray: The predicted target values.
    """
    
    model.eval()
    with torch.no_grad():
        X_test = X_test.to(device)
        y_pred_tensor = model(X_test)
    return y_pred_tensor.cpu().numpy().flatten()



def adjusted_r2_score(r2, n, p):
    return 1 - ((1 - r2) * (n - 1) / (n - p - 1))


def perform_grid_search_mlp(X, Y, n_splits, param_grid, device):
    best_r2 = -np.inf
    best_params = None
    input_dim = X.shape[1]
    scaler = StandardScaler()

    for params in ParameterGrid(param_grid):
        hidden_dim = params['hidden_dim']
        learning_rate = params['learning_rate']
        batch_size = params['batch_size']
        epochs = params['epochs']
        weight_decay = params['weight_decay']

        cv = KFold(n_splits=n_splits, random_state=42, shuffle=True)
        r2_scores = []

        for train_index, test_index in cv.split(X):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = Y[train_index], Y[test_index]

            scaler.fit(X_train)
            X_train = scaler.transform(X_train)
            X_test = scaler.transform(X_test)

            X_train_tensor = torch.tensor(X_train, dtype=torch.float32)
            y_train_tensor = torch.tensor(y_train, dtype=torch.float32).view(-1, 1)
            X_test_tensor = torch.tensor(X_test, dtype=torch.float32)

            model = train_mlp_model(X_train_tensor, y_train_tensor, input_dim, hidden_dim, learning_rate, epochs, batch_size, device, weight_decay)
            y_pred = evaluate_model(model, X_test_tensor, device)

            r2 = r2_score(y_test, y_pred)
            r2_scores.append(r2)

        mean_r2 = np.mean(r2_scores)
        if mean_r2 > best_r2:
            best_r2 = mean_r2
            best_params = params

    return best_params


def regress_one_feature(i, X, Y, regression_type, n_splits, alpha_reg, device):
    print(f"Processing feature {i}...", flush=True)
    scaler = StandardScaler()
    Y_feature = Y[:, i]

    if regression_type == "linear":
        model_reg = LinearRegression()
    elif regression_type == 'mlp':
        param_grid = {
            'hidden_dim': [100, 200],
            'learning_rate': [0.01, 0.001],
            'weight_decay': [0, 0.01],
            'batch_size': [2056],
            'epochs': [10]
        }
        best_params = perform_grid_search_mlp(X, Y_feature, n_splits, param_grid, device)
        hidden_dim = best_params['hidden_dim']
        learning_rate = best_params['learning_rate']
        weight_decay = best_params['weight_decay']
        batch_size = best_params['batch_size']
        epochs = best_params['epochs']

        model_reg = MLPRegressorPyTorch(X.shape[1], hidden_dim)
        criterion = nn.MSELoss()
        optimizer = torch.optim.Adam(model_reg.parameters(), lr=learning_rate, weight_decay=weight_decay)
    elif regression_type == 'svr':
        model_reg = SVR(kernel='rbf')
    elif regression_type == 'kernelridge':
        model_reg = KernelRidge(kernel='rbf', alpha=alpha_reg)

    else:
        raise ValueError("Model type not recognized. Choose between 'linear', 'mlp', 'svr', 'kernelridge', and 'polyfit'.")

    cv = KFold(n_splits=n_splits, random_state=42, shuffle=True)
    r2_scores = []
    adjusted_r2_scores = []

    for train_index, test_index in cv.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = Y_feature[train_index], Y_feature[test_index]

        if regression_type == 'kernelridge' and len(X_train) > 8000:
            X_train, y_train = resample(X_train, y_train, n_samples=8000, random_state=42)

        scaler.fit(X_train)
        X_train = scaler.transform(X_train)
        X_test = scaler.transform(X_test)

        if regression_type == 'mlp':
            X_train_tensor = torch.tensor(X_train, dtype=torch.float32)
            y_train_tensor = torch.tensor(y_train, dtype=torch.float32).view(-1, 1)
            X_test_tensor = torch.tensor(X_test, dtype=torch.float32)

            mlp_model = train_mlp_model(X_train_tensor, y_train_tensor, X.shape[1], hidden_dim, learning_rate, epochs, batch_size, device, weight_decay)
            y_pred = evaluate_model(mlp_model, X_test_tensor, device)
        else:
            model_reg.fit(X_train, y_train)
            y_pred = model_reg.predict(X_test)

        n = len(y_test)
        p = X.shape[1]

        r2 = r2_score(y_test, y_pred)
        r2_scores.append(r2)

        adjusted_r2 = adjusted_r2_score(r2, n, p)
        adjusted_r2_scores.append(adjusted_r2)

    mean_r2 = np.mean(r2_scores)
    mean_adjusted_r2 = np.mean(adjusted_r2_scores)

    return mean_r2, mean_adjusted_r2

def regression_handcrafted_features_prediction(embedding, df_handcrafted, regression_type="linear", n_splits=5, alpha_reg=1, n_jobs=-1, saving_folder="results", model=""):
    

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    if os.path.exists(os.path.join(saving_folder, model, "r2_scores.csv")) and os.path.exists(os.path.join(saving_folder, model, "adjusted_r2_scores.csv")):
        df_r2 = pd.read_csv(os.path.join(saving_folder, model, "r2_scores.csv"))
        df_adjusted_r2 = pd.read_csv(os.path.join(saving_folder, model, "adjusted_r2_scores.csv"))

        
    else:
        X = embedding.X
        Y = df_handcrafted.values
        handcrafted_features_names = df_handcrafted.columns

        results = Parallel(n_jobs=n_jobs, backend='loky')(delayed(regress_one_feature)(i, X, Y, regression_type, n_splits, alpha_reg, device) for i in range(Y.shape[1]))

        mean_cv_r2_scores, mean_cv_adjusted_r2_scores = zip(*results)

        if not os.path.exists(os.path.join(saving_folder, model)):
            os.makedirs(os.path.join(saving_folder, model))

        df_r2 = pd.DataFrame(mean_cv_r2_scores, index=handcrafted_features_names)
        df_r2.to_csv(os.path.join(saving_folder, model, "r2_scores.csv"))

        df_adjusted_r2 = pd.DataFrame(mean_cv_adjusted_r2_scores, index=handcrafted_features_names)
        df_adjusted_r2.to_csv(os.path.join(saving_folder, model, "adjusted_r2_scores.csv"))

    return df_r2, df_adjusted_r2


def main():
    parser = argparse.ArgumentParser(description="Run regression on handcrafted features.")
    parser.add_argument("--embedding_path", type=str, required=True, help="Path to the embedding .h5ad file.")
    parser.add_argument("--handcrafted_features_path", type=str, required=True, help="Path to the handcrafted features CSV file.")
    parser.add_argument("--regression_type", type=str, default="linear", choices=["linear", "mlp", "svr", "kernelridge"], help="Type of regression model to use.")
    parser.add_argument("--n_splits", type=int, default=5, help="Number of splits for cross-validation.")
    parser.add_argument("--alpha_reg", type=float, default=1.0, help="Regularization strength for kernel ridge regression.")
    parser.add_argument("--n_jobs", type=int, default=-1, help="Number of parallel jobs to run.")
    parser.add_argument("--saving_folder", type=str, default="results", help="Folder to save the results.")
    parser.add_argument("--model", type=str, default="", help="Model name for saving results.")
    parser.add_argument("--invasive_path", type=str, default="", help="Path to the invasive CSV file.")
    parser.add_argument("--on_invasive", type="store_true", help="Whether to run the regression on the invasive features.")

    args = parser.parse_args()

    embedding = ad.read_h5ad(args.embedding_path)
    df_handcrafted = pd.read_csv(args.handcrafted_features_path, index_col=0)
    
    if args.on_invasive:
        df_invasive = pd.read_csv(args.invasive_path, index_col=0)
        embedding = embedding[df_invasive.index]
        df_handcrafted = df_handcrafted.merge(df_invasive, rleft_index=True, right_index=True)

    df_r2, df_adjusted_r2 = regression_handcrafted_features_prediction(
        embedding, df_handcrafted, args.regression_type, args.n_splits, args.alpha_reg, args.n_jobs, args.saving_folder, args.model
    )

    print("R2 Scores:")
    print(df_r2)
    print("Adjusted R2 Scores:")
    print(df_adjusted_r2)

if __name__ == "__main__":
    main()