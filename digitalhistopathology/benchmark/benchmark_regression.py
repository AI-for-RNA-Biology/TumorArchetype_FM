#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#

import os
import glob
import numpy as np
import pandas as pd
import torch
from sklearn.model_selection import KFold, ParameterGrid
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from sklearn.kernel_ridge import KernelRidge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
from sklearn.utils import resample
from joblib import Parallel, delayed
import seaborn as sns
from numpy.polynomial.polynomial import polyfit, polyval

from digitalhistopathology.regression.regression import MLPRegressorPyTorch, adjusted_r2_score, train_mlp_model, evaluate_model
from digitalhistopathology.benchmark.benchmark_base import BenchmarkBase


import torch.nn as nn
import matplotlib.pyplot as plt



class BenchmarkRegression(BenchmarkBase):
    """

    """

    def __init__(self, 
                 path_to_pipeline=None, 
                 pipelines_list=None, 
                 dataset_name=None,
                 results_folder="../results",
                 emb_df_csv_path=None,
                 saving_folder=None,
                 image_embedding_name="image_embedding.h5ad",
                 engineered_features_saving_folder=None,
                 engineered_features_type='scMTOP',
                 extension='png',
                 group='tumor',
                 regression_type="linear",
                 n_splits=5,
                 alpha_reg=1,
                 on_invasive=False,
                 dataset="HER2",
):        
        
        super().__init__(path_to_pipeline=path_to_pipeline, 
                         pipelines_list=pipelines_list, 
                         dataset_name=dataset_name,
                         results_folder=results_folder,
                         emb_df_csv_path=emb_df_csv_path,
                         saving_folder=saving_folder,
                         image_embedding_name=image_embedding_name,
                         engineered_features_saving_folder=engineered_features_saving_folder,
                         engineered_features_type=engineered_features_type,
                         extension=extension,
                         group=group,
                         dataset=dataset)
        
        self.regression_type = regression_type
        self.n_splits = n_splits
        self.alpha_reg = alpha_reg
        self.on_invasive = on_invasive
        
    
    @staticmethod
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
    
    @staticmethod
    def perform_grid_search_polyfit(X, Y, n_splits, param_grid):
        best_r2 = -np.inf
        best_params = None
        scaler = StandardScaler()

        for params in ParameterGrid(param_grid):
            degree = params['degree']

            cv = KFold(n_splits=n_splits, random_state=42, shuffle=True)
            r2_scores = []

            for train_index, test_index in cv.split(X):
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = Y[train_index], Y[test_index]

                scaler.fit(X_train)
                X_train = scaler.transform(X_train)
                X_test = scaler.transform(X_test)

                poly_coeffs = polyfit(X_train.flatten(), y_train, degree)
                y_pred = polyval(poly_coeffs, X_test.flatten())

                r2 = r2_score(y_test, y_pred)
                r2_scores.append(r2)

            mean_r2 = np.mean(r2_scores)
            if mean_r2 > best_r2:
                best_r2 = mean_r2
                best_params = params

        return best_params
    
    @staticmethod
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
            best_params = BenchmarkRegression.perform_grid_search_mlp(X, Y_feature, n_splits, param_grid, device)
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
        elif regression_type == 'polyfit':
            param_grid = {
                'degree': [1, 2, 3, 4, 5]
            }
            best_params = BenchmarkRegression.perform_grid_search_polyfit(X, Y_feature, n_splits, param_grid)
            degree = best_params['degree']
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
            elif regression_type == 'polyfit':
                poly_coeffs = polyfit(X_train.flatten(), y_train, degree)
                y_pred = polyval(poly_coeffs, X_test.flatten())
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

    def regression_handcrafted_features_prediction(self, regression_type="linear", n_splits=5, alpha_reg=1, n_jobs=-1):
        

        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        scaler = StandardScaler()

        r2_scores_per_model = {}
        adjusted_r2_scores_per_model = {}
        for model in self.pipelines_list:
            if os.path.exists(os.path.join(self.saving_folder, model, "r2_scores.csv")) and os.path.exists(os.path.join(self.saving_folder, model, "adjusted_r2_scores.csv")):
                df_r2 = pd.read_csv(os.path.join(self.saving_folder, model, "r2_scores.csv"))
                df_adjusted_r2 = pd.read_csv(os.path.join(self.saving_folder, model, "adjusted_r2_scores.csv"))
                r2_scores_per_model[model] = df_r2
                adjusted_r2_scores_per_model[model] = df_adjusted_r2
                continue

            X = self.image_embeddings[model].emb.X
            Y = self.ef.emb.X
            print(f"Computing R^2 scores for model {model} and regression {regression_type}...", flush=True)

            results = Parallel(n_jobs=n_jobs, backend='loky')(delayed(self.regress_one_feature)(i, X, Y, regression_type, n_splits, alpha_reg, device) for i in range(Y.shape[1]))

            mean_cv_r2_scores, mean_cv_adjusted_r2_scores = zip(*results)

            r2_scores_per_model[model] = mean_cv_r2_scores
            adjusted_r2_scores_per_model[model] = mean_cv_adjusted_r2_scores

            if not os.path.exists(os.path.join(self.saving_folder, model)):
                os.makedirs(os.path.join(self.saving_folder, model))

            df_r2 = pd.DataFrame(mean_cv_r2_scores, index=self.ef.emb.var_names)
            df_r2.to_csv(os.path.join(self.saving_folder, model, "r2_scores.csv"))

            df_adjusted_r2 = pd.DataFrame(mean_cv_adjusted_r2_scores, index=self.ef.emb.var_names)
            df_adjusted_r2.to_csv(os.path.join(self.saving_folder, model, "adjusted_r2_scores.csv"))

        df_r2 = pd.DataFrame.from_dict(r2_scores_per_model)
        df_r2.index = self.ef.emb.var_names

        df_adjusted_r2 = pd.DataFrame.from_dict(adjusted_r2_scores_per_model)
        df_adjusted_r2.index = self.ef.emb.var_names

        return df_r2, df_adjusted_r2

    def regression_plots(self, df, 
                         list_base_models, 
                         pairs_to_compare=[("uni", "uni_retrained"), 
                                           ("simclr", "simclr_6_freezed_layer_01_temperature")],
                         regression_type="linear",
                         adjusted=True):
        
        if adjusted:
            adjusted_str = "adjusted_"
        else:
            adjusted_str = ""
        plt.figure()
        sns.barplot(data=df)
        sns.despine()
        plt.xticks(rotation=90)
        plt.title(f"{adjusted_str.replace('_', ' ').capitalize()} R2 scores all handcrafted features prediction \n {regression_type} models", weight='bold')
        plt.savefig(os.path.join(self.saving_folder, f"{adjusted_str}r2_scores_{regression_type}_all_models_barplot.{self.extension}"), bbox_inches='tight')


        plt.figure(figsize=(10, 15))

        plt.subplot(1+int(np.ceil(len(pairs_to_compare)/2)), 2, 1)
        sns.barplot(data=df[list_base_models])
        sns.despine()
        plt.xticks(rotation=90)
        plt.title(f"{adjusted_str.replace('_', ' ').capitalize()} R2 scores all handcrafted features prediction \n {regression_type} models", weight='bold')
        plt.ylabel(f"{adjusted_str.replace('_', ' ').capitalize()} R2 score")

        for i, pair in enumerate(pairs_to_compare):
            df_change = df[[pair[0]]].apply(lambda x: df[pair[1]] - x)
            plt.subplot(1+int(np.ceil(len(pairs_to_compare)/2)), 2, 3+i)            
            sns.barplot(data=df_change)
            sns.despine()
            plt.xticks(rotation=90)
            plt.title(f"Change in {adjusted_str.replace('_', ' ').capitalize()} scores \n {pair[1]} - {pair[0]}", weight='bold')
            plt.ylabel(r"$\Delta$R2 score")
            # plt.ylim(-0.1, 0.1)

        plt.tight_layout()
        plt.savefig(os.path.join(self.saving_folder, f"{adjusted_str}r2_scores_{regression_type}_base_models_comparison_with_retrained.{self.extension}"), bbox_inches='tight')


        plt.figure(figsize=(10, 15))
        df_melted = pd.melt(df, var_name='model', value_name=f'{adjusted_str}r2_score')
        df_melted.index = list(df.index) * len(df.columns)
        df_melted['feature_type'] = [x.split('_')[0] for x in df_melted.index]
        plt.subplot(1+int(np.ceil(len(pairs_to_compare)/2)), 2, 1)
        sns.barplot(df_melted[df_melted['model'].isin(list_base_models)], x='feature_type', y=f'{adjusted_str}r2_score', hue='model')
        sns.despine()
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.xticks(rotation=90)
        plt.title(f"{adjusted_str.replace('_', ' ').capitalize()} R2 scores for base models", weight='bold')
        plt.ylabel(f"{adjusted_str.replace('_', ' ').capitalize()} R2 score")

        for i, pair in enumerate(pairs_to_compare):
            df_change = df[[pair[0]]].apply(lambda x: df[pair[1]] - x)
            df_melted = pd.melt(df_change, var_name='model', value_name=f'delta_{adjusted_str}r2_score')
            df_melted.index = list(df_change.index) * len(df_change.columns)
            df_melted['feature_type'] = [x.split('_')[0] for x in df_melted.index]

            plt.subplot(1+int(np.ceil(len(pairs_to_compare)/2)), 2, 3+i)    
            sns.barplot(df_melted[df_melted['model'].isin(list_base_models)], x='feature_type', y=f'delta_{adjusted_str}r2_score', hue='model')
            sns.despine()
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            plt.xticks(rotation=90)
            plt.title(f"Change in {adjusted_str.replace('_', ' ').capitalize()} R2 score \n {pair[1]} - {pair[0]}", weight='bold')
            plt.ylabel(r"$\Delta$R2 score")
            # plt.ylim(-0.1, 0.1)


        plt.suptitle(f"R2 scores for handcrafted features prediction per feature type \n {regression_type} models", weight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(self.saving_folder, f"{adjusted_str}r2_scores_{regression_type}_base_models_comparison_with_retrained_per_feature_type.{self.extension}"), bbox_inches='tight')


    def regression_pipeline(self, 
                            regression_type="linear", 
                            n_splits=5, 
                            alpha_reg=1,
                            on_invasive=False,
                            ):
        
        

        if not os.path.exists(os.path.join(self.saving_folder, regression_type)):
            os.makedirs(os.path.join(self.saving_folder, regression_type))


        self.saving_folder = os.path.join(self.saving_folder, regression_type)

        if on_invasive:
            self.saving_folder = os.path.join(self.saving_folder, "invasive")

        if not os.path.exists(self.saving_folder):
            os.makedirs(self.saving_folder)

        if os.path.exists(f"{self.saving_folder}/r2_scores_{regression_type}.csv"):
            df_r2 = pd.read_csv(f"{self.saving_folder}/r2_scores_{regression_type}.csv", index_col=0)
            df_adjusted_r2 = pd.read_csv(f"{self.saving_folder}/adjusted_r2_scores_{regression_type}.csv", index_col=0)
        else:
            df_r2, df_adjusted_r2 = self.regression_handcrafted_features_prediction(regression_type=regression_type, 
                                                                                    n_splits=n_splits,
                                                                                    alpha_reg=alpha_reg)
            df_r2.to_csv(f"{self.saving_folder}/r2_scores_{regression_type}.csv")
            df_adjusted_r2.to_csv(f"{self.saving_folder}/adjusted_r2_scores_{regression_type}.csv")

    def execute_pipeline(self):

        if self.on_invasive:

            self.load_invasive_image_embeddings()
            self.image_embeddings = self.invasive_image_embeddings.copy()
        
        self.regression_pipeline(regression_type=self.regression_type,
                                 n_splits=self.n_splits,
                                 alpha_reg=self.alpha_reg,
                                 on_invasive=self.on_invasive)