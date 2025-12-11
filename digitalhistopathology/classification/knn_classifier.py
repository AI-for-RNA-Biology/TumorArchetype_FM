#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#

import os
import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score

import matplotlib.pyplot as plt

from digitalhistopathology.dimred.dimensionality_reduction import DimRed


class Classification_knn(DimRed):
    """
    KNN classifier for embedding classification.
    """

    def __init__(self, emb=None, saving_plots=False, result_saving_folder=None):

        DimRed.__init__(self, 
                                                  emb=emb, 
                                                  saving_plots=saving_plots, 
                                                  result_saving_folder=result_saving_folder)

    def get_optimal_k_for_knn_per_patient(
        self,
        svd_comp=None,
        weights="uniform",
        metric="minkowski",
        k_list=[3, 5, 7, 9, 11, 13, 15, 17, 19, 21],
        verbose=True,
    ):
        """Calculate the optimal k for knn computed on each patient by computing the median of the f1 score of each patient for each k, and select the k with the highest median. Plot
        boxplots of accuracies and f1 scores for each k across patients.

        Args:
            svd_comp (int, optional): Number of svd components on which the knn will be run. If None the entire raw embedding matrix is picked. Defaults to None.
            weights (str, optional): Weights argument of KNeighborsClassifier. Defaults to "uniform".
            metric (str, optional): Metric argument of KNeighborsClassifier. Defaults to "minkowski".
            k_list (list, optional): The list of different ks to be tested. Defaults to [3, 5, 7, 9, 11, 13, 15, 17, 19, 21].
            plot_results (bool, optional): If the accuracies and f1 scores boxplots are plotted. Defaults to True.
            verbose (bool, optional): If printed information during the computation is wanted. Defaults to True.

        Returns:
            int: the optimal number of neighbours for knn (k)
        """
        print(f"Here emb size: {self.emb.shape}", flush=True)
        label_emb = Classification_knn()
        label_emb.emb = self.emb[~self.emb.obs["label"].isna()]
        label_emb.emb = label_emb.emb[label_emb.emb.obs["label"] != "undetermined"]

        print(f"Here label_emb size: {label_emb.emb.shape}", flush=True)
        accuracies_dict = dict()
        f1_dict = dict()
        for patient in list(self.emb.obs.tumor.unique()):
            if verbose:
                print("\nPatient " + patient)
            current_emb = Classification_knn()
            current_emb.emb = label_emb.emb[(label_emb.emb.obs["tumor"] == patient), :]

            if svd_comp is None:
                X = current_emb.emb.X
            else:
                try:
                    label_emb.svd["U_df"] = self.svd["U_df"].loc[label_emb.emb.obs.index, :]
                    X = self.svd["U_df"].loc[current_emb.emb.obs.index, :].iloc[:, 0:svd_comp].to_numpy()
                except Exception as e:
                    print("Compute svd before wanting to use svd comp for knn: {}".format(e))
                    return
            print(f"Here current_emb size: {current_emb.emb.shape} for patient {patient}", flush=True)
            
            print(f"Here X size: {X.shape} for patient {patient}", flush=True)
            
            X_train, X_test, y_train, y_test = train_test_split(
                X, current_emb.emb.obs.label.values, test_size=0.2, random_state=42
            )
            current_accuracies_list = []
            current_f1_list = []
            for k in k_list:
                knn = KNeighborsClassifier(n_neighbors=k, weights=weights, metric=metric)
                knn.fit(X_train, y_train)
                y_pred = knn.predict(X_test)

                accuracy = accuracy_score(y_test, y_pred)
                f1 = f1_score(y_test, y_pred, average="weighted")
                current_accuracies_list.append(accuracy)
                current_f1_list.append(f1)

                if verbose:
                    print("k = {}".format(k))
                    print("Accuracy for patient {} and k={}: {}".format(patient, k, accuracy))
                    print("F1 for patient {} and k={}: {}".format(patient, k, f1))
            accuracies_dict[patient] = current_accuracies_list
            f1_dict[patient] = current_f1_list

        self.plot_f1_accuracies_knn(accuracies_dict, metric="Accuracy", k_list=k_list)
        self.plot_f1_accuracies_knn(f1_dict, metric="F1", k_list=k_list)

        # TODO: Is it the correct way to calculate this ?
        # df_acc = pd.DataFrame.from_records(accuracies_dict).set_axis(k_list)
        df_f1 = pd.DataFrame.from_records(f1_dict).set_axis(k_list)
        # df_total = (df_acc + df_f1) / 2
        # df_total["median"] = df_total.median(axis=1)
        # optimal_k = df_total[["median"]].idxmax().values[0]
        df_f1["median"] = df_f1.median(axis=1)
        optimal_k = df_f1[["median"]].idxmax().values[0]
        print("The optimal number of neighbours k = {}".format(optimal_k))
        return optimal_k

    def plot_f1_accuracies_knn(
        self, dict_scores, metric="Accuracy", ylim=[0.5, 1], k_list=[3, 5, 7, 9, 11, 13, 15, 17, 19, 21]
    ):
        """Plot boxplots of scores for each k across patients.

        Args:
            dict_scores (dict): Dictionnary with patients as keys and dictionnary as value. The patient dictionnary has the different ks as keys and the score as value.
            metric (str, optional): Name of the metric. Defaults to "Accuracy".
            ylim (list, optional): Y axis limits for the plot. Defaults to [0.5, 1].
            k_list (list, optional): The list of different ks to be tested. Defaults to [3, 5, 7, 9, 11, 13, 15, 17, 19, 21].
        """
        df_plot = pd.DataFrame.from_records(dict_scores)
        df_plot.index = k_list
        df_plot.T.boxplot(grid=False)
        plt.xlabel("k")
        plt.ylabel(metric)
        plt.ylim(ylim)
        plt.title("Boxplot for each k of the knn {} accross patient".format(metric.lower()))

        colors = plt.cm.get_cmap("tab10", len(df_plot.columns))

        for idx, column in enumerate(df_plot.columns):
            plt.scatter(np.arange(1, len(list(df_plot.index)) + 1), df_plot[column], color=colors(idx), label=column, alpha=0.6)

        plt.legend(title="Patient", bbox_to_anchor=(1.05, 1), loc="upper left")
        if self.saving_plots:
            plt.savefig(
                os.path.join(
                    self.result_saving_folder,
                    "knn_optimal_k_search_{}_boxplot.pdf".format(metric.lower()),
                ),
                bbox_inches="tight",
            )
            plt.close()
        else:
            plt.show()

        plt.show()

    def predict_labels_with_knn(
        self, n_neighbors, reclassify_all_data=True, svd_comp=None, weights="uniform", metric="minkowski"
    ):
        """Run knn with k=n_neighbors for each patient and predict the label of the non-labeled patches and can even reannotate the others to mitigate potential mislabelling
        problems.

        Args:
            n_neighbors (int): The number of neighbours (k) for knn.
            reclassify_all_data (bool, optional): If already labeled patches by the pathologist should be relabeled by knn. Defaults to True.
            svd_comp (int, optional): Number of svd components on which the knn will be run. If None the entire raw embedding matrix is picked. Defaults to None.
            weights (str, optional): Weights argument of KNeighborsClassifier. Defaults to "uniform".
            metric (str, optional): Metric argument of KNeighborsClassifier. Defaults to "minkowski".
        """
        self.emb.obs["knn_predicted_label"] = self.emb.obs["label"]

        for patient in list(self.emb.obs.tumor.unique()):
            print("\nPatient " + patient)
            current_emb = Classification_knn()
            current_emb.emb = self.emb[self.emb.obs["tumor"] == patient]

            knn = KNeighborsClassifier(n_neighbors=n_neighbors, weights=weights, metric=metric)

            # Fit the classifier to the training data
            print("Shape of the data: {}".format(current_emb.emb.shape))

            filtered_data = current_emb.emb[
                (~current_emb.emb.obs["label"].isna()) & (current_emb.emb.obs["label"] != "undetermined") & (current_emb.emb.obs["label"] != "nan")
            ]

            print("Shape of the filtered data: {}".format(filtered_data.shape))
            if svd_comp is None:
                X_fit = filtered_data.X
            else:
                try:
                    X_fit = self.svd["U_df"].loc[filtered_data.obs.index, :].iloc[:, 0:svd_comp].to_numpy()
                except Exception as e:
                    print("Compute svd before wanting to use svd comp for knn: {}".format(e))
                    return

            knn.fit(X_fit, filtered_data.obs.label)

            if reclassify_all_data:
                if svd_comp is None:
                    X_predict = current_emb.emb.X
                else:
                    X_predict = self.svd["U_df"].loc[current_emb.emb.obs.index, :].iloc[:, 0:svd_comp].to_numpy()

                if weights == "uniform":
                    # Predict the labeled patches by keeping them in the trained model. As it is not weighted by distance it works.
                    # More stringent on the differences between pathologist's labels and the predicted ones if the pathologist's label of
                    # the patch is included in the knn decision
                    self.emb.obs["knn_predicted_label"][self.emb.obs["tumor"] == patient] = knn.predict(X_predict)
                else:
                    # Must fit one model per pathologist's labeled patch by fitting on the labeled data without the current labeled patch
                    # Otherwise, if it current labeled patch is in the training data, the distance between him and himself will be 0 and
                    # so the weight will 1/0 = infinite
                    non_label_emb = self.emb[
                        (self.emb.obs["tumor"] == patient)
                        & (self.emb.obs["label"].isna() | (self.emb.obs["label"] == "undetermined"))
                    ]
                    non_label_pred = knn.predict(
                        non_label_emb.X
                        if svd_comp is None
                        else self.svd["U_df"].loc[non_label_emb.obs.index, :].iloc[:, 0:svd_comp].to_numpy()
                    )
                    label_pred = []
                    for i in range(X_fit.shape[0]):
                        current_X = np.delete(X_fit, (i), axis=0)
                        current_y = np.delete(filtered_data.obs.label, (i))
                        knn = KNeighborsClassifier(n_neighbors=n_neighbors, weights=weights, metric=metric)
                        knn.fit(current_X, current_y)
                        label_pred.append(knn.predict(X_fit[i, :].reshape(1, -1))[0])
                    self.emb.obs["knn_predicted_label"][
                        (self.emb.obs["tumor"] == patient)
                        & (self.emb.obs["label"].isna() | (self.emb.obs["label"] == "undetermined"))
                    ] = non_label_pred
                    self.emb.obs["knn_predicted_label"][
                        (self.emb.obs["tumor"] == patient)
                        & (~self.emb.obs["label"].isna() & (self.emb.obs["label"] != "undetermined"))
                    ] = label_pred
            else:
                emb_predict = current_emb.emb[
                    current_emb.emb.obs["label"].isna() | (current_emb.emb.obs["label"] == "undetermined") | (current_emb.emb.obs["label"] == "nan")
                ]
                print("Shape of the data to predict: {}".format(emb_predict.shape))
                if svd_comp is None:
                    X_predict = emb_predict.X
                else:
                    X_predict = self.svd["U_df"].loc[emb_predict.obs.index, :].iloc[:, 0:svd_comp].to_numpy()

                label_pred = knn.predict(X_predict)

                self.emb.obs["knn_predicted_label"][
                    (self.emb.obs["tumor"] == patient)
                    & (self.emb.obs["label"].isna() | (self.emb.obs["label"] == "undetermined"))
                ] = label_pred

    def get_enrichement_score(self, groupby="knn_predicted_label"):
        """Compute the enrichement score of labels by calculating the fraction of pathologist labeled and knn labeled (groupby) patches that are labeled in the same cluster for each label.

        Args:
            groupby (str, optional): Obs column name on which groupby the data into different labels. Defaults to "knn_predicted_label".

        Returns:
            dict: dictionnary with label as keys and enrichement score as values.
        """
        cluster_info_df = self.get_predicted_clusters_infos(
            groupby=groupby, fraction=True, predicted_clusters_to_consider=None
        ).set_index("name")
        results = cluster_info_df.apply(lambda r: r[r.name], axis=1)
        return dict(results)