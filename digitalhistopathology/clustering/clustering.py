#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
from digitalhistopathology.clustering.find_clusters import FindClusters
from digitalhistopathology.clustering.clusters_analysis import ClustersAnalysis


class Clustering(FindClusters, ClustersAnalysis):
    """
    This class is a wrapper for the FindClusters and ClustersAnalysis classes.
    It provides a simple interface for clustering analysis.
    """

    def __init__(self, emb=None, saving_plots=False, result_saving_folder=None):
        FindClusters.__init__(self, emb=emb, saving_plots=saving_plots, result_saving_folder=result_saving_folder)
        ClustersAnalysis.__init__(self, emb=emb, saving_plots=saving_plots, result_saving_folder=result_saving_folder)

        