#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import scipy
import numpy as np
from sklearn.cluster import KMeans
import ot
from scipy.spatial.distance import pdist, squareform, mahalanobis


def kl_divergence(p, q):
    """Compute the Kullback-Leibler divergence of distribution p from distribution q.

    Args:
        p (np.array): Distribution p.
        q (np.array): Distribution q.

    Returns:
        float: Kullback-Leibler divergence of p from q.
    """
    return np.sum(scipy.special.rel_entr(p, q))


def quantized_wasserstein(matrix, idx_samples_cluster1, idx_samples_cluster2, ref_space=None, k=5000):
    """
    Compute the quantized Wasserstein distance between two clusters of samples.
    Parameters:
    matrix (numpy.ndarray): The data matrix containing all samples.
    idx_samples_cluster1 (list or numpy.ndarray): Indices of samples belonging to the first cluster.
    idx_samples_cluster2 (list or numpy.ndarray): Indices of samples belonging to the second cluster.
    ref_space (numpy.ndarray, optional): Reference matrix for normalization. Default is None. If not None, the
    current space will be normalized to have the same diameter as the reference space.
    k (int, optional): Number of clusters for k-means. Default is 5000.
    Returns:
    float: The quantized Wasserstein distance between the two clusters.
    """
    
    if ref_space is not None:
        # Compute the reference space diameter
        ref_diameter = find_maximal_distance(ref_space)
        
        # Compute the current space diameter
        current_diameter = find_maximal_distance(matrix)
        
        # Normalize the matrix (elongation or dilatation of the current space)
        matrix = matrix * current_diameter / ref_diameter
    
    matrix1 = matrix[idx_samples_cluster1,:]
    matrix2 = matrix[idx_samples_cluster2,:]
    
    k = min(k, len(matrix1), len(matrix2))
    
    kmeans1 = KMeans(n_clusters=k)
    labels1 = kmeans1.fit_predict(matrix1)   
    
    kmeans2 = KMeans(n_clusters=k) 
    labels2 = kmeans2.fit_predict(matrix2)
    
    # Compute the cluster centroids
    centroids1 = kmeans1.cluster_centers_
    centroids2 = kmeans2.cluster_centers_
    
    # Compute the mass of each cluster

    # Count the number of points in each cluster using np.bincount
    counts1 = np.bincount(labels1)

    # Compute the mass of each cluster
    total_points1 = len(labels1)
    masses1 = counts1 / total_points1

    # Repeat the same process for the second set of points
    counts2 = np.bincount(labels2)
    total_points2 = len(labels2)
    masses2 = counts2 / total_points2
    
    # Compute the cost matrix using squared Euclidean distance
    cost_matrix = ot.dist(centroids1, centroids2, metric='sqeuclidean')

    # Calculate the Wasserstein distance
    print(f"Shapes for input into wasserstein distance computation: masses1: {masses1.shape}, masses2: {masses2.shape}, cost: {cost_matrix.shape}")
    try:
        wasserstein_distance = ot.emd2(masses1, masses2, cost_matrix)
    except:
        print("Error in computing Wasserstein distance. Returning NaN.")
        wasserstein_distance = np.nan
    print(f"Wasserstein distance: {wasserstein_distance}")
    
    return wasserstein_distance
    

        
    
    
def find_maximal_distance(matrix, distance='euclidean', chunk_size=1000):
    """
    Find the maximal pairwise distance in a given matrix.
    This function computes the maximal distance between any two points in the 
    given matrix. If the matrix is large, it processes the matrix in chunks 
    to avoid memory issues.
    Parameters:
    matrix (numpy.ndarray): A 2D array where each row represents a point in space.
    distance (str, optional): The distance metric to use. Default is 'euclidean'.
    chunk_size (int, optional): The size of the chunks to process at a time. Default is 1000.
    Returns:
    float: The maximal pairwise distance found in the matrix.
    """
    
    
    
    if matrix.shape[0] > chunk_size:
        n = matrix.shape[0]
        max_distance = 0
        
        for i in range(0, n, chunk_size):
            for j in range(i, n, chunk_size):
                chunk1 = matrix[i:i+chunk_size]
                chunk2 = matrix[j:j+chunk_size]
                
                if i == j:
                    distances = pdist(chunk1, distance)
                else:
                    distances = squareform(pdist(np.vstack([chunk1, chunk2]), distance))[:len(chunk1), len(chunk1):]
                
                max_distance = max(max_distance, np.max(distances))
    
    else:
        # Compute pairwise Euclidean distances
        distances = pdist(matrix, distance)
        
        # Find the maximum distance
        max_distance = np.max(distances)
        
    return max_distance

    
def calculate_mahalanobis_distance(point, data, both_sided=False):
    # Calculate the mean vector
    mean_vector = np.mean(data, axis=0)
    
    # Calculate the covariance matrix
    cov_matrix = np.cov(data, rowvar=False)
    
    # Calculate the inverse of the covariance matrix
    inv_cov_matrix = np.linalg.inv(cov_matrix)
    
    # Calculate the Mahalanobis distance using scipy's function
    mahalanobis_distance = mahalanobis(point, mean_vector, inv_cov_matrix)
    
    if both_sided:
        diff_vector = point - mean_vector
        return mahalanobis_distance * np.sign(diff_vector)
    else:   
        return mahalanobis_distance
    
    
def calculate_z_score(point, data):
    
    # Calculate the mean vector
    mean_vector = np.mean(data, axis=0)
    
    # Consider the mean vector as a distribution and compute its mean and std
    mean = np.mean(mean_vector)
    std = np.std(mean_vector)
    
    # Calculate the z-score
    z_score = (point.mean() - mean) / ((np.sqrt(std**2 + point.std()**2) / len(point))) 
    
    return z_score




