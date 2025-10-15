import random
import torch

class KMeans_gpu(object):
    def __init__(
        self,
        n_clusters=8,
        n_init=10,
        init="k-means++",
        max_iter=300,
        tol=1e-4,
        random_state=0,
        verbose=True,
        method='kmeans' # can also be kmedoids
    ):
        self.n_clusters = n_clusters
        self.n_init = n_init
        self.init = init
        self.max_iter = max_iter
        self.tol = tol
        self.random_state = random_state
        
        self.verbose = verbose
        self.labels_ = None
        self.cluster_centers_ = None
        self.method = method
        assert method in ['kmeans', 'kmedoids']
    def pairwise_dists(self, X, Y):
        
        return torch.cdist(X, Y, p=2)
    
    def init_centroids(self, X):
        if self.verbose:
            print("-- init centroids --", flush=True)
        random.seed(self.random_state)

        if self.init == "random":
            centroids = X[random.sample(range(X.shape[0]), self.n_clusters), :]

        elif self.init == "k-means++":

            centroids = torch.zeros(
                (self.n_clusters, X.shape[-1]), dtype=X.dtype, device=X.device
            )
            indices = set(range(X.shape[0]))

            for i in range(self.n_clusters):

                if i == 0:
                    idx = random.sample(range(X.shape[0]), self.n_clusters)[0]
                    centroids[i] = X[idx].clone()
                    indices.remove(idx)
                else:
                    
                    l_indices = list(indices)
                    if i == 1: # only one centroid still
                        distances = self.pairwise_dists(X[l_indices], centroids[0][None, :])
                    else:
                        distances = self.pairwise_dists(X[l_indices], centroids[:i])
                    
                    idx = distances.sum(dim=1).argmax().item()
                    centroids[i] = X[l_indices[idx]]
                    indices.remove(l_indices[idx])

        return centroids


    def fit(self, X):
        n_samples = X.shape[0]

        self.inertia = None

        for run_it in range(self.n_init):
            centroids = self.init_centroids(X)

            #for it in tqdm(
            #    range(self.max_iter), desc="fit kmeans (seed = %s)" % run_it
            #):
            
            for it in range(self.max_iter):
                distances = self.pairwise_dists(centroids, X)
                
                labels = distances.argmin(dim=0)
                new_centroids = torch.zeros(
                    (self.n_clusters, X.shape[-1]), dtype=X.dtype, device=X.device
                )
                for i in range(self.n_clusters):
                    indices = torch.where(labels == i)[0]
                    if len(indices) > 0:
                        if self.method == 'kmeans':
                            new_centroids[i, :] = X[indices].mean(dim=0)
                        elif self.method == 'kmedoids':
                            new_centroids[i, :] = X[indices].median(dim=0).values
                        
                    else:  # handle empty cluster

                        new_centroids[i, :] = X[random.sample(range(n_samples), 1), :]

                diff_distances = ((centroids - new_centroids) ** 2).sum().item()
                
                if self.verbose:
                    print(
                        f"Seed : {run_it} / step: {it} / diff_distances: {diff_distances}", flush=True
                    )
                centroids = new_centroids.clone()
                if diff_distances < self.tol:
                    break

            distances = self.pairwise_dists(centroids, X)
            
            labels = distances.argmin(dim=0)

            inertia = 0.0
            for i in range(self.n_clusters):
                indices = torch.where(labels == i)[0]
                inertia += distances[i, indices].sum().item()
            
            if (self.inertia == None) or (inertia < self.inertia):
                self.inertia = inertia
                self.labels_ = labels.clone()
                self.cluster_centers_ = centroids.clone()

            if self.verbose:
                print("Iteration: {} - Best Inertia: {}".format(run_it, self.inertia), flush=True)

    def fit_predict(self, X):
        self.fit(X)
        return self.labels_

    def fit_transform(self, X):
        self.fit(X)
        return self.transform(X)

    def predict(self, X):
        distances = self.transform(X)
        return distances.argmin(dim=0)

    def transform(self, X):
        distances = self.pairwise_dists(self.cluster_centers_, X)
        return distances
    