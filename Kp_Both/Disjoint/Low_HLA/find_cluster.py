import numpy as np
import matplotlib.pyplot as plt
import random
import matplotlib.cm as cm
import pandas as pd
import copy
def NKG_Cluster(dim,in_pxl,current_path):
    binary_intensity = np.loadtxt(f'{current_path}/NKG2D_input_file.txt').astype(np.int64)
    clusters = get_clusters(binary_intensity)
    num_new_clustor = cluster_density(clusters,binary_intensity,dim)
    return clusters, num_new_clustor, binary_intensity
def KIR_Cluster(dim,in_pxl,current_path):
    binary_intensity = np.loadtxt(f'{current_path}/KIR2DL1_input_file.txt').astype(np.int64)
    clusters = get_clusters(binary_intensity)
    num_new_clustor = cluster_density(clusters,binary_intensity,dim)
    return clusters, num_new_clustor,binary_intensity

def get_clusters(binary_matrix):
    clusters = []
    visited = np.zeros_like(binary_matrix, dtype=bool)
    for y in range(binary_matrix.shape[0]):
        for x in range(binary_matrix.shape[1]):
            if binary_matrix[y, x] == 1 and not visited[y, x]:
                # Found a new cluster
                cluster = []
                # Use BFS or DFS to find the entire cluster
                stack = [(y, x)]
                while stack:
                    cy, cx = stack.pop()
                    if visited[cy, cx]:
                        continue
                    visited[cy, cx] = True
                    cluster.append((cy, cx))
                    # Check neighbors
                    for dy, dx in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                        ny, nx = cy + dy, cx + dx
                        if 0 <= ny < binary_matrix.shape[0] and 0 <= nx < binary_matrix.shape[1]:
                            if binary_matrix[ny, nx] == 1 and not visited[ny, nx]:
                                stack.append((ny, nx))
                clusters.append(cluster)
    return clusters

def cluster_density(clusters,binary_intensity,dim):
    clustor_density = len(clusters)/(binary_intensity.shape[1]*binary_intensity.shape[0])
    num_new_clustor = int(clustor_density*dim[0]*dim[1])
    return num_new_clustor
