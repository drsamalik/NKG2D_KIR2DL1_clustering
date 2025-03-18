import numpy as np
import random
import copy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from find_cluster import get_clusters
# Rotate a cluster using your method
def rotate_cluster(cluster, angle):
    min_y, min_x = min(y for y, x in cluster), min(x for y, x in cluster)
    cluster_shifted = [(y - min_y, x - min_x) for y, x in cluster]
    max_dim = max(max(y for y, x in cluster_shifted), max(x for y, x in cluster_shifted))
    if angle == 0:
        return cluster
    elif angle == 90:
        return [(x, max_dim - y) for y, x in cluster_shifted]
    elif angle == 180:
        return [(max_dim-y, max_dim-x) for y, x in cluster_shifted]
    elif angle == 270:
        return [(max_dim-x, y) for y, x in cluster_shifted]
# Calculate center of mass for a cluster
def calculate_com(cluster):
    ys, xs = zip(*cluster)
    return round(np.mean(ys)), round(np.mean(xs))
# Translate a cluster based on CoM
def translate_cluster(cluster, target_com):
    current_com = calculate_com(cluster)
    dy, dx = target_com[0] - current_com[0], target_com[1] - current_com[1]
    return [(y + dy, x + dx) for y, x in cluster]
# Check placement validity (side-wise non-touching)
def can_place_cluster(new_matrix, cluster):
    rows, cols = new_matrix.shape
    for y, x in cluster:
        if not (0 <= y < rows and 0 <= x < cols):
            return False
        if new_matrix[y, x] != 0:
            return False
        # Check neighbors (side-wise only)
        for dy, dx in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            ny, nx = y + dy, x + dx
            if 0 <= ny < rows and 0 <= nx < cols and new_matrix[ny, nx] != 0:
                return False
    return True
def random_relocation(clusters, total_clusters, new_size, nucleated_matrix=None, init=False):
    if nucleated_matrix is None:
        new_matrix = np.zeros(new_size)
        cluster_id, num_of_placed_clustor = 1, 0
    else:
        new_matrix = nucleated_matrix
        cluster_id, num_of_placed_clustor = 2, len(get_clusters(new_matrix))
    rows, cols = new_matrix.shape
    max_iter = 1000
    iter = 0
    while  num_of_placed_clustor < total_clusters and iter < max_iter:
        i = random.randint(0, len(clusters)-1) if not init else iter%len(clusters)
        cluster = clusters[i]
        placed = False
        attempts = 0
        max_attempts = 100
        while not placed and attempts < max_attempts:
            attempts += 1
            angle = random.choice([0, 90, 180, 270])
            rotated_cluster = rotate_cluster(cluster, angle)
            target_com = (random.randint(0, rows - 1), random.randint(0, cols - 1))
            translated_cluster = translate_cluster(rotated_cluster, target_com)
            if can_place_cluster(new_matrix, translated_cluster):
                for y, x in translated_cluster:
                    new_matrix[y, x] = cluster_id
                placed = True
        iter +=1
        bi_n_matrix = copy.deepcopy(np.array(new_matrix))
        bi_n_matrix[bi_n_matrix>=1] = 1
        num_of_placed_clustor = len(get_clusters(bi_n_matrix))
    return new_matrix