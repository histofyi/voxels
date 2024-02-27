from typing import List
from hdbscan import HDBSCAN


import json


def tensorize(voxel_labels:List[str]) -> List[List[int]]:
    tensorized = [[int(value) for value in voxel_label.split('_')] for voxel_label in voxel_labels]
    return tensorized

# Load the file containing the position of the used voxels for each structure

voxel_map_hash = 'e91d9bdc62da8457549cfbeed4c2b0aa'

used_voxels = json.load(open(f"output/voxel_sets/{voxel_map_hash}/used_voxels.json", 'r'))

labels = list(used_voxels.keys())

all_positions_data = [tensorize(voxel_labels) for voxel_labels in used_voxels.values()]

positions_for_clustering = [5, 6]
min_cluster_size = 3

clustering_data = []

for row in all_positions_data:
    p = 1
    processed_row = []
    for position in row:
        if p in positions_for_clustering:

            for item in position:
                processed_row.append(item)
            
        p += 1
    clustering_data.append(processed_row)

hdb = HDBSCAN(min_cluster_size=min_cluster_size)
hdb.fit(clustering_data)

print (f"Clustering at positions {positions_for_clustering} with a minimum cluster size of {min_cluster_size} gives {len(set(hdb.labels_)) - 1} meaningful clusters and a group of outliers.")

clusters = {}

for i in range(len(hdb.labels_)):
    cluster = int(hdb.labels_[i]) + 1
    if cluster not in clusters:
        clusters[cluster] = []
    clusters[cluster].append(labels[i])


for cluster in clusters:
    print (cluster, len(clusters[cluster]))


positions_for_clustering = [str(p) for p in positions_for_clustering]

json.dump(clusters, open(f"output/clusters/{voxel_map_hash}__{min_cluster_size}__{'_'.join(positions_for_clustering)}.json", 'w'), indent=4)
