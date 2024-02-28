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

positions_for_clustering = [4, 5, 6]
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


clusters = {}

for i in range(len(hdb.labels_)):
    cluster = int(hdb.labels_[i]) + 1
    if cluster not in clusters:
        clusters[cluster] = []
    clusters[cluster].append(labels[i])

print (f"Clustering at positions {positions_for_clustering} with a minimum cluster size of {min_cluster_size} gives {len(set(hdb.labels_)) - 1} meaningful clusters and a group of outliers with {len(clusters[0])} members.")

print ("The clusters are as follows:")

structure_info = json.load(open('output/structure_information/all.json', 'r'))

alleles = {}

for key, value in sorted(clusters.items(), key=lambda item: int(item[0])):
    print (f"\nCluster {key}: {len(value)} members")

    

    if key != 0:
        cluster_alleles = {}
        cluster_used_voxels = {}
        for i in range(1,10):
            cluster_used_voxels[i] = {}
        print (value)
        for member in value:
            
            this_used_voxels = used_voxels[member]

            i = 1
            for voxel in this_used_voxels:
                if not voxel in cluster_used_voxels:
                    cluster_used_voxels[i][voxel] = 0
                cluster_used_voxels[i][voxel] += 1
                i += 1

            this_allele = structure_info['structures'][member]['allele']
            if this_allele not in alleles:
                alleles[this_allele] = []
            if key not in alleles[this_allele]:
                alleles[this_allele].append(key)
            if this_allele not in cluster_alleles:
                cluster_alleles[this_allele] = 0
            cluster_alleles[this_allele] += 1
        print (f"Allele: {cluster_alleles}")
        print (f"\nVoxels:")
        for i in range(1,10):
            print (f"Position {i}: {len(cluster_used_voxels[i])} voxels")

print ("\n\nAllele cluster distribution:") 
for allele in alleles:
    print (f"{allele}: {len(alleles[allele])} clusters {alleles[allele]}")     



positions_for_clustering = [str(p) for p in positions_for_clustering]

filename = f"output/clusters/{voxel_map_hash}__{min_cluster_size}__{'_'.join(positions_for_clustering)}.json"
json.dump(clusters, open(filename, 'w'), indent=4)

print (f"The data for this run can be found in {filename}")

