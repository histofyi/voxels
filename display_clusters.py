

import json
from functions.structures import load_pdb_file_to_pymol

structure_info = json.load(open('output/structure_information/all.json', 'r'))



voxel_map_hash = 'e91d9bdc62da8457549cfbeed4c2b0aa'

used_voxels = json.load(open(f"output/voxel_sets/{voxel_map_hash}/used_voxels.json", 'r'))

clusters = json.load(open(f"output/clusters/{voxel_map_hash}__3__5_6.json", 'r'))

cluster_numbers = sorted([int(cluster) for cluster in clusters.keys() if cluster != '0'])

print (cluster_numbers)
structure_info = json.load(open('output/structure_information/all.json', 'r'))

for cluster_number in cluster_numbers:
    pdb_code = clusters[str(cluster_number)][0]
    print (f"Cluster {cluster_number}")
    load_pdb_file_to_pymol(pdb_code, structure_info, name=f"cluster_{cluster_number}_{pdb_code}") 