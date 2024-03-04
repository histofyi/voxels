from pymol import cmd

import json
import sys

print (sys.argv)

cluster_numbers = []

for i in sys.argv[1:]: 
    if i.isdigit():
        cluster_numbers.append(int(i))

if len(cluster_numbers) == 0:
    cluster_numbers = [1]

structure_info = json.load(open('output/structure_information/all.json', 'r'))

def load_structure(pdb_code):
    print (f"{pdb_code}: {structure_info['structures'][pdb_code]['allele']} binding {structure_info['structures'][pdb_code]['peptide']} at {structure_info['structures'][pdb_code]['resolution']}Å resolution")
    cmd.load(f"structures/{pdb_code}_peptide.pdb")

voxel_map_hash = 'e91d9bdc62da8457549cfbeed4c2b0aa'
min_cluster_size = 3

#used_voxels = json.load(open(f"output/voxel_sets/{voxel_map_hash}/used_voxels.json", 'r'))

clusters = json.load(open(f"output/clusters/{voxel_map_hash}/cluster_size_{min_cluster_size}/clustering_4_5_6.json", 'r'))


#filter = 'hla_b_07_02'
filter = None

structure_info = json.load(open('output/structure_information/all.json', 'r'))


print ("Showing structures for cluster(s):", cluster_numbers)

for cluster_number in cluster_numbers:
    print (f"Cluster {cluster_number}")
    i = 0
    for pdb_code in clusters['clusters'][str(cluster_number)]:
        if not filter:
            load_structure(pdb_code)
        elif filter in structure_info['structures'][pdb_code]['allele_slug']:
            load_structure(pdb_code)


