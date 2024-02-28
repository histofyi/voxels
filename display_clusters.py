from pymol import cmd

import json


structure_info = json.load(open('output/structure_information/all.json', 'r'))

def load_structure(pdb_code):
    print (f"{pdb_code}: {structure_info['structures'][pdb_code]['allele']} binding {structure_info['structures'][pdb_code]['peptide']} at {structure_info['structures'][pdb_code]['resolution']}Å resolution")
    cmd.load(f"structures/{pdb_code}_peptide.pdb")

voxel_map_hash = 'e91d9bdc62da8457549cfbeed4c2b0aa'

used_voxels = json.load(open(f"output/voxel_sets/{voxel_map_hash}/used_voxels.json", 'r'))

clusters = json.load(open(f"output/clusters/{voxel_map_hash}__3__4_5_6.json", 'r'))

cluster_numbers = [35]

#filter = 'hla_b_07_02'
filter = None

structure_info = json.load(open('output/structure_information/all.json', 'r'))


print ("Showing structures for cluster(s):", cluster_numbers)

for cluster_number in cluster_numbers:
    for pdb_code in clusters[str(cluster_number)]:
        if not filter:
            load_structure(pdb_code)
        elif filter in structure_info['structures'][pdb_code]['allele_slug']:
            load_structure(pdb_code)
