from pymol import cmd

import json
import sys

from time import sleep

print (sys.argv)


structure_info = json.load(open('output/structure_information/all.json', 'r'))

def load_structure(pdb_code):
    print (f"{pdb_code}: {structure_info['structures'][pdb_code]['allele']} binding {structure_info['structures'][pdb_code]['peptide']} at {structure_info['structures'][pdb_code]['resolution']}Å resolution")
    cmd.load(f"structures/{pdb_code}_peptide.pdb")

voxel_map_hash = 'e91d9bdc62da8457549cfbeed4c2b0aa'
min_cluster_size = 3


clusters = json.load(open(f"output/clusters/{voxel_map_hash}/cluster_size_{min_cluster_size}/clustering_4_5_6.json", 'r'))


structure_info = json.load(open('output/structure_information/all.json', 'r'))

cluster_numbers = sorted([int(cluster) for cluster in clusters['clusters'].keys()])

for cluster_number in cluster_numbers:

    i = 0
    for pdb_code in clusters['clusters'][str(cluster_number)]:
        load_structure(pdb_code)
        
    cmd.remove('hydro')            
    cmd.show('sticks', 'resi 2 and not (name c,n)')

    cmd.show('sticks', 'resi 9 and not (name c,n)')

    for position in [5,6]:
        residue_selection = f"resi {position} and not (name c,n)"
        cmd.show('sticks', residue_selection)
        cmd.png(f"output/images/{voxel_map_hash}/cluster_size_{min_cluster_size}/clustering_4_5_6/cluster_{cluster_number}__{position}.png", width=3000, height=3000, ray=1)

        cmd.hide('sticks', residue_selection)

    cmd.delete('all')