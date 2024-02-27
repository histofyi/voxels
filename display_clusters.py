from pymol import cmd

import json

voxel_map_hash = 'e91d9bdc62da8457549cfbeed4c2b0aa'

used_voxels = json.load(open(f"output/voxel_sets/{voxel_map_hash}/used_voxels.json", 'r'))

clusters = json.load(open(f"output/clusters/{voxel_map_hash}__3__5_6.json", 'r'))

cluster_number = 15






for pdb_code in clusters[str(cluster_number)]:
    cmd.load(f"structures/{pdb_code}_peptide.pdb")
