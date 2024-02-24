import json

voxel_grid_hash = 'e91d9bdc62da8457549cfbeed4c2b0aa'

voxel_set_filepath = f"output/voxel_sets/{voxel_grid_hash}"

structure_info = json.load(open("output/structure_information/all.json", 'r'))

used_voxel_set = {}
position_voxel_set = {}

for pdb_code in structure_info['structures']:
    filename = f"{voxel_set_filepath}/{pdb_code}.json"
    voxels = json.load(open(filename, 'r'))['structure_voxels']


    used_voxels = [voxels[position]['voxel_label'] for position in voxels]
    used_voxel_set[pdb_code] = used_voxels

    position = 1
    for voxel_label in used_voxels:
        if position not in position_voxel_set:
            position_voxel_set[position] = {}
        if voxel_label not in position_voxel_set[position]:
            position_voxel_set[position][voxel_label] = {'count':0,'members':[]}
        if pdb_code not in position_voxel_set[position][voxel_label]['members']:
            position_voxel_set[position][voxel_label]['members'].append(pdb_code)
        position_voxel_set[position][voxel_label]['count'] += 1
        position += 1

with open(f"{voxel_set_filepath}/used_voxels.json", 'w') as filehandle:
    json.dump(used_voxel_set, filehandle, indent=4)

with open(f"{voxel_set_filepath}/position_voxels.json", 'w') as filehandle:
    json.dump(position_voxel_set, filehandle, indent=4)

print (f"Voxel usage for {voxel_grid_hash} has been written to {voxel_set_filepath}/used_voxels.json and {voxel_set_filepath}/position_voxels.json")