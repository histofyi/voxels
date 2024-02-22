from typing import Dict
import os
import requests
import json

from biopandas.pdb import PandasPdb


import constants
import functions

pdb_code = '1hhk'


def load_pdb_file_to_dataframe(pdb_code:str, domain:str) -> Dict:
    filename = f"structures/{pdb_code}_{domain}.pdb"
    url = f"{constants.base_url}/{pdb_code}_1_{domain}.pdb"

    if os.path.exists(filename):
        structure_df = PandasPdb().read_pdb(filename)
    else:
    
        r = requests.get(url)
        if r.status_code == 200:
            structure_data = r.text

            with open(filename, 'w') as filehandle:
                filehandle.write(structure_data)

            structure_df = PandasPdb().read_pdb(filename)
        else:
            print (f"PDB code {pdb_code} has failed with status {r.status_code}")
            structure_df = None
    return structure_df.df


def find_voxels_for_structure(peptide_df, voxels, voxel_size:int=1):

    coordinates_to_check = []
    sequence = []
    positions = []
    for index, row in peptide_df['ATOM'].iterrows():
        atom_name = row['atom_name']
        residue_name = row['residue_name']
        chain_id = row['chain_id']
        x, y, z = row['x_coord'], row['y_coord'], row['z_coord']
        residue_number = row['residue_number']
        if atom_name == 'CA':
            if residue_number not in positions:
                sequence.append(row['residue_name'])
                coordinates_to_check.append((x,y,z))
                positions.append(residue_number)
    p = 1
    structure_voxels = {}

    for coordinate in coordinates_to_check:
        for voxel in voxels:
            if functions.is_coordinate_inside_voxel(coordinate[0],coordinate[1],coordinate[2], voxels[voxel]['start'], voxel_size):
                structure_voxels[str(p)] = {
                    'position': p,
                    'voxel_label': voxel,
                    'atom_coordinates': coordinate,
                    'voxel_start': voxels[voxel]['start'],
                    'voxel_centre': voxels[voxel]['centre'],
                    'residue': sequence[p-1]
                }
        p += 1
    return structure_voxels


dataframe = load_pdb_file_to_dataframe(pdb_code, 'peptide')
voxel_grid = functions.create_voxel_grid(constants.centre_of_mass, constants.box_xyz, constants.voxel_size, range_offset=1, y_offset=8)


structure_voxels = find_voxels_for_structure(dataframe, voxel_grid['voxels'], voxel_size=constants.voxel_size)

filename = f"output/voxels/{pdb_code}.json"

with open(filename, 'w') as filehandle:
    json.dump(structure_voxels, filehandle, indent=4)
print (f"Voxels for {pdb_code} have been written to {filename}")

