from typing import Dict
import os
import requests
import json

from biopandas.pdb import PandasPdb


import constants
import functions


### Functions used by this script ###

def load_pdb_file_to_dataframe(pdb_code:str, domain:str) -> Dict:
    """
    Loads a PDB file from the Histo website and returns a dataframe containing the data.

    Args:
        pdb_code (str): The PDB code of the structure to be loaded, e.g. '1hhk'.
        domain (str): The domain of the structure to be loaded, e.g. 'peptide'.

    """
    # lowercase the pdb code for both the filename and the URL
    pdb_code = pdb_code.lower()

    # create the filename and the URL
    filename = f"structures/{pdb_code}_{domain}.pdb"
    url = f"{constants.base_url}/{pdb_code}_1_{domain}.pdb"

    # if the file exists, read it into a dataframe
    if os.path.exists(filename):
        structure_df = PandasPdb().read_pdb(filename)

    # otherwise, download the file and then read it into a dataframe
    else:
    
        r = requests.get(url)
        if r.status_code == 200:
            structure_data = r.text

            with open(filename, 'w') as filehandle:
                filehandle.write(structure_data)

            structure_df = PandasPdb().read_pdb(filename)
        # if the request fails, set structure_df None and print an error message
        else:
            print (f"PDB code {pdb_code} has failed with status {r.status_code}")
            structure_df = None

    # if the structure_df is not None, return the dataframe
    if structure_df:        
        return structure_df.df
    else:
        return None


def find_voxels_for_structure(peptide_df, voxels, voxel_size:int=1):

    # initialise the lists and dictionaries
    coordinates_to_check = []
    sequence = []
    positions = []
    structure_voxels = {}

    target_atom_name = 'CA', # this is hard-coded for now, but should be a parameter in the future    

    # iterate through the rows of the dataframe to generate a list of coordinates to check
    for index, row in peptide_df['ATOM'].iterrows():

        # first, extract the relevant data from the row
        atom_name = row['atom_name']
        residue_name = row['residue_name']
        chain_id = row['chain_id']
        x, y, z = row['x_coord'], row['y_coord'], row['z_coord']
        residue_number = row['residue_number']
        
        # if the atom name is the target atom name, add the residue name to the sequence list and the coordinates to the coordinates_to_check list
        if atom_name == target_atom_name:
            if residue_number not in positions:
                sequence.append(row['residue_name'])
                coordinates_to_check.append((x,y,z))
                positions.append(residue_number)
    
    # set the variable for the position number to 1
    p = 1

    # iterate through the coordinates to check and the voxels to find the voxels that contain the coordinates
    for coordinate in coordinates_to_check:
        for voxel in voxels:
            # if the atom is inside a specific voxel, add the voxel and the atom information to the structure_voxels dictionary
            if functions.is_coordinate_inside_voxel(coordinate[0],coordinate[1],coordinate[2], voxels[voxel]['start'], voxel_size):
                structure_voxels[str(p)] = {
                    'position': p,
                    'voxel_label': voxel,
                    'atom_name': target_atom_name,
                    'atom_coordinates': coordinate,
                    'voxel_start': voxels[voxel]['start'],
                    'voxel_centre': voxels[voxel]['centre'],
                    'residue': sequence[p-1]
                }
        p += 1
    # return the structure_voxels dictionary
    return structure_voxels



### Body of the script ###

# create the voxel grid, we only need to do this once, so it's outside the loop
voxel_grid = functions.create_voxel_grid(constants.centre_of_mass, constants.box_xyz, constants.voxel_size, range_offset=1, y_offset=8)

# iterate through the test PDB codes
for pdb_code in constants.test_pdb_codes:

    # load the peptide dataframe for the PDB code
    dataframe = load_pdb_file_to_dataframe(pdb_code, 'peptide')
    
    if dataframe:
    # if the dataframe is not None, find the voxels for the structure and write them to a file
        structure_voxels = find_voxels_for_structure(dataframe, voxel_grid['voxels'], voxel_size=constants.voxel_size)

        filename = f"output/voxels/{pdb_code}.json"

        with open(filename, 'w') as filehandle:
            json.dump(structure_voxels, filehandle, indent=4)
        print (f"Voxels for {pdb_code} have been written to {filename}")

