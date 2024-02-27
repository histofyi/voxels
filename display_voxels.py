from typing import List, Dict
import copy

import functions
import constants


import json




def load_mhc_abd(pdb_code:str):
    abd_id = f"{pdb_code}_1_abd"
    histo_abd_url = f"{constants.base_url}/{abd_id}.cif"
    cmd.load(histo_abd_url)
    cmd.color('gray40', abd_id)


def load_peptide(pdb_code: str):
    peptide_id = f"{pdb_code}_1_peptide"
    histo_peptide_url = f"{constants.base_url}/{peptide_id}.cif"
    cmd.load(histo_peptide_url)


def display_pseudoatom(name:str, pos:List[float], color:str='gray80', sphere_scale:float=0.3, opacity:int=1):
    cmd.pseudoatom(name, pos=pos)
    cmd.show('spheres',name)
    cmd.color(color,name)
    cmd.set ("sphere_scale", sphere_scale, name)
    cmd.set ("sphere_transparency", opacity, name)


def show_axis_voxel(voxel_coords:List, decimation_factor:int):
    full_voxel_coords = copy.deepcopy(voxel_coords)
    show_voxel_bool = False
    axis_color = None
    axis = None
    if 0 in voxel_coords:
        if voxel_coords.index(0) == 0:
            axis_color = 'deepteal'
            axis = 'x'
        elif voxel_coords.index(0) == 1:
            axis_color = 'deeppurple'
            axis = 'y'
        elif voxel_coords.index(0) == 2:
            axis_color = 'chocolate'
            axis = 'z'
        voxel_coords.remove(0)
        if voxel_coords[0] % decimation_factor == 0 and voxel_coords[1] % decimation_factor == 0:
            show_voxel_bool = True
    return show_voxel_bool, axis_color, axis
    

def display_voxel_box(centre_of_mass:List[float], box_xyz:List[float], voxel_dict:Dict):

    display_pseudoatom('first_voxel', voxel_dict['voxels'][voxel_dict['labels'][0]]['start'], color='white')
    display_pseudoatom('last_voxel', voxel_dict['voxels'][voxel_dict['labels'][voxel_dict['params']['voxel_count'] - 1]]['start'], color='gray80')

    for key in voxel_dict['voxels']:
        x, y, z = key.split('_')
        x, y, z = int(x), int(y), int(z)
        show_voxel_bool, axis_color, axis = show_axis_voxel([x, y, z], 4)

        if show_voxel_bool:
            display_pseudoatom(f'voxel_grid_{axis}', voxel_dict['voxels'][key]['start'], color=axis_color)
    for axis in ['x', 'y', 'z']:
        cmd.set ('sphere_transparency', 1, f'voxel_grid_{axis}')

voxel_map_hash = 'e91d9bdc62da8457549cfbeed4c2b0aa'

voxel_grid = json.load(open(f"output/voxel_sets/{voxel_map_hash}/voxel_set.json", 'r'))


# Load the MHC and peptide
load_mhc_abd(constants.canonical_pdb_code)
load_peptide(constants.canonical_pdb_code)


# Display the centre of mass for the canonical class I molecule and the centre of the voxel grid
display_pseudoatom('centre_of_box', constants.centre_of_mass, color='white')

# Display the voxel grid
display_voxel_box(constants.centre_of_mass, constants.box_xyz, voxel_grid)

# Load the file containing the position of the used voxels and their frequency
position_voxels = json.load(open(f"output/voxel_sets/{voxel_map_hash}/position_voxels.json", 'r'))


for position in position_voxels:

    voxel_labels = position_voxels[position].keys()

    voxel_counts = [position_voxels[position][label]['count'] for label in voxel_labels]

    voxel_total = sum(voxel_counts)

    voxel_percentages = [round(count / voxel_total *100, 2) for count in voxel_counts]

    voxel_colours = [f"s{round(percentage *10):03}" for percentage in voxel_percentages]

    voxel_opacity = [0.6 - (int((percentage * 0.1) + 1)) / 10  for percentage in voxel_percentages]
    print (voxel_colours)
    
    i = 0
    for voxel_label in voxel_labels:
        print (voxel_label)
        print (voxel_opacity[i])
        display_pseudoatom(f"p{position}_{voxel_label}", voxel_grid['voxels'][voxel_label]['centre'], color=voxel_colours[i], sphere_scale=0.7, opacity=voxel_opacity[i])


        i += 1

        

#    for position in structure_voxels:
#        display_pseudoatom(f"{pdb_code}_voxels", structure_voxels[position]['voxel_centre'], color='green', sphere_scale=1)
#        display_pseudoatom(f"p{position}_voxels", structure_voxels[position]['voxel_centre'], color='red', sphere_scale=1)

cmd.reset()