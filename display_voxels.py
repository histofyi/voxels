from typing import List, Dict
import copy

import functions
import constants







def load_mhc_abd(pdb_code:str):
    abd_id = f"{pdb_code}_1_abd"
    histo_abd_url = f"{constants.base_url}/{abd_id}.cif"
    cmd.load(histo_abd_url)
    cmd.color('gray40', abd_id)


def load_peptide(pdb_code: str):
    peptide_id = f"{pdb_code}_1_peptide"
    histo_peptide_url = f"{constants.base_url}/{peptide_id}.cif"
    cmd.load(histo_peptide_url)


def display_pseudoatom(name:str, pos:List[float], color:str='gray80', sphere_scale:float=0.3):
    cmd.pseudoatom(name, pos=pos)
    cmd.show('spheres',name)
    cmd.color(color,name)
    cmd.set ("sphere_scale", sphere_scale, name)


def show_voxel(voxel_coords:List, decimation_factor:int):
    full_voxel_coords = copy.deepcopy(voxel_coords)
    show_voxel_bool = False
    if 0 in voxel_coords:
        voxel_coords.remove(0)
        if voxel_coords[0] % decimation_factor == 0 and voxel_coords[1] % decimation_factor == 0:
                show_voxel_bool = True
    return show_voxel_bool
    

def display_voxel_box(centre_of_mass:List[float], box_xyz:List[float], voxel_dict:Dict):

    display_pseudoatom('first_voxel', voxel_dict['voxels'][voxel_dict['labels'][0]]['start'], color='red')
    display_pseudoatom('last_voxel', voxel_dict['voxels'][voxel_dict['labels'][voxel_dict['params']['count'] - 1]]['start'], color='blue')

    for key in voxel_dict['voxels']:
        x, y, z = key.split('_')
        x, y, z = int(x), int(y), int(z)

        if show_voxel([x, y, z], 4):
            display_pseudoatom(f'voxel_grid', voxel_dict['voxels'][key]['start'], color='gray50')
    
    cmd.set ('sphere_transparency', 1, 'voxel_grid')


load_mhc_abd(constants.canonical_pdb_code)
load_peptide(constants.canonical_pdb_code)

voxel_grid = functions.create_voxel_grid(constants.centre_of_mass, constants.box_xyz, constants.voxel_size, range_offset=1, y_offset=8)

display_pseudoatom('centre_of_box', constants.centre_of_mass, color='white')

display_voxel_box(constants.centre_of_mass, constants.box_xyz, voxel_grid)

cmd.reset()