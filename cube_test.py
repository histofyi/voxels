from typing import List, Dict
import copy

import helpers
import constants


import json

from pymol import cmd, cgo
from colorsys import hls_to_rgb

def cgo_cube(x, y, z, r):
    r *= 3 ** -.5
    return [
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., 0., 1.,
        cgo.VERTEX, x + r, y + r, z + r,
        cgo.VERTEX, x + r, y - r, z + r,
        cgo.VERTEX, x - r, y + r, z + r,
        cgo.VERTEX, x - r, y - r, z + r,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 1., 0., 0.,
        cgo.VERTEX, x + r, y - r, z - r,
        cgo.VERTEX, x + r, y + r, z - r,
        cgo.VERTEX, x + r, y - r, z + r,
        cgo.VERTEX, x + r, y + r, z + r,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., 1., 0.,
        cgo.VERTEX, x + r, y + r, z - r,
        cgo.VERTEX, x - r, y + r, z - r,
        cgo.VERTEX, x + r, y + r, z + r,
        cgo.VERTEX, x - r, y + r, z + r,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., 0., -1.,
        cgo.VERTEX, x - r, y - r, z - r,
        cgo.VERTEX, x - r, y + r, z - r,
        cgo.VERTEX, x + r, y - r, z - r,
        cgo.VERTEX, x + r, y + r, z - r,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, -1., 0., 0.,
        cgo.VERTEX, x - r, y + r, z + r,
        cgo.VERTEX, x - r, y - r, z + r,
        cgo.VERTEX, x - r, y + r, z - r,
        cgo.VERTEX, x - r, y - r, z - r,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., -1., 0.,
        cgo.VERTEX, x - r, y - r, z + r,
        cgo.VERTEX, x + r, y - r, z + r,
        cgo.VERTEX, x - r, y - r, z - r,
        cgo.VERTEX, x + r, y - r, z - r,
        cgo.END,
    ]



def display_pseudoatom(name:str, pos:List[float], color:str='gray80', sphere_scale:float=0.3, opacity:int=1):
    cmd.pseudoatom(name, pos=pos)
    cmd.show('spheres',name)
    cmd.color(color,name)
    cmd.set ("sphere_scale", sphere_scale, name)
    cmd.set ("sphere_transparency", opacity, name)


def load_mhc_abd(pdb_code:str):
    abd_id = f"{pdb_code}_1_abd"
    histo_abd_url = f"{constants.base_url}/{abd_id}.cif"
    cmd.load(histo_abd_url)
    cmd.color('gray40', abd_id)


def load_peptide(pdb_code: str):
    peptide_id = f"{pdb_code}_1_peptide"
    histo_peptide_url = f"{constants.base_url}/{peptide_id}.cif"
    cmd.load(histo_peptide_url)


voxel_map_hash = 'e91d9bdc62da8457549cfbeed4c2b0aa'

voxel_grid = json.load(open(f"output/voxel_sets/{voxel_map_hash}/voxel_set.json", 'r'))


# Load the MHC and peptide
load_mhc_abd(constants.canonical_pdb_code)
load_peptide(constants.canonical_pdb_code)


# Display the centre of mass for the canonical class I molecule and the centre of the voxel grid
display_pseudoatom('centre_of_box', constants.centre_of_mass, color='white', opacity=0)

cgo_cube(*constants.centre_of_mass, 10)

# Load the file containing the position of the used voxels and their frequency
position_voxels = json.load(open(f"output/voxel_sets/{voxel_map_hash}/position_voxels.json", 'r'))
voxel_grid = json.load(open(f"output/voxel_sets/{voxel_map_hash}/voxel_set.json", 'r'))

def cubes(selection='all', name='', state=0, scale=0.5, color='green', opacity=0, _func=cgo_cube):
    '''
    DESCRIPTION

    Create a cube representation CGO for all atoms in selection.

    ARGUMENTS

    selection = string: atom selection {default: all}

    name = string: name of CGO object to create

    state = int: object state {default: 0 = all states}

    scale = float: scaling factor. If scale=1.0, the corners of the cube will
    be on the VDW surface of the atom {default: 0.5}

    atomcolors = 0/1: use atom colors (cannot be changed), otherwise
    apply one color to the object (can be changed with color command)
    {default: 1}

    SEE ALSO

    tetrahedra
    '''
    if not name:
        name = 'cube_me'
    state, scale, color, opacity = int(state), float(scale), str(color), float(opacity)
    if state < 0:
        state = cmd.get_setting_int('state')
    states = [state] if state else list(range(1, cmd.count_states(selection) + 1))

    def callback(x, y, z, vdw, color):
        vdw = 1.0
        if color:
            # this is run per atom
            obj.append(cgo.COLOR)
            # This is the color of the cube
            obj.extend((1.0, 1.0, 1.0))
        obj.extend(_func(x, y, z, vdw * scale))

    space = {'xcb': callback}
    for state in states:
        obj = []
        cmd.iterate_state(state, selection,
                          'xcb(x, y, z, vdw, color)', space=space)
        cmd.load_cgo(obj, name, state)
    cmd.color(color, name)

    cmd.set("cgo_transparency", opacity, name)
    cmd.delete(selection)


cmd.extend('cubes', cubes)

def rainbow_color_steps(steps:int, end=6/8):
    raw_colors = [hls_to_rgb(end * i/(steps-1), 0.5, 1) for i in range(steps)]
    return [[int(175 * x) for x in color] for color in raw_colors]

colour_steps = rainbow_color_steps(101)

for i, colour in enumerate(colour_steps):
    cmd.set_color(f'occupancy_{i}', colour)


voxels = []
occupancy_counts = []
centres = []

for position in position_voxels:
    structure_count = 0
    for voxel in position_voxels[position]:
        voxels.append(voxel)
        occupancy_counts.append(position_voxels[position][voxel]['count'])
        structure_count += position_voxels[position][voxel]['count']
        centres.append(voxel_grid['voxels'][voxel]['centre'])

max_occupancy = max(occupancy_counts)
min_occupancy = min(occupancy_counts)

print (max_occupancy, min_occupancy)

for i, voxel in enumerate(voxels):
        print (voxel)
        print (occupancy_counts[i])
        normalised_count = (occupancy_counts[i] - min_occupancy) / (max_occupancy - min_occupancy) * structure_count
        print (f"Normalised count {normalised_count}")
        occupancy = normalised_count / structure_count
        percentage_occupancy = round(occupancy * 100)
        opacity = (1.0 - round(occupancy, 1)) - 0.4
        print (f"Percentage occupancy {percentage_occupancy}")
        voxel_centre = centres[i]
        print (f"Opacity {opacity}")
        
        
        display_pseudoatom(voxel, voxel_centre, color='white', opacity=opacity)
        cubes(selection=voxel, name=f'voxel_{voxel}', color=f'occupancy_{percentage_occupancy}', opacity=opacity)

cmd.reset()





