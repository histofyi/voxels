from typing import List, Dict
import copy

import helpers
import constants


import json

from pymol import cmd, cgo
from colorsys import hls_to_rgb

'''
Includes the cgo_cube and ammended cubes function from the cubes pymol extension by Thomas Holder

cgo_cube has not been changed
cubes has been changed to include the colouring/transparency of the cubes based on the occupancy of the voxel

Square and Tetrahedra representations

(c) 2013 Thomas Holder

License: BSD-2-Clause
'''


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
    abd_id = "antigen_binding_domain"
    print (f"{pdb_code} is the antigen binding domain")
    histo_abd_url = f"{constants.base_url}/{pdb_code}_1_abd.cif"
    cmd.load(histo_abd_url, abd_id)
    cmd.color('gray40', abd_id)


def load_peptide(pdb_code: str):
    peptide_id = "peptide"
    print (f"{pdb_code} is the peptide")
    histo_peptide_url = f"{constants.base_url}/{pdb_code}_1_peptide.cif"
    cmd.load(histo_peptide_url, peptide_id)
    cmd.hide('everything', peptide_id)
    cmd.color('gray80', peptide_id)
    cmd.show('ribbon', peptide_id)



def display_cleft_bounding_box(com:List, xyz:List, x_offset:int=0, y_offset:int=0, z_offset:int=0, show_spheres=False):

    #TODO automate all of this
    # combinations, this and the bonding feels automateable
    """
        ---
        -+-
        +--
        --+
        ++-
        -++
        +-+
        +++
    """
    xyz = [item/2 for item in xyz]

    print (y_offset)

    lower_right_a1 = [(com[0] - xyz[0]) + x_offset, (com[1] - xyz[1]) + y_offset, (com[2] - xyz[2]) + z_offset]
    upper_right_a1 = [(com[0] - xyz[0]) + x_offset, (com[1] + xyz[1]) + y_offset, (com[2] - xyz[2]) + z_offset]
    lower_left_a1 = [(com[0] + xyz[0]) + x_offset, (com[1] - xyz[1]) + y_offset, (com[2] - xyz[2]) + z_offset]
    upper_left_a1 = [(com[0] + xyz[0]) + x_offset, (com[1] + xyz[1]) + y_offset, (com[2] - xyz[2]) + z_offset]


    lower_right_a2 = [(com[0] - xyz[0]) + x_offset, (com[1] - xyz[1]) + y_offset, (com[2] + xyz[2]) + z_offset]
    upper_right_a2 = [(com[0] - xyz[0]) + x_offset, (com[1] + xyz[1]) + y_offset, (com[2] + xyz[2]) + z_offset]
    lower_left_a2 = [(com[0] + xyz[0]) + x_offset, (com[1] - xyz[1]) + y_offset, (com[2] + xyz[2]) + z_offset]
    upper_left_a2 = [(com[0] + xyz[0]) + x_offset, (com[1] + xyz[1]) + y_offset, (com[2] + xyz[2]) + z_offset]

    cmd.pseudoatom('box', pos=lower_right_a1, name='LRA1')
    cmd.pseudoatom('box', pos=upper_right_a1, name='URA1')
    cmd.pseudoatom('box', pos=lower_left_a1, name='LLA1')
    cmd.pseudoatom('box', pos=upper_left_a1, name='ULA1')
    cmd.pseudoatom('box', pos=lower_right_a2, name='LRA2')
    cmd.pseudoatom('box', pos=upper_right_a2, name='URA2')
    cmd.pseudoatom('box', pos=lower_left_a2, name='LLA2')
    cmd.pseudoatom('box', pos=upper_left_a2, name='ULA2')

    if show_spheres:
        cmd.show('spheres','box')
    cmd.color('gray40','box')

    cmd.bond("box////LRA1", "box////LLA1")
    cmd.bond("box////LLA1", "box////LLA2")
    cmd.bond("box////LLA2", "box////LRA2")
    cmd.bond("box////LRA2", "box////LRA1")

    cmd.bond("box////URA1", "box////ULA1")
    cmd.bond("box////ULA1", "box////ULA2")
    cmd.bond("box////ULA2", "box////URA2")
    cmd.bond("box////URA2", "box////URA1")

    cmd.bond("box////URA1", "box////LRA1")
    cmd.bond("box////ULA1", "box////LLA1")
    cmd.bond("box////URA2", "box////LRA2")
    cmd.bond("box////ULA2", "box////LLA2")




def build_plane(plane:str, corners:List):
    plane_colours = {
        'x':'cyan',
        'y':'magenta',
        'z':'orange'
    }

    plane_name = f"{plane}_plane"

    cmd.pseudoatom(plane_name, pos=corners[0], name='1')
    cmd.pseudoatom(plane_name, pos=corners[1], name='2')
    cmd.pseudoatom(plane_name, pos=corners[2], name='3')
    cmd.pseudoatom(plane_name, pos=corners[3], name='4')

    cmd.bond(f"{plane_name}////1", f"{plane_name}////2")
    cmd.bond(f"{plane_name}////1", f"{plane_name}////3")
    cmd.bond(f"{plane_name}////2", f"{plane_name}////4")
    cmd.bond(f"{plane_name}////3", f"{plane_name}////4")
    cmd.color(plane_colours[plane], plane_name)



def display_reference_planes(com:List):
    xyz = [27, 30, 25]

    y_offset = 5

    x1 = [(com[0]) - xyz[0], (com[1] - y_offset), com[2] - xyz[2] ]
    x2 = [(com[0]) + xyz[0], (com[1] - y_offset), com[2] - xyz[2] ]
    x3 = [(com[0]) - xyz[0], (com[1] - y_offset), com[2] + xyz[2] ]
    x4 = [(com[0]) + xyz[0], (com[1] - y_offset), com[2] + xyz[2] ]
    
    build_plane('x', [x1,x2,x3,x4])

    y1 = [(com[0]) - xyz[0], (com[1] - y_offset), com[2]]
    y2 = [(com[0]) + xyz[0], (com[1] - y_offset), com[2]]
    y3 = [(com[0]) - xyz[0], (com[1] - y_offset + xyz[1]), com[2]]
    y4 = [(com[0]) + xyz[0], (com[1] - y_offset + xyz[1]), com[2]]
    
    build_plane('y', [y1,y2,y3,y4])

    z1 = [(com[0]), (com[1] - y_offset), (com[2] - xyz[2])]
    z2 = [(com[0]), (com[1] - y_offset), (com[2] + xyz[2])]
    z3 = [(com[0]), (com[1] - y_offset + xyz[1]), (com[2]- xyz[2])]
    z4 = [(com[0]), (com[1] - y_offset + xyz[1]), (com[2] + xyz[2])]

    build_plane('z', [z1,z2,z3,z4])

    v1 = [(com[0]), (com[1] - y_offset), com[2] ]
    v2 = [(com[0]), (com[1] - y_offset + xyz[1]), com[2] ]

    cmd.pseudoatom('vertical', pos=v1, name='1')
    cmd.pseudoatom('vertical', pos=v2, name='2')
    cmd.bond("vertical////1", "vertical////2")
    cmd.color('gray30','vertical')
    pass


def cubes(selection='all', name='', state=0, scale=0.5, color='green', opacity=0, _func=cgo_cube):
    """
    This function is a modified version of the cubes function from the cubes.py extension by Thomas Holder

    It takes a selection and creates a cube for each atom in the selection

    Args:
        selection (str): The selection of atoms to create the cubes for
        name (str): The name of the object to create
        state (int): The state of the object
        scale (float): The scale of the cube
        color (str): The color of the cube
        opacity (float): The opacity of the cube
        _func (cgo_cube): The function to create the cube

    Returns:
        None
    """
    if not name:
        name = 'voxel_cube'
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


def rainbow_color_steps(steps:int, end=6/8):
    raw_colors = [hls_to_rgb(end * i/(steps-1), 0.5, 1) for i in range(steps)]
    return [[int(175 * x) for x in color] for color in raw_colors]

### Main function below ###

voxel_map_hash = 'e91d9bdc62da8457549cfbeed4c2b0aa'

voxel_grid = json.load(open(f"output/voxel_sets/{voxel_map_hash}/voxel_set.json", 'r'))

cmd.set_view("""\
     0.712135971,    0.638299167,   -0.292294502,\
     0.037353292,    0.381308764,    0.923692703,\
     0.701047778,   -0.668712735,    0.247700989,\
     0.000000000,    0.000000000, -231.989013672,\
   -42.237426758,   56.492870331,   63.704917908,\
   195.578369141,  268.399658203,  -20.000000000 """)

# Load the MHC and peptide
load_mhc_abd(constants.canonical_pdb_code)


# Display the centre of mass for the canonical class I molecule and the centre of the voxel grid
display_pseudoatom('centre_of_box', voxel_grid['params']['centre_of_mass'], color='white', opacity=0)

print ('image of MHC and centre of mass')
cmd.png(f"output/images/mhc_and_com.png", width=3000, height=3000, ray=1)


load_peptide(constants.canonical_pdb_code)

cmd.reset()

cmd.set_view("""\
     0.712135971,    0.638299167,   -0.292294502,\
     0.037353292,    0.381308764,    0.923692703,\
     0.701047778,   -0.668712735,    0.247700989,\
     0.000000000,    0.000000000, -231.989013672,\
   -42.237426758,   56.492870331,   63.704917908,\
   195.578369141,  268.399658203,  -20.000000000 """)

cmd.png(f"output/images/mhc_com_and_peptide.png", width=3000, height=3000, ray=1)


display_reference_planes(voxel_grid['params']['centre_of_mass'])

cmd.png(f"output/images/mhc_com_and_peptide_and_planes.png", width=3000, height=3000, ray=1)

display_cleft_bounding_box(voxel_grid['params']['centre_of_mass'], voxel_grid['params']['box_xyz'], y_offset=voxel_grid['params']['y_offset'], show_spheres=False)

cmd.png(f"output/images/mhc_com_and_peptide_and_planes_and_box.png", width=3000, height=3000, ray=1)

# Load the file containing the position of the used voxels and their frequency
position_voxels = json.load(open(f"output/voxel_sets/{voxel_map_hash}/position_voxels.json", 'r'))
voxel_grid = json.load(open(f"output/voxel_sets/{voxel_map_hash}/voxel_set.json", 'r'))

# we'll now make the cubes method available as an extension in Pymol
cmd.extend('cubes', cubes)



# first we'll set up the colours for the occupancies
colour_steps = rainbow_color_steps(101)

for i, colour in enumerate(colour_steps):
    cmd.set_color(f'occupancy_{i}', colour)

# next we'll create arrays for the voxes, the raw occupancy counts and the voxel centres
voxels = []
occupancy_counts = []
centres = []

# iterate through the position voxels dictionary
for position in position_voxels:
    # we'll create a count of the structures in the dictionary for use later in creating percentages
    structure_count = 0
    # for each voxel at a specific position
    for voxel in position_voxels[position]:
        # append the voxe address to the voxels array
        voxels.append(voxel)
        # add the occupancy count to the occupancy_counts array
        occupancy_counts.append(position_voxels[position][voxel]['count'])
        # increment the structure count by the count for that voxel
        structure_count += position_voxels[position][voxel]['count']
        # append the voxel centre to the array
        centres.append(voxel_grid['voxels'][voxel]['centre'])

# now we'll find the max and min occupancy for normalisation
max_occupancy = max(occupancy_counts)
min_occupancy = min(occupancy_counts)

# for each voxel address
for i, voxel in enumerate(voxels):
        # get a normalised count
        normalised_count = (occupancy_counts[i] - min_occupancy) / (max_occupancy - min_occupancy)
        
        occupancy = normalised_count
        percentage_occupancy = round(occupancy * 100)
        opacity = (1.0 - round(occupancy, 1)) - 0.4
        voxel_centre = centres[i]
        
        display_pseudoatom(voxel, voxel_centre, color='white', opacity=opacity)
        cubes(selection=voxel, name=f'voxel_{voxel}', color=f'occupancy_{percentage_occupancy}', opacity=opacity)


cmd.reset()

cmd.set_view("""\
     0.712135971,    0.638299167,   -0.292294502,\
     0.037353292,    0.381308764,    0.923692703,\
     0.701047778,   -0.668712735,    0.247700989,\
     0.000000000,    0.000000000, -231.989013672,\
   -42.237426758,   56.492870331,   63.704917908,\
   195.578369141,  268.399658203,  -20.000000000 """)

cmd.png(f"output/images/all.png", width=3000, height=3000, ray=1)


cmd.hide('everything', 'centre_of_box')
cmd.hide('everything', 'x_plane')
cmd.hide('everything', 'y_plane')
cmd.hide('everything', 'z_plane')
cmd.hide('everything', 'vertical')

cmd.png(f"output/images/no_planes_no_com.png", width=3000, height=3000, ray=1)

cmd.hide('everything', 'antigen_binding_domain')

cmd.png(f"output/images/peptide_and_bounding_box.png", width=3000, height=3000, ray=1)
cmd.hide('everything', 'peptide')

cmd.png(f"output/images/voxels_and_bounding_box.png", width=3000, height=3000, ray=1)







