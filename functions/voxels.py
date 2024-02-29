
from typing import List, Dict

def create_voxel_grid(centre_of_mass:List[float], box_xyz:List[float], voxel_size:float, range_offset:int=0, x_offset:int=0, y_offset:int=0, z_offset:int=0) -> Dict:
    """
    This function creates a voxel grid based on the centre of mass and the dimensions of the box, the size of voxels, and the range and x, y, and z offsets.

    Args:
        centre_of_mass (List[float]): The centre of mass of the structure.
        box_xyz (List[float]): The dimensions of the box.
        voxel_size (float): The size of the voxels.
        range_offset (int): The range offset (this is currently only used for the current PyMol display of the grid).
        x_offset (int): The x offset.
        y_offset (int): The y offset.
        z_offset (int): The z offset.

    Returns:
        Dict: A dictionary containing the voxel grid.
    """

    params = dict(zip(locals().keys(), locals().values()))

    voxel_grid = {
        "params": params,
        "voxels": {}, 
        "labels": [],
        "metadata": {}
    }

    cx, cy, cz = centre_of_mass
    length, width, height = box_xyz

    start_x = cx - (length / 2) + x_offset
    start_y = cy - (width / 2) + y_offset
    start_z = cz - (height / 2) + z_offset  

    end_x = start_x + length
    end_y = start_y + width
    end_z = start_z + height

    num_cubes_x = int(length / voxel_size)
    num_cubes_y = int(width / voxel_size)
    num_cubes_z = int(height / voxel_size)


    for z in range(num_cubes_z + range_offset):
        for y in range(num_cubes_y + range_offset):
            for x in range(num_cubes_x + range_offset):
                voxel_label = f"{x}_{y}_{z}"
                voxel_grid['voxels'][voxel_label] = {
                    'start': [start_x + x * voxel_size, start_y + y * voxel_size, start_z + z * voxel_size],
                    'end': [start_x + (x + 1) * voxel_size, start_y + (y + 1) * voxel_size, start_z + (z + 1) * voxel_size],
                    'centre': [start_x + (x + 0.5) * voxel_size, start_y + (y + 0.5) * voxel_size, start_z + (z + 0.5) * voxel_size]
                }
                voxel_grid['labels'].append(voxel_label)
    
    voxel_grid['params']['voxel_count'] = len(voxel_grid['voxels'])
    
    return voxel_grid


def find_voxels_for_structure(peptide_df, voxels, voxel_size:int=1):

    # initialise the lists and dictionaries
    coordinates_to_check = []
    sequence = []
    positions = []
    structure_voxels = {}


    # iterate through the rows of the dataframe to generate a list of coordinates to check
    for index, row in peptide_df['ATOM'].iterrows():

        # first, extract the relevant data from the row
        atom_name = row['atom_name']
        residue_name = row['residue_name']
        chain_id = row['chain_id']
        x, y, z = row['x_coord'], row['y_coord'], row['z_coord']
        residue_number = row['residue_number']
        
        # if the atom name is the target atom name, add the residue name to the sequence list and the coordinates to the coordinates_to_check list
        if atom_name == 'CA': # this is hard-coded for now, but should be a parameter in the future 
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
            if is_coordinate_inside_voxel(coordinate[0],coordinate[1],coordinate[2], voxels[voxel]['start'], voxel_size):
                structure_voxels[str(p)] = {
                    'position': p,
                    'voxel_label': voxel,
                    'atom_name': 'CA',
                    'atom_coordinates': coordinate,
                    'voxel_start': voxels[voxel]['start'],
                    'voxel_centre': voxels[voxel]['centre'],
                    'residue': sequence[p-1]
                }
                break
        p += 1
    # return the structure_voxels dictionary
    return structure_voxels


def is_coordinate_inside_voxel(x:float, y:float, z:float, voxel_corner, voxel_size) -> bool:
    """
    Detects whether a three-dimensional coordinate is contained within a specific cube.

    Parameters:
        x (float): The x-coordinate of the target point.
        y (float): The y-coordinate of the target point.
        n (float): The size of the cube (cube's side length).
        voxel_corner (tuple): The three-dimensional coordinates (x, y, z) of the cube's corner.
        voxel_size (int): The size of the voxel

    Returns:
        bool: True if the coordinate is inside the cube, False otherwise.

    To Do:
        Can this be made easier by using the delta between the coordinates and the centre of mass?
    """

    # Extract the individual components of the cube's corner position
    x_corner, y_corner, z_corner = voxel_corner

    # Calculate the maximum and minimum values for each dimension of the cube
    x_max = x_corner + float(voxel_size)
    y_max = y_corner + float(voxel_size)
    z_max = z_corner + float(voxel_size)

    # Check if the target coordinate falls within the cube's dimensions
    is_inside = (x_corner <= x < x_max) and (y_corner <= y < y_max) and (z_corner <= z < z_max)

    return is_inside
