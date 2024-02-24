from typing import List, Dict



def create_voxel_grid(centre_of_mass:List[float], box_xyz:List[float], voxel_size:float, range_offset:int=0, x_offset:int=0, y_offset:int=0, z_offset:int=0):
    

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
    is_inside = (x_corner <= x <= x_max) and (y_corner <= y <= y_max) and (z_corner <= z <= z_max)

    return is_inside



  