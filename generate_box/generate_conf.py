import subprocess
import shlex
from pathlib import Path
import numpy as np
import MDAnalysis as mda
# from MDAnalysis.transformations import rotateby


def define_region(shape, **params):
    """
    Returns a function that determines if a point is within the defined region.

    :param shape: str, type of shape ('cuboid', 'cylinder', 'sphere', 'custom')
    :param params: dict, parameters necessary to define the shape (like center, radius, height, etc.)
    :return: function, a function that takes a coordinate as input and returns True if inside the region
    """
    if shape == 'cuboid':
        # Cuboid defined by its center and dimensions (length, width, height)
        center = np.array(params['center'])
        dimensions = np.array(params['dimensions']) / 2  # half-dimensions for easier calculations

        def cuboid_region(point):
            return np.all(np.abs(point - center) <= dimensions, axis=-1)

        return cuboid_region

    elif shape == 'cylinder':
        # Cylinder defined by center, radius, height along z-axis
        center = np.array(params['center'])
        radius = params['radius']
        height = params['height'] / 2

        def cylinder_region(point):
            radial_distance = np.sqrt((point[0] - center[0]) ** 2 + (point[1] - center[1]) ** 2)
            axial_distance = abs(point[2] - center[2])
            return radial_distance <= radius and axial_distance <= height

        return cylinder_region

    elif shape == 'sphere':
        # Sphere defined by center and radius
        center = np.array(params['center'])
        radius = params['radius']

        def sphere_region(point):
            return np.sqrt(np.sum((point - center) ** 2)) <= radius

        return sphere_region

    elif shape == 'custom':
        # Custom region defined by a user-provided function
        return params['region_func']

    else:
        raise ValueError("Unsupported region shape")


def main():
    cuboid = define_region('cuboid', center=[5, 5, 5], dimensions=[10, 10, 10])
    print(cuboid)


if __name__ == "__main__":
    main()
