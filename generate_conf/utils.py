import subprocess
import shlex
from pathlib import Path
import numpy as np
import MDAnalysis as mda
# from MDAnalysis.transformations import rotateby
import pyperclip


def run_command_by_user(command_to_run: str):
    pyperclip.copy(command_to_run)
    input("Command is copied. Run command. Press Enter after finished:")


class Region:
    def __init__(self, region_func):
        self.region_func = region_func

    def __call__(self, point):
        return self.region_func(point)

    def union(self, other):
        """Return a new Region representing the union of this region with another."""
        return Region(lambda point: self(point) or other(point))

    def intersection(self, other):
        """Return a new Region representing the intersection of this region with another."""
        return Region(lambda point: self(point) and other(point))

    def difference(self, other):
        """Return a new Region representing the difference of this region from another."""
        return Region(lambda point: self(point) and not other(point))


def define_region(shape, **params):
    """
    Returns a function that determines if a point is within the defined region.

    :param shape: str, type of shape ('cuboid', 'cylinder', 'sphere', 'custom')
    :param params: dict, parameters necessary to define the shape (like center, radius, height, etc.)
    :return: function, a function that takes a coordinate as input and returns True if inside the region
    """
    if shape == 'cuboid':
        # Cuboid defined by its center and dimensions (length, width, height)
        x_range = params["x_range"]
        y_range = params["y_range"]
        z_range = params["z_range"]
        center = np.array([np.mean(x_range), np.mean(y_range), np.mean(z_range)])
        dimensions = np.array([x_range[1] - x_range[0], y_range[1] - y_range[0], z_range[1] - z_range[0]]) / 2

        def cuboid_region(point):
            temp = np.abs(point - center) <= dimensions
            return np.all(temp)

        return Region(cuboid_region)

    elif shape == 'cylinder':
        # Cylinder defined by center, radius, height along z-axis
        z_range = params["z_range"]
        center = np.array(params['center'])
        radius = params['radius']
        z_center = np.mean(z_range)
        height = (z_range[1] - z_range[0]) / 2

        def cylinder_region(point):
            radial_distance = np.sqrt(np.sum(np.square(point[..., :2] - center), axis=-1))
            axial_distance = np.abs(point[..., 2] - z_center)
            return (radial_distance <= radius) & (axial_distance <= height)

        return Region(cylinder_region)

    elif shape == 'sphere':
        # Sphere defined by center and radius
        center = np.array(params['center'])
        radius = params['radius']

        def sphere_region(point):
            return np.sqrt(np.sum(np.square(point - center), axis=-1)) <= radius

        return Region(sphere_region)

    elif shape == 'custom':
        # Custom region defined by a user-provided function
        return Region(params['region_func'])

    else:
        raise ValueError("Unsupported region shape")


def select_atoms_by_region(atom_group: mda.AtomGroup, region: Region, keep_partial=False):
    """
    Selects atoms from an MDAnalysis AtomGroup based on a region function.
    Can optionally keep or discard entire molecules based on partial inclusion within the region.

    :param atom_group: MDAnalysis AtomGroup, the group of atoms to filter
    :param region: Region, a function that takes a coordinate as input and returns True if inside the region
    :param keep_partial: bool, if True, keeps entire molecules if any part is within the region;
                               if False, only keeps molecules fully within the region
    :return: MDAnalysis AtomGroup, the filtered group of atoms
    """
    # Identify all molecules in the AtomGroup
    unique_resids = set(atom.resid for atom in atom_group)

    # Container for selected residues
    selected_resids = set()

    # Iterate over each molecule by its residue ID
    for resid in unique_resids:
        # Select all atoms in the current residue
        residue_atoms = atom_group.select_atoms(f"resid {resid}")

        # Check if any or all atoms are within the defined region based on keep_partial
        if keep_partial:
            # Include the molecule if any atom is within the region
            if np.any([region(atom.position) for atom in residue_atoms]):
                selected_resids.add(resid)
        else:
            # Include the molecule only if all atoms are within the region
            if np.all([region(atom.position) for atom in residue_atoms]):
                selected_resids.add(resid)

    # Build the final AtomGroup from the selected residues
    if selected_resids:
        selected_atoms = atom_group.select_atoms(f"resid {' '.join([str(x) for x in selected_resids])}")
    else:
        selected_atoms = atom_group.select_atoms("none")  # Empty selection

    return selected_atoms


def modify_resname(atom_group: mda.AtomGroup, resname_map: dict[str, str]):
    for resname, new_resname in resname_map.items():
        atoms_to_change = atom_group.select_atoms(f"resname {resname}")
        for atom in atoms_to_change:
            atom.residue.resname = new_resname


def modify_name(atom_group: mda.AtomGroup, name_map: dict[str, str]):
    for atomname, new_atomname in name_map.items():
        atoms_to_change = atom_group.select_atoms(f"name {atomname}")
        for atom in atoms_to_change:
            atom.name = new_atomname


def merge_atom_groups(atom_groups: list[mda.AtomGroup], dimensions=np.array([0, 0, 0, 90, 90, 90])):
    """
    Merges multiple MDAnalysis AtomGroups from different universes into a single new universe.

    :param atom_groups: List of MDAnalysis AtomGroups to be merged.
    :param dimensions: Dimensions of the box
    :return: MDAnalysis Universe containing all atoms from the provided AtomGroups.
    """
    # Initialize lists to collect atom data
    atom_names = []
    positions = []
    resnames = []

    # Mappings
    atom_resindex = []

    # Counters for unique IDs
    current_residue_index = 0

    # Dictionaries to track unique residues and segments
    residue_dict = {}
    segment_dict = {}

    for i, ag in enumerate(atom_groups):
        for atom in ag:
            # Collect atom attributes
            atom_names.append(atom.name)
            positions.append(atom.position)

            # Handle residues
            resid = (i, atom.resid)
            if resid not in residue_dict:
                resnames.append(atom.resname)
                residue_dict[resid] = current_residue_index
                current_residue_index += 1
            atom_resindex.append(residue_dict[resid])

    # Create residue_segindex mapping
    # residue_segindex = [segment_dict[segid] for segid, _ in residue_dict.keys()]

    resids = list(residue_dict.values())

    n_atoms = len(atom_names)
    n_residues = current_residue_index

    # Create the new Universe
    merged_universe = mda.Universe.empty(
        n_atoms,
        n_residues=n_residues,
        n_segments=1,
        atom_resindex=atom_resindex,
        # residue_segindex=1,
        trajectory=True,
    )

    # Add topology attributes
    merged_universe.add_TopologyAttr('name', atom_names)
    merged_universe.add_TopologyAttr('resname', resnames)
    merged_universe.add_TopologyAttr('resid', resids)

    # Load positions
    merged_universe.atoms.positions = np.array(positions)
    merged_universe.dimensions = dimensions
    return merged_universe


def main():
    ...


if __name__ == "__main__":
    main()
