# Author: Mian Qin
# Date Created: 2024/12/2
import subprocess
import shlex
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from MDAnalysis.transformations import rotateby

from utils import *


file_box_of_ice = Path("pseudoice/box_of_ice.gro")


def generate_ice():
    lattice_target = np.array([7, 7, 7])
    lattice_unit = np.array([0.78228388, 0.90365185, 0.73535726])
    rep = np.ceil(lattice_target / lattice_unit).astype(int)
    rep = f"{rep[0]} {rep[2]} {rep[1]}"
    temp_file_name = "temp.gro"
    command_to_run = f"genice2 ice1h --water tip4p --format gromacs --rep {rep} > {temp_file_name}"
    run_command_by_user(command_to_run)

    u = mda.Universe(temp_file_name)
    rotation = rotateby(270, [1, 0, 0], ag=u.atoms)
    u.trajectory.add_transformations(rotation)

    with mda.Writer(temp_file_name, u.atoms.n_atoms) as writer:
        for ts in u.trajectory:
            writer.write(u)

    dimensions = u.dimensions
    command_to_run = f"gmx editconf -f {temp_file_name} -o {file_box_of_ice} -c -box {dimensions[0] / 10} {dimensions[2] / 10} {dimensions[1] / 10}"
    run_command_by_user(command_to_run)
    Path(temp_file_name).unlink()


def generate_water():
    u_box_of_ice = mda.Universe(file_box_of_ice)
    u_water = mda.Universe("box_of_water.gro")

    x_max, y_max = u_box_of_ice.dimensions[:2]
    z_max = 30

    small_box_region = define_region("cuboid", x_range=[0, x_max], y_range=[0, y_max], z_range=[0, z_max])
    small_box_ag = select_atoms_by_region(u_water.select_atoms("all"), small_box_region)

    merged_universe = merge_atom_groups([small_box_ag], dimensions=np.array([x_max, y_max, z_max, 90, 90, 90]))
    all_atoms = merged_universe.select_atoms("all")

    with mda.Writer("pseudoice/small_box_of_water.gro", all_atoms.n_atoms) as w:
        for ts in merged_universe.trajectory:
            w.write(all_atoms)


def main():
    # generate_ice()
    generate_water()


if __name__ == "__main__":
    main()
