# Author: Mian Qin
# Date Created: 2024/12/2
import subprocess
import shlex
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from MDAnalysis.transformations import rotateby


def generate_ice():
    lattice_target = np.array([8, 8, 8])
    lattice_unit = np.array([0.78228388, 0.90365185, 0.73535726])
    rep = np.ceil(lattice_target / lattice_unit).astype(int)
    rep = f"{rep[0]} {rep[2]} {rep[1]}"
    command = f"genice2 ice1h --water tip4p --format gromacs --rep {rep}"

    temp_file_name = ".temp.gro"
    with open(temp_file_name, "w") as file:
        subprocess.run(shlex.split(command), stdout=file)

    u = mda.Universe(temp_file_name)
    rotation = rotateby(270, [1, 0, 0], ag=u.atoms)
    u.trajectory.add_transformations(rotation)

    with mda.Writer(temp_file_name, u.atoms.n_atoms) as writer:
        for ts in u.trajectory:
            writer.write(u)

    dimensions = u.dimensions
    command = f"gmx editconf -f {temp_file_name} -o box_of_ice.gro -c -box {dimensions[0] / 10} {dimensions[2] / 10} {dimensions[1] / 10}"
    print(shlex.split(command))
    subprocess.run(shlex.split(command))
    Path(temp_file_name).unlink()


def main():
    generate_ice()


if __name__ == "__main__":
    main()
