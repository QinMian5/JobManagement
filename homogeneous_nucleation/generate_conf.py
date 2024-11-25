import subprocess
import shlex
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from MDAnalysis.transformations import rotateby


def _1_generate_ice():
    lattice_target = np.array([6, 6, 4])
    lattice_unit = np.array([0.78228388, 0.90365185, 0.73535726])
    rep = np.ceil(lattice_target / lattice_unit).astype(int)
    rep = f"{rep[0]} {rep[2]} {rep[1]}"
    command = f"genice2 ice1h --water tip4p --format gromacs --rep {rep}"
    with open("temp.gro", "w") as file:
        subprocess.run(shlex.split(command), stdout=file)


def _2_rotate_ice():
    u = mda.Universe("temp.gro")
    rotation = rotateby(270, [1, 0, 0], ag=u.atoms)
    u.trajectory.add_transformations(rotation)

    with mda.Writer("temp.gro", u.atoms.n_atoms) as writer:
        for ts in u.trajectory:
            writer.write(u)

    dimensions = u.dimensions
    command = f"gmx editconf -f temp.gro -o conf.gro -c -box {dimensions[0]/10} {dimensions[2]/10} {dimensions[1]/10}"
    print(shlex.split(command))
    subprocess.run(shlex.split(command))
    Path("temp.gro").unlink()


def _3_select_pseudo_ice():
    u = mda.Universe("conf.gro")

    x_size, y_size, z_size = u.dimensions[:3]

    all_o_atoms = u.select_atoms("name OW")
    SOL_o_atoms = u.atoms[[]]
    PI_o_atoms = u.atoms[[]]
    for atom in all_o_atoms:
        x_pos, y_pos, z_pos = atom.position
        assert x_pos >= 0 and y_pos >= 0 and z_pos >= 0, "All coordinates must be positive."
        # Base
        if z_pos <= 10:
            PI_o_atoms += atom
            continue
        # Others are solvant
        SOL_o_atoms += atom
    SOL_resids = " ".join(map(str, SOL_o_atoms.resids))
    SOL_water_molecules = u.select_atoms(f"resid {SOL_resids}")
    PI_resids = " ".join(map(str, PI_o_atoms.resids))
    PI_water_molecules = u.select_atoms(f"resid {PI_resids}")
    for atom in PI_water_molecules:
        atom.name = atom.name.replace("W", "_PI")
    print(f"SOL     {len(SOL_o_atoms)}")
    print(f"PI      {len(PI_o_atoms)}")

    # Write .gro file
    with open("conf.gro", "w") as f:
        f.write("Pseudo Ice\n")
        f.write(f"{len(SOL_water_molecules) + len(PI_water_molecules):5d}\n")

        for i, atom in enumerate(SOL_water_molecules, start=0):
            line = f"{i // 4 + 1:5d}{'SOL':<5s}{atom.name:>5s}{i + 1:5d}"
            line += f"{atom.position[0] / 10:8.3f}{atom.position[1] / 10:8.3f}{atom.position[2] / 10:8.3f}\n"
            f.write(line)
        for i, atom in enumerate(PI_water_molecules, start=len(SOL_water_molecules)):
            line = f"{i // 4 + 1:5d}{'PI':<5s}{atom.name:>5s}{i + 1:5d}"
            line += f"{atom.position[0] / 10:8.3f}{atom.position[1] / 10:8.3f}{atom.position[2] / 10:8.3f}\n"
            f.write(line)

        f.write(f"   {u.dimensions[0] / 10:10.5f}{u.dimensions[1] / 10:10.5f}{u.dimensions[2] / 10:10.5f}\n")


def _4_create_box_of_water():
    h_z = 40.0
    u_ice = mda.Universe("conf.gro")
    u_water = mda.Universe("box_of_water.gro")

    x_size, y_size, z_size = u_ice.dimensions[:3]
    water_to_add_o = u_water.select_atoms(
        f"name OW"
        f" and (prop x >= 0.0) and (prop x <= {x_size})"
        f" and (prop y >= 0.0) and (prop y <= {y_size})"
        f" and (prop z >= 0.0) and (prop z <= {h_z})")
    water_to_add_resids_o = {str(x) for x in water_to_add_o.resids}
    water_to_add_h1 = u_water.select_atoms(
        f"name HW1"
        f" and (prop x >= 0.0) and (prop x <= {x_size})"
        f" and (prop y >= 0.0) and (prop y <= {y_size})"
        f" and (prop z >= 0.0) and (prop z <= {h_z})")
    water_to_add_resids_h1 = {str(x) for x in water_to_add_h1.resids}
    PW_to_add_h2 = u_water.select_atoms(
        f"name HW2"
        f" and (prop x >= 0.0) and (prop x <= {x_size})"
        f" and (prop y >= 0.0) and (prop y <= {y_size})"
        f" and (prop z >= 0.0) and (prop z <= {h_z})")
    water_to_add_resids_h2 = {str(x) for x in PW_to_add_h2.resids}
    water_to_add_resids = water_to_add_resids_o.intersection(water_to_add_resids_h1.intersection(water_to_add_resids_h2))
    water_to_add = u_water.select_atoms(f"resid {' '.join(water_to_add_resids)}")

    n_SOL = len(water_to_add) // 4
    print(f"SOL     {n_SOL}")

    with open("small_box_of_water.gro", 'w') as f:
        f.write("Small Box of Water\n")
        f.write(f"{n_SOL * 4:5d}\n")

        i = 0
        for atom in water_to_add:
            assert atom.resname == "SOL"
            line = f"{i // 4 + 1:5d}{'SOL':<5s}{atom.name:>5s}{i + 1:5d}{atom.position[0] / 10:8.3f}{atom.position[1] / 10:8.3f}{atom.position[2] / 10:8.3f}\n"
            f.write(line)
            i += 1
        f.write(f"   {x_size / 10:10.5f}{y_size / 10:10.5f}{h_z / 10:10.5f}\n")


def _5_insert_PW():
    h_PW = 10.0
    h_water = 10.0
    u_ice = mda.Universe("conf.gro")
    u_water = mda.Universe("small_box_of_water_out.gro")

    # TODO: 选出来水分子进行EM，NVT，NPT，NPT固定x、y的大小
    x_size, y_size, z_size = u_ice.dimensions[:3]
    PI_origin = u_ice.select_atoms("resname PI")
    SOL_origin = u_ice.select_atoms("resname SOL")
    PW_to_add_o = u_water.select_atoms(
        f"name OW"
        # f" and (prop x >= 0.0) and (prop x <= {x_size})"
        # f" and (prop y >= 0.0) and (prop y <= {y_size})"
        f" and (prop z >= 0.0) and (prop z <= {h_PW})")
    PW_to_add_resids_o = {str(x) for x in PW_to_add_o.resids}
    PW_to_add_h1 = u_water.select_atoms(
        f"name HW1"
        # f" and (prop x >= 0.0) and (prop x <= {x_size})"
        # f" and (prop y >= 0.0) and (prop y <= {y_size})"
        f" and (prop z >= 0.0) and (prop z <= {h_PW})")
    PW_to_add_resids_h1 = {str(x) for x in PW_to_add_h1.resids}
    PW_to_add_h2 = u_water.select_atoms(
        f"name HW2"
        # f" and (prop x >= 0.0) and (prop x <= {x_size})"
        # f" and (prop y >= 0.0) and (prop y <= {y_size})"
        f" and (prop z >= 0.0) and (prop z <= {h_PW})")
    PW_to_add_resids_h2 = {str(x) for x in PW_to_add_h2.resids}
    PW_to_add_resids = PW_to_add_resids_o.intersection(PW_to_add_resids_h1.intersection(PW_to_add_resids_h2))
    PW_to_add = u_water.select_atoms(f"resid {' '.join(PW_to_add_resids)}")

    water_to_add_o = u_water.select_atoms(
        f"name OW"
        # f" and (prop x >= 0.0) and (prop x <= {x_size})"
        # f" and (prop y >= 0.0) and (prop y <= {y_size})"
        f" and (prop z >= {h_PW + 5.0}) and (prop z <= {h_PW + h_water + 5.0})")
    water_to_add_resids_o = {str(x) for x in water_to_add_o.resids}
    water_to_add_h1 = u_water.select_atoms(
        f"name HW1"
        # f" and (prop x >= 0.0) and (prop x <= {x_size})"
        # f" and (prop y >= 0.0) and (prop y <= {y_size})"
        f" and (prop z >= {h_PW + 5.0}) and (prop z <= {h_PW + h_water + 5.0})")
    water_to_add_resids_h1 = {str(x) for x in water_to_add_h1.resids}
    water_to_add_h2 = u_water.select_atoms(
        f"name HW2"
        # f" and (prop x >= 0.0) and (prop x <= {x_size})"
        # f" and (prop y >= 0.0) and (prop y <= {y_size})"
        f" and (prop z >= {h_PW + 5.0}) and (prop z <= {h_PW + h_water + 5.0})")
    water_to_add_resids_h2 = {str(x) for x in water_to_add_h2.resids}
    water_to_add_resids = water_to_add_resids_o.intersection(water_to_add_resids_h1.intersection(water_to_add_resids_h2))
    water_to_add = u_water.select_atoms(f"resid {' '.join(water_to_add_resids)}")

    interval = 1.0
    PI_origin.translate([0, 0, h_PW + interval])
    SOL_origin.translate([0, 0, h_PW + interval])
    water_to_add.translate([0, 0, -(h_PW + 5.0) + h_PW + z_size + 2 * interval])

    for atom in PW_to_add:
        atom.name = atom.name.replace("W", "_PW")

    n_SOL = len(SOL_origin) // 4 + len(water_to_add) // 4
    n_PI = len(PI_origin) // 4
    n_PW = len(PW_to_add) // 4

    print(f"SOL     {n_SOL}")
    print(f"PI      {n_PI}")
    print(f"PW      {n_PW}")

    # Write .gro file
    with open("confout.gro", "w") as f:
        f.write("Pseudo Ice\n")
        f.write(f"{(n_SOL + n_PI + n_PW) * 4:5d}\n")

        i = 0
        for atom in SOL_origin:
            assert atom.resname == "SOL"
            line = f"{i // 4 + 1:5d}{'SOL':<5s}{atom.name:>5s}{i + 1:5d}{atom.position[0] / 10:8.3f}{atom.position[1] / 10:8.3f}{atom.position[2] / 10:8.3f}\n"
            f.write(line)
            i += 1
        for atom in water_to_add:
            assert atom.resname == "SOL"
            line = F"{i // 4 + 1:5d}{'SOL':<5s}{atom.name:>5s}{i + 1:5d}{atom.position[0] / 10:8.3f}{atom.position[1] / 10:8.3f}{atom.position[2] / 10:8.3f}\n"
            f.write(line)
            i += 1
        for atom in PI_origin:
            assert atom.resname == "PI"
            line = f"{i // 4 + 1:5d}{'PI':<5s}{atom.name:>5s}{i + 1:5d}{atom.position[0] / 10:8.3f}{atom.position[1] / 10:8.3f}{atom.position[2] / 10:8.3f}\n"
            f.write(line)
            i += 1
        for atom in PW_to_add:
            assert atom.resname == "SOL"
            line = f"{i // 4 + 1:5d}{'PW':<5s}{atom.name:>5s}{i + 1:5d}{atom.position[0] / 10:8.3f}{atom.position[1] / 10:8.3f}{atom.position[2] / 10:8.3f}\n"
            f.write(line)
            i += 1

        f.write(f"   {x_size / 10:10.5f}{y_size / 10:10.5f}{(z_size + h_PW + h_water + 3 * interval) / 10:10.5f}\n")


def main():
    _1_generate_ice()
    _2_rotate_ice()
    _3_select_pseudo_ice()
    # _4_create_box_of_water()
    _5_insert_PW()


if __name__ == "__main__":
    main()
