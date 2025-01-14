# Author: Mian Qin
# Date Created: 12/2/24
import MDAnalysis as mda
from utils import *


def generate_hom():
    margin = 7

    u_box_of_ice = mda.Universe("box_of_ice.gro")
    u_large_box_of_water = mda.Universe("box_of_water.gro")

    x_max, y_max, z_max = u_box_of_ice.dimensions[:3]
    SOL_ice_region = define_region("cuboid", x_range=[margin, x_max - margin], y_range=[margin, y_max - margin],
                                   z_range=[margin, z_max - margin])
    SOL_ice_ag = select_atoms_by_region(u_box_of_ice.select_atoms("all"), SOL_ice_region)
    SOL_water_region = define_region("cuboid", x_range=[0, x_max], y_range=[0, y_max],
                                                        z_range=[0, z_max]).difference(SOL_ice_region)
    SOL_water_ag = select_atoms_by_region(u_large_box_of_water.select_atoms("all"), SOL_water_region)

    modify_resname(SOL_ice_ag, {"ICE": "SOL"})

    merged_universe = merge_atom_groups([SOL_ice_ag, SOL_water_ag],
                                        dimensions=np.array([x_max, y_max, z_max, 90, 90, 90]))
    all_atoms = merged_universe.select_atoms("all")

    with mda.Writer("../homogeneous_nucleation/conf.gro", all_atoms.n_atoms) as w:
        for ts in merged_universe.trajectory:
            w.write(all_atoms)


def generate_het():
    # Unit: Angstrom
    h_PI = 10
    h_PW = 10
    interval_between_PW_PI = 20
    h_SOL_ice = 30
    h_SOL_water = 15
    margin = 7

    u_box_of_ice = mda.Universe("box_of_ice.gro")
    u_large_box_of_water = mda.Universe("box_of_water.gro")
    u_small_box_of_water = mda.Universe("small_box_of_water.gro")

    x_max, y_max = u_box_of_ice.dimensions[:2]
    z_max = h_PW + interval_between_PW_PI + h_PI + h_SOL_ice + h_SOL_water

    PW_region = define_region("cuboid", x_range=[0, x_max], y_range=[0, y_max], z_range=[0, h_PW])
    PW_ag = select_atoms_by_region(u_small_box_of_water.select_atoms("all"), PW_region)
    PI_region = define_region("cuboid", x_range=[0, x_max], y_range=[0, y_max], z_range=[0, h_PI])
    PI_ag = select_atoms_by_region(u_box_of_ice.select_atoms("all"), PI_region)
    SOL_ice_region = define_region("cuboid", x_range=[margin, x_max - margin], y_range=[margin, y_max - margin],
                                   z_range=[h_PI, h_PI + h_SOL_ice])
    SOL_ice_ag = select_atoms_by_region(u_box_of_ice.select_atoms("all"), SOL_ice_region)
    SOL_water_region = define_region("cuboid", x_range=[0, x_max], y_range=[0, y_max],
                                     z_range=[h_PI, h_PI + h_SOL_ice + h_SOL_water]).difference(SOL_ice_region)
    SOL_water_ag = select_atoms_by_region(u_large_box_of_water.select_atoms("all"), SOL_water_region)
    ice_translation = [0, 0, 10 + interval_between_PW_PI]  # 1 Angstrom spacing between PW and PI
    PI_ag.translate(ice_translation)
    SOL_ice_ag.translate(ice_translation)
    SOL_water_ag.translate(ice_translation)
    print("Renaming atoms")
    modify_resname(SOL_ice_ag, {"ICE": "SOL"})
    modify_resname(PW_ag, {"SOL": "PW"})
    modify_name(PW_ag, {"OW": "O_PW", "HW1": "H_PW1", "HW2": "H_PW2", "MW": "M_PW"})
    modify_resname(PI_ag, {"ICE": "PI"})
    modify_name(PI_ag, {"OW": "O_PI", "HW1": "H_PI1", "HW2": "H_PI2", "MW": "M_PI"})
    print("Finished renaming atoms")
    merged_universe = merge_atom_groups([SOL_ice_ag, SOL_water_ag, PI_ag, PW_ag],
                                        dimensions=np.array([x_max, y_max, z_max, 90, 90, 90]))
    all_atoms = merged_universe.select_atoms("all")

    with mda.Writer("../heterogeneous_nucleation/conf.gro", all_atoms.n_atoms) as w:
        for ts in merged_universe.trajectory:
            w.write(all_atoms)


def main():
    generate_hom()
    # generate_het()


if __name__ == "__main__":
    main()
