# Author: Mian Qin
# Date Created: 12/2/24
from utils import *


def main():
    h_PI = 10
    h_PW = 10
    interval_between_PW_PI = 1
    h_SOL_ice = 60
    h_SOL_water = 10
    margin = 10

    u_box_of_ice = mda.Universe("box_of_ice.gro")
    u_large_box_of_water = mda.Universe("large_box_of_water.gro")
    u_small_box_of_water = mda.Universe("small_box_of_water.gro")

    x_max, y_max = u_box_of_ice.dimensions[:2]

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
    ice_translation = [0, 0, 10 + 1]  # 1 Angstrom spacing between PW and PI
    PI_ag.translate(ice_translation)
    SOL_ice_ag.translate(ice_translation)
    SOL_water_ag.translate(ice_translation)
    print("Renaming atoms")
    modify_resname(PW_ag, {"SOL": "PW"})
    modify_name(PW_ag, {"OW": "O_PW", "HW1": "H_PW1", "HW2": "H_PW2", "MW": "M_PW"})
    modify_resname(PI_ag, {"ICE": "PI"})
    modify_name(PI_ag, {"OW": "O_PI", "HW1": "H_PI1", "HW2": "H_PI2", "MW": "M_PI"})
    merge_atom_groups([SOL_ice_ag, SOL_water_ag, PI_ag, PW_ag])


if __name__ == "__main__":
    main()
