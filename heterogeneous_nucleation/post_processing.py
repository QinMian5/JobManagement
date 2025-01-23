# Author: Mian Qin
# Date Created: 9/24/24
from pathlib import Path
import argparse
import json
import pickle
from collections import OrderedDict
from itertools import product

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
import MDAnalysis as mda
from skimage import measure
import freud


def read_index_dict(file_path: Path) -> OrderedDict[str, list[str]]:
    index_dict = OrderedDict()
    with open(file_path) as file:
        for line in file:
            line = line.strip().split()
            if len(line) == 0:  # End of file
                break
            t = float(line[0])
            indices = [str(x) for x in line[1:]]
            index_dict[f"{t:.1f}"] = indices
    return index_dict


def get_water_indices(indices: list[str]):
    indices = indices.copy()
    indices.append(f"{3946 * 4 + 1}")  # num_SOL * 4 + 1
    water_indices = []
    for i in range(len(indices) - 1):
        first = int(indices[i])
        second = int(indices[i + 1])
        for index in range(first + 4, second, 4):
            water_indices.append(str(index))
    return water_indices


def write_index_dict(index_dict: OrderedDict[str, list[str]], file_path: Path):
    with open(file_path, 'w') as file:
        for t, indices in index_dict.items():
            file.write(f"{t} {' '.join(indices)}\n")
    print(f"Write indices to {file_path}")


def combine_indices(indices_list: list[OrderedDict], allow_duplicate=False):
    new_index_dict = OrderedDict()
    if indices_list:
        for t in indices_list[0].keys():
            combined_indices = []
            for indices_dict in indices_list:
                indices = indices_dict[t]
                combined_indices.extend(indices)
            if not allow_duplicate and len(set(combined_indices)) != len(combined_indices):
                raise RuntimeError("Duplicate indices")
            combined_indices.sort(key=lambda x: int(x))
            new_index_dict[t] = combined_indices
    return new_index_dict


def filter_solid_like_atoms(solid_like_atoms_dict: OrderedDict):
    for k, v, in solid_like_atoms_dict.items():
        for i in range(len(v)):
            if int(v[i]) > 3946 * 4:
                solid_like_atoms_dict[k] = v[:i]
                break
    return solid_like_atoms_dict


def calc_scale_offset(x_range: tuple[float, float], n_x: int) \
        -> tuple[float, float]:
    start, end = x_range
    scale = (end - start) / (n_x - 1)
    offset = start
    return scale, offset


def generate_grid(x_range: tuple[float, float], y_range: tuple[float, float], z_range: tuple[float, float],
                  step: float):
    n_x = int(np.ceil((x_range[1] - x_range[0]) / step))
    n_y = int(np.ceil((y_range[1] - y_range[0]) / step))
    n_z = int(np.ceil((z_range[1] - z_range[0]) / step))
    L_x = n_x * step
    L_y = n_y * step
    L_z = n_z * step
    new_x_range = x_range[0] - (L_x - (x_range[1] - x_range[0])) / 2, x_range[1] + (L_x - (x_range[1] - x_range[0])) / 2
    new_y_range = y_range[0] - (L_y - (y_range[1] - y_range[0])) / 2, y_range[1] + (L_y - (y_range[1] - y_range[0])) / 2
    new_z_range = z_range[0] - (L_z - (z_range[1] - z_range[0])) / 2, z_range[1] + (L_z - (z_range[1] - z_range[0])) / 2
    x_scale, x_offset = calc_scale_offset(new_x_range, n_x)
    y_scale, y_offset = calc_scale_offset(new_y_range, n_y)
    z_scale, z_offset = calc_scale_offset(new_z_range, n_z)
    scale = np.array([x_scale, y_scale, z_scale]).reshape(1, 3)
    offset = np.array([x_offset, y_offset, z_offset]).reshape(1, 3)
    x, y, z = np.linspace(*new_x_range, n_x), np.linspace(*new_y_range, n_y), np.linspace(*new_z_range, n_z)
    X, Y, Z = np.meshgrid(x, y, z)
    pos_grid = np.stack((X, Y, Z), axis=3)
    return pos_grid, scale, offset


def wrap_pos_vec_array(pos_vec_array: np.ndarray, box_size: np.ndarray, mode: str = "closest"):
    if mode == "closest":
        return pos_vec_array - np.round(pos_vec_array / box_size) * box_size
    elif mode == "positive":
        return pos_vec_array - np.floor(pos_vec_array / box_size) * box_size
    elif mode == "negative":
        return pos_vec_array + np.ceil(pos_vec_array / box_size) * box_size
    else:
        raise ValueError("Invalid mode. Choose from 'closest', 'positive', or 'negative'.")


def calc_concentration_at_point(pos_atoms: np.ndarray, pos_points: np.ndarray, box_size: np.ndarray,
                                ksi: float) -> np.ndarray:
    """

    :param pos_atoms: Positions of O atoms. Shape: (N, 3)
    :param pos_points: Positions of points to calculate. Shape: (N, 3)
    :param box_size: Size of the box. Shape: (3)
    :param ksi: Variance of Gaussian
    :return:
    """
    box = freud.box.Box.from_box(box_size)
    aq = freud.locality.AABBQuery(box, pos_atoms)
    r_max = 3 * ksi
    nlist = aq.query(pos_points, {"r_max": r_max}).toNeighborList()
    distances = nlist.distances
    concentrations = 1 / ((np.sqrt(2 * np.pi) * ksi) ** 3) * np.exp(-distances ** 2 / (2 * ksi ** 2))
    total_concentrations = np.bincount(nlist.query_point_indices, weights=concentrations, minlength=len(pos_points))
    return total_concentrations


def calc_concentration_on_grid(pos_atom: np.ndarray, pos_grid: np.ndarray,
                               box_size: np.ndarray, ksi: float) -> np.ndarray:
    """

    :param pos_atom: Positions of O atoms. Shape: (N, 3)
    :param pos_grid: Positions of grids. Shape: (X, Y, Z, 3)
    :param box_size: Size of the box. Shape: (3)
    :param ksi: Variance of Gaussian
    :return: Shape: (X, Y ,Z)
    """
    shape = pos_grid.shape[:3]
    concentrations = calc_concentration_at_point(pos_atom, pos_grid.reshape(-1, 3), box_size, ksi).reshape(shape)
    return concentrations


def generate_interface_from_concentration(concentration: np.ndarray, level=0.015) \
        -> tuple[np.ndarray, np.ndarray]:
    if level > concentration.max() or level < concentration.min():
        nodes, faces = np.zeros((0, 3)), np.zeros((0, 3), dtype=int)
    else:
        nodes, faces, normals, values = measure.marching_cubes(concentration, level=level)
    return nodes, faces


def calc_interface_in_t_range(u: mda.Universe, solid_like_atoms_dict: dict[str, list[str]], pos_grid, scale,
                              offset, t_range: tuple[float, float], ksi=3.5 / 2):
    instantaneous_interface_dict = {}
    acc_concentration = {"ice": np.zeros(pos_grid.shape[:3]),
                         "water": np.zeros(pos_grid.shape[:3]),
                         "surface": np.zeros(pos_grid.shape[:3])}
    acc_n = 0
    for ts in u.trajectory[::10]:
        t = ts.time
        box_size = u.dimensions[:3]
        ice_atoms_indices = solid_like_atoms_dict[f"{t:.1f}"]
        ice_atoms = u.select_atoms(f"bynum {' '.join(ice_atoms_indices)}") if ice_atoms_indices else u.select_atoms("")
        water_indices = get_water_indices(ice_atoms_indices)
        water_atoms = u.select_atoms(f"bynum {' '.join(water_indices)}") if water_indices else u.select_atoms("")
        surface_atoms = u.select_atoms("resname PI and name O_PI")
        atoms_dict = {"ice": ice_atoms, "water": water_atoms, "surface": surface_atoms}
        pos_dict = {k: v.positions for k, v in atoms_dict.items()}
        concentration_dict = {k: calc_concentration_on_grid(v, pos_grid, box_size, ksi) for k, v in pos_dict.items()}
        nodes, faces = generate_interface_from_concentration(concentration_dict["ice"])
        nodes = nodes * scale + offset

        # identify ice-water or ice-surface interface
        triangles = nodes[faces]  # [triangle_index, node_index, pos]
        triangles_center = np.mean(triangles, axis=1)
        shape = concentration_dict["ice"].shape
        x = np.arange(shape[0]) * scale[0, 0] + offset[0, 0]
        y = np.arange(shape[1]) * scale[0, 1] + offset[0, 1]
        z = np.arange(shape[2]) * scale[0, 2] + offset[0, 2]
        mean_concentration_water = RegularGridInterpolator((x, y, z), concentration_dict["water"])(triangles_center)
        mean_concentration_surface = RegularGridInterpolator((x, y, z), concentration_dict["surface"])(triangles_center)
        # print(mean_concentration_water.min(), mean_concentration_surface.min())
        interface_type = np.where(mean_concentration_water > mean_concentration_surface, 0,
                                  1)  # 0: ice-water, 1: ice-surface

        instantaneous_interface_dict[f"{t:.1f}"] = [nodes, faces, interface_type]
        if t_range[0] <= t <= t_range[1]:
            for k, v in concentration_dict.items():
                acc_concentration[k] += v
            acc_n += 1
    mean_concentration_dict = {k: v / acc_n for k, v in acc_concentration.items()}
    nodes, faces = generate_interface_from_concentration(mean_concentration_dict["ice"])
    nodes = nodes * scale + offset

    # identify ice-water or ice-surface interface
    triangles = nodes[faces]
    triangles_center = np.mean(triangles, axis=1)
    x = np.arange(mean_concentration_dict["ice"].shape[0]) * scale[0, 0] + offset[0, 0]
    y = np.arange(mean_concentration_dict["ice"].shape[1]) * scale[0, 1] + offset[0, 1]
    z = np.arange(mean_concentration_dict["ice"].shape[2]) * scale[0, 2] + offset[0, 2]
    mean_concentration_water = RegularGridInterpolator((x, y, z), mean_concentration_dict["water"])(triangles_center)
    mean_concentration_surface = RegularGridInterpolator((x, y, z), mean_concentration_dict["surface"])(
        triangles_center)
    interface_type = np.where(mean_concentration_water > mean_concentration_surface, 0,
                              1)  # 0: ice-water, 1: ice-surface
    return (nodes, faces, interface_type), instantaneous_interface_dict


def post_processing_chillplus():
    filename_list = ["Hex.index", "Cubic.index", "IntIce.index"]
    indices_list = []
    for filename in filename_list:
        index_dict = read_index_dict(Path(filename))
        indices_list.append(index_dict)
    new_index_dict = combine_indices(indices_list)
    write_index_dict(new_index_dict, Path("solid_like_atoms.index"))


# def post_processing_with_PI():
#     file_path = Path("solid_like_atoms.index")
#     index_dict = read_index_dict(file_path)
#     filtered_index_dict = filter_solid_like_atoms(index_dict)
#     write_index_dict(filtered_index_dict, file_path)


def post_processing_combine_op():
    file_op = Path("op.out")
    method_list = ["chillplus"]
    content = file_op.read_text()
    first_line = content.splitlines()[0]
    columns = first_line.split()[1:]  # Skip '#'
    df = pd.read_csv("op.out", sep=r'\s+', header=None, names=columns, comment='#')
    df_method_list = []
    for method in method_list:
        path_index_file = Path(f"post_processing_{method}/solid_like_atoms.index")
        index_dict = read_index_dict(path_index_file)
        t_list = []
        n_list = []
        for t, indices in index_dict.items():
            t_list.append(float(t))
            n_list.append(len(indices))
        df_method = pd.DataFrame({"t[ps]": t_list, f"lambda_{method}": n_list})
        df_method_list.append(df_method)
    for df_method in df_method_list:
        df = df.merge(df_method, on="t[ps]", how="inner")
    df.to_csv("op_combined.csv", index=False)


def post_processing_interface():
    current_path = Path.cwd()
    u = mda.Universe("../../conf.gro", "trajout.xtc")
    x_max, y_max, z_max = u.dimensions[:3]
    step = 1.0
    x_range, y_range, z_range = (7 - 2, x_max - 7 + 2), (7 - 2, y_max - 7 + 2), (40 + 2, 70 + 2)
    pos_grid, scale, offset = generate_grid(x_range, y_range, z_range, step=step)
    index_path = Path("post_processing_chillplus/solid_like_atoms.index")
    index_dict = filter_solid_like_atoms(read_index_dict(index_path))

    with open(current_path.parent.parent / "job_params.json", 'r') as file:
        job_params = json.load(file)
    job_name = current_path.name
    params = job_params[job_name]

    t_start = params["RAMP_TIME"] + 200
    t_end = params["RAMP_TIME"] + params["PRD_TIME"]
    t_range = (t_start, t_end)
    mean_interface, instantaneous_interface_dict = calc_interface_in_t_range(u, index_dict, pos_grid,
                                                                             scale, offset, t_range)
    with open("interface.pickle", "wb") as file:
        pickle.dump(mean_interface, file)
    with open("instantaneous_interface.pickle", "wb") as file:
        pickle.dump(instantaneous_interface_dict, file)


def main():
    parser = argparse.ArgumentParser(description="Post processing.")
    parser.add_argument("--chillplus", action="store_true", help="Post processing for chillplus.")
    parser.add_argument("--combine_op", action="store_true", help="Combine op.")
    parser.add_argument("--interface", action="store_true", help="Calculate mean interface.")

    args = parser.parse_args()
    if args.chillplus:
        print("Processing chillplus...")
        post_processing_chillplus()
        print("Done.")
    if args.combine_op:
        print("Combining op...")
        post_processing_combine_op()
        print("Done.")
    if args.interface:
        print("Calculating mean interface...")
        post_processing_interface()
        print("Done.")


if __name__ == "__main__":
    main()
