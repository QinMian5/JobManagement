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
import MDAnalysis as mda
from skimage import measure


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
            if int(v[i]) > 15784:
                solid_like_atoms_dict[k] = v[:i]
                break
    return solid_like_atoms_dict


def calc_scale_offset(x_range: tuple[float, float], n_x: int) \
        -> tuple[float, float]:
    start, end = x_range
    scale = (end - start) / (n_x - 1)
    offset = start
    return scale, offset


def generate_grid(x_range: tuple[float, float], y_range: tuple[float, float], z_range: tuple[float, float], step: float):
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


def calc_concentration(pos_atom: np.ndarray, pos_grid: np.ndarray,
                       box_size: np.ndarray, ksi: float) -> np.ndarray:
    """

    :param pos_atom: Positions of O atoms. Shape: (N, 3)
    :param pos_grid: Positions of grids. Shape: (X, Y, Z, 3)
    :param box_size: Size of the box. Shape: (3)
    :param ksi:
    :return: Shape: (X, Y ,Z)
    """
    cutoff = 3 * ksi
    grid_size = pos_grid[-1, -1, -1] - pos_grid[0, 0, 0]
    n_x, n_y, n_z = np.ceil(grid_size / (cutoff * 2)).astype(int)
    concentrations = [[[None for _ in range(n_z)] for _ in range(n_y)] for _ in range(n_x)]
    sub_grids = np.array_split(pos_grid, n_x, axis=0)  # 沿 x 轴划分
    sub_grids = [np.array_split(sub_grid, n_y, axis=1) for sub_grid in sub_grids]  # 沿 y 轴划分
    sub_grids = [[np.array_split(sub_sub_grid, n_z, axis=2) for sub_sub_grid in sub_grid] for sub_grid in sub_grids]

    for i, j, k in product(range(n_x), range(n_y), range(n_z)):
        sub_grid = sub_grids[i][j][k]
        corner1 = sub_grid[0, 0, 0] - cutoff
        corner2 = sub_grid[-1, -1, -1] + cutoff
        pos_vec_to_corner1 = wrap_pos_vec_array(pos_atom - corner1, box_size, mode="positive")
        pos_vec_to_corner2 = wrap_pos_vec_array(pos_atom - corner2, box_size, mode="negative")
        in_box1 = np.all(pos_vec_to_corner1 < corner2 - corner1, axis=1)
        in_box2 = np.all(pos_vec_to_corner2 > -(corner2 - corner1), axis=1)
        in_box = np.logical_and(in_box1, in_box2)
        sub_atom = pos_atom[in_box]

        pos_vec_array = np.expand_dims(sub_grid, axis=3) - np.expand_dims(sub_atom, axis=(0, 1, 2))
        wrapped_pos_vec_array = wrap_pos_vec_array(pos_vec_array, box_size)
        distance_array = np.linalg.norm(wrapped_pos_vec_array, axis=4)
        concentration = np.sum(1 / ((np.sqrt(2 * np.pi) * ksi) ** 3) * np.exp(-distance_array ** 2 / (2 * ksi ** 2)), axis=3)
        concentrations[i][j][k] = concentration
    concentrations = np.concatenate([np.concatenate([np.concatenate(sub_sub_con, axis=2) for sub_sub_con in sub_con], axis=1) for sub_con in concentrations], axis=0)
    return concentrations


def generate_interface_from_concentration(concentration: np.ndarray, level=0.015) \
        -> tuple[np.ndarray, np.ndarray]:
    if level > concentration.max() or level < concentration.min():
        nodes, faces = np.zeros((0, 3)), np.zeros((0, 3))
    else:
        nodes, faces, normals, values = measure.marching_cubes(concentration, level=level)
    return nodes, faces


def calc_interface_in_t_range(u: mda.Universe, solid_like_atoms_dict: dict[str, list[str]], pos_grid, scale,
                              offset, t_range: tuple[float, float], ksi=3.5 / 2):
    instantaneous_interface_dict = {}
    acc_concentration = np.zeros(pos_grid.shape[:3])
    acc_n = 0
    for ts in u.trajectory[::10]:  # 1 frame / 100 ps
        t = ts.time
        solid_like_atoms_id = solid_like_atoms_dict[f"{t:.1f}"]
        if solid_like_atoms_id:
            solid_like_atoms = u.select_atoms(f"bynum {' '.join(solid_like_atoms_id)}")
            pos_ice = solid_like_atoms.positions
            concentration = calc_concentration(pos_ice, pos_grid, u.dimensions[:3], ksi)
            # TODO
            nodes, faces = generate_interface_from_concentration(concentration)
            nodes = nodes * scale + offset
        else:
            concentration = 0
            nodes, faces = [], []
        instantaneous_interface_dict[f"{t:.1f}"] = [nodes, faces]
        if t_range[0] <= t <= t_range[1]:
            acc_concentration += concentration
            acc_n += 1
    mean_concentration = acc_concentration / acc_n
    nodes, faces = generate_interface_from_concentration(mean_concentration)
    nodes = nodes * scale + offset
    return nodes, faces, instantaneous_interface_dict


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
    step = 0.5
    x_range, y_range, z_range = (7 - 2, x_max - 7 + 2), (7 - 2, y_max -7 + 2), (40 - 2, 70 + 2)
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
    nodes, faces, instantaneous_interface_dict = calc_interface_in_t_range(u, index_dict, pos_grid, scale, offset, t_range)
    with open("interface.pickle", "wb") as file:
        pickle.dump([nodes, faces], file)
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
