# Author: Mian Qin
# Date Created: 9/24/24
from pathlib import Path
import argparse
import json
import pickle
from collections import OrderedDict

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
            if int(v[i]) > 32768:
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
                  n_x: int = 50, n_y: int = 50, n_z: int = 50):
    x_scale, x_offset = calc_scale_offset(x_range, n_x)
    y_scale, y_offset = calc_scale_offset(y_range, n_y)
    z_scale, z_offset = calc_scale_offset(z_range, n_z)
    scale = np.array([x_scale, y_scale, z_scale]).reshape(1, 3)
    offset = np.array([x_offset, y_offset, z_offset]).reshape(1, 3)
    x, y, z = np.linspace(*x_range, n_x), np.linspace(*y_range, n_y), np.linspace(*z_range, n_z)
    X, Y, Z = np.meshgrid(x, y, z)
    pos_grid = np.stack((X, Y, Z), axis=3)
    return pos_grid, scale, offset


def wrap_position_vector_array(position_vector_array: np.ndarray, box_size: np.ndarray):
    half_box_size = box_size / 2
    lt = position_vector_array < -half_box_size
    gt = position_vector_array > half_box_size
    position_vector_array += lt * box_size
    position_vector_array -= gt * box_size
    return position_vector_array


def calc_concentration(pos_ice: np.ndarray, pos_grid: np.ndarray,
                       box_size: np.ndarray, ksi: float) -> np.ndarray:
    """

    :param pos_ice: Positions of O atoms in ice-like molecules. Shape: (N, 3)
    :param pos_grid: Positions of grids. Shape: (X, Y, Z, 3)
    :param box_size: Size of the box. Shape: (3)
    :param ksi:
    :return: Shape: (X, Y ,Z)
    """
    position_vector_array = np.expand_dims(pos_grid, axis=3) - np.expand_dims(pos_ice, axis=(0, 1, 2))
    wrapped_position_vector_array = wrap_position_vector_array(position_vector_array, box_size)
    distance_array = np.linalg.norm(wrapped_position_vector_array, axis=4)
    concentration = np.sum(1 / ((np.sqrt(2 * np.pi) * ksi) ** 3) * np.exp(-distance_array ** 2 / (2 * ksi ** 2)),
                           axis=3)
    return concentration


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
    x_range, y_range, z_range = (0, x_max), (0, y_max), (0, z_max)
    pos_grid, scale, offset = generate_grid(x_range, y_range, z_range, n_x=70, n_y=70, n_z=70)
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
