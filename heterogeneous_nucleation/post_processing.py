# Author: Mian Qin
# Date Created: 9/24/24
from pathlib import Path
import argparse
import json
import pickle
from collections import OrderedDict
import inspect

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
import MDAnalysis as mda
from skimage import measure
import freud
import networkx as nx


def read_index_dict(file_path: Path) -> OrderedDict[str, list[str]]:
    """
    Read index file and return ordered dictionary with timestamps as keys.

    :param file_path: Path to index file containing timestamps and atom indices.
    :return: OrderedDict where keys are timestamps (str) and values are lists of atom indices (str).
    """
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
    """
    Generate water molecule indices from given ice-like atom indices.

    :param indices: List of ice-like atom indices.
    :return: List of water molecule indices not present in the input solid-like indices.
    """
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
    """
    Write index dictionary to file.

    :param index_dict: Dictionary with timestamps as keys and atom index lists as values.
    :param file_path: Output file path for writing the index data.
    :return:
    """
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


def generate_graph(pos_atoms: np.ndarray, ids: np.ndarray, box_size, type2indices: dict[str, list[int]],
                   threshold_edge=3.5) -> nx.Graph:
    """
    Generate molecular graph based on spatial proximity.

    :param pos_atoms: Positions of atoms. Shape: (N, 3)
    :param ids: Indices of atoms. Shape: (N)
    :param box_size:
    :param type2indices: Ice/Water/Surface
    :param threshold_edge: In Angstrom
    :return:
    """
    box = freud.box.Box.from_box(box_size)
    aq = freud.locality.AABBQuery(box, pos_atoms)
    r_max = threshold_edge
    nlist = aq.query(pos_atoms, {"r_max": r_max, "r_min": 0.1}).toNeighborList()
    i_array = nlist.point_indices
    j_array = nlist.query_point_indices
    non_duplicate_index = np.where(i_array < j_array)
    i_array = i_array[non_duplicate_index]
    j_array = j_array[non_duplicate_index]

    i_ids_array = ids[i_array]
    j_ids_array = ids[j_array]

    G = nx.Graph()
    for t, indices in type2indices.items():
        G.add_nodes_from(indices, type=t)
    G.add_edges_from(list(zip(i_ids_array, j_ids_array)))
    return G


def update_graph(G: nx.Graph, type2indices: dict):
    for t, indices in type2indices.items():
        for index in indices:
            G.nodes[index]["type"] = t
    return G


def _correct_bulk(G: nx.Graph, max_iter=10) -> tuple[np.ndarray, np.ndarray]:
    _type_to_index = {"ice": 0, "water": 1, "surface": 2}

    new_ice_array = np.zeros(len(G.nodes))
    new_water_array = np.zeros(len(G.nodes))
    node_type_array = np.zeros((len(G.nodes), len(_type_to_index)), dtype=int)
    for i, node in enumerate(G.nodes):
        node_type = G.nodes[node]["type"]
        node_type_array[i][_type_to_index[node_type]] = 1

    adj_matrix_first_shell = nx.to_scipy_sparse_array(G)
    adj_matrix2 = adj_matrix_first_shell @ adj_matrix_first_shell
    adj_matrix_two_shell = adj_matrix_first_shell + adj_matrix2
    adj_matrix_two_shell.setdiag(0)  # Exclude self
    adj_matrix_two_shell.data = np.minimum(adj_matrix_two_shell.data, 1)
    adj_matrix_two_shell.eliminate_zeros()
    is_surface = node_type_array[:, _type_to_index["surface"]]
    n_surface_two_shell = adj_matrix_two_shell @ is_surface
    cond_surface = n_surface_two_shell >= 1
    for i in range(max_iter):
        # adj_matrix_second_shell = adj_matrix_two_shell - adj_matrix_first_shell
        # adj_matrix_second_shell.eliminate_zeros()
        is_ice = node_type_array[:, _type_to_index["ice"]]
        is_water = node_type_array[:, _type_to_index["water"]]
        n_ice_two_shell = adj_matrix_two_shell @ is_ice
        n_water_two_shell = adj_matrix_two_shell @ is_water
        n_mobile_two_shell = n_ice_two_shell + n_water_two_shell
        cond_bulk = ~cond_surface & (n_mobile_two_shell >= 8)

        water_to_ice = is_water & cond_bulk & (n_water_two_shell <= 1)
        ice_to_water = is_ice & cond_bulk & (n_ice_two_shell <= 1)

        if water_to_ice.sum() + ice_to_water.sum() == 0:  # converge
            break
        node_type_array[:, _type_to_index["ice"]] += water_to_ice
        node_type_array[:, _type_to_index["water"]] -= water_to_ice
        node_type_array[:, _type_to_index["water"]] += ice_to_water
        node_type_array[:, _type_to_index["ice"]] -= ice_to_water
        new_ice_array += water_to_ice
        new_water_array += ice_to_water
    else:
        print(f"{inspect.currentframe().f_code.co_name}: Failed to converge after {max_iter} iterations.")
    node_index_array = np.array(G.nodes())
    water_to_ice_array = node_index_array[new_ice_array.astype(bool)]
    ice_to_water_array = node_index_array[new_water_array.astype(bool)]
    return water_to_ice_array, ice_to_water_array


def _correct_surface_water2ice(G: nx.Graph, max_iter=20) -> np.ndarray:
    _type_to_index = {"ice": 0, "water": 1, "surface": 2, "maybe_ice": 3}

    new_ice_array = np.zeros(len(G.nodes)).astype(bool)
    node_type_array = np.zeros((len(G.nodes), len(_type_to_index)), dtype=bool)
    for i, node in enumerate(G.nodes):
        node_type = G.nodes[node]["type"]
        node_type_array[i][_type_to_index[node_type]] = True

    adj_matrix_first_shell = nx.to_scipy_sparse_array(G)
    adj_matrix2 = adj_matrix_first_shell @ adj_matrix_first_shell
    adj_matrix3 = adj_matrix2 @ adj_matrix_first_shell
    adj_matrix_two_shell = adj_matrix_first_shell + adj_matrix2
    adj_matrix_two_shell.setdiag(0)  # Exclude self
    adj_matrix_two_shell.data = np.minimum(adj_matrix_two_shell.data, 1)
    adj_matrix_two_shell.eliminate_zeros()
    adj_matrix_three_shell = adj_matrix_two_shell + adj_matrix3
    adj_matrix_three_shell.setdiag(0)
    adj_matrix_three_shell.data = np.minimum(adj_matrix_three_shell.data, 1)
    adj_matrix_three_shell.eliminate_zeros()
    is_surface = node_type_array[:, _type_to_index["surface"]]
    n_surface_two_shell = adj_matrix_two_shell @ is_surface
    near_surface = n_surface_two_shell >= 1

    update_n_list = []
    for i in range(max_iter):
        is_ice = node_type_array[:, _type_to_index["ice"]]
        is_water = node_type_array[:, _type_to_index["water"]]
        n_ice_two_shell = adj_matrix_two_shell @ is_ice
        water_to_ice = is_water & near_surface & (n_ice_two_shell >= 8)
        update_n = np.sum(water_to_ice)
        if update_n == 0:
            break
        update_n_list.append(update_n)
        new_ice_array |= water_to_ice
        node_type_array[:, _type_to_index["ice"]] |= water_to_ice
        node_type_array[:, _type_to_index["water"]] &= ~water_to_ice
    else:
        print(f"{inspect.currentframe().f_code.co_name}: Failed to converge after {max_iter} iterations.")
        print(update_n_list)
    node_index_array = np.array(G.nodes())
    water_to_ice_array = node_index_array[new_ice_array.astype(bool)]
    return water_to_ice_array


def correct_ice_index(path_to_conf: Path, path_to_traj: Path, path_to_index: Path, correct_surface=True,
                      correct_bulk=True) -> tuple[OrderedDict, OrderedDict]:
    u = mda.Universe(path_to_conf, path_to_traj)
    ice_index_dict = filter_solid_like_atoms(read_index_dict(path_to_index))
    new_ice_index_dict = OrderedDict()
    new_water_index_dict = OrderedDict()
    for ts in u.trajectory:
        if int(float(ts.time)) % 500 == 0:
            print(ts.time)
        ice_index = ice_index_dict[f"{ts.time:.1f}"]
        O_atoms_all = u.select_atoms("name OW")
        if len(ice_index) == 0:
            O_atoms_water = O_atoms_all
            O_atoms_ice = O_atoms_all - O_atoms_water
        else:
            O_atoms_ice = u.select_atoms(f"bynum {' '.join(ice_index)}")
            O_atoms_water = O_atoms_all - O_atoms_ice
        O_atoms_surface = u.select_atoms("name O_PI")
        pos_ice = O_atoms_ice.positions
        pos_water = O_atoms_water.positions
        pos_surface = O_atoms_surface.positions
        box_size = u.dimensions[:3]
        pos_atoms = np.concatenate((pos_ice, pos_water, pos_surface))
        ids = np.concatenate((O_atoms_ice.ids, O_atoms_water.ids, O_atoms_surface.ids))
        type2indices = {"ice": O_atoms_ice.ids,
                        "water": O_atoms_water.ids,
                        "surface": O_atoms_surface.ids}
        G = generate_graph(pos_atoms, ids, box_size, type2indices)

        new_water_set = set()
        new_ice_set = set()
        if correct_bulk:
            bulk_water_to_ice, bulk_ice_to_water = _correct_bulk(G)
            # Update sets with bulk corrections
            new_ice_set.update(bulk_water_to_ice)
            new_water_set.update(bulk_ice_to_water)
            # Remove any conflicts (water->ice should not be in water set)
            new_water_set -= set(bulk_water_to_ice)
            new_ice_set -= set(bulk_ice_to_water)
            # Update graph with current changes
            G = update_graph(G, {
                "ice": list(bulk_water_to_ice),
                "water": list(bulk_ice_to_water)
            })
        if correct_surface:
            surface_water_to_ice = _correct_surface_water2ice(G)
            # Update sets with surface corrections
            new_ice_set.update(surface_water_to_ice)
            # Remove converted waters from water set
            new_water_set -= set(surface_water_to_ice)
            G = update_graph(G, {
                "ice": list(surface_water_to_ice),
            })

        new_ice_index = [str(x) for x in sorted(new_ice_set, key=int)]
        new_water_index = [str(x) for x in sorted(new_water_set, key=int)]

        new_ice_index_dict[f"{ts.time:.1f}"] = new_ice_index
        new_water_index_dict[f"{ts.time:.1f}"] = new_water_index
    return new_ice_index_dict, new_water_index_dict


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
    assert len(pos_atoms) != 0
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
    if len(pos_atom) == 0:
        concentrations = np.zeros(shape)
    else:
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
    avg_centroid = u.dimensions[:3] / 2

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
        surface_atoms = u.select_atoms("name O_PI")
        atoms_dict = {"ice": ice_atoms, "water": water_atoms, "surface": surface_atoms}
        pos_dict = {k: v.positions for k, v in atoms_dict.items()}
        if pos_dict["ice"].shape[0] > 0:
            current_centroid = ice_atoms.center_of_geometry()
            shift = avg_centroid - current_centroid
            shift[2] = 0  # do not shift z
            for k, v in pos_dict.items():
                pos_dict[k] = v + shift
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


def post_processing_with_PI():
    file_path = Path("solid_like_atoms.index")
    index_dict = read_index_dict(file_path)
    filtered_index_dict = filter_solid_like_atoms(index_dict)
    write_index_dict(filtered_index_dict, file_path)


def post_processing_correct_ice():
    path_to_conf = Path("../../conf.gro")
    path_to_traj = Path("./trajout.xtc")
    method_list = ["with_PI"]
    for method in method_list:
        method_folder = Path(f"./post_processing_{method}")
        path_to_index = method_folder / "solid_like_atoms.index"
        new_ice_save_path = method_folder / "new_ice.index"
        new_water_save_path = method_folder / "new_water.index"
        ice_save_path = method_folder / "corrected_ice.index"
        water_save_path = method_folder / "corrected_water.index"
        new_ice_index_dict, new_water_index_dict = correct_ice_index(path_to_conf, path_to_traj, path_to_index)

        write_index_dict(new_ice_index_dict, new_ice_save_path)
        write_index_dict(new_water_index_dict, new_water_save_path)

        ice_index_dict = filter_solid_like_atoms(read_index_dict(path_to_index))
        corrected_ice_dict = combine_indices([ice_index_dict, new_ice_index_dict])
        corrected_water_index_dict = OrderedDict(
            (k, get_water_indices(ice_index)) for k, ice_index in corrected_ice_dict.items())
        write_index_dict(corrected_ice_dict, ice_save_path)
        write_index_dict(corrected_water_index_dict, water_save_path)


def post_processing_combine_op():
    file_op = Path("op.out")
    method_list = ["chillplus", "with_PI"]
    content = file_op.read_text()
    first_line = content.splitlines()[0]
    columns = first_line.split()[1:]  # Skip '#'
    df = pd.read_csv("op.out", sep=r'\s+', header=None, names=columns, comment='#')
    df_method_list = []
    for method in method_list:
        path_index_file = Path(f"post_processing_{method}/corrected_ice.index")
        column_name = f"lambda_{method}"
        if not path_index_file.exists():
            path_index_file = Path(f"post_processing_{method}/solid_like_atoms.index")
            column_name = f"{method}"
        index_dict = read_index_dict(path_index_file)
        t_list = []
        n_list = []
        for t, indices in index_dict.items():
            t_list.append(float(t))
            n_list.append(len(indices))
        df_method = pd.DataFrame({"t[ps]": t_list, column_name: n_list})
        df_method_list.append(df_method)
    for df_method in df_method_list:
        df = df.merge(df_method, on="t[ps]", how="inner")
    df.to_csv("op_combined.csv", index=False)


def post_processing_interface():
    current_path = Path.cwd()
    u = mda.Universe("../../conf.gro", "trajout.xtc")
    x_max, y_max, z_max = u.dimensions[:3]
    step = 1.0
    x_range, y_range, z_range = (0, x_max), (0, y_max), (40 - 5, 70 + 2)
    pos_grid, scale, offset = generate_grid(x_range, y_range, z_range, step=step)
    index_path = Path("post_processing_with_PI/corrected_ice.index")
    if not index_path.exists():
        index_path = Path("post_processing_with_PI/solid_like_atoms.index")
    index_dict = filter_solid_like_atoms(read_index_dict(index_path))
    print(f"Read index from {index_path}")

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
    parser.add_argument("--with_PI", action="store_true", help="Post processing for with_PI.")
    parser.add_argument("--correct_ice", action="store_true", help="Correct ice.")
    parser.add_argument("--combine_op", action="store_true", help="Combine op.")
    parser.add_argument("--interface", action="store_true", help="Calculate mean interface.")

    args = parser.parse_args()
    if args.chillplus:
        print("Processing chillplus...")
        post_processing_chillplus()
        print("Done.")
    if args.with_PI:
        print("Processing with_PI...")
        post_processing_with_PI()
        print("Done.")
    if args.correct_ice:
        print("Correcting ice...")
        post_processing_correct_ice()
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
