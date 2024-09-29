# Author: Mian Qin
# Date Created: 9/24/24
from pathlib import Path
import argparse


def read_index_dict(file_path: Path) -> dict[str, list[str]]:
    index_dict = {}
    with open(file_path) as file:
        for line in file:
            line = line.strip().split()
            if len(line) == 0:  # End of file
                break
            t = float(line[0])
            indices = [str(x) for x in line[1:]]
            index_dict[f"{t:.1f}"] = indices
    return index_dict


def write_index_dict(index_dict: dict[str, list[str]], file_path: Path):
    with open(file_path, 'w') as file:
        for t, indices in index_dict.items():
            file.write(f"{t} {' '.join(indices)}\n")
    print(f"Write indices to {file_path}")


def combine_indices(indices_list: list[dict], allow_duplicate=False):
    new_index_dict = {}
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


def filter_solid_like_atoms(solid_like_atoms_dict: dict):
    for k, v, in solid_like_atoms_dict.items():
        for i in range(len(v)):
            if int(v[i]) > 11892:
                solid_like_atoms_dict[k] = v[:i]
                break
    return solid_like_atoms_dict


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


def main():
    parser = argparse.ArgumentParser(description="Post processing.")
    parser.add_argument("--chillplus", action="store_true", help="Post processing for chillplus.")
    parser.add_argument("--with_PI", action="store_true", help="Post processing for with_PI.")

    args = parser.parse_args()
    if args.chillplus:
        print("Processing chillplus...")
        post_processing_chillplus()
        print("Done.")
    if args.with_PI:
        print("Processing with_PI...")
        post_processing_with_PI()
        print("Done.")


if __name__ == "__main__":
    main()
