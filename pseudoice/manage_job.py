# Author: Mian Qin
# Date Created: 5/29/24
from pathlib import Path
import shutil
from itertools import count
import subprocess
import shlex
import json
import argparse


_templates = {
    "job.sh": "job.sh",
    "grompp.mdp": "grompp.mdp",
    "order_parameters.dat": "order_parameters.dat",
    "post_processing_qbar.dat": "post_processing_qbar/post_processing_qbar.dat",
    "post_processing_with_PI.dat": "post_processing_with_PI/post_processing_with_PI.dat",
    "post_processing_chillplus.dat": "post_processing_chillplus/post_processing_chillplus.dat"
}
_files_to_gather = {
    "op.out": "op.out",
    "trajout.xtc": "trajout.xtc",
    "post_processing_qbar/solid_like_atoms_qbar.index": "solid_like_atoms_qbar.index",
    "post_processing_with_PI/solid_like_atoms_with_PI.index": "solid_like_atoms_with_PI.index",
    "post_processing_chillplus/solid_like_atoms_chillplus.index": "solid_like_atoms_chillplus.index"
}
_filename_job = "job.sh"
_path_data_dir = Path("./data")
_path_result_dir = Path("./result")
_path_template_dir = Path("./template")


def _load_params() -> dict[str, dict]:
    with open("job_params.json", 'r') as file:
        job_params = json.load(file)
    return job_params


def _backup(path: Path) -> None:
    if path.exists():
        for i in count(start=1):
            filename = path.name
            backup_name = f"#{filename}.{i}#"
            backup_path = path.parent / backup_name
            if not backup_path.exists():
                shutil.move(path, backup_path)
                print(f"Backed up {path} to {backup_path}")
                return


def _flatten_dict(d: dict, parent_key="", sep='.') -> dict:
    flattened_dict = {}
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, dict):
            flattened_dict.update(_flatten_dict(v, new_key))
        else:
            flattened_dict[new_key] = v
    return flattened_dict


def _replace_params(content: str, params: dict) -> str:
    flattened_dict = _flatten_dict(params)
    for k, v in flattened_dict.items():
        content = content.replace(k, str(v))
    return content


def create(backup=True) -> None:
    _write_params()

    if backup:
        _backup(_path_data_dir)
        _path_data_dir.mkdir()

    path_template_op = _path_template_dir / _filename_op
    template_content_op = path_template_op.read_text()
    path_template_pp = _path_template_dir / _filename_pp
    template_content_pp = path_template_pp.read_text()
    path_template_pp2 = _path_template_dir / _filename_pp2
    template_content_pp2 = path_template_pp2.read_text()
    path_template_pp3 = _path_template_dir / _filename_pp3
    template_content_pp3 = path_template_pp3.read_text()
    path_template_job = _path_template_dir / _filename_job
    template_content_job = path_template_job.read_text()
    path_template_job_post = _path_template_dir / _filename_job_post
    template_content_job_post = path_template_job_post.read_text()

    job_params = _load_params()
    for name, params in job_params.items():
        path_dir = _path_data_dir / f"{name}"
        path_dir.mkdir(exist_ok=True)

        path_op = path_dir / _filename_op
        path_op.write_text(_replace_params(template_content_op, params))
        path_pp = path_dir / _folder_name_pp / _filename_pp
        path_pp.parent.mkdir(exist_ok=True)
        path_pp.write_text(_replace_params(template_content_pp, params))
        path_pp2 = path_dir / _folder_name_pp2 / _filename_pp2
        path_pp2.parent.mkdir(exist_ok=True)
        path_pp2.write_text(_replace_params(template_content_pp2, params))
        path_pp3 = path_dir / _folder_name_pp3 / _filename_pp3
        path_pp3.parent.mkdir(exist_ok=True)
        path_pp3.write_text(_replace_params(template_content_pp3, params))
        path_job = path_dir / _filename_job
        path_job.write_text(_replace_params(template_content_job, params))
        path_job_post = path_dir / _filename_job_post
        path_job_post.write_text(_replace_params(template_content_job_post, params))


def submit() -> None:
    job_params = _load_params()
    for name in job_params.keys():
        working_dir = _path_data_dir / f"{name}"
        command = f"sbatch {_filename_job}"
        subprocess.run(shlex.split(command), cwd=working_dir)


def submit_post() -> None:
    job_params = _load_params()
    for name in job_params.keys():
        working_dir = _path_data_dir / f"{name}"
        command = f"sbatch {_filename_job_post}"
        subprocess.run(shlex.split(command), cwd=working_dir)


def gather() -> None:
    _backup(_path_result_dir)
    _path_result_dir.mkdir()
    shutil.copyfile("./conf.gro", _path_result_dir / "conf.gro")
    shutil.copyfile(_filename_job_params, _path_result_dir / _filename_job_params)

    job_params = _load_params()
    for name in job_params.keys():
        path_src_dir = _path_data_dir / f"{name}"
        path_dst_dir = _path_result_dir / f"{name}"
        path_dst_dir.mkdir()
        # Copy result
        shutil.copyfile(path_src_dir / _filename_op_out, path_dst_dir / _filename_op_out)
        shutil.copyfile(path_src_dir / _filename_trj, path_dst_dir / _filename_trj)
        shutil.copyfile(path_src_dir / _folder_name_pp / _filename_index, path_dst_dir / _filename_index)
        shutil.copyfile(path_src_dir / _folder_name_pp2 / _filename_index, path_dst_dir / _filename_index2)
        shutil.copyfile(path_src_dir / _folder_name_pp3 / _filename_index, path_dst_dir / _filename_index3)


def clean() -> None:
    job_params = _load_params()
    for backup_folder in Path("..").glob("#*"):
        shutil.rmtree(backup_folder)
        print(f"Deleted {backup_folder.resolve()}")
    for job_name in job_params.keys():
        working_dir = _path_data_dir / job_name
        for backup_file in working_dir.glob("#*"):
            backup_file.unlink()
            print(f"Deleted {backup_file.resolve()}")
        for old_log_file in working_dir.glob(r"slurm*"):
            old_log_file.unlink()
            print(f"Deleted {old_log_file.resolve()}")
        old_pp_file_list = ["F_q.out", "post_processing_qbar.dat", "post_processing_with_PI.dat", "solid_like_atoms.index",
                            "solid_like_atoms_with_PI.index", "time_samples_q.out"]
        for old_pp_file in old_pp_file_list:
            old_file = working_dir / old_pp_file
            if old_file.exists():
                old_file.unlink()
                print(f"Deleted {old_file.resolve()}")


def main():
    parser = argparse.ArgumentParser(description="Control the program functions.")
    parser.add_argument("--create", action="store_true", help="Execute the create function")
    parser.add_argument("--no_backup", action="store_true", help="Whether to backup during creating")
    parser.add_argument("--submit", action="store_true", help="Execute the submit function")
    parser.add_argument("--submit_post", action="store_true", help="Execute the submit_post function")
    parser.add_argument("--gather", action="store_true", help="Execute the gather function")
    parser.add_argument("--clean", action="store_true", help="Execute the clean function")

    args = parser.parse_args()
    if args.create:
        print("Creating job folders...")
        backup = not args.no_backup
        create(backup=backup)
        print("Done.")
    if args.submit:
        print("Submitting job.sh...")
        submit()
        print("Done.")
    if args.submit_post:
        print("Submitting job_post.sh...")
        submit_post()
        print("Done.")
    if args.gather:
        print("Gathering results...")
        gather()
        print("Done.")
    if args.clean:
        print("Cleaning backup files...")
        clean()
        print("Done.")


if __name__ == "__main__":
    main()
