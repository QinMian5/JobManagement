# Author: Mian Qin
# Date Created: 9/23/24
from pathlib import Path
import shutil
import subprocess
import shlex
import json
import argparse
import logging
from itertools import count

logging.basicConfig(level=logging.INFO)


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
        content = content.replace(f"${{{k}}}", str(v))
    return content


class GromacsJob:
    templates: dict[str, str] = {}  # {src: dst}
    files_to_gather = {}  # {src: dst}
    data_dir = Path("./data")
    result_dir = Path("./result")
    template_dir = Path("./template")

    def __init__(self, name, params: dict):
        self.name = name
        self.params = params
        self.root_dir = self.data_dir / name

    def create(self):
        logging.info(f"Creating job folder for {self.name}")
        self.root_dir.mkdir(exist_ok=True)
        for filename_src, filename_dst in self.templates.items():
            template_content = (self.template_dir / filename_src).read_text()
            content = _replace_params(template_content, self.params)
            destination = self.root_dir / filename_dst
            destination.parent.mkdir(exist_ok=True)
            destination.write_text(content)
            logging.debug(f"Written {destination}")

    def submit(self):
        logging.info(f"Submitting job.sh for {self.name}")
        command = f"sbatch job.sh"
        subprocess.run(shlex.split(command), cwd=self.root_dir)

    def gather(self):
        logging.info(f"Gathering results for {self.name}")
        destination_dir = self.result_dir / self.name
        destination_dir.mkdir(exist_ok=True)
        # Copy necessary files
        for filename_src, filename_dst in self.files_to_gather:
            src = self.root_dir / filename_src
            dst = destination_dir / filename_dst
            shutil.copyfile(src, dst)
            logging.debug(f"Copied {src} to {dst}")

    def clean(self):
        logging.info(f"Cleaning up for {self.name}")
        backup_files = self.root_dir.glob("#*")
        for backup_file in backup_files:
            backup_file.unlink()
            logging.debug(f"Deleted {backup_file}")


def main():
    ...


if __name__ == "__main__":
    main()
