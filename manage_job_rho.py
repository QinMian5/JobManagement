# Author: Mian Qin
# Date Created: 6/28/24
from pathlib import Path
import shutil
from itertools import count
import subprocess
import shlex
import json


ROOT_DIR = Path("/home/qinmian/data/gromacs/pseudoice/data")


def copy_submit(rho):
    dst_dir = ROOT_DIR / str(rho)
    shutil.copytree(ROOT_DIR / "template", dst_dir)
    shutil.copyfile(ROOT_DIR.parent / "topol" / f"topol_{rho}.top", dst_dir / "topol.top")
    command = "sbatch job.sh"
    working_dir = dst_dir / "equilibrium"
    subprocess.run(shlex.split(command), cwd=working_dir)


def main():
    for rho in [0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9]:
        copy_submit(rho)


if __name__ == "__main__":
    main()
