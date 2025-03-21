# Author: Mian Qin
# Date Created: 6/8/24
import os
from pathlib import Path
import shutil
from itertools import count
import subprocess
import shlex
import json


run_env = os.environ.get("RUN_ENVIRONMENT")
if run_env == "wsl":
    home_path = Path("/home/qinmian")
elif run_env == "chestnut":
    home_path = Path("/home/mianqin")
elif run_env == "mianqin_PC":
    home_path = Path("//wsl.localhost/Ubuntu-22.04/home/qinmian")
else:
    raise RuntimeError(f"Unknown environment: {run_env}")


def write(rho: float):
    root_dir = home_path / "data/gromacs/pseudoice/topol"
    path_template = root_dir / "topol_template.top"
    path_topol = root_dir / f"topol_{rho}.top"
    template_content = path_template.read_text()

    H_charge = 0.5897
    H_charge_PI = f"{H_charge * rho:.4f}"
    O_charge_PI = f"{-2 * float(H_charge_PI):.4f}"

    content = template_content
    content = content.replace("PI.RHO", str(rho))
    content = content.replace("PI.H.CHARGE", H_charge_PI)
    content = content.replace("PI.O.CHARGE", O_charge_PI)
    path_topol.write_text(content)


def main():
    for rho in [0.0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1.0]:
        write(rho)


if __name__ == "__main__":
    main()
