# Author: Mian Qin
# Date Created: 9/2/24
from pathlib import Path
import subprocess
import shlex

import scipy.constants as c


class PseudoWaterGenerator:
    _root_dir = Path("/home/qinmian/data/gromacs/pseudoice/generate_conf")
    _filepath_packmol_template = _root_dir / "template" / "packmol_template.txt"
    _filepath_tip4pice_pdb = _root_dir / "template" / "tip4pice.pdb"

    _default_working_dir = Path("/home/qinmian/data/gromacs/pseudoice/generate_conf/pseudo_water")
    _folder_name_preparation = "preparation"
    _filename_packmol_input = "packmol.txt"

    def __init__(self, working_dir=None):
        if working_dir is None:
            self.working_dir = self._default_working_dir
        self.working_dir.mkdir(parents=True, exist_ok=True)

    def generate_box_of_water(self, x, y, z):
        """

        :param x: In Angstrom.
        :param y: In Angstrom.
        :param z: In Angstrom.
        """
        volume = x * y * z * 1e-30  # Convert to m^3
        rho = 5e4 * c.N_A  # $5\times10^4$ mol/m^3
        N_water = int(rho * volume)
        packmol_content = self._filepath_packmol_template.read_text()
        packmol_content = packmol_content.replace("N_WATER", f"{N_water}")
        packmol_content = packmol_content.replace("BOX_X", f"{x}")
        packmol_content = packmol_content.replace("BOX_Y", f"{y}")
        packmol_content = packmol_content.replace("BOX_Z", f"{z}")

        save_dir = self.working_dir / self._folder_name_preparation
        save_dir.mkdir(exist_ok=True)
        save_path = save_dir / self._filename_packmol_input
        save_path.write_text(packmol_content)

        with open(save_path) as file:
            result = subprocess.run("packmol", stdin=file, cwd=save_dir)


def main():
    generator = PseudoWaterGenerator()
    generator.generate_box_of_water(40, 40, 40)


if __name__ == "__main__":
    main()
