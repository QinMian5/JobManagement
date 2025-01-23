# Author: Mian Qin
# Date Created: 1/22/25
from pathlib import Path

import numpy as np
import MDAnalysis as mda
import freud

from .post_processing import generate_graph


def main():
    u = mda.Universe("./data/conf.gro")
    dimensions = u.dimensions[:3]
    box = freud.box.Box.from_box(dimensions)
    O_water = u.select_atoms("resname SOL and name OW")
    pos = O_water.positions
    aq = freud.locality.AABBQuery(box, pos)
    r_max = 3.5
    nlist = aq.query(pos, {"r_max": r_max, "r_min": 0.1}).toNeighborList()
    print()


if __name__ == "__main__":
    main()
