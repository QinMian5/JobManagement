# Author: Mian Qin
# Date Created: 1/22/25
from pathlib import Path

import numpy as np
import MDAnalysis as mda
import freud


def main():
    u = mda.Universe("./data/conf.gro")
    dimensions = u.dimensions[:3]
    box = freud.box.Box.from_box(dimensions)
    O_ag = u.select_atoms("resname SOL and name OW")
    pos = O_ag.positions
    aq = freud.locality.AABBQuery(box, pos)
    r_max = 3.0
    nlist = aq.query(pos, {"r_max": r_max}).toNeighborList()
    print(nlist)


if __name__ == "__main__":
    main()
