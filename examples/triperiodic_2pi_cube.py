"""Example: droplets in a tri-periodic (2*pi)^3 cube.

Grid 256^3 (dx = 2*pi/256) with 32 points per diameter -> d = pi/4 (about 0.785).
Droplets stay fully inside the box on every axis (they never straddle a periodic face),
and droplet spacing is checked across periodic faces as well.

Writes spheres.in to the current directory, so run it from the case directory,
or pass an output path:
    python /path/to/generateDroplets/examples/triperiodic_2pi_cube.py [outfile]
Equivalent CLI call:
    python main.py --bc triperiodic --lx 6.283185307179586 --ly 6.283185307179586 \
        --lz 6.283185307179586 --nx 256 --ny 256 --nz 256 --ppd 32 --alpha 0.1 \
        --method boxed --seed 1 --outfile spheres.in
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import droplets  # noqa: E402
from config import Config  # noqa: E402

L = 2.0 * np.pi
cfg = Config(
    bc="triperiodic",
    lx=L, ly=L, lz=L,
    nx=256, ny=256, nz=256,
    ppd=32,
    alpha=0.1,
    method="boxed",
    epsf=0.5,
    gap_x=0.12, gap_y=0.12, gap_z=0.12,
    drop_gap=0.4,
    seed=1,
    outfile=sys.argv[1] if len(sys.argv) > 1 else "spheres.in",
).validate()

print(cfg.summary())
pos = droplets.generate(cfg, np.random.default_rng(cfg.seed))
droplets.write_spheres(cfg.outfile, pos, cfg.R)
print(f"written {len(pos)} droplets to {cfg.outfile}")
pos_file, radii = droplets.read_spheres(cfg.outfile)  # validate what was written
print(droplets.check_config(pos_file, cfg, radii))
