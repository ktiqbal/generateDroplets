"""Example: droplets in a 6 x 3 x 2 channel (periodic in x, y; walls at z = 0 and z = 2).

Grid 384 x 192 x 128 (dx = 1/64) with 16 points per diameter -> d = 0.25, R = 0.125.
Every droplet keeps a surface gap of at least 0.5 R (0.0625) to the walls.

Writes spheres.in to the current directory, so run it from the case directory,
or pass an output path:
    python /path/to/generateDroplets/examples/channel_6x3x2.py [outfile]
Equivalent CLI call:
    python main.py --bc channel --lx 6 --ly 3 --lz 2 --nx 384 --ny 192 --nz 128 \
        --ppd 16 --alpha 0.05 --method monte-carlo --gap-z 0.5 --seed 1 \
        --outfile spheres.in
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import droplets  # noqa: E402
from config import Config  # noqa: E402

cfg = Config(
    bc="channel",
    lx=6.0, ly=3.0, lz=2.0,
    nx=384, ny=192, nz=128,
    ppd=16,
    alpha=0.05,
    method="monte-carlo",
    gap_x=0.12, gap_y=0.12,
    gap_z=0.5,          # minimum droplet-surface to wall gap, in units of R
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
