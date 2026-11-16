"""Generate droplet positions.  Defaults come from param.py; any value can be overridden:

    python main.py --bc channel --method monte-carlo --gap-z 0.3 --alpha 0.05 --seed 1
"""

import sys

import numpy as np

import droplets
from config import config_from_cli


def plot(pos, cfg):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.scatter(*pos.T, s=20)
    ax.set_xlim(0, cfg.lx); ax.set_ylim(0, cfg.ly); ax.set_zlim(0, cfg.lz)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    plt.show()


def main(argv=None):
    cfg = config_from_cli(argv).validate()
    print(cfg.summary())
    rng = np.random.default_rng(cfg.seed)
    try:
        pos = droplets.generate(cfg, rng)
    except droplets.GenerationError as e:
        print("ERROR:", e, file=sys.stderr)
        return 1

    if cfg.do_write:
        droplets.write_spheres(cfg.outfile, pos, cfg.R)
        print(f"written {len(pos)} droplets to {cfg.outfile}")
        pos_chk, radii = droplets.read_spheres(cfg.outfile)  # validate what was written
    else:
        pos_chk, radii = pos, None
    rep = droplets.check_config(pos_chk, cfg, radii)
    print(rep)
    if cfg.do_plot:
        plot(pos, cfg)
    return 0 if rep.ok else 2


if __name__ == "__main__":
    sys.exit(main())
