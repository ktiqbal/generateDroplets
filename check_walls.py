"""Validate an existing sphere file against the configured domain and clearances.

    python check_walls.py spheres.in [--bc channel --lz 10 --gap-z 0.12 ...]
"""

import sys

import droplets
from config import Config, build_parser


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0].startswith("-"):
        print(__doc__)
        return 1
    fname, rest = argv[0], argv[1:]
    cfg = Config(**vars(build_parser(Config.from_param()).parse_args(rest)))
    pos, radii = droplets.read_spheres(fname)
    rep = droplets.check_config(pos, cfg, radii)
    print(rep)
    return 0 if rep.ok else 2


if __name__ == "__main__":
    sys.exit(main())
