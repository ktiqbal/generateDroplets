"""Configuration for droplet generation: defaults from param.py, overridable from the CLI."""

import argparse
import dataclasses
import math
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

BCS = ("triperiodic", "channel")
METHODS = ("monte-carlo", "boxed")


@dataclass
class Config:
    bc: str = "channel"
    lx: float = 16.0
    ly: float = 16.0
    lz: float = 10.0
    nx: int = 256
    ny: int = 256
    nz: int = 160
    ppd: float = 32
    alpha: float = 0.1
    method: str = "boxed"
    epsf: float = 0.5
    gap_x: float = 0.12
    gap_y: float = 0.12
    gap_z: float = 0.12
    drop_gap: float = 0.4
    max_attempts: int = 100000
    n_restarts: int = 20
    seed: Optional[int] = None
    outfile: str = "bub.in"
    do_write: bool = True
    do_plot: bool = False

    # ---- construction -------------------------------------------------
    @classmethod
    def from_param(cls, module=None):
        """Build a Config from the module-level values in param.py."""
        if module is None:
            import param as module
        kw = {f.name: getattr(module, f.name) for f in dataclasses.fields(cls)
              if hasattr(module, f.name)}
        return cls(**kw)

    def replace(self, **kw):
        return dataclasses.replace(self, **kw)

    # ---- derived quantities -------------------------------------------
    @property
    def L(self):
        return np.array([self.lx, self.ly, self.lz], dtype=float)

    @property
    def dxyz(self):
        return self.L / np.array([self.nx, self.ny, self.nz], dtype=float)

    @property
    def d(self):
        """Droplet diameter = ppd * dx."""
        return self.ppd * self.dxyz[0]

    @property
    def R(self):
        return 0.5 * self.d

    @property
    def volume(self):
        return float(np.prod(self.L))

    @property
    def drop_volume(self):
        return 4.0 / 3.0 * math.pi * self.R ** 3

    @property
    def nb(self):
        """Number of droplets for the requested volume fraction."""
        return int(math.floor(self.alpha * self.volume / self.drop_volume + 1e-12))

    @property
    def actual_alpha(self):
        return self.nb * self.drop_volume / self.volume

    @property
    def periodic(self):
        return np.array([True, True, self.bc == "triperiodic"])

    @property
    def gaps(self):
        return np.array([self.gap_x, self.gap_y, self.gap_z], dtype=float)

    @property
    def lo(self):
        """Lowest admissible centre coordinate per axis."""
        return self.R * (1.0 + self.gaps)

    @property
    def hi(self):
        """Highest admissible centre coordinate per axis."""
        return self.L - self.R * (1.0 + self.gaps)

    @property
    def min_spacing(self):
        """Minimum allowed centre-to-centre distance."""
        return self.d + self.drop_gap * self.R

    # ---- checks -------------------------------------------------------
    def validate(self):
        if self.bc not in BCS:
            raise ValueError(f"bc must be one of {BCS}, got {self.bc!r}")
        if self.method not in METHODS:
            raise ValueError(f"method must be one of {METHODS}, got {self.method!r}")
        if np.any(self.L <= 0):
            raise ValueError("lx, ly, lz must be positive")
        if min(self.nx, self.ny, self.nz) <= 0 or self.ppd <= 0:
            raise ValueError("nx, ny, nz and ppd must be positive")
        if not 0.0 <= self.alpha < 1.0:
            raise ValueError("alpha must be in [0, 1)")
        if not 0.0 <= self.epsf <= 1.0:
            raise ValueError("epsf must be in [0, 1]")
        if np.any(self.gaps < 0) or self.drop_gap < 0:
            raise ValueError("gap_x, gap_y, gap_z and drop_gap must be >= 0")
        if self.max_attempts < 1 or self.n_restarts < 0:
            raise ValueError("max_attempts must be >= 1 and n_restarts >= 0")
        dx = self.dxyz
        if np.max(np.abs(dx - dx[0])) > 1e-12 * dx[0]:
            warnings.warn(f"anisotropic grid (dx, dy, dz) = {tuple(dx)}; "
                          "droplet diameter is computed from dx only", stacklevel=2)
        if np.any(self.hi < self.lo):
            bad = "xyz"[int(np.argmax(self.lo - self.hi))]
            raise ValueError(f"droplet (d={self.d:g}) plus face gaps do not fit in the "
                             f"{bad}-direction of the domain")
        if self.nb == 0:
            warnings.warn("alpha is too small for a single droplet; nb = 0", stacklevel=2)
        return self

    def summary(self):
        return (f"bc={self.bc} method={self.method} domain=({self.lx:g},{self.ly:g},{self.lz:g}) "
                f"grid=({self.nx},{self.ny},{self.nz}) ppd={self.ppd:g}\n"
                f"d={self.d:g} R={self.R:g} nb={self.nb} alpha(target)={self.alpha:g} "
                f"alpha(actual)={self.actual_alpha:.5f}\n"
                f"face gaps/R=({self.gap_x:g},{self.gap_y:g},{self.gap_z:g}) "
                f"drop_gap/R={self.drop_gap:g}")


def _bool(s):
    if isinstance(s, bool):
        return s
    if s.lower() in ("1", "true", "yes", "y", "on"):
        return True
    if s.lower() in ("0", "false", "no", "n", "off"):
        return False
    raise argparse.ArgumentTypeError(f"not a boolean: {s!r}")


def _seed(s):
    return None if s.lower() == "none" else int(s)


_TYPES = {"bc": str, "method": str, "outfile": str, "seed": _seed,
          "do_write": _bool, "do_plot": _bool,
          "nx": int, "ny": int, "nz": int, "max_attempts": int, "n_restarts": int}

_HELP = {
    "bc": "triperiodic | channel (walls at z=0, z=lz)",
    "method": "monte-carlo | boxed",
    "ppd": "grid points per droplet diameter (d = ppd*lx/nx)",
    "alpha": "target volume fraction",
    "epsf": "boxed jitter fraction in [0,1]",
    "gap_z": "min droplet-surface to z-face gap, in units of R",
    "gap_x": "min droplet-surface to x-face gap, in units of R",
    "gap_y": "min droplet-surface to y-face gap, in units of R",
    "drop_gap": "min droplet-droplet surface gap, in units of R",
    "seed": "RNG seed (int or none)",
}


def build_parser(defaults: Config):
    p = argparse.ArgumentParser(description="Generate non-overlapping droplet positions.")
    for f in dataclasses.fields(Config):
        typ = _TYPES.get(f.name, float)
        kw = dict(type=typ, default=getattr(defaults, f.name), help=_HELP.get(f.name))
        if f.name == "bc":
            kw["choices"] = BCS
        if f.name == "method":
            kw["choices"] = METHODS
        p.add_argument("--" + f.name.replace("_", "-"), dest=f.name, **kw)
    return p


def config_from_cli(argv=None, defaults: Optional[Config] = None):
    """param.py defaults, overridden by command-line flags."""
    defaults = defaults if defaults is not None else Config.from_param()
    args = build_parser(defaults).parse_args(argv)
    return Config(**vars(args))
