"""
Random, non-overlapping droplet positions for tri-periodic or channel domains.

Original: April 2024, ktiqbal (modified from code by marcres). Refactored Nov 2026.

All constraints are enforced while sampling, not by rejecting whole sets afterwards:
  * every droplet lies strictly inside the box on every axis, with surface-to-face
    clearance >= gap_i * R  (i = x, y, z)
  * every pair of droplets has surface-to-surface clearance >= drop_gap * R, measured
    with the minimum-image convention along periodic axes (x, y always; z if triperiodic)
"""

from dataclasses import dataclass

import numpy as np

from config import Config

# relative tolerance (in units of R) used when validating; covers output rounding
CHECK_TOL = 1e-6


class GenerationError(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------
def min_image_delta(p, P, L, periodic):
    """Separation vectors P - p using the minimum-image convention on periodic axes."""
    delta = np.asarray(P, dtype=float) - np.asarray(p, dtype=float)
    L = np.asarray(L, dtype=float)
    per = np.asarray(periodic, dtype=bool)
    delta[..., per] -= L[per] * np.round(delta[..., per] / L[per])
    return delta


def min_image_dist(p, P, L, periodic):
    """Minimum-image distance(s) between point p and point(s) P."""
    return np.linalg.norm(min_image_delta(p, P, L, periodic), axis=-1)


# ---------------------------------------------------------------------------
# generators
# ---------------------------------------------------------------------------
def generate_monte_carlo(cfg: Config, rng: np.random.Generator, batch: int = 256):
    """Random sequential addition with per-droplet attempt cap and full restarts."""
    nb, lo, hi, L, per = cfg.nb, cfg.lo, cfg.hi, cfg.L, cfg.periodic
    smin = cfg.min_spacing
    best = 0
    for restart in range(cfg.n_restarts + 1):
        pos = np.empty((nb, 3))
        n = 0
        while n < nb:
            placed = False
            tries = 0
            while tries < cfg.max_attempts:
                m = min(batch, cfg.max_attempts - tries)
                cand = lo + (hi - lo) * rng.random((m, 3))
                tries += m
                if n == 0:
                    ok = np.ones(m, dtype=bool)
                else:
                    dist = np.linalg.norm(
                        min_image_delta(cand[:, None, :], pos[None, :n, :], L, per), axis=-1)
                    ok = np.all(dist >= smin, axis=1)
                idx = np.flatnonzero(ok)
                if idx.size:
                    pos[n] = cand[idx[0]]
                    n += 1
                    placed = True
                    break
            if not placed:
                break
        if n == nb:
            if restart:
                print(f"monte-carlo: succeeded after {restart} restart(s)")
            return pos
        best = max(best, n)
    raise GenerationError(
        f"monte-carlo jammed: placed at most {best}/{nb} droplets after "
        f"{cfg.n_restarts} restarts (max_attempts={cfg.max_attempts}). "
        "Lower alpha/drop_gap/gaps, or use method='boxed'.")


def _lattice_axis(W, L, s, periodic):
    """Base lattice along one axis inside an admissible interval of width W.

    Returns (n, pitch, jitter_amplitude, offset) such that neighbours (and, on a
    periodic axis, the pair across the periodic face) stay >= s apart after jitter.
    """
    Weff = min(W, L - s) if periodic else W
    if Weff < 0:
        n = 1
    else:
        n = int(np.floor(Weff / s + 1e-12)) + 1
    if n == 1:
        return 1, 0.0, 0.5 * W, 0.5 * W
    p = Weff / (n - 1)
    slack = p - s
    if periodic:
        slack = min(slack, (L - (n - 1) * p) - s)
    amp = 0.5 * max(slack, 0.0)
    offset = 0.5 * (W - (n - 1) * p)
    return n, p, amp, offset


def boxed_capacity(cfg: Config):
    W = cfg.hi - cfg.lo
    s = cfg.min_spacing
    return [_lattice_axis(W[a], cfg.L[a], s, cfg.periodic[a]) for a in range(3)]


def generate_boxed(cfg: Config, rng: np.random.Generator):
    """Jittered lattice: pick nb distinct lattice sites at random, jitter each uniformly."""
    nb = cfg.nb
    axes = boxed_capacity(cfg)
    ncell = [a[0] for a in axes]
    ntot = int(np.prod(ncell))
    max_alpha = ntot * cfg.drop_volume / cfg.volume
    print(f"boxed: lattice {ncell[0]} x {ncell[1]} x {ncell[2]} = {ntot} sites "
          f"(max alpha with these gaps = {max_alpha * 100:.2f}%)")
    if nb > ntot:
        raise GenerationError(
            f"boxed: need {nb} droplets but only {ntot} lattice sites fit "
            f"(max alpha = {max_alpha:.4f}). Lower alpha/drop_gap/gaps or use monte-carlo.")
    sites = rng.choice(ntot, size=nb, replace=False)
    ijk = np.stack(np.unravel_index(sites, ncell), axis=1)  # correct for any nx != ny != nz
    pos = np.empty((nb, 3))
    for a, (n, p, amp, off) in enumerate(axes):
        base = cfg.lo[a] + off + ijk[:, a] * p
        jitter = cfg.epsf * amp * rng.uniform(-1.0, 1.0, nb)
        pos[:, a] = np.clip(base + jitter, cfg.lo[a], cfg.hi[a])
    return pos


def generate(cfg: Config, rng=None):
    cfg.validate()
    if rng is None or isinstance(rng, (int, np.integer)):
        rng = np.random.default_rng(cfg.seed if rng is None else rng)
    if cfg.nb == 0:
        return np.empty((0, 3))
    if cfg.method == "monte-carlo":
        pos = generate_monte_carlo(cfg, rng)
    else:
        pos = generate_boxed(cfg, rng)
    rep = check_config(pos, cfg)
    if not rep.ok:  # should be impossible by construction
        raise GenerationError("internal error, generated set fails validation:\n" + str(rep))
    return pos


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------
@dataclass
class Report:
    ok: bool
    n: int
    min_face_gap: np.ndarray      # per axis, surface-to-face, units of R
    min_drop_gap: float           # surface-to-surface, units of R (inf if n < 2)
    messages: list

    def __str__(self):
        fg = ", ".join(f"{c}={g:.4f}" for c, g in zip("xyz", self.min_face_gap))
        lines = [f"check {'PASSED' if self.ok else 'FAILED'}: n={self.n}",
                 f"  min face gap / R : {fg}",
                 f"  min drop gap / R : {self.min_drop_gap:.4f}"]
        return "\n".join(lines + ["  " + m for m in self.messages])


def _nearest_dist(pos, L, per):
    """Smallest minimum-image centre distance among all pairs."""
    n = len(pos)
    try:
        from scipy.spatial import cKDTree
        # non-periodic axes get a box so large that wrapping never matters
        box = np.where(per, L, 3.0 * L)
        shifted = np.where(per, np.mod(pos, L), pos + L)
        shifted = np.where(per & (shifted >= L), 0.0, shifted)
        tree = cKDTree(shifted, boxsize=box)
        dd, _ = tree.query(shifted, k=2)
        return float(dd[:, 1].min())
    except ImportError:
        best = np.inf
        for i in range(n - 1):
            best = min(best, float(min_image_dist(pos[i], pos[i + 1:], L, per).min()))
        return best


def check_config(pos, cfg: Config, radii=None):
    """Validate face clearance on all axes and pairwise min-image clearance."""
    pos = np.asarray(pos, dtype=float).reshape(-1, 3)
    n = len(pos)
    R = cfg.R
    r = np.full(n, R) if radii is None else np.asarray(radii, dtype=float)
    msgs = []
    tol = CHECK_TOL
    if n == 0:
        return Report(True, 0, np.full(3, np.inf), np.inf, ["no droplets"])
    if radii is not None and np.max(np.abs(r - R)) > tol * R:
        msgs.append(f"radius in file ({r.min():g}..{r.max():g}) differs from config R={R:g}")
    lower = (pos - r[:, None]) / R
    upper = (cfg.L - pos - r[:, None]) / R
    face = np.minimum(lower, upper).min(axis=0)
    ok = True
    for a, c in enumerate("xyz"):
        if face[a] < cfg.gaps[a] - tol:
            ok = False
            kind = "wall" if (c == "z" and cfg.bc == "channel") else "face"
            msgs.append(f"{c}-{kind} clearance {face[a]:.6f} R < required {cfg.gaps[a]:g} R")
    if n >= 2:
        dmin = _nearest_dist(pos, cfg.L, cfg.periodic)
        dgap = (dmin - 2.0 * r.max()) / R
        if dgap < cfg.drop_gap - tol:
            ok = False
            msgs.append(f"droplet clearance {dgap:.6f} R < required {cfg.drop_gap:g} R"
                        + (" (overlap!)" if dgap < 0 else ""))
    else:
        dgap = np.inf
    return Report(ok, n, face, dgap, msgs)


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def write_spheres(path, pos, R):
    """One line per droplet: x y z r  (format read by CaNS-Fizzy read_sphere_file)."""
    with open(path, "w") as f:
        for x, y, z in pos:
            f.write(f"{x:.8f} {y:.8f} {z:.8f} {R:.8f}\n")


def read_spheres(path):
    data = np.loadtxt(path, ndmin=2)
    if data.size == 0:
        return np.empty((0, 3)), np.empty(0)
    return data[:, :3], data[:, 3]
