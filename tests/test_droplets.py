import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import droplets  # noqa: E402
from config import Config, config_from_cli  # noqa: E402

# small, fast base case: d = 32 * 8/128 = 2  ->  R = 1
BASE = dict(lx=8.0, ly=8.0, lz=6.0, nx=128, ny=128, nz=96, ppd=32, alpha=0.1,
            gap_x=0.12, gap_y=0.12, gap_z=0.12, drop_gap=0.4, do_write=False)


def cfg_of(**kw):
    return Config(**{**BASE, **kw})


def brute_min_gap(pos, cfg):
    """Independent O(N^2) check with explicit periodic images."""
    shifts = [np.array([i, j, k]) * cfg.L
              for i in (-1, 0, 1) for j in (-1, 0, 1)
              for k in ((-1, 0, 1) if cfg.periodic[2] else (0,))]
    best = np.inf
    for a in range(len(pos)):
        for b in range(a + 1, len(pos)):
            for s in shifts:
                best = min(best, np.linalg.norm(pos[a] - pos[b] + s))
    return (best - cfg.d) / cfg.R


def assert_valid(pos, cfg):
    assert pos.shape == (cfg.nb, 3)
    R = cfg.R
    # strictly inside on every axis with required surface gap
    lower = (pos - R) / R
    upper = (cfg.L - pos - R) / R
    assert np.all(lower >= cfg.gaps - 1e-9)
    assert np.all(upper >= cfg.gaps - 1e-9)
    assert np.all(pos - R > 0) and np.all(pos + R < cfg.L)
    if len(pos) > 1:
        assert brute_min_gap(pos, cfg) >= cfg.drop_gap - 1e-9
    assert droplets.check_config(pos, cfg).ok


@pytest.mark.parametrize("method", ["monte-carlo", "boxed"])
@pytest.mark.parametrize("bc", ["triperiodic", "channel"])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_constraints_hold(method, bc, seed):
    cfg = cfg_of(method=method, bc=bc)
    pos = droplets.generate(cfg, np.random.default_rng(seed))
    assert len(pos) == cfg.nb > 5
    assert_valid(pos, cfg)


@pytest.mark.parametrize("method", ["monte-carlo", "boxed"])
def test_non_cubic_domain(method):
    # nbox_x != nbox_y != nbox_z: regression for the old j = rem // nbox[1] bug
    cfg = cfg_of(method=method, lx=12.0, ly=9.0, lz=5.0, nx=192, ny=144, nz=80, alpha=0.12)
    if method == "boxed":
        ncell = [a[0] for a in droplets.boxed_capacity(cfg)]
        assert len(set(ncell)) == 3
    for seed in range(5):
        assert_valid(droplets.generate(cfg, np.random.default_rng(seed)), cfg)


@pytest.mark.parametrize("method", ["monte-carlo", "boxed"])
def test_large_wall_gap_respected(method):
    cfg = cfg_of(method=method, bc="channel", gap_z=1.0, alpha=0.05)
    pos = droplets.generate(cfg, np.random.default_rng(3))
    assert np.all(pos[:, 2] >= 2.0 * cfg.R - 1e-9)
    assert np.all(pos[:, 2] <= cfg.lz - 2.0 * cfg.R + 1e-9)
    assert_valid(pos, cfg)


def test_infeasible_gap_raises():
    with pytest.raises(ValueError):
        cfg_of(gap_z=2.5).validate()   # 2R(1+2.5)=7 > lz=6


@pytest.mark.parametrize("method", ["monte-carlo", "boxed"])
def test_infeasible_alpha_raises_without_hang(method):
    cfg = cfg_of(method=method, alpha=0.45, max_attempts=2000, n_restarts=1)
    with pytest.raises(droplets.GenerationError):
        droplets.generate(cfg, np.random.default_rng(0))


def test_min_image_periodic_and_channel():
    L = np.array([10.0, 10.0, 10.0])
    a, b = np.array([0.5, 5, 5]), np.array([9.5, 5, 5])
    assert droplets.min_image_dist(a, b, L, [True, True, True]) == pytest.approx(1.0)
    a, b = np.array([5, 5, 0.5]), np.array([5, 5, 9.5])
    assert droplets.min_image_dist(a, b, L, [True, True, False]) == pytest.approx(9.0)
    assert droplets.min_image_dist(a, b, L, [True, True, True]) == pytest.approx(1.0)


def test_check_detects_cross_periodic_touching():
    # old code accepted this: two drops touching across the x periodic face
    cfg = cfg_of(gap_x=0.0, gap_y=0.0, gap_z=0.0)
    R = cfg.R
    pos = np.array([[R, 4.0, 3.0], [cfg.lx - R, 4.0, 3.0]])
    rep = droplets.check_config(pos, cfg)
    assert not rep.ok and rep.min_drop_gap == pytest.approx(0.0, abs=1e-9)


def test_check_detects_wall_violation():
    cfg = cfg_of(bc="channel")
    pos = np.array([[4.0, 4.0, cfg.R * 1.05]])  # gap 0.05R < 0.12R
    rep = droplets.check_config(pos, cfg)
    assert not rep.ok and any("z-wall" in m for m in rep.messages)


def test_periodic_boundaries_never_crossed_boxed_tight():
    # zero face gap, tight packing: lattice must still respect cross-face spacing
    cfg = cfg_of(method="boxed", bc="triperiodic", gap_x=0.0, gap_y=0.0, gap_z=0.0,
                 alpha=0.18, drop_gap=0.1)
    for seed in range(5):
        assert_valid(droplets.generate(cfg, np.random.default_rng(seed)), cfg)


def test_diameter_follows_grid():
    assert cfg_of().d == pytest.approx(2.0)
    assert cfg_of(lx=16.0, nx=128).d == pytest.approx(4.0)
    with pytest.warns(UserWarning, match="anisotropic"):
        cfg_of(ny=100).validate()


def test_reproducible_with_seed():
    for m in ("monte-carlo", "boxed"):
        cfg = cfg_of(method=m)
        p1 = droplets.generate(cfg, np.random.default_rng(42))
        p2 = droplets.generate(cfg, np.random.default_rng(42))
        np.testing.assert_array_equal(p1, p2)


def test_roundtrip_file(tmp_path):
    cfg = cfg_of(method="boxed", gap_x=0.0, gap_y=0.0, gap_z=0.0)
    pos = droplets.generate(cfg, np.random.default_rng(0))
    f = tmp_path / "bub.in"
    droplets.write_spheres(f, pos, cfg.R)
    pos2, r2 = droplets.read_spheres(f)
    np.testing.assert_allclose(pos2, pos, atol=1e-8)
    assert droplets.check_config(pos2, cfg, r2).ok


def test_cli_overrides_param():
    cfg = config_from_cli(["--bc", "triperiodic", "--gap-z", "0.3", "--lx", "8",
                           "--nx", "128", "--seed", "7", "--method", "monte-carlo"])
    assert (cfg.bc, cfg.gap_z, cfg.lx, cfg.nx, cfg.seed, cfg.method) == \
        ("triperiodic", 0.3, 8.0, 128, 7, "monte-carlo")
    assert cfg.d == pytest.approx(2.0)
