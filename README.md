# generateDroplets

Generates random, non-overlapping positions for equal-sized spheres (droplets, bubbles or particles) in a rectangular box, for initialising dispersed-phase simulations. The output is a plain text file with one line per sphere: `x y z r`.

Supported domains:
- `triperiodic`: periodic in x, y and z.
- `channel`: periodic in x and y, with walls at `z = 0` and `z = lz`.

Every sphere is guaranteed to:
- stay fully inside the box on all three axes (spheres never cross or touch any face), with a surface-to-face gap of at least `gap_x`, `gap_y`, `gap_z` times R;
- keep a surface-to-surface gap of at least `drop_gap` times R from every other sphere, measured across periodic faces as well.

These constraints are enforced while the spheres are placed, and the written file is checked again afterwards.

## Quick start

```bash
python3 -m venv .venv
.venv/bin/pip install numpy scipy pytest      # scipy is optional (it speeds up the check)

.venv/bin/python main.py                      # uses the defaults in param.py
.venv/bin/python main.py -h                   # lists all options
```

`main.py` is the entry point. It reads the defaults from `param.py`, and any value can be overridden on the command line:

```bash
.venv/bin/python main.py --bc channel --method monte-carlo \
    --lx 6 --ly 3 --lz 2 --nx 384 --ny 192 --nz 128 --ppd 16 \
    --alpha 0.05 --gap-z 0.5 --seed 1            # writes spheres.in
```

It prints a summary and writes the file, then re-reads and checks it. The exit code is `0` if the check passes and non-zero if it fails.

To check an existing file:

```bash
.venv/bin/python check_walls.py spheres.in --bc channel --gap-z 0.5   # plus the same domain flags
```

## Parameters (`param.py`)

| Name | Meaning |
|---|---|
| `bc` | `channel` or `triperiodic` |
| `lx, ly, lz` / `nx, ny, nz` | Domain lengths and grid points |
| `ppd` | Grid points per diameter: `d = ppd * lx / nx`. You get a warning if the grid is anisotropic. |
| `alpha` | Target volume fraction. The number of spheres is `floor(alpha V / V_sphere)`. |
| `method` | `monte-carlo` (random sequential addition) or `boxed` (randomly chosen sites on a jittered lattice) |
| `epsf` | `boxed` only: jitter amplitude as a fraction of the free space, in [0, 1] |
| `gap_x, gap_y, gap_z` | Minimum gap between the sphere surface and each face, in units of R. `gap_z` is the wall offset. |
| `drop_gap` | Minimum gap between sphere surfaces, in units of R. The default of 0.4 gives a centre spacing of 1.2 d. |
| `max_attempts, n_restarts` | Monte-carlo limits. The run fails with a message instead of hanging. |
| `seed` | Random seed, for reproducible output (`None` gives a different result each run) |
| `outfile, do_write, do_plot` | Output file name (default `spheres.in`), whether to write it, and an optional 3D scatter plot (needs matplotlib) |

On the command line, underscores become dashes, e.g. `--gap-z`, `--drop-gap`, `--max-attempts`.

**Choosing a method:** `boxed` is fast and reliable up to its lattice limit, which it prints as "max alpha". `monte-carlo` gives more disordered configurations, but it jams at lower `alpha`. In the default 16 x 16 x 10 box with d = 2 and the default gaps, monte-carlo succeeds up to about 17%, while the boxed lattice reaches about 23.6%. Smaller gaps or a larger domain relative to d raise both limits. If monte-carlo jams, it reports how many spheres it managed to place.

## Examples

| Script | Case |
|---|---|
| `examples/channel_6x3x2.py` | 6 x 3 x 2 channel, d = 0.25, alpha = 5%, wall gap 0.5 R, monte-carlo |
| `examples/triperiodic_2pi_cube.py` | (2 pi)^3 tri-periodic cube, d = pi/4, alpha = 10%, boxed |

Each script writes `spheres.in` to the current directory, or to the path you pass as the first argument. Each one also shows how to use `Config` and `droplets.generate` from Python, and lists the equivalent `main.py` command in its docstring.

### Using an example to set up a simulation

Generate the file in the case directory and check it (`G` is the path to this directory):

```bash
G=/path/to/generateDroplets
mkdir -p run_channel && cd run_channel
$G/.venv/bin/python $G/examples/channel_6x3x2.py     # writes ./spheres.in
$G/.venv/bin/python $G/check_walls.py spheres.in --bc channel \
    --lx 6 --ly 3 --lz 2 --nx 384 --ny 192 --nz 128 --ppd 16 --alpha 0.05 --gap-z 0.5
```

Then rename or point your solver at the file if it expects a different name or format. Keep these settings consistent with the simulation:
- **Domain and grid:** lengths `lx, ly, lz` and grid points `nx, ny, nz` must match the simulation, otherwise the resolution per diameter (`ppd`) will be wrong. For the examples: 6 x 3 x 2 on 384 x 192 x 128, and (2 pi)^3 on 256^3.
- **Boundary conditions:** `bc` must match the simulation: walls in z for `channel`, periodic everywhere for `triperiodic`.
- **Grid stretching:** the generator assumes a uniform grid. With a stretched z grid, gaps are still enforced in physical length units, but `ppd` refers only to dx.

## Layout

```
main.py         entry point (command line)
param.py        default parameters
config.py       Config dataclass: derived quantities, input checks, command-line parser
droplets.py     generators, check_config, file read/write
check_walls.py  checks an existing sphere file
examples/       example scripts
tests/          pytest suite
```

## Tests

```bash
.venv/bin/python -m pytest -q tests/
```

## Notes

- Running the generator again overwrites the output file.
- Output files (`*.in`) are git-ignored.
