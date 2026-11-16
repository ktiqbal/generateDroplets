# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Coding guidelines

- Don't write comments unless they're absolutely necessary to explain genuinely non-obvious logic (in which case, keep them concise and clear), or a TODO/FIXME with a concrete reason. Don't write comments that just restate what the code does — that's AI slop, and it's a DRY violation: it's now a second place to update when the logic changes.
- Don't assume existing code is correct or well-structured — call out bad patterns and propose better ones rather than copying them forward.
- Prefer readable, intuitive code over clever abstractions. No speculative flexibility/configurability that wasn't asked for.
- When adding a feature or fixing a bug, add tests that actually fail if the behavior regresses — not tests that only exercise structure or pad coverage.
- Keep commit messages concise: a brief one-line summary, then one paragraph (max) on the additions.
- Never put a Claude Code session id/link (e.g. a `Claude-Session:` trailer) in a commit message or PR description — strip it if the harness auto-appends one. A `Co-Authored-By: Claude ...` trailer is fine to keep.

## Think before coding

Don't assume, don't hide confusion, surface tradeoffs. Before implementing:

- State your assumptions explicitly. If uncertain, ask.
- If multiple reasonable interpretations exist, present them — don't pick one silently.
- If a simpler approach exists than the one implied by the request, say so.
- If something is unclear, stop and name what's confusing rather than guessing.

## What this is

Generates random, non-overlapping equal-sized sphere positions (droplets/bubbles) in a box, for initialising dispersed-phase simulations. Output is a text file, one line per sphere: `x y z r`. See `README.md` for the full parameter table and usage.

## Commands

Use the project venv (`.venv`, git-ignored; create with `python3 -m venv .venv` and install `numpy scipy pytest`; matplotlib only for `--do-plot`).

```bash
.venv/bin/python main.py -h                       # all options (flags = param names with dashes)
.venv/bin/python main.py --bc channel --method boxed --alpha 0.1 --seed 1
.venv/bin/python check_walls.py spheres.in --bc channel --gap-z 0.5   # re-check an existing file (pass same domain flags)
.venv/bin/python -m pytest -q tests/              # full suite
.venv/bin/python -m pytest -q tests/test_droplets.py -k <name>   # single test
.venv/bin/python examples/channel_6x3x2.py [outpath]
```

`main.py` exit codes: 0 ok, 1 generation failed (`GenerationError`), 2 post-write check failed. No linter is configured.

## Architecture

- `param.py` holds module-level defaults. `Config.from_param()` (in `config.py`) loads them via the dataclass fields; `config_from_cli` then layers argparse overrides on top. The parser is auto-generated from `dataclasses.fields(Config)`, so adding a parameter means: add the field to `Config`, add a default in `param.py`, and (if not float) register its type in `_TYPES` in `config.py`. Check it in `Config.validate()`.
- `Config` also carries all derived geometry: `d = ppd * lx/nx` (diameter derives from `dx` only), `nb` from `alpha`, admissible centre bounds `lo`/`hi` (face gaps), `min_spacing = d + drop_gap*R`, and `periodic` (x, y always periodic; z periodic only for `triperiodic`).
- `droplets.py` has the two generators (`generate_monte_carlo` = random sequential addition with restarts; `generate_boxed` = random sites on a jittered lattice with a computed capacity), dispatched by `generate(cfg, rng)`. Distances use minimum-image convention per periodic axis. `check_config` independently verifies bounds and pair spacing and returns a `Report`; `main.py` runs it on the re-read output file, so write/read round-trip is validated.
- `legacy/` is a git-ignored local reference copy of the original code; do not edit or depend on it.
- Output `*.in` files are git-ignored.
