# warburg-sim

Do fermenting ("Warburg") tumor cells win by natural selection? This is a grid
simulation where every cell runs a small flux balance analysis (FBA) to choose how to
make ATP, and the lactate it makes spreads as acid that can kill its neighbors.

**Try it in the browser:** https://warburg-sim.vercel.app/

## What's here

| Path | What it is |
|---|---|
| `main.py` | Entry point: runs a simulation and writes figures to `output/`. |
| `fba.py`, `fba_lp.py`, `linalg.py` | Per-cell metabolism: closed-form FBA, with an LP cross-check (scipy). |
| `fields.py`, `rules.py`, `scenarios.py` | Diffusion of glucose, oxygen and acid; division, death and mutation rules; preset scenarios. |
| `experiments.py` | Parameter sweeps over seeds. Results in `output/experiments/results.md`. |
| `FINDINGS.md` | What the model says so far. |
| `PLAN.md` | The build plan, step by step. |
| `web/` | The plain-JS browser version, checked against the Python model with `golden.json`. |

## Run it

The model is standard-library Python. scipy is only needed for the LP cross-check tests.

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python -m unittest -v
python3 main.py
```

Web version: `cd web && node model.test.js`, then serve `web/` with any static server.
