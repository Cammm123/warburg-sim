# Warburg simulation: web version (plan)

## Context

The Python `warburg-sim/` (steps 1–9, 178 tests, findings in `warburg-sim/FINDINGS.md`) is
done. Cameron wants an interactive web version, and then to make it the November blog
post (due 2026-11-30, written with `/blog`). Decisions made by Cameron:
- **Hosting:** a separate Vercel site, like spatial-game.vercel.app. Cameron deploys it;
  the blog links to it.
- **Features:** all four:
  - scenario presets
  - live sliders
  - click a cell to see its FBA / null-space point
  - live charts
- **Build mode:** Claude builds each step with node tests and shows it running; Cameron
  reviews and steers.

## Approach

Follow the spatial-games web pattern (`02-120/spatial-games/web/`): plain JS with no
build step, a `rules.js`-style model module that works in both node and the browser,
`*.test.js` files run with `node` using the same `check` / `checkClose` helpers as
`web/rules.test.js`, and `index.html` + `app.js`. The JS is a faithful port of the Python
model, so the page shows exactly what the experiments measured.

**Cross-language check (the key verification):** `warburg-sim/web/make_golden.py` runs
the Python functions on fixed inputs and writes `golden.json`. `model.test.js` loads it
and requires the JS to match to 1e-9. Inputs covered:
- `solve_cell`, `best_budget_amounts`, `pathway_amounts`, the null-space basis;
- `diffuse` on spikes, `refill`, `metabolism` on small grids;
- one full `step` on grids with no divisions, so no randomness is involved.

Randomness can't match Python's generator, so dice and mutation are checked
statistically, as in Python (tumor wins about 75%, mutation sd 0.05), using a seeded
JS RNG (mulberry32 + Box–Muller).

## Files (new folder `warburg-sim/web/`)

| File | Holds |
|---|---|
| `SPEC.md` | The build spec (pattern: `spatial-games/spec.txt`). |
| `model.js` | Port of `linalg` (rref, null space, gcd → integer basis), `fba` (S, pathways computed from the null space, `pathwayAmounts`, `solveCell`, `budgetTraits`, `bestBudgetAmounts` corner method), `fields` (`diffuse`, refill, applyChange), `rules` (`neighbors`, `cellFate`, `cellFluxes`, `chooseDivisions` with weighted dice, `makeDaughter` with marker, `step` with substeps), `metrics`, `scenarios` (calibrated params, `healthyTissue`, presets). |
| `model.test.js` | Hand checks ported from the Python tests, plus golden comparisons and statistical checks. |
| `make_golden.py`, `golden.json` | Python → JSON reference cases. |
| `index.html`, `style.css`, `app.js` | The page. |

## Page design

- **Simulation card:** three canvases (cells, O₂, lactate), Play/Pause/Step/Reset, a
  speed control, and step / healthy / tumor counts.
- **Presets** (each sets params and the starting grid, as in `main.py` / `experiments.py`):
  - Warburg invasion (run 1)
  - Free mitochondria evolving (run 2)
  - Budget evolution (run 3)
  - Plastic cells poison themselves (run 4)
  - Acid off (no invasion)
- **Sliders** (take effect on the next step):
  - healthy acid tolerance (0.5–4)
  - blood glucose (2–10)
  - mutation on/off
  - tumor dice weight (1–5)
- **Live charts:** healthy vs tumor counts; mean mitochondria share (or M) vs mean
  Warburg index. Canvas line charts that update while playing.
- **Click a cell:** a panel showing that cell's local glucose and O₂, its a (fermented)
  and b (respired), ATP, lactate and Warburg index. A small plot shows the cell as a
  point in the (a, b) null-space plane, with the constraint lines.
- **Explainer** (right column): the model in plain words, Sv = 0 and the two pathways,
  the budget trade-off, and what the presets show (from FINDINGS.md, with caveats).
- Mobile layout and dark-mode-safe colors, as in spatial-games.
- Default 60×60 grid with 10 substeps; several steps per animation frame if fast enough.

## Build steps (each runs and is tested before the next)

1. Golden file + model core (linalg, fba, budget): node tests, including golden matches.
2. Fields, rules and `step` with the seeded RNG: hand, golden and statistical tests.
3. Page skeleton: canvases, Play/Step/Reset, run-1 preset. Check it in the browser pane.
4. Presets and sliders.
5. Live charts.
6. Click-a-cell inspector with the null-space plot.
7. Explainer text, mobile layout and polish. Write the SPEC status, memory note and repo map.
8. Cameron deploys on Vercel (root directory `warburg-sim/web`); then `/blog` for the
   November post.

## Verification

- `cd warburg-sim/web && node model.test.js`: every check passes, including the
  Python-golden ones.
- The Python suite stays green: `.venv/bin/python -m unittest`.
- In the browser pane (a new launch config serving `warburg-sim/web`):
  - each preset runs and looks like its Python run (run 1: healthy cleared by acid,
    tumor around capillaries; acid off: no growth);
  - sliders change behavior;
  - clicking a cell shows numbers consistent with `solveCell`;
  - layout works at phone width.
- Commit and push after each step.
