# Warburg-effect simulation: plan

## Context

This grew out of a 15-151/21-241 tutor session. The thread ran: solution sets are a particular solution plus the kernel, then metabolic networks (Sv = 0), then flux balance analysis (FBA), then the Warburg effect. Cameron wants a spatial simulation in which **each cell runs a tiny FBA** and lactate (acid) spreads across a grid. The question it should answer: **do fermenting ("Warburg") tumor cells win by selection when acid harms their neighbors?** It is a passion project: it supports goal 3 in `claude-config/strategy/fall-2026.md` and could become a blog post or web page later.

Decisions already made by Cameron:
- **Language:** Python first, web version later.
- **Pace:** built in the style of the vibe-coding series, one function at a time, with tests.
- **Cell model:** closed-form FBA in the sim, plus an LP version used in tests to check it.
- **Location:** a new top-level `warburg-sim/`.

## The cell model (the linear algebra core)

We use the 6-reaction network from the session. Columns are v_in, v_g, v_e, v_l, v_m, v_out. Rows are glucose, pyruvate, NAD⁺, NADH, lactate.

Sv = 0 forces two things:
- **v_e = v_m = r**, the pyruvate that is respired.
- **v_l = 2·v_g − r**, the pyruvate that is fermented.

Per cell, the bounds are:
- glucose: v_g ≤ G (local glucose)
- O₂: 3r ≤ O (local O₂). This assumes ½ O₂ per v_e and 2.5 O₂ per v_m, which gives 6 O₂ per glucose and matches respiration.
- mitochondrial capacity: r ≤ M, a **heritable trait**.

Maximizing ATP = 2·v_g + 15·r gives the **closed form**:
- v_g = G
- r = min(2G, O/3, M)
- lactate out = 2G − r

When glucose is plentiful but M is small, cells ferment **even with O₂ present**. That's the Warburg effect, produced by a capacity limit rather than hard-coded.

## Files (new folder `warburg-sim/`)

| File | Holds |
|---|---|
| `CLAUDE.md` | Copy of the working style in `02-120/spatial-games/CLAUDE.md` ("one idea at a time", tests, type hints, `main.py`). |
| `fba.py` | `stoichiometric_matrix()`, `mat_vec(S, v)` (plain lists), `solve_cell(G, O, M) -> dict` (closed form), `atp(fluxes)`. |
| `fba_lp.py` | `solve_cell_lp(G, O, M)`: the same problem through `scipy.optimize.linprog`. Used only by tests. |
| `datatypes.py` | `Cell` (kind: EMPTY / HEALTHY / TUMOR, `mito_capacity`, `acid_tolerance`) and `Grid` (cells plus the O₂, glucose and lactate fields as 2D lists). |
| `fields.py` | `diffuse(field, rate)`, which is **the same scheme as the 02-120 Gray–Scott model**. It uses the 3×3 kernel from `Starter Code/gray-scott/main.py` ([[0.05, 0.2, 0.05], [0.2, −1, 0.2], [0.05, 0.2, 0.05]], whose weights sum to 0) times a per-chemical diffusion rate, with edges fixed. Also `add_sources` (blood-vessel squares) and `consume` (applies each cell's fluxes). The FBA output plays the role of Gray–Scott's "reaction" term. |
| `rules.py` | `neighbors(row, col, n_rows, n_cols)` (Moore neighborhood, same idea as `neighbors_fixed` in `02-120/spatial-games/functions.py`), `cell_fate(...)` (divide, survive or die from ATP and local acid), `step(grid, params)`. |
| `render.py` | Grid to SVG plus a stdlib HTML frame viewer. The pattern comes from `02-120/spatial-games/make_test_figures.py`. Also a line chart of tumor % and mean M over time. |
| `main.py` | Runs a scenario and writes `output/`. |
| `test_*.py` | One `unittest` file per module. |
| `requirements.txt` | `scipy` (LP tests only). Everything else is the standard library, matching 02-120. |

Environment: Homebrew Python 3.14 has no numpy or scipy, so we make `warburg-sim/.venv` (git-ignored) and install scipy there. If scipy is missing, the LP tests are skipped, not failed.

## Build order (each step runs and is tested before the next)

1. **Matrix check:** `stoichiometric_matrix` and `mat_vec`. Test: the session's hand-worked flux vectors give Sv = 0, and a wrong vector does not.
2. **Closed-form cell:** `solve_cell`. Tests:
   - no O₂ gives r = 0 and 2G lactate;
   - plenty of everything gives r = 2G, no lactate and 32 ATP per glucose;
   - small M gives fermentation even with O₂;
   - the output always satisfies Sv = 0 and every bound.
3. **LP cross-check:** `solve_cell_lp` agrees with `solve_cell` on a grid of (G, O, M) values. The formula is derived by hand; scipy's general linear-programming solver knows nothing about our derivation. If both give the same optimum on hundreds of inputs, the derivation is right. If they disagree, the test shows exactly which (G, O, M) breaks it.
4. **Diffusion (Gray–Scott style):** `diffuse` on a single spike spreads it symmetrically and conserves the total exactly (closed-box edges: a square only sends to neighbors that exist). `add_sources` and `consume` get their own tests.
5. **Cell rules:** `neighbors` and `cell_fate`, with each threshold tested at, just above and just below it.
6. **One full step** on a tiny 5×5 grid, checked by hand.
7. **Render:** SVG frames and an HTML viewer. Then the first real run, about 100×100 with a vessel row along one edge.
8a. **Null space made explicit (added 2026-10-10):** `linalg.py` (exact rref, rank, null space, left null space, integer basis via Euclid's gcd). `solve_cell` is built on the computed basis (a·FERMENT + b·RESPIRE), plus a per-cell Warburg index a/(a+b).
8. **Evolution (8b):** daughter cells inherit M with a small random mutation, under a fixed seed. Plot mean M over time.
9. **Experiment:** run with acid harm on and off. Does low M (Warburg) take over only when acid hurts healthy cells? This could become a blog post.

Later (v2): immune cells killed by lactate, an evolving acid tolerance, the JS web version (`web/`, the spatial-games pattern), and a blog post.

## Memory and bookkeeping

- Add `claude-config/memory/research/warburg-sim.md` with a row in `research/INDEX.md`, and move the idea line out of `spatial-games.md`.
- Add `warburg-sim/` to the repo map in the root `CLAUDE.md`.
- Commit and push after each completed step.

## Verification

- At each step: `python3 -m unittest -v` from `warburg-sim/`, with the output shown to Cameron.
- `test_fba.py` checks Sv = 0 on every `solve_cell` output, so the linear algebra is checked every time it runs.
- LP agreement (step 3) confirms the closed form is the true FBA optimum.
- Sanity scenarios in step 7: a tumor-free grid stays stable, and a grid with no O₂ produces only lactate.
- Pictures: open the HTML viewer in the browser pane and look at the frames.
