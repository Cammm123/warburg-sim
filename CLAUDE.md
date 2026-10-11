# Warburg-effect simulation — working with Cameron

A grid simulation where every cell runs a tiny flux balance analysis (FBA) and
lactate spreads as acid. Question: do fermenting ("Warburg") tumor cells win by
selection? Full plan and build order: `PLAN.md`. Status: `claude-config/memory/research/warburg-sim.md`.

## How I like to work

Same as `02-120/spatial-games/CLAUDE.md`:

- **One idea at a time.** One function or change per step; I understand it before we move on.
- **Never hand me code I can't verify.** Each piece comes with a plain-English
  explanation, tests that define "correct" (use `test-every-function`), and the
  command to run them with the output shown.
- **Small, runnable steps.** Every step leaves the program runnable.
- **Ask before assuming** when a design choice is mine.

## Code conventions

- Python, plain functions, type hints everywhere, short docstrings, plain lists
  (standard library only, like 02-120). `main.py` is the entry point.
- scipy is used only by the LP cross-check tests, from `.venv/` (git-ignored).
  Setup: `python3 -m venv .venv && .venv/bin/pip install -r requirements.txt`.
- Run everything from inside this folder: `.venv/bin/python -m unittest -v`
  (plain `python3` skips the LP tests). `python3 explore.py` prints and draws what fba.py does.
