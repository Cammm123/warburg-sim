# Findings so far (2026-10-10)

The numbers are in `output/experiments/results.md`, made by `python3 experiments.py`
(40x40 grid, 200 steps, 3 seeds per setting, reported as mean ± sd). "Marker" is a neutral
trait that mutates like the real one but does nothing. A trait that moves much more than
the marker is being selected; one that doesn't could just be chance.

## What the model says

1. **The tumor can only invade if acid kills healthy cells.** With acid harm off, the
   tumor stays at its seed size (8 cells) and all healthy tissue survives. This follows
   from the rules: cells divide only into empty squares, and calibrated healthy tissue
   never dies on its own, so acid is the only way space opens up. It is a consequence of
   our assumptions (no pushing or crowding out), not a discovery.
2. **Healthy acid tolerance decides the outcome.** At tolerance 1 the tumor kills every
   healthy cell and takes over, at every blood-glucose level tested (2, 3, 5). At
   tolerance 2 it is mostly contained (6–20 cells). It only partly invades at glucose 5
   (138 ± 102 cells), where the outcome varies a lot from seed to seed.
3. **Evolution never favored more fermentation.** In the budget model the mitochondria
   share rose 0.26–0.41 every time the tumor grew, while the marker moved at most 0.1
   either way. That is selection toward respiration. The Warburg index fell from about
   0.9 to about 0.3–0.35.
4. **The tumor ferments to invade, then evolves away from fermenting.** The early, highly
   fermenting tumor makes the acid that clears out healthy tissue. Once it has space, its
   descendants shift to respiring.
5. **Free mitochondria (no budget): no clear evolution.** M changed by 0.27 ± 0.09 and the
   marker by 0.19 ± 0.30, so the change can't be told apart from chance. This corrects
   an earlier one-seed reading of run 2 as "selection against Warburg".
6. **Cells that chase the most ATP right now poison themselves.** With a plastic budget,
   cells pick pure fermentation (40 ATP vs 35) and die of their own acid in one step.

## Why this doesn't match real tumors (and what to try next)

Real tumors stay highly fermenting. The model points to the missing ingredients:

- **Acid tolerance that evolves.** Real tumors evolve acid resistance. Here it is fixed at
  3, so fermenting is always partly self-harming.
- **Carbon for building new cells.** Dividing cells need glucose carbon as building
  material, not just ATP (the Vander Heiden "biosynthesis" explanation). Our objective
  counts ATP only.
- **Pushing and crowding out.** Real tumors also grow by displacing healthy cells. Adding
  that would change finding 1.
- **Uneven blood supply.** Real tumors outgrow their vessels, leaving low-O2 cores.

## Caveats

Budget B = 20, mitochondria cost k = 10 and both acid tolerances are placeholder
numbers. There are only 3 seeds per setting, on a 40x40 grid for 200 steps. The ATP
yields (2 and 32) are standard textbook values.
