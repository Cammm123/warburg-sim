"""
Calibration: with no tumor at all, does healthy tissue stay alive?

Fills a grid with healthy cells (capillaries every `spacing` squares), runs it,
and prints how many healthy cells are left plus the lowest O2 and glucose.
- uptake_fraction (physics) is the knob we calibrate.
- diffusion_substeps (numerics) is "converged" once raising it stops changing
  the answer.
Run from this folder:  python3 calibrate.py
Results are recorded in scenarios.calibrated_params().
"""
import random

from datatypes import HEALTHY, Params
from rules import step
from scenarios import healthy_tissue


def run(size: int, spacing: int, substeps: int, fraction: float, uptake: float,
        steps: int) -> tuple[int, int, float, float]:
    """Returns (healthy cells at start, at the end, lowest O2, lowest glucose at the end)."""
    params = Params(diffusion_substeps=substeps, uptake_fraction=fraction)
    grid = healthy_tissue(size, size, spacing, params)
    start = grid.count(HEALTHY)
    rng = random.Random(0)
    for _ in range(steps):
        grid = step(grid, params, rng)
    low_o2 = min(min(row) for row in grid.oxygen)
    low_glucose = min(min(row) for row in grid.glucose)
    return start, grid.count(HEALTHY), low_o2, low_glucose


def main() -> None:
    size, steps, spacing, uptake = 30, 60, 6, 1.0
    print(f"{size}x{size} healthy tissue, capillaries every {spacing}, {steps} steps")
    for fraction in (0.2, 0.1, 0.05):
        print(f"\nuptake_fraction {fraction}:")
        for substeps in (5, 10, 20, 40):
            start, end, low_o2, low_glc = run(size, spacing, substeps, fraction, uptake, steps)
            print(f"  substeps {substeps:>3}: healthy {start} -> {end:>4} ({100 * end / start:5.1f}%),"
                  f" lowest O2 {low_o2:5.2f}, lowest glucose {low_glc:5.3f}")


if __name__ == "__main__":
    main()
