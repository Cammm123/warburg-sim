"""
Run the Warburg simulation and save pictures.

Run from this folder:  python3 main.py
Writes output/run1/: viewer.html (open in a browser; Play or drag the slider),
snapshots at a few steps, and counts.svg (healthy and tumor cells over time).
"""
import os
import random
import time

from datatypes import HEALTHY, TUMOR, Cell, Grid, Params
from render import frame_data, grid_svg, line_chart_svg, viewer_html
from rules import step
from scenarios import CAPILLARY_SPACING, calibrated_params, healthy_tissue


def warburg_tumor_cell() -> Cell:
    """Small mitochondria, high glucose uptake (ATP 2(8) + 7.5 = 23.5), acid tolerant."""
    return Cell(TUMOR, mito_capacity=0.5, acid_tolerance=3.0, glucose_uptake=8.0)


def seed_tumor(grid: Grid, center_row: int, center_col: int, radius: int, cell_maker) -> None:
    """Replace the cells in a (2 radius + 1) square around the center with tumor cells (skips vessels)."""
    vessels = set(grid.vessels)
    for r in range(center_row - radius, center_row + radius + 1):
        for c in range(center_col - radius, center_col + radius + 1):
            if (r, c) not in vessels:
                grid.cells[r][c] = cell_maker()


def run(grid: Grid, params: Params, steps: int, record_every: int, seed: int,
        oxygen_top: float, lactate_top: float) -> tuple[list[dict], list[Grid]]:
    """Run the simulation. Returns (viewer frames, the grid at every recorded step)."""
    rng = random.Random(seed)
    frames, grids = [frame_data(grid, 0, oxygen_top, lactate_top)], [grid]
    for n in range(1, steps + 1):
        grid = step(grid, params, rng)
        if n % record_every == 0:
            frames.append(frame_data(grid, n, oxygen_top, lactate_top))
            grids.append(grid)
    return frames, grids


def main() -> None:
    size, steps, record_every, seed = 60, 150, 2, 1
    params = calibrated_params()
    grid = healthy_tissue(size, size, CAPILLARY_SPACING, params)
    seed_tumor(grid, size // 2, size // 2, 1, warburg_tumor_cell)
    oxygen_top, lactate_top = params.vessel_oxygen, 4.0

    out_dir = os.path.join("output", "run1")
    os.makedirs(out_dir, exist_ok=True)
    started = time.time()
    frames, grids = run(grid, params, steps, record_every, seed, oxygen_top, lactate_top)
    print(f"Ran {steps} steps on {size}x{size} in {time.time() - started:.0f}s")

    with open(os.path.join(out_dir, "viewer.html"), "w") as file:
        file.write(viewer_html(frames, size, size, "Warburg simulation, run 1"))
    for frame, snapshot in zip(frames, grids):
        if frame["step"] in (0, steps // 3, 2 * steps // 3, steps):
            svg = grid_svg(snapshot, f"step {frame['step']}: {frame['healthy']} healthy, "
                                     f"{frame['tumor']} tumor", oxygen_top, lactate_top, square=5)
            with open(os.path.join(out_dir, f"step_{frame['step']:03d}.svg"), "w") as file:
                file.write(svg)
    chart = line_chart_svg({"healthy": [f["healthy"] for f in frames], "tumor": [f["tumor"] for f in frames]},
                           {"healthy": "#1D9E75", "tumor": "#D85A30"},
                           "Cells over time (run 1)", "cells", x_step=record_every)
    with open(os.path.join(out_dir, "counts.svg"), "w") as file:
        file.write(chart)

    for frame in frames[:: max(1, len(frames) // 10)] + [frames[-1]]:
        print(f"  step {frame['step']:>3}: {frame['healthy']:>5} healthy, {frame['tumor']:>4} tumor")
    print(f"Wrote {out_dir}/viewer.html, snapshots and counts.svg")


if __name__ == "__main__":
    main()
