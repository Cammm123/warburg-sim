"""
Run the Warburg simulation and save pictures.

Run from this folder:  python3 main.py
Each scenario writes output/<name>/: viewer.html (open in a browser; Play or drag
the slider), snapshots at a few steps, counts.svg (healthy and tumor cells over
time) and evolution.svg (mean tumor mito_capacity and Warburg index over time).
  run1: Warburg tumor seed, no mutation.
  run2: same, but tumor daughters mutate mito_capacity (sd 0.1): the control
        experiment for evolution, where mitochondria cost nothing.
  run3: budget tumor (transporters vs mitochondria share one budget), starting
        at a 50% mitochondria share that evolves (sd 0.05).
  run4: budget tumor, plastic: every cell picks its best split each substep.
"""
import os
import random
import time

from datatypes import HEALTHY, TUMOR, Cell, Grid, Params
from metrics import summary
from render import frame_data, grid_svg, line_chart_svg, viewer_html
from rules import step
from scenarios import CAPILLARY_SPACING, calibrated_params, healthy_tissue


def warburg_tumor_cell() -> Cell:
    """Small mitochondria, high glucose uptake (ATP 2(8) + 7.5 = 23.5), acid tolerant."""
    return Cell(TUMOR, mito_capacity=0.5, acid_tolerance=3.0, glucose_uptake=8.0)


def budget_tumor_cell(share: float = 0.5) -> Cell:
    """A budget tumor cell: its traits come from splitting params.budget (see fba.budget_traits)."""
    return Cell(TUMOR, acid_tolerance=3.0, mito_share=share)


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


def run_scenario(name: str, params: Params, steps: int, tumor_maker=warburg_tumor_cell,
                 trait: str = "mean_mito", seed: int = 1, size: int = 60,
                 record_every: int = 2) -> list[dict]:
    """
    Healthy tissue + 3x3 tumor seed in the middle; run, save pictures, return
    frames with metrics. trait is the metric drawn in evolution.svg
    ("mean_mito" or "mean_share").
    """
    grid = healthy_tissue(size, size, CAPILLARY_SPACING, params)
    seed_tumor(grid, size // 2, size // 2, 1, tumor_maker)
    oxygen_top, lactate_top = params.vessel_oxygen, 4.0
    out_dir = os.path.join("output", name)
    os.makedirs(out_dir, exist_ok=True)

    started = time.time()
    frames, grids = run(grid, params, steps, record_every, seed, oxygen_top, lactate_top)
    for frame, snapshot in zip(frames, grids):
        frame.update(summary(snapshot, params))
    print(f"\n{name}: {steps} steps on {size}x{size} in {time.time() - started:.0f}s")

    with open(os.path.join(out_dir, "viewer.html"), "w") as file:
        file.write(viewer_html(frames, size, size, f"Warburg simulation, {name}"))
    for frame, snapshot in zip(frames, grids):
        if frame["step"] in (0, steps // 3, 2 * steps // 3, steps):
            svg = grid_svg(snapshot, f"step {frame['step']}: {frame['healthy']} healthy, "
                                     f"{frame['tumor']} tumor", oxygen_top, lactate_top, square=5)
            with open(os.path.join(out_dir, f"step_{frame['step']:03d}.svg"), "w") as file:
                file.write(svg)
    trait_label = {"mean_mito": "mean mito capacity", "mean_share": "mean mito share"}[trait]
    charts = {
        "counts.svg": line_chart_svg({"healthy": [f["healthy"] for f in frames],
                                      "tumor": [f["tumor"] for f in frames]},
                                     {"healthy": "#1D9E75", "tumor": "#D85A30"},
                                     f"Cells over time ({name})", "cells", x_step=record_every),
        "evolution.svg": line_chart_svg({trait_label: [f[trait] for f in frames],
                                         "mean Warburg index": [f["mean_warburg"] for f in frames]},
                                        {trait_label: "#534AB7", "mean Warburg index": "#D85A30"},
                                        f"Tumor traits over time ({name})", "", x_step=record_every),
    }
    for filename, chart in charts.items():
        with open(os.path.join(out_dir, filename), "w") as file:
            file.write(chart)

    shown = frames[:: max(1, len(frames) // 10)]
    if shown[-1] is not frames[-1]:
        shown.append(frames[-1])
    for frame in shown:
        print(f"  step {frame['step']:>3}: {frame['healthy']:>5} healthy, {frame['tumor']:>4} tumor, "
              f"{trait_label} {frame[trait]:.2f}, mean Warburg index {frame['mean_warburg']:.2f}")
    print(f"  wrote {out_dir}/")
    return frames


def main() -> None:
    run_scenario("run1", calibrated_params(), steps=150)
    run_scenario("run2", calibrated_params(mutation_sd=0.1), steps=300)
    run_scenario("run3", calibrated_params(mutation_sd=0.05), steps=300,
                 tumor_maker=budget_tumor_cell, trait="mean_share")
    run_scenario("run4", calibrated_params(plastic=True), steps=150, tumor_maker=budget_tumor_cell,
                 trait="mean_share")


if __name__ == "__main__":
    main()
