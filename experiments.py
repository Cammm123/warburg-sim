"""
Step 9: experiments. No pictures, several seeds per setting.

E1  acid harm on vs off, for the free-mitochondria model (trait M) and the
    budget model (trait mito_share). Off = lactate never kills anyone.
E2  budget model: blood glucose x healthy acid tolerance, to map where
    fermentation (Warburg) wins. (A first sweep at glucose 5/10/20 and tolerance
    1/3/inf was flat: glucose >= 5 is never scarce per substep, and tolerance 3
    equals the tumor's, so it acts like acid off.)
Every evolving run also tracks the neutral marker: if the real trait moves much
more than the marker, that is selection, not chance.

Run from this folder:  python3 experiments.py
Writes output/experiments/results.md (tables) and prints them.
"""
import os
import random
import statistics
import time

from datatypes import HEALTHY, TUMOR, Cell, Params
from metrics import summary
from rules import step
from scenarios import CAPILLARY_SPACING, calibrated_params, healthy_tissue

SIZE, STEPS, SEEDS = 40, 200, (1, 2, 3)


def free_mito_tumor() -> Cell:
    """Fixed-trait Warburg tumor cell (as in run 1 and 2); M evolves."""
    return Cell(TUMOR, mito_capacity=0.5, acid_tolerance=3.0, glucose_uptake=8.0)


def budget_tumor() -> Cell:
    """Budget tumor cell starting at a 50% mitochondria share (as in run 3)."""
    return Cell(TUMOR, acid_tolerance=3.0, mito_share=0.5)


def run_once(params: Params, tumor_maker, seed: int, healthy_tolerance: float = 1.0,
             size: int = SIZE, steps: int = STEPS) -> dict[str, float]:
    """One run: healthy tissue + 3x3 tumor seed in the middle. Returns final counts and means."""
    grid = healthy_tissue(size, size, CAPILLARY_SPACING, params)
    vessels = set(grid.vessels)
    for row in grid.cells:
        for cell in row:
            if cell.kind == HEALTHY:
                cell.acid_tolerance = healthy_tolerance
    mid = size // 2
    for r in range(mid - 1, mid + 2):
        for c in range(mid - 1, mid + 2):
            if (r, c) not in vessels:
                grid.cells[r][c] = tumor_maker()
    start = summary(grid, params)
    rng = random.Random(seed)
    for _ in range(steps):
        grid = step(grid, params, rng)
    end = summary(grid, params)
    return {"healthy": grid.count(HEALTHY), "tumor": grid.count(TUMOR),
            "mito_change": end["mean_mito"] - start["mean_mito"],
            "share_change": end["mean_share"] - start["mean_share"],
            "marker_change": end["mean_marker"] - start["mean_marker"],
            "warburg": end["mean_warburg"]}


def run_seeds(params: Params, tumor_maker, healthy_tolerance: float = 1.0,
              seeds: tuple[int, ...] = SEEDS, size: int = SIZE, steps: int = STEPS) -> dict[str, tuple[float, float]]:
    """run_once over several seeds. Returns {metric: (mean, standard deviation)}."""
    runs = [run_once(params, tumor_maker, seed, healthy_tolerance, size, steps) for seed in seeds]
    return {key: (statistics.mean(r[key] for r in runs),
                  statistics.stdev(r[key] for r in runs) if len(runs) > 1 else 0.0) for key in runs[0]}


def fmt(pair: tuple[float, float], digits: int = 2) -> str:
    """'mean ± sd' with the given decimals."""
    return f"{pair[0]:.{digits}f} ± {pair[1]:.{digits}f}"


def experiment_1() -> list[str]:
    """Acid on vs off, for both trait models. Returns markdown table lines."""
    lines = ["## E1: acid harm on vs off", "",
             "| model | acid | healthy left | tumor | trait change | marker change | Warburg index |",
             "|---|---|---|---|---|---|---|"]
    for model, maker, trait, sd in [("free mitochondria (M)", free_mito_tumor, "mito_change", 0.1),
                                    ("budget (share)", budget_tumor, "share_change", 0.05)]:
        for acid in (True, False):
            r = run_seeds(calibrated_params(mutation_sd=sd, acid_harm=acid), maker)
            lines.append(f"| {model} | {'on' if acid else 'off'} | {fmt(r['healthy'], 0)} | {fmt(r['tumor'], 0)} "
                         f"| {fmt(r[trait])} | {fmt(r['marker_change'])} | {fmt(r['warburg'])} |")
            print(lines[-1])
    return lines


def experiment_2() -> list[str]:
    """Budget model: blood glucose x healthy acid tolerance. Returns markdown table lines."""
    lines = ["", "## E2: budget model, blood glucose x healthy acid tolerance", "",
             "| blood glucose | healthy acid tolerance | healthy left | tumor | share change | marker change | Warburg index |",
             "|---|---|---|---|---|---|---|"]
    for glucose in (2.0, 3.0, 5.0):   # below 2, healthy tissue starves on its own
        for tolerance in (1.0, 2.0):
            params = calibrated_params(mutation_sd=0.05, vessel_glucose=glucose)
            r = run_seeds(params, budget_tumor, healthy_tolerance=tolerance)
            lines.append(f"| {glucose:g} | {tolerance:g} | {fmt(r['healthy'], 0)} | {fmt(r['tumor'], 0)} "
                         f"| {fmt(r['share_change'])} | {fmt(r['marker_change'])} | {fmt(r['warburg'])} |")
            print(lines[-1])
    return lines


def main() -> None:
    started = time.time()
    header = [f"# Warburg experiments ({SIZE}x{SIZE}, {STEPS} steps, seeds {list(SEEDS)}; mean ± sd)", "",
              "Trait change = final mean minus starting mean of the evolving trait (M or mito share). "
              "Marker change = the same for the neutral marker (chance alone). "
              "Warburg index = mean share of glucose fermented by tumor cells at the end.", ""]
    lines = header + experiment_1() + experiment_2()
    lines += ["", f"Total time: {time.time() - started:.0f}s"]
    os.makedirs(os.path.join("output", "experiments"), exist_ok=True)
    with open(os.path.join("output", "experiments", "results.md"), "w") as file:
        file.write("\n".join(lines) + "\n")
    print(f"\nWrote output/experiments/results.md ({time.time() - started:.0f}s)")


if __name__ == "__main__":
    main()
