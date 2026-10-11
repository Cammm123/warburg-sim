"""
Write golden.json: reference answers from the Python model, so the JavaScript
port can be checked number for number (model.test.js).

Run from warburg-sim/web:  python3 make_golden.py
Only deterministic cases go in (no random choices), so the two languages must
agree to about 1e-9.
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from datatypes import HEALTHY, TUMOR, Cell, Grid  # noqa: E402
from fba import (best_budget_amounts, budget_traits, pathway_amounts, solve_cell,  # noqa: E402
                 stoichiometric_matrix)
from fields import diffuse, make_field  # noqa: E402
from linalg import integer_vector, left_null_space, null_space, rank  # noqa: E402
from rules import metabolism, step  # noqa: E402
from scenarios import CAPILLARY_SPACING, calibrated_params, healthy_tissue  # noqa: E402


def cell_json(cell: Cell) -> dict:
    """A Cell as a JSON object with the JS field names."""
    return {"kind": cell.kind, "mitoCapacity": cell.mito_capacity, "acidTolerance": cell.acid_tolerance,
            "glucoseUptake": cell.glucose_uptake, "mitoShare": cell.mito_share, "marker": cell.marker}


def grid_json(grid: Grid) -> dict:
    """A Grid as JSON (cells, three fields, vessels)."""
    return {"cells": [[cell_json(c) for c in row] for row in grid.cells], "oxygen": grid.oxygen,
            "glucose": grid.glucose, "lactate": grid.lactate, "vessels": [list(v) for v in grid.vessels]}


def params_json(params) -> dict:
    """Params with the JS field names."""
    return {"surviveAtp": params.survive_atp, "divideAtp": params.divide_atp,
            "oxygenRate": params.oxygen_rate, "glucoseRate": params.glucose_rate,
            "lactateRate": params.lactate_rate, "vesselOxygen": params.vessel_oxygen,
            "vesselGlucose": params.vessel_glucose, "tumorWeight": params.tumor_weight,
            "healthyWeight": params.healthy_weight, "diffusionSubsteps": params.diffusion_substeps,
            "uptakeFraction": params.uptake_fraction, "mutationSd": params.mutation_sd,
            "budget": params.budget, "mitoCost": params.mito_cost, "plastic": params.plastic,
            "acidHarm": params.acid_harm}


def main() -> None:
    rng = random.Random(0)
    golden = {}

    s = stoichiometric_matrix()
    golden["rank"] = rank(s)
    golden["nullSpace"] = sorted(integer_vector(v) for v in null_space(s))
    golden["leftNullSpace"] = [integer_vector(v) for v in left_null_space(s)]

    values = [0, 0.1, 0.5, 1, 2.5, 6, 10]
    golden["solveCell"] = [{"in": [g, o, m], "out": solve_cell(g, o, m)}
                           for g in values for o in values for m in values]
    golden["pathwayAmounts"] = [{"in": [g, o, m], "out": list(pathway_amounts(g, o, m))}
                                for g in values for o in values for m in values]
    golden["budgetTraits"] = [{"in": [share, 20, 10], "out": list(budget_traits(share, 20, 10))}
                              for share in (0, 0.25, 0.5, 0.9, 1)]
    cases = []
    for _ in range(300):
        args = [rng.uniform(0, 12), rng.uniform(0, 12), rng.uniform(0.5, 30), rng.uniform(0.1, 20)]
        cases.append({"in": args, "out": list(best_budget_amounts(*args))})
    golden["bestBudgetAmounts"] = cases

    spike = make_field(5, 5)
    spike[2][2] = 3.0
    random_field = [[rng.uniform(0, 5) for _ in range(6)] for _ in range(4)]
    golden["diffuse"] = [{"field": spike, "rate": 0.5, "steps": 3},
                         {"field": random_field, "rate": 0.7, "steps": 4}]
    for case in golden["diffuse"]:
        f = case["field"]
        for _ in range(case["steps"]):
            f = diffuse(f, case["rate"])
        case["out"] = f

    # Metabolism on a small mixed grid: fixed-trait, evolving budget and plastic.
    cells = [[Cell(HEALTHY, 10, 1.0), Cell(TUMOR, 0.5, 3.0, 8.0), Cell(TUMOR, acid_tolerance=3.0, mito_share=0.3)],
             [Cell(TUMOR, acid_tolerance=3.0, mito_share=0.8), Cell(HEALTHY, 10, 1.0, 0.5), Cell(TUMOR, 2.0, 3.0, 4.0)]]
    oxygen = [[rng.uniform(0, 10) for _ in range(3)] for _ in range(2)]
    glucose = [[rng.uniform(0, 5) for _ in range(3)] for _ in range(2)]
    grid = Grid(cells, oxygen, glucose, make_field(2, 3), [])
    golden["metabolism"] = []
    for params, fraction in [(calibrated_params(), 1.0), (calibrated_params(), 0.02),
                             (calibrated_params(plastic=True), 0.02)]:
        out = metabolism(grid, fraction, params)
        golden["metabolism"].append({"grid": grid_json(grid), "params": params_json(params),
                                     "fraction": fraction, "out": [list(map(list, f)) for f in out]})

    # Full steps with no deaths, so no random choices are made.
    params = calibrated_params()
    tissue = healthy_tissue(12, 12, CAPILLARY_SPACING, params)
    tissue.cells[5][5] = Cell(TUMOR, acid_tolerance=3.0, mito_share=0.6)
    tissue.cells[6][6] = Cell(TUMOR, 0.5, 3.0, 2.0)
    start = grid_json(tissue)
    g = tissue
    for _ in range(4):
        before = [[c.kind for c in row] for row in g.cells]
        g = step(g, params, random.Random(0))
        if [[c.kind for c in row] for row in g.cells] != before:
            raise RuntimeError("a cell died or divided: this case would need random choices")
    golden["step"] = {"grid": start, "params": params_json(params), "steps": 4, "out": grid_json(g)}

    with open("golden.json", "w") as file:
        json.dump(golden, file)
    print(f"wrote golden.json ({os.path.getsize('golden.json') // 1024} KB)")


if __name__ == "__main__":
    main()
