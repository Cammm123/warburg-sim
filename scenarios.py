"""
Ready-made settings and starting grids.

calibrated_params() comes from calibrate.py (2026-10-10): with these numbers a
30x30 grid of healthy cells, capillaries every 6 squares, stays 100% alive for
200 steps, and the results are identical at 10, 20 and 40 substeps (converged).
Glucose and O2 both reach every cell, but O2 is only just enough (lowest about
2.4), so extra demand from a tumor can push its neighborhood into hypoxia.
"""
from datatypes import HEALTHY, Cell, Grid, Params, empty_cell
from fields import make_field
from rules import capillaries

CAPILLARY_SPACING = 6


def calibrated_params(**overrides: float) -> Params:
    """Params tuned so healthy tissue alone is stable. Keyword args override any value."""
    settings = {"diffusion_substeps": 10, "uptake_fraction": 0.02,
                "vessel_oxygen": 10.0, "vessel_glucose": 5.0}
    settings.update(overrides)
    return Params(**settings)


def healthy_cell() -> Cell:
    """The standard healthy cell: big mitochondria, low acid tolerance, normal glucose uptake."""
    return Cell(HEALTHY, mito_capacity=10.0, acid_tolerance=1.0, glucose_uptake=1.0)


def healthy_tissue(num_rows: int, num_cols: int, spacing: int, params: Params) -> Grid:
    """Healthy cells everywhere except capillaries; fields start at blood levels."""
    vessels = capillaries(num_rows, num_cols, spacing)
    vessel_set = set(vessels)
    cells = [[empty_cell() if (r, c) in vessel_set else healthy_cell() for c in range(num_cols)]
             for r in range(num_rows)]
    return Grid(cells, make_field(num_rows, num_cols, params.vessel_oxygen),
                make_field(num_rows, num_cols, params.vessel_glucose),
                make_field(num_rows, num_cols), vessels)
