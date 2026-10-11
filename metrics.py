"""Numbers that summarize a grid, for tracking evolution over time."""
from datatypes import TUMOR, Grid, Params
from fba import warburg_index
from rules import cell_fluxes


def tumor_traits(grid: Grid) -> list[float]:
    """mito_capacity of every tumor cell."""
    return [cell.mito_capacity for row in grid.cells for cell in row if cell.kind == TUMOR]


def mean(values: list[float]) -> float:
    """Average, or 0.0 for an empty list."""
    return sum(values) / len(values) if values else 0.0


def tumor_shares(grid: Grid) -> list[float]:
    """mito_share of every tumor budget cell."""
    return [cell.mito_share for row in grid.cells for cell in row
            if cell.kind == TUMOR and cell.mito_share is not None]


def tumor_warburg_indices(grid: Grid, params: Params | None = None) -> list[float]:
    """
    Each tumor cell's Warburg index (share of glucose fermented) under the
    glucose and O2 on its square right now (same rule as metabolism).
    """
    params = params or Params()
    indices = []
    for r, row in enumerate(grid.cells):
        for c, cell in enumerate(row):
            if cell.kind == TUMOR:
                indices.append(warburg_index(cell_fluxes(cell, grid.glucose[r][c], grid.oxygen[r][c], params)))
    return indices


def summary(grid: Grid, params: Params | None = None) -> dict[str, float]:
    """Means over tumor cells (0 if none): mito_capacity, mito_share, Warburg index."""
    return {"mean_mito": mean(tumor_traits(grid)), "mean_share": mean(tumor_shares(grid)),
            "mean_warburg": mean(tumor_warburg_indices(grid, params))}
