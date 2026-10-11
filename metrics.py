"""Numbers that summarize a grid, for tracking evolution over time."""
from datatypes import TUMOR, Grid
from fba import solve_cell, warburg_index


def tumor_traits(grid: Grid) -> list[float]:
    """mito_capacity of every tumor cell."""
    return [cell.mito_capacity for row in grid.cells for cell in row if cell.kind == TUMOR]


def mean(values: list[float]) -> float:
    """Average, or 0.0 for an empty list."""
    return sum(values) / len(values) if values else 0.0


def tumor_warburg_indices(grid: Grid) -> list[float]:
    """
    Each tumor cell's Warburg index (share of glucose fermented) under the
    glucose and O2 on its square right now.
    """
    indices = []
    for r, row in enumerate(grid.cells):
        for c, cell in enumerate(row):
            if cell.kind == TUMOR:
                glucose_in = min(cell.glucose_uptake, grid.glucose[r][c])
                indices.append(warburg_index(solve_cell(glucose_in, grid.oxygen[r][c], cell.mito_capacity)))
    return indices


def summary(grid: Grid) -> dict[str, float]:
    """Mean tumor mito_capacity and mean tumor Warburg index (0 if no tumor)."""
    return {"mean_mito": mean(tumor_traits(grid)), "mean_warburg": mean(tumor_warburg_indices(grid))}
