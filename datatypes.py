"""The things that live on the grid."""

EMPTY = "empty"
HEALTHY = "healthy"
TUMOR = "tumor"
KINDS = (EMPTY, HEALTHY, TUMOR)


class Cell:
    """
    One square of tissue.
    - kind: EMPTY, HEALTHY or TUMOR
    - mito_capacity: most pyruvate per step its mitochondria can burn (M in solve_cell)
    - acid_tolerance: dies if local lactate is above this
    """
    def __init__(self, kind: str, mito_capacity: float = 0.0, acid_tolerance: float = 0.0) -> None:
        if kind not in KINDS:
            raise ValueError(f"kind must be one of {KINDS}, not {kind!r}")
        self.kind = kind
        self.mito_capacity = mito_capacity
        self.acid_tolerance = acid_tolerance

    def __repr__(self) -> str:
        return f"Cell({self.kind!r}, mito={self.mito_capacity:g}, acid_tol={self.acid_tolerance:g})"


def empty_cell() -> Cell:
    """A square with no cell in it."""
    return Cell(EMPTY)
