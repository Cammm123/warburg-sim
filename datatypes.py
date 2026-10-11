"""The things that live on the grid, and the simulation's settings."""

EMPTY = "empty"
HEALTHY = "healthy"
TUMOR = "tumor"
KINDS = (EMPTY, HEALTHY, TUMOR)

Field = list[list[float]]


class Cell:
    """
    One square of tissue.
    - kind: EMPTY, HEALTHY or TUMOR
    - mito_capacity: most pyruvate per step its mitochondria can burn (M in solve_cell)
    - acid_tolerance: dies if local lactate is above this
    - glucose_uptake: most glucose per step it can take in (tumors make more transporters)
    """
    def __init__(self, kind: str, mito_capacity: float = 0.0, acid_tolerance: float = 0.0,
                 glucose_uptake: float = 1.0) -> None:
        if kind not in KINDS:
            raise ValueError(f"kind must be one of {KINDS}, not {kind!r}")
        self.kind = kind
        self.mito_capacity = mito_capacity
        self.acid_tolerance = acid_tolerance
        self.glucose_uptake = glucose_uptake

    def __repr__(self) -> str:
        return (f"Cell({self.kind!r}, mito={self.mito_capacity:g}, "
                f"acid_tol={self.acid_tolerance:g}, uptake={self.glucose_uptake:g})")


def empty_cell() -> Cell:
    """A square with no cell in it."""
    return Cell(EMPTY)


def daughter_of(parent: Cell) -> Cell:
    """A new cell with the same traits as its parent (mutation comes in step 8)."""
    return Cell(parent.kind, parent.mito_capacity, parent.acid_tolerance, parent.glucose_uptake)


class Params:
    """Every number the simulation uses, in one place, with defaults."""
    def __init__(self, survive_atp: float = 5.0, divide_atp: float = 20.0,
                 oxygen_rate: float = 1.0, glucose_rate: float = 0.5, lactate_rate: float = 0.5,
                 vessel_oxygen: float = 10.0, vessel_glucose: float = 2.0,
                 tumor_weight: float = 3.0, healthy_weight: float = 1.0) -> None:
        self.survive_atp = survive_atp          # ATP needed to stay alive
        self.divide_atp = divide_atp            # ATP needed to divide
        self.oxygen_rate = oxygen_rate          # diffusion rates (0 to 1)
        self.glucose_rate = glucose_rate
        self.lactate_rate = lactate_rate
        self.vessel_oxygen = vessel_oxygen      # levels in the blood
        self.vessel_glucose = vessel_glucose
        self.tumor_weight = tumor_weight        # dice weights when cells compete for a square
        self.healthy_weight = healthy_weight


class Grid:
    """
    The whole tissue: a 2D list of Cells, three chemical fields of the same
    shape, and the blood-vessel squares (which never hold a cell).
    Raises ValueError if shapes differ or a vessel square holds a cell.
    """
    def __init__(self, cells: list[list[Cell]], oxygen: Field, glucose: Field, lactate: Field,
                 vessels: list[tuple[int, int]]) -> None:
        num_rows, num_cols = len(cells), len(cells[0])
        for field in (oxygen, glucose, lactate):
            if len(field) != num_rows or any(len(row) != num_cols for row in field):
                raise ValueError("every field must have the same shape as cells")
        for r, c in vessels:
            if cells[r][c].kind != EMPTY:
                raise ValueError(f"vessel square {(r, c)} holds a cell")
        self.cells = cells
        self.oxygen = oxygen
        self.glucose = glucose
        self.lactate = lactate
        self.vessels = list(vessels)

    def num_rows(self) -> int:
        return len(self.cells)

    def num_cols(self) -> int:
        return len(self.cells[0])

    def count(self, kind: str) -> int:
        """How many squares hold this kind of cell."""
        return sum(1 for row in self.cells for cell in row if cell.kind == kind)
