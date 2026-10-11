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
    - mito_share: None for a cell with fixed traits (above). A number from 0 to 1
      makes it a budget cell: it splits a protein budget between transporters
      and mitochondria, spending this share on mitochondria (see fba.budget_traits).
    - marker: a neutral trait that mutates like the real one but does nothing.
      How far it drifts shows how much a trait changes by chance alone.
    """
    def __init__(self, kind: str, mito_capacity: float = 0.0, acid_tolerance: float = 0.0,
                 glucose_uptake: float = 1.0, mito_share: float | None = None,
                 marker: float = 0.0) -> None:
        if kind not in KINDS:
            raise ValueError(f"kind must be one of {KINDS}, not {kind!r}")
        self.kind = kind
        self.mito_capacity = mito_capacity
        self.acid_tolerance = acid_tolerance
        self.glucose_uptake = glucose_uptake
        if mito_share is not None and not 0 <= mito_share <= 1:
            raise ValueError("mito_share must be between 0 and 1 (or None)")
        self.mito_share = mito_share
        self.marker = marker

    def __repr__(self) -> str:
        share = "" if self.mito_share is None else f", share={self.mito_share:g}"
        return (f"Cell({self.kind!r}, mito={self.mito_capacity:g}, "
                f"acid_tol={self.acid_tolerance:g}, uptake={self.glucose_uptake:g}{share})")


def empty_cell() -> Cell:
    """A square with no cell in it."""
    return Cell(EMPTY)


def daughter_of(parent: Cell) -> Cell:
    """A new cell with the same traits as its parent (mutation is rules.make_daughter)."""
    return Cell(parent.kind, parent.mito_capacity, parent.acid_tolerance, parent.glucose_uptake,
                parent.mito_share, parent.marker)


class Params:
    """Every number the simulation uses, in one place, with defaults."""
    def __init__(self, survive_atp: float = 5.0, divide_atp: float = 20.0,
                 oxygen_rate: float = 1.0, glucose_rate: float = 0.5, lactate_rate: float = 0.5,
                 vessel_oxygen: float = 10.0, vessel_glucose: float = 2.0,
                 tumor_weight: float = 3.0, healthy_weight: float = 1.0,
                 diffusion_substeps: int = 1, uptake_fraction: float = 1.0,
                 mutation_sd: float = 0.0, budget: float = 20.0, mito_cost: float = 10.0,
                 plastic: bool = False, acid_harm: bool = True) -> None:
        self.survive_atp = survive_atp          # ATP needed to stay alive
        self.divide_atp = divide_atp            # ATP needed to divide
        self.oxygen_rate = oxygen_rate          # diffusion rates (0 to 1)
        self.glucose_rate = glucose_rate
        self.lactate_rate = lactate_rate
        self.vessel_oxygen = vessel_oxygen      # levels in the blood
        self.vessel_glucose = vessel_glucose
        self.tumor_weight = tumor_weight        # dice weights when cells compete for a square
        self.healthy_weight = healthy_weight
        # Numerics: each cell step runs this many substeps (refill, eat, diffuse),
        # and a cell's ATP is its average rate over them. Raise it until results
        # stop changing (convergence).
        if diffusion_substeps < 1:
            raise ValueError("diffusion_substeps must be at least 1")
        self.diffusion_substeps = diffusion_substeps
        # Physics: share of a full step's appetite a cell eats per substep, i.e.
        # how fast cells consume compared with how fast chemicals spread.
        # Smaller = chemicals reach further before running out.
        if not 0 < uptake_fraction <= 1:
            raise ValueError("uptake_fraction must be in (0, 1]")
        self.uptake_fraction = uptake_fraction
        # Evolution: when a tumor cell divides, its daughter's mito_capacity is
        # the parent's plus random noise with this standard deviation (0 = exact copy).
        if mutation_sd < 0:
            raise ValueError("mutation_sd must be >= 0")
        self.mutation_sd = mutation_sd
        # Budget cells (mito_share not None): uptake + mito_cost x M = budget.
        # plastic=True: they pick the best split every substep instead of using
        # their inherited mito_share.
        if budget <= 0 or mito_cost <= 0:
            raise ValueError("budget and mito_cost must be positive")
        self.budget = budget
        self.mito_cost = mito_cost
        self.plastic = plastic
        # Experiment switch: False means lactate never kills anyone.
        self.acid_harm = acid_harm


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
