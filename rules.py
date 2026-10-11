"""Rules for what each cell does on a step, and one full step of the grid."""
import random

from datatypes import EMPTY, TUMOR, Cell, Field, Grid, Params, daughter_of, empty_cell
from fba import (atp, best_budget_amounts, budget_traits, fluxes_from_amounts, oxygen_used,
                 solve_cell)
from fields import apply_change, diffuse, make_field, refill_sources

DIE = "die"
DIVIDE = "divide"
STAY = "stay"
NOTHING = "nothing"   # an empty square


def neighbors(row: int, col: int, num_rows: int, num_cols: int) -> list[tuple[int, int]]:
    """
    The up-to-8 squares around (row, col) that are on the grid (Moore
    neighborhood, fixed edges), in reading order: top row left to right, then
    the middle row, then the bottom row. Same idea as neighbors_fixed in
    02-120/spatial-games/functions.py, without the square itself.
    """
    result = []
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            r, c = row + dr, col + dc
            if 0 <= r < num_rows and 0 <= c < num_cols:
                result.append((r, c))
    return result


def cell_fate(cell: Cell, cell_atp: float, lactate: float,
              survive_atp: float, divide_atp: float) -> str:
    """
    Decide what a cell does this step, checking in this order:
      1. empty square                      -> NOTHING
      2. lactate above acid_tolerance      -> DIE   (acid kills, whatever the ATP)
      3. ATP below survive_atp             -> DIE   (starves)
      4. ATP at or above divide_atp        -> DIVIDE (whether there's room is
                                              decided later, by the grid step)
      5. otherwise                         -> STAY
    Lactate exactly at the tolerance survives. Raises ValueError if
    divide_atp < survive_atp (a cell could divide but not survive).
    """
    if divide_atp < survive_atp:
        raise ValueError("divide_atp must be at least survive_atp")
    if cell.kind == EMPTY:
        return NOTHING
    if lactate > cell.acid_tolerance:
        return DIE
    if cell_atp < survive_atp:
        return DIE
    if cell_atp >= divide_atp:
        return DIVIDE
    return STAY


def cell_fluxes(cell: Cell, glucose_here: float, oxygen_here: float, params: Params,
                fraction: float = 1.0) -> dict[str, float]:
    """
    How one living cell metabolizes on its square, for `fraction` of a step:
    - fixed-trait cell (mito_share None): solve_cell with its uptake and M;
    - budget cell, evolving: uptake and M come from its inherited mito_share;
    - budget cell, plastic (params.plastic): best split for the local glucose and O2.
    Uptake, M and budget all scale with fraction.
    """
    if cell.mito_share is None:
        return solve_cell(min(cell.glucose_uptake * fraction, glucose_here), oxygen_here,
                          cell.mito_capacity * fraction)
    if params.plastic:
        a, b = best_budget_amounts(glucose_here, oxygen_here, params.budget * fraction, params.mito_cost)
        return fluxes_from_amounts(a, b)
    uptake, mito = budget_traits(cell.mito_share, params.budget, params.mito_cost)
    return solve_cell(min(uptake * fraction, glucose_here), oxygen_here, mito * fraction)


def metabolism(grid: Grid, fraction: float = 1.0, params: Params | None = None) -> tuple[Field, Field, Field, Field]:
    """
    Every living cell runs solve_cell on its own square's glucose and O2.
    Glucose taken in is also capped by the cell's glucose_uptake.
    fraction (0 to 1] is how much of a full step this covers: uptake and mito
    capacity are scaled by it, so k calls with fraction 1/k add up to one step
    when supply is plentiful (FBA scales: 1/k of the inputs -> 1/k of the fluxes).
    params is needed only for budget cells (default Params()).
    Returns (atp_made, glucose_change, oxygen_change, lactate_change), each a
    field: ATP per cell, and how each square's chemicals change (used up is
    negative, released is positive). Empty squares make and change nothing.
    """
    params = params or Params()
    num_rows, num_cols = grid.num_rows(), grid.num_cols()
    atp_made = make_field(num_rows, num_cols)
    glucose_change = make_field(num_rows, num_cols)
    oxygen_change = make_field(num_rows, num_cols)
    lactate_change = make_field(num_rows, num_cols)
    for r in range(num_rows):
        for c in range(num_cols):
            cell = grid.cells[r][c]
            if cell.kind == EMPTY:
                continue
            fluxes = cell_fluxes(cell, grid.glucose[r][c], grid.oxygen[r][c], params, fraction)
            atp_made[r][c] = atp(fluxes)
            glucose_change[r][c] = -fluxes["v_in"]
            oxygen_change[r][c] = -oxygen_used(fluxes)
            lactate_change[r][c] = fluxes["v_out"]
    return atp_made, glucose_change, oxygen_change, lactate_change


def make_daughter(parent: Cell, params: Params, rng: random.Random) -> Cell:
    """
    The new cell from a division. Healthy cells copy exactly. Tumor cells
    (genetically unstable) mutate one trait with Gaussian noise of standard
    deviation params.mutation_sd: mito_share (kept in [0, 1]) for budget cells,
    otherwise mito_capacity (kept >= 0). The neutral marker gets its own noise
    of the same size. Other traits copy exactly.
    """
    daughter = daughter_of(parent)
    if parent.kind == TUMOR and params.mutation_sd > 0:
        noise = rng.gauss(0.0, params.mutation_sd)
        if parent.mito_share is not None:
            daughter.mito_share = min(1.0, max(0.0, parent.mito_share + noise))
        else:
            daughter.mito_capacity = max(0.0, parent.mito_capacity + noise)
        daughter.marker = parent.marker + rng.gauss(0.0, params.mutation_sd)
    return daughter


def claim_weight(cell: Cell, params: Params) -> float:
    """A cell's weight on the dice when competing for an empty square."""
    return params.tumor_weight if cell.kind == TUMOR else params.healthy_weight


def choose_divisions(cells: list[list[Cell]], dividers: list[tuple[int, int]],
                     vessels: list[tuple[int, int]], params: Params,
                     rng: random.Random) -> dict[tuple[int, int], tuple[int, int]]:
    """
    Decide which divisions happen. Returns {empty square: parent square}.
      1. Each divider claims one random empty, non-vessel neighbor (if any).
      2. Each claimed square rolls weighted dice among its claimants:
         a cell's chance = its weight / total weight (tumor_weight vs
         healthy_weight; same kind = even odds).
    Losers don't divide this step. Ties are fair whatever order dividers are listed in.
    """
    num_rows, num_cols = len(cells), len(cells[0])
    vessel_set = set(vessels)
    claims: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for r, c in dividers:
        options = [(nr, nc) for nr, nc in neighbors(r, c, num_rows, num_cols)
                   if cells[nr][nc].kind == EMPTY and (nr, nc) not in vessel_set]
        if options:
            target = rng.choice(options)
            claims.setdefault(target, []).append((r, c))
    winners = {}
    for target in sorted(claims):
        claimants = claims[target]
        weights = [claim_weight(cells[r][c], params) for r, c in claimants]
        winners[target] = rng.choices(claimants, weights=weights)[0]
    return winners


def refill_vessels(oxygen: Field, glucose: Field, lactate: Field, vessels: list[tuple[int, int]],
                   params: Params) -> tuple[Field, Field, Field]:
    """Vessels reset O2 and glucose to blood levels and drain lactate to 0."""
    return (refill_sources(oxygen, vessels, params.vessel_oxygen),
            refill_sources(glucose, vessels, params.vessel_glucose),
            refill_sources(lactate, vessels, 0.0))


def capillaries(num_rows: int, num_cols: int, spacing: int) -> list[tuple[int, int]]:
    """
    Vessel squares spread evenly, like capillaries: every spacing-th row and
    column, starting spacing // 2 in from the edge so they sit inside the tissue.
    """
    if spacing < 1:
        raise ValueError("spacing must be at least 1")
    start = spacing // 2
    return [(r, c) for r in range(start, num_rows, spacing) for c in range(start, num_cols, spacing)]


def step(grid: Grid, params: Params, rng: random.Random) -> Grid:
    """
    One step of the simulation. Returns a new Grid (the input is not changed).
      1-3. Repeated params.diffusion_substeps (k) times, because chemicals
         move and get used much faster than cells divide:
           1. Vessels refill O2 and glucose and drain lactate.
           2. Metabolism: each cell eats uptake_fraction (f) of a full step's
              appetite (solve_cell on its square); fields change.
           3. O2, glucose and lactate diffuse.
         A cell's ATP for its fate = its average ATP rate over the substeps,
         i.e. (ATP summed over substeps) / (k x f).
      4. Every cell's fate is decided at once (cell_fate, with this step's ATP
         and its square's lactate after diffusion).
      5. Deaths happen, then divisions into empty neighbors (choose_divisions);
         tumor daughters mutate (make_daughter).
    """
    num_rows, num_cols = grid.num_rows(), grid.num_cols()
    k, f = params.diffusion_substeps, params.uptake_fraction
    oxygen, glucose, lactate = grid.oxygen, grid.glucose, grid.lactate
    atp_made = make_field(num_rows, num_cols)
    for _ in range(k):
        oxygen, glucose, lactate = refill_vessels(oxygen, glucose, lactate, grid.vessels, params)
        fed = Grid(grid.cells, oxygen, glucose, lactate, grid.vessels)
        atp_part, glucose_change, oxygen_change, lactate_change = metabolism(fed, f, params)
        atp_made = [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(atp_made, atp_part)]
        oxygen = apply_change(oxygen, oxygen_change)
        glucose = apply_change(glucose, glucose_change)
        lactate = apply_change(lactate, lactate_change)
        oxygen = diffuse(oxygen, params.oxygen_rate)
        glucose = diffuse(glucose, params.glucose_rate)
        lactate = diffuse(lactate, params.lactate_rate)

    atp_made = [[value / (k * f) for value in row] for row in atp_made]

    cells = [row[:] for row in grid.cells]
    dividers = []
    for r in range(num_rows):
        for c in range(num_cols):
            acid = lactate[r][c] if params.acid_harm else 0.0
            fate = cell_fate(grid.cells[r][c], atp_made[r][c], acid,
                             params.survive_atp, params.divide_atp)
            if fate == DIE:
                cells[r][c] = empty_cell()
            elif fate == DIVIDE:
                dividers.append((r, c))

    winners = choose_divisions(cells, dividers, grid.vessels, params, rng)
    for (tr, tc), (pr, pc) in sorted(winners.items()):
        cells[tr][tc] = make_daughter(cells[pr][pc], params, rng)

    return Grid(cells, oxygen, glucose, lactate, grid.vessels)
