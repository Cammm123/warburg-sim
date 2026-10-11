"""Rules for what each cell does on a step, and one full step of the grid."""
import random

from datatypes import EMPTY, TUMOR, Cell, Field, Grid, Params, daughter_of, empty_cell
from fba import atp, oxygen_used, solve_cell
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


def metabolism(grid: Grid) -> tuple[Field, Field, Field, Field]:
    """
    Every living cell runs solve_cell on its own square's glucose and O2.
    Glucose taken in is also capped by the cell's glucose_uptake.
    Returns (atp_made, glucose_change, oxygen_change, lactate_change), each a
    field: ATP per cell, and how each square's chemicals change (used up is
    negative, released is positive). Empty squares make and change nothing.
    """
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
            glucose_in = min(cell.glucose_uptake, grid.glucose[r][c])
            fluxes = solve_cell(glucose_in, grid.oxygen[r][c], cell.mito_capacity)
            atp_made[r][c] = atp(fluxes)
            glucose_change[r][c] = -fluxes["v_in"]
            oxygen_change[r][c] = -oxygen_used(fluxes)
            lactate_change[r][c] = fluxes["v_out"]
    return atp_made, glucose_change, oxygen_change, lactate_change


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


def step(grid: Grid, params: Params, rng: random.Random) -> Grid:
    """
    One step of the simulation. Returns a new Grid (the input is not changed).
      1. Vessels refill O2 and glucose and drain lactate.
      2. Metabolism: each cell runs solve_cell on its square; fields change.
      3. O2, glucose and lactate diffuse.
      4. Every cell's fate is decided at once (cell_fate, with this step's ATP
         and its square's lactate after diffusion).
      5. Deaths happen, then divisions into empty neighbors (choose_divisions).
    """
    oxygen = refill_sources(grid.oxygen, grid.vessels, params.vessel_oxygen)
    glucose = refill_sources(grid.glucose, grid.vessels, params.vessel_glucose)
    lactate = refill_sources(grid.lactate, grid.vessels, 0.0)

    fed = Grid(grid.cells, oxygen, glucose, lactate, grid.vessels)
    atp_made, glucose_change, oxygen_change, lactate_change = metabolism(fed)
    oxygen = apply_change(oxygen, oxygen_change)
    glucose = apply_change(glucose, glucose_change)
    lactate = apply_change(lactate, lactate_change)

    oxygen = diffuse(oxygen, params.oxygen_rate)
    glucose = diffuse(glucose, params.glucose_rate)
    lactate = diffuse(lactate, params.lactate_rate)

    num_rows, num_cols = grid.num_rows(), grid.num_cols()
    cells = [row[:] for row in grid.cells]
    dividers = []
    for r in range(num_rows):
        for c in range(num_cols):
            fate = cell_fate(grid.cells[r][c], atp_made[r][c], lactate[r][c],
                             params.survive_atp, params.divide_atp)
            if fate == DIE:
                cells[r][c] = empty_cell()
            elif fate == DIVIDE:
                dividers.append((r, c))

    for (tr, tc), (pr, pc) in choose_divisions(cells, dividers, grid.vessels, params, rng).items():
        cells[tr][tc] = daughter_of(cells[pr][pc])

    return Grid(cells, oxygen, glucose, lactate, grid.vessels)
