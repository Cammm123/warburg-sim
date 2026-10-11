"""Rules for what each cell does on a step."""
from datatypes import EMPTY, Cell

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
