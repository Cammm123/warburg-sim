"""
Chemical fields: O2, glucose and lactate, each stored as a 2D grid of amounts.

Diffusion uses the same 3x3 kernel as the 02-120 Gray-Scott model:
    0.05  0.2  0.05
    0.2  -1.0  0.2
    0.05  0.2  0.05
Read it as "each step, a square sends 20% of (rate x its amount) to each side
neighbor and 5% to each corner neighbor". The -1 is everything it sent away.
"""

Field = list[list[float]]

KERNEL = [
    [0.05, 0.2, 0.05],
    [0.2, -1.0, 0.2],
    [0.05, 0.2, 0.05],
]


def make_field(num_rows: int, num_cols: int, value: float = 0.0) -> Field:
    """Return a num_rows x num_cols field with every square set to value."""
    return [[value for _ in range(num_cols)] for _ in range(num_rows)]


def diffuse(field: Field, rate: float) -> Field:
    """
    One diffusion step. Returns a new field (the input is not changed).

    Each square sends rate x KERNEL weight x its amount to each neighbor that
    exists and keeps the rest. At the edges there is no neighbor to send to, so
    nothing leaves the grid: the total amount is conserved exactly (a closed
    box). In the interior this is the same as Gray-Scott's
        new = old + rate x (kernel applied to the 3x3 neighborhood).
    rate must be between 0 and 1, otherwise a square could send away more than
    it has. Raises ValueError otherwise.
    The math can't go below 0, but float rounding can leave crumbs like -3e-18
    when a square sends away everything; those are snapped to 0.
    """
    if not 0 <= rate <= 1:
        raise ValueError("rate must be between 0 and 1")
    num_rows, num_cols = len(field), len(field[0])
    new = [row[:] for row in field]
    for r in range(num_rows):
        for c in range(num_cols):
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    if dr == 0 and dc == 0:
                        continue
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < num_rows and 0 <= nc < num_cols:
                        moved = rate * KERNEL[dr + 1][dc + 1] * field[r][c]
                        new[r][c] -= moved
                        new[nr][nc] += moved
    for row in new:
        for c, value in enumerate(row):
            if value < 0:
                if value < -1e-9:
                    raise ValueError(f"diffusion produced {value}; this should be impossible")
                row[c] = 0.0
    return new


def refill_sources(field: Field, sources: list[tuple[int, int]], level: float) -> Field:
    """
    Blood vessels: set each (row, col) in sources back to level. Returns a new
    field. A vessel square always holds the blood's concentration, however much
    the nearby cells used.
    """
    new = [row[:] for row in field]
    for r, c in sources:
        new[r][c] = level
    return new


def apply_change(field: Field, change: Field) -> Field:
    """
    Add change to field square by square and return the new field. Cells use
    things up (negative change) and release things (positive change). Tiny
    negatives from float rounding become 0. Raises ValueError if a square would
    go clearly negative (something used more than was there) or if the shapes
    differ.
    """
    if len(field) != len(change) or any(len(a) != len(b) for a, b in zip(field, change)):
        raise ValueError("field and change must have the same shape")
    new = []
    for field_row, change_row in zip(field, change):
        new_row = []
        for amount, delta in zip(field_row, change_row):
            value = amount + delta
            if value < -1e-9:
                raise ValueError(f"amount went negative ({value}); something used more than was there")
            new_row.append(max(value, 0.0))
        new.append(new_row)
    return new


def total(field: Field) -> float:
    """Sum of every square."""
    return sum(sum(row) for row in field)
