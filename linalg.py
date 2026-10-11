"""
Exact linear algebra from scratch (21-241), with fractions so nothing is rounded.

- rref, rank, null_space, left_null_space: row reduction, as in class.
- integer_vector: scale a rational vector to the smallest integer vector, using
  Euclid's algorithm for gcd (15-151). This is the Diophantine step: the null
  space over the reals is a line or plane; the integer points on it are a lattice.
"""
from fractions import Fraction

Matrix = list[list[Fraction]]


def gcd(a: int, b: int) -> int:
    """Euclid's algorithm: gcd(a, b) = gcd(b, a mod b), and gcd(a, 0) = |a|."""
    a, b = abs(a), abs(b)
    while b != 0:
        a, b = b, a % b
    return a


def lcm(a: int, b: int) -> int:
    """Least common multiple of two positive integers."""
    return a * b // gcd(a, b)


def to_fractions(matrix: list[list[float]]) -> Matrix:
    """Copy a matrix of ints (or exact floats) into Fractions."""
    return [[Fraction(entry) for entry in row] for row in matrix]


def transpose(matrix: Matrix) -> Matrix:
    """Rows become columns."""
    return [list(column) for column in zip(*matrix)]


def rref(matrix: list[list[float]]) -> tuple[Matrix, list[int]]:
    """
    Reduced row echelon form by row reduction (swap rows, scale a row, add a
    multiple of one row to another). Returns (the rref matrix, the pivot columns).
    """
    m = to_fractions(matrix)
    num_rows, num_cols = len(m), len(m[0])
    pivots = []
    row = 0
    for col in range(num_cols):
        if row == num_rows:
            break
        pivot_row = next((r for r in range(row, num_rows) if m[r][col] != 0), None)
        if pivot_row is None:
            continue                                     # no pivot here: a free column
        m[row], m[pivot_row] = m[pivot_row], m[row]      # swap
        scale = m[row][col]
        m[row] = [entry / scale for entry in m[row]]     # make the pivot 1
        for r in range(num_rows):                        # clear the rest of the column
            if r != row and m[r][col] != 0:
                factor = m[r][col]
                m[r] = [a - factor * b for a, b in zip(m[r], m[row])]
        pivots.append(col)
        row += 1
    return m, pivots


def rank(matrix: list[list[float]]) -> int:
    """Number of pivots = dimension of the column space."""
    return len(rref(matrix)[1])


def null_space(matrix: list[list[float]]) -> Matrix:
    """
    A basis for {x : A x = 0}: one vector per free column. Set that free
    variable to 1, the other free variables to 0, and read the pivot variables
    off the rref. Rank-nullity: len(result) = columns - rank.
    """
    reduced, pivots = rref(matrix)
    num_cols = len(matrix[0])
    free = [c for c in range(num_cols) if c not in pivots]
    basis = []
    for f in free:
        x = [Fraction(0)] * num_cols
        x[f] = Fraction(1)
        for row, p in enumerate(pivots):
            x[p] = -reduced[row][f]
        basis.append(x)
    return basis


def left_null_space(matrix: list[list[float]]) -> Matrix:
    """A basis for {y : y A = 0}, the null space of the transpose."""
    return null_space(transpose(to_fractions(matrix)))


def integer_vector(vector: list[Fraction]) -> list[int]:
    """
    Smallest integer vector pointing the same way: multiply by the lcm of the
    denominators, then divide by the gcd of the entries. The first nonzero entry
    is made positive. Raises ValueError for the zero vector.
    """
    fractions = [Fraction(entry) for entry in vector]
    if all(entry == 0 for entry in fractions):
        raise ValueError("the zero vector has no direction")
    common = 1
    for entry in fractions:
        common = lcm(common, entry.denominator)
    integers = [int(entry * common) for entry in fractions]
    divisor = 0
    for entry in integers:
        divisor = gcd(divisor, entry)
    integers = [entry // divisor for entry in integers]
    first = next(entry for entry in integers if entry != 0)
    return [-entry for entry in integers] if first < 0 else integers


def mat_vec_exact(matrix: list[list[float]], vector: list[float]) -> list[Fraction]:
    """Matrix times vector with exact fractions."""
    return [sum((Fraction(a) * Fraction(b) for a, b in zip(row, vector)), Fraction(0)) for row in matrix]
