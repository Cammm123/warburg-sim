"""
The cell's chemistry as linear algebra.

Rows of S are metabolites, columns are reactions. Entry S[i][j] is how many of
metabolite i reaction j makes (+) or uses up (-). For a flux vector v (one rate
per reaction), S times v is how fast each metabolite's amount is changing, so
steady state means S v = 0.
"""

METABOLITES = ["glucose", "pyruvate", "NAD+", "NADH", "lactate"]

# v_in:  glucose enters the cell
# v_g:   glycolysis, glucose + 2 NAD+ -> 2 pyruvate + 2 NADH (+ 2 ATP)
# v_e:   electron transport (needs O2), NADH -> NAD+
# v_l:   lactate fermentation, pyruvate + NADH -> lactate + NAD+
# v_m:   pyruvate goes into the mitochondria and is burned
# v_out: lactate leaves the cell
REACTIONS = ["v_in", "v_g", "v_e", "v_l", "v_m", "v_out"]


def stoichiometric_matrix() -> list[list[int]]:
    """Return S (5 metabolites x 6 reactions) for the simplified cell."""
    return [
        # v_in  v_g  v_e  v_l  v_m  v_out
        [1, -1, 0, 0, 0, 0],    # glucose
        [0, 2, 0, -1, -1, 0],   # pyruvate
        [0, -2, 1, 1, 0, 0],    # NAD+
        [0, 2, -1, -1, 0, 0],   # NADH
        [0, 0, 0, 1, 0, -1],    # lactate
    ]


def mat_vec(matrix: list[list[float]], vector: list[float]) -> list[float]:
    """
    Multiply a matrix by a vector: entry i is row i dotted with the vector.
    Raises ValueError if a row's length doesn't match the vector's length.
    """
    result = []
    for row in matrix:
        if len(row) != len(vector):
            raise ValueError(f"row has {len(row)} entries but vector has {len(vector)}")
        total = 0
        for entry, value in zip(row, vector):
            total += entry * value
        result.append(total)
    return result
