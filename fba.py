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


# ATP made per unit of flux (usual textbook yields).
ATP_PER_GLYCOLYSIS = 2.0     # per glucose through glycolysis
ATP_PER_NADH = 2.5           # per NADH through electron transport
ATP_PER_PYRUVATE = 12.5      # per pyruvate burned in the mitochondria

# O2 used per unit of flux: 1/2 per NADH, 2.5 per pyruvate burned.
# One fully respired glucose (v_e = v_m = 2) uses 2(0.5) + 2(2.5) = 6 O2.
O2_PER_NADH = 0.5
O2_PER_PYRUVATE = 2.5


def solve_cell(glucose: float, oxygen: float, mito_capacity: float) -> dict[str, float]:
    """
    Closed-form FBA for one cell: the steady-state flux vector (S v = 0) that
    makes the most ATP, given how much glucose and O2 the cell can take in and
    how much pyruvate its mitochondria can burn.

    S v = 0 forces v_e = v_m = r (pyruvate respired) and v_l = 2 v_g - r
    (pyruvate fermented). ATP = 2 v_g + 15 r, so the best choice is to take in
    all the glucose and make r as large as every limit allows:
        r = min(2 * glucose, oxygen / 3, mito_capacity)
    Raises ValueError if any input is negative.
    """
    if glucose < 0 or oxygen < 0 or mito_capacity < 0:
        raise ValueError("glucose, oxygen and mito_capacity must be >= 0")
    v_g = glucose
    oxygen_per_r = O2_PER_NADH + O2_PER_PYRUVATE   # = 3
    r = min(2 * v_g, oxygen / oxygen_per_r, mito_capacity)
    v_l = 2 * v_g - r
    return {"v_in": v_g, "v_g": v_g, "v_e": r, "v_l": v_l, "v_m": r, "v_out": v_l}


def flux_vector(fluxes: dict[str, float]) -> list[float]:
    """Put a flux dict in the column order of S, so it can be multiplied by S."""
    return [fluxes[name] for name in REACTIONS]


def atp(fluxes: dict[str, float]) -> float:
    """ATP made per unit time by a flux vector."""
    return (ATP_PER_GLYCOLYSIS * fluxes["v_g"]
            + ATP_PER_NADH * fluxes["v_e"]
            + ATP_PER_PYRUVATE * fluxes["v_m"])


def oxygen_used(fluxes: dict[str, float]) -> float:
    """O2 used per unit time by a flux vector."""
    return O2_PER_NADH * fluxes["v_e"] + O2_PER_PYRUVATE * fluxes["v_m"]
