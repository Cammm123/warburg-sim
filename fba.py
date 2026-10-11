"""
The cell's chemistry as linear algebra.

Rows of S are metabolites, columns are reactions. Entry S[i][j] is how many of
metabolite i reaction j makes (+) or uses up (-). For a flux vector v (one rate
per reaction), S times v is how fast each metabolite's amount is changing, so
steady state means S v = 0.

The steady states are the null space of S. linalg.py computes it: it is
2-dimensional, and its integer basis is two pathways, RESPIRE and FERMENT, each
processing one glucose. Every cell's flux vector is a FERMENT + b RESPIRE.
"""
from linalg import integer_vector, null_space

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


def steady_state_pathways() -> tuple[list[int], list[int]]:
    """
    Compute the null space of S and return its integer basis as
    (respire, ferment). Respire is the one that uses electron transport.
    Raises ValueError if the null space isn't the expected two pathways.
    """
    basis = [integer_vector(v) for v in null_space(stoichiometric_matrix())]
    if len(basis) != 2:
        raise ValueError(f"expected a 2-dimensional null space, got {len(basis)}")
    v_e = REACTIONS.index("v_e")
    respire = next(v for v in basis if v[v_e] > 0)
    ferment = next(v for v in basis if v[v_e] == 0)
    return respire, ferment


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


RESPIRE, FERMENT = steady_state_pathways()

# What one glucose through each pathway costs and makes, read off the vectors.
O2_PER_RESPIRED_GLUCOSE = (O2_PER_NADH * RESPIRE[REACTIONS.index("v_e")]
                           + O2_PER_PYRUVATE * RESPIRE[REACTIONS.index("v_m")])   # 6
PYRUVATE_BURNED_PER_RESPIRED_GLUCOSE = RESPIRE[REACTIONS.index("v_m")]           # 2


def pathway_amounts(glucose: float, oxygen: float, mito_capacity: float) -> tuple[float, float]:
    """
    The best steady state in null-space coordinates: (a, b) = (glucose
    fermented, glucose respired). Respiring makes 32 ATP per glucose and
    fermenting 2, so respire as much as every limit allows and ferment the rest:
        b = min(glucose, oxygen / 6, mito_capacity / 2),   a = glucose - b
    Raises ValueError if any input is negative.
    """
    if glucose < 0 or oxygen < 0 or mito_capacity < 0:
        raise ValueError("glucose, oxygen and mito_capacity must be >= 0")
    b = min(glucose, oxygen / O2_PER_RESPIRED_GLUCOSE,
            mito_capacity / PYRUVATE_BURNED_PER_RESPIRED_GLUCOSE)
    return glucose - b, b


def fluxes_from_amounts(a: float, b: float) -> dict[str, float]:
    """The flux vector a FERMENT + b RESPIRE, by reaction name."""
    return {name: a * f + b * r for name, f, r in zip(REACTIONS, FERMENT, RESPIRE)}


def solve_cell(glucose: float, oxygen: float, mito_capacity: float) -> dict[str, float]:
    """
    Closed-form FBA for one cell: the steady-state flux vector (S v = 0) that
    makes the most ATP, given how much glucose and O2 the cell can take in and
    how much pyruvate its mitochondria can burn. Built from the null space:
        v = a FERMENT + b RESPIRE, with (a, b) from pathway_amounts.
    Raises ValueError if any input is negative.
    """
    a, b = pathway_amounts(glucose, oxygen, mito_capacity)
    return fluxes_from_amounts(a, b)


def warburg_index(fluxes: dict[str, float]) -> float:
    """
    Share of glucose that is fermented, a / (a + b): 0 = fully respiring,
    1 = fully fermenting. A cell taking in no glucose gets 0.
    """
    a = fluxes["v_l"] / FERMENT[REACTIONS.index("v_l")]
    b = fluxes["v_m"] / RESPIRE[REACTIONS.index("v_m")]
    return a / (a + b) if a + b > 0 else 0.0


# ---- The shared protein budget (step 8c) -----------------------------------
# A budget cell splits a budget B between transporters (uptake u, cost 1 per
# unit) and mitochondria (capacity M, cost mito_cost = k per unit):
#     u + k M = B
# In flux terms, taking in g glucose needs u >= g, and respiring b glucose needs
# M >= 2b, so the budget is one more linear inequality:  g + 2k b <= B.

ATP_FERMENT = 2.0    # ATP per glucose through FERMENT (= atp of the FERMENT vector)
ATP_RESPIRE = 32.0   # ATP per glucose through RESPIRE


def budget_traits(mito_share: float, budget: float, mito_cost: float) -> tuple[float, float]:
    """(glucose_uptake, mito_capacity) for a cell spending mito_share of its budget on mitochondria."""
    if not 0 <= mito_share <= 1:
        raise ValueError("mito_share must be between 0 and 1")
    return (1 - mito_share) * budget, mito_share * budget / mito_cost


def solve_2x2(p1: float, q1: float, r1: float, p2: float, q2: float, r2: float) -> tuple[float, float] | None:
    """
    Solve p1 x + q1 y = r1 and p2 x + q2 y = r2 by Cramer's rule.
    Returns None when the lines are parallel (determinant 0).
    """
    det = p1 * q2 - p2 * q1
    if abs(det) < 1e-12:
        return None
    return (r1 * q2 - r2 * q1) / det, (p1 * r2 - p2 * r1) / det


def best_budget_amounts(glucose: float, oxygen: float, budget: float,
                        mito_cost: float) -> tuple[float, float]:
    """
    A plastic cell's best split, as null-space amounts (a, b) = (fermented, respired).
    Choose glucose taken in g and respired b to maximize ATP = 2(g - b) + 32 b
    subject to:   g + 2k b <= budget,  g <= glucose,  b <= oxygen / 6,  b <= g,
                  g >= 0,  b >= 0.
    A linear objective on a polygon is largest at a corner, so: intersect every
    pair of constraint lines, keep the corners that satisfy every constraint,
    and return the corner with the most ATP.
    Raises ValueError if any input is negative or budget/mito_cost isn't positive.
    """
    if glucose < 0 or oxygen < 0:
        raise ValueError("glucose and oxygen must be >= 0")
    if budget <= 0 or mito_cost <= 0:
        raise ValueError("budget and mito_cost must be positive")
    # Each constraint p g + q b <= r.
    # Respiring one glucose needs PYRUVATE_BURNED_PER_RESPIRED_GLUCOSE (= 2) units
    # of M, each costing mito_cost, hence the 2k in g + 2k b <= budget.
    constraints = [(1.0, mito_cost * PYRUVATE_BURNED_PER_RESPIRED_GLUCOSE, budget),
                   (1.0, 0.0, glucose),
                   (0.0, 1.0, oxygen / O2_PER_RESPIRED_GLUCOSE),
                   (-1.0, 1.0, 0.0),
                   (-1.0, 0.0, 0.0),
                   (0.0, -1.0, 0.0)]
    best_g, best_b, best_atp = 0.0, 0.0, 0.0
    for i in range(len(constraints)):
        for j in range(i + 1, len(constraints)):
            corner = solve_2x2(*constraints[i], *constraints[j])
            if corner is None:
                continue
            g, b = corner
            if all(p * g + q * b <= r + 1e-9 for p, q, r in constraints):
                value = ATP_FERMENT * (g - b) + ATP_RESPIRE * b
                if value > best_atp + 1e-12:
                    best_g, best_b, best_atp = g, b, value
    a, b = best_g - best_b, best_b
    return (a if a > 0 else 0.0), (b if b > 0 else 0.0)   # also turns -0.0 into 0.0


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
