"""
The same one-cell FBA as fba.solve_cell, but solved the textbook way: hand the
whole problem (S v = 0, bounds, maximize ATP) to a general linear-programming
solver. It knows nothing about our closed-form shortcut, so agreement between
the two is an independent check of the shortcut. Used only by the tests.

Needs scipy (installed in .venv/).
"""
from scipy.optimize import linprog

from fba import (ATP_PER_GLYCOLYSIS, ATP_PER_NADH, ATP_PER_PYRUVATE,
                 O2_PER_NADH, O2_PER_PYRUVATE, REACTIONS, stoichiometric_matrix)


def solve_cell_lp(glucose: float, oxygen: float, mito_capacity: float) -> dict[str, float]:
    """
    Maximize ATP subject to:
        S v = 0                                   (steady state)
        0.5 v_e + 2.5 v_m <= oxygen               (O2 limit)
        0 <= v_in <= glucose, 0 <= v_m <= mito_capacity, every other v >= 0
    Returns the optimal fluxes by reaction name. Raises RuntimeError if the
    solver fails.
    """
    # linprog minimizes, so minimize -ATP to maximize ATP.
    atp_weights = {"v_g": ATP_PER_GLYCOLYSIS, "v_e": ATP_PER_NADH, "v_m": ATP_PER_PYRUVATE}
    cost = [-atp_weights.get(name, 0.0) for name in REACTIONS]

    s = stoichiometric_matrix()
    zeros = [0.0] * len(s)

    o2_row = [0.0] * len(REACTIONS)
    o2_row[REACTIONS.index("v_e")] = O2_PER_NADH
    o2_row[REACTIONS.index("v_m")] = O2_PER_PYRUVATE

    bounds = []
    for name in REACTIONS:
        if name == "v_in":
            bounds.append((0, glucose))
        elif name == "v_m":
            bounds.append((0, mito_capacity))
        else:
            bounds.append((0, None))

    result = linprog(cost, A_ub=[o2_row], b_ub=[oxygen], A_eq=s, b_eq=zeros,
                     bounds=bounds, method="highs")
    if not result.success:
        raise RuntimeError(f"LP solver failed: {result.message}")
    return {name: float(value) for name, value in zip(REACTIONS, result.x)}
