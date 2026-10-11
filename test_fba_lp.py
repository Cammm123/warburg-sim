"""
LP cross-check: the closed-form solve_cell must match a general LP solver.
Run with the venv's Python (it has scipy):  .venv/bin/python -m unittest -v
With plain python3 these tests are skipped, not failed.
"""
import unittest

from fba import atp, solve_cell

try:
    from fba_lp import solve_cell_lp
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

VALUES = [0, 0.1, 0.25, 0.5, 1, 1.5, 2, 3, 4.5, 6, 10]


@unittest.skipUnless(HAVE_SCIPY, "scipy not installed; use .venv/bin/python")
class TestClosedFormMatchesLP(unittest.TestCase):
    def test_hand_worked_cases(self) -> None:
        # No O2, plenty of everything, and the Warburg case from test_fba.py.
        for g, o, m in [(1, 0, 10), (1, 100, 100), (1, 100, 0.5)]:
            lp, closed = solve_cell_lp(g, o, m), solve_cell(g, o, m)
            for name in closed:
                self.assertAlmostEqual(lp[name], closed[name], places=6, msg=f"{name} at {(g, o, m)}")

    def test_same_best_atp_everywhere(self) -> None:
        # 11^3 = 1331 inputs: the best ATP must agree.
        for g in VALUES:
            for o in VALUES:
                for m in VALUES:
                    self.assertAlmostEqual(atp(solve_cell_lp(g, o, m)), atp(solve_cell(g, o, m)),
                                           places=6, msg=f"at glucose={g}, oxygen={o}, mito={m}")

    def test_same_fluxes_everywhere(self) -> None:
        # The best flux vector is unique here, so every flux must agree too.
        for g in VALUES:
            for o in VALUES:
                for m in VALUES:
                    lp, closed = solve_cell_lp(g, o, m), solve_cell(g, o, m)
                    for name in closed:
                        self.assertAlmostEqual(lp[name], closed[name], places=6,
                                               msg=f"{name} at glucose={g}, oxygen={o}, mito={m}")


if __name__ == "__main__":
    unittest.main()
