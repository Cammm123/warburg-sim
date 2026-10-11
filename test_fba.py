"""Tests for fba.py. Run from this folder: python3 -m unittest -v"""
import unittest

from fba import (METABOLITES, REACTIONS, atp, flux_vector, mat_vec,
                 oxygen_used, solve_cell, stoichiometric_matrix)


class TestMatVec(unittest.TestCase):
    def test_small_example_by_hand(self) -> None:
        # [1 2; 3 4] times (5, 6) = (1*5 + 2*6, 3*5 + 4*6) = (17, 39)
        self.assertEqual(mat_vec([[1, 2], [3, 4]], [5, 6]), [17, 39])

    def test_identity_returns_same_vector(self) -> None:
        identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        self.assertEqual(mat_vec(identity, [7, -2, 0.5]), [7, -2, 0.5])

    def test_zero_vector_gives_zero(self) -> None:
        self.assertEqual(mat_vec([[1, 2], [3, 4]], [0, 0]), [0, 0])

    def test_unit_vector_picks_out_a_column(self) -> None:
        # The session's point: S times (0, 1) is just S's second column.
        self.assertEqual(mat_vec([[1, -1], [6, -6]], [0, 1]), [-1, -6])

    def test_non_square_matrix(self) -> None:
        self.assertEqual(mat_vec([[1, 0, 2]], [3, 4, 5]), [13])

    def test_length_mismatch_raises(self) -> None:
        with self.assertRaises(ValueError):
            mat_vec([[1, 2, 3]], [1, 2])


class TestStoichiometricMatrix(unittest.TestCase):
    def test_shape(self) -> None:
        s = stoichiometric_matrix()
        self.assertEqual(len(s), len(METABOLITES))
        for row in s:
            self.assertEqual(len(row), len(REACTIONS))

    def test_no_oxygen_steady_state(self) -> None:
        # Worked in the session: no O2, one glucose -> 2 lactate.
        v = [1, 1, 0, 2, 0, 2]  # v_in, v_g, v_e, v_l, v_m, v_out
        self.assertEqual(mat_vec(stoichiometric_matrix(), v), [0, 0, 0, 0, 0])

    def test_with_oxygen_steady_state(self) -> None:
        # Worked in the session: O2 available, no fermentation.
        v = [1, 1, 2, 0, 2, 0]
        self.assertEqual(mat_vec(stoichiometric_matrix(), v), [0, 0, 0, 0, 0])

    def test_glycolysis_alone_is_not_steady(self) -> None:
        # Glycolysis with nothing to clear its products: pyruvate and NADH
        # pile up (+2 each) and NAD+ runs out (-2).
        v = [1, 1, 0, 0, 0, 0]
        self.assertEqual(mat_vec(stoichiometric_matrix(), v), [0, 2, -2, 2, 0])

    def test_every_reaction_conserves_the_nad_pool(self) -> None:
        # NAD+ row + NADH row = 0 in every column: reactions swap one form for
        # the other, so the total NAD+ + NADH never changes.
        s = stoichiometric_matrix()
        nad, nadh = METABOLITES.index("NAD+"), METABOLITES.index("NADH")
        for j in range(len(REACTIONS)):
            self.assertEqual(s[nad][j] + s[nadh][j], 0)


class TestSolveCell(unittest.TestCase):
    def assert_fluxes(self, fluxes: dict[str, float], expected: list[float]) -> None:
        """Compare against [v_in, v_g, v_e, v_l, v_m, v_out]."""
        for got, want in zip(flux_vector(fluxes), expected):
            self.assertAlmostEqual(got, want)

    def test_no_oxygen_ferments_everything(self) -> None:
        # Same answer as the session's hand calculation: 2 lactate, 2 ATP.
        f = solve_cell(glucose=1, oxygen=0, mito_capacity=10)
        self.assert_fluxes(f, [1, 1, 0, 2, 0, 2])
        self.assertAlmostEqual(atp(f), 2)

    def test_plenty_of_everything_respires_fully(self) -> None:
        # No lactate, 32 ATP per glucose, 6 O2 per glucose.
        f = solve_cell(glucose=1, oxygen=100, mito_capacity=100)
        self.assert_fluxes(f, [1, 1, 2, 0, 2, 0])
        self.assertAlmostEqual(atp(f), 32)
        self.assertAlmostEqual(oxygen_used(f), 6)

    def test_small_mito_capacity_ferments_despite_oxygen(self) -> None:
        # The Warburg case: lots of O2, but mitochondria can only burn 0.5.
        f = solve_cell(glucose=1, oxygen=100, mito_capacity=0.5)
        self.assert_fluxes(f, [1, 1, 0.5, 1.5, 0.5, 1.5])
        self.assertAlmostEqual(atp(f), 2 + 15 * 0.5)  # 9.5

    def test_oxygen_is_the_limit(self) -> None:
        # O2 = 3 allows r = 1: half the pyruvate respired, half fermented.
        f = solve_cell(glucose=1, oxygen=3, mito_capacity=100)
        self.assert_fluxes(f, [1, 1, 1, 1, 1, 1])
        self.assertAlmostEqual(oxygen_used(f), 3)

    def test_exactly_enough_oxygen(self) -> None:
        # Boundary: O2 = 6 is exactly enough for one glucose.
        f = solve_cell(glucose=1, oxygen=6, mito_capacity=100)
        self.assertAlmostEqual(f["v_l"], 0)
        self.assertAlmostEqual(f["v_e"], 2)

    def test_no_glucose_does_nothing(self) -> None:
        f = solve_cell(glucose=0, oxygen=100, mito_capacity=100)
        self.assert_fluxes(f, [0, 0, 0, 0, 0, 0])
        self.assertAlmostEqual(atp(f), 0)

    def test_scales_with_glucose(self) -> None:
        f = solve_cell(glucose=3, oxygen=0, mito_capacity=10)
        self.assert_fluxes(f, [3, 3, 0, 6, 0, 6])

    def test_negative_input_raises(self) -> None:
        for args in [(-1, 1, 1), (1, -1, 1), (1, 1, -1)]:
            with self.assertRaises(ValueError):
                solve_cell(*args)

    def test_every_answer_is_a_valid_steady_state(self) -> None:
        # Over many inputs: S v = 0, no negative flux, and every bound holds.
        s = stoichiometric_matrix()
        values = [0, 0.25, 0.5, 1, 2, 3, 6, 10]
        for g in values:
            for o in values:
                for m in values:
                    f = solve_cell(g, o, m)
                    for entry in mat_vec(s, flux_vector(f)):
                        self.assertAlmostEqual(entry, 0)
                    for value in f.values():
                        self.assertGreaterEqual(value, 0)
                    self.assertLessEqual(f["v_g"], g + 1e-12)
                    self.assertLessEqual(oxygen_used(f), o + 1e-12)
                    self.assertLessEqual(f["v_m"], m + 1e-12)

    def test_more_oxygen_or_capacity_never_lowers_atp(self) -> None:
        # Loosening a limit can only help the best answer.
        values = [0, 0.5, 1, 3, 6, 10]
        for g in [0.5, 1, 2]:
            for i in range(len(values) - 1):
                low, high = values[i], values[i + 1]
                self.assertLessEqual(atp(solve_cell(g, low, 5)), atp(solve_cell(g, high, 5)))
                self.assertLessEqual(atp(solve_cell(g, 10, low)), atp(solve_cell(g, 10, high)))


if __name__ == "__main__":
    unittest.main()
