"""Tests for fba.py. Run from this folder: python3 -m unittest -v"""
import unittest

from fba import METABOLITES, REACTIONS, mat_vec, stoichiometric_matrix


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


if __name__ == "__main__":
    unittest.main()
