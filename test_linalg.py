"""Tests for linalg.py. Run from this folder: python3 -m unittest -v"""
import random
import unittest
from fractions import Fraction

from fba import stoichiometric_matrix
from linalg import (gcd, integer_vector, lcm, left_null_space, mat_vec_exact, null_space, rank,
                    rref, transpose)


class TestGcd(unittest.TestCase):
    def test_textbook_cases(self) -> None:
        self.assertEqual(gcd(6, 10), 2)       # the 6x + 10y = 8 example
        self.assertEqual(gcd(4, 6), 2)        # quiz problem 2
        self.assertEqual(gcd(3, 5), 1)
        self.assertEqual(gcd(252, 198), 18)

    def test_zero_and_signs(self) -> None:
        self.assertEqual(gcd(0, 7), 7)
        self.assertEqual(gcd(7, 0), 7)
        self.assertEqual(gcd(0, 0), 0)
        self.assertEqual(gcd(-12, 18), 6)

    def test_lcm(self) -> None:
        self.assertEqual(lcm(4, 6), 12)
        self.assertEqual(lcm(1, 9), 9)


class TestRref(unittest.TestCase):
    def test_small_example_by_hand(self) -> None:
        # [1 2; 3 4] reduces to the identity (it's invertible).
        reduced, pivots = rref([[1, 2], [3, 4]])
        self.assertEqual(reduced, [[1, 0], [0, 1]])
        self.assertEqual(pivots, [0, 1])

    def test_dependent_rows(self) -> None:
        reduced, pivots = rref([[1, 2, 3], [2, 4, 6]])
        self.assertEqual(reduced, [[1, 2, 3], [0, 0, 0]])
        self.assertEqual(pivots, [0])

    def test_needs_a_row_swap(self) -> None:
        reduced, pivots = rref([[0, 1], [1, 0]])
        self.assertEqual(reduced, [[1, 0], [0, 1]])

    def test_exact_fractions(self) -> None:
        reduced, _ = rref([[3, 1]])
        self.assertEqual(reduced, [[1, Fraction(1, 3)]])

    def test_input_not_changed(self) -> None:
        a = [[2, 4], [1, 3]]
        rref(a)
        self.assertEqual(a, [[2, 4], [1, 3]])


class TestRankAndNullSpace(unittest.TestCase):
    def test_ranks(self) -> None:
        self.assertEqual(rank([[1, 0], [0, 1]]), 2)
        self.assertEqual(rank([[0, 0], [0, 0]]), 0)
        self.assertEqual(rank([[1, 2, 3], [2, 4, 6]]), 1)

    def test_null_space_by_hand(self) -> None:
        # x + 2y + 3z = 0: free y and z give (-2, 1, 0) and (-3, 0, 1).
        self.assertEqual(null_space([[1, 2, 3]]), [[-2, 1, 0], [-3, 0, 1]])

    def test_invertible_matrix_has_trivial_null_space(self) -> None:
        self.assertEqual(null_space([[1, 2], [3, 4]]), [])

    def test_rank_nullity_and_av_zero_on_random_matrices(self) -> None:
        # Rank-nullity theorem: rank + dim(null space) = number of columns.
        rng = random.Random(0)
        for _ in range(200):
            rows, cols = rng.randint(1, 5), rng.randint(1, 6)
            a = [[rng.randint(-3, 3) for _ in range(cols)] for _ in range(rows)]
            basis = null_space(a)
            self.assertEqual(rank(a) + len(basis), cols)
            for v in basis:
                self.assertEqual(mat_vec_exact(a, v), [0] * rows)
            if basis:
                self.assertEqual(rank(basis), len(basis))     # basis is independent

    def test_transpose(self) -> None:
        self.assertEqual(transpose([[1, 2, 3], [4, 5, 6]]), [[1, 4], [2, 5], [3, 6]])


class TestIntegerVector(unittest.TestCase):
    def test_clears_denominators(self) -> None:
        self.assertEqual(integer_vector([Fraction(1, 2), Fraction(1, 3)]), [3, 2])

    def test_divides_out_common_factor(self) -> None:
        # The Diophantine step from the session: (6, -10) points the same way as (3, -5).
        self.assertEqual(integer_vector([6, -10]), [3, -5])

    def test_first_nonzero_made_positive(self) -> None:
        self.assertEqual(integer_vector([0, -2, 4]), [0, 1, -2])

    def test_zero_vector_raises(self) -> None:
        with self.assertRaises(ValueError):
            integer_vector([0, 0])


class TestStoichiometricMatrixSpaces(unittest.TestCase):
    def test_rank_four_so_null_space_is_two_dimensional(self) -> None:
        s = stoichiometric_matrix()
        self.assertEqual(rank(s), 4)
        self.assertEqual(len(null_space(s)), 6 - 4)

    def test_integer_null_space_basis_is_respire_and_ferment(self) -> None:
        basis = sorted(integer_vector(v) for v in null_space(stoichiometric_matrix()))
        self.assertEqual(basis, [[1, 1, 0, 2, 0, 2], [1, 1, 2, 0, 2, 0]])

    def test_left_null_space_is_nad_conservation(self) -> None:
        # y S = 0 for y = (0, 0, 1, 1, 0): NAD+ + NADH never changes.
        self.assertEqual([integer_vector(v) for v in left_null_space(stoichiometric_matrix())],
                         [[0, 0, 1, 1, 0]])


if __name__ == "__main__":
    unittest.main()
