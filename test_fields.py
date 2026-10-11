"""Tests for fields.py. Run from this folder: python3 -m unittest -v"""
import unittest

from fields import apply_change, diffuse, make_field, refill_sources, total


def spike(num_rows: int, num_cols: int, row: int, col: int, amount: float = 1.0) -> list[list[float]]:
    """A field that is zero everywhere except one square."""
    field = make_field(num_rows, num_cols)
    field[row][col] = amount
    return field


class TestMakeField(unittest.TestCase):
    def test_shape_and_value(self) -> None:
        f = make_field(2, 3, 7.0)
        self.assertEqual(f, [[7.0, 7.0, 7.0], [7.0, 7.0, 7.0]])

    def test_rows_are_independent(self) -> None:
        # Changing one row must not change another (a classic list bug).
        f = make_field(2, 2)
        f[0][0] = 5
        self.assertEqual(f[1][0], 0)


class TestDiffuse(unittest.TestCase):
    def test_spike_in_middle_matches_kernel(self) -> None:
        # One step from a single 1.0: the neighbors get rate x kernel weight.
        rate = 0.5
        new = diffuse(spike(5, 5, 2, 2), rate)
        self.assertAlmostEqual(new[2][2], 1 - rate)          # kept
        for r, c in [(1, 2), (3, 2), (2, 1), (2, 3)]:
            self.assertAlmostEqual(new[r][c], 0.2 * rate)    # sides
        for r, c in [(1, 1), (1, 3), (3, 1), (3, 3)]:
            self.assertAlmostEqual(new[r][c], 0.05 * rate)   # corners
        self.assertAlmostEqual(new[0][0], 0)                 # too far away

    def test_spreads_symmetrically(self) -> None:
        f = spike(7, 7, 3, 3)
        for _ in range(5):
            f = diffuse(f, 0.8)
        for r in range(7):
            for c in range(7):
                self.assertAlmostEqual(f[r][c], f[6 - r][c])   # top/bottom mirror
                self.assertAlmostEqual(f[r][c], f[r][6 - c])   # left/right mirror
                self.assertAlmostEqual(f[r][c], f[c][r])       # diagonal mirror

    def test_total_conserved_in_middle_edge_and_corner(self) -> None:
        for row, col in [(2, 2), (0, 2), (0, 0), (4, 4)]:
            f = spike(5, 5, row, col, 3.0)
            for _ in range(20):
                f = diffuse(f, 1.0)
            self.assertAlmostEqual(total(f), 3.0)

    def test_uniform_field_stays_uniform(self) -> None:
        # Nothing to even out: every square sends and receives the same.
        new = diffuse(make_field(4, 6, 2.5), 0.7)
        for row in new:
            for value in row:
                self.assertAlmostEqual(value, 2.5)

    def test_rate_zero_changes_nothing(self) -> None:
        f = spike(3, 3, 1, 1)
        self.assertEqual(diffuse(f, 0), f)

    def test_never_negative(self) -> None:
        f = spike(5, 5, 0, 0)
        for _ in range(30):
            f = diffuse(f, 1.0)
            for row in f:
                for value in row:
                    self.assertGreaterEqual(value, 0)

    def test_rounding_crumbs_snap_to_zero(self) -> None:
        # Found by the calibration run: at rate 1 the middle square sends away
        # everything, and float rounding used to leave about -3.5e-18 behind.
        f = [[9.196592523769888e-31, 0.0, 0.0],
             [5.067951222763686e-31, 0.26869184399515444, 4.902799617558136e-31],
             [2.653069509641705e-31, 0.0, 4.569209446671899e-31]]
        new = diffuse(f, 1.0)
        for row in new:
            for value in row:
                self.assertGreaterEqual(value, 0.0)
        self.assertAlmostEqual(total(new), total(f))

    def test_evens_out_over_time(self) -> None:
        # After many steps, the spike becomes nearly uniform: 1/25 each.
        f = spike(5, 5, 0, 0)
        for _ in range(500):
            f = diffuse(f, 1.0)
        for row in f:
            for value in row:
                self.assertAlmostEqual(value, 1 / 25, places=4)

    def test_one_by_one_grid(self) -> None:
        self.assertEqual(diffuse([[4.0]], 1.0), [[4.0]])

    def test_input_not_changed(self) -> None:
        f = spike(3, 3, 1, 1)
        diffuse(f, 1.0)
        self.assertEqual(f, spike(3, 3, 1, 1))

    def test_bad_rate_raises(self) -> None:
        for rate in [-0.1, 1.1]:
            with self.assertRaises(ValueError):
                diffuse(make_field(2, 2), rate)


class TestRefillSources(unittest.TestCase):
    def test_sets_only_the_sources(self) -> None:
        f = make_field(2, 3, 0.5)
        new = refill_sources(f, [(0, 0), (1, 2)], 10.0)
        self.assertEqual(new, [[10.0, 0.5, 0.5], [0.5, 0.5, 10.0]])
        self.assertEqual(f, make_field(2, 3, 0.5))  # input not changed

    def test_no_sources(self) -> None:
        f = make_field(2, 2, 1.0)
        self.assertEqual(refill_sources(f, [], 9.0), f)


class TestApplyChange(unittest.TestCase):
    def test_use_and_release(self) -> None:
        new = apply_change([[5.0, 1.0]], [[-2.0, 3.0]])
        self.assertEqual(new, [[3.0, 4.0]])

    def test_use_exactly_everything(self) -> None:
        self.assertEqual(apply_change([[2.0]], [[-2.0]]), [[0.0]])

    def test_tiny_rounding_negative_becomes_zero(self) -> None:
        # 0.3 - (0.1 + 0.2) is about -5.6e-17 in floating point, not exactly 0.
        self.assertLess(0.3 - (0.1 + 0.2), 0)
        self.assertEqual(apply_change([[0.3]], [[-(0.1 + 0.2)]]), [[0.0]])

    def test_using_more_than_there_raises(self) -> None:
        with self.assertRaises(ValueError):
            apply_change([[1.0]], [[-1.5]])

    def test_shape_mismatch_raises(self) -> None:
        with self.assertRaises(ValueError):
            apply_change([[1.0, 2.0]], [[1.0]])


if __name__ == "__main__":
    unittest.main()
