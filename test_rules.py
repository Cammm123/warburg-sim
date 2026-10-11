"""Tests for datatypes.py and rules.py. Run from this folder: python3 -m unittest -v"""
import unittest

from datatypes import EMPTY, HEALTHY, TUMOR, Cell, empty_cell
from rules import DIE, DIVIDE, NOTHING, STAY, cell_fate, neighbors

SURVIVE, DIVIDE_AT = 5.0, 20.0


def healthy() -> Cell:
    return Cell(HEALTHY, mito_capacity=10, acid_tolerance=1.0)


def tumor() -> Cell:
    return Cell(TUMOR, mito_capacity=0.5, acid_tolerance=3.0)


class TestCell(unittest.TestCase):
    def test_fields(self) -> None:
        c = tumor()
        self.assertEqual((c.kind, c.mito_capacity, c.acid_tolerance), (TUMOR, 0.5, 3.0))

    def test_empty_cell(self) -> None:
        self.assertEqual(empty_cell().kind, EMPTY)

    def test_bad_kind_raises(self) -> None:
        with self.assertRaises(ValueError):
            Cell("zombie")


class TestNeighbors(unittest.TestCase):
    def test_middle_has_eight_in_reading_order(self) -> None:
        self.assertEqual(neighbors(1, 1, 3, 3),
                         [(0, 0), (0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1), (2, 2)])

    def test_corner_has_three(self) -> None:
        self.assertEqual(neighbors(0, 0, 3, 3), [(0, 1), (1, 0), (1, 1)])

    def test_other_corner(self) -> None:
        self.assertEqual(neighbors(2, 3, 3, 4), [(1, 2), (1, 3), (2, 2)])

    def test_edge_has_five(self) -> None:
        self.assertEqual(len(neighbors(0, 2, 5, 5)), 5)

    def test_one_by_one_has_none(self) -> None:
        self.assertEqual(neighbors(0, 0, 1, 1), [])

    def test_single_row(self) -> None:
        self.assertEqual(neighbors(0, 1, 1, 3), [(0, 0), (0, 2)])

    def test_never_includes_itself_and_all_on_grid(self) -> None:
        for r in range(4):
            for c in range(5):
                for nr, nc in neighbors(r, c, 4, 5):
                    self.assertNotEqual((nr, nc), (r, c))
                    self.assertTrue(0 <= nr < 4 and 0 <= nc < 5)

    def test_neighbor_relation_is_symmetric(self) -> None:
        # If b is next to a, then a is next to b.
        for r in range(4):
            for c in range(5):
                for nr, nc in neighbors(r, c, 4, 5):
                    self.assertIn((r, c), neighbors(nr, nc, 4, 5))


class TestCellFate(unittest.TestCase):
    def fate(self, cell: Cell, cell_atp: float, lactate: float) -> str:
        return cell_fate(cell, cell_atp, lactate, SURVIVE, DIVIDE_AT)

    def test_empty_square_does_nothing(self) -> None:
        self.assertEqual(self.fate(empty_cell(), 100, 100), NOTHING)

    def test_acid_threshold_at_above_below(self) -> None:
        # healthy tolerance 1.0; ATP 10 is fine on its own (stay).
        self.assertEqual(self.fate(healthy(), 10, 0.99), STAY)
        self.assertEqual(self.fate(healthy(), 10, 1.0), STAY)    # exactly at: survives
        self.assertEqual(self.fate(healthy(), 10, 1.01), DIE)

    def test_acid_kills_even_a_cell_with_lots_of_atp(self) -> None:
        self.assertEqual(self.fate(healthy(), 1000, 2.0), DIE)

    def test_tumor_survives_acid_that_kills_healthy(self) -> None:
        self.assertEqual(self.fate(healthy(), 10, 2.0), DIE)
        self.assertEqual(self.fate(tumor(), 10, 2.0), STAY)

    def test_tumor_dies_above_its_own_tolerance(self) -> None:
        self.assertEqual(self.fate(tumor(), 10, 3.5), DIE)

    def test_starvation_threshold_at_above_below(self) -> None:
        self.assertEqual(self.fate(healthy(), 4.99, 0), DIE)
        self.assertEqual(self.fate(healthy(), 5.0, 0), STAY)     # exactly at: survives
        self.assertEqual(self.fate(healthy(), 5.01, 0), STAY)

    def test_divide_threshold_at_above_below(self) -> None:
        self.assertEqual(self.fate(tumor(), 19.99, 0), STAY)
        self.assertEqual(self.fate(tumor(), 20.0, 0), DIVIDE)    # exactly at: divides
        self.assertEqual(self.fate(tumor(), 32, 0), DIVIDE)

    def test_zero_atp_dies(self) -> None:
        self.assertEqual(self.fate(tumor(), 0, 0), DIE)

    def test_equal_thresholds_allowed(self) -> None:
        self.assertEqual(cell_fate(healthy(), 5, 0, 5, 5), DIVIDE)

    def test_divide_below_survive_raises(self) -> None:
        with self.assertRaises(ValueError):
            cell_fate(healthy(), 10, 0, survive_atp=20, divide_atp=5)


if __name__ == "__main__":
    unittest.main()
