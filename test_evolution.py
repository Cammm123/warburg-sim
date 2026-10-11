"""Tests for mutation (rules.make_daughter) and metrics.py. Run: python3 -m unittest -v"""
import random
import statistics
import unittest

from datatypes import HEALTHY, TUMOR, Cell, Grid, Params, empty_cell
from fields import make_field
from metrics import mean, summary, tumor_traits, tumor_warburg_indices
from rules import make_daughter, step


def tumor(mito: float = 0.5) -> Cell:
    return Cell(TUMOR, mito_capacity=mito, acid_tolerance=3.0, glucose_uptake=8.0)


class TestMakeDaughter(unittest.TestCase):
    def test_healthy_copies_exactly_even_with_mutation_on(self) -> None:
        parent = Cell(HEALTHY, 10, 1.0)
        rng = random.Random(0)
        for _ in range(100):
            self.assertEqual(make_daughter(parent, Params(mutation_sd=0.5), rng).mito_capacity, 10)

    def test_tumor_copies_exactly_when_mutation_off(self) -> None:
        self.assertEqual(make_daughter(tumor(0.7), Params(), random.Random(0)).mito_capacity, 0.7)

    def test_tumor_mutation_has_right_mean_and_spread(self) -> None:
        rng = random.Random(0)
        values = [make_daughter(tumor(2.0), Params(mutation_sd=0.1), rng).mito_capacity for _ in range(5000)]
        self.assertAlmostEqual(statistics.mean(values), 2.0, delta=0.01)
        self.assertAlmostEqual(statistics.stdev(values), 0.1, delta=0.01)

    def test_never_negative(self) -> None:
        rng = random.Random(0)
        for _ in range(2000):
            self.assertGreaterEqual(make_daughter(tumor(0.05), Params(mutation_sd=1.0), rng).mito_capacity, 0.0)

    def test_other_traits_and_parent_unchanged(self) -> None:
        parent = tumor(0.5)
        child = make_daughter(parent, Params(mutation_sd=0.3), random.Random(1))
        self.assertEqual((child.kind, child.acid_tolerance, child.glucose_uptake), (TUMOR, 3.0, 8.0))
        self.assertEqual(parent.mito_capacity, 0.5)
        self.assertIsNot(child, parent)

    def test_same_seed_same_daughter(self) -> None:
        a = make_daughter(tumor(), Params(mutation_sd=0.1), random.Random(7)).mito_capacity
        b = make_daughter(tumor(), Params(mutation_sd=0.1), random.Random(7)).mito_capacity
        self.assertEqual(a, b)

    def test_negative_sd_raises(self) -> None:
        with self.assertRaises(ValueError):
            Params(mutation_sd=-0.1)


class TestMutationInStep(unittest.TestCase):
    def test_dividing_tumor_gets_a_mutated_daughter(self) -> None:
        # Well-fed tumor in a 1x2 grid divides into the empty square. Acid
        # tolerance is high because in a closed box its 14 lactate would kill it.
        cells = [[Cell(TUMOR, mito_capacity=2.0, acid_tolerance=100.0, glucose_uptake=8.0), empty_cell()]]
        g = Grid(cells, make_field(1, 2, 100.0), make_field(1, 2, 100.0), make_field(1, 2), [])
        new = step(g, Params(mutation_sd=0.2), random.Random(3))
        self.assertEqual([c.kind for c in new.cells[0]], [TUMOR, TUMOR])
        self.assertEqual(new.cells[0][0].mito_capacity, 2.0)          # parent unchanged
        self.assertNotEqual(new.cells[0][1].mito_capacity, 2.0)       # daughter mutated


class TestMetrics(unittest.TestCase):
    def grid(self) -> Grid:
        cells = [[tumor(0.5), Cell(HEALTHY, 10, 1.0), tumor(10.0)]]
        g = Grid(cells, make_field(1, 3, 10.0), make_field(1, 3, 1.0), make_field(1, 3), [])
        return g

    def test_tumor_traits_ignore_healthy(self) -> None:
        self.assertEqual(tumor_traits(self.grid()), [0.5, 10.0])

    def test_warburg_indices_by_hand(self) -> None:
        # Square glucose 1 caps uptake 8 at 1. M 0.5 -> b = 0.25, index 0.75.
        # M 10 with O2 10 -> b = min(1, 10/6, 5) = 1, index 0.
        indices = tumor_warburg_indices(self.grid())
        self.assertAlmostEqual(indices[0], 0.75)
        self.assertAlmostEqual(indices[1], 0.0)

    def test_summary(self) -> None:
        s = summary(self.grid())
        self.assertAlmostEqual(s["mean_mito"], 5.25)
        self.assertAlmostEqual(s["mean_warburg"], 0.375)

    def test_no_tumor_gives_zeros(self) -> None:
        g = Grid([[empty_cell()]], [[0.0]], [[0.0]], [[0.0]], [])
        self.assertEqual(summary(g), {"mean_mito": 0.0, "mean_warburg": 0.0})
        self.assertEqual(mean([]), 0.0)


if __name__ == "__main__":
    unittest.main()
