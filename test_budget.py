"""Tests for the shared protein budget (step 8c). Run: .venv/bin/python -m unittest -v"""
import random
import statistics
import unittest

from datatypes import HEALTHY, TUMOR, Cell, Grid, Params, empty_cell
from fba import atp, best_budget_amounts, budget_traits, fluxes_from_amounts, solve_2x2, solve_cell
from fields import make_field
from metrics import summary, tumor_shares
from rules import cell_fluxes, make_daughter, step

try:
    from scipy.optimize import linprog
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


def budget_tumor(share: float = 0.5) -> Cell:
    return Cell(TUMOR, acid_tolerance=3.0, mito_share=share)


class TestBudgetBasics(unittest.TestCase):
    def test_budget_traits(self) -> None:
        self.assertEqual(budget_traits(0.0, 20, 10), (20.0, 0.0))     # all transporters
        self.assertEqual(budget_traits(1.0, 20, 10), (0.0, 2.0))      # all mitochondria
        self.assertEqual(budget_traits(0.5, 20, 10), (10.0, 1.0))

    def test_budget_traits_spend_exactly_the_budget(self) -> None:
        for share in (0, 0.1, 0.37, 0.9, 1):
            uptake, mito = budget_traits(share, 20, 10)
            self.assertAlmostEqual(uptake + 10 * mito, 20)

    def test_bad_share_raises(self) -> None:
        with self.assertRaises(ValueError):
            budget_traits(1.2, 20, 10)
        with self.assertRaises(ValueError):
            Cell(TUMOR, mito_share=-0.1)

    def test_solve_2x2(self) -> None:
        self.assertEqual(solve_2x2(1, 1, 3, 1, -1, 1), (2.0, 1.0))   # x + y = 3, x - y = 1
        self.assertIsNone(solve_2x2(1, 2, 3, 2, 4, 5))                # parallel lines


class TestBestBudgetAmounts(unittest.TestCase):
    def check(self, args: tuple, a: float, b: float) -> None:
        got_a, got_b = best_budget_amounts(*args)
        self.assertAlmostEqual(got_a, a, msg=str(args))
        self.assertAlmostEqual(got_b, b, msg=str(args))

    def test_cases_from_the_picture(self) -> None:
        # (glucose, O2, budget, k): plenty of glucose -> pure fermentation, 20 ATP.
        self.check((10, 10, 10, 10), 10, 0)
        # Scarce glucose -> leftover budget buys mitochondria: b = (10 - 1) / 20 = 0.45, 15.5 ATP.
        self.check((1, 10, 10, 10), 0.55, 0.45)
        # No O2 -> mitochondria are worthless.
        self.check((10, 0, 10, 10), 10, 0)

    def test_simulation_budget(self) -> None:
        # B = 20: glucose 5 caps uptake, leftover 15 buys M = 1.5 -> b = 0.75, 32.5 ATP.
        self.check((5, 5, 20, 10), 4.25, 0.75)
        self.check((10, 10, 20, 10), 9.5, 0.5)

    def test_cheap_mitochondria_respire_everything(self) -> None:
        # k = 0.1: respiring one glucose costs 1.2 budget; O2 6 allows 1 glucose.
        self.check((1, 6, 20, 0.1), 0, 1)

    def test_no_glucose(self) -> None:
        self.check((0, 10, 20, 10), 0, 0)

    def test_answer_obeys_every_constraint(self) -> None:
        rng = random.Random(0)
        for _ in range(500):
            g, o, budget, k = rng.uniform(0, 12), rng.uniform(0, 12), rng.uniform(0.5, 30), rng.uniform(0.1, 20)
            a, b = best_budget_amounts(g, o, budget, k)
            self.assertGreaterEqual(a, 0)
            self.assertGreaterEqual(b, 0)
            self.assertLessEqual(a + b, g + 1e-9)
            self.assertLessEqual(6 * b, o + 1e-9)
            self.assertLessEqual(a + b + 2 * k * b, budget + 1e-9)

    def test_scales_with_budget_and_supply(self) -> None:
        # Halving glucose, O2 and budget halves the answer (what substeps rely on).
        for args in [(10, 10, 20, 10), (1, 10, 10, 10), (5, 3, 20, 4)]:
            a, b = best_budget_amounts(*args)
            ha, hb = best_budget_amounts(args[0] / 2, args[1] / 2, args[2] / 2, args[3])
            self.assertAlmostEqual(ha, a / 2)
            self.assertAlmostEqual(hb, b / 2)

    def test_negative_input_raises(self) -> None:
        with self.assertRaises(ValueError):
            best_budget_amounts(-1, 1, 1, 1)
        with self.assertRaises(ValueError):
            best_budget_amounts(1, 1, 0, 1)


@unittest.skipUnless(HAVE_SCIPY, "scipy not installed; use .venv/bin/python")
class TestBudgetMatchesLP(unittest.TestCase):
    def test_corner_method_matches_linprog(self) -> None:
        # Same problem in variables (g, b): maximize 2(g - b) + 32 b.
        rng = random.Random(1)
        for _ in range(400):
            g, o, budget, k = rng.uniform(0, 12), rng.uniform(0, 12), rng.uniform(0.5, 30), rng.uniform(0.1, 20)
            result = linprog([-2, -30], A_ub=[[1, 2 * k], [1, 0], [0, 1], [-1, 1]],
                             b_ub=[budget, g, o / 6, 0], bounds=[(0, None), (0, None)], method="highs")
            best_lp = -result.fun
            a, b = best_budget_amounts(g, o, budget, k)
            self.assertAlmostEqual(2 * a + 32 * b, best_lp, places=6, msg=f"{(g, o, budget, k)}")


class TestBudgetCellsOnTheGrid(unittest.TestCase):
    def test_fixed_trait_cells_unchanged(self) -> None:
        cell = Cell(HEALTHY, 10, 1.0)
        self.assertEqual(cell_fluxes(cell, 1, 10, Params()), solve_cell(1, 10, 10))

    def test_evolving_budget_cell_uses_its_share(self) -> None:
        # share 0.5, B 20, k 10 -> uptake 10, M 1: glucose 5 caps uptake.
        self.assertEqual(cell_fluxes(budget_tumor(0.5), 5, 10, Params()), solve_cell(5, 10, 1.0))

    def test_plastic_budget_cell_ignores_its_share(self) -> None:
        plastic = Params(plastic=True)
        best = fluxes_from_amounts(*best_budget_amounts(5, 5, 20, 10))
        for share in (0.0, 0.5, 1.0):
            self.assertEqual(cell_fluxes(budget_tumor(share), 5, 5, plastic), best)

    def test_plastic_never_worse_than_any_fixed_share(self) -> None:
        plastic = Params(plastic=True)
        for g, o in [(5, 5), (1, 10), (10, 0), (3, 8)]:
            best = atp(cell_fluxes(budget_tumor(), g, o, plastic))
            for share in [i / 20 for i in range(21)]:
                self.assertLessEqual(atp(cell_fluxes(budget_tumor(share), g, o, Params())), best + 1e-9)

    def test_share_mutates_and_stays_in_range(self) -> None:
        rng = random.Random(0)
        values = [make_daughter(budget_tumor(0.5), Params(mutation_sd=0.05), rng).mito_share for _ in range(5000)]
        self.assertAlmostEqual(statistics.mean(values), 0.5, delta=0.01)
        self.assertAlmostEqual(statistics.stdev(values), 0.05, delta=0.005)
        edge = [make_daughter(budget_tumor(0.99), Params(mutation_sd=0.5), rng).mito_share for _ in range(500)]
        self.assertTrue(all(0 <= v <= 1 for v in edge))

    def test_share_mutation_leaves_mito_capacity_alone(self) -> None:
        child = make_daughter(budget_tumor(0.5), Params(mutation_sd=0.2), random.Random(2))
        self.assertEqual(child.mito_capacity, 0.0)

    def test_budget_tumor_grows_in_a_step(self) -> None:
        # Plenty of food, acid-proof budget tumor divides (ATP 35 >= 20).
        cells = [[Cell(TUMOR, acid_tolerance=100, mito_share=0.5), empty_cell()]]
        g = Grid(cells, make_field(1, 2, 10.0), make_field(1, 2, 10.0), make_field(1, 2), [])
        new = step(g, Params(mutation_sd=0.05), random.Random(0))
        self.assertEqual(new.count(TUMOR), 2)
        self.assertNotEqual(new.cells[0][1].mito_share, 0.5)

    def test_metrics_report_shares(self) -> None:
        cells = [[budget_tumor(0.2), budget_tumor(0.6), Cell(TUMOR, 0.5, 3.0)]]
        g = Grid(cells, make_field(1, 3, 10.0), make_field(1, 3, 5.0), make_field(1, 3), [])
        self.assertEqual(tumor_shares(g), [0.2, 0.6])
        self.assertAlmostEqual(summary(g)["mean_share"], 0.4)


if __name__ == "__main__":
    unittest.main()
