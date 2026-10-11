"""Tests for Grid, metabolism, choose_divisions and step. Run: python3 -m unittest -v"""
import random
import unittest

from datatypes import EMPTY, HEALTHY, TUMOR, Cell, Grid, Params, daughter_of, empty_cell
from fba import solve_cell
from fields import make_field, total
from rules import capillaries, choose_divisions, metabolism, neighbors, step


def blank_grid(num_rows: int, num_cols: int, vessels: list[tuple[int, int]] | None = None) -> Grid:
    """An empty grid with every field at 0."""
    return Grid([[empty_cell() for _ in range(num_cols)] for _ in range(num_rows)],
                make_field(num_rows, num_cols), make_field(num_rows, num_cols),
                make_field(num_rows, num_cols), vessels or [])


def healthy() -> Cell:
    return Cell(HEALTHY, mito_capacity=10, acid_tolerance=1.0)


def tumor(mito: float = 0.5, acid_tol: float = 3.0) -> Cell:
    return Cell(TUMOR, mito_capacity=mito, acid_tolerance=acid_tol)


class TestGrid(unittest.TestCase):
    def test_count(self) -> None:
        g = blank_grid(2, 2)
        g.cells[0][0], g.cells[1][1] = healthy(), tumor()
        self.assertEqual((g.count(HEALTHY), g.count(TUMOR), g.count(EMPTY)), (1, 1, 2))

    def test_field_shape_mismatch_raises(self) -> None:
        with self.assertRaises(ValueError):
            Grid([[empty_cell()]], [[0.0, 0.0]], [[0.0]], [[0.0]], [])

    def test_cell_on_vessel_raises(self) -> None:
        with self.assertRaises(ValueError):
            Grid([[healthy()]], [[0.0]], [[0.0]], [[0.0]], [(0, 0)])

    def test_daughter_copies_traits_but_is_a_new_object(self) -> None:
        parent = Cell(TUMOR, 0.7, 2.5, 3.0)
        child = daughter_of(parent)
        self.assertIsNot(child, parent)
        self.assertEqual((child.kind, child.mito_capacity, child.acid_tolerance, child.glucose_uptake),
                         (TUMOR, 0.7, 2.5, 3.0))


class TestMetabolism(unittest.TestCase):
    def test_matches_solve_cell_and_empty_squares_do_nothing(self) -> None:
        g = blank_grid(1, 2)
        g.cells[0][0] = tumor()
        g.glucose[0][0], g.oxygen[0][0] = 1.0, 10.0
        g.glucose[0][1], g.oxygen[0][1] = 5.0, 5.0      # empty square: ignored
        atp_made, dg, do, dl = metabolism(g)
        self.assertAlmostEqual(atp_made[0][0], 9.5)       # the Warburg cell from step 2
        self.assertAlmostEqual(dg[0][0], -1.0)
        self.assertAlmostEqual(do[0][0], -1.5)
        self.assertAlmostEqual(dl[0][0], 1.5)
        self.assertEqual((atp_made[0][1], dg[0][1], do[0][1], dl[0][1]), (0, 0, 0, 0))

    def test_glucose_uptake_caps_intake(self) -> None:
        g = blank_grid(1, 1)
        g.cells[0][0] = Cell(HEALTHY, 10, 1.0, glucose_uptake=0.5)
        g.glucose[0][0], g.oxygen[0][0] = 2.0, 100.0
        atp_made, dg, _, _ = metabolism(g)
        self.assertAlmostEqual(dg[0][0], -0.5)
        self.assertAlmostEqual(atp_made[0][0], 16.0)      # half of 32

    def test_fraction_scales_everything(self) -> None:
        g = blank_grid(1, 1)
        g.cells[0][0] = tumor(mito=0.5)
        g.glucose[0][0], g.oxygen[0][0] = 10.0, 100.0
        full = metabolism(g)
        quarter = metabolism(g, 0.25)
        for f_field, q_field in zip(full, quarter):
            self.assertAlmostEqual(q_field[0][0], f_field[0][0] / 4)

    def test_local_glucose_caps_intake(self) -> None:
        g = blank_grid(1, 1)
        g.cells[0][0] = Cell(HEALTHY, 10, 1.0, glucose_uptake=5.0)
        g.glucose[0][0], g.oxygen[0][0] = 0.25, 100.0
        _, dg, _, _ = metabolism(g)
        self.assertAlmostEqual(dg[0][0], -0.25)


class TestChooseDivisions(unittest.TestCase):
    def contest(self, left: Cell, right: Cell, params: Params, trials: int = 4000,
                reverse: bool = False) -> float:
        """1x3 grid [left, empty, right]: both must claim the middle. Fraction left wins."""
        cells = [[left, empty_cell(), right]]
        dividers = [(0, 2), (0, 0)] if reverse else [(0, 0), (0, 2)]
        wins = 0
        for seed in range(trials):
            result = choose_divisions(cells, dividers, [], params, random.Random(seed))
            wins += result[(0, 1)] == (0, 0)
        return wins / trials

    def test_tumor_beats_healthy_three_to_one(self) -> None:
        self.assertAlmostEqual(self.contest(tumor(), healthy(), Params()), 0.75, delta=0.03)

    def test_order_of_dividers_does_not_matter(self) -> None:
        self.assertAlmostEqual(self.contest(tumor(), healthy(), Params(), reverse=True), 0.75, delta=0.03)

    def test_same_kind_is_even(self) -> None:
        self.assertAlmostEqual(self.contest(tumor(), tumor(), Params()), 0.5, delta=0.03)
        self.assertAlmostEqual(self.contest(healthy(), healthy(), Params()), 0.5, delta=0.03)

    def test_weight_one_removes_the_advantage(self) -> None:
        self.assertAlmostEqual(self.contest(tumor(), healthy(), Params(tumor_weight=1.0)), 0.5, delta=0.03)

    def test_someone_always_wins_a_claimed_square(self) -> None:
        cells = [[tumor(), empty_cell(), healthy()]]
        for seed in range(50):
            result = choose_divisions(cells, [(0, 0), (0, 2)], [], Params(), random.Random(seed))
            self.assertEqual(list(result), [(0, 1)])

    def test_no_empty_neighbor_no_division(self) -> None:
        cells = [[healthy(), healthy()]]
        self.assertEqual(choose_divisions(cells, [(0, 0), (0, 1)], [], Params(), random.Random(0)), {})

    def test_never_divides_into_a_vessel(self) -> None:
        cells = [[healthy(), empty_cell()]]
        self.assertEqual(choose_divisions(cells, [(0, 0)], [(0, 1)], Params(), random.Random(0)), {})

    def test_target_is_always_an_empty_neighbor(self) -> None:
        cells = [[empty_cell() for _ in range(5)] for _ in range(5)]
        cells[2][2] = healthy()
        for seed in range(100):
            (target, parent), = choose_divisions(cells, [(2, 2)], [], Params(), random.Random(seed)).items()
            self.assertEqual(parent, (2, 2))
            self.assertIn(target, neighbors(2, 2, 5, 5))

    def test_same_seed_same_answer(self) -> None:
        cells = [[tumor(), empty_cell(), healthy()]]
        a = choose_divisions(cells, [(0, 0), (0, 2)], [], Params(), random.Random(42))
        b = choose_divisions(cells, [(0, 0), (0, 2)], [], Params(), random.Random(42))
        self.assertEqual(a, b)


class TestStepByHand(unittest.TestCase):
    """Small grids whose result can be worked out on paper."""

    def test_starving_tumor_dies_and_leaves_its_lactate(self) -> None:
        # 5x5, no vessels, one tumor in the middle with 1 glucose and no O2.
        # It ferments: 2 lactate, 2 ATP < 5, so it starves. Lactate then
        # diffuses at rate 0.5: middle 2(1 - 0.5) = 1.0, sides 2(0.2)(0.5) = 0.2,
        # corners 2(0.05)(0.5) = 0.05. Closed box, so total lactate stays 2.
        g = blank_grid(5, 5)
        g.cells[2][2] = tumor()
        g.glucose[2][2] = 1.0
        new = step(g, Params(), random.Random(0))
        self.assertEqual(new.count(TUMOR), 0)
        self.assertAlmostEqual(total(new.lactate), 2.0)
        self.assertAlmostEqual(new.lactate[2][2], 1.0)
        self.assertAlmostEqual(new.lactate[1][2], 0.2)
        self.assertAlmostEqual(new.lactate[1][1], 0.05)
        self.assertAlmostEqual(total(new.glucose), 0.0)

    def test_well_fed_healthy_cell_divides_next_to_itself(self) -> None:
        # 1 glucose + 6 O2: full respiration, 32 ATP >= 20, so it divides.
        g = blank_grid(5, 5)
        g.cells[2][2] = healthy()
        g.glucose[2][2], g.oxygen[2][2] = 1.0, 6.0
        new = step(g, Params(), random.Random(0))
        self.assertEqual(new.count(HEALTHY), 2)
        daughters = [(r, c) for r in range(5) for c in range(5)
                     if new.cells[r][c].kind == HEALTHY and (r, c) != (2, 2)]
        self.assertEqual(len(daughters), 1)
        self.assertIn(daughters[0], neighbors(2, 2, 5, 5))
        self.assertAlmostEqual(total(new.oxygen), 0.0)     # all 6 O2 used
        self.assertAlmostEqual(total(new.lactate), 0.0)    # no fermentation

    def test_vessels_refill_and_drain(self) -> None:
        # Empty 5x5 with lactate 5 everywhere and a vessel row on top.
        # Vessels drain 5 squares of lactate (25 - 5 = 20 squares x 5 = 100 left)
        # and add O2 10 and glucose 2 on 5 squares. Diffusion keeps totals.
        g = blank_grid(5, 5, vessels=[(0, c) for c in range(5)])
        g.lactate = make_field(5, 5, 5.0)
        new = step(g, Params(), random.Random(0))
        self.assertAlmostEqual(total(new.lactate), 100.0)
        self.assertAlmostEqual(total(new.oxygen), 50.0)
        self.assertAlmostEqual(total(new.glucose), 10.0)

    def test_deaths_free_space_before_divisions(self) -> None:
        # [tumor, healthy]: the healthy cell has no food and dies; the tumor
        # (well fed, acid-proof) divides into the square that just opened.
        g = blank_grid(1, 2)
        g.cells[0][0] = tumor(mito=10, acid_tol=100)
        g.cells[0][1] = healthy()
        g.glucose[0][0], g.oxygen[0][0] = 1.0, 6.0
        new = step(g, Params(), random.Random(0))
        self.assertEqual([c.kind for c in new.cells[0]], [TUMOR, TUMOR])

    def test_acid_kills_healthy_but_not_tumor(self) -> None:
        # Both well fed, lactate 2 everywhere: above healthy's 1.0, below tumor's 3.0.
        g = blank_grid(1, 2)
        g.cells[0][0], g.cells[0][1] = tumor(mito=10), healthy()
        g.glucose = make_field(1, 2, 1.0)
        g.oxygen = make_field(1, 2, 6.0)
        g.lactate = make_field(1, 2, 2.0)
        new = step(g, Params(), random.Random(0))
        self.assertEqual([c.kind for c in new.cells[0]], [TUMOR, TUMOR])

    def test_full_grid_cannot_divide(self) -> None:
        g = blank_grid(1, 2)
        g.cells[0][0], g.cells[0][1] = healthy(), healthy()
        g.glucose = make_field(1, 2, 1.0)
        g.oxygen = make_field(1, 2, 6.0)
        new = step(g, Params(), random.Random(0))
        self.assertEqual(new.count(HEALTHY), 2)

    def test_input_grid_not_changed(self) -> None:
        g = blank_grid(3, 3)
        g.cells[1][1] = healthy()
        g.glucose[1][1], g.oxygen[1][1] = 1.0, 6.0
        before = ([[c.kind for c in row] for row in g.cells],
                  [row[:] for row in g.oxygen], [row[:] for row in g.glucose])
        step(g, Params(), random.Random(0))
        after = ([[c.kind for c in row] for row in g.cells],
                 [row[:] for row in g.oxygen], [row[:] for row in g.glucose])
        self.assertEqual(before, after)

    def test_lactate_made_matches_fba(self) -> None:
        # Closed box, Warburg tumor: lactate after one step = what solve_cell says.
        g = blank_grid(3, 3)
        g.cells[1][1] = tumor(mito=0.5)
        g.glucose[1][1], g.oxygen[1][1] = 1.0, 10.0
        new = step(g, Params(), random.Random(0))
        self.assertAlmostEqual(total(new.lactate), solve_cell(1.0, 10.0, 0.5)["v_out"])



class TestSubstepsAndCapillaries(unittest.TestCase):
    def test_substeps_one_matches_old_behavior(self) -> None:
        # Default Params has 1 substep: the hand-worked tests above use it.
        self.assertEqual(Params().diffusion_substeps, 1)

    def test_substeps_with_plenty_scale_chemicals_by_k_times_f(self) -> None:
        # Plenty everywhere (closed box): a Warburg cell makes 1.5 lactate per
        # full step, so k substeps eating f each make 1.5 x k x f.
        for k, f in [(1, 1.0), (4, 1.0), (4, 0.25), (10, 0.1), (10, 0.5)]:
            g = blank_grid(3, 3)
            g.cells[1][1] = tumor(mito=0.5)
            g.glucose = make_field(3, 3, 100.0)
            g.oxygen = make_field(3, 3, 1000.0)
            new = step(g, Params(diffusion_substeps=k, uptake_fraction=f), random.Random(0))
            self.assertAlmostEqual(total(new.lactate), 1.5 * k * f, msg=f"k={k}, f={f}")

    def test_atp_is_an_average_rate(self) -> None:
        # Healthy cell, plenty of everything: its ATP rate is 32 whatever k and
        # f are, so it divides every time.
        for k, f in [(1, 1.0), (3, 1.0), (8, 0.125), (8, 0.5)]:
            g = blank_grid(3, 3)
            g.cells[1][1] = healthy()
            g.glucose = make_field(3, 3, 100.0)
            g.oxygen = make_field(3, 3, 1000.0)
            new = step(g, Params(diffusion_substeps=k, uptake_fraction=f), random.Random(0))
            self.assertEqual(new.count(HEALTHY), 2, msg=f"k={k}, f={f}")

    def test_bad_uptake_fraction_raises(self) -> None:
        for f in (0, -0.5, 1.5):
            with self.assertRaises(ValueError):
                Params(uptake_fraction=f)

    def test_starved_square_limits_each_substep(self) -> None:
        # Only 0.1 glucose and no vessel: the cell can never take more than
        # what's there, so total glucose taken in is at most 0.1.
        g = blank_grid(1, 1)
        g.cells[0][0] = Cell(HEALTHY, 10, 1.0, glucose_uptake=1.0)
        g.glucose[0][0], g.oxygen[0][0] = 0.1, 100.0
        new = step(g, Params(diffusion_substeps=5, uptake_fraction=0.2), random.Random(0))
        self.assertAlmostEqual(total(new.glucose), 0.0)
        self.assertEqual(new.count(HEALTHY), 0)            # 3.2 ATP < 5: starves

    def test_vessels_refill_between_substeps(self) -> None:
        # Empty 1x3, vessel on the left: more substeps -> more O2 delivered.
        g = blank_grid(1, 3, vessels=[(0, 0)])
        one = step(g, Params(diffusion_substeps=1), random.Random(0))
        five = step(g, Params(diffusion_substeps=5), random.Random(0))
        self.assertGreater(five.oxygen[0][2], one.oxygen[0][2])

    def test_bad_substeps_raises(self) -> None:
        with self.assertRaises(ValueError):
            Params(diffusion_substeps=0)

    def test_capillaries_layout(self) -> None:
        self.assertEqual(capillaries(6, 6, 3), [(1, 1), (1, 4), (4, 1), (4, 4)])
        self.assertEqual(capillaries(3, 3, 1), [(r, c) for r in range(3) for c in range(3)])
        # Spacing 5 starts 2 squares in, which is off a 2x2 grid: no vessels.
        self.assertEqual(capillaries(2, 2, 5), [])

    def test_capillaries_bad_spacing(self) -> None:
        with self.assertRaises(ValueError):
            capillaries(5, 5, 0)


if __name__ == "__main__":
    unittest.main()
