"""Tests for the step 9 switches (acid_harm, neutral marker) and experiments.py."""
import random
import statistics
import unittest

from datatypes import HEALTHY, TUMOR, Cell, Grid, Params, empty_cell
from experiments import budget_tumor, run_once, run_seeds
from fields import make_field
from rules import cell_fluxes, make_daughter, step


class TestAcidHarmSwitch(unittest.TestCase):
    def acid_grid(self) -> Grid:
        # Well-fed healthy cell sitting in lactate 5 (tolerance 1).
        cells = [[Cell(HEALTHY, 10, 1.0)]]
        return Grid(cells, make_field(1, 1, 100.0), make_field(1, 1, 10.0), make_field(1, 1, 5.0), [])

    def test_acid_kills_when_on(self) -> None:
        self.assertEqual(step(self.acid_grid(), Params(), random.Random(0)).count(HEALTHY), 0)

    def test_acid_harmless_when_off(self) -> None:
        self.assertEqual(step(self.acid_grid(), Params(acid_harm=False), random.Random(0)).count(HEALTHY), 1)


class TestNeutralMarker(unittest.TestCase):
    def test_marker_drifts_with_the_mutation_size(self) -> None:
        rng = random.Random(0)
        parent = Cell(TUMOR, 0.5, 3.0, 8.0)
        values = [make_daughter(parent, Params(mutation_sd=0.1), rng).marker for _ in range(5000)]
        self.assertAlmostEqual(statistics.mean(values), 0.0, delta=0.01)
        self.assertAlmostEqual(statistics.stdev(values), 0.1, delta=0.01)

    def test_marker_copies_for_healthy_and_without_mutation(self) -> None:
        healthy = Cell(HEALTHY, 10, 1.0, marker=0.3)
        self.assertEqual(make_daughter(healthy, Params(mutation_sd=0.1), random.Random(0)).marker, 0.3)
        tumor = Cell(TUMOR, 0.5, 3.0, marker=0.3)
        self.assertEqual(make_daughter(tumor, Params(), random.Random(0)).marker, 0.3)

    def test_marker_has_no_effect_on_metabolism(self) -> None:
        a = Cell(TUMOR, 0.5, 3.0, 8.0, marker=0.0)
        b = Cell(TUMOR, 0.5, 3.0, 8.0, marker=99.0)
        self.assertEqual(cell_fluxes(a, 5, 5, Params()), cell_fluxes(b, 5, 5, Params()))


class TestExperimentRunner(unittest.TestCase):
    def test_run_once_returns_every_metric(self) -> None:
        result = run_once(Params(diffusion_substeps=2, uptake_fraction=0.02, vessel_glucose=5.0,
                                 mutation_sd=0.05), budget_tumor, seed=0, size=12, steps=3)
        self.assertEqual(set(result), {"healthy", "tumor", "mito_change", "share_change",
                                       "marker_change", "warburg"})
        self.assertEqual(result["mito_change"], 0.0)   # budget cells don't use M

    def test_same_seed_same_result(self) -> None:
        params = Params(diffusion_substeps=2, uptake_fraction=0.02, vessel_glucose=5.0, mutation_sd=0.05)
        self.assertEqual(run_once(params, budget_tumor, seed=4, size=12, steps=5),
                         run_once(params, budget_tumor, seed=4, size=12, steps=5))

    def test_run_seeds_gives_mean_and_spread(self) -> None:
        params = Params(diffusion_substeps=2, uptake_fraction=0.02, vessel_glucose=5.0, mutation_sd=0.05)
        result = run_seeds(params, budget_tumor, seeds=(1, 2), size=12, steps=3)
        mean, sd = result["tumor"]
        self.assertGreaterEqual(mean, 0)
        self.assertGreaterEqual(sd, 0)


if __name__ == "__main__":
    unittest.main()
