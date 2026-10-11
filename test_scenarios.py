"""Tests for scenarios.py. Run from this folder: python3 -m unittest -v"""
import random
import unittest

from datatypes import EMPTY, HEALTHY
from rules import step
from scenarios import CAPILLARY_SPACING, calibrated_params, healthy_tissue


class TestScenarios(unittest.TestCase):
    def test_calibrated_values_and_overrides(self) -> None:
        p = calibrated_params()
        self.assertEqual((p.diffusion_substeps, p.uptake_fraction, p.vessel_oxygen, p.vessel_glucose),
                         (10, 0.02, 10.0, 5.0))
        self.assertEqual(calibrated_params(tumor_weight=1.0).tumor_weight, 1.0)

    def test_healthy_tissue_layout(self) -> None:
        p = calibrated_params()
        g = healthy_tissue(12, 12, CAPILLARY_SPACING, p)
        self.assertEqual(g.vessels, [(3, 3), (3, 9), (9, 3), (9, 9)])
        for r, c in g.vessels:
            self.assertEqual(g.cells[r][c].kind, EMPTY)
        self.assertEqual(g.count(HEALTHY), 144 - 4)
        self.assertEqual(g.oxygen[0][0], p.vessel_oxygen)
        self.assertEqual(g.glucose[5][5], p.vessel_glucose)

    def test_healthy_tissue_stays_alive(self) -> None:
        # The calibration result, as a guard: no tumor, nobody dies.
        p = calibrated_params()
        g = healthy_tissue(12, 12, CAPILLARY_SPACING, p)
        rng = random.Random(0)
        for _ in range(40):
            g = step(g, p, rng)
        self.assertEqual(g.count(HEALTHY), 140)


if __name__ == "__main__":
    unittest.main()
