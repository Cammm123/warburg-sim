"""Tests for render.py. Run from this folder: python3 -m unittest -v"""
import json
import unittest

from datatypes import HEALTHY, TUMOR, Cell, Grid, empty_cell
from fields import make_field
from render import encode_cells, frame_data, grid_svg, line_chart_svg, quantize, shade, viewer_html


def tiny_grid() -> Grid:
    """2x3: vessel, healthy, tumor / empty, empty, healthy."""
    cells = [[empty_cell(), Cell(HEALTHY), Cell(TUMOR)], [empty_cell(), empty_cell(), Cell(HEALTHY)]]
    return Grid(cells, make_field(2, 3, 5.0), make_field(2, 3, 1.0), make_field(2, 3, 0.5), [(0, 0)])


class TestEncoding(unittest.TestCase):
    def test_encode_cells(self) -> None:
        self.assertEqual(encode_cells(tiny_grid()), "VHT..H")

    def test_quantize(self) -> None:
        self.assertEqual(quantize([[0.0, 5.0], [10.0, 20.0]], 10.0), [0, 128, 255, 255])

    def test_quantize_bad_top_raises(self) -> None:
        with self.assertRaises(ValueError):
            quantize([[1.0]], 0)

    def test_shade_endpoints(self) -> None:
        self.assertEqual(shade(0, 10, (24, 95, 165)), "rgb(255,255,255)")
        self.assertEqual(shade(10, 10, (24, 95, 165)), "rgb(24,95,165)")
        self.assertEqual(shade(50, 10, (24, 95, 165)), "rgb(24,95,165)")   # capped

    def test_frame_data_counts(self) -> None:
        f = frame_data(tiny_grid(), 7, oxygen_top=10.0, lactate_top=1.0)
        self.assertEqual((f["step"], f["healthy"], f["tumor"]), (7, 2, 1))
        self.assertEqual(f["oxygen"], [128] * 6)
        self.assertEqual(len(f["lactate"]), 6)


class TestPictures(unittest.TestCase):
    def test_grid_svg_has_three_panels_of_squares(self) -> None:
        svg = grid_svg(tiny_grid(), "test", 10.0, 1.0)
        self.assertTrue(svg.startswith("<svg") and svg.endswith("</svg>"))
        self.assertEqual(svg.count('width="6" height="6"'), 3 * 6)

    def test_line_chart_has_one_line_per_series(self) -> None:
        svg = line_chart_svg({"healthy": [3, 2, 1], "tumor": [0, 1, 2]},
                             {"healthy": "#1D9E75", "tumor": "#D85A30"}, "t", "cells")
        self.assertEqual(svg.count("<polyline"), 2)

    def test_viewer_embeds_valid_json(self) -> None:
        frames = [frame_data(tiny_grid(), s, 10.0, 1.0) for s in (0, 1)]
        html = viewer_html(frames, 2, 3, "Run")
        self.assertIn("<title>Run</title>", html)
        start = html.index("const DATA = ") + len("const DATA = ")
        data = json.loads(html[start:html.index(";\n", start)])
        self.assertEqual((data["rows"], data["cols"], len(data["frames"])), (2, 3, 2))


if __name__ == "__main__":
    unittest.main()
