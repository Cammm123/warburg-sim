"""
Pictures of a run, standard library only.

- grid_svg: one frame as an SVG (cells, oxygen, lactate side by side).
- viewer_html: every frame in one HTML page with a slider and Play button. Frames
  are stored compactly (one letter per cell, fields as 0-255 numbers) and drawn
  on a canvas, so long runs stay small.
- line_chart_svg: how many healthy and tumor cells there are over time.
"""
import json

from datatypes import EMPTY, HEALTHY, TUMOR, Field, Grid

CELL_COLORS = {EMPTY: "#F1EFE8", HEALTHY: "#1D9E75", TUMOR: "#D85A30"}
VESSEL_COLOR = "#A32D2D"
KIND_LETTERS = {EMPTY: ".", HEALTHY: "H", TUMOR: "T"}


def encode_cells(grid: Grid) -> str:
    """One letter per square, row by row: V vessel, H healthy, T tumor, . empty."""
    vessels = set(grid.vessels)
    return "".join("V" if (r, c) in vessels else KIND_LETTERS[cell.kind]
                   for r, row in enumerate(grid.cells) for c, cell in enumerate(row))


def quantize(field: Field, top_value: float) -> list[int]:
    """Field row by row as integers 0-255 (top_value or more -> 255)."""
    if top_value <= 0:
        raise ValueError("top_value must be positive")
    return [min(255, max(0, round(255 * value / top_value))) for row in field for value in row]


def shade(value: float, top_value: float, dark: tuple[int, int, int]) -> str:
    """Color from white (0) to dark (top_value or more)."""
    t = min(max(value / top_value, 0.0), 1.0)
    r, g, b = (round(255 + (d - 255) * t) for d in dark)
    return f"rgb({r},{g},{b})"


def grid_svg(grid: Grid, title: str, oxygen_top: float, lactate_top: float, square: int = 6) -> str:
    """One frame: cells | oxygen | lactate, side by side."""
    n_rows, n_cols = grid.num_rows(), grid.num_cols()
    panel_w, gap, top = n_cols * square, 24, 44
    width, height = 3 * panel_w + 2 * gap + 20, top + n_rows * square + 30
    vessels = set(grid.vessels)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="100%" viewBox="0 0 {width} {height}" '
             f'font-family="sans-serif" font-size="12">',
             f'<rect width="100%" height="100%" fill="white"/>',
             f'<text x="10" y="16" font-size="14">{title}</text>']
    panels = [("cells (green healthy, orange tumor, red vessel)", None),
              ("oxygen", (grid.oxygen, oxygen_top, (24, 95, 165))),
              ("lactate (acid)", (grid.lactate, lactate_top, (153, 60, 29)))]
    for p, (label, field_info) in enumerate(panels):
        x0 = 10 + p * (panel_w + gap)
        parts.append(f'<text x="{x0}" y="{top - 8}">{label}</text>')
        for r in range(n_rows):
            for c in range(n_cols):
                if field_info is None:
                    color = VESSEL_COLOR if (r, c) in vessels else CELL_COLORS[grid.cells[r][c].kind]
                else:
                    field, top_value, dark = field_info
                    color = shade(field[r][c], top_value, dark)
                parts.append(f'<rect x="{x0 + c * square}" y="{top + r * square}" '
                             f'width="{square}" height="{square}" fill="{color}"/>')
    parts.append("</svg>")
    return "\n".join(parts)


def line_chart_svg(series: dict[str, list[float]], colors: dict[str, str], title: str,
                   y_label: str, x_step: int = 1) -> str:
    """Lines over time. series maps a name to its values; x is the frame number x x_step."""
    width, height, left, top, plot_w, plot_h = 640, 310, 60, 36, 540, 210
    length = max(len(values) for values in series.values())
    top_value = max(1.0, max(max(values) for values in series.values()))
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="100%" viewBox="0 0 {width} {height}" '
             f'font-family="sans-serif" font-size="12">', '<rect width="100%" height="100%" fill="white"/>',
             f'<text x="{left}" y="20" font-size="14">{title}</text>',
             f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#888"/>',
             f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#888"/>',
             f'<text x="{left + 20}" y="{top + plot_h + 34}" text-anchor="middle">step</text>',
             f'<text x="{left - 8}" y="{top + 4}" text-anchor="end">{top_value:.3g}</text>',
             f'<text x="{left - 8}" y="{top + plot_h + 4}" text-anchor="end">0</text>',
             f'<text x="{left - 8}" y="{top - 10}" text-anchor="start">{y_label}</text>']
    for i in range(5):
        frame = round(i * (length - 1) / 4) if length > 1 else 0
        x = left + plot_w * frame / max(length - 1, 1)
        parts.append(f'<text x="{x:.1f}" y="{top + plot_h + 16}" text-anchor="middle">{frame * x_step}</text>')
    for k, (name, values) in enumerate(series.items()):
        points = " ".join(f"{left + plot_w * i / max(length - 1, 1):.1f},{top + plot_h * (1 - v / top_value):.1f}"
                          for i, v in enumerate(values))
        parts.append(f'<polyline points="{points}" fill="none" stroke="{colors[name]}" stroke-width="2.5"/>')
        parts.append(f'<rect x="{left + 200 + 170 * k}" y="{height - 14}" width="12" height="3" fill="{colors[name]}"/>')
        parts.append(f'<text x="{left + 218 + 170 * k}" y="{height - 9}">{name}</text>')
    parts.append("</svg>")
    return "\n".join(parts)


VIEWER_TEMPLATE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
body { font-family: sans-serif; margin: 16px; background: #fff; color: #222; }
.panels { display: flex; gap: 16px; flex-wrap: wrap; }
.panels div { font-size: 13px; }
canvas { image-rendering: pixelated; width: 300px; max-width: 90vw; border: 1px solid #ccc; display: block; }
.controls { display: flex; gap: 12px; align-items: center; margin: 12px 0; flex-wrap: wrap; }
input[type=range] { width: 320px; max-width: 80vw; }
</style></head><body>
<h2 style="font-weight:500">__TITLE__</h2>
<div class="controls"><button id="play">Play</button>
<input type="range" id="frame" min="0" value="0"><span id="label"></span></div>
<div class="panels">
<div>cells (green healthy, orange tumor, red vessel)<canvas id="cells"></canvas></div>
<div>oxygen<canvas id="oxygen"></canvas></div>
<div>lactate (acid)<canvas id="lactate"></canvas></div>
</div>
<script>
const DATA = __DATA__;
const colors = {".": [241,239,232], "H": [29,158,117], "T": [216,90,48], "V": [163,45,45]};
const darks = {oxygen: [24,95,165], lactate: [153,60,29]};
const slider = document.getElementById("frame"), label = document.getElementById("label");
slider.max = DATA.frames.length - 1;
for (const id of ["cells", "oxygen", "lactate"]) {
  const c = document.getElementById(id); c.width = DATA.cols; c.height = DATA.rows;
}
function paint(id, rgbAt) {
  const c = document.getElementById(id), ctx = c.getContext("2d");
  const img = ctx.createImageData(DATA.cols, DATA.rows);
  for (let i = 0; i < DATA.rows * DATA.cols; i++) {
    const [r, g, b] = rgbAt(i); img.data.set([r, g, b, 255], 4 * i);
  }
  ctx.putImageData(img, 0, 0);
}
function shade(v, dark) { const t = v / 255; return dark.map(d => Math.round(255 + (d - 255) * t)); }
function show(n) {
  const f = DATA.frames[n];
  paint("cells", i => colors[f.cells[i]]);
  paint("oxygen", i => shade(f.oxygen[i], darks.oxygen));
  paint("lactate", i => shade(f.lactate[i], darks.lactate));
  label.textContent = "step " + f.step + ": " + f.healthy + " healthy, " + f.tumor + " tumor";
}
slider.oninput = () => show(+slider.value);
let timer = null;
document.getElementById("play").onclick = (e) => {
  if (timer) { clearInterval(timer); timer = null; e.target.textContent = "Play"; return; }
  e.target.textContent = "Pause";
  timer = setInterval(() => {
    slider.value = (+slider.value + 1) % DATA.frames.length; show(+slider.value);
  }, 120);
};
show(0);
</script></body></html>
"""


def frame_data(grid: Grid, step_number: int, oxygen_top: float, lactate_top: float) -> dict:
    """One frame for the viewer: compact cells string, quantized fields, counts."""
    return {"step": step_number, "cells": encode_cells(grid),
            "oxygen": quantize(grid.oxygen, oxygen_top), "lactate": quantize(grid.lactate, lactate_top),
            "healthy": grid.count(HEALTHY), "tumor": grid.count(TUMOR)}


def viewer_html(frames: list[dict], num_rows: int, num_cols: int, title: str) -> str:
    """The whole HTML page for a list of frame_data frames."""
    data = json.dumps({"rows": num_rows, "cols": num_cols, "frames": frames}, separators=(",", ":"))
    return VIEWER_TEMPLATE.replace("__TITLE__", title).replace("__DATA__", data)
