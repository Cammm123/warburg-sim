"""
Look at what fba.py does, instead of only trusting the tests.

Run from this folder:  python3 explore.py
Prints S, S v for a few flux vectors, and solve_cell on a few cells, then writes
output/cell_response.svg: how one cell's ATP and lactate change with oxygen.
Standard library only.
"""
import os

from fba import (METABOLITES, REACTIONS, atp, flux_vector, mat_vec, oxygen_used,
                 solve_cell, stoichiometric_matrix)


def print_matrix(matrix: list[list[float]]) -> None:
    """Print S with metabolite names on the rows and reaction names on top."""
    print(" " * 10 + "".join(f"{name:>7}" for name in REACTIONS))
    for name, row in zip(METABOLITES, matrix):
        print(f"{name:<10}" + "".join(f"{entry:>7}" for entry in row))


def print_sv(label: str, v: list[float]) -> None:
    """Print v and S v, naming which metabolites pile up or run out."""
    sv = mat_vec(stoichiometric_matrix(), v)
    changes = [f"{m} {'+' if x > 0 else ''}{x:g}" for m, x in zip(METABOLITES, sv) if x != 0]
    verdict = "steady state (S v = 0)" if not changes else "NOT steady: " + ", ".join(changes)
    print(f"  {label:<26} v = {v}  ->  S v = {sv}   {verdict}")


def print_cells(cases: list[tuple[str, float, float, float]]) -> None:
    """Run solve_cell on each (label, glucose, oxygen, mito_capacity) and print a table."""
    print(f"  {'cell':<30}{'G':>5}{'O2':>6}{'M':>6} | {'respired r':>10}{'lactate':>9}{'ATP':>7}{'O2 used':>9}")
    for label, g, o, m in cases:
        f = solve_cell(g, o, m)
        print(f"  {label:<30}{g:>5g}{o:>6g}{m:>6g} | {f['v_m']:>10.2f}{f['v_out']:>9.2f}"
              f"{atp(f):>7.1f}{oxygen_used(f):>9.2f}")


def response_svg(path: str) -> None:
    """Draw ATP and lactate vs oxygen for one glucose, at three mito capacities."""
    width, height, left, top, plot_w, plot_h = 640, 320, 60, 30, 520, 220
    max_o, max_atp = 8.0, 34.0
    capacities = [(10.0, "#1D9E75", "M = 10 (healthy)"), (1.0, "#BA7517", "M = 1"),
                  (0.5, "#D85A30", "M = 0.5 (Warburg)")]

    def x(o: float) -> float:
        return left + plot_w * o / max_o

    def y(value: float, top_value: float) -> float:
        return top + plot_h * (1 - value / top_value)

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height * 2}" '
             f'font-family="sans-serif" font-size="12">',
             f'<rect width="100%" height="100%" fill="white"/>']
    for panel, (title, key, top_value) in enumerate([("ATP made per glucose", "atp", max_atp),
                                                     ("Lactate made per glucose", "lactate", 2.2)]):
        dy = panel * height
        parts.append(f'<text x="{left}" y="{top - 10 + dy}" font-size="14">{title}</text>')
        parts.append(f'<line x1="{left}" y1="{top + plot_h + dy}" x2="{left + plot_w}" '
                     f'y2="{top + plot_h + dy}" stroke="#888"/>')
        parts.append(f'<line x1="{left}" y1="{top + dy}" x2="{left}" y2="{top + plot_h + dy}" stroke="#888"/>')
        for o in range(0, 9, 2):
            parts.append(f'<text x="{x(o)}" y="{top + plot_h + 16 + dy}" text-anchor="middle">{o}</text>')
        parts.append(f'<text x="{left + plot_w / 2}" y="{top + plot_h + 34 + dy}" '
                     f'text-anchor="middle">oxygen available (6 = enough to fully respire)</text>')
        for tick in ([0, 8, 16, 24, 32] if key == "atp" else [0, 0.5, 1, 1.5, 2]):
            parts.append(f'<text x="{left - 8}" y="{y(tick, top_value) + 4 + dy}" text-anchor="end">{tick:g}</text>')
        for m, color, label in capacities:
            points = []
            for step in range(81):
                o = max_o * step / 80
                f = solve_cell(1.0, o, m)
                value = atp(f) if key == "atp" else f["v_out"]
                points.append(f"{x(o):.1f},{y(value, top_value) + dy:.1f}")
            parts.append(f'<polyline points="{" ".join(points)}" fill="none" stroke="{color}" stroke-width="2.5"/>')
        if panel == 0:
            for i, (m, color, label) in enumerate(capacities):
                parts.append(f'<rect x="{left + 20}" y="{top + 8 + 18 * i}" width="12" height="3" fill="{color}"/>')
                parts.append(f'<text x="{left + 38}" y="{top + 13 + 18 * i}">{label}</text>')
    parts.append("</svg>")
    with open(path, "w") as file:
        file.write("\n".join(parts))


def main() -> None:
    print("STEP 1: the stoichiometric matrix S (rows = metabolites, columns = reactions)\n")
    print_matrix(stoichiometric_matrix())

    print("\nS v for some flux vectors (order: v_in, v_g, v_e, v_l, v_m, v_out)\n")
    print_sv("no O2, ferment", [1, 1, 0, 2, 0, 2])
    print_sv("O2, respire", [1, 1, 2, 0, 2, 0])
    print_sv("glycolysis alone", [1, 1, 0, 0, 0, 0])
    print_sv("ferment, lactate trapped", [1, 1, 0, 2, 0, 0])

    print("\nSTEP 2: solve_cell picks the best steady state for each cell\n")
    print_cells([
        ("healthy, plenty of O2", 1, 10, 10),
        ("healthy, no O2 (sprinting)", 1, 0, 10),
        ("healthy, half the O2 needed", 1, 3, 10),
        ("Warburg: O2 but small mitos", 1, 10, 0.5),
        ("Warburg, double glucose", 2, 10, 0.5),
    ])
    first = solve_cell(1, 10, 10)
    print(f"\n  Full flux vector for the first cell: {dict(zip(REACTIONS, flux_vector(first)))}")

    os.makedirs("output", exist_ok=True)
    response_svg(os.path.join("output", "cell_response.svg"))
    print("\nWrote output/cell_response.svg")


if __name__ == "__main__":
    main()
