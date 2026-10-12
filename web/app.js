// app.js - the page: three questions, a map with three views, charts, a
// "how it works" section, and folded-away settings and cell inspector.
// All the science is in model.js; this file only runs and draws it.

const M = window.WarburgModel;
const SIZE = 60;
const CELL_PX = 8;   // the map is drawn at 8 px per grid square


const VIEWS = {
  cells: "Each square is one cell. Healthy cells are green. Tumor cells are colored by how they make energy: orange = fermenting, purple = respiring.",
  oxygen: "Oxygen comes out of the capillaries (dark dots) and spreads. Dark blue = plenty, white = none. Hungry cells use it up, so squares far from capillaries run low.",
  acid: "Lactate is the acid fermenting cells give off. Capillaries wash it away. Healthy cells die where it passes their tolerance (1); tumor cells survive up to 3. Dark brown = 3 or more.",
};

const state = {
  preset: "invasion", grid: null, params: null, rng: null, seed: 1, stepCount: 0,
  playing: false, selected: null, view: "cells", acidHarmless: false, stepsPerSecond: 6, owed: 0, lastTime: null,
  history: { healthy: [], tumor: [], trait: [], warburg: [] },
};

const $ = (id) => document.getElementById(id);

// ---- Settings ------------------------------------------------------------------
function readSettings() {
  return { tolerance: +$("tolerance").value, glucose: +$("glucose").value, weight: +$("weight").value };
}

/** Params overrides from the settings, on top of the question's own (which set evolution). */
function overrides(settings) {
  return { vesselGlucose: settings.glucose, tumorWeight: settings.weight, acidHarm: !state.acidHarmless };
}

function showSettings() {
  const s = readSettings();
  $("toleranceOut").textContent = s.tolerance.toFixed(1);
  $("glucoseOut").textContent = s.glucose.toFixed(1);
  $("weightOut").textContent = `${s.weight % 1 ? s.weight.toFixed(1) : s.weight} : 1`;
  $("harmless").classList.toggle("active", state.preset === "invasion" && state.acidHarmless);
}

/** Apply settings to the running simulation (the next step uses them). */
function applySettings() {
  const s = readSettings();
  Object.assign(state.params, overrides(s));
  for (const row of state.grid.cells) for (const cell of row) if (cell.kind === M.HEALTHY) cell.acidTolerance = s.tolerance;
  showSettings();
  draw();
}

// ---- Starting a run -------------------------------------------------------------------
function startRun() {
  const s = readSettings();
  const { grid, params } = M.presetGrid(state.preset, SIZE, overrides(s), s.tolerance);
  Object.assign(state, { grid, params, rng: M.makeRng(state.seed), stepCount: 0 });
  state.history = { healthy: [], tumor: [], trait: [], warburg: [] };
  record();
  showSettings();
  draw();
}

/** Choose a question. acidHarmless: question 1's comparison run. */
function chooseQuestion(name, acidHarmless = false) {
  state.preset = name;
  const preset = M.PRESETS[name];
  state.acidHarmless = acidHarmless;
  $("tolerance").value = 1; $("glucose").value = 5; $("weight").value = 3;
  document.querySelectorAll(".question").forEach((q) => q.classList.toggle("active", q.dataset.name === name));
  document.querySelectorAll(".run").forEach((b) => b.classList.toggle("active", b.dataset.name === name && !acidHarmless));
  $("traitCaption").innerHTML = preset.trait
    ? '<span class="swatch" style="background:#d85a30"></span>share of sugar fermented (0 = all respired, 1 = all fermented) &nbsp; <span class="swatch" style="background:#534ab7"></span>share of the budget spent on mitochondria'
    : '<span class="swatch" style="background:#d85a30"></span>share of sugar the tumor ferments (0 = all respired, 1 = all fermented)';
  state.seed = 1;
  setPlaying(false);
  startRun();
}

function record() {
  const s = M.summary(state.grid, state.params);
  const h = state.history, trait = M.PRESETS[state.preset].trait;
  h.healthy.push(M.count(state.grid, M.HEALTHY));
  h.tumor.push(M.count(state.grid, M.TUMOR));
  h.trait.push(trait ? s[trait] : null);
  h.warburg.push(M.count(state.grid, M.TUMOR) ? s.meanWarburg : null);
}

function doStep() {
  state.grid = M.step(state.grid, state.params, state.rng);
  state.stepCount++;
  record();
}

// ---- Drawing ----------------------------------------------------------------------------
function hexToRgb(hex) { const n = parseInt(hex.slice(1), 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; }
const COLORS = { empty: hexToRgb("#f1efe8"), healthy: hexToRgb("#1d9e75"), vessel: "#2c2c2a",
                 ferment: hexToRgb("#d85a30"), respire: hexToRgb("#534ab7") };
const mix = (x, y, t) => x.map((v, i) => Math.round(v + (y[i] - v) * t));
const rgb = (c) => `rgb(${c[0]},${c[1]},${c[2]})`;

/** Warburg index of the tumor cell at (r, c), judged as the simulation eats. */
function cellIndex(r, c) {
  const g = state.grid;
  return M.warburgIndex(M.cellFluxes(g.cells[r][c], g.glucose[r][c], g.oxygen[r][c], state.params, state.params.uptakeFraction));
}

/** Color of square (r, c) in the current view. */
function squareColor(r, c) {
  const g = state.grid, cell = g.cells[r][c];
  if (state.view === "oxygen") return mix([255, 255, 255], [24, 95, 165], Math.min(g.oxygen[r][c] / state.params.vesselOxygen, 1));
  if (state.view === "acid") return mix([255, 255, 255], [113, 43, 19], Math.min(g.lactate[r][c] / 3, 1));
  if (cell.kind === M.HEALTHY) return COLORS.healthy;
  if (cell.kind === M.TUMOR) return mix(COLORS.respire, COLORS.ferment, cellIndex(r, c));
  return COLORS.empty;
}

function drawMap() {
  const canvas = $("map"), ctx = canvas.getContext("2d"), g = state.grid;
  canvas.width = SIZE * CELL_PX; canvas.height = SIZE * CELL_PX;
  for (let r = 0; r < SIZE; r++) {
    for (let c = 0; c < SIZE; c++) {
      ctx.fillStyle = rgb(squareColor(r, c));
      ctx.fillRect(c * CELL_PX, r * CELL_PX, CELL_PX, CELL_PX);
    }
  }
  ctx.fillStyle = COLORS.vessel;          // capillaries: dark dots
  for (const [r, c] of g.vessels) {
    ctx.beginPath(); ctx.arc(c * CELL_PX + CELL_PX / 2, r * CELL_PX + CELL_PX / 2, CELL_PX * 0.42, 0, 2 * Math.PI); ctx.fill();
  }
  if (state.view !== "cells") {           // outline the tumor so you can see where it is
    ctx.strokeStyle = "rgba(216, 90, 48, 0.9)"; ctx.lineWidth = 1;
    for (let r = 0; r < SIZE; r++) for (let c = 0; c < SIZE; c++)
      if (g.cells[r][c].kind === M.TUMOR) ctx.strokeRect(c * CELL_PX + 0.5, r * CELL_PX + 0.5, CELL_PX - 1, CELL_PX - 1);
  }
  if (state.selected && $("inspectBox").open) {
    const [r, c] = state.selected;
    ctx.lineWidth = 2; ctx.strokeStyle = "#111";
    ctx.strokeRect(c * CELL_PX + 1, r * CELL_PX + 1, CELL_PX - 2, CELL_PX - 2);
  }
}

function drawLegend() {
  const legends = {
    cells: '<span><span class="swatch" style="background:#1d9e75"></span>Healthy</span>' +
           '<span>Tumor: fermenting<span class="ramp"></span>respiring</span>' +
           '<span><span class="dot"></span>Capillary</span>' +
           '<span><span class="swatch" style="background:#f1efe8;border:1px solid #bbb"></span>Empty (a cell died)</span>',
    oxygen: '<span>No oxygen<span class="ramp" style="background:linear-gradient(90deg,#fff,#185fa5);border:1px solid #bbb"></span>blood level</span>' +
            '<span><span class="dot"></span>Capillary</span><span><span class="swatch outline"></span>Tumor cell</span>',
    acid: '<span>No acid<span class="ramp" style="background:linear-gradient(90deg,#fff,#712b13);border:1px solid #bbb"></span>3 or more</span>' +
          '<span><span class="dot"></span>Capillary</span><span><span class="swatch outline"></span>Tumor cell</span>',
  };
  $("legend").innerHTML = legends[state.view];
  $("viewCaption").textContent = VIEWS[state.view];
}

/** A crisp canvas at the device's pixel ratio; returns its context and CSS size. */
function setupCanvas(canvas) {
  const ratio = window.devicePixelRatio || 1, w = canvas.clientWidth, h = canvas.clientHeight;
  canvas.width = Math.round(w * ratio); canvas.height = Math.round(h * ratio);
  const ctx = canvas.getContext("2d");
  ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
  return { ctx, w, h };
}

const cssVar = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

/** Line chart; values may contain null (gap). yLabel goes at the top of the axis. */
function drawChart(id, series, top, yLabel) {
  const { ctx, w, h } = setupCanvas($(id));
  const left = 42, right = 8, upper = 14, lower = 18, pw = w - left - right, ph = h - upper - lower;
  ctx.clearRect(0, 0, w, h);
  ctx.strokeStyle = cssVar("--line"); ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(left, upper); ctx.lineTo(left, upper + ph); ctx.lineTo(left + pw, upper + ph); ctx.stroke();
  ctx.fillStyle = cssVar("--muted"); ctx.font = "11px system-ui, sans-serif"; ctx.textAlign = "right";
  ctx.fillText(top >= 10 ? Math.round(top) : top.toFixed(1), left - 4, upper + 8);
  ctx.fillText("0", left - 4, upper + ph);
  ctx.textAlign = "left"; ctx.fillText(yLabel, left + 4, upper - 3);
  const n = Math.max(...series.map((s) => s.values.length));
  ctx.fillText("step 0", left, h - 4);
  ctx.textAlign = "right"; ctx.fillText(`step ${n - 1}`, left + pw, h - 4);
  for (const s of series) {
    ctx.strokeStyle = s.color; ctx.lineWidth = 2; ctx.beginPath();
    let pen = false;
    s.values.forEach((v, i) => {
      if (v === null) { pen = false; return; }
      const x = left + (n > 1 ? (pw * i) / (n - 1) : 0), y = upper + ph * (1 - v / top);
      if (pen) ctx.lineTo(x, y); else ctx.moveTo(x, y);
      pen = true;
    });
    ctx.stroke();
  }
}

function drawCharts() {
  const h = state.history, living = SIZE * SIZE - state.grid.vessels.length;
  drawChart("countsChart", [{ values: h.healthy, color: "#1d9e75" }, { values: h.tumor, color: "#d85a30" }], living, "cells");
  const series = [{ values: h.warburg, color: "#d85a30" }];
  if (M.PRESETS[state.preset].trait) series.push({ values: h.trait, color: "#534ab7" });
  drawChart("traitChart", series, 1, "share (0 to 1)");
}

function draw() {
  drawMap();
  $("step").textContent = state.stepCount;
  $("healthy").textContent = M.count(state.grid, M.HEALTHY);
  $("tumor").textContent = M.count(state.grid, M.TUMOR);
  drawCharts();
  if ($("inspectBox").open) drawInspector();
}

// ---- Cell inspector ------------------------------------------------------------------------
/** Clip a polygon to the half-plane p a + q b <= r (Sutherland-Hodgman, one edge). */
function clip(poly, p, q, r) {
  const out = [], inside = ([a, b]) => p * a + q * b <= r + 1e-12;
  poly.forEach((cur, i) => {
    const prev = poly[(i + poly.length - 1) % poly.length];
    const cross = () => {
      const [a1, b1] = prev, [a2, b2] = cur, d1 = p * a1 + q * b1 - r, d2 = p * a2 + q * b2 - r, t = d1 / (d1 - d2);
      return [a1 + t * (a2 - a1), b1 + t * (b2 - b1)];
    };
    if (inside(cur)) { if (!inside(prev)) out.push(cross()); out.push(cur); }
    else if (inside(prev)) out.push(cross());
  });
  return out;
}

/** The selected cell's limits as lines p a + q b = r (per-step rates), with labels. */
function cellLimits(cell, glucose, oxygen) {
  const f = state.params.uptakeFraction, limits = [];
  limits.push({ p: 0, q: 1, r: oxygen / (6 * f), label: "oxygen" });
  if (cell.mitoShare !== null && state.params.plastic) {
    limits.push({ p: 1, q: 1, r: glucose / f, label: "glucose here" });
    limits.push({ p: 1, q: 1 + 2 * state.params.mitoCost, r: state.params.budget, label: "budget" });
  } else {
    const [uptake, mito] = cell.mitoShare !== null
      ? M.budgetTraits(cell.mitoShare, state.params.budget, state.params.mitoCost)
      : [cell.glucoseUptake, cell.mitoCapacity];
    limits.push({ p: 1, q: 1, r: Math.min(uptake, glucose / f), label: "glucose taken in" });
    limits.push({ p: 0, q: 1, r: mito / 2, label: "mitochondria" });
  }
  return limits;
}

function drawPlane(a, b, limits) {
  const { ctx, w, h } = setupCanvas($("plane"));
  const pad = 30, size = Math.min(w, h) - pad - 8;
  const span = Math.max(1, ...limits.filter((l) => l.p !== 0).map((l) => l.r / l.p)) * 1.15;
  const X = (v) => pad + (size * v) / span, Y = (v) => 8 + size - (size * v) / span;
  ctx.clearRect(0, 0, w, h);
  let region = [[0, 0], [span, 0], [span, span], [0, span]];
  for (const l of limits) region = clip(region, l.p, l.q, l.r);
  if (region.length) {
    ctx.fillStyle = "rgba(83, 74, 183, 0.14)"; ctx.beginPath();
    region.forEach(([x, y], i) => (i ? ctx.lineTo(X(x), Y(y)) : ctx.moveTo(X(x), Y(y))));
    ctx.closePath(); ctx.fill();
  }
  ctx.strokeStyle = cssVar("--muted"); ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(X(0), Y(span)); ctx.lineTo(X(0), Y(0)); ctx.lineTo(X(span), Y(0)); ctx.stroke();
  ctx.fillStyle = cssVar("--muted"); ctx.font = "11px system-ui, sans-serif";
  ctx.textAlign = "right"; ctx.fillText("a (fermented)", X(span), Y(0) + 16);
  ctx.save(); ctx.translate(12, Y(span / 2)); ctx.rotate(-Math.PI / 2); ctx.textAlign = "center";
  ctx.fillText("b (respired)", 0, 0); ctx.restore();
  ctx.setLineDash([4, 3]); ctx.textAlign = "left";
  for (const l of limits) {
    let pts;
    if (l.p === 0) pts = [[0, l.r / l.q], [span, l.r / l.q]];
    else pts = [[l.r / l.p, 0], [Math.max(0, (l.r - l.q * span) / l.p), Math.min(span, l.r / l.q)]];
    if (l.p === 0 && l.r / l.q > span) {
      ctx.fillText(`${l.label} limit is off the chart`, X(0) + 4, Y(span) + 10 + (l.label === "oxygen" ? 0 : 12));
      continue;
    }
    ctx.strokeStyle = cssVar("--muted"); ctx.beginPath();
    ctx.moveTo(X(pts[0][0]), Y(pts[0][1])); ctx.lineTo(X(pts[1][0]), Y(pts[1][1])); ctx.stroke();
    const [lx, ly] = l.p === 0 ? [span * 0.55, l.r / l.q] : [pts[0][0] * 0.55, (l.r - l.p * pts[0][0] * 0.55) / l.q];
    ctx.fillText(l.label, Math.min(X(lx) + 3, w - 80), Y(ly) - 4);
  }
  ctx.setLineDash([]);
  ctx.fillStyle = "#d85a30"; ctx.beginPath(); ctx.arc(X(a), Y(b), 5, 0, 2 * Math.PI); ctx.fill();
}

function drawInspector() {
  if (!state.selected) return;
  const [r, c] = state.selected, g = state.grid, cell = g.cells[r][c];
  const glucose = g.glucose[r][c], oxygen = g.oxygen[r][c], lactate = g.lactate[r][c];
  const rows = [["Square", `row ${r}, column ${c}`]];
  const fmt = (v, d = 2) => (Math.abs(v) >= 100 ? v.toFixed(0) : v.toFixed(d));
  rows.push(["Glucose here", fmt(glucose)], ["Oxygen here", fmt(oxygen)], ["Lactate here", fmt(lactate)]);
  if (g.vesselSet.has(r * SIZE + c) || cell.kind === M.EMPTY) {
    rows.unshift(["Kind", g.vesselSet.has(r * SIZE + c) ? "blood vessel" : "empty"]);
    $("inspect").hidden = false; $("inspectHint").hidden = true;
    $("inspectValues").innerHTML = rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join("");
    const { ctx, w, h } = setupCanvas($("plane")); ctx.clearRect(0, 0, w, h);
    return;
  }
  const f = state.params.uptakeFraction;
  const fluxes = M.cellFluxes(cell, glucose, oxygen, state.params, f);
  const a = fluxes.v_l / 2 / f, b = fluxes.v_m / 2 / f;
  rows.unshift(["Kind", cell.kind]);
  if (cell.mitoShare !== null) rows.push(["Budget on mitochondria", state.params.plastic ? "chosen each step" : `${(100 * cell.mitoShare).toFixed(0)}%`]);
  else rows.push(["Mitochondria M", fmt(cell.mitoCapacity)], ["Glucose uptake", fmt(cell.glucoseUptake)]);
  rows.push(["a: sugar fermented / step", fmt(a)], ["b: sugar respired / step", fmt(b)],
            ["ATP / step", fmt(M.atp(fluxes) / f, 1)], ["Lactate made / step", fmt(fluxes.v_out / f)],
            ["Warburg index a/(a+b)", M.warburgIndex(fluxes).toFixed(2)],
            ["Acid tolerance", fmt(cell.acidTolerance, 1)]);
  $("inspect").hidden = false; $("inspectHint").hidden = true;
  $("inspectValues").innerHTML = rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join("");
  drawPlane(a, b, cellLimits(cell, glucose, oxygen));
}

/** Switch the map to a view (cells, oxygen, acid). */
function setView(view) {
  state.view = view;
  document.querySelectorAll(".tabs button").forEach((x) => x.setAttribute("aria-selected", String(x.dataset.view === view)));
  drawLegend(); draw();
}

// ---- Wiring --------------------------------------------------------------------------------
function setPlaying(on) {
  state.playing = on;
  state.lastTime = null; state.owed = 0;
  $("play").textContent = on ? "Pause" : "Play";
}

/** Runs about stepsPerSecond steps each second, however fast the screen refreshes. */
function frame(time) {
  if (state.playing) {
    if (state.lastTime !== null) state.owed += ((time - state.lastTime) / 1000) * state.stepsPerSecond;
    state.lastTime = time;
    const steps = Math.min(Math.floor(state.owed), 4);
    if (steps > 0) {
      state.owed -= steps;
      for (let i = 0; i < steps; i++) doStep();
      draw();
      if (M.count(state.grid, M.TUMOR) === 0) setPlaying(false);
    }
  }
  requestAnimationFrame(frame);
}

function init() {
  document.querySelectorAll(".run").forEach((b) => b.addEventListener("click", () => chooseQuestion(b.dataset.name)));
  $("harmless").addEventListener("click", () => chooseQuestion("invasion", true));
  $("play").addEventListener("click", () => setPlaying(!state.playing));
  $("stepBtn").addEventListener("click", () => { setPlaying(false); doStep(); draw(); });
  $("reset").addEventListener("click", () => { state.seed++; setPlaying(false); startRun(); });
  ["tolerance", "glucose", "weight"].forEach((id) => $(id).addEventListener("input", applySettings));
  document.querySelectorAll(".speed button").forEach((b) => b.addEventListener("click", () => {
    state.stepsPerSecond = +b.dataset.speed;
    document.querySelectorAll(".speed button").forEach((x) => x.classList.toggle("active", x === b));
  }));
  document.querySelectorAll(".tabs button, .viewlink").forEach((b) => b.addEventListener("click", () => setView(b.dataset.view)));
  $("map").addEventListener("click", (e) => {
    const rect = e.target.getBoundingClientRect();
    const c = Math.floor(((e.clientX - rect.left) / rect.width) * SIZE), r = Math.floor(((e.clientY - rect.top) / rect.height) * SIZE);
    if (r >= 0 && r < SIZE && c >= 0 && c < SIZE) { state.selected = [r, c]; $("inspectBox").open = true; draw(); }
  });
  $("inspectBox").addEventListener("toggle", draw);
  window.addEventListener("resize", draw);
  state.selected = [SIZE / 2, SIZE / 2];
  drawLegend();
  chooseQuestion("invasion");
  requestAnimationFrame(frame);
}

init();
