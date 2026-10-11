// app.js - the page: presets, settings, animation loop, charts and the cell
// inspector. All the science is in model.js; this file only runs and draws it.

const M = window.WarburgModel;
const SIZE = 60;
const CELL_PX = 8;   // the cells canvas is drawn at 8 px per grid square

const PRESET_NOTES = {
  invasion: "Fermenting tumor, no evolution. Its acid clears the healthy tissue.",
  freeMito: "Tumor mitochondria mutate and cost nothing. Does fermenting survive?",
  budget: "Tumors split a protein budget between sugar transporters and mitochondria; the split evolves.",
  plastic: "Each cell picks its best split every step. Greedy: watch what its acid does.",
  acidOff: "Budget evolution, but acid is harmless. Can the tumor still invade?",
};
const MUTATION_SD = { meanMito: 0.1, meanShare: 0.05 };

const state = {
  preset: "invasion", grid: null, params: null, rng: null, seed: 1, stepCount: 0,
  playing: false, selected: null,
  history: { healthy: [], tumor: [], trait: [], warburg: [] },
};

const $ = (id) => document.getElementById(id);

// ---- Settings ------------------------------------------------------------------
function readSettings() {
  return {
    tolerance: +$("tolerance").value, glucose: +$("glucose").value, weight: +$("weight").value,
    mutation: $("mutation").checked, acid: $("acid").checked,
  };
}

/** Params overrides from the settings panel, on top of the preset's own. */
function overrides(settings) {
  const trait = M.PRESETS[state.preset].trait;
  return {
    vesselGlucose: settings.glucose, tumorWeight: settings.weight, acidHarm: settings.acid,
    mutationSd: settings.mutation ? MUTATION_SD[trait] : 0,
  };
}

function showSettings() {
  const s = readSettings();
  $("toleranceOut").textContent = s.tolerance.toFixed(1);
  $("glucoseOut").textContent = s.glucose.toFixed(1);
  $("weightOut").textContent = `${s.weight % 1 ? s.weight.toFixed(1) : s.weight} : 1`;
  $("speedOut").textContent = `${$("speed").value} step${$("speed").value === "1" ? "" : "s"} / frame`;
}

/** Apply settings to the running simulation (next step uses them). */
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

/** Choose a preset: its own defaults for mutation and acid, other settings back to defaults. */
function choosePreset(name) {
  state.preset = name;
  const preset = M.PRESETS[name];
  $("mutation").checked = (preset.params.mutationSd || 0) > 0;
  $("acid").checked = preset.params.acidHarm !== false;
  $("tolerance").value = 1; $("glucose").value = 5; $("weight").value = 3;
  $("traitName").textContent = preset.trait === "meanShare" ? "mean share of budget on mitochondria" : "mean mitochondrial capacity M";
  document.querySelectorAll("#presets button").forEach((b) => b.classList.toggle("active", b.dataset.name === name));
  state.seed = 1;
  setPlaying(false);
  startRun();
}

function record() {
  const s = M.summary(state.grid, state.params);
  const h = state.history;
  h.healthy.push(M.count(state.grid, M.HEALTHY));
  h.tumor.push(M.count(state.grid, M.TUMOR));
  h.trait.push(s[M.PRESETS[state.preset].trait]);
  h.warburg.push(s.meanWarburg);
}

function doStep() {
  state.grid = M.step(state.grid, state.params, state.rng);
  state.stepCount++;
  record();
}

// ---- Drawing ----------------------------------------------------------------------------
function hexToRgb(hex) { const n = parseInt(hex.slice(1), 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; }
const COLORS = { empty: hexToRgb("#f1efe8"), healthy: hexToRgb("#1d9e75"), vessel: hexToRgb("#a32d2d"),
                 ferment: hexToRgb("#d85a30"), respire: hexToRgb("#534ab7") };
const mix = (x, y, t) => x.map((v, i) => Math.round(v + (y[i] - v) * t));
const rgb = (c) => `rgb(${c[0]},${c[1]},${c[2]})`;

/** Warburg index of the tumor cell at (r, c), judged as the simulation eats. */
function cellIndex(r, c) {
  const g = state.grid;
  return M.warburgIndex(M.cellFluxes(g.cells[r][c], g.glucose[r][c], g.oxygen[r][c], state.params, state.params.uptakeFraction));
}

function drawCells() {
  const canvas = $("cells"), ctx = canvas.getContext("2d"), g = state.grid;
  canvas.width = SIZE * CELL_PX; canvas.height = SIZE * CELL_PX;
  for (let r = 0; r < SIZE; r++) {
    for (let c = 0; c < SIZE; c++) {
      const cell = g.cells[r][c];
      let color;
      if (g.vesselSet.has(r * SIZE + c)) color = COLORS.vessel;
      else if (cell.kind === M.HEALTHY) color = COLORS.healthy;
      else if (cell.kind === M.TUMOR) color = mix(COLORS.respire, COLORS.ferment, cellIndex(r, c));
      else color = COLORS.empty;
      ctx.fillStyle = rgb(color);
      ctx.fillRect(c * CELL_PX, r * CELL_PX, CELL_PX, CELL_PX);
    }
  }
  if (state.selected) {
    const [r, c] = state.selected;
    ctx.lineWidth = 2; ctx.strokeStyle = "#111";   // a tight box on the square's own edge
    ctx.strokeRect(c * CELL_PX + 1, r * CELL_PX + 1, CELL_PX - 2, CELL_PX - 2);
  }
}

function drawField(id, field, top, dark) {
  const ctx = $(id).getContext("2d"), img = ctx.createImageData(SIZE, SIZE);
  for (let r = 0; r < SIZE; r++) {
    for (let c = 0; c < SIZE; c++) {
      const t = Math.min(Math.max(field[r][c] / top, 0), 1), col = mix([255, 255, 255], dark, t);
      img.data.set([col[0], col[1], col[2], 255], 4 * (r * SIZE + c));
    }
  }
  ctx.putImageData(img, 0, 0);
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

function drawChart(id, series, top) {
  const { ctx, w, h } = setupCanvas($(id));
  const left = 38, right = 8, upper = 8, lower = 18, pw = w - left - right, ph = h - upper - lower;
  ctx.clearRect(0, 0, w, h);
  ctx.strokeStyle = cssVar("--line"); ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(left, upper); ctx.lineTo(left, upper + ph); ctx.lineTo(left + pw, upper + ph); ctx.stroke();
  ctx.fillStyle = cssVar("--muted"); ctx.font = "11px system-ui, sans-serif"; ctx.textAlign = "right";
  ctx.fillText(top >= 10 ? Math.round(top) : top.toFixed(2), left - 4, upper + 8);
  ctx.fillText("0", left - 4, upper + ph);
  const n = Math.max(...series.map((s) => s.values.length));
  ctx.textAlign = "left"; ctx.fillText("step 0", left, h - 4);
  ctx.textAlign = "right"; ctx.fillText(`step ${n - 1}`, left + pw, h - 4);
  for (const s of series) {
    ctx.strokeStyle = s.color; ctx.lineWidth = 2; ctx.beginPath();
    s.values.forEach((v, i) => {
      const x = left + (n > 1 ? (pw * i) / (n - 1) : 0), y = upper + ph * (1 - v / top);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    });
    ctx.stroke();
  }
}

function drawCharts() {
  const h = state.history, living = SIZE * SIZE - state.grid.vessels.length;
  drawChart("countsChart", [{ values: h.healthy, color: "#1d9e75" }, { values: h.tumor, color: "#d85a30" }], living);
  const top = Math.max(1, ...h.trait);
  drawChart("traitChart", [{ values: h.warburg, color: "#d85a30" }, { values: h.trait, color: "#534ab7" }], top);
}

function draw() {
  drawCells();
  drawField("oxygen", state.grid.oxygen, state.params.vesselOxygen, [24, 95, 165]);
  drawField("lactate", state.grid.lactate, 4, [153, 60, 29]);
  $("step").textContent = state.stepCount;
  $("healthy").textContent = M.count(state.grid, M.HEALTHY);
  $("tumor").textContent = M.count(state.grid, M.TUMOR);
  drawCharts();
  drawInspector();
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

// ---- Wiring --------------------------------------------------------------------------------
function setPlaying(on) {
  state.playing = on;
  $("play").textContent = on ? "Pause" : "Play";
}

function frame() {
  if (state.playing) {
    for (let i = 0; i < +$("speed").value; i++) doStep();
    draw();
    if (M.count(state.grid, M.TUMOR) === 0) setPlaying(false);
  }
  requestAnimationFrame(frame);
}

function init() {
  $("presets").innerHTML = Object.entries(M.PRESETS).map(([name, p]) =>
    `<div class="preset"><button data-name="${name}">${p.label}</button><span>${PRESET_NOTES[name]}</span></div>`).join("");
  document.querySelectorAll("#presets button").forEach((b) => b.addEventListener("click", () => choosePreset(b.dataset.name)));
  $("play").addEventListener("click", () => setPlaying(!state.playing));
  $("stepBtn").addEventListener("click", () => { setPlaying(false); doStep(); draw(); });
  $("reset").addEventListener("click", () => { state.seed++; setPlaying(false); startRun(); });
  ["tolerance", "glucose", "weight", "mutation", "acid"].forEach((id) => $(id).addEventListener("input", applySettings));
  $("speed").addEventListener("input", showSettings);
  $("cells").addEventListener("click", (e) => {
    const rect = e.target.getBoundingClientRect();
    const c = Math.floor(((e.clientX - rect.left) / rect.width) * SIZE), r = Math.floor(((e.clientY - rect.top) / rect.height) * SIZE);
    if (r >= 0 && r < SIZE && c >= 0 && c < SIZE) { state.selected = [r, c]; draw(); }
  });
  window.addEventListener("resize", () => { drawCharts(); drawInspector(); });
  choosePreset("invasion");
  state.selected = [SIZE / 2, SIZE / 2];
  draw();
  requestAnimationFrame(frame);
}

init();
