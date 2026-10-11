// model.js - the Warburg simulation, ported line by line from the Python in
// warburg-sim/ (linalg.py, fba.py, fields.py, rules.py, metrics.py,
// scenarios.py). Same order of operations, so results match the Python to
// floating-point precision (checked against golden.json in model.test.js).
// Works in node (require) and in the browser (window.WarburgModel).

// ---- Cells -----------------------------------------------------------------
const EMPTY = "empty";
const HEALTHY = "healthy";
const TUMOR = "tumor";

/**
 * A cell. mitoShare null = fixed traits; a number in [0, 1] = budget cell.
 * marker = neutral trait that mutates but does nothing (control for chance).
 */
function makeCell(kind, mitoCapacity = 0, acidTolerance = 0, glucoseUptake = 1, mitoShare = null, marker = 0) {
  if (kind !== EMPTY && kind !== HEALTHY && kind !== TUMOR) throw new Error(`Unknown kind: ${kind}`);
  if (mitoShare !== null && !(mitoShare >= 0 && mitoShare <= 1)) throw new Error("mitoShare must be in [0, 1]");
  return { kind, mitoCapacity, acidTolerance, glucoseUptake, mitoShare, marker };
}

function emptyCell() { return makeCell(EMPTY); }

/** Same traits as the parent (mutation is makeDaughter). */
function daughterOf(p) {
  return makeCell(p.kind, p.mitoCapacity, p.acidTolerance, p.glucoseUptake, p.mitoShare, p.marker);
}

/** Python Params defaults; overrides replace any field. */
function makeParams(overrides = {}) {
  const p = {
    surviveAtp: 5, divideAtp: 20, oxygenRate: 1, glucoseRate: 0.5, lactateRate: 0.5,
    vesselOxygen: 10, vesselGlucose: 2, tumorWeight: 3, healthyWeight: 1,
    diffusionSubsteps: 1, uptakeFraction: 1, mutationSd: 0, budget: 20, mitoCost: 10,
    plastic: false, acidHarm: true, ...overrides,
  };
  if (p.diffusionSubsteps < 1) throw new Error("diffusionSubsteps must be at least 1");
  if (!(p.uptakeFraction > 0 && p.uptakeFraction <= 1)) throw new Error("uptakeFraction must be in (0, 1]");
  if (p.mutationSd < 0) throw new Error("mutationSd must be >= 0");
  return p;
}

// ---- Exact linear algebra with fractions (linalg.py) -----------------------
function gcd(a, b) {
  a = Math.abs(a); b = Math.abs(b);
  while (b !== 0) [a, b] = [b, a % b];
  return a;
}
function lcm(a, b) { return (a * b) / gcd(a, b); }

/** A fraction {n, d} in lowest terms with d > 0. */
function frac(n, d = 1) {
  if (d === 0) throw new Error("zero denominator");
  if (d < 0) { n = -n; d = -d; }
  const g = gcd(n, d) || 1;
  return { n: n / g, d: d / g };
}
const fsub = (x, y) => frac(x.n * y.d - y.n * x.d, x.d * y.d);
const fmul = (x, y) => frac(x.n * y.n, x.d * y.d);
const fdiv = (x, y) => frac(x.n * y.d, x.d * y.n);
const isZero = (x) => x.n === 0;

/** Reduced row echelon form. Returns {reduced, pivots}. */
function rref(matrix) {
  const m = matrix.map((row) => row.map((v) => frac(v)));
  const rows = m.length, cols = m[0].length, pivots = [];
  let row = 0;
  for (let col = 0; col < cols && row < rows; col++) {
    let p = -1;
    for (let r = row; r < rows; r++) if (!isZero(m[r][col])) { p = r; break; }
    if (p < 0) continue;
    [m[row], m[p]] = [m[p], m[row]];
    const scale = m[row][col];
    m[row] = m[row].map((x) => fdiv(x, scale));
    for (let r = 0; r < rows; r++) {
      if (r !== row && !isZero(m[r][col])) {
        const factor = m[r][col];
        m[r] = m[r].map((x, j) => fsub(x, fmul(factor, m[row][j])));
      }
    }
    pivots.push(col);
    row++;
  }
  return { reduced: m, pivots };
}

function rank(matrix) { return rref(matrix).pivots.length; }

/** Basis of {x : A x = 0}, one vector (of fractions) per free column. */
function nullSpace(matrix) {
  const { reduced, pivots } = rref(matrix);
  const cols = matrix[0].length, basis = [];
  for (let f = 0; f < cols; f++) {
    if (pivots.includes(f)) continue;
    const x = Array.from({ length: cols }, () => frac(0));
    x[f] = frac(1);
    pivots.forEach((p, r) => { x[p] = frac(-reduced[r][f].n, reduced[r][f].d); });
    basis.push(x);
  }
  return basis;
}

function transpose(matrix) { return matrix[0].map((_, j) => matrix.map((row) => row[j])); }
function leftNullSpace(matrix) { return nullSpace(transpose(matrix)); }

/** Smallest integer vector in the same direction (lcm of denominators, then gcd). */
function integerVector(fracs) {
  if (fracs.every(isZero)) throw new Error("the zero vector has no direction");
  const common = fracs.reduce((acc, x) => lcm(acc, x.d), 1);
  let ints = fracs.map((x) => (x.n * common) / x.d);
  const divisor = ints.reduce((acc, v) => gcd(acc, v), 0);
  ints = ints.map((v) => v / divisor);
  const first = ints.find((v) => v !== 0);
  return first < 0 ? ints.map((v) => -v) : ints;
}

// ---- The cell's chemistry (fba.py) ------------------------------------------
const METABOLITES = ["glucose", "pyruvate", "NAD+", "NADH", "lactate"];
const REACTIONS = ["v_in", "v_g", "v_e", "v_l", "v_m", "v_out"];

function stoichiometricMatrix() {
  return [
    [1, -1, 0, 0, 0, 0],   // glucose
    [0, 2, 0, -1, -1, 0],  // pyruvate
    [0, -2, 1, 1, 0, 0],   // NAD+
    [0, 2, -1, -1, 0, 0],  // NADH
    [0, 0, 0, 1, 0, -1],   // lactate
  ];
}

const ATP_PER_GLYCOLYSIS = 2.0, ATP_PER_NADH = 2.5, ATP_PER_PYRUVATE = 12.5;
const O2_PER_NADH = 0.5, O2_PER_PYRUVATE = 2.5;
const iv_e = REACTIONS.indexOf("v_e"), iv_l = REACTIONS.indexOf("v_l"), iv_m = REACTIONS.indexOf("v_m");

/** The null space's integer basis as [respire, ferment]. */
function steadyStatePathways() {
  const basis = nullSpace(stoichiometricMatrix()).map(integerVector);
  if (basis.length !== 2) throw new Error(`expected a 2-dimensional null space, got ${basis.length}`);
  return [basis.find((v) => v[iv_e] > 0), basis.find((v) => v[iv_e] === 0)];
}
const [RESPIRE, FERMENT] = steadyStatePathways();
const O2_PER_RESPIRED_GLUCOSE = O2_PER_NADH * RESPIRE[iv_e] + O2_PER_PYRUVATE * RESPIRE[iv_m];   // 6
const PYRUVATE_BURNED_PER_RESPIRED_GLUCOSE = RESPIRE[iv_m];                                      // 2
const ATP_FERMENT = 2.0, ATP_RESPIRE = 32.0;

/** [a, b] = [glucose fermented, glucose respired]; b = min(G, O/6, M/2). */
function pathwayAmounts(glucose, oxygen, mitoCapacity) {
  if (glucose < 0 || oxygen < 0 || mitoCapacity < 0) throw new Error("inputs must be >= 0");
  const b = Math.min(glucose, oxygen / O2_PER_RESPIRED_GLUCOSE, mitoCapacity / PYRUVATE_BURNED_PER_RESPIRED_GLUCOSE);
  return [glucose - b, b];
}

function fluxesFromAmounts(a, b) {
  const out = {};
  REACTIONS.forEach((name, i) => { out[name] = a * FERMENT[i] + b * RESPIRE[i]; });
  return out;
}

function solveCell(glucose, oxygen, mitoCapacity) {
  const [a, b] = pathwayAmounts(glucose, oxygen, mitoCapacity);
  return fluxesFromAmounts(a, b);
}

function atp(f) { return ATP_PER_GLYCOLYSIS * f.v_g + ATP_PER_NADH * f.v_e + ATP_PER_PYRUVATE * f.v_m; }
function oxygenUsed(f) { return O2_PER_NADH * f.v_e + O2_PER_PYRUVATE * f.v_m; }

/** Share of glucose fermented, a / (a + b); 0 if no glucose. */
function warburgIndex(f) {
  const a = f.v_l / FERMENT[iv_l], b = f.v_m / RESPIRE[iv_m];
  return a + b > 0 ? a / (a + b) : 0.0;
}

/** [glucoseUptake, mitoCapacity] for a cell spending mitoShare of its budget on mitochondria. */
function budgetTraits(mitoShare, budget, mitoCost) {
  if (!(mitoShare >= 0 && mitoShare <= 1)) throw new Error("mitoShare must be in [0, 1]");
  return [(1 - mitoShare) * budget, (mitoShare * budget) / mitoCost];
}

/** Cramer's rule for p1 x + q1 y = r1, p2 x + q2 y = r2; null if parallel. */
function solve2x2(p1, q1, r1, p2, q2, r2) {
  const det = p1 * q2 - p2 * q1;
  if (Math.abs(det) < 1e-12) return null;
  return [(r1 * q2 - r2 * q1) / det, (p1 * r2 - p2 * r1) / det];
}

/** A plastic cell's best split: the best corner of the feasible polygon, as [a, b]. */
function bestBudgetAmounts(glucose, oxygen, budget, mitoCost) {
  if (glucose < 0 || oxygen < 0) throw new Error("glucose and oxygen must be >= 0");
  if (budget <= 0 || mitoCost <= 0) throw new Error("budget and mitoCost must be positive");
  const cons = [[1.0, mitoCost * PYRUVATE_BURNED_PER_RESPIRED_GLUCOSE, budget], [1.0, 0.0, glucose],
                [0.0, 1.0, oxygen / O2_PER_RESPIRED_GLUCOSE], [-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]];
  let bestG = 0.0, bestB = 0.0, bestAtp = 0.0;
  for (let i = 0; i < cons.length; i++) {
    for (let j = i + 1; j < cons.length; j++) {
      const corner = solve2x2(...cons[i], ...cons[j]);
      if (corner === null) continue;
      const [g, b] = corner;
      if (cons.every(([p, q, r]) => p * g + q * b <= r + 1e-9)) {
        const value = ATP_FERMENT * (g - b) + ATP_RESPIRE * b;
        if (value > bestAtp + 1e-12) { bestG = g; bestB = b; bestAtp = value; }
      }
    }
  }
  const a = bestG - bestB, b = bestB;
  return [a > 0 ? a : 0.0, b > 0 ? b : 0.0];
}

// ---- Fields (fields.py) -----------------------------------------------------
const KERNEL = [[0.05, 0.2, 0.05], [0.2, -1.0, 0.2], [0.05, 0.2, 0.05]];

function makeField(rows, cols, value = 0.0) {
  return Array.from({ length: rows }, () => new Array(cols).fill(value));
}

/** One diffusion step: each square sends rate x weight x amount to each existing neighbor. */
function diffuse(field, rate) {
  if (!(rate >= 0 && rate <= 1)) throw new Error("rate must be between 0 and 1");
  const rows = field.length, cols = field[0].length;
  const out = field.map((row) => row.slice());
  for (let r = 0; r < rows; r++) {
    for (let c = 0; c < cols; c++) {
      for (let dr = -1; dr <= 1; dr++) {
        for (let dc = -1; dc <= 1; dc++) {
          if (dr === 0 && dc === 0) continue;
          const nr = r + dr, nc = c + dc;
          if (nr >= 0 && nr < rows && nc >= 0 && nc < cols) {
            const moved = rate * KERNEL[dr + 1][dc + 1] * field[r][c];
            out[r][c] -= moved;
            out[nr][nc] += moved;
          }
        }
      }
    }
  }
  for (const row of out) {
    for (let c = 0; c < row.length; c++) {
      if (row[c] < 0) {
        if (row[c] < -1e-9) throw new Error(`diffusion produced ${row[c]}`);
        row[c] = 0.0;
      }
    }
  }
  return out;
}

function refillSources(field, sources, level) {
  const out = field.map((row) => row.slice());
  for (const [r, c] of sources) out[r][c] = level;
  return out;
}

function applyChange(field, change) {
  return field.map((row, r) => row.map((amount, c) => {
    const value = amount + change[r][c];
    if (value < -1e-9) throw new Error(`amount went negative (${value})`);
    return Math.max(value, 0.0);
  }));
}

function total(field) { return field.reduce((s, row) => s + row.reduce((t, v) => t + v, 0), 0); }

// ---- Random numbers (seeded) -------------------------------------------------
/** mulberry32: a small seeded generator. rng.random() in [0, 1). */
function makeRng(seed) {
  let t = seed >>> 0;
  const random = () => {
    t = (t + 0x6d2b79f5) >>> 0;
    let x = Math.imul(t ^ (t >>> 15), 1 | t);
    x = (x + Math.imul(x ^ (x >>> 7), 61 | x)) ^ x;
    return ((x ^ (x >>> 14)) >>> 0) / 4294967296;
  };
  return {
    random,
    choice: (list) => list[Math.floor(random() * list.length)],
    weightedChoice: (list, weights) => {
      const totalWeight = weights.reduce((s, w) => s + w, 0);
      let x = random() * totalWeight;
      for (let i = 0; i < list.length; i++) { x -= weights[i]; if (x < 0) return list[i]; }
      return list[list.length - 1];
    },
    gauss: (mu, sd) => {   // Box-Muller
      const u = 1 - random(), v = random();
      return mu + sd * Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
    },
  };
}

// ---- Rules and one step (rules.py) -------------------------------------------
const DIE = "die", DIVIDE = "divide", STAY = "stay", NOTHING = "nothing";

function neighbors(row, col, rows, cols) {
  const out = [];
  for (let dr = -1; dr <= 1; dr++) {
    for (let dc = -1; dc <= 1; dc++) {
      if (dr === 0 && dc === 0) continue;
      const r = row + dr, c = col + dc;
      if (r >= 0 && r < rows && c >= 0 && c < cols) out.push([r, c]);
    }
  }
  return out;
}

function cellFate(cell, cellAtp, lactate, surviveAtp, divideAtp) {
  if (divideAtp < surviveAtp) throw new Error("divideAtp must be at least surviveAtp");
  if (cell.kind === EMPTY) return NOTHING;
  if (lactate > cell.acidTolerance) return DIE;
  if (cellAtp < surviveAtp) return DIE;
  if (cellAtp >= divideAtp) return DIVIDE;
  return STAY;
}

/** How one living cell metabolizes for `fraction` of a step (fixed / evolving share / plastic). */
function cellFluxes(cell, glucoseHere, oxygenHere, params, fraction = 1.0) {
  if (cell.mitoShare === null) {
    return solveCell(Math.min(cell.glucoseUptake * fraction, glucoseHere), oxygenHere, cell.mitoCapacity * fraction);
  }
  if (params.plastic) {
    const [a, b] = bestBudgetAmounts(glucoseHere, oxygenHere, params.budget * fraction, params.mitoCost);
    return fluxesFromAmounts(a, b);
  }
  const [uptake, mito] = budgetTraits(cell.mitoShare, params.budget, params.mitoCost);
  return solveCell(Math.min(uptake * fraction, glucoseHere), oxygenHere, mito * fraction);
}

/** [atpMade, glucoseChange, oxygenChange, lactateChange] fields. */
function metabolism(grid, fraction, params) {
  const rows = grid.cells.length, cols = grid.cells[0].length;
  const atpMade = makeField(rows, cols), dg = makeField(rows, cols), dox = makeField(rows, cols), dl = makeField(rows, cols);
  for (let r = 0; r < rows; r++) {
    for (let c = 0; c < cols; c++) {
      const cell = grid.cells[r][c];
      if (cell.kind === EMPTY) continue;
      const f = cellFluxes(cell, grid.glucose[r][c], grid.oxygen[r][c], params, fraction);
      atpMade[r][c] = atp(f);
      dg[r][c] = -f.v_in;
      dox[r][c] = -oxygenUsed(f);
      dl[r][c] = f.v_out;
    }
  }
  return [atpMade, dg, dox, dl];
}

function claimWeight(cell, params) { return cell.kind === TUMOR ? params.tumorWeight : params.healthyWeight; }

/** {target "r,c": parent [r, c]}: each divider claims a random empty neighbor; weighted dice per square. */
function chooseDivisions(cells, dividers, vesselSet, params, rng) {
  const rows = cells.length, cols = cells[0].length, claims = new Map();
  for (const [r, c] of dividers) {
    const options = neighbors(r, c, rows, cols).filter(([nr, nc]) =>
      cells[nr][nc].kind === EMPTY && !vesselSet.has(nr * cols + nc));
    if (options.length) {
      const [tr, tc] = rng.choice(options);
      const key = tr * cols + tc;
      if (!claims.has(key)) claims.set(key, []);
      claims.get(key).push([r, c]);
    }
  }
  const winners = new Map();
  for (const key of [...claims.keys()].sort((x, y) => x - y)) {
    const claimants = claims.get(key);
    winners.set(key, rng.weightedChoice(claimants, claimants.map(([r, c]) => claimWeight(cells[r][c], params))));
  }
  return winners;
}

/** Daughter of a division: tumor cells mutate (share or M) plus the neutral marker. */
function makeDaughter(parent, params, rng) {
  const d = daughterOf(parent);
  if (parent.kind === TUMOR && params.mutationSd > 0) {
    const noise = rng.gauss(0, params.mutationSd);
    if (parent.mitoShare !== null) d.mitoShare = Math.min(1, Math.max(0, parent.mitoShare + noise));
    else d.mitoCapacity = Math.max(0, parent.mitoCapacity + noise);
    d.marker = parent.marker + rng.gauss(0, params.mutationSd);
  }
  return d;
}

/** A grid: cells (rows of cells), oxygen / glucose / lactate fields, vessels [[r, c], ...]. */
function makeGrid(cells, oxygen, glucose, lactate, vessels) {
  const cols = cells[0].length;
  for (const [r, c] of vessels) if (cells[r][c].kind !== EMPTY) throw new Error(`vessel ${r},${c} holds a cell`);
  return { cells, oxygen, glucose, lactate, vessels, vesselSet: new Set(vessels.map(([r, c]) => r * cols + c)) };
}

function count(grid, kind) {
  let n = 0;
  for (const row of grid.cells) for (const cell of row) if (cell.kind === kind) n++;
  return n;
}

/**
 * One step: k substeps of (refill vessels, eat 1/k-sized share, diffuse); ATP = average
 * rate; fates decided at once; deaths, then divisions (with mutation).
 * Also returns the step's ATP per square as grid.atp (for the inspector).
 */
function step(grid, params, rng) {
  const rows = grid.cells.length, cols = grid.cells[0].length;
  const k = params.diffusionSubsteps, f = params.uptakeFraction;
  let oxygen = grid.oxygen, glucose = grid.glucose, lactate = grid.lactate;
  let atpMade = makeField(rows, cols);
  for (let s = 0; s < k; s++) {
    oxygen = refillSources(oxygen, grid.vessels, params.vesselOxygen);
    glucose = refillSources(glucose, grid.vessels, params.vesselGlucose);
    lactate = refillSources(lactate, grid.vessels, 0.0);
    const [part, dg, dox, dl] = metabolism({ cells: grid.cells, oxygen, glucose, lactate }, f, params);
    atpMade = atpMade.map((row, r) => row.map((v, c) => v + part[r][c]));
    oxygen = applyChange(oxygen, dox);
    glucose = applyChange(glucose, dg);
    lactate = applyChange(lactate, dl);
    oxygen = diffuse(oxygen, params.oxygenRate);
    glucose = diffuse(glucose, params.glucoseRate);
    lactate = diffuse(lactate, params.lactateRate);
  }
  atpMade = atpMade.map((row) => row.map((v) => v / (k * f)));

  const cells = grid.cells.map((row) => row.slice());
  const dividers = [];
  for (let r = 0; r < rows; r++) {
    for (let c = 0; c < cols; c++) {
      const acid = params.acidHarm ? lactate[r][c] : 0.0;
      const fate = cellFate(grid.cells[r][c], atpMade[r][c], acid, params.surviveAtp, params.divideAtp);
      if (fate === DIE) cells[r][c] = emptyCell();
      else if (fate === DIVIDE) dividers.push([r, c]);
    }
  }
  const winners = chooseDivisions(cells, dividers, grid.vesselSet, params, rng);
  for (const [key, [pr, pc]] of winners) {
    cells[Math.floor(key / cols)][key % cols] = makeDaughter(cells[pr][pc], params, rng);
  }
  const next = makeGrid(cells, oxygen, glucose, lactate, grid.vessels);
  next.atp = atpMade;
  return next;
}

// ---- Metrics (metrics.py) -------------------------------------------------------
function mean(values) { return values.length ? values.reduce((s, v) => s + v, 0) / values.length : 0.0; }

function summary(grid, params) {
  const mito = [], share = [], marker = [], index = [];
  grid.cells.forEach((row, r) => row.forEach((cell, c) => {
    if (cell.kind !== TUMOR) return;
    mito.push(cell.mitoCapacity);
    if (cell.mitoShare !== null) share.push(cell.mitoShare);
    marker.push(cell.marker);
    index.push(warburgIndex(cellFluxes(cell, grid.glucose[r][c], grid.oxygen[r][c], params, params.uptakeFraction)));
  }));
  return { meanMito: mean(mito), meanShare: mean(share), meanMarker: mean(marker), meanWarburg: mean(index) };
}

// ---- Scenarios (scenarios.py, main.py) ---------------------------------------------
const CAPILLARY_SPACING = 6;

function calibratedParams(overrides = {}) {
  return makeParams({ diffusionSubsteps: 10, uptakeFraction: 0.02, vesselOxygen: 10, vesselGlucose: 5, ...overrides });
}

function capillaries(rows, cols, spacing) {
  if (spacing < 1) throw new Error("spacing must be at least 1");
  const out = [], start = Math.floor(spacing / 2);
  for (let r = start; r < rows; r += spacing) for (let c = start; c < cols; c += spacing) out.push([r, c]);
  return out;
}

function healthyCell(acidTolerance = 1.0) { return makeCell(HEALTHY, 10, acidTolerance, 1); }
function warburgTumorCell() { return makeCell(TUMOR, 0.5, 3.0, 8.0); }
function budgetTumorCell(share = 0.5) { return makeCell(TUMOR, 0, 3.0, 1, share); }

function healthyTissue(rows, cols, spacing, params, healthyTolerance = 1.0) {
  const vessels = capillaries(rows, cols, spacing);
  const vset = new Set(vessels.map(([r, c]) => r * cols + c));
  const cells = Array.from({ length: rows }, (_, r) => Array.from({ length: cols }, (_, c) =>
    (vset.has(r * cols + c) ? emptyCell() : healthyCell(healthyTolerance))));
  return makeGrid(cells, makeField(rows, cols, params.vesselOxygen), makeField(rows, cols, params.vesselGlucose),
                  makeField(rows, cols), vessels);
}

/** Replace a (2 radius + 1) square around (row, col) with tumor cells, skipping vessels. */
function seedTumor(grid, row, col, radius, maker) {
  const cols = grid.cells[0].length;
  for (let r = row - radius; r <= row + radius; r++) {
    for (let c = col - radius; c <= col + radius; c++) {
      if (!grid.vesselSet.has(r * cols + c)) grid.cells[r][c] = maker();
    }
  }
}

/** The scenarios from main.py / experiments.py. trait: which mean to chart. */
const PRESETS = {
  invasion: { label: "Warburg invasion", params: {}, maker: warburgTumorCell, trait: "meanMito" },
  freeMito: { label: "Free mitochondria evolve", params: { mutationSd: 0.1 }, maker: warburgTumorCell, trait: "meanMito" },
  budget: { label: "Budget evolution", params: { mutationSd: 0.05 }, maker: () => budgetTumorCell(0.5), trait: "meanShare" },
  plastic: { label: "Plastic cells", params: { plastic: true }, maker: () => budgetTumorCell(0.5), trait: "meanShare" },
  acidOff: { label: "Acid off", params: { mutationSd: 0.05, acidHarm: false }, maker: () => budgetTumorCell(0.5), trait: "meanShare" },
};

/** Healthy tissue + 3x3 tumor seed in the middle, for a preset (plus slider overrides). */
function presetGrid(name, size = 60, overrides = {}, healthyTolerance = 1.0) {
  const preset = PRESETS[name];
  const params = calibratedParams({ ...preset.params, ...overrides });
  const grid = healthyTissue(size, size, CAPILLARY_SPACING, params, healthyTolerance);
  seedTumor(grid, Math.floor(size / 2), Math.floor(size / 2), 1, preset.maker);
  return { grid, params };
}

const WarburgModel = {
  EMPTY, HEALTHY, TUMOR, DIE, DIVIDE, STAY, NOTHING, METABOLITES, REACTIONS, RESPIRE, FERMENT, KERNEL,
  makeCell, emptyCell, daughterOf, makeParams, gcd, lcm, frac, rref, rank, nullSpace, leftNullSpace, transpose,
  integerVector, stoichiometricMatrix, pathwayAmounts, fluxesFromAmounts, solveCell, atp, oxygenUsed,
  warburgIndex, budgetTraits, solve2x2, bestBudgetAmounts, makeField, diffuse, refillSources, applyChange,
  total, makeRng, neighbors, cellFate, cellFluxes, metabolism, chooseDivisions, makeDaughter, makeGrid, count,
  step, mean, summary, calibratedParams, capillaries, healthyCell, warburgTumorCell, budgetTumorCell,
  healthyTissue, seedTumor, PRESETS, presetGrid, CAPILLARY_SPACING,
};
if (typeof module !== "undefined" && module.exports) module.exports = WarburgModel;
else window.WarburgModel = WarburgModel;
