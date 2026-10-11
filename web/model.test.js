// model.test.js - checks for model.js. Run from the web/ folder:  node model.test.js
// Three kinds of checks:
//   1. hand-worked numbers (the same ones as the Python tests),
//   2. golden.json: answers computed by the Python model (python3 make_golden.py),
//      which the JS must match to 1e-9,
//   3. statistics for the random parts (dice, mutation), which can't match Python's
//      generator number for number.

const M = require("./model.js");
const golden = require("./golden.json");

let passed = 0;
let failed = 0;

/** Record one check: print PASS/FAIL with its name. */
function check(name, actual, expected) {
  const ok = JSON.stringify(actual) === JSON.stringify(expected);
  if (ok) passed++;
  else failed++;
  if (!ok) console.log(`FAIL  ${name}  (got ${JSON.stringify(actual)}, expected ${JSON.stringify(expected)})`);
}

/** Like check, but for decimals: within tol of expected. */
function checkClose(name, actual, expected, tol = 1e-9) {
  const ok = Math.abs(actual - expected) <= tol;
  if (ok) passed++;
  else failed++;
  if (!ok) console.log(`FAIL  ${name}  (got ${actual}, expected ${expected})`);
}

/** Every number in two nested arrays/objects agrees to tol (one check). */
function checkDeepClose(name, actual, expected, tol = 1e-9) {
  let worst = 0, where = "";
  const walk = (a, e, path) => {
    if (typeof e === "number") {
      const diff = Math.abs(a - e);
      if (!(diff <= worst)) { worst = diff; where = path; }
    } else if (e !== null && typeof e === "object") {
      for (const key of Object.keys(e)) walk(a == null ? undefined : a[key], e[key], `${path}.${key}`);
    } else if (a !== e) { worst = Infinity; where = path; }
  };
  walk(actual, expected, "");
  const ok = worst <= tol;
  if (ok) passed++;
  else failed++;
  if (!ok) console.log(`FAIL  ${name}  (largest difference ${worst} at ${where})`);
}

function checkThrows(name, fn) {
  let threw = false;
  try { fn(); } catch (e) { threw = true; }
  check(name, threw, true);
}

const section = (title) => console.log(`-- ${title}`);

// ---- 1. Linear algebra ------------------------------------------------------
section("linear algebra");
check("gcd(6, 10) = 2", M.gcd(6, 10), 2);
check("gcd(252, 198) = 18", M.gcd(252, 198), 18);
check("gcd(-12, 18) = 6", M.gcd(-12, 18), 6);
check("gcd(0, 0) = 0", M.gcd(0, 0), 0);
check("rank [[1,2],[3,4]] = 2", M.rank([[1, 2], [3, 4]]), 2);
check("rank [[1,2,3],[2,4,6]] = 1", M.rank([[1, 2, 3], [2, 4, 6]]), 1);
// integerVector makes the first nonzero entry positive: (-2, 1, 0) -> (2, -1, 0).
check("null space of x + 2y + 3z = 0", M.nullSpace([[1, 2, 3]]).map(M.integerVector), [[2, -1, 0], [3, 0, -1]]);
check("integer vector (1/2, 1/3) -> (3, 2)", M.integerVector([M.frac(1, 2), M.frac(1, 3)]), [3, 2]);
check("integer vector (6, -10) -> (3, -5)", M.integerVector([M.frac(6), M.frac(-10)]), [3, -5]);
checkThrows("zero vector has no direction", () => M.integerVector([M.frac(0), M.frac(0)]));
check("golden: rank of S", M.rank(M.stoichiometricMatrix()), golden.rank);
check("golden: null space of S", M.nullSpace(M.stoichiometricMatrix()).map(M.integerVector).sort(), golden.nullSpace);
check("golden: left null space of S", M.leftNullSpace(M.stoichiometricMatrix()).map(M.integerVector), golden.leftNullSpace);
check("RESPIRE from the null space", M.RESPIRE, [1, 1, 2, 0, 2, 0]);
check("FERMENT from the null space", M.FERMENT, [1, 1, 0, 2, 0, 2]);
{
  // Rank-nullity on random small integer matrices.
  const rng = M.makeRng(5);
  let ok = true;
  for (let t = 0; t < 200; t++) {
    const rows = 1 + Math.floor(rng.random() * 5), cols = 1 + Math.floor(rng.random() * 6);
    const a = Array.from({ length: rows }, () => Array.from({ length: cols }, () => Math.floor(rng.random() * 7) - 3));
    const basis = M.nullSpace(a);
    if (M.rank(a) + basis.length !== cols) ok = false;
    for (const v of basis) {
      const vec = M.integerVector(v);
      for (const row of a) if (row.reduce((s, x, j) => s + x * vec[j], 0) !== 0) ok = false;
    }
  }
  check("rank-nullity and A v = 0 on 200 random matrices", ok, true);
}

// ---- 2. Cell chemistry ------------------------------------------------------
section("cell chemistry");
{
  const f = M.solveCell(1, 10, 0.5);
  checkClose("Warburg cell ATP 9.5", M.atp(f), 9.5);
  checkClose("Warburg cell lactate 1.5", f.v_out, 1.5);
  checkClose("Warburg index 0.75", M.warburgIndex(f), 0.75);
  checkClose("full respiration 32 ATP", M.atp(M.solveCell(1, 100, 100)), 32);
  checkClose("full respiration uses 6 O2", M.oxygenUsed(M.solveCell(1, 100, 100)), 6);
  checkClose("no O2: 2 ATP", M.atp(M.solveCell(1, 0, 10)), 2);
  checkThrows("negative glucose throws", () => M.solveCell(-1, 1, 1));
  const s = M.stoichiometricMatrix();
  let steady = true;
  for (const c of golden.solveCell) {
    const f2 = M.solveCell(...c.in);
    for (const row of s) if (Math.abs(row.reduce((acc, x, j) => acc + x * f2[M.REACTIONS[j]], 0)) > 1e-12) steady = false;
  }
  check("every solveCell answer satisfies S v = 0", steady, true);
}
checkDeepClose(`golden: solveCell (${golden.solveCell.length} cases)`,
  golden.solveCell.map((c) => M.solveCell(...c.in)), golden.solveCell.map((c) => c.out));
checkDeepClose(`golden: pathwayAmounts (${golden.pathwayAmounts.length} cases)`,
  golden.pathwayAmounts.map((c) => M.pathwayAmounts(...c.in)), golden.pathwayAmounts.map((c) => c.out));
checkDeepClose("golden: budgetTraits", golden.budgetTraits.map((c) => M.budgetTraits(...c.in)), golden.budgetTraits.map((c) => c.out));
checkDeepClose(`golden: bestBudgetAmounts (${golden.bestBudgetAmounts.length} random cases)`,
  golden.bestBudgetAmounts.map((c) => M.bestBudgetAmounts(...c.in)), golden.bestBudgetAmounts.map((c) => c.out));
checkDeepClose("budget: plenty of glucose -> pure fermentation", M.bestBudgetAmounts(10, 10, 10, 10), [10, 0]);
checkDeepClose("budget: scarce glucose -> b = 0.45", M.bestBudgetAmounts(1, 10, 10, 10), [0.55, 0.45]);
check("solve2x2 parallel -> null", M.solve2x2(1, 2, 3, 2, 4, 5), null);

// ---- 3. Fields --------------------------------------------------------------
section("fields");
{
  const spike = M.makeField(5, 5);
  spike[2][2] = 1;
  const one = M.diffuse(spike, 0.5);
  checkClose("spike: middle keeps 1 - rate", one[2][2], 0.5);
  checkClose("spike: side gets 0.2 x rate", one[1][2], 0.1);
  checkClose("spike: corner gets 0.05 x rate", one[1][1], 0.025);
  let f = M.makeField(5, 5); f[0][0] = 3;
  for (let i = 0; i < 20; i++) f = M.diffuse(f, 1);
  checkClose("corner spike total conserved", M.total(f), 3);
  checkThrows("rate above 1 throws", () => M.diffuse(spike, 1.1));
  checkThrows("using more than there throws", () => M.applyChange([[1]], [[-1.5]]));
  check("refill sets only the sources", M.refillSources([[0.5, 0.5]], [[0, 1]], 9), [[0.5, 9]]);
}
for (const [i, c] of golden.diffuse.entries()) {
  let f = c.field;
  for (let s = 0; s < c.steps; s++) f = M.diffuse(f, c.rate);
  checkDeepClose(`golden: diffuse case ${i}`, f, c.out);
}

// ---- 4. Rules, metabolism, step ------------------------------------------------
section("rules and step");
check("neighbors of a middle square", M.neighbors(1, 1, 3, 3),
  [[0, 0], [0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1], [2, 2]]);
check("neighbors of a corner", M.neighbors(0, 0, 3, 3), [[0, 1], [1, 0], [1, 1]]);
{
  const h = M.makeCell(M.HEALTHY, 10, 1.0);
  check("acid exactly at tolerance survives", M.cellFate(h, 10, 1.0, 5, 20), M.STAY);
  check("acid above tolerance dies", M.cellFate(h, 1000, 1.01, 5, 20), M.DIE);
  check("ATP below 5 starves", M.cellFate(h, 4.99, 0, 5, 20), M.DIE);
  check("ATP at 20 divides", M.cellFate(h, 20, 0, 5, 20), M.DIVIDE);
  check("empty square does nothing", M.cellFate(M.emptyCell(), 100, 100, 5, 20), M.NOTHING);
}

/** Rebuild a grid from its golden JSON form. */
function gridFrom(json) {
  const cells = json.cells.map((row) => row.map((c) =>
    M.makeCell(c.kind, c.mitoCapacity, c.acidTolerance, c.glucoseUptake, c.mitoShare, c.marker)));
  return M.makeGrid(cells, json.oxygen, json.glucose, json.lactate, json.vessels);
}

for (const [i, c] of golden.metabolism.entries()) {
  const out = M.metabolism(gridFrom(c.grid), c.fraction, M.makeParams(c.params));
  checkDeepClose(`golden: metabolism case ${i} (${c.params.plastic ? "plastic" : "fixed/evolving"}, fraction ${c.fraction})`, out, c.out);
}
{
  let g = gridFrom(golden.step.grid);
  const params = M.makeParams(golden.step.params);
  const rng = M.makeRng(0);
  for (let s = 0; s < golden.step.steps; s++) g = M.step(g, params, rng);
  checkDeepClose("golden: 4 full steps on a 12x12 tissue (fields)", [g.oxygen, g.glucose, g.lactate],
    [golden.step.out.oxygen, golden.step.out.glucose, golden.step.out.lactate]);
  check("golden: 4 full steps (cells)", g.cells.map((row) => row.map((x) => x.kind)),
    golden.step.out.cells.map((row) => row.map((x) => x.kind)));
}

// ---- 5. Random parts (statistics) ----------------------------------------------
section("random parts");
{
  const rng = M.makeRng(1);
  const xs = Array.from({ length: 20000 }, () => rng.random());
  checkClose("uniform mean about 0.5", M.mean(xs), 0.5, 0.01);
  const gs = Array.from({ length: 20000 }, () => rng.gauss(0, 1));
  const m = M.mean(gs), sd = Math.sqrt(M.mean(gs.map((v) => (v - m) ** 2)));
  checkClose("gauss mean about 0", m, 0, 0.03);
  checkClose("gauss sd about 1", sd, 1, 0.03);
  const a = M.makeRng(42), b = M.makeRng(42);
  check("same seed, same numbers", [a.random(), a.random()], [b.random(), b.random()]);
}
/** 1x3 [left, empty, right], both dividing: fraction of trials the left one wins the middle. */
function contest(left, right, params, trials = 6000) {
  const cells = [[left, M.emptyCell(), right]];
  const rng = M.makeRng(7);
  let wins = 0;
  for (let t = 0; t < trials; t++) {
    const w = M.chooseDivisions(cells, [[0, 0], [0, 2]], new Set(), params, rng);
    if (w.get(1)[1] === 0) wins++;
  }
  return wins / trials;
}
{
  const tumor = M.warburgTumorCell(), healthy = M.healthyCell();
  checkClose("tumor beats healthy about 75%", contest(tumor, healthy, M.makeParams()), 0.75, 0.025);
  checkClose("same kind about 50%", contest(tumor, M.warburgTumorCell(), M.makeParams()), 0.5, 0.025);
  checkClose("tumor weight 1 -> about 50%", contest(tumor, healthy, M.makeParams({ tumorWeight: 1 })), 0.5, 0.025);
  check("never divides into a vessel",
    M.chooseDivisions([[healthy, M.emptyCell()]], [[0, 0]], new Set([1]), M.makeParams(), M.makeRng(0)).size, 0);
  const rng = M.makeRng(3);
  const shares = Array.from({ length: 6000 }, () => M.makeDaughter(M.budgetTumorCell(0.5), M.makeParams({ mutationSd: 0.05 }), rng).mitoShare);
  const sm = M.mean(shares), ssd = Math.sqrt(M.mean(shares.map((v) => (v - sm) ** 2)));
  checkClose("share mutation mean 0.5", sm, 0.5, 0.005);
  checkClose("share mutation sd 0.05", ssd, 0.05, 0.004);
  check("healthy cells never mutate",
    M.makeDaughter(M.healthyCell(), M.makeParams({ mutationSd: 1 }), rng).mitoCapacity, 10);
}

// ---- 6. Scenarios ---------------------------------------------------------------
section("scenarios");
{
  check("capillaries every 3 on 6x6", M.capillaries(6, 6, 3), [[1, 1], [1, 4], [4, 1], [4, 4]]);
  const { grid, params } = M.presetGrid("invasion", 60);
  check("invasion preset: 9 tumor seed cells", M.count(grid, M.TUMOR), 9);
  // On 30x30 the center (15, 15) is a capillary, so the seed skips it.
  check("30x30 seed skips the central vessel", M.count(M.presetGrid("invasion", 30).grid, M.TUMOR), 8);
  check("preset uses calibrated substeps", params.diffusionSubsteps, 10);
  const off = M.presetGrid("acidOff", 30);
  check("acid-off preset turns acid harm off", off.params.acidHarm, false);
  // Healthy tissue alone stays alive (the calibration guard from test_scenarios.py).
  const p = M.calibratedParams();
  let g = M.healthyTissue(12, 12, M.CAPILLARY_SPACING, p);
  const rng = M.makeRng(0);
  for (let s = 0; s < 40; s++) g = M.step(g, p, rng);
  check("healthy tissue alone stays alive for 40 steps", M.count(g, M.HEALTHY), 140);
}

console.log(`\n${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
