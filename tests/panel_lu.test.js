import assert from 'node:assert/strict';
import '../js/naca.js';
import { ludcmp, baksub } from '../js/xsolve.js';
import { createCachedLudcmp } from '../js/panel_lu.js';
import { createMatrix } from '../js/arrays.js';
import { scalc, segspl } from '../js/spline.js';
import { ncalc, apcalc, ggcalc } from '../js/xpanel.js';
import { tecalc } from '../js/xfoil.js';

let seed = 24120012;
function random() {
  seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
  return seed / 0x100000000 - 0.5;
}
function system(n, Row = Float64Array) {
  return { a: Array.from({ length: n + 2 }, () => Row.from({ length: n + 3 }, random)),
    indx: Int32Array.from({ length: n + 2 }, () => -17) };
}

let comparisons = 0;
function compare(factorize, n, actual, expectedHit) {
  const expected = structuredClone(actual);
  const before = actual.a.slice();
  const expectedBefore = expected.a.slice();
  ludcmp(n, expected.a, expected.indx);
  const hit = factorize(n, actual.a, actual.indx);
  if (expectedHit !== undefined) assert.equal(hit, expectedHit, `cache, case ${comparisons}`);
  assert.deepStrictEqual(actual, expected, `factorization and padding, case ${comparisons}`);
  assert.deepStrictEqual(actual.a.map(row => before.indexOf(row)),
    expected.a.map(row => expectedBefore.indexOf(row)), `row objects, case ${comparisons}`);
  for (const rhs of [new Float64Array(n), Float64Array.from({ length: n }, random)]) {
    const reference = rhs.slice();
    baksub(n, expected.a, expected.indx, reference);
    baksub(n, actual.a, actual.indx, rhs);
    assert.deepStrictEqual(rhs, reference, `back substitution, case ${comparisons}`);
  }
  comparisons += 1;
  return hit;
}

const factorize = createCachedLudcmp();
for (const n of [1, 2, 3, 4, 8, 31, 241]) {
  for (let sample = 0; sample < (n === 241 ? 2 : 8); sample += 1) {
    const input = system(n);
    const output = structuredClone(input);
    compare(factorize, n, output, false);
    // Cache storage must be independent of both the factors and pivot output.
    for (const row of output.a) row.fill(NaN);
    output.indx.fill(-1);
    const padded = structuredClone(input);
    for (let i = 0; i < n; i += 1) padded.a[i][n] = i + 123;
    padded.indx[n] = 987;
    compare(factorize, n, padded, true);
    compare(factorize, n, structuredClone(input), true);
    const changed = structuredClone(input);
    changed.a[n - 1][n - 1] += 0.000001;
    compare(factorize, n, structuredClone(changed), false);
    compare(factorize, n, structuredClone(changed), true);
    compare(factorize, n, structuredClone(input), false); // Only one entry is retained.
  }
}

// Distinguish signed zero even if it would produce the same factors here.
const zero = system(3);
zero.a[0][1] = 0;
compare(factorize, 3, structuredClone(zero), false);
zero.a[0][1] = -0;
compare(factorize, 3, structuredClone(zero), false);
compare(factorize, 3, structuredClone(zero), true);

for (const value of [0, -0, NaN, Infinity, -Infinity, Number.MIN_VALUE, Number.MAX_VALUE]) {
  for (const [i, j] of [[0, 0], [1, 3], [3, 1]]) {
    const input = system(4);
    input.a[i][j] = value;
    compare(factorize, 4, structuredClone(input), false);
    compare(factorize, 4, structuredClone(input));
  }
}
for (const column of [[2, -2, 2, -2], [0, 0, 0, 2], [0, 0, 0, 0]]) {
  const input = system(4);
  for (let i = 0; i < 4; i += 1) input.a[i][0] = column[i];
  compare(factorize, 4, structuredClone(input), false);
  compare(factorize, 4, structuredClone(input));
}

// Other layouts use LUDCMP directly, retaining Float32 rounding and aliases.
const fallbacks = [system(4, Float32Array), system(4, Array)];
const shared = system(4);
const flat = Float64Array.from({ length: 6 * 7 }, random);
shared.a = Array.from({ length: 6 }, (_, i) => flat.subarray(i * 7, (i + 1) * 7));
fallbacks.push(shared);
const overlap = system(4);
overlap.a[1] = overlap.a[0];
fallbacks.push(overlap);
const pivotAlias = system(4);
pivotAlias.indx = new Int32Array(pivotAlias.a[3].buffer, 0, 6);
fallbacks.push(pivotAlias);
for (const input of fallbacks) {
  compare(factorize, 4, structuredClone(input), false);
  compare(factorize, 4, structuredClone(input), false);
}
assert.throws(() => factorize(516, [], new Int32Array(516)),
  { message: 'LUDCMP: Array overflow. Increase NVX.' });

// Exercise the optional GGCALC hook with the real panel assembly, including
// reflected geometry (whose influence matrix also depends on angle).
function geometry(ctx) {
  const n = ctx.N;
  scalc(ctx.X, ctx.Y, ctx.S, n);
  segspl(ctx.X, ctx.XP, ctx.S, n);
  segspl(ctx.Y, ctx.YP, ctx.S, n);
  ncalc(ctx.X, ctx.Y, ctx.S, n, ctx.NX, ctx.NY);
  ctx.XTE = 0.5 * (ctx.X[0] + ctx.X[n - 1]);
  ctx.YTE = 0.5 * (ctx.Y[0] + ctx.Y[n - 1]);
  tecalc(ctx);
  apcalc(ctx);
}
function panel(ides, sharp, image) {
  const nside = 16;
  const x = new Float64Array(2 * nside), y = new Float64Array(2 * nside);
  const generate = ides < 10000 ? globalThis.Naca.naca4 : globalThis.Naca.naca5;
  const { nb: n } = generate(ides, new Float64Array(nside), new Float64Array(nside),
    new Float64Array(nside), nside, x, y, sharp);
  const ctx = { N: n, X: x, Y: y, PI: Math.PI, QOPI: 1 / (4 * Math.PI),
    HOPI: 1 / (2 * Math.PI), QINF: 1, ALFA: 0, LIMAGE: image, YIMAGE: 2,
    GAMU: createMatrix(n + 1, 2), QINVU: createMatrix(n + 1, 2),
    AIJ: createMatrix(n + 1, n + 1), BIJ: createMatrix(n + 1, n + 1),
    AIJPIV: new Int32Array(n + 1) };
  for (const name of ['S', 'XP', 'YP', 'NX', 'NY', 'APANEL', 'GAM', 'SIG',
    'DZDG', 'DZDN', 'DQDG', 'DZDM', 'DQDM']) ctx[name] = new Float64Array(n + 1);
  geometry(ctx);
  return ctx;
}
let panelComparisons = 0;
for (const ides of [12, 2412, 23012]) {
  for (const sharp of [false, true]) {
    for (const image of [false, true]) {
      const actual = panel(ides, sharp, image);
      const expected = structuredClone(actual);
      const cached = createCachedLudcmp();
      function check(expectedHit) {
        ggcalc(expected);
        let hit;
        ggcalc(actual, (...args) => { hit = cached(...args); });
        assert.equal(hit, expectedHit, `panel cache ${ides}, sharp=${sharp}, image=${image}`);
        assert.deepStrictEqual(actual, expected, `entire panel state ${panelComparisons}`);
        panelComparisons += 1;
      }
      check(false);
      check(true);
      actual.ALFA = expected.ALFA = 0.23;
      check(!image);
      actual.QINF = expected.QINF = 1.3;
      for (let i = 0; i < actual.N; i += 1) actual.SIG[i] = expected.SIG[i] = random();
      check(true);
      actual.Y[10] += 0.001;
      expected.Y[10] += 0.001;
      geometry(actual);
      geometry(expected);
      check(false);
      check(true);
    }
  }
}
console.log(`panel_lu: ${comparisons} exact matrix and ${panelComparisons} entire panel-state comparisons passed`);
