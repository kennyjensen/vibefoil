import assert from 'node:assert/strict';
import { gauss, gauss1 } from '../js/xsolve.js';

// Reproduce the previous copying wrapper as an independent reference. Compare
// the destroyed matrix as well as the RHS, including all unused padding.
function reference(n, z, r) {
  const z0 = Array.from({ length: n }, (_, i) =>
    Float64Array.from({ length: n }, (_, j) => z[i + 1][j + 1]));
  const r0 = Float64Array.from({ length: n }, (_, i) => r[i + 1]);
  gauss(n, z0, r0);
  for (let i = 0; i < n; i += 1) {
    r[i + 1] = r0[i];
    for (let j = 0; j < n; j += 1) z[i + 1][j + 1] = z0[i][j];
  }
}

let seed = 0x23012412;
function random() {
  seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
  return seed / 0x100000000 - 0.5;
}

function makeSystem(n) {
  // Rows share storage, as the flat VS2 views do in the real solver.
  const size = n + 2;
  const flat = Float64Array.from({ length: size * size }, random);
  const z = Array.from({ length: size }, (_, i) => flat.subarray(i * size, (i + 1) * size));
  const r = Float64Array.from({ length: size }, random);
  return { z, r };
}

let cases = 0;
function compare(n, { z, r }) {
  const expected = structuredClone({ z, r });
  reference(n, expected.z, expected.r);
  gauss1(n, z, r);
  assert.deepStrictEqual(z, expected.z, `matrix, n=${n}, case=${cases}`);
  assert.deepStrictEqual(r, expected.r, `RHS, n=${n}, case=${cases}`);
  cases += 1;
}

for (const n of [1, 2, 3, 4, 8, 16, 31]) {
  for (let sample = 0; sample < 20; sample += 1) {
    const system = makeSystem(n);
    compare(n, system);
    compare(n, system); // Repeated solves must preserve the public layout.
  }
}

// Equal-magnitude pivots, an exchange with the last row, zero/singular systems,
// signed zeros, and nonfinite/extreme values must follow the original routine.
const specials = [0, -0, NaN, Infinity, -Infinity, Number.MIN_VALUE, Number.MAX_VALUE];
for (const value of specials) {
  for (const row of [1, 2, 4]) {
    for (const column of [1, 2, 4]) {
      const system = makeSystem(4);
      system.z[row][column] = value;
      compare(4, system);
    }
    const system = makeSystem(4);
    system.r[row] = value;
    compare(4, system);
  }
}
for (const column of [[2, -2, 2, -2], [0, 0, 0, 2], [0, 0, 0, 0]]) {
  const system = makeSystem(4);
  for (let i = 1; i <= 4; i += 1) system.z[i][1] = column[i - 1];
  compare(4, system);
}

console.log(`gauss1: ${cases} exact in-place/reference comparisons passed`);
