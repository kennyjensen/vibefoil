import assert from 'node:assert/strict';
import { ludcmp, baksub, baksubMany } from '../js/xsolve.js';

let seed = 0x24122301;
function random() {
  seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
  return seed / 0x100000000 - 0.5;
}
function system(n, count) {
  const a = Array.from({ length: n + 2 }, () => Float64Array.from({ length: n + 3 }, random));
  const indx = new Int32Array(n + 2).fill(-123);
  ludcmp(n, a, indx);
  const columns = Array.from({ length: count + 2 }, () => Float64Array.from({ length: n + 3 }, random));
  return { a, indx, columns };
}
let comparisons = 0;
function compare(n, count, actual) {
  const expected = structuredClone(actual);
  const rows = actual.a.slice();
  const columns = actual.columns.slice();
  for (let k = 0; k < count; k += 1) baksub(n, expected.a, expected.indx, expected.columns[k]);
  baksubMany(n, actual.a, actual.indx, actual.columns, count);
  assert.deepStrictEqual(actual, expected, `n=${n}, columns=${count}, case=${comparisons}`);
  // The coefficient matrix, pivots, padding, and inactive columns must retain
  // their values, and neither matrix nor RHS row identities may change.
  actual.a.forEach((row, i) => assert.equal(row, rows[i]));
  actual.columns.forEach((column, i) => assert.equal(column, columns[i]));
  comparisons += 1;
}
for (const n of [0, 1, 2, 3, 4, 8, 31, 241]) {
  for (const count of [0, 1, 2, 3, 4, 5, 7, 8, 9]) {
    for (let sample = 0; sample < 3; sample += 1) {
      const input = system(n, count);
      compare(n, count, input);
      compare(n, count, input); // Repeated calls with the previous solution.
      for (let k = 0; k < count; k += 1) {
        // Different leading-zero lengths, including an entirely zero RHS.
        input.columns[k].fill(k % 2 ? -0 : 0, 0, Math.floor(n * (k % 5) / 4));
      }
      compare(n, count, input);
    }
  }
}

// Control the pivot sequence to cover common and different first-nonzero
// positions, including the all-zero first row used by the Kutta condition.
for (const pivot of ['identity', 'last']) {
  for (const starts of [[0, 0, 0, 0], [1, 1, 1, 1], [0, 1, 2, 3], [8, 8, 8, 8], [8, 2, 0, 5]]) {
    const input = system(8, 8);
    for (let i = 0; i < 8; i += 1) input.indx[i] = pivot === 'last' ? 7 : i;
    for (let k = 0; k < 8; k += 1) {
      input.columns[k].fill(k % 2 ? -0 : 0, 0, starts[k % 4]);
    }
    compare(8, 8, input);
  }
}

// Subtractions involving zero, nonfinite, and extreme values cannot be
// dropped or reassociated. Check each lane of both batched and leftover RHSs.
for (const value of [0, -0, NaN, Infinity, -Infinity, Number.MIN_VALUE, Number.MAX_VALUE]) {
  for (const column of [0, 1, 2, 3, 4]) {
    for (const row of [0, 1, 7]) {
      const input = system(8, 5);
      input.columns[column][row] = value;
      compare(8, 5, input);
    }
  }
  for (const [row, column] of [[0, 0], [3, 0], [7, 7], [2, 6]]) {
    const input = system(8, 5);
    input.a[row][column] = value;
    compare(8, 5, input);
  }
}
for (const diagonal of [0, -0, Number.MIN_VALUE]) {
  const input = system(8, 8);
  for (let i = 0; i < 8; i += 1) input.a[i][i] = diagonal;
  compare(8, 8, input);
}

// Independent subarray views in one buffer are supported without changing
// inactive padding. Only overlapping views are excluded by the API contract.
const packed = system(8, 9);
const flat = Float64Array.from({ length: 11 * 11 }, random);
packed.columns = Array.from({ length: 11 }, (_, i) => flat.subarray(i * 11, (i + 1) * 11));
compare(8, 9, packed);
const defaultCount = system(8, 4);
const expected = structuredClone(defaultCount);
for (const rhs of expected.columns) baksub(8, expected.a, expected.indx, rhs);
baksubMany(8, defaultCount.a, defaultCount.indx, defaultCount.columns);
assert.deepStrictEqual(defaultCount, expected);
comparisons += 1;

console.log(`baksub_many: ${comparisons} exact batched/single-column comparisons passed`);
