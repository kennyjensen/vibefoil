import assert from 'node:assert/strict';
import { makeContext } from './blsolv_fixture.js';
import { createWasmBlsolv } from '../js/blsolv_wasm.js';
import { blsolv } from '../js/xsolve.js';

// The nested-array path is an independent reference for the flat-array path
// used by the browser. Check all mutated coefficients, not only the solution.

const accelerated = createWasmBlsolv();
assert.equal(typeof accelerated, 'function', 'Tests must exercise the WASM kernel');
let acceleratedCases = 0;
let cases = 0;
function compare(ctx, expectedSuccess = true) {
  const wasm = structuredClone(ctx);
  const wasmSuccess = accelerated(wasm);
  assert.equal(wasmSuccess, expectedSuccess);
  if (wasmSuccess) acceleratedCases += 1;
  else assert.deepStrictEqual(wasm, ctx, 'Failed acceleration must not alter its input');
  const reference = structuredClone(ctx);
  for (const name of ['VA', 'VB', 'VM', 'VDEL']) delete reference[`${name}F`];
  const success = blsolv(ctx);
  assert.equal(success, expectedSuccess);
  if (success) assert.ok(ctx.VDELF.every(Number.isFinite), 'Solution must remain finite');
  assert.equal(success, blsolv(reference));
  for (const name of ['VA', 'VB', 'VM', 'VDEL', 'VZ']) {
    assert.deepStrictEqual(ctx[name], reference[name], `${name}, NSYS=${ctx.NSYS}`);
  }
  if (wasmSuccess) assert.deepStrictEqual(wasm, ctx, 'WASM must match every context field');
  cases += 1;
}

for (const nsys of [1, 2, 3, 4, 5, 6, 7, 8, 31, 32, 33, 34, 35, 48, 120, 240]) {
  for (const vaccel of [0, 0.0001, 0.01, 1]) {
    const ctx = makeContext(nsys, vaccel);
    compare(ctx);
    // Repeated calls must still expose the original public matrix layout.
    compare(ctx);
  }
}

// Exercise all eight acceleration masks at iv=1, kv=3. Values exactly at
// the cutoff must remain inactive; nonzero source terms expose any mismatch.
// Cover no four-entry block, exact blocks, and all three tail lengths.
for (const nsys of [4, 5, 6, 7, 8, 9, 32, 33, 34, 35]) {
  for (let mask = 0; mask < 8; mask += 1) {
    const ctx = makeContext(nsys, 0.125);
    ctx.S[1] = 0;
    ctx.S[2] = 2;
    for (const name of ['VAF', 'VBF', 'VMF', 'VDELF']) ctx[name].fill(0);
    for (const row of ctx.VZ) row.fill(0);
    for (let iv = 1; iv <= nsys; iv += 1) {
      ctx.VA[1][1][iv] = 1;
      ctx.VA[2][2][iv] = 1;
      ctx.VM[3][iv][iv] = 1;
    }
    for (let k = 1; k <= 3; k += 1) {
      ctx.VM[k][1][3] = mask & (1 << (k - 1)) ? 0.25 : 0.125;
      ctx.VDEL[k][1][1] = 0.2 * k;
      ctx.VDEL[k][2][1] = -0.3 * k;
    }
    for (let l = 2; l <= nsys; l += 1) {
      ctx.VM[3][l][1] = (l % 2 ? -1 : 1) * (0.01 + 0.01 * l);
    }
    compare(ctx);
  }
}

const warn = console.warn;
try {
  console.warn = () => {};
  for (const iv of [1, 4, 8]) {
    for (const value of [0, NaN, Infinity]) {
      for (const name of ['VA11', 'VA22', 'VM33']) {
        const ctx = makeContext(8, 0);
        if (name === 'VA11') ctx.VA[1][1][iv] = value;
        if (name === 'VA22') ctx.VA[2][2][iv] = value;
        if (name === 'VM33') ctx.VM[3][iv][iv] = value;
        // Earlier elimination changes VA22/VM33, so force exact zero pivots
        // only in the first block; nonfinite values remain singular later.
        if (value === 0 && iv > 1 && name !== 'VA11') continue;
        compare(ctx, false);
      }
    }
  }
} finally {
  console.warn = warn;
}
console.log(`blsolv: ${cases} exact flat/reference comparisons passed`);
console.log(`blsolv: ${acceleratedCases} successful WASM comparisons passed`);
