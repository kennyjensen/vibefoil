import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { createWasmPanel } from '../js/panel_wasm.js';
import { panelBase64, sourceHashes } from '../js/panel_wasm_data.js';
import { ludcmp, baksub } from '../js/xsolve.js';

for (const [path, digest] of Object.entries(sourceHashes)) {
  assert.equal(createHash('sha256').update(readFileSync(new URL(`../${path}`, import.meta.url))).digest('hex'), digest,
    `Rebuild panel WASM after changing ${path}`);
}
const module = new WebAssembly.Module(Buffer.from(panelBase64, 'base64'));
assert.deepEqual(WebAssembly.Module.imports(module), []);
const accelerator = createWasmPanel();
assert.ok(accelerator);
let seed = 0x24122301;
const random = () => {
  seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
  return seed / 0x100000000 - 0.5;
};
const system = (n, count = 9) => ({
  a: Array.from({length: n+2}, () => Float64Array.from({length: n+3}, random)),
  indx: new Int32Array(n+2).fill(-17),
  b: Array.from({length: n+2}, () => Float64Array.from({length: count+4}, random)),
});
let comparisons = 0;
function factorize(n, input, required = true) {
  const expected = structuredClone(input);
  const rows = input.a.slice(), expectedRows = expected.a.slice();
  ludcmp(n, expected.a, expected.indx);
  const success = accelerator.factorize(n, input.a, input.indx);
  if (required) assert.equal(success, true);
  if (!success) ludcmp(n, input.a, input.indx);
  assert.deepStrictEqual(input, expected);
  assert.deepEqual(input.a.map(row => rows.indexOf(row)), expected.a.map(row => expectedRows.indexOf(row)));
  comparisons++;
}
function substitute(n, input, first, count) {
  const before = structuredClone(input);
  const expected = new Float64Array(n*count);
  for (let k = 0; k < count; k++) {
    const column = expected.subarray(k*n, (k+1)*n);
    for (let i = 0; i < n; i++) column[i] = input.b[i][first+k];
    baksub(n, input.a, input.indx, column);
  }
  const result = accelerator.solveColumns(n, input.a, input.indx, input.b, first, count);
  assert.deepStrictEqual(result, expected);
  assert.deepStrictEqual(input, before, 'Input matrix, pivots and RHS remain unchanged');
  // The result must survive another call and scratch-memory growth.
  const other = system(n);
  factorize(n, other);
  assert.deepStrictEqual(result, expected);
  comparisons++;
}
for (const n of [1, 2, 3, 8, 31, 64, 241, 515]) {
  const input = system(n, 17);
  factorize(n, input);
  for (const count of [1, 3, 4, 5, 17]) substitute(n, input, 2, count);
}
const large = system(241, 240);
factorize(241, large);
substitute(241, large, 1, 240);

// Pivot ties, late pivots and signed zeros must preserve JS row identities.
for (const column of [[2, -2, 2, -2], [0, 0, 0, 2], [0, -0, 0, -0]]) {
  const input = system(4);
  for (let i = 0; i < 4; i++) input.a[i][0] = column[i];
  factorize(4, input, false);
}
for (const pivot of ['identity', 'last']) {
  for (const starts of [[0, 0, 0, 0], [1, 1, 1, 1], [0, 1, 2, 3], [8, 8, 8, 8], [8, 2, 0, 5]]) {
    const input = system(8);
    factorize(8, input);
    for (let i = 0; i < 8; i++) {
      input.indx[i] = pivot === 'last' ? 7 : i;
      for (let k = 0; k < 8; k++) if (i < starts[k%4]) input.b[i][k+1] = k%2 ? -0 : 0;
    }
    substitute(8, input, 1, 8);
  }
}

// Rejection must preserve all bytes, including inactive padding and NaN bits.
function snapshot(input) {
  const buffers = new Map();
  for (const array of [...input.a, input.indx, ...input.b]) {
    if (ArrayBuffer.isView(array) && array.buffer.byteLength) buffers.set(array.buffer, new Uint8Array(array.buffer).slice());
  }
  return () => { for (const [buffer, bytes] of buffers) assert.deepEqual(new Uint8Array(buffer), bytes); };
}
function rejected(input, method = 'factorize', n = 8, first = 0, count = 4, instance = accelerator) {
  const check = snapshot(input), rows = input.a.slice();
  const result = method === 'factorize' ? instance.factorize(n, input.a, input.indx)
    : instance.solveColumns(n, input.a, input.indx, input.b, first, count);
  assert.equal(result, method === 'factorize' ? false : null);
  check();
  input.a.forEach((row, i) => assert.equal(row, rows[i]));
  comparisons++;
}
for (const mutate of [
  v => { v.a[0] = Array.from(v.a[0]); },
  v => { v.a[0] = Float32Array.from(v.a[0]); },
  v => { v.a[0] = new (class extends Float64Array {})(v.a[0]); },
  v => { v.a[0] = new Float64Array(new SharedArrayBuffer(88)); },
  v => { v.a[0] = new Float64Array(new ArrayBuffer(88, {maxByteLength: 176})); },
  v => { v.a[0] = new Float64Array(7); },
  v => { v.a[1] = v.a[0]; },
  v => { v.indx = new Int32Array(v.a[0].buffer, 0, 10); },
  v => { v.indx = new Float64Array(10); },
  v => { Object.defineProperty(v.a[0], 'length', {value: 88}); },
  v => { Object.freeze(v.a); },
  v => { v.a[0].fill(0); },
  v => { new DataView(v.a[0].buffer).setBigUint64(0, 0x7ff8000000000001n, true); },
  v => { v.a[4][5] = Infinity; },
  v => { v.a[4][5] = -Infinity; },
  v => { v.a[0].fill(Number.MIN_VALUE); },
]) {
  const input = system(8);
  mutate(input);
  rejected(input);
}
for (const mutate of [
  v => { v.indx[1] = -1; },
  v => { v.indx[1] = 8; },
  v => { v.b[1] = v.b[0]; },
  v => { v.b[1] = v.a[0]; },
  v => { v.b[0] = new Float64Array(3); },
  v => { v.b[0][2] = NaN; },
  v => { v.b[0][2] = Infinity; },
  v => { v.a[7][7] = 0; },
  v => { v.a[7][7] = Number.MIN_VALUE;v.b[7][0] = Number.MAX_VALUE; },
]) {
  const input = system(8);
  ludcmp(8, input.a, input.indx);
  mutate(input);
  rejected(input, 'solveColumns');
}
for (const n of [0, -1, 8.5, 516, NaN, Infinity]) rejected(system(8), 'factorize', n);
for (const [first, count] of [[-1, 4], [0.5, 4], [0, 0], [0, 1025], [2048, 4]]) {
  rejected(system(8), 'solveColumns', 8, first, count);
}
let getterCalls = 0;
for (const method of ['factorize', 'solveColumns']) {
  const input = system(8);
  Object.defineProperty(input.a, '0', {get() { getterCalls++;throw new Error('Unexpected getter'); }});
  const before = getterCalls;
  assert.equal(method === 'factorize' ? accelerator.factorize(8, input.a, input.indx)
    : accelerator.solveColumns(8, input.a, input.indx, input.b, 0, 4), method === 'factorize' ? false : null);
  assert.equal(getterCalls, before);
  comparisons++;
}

// A trap or growth failure discards private work. No half-published rows/pivots.
const originalWasm = globalThis.WebAssembly;
try {
  for (const failure of ['kernel', 'grow']) {
    globalThis.WebAssembly = Object.assign(Object.create(originalWasm), {Instance: class {
      constructor(module) {
        const original = new originalWasm.Instance(module).exports;
        const fail = () => { throw new Error('injected failure'); };
        this.exports = {...original,
          ...(failure === 'kernel' ? {panel_lu: fail, panel_solve: fail}
            : {memory: {buffer: new ArrayBuffer(0), grow: fail}})};
      }
    }});
    const broken = createWasmPanel();
    assert.ok(broken);
    rejected(system(8), 'factorize', 8, 0, 4, broken);
    const input = system(8);
    ludcmp(8, input.a, input.indx);
    rejected(input, 'solveColumns', 8, 0, 4, broken);
  }
  globalThis.WebAssembly = undefined;
  assert.equal(createWasmPanel(), null);
  globalThis.WebAssembly = Object.assign(Object.create(originalWasm), {Module: class { constructor() { throw new Error('blocked'); } }});
  assert.equal(createWasmPanel(), null);
} finally { globalThis.WebAssembly = originalWasm; }
factorize(515, system(515));
assert.ok(accelerator.stats.factorizations > 0 && accelerator.stats.solves > 0 && accelerator.stats.fallback > 0);
console.log(`panel_wasm: ${comparisons} exact matrix/solution and failure checks passed`);
