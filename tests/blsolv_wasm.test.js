import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { once } from 'node:events';
import { readFileSync } from 'node:fs';
import { Worker } from 'node:worker_threads';
import { createWasmBlsolv } from '../js/blsolv_wasm.js';
import { blsolvBase64, sourceSha256 } from '../js/blsolv_wasm_data.js';
import { blsolv } from '../js/xsolve.js';
import { makeContext } from './blsolv_fixture.js';

const source = readFileSync(new URL('../wasm/blsolv_kernel.cpp', import.meta.url));
assert.equal(createHash('sha256').update(source).digest('hex'), sourceSha256,
  'Rebuild the generated module after changing the kernel');
const module = new WebAssembly.Module(Buffer.from(blsolvBase64, 'base64'));
assert.deepEqual(WebAssembly.Module.imports(module), [], 'Kernel needs no external runtime');
const solve = createWasmBlsolv();
assert.equal(typeof solve, 'function');
let cases = 0;

function unchanged(ctx) {
  const before = structuredClone(ctx);
  assert.equal(solve(ctx), false);
  assert.deepStrictEqual(ctx, before);
  for (const name of ['VAF', 'VBF', 'VMF', 'VDELF']) {
    const bytes = value => new Uint8Array(value.buffer, value.byteOffset, value.byteLength);
    if (ArrayBuffer.isView(ctx[name])) assert.deepStrictEqual(bytes(ctx[name]), bytes(before[name]));
  }
  cases += 1;
}

// Reject layouts for which snapshotting would change aliasing or rounding.
for (const change of [
  ctx => { delete ctx.VAF; },
  ctx => { ctx.VAF = Array.from(ctx.VAF); },
  ctx => { ctx.VAF = Float32Array.from(ctx.VAF); },
  ctx => { ctx.VAF = new (class extends Float64Array {})(ctx.VAF); },
  ctx => { ctx.VAF = ctx.VAF.subarray(1); },
  ctx => { ctx.VBF = ctx.VAF; },
  ctx => { ctx.VAF = new Float64Array(new SharedArrayBuffer(ctx.VAF.byteLength)); },
  ctx => { ctx.S = Array.from(ctx.S); },
  ctx => { ctx.S = ctx.VAF.subarray(0, 3); },
  ctx => { ctx.IBLTE = new Int32Array(ctx.VAF.buffer, 0, 3); },
  ctx => { ctx.ISYS[1] = ctx.VMF.subarray(0, 3); },
  ctx => { ctx.ISYS[3] = ctx.VAF.subarray(0, 3); },
  ctx => { ctx.VZ[1] = ctx.VBF.subarray(0, 3); },
  ctx => { ctx.VZ[2] = Float32Array.from(ctx.VZ[2]); },
  ctx => { ctx.VZ[3] = new Float64Array(2); },
  ctx => { ctx.ISYS[1][1] = 1.5; },
  ctx => { ctx.ISYS[3][2] = 0; },
  ctx => { ctx.ISYS[3][2] = 41; },
  ctx => { ctx.N = 100; },
  ctx => { ctx.NSYS = 0; },
  ctx => { ctx.NSYS = 1025; },
  ctx => { ctx.NSYS = 40.5; },
  ctx => { ctx.VACCEL = NaN; },
  ctx => { ctx.S[2] = Infinity; },
]) {
  const ctx = makeContext(40, 0.001);
  change(ctx);
  // structuredClone drops subclass prototypes, which are irrelevant to bytes.
  if (ctx.VAF instanceof Float64Array && Object.getPrototypeOf(ctx.VAF) !== Float64Array.prototype) {
    const before = ctx.VAF.slice();
    assert.equal(solve(ctx), false);
    assert.deepEqual(Array.from(ctx.VAF), Array.from(before));
    cases += 1;
  } else unchanged(ctx);
}

// Failed solves must restart JS with the original inputs, including padding
// and NaN payloads. Compare its partial writes with the independent nested path.
const warn = console.warn;
try {
  console.warn = () => {};
  for (const iv of [1, 20, 40]) {
    for (const value of [0, NaN, Infinity, -Infinity]) {
      for (const name of ['VA11', 'VA22', 'VM33']) {
        if (value === 0 && iv > 1 && name !== 'VA11') continue;
        const ctx = makeContext(40, 0);
        if (name === 'VA11') ctx.VA[1][1][iv] = value;
        if (name === 'VA22') ctx.VA[2][2][iv] = value;
        if (name === 'VM33') ctx.VM[3][iv][iv] = value;
        unchanged(ctx);
        const reference = structuredClone(ctx);
        for (const name of ['VA', 'VB', 'VM', 'VDEL']) delete reference[`${name}F`];
        assert.equal(blsolv(ctx), blsolv(reference));
        for (const name of ['VA', 'VB', 'VM', 'VDEL', 'VZ']) {
          assert.deepStrictEqual(ctx[name], reference[name]);
        }
      }
    }
  }
  const overflow = makeContext(40, 0);
  overflow.VA[1][1][1] = Number.MIN_VALUE;
  unchanged(overflow);
  const payload = makeContext(40, 0);
  new BigUint64Array(payload.VAF.buffer)[0] = 0x7ff0000000000123n;
  unchanged(payload);
} finally {
  console.warn = warn;
}

// Vary dimensions in one instance to exercise memory growth and view renewal.
// Include shifted views, signed zeros, tiny values, and unused wake metadata.
for (const n of [32, 240, 33, 320, 1, 48]) {
  const ctx = makeContext(n, 0.0001);
  const backing = new Float64Array(ctx.VAF.length + 8);
  const shifted = backing.subarray(4, 4 + ctx.VAF.length);
  shifted.set(ctx.VAF);
  ctx.VAF = shifted;
  const reference = structuredClone(ctx);
  // Flat inputs are authoritative; nested VA now deliberately retains its old
  // backing array. Both implementations must leave those old views untouched.
  for (const name of ['VAF', 'VMF', 'VDELF']) {
    ctx[name][0] = reference[name][0] = -0;
    ctx[name][1] = reference[name][1] = Number.MIN_VALUE;
  }
  ctx.ISYS[1][1] = reference.ISYS[1][1] = n;
  delete ctx.ISYS[3];
  delete reference.ISYS[3];
  // Compute the reference in a worker with WASM disabled below.
  const worker = new Worker(`
    import { parentPort, workerData } from 'node:worker_threads';
    globalThis.WebAssembly = undefined;
    const { blsolv } = await import(workerData.module);
    const success = blsolv(workerData.ctx);
    parentPort.postMessage({ success, ctx: workerData.ctx });
  `, { eval: true, type: 'module', workerData: {
    module: new URL('../js/xsolve.js', import.meta.url).href, ctx: reference,
  } });
  try {
    const [result] = await once(worker, 'message');
    assert.equal(result.success, true);
    assert.equal(solve(ctx), true);
    assert.deepStrictEqual(ctx, result.ctx);
    for (const name of ['VAF', 'VBF', 'VMF', 'VDELF']) {
      assert.deepStrictEqual(new Uint8Array(ctx[name].buffer), new Uint8Array(result.ctx[name].buffer));
    }
    cases += 1;
  } finally {
    await worker.terminate();
  }
}

const wasm = globalThis.WebAssembly;
try {
  globalThis.WebAssembly = undefined;
  assert.equal(createWasmBlsolv(), null);
  globalThis.WebAssembly = { Module() { throw new Error('Compilation blocked'); } };
  assert.equal(createWasmBlsolv(), null);
  globalThis.WebAssembly = {
    Module: wasm.Module,
    Instance: function () {
      return { exports: {
        memory: { buffer: new ArrayBuffer(0), grow() { throw new RangeError('Memory limit'); } },
        __heap_base: { value: 0 },
        blsolv_flat() { assert.fail('Must not solve after failed memory growth'); },
      } };
    },
  };
  const limited = createWasmBlsolv();
  const ctx = makeContext(40, 0);
  const before = structuredClone(ctx);
  assert.equal(limited(ctx), false);
  assert.equal(limited(ctx), false); // Retrying must not reuse invalid views.
  assert.deepStrictEqual(ctx, before);
} finally {
  globalThis.WebAssembly = wasm;
}

// The site's complete solver must also work without WebAssembly, including
// subsequent requests and saved solution reuse in the same worker.
const workerSource = `
  import { parentPort, workerData } from 'node:worker_threads';
  console.log = () => {};
  console.warn = () => {};
  if (workerData.disabled) globalThis.WebAssembly = undefined;
  globalThis.self = { postMessage(payload, transfers = []) {
    parentPort.postMessage(payload, transfers);
  } };
  await import(workerData.module);
  parentPort.on('message', data => self.onmessage({ data }));
  parentPort.postMessage({ ready: true });
`;
const workers = [false, true].map(disabled => new Worker(workerSource, {
  eval: true, type: 'module', workerData: {
    disabled, module: new URL('../js/solver_worker.js', import.meta.url).href,
  },
}));
try {
  await Promise.all(workers.map(worker => once(worker, 'message')));
  let states = [null, null];
  for (let requestId = 1; requestId <= 3; requestId += 1) {
    const results = [];
    for (const i of [0, 1]) {
      const response = once(workers[i], 'message');
      const alphaDeg = requestId + 1;
      workers[i].postMessage({ requestId, settings: {
        source: 'naca', mode: '4', m: 2, p: 4, t: 12,
        geometryKey: 'wasm-fallback-naca2412', advancedMode: false,
        reusePanel: true, reuseSolution: false, viscous: true,
        reuseState: states[i], alphaDeg, alphaRad: alphaDeg * Math.PI / 180,
        mach: 0.1, reynolds: 3e6, ncr: 9, nIter: 30,
        flap: { x: 0.75, yrel: 0.5, deflection: 0 }, canvasWidth: 800, canvasHeight: 500,
      } });
      [results[i]] = await response;
      assert.equal(results[i].ok, true, results[i].error);
      assert.equal(results[i].converged, true);
      states[i] = results[i].reuseState;
    }
    assert.deepStrictEqual(results[0], results[1], 'Complete WASM/JS worker results');
    cases += 1;
  }
} finally {
  await Promise.all(workers.map(worker => worker.terminate()));
}
console.log(`blsolv_wasm: ${cases} storage/failure/byte-parity cases passed`);
