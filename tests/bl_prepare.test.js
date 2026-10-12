import assert from 'node:assert/strict';
import { blpini, mrchue, mrchdu } from '../js/xbl.js';
import { NCOM, ensureCtx, blprv, blkin, blprvKin } from '../js/xblsys.js';

function context(mach = 0.1, Type = Float64Array) {
  const ctx = {
    GM1BL: 0.4, QINFBL: 1, TKBL: 0.01, TKBL_MS: 0.05,
    HSTINV: 0.4 * mach * mach / (1 + 0.2 * mach * mach),
    HSTINV_MS: 0.4 / (1 + 0.2 * mach * mach) ** 2,
    RSTBL: 1, RSTBL_MS: 0.02, REYBL: 3e6, REYBL_RE: 1,
    REYBL_MS: 0.03, HVRAT: 0.35,
    COM1: Type.from({ length: NCOM + 4 }, (_, i) => i + 100),
    COM2: Type.from({ length: NCOM + 4 }, (_, i) => -i - 100),
  };
  blpini(ctx);
  ensureCtx(ctx);
  // Fields outside the predictor/kinematic subset must also be published.
  ctx.HS2 = 1.6;
  ctx.CF2 = 0.004;
  ctx.DI2 = 0.003;
  return ctx;
}

function countStores(ctx, run) {
  const data = ctx.COM2;
  let writes = 0;
  ctx.COM2 = new Proxy(data, {
    get(target, key) { return Reflect.get(target, key, target); },
    set(target, key, value) {
      if (Number(key) >= 1 && Number(key) <= NCOM) writes += 1;
      return Reflect.set(target, key, value, target);
    },
  });
  try { run(); } finally { ctx.COM2 = data; }
  return writes;
}

let comparisons = 0;
function compare(input, args) {
  const expected = structuredClone(input);
  const actual = structuredClone(input);
  for (let repeat = 0; repeat < 2; repeat += 1) {
    const beforeStores = countStores(expected, () => {
      blprv(...args, expected);
      blkin(expected);
    });
    const afterStores = countStores(actual, () => blprvKin(...args, actual));
    assert.equal(beforeStores, 2 * NCOM);
    assert.equal(afterStores, NCOM);
    assert.deepStrictEqual(actual, expected, `entire station state ${comparisons}`);
    comparisons += 1;
  }
}

const baseArgs = [0.13, 2.3, 0.021, 0.0012, 0.0028, 0, 1.03];
for (const Type of [Array, Float64Array, Float32Array]) {
  for (const mach of [0, 0.1, 0.4]) {
    const input = context(mach, Type);
    compare(input, baseArgs);
    for (let index = 0; index < baseArgs.length; index += 1) {
      for (const value of [0, -0, Number.MIN_VALUE, 1e-155, 1e155, Infinity, NaN]) {
        const args = [...baseArgs];
        args[index] = value;
        compare(input, args);
      }
    }
  }
}

// The public predictor still publishes a complete snapshot on its own, and
// the combined routine supports the same explicit-context and method forms.
{
  const a = context(), b = context();
  assert.equal(countStores(a, () => blprv.call(a, ...baseArgs)), NCOM);
  blkin.call(a);
  assert.equal(countStores(b, () => blprvKin.call(b, ...baseArgs)), NCOM);
  assert.deepStrictEqual(a, b);
  comparisons += 1;
}

// Stop at the first assembled station so hook behavior can be checked without
// substituting any BL equations. Both marching entry points use the real
// predictor and kinematics, with their original callback receiver conventions.
function firstStation(march, customPredictor, customKinematics, replaceKinematics = false) {
  const ctx = context();
  const row = () => new Float64Array(3);
  for (const key of ['XSSI', 'UEDG', 'THET', 'DSTR', 'CTAU', 'TSTR']) {
    ctx[key] = Array.from({ length: 4 }, row);
  }
  ctx.XSSI[2].fill(0.1); ctx.XSSI[3].fill(0.2);
  ctx.UEDG[2].fill(0.8);
  ctx.THET[2].fill(0.001);
  ctx.DSTR[2].fill(0.0025);
  ctx.NBL = new Int32Array([0, 3, 3]);
  ctx.IBLTE = new Int32Array([0, 3, 3]);
  ctx.ITRAN = new Int32Array([0, 3, 3]);
  ctx.ACRIT = new Float64Array([0, 9, 9]);
  ctx.XSTRIP = new Float64Array([0, 1, 1]);
  ctx.blprv = blprv;
  ctx.blkin = blkin;
  const receiver = march === mrchue ? ctx : undefined;
  const calls = [];
  function customKin(value) {
    assert.equal(this, receiver);
    assert.equal(arguments.length, 1);
    assert.equal(value, ctx);
    // A custom kinematics callback must see the predictor's snapshot.
    assert.equal(ctx.COM2[1], ctx.X2);
    assert.equal(ctx.COM2[2], ctx.U2);
    assert.equal(ctx.COM2[3], ctx.T2);
    assert.equal(ctx.COM2[4], ctx.D2);
    calls.push('kin');
    blkin(ctx);
  }
  if (customKinematics) ctx.blkin = customKin;
  if (customPredictor) ctx.blprv = function (...args) {
    assert.equal(this, receiver);
    assert.equal(args.length, 8);
    assert.equal(args[7], ctx);
    calls.push('predictor');
    blprv(...args);
    if (replaceKinematics) ctx.blkin = customKin;
  };
  const done = new Error('Reached assembled station');
  ctx.blsys = function (value) {
    assert.equal(this, receiver);
    assert.equal(value, ctx);
    throw done;
  };
  const writes = countStores(ctx, () => {
    assert.throws(() => march(ctx), error => error === done);
  });
  const state = Object.fromEntries(Object.entries(ctx).filter(([, value]) => typeof value !== 'function'));
  return { state, calls, writes };
}

for (const march of [mrchue, mrchdu]) {
  const standard = firstStation(march, false, false);
  assert.equal(standard.writes, NCOM);
  for (const [predictor, kinematics] of [[true, false], [false, true], [true, true]]) {
    const custom = firstStation(march, predictor, kinematics);
    assert.deepStrictEqual(custom.state, standard.state);
    assert.deepStrictEqual(custom.calls, [predictor && 'predictor', kinematics && 'kin'].filter(Boolean));
    assert.equal(custom.writes, 2 * NCOM);
    comparisons += 1;
  }
  const replaced = firstStation(march, true, false, true);
  assert.deepStrictEqual(replaced.state, standard.state);
  // MRCHUE reads the next method after the predictor; MRCHDU caches both.
  assert.deepStrictEqual(replaced.calls, march === mrchue ? ['predictor', 'kin'] : ['predictor']);
  assert.equal(replaced.writes, 2 * NCOM);
  comparisons += 1;
}

console.log(`bl_prepare: ${comparisons} exact station-state, snapshot-count, and callback comparisons passed`);
