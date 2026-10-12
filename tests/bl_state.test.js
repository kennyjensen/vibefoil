import assert from 'node:assert/strict';
import { blpini } from '../js/xbl.js';
import {
  NCOM, ensureCtx, syncVarsToCom, copyCom,
  blprv, blkin, blvar, blmid, blsys, trchek, trdif,
} from '../js/xblsys.js';

function state(mach, type, similarity = false) {
  const ctx = {
    GM1BL: 0.4, QINFBL: 1, TKBL: 0.01, TKBL_MS: 0.05,
    HSTINV: 0.4 * mach * mach / (1 + 0.2 * mach * mach),
    HSTINV_MS: 0.4 / (1 + 0.2 * mach * mach) ** 2,
    RSTBL: 1, RSTBL_MS: 0.02, REYBL: 3e6, REYBL_RE: 1,
    REYBL_MS: 0.03, HVRAT: 0.35, BULE: 1, IDAMPV: 0, AMCRIT: 9, XIFORC: 1,
    SIMI: similarity, TRAN: false, TURB: type !== 1, WAKE: type === 3,
    COM1: new Float64Array(NCOM + 3), COM2: new Float64Array(NCOM + 3),
  };
  // Nonzero padding detects any expansion of the active snapshot range.
  ctx.COM1[0] = ctx.COM2[0] = -123;
  ctx.COM1[NCOM + 1] = ctx.COM2[NCOM + 1] = 456;
  ctx.COM1[NCOM + 2] = ctx.COM2[NCOM + 2] = -789;
  blpini(ctx);
  ensureCtx(ctx);
  blprv(0.1, 2, 0.02, 0.001, 0.0024, 0, 1.04, ctx);
  blkin(ctx);
  blvar(type, ctx);
  copyCom(ctx, 2, 1);
  blprv(0.13, 2.3, 0.021, 0.0012, 0.0028,
    type === 3 ? 0.0001 : 0, 1.03, ctx);
  blkin(ctx);
  return ctx;
}

// Count actual snapshot stores, without changing their Float64 conversion.
function stores(ctx, run) {
  const data = ctx.COM2;
  let writes = 0;
  ctx.COM2 = new Proxy(data, {
    get(target, key) { return Reflect.get(target, key, target); },
    set(target, key, value) {
      if (Number(key) >= 1 && Number(key) <= NCOM) writes += 1;
      return Reflect.set(target, key, value, target);
    },
  });
  try {
    run();
  } finally {
    ctx.COM2 = data;
  }
  return writes;
}

let comparisons = 0;
for (const mach of [0, 0.1, 0.4]) {
  for (const type of [1, 2, 3]) {
    for (const similarity of [false, true]) {
      const ctx = state(mach, type, similarity);

      // BLKIN consumes station variables, even when COM2 contains stale data.
      // Pre-copying those variables must not affect any resulting state.
      ctx.U2 *= 1.03;
      ctx.T2 *= 1.02;
      ctx.D2 *= 0.99;
      ctx.COM2.fill(2026, 1, NCOM + 1);
      const expected = structuredClone(ctx);
      syncVarsToCom(expected, 2);
      blkin(expected);
      assert.equal(stores(ctx, () => blkin(ctx)), NCOM);
      assert.deepStrictEqual(ctx, expected);
      comparisons += 1;

      // BLVAR publishes the complete station-2 snapshot. BLMID must leave
      // that snapshot valid, including when it updates similarity station 1.
      blvar(type, ctx);
      const snapshot = ctx.COM2.slice();
      blmid(type, ctx);
      syncVarsToCom(ctx, 2);
      assert.deepStrictEqual(ctx.COM2, snapshot);
      assert.ok(snapshot.slice(1, NCOM + 1).every(Number.isFinite));
      comparisons += 1;

      // Normal BLSYS assembly needs only BLVAR's snapshot, and must expose it
      // unchanged after assembling residuals and derivatives.
      assert.equal(stores(ctx, () => blsys(ctx)), NCOM);
      assert.deepStrictEqual(ctx.COM2, snapshot);
      assert.ok(ctx.VSREZ.slice(1, 4).every(Number.isFinite));
      comparisons += 1;
    }
  }

  for (const damping of [0, 1]) {
    for (const force of [0.115, 0.13 - 1e-13, 0.13, 0.13 + 1e-13, 1]) {
      const ctx = state(mach, 1);
      ctx.IDAMPV = damping;
      ctx.XIFORC = force;
      const initial = ctx.COM2.slice();
      const writes = stores(ctx, () => trchek(ctx));
      // Each transition iteration saves the BLKIN result, then restores COM2.
      assert.ok(writes >= 2 * NCOM);
      assert.equal(writes % (2 * NCOM), 0);
      assert.deepStrictEqual(ctx.COM2, initial);
      comparisons += 1;

      if (!ctx.TRAN) continue;
      const expected = structuredClone(ctx);
      blvar(2, expected);
      blmid(2, expected);
      const station1 = expected.COM1.slice();
      const station2 = expected.COM2.slice();
      assert.equal(stores(expected, () => trdif(expected)), 5 * NCOM);
      assert.deepStrictEqual(expected.COM1, station1);
      assert.deepStrictEqual(expected.COM2, station2);
      assert.ok(expected.VSREZ.slice(1, 4).every(Number.isFinite));
      comparisons += 1;

      assert.equal(stores(ctx, () => blsys(ctx)), 6 * NCOM);
      assert.deepStrictEqual(ctx.COM2, station2);
      assert.ok(ctx.VSREZ.slice(1, 4).every(Number.isFinite));
      comparisons += 1;
    }
  }
}

console.log(`bl_state: ${comparisons} snapshot, transition-restoration, and work-count checks passed`);
