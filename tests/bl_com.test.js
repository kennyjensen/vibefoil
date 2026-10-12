import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import {
  NCOM, ensureCtx, syncComToVars, syncVarsToCom, copyCom,
} from '../js/xblsys.js';

// Use the original COMMON-block declarations as the independent field order.
// This checks every offset without duplicating the optimized JS assignments.
const source = readFileSync(new URL('../third_party/Xfoil/src/XBL.INC', import.meta.url), 'utf8');
const fields = [null, ...[1, 2].map(station => {
  const block = source.match(new RegExp(
    `COMMON/V_VAR${station}/([\\s\\S]*?)(?=\\n {6}(?:COMMON|EQUIVALENCE))`,
  ));
  assert.ok(block, `COMMON block ${station}`);
  const names = block[1].replace(/\n\s*&/g, ' ').split(',').map(name => name.trim());
  assert.equal(names.length, NCOM);
  return [null, ...names];
})];

function referenceSync(ctx, which, toArray) {
  ensureCtx(ctx);
  const station = which === 1 ? 1 : 2;
  const data = ctx[`COM${station}`];
  for (let i = 1; i <= NCOM; i += 1) {
    if (toArray) data[i] = ctx[fields[station][i]];
    else ctx[fields[station][i]] = data[i];
  }
}

function referenceCopy(ctx, from, to) {
  ensureCtx(ctx);
  const src = ctx[from === 1 ? 'COM1' : 'COM2'];
  const dst = ctx[to === 1 ? 'COM1' : 'COM2'];
  for (let i = 1; i <= NCOM; i += 1) dst[i] = src[i];
  referenceSync(ctx, to, false);
}

function state(Type, length, initialized, special) {
  const ctx = { unrelated: 1234 };
  for (const station of [1, 2]) {
    ctx[`COM${station}`] = Type.from({ length }, (_, i) => special && i % 4 === 0
      ? special[Math.floor(i / 4) % special.length] : 1000 * station + i);
  }
  if (initialized) ensureCtx(ctx);
  for (const station of [1, 2]) {
    for (let i = 1; i <= NCOM; i += 1) {
      ctx[fields[station][i]] = special && i % 3 === 0
        ? special[Math.floor(i / 3) % special.length] : -1000 * station - i;
    }
  }
  return ctx;
}

let comparisons = 0;
function compare(input, actual, expected) {
  const a = structuredClone(input);
  const b = structuredClone(input);
  actual(a);
  expected(b);
  assert.deepStrictEqual(a, b, `entire state, case ${comparisons}`);
  comparisons += 1;
}

const special = [0, -0, Number.MIN_VALUE, 1e-155, 1e155, Infinity, -Infinity, NaN, undefined, null];
for (const Type of [Array, Float64Array, Float32Array]) {
  for (const length of [0, 1, 10, NCOM + 1, NCOM + 4]) {
    for (const initialized of [false, true]) {
      for (const values of [null, special]) {
        const input = state(Type, length, initialized, values);
        for (const which of [1, 2, 0, '1', undefined]) {
          compare(input, ctx => syncComToVars(ctx, which), ctx => referenceSync(ctx, which, false));
          compare(input, ctx => syncVarsToCom(ctx, which), ctx => referenceSync(ctx, which, true));
        }
        for (const from of [1, 2]) {
          for (const to of [1, 2]) {
            compare(input, ctx => copyCom(ctx, from, to), ctx => referenceCopy(ctx, from, to));
          }
        }
      }
    }
  }
}

// Preserve copying order when the two snapshots share or overlap storage.
for (const offset of [0, 1, 5, NCOM]) {
  const input = state(Float64Array, NCOM + 4, true, special);
  const shared = Float64Array.from({ length: 2 * NCOM + 4 }, (_, i) => i - 20.5);
  input.COM1 = shared.subarray(0, NCOM + 4);
  input.COM2 = shared.subarray(offset, offset + NCOM + 4);
  for (const from of [1, 2]) {
    for (const to of [1, 2]) {
      compare(input, ctx => copyCom(ctx, from, to), ctx => referenceCopy(ctx, from, to));
    }
  }
}

// Numeric arrays retain their conversion behavior; field reads and writes must
// also occur in the original order, including with observable property access.
function trace(run) {
  const events = [];
  const input = state(Float64Array, NCOM + 4, true, special);
  for (const key of ['COM1', 'COM2']) {
    input[key] = new Proxy(input[key], {
      get(target, index) {
        events.push(['get', key, index]);
        return Reflect.get(target, index, target);
      },
      set(target, index, value) {
        events.push(['set', key, index, value]);
        return Reflect.set(target, index, value, target);
      },
    });
  }
  run(new Proxy(input, {
    get(target, key) {
      events.push(['get', key]);
      return Reflect.get(target, key, target);
    },
    set(target, key, value) {
      events.push(['set', key, value]);
      return Reflect.set(target, key, value, target);
    },
  }));
  return events;
}
for (const which of [1, 2]) {
  assert.deepStrictEqual(trace(ctx => syncComToVars(ctx, which)),
    trace(ctx => referenceSync(ctx, which, false)));
  assert.deepStrictEqual(trace(ctx => syncVarsToCom(ctx, which)),
    trace(ctx => referenceSync(ctx, which, true)));
  comparisons += 2;
}

console.log(`bl_com: ${comparisons} exact COMMON-layout, conversion, aliasing, and access-order comparisons passed`);
