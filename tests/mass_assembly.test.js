import assert from 'node:assert/strict';
import { createTensor3Flat } from '../js/arrays.js';
import { assembleMassRows } from '../js/xbl.js';

// The previous three-pass assembly is the exact reference. Compare the entire
// context, including padding, unrelated columns, and read-only source arrays.
function reference(ctx, iv, d1M, u1M, d2M, u2M, ule1M, ule2M, xiUle1, xiUle2) {
  for (let k = 1; k <= 3; k += 1) {
    const offset = 6 * k;
    const vs1 = ctx.VS1F;
    const vs2 = ctx.VS2F;
    const v1d = vs1[offset + 3];
    const v1u = vs1[offset + 4];
    const v2d = vs2[offset + 3];
    const v2u = vs2[offset + 4];
    const vx = vs1[offset + 5] + vs2[offset + 5] + ctx.VSX[k];
    for (let jv = 1; jv <= ctx.NSYS; jv += 1) {
      ctx.VM[k][jv][iv] = v1d * d1M[jv] + v1u * u1M[jv]
        + v2d * d2M[jv] + v2u * u2M[jv]
        + vx * (xiUle1 * ule1M[jv] + xiUle2 * ule2M[jv]);
    }
  }
}

let seed = 0x24122301;
function random() {
  seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
  return seed / 0x100000000 - 0.5;
}
function makeCase(nsys) {
  const vm = createTensor3Flat(3, nsys, nsys);
  for (let i = 0; i < vm.flat.length; i += 1) vm.flat[i] = random();
  return {
    ctx: { NSYS: nsys, VM: vm.view, VMF: vm.flat,
      VS1F: Float64Array.from({ length: 30 }, random),
      VS2F: Float64Array.from({ length: 30 }, random),
      VSX: Float64Array.from({ length: 5 }, random) },
    sources: Array.from({ length: 6 }, () => Float64Array.from({ length: nsys + 1 }, random)),
    xi: [random(), random()],
  };
}
let comparisons = 0;
function compare(data, iv) {
  const expected = structuredClone(data);
  reference(expected.ctx, iv, ...expected.sources, ...expected.xi);
  const nested = structuredClone(data);
  delete nested.ctx.VMF;
  const expectedNested = structuredClone(expected);
  delete expectedNested.ctx.VMF;
  assembleMassRows(nested.ctx, iv, ...nested.sources, ...nested.xi);
  assert.deepStrictEqual(nested, expectedNested, `nested NSYS=${data.ctx.NSYS}, iv=${iv}, case=${comparisons}`);
  assembleMassRows(data.ctx, iv, ...data.sources, ...data.xi);
  assert.deepStrictEqual(data, expected, `NSYS=${data.ctx.NSYS}, iv=${iv}, case=${comparisons}`);
  comparisons += 2;
}

for (const nsys of [1, 2, 8, 31, 240]) {
  for (let sample = 0; sample < 4; sample += 1) {
    const data = makeCase(nsys);
    for (const iv of new Set([1, Math.ceil(nsys / 2), nsys])) compare(data, iv);
    compare(data, nsys); // Repeated assembly overwrites prior results exactly.
  }
}
const targets = [
  data => [data.ctx.VS1F, 9], data => [data.ctx.VS1F, 17],
  data => [data.ctx.VS2F, 22], data => [data.ctx.VSX, 1],
  data => [data.sources[0], 1], data => [data.sources[5], 2],
  data => [data.xi, 0], data => [data.xi, 1],
];
for (const value of [0, -0, Number.MIN_VALUE, Number.MAX_VALUE, Infinity, -Infinity, NaN]) {
  for (const target of targets) {
    const data = makeCase(3);
    const [array, index] = target(data);
    array[index] = value;
    compare(data, 2);
  }
}
// Existing storage may be larger than the active system. Preserve its stride
// and all inactive entries rather than assuming a tightly sized flat tensor.
const oversized = makeCase(8);
oversized.ctx.NSYS = 3;
for (const iv of [1, 2, 3]) compare(oversized, iv);
console.log(`mass_assembly: ${comparisons} exact fused/reference comparisons passed`);
