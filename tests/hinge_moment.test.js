import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { computeHingeMoment } from '../js/hinge_moment.js';
import { makeHingeFixture } from './hinge_moment_fixture.js';

const reference = JSON.parse(readFileSync(new URL('./hinge_moment_reference.json', import.meta.url)));
assert.equal(createHash('sha256').update(readFileSync(new URL('../third_party/Xfoil/src/xoper.f', import.meta.url))).digest('hex'),
  reference.sourceSha256, 'Regenerate fixtures if original MHINGE changes');
for (const input of reference.cases) {
  const fixture = makeHingeFixture(input);
  const before = structuredClone(fixture);
  const value = computeHingeMoment(fixture.panel, fixture.cp, fixture.hinge);
  assert.ok(Number.isFinite(value) && Math.abs(value - input.ch) < 1e-13,
    `NACA ${input.foil}, flap ${input.deflection}: ${value} vs ${input.ch}`);
  assert.deepStrictEqual(fixture, before, 'Hinge postprocessing must not change solver state');
}
const {panel, cp, hinge} = makeHingeFixture(reference.cases[2]);
for (const invalid of [null, {x: NaN, y: 0}, {x: -1, y: 0}, {x: 2, y: 0}]) {
  assert.ok(Number.isNaN(computeHingeMoment(panel, cp, invalid)));
}
assert.ok(Number.isNaN(computeHingeMoment(panel, new Float64Array(2), hinge)));
cp[1] = NaN;
assert.ok(Number.isNaN(computeHingeMoment(panel, cp, hinge)));
console.log(`hinge_moment: ${reference.cases.length} original-XFOIL comparisons and read-only/invalid-input checks passed`);
