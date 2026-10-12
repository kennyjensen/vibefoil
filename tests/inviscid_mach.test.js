// Regression for using the requested Mach in inviscid coefficient integration.
// Values came from original double-precision XFOIL, not the JS integrator.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

const fixture = JSON.parse(readFileSync(new URL('./canonical/inviscid_mach_reference.json', import.meta.url)));
const worker = new URL('./canonical/worker.mjs', import.meta.url);
const result = spawnSync(process.execPath, [fileURLToPath(worker), 'javascript'], {
  input: fixture.records.map(record => JSON.stringify(record.case)).join('\n') + '\n',
  cwd: fileURLToPath(new URL('../', import.meta.url)),
  encoding: 'utf8', timeout: 30000, maxBuffer: 8 * 1024 * 1024,
});
assert.equal(result.status, 0, result.error?.message ?? result.stderr);
const lines = result.stdout.trim().split('\n').map(line => JSON.parse(line));
assert.equal(lines.length, fixture.records.length);
for (let i = 0; i < lines.length; i += 1) {
  const expected = fixture.records[i];
  const actual = lines[i];
  assert.equal(actual.id, expected.case.id, actual.error);
  for (const [name, value] of Object.entries(expected.coefficients)) {
    assert.ok(Math.abs(actual.result.coefficients[name] - value) <= 1e-10,
      `${actual.id} ${name}: ${actual.result.coefficients[name]} vs ${value}`);
  }
}
console.log(`inviscid_mach: ${lines.length} website cases match original XFOIL coefficients`);
