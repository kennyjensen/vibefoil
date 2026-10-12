// Time the actual browser worker's synchronous compute path in Node. Importing
// a separate checkout with --root enables before/after comparisons on one host.
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { performance } from 'node:perf_hooks';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { isDeepStrictEqual, parseArgs } from 'node:util';
import { deserialize, serialize } from 'node:v8';
import { parentPort, workerData } from 'node:worker_threads';

const { values } = parseArgs({ options: {
  root: { type: 'string', default: fileURLToPath(new URL('../', import.meta.url)) },
  runs: { type: 'string', default: '5' },
  save: { type: 'string' },
  compare: { type: 'string' },
} });
const runs = Number(values.runs);
assert.ok(Number.isInteger(runs) && runs > 0, '--runs must be a positive integer');
const root = resolve(workerData?.root ?? values.root);
let response;
globalThis.self = { postMessage(message) { response = message; } };
// Logging is deliberately excluded from the numerical timing.
console.log = () => {};
console.warn = () => {};
await import(pathToFileURL(resolve(root, 'js/solver_worker.js')));

const base = {
  mode: '4', source: 'naca', m: 2, p: 4, t: 12,
  flap: { x: 0.75, yrel: 0.5, deflection: 0 },
  geometryKey: 'benchmark-naca2412', advancedMode: true,
  reusePanel: false, reuseSolution: false, viscous: true,
  mach: 0.1, reynolds: 3e6, ncr: 9, nIter: 30,
  canvasWidth: 800, canvasHeight: 500,
};
const scenarios = [
  { name: 'inviscid', settings: { viscous: false }, angles: [4] },
  { name: 'viscous', settings: {}, angles: [4] },
  { name: 'flap', settings: { geometryKey: 'benchmark-flap',
    flap: { x: 0.75, yrel: 0.5, deflection: 10 } }, angles: [4] },
  { name: 'polar-sweep', settings: {}, angles: [2, 4, 6] },
];

function runScenario(scenario) {
  const results = [];
  let reuseState = null;
  for (const [index, alphaDeg] of scenario.angles.entries()) {
    const settings = { ...base, ...scenario.settings,
      alphaDeg, alphaRad: alphaDeg * Math.PI / 180,
      reusePanel: index > 0, reuseSolution: index > 0, reuseState };
    response = null;
    self.onmessage({ data: { requestId: index, settings } });
    assert.ok(response?.ok, response?.error ?? 'No worker response');
    assert.ok(response.converged, `${scenario.name}: unconverged at ${alphaDeg} degrees`);
    results.push(response);
    reuseState = response.reuseState;
  }
  return results;
}

// The paired runner keeps each revision in a separate worker and alternates
// measurements, avoiding shared module state and long gaps between samples.
function readThreadCpu(previous) {
  try {
    return process.threadCpuUsage?.(previous);
  } catch {
    // Optional counters must not prevent comparisons on unsupported runtimes.
    return undefined;
  }
}

if (workerData?.pairedBenchmark) {
  parentPort.on('message', index => {
    assert.ok(Number.isInteger(index) && scenarios[index], 'Unknown scenario');
    const cpuStart = readThreadCpu();
    const start = performance.now();
    const outputs = runScenario(scenarios[index]);
    const elapsedMs = performance.now() - start;
    const cpu = cpuStart ? readThreadCpu(cpuStart) : null;
    const threadCpuMs = cpu ? (cpu.user + cpu.system) / 1000 : undefined;
    parentPort.postMessage({ elapsedMs, threadCpuMs, outputs });
  });
  parentPort.postMessage({ scenarios: scenarios.map(({ name }) => name) });
} else {
  const timings = Object.fromEntries(scenarios.map(({ name }) => [name, { wallMs: [], cpuMs: [] }]));
  let outputs;
  for (let round = 0; round < runs + 2; round += 1) {
    outputs = [];
    for (const scenario of scenarios) {
      const cpuStart = process.cpuUsage();
      const start = performance.now();
      const results = runScenario(scenario);
      const wallMs = performance.now() - start;
      const cpu = process.cpuUsage(cpuStart);
      if (round >= 2) {
        timings[scenario.name].wallMs.push(wallMs);
        timings[scenario.name].cpuMs.push((cpu.user + cpu.system) / 1000);
      }
      outputs.push(results);
    }
  }
  // Preserve all numbers and typed arrays for exact same-runtime comparisons.
  if (values.compare) {
    const reference = deserialize(readFileSync(values.compare));
    assert.ok(isDeepStrictEqual(outputs, reference), 'Complete worker results differ from the reference');
  }
  if (values.save) writeFileSync(values.save, serialize(outputs));
  const median = (samples) => [...samples].sort((a, b) => a - b)[Math.floor(samples.length / 2)];
  process.stdout.write(`${JSON.stringify({
    node: process.version, arch: process.arch, root, runs, warmupRounds: 2,
    matchesReference: values.compare ? true : undefined,
    timings: Object.fromEntries(Object.entries(timings).map(([name, samples]) => [name, {
      medianMs: median(samples.wallMs), medianCpuMs: median(samples.cpuMs), ...samples,
    }])),
  }, null, 2)}\n`);
}
