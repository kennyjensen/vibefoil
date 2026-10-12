// Alternate isolated revisions to reduce drift from CPU frequency/system load.
import assert from 'node:assert/strict';
import { once } from 'node:events';
import { resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { isDeepStrictEqual, parseArgs } from 'node:util';
import { Worker } from 'node:worker_threads';

const { values } = parseArgs({ options: {
  baseline: { type: 'string' },
  root: { type: 'string', default: fileURLToPath(new URL('../', import.meta.url)) },
  runs: { type: 'string', default: '9' },
  warmup: { type: 'string', default: '4' },
} });
assert.ok(values.baseline, '--baseline must name a checkout to compare');
const runs = Number(values.runs);
const warmup = Number(values.warmup);
assert.ok(Number.isInteger(runs) && runs > 0, '--runs must be a positive integer');
assert.ok(Number.isInteger(warmup) && warmup >= 0, '--warmup must be a nonnegative integer');
const roots = [resolve(values.baseline), resolve(values.root)];
const workers = roots.map(root => new Worker(new URL('./solver.mjs', import.meta.url), {
  workerData: { root, pairedBenchmark: true }, argv: [],
}));

try {
  const ready = await Promise.all(workers.map(worker => once(worker, 'message')));
  const names = ready[0][0].scenarios;
  assert.deepEqual(names, ready[1][0].scenarios);
  const samples = names.map(() => []);
  const cpuSamples = names.map(() => []);
  for (let round = 0; round < runs + warmup; round += 1) {
    for (let index = 0; index < names.length; index += 1) {
      const pair = [];
      for (const variant of round % 2 ? [1, 0] : [0, 1]) {
        const response = once(workers[variant], 'message');
        workers[variant].postMessage(index);
        [pair[variant]] = await response;
      }
      assert.ok(isDeepStrictEqual(pair[0].outputs, pair[1].outputs),
        `${names[index]}, round ${round}: complete worker outputs differ`);
      if (round >= warmup) {
        samples[index].push(pair.map(result => result.elapsedMs));
        if (pair.every(result => Number.isFinite(result.threadCpuMs))) {
          cpuSamples[index].push(pair.map(result => result.threadCpuMs));
        }
      }
    }
  }
  const median = array => [...array].sort((a, b) => a - b)[Math.floor(array.length / 2)];
  process.stdout.write(`${JSON.stringify({
    node: process.version, arch: process.arch, roots, runs, warmupRounds: warmup,
    matchesReference: true,
    timings: Object.fromEntries(names.map((name, i) => [name, {
      beforeMs: median(samples[i].map(sample => sample[0])),
      afterMs: median(samples[i].map(sample => sample[1])),
      medianReductionPct: median(samples[i].map(sample => 100 * (1 - sample[1] / sample[0]))),
      samples: samples[i],
      threadCpu: cpuSamples[i].length === runs ? {
        beforeMs: median(cpuSamples[i].map(sample => sample[0])),
        afterMs: median(cpuSamples[i].map(sample => sample[1])),
        medianReductionPct: median(cpuSamples[i].map(sample => 100 * (1 - sample[1] / sample[0]))),
        samples: cpuSamples[i],
      } : undefined,
    }])),
  }, null, 2)}\n`);
} finally {
  await Promise.all(workers.map(worker => worker.terminate()));
}
