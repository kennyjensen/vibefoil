// Measure actual worker message serialization/transfer, separately from solving.
import assert from 'node:assert/strict';
import { once } from 'node:events';
import { resolve } from 'node:path';
import { performance } from 'node:perf_hooks';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { isDeepStrictEqual, parseArgs } from 'node:util';
import { Worker, parentPort, workerData } from 'node:worker_threads';

function bufferStats(payload) {
  const buffers = new Set();
  function collect(value) {
    if (ArrayBuffer.isView(value)) buffers.add(value.buffer);
    else if (value && typeof value === 'object') Object.values(value).forEach(collect);
  }
  collect(payload);
  return { buffers: buffers.size,
    bytes: [...buffers].reduce((total, buffer) => total + buffer.byteLength, 0) };
}

if (workerData) {
  console.log = () => {};
  console.warn = () => {};
  globalThis.self = { postMessage(payload, transfers = []) {
    const stats = bufferStats(payload);
    const start = performance.now();
    parentPort.postMessage({ payload }, transfers);
    const sendMs = performance.now() - start;
    parentPort.postMessage({ sendMs, ...stats, transferred: transfers.length });
  } };
  await import(pathToFileURL(resolve(workerData.root, 'js/solver_worker.js')));
  parentPort.on('message', data => self.onmessage({ data }));
  parentPort.postMessage({ ready: true });
} else {
  const { values } = parseArgs({ options: {
    baseline: { type: 'string' },
    root: { type: 'string', default: fileURLToPath(new URL('../', import.meta.url)) },
    runs: { type: 'string', default: '11' },
    warmup: { type: 'string', default: '4' },
    direction: { type: 'string', default: 'result' },
  } });
  assert.ok(values.baseline, '--baseline must name a checkout to compare');
  const runs = Number(values.runs), warmup = Number(values.warmup);
  assert.ok(Number.isInteger(runs) && runs > 0, '--runs must be a positive integer');
  assert.ok(Number.isInteger(warmup) && warmup >= 0, '--warmup must be a nonnegative integer');
  assert.ok(['request', 'result'].includes(values.direction), '--direction must be request or result');
  const roots = [resolve(values.baseline), resolve(values.root)];
  const reuse = values.direction === 'request'
    ? (await import(pathToFileURL(resolve(roots[1], 'js/worker_reuse.js')))).createWorkerReuse()
    : null;
  const workers = roots.map(root => new Worker(new URL(import.meta.url), {
    workerData: { root }, argv: [],
  }));
  const base = {
    source: 'naca', mode: '4', m: 2, p: 4, t: 12,
    geometryKey: 'message-benchmark', advancedMode: false,
    reusePanel: true, reuseSolution: false, viscous: true,
    mach: 0.1, reynolds: 3e6, ncr: 9, nIter: 30,
    flap: { x: 0.75, yrel: 0.5, deflection: 0 },
    canvasWidth: 800, canvasHeight: 500,
  };
  function request(worker, message) {
    return new Promise((resolve, reject) => {
      let payload, requestTiming;
      const cleanup = () => {
        worker.off('message', onMessage);
        worker.off('error', onError);
      };
      const onError = error => { cleanup(); reject(error); };
      const onMessage = data => {
        if ('payload' in data) payload = data.payload;
        else { cleanup(); resolve({ payload, timing: reuse ? requestTiming : data }); }
      };
      worker.on('message', onMessage);
      worker.once('error', onError);
      const stats = reuse ? bufferStats(message) : null;
      const start = performance.now();
      worker.postMessage(message);
      requestTiming = { sendMs: performance.now() - start, ...stats, transferred: 0 };
    });
  }
  try {
    await Promise.all(workers.map(worker => once(worker, 'message')));
    let states = [null, null];
    const samples = [];
    for (let round = 0; round < runs + warmup; round += 1) {
      const alphaDeg = round % 5;
      const pair = [];
      for (const variant of round % 2 ? [1, 0] : [0, 1]) {
        let settings = {
          ...base, alphaDeg, alphaRad: alphaDeg * Math.PI / 180, reuseState: states[variant],
        };
        if (variant === 1 && reuse) settings = reuse.prepare(round, settings);
        pair[variant] = await request(workers[variant], { requestId: round, settings });
        if (variant === 1 && reuse) reuse.accept(pair[variant].payload);
        assert.ok(pair[variant].payload.ok, pair[variant].payload.error);
        assert.ok(pair[variant].payload.converged, `Round ${round}: unconverged`);
      }
      const after = { ...pair[1].payload };
      if (reuse) delete after.reuseCacheReady; // Only the opt-in transport ack differs.
      assert.ok(isDeepStrictEqual(pair[0].payload, after), `Round ${round}: outputs differ`);
      states = pair.map(result => result.payload.reuseState);
      if (round >= warmup) samples.push(pair.map(result => result.timing));
    }
    const median = values => [...values].sort((a, b) => a - b)[Math.floor(values.length / 2)];
    const summary = variant => ({
      medianSendMs: median(samples.map(sample => sample[variant].sendMs)),
      buffers: samples[0][variant].buffers, bytes: samples[0][variant].bytes,
      transferred: samples[0][variant].transferred,
    });
    process.stdout.write(`${JSON.stringify({
      node: process.version, arch: process.arch, roots, runs, warmupRounds: warmup,
      direction: values.direction,
      matchesReference: true, before: summary(0), after: summary(1), samples,
    }, null, 2)}\n`);
  } finally {
    await Promise.all(workers.map(worker => worker.terminate()));
  }
}
