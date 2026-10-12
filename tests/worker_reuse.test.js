import assert from 'node:assert/strict';
import { once } from 'node:events';
import { isDeepStrictEqual } from 'node:util';
import { Worker } from 'node:worker_threads';
import { createWorkerReuse } from '../js/worker_reuse.js';

const saved = { geometryKey: 'foil', bl: { LVCONV: true } };
const settings = { geometryKey: 'foil', viscous: true, reusePanel: true, reuseState: saved };
const ack = { requestId: 1, ok: true, converged: true, reuseState: saved, reuseCacheReady: true };
function acknowledged() {
  const tracker = createWorkerReuse();
  tracker.prepare(1, settings);
  tracker.accept(ack);
  return tracker;
}
assert.equal(acknowledged().prepare(2, settings).reuseState, null);
for (const change of [
  { reuseState: structuredClone(saved) }, { geometryKey: 'other' },
  { reusePanel: false },
]) {
  const input = { ...settings, ...change };
  assert.equal(acknowledged().prepare(2, input).reuseState, input.reuseState);
}
for (const change of [{ viscous: false }, { advancedMode: true, reuseSolution: false }]) {
  // These modes ignore the snapshot even without a prior cache acknowledgment.
  const tracker = createWorkerReuse();
  assert.equal(tracker.prepare(1, { ...settings, ...change }).reuseState, null);
  tracker.accept({ requestId: 1, ok: true, converged: true, reuseState: null });
  assert.equal(tracker.prepare(2, settings).reuseState, saved);
}
for (const change of [
  { ok: false }, { converged: false }, { reuseCacheReady: false },
  { reuseCacheReady: undefined }, { requestId: 0 }, { type: 'dump' },
  { reuseState: { ...saved, bl: { LVCONV: false } } },
]) {
  const tracker = createWorkerReuse();
  tracker.prepare(1, settings);
  tracker.accept({ ...ack, ...change });
  assert.equal(tracker.prepare(2, settings).reuseState, saved);
}
{
  const tracker = acknowledged();
  tracker.clear(); // Restart, error, or inverse design invalidates live state.
  assert.equal(tracker.prepare(2, settings).reuseState, saved);
  tracker.clear(); // QDES queued before the preceding solve's response arrives.
  tracker.accept({ ...ack, requestId: 2 });
  assert.equal(tracker.prepare(3, settings).reuseState, saved);
  tracker.accept({ ...ack, requestId: 2 }); // A stale response cannot restore it.
  tracker.accept({ ...ack, requestId: 3, type: 'qdes' });
  tracker.accept({ ...ack, requestId: 3 });
  assert.equal(tracker.prepare(4, settings).reuseState, null);
}

// Real workers compare sending every saved snapshot with the application's
// omission/acknowledgment protocol. Both honor output transfer/detachment.
const workerSource = `
  import { parentPort, workerData } from 'node:worker_threads';
  console.log = () => {};
  console.warn = () => {};
  globalThis.self = { postMessage(payload, transfers = []) {
    parentPort.postMessage(payload, transfers);
  } };
  await import(workerData.module);
  parentPort.on('message', data => self.onmessage({ data }));
  parentPort.postMessage({ ready: true });
`;
const module = new URL('../js/solver_worker.js', import.meta.url).href;
const tracker = createWorkerReuse();
let workers = [], states = [null, null], requestId = 0, comparisons = 0, omissions = 0;
async function restart() {
  tracker.clear();
  await Promise.all(workers.map(worker => worker.terminate()));
  workers = [0, 1].map(() => new Worker(workerSource, { eval: true, type: 'module', workerData: { module } }));
  await Promise.all(workers.map(worker => once(worker, 'message')));
}

async function pair(messages, expected = {}) {
  const id = ++requestId;
  const results = [];
  for (const i of [0, 1]) {
    let message = { requestId: id, ...messages[i] };
    if (i === 1 && message.settings) {
      const prepared = tracker.prepare(id, message.settings);
      const omitted = !!message.settings.reuseState && prepared.reuseState === null;
      if (expected.omitted !== undefined) assert.equal(omitted, expected.omitted, `Request ${id}: omission`);
      if (omitted) omissions += 1;
      assert.equal(message.settings.reuseState, states[i]); // Inputs remain intact.
      message = { ...message, settings: prepared };
    } else if (i === 1 && message.action === 'qdes') {
      tracker.clear();
    }
    const response = once(workers[i], 'message');
    workers[i].postMessage(message);
    [results[i]] = await response;
    if (i === 1) tracker.accept(results[i]);
  }
  const { reuseCacheReady, ...after } = results[1];
  assert.ok(isDeepStrictEqual(results[0], after), `Request ${id}: complete results differ`);
  assert.equal(after.ok, expected.ok ?? true, after.error);
  if (expected.ready !== undefined) assert.equal(reuseCacheReady, expected.ready, `Request ${id}: cache readiness`);
  if (expected.converged !== undefined) assert.equal(after.converged, expected.converged, `Request ${id}: convergence`);
  comparisons += 1;
  return results;
}

const base = {
  source: 'naca', mode: '4', m: 2, p: 4, t: 12,
  geometryKey: 'reuse-naca2412', advancedMode: false,
  reusePanel: true, reuseSolution: false, viscous: true,
  mach: 0.1, reynolds: 3e6, ncr: 9, nIter: 30,
  alphaDeg: 2, alphaRad: 2 * Math.PI / 180,
  flap: { x: 0.75, yrel: 0.5, deflection: 0 },
  canvasWidth: 800, canvasHeight: 500,
};
async function solve(changes = {}, expected = {}) {
  const results = await pair(states.map(reuseState => ({ settings: { ...base, ...changes, reuseState } })), expected);
  // Match the page: a result without a snapshot leaves the saved backup intact.
  results.forEach((result, i) => { if (result.reuseState) states[i] = result.reuseState; });
  return results;
}
try {
  await restart();
  await solve({}, { omitted: false, ready: true, converged: true });
  let solved = await solve({ alphaDeg: 4, alphaRad: 4 * Math.PI / 180, ncr: 7 },
    { omitted: true, ready: true, converged: true });
  for (const kind of ['cp', 'bl']) await pair([0, 1].map(() => ({ action: 'dump', kind })));
  await solve({ ncr: 7, reynolds: 2e6, mach: 0.15 }, { omitted: true, ready: true, converged: true });

  await restart();
  await solve({ ncr: 7, advancedMode: true, reuseSolution: true }, { omitted: false, ready: true, converged: true });
  await solve({}, { omitted: true, ready: true, converged: true });

  // A fresh, unconverged advanced solve returns no snapshot. The page retains
  // its older converged backup, which the next request must send in full.
  await solve({ advancedMode: true, reuseSolution: false, nIter: 1,
    alphaDeg: 18, alphaRad: 18 * Math.PI / 180 }, { omitted: true, ready: false, converged: false });
  await solve({}, { omitted: false, ready: true, converged: true });
  await solve({ viscous: false }, { omitted: true, ready: false });
  await solve({}, { omitted: false, ready: true, converged: true });
  await solve({ source: 'custom', custom: null }, { omitted: true, ok: false });
  solved = await solve({}, { omitted: false, ready: true, converged: true });

  await pair(solved.map(result => ({ action: 'qdes', niter: 1,
    cpSpec: result.cpData.invAll.map(point => point.cp) })));
  await pair([0, 1].map(() => ({ action: 'dump', kind: 'cp' })));
  await solve({}, { omitted: false, ready: true });
  await solve({}, { omitted: true, ready: true });

  // A snapshot can be returned even if an incompatible restore left the live
  // BL cache empty. Only the worker can confirm this distinction reliably.
  states = states.map(state => {
    const copy = structuredClone(state);
    copy.bl.N += 1;
    return copy;
  });
  await restart();
  await solve({ advancedMode: true, reuseSolution: true }, { omitted: false, ready: false, converged: true });
  await solve({ advancedMode: true, reuseSolution: true }, { omitted: false, ready: true, converged: true });
  await solve({ advancedMode: true, reuseSolution: true }, { omitted: true, ready: true, converged: true });
  console.log(`worker_reuse.test.js: OK (${comparisons} exact worker comparisons, ${omissions} omitted snapshots)`);
} finally {
  await Promise.all(workers.map(worker => worker.terminate()));
}
