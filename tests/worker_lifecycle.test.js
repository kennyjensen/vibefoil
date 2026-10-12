import assert from 'node:assert/strict';
import { once } from 'node:events';
import { isDeepStrictEqual } from 'node:util';
import { Worker } from 'node:worker_threads';

// Run the actual solver in isolated workers. One copies messages; the other
// honors transfer lists and detaches sender buffers, as the browser does.
const workerSource = `
  import { parentPort, workerData } from 'node:worker_threads';
  console.log = () => {};
  console.warn = () => {};
  globalThis.self = { postMessage(message, transfers = []) {
    parentPort.postMessage(message, workerData.transfer ? transfers : []);
  } };
  await import(workerData.module);
  parentPort.on('message', data => self.onmessage({ data }));
  parentPort.postMessage({ ready: true });
`;
const module = new URL('../js/solver_worker.js', import.meta.url).href;
let workers = [];
let requestId = 0;
let comparisons = 0;
async function restart() {
  await Promise.all(workers.map(worker => worker.terminate()));
  workers = [false, true].map(transfer => new Worker(workerSource, {
    eval: true, type: 'module', workerData: { module, transfer },
  }));
  await Promise.all(workers.map(worker => once(worker, 'message')));
}

async function sendPair(messages) {
  const id = requestId++;
  const results = await Promise.all(workers.map(async (worker, index) => {
    const response = once(worker, 'message');
    worker.postMessage({ requestId: id, ...messages[index] });
    const [result] = await response;
    assert.ok(result.ok, result.error);
    return result;
  }));
  assert.ok(isDeepStrictEqual(results[0], results[1]), `Message ${id}: copy/transfer outputs differ`);
  comparisons += 1;
  return results;
}

const base = {
  source: 'naca', mode: '4', m: 2, p: 4, t: 12,
  geometryKey: 'ownership-naca2412', advancedMode: false,
  reusePanel: true, reuseSolution: false, viscous: true,
  mach: 0.1, reynolds: 3e6, ncr: 9, nIter: 30,
  alphaDeg: 2, alphaRad: 2 * Math.PI / 180,
  flap: { x: 0.75, yrel: 0.5, deflection: 0 },
  canvasWidth: 800, canvasHeight: 500,
};
try {
  await restart();
  let solved = await sendPair([0, 1].map(() => ({ settings: base })));
  assert.ok(solved.every(result => result.converged && result.reuseState));
  // ACRIT in the displayed coefficients references live solver state.
  solved = await sendPair(solved.map(result => ({ settings: {
    ...base, alphaDeg: 4, alphaRad: 4 * Math.PI / 180, ncr: 7,
    reuseState: result.reuseState,
  } })));
  assert.ok(solved.every(result => result.converged && result.coeffs.ACRIT[1] === 7));
  for (const kind of ['cp', 'bl']) {
    await sendPair([0, 1].map(() => ({ action: 'dump', kind })));
  }

  // Restoring a saved packed snapshot must work after losing worker caches.
  await restart();
  solved = await sendPair(solved.map(result => ({ settings: {
    ...base, advancedMode: true, reuseSolution: true, ncr: 7,
    alphaDeg: 3, alphaRad: 3 * Math.PI / 180, reuseState: result.reuseState,
  } })));
  assert.ok(solved.every(result => result.converged));

  let inverse;
  for (let repeat = 0; repeat < 2; repeat += 1) {
    inverse = await sendPair(solved.map(result => ({ action: 'qdes', niter: 1,
      cpSpec: result.cpData.invAll.map(point => point.cp),
    })));
    assert.ok(inverse.every(result => result.debug.nanCount === 0));
  }
  await sendPair([0, 1].map(() => ({ action: 'dump', kind: 'cp' })));

  // The transferred inverse-design coordinates remain valid solver inputs.
  const edited = await sendPair(inverse.map(result => ({ settings: {
    ...base, source: 'custom', geometryKey: 'ownership-edited',
    custom: { name: 'Edited', nb: result.nb, x: result.xb, y: result.yb },
    advancedMode: true, reusePanel: false, reuseSolution: false,
  } })));
  assert.ok(edited.every(result => result.converged));
  const inviscid = await sendPair([0, 1].map(() => ({ settings: { ...base, viscous: false } })));
  assert.ok(inviscid.every(result => result.converged && result.reuseState === null));
  console.log(`worker_lifecycle.test.js: OK (${comparisons} exact copy/transfer comparisons)`);
} finally {
  await Promise.all(workers.map(worker => worker.terminate()));
}
