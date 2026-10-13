// JSON-lines adapter around the actual website worker. No solver substitutions.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createInterface } from 'node:readline';

const backend = process.argv[2] ?? 'auto';
assert.ok(['auto', 'javascript'].includes(backend), 'Unknown backend');
if (backend === 'javascript') globalThis.WebAssembly = undefined;
else {
  const { createWasmBlsolv } = await import('../../js/blsolv_wasm.js');
  assert.equal(typeof createWasmBlsolv(), 'function', 'Auto-backend audit requires WebAssembly/SIMD');
  const { createWasmBoundaryLayer } = await import('../../js/boundary_layer_wasm.js');
  assert.equal(typeof createWasmBoundaryLayer(), 'function', 'Auto-backend audit requires the boundary-layer kernel');
  const { createWasmPanel } = await import('../../js/panel_wasm.js');
  assert.ok(createWasmPanel(), 'Auto-backend audit requires the panel kernel');
}
const { getBoundaryLayerWasmStats } = await import('../../js/boundary_layer_wasm.js');
const { getPanelWasmStats } = await import('../../js/panel_wasm.js');
const { cpcalc } = await import('../../js/xfoil.js');
let response;
let logs = [];
console.log = (...args) => logs.push(args.join(' '));
console.warn = (...args) => logs.push(args.join(' '));
globalThis.self = { postMessage(value) { response = value; } };
await import('../../js/solver_worker.js');

const blFields = ['XSSI', 'UEDG', 'THET', 'DSTR', 'CTAU', 'MASS', 'TSTR', 'TAU', 'DIS', 'CTQ', 'DELT'];
function geometry(foil) {
  if (foil.coordinates) {
    const rows = readFileSync(foil.coordinates, 'utf8').trim().split(/\r?\n/).slice(1)
      .filter(line => line.trim()).map(line => line.trim().split(/\s+/).map(Number));
    assert.ok(rows.every(row => row.length === 2 && row.every(Number.isFinite)));
    return { source: 'database', custom: {
      name: foil.id, nb: rows.length, x: rows.map(row => row[0]), y: rows.map(row => row[1]),
    } };
  }
  const designation = String(foil.naca);
  if (designation.length === 5) return {
    source: 'naca', mode: '5', series5: designation.slice(0, 3), t5: foil.naca % 100,
  };
  return { source: 'naca', mode: '4', m: Math.floor(foil.naca / 1000),
    p: Math.floor(foil.naca / 100) % 10, t: foil.naca % 100 };
}

function run(input) {
  logs = [];
  response = null;
  const settings = { ...geometry(input.foil),
    // Unique geometry keys prevent automatic reuse between independent points.
    geometryKey: input.id, advancedMode: false, reusePanel: false,
    reuseSolution: false, reuseState: null,
    alphaDeg: input.alpha, alphaRad: input.alpha * Math.PI / 180,
    viscous: input.viscous, mach: input.mach, reynolds: input.reynolds,
    ncr: input.ncrit, nIter: input.iterations,
    residualTolerance: input.residual_tolerance,
    flap: { x: 0.75, yrel: 0.5, deflection: 0 }, canvasWidth: 800, canvasHeight: 500,
  };
  const before = getBoundaryLayerWasmStats();
  const beforePanel = getPanelWasmStats();
  self.onmessage({ data: { requestId: input.id, settings } });
  const after = getBoundaryLayerWasmStats();
  const boundaryLayer = Object.fromEntries(Object.keys(after).map(key => [key, after[key] - before[key]]));
  const panel = Object.fromEntries(Object.entries(getPanelWasmStats()).map(([key, value]) => [key, value - beforePanel[key]]));
  assert.equal(response?.ok, true, response?.error ?? 'No worker result');
  const c = response.coeffs;
  const result = { N: response.nb, NW: 0, converged: response.converged,
    coefficients: { CL: c.CL, CM: c.CM, CDP_SURFACE: c.CDP },
    panel: response.cpData.invAll.map(point => [point.x, point.y, point.cp]),
  };
  if (input.viscous) {
    const { panel, bl } = response.reuseState;
    assert.equal(panel.N, response.nb);
    assert.equal(bl.LVCONV, response.converged);
    result.NW = panel.NW;
    Object.assign(result.coefficients, { CD: c.CD, CDF: c.CDF, CDP: c.CD - c.CDF });
    result.operating = { mach: bl.MINF, reynolds: bl.REINF, ncrit: Array.from(bl.ACRIT).slice(1, 3) };
    result.transition = { x: [c.XTRT, c.XTRB], index: Array.from(bl.ITRAN).slice(1, 3),
      forced: bl.TFORCE.slice(1, 3) };
    result.NBL = Array.from(bl.NBL).slice(1, 3);
    result.IBLTE = Array.from(bl.IBLTE).slice(1, 3);
    const cp = [...response.cpData.upper, ...response.cpData.lower, ...response.cpData.wake];
    assert.equal(cp.length, panel.N + panel.NW);
    const cpi = cpcalc(panel.QINV, panel.QINF, bl.MINF);
    result.panel = cp.map((point, i) => {
      assert.equal(point.x, panel.X[i]);
      assert.equal(point.y, panel.Y[i]);
      return [point.x, point.y, panel.QINV[i + 1], panel.QVIS[i + 1], cpi[i + 1], point.cp];
    });
    result.bl = [1, 2].map(side => Array.from({ length: bl.NBL[side] }, (_, row) => [
      bl.IPAN[row + 1][side], ...blFields.map(field => bl[field][row + 1][side]),
    ]));
  }
  if (backend === 'auto' && input.viscous) assert.ok(boundaryLayer.accelerated > 0,
    'Viscous comparison did not exercise the boundary-layer kernel');
  if (backend === 'auto' && input.viscous) assert.ok(panel.solves > 0,
    'Viscous comparison did not exercise the panel kernel');
  return { id: input.id, backend, result, logs, boundaryLayer, panel };
}

for await (const line of createInterface({ input: process.stdin, crlfDelay: Infinity })) {
  try {
    const output = run(JSON.parse(line));
    process.stdout.write(`${JSON.stringify(output)}\n`);
  } catch (error) {
    process.stdout.write(`${JSON.stringify({ error: error.stack, logs })}\n`);
  }
}
