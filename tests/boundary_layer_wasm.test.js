import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { cpSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { layout } from '../js/boundary_layer_layout.js';
import { sourceHashes, boundaryLayerBase64 } from '../js/boundary_layer_wasm_data.js';
import { createMatrix1, createMatrix1Flat } from '../js/arrays.js';

for (const [file, hash] of Object.entries({...layout.hashes, ...sourceHashes})) {
  const bytes = readFileSync(new URL(`../${file}`, import.meta.url));
  assert.equal(createHash('sha256').update(bytes).digest('hex'), hash, `${file}: regenerate the boundary-layer kernel`);
}
const module = new WebAssembly.Module(Buffer.from(boundaryLayerBase64, 'base64'));
assert.ok(WebAssembly.Module.imports(module).every(i => i.module === 'math' || i.module === 'diagnostics'));
for (const [rows, cols] of [[0,0],[3,2],[271,2]]) {
  const flat = createMatrix1Flat(rows, cols), original = createMatrix1(rows, cols);
  for (let i = 0; i <= rows; i++) for (let j = 0; j <= cols; j++) flat[i][j] = original[i][j] = i * 100 + j;
  assert.deepStrictEqual(flat, original);
  assert.equal(flat[0].buffer, flat[rows].buffer);
}

// Instrument a disposable copy of the real worker only to intercept SETBL's
// input. The expected result comes from the maintained pure JS equations.
const directory = mkdtempSync(join(tmpdir(), 'vibefoil-bl-parity-'));
const output = console.log;
const warn = console.warn;
const clean = ctx => Object.fromEntries(Object.entries(ctx).filter(([,v]) => typeof v !== 'function'));
const clone = ctx => Object.assign(structuredClone(clean(ctx)), Object.fromEntries(Object.entries(ctx).filter(([,v]) => typeof v === 'function')));
let capture, response, seed, setblReference, blsolvReference, combinedComparisons=0, combinedAccelerated=0, comparisons = 0, accelerated = 0, rejected = 0, diagnosticComparisons = 0;
let forced = false;
try {
  cpSync(new URL('../js/', import.meta.url), join(directory, 'js'), {recursive: true});
  writeFileSync(join(directory, 'package.json'), '{"type":"module"}\n');
  const file = join(directory, 'js/xoper.js');
  const code = readFileSync(file, 'utf8');
  const marker = '    const solved = tryWasmBoundaryLayer(blCtx, true);\n    if (!solved) setbl(blCtx);\n    const ok = solved || blsolv(blCtx);';
  assert.equal(code.split(marker).length, 2);
  writeFileSync(file, code.replace(marker, '    globalThis.__checkBoundaryLayer(blCtx, setbl, blsolv);\n    const ok = blsolv(blCtx);'));
  const {createWasmBoundaryLayer} = await import(pathToFileURL(join(directory, 'js/boundary_layer_wasm.js')));
  const wasm = createWasmBoundaryLayer();
  assert.equal(typeof wasm, 'function');
  console.log = (...args) => capture?.push(args.join(' '));
  console.warn = (...args) => capture?.push(args.join(' '));
  globalThis.__checkBoundaryLayer = (ctx, reference, matrixReference = blsolvReference) => {
    blsolvReference = matrixReference;
    if (forced) { ctx.XSTRIP[1] = 0.3; ctx.XSTRIP[2] = 0.6; }
    seed ??= clone(ctx);
    setblReference = reference;
    const before = clone(ctx), expected = clone(ctx);
    capture = [];
    reference(expected);
    const expectedLogs = capture;
    capture = [];
    if (wasm(ctx)) accelerated++;
    else {
      rejected++;
      assert.deepStrictEqual(clean(ctx), clean(before), 'Rejected kernel changed live state');
      assert.deepStrictEqual(capture, [], 'Rejected kernel leaked diagnostics');
      reference(ctx);
    }
    assert.deepStrictEqual(clean(ctx), clean(expected), `Entire BL state at iteration ${comparisons}`);
    assert.deepStrictEqual(capture, expectedLogs, `Diagnostics at iteration ${comparisons}`);
    comparisons++; diagnosticComparisons++;
    const expectedStep = clone(expected), actualStep = clone(before);
    capture=[];
    const saved = globalThis.WebAssembly;
    let ok;
    try {globalThis.WebAssembly=undefined;ok=matrixReference(expectedStep);} finally {globalThis.WebAssembly=saved;}
    const stepLogs=expectedLogs.concat(capture);
    capture=[];
    if (wasm(actualStep,true)) {assert.equal(ok,true);combinedAccelerated++;}
    else {
      assert.deepStrictEqual(clean(actualStep),clean(before),'Rejected combined step changed live state');
      assert.deepStrictEqual(capture,[]);reference(actualStep);assert.equal(matrixReference(actualStep),ok);
    }
    assert.deepStrictEqual(clean(actualStep),clean(expectedStep),`Complete assembled/solved state at iteration ${combinedComparisons}`);
    assert.deepStrictEqual(capture,stepLogs);combinedComparisons++;

  };
  globalThis.self = {postMessage(value) { response = value; }};
  await import(pathToFileURL(join(directory, 'js/solver_worker.js')));
  const base = {source:'naca',mode:'4',m:2,p:4,t:12,advancedMode:true,reusePanel:false,reuseSolution:false,
    viscous:true,mach:0.1,reynolds:3e6,ncr:9,nIter:40,residualTolerance:1e-8,
    alphaDeg:4,flap:{x:0.75,yrel:0.5,deflection:0},canvasWidth:800,canvasHeight:500};
  const manifest = JSON.parse(readFileSync(new URL('./canonical/cases.json', import.meta.url), 'utf8'));
  const coordinatePath = manifest.airfoils.find(foil => foil.id === 'e387').coordinates;
  const coordinates = readFileSync(new URL(`../${coordinatePath}`, import.meta.url), 'utf8').trim().split(/\r?\n/).slice(1).filter(line => line.trim()).map(line=>line.trim().split(/\s+/).map(Number));
  const e387 = {source:'database',custom:{name:'E387',nb:coordinates.length,x:coordinates.map(p=>p[0]),y:coordinates.map(p=>p[1])}};
  const cases = [
    {}, {m:0,p:0,alphaDeg:8,mach:0.3,reynolds:5e5,ncr:12},
    {m:4,p:4,alphaDeg:-4,mach:0,reynolds:5e6,ncr:5},
    {mode:'5',series5:'230',t5:12,alphaDeg:0,mach:0.15,reynolds:1e6},
    {flap:{x:0.75,yrel:0.5,deflection:10}},
    {alphaDeg:2,forced:true},
    {...e387,alphaDeg:-4,mach:0,reynolds:1e6,ncr:12,nIter:100},
    {...e387,alphaDeg:-4,mach:0.15,reynolds:1e6,ncr:9,nIter:20},
  ];
  let reuseState;
  for (const [i, settings] of cases.entries()) {
    forced = !!settings.forced;
    const input = {...base,...settings,geometryKey:`wasm-case-${i}`};
    input.alphaRad = input.alphaDeg * Math.PI / 180;
    self.onmessage({data:{requestId:i,settings:input}});
    assert.ok(response?.ok, response?.errorStack ?? response?.error);
    if (i === 0) assert.ok(response.converged);
    reuseState = response.reuseState;
  }
  forced = false;
  // Exercise live cached state across angle/Ncrit changes and reuse attempts.
  for (const alphaDeg of [2,4,6]) {
    const input = {...base,advancedMode:false,reusePanel:true,reuseSolution:true,reuseState,
      geometryKey:'wasm-reuse',alphaDeg,alphaRad:alphaDeg*Math.PI/180,ncr:alphaDeg===6?5:9};
    self.onmessage({data:{requestId:100+alphaDeg,settings:input}});
    assert.ok(response?.ok,response?.errorStack);reuseState=response.reuseState;
  }
  assert.ok(accelerated > 30, 'Must actually exercise the compiled equations');

  // Reuse one instance across contexts and memory growth, then invalidate its
  // cached layouts through replaced views and changed outer-array lengths.
  const contexts = [clone(seed), clone(seed)];
  const padding = new Float64Array(300000).fill(-0);
  const originalCom = contexts[1].COM1;
  padding.set(originalCom, 17);
  contexts[1].COM1 = padding.subarray(17, 17 + originalCom.length);
  for (let round = 0; round < 3; round++) for (const ctx of contexts) {
    if (round === 1) ctx.THET.push(new Float64Array(3));
    if (round === 2) ctx.COM2 = new Float64Array(ctx.COM2);
    globalThis.__checkBoundaryLayer(ctx, setblReference);
  }
  assert.ok(padding.subarray(0,17).every(value => Object.is(value,-0)));
  assert.ok(padding.subarray(17+originalCom.length).every(value => Object.is(value,-0)));
  {
    const ctx = contexts[0];
    structuredClone(ctx.COM1.buffer, {transfer:[ctx.COM1.buffer]});
    const snapshot = structuredClone(Object.fromEntries(Object.entries(clean(ctx)).filter(([k])=>k!=='COM1')));
    assert.equal(wasm(ctx),false,'Detached buffer in a cached plan');
    assert.deepStrictEqual(Object.fromEntries(Object.entries(clean(ctx)).filter(([k])=>k!=='COM1')),snapshot);
    assert.equal(ctx.COM1.byteLength,0);comparisons++;
  }
  {
    const ctx = contexts[1];
    Object.defineProperty(ctx._velocitySources,'count',{writable:false});
    const before=clone(ctx);
    assert.equal(wasm(ctx),false,'Read-only output in a cached plan');
    assert.deepStrictEqual(clean(ctx),clean(before));comparisons++;
  }

  assert.ok(combinedAccelerated > 200);
  {
    const ctx = clone(seed);
    globalThis.__checkBoundaryLayer(ctx, setblReference);
    // Cached metadata must survive a rejected kernel call only by being
    // republished. Fractional and oversized indices must never be truncated
    // into valid rows by the faster array-index conversion.
    for (const invalid of [1.5, -0.5, 4294967296]) {
      const saved = ctx.ISYS[2][1];
      ctx.ISYS[2][1] = invalid;
      const before = clone(ctx);capture=[];
      assert.equal(wasm(ctx, true), false, `Invalid array index ${invalid}`);
      assert.deepStrictEqual(clean(ctx), clean(before));
      assert.deepStrictEqual(capture, []);comparisons++;
      ctx.ISYS[2][1] = saved;
      globalThis.__checkBoundaryLayer(ctx, setblReference);
    }
    // A second layout can overwrite the resident metadata before rejecting
    // a scalar input. Returning to the first context must restore its layout.
    const other = clone(seed);
    const padded = new Float64Array(400000);
    padded.set(other.COM1, 31);
    other.COM1 = padded.subarray(31, 31 + other.COM1.length);
    other.CL = NaN;
    const before = clone(other);capture=[];
    assert.equal(wasm(other, true), false);
    assert.deepStrictEqual(clean(other), clean(before));
    assert.deepStrictEqual(capture, []);comparisons++;
    globalThis.__checkBoundaryLayer(ctx, setblReference);
  }
  {
    const ctx = clone(seed);
    // The geometry buffer is also a writable scratch buffer in this case.
    ctx.W1 = ctx.X;ctx.XSTRIP[1] = 0.3;
    globalThis.__checkBoundaryLayer(ctx, setblReference);
    globalThis.__checkBoundaryLayer(ctx, setblReference);
    let reads=0;const row=ctx.THET[2];
    Object.defineProperty(ctx.THET,'2',{get(){reads++;return row;},configurable:true});
    for (const solve of [false,true]) {
      assert.equal(wasm(ctx,solve),false);assert.equal(reads,0);comparisons++;
    }
  }
  for (const mutate of [
    ctx=>{ctx.VAF=new Float64Array(ctx.VAF);},
    ctx=>{ctx.VACCEL=undefined;},
    ctx=>{ctx.S=ctx.VMF.subarray(0,ctx.S.length);},
  ]) {
    const ctx=clone(seed);mutate(ctx);const before=clone(ctx);capture=[];
    assert.equal(wasm(ctx,true),false,'Unsupported combined matrix solve');
    assert.deepStrictEqual(clean(ctx),clean(before));assert.deepStrictEqual(capture,[]);comparisons++;
  }
  // Storage failures, nonfinite inputs and traps must leave every byte/value
  // untouched. In particular a cached layout must not hide a replaced array.
  const mutations = [
    ctx => { ctx.NSYS = 2048; },
    ctx => { ctx.CL = NaN; },
    ctx => { ctx.THET[2][1] = Infinity; },
    ctx => { ctx.COM1 = new Float32Array(ctx.COM1); },
    ctx => { ctx.COM1 = new Int32Array(ctx.COM1); },
    ctx => { ctx.COM1 = Array.from(ctx.COM1);ctx.COM1[1] = false; },
    ctx => { ctx.TFORCE[1] = 0; },
    ctx => { ctx.CL = true; },
    ctx => { ctx.blkin = () => {}; },
    ctx => { ctx.ISYS[2][1] = -1; },
    ctx => { ctx.COM1 = new Float64Array(new SharedArrayBuffer(ctx.COM1.byteLength)); },
  ];
  for (const mutate of mutations) {
    const ctx = clone(seed);mutate(ctx);const before=clone(ctx);capture=[];
    assert.equal(wasm(ctx),false);
    assert.deepStrictEqual(clean(ctx),clean(before));assert.deepStrictEqual(capture,[]);
    comparisons++;
  }
  {
    const ctx = clone(seed);
    ctx.N = 2501;ctx.LBLINI = false;ctx.XSTRIP[1] = 0.3;
    for (const name of ['X','Y','S','W1','W2','W3','W4']) {
      const data = new Float64Array(ctx.N + 1);
      data.set(ctx[name]);ctx[name] = data;
    }
    const before = clone(ctx);capture=[];
    assert.equal(wasm(ctx),false,'Original spline size guard must remain effective');
    assert.deepStrictEqual(clean(ctx),clean(before));assert.deepStrictEqual(capture,[]);
    assert.throws(()=>setblReference(ctx),/SPLIND: array overflow/);comparisons++;
  }
  {
    const ctx=clone(seed);let reads=0;
    Object.defineProperty(ctx,'NSYS',{get(){reads++;return 300;},configurable:true});
    assert.equal(wasm(ctx),false);assert.equal(reads,0);comparisons++;
  }
  const savedWasm = globalThis.WebAssembly;
  try {
    globalThis.WebAssembly=undefined;assert.equal(createWasmBoundaryLayer(),null);
    globalThis.WebAssembly={Module(){throw new Error('blocked');}};assert.equal(createWasmBoundaryLayer(),null);
  } finally {globalThis.WebAssembly=savedWasm;}
  // Inject a nonfinite result inside the kernel, after valid input packing.
  const savedPow=Math.pow;
  let broken;
  try {Math.pow=()=>NaN;broken=createWasmBoundaryLayer();} finally {Math.pow=savedPow;}
  const bad=clone(seed), before=clone(bad);capture=[];
  assert.equal(broken(bad),false);assert.deepStrictEqual(clean(bad),clean(before));assert.deepStrictEqual(capture,[]);
  comparisons++;
  output(`combined: ${combinedComparisons} exact checks; ${combinedAccelerated} accelerated steps`);
  output(`boundary_layer_wasm: ${comparisons} complete-state/failure checks; ${accelerated} accelerated iterations, ${rejected} safe fallbacks; ${diagnosticComparisons} exact diagnostic comparisons`);
} finally {
  console.log=output;console.warn=warn;
  delete globalThis.__checkBoundaryLayer;
  delete globalThis.self;
  rmSync(directory,{recursive:true,force:true});
}
