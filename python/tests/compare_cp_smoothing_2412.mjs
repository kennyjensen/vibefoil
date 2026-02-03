import '../../js/naca.js';
import { pangen, cpcalc, tecalc } from '../../js/xfoil.js';
import { apcalc, ncalc, ggcalc, stfind } from '../../js/xpanel.js';
import { buildBlContext, viscal } from '../../js/xoper.js';
import { scalc, segspl } from '../../js/spline.js';
import { createMatrix } from '../../js/arrays.js';

function buildCpEditPoints(nb, cpData) {
  if (!cpData) return null;
  const invAll = Array.isArray(cpData.invAll) ? cpData.invAll : [];
  const source = invAll.length ? invAll : cpData.upper.concat(cpData.lower);
  if (!source.length) return null;
  const count = Math.min(nb, source.length);
  const lePt = cpData.le;
  const tePt = cpData.te;
  if (!lePt || !tePt) return null;
  const dxChord = tePt.x - lePt.x;
  const dyChord = tePt.y - lePt.y;
  const chord2 = dxChord * dxChord + dyChord * dyChord || 1.0;
  let leIdx = 0;
  let leDist = Infinity;
  for (let i = 0; i < count; i += 1) {
    const p = source[i];
    const dx = p.x - lePt.x;
    const dy = p.y - lePt.y;
    const dist = dx * dx + dy * dy;
    if (dist < leDist) {
      leDist = dist;
      leIdx = i;
    }
  }
  const points = new Array(count);
  for (let i = 0; i < count; i += 1) {
    const p = source[i];
    const s = ((p.x - lePt.x) * dxChord + (p.y - lePt.y) * dyChord) / chord2;
    const side = i <= leIdx ? 'upper' : 'lower';
    points[i] = { index: i, s, cp: p.cp, side };
  }
  return points;
}

function smoothCpPoints(points) {
  if (!points || points.length < 3) return points;
  const smoothSide = (side) => {
    const idxs = [];
    for (let i = 0; i < points.length; i += 1) {
      if (points[i].side === side) idxs.push(i);
    }
    if (idxs.length < 3) return;
    const nextVals = new Array(idxs.length);
    for (let i = 0; i < idxs.length; i += 1) {
      const idx = idxs[i];
      if (i === 0 || i === idxs.length - 1) {
        nextVals[i] = points[idx].cp;
        continue;
      }
      const prev = points[idxs[i - 1]].cp;
      const cur = points[idx].cp;
      const next = points[idxs[i + 1]].cp;
      nextVals[i] = (prev + 2.0 * cur + next) * 0.25;
    }
    for (let i = 0; i < idxs.length; i += 1) {
      points[idxs[i]].cp = nextVals[i];
    }
  };
  smoothSide('upper');
  smoothSide('lower');
  return points;
}

function calcLeSpike(points, side) {
  const sidePts = points.filter((p) => p.side === side);
  if (sidePts.length < 3) return { spike: null, count: sidePts.length, leIdx: -1 };
  let minIdx = 0;
  for (let i = 1; i < sidePts.length; i += 1) {
    if (sidePts[i].s < sidePts[minIdx].s) minIdx = i;
  }
  const candidates = [minIdx - 1, minIdx, minIdx + 1].filter((i) => i >= 1 && i <= sidePts.length - 2);
  let maxSpike = 0.0;
  for (const idx of candidates) {
    const prev = sidePts[idx - 1].cp;
    const cur = sidePts[idx].cp;
    const next = sidePts[idx + 1].cp;
    const spike = Math.abs(2.0 * cur - prev - next);
    if (spike > maxSpike) maxSpike = spike;
  }
  return { spike: maxSpike, count: sidePts.length, leIdx: minIdx };
}

function monotonicInfo(points, side) {
  const sidePts = points.filter((p) => p.side === side);
  if (sidePts.length < 2) return { ok: true, direction: 0 };
  let inc = true;
  let dec = true;
  for (let i = 1; i < sidePts.length; i += 1) {
    if (sidePts[i].s < sidePts[i - 1].s) inc = false;
    if (sidePts[i].s > sidePts[i - 1].s) dec = false;
  }
  return { ok: inc || dec, direction: inc ? 1 : (dec ? -1 : 0) };
}

function getSurfaceIndices(nb, ctxPanel) {
  const { ist } = stfind(ctxPanel, nb);
  const upperIdx = [];
  for (let i = ist; i >= 0; i -= 1) {
    upperIdx.push(i);
  }
  const lowerIdx = [];
  for (let i = ist + 1; i < nb; i += 1) {
    lowerIdx.push(i);
  }
  return { upperIdx, lowerIdx };
}

function buildPanelContext(xb, yb, alphaRad) {
  const nb = xb.length;
  const waklen = 1.0;
  const nw = Math.floor(nb / 12) + 10 * Math.floor(waklen);
  const total = nb + nw;
  const n1 = nb + 1;
  const ctx = {
    N: nb,
    NW: nw,
    WAKLEN: waklen,
    X: new Float64Array(total),
    Y: new Float64Array(total),
    XP: new Float64Array(total),
    YP: new Float64Array(total),
    S: new Float64Array(total),
    NX: new Float64Array(total),
    NY: new Float64Array(total),
    APANEL: new Float64Array(total),
    SHARP: true,
    PI: Math.PI,
    ANTE: 0.0,
    ASTE: 0.0,
    DSTE: 0.0,
    GAM: new Float64Array(n1),
    GAM_A: new Float64Array(n1),
    QINV: new Float64Array(total + 1),
    QINV_A: new Float64Array(total + 1),
    QVIS: new Float64Array(total + 1),
    GAMU: Array.from({ length: n1 }, () => new Float64Array(2)),
    SIG: new Float64Array(total),
    QF0: new Float64Array(nb),
    QF1: new Float64Array(nb),
    QF2: new Float64Array(nb),
    QF3: new Float64Array(nb),
    QINVU: Array.from({ length: total }, () => new Float64Array(2)),
    AIJ: createMatrix(n1, n1),
    BIJ: createMatrix(n1, total),
    AIJPIV: new Int32Array(n1),
    QOPI: 1.0 / (4.0 * Math.PI),
    HOPI: 1.0 / (2.0 * Math.PI),
    ALFA: alphaRad,
    QINF: 1.0,
    LIMAGE: false,
    YIMAGE: 0.0,
    XTE: 0.0,
    YTE: 0.0,
    DZDG: new Float64Array(nb),
    DZDN: new Float64Array(nb),
    DQDG: new Float64Array(nb),
    DZDM: new Float64Array(total),
    DQDM: new Float64Array(total),
    LWAKE: false,
    LWDIJ: false,
    LADIJ: false,
    SNEW: new Float64Array(total),
  };

  for (let i = 0; i < nb; i += 1) {
    ctx.X[i] = xb[i];
    ctx.Y[i] = yb[i];
  }

  scalc(ctx.X, ctx.Y, ctx.S, nb);
  segspl(ctx.X, ctx.XP, ctx.S, nb);
  segspl(ctx.Y, ctx.YP, ctx.S, nb);
  ncalc(ctx.X, ctx.Y, ctx.S, nb, ctx.NX, ctx.NY);
  ctx.XTE = 0.5 * (ctx.X[0] + ctx.X[nb - 1]);
  ctx.YTE = 0.5 * (ctx.Y[0] + ctx.Y[nb - 1]);
  tecalc(ctx);
  apcalc(ctx);
  ggcalc(ctx);
  return ctx;
}

function main() {
  const nside = 123;
  const xx = new Float64Array(nside);
  const yt = new Float64Array(nside);
  const yc = new Float64Array(nside);
  const xb = new Float64Array(2 * nside);
  const yb = new Float64Array(2 * nside);

  const ides = 2412;
  const res = globalThis.Naca.naca4(ides, xx, yt, yc, nside, xb, yb);
  const panelRes = pangen(xb.subarray(0, res.nb), yb.subarray(0, res.nb), res.nb);
  const nb = panelRes.n;

  const alphaRad = 3.0 * Math.PI / 180.0;
  const ctxPanel = buildPanelContext(panelRes.x, panelRes.y, alphaRad);
  const blCtx = buildBlContext(nb, ctxPanel);
  let viscalResult;
  const prevLog = console.log;
  console.log = () => {};
  try {
    viscalResult = viscal(blCtx, ctxPanel, alphaRad, 3.0e6, { maxIter: 20, logSurface: false, reuseSolution: false });
  } finally {
    console.log = prevLog;
  }
  const qvis = viscalResult?.qvis;
  const cpVis = cpcalc(qvis, ctxPanel.QINF ?? 1.0, blCtx.MINF ?? 0.0);

  const { upperIdx, lowerIdx } = getSurfaceIndices(nb, ctxPanel);
  const cpUpper = upperIdx.map((idx) => ({ x: ctxPanel.X[idx], y: ctxPanel.Y[idx], cp: cpVis[idx + 1] }));
  const cpLower = lowerIdx.map((idx) => ({ x: ctxPanel.X[idx], y: ctxPanel.Y[idx], cp: cpVis[idx + 1] }));

  const lePt = { x: ctxPanel.XLE, y: ctxPanel.YLE };
  const tePt = { x: ctxPanel.XTE, y: ctxPanel.YTE };
  if (!Number.isFinite(lePt.x) || !Number.isFinite(lePt.y) || !Number.isFinite(tePt.x) || !Number.isFinite(tePt.y)) {
    let leIdx = 0;
    let teIdx = 0;
    for (let i = 1; i < nb; i += 1) {
      if (ctxPanel.X[i] < ctxPanel.X[leIdx]) leIdx = i;
      if (ctxPanel.X[i] > ctxPanel.X[teIdx]) teIdx = i;
    }
    lePt.x = ctxPanel.X[leIdx];
    lePt.y = ctxPanel.Y[leIdx];
    tePt.x = ctxPanel.X[teIdx];
    tePt.y = ctxPanel.Y[teIdx];
  }

  const cpData = {
    upper: cpUpper,
    lower: cpLower,
    invAll: [],
    le: lePt,
    te: tePt,
  };

  const points = buildCpEditPoints(nb, cpData);
  const beforeUpper = calcLeSpike(points, 'upper');
  const beforeLower = calcLeSpike(points, 'lower');
  const monoUpper = monotonicInfo(points, 'upper');
  const monoLower = monotonicInfo(points, 'lower');

  const smoothed = points.map((p) => ({ ...p }));
  smoothCpPoints(smoothed);
  const afterUpper = calcLeSpike(smoothed, 'upper');
  const afterLower = calcLeSpike(smoothed, 'lower');

  const results = {
    nb,
    converged: !!viscalResult?.converged,
    before: { upper: beforeUpper, lower: beforeLower },
    after: { upper: afterUpper, lower: afterLower },
    monotonic: { upper: monoUpper, lower: monoLower },
  };
  process.stdout.write(JSON.stringify({ results }));
}

main();
