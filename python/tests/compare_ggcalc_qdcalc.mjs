import fs from 'fs';

import '../../js/naca.js';
import { pangen, tecalc } from '../../js/xfoil.js';
import { scalc, segspl } from '../../js/spline.js';
import { apcalc, ggcalc, ncalc, qdcalc, xywake } from '../../js/xpanel.js';
import { createMatrix } from '../../js/arrays.js';

const noop = () => {};
console.log = noop;
console.warn = noop;

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const { ctx: ctxPayload } = payload;

let ctx = null;
let n = 0;
let total = 0;

if (ctxPayload) {
  n = ctxPayload.N;
  const nw = ctxPayload.NW;
  total = n + nw;
  ctx = {
    N: n,
    NW: nw,
    WAKLEN: ctxPayload.WAKLEN ?? 0.0,
    X: new Float64Array(total),
    Y: new Float64Array(total),
    XP: new Float64Array(total),
    YP: new Float64Array(total),
    S: new Float64Array(total),
    NX: new Float64Array(total),
    NY: new Float64Array(total),
    APANEL: new Float64Array(total),
    SHARP: ctxPayload.SHARP,
    PI: Math.PI,
    ANTE: ctxPayload.ANTE,
    ASTE: ctxPayload.ASTE,
    DSTE: ctxPayload.DSTE,
    GAM: new Float64Array(n + 1),
    GAMU: Array.from({ length: n + 1 }, () => new Float64Array(2)),
    SIG: new Float64Array(total),
    QINVU: Array.from({ length: total + 1 }, () => new Float64Array(2)),
    AIJ: createMatrix(n + 1, n + 1),
    BIJ: createMatrix(n + 1, total),
    AIJPIV: new Int32Array(n + 1),
    QOPI: 1.0 / (4.0 * Math.PI),
    HOPI: 1.0 / (2.0 * Math.PI),
    ALFA: ctxPayload.ALFA,
    QINF: ctxPayload.QINF,
    LIMAGE: false,
    YIMAGE: 0.0,
    XTE: ctxPayload.XTE,
    YTE: ctxPayload.YTE,
    DZDG: new Float64Array(n + 1),
    DZDN: new Float64Array(n + 1),
    DQDG: new Float64Array(n + 1),
    DZDM: new Float64Array(total),
    DQDM: new Float64Array(total),
    LWAKE: true,
    LWDIJ: false,
    LADIJ: false,
    SNEW: new Float64Array(total),
  };

  for (let i = 0; i < total; i += 1) {
    ctx.X[i] = ctxPayload.X[i];
    ctx.Y[i] = ctxPayload.Y[i];
    ctx.S[i] = ctxPayload.S[i];
    ctx.NX[i] = ctxPayload.NX[i];
    ctx.NY[i] = ctxPayload.NY[i];
    ctx.APANEL[i] = ctxPayload.APANEL[i];
    ctx.XP[i] = ctxPayload.XP[i];
    ctx.YP[i] = ctxPayload.YP[i];
  }

  for (let i = 0; i < n; i += 1) {
    ctx.GAM[i] = i < Math.floor(n / 2) ? 1.0 : -1.0;
  }
} else {
  const { ides, params, alphaRad, minf, waklen } = payload;
  const { naca4, naca5 } = globalThis.Naca;
  const nside = payload.nside;
  const xx = new Float64Array(nside);
  const yt = new Float64Array(nside);
  const yc = new Float64Array(nside);
  const xb = new Float64Array(2 * nside);
  const yb = new Float64Array(2 * nside);

  let nb = 0;
  if (ides <= 9999) {
    nb = naca4(ides, xx, yt, yc, nside, xb, yb).nb;
  } else {
    nb = naca5(ides, xx, yt, yc, nside, xb, yb).nb;
  }

  const panel = pangen(xb.subarray(0, nb), yb.subarray(0, nb), nb, params);
  n = panel.n;
  const nw = Math.floor(n / 12) + 10 * Math.floor(waklen);
  total = n + nw;

  ctx = {
    N: n,
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
    GAM: new Float64Array(n + 1),
    GAMU: Array.from({ length: n + 1 }, () => new Float64Array(2)),
    SIG: new Float64Array(total),
    QINVU: Array.from({ length: total + 1 }, () => new Float64Array(2)),
    AIJ: createMatrix(n + 1, n + 1),
    BIJ: createMatrix(n + 1, total),
    AIJPIV: new Int32Array(n + 1),
    QOPI: 1.0 / (4.0 * Math.PI),
    HOPI: 1.0 / (2.0 * Math.PI),
    ALFA: alphaRad,
    QINF: 1.0,
    LIMAGE: false,
    YIMAGE: 0.0,
    XTE: 0.0,
    YTE: 0.0,
    DZDG: new Float64Array(n + 1),
    DZDN: new Float64Array(n + 1),
    DQDG: new Float64Array(n + 1),
    DZDM: new Float64Array(total),
    DQDM: new Float64Array(total),
    LWAKE: false,
    LWDIJ: false,
    LADIJ: false,
    SNEW: new Float64Array(total),
  };

  for (let i = 0; i < n; i += 1) {
    ctx.X[i] = panel.x[i];
    ctx.Y[i] = panel.y[i];
  }

  scalc(ctx.X, ctx.Y, ctx.S, n);
  segspl(ctx.X, ctx.XP, ctx.S, n);
  segspl(ctx.Y, ctx.YP, ctx.S, n);
  ncalc(ctx.X, ctx.Y, ctx.S, n, ctx.NX, ctx.NY);
  ctx.XTE = 0.5 * (ctx.X[0] + ctx.X[n - 1]);
  ctx.YTE = 0.5 * (ctx.Y[0] + ctx.Y[n - 1]);
  tecalc(ctx);
  apcalc(ctx);

  for (let i = 0; i < n; i += 1) {
    ctx.GAM[i] = i < Math.floor(n / 2) ? 1.0 : -1.0;
  }

  xywake(ctx);
}
ggcalc(ctx);
qdcalc(ctx);

function metricsMatrix(mat, n1, n2, samples) {
  let sum = 0.0;
  let sumsq = 0.0;
  let maxabs = 0.0;
  for (let i = 0; i < n1; i += 1) {
    const row = mat[i];
    for (let j = 0; j < n2; j += 1) {
      const val = row[j];
      sum += val;
      sumsq += val * val;
      maxabs = Math.max(maxabs, Math.abs(val));
    }
  }
  const sampleVals = samples.map(([ri, ci]) => {
    const i = ri - 1;
    const j = ci - 1;
    return (mat[i] && mat[i][j]) ?? 0.0;
  });
  return { sum, sumsq, maxabs, samples: sampleVals };
}

function metricsMatrix1Based(mat, n1, n2, samples) {
  let sum = 0.0;
  let sumsq = 0.0;
  let maxabs = 0.0;
  for (let i = 1; i <= n1; i += 1) {
    const row = mat[i];
    for (let j = 1; j <= n2; j += 1) {
      const val = row[j];
      sum += val;
      sumsq += val * val;
      maxabs = Math.max(maxabs, Math.abs(val));
    }
  }
  const sampleVals = samples.map(([ri, ci]) => mat[ri][ci]);
  return { sum, sumsq, maxabs, samples: sampleVals };
}

const aiSamples = [
  [1, 1],
  [1, n],
  [n, 1],
  [n, n],
  [n + 1, 1],
  [1, n + 1],
  [n + 1, n + 1],
  [Math.floor((n + 1) / 2), Math.floor((n + 1) / 2)],
];
const diSamples = [
  [1, 1],
  [1, total],
  [total, 1],
  [total, total],
  [Math.floor((total + 1) / 2), Math.floor((total + 1) / 2)],
];

const aiMetrics = metricsMatrix(ctx.AIJ, n + 1, n + 1, aiSamples);
const diMetrics = metricsMatrix1Based(ctx.DIJ, total, total, diSamples);

const results = {
  n,
  total,
  ai: aiMetrics,
  di: diMetrics,
};

process.stdout.write(JSON.stringify({ results }));
