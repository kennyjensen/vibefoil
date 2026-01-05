import fs from 'fs';

import { MIXED } from '../../js/xqdes.js';
import { createMatrix } from '../../js/arrays.js';

const noop = () => {};
console.log = noop;
console.warn = noop;

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const {
  n,
  nsp,
  x,
  y,
  s,
  xp,
  yp,
  nx,
  ny,
  apanel,
  qspec,
  sspec,
  gam,
  sig,
  alfa,
  minf,
  qinf,
  xcmref,
  ycmref,
  psio,
  qdof0,
  qdof1,
  qdof2,
  qdof3,
  iq1,
  iq2,
  lcpXX,
  limage,
  yimage,
  sharp,
  ante,
  aste,
  dste,
  xte,
  yte,
  niterq,
  samples,
} = payload;

const ctx = {
  N: n,
  NSP: nsp,
  X: new Float64Array(n),
  Y: new Float64Array(n),
  S: new Float64Array(n),
  XP: new Float64Array(n),
  YP: new Float64Array(n),
  NX: new Float64Array(n),
  NY: new Float64Array(n),
  APANEL: new Float64Array(n),
  QSPEC: createMatrix(nsp, 1),
  QSPECP: createMatrix(nsp, 1),
  SSPEC: new Float64Array(nsp),
  GAM: new Float64Array(n),
  GAMU: Array.from({ length: n + 1 }, () => new Float64Array(2)),
  GAM_A: new Float64Array(n),
  SIG: new Float64Array(n),
  QF0: new Float64Array(n),
  QF1: new Float64Array(n),
  QF2: new Float64Array(n),
  QF3: new Float64Array(n),
  Q: createMatrix(n + 5, n + 5),
  DQ: new Float64Array(n + 5),
  DZDG: new Float64Array(n),
  DZDN: new Float64Array(n),
  DQDG: new Float64Array(n),
  DZDM: new Float64Array(n),
  DQDM: new Float64Array(n),
  QOPI: 1.0 / (4.0 * Math.PI),
  HOPI: 1.0 / (2.0 * Math.PI),
  PI: Math.PI,
  ALFA: alfa,
  MINF: minf,
  QINF: qinf,
  XCMREF: xcmref,
  YCMREF: ycmref,
  PSIO: psio,
  QDOF0: qdof0,
  QDOF1: qdof1,
  QDOF2: qdof2,
  QDOF3: qdof3,
  IQ1: iq1,
  IQ2: iq2,
  LCPXX: lcpXX,
  LIMAGE: limage,
  YIMAGE: yimage,
  SHARP: sharp,
  ANTE: ante,
  ASTE: aste,
  DSTE: dste,
  XTE: xte,
  YTE: yte,
  Z_QINF: 0.0,
  Z_ALFA: 0.0,
  Z_QDOF0: 0.0,
  Z_QDOF1: 0.0,
  Z_QDOF2: 0.0,
  Z_QDOF3: 0.0,
  CL: 0.0,
  CM: 0.0,
  CDP: 0.0,
  CL_ALF: 0.0,
  CL_MSQ: 0.0,
  COSA: Math.cos(alfa),
  SINA: Math.sin(alfa),
};

for (let i = 0; i < n; i += 1) {
  ctx.X[i] = x[i];
  ctx.Y[i] = y[i];
  ctx.S[i] = s[i];
  ctx.XP[i] = xp[i];
  ctx.YP[i] = yp[i];
  ctx.NX[i] = nx[i];
  ctx.NY[i] = ny[i];
  ctx.APANEL[i] = apanel[i];
  ctx.GAM[i] = gam[i];
  ctx.SIG[i] = sig[i];
}

for (let i = 0; i < nsp; i += 1) {
  ctx.SSPEC[i] = sspec[i];
  ctx.QSPEC[i][0] = qspec[i];
}

MIXED(ctx, 0, niterq);

function metricsArray(arr, n1, samples1) {
  let sum = 0.0;
  let sumsq = 0.0;
  let maxabs = 0.0;
  for (let i = 0; i < n1; i += 1) {
    const val = arr[i];
    sum += val;
    sumsq += val * val;
    maxabs = Math.max(maxabs, Math.abs(val));
  }
  const sampleVals = samples1.map((idx1) => {
    const idx = idx1 - 1;
    return arr[idx] ?? 0.0;
  });
  return { sum, sumsq, maxabs, samples: sampleVals };
}

const metrics = {
  x: metricsArray(ctx.X, n, samples),
  y: metricsArray(ctx.Y, n, samples),
  gam: metricsArray(ctx.GAM, n, samples),
};

const results = {
  metrics,
  psio: ctx.PSIO,
  qdof0: ctx.QDOF0,
  qdof1: ctx.QDOF1,
  qdof2: ctx.QDOF2,
  qdof3: ctx.QDOF3,
  cl: ctx.CL,
  cm: ctx.CM,
  cdp: ctx.CDP,
  clAlf: ctx.CL_ALF,
  clMsq: ctx.CL_MSQ,
};

process.stdout.write(JSON.stringify({ results }));
