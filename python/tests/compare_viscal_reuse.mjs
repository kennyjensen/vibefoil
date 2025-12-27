import fs from 'fs';

import { pangen, tecalc } from '../../js/xfoil.js';
import { scalc, segspl } from '../../js/spline.js';
import { apcalc, ggcalc, ncalc } from '../../js/xpanel.js';
import { createMatrix } from '../../js/arrays.js';
import { buildBlContext, viscal } from '../../js/xoper.js';

const noop = () => {};
console.log = noop;
console.warn = noop;

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const {
  xb,
  yb,
  nb,
  params,
  alphaRad,
  reinf,
  minf,
  waklen,
  ncrit,
} = payload;
const niter = 5;

const xbArr = Float64Array.from(xb);
const ybArr = Float64Array.from(yb);
const panel = pangen(xbArr, ybArr, nb, params);
const n = panel.n;
const nw = Math.floor(n / 12) + 10 * Math.floor(waklen);
const total = n + nw;

const ctxPanel = {
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
  GAM_A: new Float64Array(n + 1),
  QINV: new Float64Array(total + 1),
  QINV_A: new Float64Array(total + 1),
  QVIS: new Float64Array(total + 1),
  GAMU: Array.from({ length: n + 1 }, () => new Float64Array(2)),
  SIG: new Float64Array(total),
  QF0: new Float64Array(n),
  QF1: new Float64Array(n),
  QF2: new Float64Array(n),
  QF3: new Float64Array(n),
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
  ctxPanel.X[i] = panel.x[i];
  ctxPanel.Y[i] = panel.y[i];
}

scalc(ctxPanel.X, ctxPanel.Y, ctxPanel.S, n);
segspl(ctxPanel.X, ctxPanel.XP, ctxPanel.S, n);
segspl(ctxPanel.Y, ctxPanel.YP, ctxPanel.S, n);
ncalc(ctxPanel.X, ctxPanel.Y, ctxPanel.S, n, ctxPanel.NX, ctxPanel.NY);
ctxPanel.XTE = 0.5 * (ctxPanel.X[0] + ctxPanel.X[n - 1]);
ctxPanel.YTE = 0.5 * (ctxPanel.Y[0] + ctxPanel.Y[n - 1]);
tecalc(ctxPanel);
apcalc(ctxPanel);
ggcalc(ctxPanel);
const cosa = Math.cos(alphaRad);
const sina = Math.sin(alphaRad);
for (let i = 0; i < n; i += 1) {
  ctxPanel.GAM[i] = cosa * ctxPanel.GAMU[i][0] + sina * ctxPanel.GAMU[i][1];
}

const blCtx = buildBlContext(n, ctxPanel, ncrit);
blCtx.MINF = minf;
ctxPanel.QINF = 1.0;

viscal(blCtx, ctxPanel, alphaRad, reinf, { maxIter: niter, reuseSolution: true });
viscal(blCtx, ctxPanel, alphaRad, reinf, { maxIter: niter, reuseSolution: true });

const maxNbl = Math.max(blCtx.NBL[1], blCtx.NBL[2]);
function extractSide(mat, is) {
  const out = new Array(maxNbl + 1).fill(0.0);
  for (let ibl = 1; ibl <= blCtx.NBL[is]; ibl += 1) {
    out[ibl] = mat[ibl][is];
  }
  return out;
}

const results = {
  CL: blCtx.CL,
  CM: blCtx.CM,
  CD: blCtx.CD,
  CDF: blCtx.CDF,
  CDP: blCtx.CDP,
  NBL: Array.from(blCtx.NBL),
  IBLTE: Array.from(blCtx.IBLTE),
  ITRAN: Array.from(blCtx.ITRAN),
  XSSITR: Array.from(blCtx.XSSITR),
  THET: { 1: extractSide(blCtx.THET, 1), 2: extractSide(blCtx.THET, 2) },
  DSTR: { 1: extractSide(blCtx.DSTR, 1), 2: extractSide(blCtx.DSTR, 2) },
  CTAU: { 1: extractSide(blCtx.CTAU, 1), 2: extractSide(blCtx.CTAU, 2) },
  UEDG: { 1: extractSide(blCtx.UEDG, 1), 2: extractSide(blCtx.UEDG, 2) },
  MASS: { 1: extractSide(blCtx.MASS, 1), 2: extractSide(blCtx.MASS, 2) },
  TAU: { 1: extractSide(blCtx.TAU, 1), 2: extractSide(blCtx.TAU, 2) },
  DIS: { 1: extractSide(blCtx.DIS, 1), 2: extractSide(blCtx.DIS, 2) },
  CTQ: { 1: extractSide(blCtx.CTQ, 1), 2: extractSide(blCtx.CTQ, 2) },
  DELT: { 1: extractSide(blCtx.DELT, 1), 2: extractSide(blCtx.DELT, 2) },
  TSTR: { 1: extractSide(blCtx.TSTR, 1), 2: extractSide(blCtx.TSTR, 2) },
};

process.stdout.write(JSON.stringify({ results }));
