import fs from 'fs';

import { ggcalc, qdcalc, qwcalc, qiset, stfind, iblpan, xicalc, uicalc, qvfue, gamqv, stmove, xywake } from '../../js/xpanel.js';
import { buildBlContext } from '../../js/xoper.js';
import { setbl, iblsys } from '../../js/xbl.js';
import { createMatrix, createMatrix1, createTensor3 } from '../../js/arrays.js';

const noop = () => {};
console.log = noop;
console.warn = noop;

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const { ctx: ctxPayload } = payload;

if (!ctxPayload) {
  throw new Error('ctx payload required');
}

const n = ctxPayload.N;
const nw = ctxPayload.NW;
const total = n + nw;
const alphaRad = ctxPayload.ALFA;
const minf = ctxPayload.MINF;
const reinf = ctxPayload.REINF;
const ncrit = ctxPayload.NCRIT ?? 9.0;

const ctxPanel = {
  N: n,
  NW: nw,
  WAKLEN: ctxPayload.WAKLEN,
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
  QINF: ctxPayload.QINF,
  MINF: minf,
  LIMAGE: false,
  YIMAGE: 0.0,
  XTE: ctxPayload.XTE,
  YTE: ctxPayload.YTE,
  CHORD: ctxPayload.CHORD,
  DZDG: new Float64Array(n + 1),
  DZDN: new Float64Array(n + 1),
  DQDG: new Float64Array(n + 1),
  DZDM: new Float64Array(total),
  DQDM: new Float64Array(total),
  LWAKE: false,
  LWDIJ: false,
  LADIJ: false,
  SNEW: new Float64Array(total),
  XCMREF: ctxPayload.XCMREF,
  YCMREF: ctxPayload.YCMREF,
};

for (let i = 0; i < total; i += 1) {
  ctxPanel.X[i] = ctxPayload.X[i];
  ctxPanel.Y[i] = ctxPayload.Y[i];
  ctxPanel.S[i] = ctxPayload.S[i];
  ctxPanel.NX[i] = ctxPayload.NX[i];
  ctxPanel.NY[i] = ctxPayload.NY[i];
  ctxPanel.APANEL[i] = ctxPayload.APANEL[i];
  ctxPanel.XP[i] = ctxPayload.XP[i];
  ctxPanel.YP[i] = ctxPayload.YP[i];
}

ggcalc(ctxPanel);

for (let i = 1; i <= n; i += 1) {
  ctxPanel.GAM[i - 1] = ctxPayload.GAM[i];
}
for (let i = 1; i <= total; i += 1) {
  ctxPanel.QINVU[i - 1][0] = ctxPayload.QINVU1[i];
  ctxPanel.QINVU[i - 1][1] = ctxPayload.QINVU2[i];
}

const blCtx = buildBlContext(n, ctxPanel, ncrit);
blCtx.MINF = minf;
blCtx.REINF = reinf;
blCtx.MINF1 = minf;
blCtx.REINF1 = reinf;
blCtx.HVRAT = 0.0;
blCtx.QINF = ctxPayload.QINF;
blCtx.LVISC = true;
blCtx.ALFA = alphaRad;
blCtx.ADEG = alphaRad / blCtx.DTOR;

function ensureViscousArrays(blCtxLocal) {
  const nsys = blCtxLocal.NSYS;
  if (!Number.isFinite(nsys) || nsys <= 0) return;

  const needsTensor = !blCtxLocal.VA
    || blCtxLocal.VA.length <= 3
    || blCtxLocal.VA[1]?.length <= 2
    || blCtxLocal.VA[1]?.[1]?.length <= nsys;
  const needsVm = !blCtxLocal.VM
    || blCtxLocal.VM.length <= 3
    || blCtxLocal.VM[1]?.length <= nsys
    || blCtxLocal.VM[1]?.[1]?.length <= nsys;

  if (needsTensor) {
    blCtxLocal.VA = createTensor3(3, 2, nsys);
    blCtxLocal.VB = createTensor3(3, 2, nsys);
    blCtxLocal.VDEL = createTensor3(3, 2, nsys);
  }
  if (needsVm) {
    blCtxLocal.VM = createTensor3(3, nsys, nsys);
  }
  if (!blCtxLocal.VZ || blCtxLocal.VZ.length <= 3) {
    blCtxLocal.VZ = createMatrix1(3, 2);
  }
  if (!blCtxLocal.XOCTR || blCtxLocal.XOCTR.length < 3) {
    blCtxLocal.XOCTR = new Float64Array(3);
  }
  if (!blCtxLocal.YOCTR || blCtxLocal.YOCTR.length < 3) {
    blCtxLocal.YOCTR = new Float64Array(3);
  }
  if (!blCtxLocal.TINDEX || blCtxLocal.TINDEX.length < 3) {
    blCtxLocal.TINDEX = new Float64Array(3);
  }
}

xywake(ctxPanel);
qwcalc(ctxPanel);
const { qinv, qinvA } = qiset(ctxPanel, alphaRad);
const st = stfind(ctxPanel, n);
blCtx.IST = st.ist + 1;
blCtx.SST = st.sst;

iblpan(blCtx, n, nw);
xicalc(blCtx, ctxPanel);
iblsys(blCtx);
ensureViscousArrays(blCtx);
uicalc(blCtx, qinv, qinvA);

for (let is = 1; is <= 2; is += 1) {
  for (let ibl = 1; ibl <= blCtx.NBL[is]; ibl += 1) {
    blCtx.UEDG[ibl][is] = blCtx.UINV[ibl][is];
  }
}

qdcalc(ctxPanel);
blCtx.DIJ = ctxPanel.DIJ;

qvfue(blCtx, ctxPanel.QVIS);
gamqv(ctxPanel, ctxPanel.QVIS, qinvA);
stmove(ctxPanel, blCtx, qinv, qinvA);

setbl(blCtx);

function flattenTensor3(mat, d1, d2, d3) {
  const data = new Array(d1 * d2 * d3);
  let idx = 0;
  for (let k = 1; k <= d1; k += 1) {
    for (let j = 1; j <= d2; j += 1) {
      for (let i = 1; i <= d3; i += 1) {
        data[idx] = mat[k][j][i] ?? 0.0;
        idx += 1;
      }
    }
  }
  return data;
}

const nsys = blCtx.NSYS;
const results = {
  nsys,
  reybl: blCtx.REYBL,
  rstbl: blCtx.RSTBL,
  hstinv: blCtx.HSTINV,
  va: { dims: [3, 2, nsys], data: flattenTensor3(blCtx.VA, 3, 2, nsys) },
  vb: { dims: [3, 2, nsys], data: flattenTensor3(blCtx.VB, 3, 2, nsys) },
  vdel: { dims: [3, 2, nsys], data: flattenTensor3(blCtx.VDEL, 3, 2, nsys) },
  vm: { dims: [3, nsys, nsys], data: flattenTensor3(blCtx.VM, 3, nsys, nsys) },
};

process.stdout.write(JSON.stringify({ results }));
