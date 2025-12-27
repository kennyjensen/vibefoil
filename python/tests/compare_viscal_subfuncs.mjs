import fs from 'fs';

import '../../js/naca.js';
import { pangen, tecalc, clcalc, cdcalc, cpcalc } from '../../js/xfoil.js';
import { scalc, segspl } from '../../js/spline.js';
import { apcalc, ggcalc, ncalc, qdcalc, qwcalc, qiset, stfind, iblpan, xicalc, uicalc, qvfue, gamqv, stmove, xywake } from '../../js/xpanel.js';
import { buildBlContext } from '../../js/xoper.js';
import { setbl, update, iblsys } from '../../js/xbl.js';
import { blsolv } from '../../js/xsolve.js';
import { createMatrix, createMatrix1, createTensor3 } from '../../js/arrays.js';

const noop = () => {};
console.log = noop;
console.warn = noop;

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const { ctx: ctxPayload } = payload;

let ctxPanel = null;
let n = 0;
let nw = 0;
let total = 0;
let alphaRad = 0.0;
let minf = 0.0;
let reinf = 0.0;
let ncrit = 9.0;
let ides = 0;
let params = null;
let waklen = 0.0;
let nside = 0;

if (ctxPayload) {
  n = ctxPayload.N;
  nw = ctxPayload.NW;
  total = n + nw;
  alphaRad = ctxPayload.ALFA;
  minf = ctxPayload.MINF;
  reinf = ctxPayload.REINF;
  ncrit = ctxPayload.NCRIT ?? 9.0;

  ctxPanel = {
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
} else {
  ({ ides, params, alphaRad, reinf, minf, waklen, ncrit, nside } = payload);
  const { naca4, naca5 } = globalThis.Naca;
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
  nw = Math.floor(n / 12) + 10 * Math.floor(waklen);
  total = n + nw;

  ctxPanel = {
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
    MINF: minf,
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
  for (let i = 0; i < n; i += 1) {
    ctxPanel.GAM[i] = i < Math.floor(n / 2) ? 1.0 : -1.0;
  }
}

const blCtx = buildBlContext(n, ctxPanel, ncrit);
blCtx.MINF = minf;
blCtx.REINF = reinf;
blCtx.MINF1 = minf;
blCtx.REINF1 = reinf;
blCtx.HVRAT = 0.0;
blCtx.QINF = 1.0;
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

function sampleIndices(len) {
  if (len <= 0) return [1];
  const mid = Math.max(1, Math.floor((len + 1) / 2));
  return [1, mid, len];
}

function metrics1d(arr, len) {
  let sum = 0.0;
  let sumsq = 0.0;
  let maxabs = 0.0;
  for (let i = 1; i <= len; i += 1) {
    const val = arr[i] ?? 0.0;
    sum += val;
    sumsq += val * val;
    maxabs = Math.max(maxabs, Math.abs(val));
  }
  const samples = sampleIndices(len).map((i) => arr[i] ?? 0.0);
  return { sum, sumsq, maxabs, samples };
}

function metrics1dPanel(arr, len) {
  let sum = 0.0;
  let sumsq = 0.0;
  let maxabs = 0.0;
  for (let i = 1; i <= len; i += 1) {
    const val = arr[i - 1] ?? 0.0;
    sum += val;
    sumsq += val * val;
    maxabs = Math.max(maxabs, Math.abs(val));
  }
  const samples = sampleIndices(len).map((i) => arr[i - 1] ?? 0.0);
  return { sum, sumsq, maxabs, samples };
}

function metrics2d(mat, nbl1, nbl2) {
  let sum = 0.0;
  let sumsq = 0.0;
  let maxabs = 0.0;
  for (let is = 1; is <= 2; is += 1) {
    const nbl = is === 1 ? nbl1 : nbl2;
    for (let ibl = 1; ibl <= nbl; ibl += 1) {
      const val = mat[ibl][is] ?? 0.0;
      sum += val;
      sumsq += val * val;
      maxabs = Math.max(maxabs, Math.abs(val));
    }
  }
  const s1 = Math.min(nbl1, Math.max(1, Math.floor((nbl1 + 1) / 2)));
  const s2 = Math.min(nbl2, Math.max(1, Math.floor((nbl2 + 1) / 2)));
  const samples = [
    mat[Math.min(2, nbl1)][1] ?? 0.0,
    mat[nbl1][1] ?? 0.0,
    mat[Math.min(2, nbl2)][2] ?? 0.0,
    mat[nbl2][2] ?? 0.0,
    mat[s1][1] ?? 0.0,
    mat[s2][2] ?? 0.0,
  ];
  return { sum, sumsq, maxabs, samples };
}

function metrics3d(mat, nsys) {
  let sum = 0.0;
  let sumsq = 0.0;
  let maxabs = 0.0;
  for (let k = 1; k <= 3; k += 1) {
    for (let j = 1; j <= 2; j += 1) {
      for (let iv = 1; iv <= nsys; iv += 1) {
        const val = mat[k][j][iv] ?? 0.0;
        sum += val;
        sumsq += val * val;
        maxabs = Math.max(maxabs, Math.abs(val));
      }
    }
  }
  const mid = Math.max(1, Math.floor((nsys + 1) / 2));
  const samples = [
    mat[1][1][1] ?? 0.0,
    mat[2][1][mid] ?? 0.0,
    mat[3][2][nsys] ?? 0.0,
  ];
  return { sum, sumsq, maxabs, samples };
}

function countNonFinite3d(mat, nsys) {
  let count = 0;
  for (let k = 1; k <= 3; k += 1) {
    for (let j = 1; j <= 2; j += 1) {
      for (let iv = 1; iv <= nsys; iv += 1) {
        const val = mat[k][j][iv];
        if (!Number.isFinite(val)) count += 1;
      }
    }
  }
  return count;
}

function countNonFinite2d(mat, nbl1, nbl2) {
  let count = 0;
  for (let is = 1; is <= 2; is += 1) {
    const nbl = is === 1 ? nbl1 : nbl2;
    for (let ibl = 1; ibl <= nbl; ibl += 1) {
      const val = mat[ibl][is];
      if (!Number.isFinite(val)) count += 1;
    }
  }
  return count;
}

function countNonPositive2d(mat, nbl1, nbl2) {
  let count = 0;
  for (let is = 1; is <= 2; is += 1) {
    const nbl = is === 1 ? nbl1 : nbl2;
    for (let ibl = 1; ibl <= nbl; ibl += 1) {
      const val = mat[ibl][is];
      if (!Number.isFinite(val) || val <= 0.0) count += 1;
    }
  }
  return count;
}

function min2d(mat, nbl1, nbl2) {
  let minVal = Infinity;
  for (let is = 1; is <= 2; is += 1) {
    const nbl = is === 1 ? nbl1 : nbl2;
    for (let ibl = 2; ibl <= nbl; ibl += 1) {
      const val = mat[ibl][is];
      if (Number.isFinite(val)) {
        minVal = Math.min(minVal, val);
      }
    }
  }
  return Number.isFinite(minVal) ? minVal : null;
}

function countNonFinite1d(arr, n) {
  let count = 0;
  for (let i = 1; i <= n; i += 1) {
    if (!Number.isFinite(arr[i])) count += 1;
  }
  return count;
}

function countNonPositiveTurb(mat, itran, nbl1, nbl2) {
  let count = 0;
  for (let is = 1; is <= 2; is += 1) {
    const nbl = is === 1 ? nbl1 : nbl2;
    for (let ibl = Math.max(2, itran[is]); ibl <= nbl; ibl += 1) {
      const val = mat[ibl][is];
      if (!Number.isFinite(val) || val <= 0.0) count += 1;
    }
  }
  return count;
}

function findFirstNonFinite2d(mat, nbl1, nbl2) {
  for (let is = 1; is <= 2; is += 1) {
    const nbl = is === 1 ? nbl1 : nbl2;
    for (let ibl = 1; ibl <= nbl; ibl += 1) {
      const val = mat[ibl][is];
      if (!Number.isFinite(val)) return { is, ibl, val };
    }
  }
  return null;
}

const results = {};

xywake(ctxPanel);
results.xywake = {
  X: metrics1dPanel(ctxPanel.X, total),
  Y: metrics1dPanel(ctxPanel.Y, total),
  S: metrics1dPanel(ctxPanel.S, total),
  NX: metrics1dPanel(ctxPanel.NX, total),
  NY: metrics1dPanel(ctxPanel.NY, total),
  APANEL: metrics1dPanel(ctxPanel.APANEL, total),
};

qwcalc(ctxPanel);
const qinvu1 = new Float64Array(total + 1);
const qinvu2 = new Float64Array(total + 1);
for (let i = 1; i <= total; i += 1) {
  qinvu1[i] = ctxPanel.QINVU[i - 1][0];
  qinvu2[i] = ctxPanel.QINVU[i - 1][1];
}
results.qwcalc = {
  QINVU1: metrics1d(qinvu1, total),
  QINVU2: metrics1d(qinvu2, total),
};

const { qinv, qinvA } = qiset(ctxPanel, alphaRad);
results.qiset = {
  QINV: metrics1d(qinv, total),
  QINV_A: metrics1d(qinvA, total),
};

const st = stfind(ctxPanel, n);
blCtx.IST = st.ist + 1;
blCtx.SST = st.sst;
results.stfind = { IST: blCtx.IST, SST: blCtx.SST };

iblpan(blCtx, n, nw);
results.iblpan = {
  NBL: Array.from(blCtx.NBL),
  IBLTE: Array.from(blCtx.IBLTE),
  IPAN: metrics2d(blCtx.IPAN, blCtx.NBL[1], blCtx.NBL[2]),
  VTI: metrics2d(blCtx.VTI, blCtx.NBL[1], blCtx.NBL[2]),
};

xicalc(blCtx, ctxPanel);
results.xicalc = {
  XSSI: metrics2d(blCtx.XSSI, blCtx.NBL[1], blCtx.NBL[2]),
};

iblsys(blCtx);
ensureViscousArrays(blCtx);
results.iblsys = {
  NSYS: blCtx.NSYS,
  ISYS: metrics2d(blCtx.ISYS, blCtx.NBL[1], blCtx.NBL[2]),
};

uicalc(blCtx, qinv, qinvA);
results.uicalc = {
  UINV: metrics2d(blCtx.UINV, blCtx.NBL[1], blCtx.NBL[2]),
  UINV_A: metrics2d(blCtx.UINV_A, blCtx.NBL[1], blCtx.NBL[2]),
};

for (let is = 1; is <= 2; is += 1) {
  for (let ibl = 1; ibl <= blCtx.NBL[is]; ibl += 1) {
    blCtx.UEDG[ibl][is] = blCtx.UINV[ibl][is];
  }
}

qdcalc(ctxPanel);
blCtx.DIJ = ctxPanel.DIJ;

qvfue(blCtx, ctxPanel.QVIS);
results.qvfue = {
  QVIS: metrics1d(ctxPanel.QVIS, total),
};

gamqv(ctxPanel, ctxPanel.QVIS, qinvA);
results.gamqv = {
  GAM: metrics1dPanel(ctxPanel.GAM, n),
  GAM_A: metrics1dPanel(ctxPanel.GAM_A, n),
};

stmove(ctxPanel, blCtx, qinv, qinvA);
results.stmove = {
  IST: blCtx.IST,
  XSSI: metrics2d(blCtx.XSSI, blCtx.NBL[1], blCtx.NBL[2]),
};

setbl(blCtx);
results.setbl = {
  REYBL: blCtx.REYBL,
  RSTBL: blCtx.RSTBL,
  HSTINV: blCtx.HSTINV,
  VA: metrics3d(blCtx.VA, blCtx.NSYS),
  VB: metrics3d(blCtx.VB, blCtx.NSYS),
  VDEL: metrics3d(blCtx.VDEL, blCtx.NSYS),
  VM: metrics3d(blCtx.VM, blCtx.NSYS),
};

results.preUpdate = {
  CTAU: countNonFinite2d(blCtx.CTAU, blCtx.NBL[1], blCtx.NBL[2]),
  THET: countNonFinite2d(blCtx.THET, blCtx.NBL[1], blCtx.NBL[2]),
  DSTR: countNonFinite2d(blCtx.DSTR, blCtx.NBL[1], blCtx.NBL[2]),
  UEDG: countNonFinite2d(blCtx.UEDG, blCtx.NBL[1], blCtx.NBL[2]),
  THET_NONPOS: countNonPositive2d(blCtx.THET, blCtx.NBL[1], blCtx.NBL[2]),
  DSTR_NONPOS: countNonPositive2d(blCtx.DSTR, blCtx.NBL[1], blCtx.NBL[2]),
  UEDG_NONPOS: countNonPositive2d(blCtx.UEDG, blCtx.NBL[1], blCtx.NBL[2]),
  CTAU_TURB_NONPOS: countNonPositiveTurb(blCtx.CTAU, blCtx.ITRAN, blCtx.NBL[1], blCtx.NBL[2]),
  THET_MIN: min2d(blCtx.THET, blCtx.NBL[1], blCtx.NBL[2]),
  DSTR_MIN: min2d(blCtx.DSTR, blCtx.NBL[1], blCtx.NBL[2]),
  UEDG_MIN: min2d(blCtx.UEDG, blCtx.NBL[1], blCtx.NBL[2]),
  X_NONFIN: countNonFinite1d(blCtx.X, blCtx.N),
  Y_NONFIN: countNonFinite1d(blCtx.Y, blCtx.N),
};

const blsolvOk = blsolv(blCtx);
results.blsolv = {
  ok: blsolvOk,
  VDEL: metrics3d(blCtx.VDEL, blCtx.NSYS),
  VDEL_NAN: countNonFinite3d(blCtx.VDEL, blCtx.NSYS),
};

update(blCtx);
results.update = {
  CTAU: metrics2d(blCtx.CTAU, blCtx.NBL[1], blCtx.NBL[2]),
  THET: metrics2d(blCtx.THET, blCtx.NBL[1], blCtx.NBL[2]),
  DSTR: metrics2d(blCtx.DSTR, blCtx.NBL[1], blCtx.NBL[2]),
  UEDG: metrics2d(blCtx.UEDG, blCtx.NBL[1], blCtx.NBL[2]),
  CL: blCtx.CL,
  ALFA: blCtx.ALFA,
  RMSBL: blCtx.RMSBL,
  RMXBL: blCtx.RMXBL,
  RLX: blCtx.RLX,
  DAC: blCtx.DAC,
};

results.updateNonFinite = {
  CTAU: findFirstNonFinite2d(blCtx.CTAU, blCtx.NBL[1], blCtx.NBL[2]),
  THET: findFirstNonFinite2d(blCtx.THET, blCtx.NBL[1], blCtx.NBL[2]),
  DSTR: findFirstNonFinite2d(blCtx.DSTR, blCtx.NBL[1], blCtx.NBL[2]),
  UEDG: findFirstNonFinite2d(blCtx.UEDG, blCtx.NBL[1], blCtx.NBL[2]),
};

const cpi = cpcalc(qinv, ctxPanel.QINF, minf);
const clres = clcalc(n, ctxPanel.X, ctxPanel.Y, ctxPanel.GAM, ctxPanel.GAM_A, alphaRad, minf, ctxPanel.QINF, ctxPanel.XCMREF ?? 0.25, ctxPanel.YCMREF ?? 0.0);
const cdres = cdcalc(ctxPanel, blCtx, alphaRad, ctxPanel.QINF);
results.cpcalc = {
  CPI: metrics1d(cpi, total),
};
results.clcalc = {
  CL: clres.cl,
  CM: clres.cm,
  CDP: clres.cdp,
};
results.cdcalc = {
  CD: cdres.cd,
  CDF: cdres.cdf,
};

process.stdout.write(JSON.stringify({ results }));
