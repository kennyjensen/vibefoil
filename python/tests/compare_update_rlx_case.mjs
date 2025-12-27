import fs from 'fs';

import { ggcalc, qdcalc, qwcalc, qiset, stfind, iblpan, xicalc, uicalc, qvfue, gamqv, stmove, xywake } from '../../js/xpanel.js';
import { buildBlContext } from '../../js/xoper.js';
import { setbl, iblsys } from '../../js/xbl.js';
import { blsolv } from '../../js/xsolve.js';
import { createMatrix, createMatrix1, createTensor3 } from '../../js/arrays.js';

const noop = () => {};
console.log = noop;
console.warn = noop;

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const ctxPayload = payload.ctx ?? payload;

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
  if (!blCtxLocal.VA || blCtxLocal.VA.length <= 3 || blCtxLocal.VA[1]?.length <= 2) {
    blCtxLocal.VA = createTensor3(3, 2, nsys);
    blCtxLocal.VB = createTensor3(3, 2, nsys);
    blCtxLocal.VDEL = createTensor3(3, 2, nsys);
  }
  if (!blCtxLocal.VM || blCtxLocal.VM.length <= 3) {
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
blCtx.LBLINI = false;
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
blsolv(blCtx);

function computeUpdateMetrics(ctxLocal) {
  const unew = createMatrix1(ctxLocal.IVX, 2);
  const uAc = createMatrix1(ctxLocal.IVX, 2);
  const qnew = new Float64Array(ctxLocal.IQX + 1);
  const qAc = new Float64Array(ctxLocal.IQX + 1);

  let dclmin = -0.5;
  const dclmax = 0.5;
  if (ctxLocal.MATYP !== 1) dclmin = Math.max(-0.5, -0.9 * ctxLocal.CL);

  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= ctxLocal.NBL[is]; ibl += 1) {
      const i = ctxLocal.IPAN[ibl][is];
      let dui = 0.0;
      let duiAc = 0.0;
      for (let js = 1; js <= 2; js += 1) {
        for (let jbl = 2; jbl <= ctxLocal.NBL[js]; jbl += 1) {
          const j = ctxLocal.IPAN[jbl][js];
          const jv = ctxLocal.ISYS[jbl][js];
          const ueM = -ctxLocal.VTI[ibl][is] * ctxLocal.VTI[jbl][js] * ctxLocal.DIJ[i][j];
          dui += ueM * (ctxLocal.MASS[jbl][js] + ctxLocal.VDEL[3][1][jv]);
          duiAc += ueM * (0.0 - ctxLocal.VDEL[3][2][jv]);
        }
      }
      const uinvAc = ctxLocal.LALFA ? 0.0 : ctxLocal.UINV_A[ibl][is];
      unew[ibl][is] = ctxLocal.UINV[ibl][is] + dui;
      uAc[ibl][is] = uinvAc + duiAc;
    }
  }

  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= ctxLocal.IBLTE[is]; ibl += 1) {
      const i = ctxLocal.IPAN[ibl][is];
      qnew[i] = ctxLocal.VTI[ibl][is] * unew[ibl][is];
      qAc[i] = ctxLocal.VTI[ibl][is] * uAc[ibl][is];
    }
  }

  const sa = Math.sin(ctxLocal.ALFA);
  const ca = Math.cos(ctxLocal.ALFA);
  const beta = Math.sqrt(1.0 - ctxLocal.MINF ** 2);
  const betaMsq = -0.5 / beta;
  const bfac = 0.5 * ctxLocal.MINF ** 2 / (1.0 + beta);
  const bfacMsq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * betaMsq;

  let clnew = 0.0;
  let clA = 0.0;
  let clMs = 0.0;
  let clAc = 0.0;

  let i = 1;
  let cginc = 1.0 - (qnew[i] / ctxLocal.QINF) ** 2;
  let cpg1 = cginc / (beta + bfac * cginc);
  let cpg1Ms = -cpg1 / (beta + bfac * cginc) * (betaMsq + bfacMsq * cginc);
  let cpiQ = -2.0 * qnew[i] / ctxLocal.QINF ** 2;
  let cpcCpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
  let cpg1Ac = cpcCpi * cpiQ * qAc[i];

  for (i = 1; i <= ctxLocal.N; i += 1) {
    let ip = i + 1;
    if (i === ctxLocal.N) ip = 1;
    cginc = 1.0 - (qnew[ip] / ctxLocal.QINF) ** 2;
    const cpg2 = cginc / (beta + bfac * cginc);
    const cpg2Ms = -cpg2 / (beta + bfac * cginc) * (betaMsq + bfacMsq * cginc);
    cpiQ = -2.0 * qnew[ip] / ctxLocal.QINF ** 2;
    cpcCpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
    const cpg2Ac = cpcCpi * cpiQ * qAc[ip];
    const dx = (ctxLocal.X[ip] - ctxLocal.X[i]) * ca + (ctxLocal.Y[ip] - ctxLocal.Y[i]) * sa;
    const dxA = -(ctxLocal.X[ip] - ctxLocal.X[i]) * sa + (ctxLocal.Y[ip] - ctxLocal.Y[i]) * ca;
    const ag = 0.5 * (cpg2 + cpg1);
    const agMs = 0.5 * (cpg2Ms + cpg1Ms);
    const agAc = 0.5 * (cpg2Ac + cpg1Ac);
    clnew += dx * ag;
    clA += dxA * ag;
    clMs += dx * agMs;
    clAc += dx * agAc;
    cpg1 = cpg2;
    cpg1Ms = cpg2Ms;
    cpg1Ac = cpg2Ac;
  }

  let rlx = 1.0;
  let dac = 0.0;
  if (ctxLocal.LALFA) {
    dac = (clnew - ctxLocal.CL)
      / (1.0 - clAc - clMs * 2.0 * ctxLocal.MINF * ctxLocal.MINF_CL);
    if (rlx * dac > dclmax) rlx = dclmax / dac;
    if (rlx * dac < dclmin) rlx = dclmin / dac;
  } else {
    const dalmax = 0.5 * ctxLocal.DTOR;
    const dalmin = -0.5 * ctxLocal.DTOR;
    dac = (clnew - ctxLocal.CLSPEC) / (0.0 - clAc - clA);
    if (rlx * dac > dalmax) rlx = dalmax / dac;
    if (rlx * dac < dalmin) rlx = dalmin / dac;
  }
  const rlxDac = rlx;

  let rmsbl = 0.0;
  let rmxbl = 0.0;
  let vmxbl = ' ';
  let imxbl = 0;
  let ismxbl = 0;
  let dnMax = [0.0, 0.0, 0.0, 0.0];

  const dhi = 1.5;
  const dlo = -0.5;

  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= ctxLocal.NBL[is]; ibl += 1) {
      const iv = ctxLocal.ISYS[ibl][is];
      const dctau = ctxLocal.VDEL[1][1][iv] - dac * ctxLocal.VDEL[1][2][iv];
      const dthet = ctxLocal.VDEL[2][1][iv] - dac * ctxLocal.VDEL[2][2][iv];
      const dmass = ctxLocal.VDEL[3][1][iv] - dac * ctxLocal.VDEL[3][2][iv];
      const duedg = unew[ibl][is] + dac * uAc[ibl][is] - ctxLocal.UEDG[ibl][is];
      const ddstr = (dmass - ctxLocal.DSTR[ibl][is] * duedg) / ctxLocal.UEDG[ibl][is];
      let dn1 = 0.0;
      if (ibl < ctxLocal.ITRAN[is]) dn1 = dctau / 10.0;
      if (ibl >= ctxLocal.ITRAN[is]) dn1 = dctau / ctxLocal.CTAU[ibl][is];
      const dn2 = dthet / ctxLocal.THET[ibl][is];
      const dn3 = ddstr / ctxLocal.DSTR[ibl][is];
      const dn4 = Math.abs(duedg) / 0.25;
      rmsbl += dn1 ** 2 + dn2 ** 2 + dn3 ** 2 + dn4 ** 2;

      let rdn1 = rlx * dn1;
      if (Math.abs(dn1) > Math.abs(rmxbl)) {
        rmxbl = dn1;
        vmxbl = ibl < ctxLocal.ITRAN[is] ? 'n' : 'C';
        imxbl = ibl;
        ismxbl = is;
        dnMax = [dn1, dn2, dn3, dn4];
      }
      if (rdn1 > dhi) rlx = dhi / dn1;
      if (rdn1 < dlo) rlx = dlo / dn1;

      let rdn2 = rlx * dn2;
      if (Math.abs(dn2) > Math.abs(rmxbl)) {
        rmxbl = dn2;
        vmxbl = 'T';
        imxbl = ibl;
        ismxbl = is;
        dnMax = [dn1, dn2, dn3, dn4];
      }
      if (rdn2 > dhi) rlx = dhi / dn2;
      if (rdn2 < dlo) rlx = dlo / dn2;

      let rdn3 = rlx * dn3;
      if (Math.abs(dn3) > Math.abs(rmxbl)) {
        rmxbl = dn3;
        vmxbl = 'D';
        imxbl = ibl;
        ismxbl = is;
        dnMax = [dn1, dn2, dn3, dn4];
      }
      if (rdn3 > dhi) rlx = dhi / dn3;
      if (rdn3 < dlo) rlx = dlo / dn3;

      let rdn4 = rlx * dn4;
      if (Math.abs(dn4) > Math.abs(rmxbl)) {
        rmxbl = duedg;
        vmxbl = 'U';
        imxbl = ibl;
        ismxbl = is;
        dnMax = [dn1, dn2, dn3, dn4];
      }
      if (rdn4 > dhi) rlx = dhi / dn4;
      if (rdn4 < dlo) rlx = dlo / dn4;
    }
  }

  rmsbl = Math.sqrt(rmsbl / (4.0 * (ctxLocal.NBL[1] + ctxLocal.NBL[2])));

  return {
    clnew,
    clAc,
    clMs,
    dac,
    rlxDac,
    rlxFinal: rlx,
    rmsbl,
    rmxbl,
    vmxbl,
    imxbl,
    ismxbl,
    dnMax,
  };
}

const results = computeUpdateMetrics(blCtx);
process.stdout.write(JSON.stringify({ results }));
