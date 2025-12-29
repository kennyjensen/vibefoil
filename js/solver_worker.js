import './naca.js';
import { scalc, segspl } from './spline.js';
import {
  apcalc,
  ncalc,
  psilin,
  ggcalc,
  stfind,
} from './xpanel.js';
import { computeCoefficients, cpcalc, tecalc, pangen } from './xfoil.js';
import { buildBlContext, computeQvisFromUedg, specal, viscal } from './xoper.js';
import { createMatrix } from './arrays.js';
import { flap as applyFlap } from './xgdes.js';

const nside = 123;
const xx = new Float64Array(nside);
const yt = new Float64Array(nside);
const yc = new Float64Array(nside);
const xb = new Float64Array(2 * nside);
const yb = new Float64Array(2 * nside);
const xbBuffer = new Float64Array(2 * nside);
const ybBuffer = new Float64Array(2 * nside);

let panelCtx = null;
let panelX = null;
let panelY = null;
let panelXP = null;
let panelYP = null;
const panelCache = { ctx: null, key: null };
const blCache = { ctx: null, key: null };

function rotatePoint(x, y, angle, ox, oy) {
  const dx = x - ox;
  const dy = y - oy;
  const ca = Math.cos(angle);
  const sa = Math.sin(angle);
  return {
    x: ox + dx * ca - dy * sa,
    y: oy + dx * sa + dy * ca,
  };
}

function computeBounds(nb, angle, width, height) {
  let xmin = xb[0];
  let xmax = xb[0];
  let ymin = yb[0];
  let ymax = yb[0];

  for (let i = 1; i < nb; i += 1) {
    xmin = Math.min(xmin, xb[i]);
    xmax = Math.max(xmax, xb[i]);
    ymin = Math.min(ymin, yb[i]);
    ymax = Math.max(ymax, yb[i]);
  }

  const chord = xmax - xmin || 1.0;
  const ox = 0.5 * (xmin + xmax);
  const oy = 0.5 * (ymin + ymax);

  let rxmin = xmin;
  let rxmax = xmax;
  let rymin = ymin;
  let rymax = ymax;

  if (angle !== 0.0) {
    rxmin = Infinity;
    rxmax = -Infinity;
    rymin = Infinity;
    rymax = -Infinity;
    for (let i = 0; i < nb; i += 1) {
      const rp = rotatePoint(xb[i], yb[i], angle, ox, oy);
      rxmin = Math.min(rxmin, rp.x);
      rxmax = Math.max(rxmax, rp.x);
      rymin = Math.min(rymin, rp.y);
      rymax = Math.max(rymax, rp.y);
    }
  }

  const marginX = 0.05 * chord;
  const marginY = 0.2 * chord;
  rxmin -= marginX;
  rxmax += marginX;
  rymin -= marginY;
  rymax += marginY;

  const padding = 32;
  const spanX = rxmax - rxmin || 1.0;
  const spanY = rymax - rymin || 1.0;
  const scale = Math.min((width - padding * 2) / spanX, (height - padding * 2) / spanY);

  return {
    xmin: rxmin,
    xmax: rxmax,
    ymin: rymin,
    ymax: rymax,
    width,
    height,
    padding,
    scale,
    angle,
    ox,
    oy,
  };
}

function pointInPolygon(x, y, px, py, n) {
  let inside = false;
  for (let i = 0, j = n - 1; i < n; j = i, i += 1) {
    const xi = px[i];
    const yi = py[i];
    const xj = px[j];
    const yj = py[j];

    const intersect = ((yi > y) !== (yj > y))
      && (x < (xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi);
    if (intersect) {
      inside = !inside;
    }
  }
  return inside;
}

function buildPanelContext(nb, alphaRad, opts = {}) {
  const { reusePanel = false, geometryKey = '', reuseState = null } = opts;
  const waklen = 1.0;
  const nw = Math.floor(nb / 12) + 10 * Math.floor(waklen);
  const total = nb + nw;
  if (reusePanel && panelCache.ctx && panelCache.key === geometryKey
    && panelCache.ctx.N === nb && panelCache.ctx.NW === nw) {
    panelCtx = panelCache.ctx;
    panelCtx.ALFA = alphaRad;
    return panelCtx;
  }

  const resetViscous = !reusePanel || panelCache.key !== geometryKey;
  if (!panelCtx || panelCtx.N !== nb || panelCtx.NW !== nw) {
    panelX = new Float64Array(total);
    panelY = new Float64Array(total);
    panelXP = new Float64Array(total);
    panelYP = new Float64Array(total);
    const n1 = nb + 1;
    panelCtx = {
      N: nb,
      NW: nw,
      WAKLEN: waklen,
      X: panelX,
      Y: panelY,
      XP: panelXP,
      YP: panelYP,
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
      ALFA: 0.0,
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
  }

  for (let i = 0; i < nb; i += 1) {
    panelX[i] = xb[i];
    panelY[i] = yb[i];
  }

  scalc(panelX, panelY, panelCtx.S, nb);
  segspl(panelX, panelXP, panelCtx.S, nb);
  segspl(panelY, panelYP, panelCtx.S, nb);
  ncalc(panelX, panelY, panelCtx.S, nb, panelCtx.NX, panelCtx.NY);
  panelCtx.XTE = 0.5 * (panelX[0] + panelX[nb - 1]);
  panelCtx.YTE = 0.5 * (panelY[0] + panelY[nb - 1]);
  tecalc(panelCtx);
  apcalc(panelCtx);
  panelCtx.ALFA = alphaRad;
  ggcalc(panelCtx);

  if (resetViscous && panelCtx.QVIS) {
    panelCtx.QVIS.fill(0.0);
    panelCtx.LWAKE = false;
    panelCtx.LWDIJ = false;
    panelCtx.LADIJ = false;
  }

  panelCache.ctx = panelCtx;
  panelCache.key = geometryKey;

  if (reuseState && reuseState.panel) {
    const restored = applyPanelSnapshot(panelCtx, reuseState.panel);
    if (restored) {
      panelCache.ctx = panelCtx;
      panelCache.key = geometryKey;
    }
  }
  return panelCtx;
}

function cloneFloat64Array(arr) {
  return arr ? new Float64Array(arr) : null;
}

function cloneInt32Array(arr) {
  return arr ? new Int32Array(arr) : null;
}

function cloneArrayOfFloat64Arrays(arr) {
  if (!arr) return null;
  const out = new Array(arr.length);
  for (let i = 0; i < arr.length; i += 1) {
    const row = arr[i];
    out[i] = row instanceof Float64Array ? new Float64Array(row) : row;
  }
  return out;
}

function cloneMatrix1(mat) {
  return cloneArrayOfFloat64Arrays(mat);
}

function snapshotPanelCtx(ctx) {
  return {
    N: ctx.N,
    NW: ctx.NW,
    WAKLEN: ctx.WAKLEN,
    LWAKE: ctx.LWAKE,
    LWDIJ: ctx.LWDIJ,
    LADIJ: ctx.LADIJ,
    SHARP: ctx.SHARP,
    ANTE: ctx.ANTE,
    ASTE: ctx.ASTE,
    DSTE: ctx.DSTE,
    XTE: ctx.XTE,
    YTE: ctx.YTE,
    QINF: ctx.QINF,
    ALFA: ctx.ALFA,
    X: cloneFloat64Array(ctx.X),
    Y: cloneFloat64Array(ctx.Y),
    XP: cloneFloat64Array(ctx.XP),
    YP: cloneFloat64Array(ctx.YP),
    S: cloneFloat64Array(ctx.S),
    NX: cloneFloat64Array(ctx.NX),
    NY: cloneFloat64Array(ctx.NY),
    APANEL: cloneFloat64Array(ctx.APANEL),
    SIG: cloneFloat64Array(ctx.SIG),
    QF0: cloneFloat64Array(ctx.QF0),
    QF1: cloneFloat64Array(ctx.QF1),
    QF2: cloneFloat64Array(ctx.QF2),
    QF3: cloneFloat64Array(ctx.QF3),
    DZDG: cloneFloat64Array(ctx.DZDG),
    DZDN: cloneFloat64Array(ctx.DZDN),
    DQDG: cloneFloat64Array(ctx.DQDG),
    DZDM: cloneFloat64Array(ctx.DZDM),
    DQDM: cloneFloat64Array(ctx.DQDM),
    QINV: cloneFloat64Array(ctx.QINV),
    QINV_A: cloneFloat64Array(ctx.QINV_A),
    QVIS: cloneFloat64Array(ctx.QVIS),
    GAM: cloneFloat64Array(ctx.GAM),
    GAM_A: cloneFloat64Array(ctx.GAM_A),
    GAMU: cloneArrayOfFloat64Arrays(ctx.GAMU),
    QINVU: cloneArrayOfFloat64Arrays(ctx.QINVU),
    AIJ: cloneArrayOfFloat64Arrays(ctx.AIJ),
    BIJ: cloneArrayOfFloat64Arrays(ctx.BIJ),
    AIJPIV: cloneInt32Array(ctx.AIJPIV),
    DIJ: cloneArrayOfFloat64Arrays(ctx.DIJ),
  };
}

function applyPanelSnapshot(ctx, snap) {
  if (!snap || snap.N !== ctx.N || snap.NW !== ctx.NW) return false;
  ctx.WAKLEN = snap.WAKLEN ?? ctx.WAKLEN;
  ctx.LWAKE = !!snap.LWAKE;
  ctx.LWDIJ = !!snap.LWDIJ;
  ctx.LADIJ = !!snap.LADIJ;
  ctx.SHARP = !!snap.SHARP;
  ctx.ANTE = snap.ANTE ?? ctx.ANTE;
  ctx.ASTE = snap.ASTE ?? ctx.ASTE;
  ctx.DSTE = snap.DSTE ?? ctx.DSTE;
  ctx.XTE = snap.XTE ?? ctx.XTE;
  ctx.YTE = snap.YTE ?? ctx.YTE;
  ctx.QINF = snap.QINF ?? ctx.QINF;
  ctx.ALFA = snap.ALFA ?? ctx.ALFA;
  if (snap.X) ctx.X.set(snap.X);
  if (snap.Y) ctx.Y.set(snap.Y);
  if (snap.XP) ctx.XP.set(snap.XP);
  if (snap.YP) ctx.YP.set(snap.YP);
  if (snap.S) ctx.S.set(snap.S);
  if (snap.NX) ctx.NX.set(snap.NX);
  if (snap.NY) ctx.NY.set(snap.NY);
  if (snap.APANEL) ctx.APANEL.set(snap.APANEL);
  if (snap.SIG) ctx.SIG.set(snap.SIG);
  if (snap.QF0) ctx.QF0.set(snap.QF0);
  if (snap.QF1) ctx.QF1.set(snap.QF1);
  if (snap.QF2) ctx.QF2.set(snap.QF2);
  if (snap.QF3) ctx.QF3.set(snap.QF3);
  if (snap.DZDG) ctx.DZDG.set(snap.DZDG);
  if (snap.DZDN) ctx.DZDN.set(snap.DZDN);
  if (snap.DQDG) ctx.DQDG.set(snap.DQDG);
  if (snap.DZDM) ctx.DZDM.set(snap.DZDM);
  if (snap.DQDM) ctx.DQDM.set(snap.DQDM);
  if (snap.QINV) ctx.QINV.set(snap.QINV);
  if (snap.QINV_A) ctx.QINV_A.set(snap.QINV_A);
  if (snap.QVIS) ctx.QVIS.set(snap.QVIS);
  if (snap.GAM) ctx.GAM.set(snap.GAM);
  if (snap.GAM_A) ctx.GAM_A.set(snap.GAM_A);
  if (snap.GAMU) {
    for (let i = 0; i < snap.GAMU.length && i < ctx.GAMU.length; i += 1) {
      if (snap.GAMU[i]) ctx.GAMU[i].set(snap.GAMU[i]);
    }
  }
  if (snap.QINVU) {
    for (let i = 0; i < snap.QINVU.length && i < ctx.QINVU.length; i += 1) {
      if (snap.QINVU[i]) ctx.QINVU[i].set(snap.QINVU[i]);
    }
  }
  if (snap.AIJ) {
    for (let i = 0; i < snap.AIJ.length && i < ctx.AIJ.length; i += 1) {
      if (snap.AIJ[i]) ctx.AIJ[i].set(snap.AIJ[i]);
    }
  }
  if (snap.BIJ) {
    for (let i = 0; i < snap.BIJ.length && i < ctx.BIJ.length; i += 1) {
      if (snap.BIJ[i]) ctx.BIJ[i].set(snap.BIJ[i]);
    }
  }
  if (snap.AIJPIV) ctx.AIJPIV.set(snap.AIJPIV);
  if (snap.DIJ) {
    if (!ctx.DIJ || ctx.DIJ.length !== snap.DIJ.length) {
      ctx.DIJ = cloneArrayOfFloat64Arrays(snap.DIJ);
    } else {
      for (let i = 0; i < snap.DIJ.length; i += 1) {
        if (snap.DIJ[i]) ctx.DIJ[i].set(snap.DIJ[i]);
      }
    }
  } else {
    ctx.LWDIJ = false;
    ctx.LADIJ = false;
  }
  if ((!ctx.DIJ) && (ctx.LWDIJ || ctx.LADIJ)) {
    ctx.LWDIJ = false;
    ctx.LADIJ = false;
  }
  return true;
}

function snapshotBlCtx(ctx) {
  const maxNbl = Math.max(ctx.NBL[1], ctx.NBL[2]);
  return {
    N: ctx.N,
    NW: ctx.NW,
    LBLINI: !!ctx.LBLINI,
    LVCONV: !!ctx.LVCONV,
    LIPAN: !!ctx.LIPAN,
    NSYS: ctx.NSYS,
    IST: ctx.IST,
    SST: ctx.SST,
    SST_GO: ctx.SST_GO,
    SST_GP: ctx.SST_GP,
    SLE: ctx.SLE,
    XLE: ctx.XLE,
    YLE: ctx.YLE,
    XTE: ctx.XTE,
    YTE: ctx.YTE,
    AVISC: ctx.AVISC,
    MVISC: ctx.MVISC,
    MINF: ctx.MINF,
    MINF1: ctx.MINF1,
    REINF: ctx.REINF,
    REINF1: ctx.REINF1,
    ACRIT: cloneFloat64Array(ctx.ACRIT),
    XSTRIP: cloneFloat64Array(ctx.XSTRIP),
    NBL: cloneInt32Array(ctx.NBL),
    IBLTE: cloneInt32Array(ctx.IBLTE),
    ITRAN: cloneInt32Array(ctx.ITRAN),
    XSSITR: cloneFloat64Array(ctx.XSSITR),
    TFORCE: Array.from(ctx.TFORCE ?? []),
    WGAP: cloneFloat64Array(ctx.WGAP),
    IPAN: cloneMatrix1(ctx.IPAN),
    ISYS: cloneMatrix1(ctx.ISYS),
    VTI: cloneMatrix1(ctx.VTI),
    XSSI: cloneMatrix1(ctx.XSSI),
    UEDG: cloneMatrix1(ctx.UEDG),
    UINV: cloneMatrix1(ctx.UINV),
    UINV_A: cloneMatrix1(ctx.UINV_A),
    THET: cloneMatrix1(ctx.THET),
    DSTR: cloneMatrix1(ctx.DSTR),
    CTAU: cloneMatrix1(ctx.CTAU),
    MASS: cloneMatrix1(ctx.MASS),
    TAU: cloneMatrix1(ctx.TAU),
    DIS: cloneMatrix1(ctx.DIS),
    CTQ: cloneMatrix1(ctx.CTQ),
    DELT: cloneMatrix1(ctx.DELT),
    TSTR: cloneMatrix1(ctx.TSTR),
    maxNbl,
  };
}

function applyBlSnapshot(ctx, snap) {
  if (!snap || snap.N !== ctx.N || snap.NW !== ctx.NW) return false;
  ctx.LBLINI = !!snap.LBLINI;
  ctx.LVCONV = !!snap.LVCONV;
  ctx.LIPAN = !!snap.LIPAN;
  ctx.NSYS = snap.NSYS ?? ctx.NSYS;
  ctx.IST = snap.IST ?? ctx.IST;
  ctx.SST = snap.SST ?? ctx.SST;
  ctx.SST_GO = snap.SST_GO ?? ctx.SST_GO;
  ctx.SST_GP = snap.SST_GP ?? ctx.SST_GP;
  ctx.SLE = snap.SLE ?? ctx.SLE;
  ctx.XLE = snap.XLE ?? ctx.XLE;
  ctx.YLE = snap.YLE ?? ctx.YLE;
  ctx.XTE = snap.XTE ?? ctx.XTE;
  ctx.YTE = snap.YTE ?? ctx.YTE;
  ctx.AVISC = snap.AVISC ?? ctx.AVISC;
  ctx.MVISC = snap.MVISC ?? ctx.MVISC;
  ctx.MINF = snap.MINF ?? ctx.MINF;
  ctx.MINF1 = snap.MINF1 ?? ctx.MINF1;
  ctx.REINF = snap.REINF ?? ctx.REINF;
  ctx.REINF1 = snap.REINF1 ?? ctx.REINF1;
  if (snap.ACRIT) ctx.ACRIT.set(snap.ACRIT);
  if (snap.XSTRIP) ctx.XSTRIP.set(snap.XSTRIP);
  if (snap.NBL) ctx.NBL.set(snap.NBL);
  if (snap.IBLTE) ctx.IBLTE.set(snap.IBLTE);
  if (snap.ITRAN) ctx.ITRAN.set(snap.ITRAN);
  if (snap.XSSITR) ctx.XSSITR.set(snap.XSSITR);
  if (snap.TFORCE && ctx.TFORCE) {
    for (let i = 0; i < snap.TFORCE.length && i < ctx.TFORCE.length; i += 1) {
      ctx.TFORCE[i] = snap.TFORCE[i];
    }
  }
  if (snap.WGAP) ctx.WGAP.set(snap.WGAP);
  if (snap.IPAN) {
    for (let i = 0; i < snap.IPAN.length && i < ctx.IPAN.length; i += 1) {
      if (snap.IPAN[i]) ctx.IPAN[i].set(snap.IPAN[i]);
    }
  }
  if (snap.ISYS) {
    for (let i = 0; i < snap.ISYS.length && i < ctx.ISYS.length; i += 1) {
      if (snap.ISYS[i]) ctx.ISYS[i].set(snap.ISYS[i]);
    }
  }
  if (snap.VTI) {
    for (let i = 0; i < snap.VTI.length && i < ctx.VTI.length; i += 1) {
      if (snap.VTI[i]) ctx.VTI[i].set(snap.VTI[i]);
    }
  }
  if (snap.XSSI) {
    for (let i = 0; i < snap.XSSI.length && i < ctx.XSSI.length; i += 1) {
      if (snap.XSSI[i]) ctx.XSSI[i].set(snap.XSSI[i]);
    }
  }
  if (snap.UEDG) {
    for (let i = 0; i < snap.UEDG.length && i < ctx.UEDG.length; i += 1) {
      if (snap.UEDG[i]) ctx.UEDG[i].set(snap.UEDG[i]);
    }
  }
  if (snap.UINV) {
    for (let i = 0; i < snap.UINV.length && i < ctx.UINV.length; i += 1) {
      if (snap.UINV[i]) ctx.UINV[i].set(snap.UINV[i]);
    }
  }
  if (snap.UINV_A) {
    for (let i = 0; i < snap.UINV_A.length && i < ctx.UINV_A.length; i += 1) {
      if (snap.UINV_A[i]) ctx.UINV_A[i].set(snap.UINV_A[i]);
    }
  }
  const fields = ['THET', 'DSTR', 'CTAU', 'MASS', 'TAU', 'DIS', 'CTQ', 'DELT', 'TSTR'];
  for (let f = 0; f < fields.length; f += 1) {
    const key = fields[f];
    const snapMat = snap[key];
    const ctxMat = ctx[key];
    if (!snapMat || !ctxMat) continue;
    for (let i = 0; i < snapMat.length && i < ctxMat.length; i += 1) {
      if (snapMat[i]) ctxMat[i].set(snapMat[i]);
    }
  }
  return true;
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

function getChordPoints(nb) {
  let leIdx = 0;
  let teIdx = 0;
  for (let i = 1; i < nb; i += 1) {
    if (xb[i] < xb[leIdx]) leIdx = i;
    if (xb[i] > xb[teIdx]) teIdx = i;
  }
  return {
    le: { x: xb[leIdx], y: yb[leIdx] },
    te: { x: xb[teIdx], y: yb[teIdx] },
  };
}

function buildBoundaryLayerLines(blCtx, ctxPanel, userScale) {
  const scale = Math.max(userScale, 0.05);
  const lines = {
    upper: [],
    lower: [],
    wakeUpper: [],
    wakeLower: [],
  };
  if (!blCtx || !ctxPanel) return lines;

  for (let ibl = 2; ibl <= blCtx.IBLTE[1]; ibl += 1) {
    const i = blCtx.IPAN[ibl][1];
    const dstr = blCtx.DSTR[ibl][1];
    if (!Number.isFinite(dstr)) continue;
    const dVis = dstr * scale;
    const x = ctxPanel.X[i - 1] + ctxPanel.NX[i - 1] * dVis;
    const y = ctxPanel.Y[i - 1] + ctxPanel.NY[i - 1] * dVis;
    lines.upper.push({ x, y });
  }

  for (let ibl = 2; ibl <= blCtx.IBLTE[2]; ibl += 1) {
    const i = blCtx.IPAN[ibl][2];
    const dstr = blCtx.DSTR[ibl][2];
    if (!Number.isFinite(dstr)) continue;
    const dVis = dstr * scale;
    const x = ctxPanel.X[i - 1] + ctxPanel.NX[i - 1] * dVis;
    const y = ctxPanel.Y[i - 1] + ctxPanel.NY[i - 1] * dVis;
    lines.lower.push({ x, y });
  }

  const nw = blCtx.NW ?? 0;
  if (nw > 0) {
    const is = 2;
    const dstrTe = blCtx.DSTR[blCtx.IBLTE[is] + 1][is];
    let dsf1 = 0.5;
    let dsf2 = 0.5;
    if (dstrTe !== 0.0) {
      dsf1 = (blCtx.DSTR[blCtx.IBLTE[1]][1] + 0.5 * ctxPanel.ANTE) / dstrTe;
      dsf2 = (blCtx.DSTR[blCtx.IBLTE[2]][2] + 0.5 * ctxPanel.ANTE) / dstrTe;
    }

    for (let ibl = blCtx.IBLTE[is] + 1; ibl <= blCtx.NBL[is]; ibl += 1) {
      const i = blCtx.IPAN[ibl][is];
      const dstr = blCtx.DSTR[ibl][is];
      if (!Number.isFinite(dstr)) continue;
      const dVis = dstr * dsf1 * scale;
      const x = ctxPanel.X[i - 1] - ctxPanel.NX[i - 1] * dVis;
      const y = ctxPanel.Y[i - 1] - ctxPanel.NY[i - 1] * dVis;
      lines.wakeUpper.push({ x, y });
    }

    for (let ibl = blCtx.IBLTE[is] + 1; ibl <= blCtx.NBL[is]; ibl += 1) {
      const i = blCtx.IPAN[ibl][is];
      const dstr = blCtx.DSTR[ibl][is];
      if (!Number.isFinite(dstr)) continue;
      const dVis = dstr * dsf2 * scale;
      const x = ctxPanel.X[i - 1] + ctxPanel.NX[i - 1] * dVis;
      const y = ctxPanel.Y[i - 1] + ctxPanel.NY[i - 1] * dVis;
      lines.wakeLower.push({ x, y });
    }
  }

  return lines;
}

function buildStreamlineGrid(bounds, nb, ctxPanel) {
  if (!ctxPanel) return null;
  const gridX = 60;
  const gridY = 32;
  const xSpan = bounds.xmax - bounds.xmin;
  const ySpan = bounds.ymax - bounds.ymin;
  const dx = xSpan / (gridX - 1);
  const dy = ySpan / (gridY - 1);
  const grid = new Float64Array(gridX * gridY);

  let psiMin = Infinity;
  let psiMax = -Infinity;

  for (let j = 0; j < gridY; j += 1) {
    for (let i = 0; i < gridX; i += 1) {
      const x = bounds.xmin + i * dx;
      const y = bounds.ymin + j * dy;
      const idx = j * gridX + i;

      if (pointInPolygon(x, y, ctxPanel.X, ctxPanel.Y, nb)) {
        grid[idx] = NaN;
        continue;
      }

      const { psi } = psilin(nb, x, y, 0.0, 0.0, false, false, ctxPanel);
      grid[idx] = psi;
      if (psi < psiMin) psiMin = psi;
      if (psi > psiMax) psiMax = psi;
    }
  }

  if (!Number.isFinite(psiMin) || !Number.isFinite(psiMax)) {
    return null;
  }

  return {
    grid,
    gridX,
    gridY,
    psiMin,
    psiMax,
  };
}

function generateGeometry(settings) {
  const {
    mode,
    source,
    m,
    p,
    t,
    series5,
    t5,
    profile6,
    t6,
    cl6,
    custom,
    flap,
  } = settings;

  let nb = 0;
  let airfoilName = 'Airfoil';

  if (source === 'custom' || source === 'database') {
    if (!custom) {
      return { ok: false };
    }
    nb = Math.min(custom.nb, xb.length);
    const cx = Array.from(custom.x).slice(0, nb);
    const cy = Array.from(custom.y).slice(0, nb);
    const flapped = applyFlap(cx, cy, flap);
    nb = flapped.nb;
    xb.set(flapped.xb);
    yb.set(flapped.yb);
    airfoilName = custom.name || 'Custom Airfoil';
    return { ok: true, nb, airfoilName, hinge: flapped.hinge };
  }

  if (mode === '4') {
    const ides = m * 1000 + p * 100 + t;
    const res = globalThis.Naca.naca4(ides, xx, yt, yc, nside, xbBuffer, ybBuffer);
    const flapped = applyFlap(Array.from(xbBuffer).slice(0, res.nb), Array.from(ybBuffer).slice(0, res.nb), flap);
    const panelRes = pangen(Float64Array.from(flapped.xb), Float64Array.from(flapped.yb), flapped.nb);
    nb = panelRes.n;
    xb.set(panelRes.x);
    yb.set(panelRes.y);
    airfoilName = res.name;
    return { ok: true, nb, airfoilName, hinge: flapped.hinge };
  }

  if (mode === '5') {
    const n5 = parseInt(series5.charAt(0), 10);
    const n4 = parseInt(series5.charAt(1), 10);
    const n3 = parseInt(series5.charAt(2), 10);
    const ides = n5 * 10000 + n4 * 1000 + n3 * 100 + t5;
    const result = globalThis.Naca.naca5(ides, xx, yt, yc, nside, xbBuffer, ybBuffer);
    if (!result.ok) return { ok: false };
    const flapped = applyFlap(Array.from(xbBuffer).slice(0, result.nb), Array.from(ybBuffer).slice(0, result.nb), flap);
    const panelRes = pangen(Float64Array.from(flapped.xb), Float64Array.from(flapped.yb), flapped.nb);
    nb = panelRes.n;
    xb.set(panelRes.x);
    yb.set(panelRes.y);
    airfoilName = result.name;
    return { ok: true, nb, airfoilName, hinge: flapped.hinge };
  }

  const result = globalThis.Naca.naca6(
    {
      profile: profile6,
      toc: t6 / 100,
      camber: settings.camber6,
      cl: Number.isFinite(cl6) ? cl6 : 0.0,
      a: settings.fallbackA6,
    },
    xx,
    yt,
    yc,
    nside,
    xbBuffer,
    ybBuffer,
  );
  if (!result.ok) return { ok: false };
  const flapped = applyFlap(Array.from(xbBuffer).slice(0, result.nb), Array.from(ybBuffer).slice(0, result.nb), flap);
  const panelRes = pangen(Float64Array.from(flapped.xb), Float64Array.from(flapped.yb), flapped.nb);
  nb = panelRes.n;
  xb.set(panelRes.x);
  yb.set(panelRes.y);
  airfoilName = result.name || airfoilName;
  return { ok: true, nb, airfoilName, hinge: flapped.hinge };
}

function computeCase(settings) {
  const {
    alphaDeg,
    alphaRad,
    geometryKey,
    advancedMode = false,
    reusePanel,
    reuseSolution,
    reuseState,
    viscous,
    mach,
    reynolds,
    ncr,
    nIter,
    canvasWidth,
    canvasHeight,
  } = settings;

  const geom = generateGeometry(settings);
  if (!geom.ok) return { ok: false };

  const nb = geom.nb;
  const displayAngle = -alphaRad;
  const bounds = computeBounds(nb, displayAngle, canvasWidth, canvasHeight);

  const matchedReuseState = reuseState && reuseState.geometryKey === geometryKey
    ? reuseState
    : null;
  const prevConverged = matchedReuseState?.bl?.LVCONV === true
    || (blCache.ctx && blCache.key === geometryKey && blCache.ctx.LVCONV === true);
  const reuseAttempt = advancedMode ? reuseSolution : prevConverged;
  const panelReuse = reusePanel && (advancedMode ? reuseSolution : reuseAttempt);
  let ctxPanel = buildPanelContext(nb, alphaRad, {
    reusePanel: panelReuse,
    geometryKey,
    reuseState: reuseAttempt ? matchedReuseState : null,
  });
  let blCtx = null;
  let qinv = null;
  let qinvA = null;
  if (viscous && ctxPanel) {
    const maxIter = Number.isFinite(nIter) && nIter > 0 ? nIter : 20;
    const applyMach = () => {
      blCtx.MINF = Number.isFinite(mach) ? mach : 0.0;
      blCtx.MINF1 = blCtx.MINF;
    };

    const runViscous = (useReuse, allowRestore) => {
      const reuseCache = useReuse && blCache.ctx && blCache.key === geometryKey;
      if (reuseCache) {
        blCtx = blCache.ctx;
        const acrit = Number.isFinite(ncr) ? ncr : 9.0;
        blCtx.ACRIT[1] = acrit;
        blCtx.ACRIT[2] = acrit;
      } else {
        blCtx = buildBlContext(nb, ctxPanel, ncr);
        if (useReuse && allowRestore && matchedReuseState?.bl) {
          const restored = applyBlSnapshot(blCtx, matchedReuseState.bl);
          if (restored) {
            blCache.ctx = blCtx;
            blCache.key = geometryKey;
          }
        } else {
          blCache.ctx = blCtx;
          blCache.key = geometryKey;
        }
      }

      applyMach();
      return viscal(
        blCtx,
        ctxPanel,
        alphaRad,
        Number.isFinite(reynolds) ? reynolds : 0.0,
        {
          maxIter,
          logSurface: true,
          reuseSolution: useReuse,
        },
      );
    };

    let viscalResult = runViscous(reuseAttempt, reuseAttempt);
    qinv = viscalResult.qinv;
    qinvA = viscalResult.qinvA;
    if (!advancedMode && reuseAttempt && !viscalResult.converged) {
      ctxPanel = buildPanelContext(nb, alphaRad, {
        reusePanel: false,
        geometryKey,
        reuseState: null,
      });
      viscalResult = runViscous(false, false);
      qinv = viscalResult.qinv;
      qinvA = viscalResult.qinvA;
    }
  } else if (ctxPanel) {
    ({ qinv, qinvA } = specal(ctxPanel, alphaRad));
  }

  const streamlines = ctxPanel ? buildStreamlineGrid(bounds, nb, ctxPanel) : null;

  let cpData = null;
  let blLines = null;
  let coeffsDisplay = null;
  let converged = true;

  if (ctxPanel && qinv) {
    const qinf = ctxPanel.QINF ?? 1.0;
    const minf = viscous && blCtx ? blCtx.MINF ?? 0.0 : (Number.isFinite(mach) ? mach : 0.0);
    const total = ctxPanel.N + (ctxPanel.NW ?? 0);
    const qvis = (ctxPanel.QVIS && ctxPanel.QVIS.length === total + 1)
      ? ctxPanel.QVIS
      : computeQvisFromUedg(blCtx, nb, qinv);
    const cpInv = cpcalc(qinv, qinf, minf);
    const cpVis = cpcalc(qvis, qinf, minf);

    let cpUpper = [];
    let cpLower = [];
    let cpWake = [];
    if (viscous && blCtx) {
      const ile1 = blCtx.IPAN[2][1] || 0;
      const ile2 = blCtx.IPAN[2][2] || 0;
      for (let i = 1; i <= ile1; i += 1) {
        cpUpper.push({ x: ctxPanel.X[i - 1], y: ctxPanel.Y[i - 1], cp: cpVis[i] });
      }
      for (let i = ile2; i <= nb; i += 1) {
        cpLower.push({ x: ctxPanel.X[i - 1], y: ctxPanel.Y[i - 1], cp: cpVis[i] });
      }
      for (let i = nb + 1; i <= total; i += 1) {
        cpWake.push({ x: ctxPanel.X[i - 1], y: ctxPanel.Y[i - 1], cp: cpVis[i] });
      }
    } else {
      const { upperIdx, lowerIdx } = getSurfaceIndices(nb, ctxPanel);
      for (let k = 0; k < upperIdx.length; k += 1) {
        const i = upperIdx[k] + 1;
        cpUpper.push({ x: ctxPanel.X[i - 1], y: ctxPanel.Y[i - 1], cp: cpInv[i] });
      }
      for (let k = 0; k < lowerIdx.length; k += 1) {
        const i = lowerIdx[k] + 1;
        cpLower.push({ x: ctxPanel.X[i - 1], y: ctxPanel.Y[i - 1], cp: cpInv[i] });
      }
    }

    const cpInvAll = [];
    if (viscous && blCtx) {
      for (let i = 1; i <= total; i += 1) {
        cpInvAll.push({ x: ctxPanel.X[i - 1], y: ctxPanel.Y[i - 1], cp: cpInv[i] });
      }
    }

    let lePt = { x: ctxPanel.XLE, y: ctxPanel.YLE };
    let tePt = { x: ctxPanel.XTE, y: ctxPanel.YTE };
    if (!Number.isFinite(lePt.x) || !Number.isFinite(lePt.y)
      || !Number.isFinite(tePt.x) || !Number.isFinite(tePt.y)) {
      const chordPts = getChordPoints(nb);
      lePt = chordPts.le;
      tePt = chordPts.te;
    }

    cpData = {
      upper: cpUpper,
      lower: cpLower,
      wake: cpWake,
      invAll: cpInvAll,
      le: lePt,
      te: tePt,
    };
  }

  const coeffs = ctxPanel
    ? computeCoefficients(nb, ctxPanel, blCtx, alphaRad, qinvA, viscous)
    : null;
  if (coeffs) {
    coeffsDisplay = {
      CL: coeffs.cl,
      CM: coeffs.cm,
      CD: coeffs.cd,
      CDF: coeffs.cdf,
      CDP: coeffs.cdp,
      ACRIT: blCtx?.ACRIT,
      REINF1: blCtx?.REINF1,
    };
    if (!viscous) {
      coeffsDisplay.CD = coeffs.cdp;
    }
  }

  if (blCtx) {
    converged = blCtx.LVCONV === true;
    blLines = buildBoundaryLayerLines(blCtx, ctxPanel, 1.0);
  }

  let reuseSnapshot = null;
  if (viscous && blCtx && ctxPanel && (advancedMode ? reuseSolution : true)) {
    reuseSnapshot = {
      geometryKey,
      panel: snapshotPanelCtx(ctxPanel),
      bl: snapshotBlCtx(blCtx),
    };
  }

  return {
    ok: true,
    nb,
    airfoilName: geom.airfoilName,
    hinge: geom.hinge,
    xb: xb.slice(0, nb),
    yb: yb.slice(0, nb),
    bounds,
    streamlines,
    cpData,
    coeffs: coeffsDisplay,
    converged,
    blLines,
    alphaDeg,
    alphaRad,
    viscous,
    reuseState: reuseSnapshot,
  };
}

self.onmessage = (event) => {
  const { requestId, settings } = event.data;
  try {
    const payload = computeCase(settings);
    self.postMessage({ requestId, ...payload });
  } catch (err) {
    self.postMessage({
      requestId,
      ok: false,
      error: err instanceof Error ? err.message : String(err),
      errorStack: err instanceof Error ? err.stack : null,
    });
  }
};
