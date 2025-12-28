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
  const { reusePanel = false, geometryKey = '' } = opts;
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
  return panelCtx;
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
    reusePanel,
    reuseSolution,
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

  const ctxPanel = buildPanelContext(nb, alphaRad, { reusePanel, geometryKey });
  let blCtx = null;
  let qinv = null;
  let qinvA = null;
  if (viscous && ctxPanel) {
    const reuse = reuseSolution && blCache.ctx && blCache.key === geometryKey;
    if (reuse) {
      blCtx = blCache.ctx;
      const acrit = Number.isFinite(ncr) ? ncr : 9.0;
      blCtx.ACRIT[1] = acrit;
      blCtx.ACRIT[2] = acrit;
    } else {
      blCtx = buildBlContext(nb, ctxPanel, ncr);
      blCache.ctx = blCtx;
      blCache.key = geometryKey;
    }

    const maxIter = Number.isFinite(nIter) && nIter > 0 ? nIter : 20;
    blCtx.MINF = Number.isFinite(mach) ? mach : 0.0;
    blCtx.MINF1 = blCtx.MINF;

    ({ qinv, qinvA } = viscal(
      blCtx,
      ctxPanel,
      alphaRad,
      Number.isFinite(reynolds) ? reynolds : 0.0,
      {
        maxIter,
        logSurface: true,
        reuseSolution: reuse,
      },
    ));
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
    });
  }
};
