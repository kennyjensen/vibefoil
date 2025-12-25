// Minimal driver to exercise panel routines on a canonical NACA 2412.
import './naca.js';
import { scalc } from './spline.js';
import { apcalc, ncalc, psilin } from './xpanel.js';

// Construct a reference airfoil for sanity checks of panel influence routines.
function buildNaca2412(nside = 80) {
  const xx = new Float64Array(nside);
  const yt = new Float64Array(nside);
  const yc = new Float64Array(nside);
  const xb = new Float64Array(2 * nside);
  const yb = new Float64Array(2 * nside);

  const { nb } = globalThis.Naca.naca4(2412, xx, yt, yc, nside, xb, yb);

  const x = new Float64Array(nb);
  const y = new Float64Array(nb);
  for (let i = 0; i < nb; i += 1) {
    x[i] = xb[i];
    y[i] = yb[i];
  }

  return { x, y, n: nb };
}

// Build the minimal panel context required by psilin/apcalc.
function buildContext(x, y, n) {
  const s = new Float64Array(n);
  scalc(x, y, s, n);

  const nx = new Float64Array(n);
  const ny = new Float64Array(n);
  ncalc(x, y, s, n, nx, ny);

  const apanel = new Float64Array(n);

  const gam = new Float64Array(n);
  const gamu = Array.from({ length: n }, () => new Float64Array(2));
  const sig = new Float64Array(n);

  const dzdg = new Float64Array(n);
  const dzdn = new Float64Array(n);
  const dqdg = new Float64Array(n);
  const dzdm = new Float64Array(n);
  const dqdm = new Float64Array(n);

  const qf0 = new Float64Array(n);
  const qf1 = new Float64Array(n);
  const qf2 = new Float64Array(n);
  const qf3 = new Float64Array(n);

  const pi = Math.PI;

  return {
    N: n,
    X: x,
    Y: y,
    S: s,
    NX: nx,
    NY: ny,
    APANEL: apanel,
    SHARP: true,
    PI: pi,
    ANTE: 0.0,
    ASTE: 0.0,
    DSTE: 0.0,
    GAM: gam,
    GAMU: gamu,
    SIG: sig,
    QF0: qf0,
    QF1: qf1,
    QF2: qf2,
    QF3: qf3,
    QOPI: 1.0 / (4.0 * pi),
    HOPI: 1.0 / (2.0 * pi),
    ALFA: 0.0,
    QINF: 1.0,
    LIMAGE: false,
    YIMAGE: 0.0,
    DZDG: dzdg,
    DZDN: dzdn,
    DQDG: dqdg,
    DZDM: dzdm,
    DQDM: dqdm,
  };
}

// Demo call returning a local influence evaluation for diagnostics.
function runXpanelDemo() {
  const { x, y, n } = buildNaca2412(80);
  const ctx = buildContext(x, y, n);

  apcalc(ctx);

  const i = Math.floor(n / 4);
  const xi = x[i];
  const yi = y[i];
  const nxi = ctx.NX[i];
  const nyi = ctx.NY[i];

  const result = psilin(i, xi, yi, nxi, nyi, false, false, ctx);

  return {
    pointIndex: i,
    point: { x: xi, y: yi, nx: nxi, ny: nyi },
    psi: result.psi,
    psiNi: result.psiNi,
    zQinf: result.zQinf,
    zAlfa: result.zAlfa,
    qtan1: result.qtan1,
    qtan2: result.qtan2,
    qtanm: result.qtanm,
  };
}

export { runXpanelDemo };
