import '../../js/naca.js';
import { pangen, cpcalc, tecalc } from '../../js/xfoil.js';
import { apcalc, ncalc, ggcalc, qiset } from '../../js/xpanel.js';
import { scalc, segspl } from '../../js/spline.js';
import { createMatrix } from '../../js/arrays.js';

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
  const { qinv } = qiset(ctxPanel, alphaRad);
  const cpInv = cpcalc(qinv, ctxPanel.QINF ?? 1.0, 0.0);

  const points = [];
  for (let i = 1; i <= nb; i += 1) {
    points.push({ x: ctxPanel.X[i - 1], y: ctxPanel.Y[i - 1], cp: cpInv[i] });
  }

  let leIdx = 0;
  let teIdx = 0;
  for (let i = 1; i < nb; i += 1) {
    if (ctxPanel.X[i] < ctxPanel.X[leIdx]) leIdx = i;
    if (ctxPanel.X[i] > ctxPanel.X[teIdx]) teIdx = i;
  }
  const le = { x: ctxPanel.X[leIdx], y: ctxPanel.Y[leIdx] };
  const te = { x: ctxPanel.X[teIdx], y: ctxPanel.Y[teIdx] };

  const input = {
    nb,
    alpha: alphaRad,
    minf: 0.0,
    qinf: 1.0,
    sharp: ctxPanel.SHARP ? 1 : 0,
    ante: ctxPanel.ANTE,
    aste: ctxPanel.ASTE,
    dste: ctxPanel.DSTE,
    qopi: ctxPanel.QOPI,
    hopi: ctxPanel.HOPI,
    pi: ctxPanel.PI,
    arrays: {
      x: Array.from(ctxPanel.X.slice(0, nb)),
      y: Array.from(ctxPanel.Y.slice(0, nb)),
      s: Array.from(ctxPanel.S.slice(0, nb)),
      xp: Array.from(ctxPanel.XP.slice(0, nb)),
      yp: Array.from(ctxPanel.YP.slice(0, nb)),
      nx: Array.from(ctxPanel.NX.slice(0, nb)),
      ny: Array.from(ctxPanel.NY.slice(0, nb)),
      sig: Array.from(ctxPanel.SIG.slice(0, nb)).map(() => 0.0),
      qf0: Array.from(ctxPanel.QF0),
      qf1: Array.from(ctxPanel.QF1),
      qf2: Array.from(ctxPanel.QF2),
      qf3: Array.from(ctxPanel.QF3),
    },
    jsPoints: points,
    le,
    te,
  };

  process.stdout.write(JSON.stringify(input));
}

main();
