// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

import { scalc, splind, trisol } from './spline.js';
import { ncalc, psilin } from './xpanel.js';
import { clcalc, tecalc } from './xfoil.js';
import { gauss } from './xsolve.js';
import { atanc } from './xutils.js';

function QDES() {
  throw new Error('QDES interactive command loop is not ported.');
}

function NEWPLOTQ() {
  return null;
}

function QPLINI() {
  return null;
}

function QSPLOT() {
  return null;
}

function QSPPLT() {
  return null;
}

function IQSGET() {
  throw new Error('IQSGET cursor selection is not ported.');
}

function SPLQSP(ctx, kqsp) {
  if (ctx.NSP < 2) return;

  const nmid = ctx.NSP - 2;
  if (nmid > 0) {
    const qLocal = new Float64Array(nmid);
    const sLocal = new Float64Array(nmid);
    const qpLocal = new Float64Array(nmid);
    for (let i = 0; i < nmid; i += 1) {
      const idx = i + 1;
      qLocal[i] = ctx.QSPEC[idx][kqsp];
      sLocal[i] = ctx.SSPEC[idx];
    }
    splind(qLocal, qpLocal, sLocal, nmid, -999.0, -999.0);
    for (let i = 0; i < nmid; i += 1) {
      const idx = i + 1;
      ctx.QSPECP[idx][kqsp] = qpLocal[i];
    }
  }

  let qLocal = new Float64Array(2);
  let sLocal = new Float64Array(2);
  let qpLocal = new Float64Array(2);
  qLocal[0] = ctx.QSPEC[0][kqsp];
  qLocal[1] = ctx.QSPEC[1][kqsp];
  sLocal[0] = ctx.SSPEC[0];
  sLocal[1] = ctx.SSPEC[1];
  splind(qLocal, qpLocal, sLocal, 2, -999.0, ctx.QSPECP[1][kqsp]);
  ctx.QSPECP[0][kqsp] = qpLocal[0];

  qLocal = new Float64Array(2);
  sLocal = new Float64Array(2);
  qpLocal = new Float64Array(2);
  qLocal[0] = ctx.QSPEC[ctx.NSP - 2][kqsp];
  qLocal[1] = ctx.QSPEC[ctx.NSP - 1][kqsp];
  sLocal[0] = ctx.SSPEC[ctx.NSP - 2];
  sLocal[1] = ctx.SSPEC[ctx.NSP - 1];
  splind(qLocal, qpLocal, sLocal, 2, ctx.QSPECP[ctx.NSP - 2][kqsp], -999.0);
  ctx.QSPECP[ctx.NSP - 2][kqsp] = qpLocal[0];
}

function SMOOQ(ctx, kq1, kq2, kqsp) {
  for (let i = 0; i < ctx.NSP; i += 1) {
    ctx.W8[i] = ctx.SSPEC[i];
  }

  if (kq2 - kq1 < 2) return;

  const smool = 0.002 * (ctx.W8[ctx.NSP - 1] - ctx.W8[0]);
  const smoosq = smool * smool;

  for (let i = kq1 + 1; i <= kq2 - 1; i += 1) {
    const dsm = ctx.W8[i] - ctx.W8[i - 1];
    const dsp = ctx.W8[i + 1] - ctx.W8[i];
    const dso = 0.5 * (ctx.W8[i + 1] - ctx.W8[i - 1]);
    ctx.W1[i] = smoosq * (-1.0 / dsm) / dso;
    ctx.W2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
    ctx.W3[i] = smoosq * (-1.0 / dsp) / dso;
  }

  ctx.W2[kq1] = 1.0;
  ctx.W3[kq1] = 0.0;
  ctx.W1[kq2] = 0.0;
  ctx.W2[kq2] = 1.0;

  if (ctx.LQSLOP) {
    let i = kq1 + 1;
    let dsm = ctx.W8[i] - ctx.W8[i - 1];
    let dsp = ctx.W8[i + 1] - ctx.W8[i];
    let ds = ctx.W8[i + 1] - ctx.W8[i - 1];
    ctx.W1[i] = -1.0 / dsm - (dsm / ds) / dsm;
    ctx.W2[i] = 1.0 / dsm + (dsm / ds) / dsm + (dsm / ds) / dsp;
    ctx.W3[i] = -(dsm / ds) / dsp;
    const qsp1 = ctx.W1[i] * ctx.QSPEC[i - 1][kqsp]
      + ctx.W2[i] * ctx.QSPEC[i][kqsp]
      + ctx.W3[i] * ctx.QSPEC[i + 1][kqsp];

    i = kq2 - 1;
    dsm = ctx.W8[i] - ctx.W8[i - 1];
    dsp = ctx.W8[i + 1] - ctx.W8[i];
    ds = ctx.W8[i + 1] - ctx.W8[i - 1];
    ctx.W1[i] = (dsp / ds) / dsm;
    ctx.W2[i] = -1.0 / dsp - (dsp / ds) / dsp - (dsp / ds) / dsm;
    ctx.W3[i] = 1.0 / dsp + (dsp / ds) / dsp;
    const qsp2 = ctx.W1[i] * ctx.QSPEC[i - 1][kqsp]
      + ctx.W2[i] * ctx.QSPEC[i][kqsp]
      + ctx.W3[i] * ctx.QSPEC[i + 1][kqsp];

    ctx.QSPEC[kq1 + 1][kqsp] = qsp1;
    ctx.QSPEC[kq2 - 1][kqsp] = qsp2;
  }

  const kk = kq2 - kq1 + 1;
  const a = new Float64Array(kk);
  const b = new Float64Array(kk);
  const c = new Float64Array(kk);
  const d = new Float64Array(kk);
  for (let i = 0; i < kk; i += 1) {
    const idx = kq1 + i;
    a[i] = ctx.W2[idx];
    b[i] = ctx.W1[idx];
    c[i] = ctx.W3[idx];
    d[i] = ctx.QSPEC[idx][kqsp];
  }
  trisol(a, b, c, d, kk);
  for (let i = 0; i < kk; i += 1) {
    const idx = kq1 + i;
    ctx.QSPEC[idx][kqsp] = d[i];
  }
}

function QINCOM(qc, qinf, tklam) {
  if (tklam < 1.0e-4 || Math.abs(qc) < 1.0e-4) {
    return qc / (1.0 - tklam);
  }
  const tmp = 0.5 * (1.0 - tklam) * qinf / (qc * tklam);
  return qinf * tmp * (Math.sqrt(1.0 + 1.0 / (tklam * tmp * tmp)) - 1.0);
}

function GAMQSP(ctx, kqsp) {
  ctx.ALQSP[kqsp] = ctx.ALGAM;
  ctx.CLQSP[kqsp] = ctx.CLGAM;
  ctx.CMQSP[kqsp] = ctx.CMGAM;
  for (let i = 0; i < ctx.NSP; i += 1) {
    ctx.QSPEC[i][kqsp] = ctx.QGAMM[i];
  }
  ctx.QDOF0 = 0.0;
  ctx.QDOF1 = 0.0;
  ctx.QDOF2 = 0.0;
  ctx.QDOF3 = 0.0;
  SPLQSP(ctx, kqsp);
  if (!ctx.LIQSET) {
    ctx.IQ1 = 0;
    ctx.IQ2 = ctx.NSP - 1;
  }
}

function SYMQSP(ctx, kqsp) {
  ctx.ALQSP[kqsp] = 0.0;
  ctx.CLQSP[kqsp] = 0.0;
  ctx.CMQSP[kqsp] = 0.0;
  const sspmid = 0.5 * (ctx.SSPEC[ctx.NSP - 1] - ctx.SSPEC[0]);
  const mid = Math.floor((ctx.NSP - 1) / 2);
  for (let i = 0; i <= mid; i += 1) {
    ctx.SSPEC[i] = sspmid + 0.5 * (ctx.SSPEC[i] - ctx.SSPEC[ctx.NSP - i - 1]);
    ctx.QSPEC[i][kqsp] = 0.5 * (ctx.QSPEC[i][kqsp] - ctx.QSPEC[ctx.NSP - i - 1][kqsp]);
  }
  for (let i = mid + 1; i < ctx.NSP; i += 1) {
    ctx.SSPEC[i] = -ctx.SSPEC[ctx.NSP - i - 1] + 2.0 * sspmid;
    ctx.QSPEC[i][kqsp] = -ctx.QSPEC[ctx.NSP - i - 1][kqsp];
  }
  ctx.QDOF0 = 0.0;
  ctx.QDOF1 = 0.0;
  ctx.QDOF2 = 0.0;
  ctx.QDOF3 = 0.0;
  SPLQSP(ctx, kqsp);
}

function MIXED(ctx, kqsp, niterq) {
  const bwt = 0.1;
  const n = ctx.N;

  ctx.COSA = Math.cos(ctx.ALFA);
  ctx.SINA = Math.sin(ctx.ALFA);
  scalc(ctx.X, ctx.Y, ctx.S, n);

  for (let i = 0; i < n; i += 1) {
    ctx.QF0[i] = 0.0;
    ctx.QF1[i] = 0.0;
    ctx.QF2[i] = 0.0;
    ctx.QF3[i] = 0.0;
  }

  for (let i = ctx.IQ1; i <= ctx.IQ2; i += 1) {
    const fs = (ctx.S[i] - ctx.S[ctx.IQ1]) / (ctx.S[ctx.IQ2] - ctx.S[ctx.IQ1]);
    ctx.QF0[i] = 1.0 - fs;
    ctx.QF1[i] = fs;
    if (ctx.LCPXX) {
      ctx.QF2[i] = Math.exp(-5.0 * fs);
      ctx.QF3[i] = Math.exp(-5.0 * (1.0 - fs));
    } else {
      ctx.QF2[i] = 0.0;
      ctx.QF3[i] = 0.0;
    }
    ctx.GAM[i] = ctx.QSPEC[i][kqsp]
      + ctx.QDOF0 * ctx.QF0[i]
      + ctx.QDOF1 * ctx.QF1[i]
      + ctx.QDOF2 * ctx.QF2[i]
      + ctx.QDOF3 * ctx.QF3[i];
  }

  for (let iter = 0; iter < niterq; iter += 1) {
    const nsys = n + 5;
    for (let i = 0; i < nsys; i += 1) {
      for (let j = 0; j < nsys; j += 1) {
        ctx.Q[i][j] = 0.0;
      }
    }

    ncalc(ctx.X, ctx.Y, ctx.S, n, ctx.NX, ctx.NY);

    for (let i = 0; i < n; i += 1) {
      const { psi, psiNi } = psilin(i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i], true, false, ctx);
      ctx.DZDN[i] = ctx.DZDN[i] + psiNi;

      for (let j = 0; j <= ctx.IQ1 - 1; j += 1) {
        ctx.Q[i][j] = ctx.Q[i][j] + ctx.DZDG[j];
      }
      for (let j = ctx.IQ1; j <= ctx.IQ2; j += 1) {
        ctx.Q[i][j] = ctx.Q[i][j] + ctx.DZDN[j];
      }
      for (let j = ctx.IQ2 + 1; j < n; j += 1) {
        ctx.Q[i][j] = ctx.Q[i][j] + ctx.DZDG[j];
      }

      ctx.DQ[i] = ctx.PSIO - psi;

      ctx.Q[i][n] -= 1.0;
      ctx.Q[i][n + 1] += ctx.Z_QDOF0;
      ctx.Q[i][n + 2] += ctx.Z_QDOF1;
      ctx.Q[i][n + 3] += ctx.Z_QDOF2;
      ctx.Q[i][n + 4] += ctx.Z_QDOF3;
    }

    ctx.DQ[n] = -(ctx.GAM[0] + ctx.GAM[n - 1]);
    GAMLIN(ctx, n, 0, 1.0);
    GAMLIN(ctx, n, n - 1, 1.0);

    if (ctx.SHARP) {
      const ag1 = Math.atan2(-ctx.YP[0], -ctx.XP[0]);
      const ag2 = atanc(ctx.YP[n - 1], ctx.XP[n - 1], ag1);
      const abis = 0.5 * (ag1 + ag2);
      const cbis = Math.cos(abis);
      const sbis = Math.sin(abis);

      const ds1 = Math.hypot(ctx.X[0] - ctx.X[1], ctx.Y[0] - ctx.Y[1]);
      const ds2 = Math.hypot(ctx.X[n - 1] - ctx.X[n - 2], ctx.Y[n - 1] - ctx.Y[n - 2]);
      const dsmin = Math.min(ds1, ds2);

      const xbis = ctx.XTE - bwt * dsmin * cbis;
      const ybis = ctx.YTE - bwt * dsmin * sbis;
      const { psiNi } = psilin(-1, xbis, ybis, -sbis, cbis, false, true, ctx);
      const res = psiNi;

      for (let j = 0; j < nsys; j += 1) {
        ctx.Q[n - 1][j] = 0.0;
      }

      for (let j = 0; j < n; j += 1) {
        GAMLIN(ctx, n - 1, j, ctx.DQDG[j]);
        ctx.Q[n - 1][j] = ctx.DQDG[j];
      }

      ctx.Q[n - 1][n] = 0.0;
      ctx.DQ[n - 1] = -res;
    }

    ctx.Q[n + 1][ctx.IQ1] = 1.0;
    ctx.DQ[n + 1] = 0.0;
    ctx.Q[n + 2][ctx.IQ2] = 1.0;
    ctx.DQ[n + 2] = 0.0;

    if (ctx.IQ1 > 0 && ctx.LCPXX) {
      const res = ctx.GAM[ctx.IQ1 - 1] - 2.0 * ctx.GAM[ctx.IQ1] + ctx.GAM[ctx.IQ1 + 1]
        - (ctx.QSPEC[ctx.IQ1 - 1][kqsp] - 2.0 * ctx.QSPEC[ctx.IQ1][kqsp] + ctx.QSPEC[ctx.IQ1 + 1][kqsp]);
      GAMLIN(ctx, n + 3, ctx.IQ1 - 1, 1.0);
      GAMLIN(ctx, n + 3, ctx.IQ1, -2.0);
      GAMLIN(ctx, n + 3, ctx.IQ1 + 1, 1.0);
      ctx.DQ[n + 3] = -res;
    } else {
      ctx.Q[n + 3][n + 3] = 1.0;
      ctx.DQ[n + 3] = -ctx.QDOF2;
    }

    if (ctx.IQ2 < n - 1 && ctx.LCPXX) {
      const res = ctx.GAM[ctx.IQ2 - 1] - 2.0 * ctx.GAM[ctx.IQ2] + ctx.GAM[ctx.IQ2 + 1]
        - (ctx.QSPEC[ctx.IQ2 - 1][kqsp] - 2.0 * ctx.QSPEC[ctx.IQ2][kqsp] + ctx.QSPEC[ctx.IQ2 + 1][kqsp]);
      GAMLIN(ctx, n + 4, ctx.IQ2 - 1, 1.0);
      GAMLIN(ctx, n + 4, ctx.IQ2, -2.0);
      GAMLIN(ctx, n + 4, ctx.IQ2 + 1, 1.0);
      ctx.DQ[n + 4] = -res;
    } else {
      ctx.Q[n + 4][n + 4] = 1.0;
      ctx.DQ[n + 4] = -ctx.QDOF3;
    }

    gauss(nsys, ctx.Q, ctx.DQ, 1);

    let dnmax = 0.0;
    let dgmax = 0.0;
    let maxGeomStep = 0.0;
    for (let i = ctx.IQ1; i <= ctx.IQ2; i += 1) {
      const step = Math.abs(ctx.DQ[i]);
      if (step > maxGeomStep) {
        maxGeomStep = step;
      }
    }
    const maxStep = Number.isFinite(ctx.QDES_MAXSTEP) ? ctx.QDES_MAXSTEP : 0.02;
    let relax = 1.0;
    if (maxGeomStep > maxStep && maxGeomStep > 0.0) {
      relax = maxStep / maxGeomStep;
      for (let i = 0; i < nsys; i += 1) {
        ctx.DQ[i] *= relax;
      }
    }

    for (let i = 0; i <= ctx.IQ1 - 1; i += 1) {
      ctx.GAM[i] = ctx.GAM[i] + ctx.DQ[i];
      if (Math.abs(ctx.DQ[i]) > Math.abs(dgmax)) {
        dgmax = ctx.DQ[i];
      }
    }

    for (let i = ctx.IQ1; i <= ctx.IQ2; i += 1) {
      ctx.X[i] = ctx.X[i] + ctx.NX[i] * ctx.DQ[i];
      ctx.Y[i] = ctx.Y[i] + ctx.NY[i] * ctx.DQ[i];
      if (Math.abs(ctx.DQ[i]) > Math.abs(dnmax)) {
        dnmax = ctx.DQ[i];
      }
    }

    for (let i = ctx.IQ2 + 1; i < n; i += 1) {
      ctx.GAM[i] = ctx.GAM[i] + ctx.DQ[i];
      if (Math.abs(ctx.DQ[i]) > Math.abs(dgmax)) {
        dgmax = ctx.DQ[i];
      }
    }

    ctx.PSIO = ctx.PSIO + ctx.DQ[n];
    ctx.QDOF0 = ctx.QDOF0 + ctx.DQ[n + 1];
    ctx.QDOF1 = ctx.QDOF1 + ctx.DQ[n + 2];
    ctx.QDOF2 = ctx.QDOF2 + ctx.DQ[n + 3];
    ctx.QDOF3 = ctx.QDOF3 + ctx.DQ[n + 4];

    ctx.COSA = Math.cos(ctx.ALFA);
    ctx.SINA = Math.sin(ctx.ALFA);
    scalc(ctx.X, ctx.Y, ctx.S, n);

    for (let i = ctx.IQ1; i <= ctx.IQ2; i += 1) {
      ctx.GAM[i] = ctx.QSPEC[i][kqsp]
        + ctx.QDOF0 * ctx.QF0[i]
        + ctx.QDOF1 * ctx.QF1[i]
        + ctx.QDOF2 * ctx.QF2[i]
        + ctx.QDOF3 * ctx.QF3[i];
    }

    tecalc(ctx);
    const coeffs = clcalc(n, ctx.X, ctx.Y, ctx.GAM, ctx.GAM_A, ctx.ALFA, ctx.MINF, ctx.QINF, ctx.XCMREF, ctx.YCMREF);
    ctx.CL = coeffs.cl;
    ctx.CM = coeffs.cm;
    ctx.CDP = coeffs.cdp;
    ctx.CL_ALF = coeffs.clAlf;
    ctx.CL_MSQ = coeffs.clMsq;

    if (ctx.QDES_DEBUG) {
      console.log('[QDES] iter', {
        pass: ctx.QDES_PASS || null,
        iter: iter + 1,
        dnmax,
        dgmax,
        dpsio: ctx.DQ[n],
        dQf: [ctx.DQ[n + 1], ctx.DQ[n + 2], ctx.DQ[n + 3], ctx.DQ[n + 4]],
        psio: ctx.PSIO,
        relax,
        maxGeomStep,
      });
    }

    if (Math.abs(dnmax) < 5.0e-5 && Math.abs(dgmax) < 5.0e-4) {
      return;
    }
  }
}

function GAMLIN(ctx, i, j, coef) {
  const n = ctx.N;
  if (j >= ctx.IQ1 && j <= ctx.IQ2) {
    ctx.Q[i][n + 1] += coef * ctx.QF0[j];
    ctx.Q[i][n + 2] += coef * ctx.QF1[j];
    ctx.Q[i][n + 3] += coef * ctx.QF2[j];
    ctx.Q[i][n + 4] += coef * ctx.QF3[j];
  } else {
    ctx.Q[i][j] += coef;
  }
}

export {
  QDES,
  NEWPLOTQ,
  QPLINI,
  QSPLOT,
  QSPPLT,
  IQSGET,
  SPLQSP,
  SMOOQ,
  QINCOM,
  GAMQSP,
  SYMQSP,
  MIXED,
  GAMLIN,
};
