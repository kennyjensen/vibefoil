// High-level XFOIL utilities: paneling and coefficient integrations.

import {
  scalc,
  segspl,
  seval,
  deval,
  curv,
  trisol,
} from './spline.js';
import { lefind } from './xgeom.js';

const PANGEN_DEFAULTS = {
  npan: 160,
  cvpar: 1.0,
  cterat: 0.15,
  ctrrat: 0.2,
  xsref1: 1.0,
  xsref2: 1.0,
  xpref1: 1.0,
  xpref2: 1.0,
};

// Compressibility-corrected Cp (CPCALC), matching XFOIL's formulation.
function cpcalc(qvals, qinf, minf) {
  const n = qvals.length - 1;
  const cp = new Float64Array(n + 1);
  // Prandtl-Glauert compressibility correction.
  // Prandtl-Glauert compressibility correction.
  const beta = Math.sqrt(Math.max(1.0 - minf ** 2, 0.0));
  const bfac = 0.5 * minf ** 2 / (1.0 + beta);
  for (let i = 1; i <= n; i += 1) {
    const q = qvals[i] ?? 0.0;
    const cpinc = 1.0 - (q / qinf) ** 2;
    const den = beta + bfac * cpinc;
    cp[i] = den !== 0.0 ? cpinc / den : cpinc;
  }
  return cp;
}

// Integrate surface pressures for lift, moment, and pressure drag (CLCALC).
function clcalc(nb, x, y, gam, gamA, alfa, minf, qinf, xref = 0.25, yref = 0.0) {
  const sa = Math.sin(alfa);
  const ca = Math.cos(alfa);

  const beta = Math.sqrt(Math.max(1.0 - minf ** 2, 0.0));
  const betaMsq = beta === 0.0 ? 0.0 : -0.5 / beta;
  const bfac = 0.5 * minf ** 2 / (1.0 + beta);
  const bfacMsq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * betaMsq;

  let cl = 0.0;
  let cm = 0.0;
  let cdp = 0.0;
  let clAlf = 0.0;
  let clMsq = 0.0;

  // Initialize Cp at panel 1.
  let cginc = 1.0 - (gam[0] / qinf) ** 2;
  let cpg1 = cginc / (beta + bfac * cginc);
  let cpg1Msq = -cpg1 / (beta + bfac * cginc) * (betaMsq + bfacMsq * cginc);
  let cpiGam = -2.0 * gam[0] / qinf ** 2;
  let cpcCpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
  let cpg1Alf = cpcCpi * cpiGam * (gamA?.[0] ?? 0.0);

  for (let i = 0; i < nb; i += 1) {
    const ip = i === nb - 1 ? 0 : i + 1;
    cginc = 1.0 - (gam[ip] / qinf) ** 2;
    const cpg2 = cginc / (beta + bfac * cginc);
    const cpg2Msq = -cpg2 / (beta + bfac * cginc) * (betaMsq + bfacMsq * cginc);
    cpiGam = -2.0 * gam[ip] / qinf ** 2;
    cpcCpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
    const cpg2Alf = cpcCpi * cpiGam * (gamA?.[ip] ?? 0.0);

    // Accumulate CL/CM/CDP integrals.
    const dx = (x[ip] - x[i]) * ca + (y[ip] - y[i]) * sa;
    const dy = (y[ip] - y[i]) * ca - (x[ip] - x[i]) * sa;
    const dg = cpg2 - cpg1;

    const ax = (0.5 * (x[ip] + x[i]) - xref) * ca + (0.5 * (y[ip] + y[i]) - yref) * sa;
    const ay = (0.5 * (y[ip] + y[i]) - yref) * ca - (0.5 * (x[ip] + x[i]) - xref) * sa;
    const ag = 0.5 * (cpg2 + cpg1);

    const dxAlf = -(x[ip] - x[i]) * sa + (y[ip] - y[i]) * ca;
    const agAlf = 0.5 * (cpg2Alf + cpg1Alf);
    const agMsq = 0.5 * (cpg2Msq + cpg1Msq);

    cl += dx * ag;
    cdp -= dy * ag;
    cm -= dx * (ag * ax + dg * dx / 12.0) + dy * (ag * ay + dg * dy / 12.0);
    clAlf += dx * agAlf + ag * dxAlf;
    clMsq += dx * agMsq;

    cpg1 = cpg2;
    cpg1Alf = cpg2Alf;
    cpg1Msq = cpg2Msq;
  }

  return { cl, cm, cdp, clAlf, clMsq };
}

// Drag calculation from wake extrapolation and skin friction (CDCALC).
function cdcalc(ctxPanel, blCtx, alfa, qinf) {
  let cd = 0.0;
  let cdf = 0.0;

  if (blCtx?.LVISC && blCtx?.LBLINI) {
    // Extrapolate wake to downstream infinity using Squire-Young relation.
    const thwake = blCtx.THET[blCtx.NBL[2]][2];
    const urat = blCtx.UEDG[blCtx.NBL[2]][2] / qinf;
    const uewake = blCtx.UEDG[blCtx.NBL[2]][2]
      * (1.0 - blCtx.TKLAM) / (1.0 - blCtx.TKLAM * urat ** 2);
    const shwake = blCtx.DSTR[blCtx.NBL[2]][2] / blCtx.THET[blCtx.NBL[2]][2];
    cd = 2.0 * thwake * (uewake / qinf) ** (0.5 * (5.0 + shwake));
  }

  // Friction drag coefficient from shear stress integration.
  const sa = Math.sin(alfa);
  const ca = Math.cos(alfa);
  const x = ctxPanel.X;
  const y = ctxPanel.Y;
  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 3; ibl <= blCtx.IBLTE[is]; ibl += 1) {
      const i = blCtx.IPAN[ibl][is] - 1;
      const im = blCtx.IPAN[ibl - 1][is] - 1;
      const dx = (x[i] - x[im]) * ca + (y[i] - y[im]) * sa;
      cdf += 0.5 * (blCtx.TAU[ibl][is] + blCtx.TAU[ibl - 1][is]) * dx * 2.0 / qinf ** 2;
    }
  }

  return { cd, cdf };
}

// Trailing-edge geometry diagnostics (TECALC).
// Calculates total and projected TE gap areas and flags sharp trailing edges.
function tecalc(ctx) {
  const { N: n, X: x, Y: y, XP: xp, YP: yp } = ctx;

  const dxte = x[0] - x[n - 1];
  const dyte = y[0] - y[n - 1];
  const dxs = 0.5 * (-xp[0] + xp[n - 1]);
  const dys = 0.5 * (-yp[0] + yp[n - 1]);

  const ante = dxs * dyte - dys * dxte;
  const aste = dxs * dxte + dys * dyte;
  const dste = Math.sqrt(dxte * dxte + dyte * dyte);

  let xmin = Infinity;
  let xmax = -Infinity;
  for (let i = 0; i < n; i += 1) {
    xmin = Math.min(xmin, x[i]);
    xmax = Math.max(xmax, x[i]);
  }
  const chord = Math.max(xmax - xmin, 1.0e-12);
  const sharp = dste < 0.0001 * chord;

  ctx.ANTE = ante;
  ctx.ASTE = aste;
  ctx.DSTE = dste;
  ctx.SHARP = sharp;
}

// Panel distribution generator (PANGEN) ported from XFOIL xfoil.f.
function pangen(xb, yb, nb, params = {}) {
  if (nb < 2) {
    return { x: new Float64Array(0), y: new Float64Array(0), s: new Float64Array(0), n: 0 };
  }

  const {
    npan,
    cvpar,
    cterat,
    ctrrat,
    xsref1,
    xsref2,
    xpref1,
    xpref2,
  } = { ...PANGEN_DEFAULTS, ...params };

  // Arc-length spline setup for buffer airfoil.
  const sb = new Float64Array(nb);
  const xbp = new Float64Array(nb);
  const ybp = new Float64Array(nb);
  scalc(xb, yb, sb, nb);
  segspl(xb, xbp, sb, nb);
  segspl(yb, ybp, sb, nb);

  // Reference length for curvature scaling.
  const sbref = 0.5 * (sb[nb - 1] - sb[0]);
  const w5 = new Float64Array(nb);
  for (let i = 0; i < nb; i += 1) {
    w5[i] = Math.abs(curv(sb[i], xb, xbp, yb, ybp, sb, nb)) * sbref;
  }

  // Find LE parameter.
  const sble = lefind(xb, xbp, yb, ybp, sb, nb);
  let cvle = Math.abs(curv(sble, xb, xbp, yb, ybp, sb, nb)) * sbref;

  let ible1 = 0;
  // Detect duplicate-LE points in the buffer.
  for (let i = 0; i < nb - 1; i += 1) {
    if (sble === sb[i] && sble === sb[i + 1]) {
      ible1 = i + 1;
      break;
    }
  }

  const xble = seval(sble, xb, xbp, sb, nb);
  const yble = seval(sble, yb, ybp, sb, nb);
  const xbte = 0.5 * (xb[0] + xb[nb - 1]);
  const ybte = 0.5 * (yb[0] + yb[nb - 1]);
  const chbsq = (xbte - xble) ** 2 + (ybte - yble) ** 2;

  // Average curvature around LE to avoid singular spikes.
  const nk = 3;
  let cvsum = 0.0;
  for (let k = -nk; k <= nk; k += 1) {
    const frac = k / nk;
    const sbk = sble + frac * sbref / Math.max(cvle, 20.0);
    const cvk = Math.abs(curv(sbk, xb, xbp, yb, ybp, sb, nb)) * sbref;
    cvsum += cvk;
  }
  let cvavg = cvsum / (2 * nk + 1);
  if (ible1 !== 0) cvavg = 10.0;

  // Set up tri-diagonal system for smoothed curvature.
  const cc = 6.0 * cvpar;
  const cvte = cvavg * cterat;
  w5[0] = cvte;
  w5[nb - 1] = cvte;

  const smool = Math.max(1.0 / Math.max(cvavg, 20.0), 0.25 / (npan / 2));
  const smoosq = (smool * sbref) ** 2;

  const w1 = new Float64Array(nb);
  const w2 = new Float64Array(nb);
  const w3 = new Float64Array(nb);

  w2[0] = 1.0;
  w3[0] = 0.0;
  for (let i = 1; i <= nb - 2; i += 1) {
    const dsm = sb[i] - sb[i - 1];
    const dsp = sb[i + 1] - sb[i];
    const dso = 0.5 * (sb[i + 1] - sb[i - 1]);

    // Leave curvature at corner point unchanged.
    if (dsm === 0.0 || dsp === 0.0) {
      w1[i] = 0.0;
      w2[i] = 1.0;
      w3[i] = 0.0;
    } else {
      w1[i] = smoosq * (-1.0 / dsm) / dso;
      w2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
      w3[i] = smoosq * (-1.0 / dsp) / dso;
    }
  }
  w1[nb - 1] = 0.0;
  w2[nb - 1] = 1.0;

  // Fix curvature at LE point by modifying adjacent equations.
  for (let i = 1; i <= nb - 2; i += 1) {
    const i1 = i + 1;
    if (sb[i] === sble || i1 === ible1 || i1 === ible1 + 1) {
      w1[i] = 0.0;
      w2[i] = 1.0;
      w3[i] = 0.0;
      w5[i] = cvle;
    } else if (sb[i - 1] < sble && sb[i] > sble) {
      let dsm = sb[i - 1] - sb[i - 2];
      let dsp = sble - sb[i - 1];
      let dso = 0.5 * (sble - sb[i - 2]);
      w1[i - 1] = smoosq * (-1.0 / dsm) / dso;
      w2[i - 1] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
      w3[i - 1] = 0.0;
      w5[i - 1] += smoosq * cvle / (dsp * dso);

      dsm = sb[i] - sble;
      dsp = sb[i + 1] - sb[i];
      dso = 0.5 * (sb[i + 1] - sble);
      w1[i] = 0.0;
      w2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
      w3[i] = smoosq * (-1.0 / dsp) / dso;
      w5[i] += smoosq * cvle / (dsm * dso);
      break;
    }
  }

  // Set artificial curvature at bunching points and fix it there.
  for (let i = 1; i <= nb - 2; i += 1) {
    const xoc = ((xb[i] - xble) * (xbte - xble) + (yb[i] - yble) * (ybte - yble)) / chbsq;
    if (sb[i] < sble) {
      if (xoc > xsref1 && xoc < xsref2) {
        w1[i] = 0.0;
        w2[i] = 1.0;
        w3[i] = 0.0;
        w5[i] = cvle * ctrrat;
      }
    } else if (xoc > xpref1 && xoc < xpref2) {
      w1[i] = 0.0;
      w2[i] = 1.0;
      w3[i] = 0.0;
      w5[i] = cvle * ctrrat;
    }
  }

  // Solve for smoothed curvature array W5.
  if (ible1 === 0) {
    trisol(w2, w1, w3, w5, nb);
  } else {
    const firstLen = ible1;
    trisol(w2.subarray(0, firstLen), w1.subarray(0, firstLen), w3.subarray(0, firstLen), w5, firstLen, 0);
    const secondLen = nb - firstLen;
    if (secondLen > 0) {
      trisol(
        w2.subarray(firstLen),
        w1.subarray(firstLen),
        w3.subarray(firstLen),
        w5,
        secondLen,
        firstLen,
      );
    }
  }

  // Find max curvature.
  let cvmax = 0.0;
  for (let i = 0; i < nb; i += 1) {
    cvmax = Math.max(cvmax, Math.abs(w5[i]));
  }
  // Normalize curvature array.
  for (let i = 0; i < nb; i += 1) {
    w5[i] /= cvmax;
  }

  // Spline curvature array.
  const w6 = new Float64Array(nb);
  segspl(w5, w6, sb, nb);

  // Set initial guess for node positions uniform in s (IPFAC refinement).
  const ipfac = 5;
  const n = npan;
  const nn = ipfac * (n - 1) + 1;
  const snew = new Float64Array(nn);

  // Ratio of panel at TE to one away from TE.
  const rdste = 0.667;
  const rtf = (rdste - 1.0) * 2.0 + 1.0;

  let nn1 = 0;
  if (ible1 === 0) {
    const dsavg = (sb[nb - 1] - sb[0]) / ((nn - 3) + 2.0 * rtf);
    snew[0] = sb[0];
    for (let i = 1; i <= nn - 2; i += 1) {
      snew[i] = sb[0] + dsavg * ((i - 1) + rtf);
    }
    snew[nn - 1] = sb[nb - 1];
  } else {
    const nfrac1 = Math.floor((n * ible1) / nb);
    nn1 = ipfac * (nfrac1 - 1) + 1;
    const dsavg1 = (sble - sb[0]) / ((nn1 - 2) + rtf);
    snew[0] = sb[0];
    for (let i = 1; i <= nn1 - 1; i += 1) {
      snew[i] = sb[0] + dsavg1 * ((i - 1) + rtf);
    }

    const nn2 = nn - nn1 + 1;
    const dsavg2 = (sb[nb - 1] - sble) / ((nn2 - 2) + rtf);
    for (let i = 1; i <= nn2 - 2; i += 1) {
      snew[i - 1 + nn1] = sble + dsavg2 * ((i - 1) + rtf);
    }
    snew[nn - 1] = sb[nb - 1];
  }

  // Newton iteration loop for new node positions.
  const t1 = new Float64Array(nn);
  const t2 = new Float64Array(nn);
  const t3 = new Float64Array(nn);
  const t4 = new Float64Array(nn);

  for (let iter = 0; iter < 20; iter += 1) {
    let cv1 = seval(snew[0], w5, w6, sb, nb);
    let cv2 = seval(snew[1], w5, w6, sb, nb);
    let cvs1 = deval(snew[0], w5, w6, sb, nb);
    let cvs2 = deval(snew[1], w5, w6, sb, nb);

    let cavm = Math.sqrt(cv1 * cv1 + cv2 * cv2);
    let cavmS1 = cavm === 0.0 ? 0.0 : cvs1 * cv1 / cavm;
    let cavmS2 = cavm === 0.0 ? 0.0 : cvs2 * cv2 / cavm;

    for (let i = 1; i <= nn - 2; i += 1) {
      const dsm = snew[i] - snew[i - 1];
      const dsp = snew[i] - snew[i + 1];
      const cv3 = seval(snew[i + 1], w5, w6, sb, nb);
      const cvs3 = deval(snew[i + 1], w5, w6, sb, nb);

      const cavp = Math.sqrt(cv3 * cv3 + cv2 * cv2);
      const cavpS2 = cavp === 0.0 ? 0.0 : cvs2 * cv2 / cavp;
      const cavpS3 = cavp === 0.0 ? 0.0 : cvs3 * cv3 / cavp;

      const fm = cc * cavm + 1.0;
      const fp = cc * cavp + 1.0;
      const rez = dsp * fp + dsm * fm;

      t1[i] = -fm + cc * dsm * cavmS1;
      t2[i] = fp + fm + cc * (dsp * cavpS2 + dsm * cavmS2);
      t3[i] = -fp + cc * dsp * cavpS3;
      t4[i] = -rez;

      cv1 = cv2;
      cv2 = cv3;
      cvs1 = cvs2;
      cvs2 = cvs3;
      cavm = cavp;
      cavmS1 = cavpS2;
      cavmS2 = cavpS3;
    }

    // Apply boundary conditions at TE.
    t2[0] = 1.0;
    t3[0] = 0.0;
    t4[0] = 0.0;
    t1[nn - 1] = 0.0;
    t2[nn - 1] = 1.0;
    t4[nn - 1] = 0.0;

    // Enforce TE spacing ratio if requested.
    if (rtf !== 1.0) {
      t4[1] = -((snew[1] - snew[0]) + rtf * (snew[1] - snew[2]));
      t1[1] = -1.0;
      t2[1] = 1.0 + rtf;
      t3[1] = -rtf;

      t4[nn - 2] = -((snew[nn - 2] - snew[nn - 1]) + rtf * (snew[nn - 2] - snew[nn - 3]));
      t3[nn - 2] = -1.0;
      t2[nn - 2] = 1.0 + rtf;
      t1[nn - 2] = -rtf;
    }

    // Constrain LE station if specified.
    if (ible1 !== 0) {
      const idx = nn1 - 1;
      t1[idx] = 0.0;
      t2[idx] = 1.0;
      t3[idx] = 0.0;
      t4[idx] = sble - snew[idx];
    }

    trisol(t2, t1, t3, t4, nn);

    // Under-relaxation to keep spacing ratios reasonable.
    let rlx = 1.0;
    let dmax = 0.0;
    for (let i = 0; i < nn - 1; i += 1) {
      const ds = snew[i + 1] - snew[i];
      const dds = t4[i + 1] - t4[i];
      const dsrat = 1.0 + rlx * dds / ds;
      if (dsrat > 4.0) rlx = (4.0 - 1.0) * ds / dds;
      if (dsrat < 0.2) rlx = (0.2 - 1.0) * ds / dds;
      dmax = Math.max(Math.abs(t4[i]), dmax);
    }

    for (let i = 1; i <= nn - 2; i += 1) {
      snew[i] += rlx * t4[i];
    }

    if (Math.abs(dmax) < 1.0e-3) break;
  }

  const x = new Float64Array(n + 20);
  const y = new Float64Array(n + 20);
  const s = new Float64Array(n + 20);

  // Interpolate to final panel coordinates.
  for (let i = 0; i < n; i += 1) {
    const ind = ipfac * i;
    s[i] = snew[ind];
    x[i] = seval(snew[ind], xb, xbp, sb, nb);
    y[i] = seval(snew[ind], yb, ybp, sb, nb);
  }

  // Insert corner points from original buffer if needed.
  let nout = n;
  for (let ib = 0; ib < nb - 1; ib += 1) {
    if (sb[ib] !== sb[ib + 1]) continue;
    const sbcorn = sb[ib];
    const xbcorn = xb[ib];
    const ybcorn = yb[ib];

    for (let i = 0; i < nout; i += 1) {
      if (s[i] <= sbcorn) continue;
      for (let j = nout; j >= i; j -= 1) {
        x[j + 1] = x[j];
        y[j + 1] = y[j];
        s[j + 1] = s[j];
      }
      nout += 1;
      x[i] = xbcorn;
      y[i] = ybcorn;
      s[i] = sbcorn;

      if (i - 2 >= 0) {
        s[i - 1] = 0.5 * (s[i] + s[i - 2]);
        x[i - 1] = seval(s[i - 1], xb, xbp, sb, nb);
        y[i - 1] = seval(s[i - 1], yb, ybp, sb, nb);
      }
      if (i + 2 <= nout - 1) {
        s[i + 1] = 0.5 * (s[i] + s[i + 2]);
        x[i + 1] = seval(s[i + 1], xb, xbp, sb, nb);
        y[i + 1] = seval(s[i + 1], yb, ybp, sb, nb);
      }
      break;
    }
  }

  const xout = x.subarray(0, nout);
  const yout = y.subarray(0, nout);
  const sout = s.subarray(0, nout);

  const xp = new Float64Array(nout);
  const yp = new Float64Array(nout);
  scalc(xout, yout, sout, nout);
  segspl(xout, xp, sout, nout);
  segspl(yout, yp, sout, nout);

  return { x: xout, y: yout, s: sout, xp, yp, n: nout };
}

function computeCoefficients(nb, ctxPanel, blCtx, alphaRad, qinvA, viscous) {
  const qinf = viscous && Number.isFinite(blCtx?.QINF) ? blCtx.QINF : (ctxPanel.QINF ?? 1.0);
  const minf = viscous ? blCtx?.MINF ?? 0.0 : 0.0;
  const alfa = viscous && Number.isFinite(blCtx?.ALFA) ? blCtx.ALFA : alphaRad;

  const gamA = viscous
    ? ctxPanel.GAM_A
    : (qinvA ? qinvA.subarray(1, nb + 1) : null);

  const clRes = clcalc(
    nb,
    ctxPanel.X,
    ctxPanel.Y,
    ctxPanel.GAM,
    gamA,
    alfa,
    minf,
    qinf,
  );

  let cdRes = { cd: 0.0, cdf: 0.0 };
  if (viscous && blCtx) {
    cdRes = cdcalc(ctxPanel, blCtx, alfa, qinf);
  }

  return {
    cl: clRes.cl,
    cm: clRes.cm,
    cd: cdRes.cd,
    cdf: cdRes.cdf,
    cdp: clRes.cdp,
    clAlf: clRes.clAlf,
    clMsq: clRes.clMsq,
  };
}

export {
  PANGEN_DEFAULTS,
  cpcalc,
  clcalc,
  cdcalc,
  tecalc,
  pangen,
  computeCoefficients,
};
