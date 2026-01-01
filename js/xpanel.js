// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

// Port of src/xpanel.f from XFOIL (Mark Drela). Numerical formulas preserved.
import { segspl } from './spline.js';
import { atanc, setexp } from './xutils.js?v=1';
import { ludcmp, baksub } from './xsolve.js';
import { iblsys } from './xbl.js';

// Panel-method core: geometry processing, wake setup, and influence calculations.

// Sign helper used in corner/angle logic where XFOIL expects +/-1.
function sign1(value) {
  return value >= 0.0 ? 1.0 : -1.0;
}

// Compute panel angles with trailing-edge handling (APCALC).
function apcalc(ctx) {
  const { N: n, X: x, Y: y, NX: nx, NY: ny, APANEL: apanel, SHARP: sharp, PI: pi, ANTE: ante, ASTE: aste, DSTE: dste } = ctx;

  for (let i = 0; i < n - 1; i += 1) {
    const sx = x[i + 1] - x[i];
    const sy = y[i + 1] - y[i];
    if (sx === 0.0 && sy === 0.0) {
      apanel[i] = Math.atan2(-ny[i], -nx[i]);
    } else {
      apanel[i] = Math.atan2(sx, -sy);
    }
  }

  const i = n - 1;
  const ip = 0;
  if (sharp) {
    apanel[i] = pi;
  } else {
    const sx = x[ip] - x[i];
    const sy = y[ip] - y[i];
    apanel[i] = Math.atan2(-sx, sy) + pi;
  }
}

// Generate wake panels aligned with local streamline direction (XYWAKE).
function xywake(ctx) {
  const { N: n, X: x, Y: y, XP: xp, YP: yp, S: s, NX: nx, NY: ny, APANEL: apanel } = ctx;
  const waklen = ctx.WAKLEN ?? 1.0;
  const nw = ctx.NW ?? 0;
  if (nw <= 0) {
    ctx.NW = 0;
    return;
  }

  // First wake point spacing from average of TE-adjacent panel lengths.
  const ds1 = 0.5 * ((s[1] - s[0]) + (s[n - 1] - s[n - 2]));
  let chord = ctx.CHORD;
  if (!Number.isFinite(chord) || chord <= 0.0) {
    let xmin = Infinity;
    let xmax = -Infinity;
    for (let i = 0; i < n; i += 1) {
      xmin = Math.min(xmin, x[i]);
      xmax = Math.max(xmax, x[i]);
    }
    chord = xmax - xmin;
  }
  const snew = ctx.SNEW ?? new Float64Array(n + nw);
  if (!ctx.SNEW || ctx.SNEW.length < n + nw) {
    ctx.SNEW = snew;
  }
  const wakeS = new Float64Array(nw);
  setexp(wakeS, ds1, waklen * chord, nw);
  for (let i = 0; i < nw; i += 1) {
    snew[n + i] = wakeS[i];
  }

  const xte = 0.5 * (x[0] + x[n - 1]);
  const yte = 0.5 * (y[0] + y[n - 1]);
  ctx.XTE = xte;
  ctx.YTE = yte;

  // Set first wake point a tiny distance behind TE.
  const i0 = n;
  const sx = 0.5 * (yp[n - 1] - yp[0]);
  const sy = 0.5 * (xp[0] - xp[n - 1]);
  const smod = Math.sqrt(sx * sx + sy * sy) || 1.0;
  nx[i0] = sx / smod;
  ny[i0] = sy / smod;
  x[i0] = xte - 0.0001 * ny[i0];
  y[i0] = yte + 0.0001 * nx[i0];
  s[i0] = s[n - 1];

  // Streamfunction gradient components at the first wake point.
  const psiX = psilin(i0, x[i0], y[i0], 1.0, 0.0, false, false, ctx).psiNi;
  const psiY = psilin(i0, x[i0], y[i0], 0.0, 1.0, false, false, ctx).psiNi;
  const denom = Math.sqrt(psiX * psiX + psiY * psiY) || 1.0;
  nx[i0 + 1] = -psiX / denom;
  ny[i0 + 1] = -psiY / denom;
  apanel[i0] = Math.atan2(psiY, psiX);

  // Set rest of wake points downstream.
  for (let iw = 1; iw < nw; iw += 1) {
    const i = n + iw;
    const ds = snew[i] - snew[i - 1];
    x[i] = x[i - 1] - ds * ny[i];
    y[i] = y[i - 1] + ds * nx[i];
    s[i] = s[i - 1] + ds;
    if (i === n + nw - 1) break;

    // Normal vector for next wake point.
    const psiXw = psilin(i, x[i], y[i], 1.0, 0.0, false, false, ctx).psiNi;
    const psiYw = psilin(i, x[i], y[i], 0.0, 1.0, false, false, ctx).psiNi;
    const denomw = Math.sqrt(psiXw * psiXw + psiYw * psiYw) || 1.0;
    nx[i + 1] = -psiXw / denomw;
    ny[i + 1] = -psiYw / denomw;
    apanel[i] = Math.atan2(psiYw, psiXw);
  }

  ctx.LWAKE = true;
  ctx.AWAKE = ctx.ALFA;
  ctx.LWDIJ = false;
}

// Wake strength initialization based on Kutta condition and geometry (QWCALC).
function qwcalc(ctx) {
  const { N: n, NW: nw = 0, QINVU: qinvu, X: x, Y: y, NX: nx, NY: ny } = ctx;
  // First wake point (same as TE).
  if (nw <= 0) return;
  qinvu[n][0] = qinvu[n - 1][0];
  qinvu[n][1] = qinvu[n - 1][1];
  // Rest of wake.
  for (let i = n + 1; i < n + nw; i += 1) {
    const res = psilin(i, x[i], y[i], nx[i], ny[i], false, false, ctx);
    qinvu[i][0] = res.qtan1;
    qinvu[i][1] = res.qtan2;
  }
}

// Compute unit normals along the surface from spline parameterization (NCALC).
function ncalc(x, y, s, n, xn, yn) {
  if (n <= 1) {
    return;
  }

  segspl(x, xn, s, n);
  segspl(y, yn, s, n);

  for (let i = 0; i < n; i += 1) {
    const sx = yn[i];
    const sy = -xn[i];
    const smod = Math.sqrt(sx * sx + sy * sy);
    if (smod === 0.0) {
      xn[i] = -1.0;
      yn[i] = 0.0;
    } else {
      xn[i] = sx / smod;
      yn[i] = sy / smod;
    }
  }

  for (let i = 0; i < n - 1; i += 1) {
    if (s[i] === s[i + 1]) {
      const sx = 0.5 * (xn[i] + xn[i + 1]);
      const sy = 0.5 * (yn[i] + yn[i + 1]);
      const smod = Math.sqrt(sx * sx + sy * sy);
      if (smod === 0.0) {
        xn[i] = -1.0;
        yn[i] = 0.0;
        xn[i + 1] = -1.0;
        yn[i + 1] = 0.0;
      } else {
        xn[i] = sx / smod;
        yn[i] = sy / smod;
        xn[i + 1] = sx / smod;
        yn[i + 1] = sy / smod;
      }
    }
  }
}

// Assemble the panel influence matrix and RHS for gamma solution (GGCALC).
function ggcalc(ctx) {
  const {
    N: n,
    X: x,
    Y: y,
    XP: xp,
    YP: yp,
    NX: nx,
    NY: ny,
    APANEL: apanel,
    SHARP: sharp,
    PI: pi,
    GAM: gam,
    GAMU: gamu,
    SIG: sig,
    QINF: qinf,
    ALFA: alfa,
    AIJ: aij,
    BIJ: bij,
    AIJPIV: aijpiv,
    QINVU: qinvu,
    XTE: xte,
    YTE: yte,
  } = ctx;

  const bwt = 0.1;

  for (let i = 0; i < n; i += 1) {
    gam[i] = 0.0;
    gamu[i][0] = 0.0;
    gamu[i][1] = 0.0;
  }

  // Set up influence matrix rows for each control point.
  for (let i = 0; i < n; i += 1) {
    const { psi } = psilin(i, x[i], y[i], nx[i], ny[i], false, true, ctx);

    const res1 = qinf * y[i];
    const res2 = -qinf * x[i];

    for (let j = 0; j < n; j += 1) {
      aij[i][j] = ctx.DZDG[j];
      bij[i][j] = -ctx.DZDM[j];
    }

    aij[i][n] = -1.0;
    gamu[i][0] = -res1;
    gamu[i][1] = -res2;
  }

  // Kutta condition row.
  for (let j = 0; j < n + 1; j += 1) {
    aij[n][j] = 0.0;
  }
  aij[n][0] = 1.0;
  aij[n][n - 1] = 1.0;
  gamu[n][0] = 0.0;
  gamu[n][1] = 0.0;

  // Multiply each dPsi/dSig vector by inverse of dPsi/dGam matrix.
  for (let j = 0; j < n; j += 1) {
    bij[n][j] = 0.0;
  }

  // Sharp TE gamma extrapolation row.
  if (sharp) {
    const ag1 = Math.atan2(-yp[0], -xp[0]);
    const ag2 = atanc(yp[n - 1], xp[n - 1], ag1);
    const abis = 0.5 * (ag1 + ag2);
    const cbis = Math.cos(abis);
    const sbis = Math.sin(abis);

    const ds1 = Math.sqrt((x[0] - x[1]) ** 2 + (y[0] - y[1]) ** 2);
    const ds2 = Math.sqrt((x[n - 1] - x[n - 2]) ** 2 + (y[n - 1] - y[n - 2]) ** 2);
    const dsmin = Math.min(ds1, ds2);

    const xbis = xte - bwt * dsmin * cbis;
    const ybis = yte - bwt * dsmin * sbis;

    psilin(-1, xbis, ybis, -sbis, cbis, false, true, ctx);

    for (let j = 0; j < n; j += 1) {
      aij[n - 1][j] = ctx.DQDG[j];
      bij[n - 1][j] = -ctx.DQDM[j];
    }

    aij[n - 1][n] = 0.0;
    gamu[n - 1][0] = -cbis;
    gamu[n - 1][1] = -sbis;
  }

  // Factor matrix and solve for alpha=0 and alpha=90 distributions.
  ludcmp(n + 1, aij, aijpiv);

  const rhs1 = new Float64Array(n + 1);
  const rhs2 = new Float64Array(n + 1);
  for (let i = 0; i < n + 1; i += 1) {
    rhs1[i] = gamu[i][0];
    rhs2[i] = gamu[i][1];
  }

  baksub(n + 1, aij, aijpiv, rhs1);
  baksub(n + 1, aij, aijpiv, rhs2);

  for (let i = 0; i < n + 1; i += 1) {
    gamu[i][0] = rhs1[i];
    gamu[i][1] = rhs2[i];
  }

  // Superimpose alpha=0 and alpha=90 distributions.
  const cosa = Math.cos(alfa);
  const sina = Math.sin(alfa);
  for (let i = 0; i < n; i += 1) {
    qinvu[i][0] = gamu[i][0];
    qinvu[i][1] = gamu[i][1];
    gam[i] = cosa * gamu[i][0] + sina * gamu[i][1];
  }
}

// Streamfunction/velocity influence at (xi, yi) from a straight panel (PSILIN).
// Derived from classic source/vortex panel integrals (Drela formulation).
function psilin(i, xi, yi, nxi, nyi, geolin, siglin, ctx) {
  // Fortran comments (PSILIN) highlight:
  // - Skip zero-length TE panel.
  // - Use sign of yy for branch cuts in atan/log.
  // - Apply TE source/vortex correction for sharp/round TE.
  // - Accumulate streamfunction, normal derivative, and tangent velocities.
  const {
    N: n,
    X: x,
    Y: y,
    S: s,
    NX: nx,
    NY: ny,
    APANEL: apanel,
    SHARP: sharp,
    PI: pi,
    ANTE: ante,
    ASTE: aste,
    DSTE: dste,
    GAM: gam,
    GAMU: gamu,
    SIG: sig,
    QF0: qf0,
    QF1: qf1,
    QF2: qf2,
    QF3: qf3,
    QOPI: qopi,
    HOPI: hopi,
    ALFA: alfa,
    QINF: qinf,
    LIMAGE: limage,
    YIMAGE: yimage,
    DZDG: dzdg,
    DZDN: dzdn,
    DQDG: dqdg,
    DZDM: dzdm,
    DQDM: dqdm,
  } = ctx;

  const seps = (s[n - 1] - s[0]) * 1.0e-5;
  const io = i;
  const cosa = Math.cos(alfa);
  const sina = Math.sin(alfa);

  for (let jo = 0; jo < n; jo += 1) {
    dzdg[jo] = 0.0;
    dzdn[jo] = 0.0;
    dqdg[jo] = 0.0;
  }
  let skipTeImage = false;
  for (let jo = 0; jo < n; jo += 1) {
    dzdm[jo] = 0.0;
    dqdm[jo] = 0.0;
  }

  let zQinf = 0.0;
  let zAlfa = 0.0;
  let zQdof0 = 0.0;
  let zQdof1 = 0.0;
  let zQdof2 = 0.0;
  let zQdof3 = 0.0;

  let psi = 0.0;
  let psiNi = 0.0;

  let qtan1 = 0.0;
  let qtan2 = 0.0;
  let qtanm = 0.0;

  let scs = 1.0;
  let sds = 0.0;
  if (!sharp) {
    scs = ante / dste;
    sds = aste / dste;
  }

  let x1 = 0.0;
  let x2 = 0.0;
  let yy = 0.0;
  let g1 = 0.0;
  let g2 = 0.0;
  let t1 = 0.0;
  let t2 = 0.0;
  let apan = 0.0;
  let x1i = 0.0;
  let x2i = 0.0;
  let yyi = 0.0;
  let x1o = 0.0;
  let x1p = 0.0;
  let x2o = 0.0;
  let x2p = 0.0;
  let yyo = 0.0;
  let yyp = 0.0;
  let teJo = 0;
  let teJp = 0;
  let skipTe = false;

  for (let jo = 0; jo < n; jo += 1) {
    let jp = jo + 1;
    let jm = jo - 1;
    let jq = jp + 1;

    if (jo === 0) {
      jm = jo;
    } else if (jo === n - 2) {
      jq = jp;
    } else if (jo === n - 1) {
      jp = 0;
      const dxTe = x[jo] - x[jp];
      const dyTe = y[jo] - y[jp];
      if ((dxTe * dxTe + dyTe * dyTe) < seps * seps) {
        skipTe = true;
        break;
      }
    }

    const xjo = x[jo];
    const yjo = y[jo];
    const xjp = x[jp];
    const yjp = y[jp];
    const dxPanel = xjo - xjp;
    const dyPanel = yjo - yjp;
    const dso = Math.sqrt(dxPanel * dxPanel + dyPanel * dyPanel);
    if (dso === 0.0) {
      continue;
    }

    const dsio = 1.0 / dso;
    apan = apanel[jo];

    const rx1 = xi - xjo;
    const ry1 = yi - yjo;
    const rx2 = xi - xjp;
    const ry2 = yi - yjp;

    const sx = (xjp - xjo) * dsio;
    const sy = (yjp - yjo) * dsio;

    x1 = sx * rx1 + sy * ry1;
    x2 = sx * rx2 + sy * ry2;
    yy = sx * ry1 - sy * rx1;

    const rs1 = rx1 * rx1 + ry1 * ry1;
    const rs2 = rx2 * rx2 + ry2 * ry2;

    const sgn = (io >= 0 && io <= n - 1) ? 1.0 : sign1(yy);

    // Guard self-influence for logarithms.
    if (io !== jo && rs1 > 0.0) {
      g1 = Math.log(rs1);
      t1 = Math.atan2(sgn * x1, sgn * yy) + (0.5 - 0.5 * sgn) * pi;
    } else {
      g1 = 0.0;
      t1 = 0.0;
    }

    if (io !== jp && rs2 > 0.0) {
      g2 = Math.log(rs2);
      t2 = Math.atan2(sgn * x2, sgn * yy) + (0.5 - 0.5 * sgn) * pi;
    } else {
      g2 = 0.0;
      t2 = 0.0;
    }

    x1i = sx * nxi + sy * nyi;
    x2i = sx * nxi + sy * nyi;
    yyi = sx * nyi - sy * nxi;

    // Geometry derivative terms for linearized influence coefficients.
    if (geolin) {
      const nxo = nx[jo];
      const nyo = ny[jo];
      const nxp = nx[jp];
      const nyp = ny[jp];

      x1o = -((rx1 - x1 * sx) * nxo + (ry1 - x1 * sy) * nyo) * dsio - (sx * nxo + sy * nyo);
      x1p = ((rx1 - x1 * sx) * nxp + (ry1 - x1 * sy) * nyp) * dsio;
      x2o = -((rx2 - x2 * sx) * nxo + (ry2 - x2 * sy) * nyo) * dsio;
      x2p = ((rx2 - x2 * sx) * nxp + (ry2 - x2 * sy) * nyp) * dsio - (sx * nxp + sy * nyp);
      yyo = ((rx1 + x1 * sy) * nyo - (ry1 - x1 * sx) * nxo) * dsio - (sx * nyo - sy * nxo);
      yyp = -((rx1 - x1 * sy) * nyp - (ry1 + x1 * sx) * nxp) * dsio;
    }

    // Trailing edge panel: store indices for TE correction.
    if (jo === n - 1) {
      teJo = jo;
      teJp = jp;
      break;
    }

    // Linear source variation on panels.
    if (siglin) {
      const x0 = 0.5 * (x1 + x2);
      const rs0 = x0 * x0 + yy * yy;
      const g0 = Math.log(rs0);
      const t0 = Math.atan2(sgn * x0, sgn * yy) + (0.5 - 0.5 * sgn) * pi;

      let dxinv = 1.0 / (x1 - x0);
      let psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
      let pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan)
        + (x0 - x1) * yy) * dxinv;

      let psx1 = -(t1 - apan);
      let psx0 = t0 - apan;
      let psyy = 0.5 * (g1 - g0);

      let pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
      let pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
      let pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;

      const dxm = xjp - x[jm];
      const dym = yjp - y[jm];
      const dsm = Math.sqrt(dxm * dxm + dym * dym);
      const dsim = 1.0 / dsm;

      let ssum = (sig[jp] - sig[jo]) * dsio + (sig[jp] - sig[jm]) * dsim;
      let sdif = (sig[jp] - sig[jo]) * dsio - (sig[jp] - sig[jm]) * dsim;

      // Accumulate streamfunction and its normal derivative.
      psi += qopi * (psum * ssum + pdif * sdif);

      dzdm[jm] += qopi * (-psum * dsim + pdif * dsim);
      dzdm[jo] += qopi * (-psum * dsio - pdif * dsio);
      dzdm[jp] += qopi * (psum * (dsio + dsim) + pdif * (dsio - dsim));

      let psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
      let pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
      psiNi += qopi * (psni * ssum + pdni * sdif);

      qtanm += qopi * (psni * ssum + pdni * sdif);

      dqdm[jm] += qopi * (-psni * dsim + pdni * dsim);
      dqdm[jo] += qopi * (-psni * dsio - pdni * dsio);
      dqdm[jp] += qopi * (psni * (dsio + dsim) + pdni * (dsio - dsim));

      dxinv = 1.0 / (x0 - x2);
      psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
      pdif = ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan)
        + (x2 - x0) * yy) * dxinv;

      psx0 = -(t0 - apan);
      let psx2 = t2 - apan;
      psyy = 0.5 * (g0 - g2);

      pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
      let pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
      pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv;

      const dxp = x[jq] - xjo;
      const dyp = y[jq] - yjo;
      const dsp = Math.sqrt(dxp * dxp + dyp * dyp);
      const dsip = 1.0 / dsp;

      ssum = (sig[jq] - sig[jo]) * dsip + (sig[jp] - sig[jo]) * dsio;
      sdif = (sig[jq] - sig[jo]) * dsip - (sig[jp] - sig[jo]) * dsio;

      psi += qopi * (psum * ssum + pdif * sdif);

      dzdm[jo] += qopi * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
      dzdm[jp] += qopi * (psum * dsio - pdif * dsio);
      dzdm[jq] += qopi * (psum * dsip + pdif * dsip);

      psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
      pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
      psiNi += qopi * (psni * ssum + pdni * sdif);

      qtanm += qopi * (psni * ssum + pdni * sdif);

      dqdm[jo] += qopi * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
      dqdm[jp] += qopi * (psni * dsio - pdni * dsio);
      dqdm[jq] += qopi * (psni * dsip + pdni * dsip);
    }

    const dxinv = 1.0 / (x1 - x2);
    const psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
    const psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;

    const psx1 = 0.5 * g1;
    const psx2 = -0.5 * g2;
    const psyy = t1 - t2;

    const pdx1 = ((x1 + x2) * psx1 + psis - x1 * g1 - psid) * dxinv;
    const pdx2 = ((x1 + x2) * psx2 + psis + x2 * g2 + psid) * dxinv;
    const pdyy = ((x1 + x2) * psyy - yy * (g1 - g2)) * dxinv;

    const gsum1 = gamu[jp][0] + gamu[jo][0];
    const gsum2 = gamu[jp][1] + gamu[jo][1];
    const gdif1 = gamu[jp][0] - gamu[jo][0];
    const gdif2 = gamu[jp][1] - gamu[jo][1];

    const gsum = gam[jp] + gam[jo];
    const gdif = gam[jp] - gam[jo];

    psi += qopi * (psis * gsum + psid * gdif);

    dzdg[jo] += qopi * (psis - psid);
    dzdg[jp] += qopi * (psis + psid);

    const psni = psx1 * x1i + psx2 * x2i + psyy * yyi;
    const pdni = pdx1 * x1i + pdx2 * x2i + pdyy * yyi;
    psiNi += qopi * (gsum * psni + gdif * pdni);

    qtan1 += qopi * (gsum1 * psni + gdif1 * pdni);
    qtan2 += qopi * (gsum2 * psni + gdif2 * pdni);

    dqdg[jo] += qopi * (psni - pdni);
    dqdg[jp] += qopi * (psni + pdni);

    if (geolin) {
      dzdn[jo] += qopi * gsum * (psx1 * x1o + psx2 * x2o + psyy * yyo)
        + qopi * gdif * (pdx1 * x1o + pdx2 * x2o + pdyy * yyo);
      dzdn[jp] += qopi * gsum * (psx1 * x1p + psx2 * x2p + psyy * yyp)
        + qopi * gdif * (pdx1 * x1p + pdx2 * x2p + pdyy * yyp);

      zQdof0 += qopi * ((psis - psid) * qf0[jo] + (psis + psid) * qf0[jp]);
      zQdof1 += qopi * ((psis - psid) * qf1[jo] + (psis + psid) * qf1[jp]);
      zQdof2 += qopi * ((psis - psid) * qf2[jo] + (psis + psid) * qf2[jp]);
      zQdof3 += qopi * ((psis - psid) * qf3[jo] + (psis + psid) * qf3[jp]);
    }
  }

  if (!skipTe) {
    const psig = 0.5 * yy * (g1 - g2) + x2 * (t2 - apan) - x1 * (t1 - apan);
    const pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);

    const psigx1 = -(t1 - apan);
    const psigx2 = t2 - apan;
    const psigyy = 0.5 * (g1 - g2);
    const pgamx1 = 0.5 * g1;
    const pgamx2 = -0.5 * g2;
    const pgamyy = t1 - t2;

    const psigni = psigx1 * x1i + psigx2 * x2i + psigyy * yyi;
    const pgamni = pgamx1 * x1i + pgamx2 * x2i + pgamyy * yyi;

    const sigte1 = 0.5 * scs * (gamu[teJp][0] - gamu[teJo][0]);
    const sigte2 = 0.5 * scs * (gamu[teJp][1] - gamu[teJo][1]);
    const gamte1 = -0.5 * sds * (gamu[teJp][0] - gamu[teJo][0]);
    const gamte2 = -0.5 * sds * (gamu[teJp][1] - gamu[teJo][1]);

    const sigte = 0.5 * scs * (gam[teJp] - gam[teJo]);
    const gamte = -0.5 * sds * (gam[teJp] - gam[teJo]);

    psi += hopi * (psig * sigte + pgam * gamte);

    dzdg[teJo] -= hopi * psig * scs * 0.5;
    dzdg[teJp] += hopi * psig * scs * 0.5;

    dzdg[teJo] += hopi * pgam * sds * 0.5;
    dzdg[teJp] -= hopi * pgam * sds * 0.5;

    psiNi += hopi * (psigni * sigte + pgamni * gamte);

    qtan1 += hopi * (psigni * sigte1 + pgamni * gamte1);
    qtan2 += hopi * (psigni * sigte2 + pgamni * gamte2);

    dqdg[teJo] -= hopi * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
    dqdg[teJp] += hopi * (psigni * 0.5 * scs - pgamni * 0.5 * sds);

    if (geolin) {
      dzdn[teJo] += hopi * (psigx1 * x1o + psigx2 * x2o + psigyy * yyo) * sigte
        + hopi * (pgamx1 * x1o + pgamx2 * x2o + pgamyy * yyo) * gamte;
      dzdn[teJp] += hopi * (psigx1 * x1p + psigx2 * x2p + psigyy * yyp) * sigte
        + hopi * (pgamx1 * x1p + pgamx2 * x2p + pgamyy * yyp) * gamte;

      zQdof0 += hopi * psig * 0.5 * (qf0[teJp] - qf0[teJo]) * scs
        - hopi * pgam * 0.5 * (qf0[teJp] - qf0[teJo]) * sds;
      zQdof1 += hopi * psig * 0.5 * (qf1[teJp] - qf1[teJo]) * scs
        - hopi * pgam * 0.5 * (qf1[teJp] - qf1[teJo]) * sds;
      zQdof2 += hopi * psig * 0.5 * (qf2[teJp] - qf2[teJo]) * scs
        - hopi * pgam * 0.5 * (qf2[teJp] - qf2[teJo]) * sds;
      zQdof3 += hopi * psig * 0.5 * (qf3[teJp] - qf3[teJo]) * scs
        - hopi * pgam * 0.5 * (qf3[teJp] - qf3[teJo]) * sds;
    }
  }

  psi += qinf * (cosa * yi - sina * xi);

  psiNi += qinf * (cosa * nyi - sina * nxi);

  qtan1 += qinf * nyi;
  qtan2 -= qinf * nxi;

  zQinf += (cosa * yi - sina * xi);
  zAlfa -= qinf * (sina * yi + cosa * xi);

  ctx.QTAN1 = qtan1;
  ctx.QTAN2 = qtan2;
  ctx.QTANM = qtanm;
  ctx.Z_QINF = zQinf;
  ctx.Z_ALFA = zAlfa;
  ctx.Z_QDOF0 = zQdof0;
  ctx.Z_QDOF1 = zQdof1;
  ctx.Z_QDOF2 = zQdof2;
  ctx.Z_QDOF3 = zQdof3;

  if (!limage) {
    return {
      psi,
      psiNi,
      zQinf,
      zAlfa,
      zQdof0,
      zQdof1,
      zQdof2,
      zQdof3,
      qtan1,
      qtan2,
      qtanm,
    };
  }

  for (let jo = 0; jo < n; jo += 1) {
    let jp = jo + 1;
    let jm = jo - 1;
    let jq = jp + 1;

    if (jo === 0) {
      jm = jo;
    } else if (jo === n - 2) {
      jq = jp;
    } else if (jo === n - 1) {
      jp = 0;
      const dxTe = x[jo] - x[jp];
      const dyTe = y[jo] - y[jp];
      if ((dxTe * dxTe + dyTe * dyTe) < seps * seps) {
        skipTeImage = true;
        break;
      }
    }

    const xjo0 = x[jo];
    const yjo0 = y[jo];
    const xjp0 = x[jp];
    const yjp0 = y[jp];
    const dxPanel = xjo0 - xjp0;
    const dyPanel = yjo0 - yjp0;
    const dso = Math.sqrt(dxPanel * dxPanel + dyPanel * dyPanel);
    if (dso === 0.0) {
      continue;
    }

    const dsio = 1.0 / dso;
    apan = pi - apanel[jo] + 2.0 * alfa;

    const xjo = xjo0 + 2.0 * (yimage + yjo0) * sina;
    const yjo = yjo0 - 2.0 * (yimage + yjo0) * cosa;
    const xjp = xjp0 + 2.0 * (yimage + yjp0) * sina;
    const yjp = yjp0 - 2.0 * (yimage + yjp0) * cosa;

    const rx1 = xi - xjo;
    const ry1 = yi - yjo;
    const rx2 = xi - xjp;
    const ry2 = yi - yjp;

    const sx = (xjp - xjo) * dsio;
    const sy = (yjp - yjo) * dsio;

    x1 = sx * rx1 + sy * ry1;
    x2 = sx * rx2 + sy * ry2;
    yy = sx * ry1 - sy * rx1;

    const rs1 = rx1 * rx1 + ry1 * ry1;
    const rs2 = rx2 * rx2 + ry2 * ry2;

    const sgn = (io >= 0 && io <= n - 1) ? 1.0 : sign1(yy);

    g1 = Math.log(rs1);
    t1 = Math.atan2(sgn * x1, sgn * yy) + (0.5 - 0.5 * sgn) * pi;

    g2 = Math.log(rs2);
    t2 = Math.atan2(sgn * x2, sgn * yy) + (0.5 - 0.5 * sgn) * pi;

    x1i = sx * nxi + sy * nyi;
    x2i = sx * nxi + sy * nyi;
    yyi = sx * nyi - sy * nxi;

    if (geolin) {
      const nxo = nx[jo];
      const nyo = ny[jo];
      const nxp = nx[jp];
      const nyp = ny[jp];

      x1o = -((rx1 - x1 * sx) * nxo + (ry1 - x1 * sy) * nyo) * dsio - (sx * nxo + sy * nyo);
      x1p = ((rx1 - x1 * sx) * nxp + (ry1 - x1 * sy) * nyp) * dsio;
      x2o = -((rx2 - x2 * sx) * nxo + (ry2 - x2 * sy) * nyo) * dsio;
      x2p = ((rx2 - x2 * sx) * nxp + (ry2 - x2 * sy) * nyp) * dsio - (sx * nxp + sy * nyp);
      yyo = ((rx1 + x1 * sy) * nyo - (ry1 - x1 * sx) * nxo) * dsio - (sx * nyo - sy * nxo);
      yyp = -((rx1 - x1 * sy) * nyp - (ry1 + x1 * sx) * nxp) * dsio;
    }

    if (jo === n - 1) {
      teJo = jo;
      teJp = jp;
      break;
    }

    if (siglin) {
      const x0 = 0.5 * (x1 + x2);
      const rs0 = x0 * x0 + yy * yy;
      const g0 = Math.log(rs0);
      const t0 = Math.atan2(sgn * x0, sgn * yy) + (0.5 - 0.5 * sgn) * pi;

      let dxinv = 1.0 / (x1 - x0);
      let psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
      let pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan)
        + (x0 - x1) * yy) * dxinv;

      let psx1 = -(t1 - apan);
      let psx0 = t0 - apan;
      let psyy = 0.5 * (g1 - g0);

      let pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
      let pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
      let pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;

      const dxm = x[jp] - x[jm];
      const dym = y[jp] - y[jm];
      const dsm = Math.sqrt(dxm * dxm + dym * dym);
      const dsim = 1.0 / dsm;

      let ssum = (sig[jp] - sig[jo]) * dsio + (sig[jp] - sig[jm]) * dsim;
      let sdif = (sig[jp] - sig[jo]) * dsio - (sig[jp] - sig[jm]) * dsim;

      psi += qopi * (psum * ssum + pdif * sdif);

      dzdm[jm] += qopi * (-psum * dsim + pdif * dsim);
      dzdm[jo] += qopi * (-psum * dsio - pdif * dsio);
      dzdm[jp] += qopi * (psum * (dsio + dsim) + pdif * (dsio - dsim));

      let psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
      let pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
      psiNi += qopi * (psni * ssum + pdni * sdif);

      qtanm += qopi * (psni * ssum + pdni * sdif);

      dqdm[jm] += qopi * (-psni * dsim + pdni * dsim);
      dqdm[jo] += qopi * (-psni * dsio - pdni * dsio);
      dqdm[jp] += qopi * (psni * (dsio + dsim) + pdni * (dsio - dsim));

      dxinv = 1.0 / (x0 - x2);
      psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
      pdif = ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan)
        + (x2 - x0) * yy) * dxinv;

      psx0 = -(t0 - apan);
      let psx2 = t2 - apan;
      psyy = 0.5 * (g0 - g2);

      pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
      let pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
      pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv;

      const dxp = x[jq] - x[jo];
      const dyp = y[jq] - y[jo];
      const dsp = Math.sqrt(dxp * dxp + dyp * dyp);
      const dsip = 1.0 / dsp;

      ssum = (sig[jq] - sig[jo]) * dsip + (sig[jp] - sig[jo]) * dsio;
      sdif = (sig[jq] - sig[jo]) * dsip - (sig[jp] - sig[jo]) * dsio;

      psi += qopi * (psum * ssum + pdif * sdif);

      dzdm[jo] += qopi * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
      dzdm[jp] += qopi * (psum * dsio - pdif * dsio);
      dzdm[jq] += qopi * (psum * dsip + pdif * dsip);

      psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
      pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
      psiNi += qopi * (psni * ssum + pdni * sdif);

      qtanm += qopi * (psni * ssum + pdni * sdif);

      dqdm[jo] += qopi * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
      dqdm[jp] += qopi * (psni * dsio - pdni * dsio);
      dqdm[jq] += qopi * (psni * dsip + pdni * dsip);
    }

    const dxinv = 1.0 / (x1 - x2);
    const psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
    const psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;

    const psx1 = 0.5 * g1;
    const psx2 = -0.5 * g2;
    const psyy = t1 - t2;

    const pdx1 = ((x1 + x2) * psx1 + psis - x1 * g1 - psid) * dxinv;
    const pdx2 = ((x1 + x2) * psx2 + psis + x2 * g2 + psid) * dxinv;
    const pdyy = ((x1 + x2) * psyy - yy * (g1 - g2)) * dxinv;

    const gsum1 = gamu[jp][0] + gamu[jo][0];
    const gsum2 = gamu[jp][1] + gamu[jo][1];
    const gdif1 = gamu[jp][0] - gamu[jo][0];
    const gdif2 = gamu[jp][1] - gamu[jo][1];

    const gsum = gam[jp] + gam[jo];
    const gdif = gam[jp] - gam[jo];

    psi -= qopi * (psis * gsum + psid * gdif);

    dzdg[jo] -= qopi * (psis - psid);
    dzdg[jp] -= qopi * (psis + psid);

    const psni = psx1 * x1i + psx2 * x2i + psyy * yyi;
    const pdni = pdx1 * x1i + pdx2 * x2i + pdyy * yyi;
    psiNi -= qopi * (gsum * psni + gdif * pdni);

    qtan1 -= qopi * (gsum1 * psni + gdif1 * pdni);
    qtan2 -= qopi * (gsum2 * psni + gdif2 * pdni);

    dqdg[jo] -= qopi * (psni - pdni);
    dqdg[jp] -= qopi * (psni + pdni);

    if (geolin) {
      dzdn[jo] -= qopi * gsum * (psx1 * x1o + psx2 * x2o + psyy * yyo)
        + qopi * gdif * (pdx1 * x1o + pdx2 * x2o + pdyy * yyo);
      dzdn[jp] -= qopi * gsum * (psx1 * x1p + psx2 * x2p + psyy * yyp)
        + qopi * gdif * (pdx1 * x1p + pdx2 * x2p + pdyy * yyp);

      zQdof0 -= qopi * ((psis - psid) * qf0[jo] + (psis + psid) * qf0[jp]);
      zQdof1 -= qopi * ((psis - psid) * qf1[jo] + (psis + psid) * qf1[jp]);
      zQdof2 -= qopi * ((psis - psid) * qf2[jo] + (psis + psid) * qf2[jp]);
      zQdof3 -= qopi * ((psis - psid) * qf3[jo] + (psis + psid) * qf3[jp]);
    }
  }

  if (!skipTeImage) {
    const psig = 0.5 * yy * (g1 - g2) + x2 * (t2 - apan) - x1 * (t1 - apan);
    const pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);

    const psigx1 = -(t1 - apan);
    const psigx2 = t2 - apan;
    const psigyy = 0.5 * (g1 - g2);
    const pgamx1 = 0.5 * g1;
    const pgamx2 = -0.5 * g2;
    const pgamyy = t1 - t2;

    const psigni = psigx1 * x1i + psigx2 * x2i + psigyy * yyi;
    const pgamni = pgamx1 * x1i + pgamx2 * x2i + pgamyy * yyi;

    const sigte1 = 0.5 * scs * (gamu[teJp][0] - gamu[teJo][0]);
    const sigte2 = 0.5 * scs * (gamu[teJp][1] - gamu[teJo][1]);
    const gamte1 = -0.5 * sds * (gamu[teJp][0] - gamu[teJo][0]);
    const gamte2 = -0.5 * sds * (gamu[teJp][1] - gamu[teJo][1]);

    const sigte = 0.5 * scs * (gam[teJp] - gam[teJo]);
    const gamte = -0.5 * sds * (gam[teJp] - gam[teJo]);

    psi += hopi * (psig * sigte - pgam * gamte);

    dzdg[teJo] -= hopi * psig * scs * 0.5;
    dzdg[teJp] += hopi * psig * scs * 0.5;

    dzdg[teJo] -= hopi * pgam * sds * 0.5;
    dzdg[teJp] += hopi * pgam * sds * 0.5;

    psiNi += hopi * (psigni * sigte - pgamni * gamte);

    qtan1 += hopi * (psigni * sigte1 - pgamni * gamte1);
    qtan2 += hopi * (psigni * sigte2 - pgamni * gamte2);

    dqdg[teJo] -= hopi * (psigni * 0.5 * scs + pgamni * 0.5 * sds);
    dqdg[teJp] += hopi * (psigni * 0.5 * scs + pgamni * 0.5 * sds);

    if (geolin) {
      dzdn[teJo] += hopi * (psigx1 * x1o + psigx2 * x2o + psigyy * yyo) * sigte
        - hopi * (pgamx1 * x1o + pgamx2 * x2o + pgamyy * yyo) * gamte;
      dzdn[teJp] += hopi * (psigx1 * x1p + psigx2 * x2p + psigyy * yyp) * sigte
        - hopi * (pgamx1 * x1p + pgamx2 * x2p + pgamyy * yyp) * gamte;

      zQdof0 += hopi * psig * 0.5 * (qf0[teJp] - qf0[teJo]) * scs
        + hopi * pgam * 0.5 * (qf0[teJp] - qf0[teJo]) * sds;
      zQdof1 += hopi * psig * 0.5 * (qf1[teJp] - qf1[teJo]) * scs
        + hopi * pgam * 0.5 * (qf1[teJp] - qf1[teJo]) * sds;
      zQdof2 += hopi * psig * 0.5 * (qf2[teJp] - qf2[teJo]) * scs
        + hopi * pgam * 0.5 * (qf2[teJp] - qf2[teJo]) * sds;
      zQdof3 += hopi * psig * 0.5 * (qf3[teJp] - qf3[teJo]) * scs
        + hopi * pgam * 0.5 * (qf3[teJp] - qf3[teJo]) * sds;
    }
  }

  return {
    psi,
    psiNi,
    zQinf,
    zAlfa,
    zQdof0,
    zQdof1,
    zQdof2,
    zQdof3,
    qtan1,
    qtan2,
    qtanm,
  };
}

// Wake-panel influence (PSWLIN), used for Kutta and wake alignment.
function pswlin(i, xi, yi, nxi, nyi, ctx) {
  const {
    N: n,
    NW: nw = 0,
    X: x,
    Y: y,
    NX: nx,
    NY: ny,
    APANEL: apanel,
    SIG: sig,
    QOPI: qopi,
    PI: pi,
    DZDM: dzdm,
    DQDM: dqdm,
  } = ctx;

  if (nw <= 0) {
    return { psi: 0.0, psiNi: 0.0 };
  }

  const io = i;

  for (let jo = n; jo < n + nw; jo += 1) {
    dzdm[jo] = 0.0;
    dqdm[jo] = 0.0;
  }

  let psi = 0.0;
  let psiNi = 0.0;

  for (let jo = n; jo < n + nw - 1; jo += 1) {
    const jp = jo + 1;
    let jm = jo - 1;
    let jq = jp + 1;
    if (jo === n) {
      jm = jo;
    } else if (jo === n + nw - 2) {
      jq = jp;
    }

    const dso = Math.sqrt((x[jo] - x[jp]) ** 2 + (y[jo] - y[jp]) ** 2);
    const dsio = 1.0 / dso;

    const apan = apanel[jo];

    const rx1 = xi - x[jo];
    const ry1 = yi - y[jo];
    const rx2 = xi - x[jp];
    const ry2 = yi - y[jp];

    const sx = (x[jp] - x[jo]) * dsio;
    const sy = (y[jp] - y[jo]) * dsio;

    const x1 = sx * rx1 + sy * ry1;
    const x2 = sx * rx2 + sy * ry2;
    const yy = sx * ry1 - sy * rx1;

    const rs1 = rx1 * rx1 + ry1 * ry1;
    const rs2 = rx2 * rx2 + ry2 * ry2;

    const sgn = (io >= n && io <= n + nw - 1) ? 1.0 : sign1(yy);

    let g1 = 0.0;
    let t1 = 0.0;
    if (io !== jo && rs1 > 0.0) {
      g1 = Math.log(rs1);
      t1 = Math.atan2(sgn * x1, sgn * yy) - (0.5 - 0.5 * sgn) * pi;
    }

    let g2 = 0.0;
    let t2 = 0.0;
    if (io !== jp && rs2 > 0.0) {
      g2 = Math.log(rs2);
      t2 = Math.atan2(sgn * x2, sgn * yy) - (0.5 - 0.5 * sgn) * pi;
    }

    const x1i = sx * nxi + sy * nyi;
    const x2i = sx * nxi + sy * nyi;
    const yyi = sx * nyi - sy * nxi;

    const x0 = 0.5 * (x1 + x2);
    const rs0 = x0 * x0 + yy * yy;
    const g0 = Math.log(rs0);
    const t0 = Math.atan2(sgn * x0, sgn * yy) - (0.5 - 0.5 * sgn) * pi;

    let dxinv = 1.0 / (x1 - x0);
    let psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
    let pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan)
      + (x0 - x1) * yy) * dxinv;

    let psx1 = -(t1 - apan);
    let psx0 = t0 - apan;
    let psyy = 0.5 * (g1 - g0);

    let pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
    let pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
    let pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;

    const dsm = Math.sqrt((x[jp] - x[jm]) ** 2 + (y[jp] - y[jm]) ** 2);
    const dsim = 1.0 / dsm;

    let ssum = (sig[jp] - sig[jo]) * dsio + (sig[jp] - sig[jm]) * dsim;
    let sdif = (sig[jp] - sig[jo]) * dsio - (sig[jp] - sig[jm]) * dsim;

    psi += qopi * (psum * ssum + pdif * sdif);

    dzdm[jm] += qopi * (-psum * dsim + pdif * dsim);
    dzdm[jo] += qopi * (-psum * dsio - pdif * dsio);
    dzdm[jp] += qopi * (psum * (dsio + dsim) + pdif * (dsio - dsim));

    let psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
    let pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
    psiNi += qopi * (psni * ssum + pdni * sdif);

    dqdm[jm] += qopi * (-psni * dsim + pdni * dsim);
    dqdm[jo] += qopi * (-psni * dsio - pdni * dsio);
    dqdm[jp] += qopi * (psni * (dsio + dsim) + pdni * (dsio - dsim));

    dxinv = 1.0 / (x0 - x2);
    psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
    pdif = ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan)
      + (x2 - x0) * yy) * dxinv;

    psx0 = -(t0 - apan);
    let psx2 = t2 - apan;
    psyy = 0.5 * (g0 - g2);

    pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
    let pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
    pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv;

    const dsp = Math.sqrt((x[jq] - x[jo]) ** 2 + (y[jq] - y[jo]) ** 2);
    const dsip = 1.0 / dsp;

    ssum = (sig[jq] - sig[jo]) * dsip + (sig[jp] - sig[jo]) * dsio;
    sdif = (sig[jq] - sig[jo]) * dsip - (sig[jp] - sig[jo]) * dsio;

    psi += qopi * (psum * ssum + pdif * sdif);

    dzdm[jo] += qopi * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
    dzdm[jp] += qopi * (psum * dsio - pdif * dsio);
    dzdm[jq] += qopi * (psum * dsip + pdif * dsip);

    psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
    pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
    psiNi += qopi * (psni * ssum + pdni * sdif);

    dqdm[jo] += qopi * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
    dqdm[jp] += qopi * (psni * dsio - pdni * dsio);
    dqdm[jq] += qopi * (psni * dsip + pdni * dsip);
  }

  return { psi, psiNi };
}

// dQ/dgamma sensitivity for viscous coupling (QDCALC).
function qdcalc(ctx) {
  // Fortran comments (QDCALC) highlight:
  // - Multiply dPsi/dSig by inverse of dPsi/dGam.
  // - No direct source influence for Kutta/TE extrapolation rows.
  // - Wake-point dQtan/dGam and dQtan/dSig contributions.
  const n = ctx.N;
  const nw = ctx.NW ?? 0;
  const total = n + nw;
  const n1 = n + 1;
  const dij = new Array(total + 1);
  for (let i = 0; i <= total; i += 1) {
    dij[i] = new Float64Array(total + 1);
  }

  for (let j = 0; j < n; j += 1) {
    const rhs = new Float64Array(n1);
    for (let i = 0; i < n1; i += 1) {
      rhs[i] = ctx.BIJ[i][j];
    }
    baksub(n1, ctx.AIJ, ctx.AIJPIV, rhs);
    for (let i = 0; i < n; i += 1) {
      dij[i + 1][j + 1] = rhs[i];
    }
  }

  if (nw > 0) {
    // Source influence on airfoil surface.
    for (let i = 0; i < n; i += 1) {
      pswlin(i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i], ctx);
      for (let j = n; j < total; j += 1) {
        ctx.BIJ[i][j] = -ctx.DZDM[j];
      }
    }

    // Kutta condition: no direct source influence.
    for (let j = n; j < total; j += 1) {
      ctx.BIJ[n][j] = 0.0;
    }
    if (ctx.SHARP) {
      const bwt = 0.1;
      const ag1 = Math.atan2(-ctx.YP[0], -ctx.XP[0]);
      const ag2 = atanc(ctx.YP[n - 1], ctx.XP[n - 1], ag1);
      const abis = 0.5 * (ag1 + ag2);
      const cbis = Math.cos(abis);
      const sbis = Math.sin(abis);
      const ds1 = Math.sqrt((ctx.X[0] - ctx.X[1]) ** 2 + (ctx.Y[0] - ctx.Y[1]) ** 2);
      const ds2 = Math.sqrt((ctx.X[n - 1] - ctx.X[n - 2]) ** 2 + (ctx.Y[n - 1] - ctx.Y[n - 2]) ** 2);
      const dsmin = Math.min(ds1, ds2);
      const xbis = ctx.XTE - bwt * dsmin * cbis;
      const ybis = ctx.YTE - bwt * dsmin * sbis;
      pswlin(0, xbis, ybis, -sbis, cbis, ctx);
      for (let j = n; j < total; j += 1) {
        ctx.BIJ[n - 1][j] = -ctx.DQDM[j];
      }
    }

    for (let j = n; j < total; j += 1) {
      const rhs = new Float64Array(n1);
      for (let i = 0; i < n1; i += 1) {
        rhs[i] = ctx.BIJ[i][j];
      }
      baksub(n1, ctx.AIJ, ctx.AIJPIV, rhs);
      for (let i = 0; i < n; i += 1) {
        ctx.BIJ[i][j] = rhs[i];
        dij[i + 1][j + 1] = rhs[i];
      }
    }

    const cij = new Array(nw);
    for (let iw = 0; iw < nw; iw += 1) {
      cij[iw] = new Float64Array(n);
    }

    // Wake contribution: dQtan/dGam and dQtan/dSig at wake points.
    for (let i = n; i < total; i += 1) {
      const iw = i - n;
      psilin(i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i], false, true, ctx);
      for (let j = 0; j < n; j += 1) {
        cij[iw][j] = ctx.DQDG[j];
        dij[i + 1][j + 1] = ctx.DQDM[j];
      }
      pswlin(i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i], ctx);
      for (let j = n; j < total; j += 1) {
        dij[i + 1][j + 1] = ctx.DQDM[j];
      }
    }

    // Add effect of all sources on airfoil vorticity (wake Qtan).
    for (let i = n; i < total; i += 1) {
      const iw = i - n;
      for (let j = 0; j < n; j += 1) {
        let sum = 0.0;
        for (let k = 0; k < n; k += 1) {
          sum += cij[iw][k] * dij[k + 1][j + 1];
        }
        dij[i + 1][j + 1] += sum;
      }
      for (let j = n; j < total; j += 1) {
        let sum = 0.0;
        for (let k = 0; k < n; k += 1) {
          sum += cij[iw][k] * ctx.BIJ[k][j];
        }
        dij[i + 1][j + 1] += sum;
      }
    }

    // First wake point has same velocity as trailing edge.
    for (let j = 1; j <= total; j += 1) {
      dij[n + 1][j] = dij[n][j];
    }
  }

  ctx.DIJ = dij;
  ctx.LADIJ = true;
  ctx.LWDIJ = true;
  return dij;
}

// Stagnation point from circulation sign change; used to split upper/lower (STFIND).
function stfind(ctxPanel, nb) {
  const gam = ctxPanel.GAM;
  // Locate stagnation point where GAM changes sign.
  let ist = Math.floor(nb / 2) - 1;
  let found = false;
  for (let i = 0; i < nb - 1; i += 1) {
    if (gam[i] >= 0.0 && gam[i + 1] < 0.0) {
      ist = i;
      found = true;
      break;
    }
  }

  const s = ctxPanel.S;
  if (!found) {
    console.warn('STFIND: Stagnation point not found. Continuing ...');
    ist = Math.floor(nb / 2) - 1;
  }

  const dgam = gam[ist + 1] - gam[ist];
  const ds = s[ist + 1] - s[ist];
  let sst;
  if (dgam !== 0.0) {
    if (gam[ist] < -gam[ist + 1]) {
      sst = s[ist] - ds * (gam[ist] / dgam);
    } else {
      sst = s[ist + 1] - ds * (gam[ist + 1] / dgam);
    }
  } else {
    sst = 0.5 * (s[ist] + s[ist + 1]);
  }

  if (sst <= s[ist]) sst = s[ist] + 1.0e-7;
  if (sst >= s[ist + 1]) sst = s[ist + 1] - 1.0e-7;

  const sstGo = dgam !== 0.0 ? (sst - s[ist + 1]) / dgam : 0.0;
  const sstGp = dgam !== 0.0 ? (s[ist] - sst) / dgam : 0.0;

  return { ist, sst, sstGo, sstGp };
}

// Map boundary-layer station indices to panel indices (IBLPAN).
function iblpan(blCtx, nb, nw = 0) {
  const ist = blCtx.IST;

  // Top surface first.
  let ibl = 1;
  for (let i = ist; i >= 1; i -= 1) {
    ibl += 1;
    blCtx.IPAN[ibl][1] = i;
    blCtx.VTI[ibl][1] = 1.0;
  }
  blCtx.IBLTE[1] = ibl;
  blCtx.NBL[1] = ibl;

  // Bottom surface next.
  ibl = 1;
  for (let i = ist + 1; i <= nb; i += 1) {
    ibl += 1;
    blCtx.IPAN[ibl][2] = i;
    blCtx.VTI[ibl][2] = -1.0;
  }
  blCtx.IBLTE[2] = ibl;
  // Wake points.
  for (let iw = 1; iw <= nw; iw += 1) {
    const i = nb + iw;
    const iblw = blCtx.IBLTE[2] + iw;
    blCtx.IPAN[iblw][2] = i;
    blCtx.VTI[iblw][2] = -1.0;
  }
  blCtx.NBL[2] = blCtx.IBLTE[2] + nw;

  // Upper wake pointers (for plotting).
  for (let iw = 1; iw <= nw; iw += 1) {
    blCtx.IPAN[blCtx.IBLTE[1] + iw][1] = blCtx.IPAN[blCtx.IBLTE[2] + iw][2];
    blCtx.VTI[blCtx.IBLTE[1] + iw][1] = 1.0;
  }

  blCtx.LIPAN = true;
}

// Surface arc-length mapping for BL stations (XICALC).
function xicalc(blCtx, ctxPanel) {
  const s = ctxPanel.S;
  const nb = ctxPanel.N;
  const nw = ctxPanel.NW ?? 0;
  const x = ctxPanel.X;
  const y = ctxPanel.Y;
  const xp = ctxPanel.XP;
  const yp = ctxPanel.YP;
  const xeps = 1.0e-7 * (s[nb - 1] - s[0]);
  const sAt = (i) => s[i - 1];

  // Minimum xi node arc length near stagnation point.
  blCtx.XSSI[1][1] = 0.0;
  for (let ibl = 2; ibl <= blCtx.IBLTE[1]; ibl += 1) {
    const i = blCtx.IPAN[ibl][1];
    blCtx.XSSI[ibl][1] = Math.max(blCtx.SST - sAt(i), xeps);
  }

  blCtx.XSSI[1][2] = 0.0;
  for (let ibl = 2; ibl <= blCtx.IBLTE[2]; ibl += 1) {
    const i = blCtx.IPAN[ibl][2];
    blCtx.XSSI[ibl][2] = Math.max(sAt(i) - blCtx.SST, xeps);
  }

  if (nw <= 0) {
    return;
  }

  // Wake arc lengths (both sides).
  const ibl1 = blCtx.IBLTE[1] + 1;
  blCtx.XSSI[ibl1][1] = blCtx.XSSI[ibl1 - 1][1];

  const ibl2 = blCtx.IBLTE[2] + 1;
  blCtx.XSSI[ibl2][2] = blCtx.XSSI[ibl2 - 1][2];

  for (let ibl = blCtx.IBLTE[2] + 2; ibl <= blCtx.NBL[2]; ibl += 1) {
    const i = blCtx.IPAN[ibl][2];
    const dxssi = Math.sqrt((x[i - 1] - x[i - 2]) ** 2 + (y[i - 1] - y[i - 2]) ** 2);

    const ibl1w = blCtx.IBLTE[1] + ibl - blCtx.IBLTE[2];
    const ibl2w = blCtx.IBLTE[2] + ibl - blCtx.IBLTE[2];
    blCtx.XSSI[ibl1w][1] = blCtx.XSSI[ibl1w - 1][1] + dxssi;
    blCtx.XSSI[ibl2w][2] = blCtx.XSSI[ibl2w - 1][2] + dxssi;
  }

  const telrat = 2.5;
  const crosp = (xp[0] * yp[nb - 1] - yp[0] * xp[nb - 1])
    / Math.sqrt((xp[0] ** 2 + yp[0] ** 2) * (xp[nb - 1] ** 2 + yp[nb - 1] ** 2));
  let dwdxte = crosp / Math.sqrt(1.0 - crosp ** 2);
  dwdxte = Math.max(dwdxte, -3.0 / telrat);
  dwdxte = Math.min(dwdxte, 3.0 / telrat);

  const aa = 3.0 + telrat * dwdxte;
  const bb = -2.0 - telrat * dwdxte;

  if (ctxPanel.SHARP) {
    for (let iw = 1; iw <= nw; iw += 1) {
      blCtx.WGAP[iw] = 0.0;
    }
  } else {
    for (let iw = 1; iw <= nw; iw += 1) {
      const ibl = blCtx.IBLTE[2] + iw;
      const zn = 1.0 - (blCtx.XSSI[ibl][2] - blCtx.XSSI[blCtx.IBLTE[2]][2]) / (telrat * ctxPanel.ANTE);
      blCtx.WGAP[iw] = 0.0;
      if (zn >= 0.0) {
        blCtx.WGAP[iw] = ctxPanel.ANTE * (aa + bb * zn) * zn ** 2;
      }
    }
  }
}

// Build BL input arrays from inviscid edge velocities (XFOIL UICALC).
// UINV(IBL,IS) = VTI * QINV(I), UINV_A likewise.
function uicalc(blCtx, qinv, qinvA) {
  for (let is = 1; is <= 2; is += 1) {
    blCtx.UINV[1][is] = 0.0;
    blCtx.UINV_A[1][is] = 0.0;
    for (let ibl = 2; ibl <= blCtx.NBL[is]; ibl += 1) {
      const i = blCtx.IPAN[ibl][is];
      blCtx.UINV[ibl][is] = blCtx.VTI[ibl][is] * qinv[i];
      blCtx.UINV_A[ibl][is] = blCtx.VTI[ibl][is] * qinvA[i];
    }
  }
}

// Compute viscous edge speed array from BL variables (XFOIL QVFUE).
// QVIS(I) = VTI * UEDG(IBL,IS).
function qvfue(blCtx, qvis) {
  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= blCtx.NBL[is]; ibl += 1) {
      const i = blCtx.IPAN[ibl][is];
      qvis[i] = blCtx.VTI[ibl][is] * blCtx.UEDG[ibl][is];
    }
  }
}

// Set inviscid panel tangential velocity for current alpha (XFOIL QISET).
// QINV = cos(alpha) * QINVU0 + sin(alpha) * QINVU90.
function qiset(ctxPanel, alphaRad) {
  const total = ctxPanel.N + (ctxPanel.NW ?? 0);
  const qinv = new Float64Array(total + 1);
  const qinvA = new Float64Array(total + 1);
  // Superimpose alpha=0 and alpha=90 distributions.
  const alfa = Number.isFinite(alphaRad) ? alphaRad : (ctxPanel.ALFA ?? 0.0);
  const cosA = Math.cos(alfa);
  const sinA = Math.sin(alfa);
  for (let i = 1; i <= total; i += 1) {
    const q0 = ctxPanel.QINVU[i - 1][0];
    const q90 = ctxPanel.QINVU[i - 1][1];
    const qi = cosA * q0 + sinA * q90;
    qinv[i] = Number.isFinite(qi) ? qi : 0.0;
    const qAi = -sinA * q0 + cosA * q90;
    qinvA[i] = Number.isFinite(qAi) ? qAi : 0.0;
  }
  if (ctxPanel.QINV) {
    ctxPanel.QINV.set(qinv);
  }
  if (ctxPanel.QINV_A) {
    ctxPanel.QINV_A.set(qinvA);
  }
  return { qinv, qinvA };
}

// Update panel circulation from viscous edge velocities (GAMQV).
// GAM(I) = QVIS(I), GAM_A(I) = QINV_A(I).
function gamqv(ctxPanel, qvis, qinvA) {
  const nb = ctxPanel.N;
  for (let i = 1; i <= nb; i += 1) {
    ctxPanel.GAM[i - 1] = qvis[i];
    if (ctxPanel.GAM_A) {
      ctxPanel.GAM_A[i - 1] = qinvA[i];
    }
  }
  if (ctxPanel.QVIS) {
    ctxPanel.QVIS.set(qvis);
  }
}

// Relocate stagnation point and update BL indices (XFOIL STMOVE).
// Shifts BL arrays when IST changes so marching can resume consistently.
function stmove(ctxPanel, blCtx, qinv, qinvA) {
  // Fortran comments (STMOVE) highlight:
  // - Reset BL pointers when stagnation point moves.
  // - Shift BL arrays between top/bottom sides.
  // - Clamp Ue to small positive value at stagnation.
  const nb = ctxPanel.N;
  const istOld = blCtx.IST;
  const { ist, sst, sstGo, sstGp } = stfind(ctxPanel, nb);
  const istNew = ist + 1;

  blCtx.IST = istNew;
  blCtx.SST = sst;
  blCtx.SST_GO = sstGo;
  blCtx.SST_GP = sstGp;

  if (istNew === istOld) {
    xicalc(blCtx, ctxPanel);
    return;
  }

  // Reset BL pointers and Ue mapping for new stagnation point.
  iblpan(blCtx, nb, ctxPanel.NW ?? 0);
  uicalc(blCtx, qinv, qinvA);
  xicalc(blCtx, ctxPanel);
  iblsys(blCtx);

  if (istNew > istOld) {
    const idif = istNew - istOld;
    blCtx.ITRAN[1] += idif;
    blCtx.ITRAN[2] -= idif;

    for (let ibl = blCtx.NBL[1]; ibl >= idif + 2; ibl -= 1) {
      blCtx.CTAU[ibl][1] = blCtx.CTAU[ibl - idif][1];
      blCtx.THET[ibl][1] = blCtx.THET[ibl - idif][1];
      blCtx.DSTR[ibl][1] = blCtx.DSTR[ibl - idif][1];
      blCtx.UEDG[ibl][1] = blCtx.UEDG[ibl - idif][1];
    }

    const dudx = blCtx.UEDG[idif + 2][1] / blCtx.XSSI[idif + 2][1];
    for (let ibl = idif + 1; ibl >= 2; ibl -= 1) {
      blCtx.CTAU[ibl][1] = blCtx.CTAU[idif + 2][1];
      blCtx.THET[ibl][1] = blCtx.THET[idif + 2][1];
      blCtx.DSTR[ibl][1] = blCtx.DSTR[idif + 2][1];
      blCtx.UEDG[ibl][1] = dudx * blCtx.XSSI[ibl][1];
    }

    for (let ibl = 2; ibl <= blCtx.NBL[2]; ibl += 1) {
      blCtx.CTAU[ibl][2] = blCtx.CTAU[ibl + idif][2];
      blCtx.THET[ibl][2] = blCtx.THET[ibl + idif][2];
      blCtx.DSTR[ibl][2] = blCtx.DSTR[ibl + idif][2];
      blCtx.UEDG[ibl][2] = blCtx.UEDG[ibl + idif][2];
    }
  } else {
    const idif = istOld - istNew;
    blCtx.ITRAN[1] -= idif;
    blCtx.ITRAN[2] += idif;

    for (let ibl = blCtx.NBL[2]; ibl >= idif + 2; ibl -= 1) {
      blCtx.CTAU[ibl][2] = blCtx.CTAU[ibl - idif][2];
      blCtx.THET[ibl][2] = blCtx.THET[ibl - idif][2];
      blCtx.DSTR[ibl][2] = blCtx.DSTR[ibl - idif][2];
      blCtx.UEDG[ibl][2] = blCtx.UEDG[ibl - idif][2];
    }

    const dudx = blCtx.UEDG[idif + 2][2] / blCtx.XSSI[idif + 2][2];
    for (let ibl = idif + 1; ibl >= 2; ibl -= 1) {
      blCtx.CTAU[ibl][2] = blCtx.CTAU[idif + 2][2];
      blCtx.THET[ibl][2] = blCtx.THET[idif + 2][2];
      blCtx.DSTR[ibl][2] = blCtx.DSTR[idif + 2][2];
      blCtx.UEDG[ibl][2] = dudx * blCtx.XSSI[ibl][2];
    }

    for (let ibl = 2; ibl <= blCtx.NBL[1]; ibl += 1) {
      blCtx.CTAU[ibl][1] = blCtx.CTAU[ibl + idif][1];
      blCtx.THET[ibl][1] = blCtx.THET[ibl + idif][1];
      blCtx.DSTR[ibl][1] = blCtx.DSTR[ibl + idif][1];
      blCtx.UEDG[ibl][1] = blCtx.UEDG[ibl + idif][1];
    }
  }

  const ueps = 1.0e-7;
  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= blCtx.NBL[is]; ibl += 1) {
      if (blCtx.UEDG[ibl][is] <= ueps) {
        blCtx.UEDG[ibl][is] = ueps;
        const i = blCtx.IPAN[ibl][is];
        if (ctxPanel.QVIS && i >= 1 && i <= ctxPanel.QVIS.length - 1) {
          ctxPanel.QVIS[i] = blCtx.VTI[ibl][is] * ueps;
        }
        if (i >= 1 && i <= nb) {
          ctxPanel.GAM[i - 1] = blCtx.VTI[ibl][is] * ueps;
        }
      }
      blCtx.MASS[ibl][is] = blCtx.DSTR[ibl][is] * blCtx.UEDG[ibl][is];
    }
  }

  const ist0 = blCtx.IST - 1;
  const upperIdx = [];
  for (let i = ist0; i >= 0; i -= 1) {
    upperIdx.push(i);
  }
  const lowerIdx = [];
  for (let i = ist0 + 1; i < nb; i += 1) {
    lowerIdx.push(i);
  }
  blCtx.upperIdx = upperIdx;
  blCtx.lowerIdx = lowerIdx;
}

export {
  apcalc,
  xywake,
  qwcalc,
  ncalc,
  ggcalc,
  psilin,
  pswlin,
  qdcalc,
  stfind,
  iblpan,
  xicalc,
  uicalc,
  qvfue,
  qiset,
  gamqv,
  stmove,
};
