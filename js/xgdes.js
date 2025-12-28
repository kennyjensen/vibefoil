// Geometry design utilities (GDES) ported from XFOIL, including FLAP.

import {
  scalc,
  segspl,
  seval,
  deval,
  d2val,
  sinvrt,
} from './spline.js';

const EPS = 1.0e-5;
const SFRAC = 0.33333;

function inside(x, y, n, xf, yf) {
  let angle = 0.0;
  for (let i = 0; i < n; i += 1) {
    const ip = i === n - 1 ? 0 : i + 1;
    const xb1 = x[i] - xf;
    const yb1 = y[i] - yf;
    const xb2 = x[ip] - xf;
    const yb2 = y[ip] - yf;
    const denom = Math.sqrt((xb1 * xb1 + yb1 * yb1) * (xb2 * xb2 + yb2 * yb2));
    angle += (xb1 * yb2 - yb1 * xb2) / denom;
  }
  return Math.abs(angle) > 1.0;
}

function getxyf(x, xp, y, yp, s, n, xf, yrel) {
  if (!Number.isFinite(xf)) {
    throw new Error('GETXYF: hinge x not specified.');
  }

  let tops = s[0] + (x[0] - xf);
  let bots = s[n - 1] - (x[n - 1] - xf);
  tops = sinvrt(tops, xf, x, xp, s, n);
  bots = sinvrt(bots, xf, x, xp, s, n);
  const topy = seval(tops, y, yp, s, n);
  const boty = seval(bots, y, yp, s, n);

  const rel = Number.isFinite(yrel) ? Math.min(1.0, Math.max(0.0, yrel)) : 0.5;
  const yf = topy * rel + boty * (1.0 - rel);
  return { tops, bots, xf, yf };
}

function sss(ss, del, xbf, ybf, x, xp, y, yp, s, n, iside) {
  const stot = Math.abs(s[n - 1] - s[0]);
  const sind = Math.sin(0.5 * Math.abs(del));
  const ssgn = iside === 1 ? -1.0 : 1.0;

  const dx0 = seval(ss, x, xp, s, n) - xbf;
  const dy0 = seval(ss, y, yp, s, n) - ybf;
  const rsq = dx0 * dx0 + dy0 * dy0;
  let s1 = ss - (sind * Math.sqrt(rsq) + EPS * stot) * ssgn;
  let s2 = ss + (sind * Math.sqrt(rsq) + EPS * stot) * ssgn;

  for (let iter = 0; iter < 10; iter += 1) {
    const x1 = seval(s1, x, xp, s, n);
    const x1p = deval(s1, x, xp, s, n);
    const y1 = seval(s1, y, yp, s, n);
    const y1p = deval(s1, y, yp, s, n);

    const x2 = seval(s2, x, xp, s, n);
    const x2p = deval(s2, x, xp, s, n);
    const y2 = seval(s2, y, yp, s, n);
    const y2p = deval(s2, y, yp, s, n);

    const r1sq = (x1 - xbf) ** 2 + (y1 - ybf) ** 2;
    const r2sq = (x2 - xbf) ** 2 + (y2 - ybf) ** 2;
    const r1 = Math.sqrt(r1sq);
    const r2 = Math.sqrt(r2sq);
    const rrsq = (x1 - x2) ** 2 + (y1 - y2) ** 2;
    const rr = Math.sqrt(rrsq);

    if (r1 <= EPS * stot || r2 <= EPS * stot) {
      return { s1: ss, s2: ss };
    }

    const r1S1 = (x1p * (x1 - xbf) + y1p * (y1 - ybf)) / r1;
    const r2S2 = (x2p * (x2 - xbf) + y2p * (y2 - ybf)) / r2;

    let rs1;
    let rs2;
    let a11;
    let a12;
    let a21;
    let a22;

    if (sind > 0.01) {
      if (rr === 0.0) {
        return { s1, s2 };
      }
      const rrS1 = (x1p * (x1 - x2) + y1p * (y1 - y2)) / rr;
      const rrS2 = -(x2p * (x1 - x2) + y2p * (y1 - y2)) / rr;

      rs1 = ((xbf - x1) * (x2 - x1) + (ybf - y1) * (y2 - y1)) / rr - sind * r1;
      a11 = ((xbf - x1) * (-x1p) + (ybf - y1) * (-y1p)) / rr
        + ((-x1p) * (x2 - x1) + (-y1p) * (y2 - y1)) / rr
        - ((xbf - x1) * (x2 - x1) + (ybf - y1) * (y2 - y1)) * rrS1 / rrsq
        - sind * r1S1;
      a12 = ((xbf - x1) * (x2p) + (ybf - y1) * (y2p)) / rr
        - ((xbf - x1) * (x2 - x1) + (ybf - y1) * (y2 - y1)) * rrS2 / rrsq;

      rs2 = r1 - r2;
      a21 = r1S1;
      a22 = -r2S2;
    } else {
      rs1 = (r1 + r2) * sind + (s1 - s2) * ssgn;
      a11 = r1S1 * sind + ssgn;
      a12 = r2S2 * sind - ssgn;

      const x1pp = d2val(s1, x, xp, s, n);
      const y1pp = d2val(s1, y, yp, s, n);
      const x2pp = d2val(s2, x, xp, s, n);
      const y2pp = d2val(s2, y, yp, s, n);
      const xtot = x1 + x2 - 2.0 * xbf;
      const ytot = y1 + y2 - 2.0 * ybf;

      rs2 = xtot * (x1p + x2p) + ytot * (y1p + y2p);
      a21 = x1p * (x1p + x2p) + y1p * (y1p + y2p) + xtot * x1pp + ytot * y1pp;
      a22 = x2p * (x1p + x2p) + y2p * (y1p + y2p) + xtot * x2pp + ytot * y2pp;
    }

    const det = a11 * a22 - a12 * a21;
    let ds1 = -(rs1 * a22 - a12 * rs2) / det;
    let ds2 = -(a11 * rs2 - rs1 * a21) / det;
    const step = 0.01 * stot;
    ds1 = Math.max(-step, Math.min(step, ds1));
    ds2 = Math.max(-step, Math.min(step, ds2));

    s1 += ds1;
    s2 += ds2;
    if (Math.abs(ds1) + Math.abs(ds2) < EPS * stot) {
      break;
    }
  }

  if (del === 0.0) {
    const mid = 0.5 * (s1 + s2);
    return { s1: mid, s2: mid };
  }

  return { s1, s2 };
}

function flap(xb, yb, flapParams) {
  if (!flapParams || !Number.isFinite(flapParams.deflection)) {
    return { xb, yb, nb: xb.length };
  }
  if (Math.abs(flapParams.deflection) < 1.0e-6) {
    return { xb, yb, nb: xb.length };
  }

  const x = Array.from(xb);
  const y = Array.from(yb);
  let nb = x.length;

  const sb = new Float64Array(nb);
  const xbp = new Float64Array(nb);
  const ybp = new Float64Array(nb);
  scalc(x, y, sb, nb);
  segspl(x, xbp, sb, nb);
  segspl(y, ybp, sb, nb);

  const hinge = getxyf(x, xbp, y, ybp, sb, nb, flapParams.x, flapParams.yrel);
  const xbf = hinge.xf;
  const ybf = hinge.yf;

  const ddef = flapParams.deflection;
  const rdef = ddef * (Math.PI / 180.0);
  if (rdef === 0.0) {
    return { xb: x, yb: y, nb };
  }

  let atop;
  let abot;
  if (inside(x, y, nb, xbf, ybf)) {
    atop = Math.max(0.0, -rdef);
    abot = Math.max(0.0, rdef);
  } else {
    const chx = deval(hinge.bots, x, xbp, sb, nb) - deval(hinge.tops, x, xbp, sb, nb);
    const chy = deval(hinge.bots, y, ybp, sb, nb) - deval(hinge.tops, y, ybp, sb, nb);
    const fvx = seval(hinge.bots, x, xbp, sb, nb) + seval(hinge.tops, x, xbp, sb, nb);
    const fvy = seval(hinge.bots, y, ybp, sb, nb) + seval(hinge.tops, y, ybp, sb, nb);
    const crsp = chx * (ybf - 0.5 * fvy) - chy * (xbf - 0.5 * fvx);
    if (crsp > 0.0) {
      atop = Math.max(0.0, rdef);
      abot = Math.max(0.0, rdef);
    } else {
      atop = Math.max(0.0, -rdef);
      abot = Math.max(0.0, -rdef);
    }
  }

  const topBreaks = sss(hinge.tops, atop, xbf, ybf, x, xbp, y, ybp, sb, nb, 1);
  const botBreaks = sss(hinge.bots, abot, xbf, ybf, x, xbp, y, ybp, sb, nb, 2);
  const st1 = topBreaks.s1;
  const st2 = topBreaks.s2;
  const sb1 = botBreaks.s1;
  const sb2 = botBreaks.s2;

  const xt1 = seval(st1, x, xbp, sb, nb);
  const yt1 = seval(st1, y, ybp, sb, nb);
  const xt2 = seval(st2, x, xbp, sb, nb);
  const yt2 = seval(st2, y, ybp, sb, nb);
  const xb1 = seval(sb1, x, xbp, sb, nb);
  const yb1 = seval(sb1, y, ybp, sb, nb);
  const xb2 = seval(sb2, x, xbp, sb, nb);
  const yb2 = seval(sb2, y, ybp, sb, nb);

  let it1 = 1;
  let it2 = 1;
  let ib1 = 1;
  let ib2 = nb - 1;
  for (let i = 0; i < nb - 1; i += 1) {
    if (sb[i] <= st1 && sb[i + 1] > st1) it1 = i + 1;
    if (sb[i] < st2 && sb[i + 1] >= st2) it2 = i;
    if (sb[i] <= sb1 && sb[i + 1] > sb1) ib1 = i;
    if (sb[i] < sb2 && sb[i + 1] >= sb2) ib2 = i + 1;
  }

  const dsavg = (sb[nb - 1] - sb[0]) / Math.max(nb - 1, 1);

  let xt1new = 0.0;
  let yt1new = 0.0;
  let xt2new = 0.0;
  let yt2new = 0.0;
  let xb1new = 0.0;
  let yb1new = 0.0;
  let xb2new = 0.0;
  let yb2new = 0.0;
  let lt1new = false;
  let lt2new = false;
  let lb1new = false;
  let lb2new = false;

  if (atop !== 0.0) {
    const st1p = st1 + SFRAC * (sb[it1] - st1);
    const st1q = st1 + SFRAC * (sb[Math.min(it1 + 1, nb - 1)] - st1);
    if (sb[it1] < st1q) {
      xt1new = seval(st1q, x, xbp, sb, nb);
      yt1new = seval(st1q, y, ybp, sb, nb);
      lt1new = false;
    } else {
      xt1new = seval(st1p, x, xbp, sb, nb);
      yt1new = seval(st1p, y, ybp, sb, nb);
      lt1new = true;
    }

    const st2p = st2 + SFRAC * (sb[it2] - st2);
    const it2q = Math.max(it2 - 1, 0);
    const st2q = st2 + SFRAC * (sb[it2q] - st2);
    if (sb[it2] > st2q) {
      xt2new = seval(st2q, x, xbp, sb, nb);
      yt2new = seval(st2q, y, ybp, sb, nb);
      lt2new = false;
    } else {
      xt2new = seval(st2p, x, xbp, sb, nb);
      yt2new = seval(st2p, y, ybp, sb, nb);
      lt2new = true;
    }
  }

  if (abot !== 0.0) {
    const sb1p = sb1 + SFRAC * (sb[ib1] - sb1);
    const sb1q = sb1 + SFRAC * (sb[Math.max(ib1 - 1, 0)] - sb1);
    if (sb[ib1] > sb1q) {
      xb1new = seval(sb1q, x, xbp, sb, nb);
      yb1new = seval(sb1q, y, ybp, sb, nb);
      lb1new = false;
    } else {
      xb1new = seval(sb1p, x, xbp, sb, nb);
      yb1new = seval(sb1p, y, ybp, sb, nb);
      lb1new = true;
    }

    const sb2p = sb2 + SFRAC * (sb[ib2] - sb2);
    const ib2q = Math.min(ib2 + 1, nb - 1);
    const sb2q = sb2 + SFRAC * (sb[ib2q] - sb2);
    if (sb[ib2] < sb2q) {
      xb2new = seval(sb2q, x, xbp, sb, nb);
      yb2new = seval(sb2q, y, ybp, sb, nb);
      lb2new = false;
    } else {
      xb2new = seval(sb2p, x, xbp, sb, nb);
      yb2new = seval(sb2p, y, ybp, sb, nb);
      lb2new = true;
    }
  }

  const sind = Math.sin(rdef);
  const cosd = Math.cos(rdef);
  for (let i = 0; i < nb; i += 1) {
    if (i >= it1 && i <= ib1) continue;
    const xbar = x[i] - xbf;
    const ybar = y[i] - ybf;
    x[i] = xbf + xbar * cosd + ybar * sind;
    y[i] = ybf - xbar * sind + ybar * cosd;
  }

  let idif = it1 - it2 - 1;
  if (idif > 0) {
    x.splice(it2 + 1, idif);
    y.splice(it2 + 1, idif);
    nb -= idif;
    it1 -= idif;
    ib1 -= idif;
    ib2 -= idif;
  }

  idif = ib2 - ib1 - 1;
  if (idif > 0) {
    x.splice(ib1 + 1, idif);
    y.splice(ib1 + 1, idif);
    nb -= idif;
    ib2 -= idif;
  }

  if (atop === 0.0) {
    const dsnew = Math.abs(rdef) * Math.sqrt((xt1 - xbf) ** 2 + (yt1 - ybf) ** 2);
    const npadd = Math.floor(1.5 * dsnew / dsavg + 1.0);
    if (npadd > 0) {
      x.splice(it1, 0, ...new Array(npadd).fill(0.0));
      y.splice(it1, 0, ...new Array(npadd).fill(0.0));
      nb += npadd;
      it1 += npadd;
      ib1 += npadd;
      ib2 += npadd;
      const dang = rdef / npadd;
      const xbar = xt1 - xbf;
      const ybar = yt1 - ybf;
      for (let ip = 1; ip <= npadd; ip += 1) {
        const ang = dang * (ip - 0.5);
        const ca = Math.cos(ang);
        const sa = Math.sin(ang);
        const idx = it1 - ip;
        x[idx] = xbf + xbar * ca + ybar * sa;
        y[idx] = ybf - xbar * sa + ybar * ca;
      }
    }
  } else {
    let npadd = 1;
    if (lt2new) npadd += 1;
    if (lt1new) npadd += 1;
    x.splice(it1, 0, ...new Array(npadd).fill(0.0));
    y.splice(it1, 0, ...new Array(npadd).fill(0.0));
    nb += npadd;
    it1 += npadd;
    ib1 += npadd;
    ib2 += npadd;

    if (lt1new) {
      x[it1 - 1] = xt1new;
      y[it1 - 1] = yt1new;
      x[it1 - 2] = xt1;
      y[it1 - 2] = yt1;
    } else {
      x[it1] = xt1new;
      y[it1] = yt1new;
      x[it1 - 1] = xt1;
      y[it1 - 1] = yt1;
    }

    const xbar = xt2new - xbf;
    const ybar = yt2new - ybf;
    if (lt2new) {
      x[it2 + 1] = xbf + xbar * cosd + ybar * sind;
      y[it2 + 1] = ybf - xbar * sind + ybar * cosd;
    } else {
      x[it2] = xbf + xbar * cosd + ybar * sind;
      y[it2] = ybf - xbar * sind + ybar * cosd;
    }
  }

  if (abot === 0.0) {
    const dsnew = Math.abs(rdef) * Math.sqrt((xb1 - xbf) ** 2 + (yb1 - ybf) ** 2);
    const npadd = Math.floor(1.5 * dsnew / dsavg + 1.0);
    if (npadd > 0) {
      x.splice(ib2, 0, ...new Array(npadd).fill(0.0));
      y.splice(ib2, 0, ...new Array(npadd).fill(0.0));
      nb += npadd;
      ib2 += npadd;
      const dang = rdef / npadd;
      const xbar = xb1 - xbf;
      const ybar = yb1 - ybf;
      for (let ip = 1; ip <= npadd; ip += 1) {
        const ang = dang * (ip - 0.5);
        const ca = Math.cos(ang);
        const sa = Math.sin(ang);
        const idx = ib1 + ip;
        x[idx] = xbf + xbar * ca + ybar * sa;
        y[idx] = ybf - xbar * sa + ybar * ca;
      }
    }
  } else {
    let npadd = 1;
    if (lb2new) npadd += 1;
    if (lb1new) npadd += 1;
    x.splice(ib2, 0, ...new Array(npadd).fill(0.0));
    y.splice(ib2, 0, ...new Array(npadd).fill(0.0));
    nb += npadd;
    ib2 += npadd;

    if (lb1new) {
      x[ib1 + 1] = xb1new;
      y[ib1 + 1] = yb1new;
      x[ib1 + 2] = xb1;
      y[ib1 + 2] = yb1;
    } else {
      x[ib1] = xb1new;
      y[ib1] = yb1new;
      x[ib1 + 1] = xb1;
      y[ib1 + 1] = yb1;
    }

    const xbar = xb2new - xbf;
    const ybar = yb2new - ybf;
    if (lb2new) {
      x[ib2 - 1] = xbf + xbar * cosd + ybar * sind;
      y[ib2 - 1] = ybf - xbar * sind + ybar * cosd;
    } else {
      x[ib2] = xbf + xbar * cosd + ybar * sind;
      y[ib2] = ybf - xbar * sind + ybar * cosd;
    }
  }

  return { xb: x, yb: y, nb };
}

export {
  flap,
  getxyf,
  inside,
  sss,
};
