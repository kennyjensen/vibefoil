// Port of src/spline.f from XFOIL (Mark Drela). Numerical formulas preserved.

const NMAX = 1000;
const KMAX = 32;

// Bounds guard mirroring the original Fortran limits, preserving algorithmic
// assumptions about fixed-size work arrays and spline segment counts.
function checkN(n, name) {
  if (n > NMAX) {
    throw new Error(`${name}: array overflow, increase NMAX`);
  }
}

// Thomas algorithm for a tridiagonal system; used by spline solvers to obtain
// first derivatives under boundary constraints (natural or clamped).
function trisol(a, b, c, d, kk, dOffset = 0) {
  for (let k = 1; k < kk; k += 1) {
    const km = k - 1;
    c[km] = c[km] / a[km];
    const dk = d[dOffset + km] / a[km];
    d[dOffset + km] = dk;
    a[k] = a[k] - b[k] * c[km];
    d[dOffset + k] = d[dOffset + k] - b[k] * dk;
  }

  d[dOffset + kk - 1] = d[dOffset + kk - 1] / a[kk - 1];

  for (let k = kk - 2; k >= 0; k -= 1) {
    d[dOffset + k] = d[dOffset + k] - c[k] * d[dOffset + k + 1];
  }
}

// Cubic spline with end conditions from XFOIL's SPLINE: solves for derivatives
// in s-parametric space, consistent with NASA airfoil geometry workflows.
function spline(x, xs, s, n) {
  checkN(n, 'SPLINE');
  const a = new Float64Array(n);
  const b = new Float64Array(n);
  const c = new Float64Array(n);

  for (let i = 1; i <= n - 2; i += 1) {
    const dsm = s[i] - s[i - 1];
    const dsp = s[i + 1] - s[i];
    b[i] = dsp;
    a[i] = 2.0 * (dsm + dsp);
    c[i] = dsm;
    xs[i] = 3.0 * ((x[i + 1] - x[i]) * dsm / dsp + (x[i] - x[i - 1]) * dsp / dsm);
  }

  a[0] = 2.0;
  c[0] = 1.0;
  xs[0] = 3.0 * (x[1] - x[0]) / (s[1] - s[0]);

  b[n - 1] = 1.0;
  a[n - 1] = 2.0;
  xs[n - 1] = 3.0 * (x[n - 1] - x[n - 2]) / (s[n - 1] - s[n - 2]);

  trisol(a, b, c, xs, n);
}

// Spline derivative solver with explicit end-slope options:
//  999/-999 sentinel values reproduce XFOIL's "natural"/"not-a-knot" handling.
function splindOffset(x, xs, s, n, xs1, xs2, offset) {
  checkN(n, 'SPLIND');
  const a = new Float64Array(n);
  const b = new Float64Array(n);
  const c = new Float64Array(n);

  for (let i = 1; i <= n - 2; i += 1) {
    const im1 = offset + i - 1;
    const i0 = offset + i;
    const ip1 = offset + i + 1;
    const dsm = s[i0] - s[im1];
    const dsp = s[ip1] - s[i0];
    b[i] = dsp;
    a[i] = 2.0 * (dsm + dsp);
    c[i] = dsm;
    xs[offset + i] = 3.0 * ((x[ip1] - x[i0]) * dsm / dsp + (x[i0] - x[im1]) * dsp / dsm);
  }

  if (xs1 === 999.0) {
    a[0] = 2.0;
    c[0] = 1.0;
    xs[offset + 0] = 3.0 * (x[offset + 1] - x[offset + 0]) / (s[offset + 1] - s[offset + 0]);
  } else if (xs1 === -999.0) {
    a[0] = 1.0;
    c[0] = 1.0;
    xs[offset + 0] = 2.0 * (x[offset + 1] - x[offset + 0]) / (s[offset + 1] - s[offset + 0]);
  } else {
    a[0] = 1.0;
    c[0] = 0.0;
    xs[offset + 0] = xs1;
  }

  if (xs2 === 999.0) {
    b[n - 1] = 1.0;
    a[n - 1] = 2.0;
    xs[offset + n - 1] = 3.0 * (x[offset + n - 1] - x[offset + n - 2]) / (s[offset + n - 1] - s[offset + n - 2]);
  } else if (xs2 === -999.0) {
    b[n - 1] = 1.0;
    a[n - 1] = 1.0;
    xs[offset + n - 1] = 2.0 * (x[offset + n - 1] - x[offset + n - 2]) / (s[offset + n - 1] - s[offset + n - 2]);
  } else {
    a[n - 1] = 1.0;
    b[n - 1] = 0.0;
    xs[offset + n - 1] = xs2;
  }

  if (n === 2 && xs1 === -999.0 && xs2 === -999.0) {
    b[n - 1] = 1.0;
    a[n - 1] = 2.0;
    xs[offset + n - 1] = 3.0 * (x[offset + n - 1] - x[offset + n - 2]) / (s[offset + n - 1] - s[offset + n - 2]);
  }

  trisol(a, b, c, xs, n, offset);
}

// Convenience wrapper for a contiguous segment.
function splind(x, xs, s, n, xs1, xs2) {
  splindOffset(x, xs, s, n, xs1, xs2, 0);
}

// Piecewise linear derivative estimate (used as a fallback / initializer).
function splina(x, xs, s, n) {
  let lend = true;
  let xs1 = 0.0;
  let xs2 = 0.0;

  for (let i = 0; i < n - 1; i += 1) {
    const ds = s[i + 1] - s[i];
    if (ds === 0.0) {
      xs[i] = xs1;
      lend = true;
    } else {
      const dx = x[i + 1] - x[i];
      xs2 = dx / ds;
      if (lend) {
        xs[i] = xs2;
        lend = false;
      } else {
        xs[i] = 0.5 * (xs1 + xs2);
      }
    }
    xs1 = xs2;
  }
  xs[n - 1] = xs1;
}

// Cubic spline evaluation in s; follows Abbott & von Doenhoff style parametric
// reconstruction with Hermite-form coefficients.
function seval(ss, x, xs, s, n) {
  let ilow = 0;
  let i = n - 1;

  while (i - ilow > 1) {
    const imid = Math.floor((i + ilow) / 2);
    if (ss < s[imid]) {
      i = imid;
    } else {
      ilow = imid;
    }
  }

  const ds = s[i] - s[i - 1];
  const t = (ss - s[i - 1]) / ds;
  const cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
  const cx2 = ds * xs[i] - x[i] + x[i - 1];
  return t * x[i] + (1.0 - t) * x[i - 1] + (t - t * t) * ((1.0 - t) * cx1 - t * cx2);
}

// First derivative of spline evaluation, used for curvature and BL metrics.
function deval(ss, x, xs, s, n) {
  let ilow = 0;
  let i = n - 1;

  while (i - ilow > 1) {
    const imid = Math.floor((i + ilow) / 2);
    if (ss < s[imid]) {
      i = imid;
    } else {
      ilow = imid;
    }
  }

  const ds = s[i] - s[i - 1];
  const t = (ss - s[i - 1]) / ds;
  const cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
  const cx2 = ds * xs[i] - x[i] + x[i - 1];
  let devalValue = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
  devalValue = devalValue / ds;
  return devalValue;
}

// Second derivative of spline evaluation, used for curvature and stability checks.
function d2val(ss, x, xs, s, n) {
  let ilow = 0;
  let i = n - 1;

  while (i - ilow > 1) {
    const imid = Math.floor((i + ilow) / 2);
    if (ss < s[imid]) {
      i = imid;
    } else {
      ilow = imid;
    }
  }

  const ds = s[i] - s[i - 1];
  const t = (ss - s[i - 1]) / ds;
  const cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
  const cx2 = ds * xs[i] - x[i] + x[i - 1];
  let d2valValue = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;
  d2valValue = d2valValue / (ds ** 2);
  return d2valValue;
}

// Curvature kappa(s) from parametric cubic splines, stable near the LE by
// imposing a lower bound on |ds/dt| as done in the original solver.
function curv(ss, x, xs, y, ys, s, n) {
  let ilow = 0;
  let i = n - 1;

  while (i - ilow > 1) {
    const imid = Math.floor((i + ilow) / 2);
    if (ss < s[imid]) {
      i = imid;
    } else {
      ilow = imid;
    }
  }

  const ds = s[i] - s[i - 1];
  const t = (ss - s[i - 1]) / ds;

  const cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
  const cx2 = ds * xs[i] - x[i] + x[i - 1];
  const xd = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
  const xdd = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;

  const cy1 = ds * ys[i - 1] - y[i] + y[i - 1];
  const cy2 = ds * ys[i] - y[i] + y[i - 1];
  const yd = y[i] - y[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cy1 + t * (3.0 * t - 2.0) * cy2;
  const ydd = (6.0 * t - 4.0) * cy1 + (6.0 * t - 2.0) * cy2;

  let sd = Math.sqrt(xd * xd + yd * yd);
  sd = Math.max(sd, 0.001 * ds);

  return (xd * ydd - yd * xdd) / (sd ** 3);
}

// Derivative of curvature with respect to s; used in some smoothing/logistics.
function curvs(ss, x, xs, y, ys, s, n) {
  let ilow = 0;
  let i = n - 1;

  while (i - ilow > 1) {
    const imid = Math.floor((i + ilow) / 2);
    if (ss < s[imid]) {
      i = imid;
    } else {
      ilow = imid;
    }
  }

  const ds = s[i] - s[i - 1];
  const t = (ss - s[i - 1]) / ds;

  const cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
  const cx2 = ds * xs[i] - x[i] + x[i - 1];
  const xd = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
  const xdd = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;
  const xddd = 6.0 * cx1 + 6.0 * cx2;

  const cy1 = ds * ys[i - 1] - y[i] + y[i - 1];
  const cy2 = ds * ys[i] - y[i] + y[i - 1];
  const yd = y[i] - y[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cy1 + t * (3.0 * t - 2.0) * cy2;
  const ydd = (6.0 * t - 4.0) * cy1 + (6.0 * t - 2.0) * cy2;
  const yddd = 6.0 * cy1 + 6.0 * cy2;

  let sd = Math.sqrt(xd * xd + yd * yd);
  sd = Math.max(sd, 0.001 * ds);

  const bot = sd ** 3;
  const dbotdt = 3.0 * sd * (xd * xdd + yd * ydd);

  const top = xd * ydd - yd * xdd;
  const dtopdt = xd * yddd - yd * xddd;

  return (dtopdt * bot - dbotdt * top) / (bot ** 2);
}

// Invert s(x) by Newton iteration on the spline (XFOIL's SINVRT).
// Keeps the same termination tolerance and returns input on failure.
function sinvrt(si, xi, x, xs, s, n) {
  const sisav = si;
  let siValue = si;

  for (let iter = 0; iter < 10; iter += 1) {
    const res = seval(siValue, x, xs, s, n) - xi;
    const resp = deval(siValue, x, xs, s, n);
    const ds = -res / resp;
    siValue += ds;
    if (Math.abs(ds / (s[n - 1] - s[0])) < 1.0e-5) {
      return siValue;
    }
  }

  // Fortran writes a warning and returns input value.
  return sisav;
}

// Cumulative arc-length parameterization along (x,y), classic spline setup.
function scalc(x, y, s, n) {
  s[0] = 0.0;
  for (let i = 1; i < n; i += 1) {
    const dx = x[i] - x[i - 1];
    const dy = y[i] - y[i - 1];
    s[i] = s[i - 1] + Math.sqrt(dx * dx + dy * dy);
  }
}

// Iterative arc-length reparameterization with Simpson-like correction, used
// to reduce interpolation error in geometry spline fits.
function splnxy(x, xs, y, ys, s, n) {
  const xt = new Float64Array(KMAX + 1);
  const yt = new Float64Array(KMAX + 1);

  const kk = KMAX;
  const npass = 10;

  scalc(x, y, s, n);
  segspl(x, xs, s, n);
  segspl(y, ys, s, n);

  for (let ipass = 0; ipass < npass; ipass += 1) {
    let serr = 0.0;

    let ds = s[1] - s[0];
    for (let i = 1; i < n; i += 1) {
      const dx = x[i] - x[i - 1];
      const dy = y[i] - y[i - 1];

      const cx1 = ds * xs[i - 1] - dx;
      const cx2 = ds * xs[i] - dx;
      const cy1 = ds * ys[i - 1] - dy;
      const cy2 = ds * ys[i] - dy;

      xt[0] = 0.0;
      yt[0] = 0.0;
      for (let k = 1; k <= kk - 1; k += 1) {
        const t = k / kk;
        xt[k] = t * dx + (t - t * t) * ((1.0 - t) * cx1 - t * cx2);
        yt[k] = t * dy + (t - t * t) * ((1.0 - t) * cy1 - t * cy2);
      }
      xt[kk] = dx;
      yt[kk] = dy;

      let sint1 = 0.0;
      for (let k = 1; k <= kk; k += 1) {
        const ddx = xt[k] - xt[k - 1];
        const ddy = yt[k] - yt[k - 1];
        sint1 += Math.sqrt(ddx * ddx + ddy * ddy);
      }

      let sint2 = 0.0;
      for (let k = 2; k <= kk; k += 2) {
        const ddx = xt[k] - xt[k - 2];
        const ddy = yt[k] - yt[k - 2];
        sint2 += Math.sqrt(ddx * ddx + ddy * ddy);
      }

      const sint = (4.0 * sint1 - sint2) / 3.0;

      if (Math.abs(sint - ds) > Math.abs(serr)) {
        serr = sint - ds;
      }

      if (i < n - 1) {
        ds = s[i + 1] - s[i];
      }

      s[i] = s[i - 1] + Math.sqrt(sint);
    }

    serr = serr / (s[n - 1] - s[0]);
    // Fortran writes IPASS and SERR each pass.
    console.log(ipass + 1, serr);

    segspl(x, xs, s, n);
    segspl(y, ys, s, n);

    if (Math.abs(serr) < 1.0e-7) {
      return;
    }
  }
}

// Segment-wise spline derivative solve that handles repeated s entries (cusps).
function segspl(x, xs, s, n) {
  if (s[0] === s[1]) {
    throw new Error('SEGSPL:  First input point duplicated');
  }
  if (s[n - 1] === s[n - 2]) {
    throw new Error('SEGSPL:  Last  input point duplicated');
  }

  let iseg0 = 0;
  for (let iseg = 1; iseg <= n - 3; iseg += 1) {
    if (s[iseg] === s[iseg + 1]) {
      const nseg = iseg - iseg0 + 1;
      splindOffset(x, xs, s, nseg, -999.0, -999.0, iseg0);
      iseg0 = iseg + 1;
    }
  }

  const nseg = n - iseg0;
  splindOffset(x, xs, s, nseg, -999.0, -999.0, iseg0);
}

// Segment-wise spline derivative solve with user-specified end slopes.
function segspld(x, xs, s, n, xs1, xs2) {
  if (s[0] === s[1]) {
    throw new Error('SEGSPL:  First input point duplicated');
  }
  if (s[n - 1] === s[n - 2]) {
    throw new Error('SEGSPL:  Last  input point duplicated');
  }

  let iseg0 = 0;
  for (let iseg = 1; iseg <= n - 3; iseg += 1) {
    if (s[iseg] === s[iseg + 1]) {
      const nseg = iseg - iseg0 + 1;
      splindOffset(x, xs, s, nseg, xs1, xs2, iseg0);
      iseg0 = iseg + 1;
    }
  }

  const nseg = n - iseg0;
  splindOffset(x, xs, s, nseg, xs1, xs2, iseg0);
}

export {
  spline,
  splind,
  splina,
  trisol,
  seval,
  deval,
  d2val,
  curv,
  curvs,
  sinvrt,
  scalc,
  splnxy,
  segspl,
  segspld,
};
