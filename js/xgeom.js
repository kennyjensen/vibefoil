// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

// Geometry utilities for airfoil repaneling and analysis (XFOIL port).

import { seval, deval, d2val, curv } from './spline.js';
import { atanc } from './xutils.js';

// Locate leading edge spline parameter SLE where chord is normal to tangent (LEFIND).
function lefind(x, xp, y, yp, s, n) {
  const dseps = (s[n - 1] - s[0]) * 1.0e-5;
  const xte = 0.5 * (x[0] + x[n - 1]);
  const yte = 0.5 * (y[0] + y[n - 1]);

  let guessIdx = n - 3;
  for (let i = 2; i <= n - 3; i += 1) {
    const dxte = x[i] - xte;
    const dyte = y[i] - yte;
    const dx = x[i + 1] - x[i];
    const dy = y[i + 1] - y[i];
    const dotp = dxte * dx + dyte * dy;
    if (dotp < 0.0) {
      guessIdx = i;
      break;
    }
  }

  let sle = s[guessIdx];
  if (sle === s[guessIdx - 1]) {
    return sle;
  }

  for (let iter = 0; iter < 50; iter += 1) {
    const xle = seval(sle, x, xp, s, n);
    const yle = seval(sle, y, yp, s, n);
    const dxds = deval(sle, x, xp, s, n);
    const dyds = deval(sle, y, yp, s, n);
    const dxdd = d2val(sle, x, xp, s, n);
    const dydd = d2val(sle, y, yp, s, n);

    const xchord = xle - xte;
    const ychord = yle - yte;
    const res = xchord * dxds + ychord * dyds;
    const ress = dxds * dxds + dyds * dyds + xchord * dxdd + ychord * dydd;

    let dsle = -res / ress;
    const limit = 0.02 * Math.abs(xchord + ychord);
    dsle = Math.max(dsle, -limit);
    dsle = Math.min(dsle, limit);
    sle += dsle;
    if (Math.abs(dsle) < dseps) {
      return sle;
    }
  }

  console.warn('LEFIND: LE point not found, using initial guess');
  return s[guessIdx];
}

// Display curvature at panel nodes (CLIS).
function clis(x, xp, y, yp, s, n) {
  const pi = Math.PI;
  let cmax = 0.0;
  let imax = 1;
  let arad = Math.atan2(-yp[0], -xp[0]);

  console.log('\n  i        x         y         s       theta        curv');
  for (let i = 0; i < n; i += 1) {
    if (i > 0) {
      arad = atanc(-yp[i], -xp[i], arad);
    }
    const adeg = arad * 180.0 / pi;
    const cv = curv(s[i], x, xp, y, yp, s, n);
    console.log(
      `${String(i + 1).padStart(3)}`
      + `${x[i].toFixed(5).padStart(10)}`
      + `${y[i].toFixed(5).padStart(10)}`
      + `${s[i].toFixed(5).padStart(10)}`
      + `${adeg.toFixed(3).padStart(11)}`
      + `${cv.toFixed(3).padStart(12)}`,
    );
    if (Math.abs(cv) > Math.abs(cmax)) {
      cmax = cv;
      imax = i + 1;
    }
  }
  console.log(
    `\n Maximum curvature =${cmax.toFixed(3).padStart(14)}`
    + `   at  i,x,y,s  = ${String(imax).padStart(3)}`
    + `${x[imax - 1].toFixed(4).padStart(9)}`
    + `${y[imax - 1].toFixed(4).padStart(9)}`
    + `${s[imax - 1].toFixed(4).padStart(9)}`,
  );
}

export { lefind, clis };
