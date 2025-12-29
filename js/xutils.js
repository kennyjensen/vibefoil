// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

// Port of xutils.f ATANC helper.
// Utility math routines used by panel/wake geometry generation.

// Continuous arctangent that unwraps across 2π to avoid angular discontinuities.
function atanc(y, x, thold) {
  const pi = Math.PI;
  const tpi = 2.0 * Math.PI;
  const thnew = Math.atan2(y, x);
  const dthet = thnew - thold;
  const dtcorr = dthet - tpi * Math.trunc((dthet + Math.sign(dthet) * pi) / tpi);
  return thold + dtcorr;
}

// Exponentially stretched spacing for wake/panel grids (XFOIL SETEXP).
function setexp(s, ds1, smax, nn) {
  const sigma = smax / ds1;
  const nex = nn - 1;
  const rnex = nex;
  const rni = 1.0 / rnex;

  const aaa = rnex * (rnex - 1.0) * (rnex - 2.0) / 6.0;
  const bbb = rnex * (rnex - 1.0) / 2.0;
  const ccc = rnex - sigma;

  let disc = bbb ** 2 - 4.0 * aaa * ccc;
  disc = Math.max(0.0, disc);

  let ratio = 1.0;
  if (nex <= 1) {
    throw new Error('SETEXP: Cannot fill array. N too small.');
  } else if (nex === 2) {
    ratio = -ccc / bbb + 1.0;
  } else {
    ratio = (-bbb + Math.sqrt(disc)) / (2.0 * aaa) + 1.0;
  }

  if (ratio !== 1.0) {
    for (let iter = 0; iter < 100; iter += 1) {
      const sigman = (ratio ** nex - 1.0) / (ratio - 1.0);
      const res = sigman ** rni - sigma ** rni;
      const dresdr = rni * sigman ** rni
        * (rnex * ratio ** (nex - 1) - sigman) / (ratio ** nex - 1.0);
      const dratio = -res / dresdr;
      ratio += dratio;
      if (Math.abs(dratio) < 1.0e-5) break;
    }
  }

  s[0] = 0.0;
  let ds = ds1;
  for (let i = 1; i < nn; i += 1) {
    s[i] = s[i - 1] + ds;
    ds *= ratio;
  }
}

export { atanc, setexp };
