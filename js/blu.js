// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

// Port of blu.f (BL utilities). Exposes CFT and a profile helper.
// Auxiliary BL profile evaluation and skin-friction correlations.

const KPRX = 129;

// Assemble a velocity profile from integral quantities (turbulent/laminar).
// Mirrors BLU's profile reconstruction using PRWALL/UWALL or FS fallback.
function bluProfile(params, deps) {
  if (!params || typeof params.hk !== 'number') {
    throw new Error('bluProfile requires params.hk');
  }
  if (!deps || typeof deps.prwall !== 'function' || typeof deps.uwall !== 'function' || typeof deps.fs !== 'function') {
    throw new Error('bluProfile requires deps.prwall, deps.uwall, deps.fs');
  }

  const hk = params.hk;
  const ret = params.ret ?? 0.0;
  const th = params.th ?? 1.0;
  const ue = params.ue ?? 1.0;
  const mach = params.mach ?? 0.0;
  const nn = Math.min(params.nn ?? KPRX, KPRX);

  const msq = mach * mach;
  const uo = 1.0;
  const dk = hk * th;
  const ct = 0.0;

  const yy = new Array(nn + 1).fill(0.0);
  const xx = new Array(nn + 1).fill(0.0);
  const ffs = new Array(nn + 1).fill(0.0);
  const sfs = new Array(nn + 1).fill(0.0);

  let de = 0.0;
  let us = 0.0;
  let cf = 0.0;
  let bb = 0.0;

  if (ret > 0.0) {
    const pr = deps.prwall({
      dstar: dk,
      theta: th,
      uo,
      rt: ret,
      ms: msq,
      ct,
    });
    de = pr.de;
    us = pr.us;
    cf = pr.cf;
    bb = pr.bb;

    deps.uwall({
      th,
      uo,
      de,
      us,
      rt: ret,
      cf,
      bb,
      yy,
      xx,
      n: nn,
    });

    for (let k = 1; k <= nn; k += 1) {
      xx[k] *= ue;
    }
  } else {
    const inorm = 3;
    const ispec = 2;
    const hspec = hk;
    const etae = 1.5 * (3.15 + 1.72 / (hk - 1.0) + hk) + 2.0;
    const geo = 1.0;

    deps.fs({
      inorm,
      ispec,
      bspec: 0.0,
      hspec,
      n: nn,
      etae,
      geo,
      eta: yy,
      f: ffs,
      u: xx,
      s: sfs,
    });
    de = etae * th;

    for (let k = 1; k <= nn; k += 1) {
      xx[k] *= ue;
      yy[k] *= th;
    }
  }

  return {
    yy,
    xx,
    de,
    us,
    cf,
    bb,
    nn,
  };
}

// Turbulent skin-friction correlation with compressibility correction (CFT).
function cft(hk, rt, msq, cffac = 1.0) {
  const gam = 1.4;
  const gm1 = gam - 1.0;
  const fc = Math.sqrt(1.0 + 0.5 * gm1 * msq);
  let grt = Math.log(rt / fc);
  grt = Math.max(grt, 3.0);

  const gex = -1.74 - 0.31 * hk;

  let arg = -1.33 * hk;
  arg = Math.max(-20.0, arg);

  const thk = Math.tanh(4.0 - hk / 0.875);

  const cfo = cffac * 0.3 * Math.exp(arg) * (grt / 2.3026) ** gex;
  const cf = (cfo + 1.1e-4 * (thk - 1.0)) / fc;
  const cfHk = (-1.33 * cfo - 0.31 * Math.log(grt / 2.3026) * cfo
    - 1.1e-4 * (1.0 - thk ** 2) / 0.875) / fc;
  const cfRt = gex * cfo / (fc * grt) / rt;
  const cfMsq = gex * cfo / (fc * grt) * (-0.25 * gm1 / fc ** 2)
    - 0.25 * gm1 * cf / fc ** 2;

  return {
    cf,
    cfHk,
    cfRt,
    cfMsq,
  };
}

export { bluProfile, cft };
