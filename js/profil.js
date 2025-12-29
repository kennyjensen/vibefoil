// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

// Port of src/profil.f (profile routines).
// Velocity profile construction for BL integral methods and closures.
import { cft } from './blu.js';

// Sign helper for branch selection in velocity-profile relations.
function sign1(value) {
  return value >= 0.0 ? 1.0 : -1.0;
}

// Turbulent wall-profile solver: returns displacement thickness and skin friction.
// Implements PRWALL from profil.f.
function prwall(params) {
  // Fortran comments (PRWALL) highlight:
  // - Turbulent wall profile with compressibility corrections.
  const dstar = params.dstar;
  const theta = params.theta;
  const rt = params.rt;
  const ms = params.ms;
  const ct = params.ct ?? 0.0;
  const cffac = params.cffac ?? 1.0;

  const n = 65;
  const eta = new Array(n + 1).fill(0.0);
  const uip = new Array(n + 1).fill(0.0);
  const uipDp = new Array(n + 1).fill(0.0);
  const g = new Array(n + 1).fill(0.0);
  const gBb = new Array(n + 1).fill(0.0);

  const hpi = 1.570796327;
  const ak = 0.09;
  const vkap = 0.40;
  const vb = 5.0;

  const hk = dstar / theta;
  const uo = 1.0;
  let bb = 1.0;

  const cftOut = cft(hk, rt, ms, cffac);
  let cf = cftOut.cf;

  let sgn = sign1(cf);
  let ut = sgn * Math.sqrt(0.5 * Math.abs(cf));

  let ui = Math.min((ut / ak) * hpi, 0.90);
  let d0 = hk * theta / (1.0 - 0.5 * (uo + ui));

  const ebk = Math.exp(-vb * vkap);

  let b11 = 0.0;
  let b12 = 0.0;
  let b21 = 0.0;
  let b22 = 0.0;

  let dsn = 0.0;
  let dsnTh = 0.0;
  let dsnRt = 0.0;
  let dsnMs = 0.0;
  let thn = 0.0;
  let thnTh = 0.0;
  let thnRt = 0.0;
  let thnMs = 0.0;

  let duo = 0.0;
  let duoUt = 0.0;
  let duoDo = 0.0;
  let duoTh = 0.0;
  let duoRt = 0.0;
  let duoMs = 0.0;

  for (let iter = 1; iter <= 40; iter += 1) {
    sgn = sign1(ut);

    const dp = sgn * ut * rt * d0 / theta;
    const dpDo = sgn * ut * rt / theta;
    const dpUt = sgn * rt * d0 / theta;

    const dpTh = -dp / theta;
    const dpRt = sgn * ut * d0 / theta;
    const dpMs = 0.0;

    let upe = Math.log(dp) / vkap + vb;

    let exu = 0.0;
    let dpU = 0.0;
    for (let itup = 1; itup <= 5; itup += 1) {
      const uk = upe * vkap;
      const arg = uk - vb * vkap;
      exu = Math.exp(arg);
      const rez = upe + exu - ebk * (1.0 + uk + uk * uk / 2.0 + uk * uk * uk / 6.0) - dp;
      dpU = 1.0 + (exu - ebk * (1.0 + uk + uk * uk / 2.0)) * vkap;

      if (Math.abs(rez / dp) < 1.0e-5) {
        break;
      }

      const dupe = -rez / dpU;
      upe += dupe;
    }

    const upeDp = 1.0 / dpU;

    const uk = upe * vkap;
    const dpUu = (exu - ebk * (1.0 + uk)) * vkap ** 2;
    const dpUuu = (exu - ebk) * vkap ** 3;

    const upd = 1.0 / dpU;
    const updDp = (-1.0 / dpU ** 3) * dpUu;

    const updd = updDp;
    const upddDp = (-1.0 / dpU ** 4) * dpUuu + (3.0 / dpU ** 5) * dpUu ** 2;

    const dc2 = 0.5 * dp * dp * updd - dp * upd;
    const dc2Dp = dp * updd - upd + 0.5 * dp * dp * upddDp - dp * updDp;

    const dc3 = -(dp * dp * updd - dp * upd) / 3.0;
    const dc3Dp = -(2.0 * dp * updd - upd) / 3.0 - (dp * dp * upddDp - dp * updDp) / 3.0;

    duo = uo - ut * (upe + dc2 + dc3);
    const duoDp = -ut * (upeDp + dc2Dp + dc3Dp);

    duoUt = -(upe + dc2 + dc3) + duoDp * dpUt;
    duoDo = duoDp * dpDo;

    duoTh = duoDp * dpTh;
    duoRt = duoDp * dpRt;
    duoMs = duoDp * dpMs;

    const bb1 = 3.0 * (bb + 2.0) * (bb + 3.0) / (bb + 7.0);
    const bb1Bb = 3.0 * (bb * 2.0 + 5.0) / (bb + 7.0) - bb1 / (bb + 7.0);
    const bb2 = -5.0 * (bb + 1.0) * (bb + 3.0) / (bb + 7.0);
    const bb2Bb = -5.0 * (bb * 2.0 + 4.0) / (bb + 7.0) - bb2 / (bb + 7.0);
    const bb3 = 2.0 * (bb + 1.0) * (bb + 2.0) / (bb + 7.0);
    const bb3Bb = 2.0 * (bb * 2.0 + 3.0) / (bb + 7.0) - bb3 / (bb + 7.0);

    const exupe = exu;
    const dexu = (exupe - ebk) / (n - 1);

    eta[1] = 0.0;
    uip[1] = 0.0;
    uipDp[1] = 0.0;
    g[1] = 0.0;
    gBb[1] = 0.0;

    for (let i = 2; i <= n; i += 1) {
      const exui = ebk + (dexu - 0.75 * dexu * (n - i) / (n - 1)) * (i - 1);
      const uki = Math.log(exui) + vb * vkap;
      const up = uki / vkap;

      const yp = up + exui - ebk * (1.0 + uki + uki * uki / 2.0 + uki * uki * uki / 6.0);
      const ypU = 1.0 + (exui - ebk * (1.0 + uki + uki * uki / 2.0)) * vkap;
      const ypUu = (exui - ebk * (1.0 + uki)) * vkap ** 2;

      const et = yp / dp;

      uip[i] = up + dc2 * et ** 2 + dc3 * et ** 3;
      uipDp[i] = dc2Dp * et ** 2 + dc3Dp * et ** 3;

      const etb = et ** bb;
      const ale = Math.log(et);

      g[i] = (bb1 * et + bb2 * et ** 2 + bb3 * et ** 3) * etb;
      gBb[i] = (bb1Bb * et + bb2Bb * et ** 2 + bb3Bb * et ** 3) * etb + g[i] * ale;

      eta[i] = et;
    }

    dsn = 0.0;
    let dsnDo = 0.0;
    let dsnUt = 0.0;

    dsnTh = 0.0;
    dsnRt = 0.0;
    dsnMs = 0.0;

    thn = 0.0;
    let thnDo = 0.0;
    let thnUt = 0.0;

    thnTh = 0.0;
    thnRt = 0.0;
    thnMs = 0.0;

    for (let i = 1; i <= n - 1; i += 1) {
      const deta = eta[i + 1] - eta[i];
      const ga = 0.5 * (g[i + 1] + g[i]);
      const gaBb = 0.5 * (gBb[i + 1] + gBb[i]);

      const uipa = 0.5 * (uip[i + 1] + uip[i]);
      const uipaDp = 0.5 * (uipDp[i + 1] + uipDp[i]);

      const uval = ut * uipa + duo * ga;
      const uDp = ut * uipaDp;

      const uDo = duoDo * ga + uDp * dpDo;
      const uUt = uipa + duoUt * ga + uDp * dpUt;
      const uBb = duo * gaBb;

      const uTh = duoTh * ga + uDp * dpTh;
      const uRt = duoRt * ga + uDp * dpRt;
      const uMs = duoMs * ga + uDp * dpMs;

      dsn += (1.0 - uval) * deta;
      dsnDo -= uDo * deta;
      dsnUt -= uUt * deta;

      dsnTh -= uTh * deta;
      dsnRt -= uRt * deta;
      dsnMs -= uMs * deta;

      thn += (uval - uval * uval) * deta;
      thnDo += (1.0 - 2.0 * uval) * uDo * deta;
      thnUt += (1.0 - 2.0 * uval) * uUt * deta;

      thnTh += (1.0 - 2.0 * uval) * uTh * deta;
      thnRt += (1.0 - 2.0 * uval) * uRt * deta;
      thnMs += (1.0 - 2.0 * uval) * uMs * deta;
    }

    const rez1 = d0 * dsn - theta * hk;
    const a11 = d0 * dsnDo + dsn;
    const a12 = d0 * dsnUt;

    const rez2 = d0 * thn - theta;
    const a21 = d0 * thnDo + thn;
    const a22 = d0 * thnUt;

    if (Math.abs(rez1 / theta) < 2.0e-5 && Math.abs(rez2 / theta) < 2.0e-5) {
      break;
    }

    const det = a11 * a22 - a12 * a21;
    b11 = a22 / det;
    b12 = -a12 / det;
    b21 = -a21 / det;
    b22 = a11 / det;

    const ddo = -(b11 * rez1 + b12 * rez2);
    const dut = -(b21 * rez1 + b22 * rez2);

    const dmax = Math.max(Math.abs(ddo / d0), Math.abs(dut / 0.05));
    let rlx = 1.0;
    if (dmax > 0.5) {
      rlx = 0.5 / dmax;
    }

    d0 += rlx * ddo;
    ut += rlx * dut;
  }

  const z1Hk = -theta;
  const z1Th = d0 * dsnTh - hk;
  const z1Rt = d0 * dsnRt;
  const z1Ms = d0 * dsnMs;

  const z2Hk = 0.0;
  const z2Th = d0 * thnTh - 1.0;
  const z2Rt = d0 * thnRt;
  const z2Ms = d0 * thnMs;

  const doHk = -(b11 * z1Hk + b12 * z2Hk);
  const doTh = -(b11 * z1Th + b12 * z2Th);
  const doRt = -(b11 * z1Rt + b12 * z2Rt);
  const doMs = -(b11 * z1Ms + b12 * z2Ms);

  const utHk = -(b21 * z1Hk + b22 * z2Hk);
  const utTh = -(b21 * z1Th + b22 * z2Th);
  const utRt = -(b21 * z1Rt + b22 * z2Rt);
  const utMs = -(b21 * z1Ms + b22 * z2Ms);

  cf = sgn * 2.0 * ut ** 2;
  const cfUt = sgn * 4.0 * ut;
  const cfDo = 0.0;

  const cfHk = cfUt * utHk + cfDo * doHk;
  const cfTh = cfUt * utTh + cfDo * doTh;
  const cfRtOut = cfUt * utRt + cfDo * doRt;
  const cfMsOut = cfUt * utMs + cfDo * doMs;

  ui = uo - duo;
  const uiUt = -duoUt;
  const uiDo = -duoDo;

  const uiHk = uiUt * utHk + uiDo * doHk;
  const uiTh = uiUt * utTh + uiDo * doTh - duoTh;
  const uiRt = uiUt * utRt + uiDo * doRt - duoRt;
  const uiMs = uiUt * utMs + uiDo * doMs - duoMs;

  const doDs = doHk / theta;
  const uiDs = uiHk / theta;
  const cfDs = cfHk / theta;

  return {
    bb,
    do: d0,
    de: d0,
    doDs,
    doTh,
    doUo: 0.0,
    doRt,
    doMs,
    ui,
    us: ui,
    uiDs,
    uiTh,
    uiUo: 0.0,
    uiRt,
    uiMs,
    hs: 0.0,
    hsDs: 0.0,
    hsTh: 0.0,
    hsUo: 0.0,
    hsRt: 0.0,
    hsMs: 0.0,
    cf,
    cfDs,
    cfTh,
    cfUo: 0.0,
    cfRt: cfRtOut,
    cfMs: cfMsOut,
    cd: 0.0,
    cdDs: 0.0,
    cdTh: 0.0,
    cdUo: 0.0,
    cdRt: 0.0,
    cdMs: 0.0,
    cdCt: 0.0,
    ct,
  };
}

// Evaluate wall-normal velocity profile for given integral quantities.
// Turbulent wall-law velocity profile construction (UWALL).
function uwall(params) {
  // Fortran comments (UWALL) highlight:
  // - Construct velocity profile using wall-law integrals.
  const th = params.th;
  const uo = params.uo;
  const d0 = params.de;
  const ui = params.us;
  const rt = params.rt;
  const cf = params.cf;
  const bb = params.bb;
  const y = params.yy;
  const u = params.xx;
  const n = params.n;

  const hpi = 1.570796327;
  const ak = 0.09;
  const vkap = 0.40;
  const vb = 5.0;

  const ebk = Math.exp(-vb * vkap);

  const sgn = sign1(cf);
  const ut = sgn * Math.sqrt(0.5 * Math.abs(cf));

  const dp = sgn * ut * rt * d0 / th;

  let upe = Math.log(dp) / vkap + vb;

  let exu = 0.0;
  let dpU = 0.0;
  for (let itup = 1; itup <= 5; itup += 1) {
    const uk = upe * vkap;
    const arg = uk - vb * vkap;
    exu = Math.exp(arg);
    const rez = upe + exu - ebk * (1.0 + uk + uk * uk / 2.0 + uk * uk * uk / 6.0) - dp;
    dpU = 1.0 + (exu - ebk * (1.0 + uk + uk * uk / 2.0)) * vkap;

    if (Math.abs(rez / dp) < 1.0e-5) {
      break;
    }

    const dupe = -rez / dpU;
    upe += dupe;
  }

  const uk = upe * vkap;
  const dpUu = (exu - ebk * (1.0 + uk)) * vkap ** 2;
  const dpUuu = (exu - ebk) * vkap ** 3;

  const upd = 1.0 / dpU;
  const updd = (-1.0 / dpU ** 3) * dpUu;

  const dc2 = 0.5 * dp * dp * updd - dp * upd;
  const dc3 = -(dp * dp * updd - dp * upd) / 3.0;

  const duo = uo - ut * (upe + dc2 + dc3);

  const bb1 = 3.0 * (bb + 2.0) * (bb + 3.0) / (bb + 7.0);
  const bb2 = -5.0 * (bb + 1.0) * (bb + 3.0) / (bb + 7.0);
  const bb3 = 2.0 * (bb + 1.0) * (bb + 2.0) / (bb + 7.0);

  const ne = n;
  const exupe = exu;
  const dexu = (exupe - ebk) / (ne - 1);

  y[1] = 0.0;
  u[1] = 0.0;
  for (let i = 2; i <= ne; i += 1) {
    const exui = ebk + (dexu - 0.75 * dexu * (ne - i) / (ne - 1)) * (i - 1);
    const uki = Math.log(exui) + vb * vkap;
    const up = uki / vkap;

    const yp = up + exui - ebk * (1.0 + uki + uki * uki / 2.0 + uki * uki * uki / 6.0);
    const ypUp = 1.0 + (exui - ebk * (1.0 + uki + uki * uki / 2.0)) * vkap;

    const et = yp / dp;

    const uip = up + dc2 * et ** 2 + dc3 * et ** 3;

    const etb = et ** bb;
    const g = (bb1 * et + bb2 * et ** 2 + bb3 * et ** 3) * etb;

    y[i] = et * d0;
    u[i] = ut * uip + duo * g;
  }

  return { y, u };
}

// Allocate Fortran-style 1-based 3D arrays for profile tables.
function create3d(n1, n2, n3) {
  const arr = new Array(n1 + 1);
  for (let i = 0; i <= n1; i += 1) {
    arr[i] = new Array(n2 + 1);
    for (let j = 0; j <= n2; j += 1) {
      arr[i][j] = new Array(n3 + 1).fill(0.0);
    }
  }
  return arr;
}

// Block tridiagonal solver for profile equations (3x3 blocks).
// Solve a 3x3 block tridiagonal system (B3SOLV).
function b3solv(A, B, C, R, n, nrhs, nrmax) {
  for (let i = 1; i <= n; i += 1) {
    const im = i - 1;
    if (i !== 1) {
      for (let k = 1; k <= 3; k += 1) {
        for (let l = 1; l <= 3; l += 1) {
          A[k][l][i] = A[k][l][i]
            - (B[k][1][i] * C[1][l][im] + B[k][2][i] * C[2][l][im] + B[k][3][i] * C[3][l][im]);
        }
        for (let l = 1; l <= nrhs; l += 1) {
          R[k][l][i] = R[k][l][i]
            - (B[k][1][i] * R[1][l][im] + B[k][2][i] * R[2][l][im] + B[k][3][i] * R[3][l][im]);
        }
      }
    }

    for (let kpiv = 1; kpiv <= 2; kpiv += 1) {
      const kp1 = kpiv + 1;
      let kx = kpiv;
      for (let k = kp1; k <= 3; k += 1) {
        if (Math.abs(A[k][kpiv][i]) > Math.abs(A[kx][kpiv][i])) {
          kx = k;
        }
      }
      if (A[kx][kpiv][i] === 0.0) {
        throw new Error(`Singular A block, i = ${i}`);
      }

      const pivot = 1.0 / A[kx][kpiv][i];
      A[kx][kpiv][i] = A[kpiv][kpiv][i];

      for (let l = kp1; l <= 3; l += 1) {
        const temp = A[kx][l][i] * pivot;
        A[kx][l][i] = A[kpiv][l][i];
        A[kpiv][l][i] = temp;
      }

      for (let l = 1; l <= 3; l += 1) {
        const temp = C[kx][l][i] * pivot;
        C[kx][l][i] = C[kpiv][l][i];
        C[kpiv][l][i] = temp;
      }

      for (let l = 1; l <= nrhs; l += 1) {
        const temp = R[kx][l][i] * pivot;
        R[kx][l][i] = R[kpiv][l][i];
        R[kpiv][l][i] = temp;
      }

      for (let k = kp1; k <= 3; k += 1) {
        for (let l = kp1; l <= 3; l += 1) {
          A[k][l][i] = A[k][l][i] - A[k][kpiv][i] * A[kpiv][l][i];
        }
        C[k][1][i] = C[k][1][i] - A[k][kpiv][i] * C[kpiv][1][i];
        C[k][2][i] = C[k][2][i] - A[k][kpiv][i] * C[kpiv][2][i];
        C[k][3][i] = C[k][3][i] - A[k][kpiv][i] * C[kpiv][3][i];
        for (let l = 1; l <= nrhs; l += 1) {
          R[k][l][i] = R[k][l][i] - A[k][kpiv][i] * R[kpiv][l][i];
        }
      }
    }

    if (A[3][3][i] === 0.0) {
      throw new Error(`Singular A block, i = ${i}`);
    }
    const pivot = 1.0 / A[3][3][i];
    C[3][1][i] *= pivot;
    C[3][2][i] *= pivot;
    C[3][3][i] *= pivot;
    for (let l = 1; l <= nrhs; l += 1) {
      R[3][l][i] *= pivot;
    }

    for (let kpiv = 2; kpiv >= 1; kpiv -= 1) {
      const kp1 = kpiv + 1;
      for (let k = kp1; k <= 3; k += 1) {
        C[kpiv][1][i] -= A[kpiv][k][i] * C[k][1][i];
        C[kpiv][2][i] -= A[kpiv][k][i] * C[k][2][i];
        C[kpiv][3][i] -= A[kpiv][k][i] * C[k][3][i];
        for (let l = 1; l <= nrhs; l += 1) {
          R[kpiv][l][i] -= A[kpiv][k][i] * R[k][l][i];
        }
      }
    }
  }

  for (let i = n - 1; i >= 1; i -= 1) {
    const ip = i + 1;
    for (let l = 1; l <= nrhs; l += 1) {
      for (let k = 1; k <= 3; k += 1) {
        R[k][l][i] = R[k][l][i]
          - (R[1][l][ip] * C[k][1][i] + R[2][l][ip] * C[k][2][i] + R[3][l][ip] * C[k][3][i]);
      }
    }
  }
}

// Laminar similarity solution helper (Falkner-Skan family).
// Falkner-Skan laminar similarity profile (FS).
function fs(params) {
  // Fortran comments (FS) highlight:
  // - Falkner-Skan similarity solution for laminar profiles.
  const inorm = params.inorm;
  const ispec = params.ispec;
  let bspec = params.bspec;
  let hspec = params.hspec;
  const n = params.n;
  const etae = params.etae;
  const geo = params.geo;
  const eta = params.eta;
  const f = params.f;
  const u = params.u;
  const s = params.s;

  const nmax = 257;
  const nrmax = 3;
  const nrhs = 3;
  const itmax = 40;

  if (n > nmax) {
    throw new Error('FS: Array overflow.');
  }

  const pi = 4.0 * Math.atan(1.0);

  let h = 0.0;
  let bu = 0.0;
  if (ispec === 1) {
    h = 2.6;
    bu = bspec;
  } else {
    h = hspec;
    if (h <= 2.17) {
      h = 2.17;
    }
    bu = (0.058 * (h - 4.0) ** 2 / (h - 1.0) - 0.068) / (6.54 * h - 14.07) * h ** 2;
    if (h > 4.0) {
      bu = Math.min(bu, 0.0);
    }
  }

  let tn = 1.0;
  if (inorm === 3) {
    tn = (6.54 * h - 14.07) / h ** 2;
  }

  let deta = 1.0;
  eta[1] = 0.0;
  for (let i = 2; i <= n; i += 1) {
    eta[i] = eta[i - 1] + deta;
    deta = geo * deta;
  }
  for (let i = 1; i <= n; i += 1) {
    eta[i] = eta[i] * etae / eta[n];
  }

  let ijoin = 1;
  if (h <= 3.0) {
    let etjoin = 5.0;
    if (inorm === 3) {
      etjoin = 7.3;
    }
    const efac = 0.5 * pi / etjoin;
    for (let i = 1; i <= n; i += 1) {
      u[i] = Math.sin(efac * eta[i]);
      f[i] = 1.0 / efac * (1.0 - Math.cos(efac * eta[i]));
      s[i] = efac * Math.cos(efac * eta[i]);
      if (eta[i] > etjoin) {
        ijoin = i;
        break;
      }
      ijoin = i;
    }
    for (let i = ijoin + 1; i <= n; i += 1) {
      u[i] = 1.0;
      f[i] = f[ijoin] + eta[i] - eta[ijoin];
      s[i] = 0.0;
    }
  } else {
    let etjoin1 = 0.0;
    let etjoin2 = 8.0;
    if (inorm === 3) {
      if (h > 4.0) {
        etjoin1 = h - 4.0;
        etjoin2 = etjoin1 + 8.0;
      }
    }
    for (let i = 1; i <= n; i += 1) {
      u[i] = 0.0;
      s[i] = 0.0;
      f[i] = 0.0;
      if (eta[i] >= etjoin1) {
        ijoin = i;
        break;
      }
      ijoin = i;
    }
    const efac = 0.5 * pi / (etjoin2 - etjoin1);
    for (let i = ijoin + 1; i <= n; i += 1) {
      const ebar = eta[i] - etjoin1;
      u[i] = 0.5 - 0.5 * Math.cos(2.0 * efac * ebar);
      f[i] = 0.5 * ebar - 0.25 / efac * Math.sin(2.0 * efac * ebar);
      s[i] = efac * Math.sin(2.0 * efac * ebar);
      if (eta[i] >= etjoin2) {
        ijoin = i;
        break;
      }
      ijoin = i;
    }
    for (let i = ijoin + 1; i <= n; i += 1) {
      u[i] = 1.0;
      f[i] = f[ijoin] + eta[i] - eta[ijoin];
      s[i] = 0.0;
    }
  }

  let rms = 1.0;
  const A = create3d(3, 3, n);
  const B = create3d(3, 3, n);
  const C = create3d(3, 3, n);
  const R = create3d(3, nrmax, n);

  for (let iter = 1; iter <= itmax; iter += 1) {
    for (let i = 1; i <= n; i += 1) {
      for (let ii = 1; ii <= 3; ii += 1) {
        for (let iii = 1; iii <= 3; iii += 1) {
          A[ii][iii][i] = 0.0;
          B[ii][iii][i] = 0.0;
          C[ii][iii][i] = 0.0;
        }
        R[ii][1][i] = 0.0;
        R[ii][2][i] = 0.0;
        R[ii][3][i] = 0.0;
      }
    }

    let thi = 0.0;
    for (let i = 1; i <= n - 1; i += 1) {
      const us = u[i] + u[i + 1];
      const detai = eta[i + 1] - eta[i];
      thi += (1.0 - 0.5 * us) * 0.5 * us * detai;
    }

    A[1][1][1] = 1.0;
    A[2][2][1] = 1.0;
    A[3][2][n] = 1.0;
    R[1][1][1] = f[1];
    R[2][1][1] = u[1];
    R[3][1][n] = u[n] - 1.0;

    let betu = 0.0;
    let betuBu = 0.0;
    let betn = 0.0;
    let betnBu = 0.0;
    if (inorm === 2) {
      betu = 2.0 * bu / (bu + 1.0);
      betuBu = (2.0 - betu / (bu + 1.0)) / (bu + 1.0);
      betn = 1.0;
      betnBu = 0.0;
    } else {
      betu = bu;
      betuBu = 1.0;
      betn = 0.5 * (1.0 + bu);
      betnBu = 0.5;
    }

    for (let i = 1; i <= n - 1; i += 1) {
      const detai = eta[i + 1] - eta[i];
      R[1][1][i + 1] = f[i + 1] - f[i] - 0.5 * detai * (u[i + 1] + u[i]);
      R[2][1][i + 1] = u[i + 1] - u[i] - 0.5 * detai * (s[i + 1] + s[i]);
      R[3][1][i] = s[i + 1] - s[i]
        + tn * (betn * detai * 0.5 * (f[i + 1] * s[i + 1] + f[i] * s[i])
          + betu * detai * (1.0 - 0.5 * (u[i + 1] ** 2 + u[i] ** 2)));

      A[3][1][i] = tn * betn * 0.5 * detai * s[i];
      C[3][1][i] = tn * betn * 0.5 * detai * s[i + 1];
      A[3][2][i] = -tn * betu * detai * u[i];
      C[3][2][i] = -tn * betu * detai * u[i + 1];
      A[3][3][i] = tn * betn * 0.5 * detai * f[i] - 1.0;
      C[3][3][i] = tn * betn * 0.5 * detai * f[i + 1] + 1.0;

      B[1][1][i + 1] = -1.0;
      A[1][1][i + 1] = 1.0;
      B[1][2][i + 1] = -0.5 * detai;
      A[1][2][i + 1] = -0.5 * detai;

      B[2][2][i + 1] = -1.0;
      A[2][2][i + 1] = 1.0;
      B[2][3][i + 1] = -0.5 * detai;
      A[2][3][i + 1] = -0.5 * detai;

      R[3][2][i] = tn * (betnBu * detai * 0.5 * (f[i + 1] * s[i + 1] + f[i] * s[i])
        + betuBu * detai * (1.0 - 0.5 * (u[i + 1] ** 2 + u[i] ** 2)));
      R[3][3][i] = betn * detai * 0.5 * (f[i + 1] * s[i + 1] + f[i] * s[i])
        + betu * detai * (1.0 - 0.5 * (u[i + 1] ** 2 + u[i] ** 2));
    }

    b3solv(A, B, C, R, n, nrhs, nrmax);

    let dsi = 0.0;
    let dsi1 = 0.0;
    let dsi2 = 0.0;
    let dsi3 = 0.0;

    thi = 0.0;
    let thi1 = 0.0;
    let thi2 = 0.0;
    let thi3 = 0.0;

    for (let i = 1; i <= n - 1; i += 1) {
      const us = u[i] + u[i + 1];
      const detai = eta[i + 1] - eta[i];

      dsi += (1.0 - 0.5 * us) * detai;
      const dsiUs = -0.5 * detai;

      thi += (1.0 - 0.5 * us) * 0.5 * us * detai;
      const thiUs = (0.5 - 0.5 * us) * detai;

      dsi1 += dsiUs * (R[2][1][i] + R[2][1][i + 1]);
      dsi2 += dsiUs * (R[2][2][i] + R[2][2][i + 1]);
      dsi3 += dsiUs * (R[2][3][i] + R[2][3][i + 1]);

      thi1 += thiUs * (R[2][1][i] + R[2][1][i + 1]);
      thi2 += thiUs * (R[2][2][i] + R[2][2][i + 1]);
      thi3 += thiUs * (R[2][3][i] + R[2][3][i + 1]);
    }

    let r1 = 0.0;
    let q11 = 0.0;
    let q12 = 0.0;
    if (ispec === 1) {
      r1 = bspec - bu;
      q11 = 1.0;
      q12 = 0.0;
    } else {
      r1 = dsi - hspec * thi - dsi1 + hspec * thi1;
      q11 = -dsi2 + hspec * thi2;
      q12 = -dsi3 + hspec * thi3;
    }

    let r2 = 0.0;
    let q21 = 0.0;
    let q22 = 0.0;
    if (inorm === 3) {
      r2 = thi - 1.0 - thi1;
      q21 = -thi2;
      q22 = -thi3;
    } else {
      r2 = tn - 1.0;
      q21 = 0.0;
      q22 = 1.0;
    }

    const det = q11 * q22 - q12 * q21;
    const dbu = -(r1 * q22 - q12 * r2) / det;
    const dtn = -(q11 * r2 - r1 * q21) / det;

    let rmax = 0.0;
    rms = 0.0;
    for (let i = 1; i <= n; i += 1) {
      const df = -R[1][1][i] - dbu * R[1][2][i] - dtn * R[1][3][i];
      const du = -R[2][1][i] - dbu * R[2][2][i] - dtn * R[2][3][i];
      const ds = -R[3][1][i] - dbu * R[3][2][i] - dtn * R[3][3][i];
      rmax = Math.max(rmax, Math.abs(df), Math.abs(du), Math.abs(ds));
      rms = df ** 2 + du ** 2 + ds ** 2 + rms;
    }
    rms = Math.sqrt(rms / (3.0 * n + 3.0));

    rmax = Math.max(rmax, Math.abs(dbu / 1.0), Math.abs(dtn / tn));

    let rlx = 1.0;
    if (rmax > 0.5) {
      rlx = 0.5 / rmax;
    }

    for (let i = 1; i <= n; i += 1) {
      const df = -R[1][1][i] - dbu * R[1][2][i] - dtn * R[1][3][i];
      const du = -R[2][1][i] - dbu * R[2][2][i] - dtn * R[2][3][i];
      const ds = -R[3][1][i] - dbu * R[3][2][i] - dtn * R[3][3][i];
      f[i] += rlx * df;
      u[i] += rlx * du;
      s[i] += rlx * ds;
    }

    bu += rlx * dbu;
    tn += rlx * dtn;

    if (iter > 3 && rms < 1.0e-6) {
      break;
    }
  }

  hspec = 0.0;
  bspec = bu;
  let dsi = 0.0;
  let thi = 0.0;
  for (let i = 1; i <= n - 1; i += 1) {
    const us = u[i] + u[i + 1];
    const detai = eta[i + 1] - eta[i];
    dsi += (1.0 - 0.5 * us) * detai;
    thi += (1.0 - 0.5 * us) * 0.5 * us * detai;
  }
  hspec = dsi / thi;

  const delta = Math.sqrt(tn);
  return {
    bspec,
    hspec,
    delta,
  };
}

export { prwall, uwall, fs, b3solv };
