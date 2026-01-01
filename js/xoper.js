// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

// Port of xoper.f viscous iteration driver.
// Orchestrates the viscous Newton loop for the webapp.

import {
  setbl,
  update as blUpdate,
  mrcl,
  comset,
  blpini,
  iblsys,
} from './xbl.js';
import { blsolv } from './xsolve.js';
import { seval } from './spline.js';
import { lefind } from './xgeom.js';
import { clcalc, cdcalc, cpcalc } from './xfoil.js';
import {
  xywake,
  qwcalc,
  qdcalc,
  stfind,
  iblpan,
  xicalc,
  uicalc,
  qvfue,
  qiset,
  gamqv,
  stmove,
} from './xpanel.js';
import {
  NCOM,
  ensureCtx,
  syncComToVars,
  syncVarsToCom,
  blprv,
  blkin,
  blvar,
  blmid,
  trchek,
  tesys,
  blsys,
  hkin,
} from './xblsys.js';
import { createMatrix1, createTensor3Flat } from './arrays.js';

function formatFixed(value, width, decimals) {
  const text = Number.isFinite(value) ? value.toFixed(decimals) : String(value);
  if (width <= 0 || text.length >= width) {
    return text;
  }
  return text.padStart(width, ' ');
}

function bstrip(line) {
  return line.replace(/ /g, '');
}

function cpdump(ctxPanel, blCtx = null, opts = {}) {
  const kdelim = Number.isFinite(opts.kdelim) ? opts.kdelim : (Number.isFinite(ctxPanel?.KDELIM) ? ctxPanel.KDELIM : 1);
  let delim = ' ';
  if (kdelim === 1) {
    delim = ',';
  } else if (kdelim === 2) {
    delim = '\t';
  } else if (kdelim !== 0) {
    console.warn('? Illegal delimiter.  Using blank.');
    delim = ' ';
  }

  const n = ctxPanel?.N ?? blCtx?.N ?? 0;
  const qinf = Number.isFinite(ctxPanel?.QINF) ? ctxPanel.QINF : (blCtx?.QINF ?? 1.0);
  const minf = Number.isFinite(blCtx?.MINF) ? blCtx.MINF : (Number.isFinite(ctxPanel?.MINF) ? ctxPanel.MINF : 0.0);

  const lines = [];
  if (kdelim === 0) {
    lines.push('#      x          Cp  ');
  } else {
    lines.push(`#x${delim}Cp`);
  }

  const beta = Math.sqrt(1.0 - minf ** 2);
  const bfac = 0.5 * minf ** 2 / (1.0 + beta);
  for (let i = 1; i <= n; i += 1) {
    const gam = ctxPanel?.GAM?.[i - 1] ?? 0.0;
    const cpinc = 1.0 - (gam / qinf) ** 2;
    const den = beta + bfac * cpinc;
    const cpcom = cpinc / den;
    if (kdelim === 0) {
      const line = ` ${formatFixed(ctxPanel?.X?.[i - 1] ?? blCtx?.X?.[i] ?? 0.0, 11, 5)}${formatFixed(cpcom, 11, 5)}`;
      lines.push(line);
    } else {
      const line = ` ${formatFixed(ctxPanel?.X?.[i - 1] ?? blCtx?.X?.[i] ?? 0.0, 11, 5)}${delim}${formatFixed(cpcom, 11, 5)}${delim}`;
      lines.push(bstrip(line));
    }
  }

  return `${lines.join('\n')}\n`;
}

function bldump(ctxPanel, blCtx = null, opts = {}) {
  const kdelim = Number.isFinite(opts.kdelim) ? opts.kdelim : (Number.isFinite(ctxPanel?.KDELIM) ? ctxPanel.KDELIM : 1);
  let delim = ' ';
  if (kdelim === 1) {
    delim = ',';
  } else if (kdelim === 2) {
    delim = '\t';
  } else if (kdelim !== 0) {
    console.warn('? Illegal delimiter.  Using blank.');
    delim = ' ';
  }

  const n = ctxPanel?.N ?? blCtx?.N ?? 0;
  const nw = blCtx?.NW ?? ctxPanel?.NW ?? 0;
  const qinf = Number.isFinite(ctxPanel?.QINF) ? ctxPanel.QINF : (blCtx?.QINF ?? 1.0);
  const minf = Number.isFinite(blCtx?.MINF) ? blCtx.MINF : (Number.isFinite(ctxPanel?.MINF) ? ctxPanel.MINF : 0.0);
  const gamm1 = Number.isFinite(blCtx?.GAMM1) ? blCtx.GAMM1 : 0.4;
  const beta = Math.sqrt(1.0 - minf ** 2);
  const tklam = minf ** 2 / (1.0 + beta) ** 2;

  const lines = [];
  if (kdelim === 0) {
    lines.push(
      '#    s        x        y     Ue/Vinf    Dstar     Theta      Cf       H       H*        P         m          K          tau         Di',
    );
  } else {
    lines.push(`#s${delim}x${delim}y${delim}Ue/Vinf${delim}Dstar${delim}Theta${delim}Cf${delim}H`);
  }

  const hstinv = gamm1 * (minf / qinf) ** 2 / (1.0 + 0.5 * gamm1 * minf ** 2);

  for (let i = 1; i <= n; i += 1) {
    const gam = ctxPanel?.GAM?.[i - 1] ?? 0.0;
    const is = gam < 0.0 ? 2 : 1;
    let ds = 0.0;
    let th = 0.0;
    let ts = 0.0;
    let cf = 0.0;
    let h = 1.0;
    let hs = 2.0;
    if (blCtx?.LIPAN && blCtx?.LVISC) {
      let ibl;
      if (is === 1) {
        ibl = blCtx.IBLTE[is] - i + 1;
      } else {
        ibl = blCtx.IBLTE[is] + i - n;
      }
      ds = blCtx.DSTR[ibl][is];
      th = blCtx.THET[ibl][is];
      ts = blCtx.TSTR[ibl][is];
      cf = blCtx.TAU[ibl][is] / (0.5 * qinf ** 2);
      if (th === 0.0) {
        h = 1.0;
        hs = 1.0;
      } else {
        h = ds / th;
        hs = ts / th;
      }
    }
    const ue = (gam / qinf) * (1.0 - tklam) / (1.0 - tklam * (gam / qinf) ** 2);
    const amsq = ue * ue * hstinv / (gamm1 * (1.0 - 0.5 * ue * ue * hstinv));
    const hk = hkin(h, amsq).hk;

    if (kdelim === 0) {
      const line = ` ${formatFixed(blCtx?.S?.[i] ?? ctxPanel?.S?.[i - 1] ?? 0.0, 9, 5)}${formatFixed(blCtx?.X?.[i] ?? ctxPanel?.X?.[i - 1] ?? 0.0, 9, 5)}${formatFixed(blCtx?.Y?.[i] ?? ctxPanel?.Y?.[i - 1] ?? 0.0, 9, 5)}${formatFixed(ue, 9, 5)}${formatFixed(ds, 10, 6)}${formatFixed(th, 10, 6)}${formatFixed(cf, 10, 6)}${formatFixed(hk, 10, 4)}${formatFixed(hs, 10, 4)}${formatFixed(th * ue ** 2, 9, 5)}${formatFixed(ds * ue, 9, 5)}${formatFixed(ts * ue ** 3, 9, 5)}`;
      lines.push(line);
    } else {
      const line = ` ${formatFixed(blCtx?.S?.[i] ?? ctxPanel?.S?.[i - 1] ?? 0.0, 9, 5)}${delim}${formatFixed(blCtx?.X?.[i] ?? ctxPanel?.X?.[i - 1] ?? 0.0, 9, 5)}${delim}${formatFixed(blCtx?.Y?.[i] ?? ctxPanel?.Y?.[i - 1] ?? 0.0, 9, 5)}${delim}${formatFixed(ue, 9, 5)}${delim}${formatFixed(ds, 10, 6)}${delim}${formatFixed(th, 10, 6)}${delim}${formatFixed(cf, 10, 6)}${delim}${formatFixed(hk, 10, 4)}`;
      lines.push(bstrip(line));
    }
  }

  if (blCtx?.LWAKE) {
    const is = 2;
    for (let i = n + 1; i <= n + nw; i += 1) {
      const ibl = blCtx.IBLTE[is] + i - n;
      const ds = blCtx.DSTR[ibl][is];
      const th = blCtx.THET[ibl][is];
      const h = ds / th;
      const cf = 0.0;
      const ui = blCtx.UEDG[ibl][is];
      const ue = (ui / qinf) * (1.0 - tklam) / (1.0 - tklam * (ui / qinf) ** 2);
      const amsq = ue * ue * hstinv / (gamm1 * (1.0 - 0.5 * ue * ue * hstinv));
      const hk = hkin(h, amsq).hk;

      if (kdelim === 0) {
        const line = ` ${formatFixed(blCtx.S?.[i] ?? ctxPanel?.S?.[i - 1] ?? 0.0, 9, 5)}${formatFixed(blCtx.X?.[i] ?? ctxPanel?.X?.[i - 1] ?? 0.0, 9, 5)}${formatFixed(blCtx.Y?.[i] ?? ctxPanel?.Y?.[i - 1] ?? 0.0, 9, 5)}${formatFixed(ue, 9, 5)}${formatFixed(ds, 10, 6)}${formatFixed(th, 10, 6)}${formatFixed(cf, 10, 6)}${formatFixed(hk, 10, 4)}`;
        lines.push(line);
      } else {
        const line = ` ${formatFixed(blCtx.S?.[i] ?? ctxPanel?.S?.[i - 1] ?? 0.0, 9, 5)}${delim}${formatFixed(blCtx.X?.[i] ?? ctxPanel?.X?.[i - 1] ?? 0.0, 9, 5)}${delim}${formatFixed(blCtx.Y?.[i] ?? ctxPanel?.Y?.[i - 1] ?? 0.0, 9, 5)}${delim}${formatFixed(ue, 9, 5)}${delim}${formatFixed(ds, 10, 6)}${delim}${formatFixed(th, 10, 6)}${delim}${formatFixed(cf, 10, 6)}${delim}${formatFixed(hk, 10, 4)}`;
        lines.push(bstrip(line));
      }
    }
  }

  return `${lines.join('\n')}\n`;
}

// Initialize BL operating point (alpha, Mach, Re) to match XFOIL's OPER init.
// Mirrors OPER/SPECAL setup in xoper.f.
function applyXfoilOperInit(blCtx, alphaRad, reinf1, opts = {}) {
  const { resetConvergence = true } = opts;
  blCtx.LVISC = true;
  if (blCtx.LVCONV) {
    const da = Math.abs((blCtx.AVISC ?? NaN) - alphaRad);
    const dm = Math.abs((blCtx.MVISC ?? NaN) - (blCtx.MINF ?? NaN));
    if ((Number.isFinite(da) && da > 1.0e-5) || (Number.isFinite(dm) && dm > 1.0e-5)) {
      blCtx.LVCONV = false;
    }
  }
  if (resetConvergence) {
    blCtx.LVCONV = false;
  }
  blCtx.LWAKE = blCtx.NW > 0;
  blCtx.LALFA = true;
  blCtx.ALFA = alphaRad;
  blCtx.ADEG = alphaRad / blCtx.DTOR;
  blCtx.AWAKE = alphaRad;
  blCtx.MINF1 = blCtx.MINF;
  blCtx.REINF1 = reinf1;
  blCtx.REINF = reinf1;
  blCtx.MINF_CL = 0.0;
  blCtx.REINF_CL = 0.0;
}

// Reinitialize BL state for retry without changing the governing equations.
function resetBlState(blCtx) {
  blCtx.LBLINI = false;
  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 1; ibl <= blCtx.NBL[is]; ibl += 1) {
      blCtx.CTAU[ibl][is] = 0.0;
      blCtx.THET[ibl][is] = 0.0;
      blCtx.DSTR[ibl][is] = 0.0;
      blCtx.MASS[ibl][is] = 0.0;
      blCtx.TAU[ibl][is] = 0.0;
      blCtx.DIS[ibl][is] = 0.0;
      blCtx.CTQ[ibl][is] = 0.0;
      blCtx.DELT[ibl][is] = 0.0;
      blCtx.TSTR[ibl][is] = 0.0;
      const uinv = blCtx.UINV[ibl][is];
      const ue = Number.isFinite(uinv) && uinv > 0.0 ? uinv : 1.0e-6;
      blCtx.UEDG[ibl][is] = ue;
    }
  }
}

function ensureViscousArrays(blCtx) {
  const nsys = blCtx.NSYS;
  if (!Number.isFinite(nsys) || nsys <= 0) return;

  const tensorSize = (3 + 1) * (2 + 1) * (nsys + 1);
  const vmSize = (3 + 1) * (nsys + 1) * (nsys + 1);
  const needsTensor = !blCtx.VAF
    || blCtx.VAF.length < tensorSize
    || !blCtx.VA
    || blCtx.VA.length <= 3
    || blCtx.VA[1]?.length <= 2
    || blCtx.VA[1]?.[1]?.length <= nsys;
  const needsVm = !blCtx.VMF
    || blCtx.VMF.length < vmSize
    || !blCtx.VM
    || blCtx.VM.length <= 3
    || blCtx.VM[1]?.length <= nsys
    || blCtx.VM[1]?.[1]?.length <= nsys;

  if (needsTensor) {
    const va = createTensor3Flat(3, 2, nsys);
    const vb = createTensor3Flat(3, 2, nsys);
    const vdel = createTensor3Flat(3, 2, nsys);
    blCtx.VAF = va.flat;
    blCtx.VBF = vb.flat;
    blCtx.VDELF = vdel.flat;
    blCtx.VA = va.view;
    blCtx.VB = vb.view;
    blCtx.VDEL = vdel.view;
  }
  if (needsVm) {
    const vm = createTensor3Flat(3, nsys, nsys);
    blCtx.VMF = vm.flat;
    blCtx.VM = vm.view;
  }
  if (!blCtx.VZ || blCtx.VZ.length <= 3) {
    blCtx.VZ = createMatrix1(3, 2);
  }
  if (!blCtx.XOCTR || blCtx.XOCTR.length < 3) {
    blCtx.XOCTR = new Float64Array(3);
  }
  if (!blCtx.YOCTR || blCtx.YOCTR.length < 3) {
    blCtx.YOCTR = new Float64Array(3);
  }
  if (!blCtx.TINDEX || blCtx.TINDEX.length < 3) {
    blCtx.TINDEX = new Float64Array(3);
  }
}

// Allocate/initialize BL arrays and geometry mapping used by viscous solver.
function buildBlContext(nb, ctxPanel, ncr) {
  const nw = ctxPanel.NW ?? 0;
  const total = nb + nw;
  const nbl = total + 1;
  const ctxBl = {
    N: nb,
    NW: nw,
    IQX: nb,
    X: new Float64Array(total + 1),
    Y: new Float64Array(total + 1),
    S: new Float64Array(total + 1),
    NBL: new Int32Array(3),
    IBLTE: new Int32Array(3),
    ISYS: createMatrix1(nbl, 2),
    IPAN: createMatrix1(nbl, 2),
    VTI: createMatrix1(nbl, 2),
    XSSI: createMatrix1(nbl, 2),
    UEDG: createMatrix1(nbl, 2),
    UINV: createMatrix1(nbl, 2),
    UINV_A: createMatrix1(nbl, 2),
    THET: createMatrix1(nbl, 2),
    DSTR: createMatrix1(nbl, 2),
    CTAU: createMatrix1(nbl, 2),
    MASS: createMatrix1(nbl, 2),
    TAU: createMatrix1(nbl, 2),
    DIS: createMatrix1(nbl, 2),
    CTQ: createMatrix1(nbl, 2),
    DELT: createMatrix1(nbl, 2),
    TSTR: createMatrix1(nbl, 2),
    WGAP: new Float64Array(nw + 1),
    ACRIT: new Float64Array(3),
    XSTRIP: new Float64Array(3),
    ITRAN: new Int32Array(3),
    TFORCE: new Array(3).fill(false),
    XSSITR: new Float64Array(3),
    W1: new Float64Array(nb + 1),
    W2: new Float64Array(nb + 1),
    W3: new Float64Array(nb + 1),
    W4: new Float64Array(nb + 1),
    DTOR: Math.PI / 180.0,
    ANTE: ctxPanel.ANTE ?? 0.0,
    ASTE: ctxPanel.ASTE ?? 0.0,
    DSTE: ctxPanel.DSTE ?? 0.0,
    GAMBL: 1.4,
    GM1BL: 0.4,
    GAMM1: 0.4,
    MINF: 0.0,
    QINF: 1.0,
    QINFBL: 1.0,
    TKBL: 0.0,
    TKBL_MS: 0.0,
    RSTBL: 1.0,
    RSTBL_MS: 0.0,
    REYBL: 1.0e6,
    REYBL_RE: 1.0,
    REYBL_MS: 0.0,
    HSTINV: 0.0,
    HSTINV_MS: 0.0,
    HVRAT: 0.0,
    IDAMPV: 0,
    NCOM,
    blprv,
    blkin,
    blvar,
    blmid,
    trchek,
    tesys,
    blsys,
    hkin,
    syncComToVars,
    syncVarsToCom,
    ensureViscousArrays,
    LIPAN: false,
  };

  for (let i = 1; i <= total; i += 1) {
    ctxBl.X[i] = ctxPanel.X[i - 1];
    ctxBl.Y[i] = ctxPanel.Y[i - 1];
    ctxBl.S[i] = ctxPanel.S[i - 1];
  }

  const sLE = lefind(ctxPanel.X, ctxPanel.XP, ctxPanel.Y, ctxPanel.YP, ctxPanel.S, nb);
  ctxBl.SLE = sLE;
  ctxBl.XLE = seval(sLE, ctxPanel.X, ctxPanel.XP, ctxPanel.S, nb);
  ctxBl.YLE = seval(sLE, ctxPanel.Y, ctxPanel.YP, ctxPanel.S, nb);
  ctxBl.XTE = 0.5 * (ctxBl.X[1] + ctxBl.X[nb]);
  ctxBl.YTE = 0.5 * (ctxBl.Y[1] + ctxBl.Y[nb]);

  ctxBl.IVX = nbl;

  const acrit = Number.isFinite(ncr) ? ncr : 9.0;
  ctxBl.ACRIT[1] = acrit;
  ctxBl.ACRIT[2] = acrit;
  ctxBl.XSTRIP[1] = 1.0;
  ctxBl.XSTRIP[2] = 1.0;

  ensureCtx(ctxBl);
  blpini(ctxBl);

  ctxBl.XP = new Float64Array(nb + 1);
  ctxBl.YP = new Float64Array(nb + 1);
  for (let i = 1; i <= nb; i += 1) {
    ctxBl.XP[i] = ctxPanel.XP[i - 1];
    ctxBl.YP[i] = ctxPanel.YP[i - 1];
  }

  ctxBl.VACCEL = 0.01;
  ctxBl.XOCTR = new Float64Array(3);
  ctxBl.YOCTR = new Float64Array(3);
  ctxBl.TINDEX = new Float64Array(3);
  ctxBl.XOCTR[1] = 1.0;
  ctxBl.XOCTR[2] = 1.0;
  ctxBl.YOCTR[1] = 0.0;
  ctxBl.YOCTR[2] = 0.0;
  ctxBl.CL = 0.0;
  ctxBl.CLSPEC = 0.0;
  ctxBl.LALFA = true;
  ctxBl.IDAMP = 0;
  ctxBl.LBLINI = false;
  ctxBl.GAMMA = 1.4;
  ctxBl.GAMM1 = 0.4;
  ctxBl.MINF1 = ctxBl.MINF;
  ctxBl.REINF1 = ctxBl.REYBL;
  ctxBl.MATYP = 1;
  ctxBl.RETYP = 1;
  ctxBl.MINF_CL = 0.0;
  ctxBl.REINF_CL = 0.0;
  return ctxBl;
}

function transitionXc(blCtx, ctxPanel, is) {
  const xoc = blCtx.XOCTR?.[is];
  if (Number.isFinite(xoc) && xoc > 0.0) {
    return xoc;
  }
  const xt = blCtx.XSSITR?.[is] ?? 0.0;
  if (!Number.isFinite(xt) || xt <= 0.0) return NaN;
  const str = is === 1 ? blCtx.SST - xt : blCtx.SST + xt;
  const xtr = seval(str, ctxPanel.X, ctxPanel.XP, ctxPanel.S, ctxPanel.N);
  const ytr = seval(str, ctxPanel.Y, ctxPanel.YP, ctxPanel.S, ctxPanel.N);
  const chx = blCtx.XTE - blCtx.XLE;
  const chy = blCtx.YTE - blCtx.YLE;
  const chsq = chx * chx + chy * chy || 1.0;
  return ((xtr - blCtx.XLE) * chx + (ytr - blCtx.YLE) * chy) / chsq;
}

// VISCAL prelude: wake setup and QISET.
function initViscousPanel(ctxPanel) {
  if (!ctxPanel.LWAKE) {
    xywake(ctxPanel);
  }

  qwcalc(ctxPanel);
  return qiset(ctxPanel);
}

// VISCAL prelude: BL geometry mapping, Ue init, and DIJ build.
function initViscousBl(blCtx, ctxPanel, qinv, qinvA) {
  if (!blCtx.LIPAN) {
    if (blCtx.LBLINI) {
      const total = ctxPanel.N + (ctxPanel.NW ?? 0);
      if (ctxPanel.QVIS && ctxPanel.QVIS.length === total + 1) {
        gamqv(ctxPanel, ctxPanel.QVIS, qinvA);
      }
    }

    const { ist, sst, sstGo, sstGp } = stfind(ctxPanel, ctxPanel.N);
    blCtx.IST = ist + 1;
    blCtx.SST = sst;
    blCtx.SST_GO = sstGo;
    blCtx.SST_GP = sstGp;

    iblpan(blCtx, ctxPanel.N, ctxPanel.NW ?? 0);
    xicalc(blCtx, ctxPanel);
    iblsys(blCtx);
    ensureViscousArrays(blCtx);
  }

  if (!blCtx.VA || !blCtx.VB || !blCtx.VDEL || !blCtx.VM) {
    if (!Number.isFinite(blCtx.NSYS) || blCtx.NSYS <= 0) {
      iblsys(blCtx);
    }
    ensureViscousArrays(blCtx);
  }

  uicalc(blCtx, qinv, qinvA);

  if (!blCtx.LBLINI) {
    for (let is = 1; is <= 2; is += 1) {
      for (let ibl = 1; ibl <= blCtx.NBL[is]; ibl += 1) {
        blCtx.UEDG[ibl][is] = blCtx.UINV[ibl][is];
      }
    }
  }

}

// Iterative viscous solve with optional re-init (VISCAL).
function solveViscous(
  blCtx,
  ctxPanel,
  qinv,
  qinvA,
  maxIter = 20,
) {
  // Fortran comments (VISCAL) highlight:
  // - Iterate BL/panel coupling to convergence.
  const total = ctxPanel.N + (ctxPanel.NW ?? 0);
  const qvis = new Float64Array(total + 1);
  const eps1 = 1.0e-4;

  let converged = false;
  console.log('');
  console.log('Solving BL system ...');
  console.log(' iter        rms        max  var at   rlx        a        CL        Cm        CD       CDf       CDp     tr1     tr2');
  if (!blCtx.LBLINI) {
    for (let ibl = 1; ibl <= blCtx.NBL[1]; ibl += 1) {
      blCtx.UEDG[ibl][1] = blCtx.UINV[ibl][1];
    }
    for (let ibl = 1; ibl <= blCtx.NBL[2]; ibl += 1) {
      blCtx.UEDG[ibl][2] = blCtx.UINV[ibl][2];
    }
  }

  for (let iter = 0; iter < maxIter; iter += 1) {
    // Assemble BL system and solve for Newton deltas.
    setbl(blCtx);
    const ok = blsolv(blCtx);
    if (!ok) break;
    // Update BL state with Newton step.
    blUpdate(blCtx);

    if (blCtx.LALFA) {
      mrcl(blCtx, blCtx.CL);
      comset(blCtx);
    } else {
      const qisetRes = qiset(ctxPanel);
      qinv = qisetRes.qinv;
      qinvA = qisetRes.qinvA;
      uicalc(blCtx, qinv, qinvA);
    }

    qvfue(blCtx, qvis);
    gamqv(ctxPanel, qvis, qinvA);
    stmove(ctxPanel, blCtx, qinv, qinvA);

    const coeffs = clcalc(
      ctxPanel.N,
      ctxPanel.X,
      ctxPanel.Y,
      ctxPanel.GAM,
      ctxPanel.GAM_A,
      blCtx.ALFA,
      blCtx.MINF,
      ctxPanel.QINF ?? 1.0,
    );
    const drag = cdcalc(ctxPanel, blCtx, blCtx.ALFA, ctxPanel.QINF ?? 1.0);
    blCtx.CL = coeffs.cl;
    blCtx.CM = coeffs.cm;
    blCtx.CDP = coeffs.cdp;
    blCtx.CD = drag.cd;
    blCtx.CDF = drag.cdf;

    const tr1 = transitionXc(blCtx, ctxPanel, 1);
    const tr2 = transitionXc(blCtx, ctxPanel, 2);
    const rlx = Number.isFinite(blCtx.RLX) ? blCtx.RLX : 1.0;
    const cdp = drag.cd - drag.cdf;
    const rawVmx = blCtx.VMXBL ?? ' ';
    const vmx = rawVmx.toString().trim().length ? rawVmx : '?';
    const imx = blCtx.IMXBL ?? 0;
    const ismx = blCtx.ISMXBL ?? 0;
    const iterLine = `${String(iter + 1).padStart(5)}`
      + `${Number.isFinite(blCtx.RMSBL) ? blCtx.RMSBL.toExponential(3).padStart(12) : '     NaN'.padStart(12)}`
      + `${Number.isFinite(blCtx.RMXBL) ? blCtx.RMXBL.toExponential(3).padStart(12) : '     NaN'.padStart(12)}`
      + `${vmx.toString().padStart(4)} at ${String(imx).padStart(4)}${String(ismx).padStart(3)}`
      + `${rlx.toFixed(2).padStart(7)}`
      + `${(blCtx.ALFA / blCtx.DTOR).toFixed(3).padStart(9)}`
      + `${coeffs.cl.toFixed(5).padStart(10)}`
      + `${coeffs.cm.toFixed(5).padStart(10)}`
      + `${drag.cd.toFixed(6).padStart(10)}`
      + `${drag.cdf.toFixed(6).padStart(10)}`
      + `${cdp.toFixed(6).padStart(10)}`
      + `${Number.isFinite(tr1) ? tr1.toFixed(4).padStart(8) : '   NaN'.padStart(8)}`
      + `${Number.isFinite(tr2) ? tr2.toFixed(4).padStart(8) : '   NaN'.padStart(8)}`;
    console.log(iterLine);

    // Convergence check (RMS residual).
    if (Number.isFinite(blCtx.RMSBL) && blCtx.RMSBL < eps1) {
      blCtx.LVCONV = true;
      blCtx.AVISC = blCtx.ALFA;
      blCtx.MVISC = blCtx.MINF;
      converged = true;
      break;
    }
  }

  if (!converged) {
    console.log('VISCAL:  Convergence failed');
  }

  return { qvis, converged };
}

// Converge to specified alpha for inviscid solutions (SPECAL).
function specal(ctxPanel, alphaRad) {
  return qiset(ctxPanel, alphaRad);
}

// Converge viscous solution for current alpha (VISCAL).
function viscal(
  blCtx,
  ctxPanel,
  alphaRad,
  reinf,
  opts = {},
) {
  const {
    maxIter = 20,
    logSurface = false,
    reuseSolution = false,
  } = opts;

  const { qinv, qinvA } = initViscousPanel(ctxPanel);
  applyXfoilOperInit(blCtx, alphaRad, Number.isFinite(reinf) ? reinf : 1.0e6, {
    resetConvergence: !reuseSolution,
  });
  initViscousBl(blCtx, ctxPanel, qinv, qinvA);
  if (blCtx.LALFA) {
    const coeffs = clcalc(
      ctxPanel.N,
      ctxPanel.X,
      ctxPanel.Y,
      ctxPanel.GAM,
      ctxPanel.GAM_A,
      blCtx.ALFA,
      blCtx.MINF,
      ctxPanel.QINF ?? 1.0,
    );
    blCtx.CL = coeffs.cl;
    blCtx.CM = coeffs.cm;
    blCtx.CDP = coeffs.cdp;
    blCtx.CL_ALF = coeffs.clAlf;
    blCtx.CL_MSQ = coeffs.clMsq;
  }

  if (blCtx.LVCONV) {
    const total = ctxPanel.N + (ctxPanel.NW ?? 0);
    const qvis = ctxPanel.QVIS && ctxPanel.QVIS.length === total + 1
      ? ctxPanel.QVIS
      : new Float64Array(total + 1);
    qvfue(blCtx, qvis);
    cpcalc(qvis, ctxPanel.QINF ?? 1.0, blCtx.MINF);
    cpcalc(qinv, ctxPanel.QINF ?? 1.0, blCtx.MINF);
    gamqv(ctxPanel, qvis, qinvA);
    const coeffs = clcalc(
      ctxPanel.N,
      ctxPanel.X,
      ctxPanel.Y,
      ctxPanel.GAM,
      ctxPanel.GAM_A,
      blCtx.ALFA,
      blCtx.MINF,
      ctxPanel.QINF ?? 1.0,
    );
    const drag = cdcalc(ctxPanel, blCtx, blCtx.ALFA, ctxPanel.QINF ?? 1.0);
    blCtx.CL = coeffs.cl;
    blCtx.CM = coeffs.cm;
    blCtx.CDP = coeffs.cdp;
    blCtx.CD = drag.cd;
    blCtx.CDF = drag.cdf;
  }

  if (!ctxPanel.LWDIJ || !ctxPanel.LADIJ) {
    qdcalc(ctxPanel);
  }
  blCtx.DIJ = ctxPanel.DIJ;

  const turbTop = 100.0 * Math.exp(-(blCtx.ACRIT[1] + 8.43) / 2.4);
  const turbBot = 100.0 * Math.exp(-(blCtx.ACRIT[2] + 8.43) / 2.4);
  const waklen = ctxPanel.WAKLEN ?? 1.0;
  console.log('');
  console.log('.OPERv   c>  vpar');
  console.log('');
  console.log(` Xtr/c     =${blCtx.XSTRIP[1].toFixed(4).padStart(8)}    top    side`);
  console.log(` Xtr/c     =${blCtx.XSTRIP[2].toFixed(4).padStart(8)}    bottom side`);
  console.log(` Ncrit     =${blCtx.ACRIT[1].toFixed(2).padStart(8)}   (${turbTop.toFixed(3).padStart(6)} % turb. level )`);
  console.log(` Ncrit     =${blCtx.ACRIT[2].toFixed(2).padStart(8)}   (${turbBot.toFixed(3).padStart(6)} % turb. level )`);
  console.log(` Vacc      =${blCtx.VACCEL.toFixed(4).padStart(8)}`);
  console.log(` WakeL/c   =${waklen.toFixed(3).padStart(8)}`);
  console.log('');
  console.log(` Klag  =${blCtx.SCCON.toFixed(4).padStart(8)}     Uxwt  =${blCtx.DUXCON.toFixed(2).padStart(8)}       Kdl =${blCtx.DLCON.toFixed(4).padStart(8)}`);
  console.log(` A     =${blCtx.GACON.toFixed(4).padStart(8)}     B     =${blCtx.GBCON.toFixed(4).padStart(8)}       KCt =${blCtx.CTCON.toFixed(5).padStart(8)}`);
  console.log(` CtiniK=${blCtx.CTRCON.toFixed(4).padStart(8)}     CtiniX=${blCtx.CTRCEX.toFixed(4).padStart(8)}`);

  const { qvis, converged } = solveViscous(blCtx, ctxPanel, qinv, qinvA, maxIter);

  for (let i = 0; i < ctxPanel.N; i += 1) {
    if (!Number.isFinite(ctxPanel.GAM[i])) {
      ctxPanel.GAM[i] = qvis[i + 1] ?? qinv[i + 1] ?? 0.0;
    }
  }

  if (converged) {
    blCtx.LVCONV = true;
  }

  const stats = null;
  return {
    qinv,
    qinvA,
    qvis,
    converged,
    stats,
  };
}

// Reconstruct viscous edge velocities from BL state when QVIS is absent.
function computeQvisFromUedg(blCtx, nb, qinv) {
  const qvis = new Float64Array(nb + 1);
  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= blCtx.NBL[is]; ibl += 1) {
      const i = blCtx.IPAN[ibl][is];
      const qv = blCtx.VTI[ibl][is] * blCtx.UEDG[ibl][is];
      if (i >= 1 && i <= nb) qvis[i] = qv;
    }
  }
  for (let i = 1; i <= nb; i += 1) {
    if (!Number.isFinite(qvis[i]) || qvis[i] === 0.0) {
      qvis[i] = qinv[i] ?? 0.0;
    }
  }
  return qvis;
}

export {
  applyXfoilOperInit,
  buildBlContext,
  bldump,
  cpdump,
  initViscousPanel,
  initViscousBl,
  resetBlState,
  solveViscous,
  specal,
  viscal,
  computeQvisFromUedg,
};
