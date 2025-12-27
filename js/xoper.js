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
import { clcalc, cdcalc } from './xfoil.js';
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
import { createMatrix1, createTensor3 } from './arrays.js';

// Initialize BL operating point (alpha, Mach, Re) to match XFOIL's OPER init.
// Mirrors OPER/SPECAL setup in xoper.f.
function applyXfoilOperInit(blCtx, alphaRad, reinf1) {
  blCtx.LVISC = true;
  blCtx.LVCONV = false;
  blCtx.LWAKE = blCtx.NW > 0;
  blCtx.LALFA = true;
  blCtx.ALFA = alphaRad;
  blCtx.ADEG = alphaRad / blCtx.DTOR;
  blCtx.AWAKE = alphaRad;
  blCtx.AVISC = alphaRad;
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

  const needsTensor = !blCtx.VA
    || blCtx.VA.length <= 3
    || blCtx.VA[1]?.length <= 2
    || blCtx.VA[1]?.[1]?.length <= nsys;
  const needsVm = !blCtx.VM
    || blCtx.VM.length <= 3
    || blCtx.VM[1]?.length <= nsys
    || blCtx.VM[1]?.[1]?.length <= nsys;

  if (needsTensor) {
    blCtx.VA = createTensor3(3, 2, nsys);
    blCtx.VB = createTensor3(3, 2, nsys);
    blCtx.VDEL = createTensor3(3, 2, nsys);
  }
  if (needsVm) {
    blCtx.VM = createTensor3(3, nsys, nsys);
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

  uicalc(blCtx, qinv, qinvA);

  if (!blCtx.LBLINI) {
    for (let is = 1; is <= 2; is += 1) {
      for (let ibl = 1; ibl <= blCtx.NBL[is]; ibl += 1) {
        blCtx.UEDG[ibl][is] = blCtx.UINV[ibl][is];
      }
    }
  }

  if (!ctxPanel.LWDIJ || !ctxPanel.LADIJ) {
    qdcalc(ctxPanel);
  }
  blCtx.DIJ = ctxPanel.DIJ;
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
  console.log(' iter        rms        max    rlx        a        CL        Cm        CD       CDf       CDp     tr1     tr2');
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
    const iterLine = `${String(iter + 1).padStart(5)}`
      + `${Number.isFinite(blCtx.RMSBL) ? blCtx.RMSBL.toExponential(3).padStart(12) : '     NaN'.padStart(12)}`
      + `${Number.isFinite(blCtx.RMXBL) ? blCtx.RMXBL.toExponential(3).padStart(12) : '     NaN'.padStart(12)}`
      + `${rlx.toFixed(2).padStart(7)}`
      + `${(blCtx.ALFA / blCtx.DTOR).toFixed(3).padStart(9)}`
      + `${coeffs.cl.toFixed(5).padStart(10)}`
      + `${coeffs.cm.toFixed(5).padStart(10)}`
      + `${drag.cd.toFixed(6).padStart(10)}`
      + `${drag.cdf.toFixed(6).padStart(10)}`
      + `${coeffs.cdp.toFixed(6).padStart(10)}`
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
  } = opts;

  const { qinv, qinvA } = initViscousPanel(ctxPanel);
  applyXfoilOperInit(blCtx, alphaRad, Number.isFinite(reinf) ? reinf : 1.0e6);
  initViscousBl(blCtx, ctxPanel, qinv, qinvA);

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

  blCtx.LVCONV = converged;

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
  initViscousPanel,
  initViscousBl,
  resetBlState,
  solveViscous,
  specal,
  viscal,
  computeQvisFromUedg,
};
