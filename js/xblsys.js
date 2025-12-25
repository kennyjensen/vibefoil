// Port of src/xblsys.f (boundary-layer system routines).
// Implements XFOIL's integral BL equations, transition logic, and closures.
const NCOM = 73;

const COM1_NAMES = [
  null,
  'X1', 'U1', 'T1', 'D1', 'S1', 'AMPL1', 'U1_UEI', 'U1_MS', 'DW1',
  'H1', 'H1_T1', 'H1_D1',
  'M1', 'M1_U1', 'M1_MS',
  'R1', 'R1_U1', 'R1_MS',
  'V1', 'V1_U1', 'V1_MS', 'V1_RE',
  'HK1', 'HK1_U1', 'HK1_T1', 'HK1_D1', 'HK1_MS',
  'HS1', 'HS1_U1', 'HS1_T1', 'HS1_D1', 'HS1_MS', 'HS1_RE',
  'HC1', 'HC1_U1', 'HC1_T1', 'HC1_D1', 'HC1_MS',
  'RT1', 'RT1_U1', 'RT1_T1', 'RT1_MS', 'RT1_RE',
  'CF1', 'CF1_U1', 'CF1_T1', 'CF1_D1', 'CF1_MS', 'CF1_RE',
  'DI1', 'DI1_U1', 'DI1_T1', 'DI1_D1', 'DI1_S1', 'DI1_MS', 'DI1_RE',
  'US1', 'US1_U1', 'US1_T1', 'US1_D1', 'US1_MS', 'US1_RE',
  'CQ1', 'CQ1_U1', 'CQ1_T1', 'CQ1_D1', 'CQ1_MS', 'CQ1_RE',
  'DE1', 'DE1_U1', 'DE1_T1', 'DE1_D1', 'DE1_MS',
];

const COM2_NAMES = [
  null,
  'X2', 'U2', 'T2', 'D2', 'S2', 'AMPL2', 'U2_UEI', 'U2_MS', 'DW2',
  'H2', 'H2_T2', 'H2_D2',
  'M2', 'M2_U2', 'M2_MS',
  'R2', 'R2_U2', 'R2_MS',
  'V2', 'V2_U2', 'V2_MS', 'V2_RE',
  'HK2', 'HK2_U2', 'HK2_T2', 'HK2_D2', 'HK2_MS',
  'HS2', 'HS2_U2', 'HS2_T2', 'HS2_D2', 'HS2_MS', 'HS2_RE',
  'HC2', 'HC2_U2', 'HC2_T2', 'HC2_D2', 'HC2_MS',
  'RT2', 'RT2_U2', 'RT2_T2', 'RT2_MS', 'RT2_RE',
  'CF2', 'CF2_U2', 'CF2_T2', 'CF2_D2', 'CF2_MS', 'CF2_RE',
  'DI2', 'DI2_U2', 'DI2_T2', 'DI2_D2', 'DI2_S2', 'DI2_MS', 'DI2_RE',
  'US2', 'US2_U2', 'US2_T2', 'US2_D2', 'US2_MS', 'US2_RE',
  'CQ2', 'CQ2_U2', 'CQ2_T2', 'CQ2_D2', 'CQ2_MS', 'CQ2_RE',
  'DE2', 'DE2_U2', 'DE2_T2', 'DE2_D2', 'DE2_MS',
];

// Allocate Fortran-style 1-based 2D arrays for BL state.
function create2d(rows, cols) {
  const arr = new Array(rows + 1);
  for (let i = 0; i <= rows; i += 1) {
    arr[i] = new Array(cols + 1).fill(0.0);
  }
  return arr;
}

// Optional consistency check for BL arrays; mirrors debug outputs in XFOIL.
function blCheck(ctx, label, extra) {
  if (!ctx.DEBUG_BL) return;
  const fields = {
    x1: ctx.X1,
    x2: ctx.X2,
    u1: ctx.U1,
    u2: ctx.U2,
    t1: ctx.T1,
    t2: ctx.T2,
    d1: ctx.D1,
    d2: ctx.D2,
    hk1: ctx.HK1,
    hk2: ctx.HK2,
    rt1: ctx.RT1,
    rt2: ctx.RT2,
    ampl1: ctx.AMPL1,
    ampl2: ctx.AMPL2,
  };
  const bad = Object.entries(fields).filter(([, v]) => !Number.isFinite(v));
  if (bad.length === 0) return;
  const record = { label, bad: bad.map(([k]) => k), ...fields, ...extra };
  if (!ctx.BL_TRACE) ctx.BL_TRACE = [];
  if (ctx.BL_TRACE.length < 50) ctx.BL_TRACE.push(record);
  console.warn('BLCHK:', record);
  if (ctx.DEBUG_BL_FAILFAST) {
    throw new Error(`BLCHK ${label}: non-finite ${bad.map(([k]) => k).join(', ')}`);
  }
}

// Ensure BL context has all arrays; used to reuse or initialize state.
function ensureCtx(ctx) {
  if (ctx.NCOM == null) ctx.NCOM = NCOM;
  if (!ctx.COM1) ctx.COM1 = new Array(NCOM + 1).fill(0.0);
  if (!ctx.COM2) ctx.COM2 = new Array(NCOM + 1).fill(0.0);
  if (!ctx.C1SAV) ctx.C1SAV = new Array(NCOM + 1).fill(0.0);
  if (!ctx.C2SAV) ctx.C2SAV = new Array(NCOM + 1).fill(0.0);
  if (!ctx.VS1) ctx.VS1 = create2d(4, 5);
  if (!ctx.VS2) ctx.VS2 = create2d(4, 5);
  if (!ctx.VSREZ) ctx.VSREZ = new Array(5).fill(0.0);
  if (!ctx.VSR) ctx.VSR = new Array(5).fill(0.0);
  if (!ctx.VSM) ctx.VSM = new Array(5).fill(0.0);
  if (!ctx.VSX) ctx.VSX = new Array(5).fill(0.0);
}

// Copy COM common-block data into local variables for computations.
function syncComToVars(ctx, which) {
  ensureCtx(ctx);
  const names = which === 1 ? COM1_NAMES : COM2_NAMES;
  const src = which === 1 ? ctx.COM1 : ctx.COM2;
  for (let i = 1; i <= NCOM; i += 1) {
    ctx[names[i]] = src[i];
  }
}

// Copy local variables back into COM common-block storage.
function syncVarsToCom(ctx, which) {
  ensureCtx(ctx);
  const names = which === 1 ? COM1_NAMES : COM2_NAMES;
  const dst = which === 1 ? ctx.COM1 : ctx.COM2;
  for (let i = 1; i <= NCOM; i += 1) {
    dst[i] = ctx[names[i]] ?? 0.0;
  }
}

// Snapshot/restore of COM state for solver staging and retries.
function copyCom(ctx, from, to) {
  ensureCtx(ctx);
  const src = from === 1 ? ctx.COM1 : ctx.COM2;
  const dst = to === 1 ? ctx.COM1 : ctx.COM2;
  for (let i = 1; i <= NCOM; i += 1) {
    dst[i] = src[i];
  }
  syncComToVars(ctx, to);
}

// Transition check wrapper (delegates to the main routine).
function trchek(ctxIn) {
  const ctx = ctxIn || this;
  trchek2(ctx);
}

// Auxiliary quantities for transition/separation damping logic.
function axset(hk1, t1, rt1, a1, hk2, t2, rt2, a2, acrit, idampv) {
  let ax1;
  let ax1Hk1;
  let ax1T1;
  let ax1Rt1;
  let ax2;
  let ax2Hk2;
  let ax2T2;
  let ax2Rt2;

  if (idampv === 0) {
    ({ ax: ax1, axHk: ax1Hk1, axTh: ax1T1, axRt: ax1Rt1 } = dampl(hk1, t1, rt1));
    ({ ax: ax2, axHk: ax2Hk2, axTh: ax2T2, axRt: ax2Rt2 } = dampl(hk2, t2, rt2));
  } else {
    ({ ax: ax1, axHk: ax1Hk1, axTh: ax1T1, axRt: ax1Rt1 } = dampl2(hk1, t1, rt1));
    ({ ax: ax2, axHk: ax2Hk2, axTh: ax2T2, axRt: ax2Rt2 } = dampl2(hk2, t2, rt2));
  }

  const axsq = 0.5 * (ax1 ** 2 + ax2 ** 2);
  let axa;
  let axaAx1;
  let axaAx2;
  if (axsq <= 0.0) {
    axa = 0.0;
    axaAx1 = 0.0;
    axaAx2 = 0.0;
  } else {
    axa = Math.sqrt(axsq);
    axaAx1 = 0.5 * ax1 / axa;
    axaAx2 = 0.5 * ax2 / axa;
  }

  const arg = Math.min(20.0 * (acrit - 0.5 * (a1 + a2)), 20.0);
  let exn;
  let exnA1;
  let exnA2;
  if (arg <= 0.0) {
    exn = 1.0;
    exnA1 = 0.0;
    exnA2 = 0.0;
  } else {
    exn = Math.exp(-arg);
    exnA1 = 20.0 * 0.5 * exn;
    exnA2 = 20.0 * 0.5 * exn;
  }

  const dax = exn * 0.002 / (t1 + t2);
  const daxA1 = exnA1 * 0.002 / (t1 + t2);
  const daxA2 = exnA2 * 0.002 / (t1 + t2);
  const daxT1 = -dax / (t1 + t2);
  const daxT2 = -dax / (t1 + t2);

  const ax = axa + dax;
  const axHk1 = axaAx1 * ax1Hk1;
  const axT1 = axaAx1 * ax1T1 + daxT1;
  const axRt1 = axaAx1 * ax1Rt1;
  const axA1 = daxA1;
  const axHk2 = axaAx2 * ax2Hk2;
  const axT2 = axaAx2 * ax2T2 + daxT2;
  const axRt2 = axaAx2 * ax2Rt2;
  const axA2 = daxA2;

  return {
    ax,
    axHk1,
    axT1,
    axRt1,
    axA1,
    axHk2,
    axT2,
    axRt2,
    axA2,
  };
}

// Main transition check pass (e^N criterion, separation adjustments).
function trchek2(ctx) {
  ensureCtx(ctx);
  const daeps = 5.0e-5;

  for (let icom = 1; icom <= NCOM; icom += 1) {
    ctx.C2SAV[icom] = ctx.COM2[icom];
  }

  let axOut = axset(ctx.HK1, ctx.T1, ctx.RT1, ctx.AMPL1,
    ctx.HK2, ctx.T2, ctx.RT2, ctx.AMPL2, ctx.AMCRIT, ctx.IDAMPV);
  let ax = axOut.ax;
  ctx.AMPL2 = ctx.AMPL1 + ax * (ctx.X2 - ctx.X1);

  let wf2 = 0.0;
  let wf2A1 = 0.0;
  let wf2A2 = 0.0;
  let wf2X1 = 0.0;
  let wf2X2 = 0.0;
  let wf2Xf = 0.0;
  let wf1 = 0.0;
  let wf1A1 = 0.0;
  let wf1A2 = 0.0;
  let wf1X1 = 0.0;
  let wf1X2 = 0.0;
  let wf1Xf = 0.0;

  let xt = 0.0;
  let tt = 0.0;
  let dt = 0.0;
  let ut = 0.0;
  let xtA2 = 0.0;
  let ttA2 = 0.0;
  let dtA2 = 0.0;
  let utA2 = 0.0;

  let hkt = 0.0;
  let hktTt = 0.0;
  let hktDt = 0.0;
  let hktUt = 0.0;
  let hktMs = 0.0;
  let rtt = 0.0;
  let rttTt = 0.0;
  let rttUt = 0.0;
  let rttMs = 0.0;
  let rttRe = 0.0;
  let amplt = 0.0;
  let ampltA2 = 0.0;

  for (let itam = 1; itam <= 30; itam += 1) {
    let sfa;
    let sfaA1;
    let sfaA2;
    if (ctx.AMPL2 <= ctx.AMCRIT) {
      amplt = ctx.AMPL2;
      ampltA2 = 1.0;
      sfa = 1.0;
      sfaA1 = 0.0;
      sfaA2 = 0.0;
    } else {
      amplt = ctx.AMCRIT;
      ampltA2 = 0.0;
      sfa = (amplt - ctx.AMPL1) / (ctx.AMPL2 - ctx.AMPL1);
      sfaA1 = (sfa - 1.0) / (ctx.AMPL2 - ctx.AMPL1);
      sfaA2 = (-sfa) / (ctx.AMPL2 - ctx.AMPL1);
    }

    let sfx;
    let sfxX1;
    let sfxX2;
    let sfxXf;
    if (ctx.XIFORC < ctx.X2) {
      sfx = (ctx.XIFORC - ctx.X1) / (ctx.X2 - ctx.X1);
      sfxX1 = (sfx - 1.0) / (ctx.X2 - ctx.X1);
      sfxX2 = (-sfx) / (ctx.X2 - ctx.X1);
      sfxXf = 1.0 / (ctx.X2 - ctx.X1);
    } else {
      sfx = 1.0;
      sfxX1 = 0.0;
      sfxX2 = 0.0;
      sfxXf = 0.0;
    }

    if (sfa < sfx) {
      wf2 = sfa;
      wf2A1 = sfaA1;
      wf2A2 = sfaA2;
      wf2X1 = 0.0;
      wf2X2 = 0.0;
      wf2Xf = 0.0;
    } else {
      wf2 = sfx;
      wf2A1 = 0.0;
      wf2A2 = 0.0;
      wf2X1 = sfxX1;
      wf2X2 = sfxX2;
      wf2Xf = sfxXf;
    }

    wf1 = 1.0 - wf2;
    wf1A1 = -wf2A1;
    wf1A2 = -wf2A2;
    wf1X1 = -wf2X1;
    wf1X2 = -wf2X2;
    wf1Xf = -wf2Xf;

    xt = ctx.X1 * wf1 + ctx.X2 * wf2;
    tt = ctx.T1 * wf1 + ctx.T2 * wf2;
    dt = ctx.D1 * wf1 + ctx.D2 * wf2;
    ut = ctx.U1 * wf1 + ctx.U2 * wf2;

    xtA2 = ctx.X1 * wf1A2 + ctx.X2 * wf2A2;
    ttA2 = ctx.T1 * wf1A2 + ctx.T2 * wf2A2;
    dtA2 = ctx.D1 * wf1A2 + ctx.D2 * wf2A2;
    utA2 = ctx.U1 * wf1A2 + ctx.U2 * wf2A2;

    ctx.X2 = xt;
    ctx.T2 = tt;
    ctx.D2 = dt;
    ctx.U2 = ut;
    syncVarsToCom(ctx, 2);

    blkin(ctx);

    hkt = ctx.HK2;
    hktTt = ctx.HK2_T2;
    hktDt = ctx.HK2_D2;
    hktUt = ctx.HK2_U2;
    hktMs = ctx.HK2_MS;

    rtt = ctx.RT2;
    rttTt = ctx.RT2_T2;
    rttUt = ctx.RT2_U2;
    rttMs = ctx.RT2_MS;
    rttRe = ctx.RT2_RE;

    const amsave = ctx.AMPL2;
    for (let icom = 1; icom <= NCOM; icom += 1) {
      ctx.COM2[icom] = ctx.C2SAV[icom];
    }
    syncComToVars(ctx, 2);
    ctx.AMPL2 = amsave;

    axOut = axset(ctx.HK1, ctx.T1, ctx.RT1, ctx.AMPL1,
      hkt, tt, rtt, amplt, ctx.AMCRIT, ctx.IDAMPV);
    ax = axOut.ax;

    if (ax <= 0.0) break;

    const axA2 = (axOut.axHk2 * hktTt + axOut.axT2 + axOut.axRt2 * rttTt) * ttA2
      + (axOut.axHk2 * hktDt) * dtA2
      + (axOut.axHk2 * hktUt + axOut.axRt2 * rttUt) * utA2
      + axOut.axA2 * ampltA2;

    const res = ctx.AMPL2 - ctx.AMPL1 - ax * (ctx.X2 - ctx.X1);
    const resA2 = 1.0 - axA2 * (ctx.X2 - ctx.X1);
    let da2 = -res / resA2;

    let rlx = 1.0;
    const dxt = xtA2 * da2;
    if (rlx * Math.abs(dxt / (ctx.X2 - ctx.X1)) > 0.05) {
      rlx = 0.05 * Math.abs((ctx.X2 - ctx.X1) / dxt);
    }
    if (rlx * Math.abs(da2) > 1.0) {
      rlx = 1.0 * Math.abs(1.0 / da2);
    }

    if (Math.abs(da2) < daeps) break;

    if ((ctx.AMPL2 > ctx.AMCRIT && ctx.AMPL2 + rlx * da2 < ctx.AMCRIT)
      || (ctx.AMPL2 < ctx.AMCRIT && ctx.AMPL2 + rlx * da2 > ctx.AMCRIT)) {
      ctx.AMPL2 = ctx.AMCRIT;
    } else {
      ctx.AMPL2 = ctx.AMPL2 + rlx * da2;
    }
  }

  ctx.XT = xt;
  ctx.XT_A2 = xtA2;

  ctx.TRFREE = ctx.AMPL2 >= ctx.AMCRIT;
  ctx.TRFORC = ctx.XIFORC > ctx.X1 && ctx.XIFORC <= ctx.X2;
  ctx.TRAN = ctx.TRFORC || ctx.TRFREE;
  if (!ctx.TRAN) return;

  if (ctx.TRFREE && ctx.TRFORC) {
    ctx.TRFORC = ctx.XIFORC < xt;
    ctx.TRFREE = ctx.XIFORC >= xt;
  }

  if (ctx.TRFORC) {
    ctx.XT = ctx.XIFORC;
    ctx.XT_A1 = 0.0;
    ctx.XT_X1 = 0.0;
    ctx.XT_T1 = 0.0;
    ctx.XT_D1 = 0.0;
    ctx.XT_U1 = 0.0;
    ctx.XT_X2 = 0.0;
    ctx.XT_T2 = 0.0;
    ctx.XT_D2 = 0.0;
    ctx.XT_U2 = 0.0;
    ctx.XT_MS = 0.0;
    ctx.XT_RE = 0.0;
    ctx.XT_XF = 1.0;
    return;
  }

  ctx.XT = ctx.X1 * wf1 + ctx.X2 * wf2;
  ctx.XT_A2 = xtA2;
  tt = ctx.T1 * wf1 + ctx.T2 * wf2;
  dt = ctx.D1 * wf1 + ctx.D2 * wf2;
  ut = ctx.U1 * wf1 + ctx.U2 * wf2;

  ctx.XT_X1 = wf1;
  let ttT1 = wf1;
  let dtD1 = wf1;
  let utU1 = wf1;

  ctx.XT_X2 = wf2;
  let ttT2 = wf2;
  let dtD2 = wf2;
  let utU2 = wf2;

  ctx.XT_A1 = ctx.X1 * wf1A1 + ctx.X2 * wf2A1;
  let ttA1 = ctx.T1 * wf1A1 + ctx.T2 * wf2A1;
  let dtA1 = ctx.D1 * wf1A1 + ctx.D2 * wf2A1;
  let utA1 = ctx.U1 * wf1A1 + ctx.U2 * wf2A1;

  ctx.XT_X1 = ctx.X1 * wf1X1 + ctx.X2 * wf2X1 + ctx.XT_X1;
  let ttX1 = ctx.T1 * wf1X1 + ctx.T2 * wf2X1;
  let dtX1 = ctx.D1 * wf1X1 + ctx.D2 * wf2X1;
  let utX1 = ctx.U1 * wf1X1 + ctx.U2 * wf2X1;

  ctx.XT_X2 = ctx.X1 * wf1X2 + ctx.X2 * wf2X2 + ctx.XT_X2;
  let ttX2 = ctx.T1 * wf1X2 + ctx.T2 * wf2X2;
  let dtX2 = ctx.D1 * wf1X2 + ctx.D2 * wf2X2;
  let utX2 = ctx.U1 * wf1X2 + ctx.U2 * wf2X2;

  ctx.XT_XF = ctx.X1 * wf1Xf + ctx.X2 * wf2Xf;
  let ttXf = ctx.T1 * wf1Xf + ctx.T2 * wf2Xf;
  let dtXf = ctx.D1 * wf1Xf + ctx.D2 * wf2Xf;
  let utXf = ctx.U1 * wf1Xf + ctx.U2 * wf2Xf;

  let axT1 = axOut.axHk1 * ctx.HK1_T1 + axOut.axT1 + axOut.axRt1 * ctx.RT1_T1
    + (axOut.axHk2 * hktTt + axOut.axT2 + axOut.axRt2 * rttTt) * ttT1;
  let axD1 = axOut.axHk1 * ctx.HK1_D1 + (axOut.axHk2 * hktDt) * dtD1;
  let axU1 = axOut.axHk1 * ctx.HK1_U1 + axOut.axRt1 * ctx.RT1_U1
    + (axOut.axHk2 * hktUt + axOut.axRt2 * rttUt) * utU1;
  let axA1 = axOut.axA1
    + (axOut.axHk2 * hktTt + axOut.axT2 + axOut.axRt2 * rttTt) * ttA1
    + (axOut.axHk2 * hktDt) * dtA1
    + (axOut.axHk2 * hktUt + axOut.axRt2 * rttUt) * utA1;
  let axX1 = (axOut.axHk2 * hktTt + axOut.axT2 + axOut.axRt2 * rttTt) * ttX1
    + (axOut.axHk2 * hktDt) * dtX1
    + (axOut.axHk2 * hktUt + axOut.axRt2 * rttUt) * utX1;

  let axT2 = (axOut.axHk2 * hktTt + axOut.axT2 + axOut.axRt2 * rttTt) * ttT2;
  let axD2 = (axOut.axHk2 * hktDt) * dtD2;
  let axU2 = (axOut.axHk2 * hktUt + axOut.axRt2 * rttUt) * utU2;
  let axA2 = axOut.axA2 * ampltA2
    + (axOut.axHk2 * hktTt + axOut.axT2 + axOut.axRt2 * rttTt) * ttA2
    + (axOut.axHk2 * hktDt) * dtA2
    + (axOut.axHk2 * hktUt + axOut.axRt2 * rttUt) * utA2;
  let axX2 = (axOut.axHk2 * hktTt + axOut.axT2 + axOut.axRt2 * rttTt) * ttX2
    + (axOut.axHk2 * hktDt) * dtX2
    + (axOut.axHk2 * hktUt + axOut.axRt2 * rttUt) * utX2;

  let axXf = (axOut.axHk2 * hktTt + axOut.axT2 + axOut.axRt2 * rttTt) * ttXf
    + (axOut.axHk2 * hktDt) * dtXf
    + (axOut.axHk2 * hktUt + axOut.axRt2 * rttUt) * utXf;

  let axMs = axOut.axHk2 * hktMs + axOut.axRt2 * rttMs + axOut.axHk1 * ctx.HK1_MS + axOut.axRt1 * ctx.RT1_MS;
  let axRe = axOut.axRt2 * rttRe + axOut.axRt1 * ctx.RT1_RE;

  const zAx = -(ctx.X2 - ctx.X1);
  const zA1 = zAx * axA1 - 1.0;
  const zT1 = zAx * axT1;
  const zD1 = zAx * axD1;
  const zU1 = zAx * axU1;
  const zX1 = zAx * axX1 + ax;
  const zA2 = zAx * axA2 + 1.0;
  const zT2 = zAx * axT2;
  const zD2 = zAx * axD2;
  const zU2 = zAx * axU2;
  const zX2 = zAx * axX2 - ax;
  const zXf = zAx * axXf;
  const zMs = zAx * axMs;
  const zRe = zAx * axRe;

  ctx.XT_A1 = ctx.XT_A1 - (ctx.XT_A2 / zA2) * zA1;
  ctx.XT_T1 = -(ctx.XT_A2 / zA2) * zT1;
  ctx.XT_D1 = -(ctx.XT_A2 / zA2) * zD1;
  ctx.XT_U1 = -(ctx.XT_A2 / zA2) * zU1;
  ctx.XT_X1 = ctx.XT_X1 - (ctx.XT_A2 / zA2) * zX1;
  ctx.XT_T2 = -(ctx.XT_A2 / zA2) * zT2;
  ctx.XT_D2 = -(ctx.XT_A2 / zA2) * zD2;
  ctx.XT_U2 = -(ctx.XT_A2 / zA2) * zU2;
  ctx.XT_X2 = ctx.XT_X2 - (ctx.XT_A2 / zA2) * zX2;
  ctx.XT_MS = -(ctx.XT_A2 / zA2) * zMs;
  ctx.XT_RE = -(ctx.XT_A2 / zA2) * zRe;
  ctx.XT_XF = 0.0;
}

// Assemble BL system matrices and residuals (two-equation integral method).
function blsys(ctxIn) {
  const ctx = ctxIn || this;
  ensureCtx(ctx);

  if (ctx.WAKE) {
    blvar(3, ctx);
    blmid(3, ctx);
  } else if (ctx.TURB || ctx.TRAN) {
    blvar(2, ctx);
    blmid(2, ctx);
  } else {
    blvar(1, ctx);
    blmid(1, ctx);
  }

  syncVarsToCom(ctx, 2);

  if (ctx.SIMI) {
    for (let icom = 1; icom <= NCOM; icom += 1) {
      ctx.COM1[icom] = ctx.COM2[icom];
    }
    syncComToVars(ctx, 1);
  }

  if (ctx.TRAN) {
    trdif(ctx);
  } else if (ctx.SIMI) {
    bldif(0, ctx);
  } else if (!ctx.TURB) {
    bldif(1, ctx);
  } else if (ctx.WAKE) {
    bldif(3, ctx);
  } else {
    bldif(2, ctx);
  }

  blCheck(ctx, 'BLSYS:after-bldif', {
    simi: ctx.SIMI,
    tran: ctx.TRAN,
    turb: ctx.TURB,
    wake: ctx.WAKE,
  });

  if (!Number.isFinite(ctx.VS1[1][1]) || !Number.isFinite(ctx.VS2[1][1])) {
    console.warn('BLSYS: non-finite VS(1,1)', {
      simi: ctx.SIMI,
      tran: ctx.TRAN,
      turb: ctx.TURB,
      wake: ctx.WAKE,
      x1: ctx.X1,
      x2: ctx.X2,
      t1: ctx.T1,
      t2: ctx.T2,
      d1: ctx.D1,
      d2: ctx.D2,
      u1: ctx.U1,
      u2: ctx.U2,
      hk1: ctx.HK1,
      hk2: ctx.HK2,
      rt1: ctx.RT1,
      rt2: ctx.RT2,
      ampl1: ctx.AMPL1,
      ampl2: ctx.AMPL2,
      s1: ctx.S1,
      s2: ctx.S2,
      vs11: ctx.VS1[1][1],
      vs21: ctx.VS2[1][1],
    });
  }

  if (ctx.SIMI) {
    for (let k = 1; k <= 4; k += 1) {
      for (let l = 1; l <= 5; l += 1) {
        ctx.VS2[k][l] = ctx.VS1[k][l] + ctx.VS2[k][l];
        ctx.VS1[k][l] = 0.0;
      }
    }
  }

  for (let k = 1; k <= 4; k += 1) {
    const resU1 = ctx.VS1[k][4];
    const resU2 = ctx.VS2[k][4];
    const resMs = ctx.VSM[k];

    ctx.VS1[k][4] = resU1 * ctx.U1_UEI;
    ctx.VS2[k][4] = resU2 * ctx.U2_UEI;
    ctx.VSM[k] = resU1 * ctx.U1_MS + resU2 * ctx.U2_MS + resMs;
  }
}

// Trailing-edge system coupling for displacement thickness in the wake.
function tesys(cte, tte, dte, ctxIn) {
  const ctx = ctxIn || this;
  ensureCtx(ctx);

  for (let k = 1; k <= 4; k += 1) {
    ctx.VSREZ[k] = 0.0;
    ctx.VSM[k] = 0.0;
    ctx.VSR[k] = 0.0;
    ctx.VSX[k] = 0.0;
    for (let l = 1; l <= 5; l += 1) {
      ctx.VS1[k][l] = 0.0;
      ctx.VS2[k][l] = 0.0;
    }
  }

  blvar(3, ctx);

  ctx.VS1[1][1] = -1.0;
  ctx.VS2[1][1] = 1.0;
  ctx.VSREZ[1] = cte - ctx.S2;

  ctx.VS1[2][2] = -1.0;
  ctx.VS2[2][2] = 1.0;
  ctx.VSREZ[2] = tte - ctx.T2;

  ctx.VS1[3][3] = -1.0;
  ctx.VS2[3][3] = 1.0;
  ctx.VSREZ[3] = dte - ctx.D2 - ctx.DW2;
}

// Predictor step for BL marching (advance with previous station values).
function blprv(xsi, ami, cti, thi, dsi, dswaki, uei, ctxIn) {
  const ctx = ctxIn || this;
  ensureCtx(ctx);

  ctx.X2 = xsi;
  ctx.AMPL2 = ami;
  ctx.S2 = cti;
  ctx.T2 = thi;
  ctx.D2 = dsi - dswaki;
  ctx.DW2 = dswaki;

  ctx.U2 = uei * (1.0 - ctx.TKBL) / (1.0 - ctx.TKBL * (uei / ctx.QINFBL) ** 2);
  ctx.U2_UEI = (1.0 + ctx.TKBL * (2.0 * ctx.U2 * uei / ctx.QINFBL ** 2 - 1.0))
    / (1.0 - ctx.TKBL * (uei / ctx.QINFBL) ** 2);
  ctx.U2_MS = (ctx.U2 * (uei / ctx.QINFBL) ** 2 - uei) * ctx.TKBL_MS
    / (1.0 - ctx.TKBL * (uei / ctx.QINFBL) ** 2);

  syncVarsToCom(ctx, 2);
}

// Initialize BL profiles at the leading edge (laminar start).
function blkin(ctxIn) {
  const ctx = ctxIn || this;
  ensureCtx(ctx);

  ctx.M2 = ctx.U2 * ctx.U2 * ctx.HSTINV / (ctx.GM1BL * (1.0 - 0.5 * ctx.U2 * ctx.U2 * ctx.HSTINV));
  const tr2 = 1.0 + 0.5 * ctx.GM1BL * ctx.M2;
  ctx.M2_U2 = 2.0 * ctx.M2 * tr2 / ctx.U2;
  ctx.M2_MS = ctx.U2 * ctx.U2 * tr2 / (ctx.GM1BL * (1.0 - 0.5 * ctx.U2 * ctx.U2 * ctx.HSTINV))
    * ctx.HSTINV_MS;

  ctx.R2 = ctx.RSTBL * tr2 ** (-1.0 / ctx.GM1BL);
  ctx.R2_U2 = -ctx.R2 / tr2 * 0.5 * ctx.M2_U2;
  ctx.R2_MS = -ctx.R2 / tr2 * 0.5 * ctx.M2_MS + ctx.RSTBL_MS * tr2 ** (-1.0 / ctx.GM1BL);

  ctx.H2 = ctx.D2 / ctx.T2;
  ctx.H2_D2 = 1.0 / ctx.T2;
  ctx.H2_T2 = -ctx.H2 / ctx.T2;

  const herat = 1.0 - 0.5 * ctx.U2 * ctx.U2 * ctx.HSTINV;
  const heU2 = -ctx.U2 * ctx.HSTINV;
  const heMs = -0.5 * ctx.U2 * ctx.U2 * ctx.HSTINV_MS;

  ctx.V2 = Math.sqrt(herat ** 3) * (1.0 + ctx.HVRAT) / (herat + ctx.HVRAT) / ctx.REYBL;
  const v2He = ctx.V2 * (1.5 / herat - 1.0 / (herat + ctx.HVRAT));

  ctx.V2_U2 = v2He * heU2;
  ctx.V2_MS = -ctx.V2 / ctx.REYBL * ctx.REYBL_MS + v2He * heMs;
  ctx.V2_RE = -ctx.V2 / ctx.REYBL * ctx.REYBL_RE;

  const hkinOut = hkin(ctx.H2, ctx.M2);
  ctx.HK2 = hkinOut.hk;
  const hk2H2 = hkinOut.hkH;
  const hk2M2 = hkinOut.hkMsq;

  ctx.HK2_U2 = hk2M2 * ctx.M2_U2;
  ctx.HK2_T2 = hk2H2 * ctx.H2_T2;
  ctx.HK2_D2 = hk2H2 * ctx.H2_D2;
  ctx.HK2_MS = hk2M2 * ctx.M2_MS;

  ctx.RT2 = ctx.R2 * ctx.U2 * ctx.T2 / ctx.V2;
  ctx.RT2_U2 = ctx.RT2 * (1.0 / ctx.U2 + ctx.R2_U2 / ctx.R2 - ctx.V2_U2 / ctx.V2);
  ctx.RT2_T2 = ctx.RT2 / ctx.T2;
  ctx.RT2_MS = ctx.RT2 * (ctx.R2_MS / ctx.R2 - ctx.V2_MS / ctx.V2);
  ctx.RT2_RE = ctx.RT2 * (-ctx.V2_RE / ctx.V2);

  syncVarsToCom(ctx, 2);
}

// Variable updates from BL system solution (laminar/turbulent branches).
function blvar(ityp, ctxIn) {
  const ctx = ctxIn || this;
  ensureCtx(ctx);

  if (ityp === 3) ctx.HK2 = Math.max(ctx.HK2, 1.00005);
  if (ityp !== 3) ctx.HK2 = Math.max(ctx.HK2, 1.05000);

  const hcOut = hct(ctx.HK2, ctx.M2);
  ctx.HC2 = hcOut.hc;
  const hc2Hk2 = hcOut.hcHk;
  const hc2M2 = hcOut.hcMsq;
  ctx.HC2_U2 = hc2Hk2 * ctx.HK2_U2 + hc2M2 * ctx.M2_U2;
  ctx.HC2_T2 = hc2Hk2 * ctx.HK2_T2;
  ctx.HC2_D2 = hc2Hk2 * ctx.HK2_D2;
  ctx.HC2_MS = hc2Hk2 * ctx.HK2_MS + hc2M2 * ctx.M2_MS;

  let hsOut;
  if (ityp === 1) {
    hsOut = hsl(ctx.HK2, ctx.RT2, ctx.M2);
  } else {
    hsOut = hst(ctx.HK2, ctx.RT2, ctx.M2);
  }
  ctx.HS2 = hsOut.hs;
  const hs2Hk2 = hsOut.hsHk;
  const hs2Rt2 = hsOut.hsRt;
  const hs2M2 = hsOut.hsMsq;

  ctx.HS2_U2 = hs2Hk2 * ctx.HK2_U2 + hs2Rt2 * ctx.RT2_U2 + hs2M2 * ctx.M2_U2;
  ctx.HS2_T2 = hs2Hk2 * ctx.HK2_T2 + hs2Rt2 * ctx.RT2_T2;
  ctx.HS2_D2 = hs2Hk2 * ctx.HK2_D2;
  ctx.HS2_MS = hs2Hk2 * ctx.HK2_MS + hs2Rt2 * ctx.RT2_MS + hs2M2 * ctx.M2_MS;
  ctx.HS2_RE = hs2Rt2 * ctx.RT2_RE;

  ctx.US2 = 0.5 * ctx.HS2 * (1.0 - (ctx.HK2 - 1.0) / (ctx.GBCON * ctx.H2));
  const us2Hs2 = 0.5 * (1.0 - (ctx.HK2 - 1.0) / (ctx.GBCON * ctx.H2));
  const us2Hk2 = 0.5 * ctx.HS2 * (-1.0 / (ctx.GBCON * ctx.H2));
  const us2H2 = 0.5 * ctx.HS2 * (ctx.HK2 - 1.0) / (ctx.GBCON * ctx.H2 ** 2);

  ctx.US2_U2 = us2Hs2 * ctx.HS2_U2 + us2Hk2 * ctx.HK2_U2;
  ctx.US2_T2 = us2Hs2 * ctx.HS2_T2 + us2Hk2 * ctx.HK2_T2 + us2H2 * ctx.H2_T2;
  ctx.US2_D2 = us2Hs2 * ctx.HS2_D2 + us2Hk2 * ctx.HK2_D2 + us2H2 * ctx.H2_D2;
  ctx.US2_MS = us2Hs2 * ctx.HS2_MS + us2Hk2 * ctx.HK2_MS;
  ctx.US2_RE = us2Hs2 * ctx.HS2_RE;

  if (ityp <= 2 && ctx.US2 > 0.95) {
    ctx.US2 = 0.98;
    ctx.US2_U2 = 0.0;
    ctx.US2_T2 = 0.0;
    ctx.US2_D2 = 0.0;
    ctx.US2_MS = 0.0;
    ctx.US2_RE = 0.0;
  }
  if (ityp === 3 && ctx.US2 > 0.99995) {
    ctx.US2 = 0.99995;
    ctx.US2_U2 = 0.0;
    ctx.US2_T2 = 0.0;
    ctx.US2_D2 = 0.0;
    ctx.US2_MS = 0.0;
    ctx.US2_RE = 0.0;
  }

  let gcc = 0.0;
  let hkc = ctx.HK2 - 1.0;
  let hkcHk2 = 1.0;
  let hkcRt2 = 0.0;
  if (ityp === 2) {
    gcc = ctx.GCCON;
    hkc = ctx.HK2 - 1.0 - gcc / ctx.RT2;
    hkcHk2 = 1.0;
    hkcRt2 = gcc / ctx.RT2 ** 2;
    if (hkc < 0.01) {
      hkc = 0.01;
      hkcHk2 = 0.0;
      hkcRt2 = 0.0;
    }
  }

  const hkb = ctx.HK2 - 1.0;
  const usb = 1.0 - ctx.US2;
  ctx.CQ2 = Math.sqrt(ctx.CTCON * ctx.HS2 * hkb * hkc ** 2 / (usb * ctx.H2 * ctx.HK2 ** 2));
  const cq2Hs2 = ctx.CTCON * hkb * hkc ** 2 / (usb * ctx.H2 * ctx.HK2 ** 2) * 0.5 / ctx.CQ2;
  const cq2Us2 = ctx.CTCON * ctx.HS2 * hkb * hkc ** 2 / (usb * ctx.H2 * ctx.HK2 ** 2) / usb * 0.5 / ctx.CQ2;
  const cq2Hk2 = ctx.CTCON * ctx.HS2 * hkc ** 2 / (usb * ctx.H2 * ctx.HK2 ** 2) * 0.5 / ctx.CQ2
    - ctx.CTCON * ctx.HS2 * hkb * hkc ** 2 / (usb * ctx.H2 * ctx.HK2 ** 3) * 2.0 * 0.5 / ctx.CQ2
    + ctx.CTCON * ctx.HS2 * hkb * hkc / (usb * ctx.H2 * ctx.HK2 ** 2) * 2.0 * 0.5 / ctx.CQ2 * hkcHk2;
  const cq2Rt2 = ctx.CTCON * ctx.HS2 * hkb * hkc / (usb * ctx.H2 * ctx.HK2 ** 2) * 2.0 * 0.5 / ctx.CQ2 * hkcRt2;
  const cq2H2 = -ctx.CTCON * ctx.HS2 * hkb * hkc ** 2 / (usb * ctx.H2 * ctx.HK2 ** 2) / ctx.H2 * 0.5 / ctx.CQ2;

  ctx.CQ2_U2 = cq2Hs2 * ctx.HS2_U2 + cq2Us2 * ctx.US2_U2 + cq2Hk2 * ctx.HK2_U2;
  ctx.CQ2_T2 = cq2Hs2 * ctx.HS2_T2 + cq2Us2 * ctx.US2_T2 + cq2Hk2 * ctx.HK2_T2;
  ctx.CQ2_D2 = cq2Hs2 * ctx.HS2_D2 + cq2Us2 * ctx.US2_D2 + cq2Hk2 * ctx.HK2_D2;
  ctx.CQ2_MS = cq2Hs2 * ctx.HS2_MS + cq2Us2 * ctx.US2_MS + cq2Hk2 * ctx.HK2_MS;
  ctx.CQ2_RE = cq2Hs2 * ctx.HS2_RE + cq2Us2 * ctx.US2_RE;

  ctx.CQ2_U2 += cq2Rt2 * ctx.RT2_U2;
  ctx.CQ2_T2 += cq2H2 * ctx.H2_T2 + cq2Rt2 * ctx.RT2_T2;
  ctx.CQ2_D2 += cq2H2 * ctx.H2_D2;
  ctx.CQ2_MS += cq2Rt2 * ctx.RT2_MS;
  ctx.CQ2_RE += cq2Rt2 * ctx.RT2_RE;

  let cfOut;
  if (ityp === 3) {
    ctx.CF2 = 0.0;
    ctx.CF2_HK2 = 0.0;
    ctx.CF2_RT2 = 0.0;
    ctx.CF2_M2 = 0.0;
  } else if (ityp === 1) {
    cfOut = cfl(ctx.HK2, ctx.RT2, ctx.M2);
    ctx.CF2 = cfOut.cf;
    ctx.CF2_HK2 = cfOut.cfHk;
    ctx.CF2_RT2 = cfOut.cfRt;
    ctx.CF2_M2 = cfOut.cfMsq;
  } else {
    cfOut = cft(ctx.HK2, ctx.RT2, ctx.M2, ctx.CFFAC ?? 1.0);
    ctx.CF2 = cfOut.cf;
    ctx.CF2_HK2 = cfOut.cfHk;
    ctx.CF2_RT2 = cfOut.cfRt;
    ctx.CF2_M2 = cfOut.cfMsq;
    const cfL = cfl(ctx.HK2, ctx.RT2, ctx.M2);
    if (cfL.cf > ctx.CF2) {
      ctx.CF2 = cfL.cf;
      ctx.CF2_HK2 = cfL.cfHk;
      ctx.CF2_RT2 = cfL.cfRt;
      ctx.CF2_M2 = cfL.cfMsq;
    }
  }

  ctx.CF2_U2 = ctx.CF2_HK2 * ctx.HK2_U2 + ctx.CF2_RT2 * ctx.RT2_U2 + ctx.CF2_M2 * ctx.M2_U2;
  ctx.CF2_T2 = ctx.CF2_HK2 * ctx.HK2_T2 + ctx.CF2_RT2 * ctx.RT2_T2;
  ctx.CF2_D2 = ctx.CF2_HK2 * ctx.HK2_D2;
  ctx.CF2_MS = ctx.CF2_HK2 * ctx.HK2_MS + ctx.CF2_RT2 * ctx.RT2_MS + ctx.CF2_M2 * ctx.M2_MS;
  ctx.CF2_RE = ctx.CF2_RT2 * ctx.RT2_RE;

  if (ityp === 1) {
    const diOut = dil(ctx.HK2, ctx.RT2);
    ctx.DI2 = diOut.di;
    const di2Hk2 = diOut.diHk;
    const di2Rt2 = diOut.diRt;

    ctx.DI2_U2 = di2Hk2 * ctx.HK2_U2 + di2Rt2 * ctx.RT2_U2;
    ctx.DI2_T2 = di2Hk2 * ctx.HK2_T2 + di2Rt2 * ctx.RT2_T2;
    ctx.DI2_D2 = di2Hk2 * ctx.HK2_D2;
    ctx.DI2_S2 = 0.0;
    ctx.DI2_MS = di2Hk2 * ctx.HK2_MS + di2Rt2 * ctx.RT2_MS;
    ctx.DI2_RE = di2Rt2 * ctx.RT2_RE;
  } else if (ityp === 2) {
    const cf2t = cft(ctx.HK2, ctx.RT2, ctx.M2, ctx.CFFAC ?? 1.0);
    const cf2tU2 = cf2t.cfHk * ctx.HK2_U2 + cf2t.cfRt * ctx.RT2_U2 + cf2t.cfMsq * ctx.M2_U2;
    const cf2tT2 = cf2t.cfHk * ctx.HK2_T2 + cf2t.cfRt * ctx.RT2_T2;
    const cf2tD2 = cf2t.cfHk * ctx.HK2_D2;
    const cf2tMs = cf2t.cfHk * ctx.HK2_MS + cf2t.cfRt * ctx.RT2_MS + cf2t.cfMsq * ctx.M2_MS;
    const cf2tRe = cf2t.cfRt * ctx.RT2_RE;

    ctx.DI2 = (0.5 * cf2t.cf * ctx.US2) * 2.0 / ctx.HS2;
    const di2Hs2 = -(0.5 * cf2t.cf * ctx.US2) * 2.0 / ctx.HS2 ** 2;
    const di2Us2 = (0.5 * cf2t.cf) * 2.0 / ctx.HS2;
    const di2Cf2t = (0.5 * ctx.US2) * 2.0 / ctx.HS2;

    ctx.DI2_S2 = 0.0;
    ctx.DI2_U2 = di2Hs2 * ctx.HS2_U2 + di2Us2 * ctx.US2_U2 + di2Cf2t * cf2tU2;
    ctx.DI2_T2 = di2Hs2 * ctx.HS2_T2 + di2Us2 * ctx.US2_T2 + di2Cf2t * cf2tT2;
    ctx.DI2_D2 = di2Hs2 * ctx.HS2_D2 + di2Us2 * ctx.US2_D2 + di2Cf2t * cf2tD2;
    ctx.DI2_MS = di2Hs2 * ctx.HS2_MS + di2Us2 * ctx.US2_MS + di2Cf2t * cf2tMs;
    ctx.DI2_RE = di2Hs2 * ctx.HS2_RE + di2Us2 * ctx.US2_RE + di2Cf2t * cf2tRe;

    const grt = Math.log(ctx.RT2);
    const hmin = 1.0 + 2.1 / grt;
    const hmRt2 = -(2.1 / grt ** 2) / ctx.RT2;

    const fl = (ctx.HK2 - 1.0) / (hmin - 1.0);
    const flHk2 = 1.0 / (hmin - 1.0);
    const flRt2 = (-fl / (hmin - 1.0)) * hmRt2;

    const tfl = Math.tanh(fl);
    const dfac = 0.5 + 0.5 * tfl;
    const dfFl = 0.5 * (1.0 - tfl ** 2);

    const dfHk2 = dfFl * flHk2;
    const dfRt2 = dfFl * flRt2;

    ctx.DI2_S2 = ctx.DI2_S2 * dfac;
    ctx.DI2_U2 = ctx.DI2_U2 * dfac + ctx.DI2 * (dfHk2 * ctx.HK2_U2 + dfRt2 * ctx.RT2_U2);
    ctx.DI2_T2 = ctx.DI2_T2 * dfac + ctx.DI2 * (dfHk2 * ctx.HK2_T2 + dfRt2 * ctx.RT2_T2);
    ctx.DI2_D2 = ctx.DI2_D2 * dfac + ctx.DI2 * (dfHk2 * ctx.HK2_D2);
    ctx.DI2_MS = ctx.DI2_MS * dfac + ctx.DI2 * (dfHk2 * ctx.HK2_MS + dfRt2 * ctx.RT2_MS);
    ctx.DI2_RE = ctx.DI2_RE * dfac + ctx.DI2 * (dfRt2 * ctx.RT2_RE);
    ctx.DI2 = ctx.DI2 * dfac;
  } else {
    ctx.DI2 = 0.0;
    ctx.DI2_S2 = 0.0;
    ctx.DI2_U2 = 0.0;
    ctx.DI2_T2 = 0.0;
    ctx.DI2_D2 = 0.0;
    ctx.DI2_MS = 0.0;
    ctx.DI2_RE = 0.0;
  }

  if (ityp !== 1) {
    let dd = ctx.S2 ** 2 * (0.995 - ctx.US2) * 2.0 / ctx.HS2;
  let ddHs2 = -(ctx.S2 ** 2) * (0.995 - ctx.US2) * 2.0 / ctx.HS2 ** 2;
  let ddUs2 = -(ctx.S2 ** 2) * 2.0 / ctx.HS2;
    let ddS2 = ctx.S2 * 2.0 * (0.995 - ctx.US2) * 2.0 / ctx.HS2;

    ctx.DI2 += dd;
    ctx.DI2_S2 = ddS2;
    ctx.DI2_U2 += ddHs2 * ctx.HS2_U2 + ddUs2 * ctx.US2_U2;
    ctx.DI2_T2 += ddHs2 * ctx.HS2_T2 + ddUs2 * ctx.US2_T2;
    ctx.DI2_D2 += ddHs2 * ctx.HS2_D2 + ddUs2 * ctx.US2_D2;
    ctx.DI2_MS += ddHs2 * ctx.HS2_MS + ddUs2 * ctx.US2_MS;
    ctx.DI2_RE += ddHs2 * ctx.HS2_RE + ddUs2 * ctx.US2_RE;

    dd = 0.15 * (0.995 - ctx.US2) ** 2 / ctx.RT2 * 2.0 / ctx.HS2;
    const ddUs2b = -0.15 * (0.995 - ctx.US2) * 2.0 / ctx.RT2 * 2.0 / ctx.HS2;
    const ddHs2b = -dd / ctx.HS2;
    const ddRt2 = -dd / ctx.RT2;

    ctx.DI2 += dd;
    ctx.DI2_U2 += ddHs2b * ctx.HS2_U2 + ddUs2b * ctx.US2_U2 + ddRt2 * ctx.RT2_U2;
    ctx.DI2_T2 += ddHs2b * ctx.HS2_T2 + ddUs2b * ctx.US2_T2 + ddRt2 * ctx.RT2_T2;
    ctx.DI2_D2 += ddHs2b * ctx.HS2_D2 + ddUs2b * ctx.US2_D2;
    ctx.DI2_MS += ddHs2b * ctx.HS2_MS + ddUs2b * ctx.US2_MS + ddRt2 * ctx.RT2_MS;
    ctx.DI2_RE += ddHs2b * ctx.HS2_RE + ddUs2b * ctx.US2_RE + ddRt2 * ctx.RT2_RE;
  }

  if (ityp === 2) {
    const di2l = dil(ctx.HK2, ctx.RT2);
    if (di2l.di > ctx.DI2) {
      ctx.DI2 = di2l.di;
      ctx.DI2_S2 = 0.0;
      ctx.DI2_U2 = di2l.diHk * ctx.HK2_U2 + di2l.diRt * ctx.RT2_U2;
      ctx.DI2_T2 = di2l.diHk * ctx.HK2_T2 + di2l.diRt * ctx.RT2_T2;
      ctx.DI2_D2 = di2l.diHk * ctx.HK2_D2;
      ctx.DI2_MS = di2l.diHk * ctx.HK2_MS + di2l.diRt * ctx.RT2_MS;
      ctx.DI2_RE = di2l.diRt * ctx.RT2_RE;
    }
  }

  if (ityp === 3) {
    const di2l = dilw(ctx.HK2, ctx.RT2);
    if (di2l.di > ctx.DI2) {
      ctx.DI2 = di2l.di;
      ctx.DI2_S2 = 0.0;
      ctx.DI2_U2 = di2l.diHk * ctx.HK2_U2 + di2l.diRt * ctx.RT2_U2;
      ctx.DI2_T2 = di2l.diHk * ctx.HK2_T2 + di2l.diRt * ctx.RT2_T2;
      ctx.DI2_D2 = di2l.diHk * ctx.HK2_D2;
      ctx.DI2_MS = di2l.diHk * ctx.HK2_MS + di2l.diRt * ctx.RT2_MS;
      ctx.DI2_RE = di2l.diRt * ctx.RT2_RE;
    }
  }

  if (ityp === 3) {
    ctx.DI2 *= 2.0;
    ctx.DI2_S2 *= 2.0;
    ctx.DI2_U2 *= 2.0;
    ctx.DI2_T2 *= 2.0;
    ctx.DI2_D2 *= 2.0;
    ctx.DI2_MS *= 2.0;
    ctx.DI2_RE *= 2.0;
  }

  ctx.DE2 = (3.15 + 1.72 / (ctx.HK2 - 1.0)) * ctx.T2 + ctx.D2;
  const de2Hk2 = (-1.72 / (ctx.HK2 - 1.0) ** 2) * ctx.T2;

  ctx.DE2_U2 = de2Hk2 * ctx.HK2_U2;
  ctx.DE2_T2 = de2Hk2 * ctx.HK2_T2 + (3.15 + 1.72 / (ctx.HK2 - 1.0));
  ctx.DE2_D2 = de2Hk2 * ctx.HK2_D2 + 1.0;
  ctx.DE2_MS = de2Hk2 * ctx.HK2_MS;

  const hdmax = 12.0;
  if (ctx.DE2 > hdmax * ctx.T2) {
    ctx.DE2 = hdmax * ctx.T2;
    ctx.DE2_U2 = 0.0;
    ctx.DE2_T2 = hdmax;
    ctx.DE2_D2 = 0.0;
    ctx.DE2_MS = 0.0;
  }

  syncVarsToCom(ctx, 2);
}

// Midpoint evaluation for BL integrals; stabilizes marching scheme.
function blmid(ityp, ctxIn) {
  const ctx = ctxIn || this;
  ensureCtx(ctx);

  if (ctx.SIMI) {
    ctx.HK1 = ctx.HK2;
    ctx.HK1_T1 = ctx.HK2_T2;
    ctx.HK1_D1 = ctx.HK2_D2;
    ctx.HK1_U1 = ctx.HK2_U2;
    ctx.HK1_MS = ctx.HK2_MS;
    ctx.RT1 = ctx.RT2;
    ctx.RT1_T1 = ctx.RT2_T2;
    ctx.RT1_U1 = ctx.RT2_U2;
    ctx.RT1_MS = ctx.RT2_MS;
    ctx.RT1_RE = ctx.RT2_RE;
    ctx.M1 = ctx.M2;
    ctx.M1_U1 = ctx.M2_U2;
    ctx.M1_MS = ctx.M2_MS;
  }

  const hka = 0.5 * (ctx.HK1 + ctx.HK2);
  const rta = 0.5 * (ctx.RT1 + ctx.RT2);
  const ma = 0.5 * (ctx.M1 + ctx.M2);

  if (ityp === 3) {
    ctx.CFM = 0.0;
    ctx.CFM_HKA = 0.0;
    ctx.CFM_RTA = 0.0;
    ctx.CFM_MA = 0.0;
    ctx.CFM_MS = 0.0;
  } else if (ityp === 1) {
    const cfOut = cfl(hka, rta, ma);
    ctx.CFM = cfOut.cf;
    ctx.CFM_HKA = cfOut.cfHk;
    ctx.CFM_RTA = cfOut.cfRt;
    ctx.CFM_MA = cfOut.cfMsq;
  } else {
    let cfOut = cft(hka, rta, ma, ctx.CFFAC ?? 1.0);
    ctx.CFM = cfOut.cf;
    ctx.CFM_HKA = cfOut.cfHk;
    ctx.CFM_RTA = cfOut.cfRt;
    ctx.CFM_MA = cfOut.cfMsq;
    const cfml = cfl(hka, rta, ma);
    if (cfml.cf > ctx.CFM) {
      ctx.CFM = cfml.cf;
      ctx.CFM_HKA = cfml.cfHk;
      ctx.CFM_RTA = cfml.cfRt;
      ctx.CFM_MA = cfml.cfMsq;
    }
  }

  ctx.CFM_U1 = 0.5 * (ctx.CFM_HKA * ctx.HK1_U1 + ctx.CFM_MA * ctx.M1_U1 + ctx.CFM_RTA * ctx.RT1_U1);
  ctx.CFM_T1 = 0.5 * (ctx.CFM_HKA * ctx.HK1_T1 + ctx.CFM_RTA * ctx.RT1_T1);
  ctx.CFM_D1 = 0.5 * (ctx.CFM_HKA * ctx.HK1_D1);

  ctx.CFM_U2 = 0.5 * (ctx.CFM_HKA * ctx.HK2_U2 + ctx.CFM_MA * ctx.M2_U2 + ctx.CFM_RTA * ctx.RT2_U2);
  ctx.CFM_T2 = 0.5 * (ctx.CFM_HKA * ctx.HK2_T2 + ctx.CFM_RTA * ctx.RT2_T2);
  ctx.CFM_D2 = 0.5 * (ctx.CFM_HKA * ctx.HK2_D2);

  ctx.CFM_MS = 0.5 * (ctx.CFM_HKA * ctx.HK1_MS + ctx.CFM_MA * ctx.M1_MS + ctx.CFM_RTA * ctx.RT1_MS
    + ctx.CFM_HKA * ctx.HK2_MS + ctx.CFM_MA * ctx.M2_MS + ctx.CFM_RTA * ctx.RT2_MS);
  ctx.CFM_RE = 0.5 * (ctx.CFM_RTA * ctx.RT1_RE + ctx.CFM_RTA * ctx.RT2_RE);
}

// Transition differential system used in Newton updates.
function trdif(ctxIn) {
  const ctx = ctxIn || this;
  ensureCtx(ctx);
  const bl1 = create2d(4, 5);
  const bl2 = create2d(4, 5);
  const blrez = new Array(5).fill(0.0);
  const blm = new Array(5).fill(0.0);
  const blr = new Array(5).fill(0.0);
  const blx = new Array(5).fill(0.0);
  const bt1 = create2d(4, 5);
  const bt2 = create2d(4, 5);
  const btrez = new Array(5).fill(0.0);
  const btm = new Array(5).fill(0.0);
  const btr = new Array(5).fill(0.0);
  const btx = new Array(5).fill(0.0);

  for (let icom = 1; icom <= NCOM; icom += 1) {
    ctx.C1SAV[icom] = ctx.COM1[icom];
    ctx.C2SAV[icom] = ctx.COM2[icom];
  }

  let wf2 = (ctx.XT - ctx.X1) / (ctx.X2 - ctx.X1);
  const wf2Xt = 1.0 / (ctx.X2 - ctx.X1);

  const wf2A1 = wf2Xt * ctx.XT_A1;
  const wf2X1 = wf2Xt * ctx.XT_X1 + (wf2 - 1.0) / (ctx.X2 - ctx.X1);
  const wf2X2 = wf2Xt * ctx.XT_X2 - wf2 / (ctx.X2 - ctx.X1);
  const wf2T1 = wf2Xt * ctx.XT_T1;
  const wf2T2 = wf2Xt * ctx.XT_T2;
  const wf2D1 = wf2Xt * ctx.XT_D1;
  const wf2D2 = wf2Xt * ctx.XT_D2;
  const wf2U1 = wf2Xt * ctx.XT_U1;
  const wf2U2 = wf2Xt * ctx.XT_U2;
  const wf2Ms = wf2Xt * ctx.XT_MS;
  const wf2Re = wf2Xt * ctx.XT_RE;
  const wf2Xf = wf2Xt * ctx.XT_XF;

  const wf1 = 1.0 - wf2;
  const wf1A1 = -wf2A1;
  const wf1X1 = -wf2X1;
  const wf1X2 = -wf2X2;
  const wf1T1 = -wf2T1;
  const wf1T2 = -wf2T2;
  const wf1D1 = -wf2D1;
  const wf1D2 = -wf2D2;
  const wf1U1 = -wf2U1;
  const wf1U2 = -wf2U2;
  const wf1Ms = -wf2Ms;
  const wf1Re = -wf2Re;
  const wf1Xf = -wf2Xf;

  const tt = ctx.T1 * wf1 + ctx.T2 * wf2;
  const ttA1 = ctx.T1 * wf1A1 + ctx.T2 * wf2A1;
  const ttX1 = ctx.T1 * wf1X1 + ctx.T2 * wf2X1;
  const ttX2 = ctx.T1 * wf1X2 + ctx.T2 * wf2X2;
  const ttT1 = ctx.T1 * wf1T1 + ctx.T2 * wf2T1 + wf1;
  const ttT2 = ctx.T1 * wf1T2 + ctx.T2 * wf2T2 + wf2;
  const ttD1 = ctx.T1 * wf1D1 + ctx.T2 * wf2D1;
  const ttD2 = ctx.T1 * wf1D2 + ctx.T2 * wf2D2;
  const ttU1 = ctx.T1 * wf1U1 + ctx.T2 * wf2U1;
  const ttU2 = ctx.T1 * wf1U2 + ctx.T2 * wf2U2;
  const ttMs = ctx.T1 * wf1Ms + ctx.T2 * wf2Ms;
  const ttRe = ctx.T1 * wf1Re + ctx.T2 * wf2Re;
  const ttXf = ctx.T1 * wf1Xf + ctx.T2 * wf2Xf;

  const dt = ctx.D1 * wf1 + ctx.D2 * wf2;
  const dtA1 = ctx.D1 * wf1A1 + ctx.D2 * wf2A1;
  const dtX1 = ctx.D1 * wf1X1 + ctx.D2 * wf2X1;
  const dtX2 = ctx.D1 * wf1X2 + ctx.D2 * wf2X2;
  const dtT1 = ctx.D1 * wf1T1 + ctx.D2 * wf2T1;
  const dtT2 = ctx.D1 * wf1T2 + ctx.D2 * wf2T2;
  const dtD1 = ctx.D1 * wf1D1 + ctx.D2 * wf2D1 + wf1;
  const dtD2 = ctx.D1 * wf1D2 + ctx.D2 * wf2D2 + wf2;
  const dtU1 = ctx.D1 * wf1U1 + ctx.D2 * wf2U1;
  const dtU2 = ctx.D1 * wf1U2 + ctx.D2 * wf2U2;
  const dtMs = ctx.D1 * wf1Ms + ctx.D2 * wf2Ms;
  const dtRe = ctx.D1 * wf1Re + ctx.D2 * wf2Re;
  const dtXf = ctx.D1 * wf1Xf + ctx.D2 * wf2Xf;

  const ut = ctx.U1 * wf1 + ctx.U2 * wf2;
  const utA1 = ctx.U1 * wf1A1 + ctx.U2 * wf2A1;
  const utX1 = ctx.U1 * wf1X1 + ctx.U2 * wf2X1;
  const utX2 = ctx.U1 * wf1X2 + ctx.U2 * wf2X2;
  const utT1 = ctx.U1 * wf1T1 + ctx.U2 * wf2T1;
  const utT2 = ctx.U1 * wf1T2 + ctx.U2 * wf2T2;
  const utD1 = ctx.U1 * wf1D1 + ctx.U2 * wf2D1;
  const utD2 = ctx.U1 * wf1D2 + ctx.U2 * wf2D2;
  const utU1 = ctx.U1 * wf1U1 + ctx.U2 * wf2U1 + wf1;
  const utU2 = ctx.U1 * wf1U2 + ctx.U2 * wf2U2 + wf2;
  const utMs = ctx.U1 * wf1Ms + ctx.U2 * wf2Ms;
  const utRe = ctx.U1 * wf1Re + ctx.U2 * wf2Re;
  const utXf = ctx.U1 * wf1Xf + ctx.U2 * wf2Xf;

  ctx.X2 = ctx.XT;
  ctx.T2 = tt;
  ctx.D2 = dt;
  ctx.U2 = ut;

  ctx.AMPL2 = ctx.AMCRIT;
  ctx.S2 = 0.0;
  syncVarsToCom(ctx, 2);

  blkin(ctx);
  blvar(1, ctx);
  blmid(1, ctx);
  bldif(1, ctx);

  for (let k = 2; k <= 3; k += 1) {
    blrez[k] = ctx.VSREZ[k];
    blm[k] = ctx.VSM[k]
      + ctx.VS2[k][2] * ttMs
      + ctx.VS2[k][3] * dtMs
      + ctx.VS2[k][4] * utMs
      + ctx.VS2[k][5] * ctx.XT_MS;
    blr[k] = ctx.VSR[k]
      + ctx.VS2[k][2] * ttRe
      + ctx.VS2[k][3] * dtRe
      + ctx.VS2[k][4] * utRe
      + ctx.VS2[k][5] * ctx.XT_RE;
    blx[k] = ctx.VSX[k]
      + ctx.VS2[k][2] * ttXf
      + ctx.VS2[k][3] * dtXf
      + ctx.VS2[k][4] * utXf
      + ctx.VS2[k][5] * ctx.XT_XF;

    bl1[k][1] = ctx.VS1[k][1]
      + ctx.VS2[k][2] * ttA1
      + ctx.VS2[k][3] * dtA1
      + ctx.VS2[k][4] * utA1
      + ctx.VS2[k][5] * ctx.XT_A1;
    bl1[k][2] = ctx.VS1[k][2]
      + ctx.VS2[k][2] * ttT1
      + ctx.VS2[k][3] * dtT1
      + ctx.VS2[k][4] * utT1
      + ctx.VS2[k][5] * ctx.XT_T1;
    bl1[k][3] = ctx.VS1[k][3]
      + ctx.VS2[k][2] * ttD1
      + ctx.VS2[k][3] * dtD1
      + ctx.VS2[k][4] * utD1
      + ctx.VS2[k][5] * ctx.XT_D1;
    bl1[k][4] = ctx.VS1[k][4]
      + ctx.VS2[k][2] * ttU1
      + ctx.VS2[k][3] * dtU1
      + ctx.VS2[k][4] * utU1
      + ctx.VS2[k][5] * ctx.XT_U1;
    bl1[k][5] = ctx.VS1[k][5]
      + ctx.VS2[k][2] * ttX1
      + ctx.VS2[k][3] * dtX1
      + ctx.VS2[k][4] * utX1
      + ctx.VS2[k][5] * ctx.XT_X1;

    bl2[k][1] = 0.0;
    bl2[k][2] = ctx.VS2[k][2] * ttT2
      + ctx.VS2[k][3] * dtT2
      + ctx.VS2[k][4] * utT2
      + ctx.VS2[k][5] * ctx.XT_T2;
    bl2[k][3] = ctx.VS2[k][2] * ttD2
      + ctx.VS2[k][3] * dtD2
      + ctx.VS2[k][4] * utD2
      + ctx.VS2[k][5] * ctx.XT_D2;
    bl2[k][4] = ctx.VS2[k][2] * ttU2
      + ctx.VS2[k][3] * dtU2
      + ctx.VS2[k][4] * utU2
      + ctx.VS2[k][5] * ctx.XT_U2;
    bl2[k][5] = ctx.VS2[k][2] * ttX2
      + ctx.VS2[k][3] * dtX2
      + ctx.VS2[k][4] * utX2
      + ctx.VS2[k][5] * ctx.XT_X2;
  }

  blvar(2, ctx);
  const ctr = ctx.CTRCON * Math.exp(-ctx.CTRCEX / (ctx.HK2 - 1.0));
  const ctrHk2 = ctr * ctx.CTRCEX / (ctx.HK2 - 1.0) ** 2;

  const st = ctr * ctx.CQ2;
  const stTt = ctr * ctx.CQ2_T2 + ctx.CQ2 * ctrHk2 * ctx.HK2_T2;
  const stDt = ctr * ctx.CQ2_D2 + ctx.CQ2 * ctrHk2 * ctx.HK2_D2;
  const stUt = ctr * ctx.CQ2_U2 + ctx.CQ2 * ctrHk2 * ctx.HK2_U2;
  const stMs = ctr * ctx.CQ2_MS + ctx.CQ2 * ctrHk2 * ctx.HK2_MS;
  const stRe = ctr * ctx.CQ2_RE;

  const stA1 = stTt * ttA1 + stDt * dtA1 + stUt * utA1;
  const stX1 = stTt * ttX1 + stDt * dtX1 + stUt * utX1;
  const stX2 = stTt * ttX2 + stDt * dtX2 + stUt * utX2;
  const stT1 = stTt * ttT1 + stDt * dtT1 + stUt * utT1;
  const stT2 = stTt * ttT2 + stDt * dtT2 + stUt * utT2;
  const stD1 = stTt * ttD1 + stDt * dtD1 + stUt * utD1;
  const stD2 = stTt * ttD2 + stDt * dtD2 + stUt * utD2;
  const stU1 = stTt * ttU1 + stDt * dtU1 + stUt * utU1;
  const stU2 = stTt * ttU2 + stDt * dtU2 + stUt * utU2;
  const stMs2 = stTt * ttMs + stDt * dtMs + stUt * utMs + stMs;
  const stRe2 = stTt * ttRe + stDt * dtRe + stUt * utRe + stRe;
  const stXf = stTt * ttXf + stDt * dtXf + stUt * utXf;

  ctx.AMPL2 = 0.0;
  ctx.S2 = st;

  blvar(2, ctx);

  for (let icom = 1; icom <= NCOM; icom += 1) {
    ctx.COM1[icom] = ctx.COM2[icom];
    ctx.COM2[icom] = ctx.C2SAV[icom];
  }
  syncComToVars(ctx, 1);
  syncComToVars(ctx, 2);

  blmid(2, ctx);
  bldif(2, ctx);

  for (let k = 1; k <= 3; k += 1) {
    btrez[k] = ctx.VSREZ[k];
    btm[k] = ctx.VSM[k]
      + ctx.VS1[k][1] * stMs2
      + ctx.VS1[k][2] * ttMs
      + ctx.VS1[k][3] * dtMs
      + ctx.VS1[k][4] * utMs
      + ctx.VS1[k][5] * ctx.XT_MS;
    btr[k] = ctx.VSR[k]
      + ctx.VS1[k][1] * stRe2
      + ctx.VS1[k][2] * ttRe
      + ctx.VS1[k][3] * dtRe
      + ctx.VS1[k][4] * utRe
      + ctx.VS1[k][5] * ctx.XT_RE;
    btx[k] = ctx.VSX[k]
      + ctx.VS1[k][1] * stXf
      + ctx.VS1[k][2] * ttXf
      + ctx.VS1[k][3] * dtXf
      + ctx.VS1[k][4] * utXf
      + ctx.VS1[k][5] * ctx.XT_XF;

    bt1[k][1] = ctx.VS1[k][1] * stA1
      + ctx.VS1[k][2] * ttA1
      + ctx.VS1[k][3] * dtA1
      + ctx.VS1[k][4] * utA1
      + ctx.VS1[k][5] * ctx.XT_A1;
    bt1[k][2] = ctx.VS1[k][1] * stT1
      + ctx.VS1[k][2] * ttT1
      + ctx.VS1[k][3] * dtT1
      + ctx.VS1[k][4] * utT1
      + ctx.VS1[k][5] * ctx.XT_T1;
    bt1[k][3] = ctx.VS1[k][1] * stD1
      + ctx.VS1[k][2] * ttD1
      + ctx.VS1[k][3] * dtD1
      + ctx.VS1[k][4] * utD1
      + ctx.VS1[k][5] * ctx.XT_D1;
    bt1[k][4] = ctx.VS1[k][1] * stU1
      + ctx.VS1[k][2] * ttU1
      + ctx.VS1[k][3] * dtU1
      + ctx.VS1[k][4] * utU1
      + ctx.VS1[k][5] * ctx.XT_U1;
    bt1[k][5] = ctx.VS1[k][1] * stX1
      + ctx.VS1[k][2] * ttX1
      + ctx.VS1[k][3] * dtX1
      + ctx.VS1[k][4] * utX1
      + ctx.VS1[k][5] * ctx.XT_X1;

    bt2[k][1] = ctx.VS2[k][1];
    bt2[k][2] = ctx.VS2[k][2]
      + ctx.VS1[k][1] * stT2
      + ctx.VS1[k][2] * ttT2
      + ctx.VS1[k][3] * dtT2
      + ctx.VS1[k][4] * utT2
      + ctx.VS1[k][5] * ctx.XT_T2;
    bt2[k][3] = ctx.VS2[k][3]
      + ctx.VS1[k][1] * stD2
      + ctx.VS1[k][2] * ttD2
      + ctx.VS1[k][3] * dtD2
      + ctx.VS1[k][4] * utD2
      + ctx.VS1[k][5] * ctx.XT_D2;
    bt2[k][4] = ctx.VS2[k][4]
      + ctx.VS1[k][1] * stU2
      + ctx.VS1[k][2] * ttU2
      + ctx.VS1[k][3] * dtU2
      + ctx.VS1[k][4] * utU2
      + ctx.VS1[k][5] * ctx.XT_U2;
    bt2[k][5] = ctx.VS2[k][5]
      + ctx.VS1[k][1] * stX2
      + ctx.VS1[k][2] * ttX2
      + ctx.VS1[k][3] * dtX2
      + ctx.VS1[k][4] * utX2
      + ctx.VS1[k][5] * ctx.XT_X2;
  }

  ctx.VSREZ[1] = btrez[1];
  ctx.VSREZ[2] = blrez[2] + btrez[2];
  ctx.VSREZ[3] = blrez[3] + btrez[3];
  ctx.VSM[1] = btm[1];
  ctx.VSM[2] = blm[2] + btm[2];
  ctx.VSM[3] = blm[3] + btm[3];
  ctx.VSR[1] = btr[1];
  ctx.VSR[2] = blr[2] + btr[2];
  ctx.VSR[3] = blr[3] + btr[3];
  ctx.VSX[1] = btx[1];
  ctx.VSX[2] = blx[2] + btx[2];
  ctx.VSX[3] = blx[3] + btx[3];
  for (let l = 1; l <= 5; l += 1) {
    ctx.VS1[1][l] = bt1[1][l];
    ctx.VS2[1][l] = bt2[1][l];
    ctx.VS1[2][l] = bl1[2][l] + bt1[2][l];
    ctx.VS2[2][l] = bl2[2][l] + bt2[2][l];
    ctx.VS1[3][l] = bl1[3][l] + bt1[3][l];
    ctx.VS2[3][l] = bl2[3][l] + bt2[3][l];
  }

  for (let icom = 1; icom <= NCOM; icom += 1) {
    ctx.COM1[icom] = ctx.C1SAV[icom];
  }
  syncComToVars(ctx, 1);
}

// Differential equations for BL integrals (laminar/turbulent).
function bldif(ityp, ctxIn) {
  const ctx = ctxIn || this;
  ensureCtx(ctx);

  let xlog;
  let ulog;
  let tlog;
  let hlog;
  let ddlog;

  if (ityp === 0) {
    xlog = 1.0;
    ulog = ctx.BULE;
    tlog = 0.5 * (1.0 - ctx.BULE);
    hlog = 0.0;
    ddlog = 0.0;
  } else {
    xlog = Math.log(ctx.X2 / ctx.X1);
    ulog = Math.log(ctx.U2 / ctx.U1);
    tlog = Math.log(ctx.T2 / ctx.T1);
    hlog = Math.log(ctx.HS2 / ctx.HS1);
    ddlog = 1.0;
  }

  for (let k = 1; k <= 4; k += 1) {
    ctx.VSREZ[k] = 0.0;
    ctx.VSM[k] = 0.0;
    ctx.VSR[k] = 0.0;
    ctx.VSX[k] = 0.0;
    for (let l = 1; l <= 5; l += 1) {
      ctx.VS1[k][l] = 0.0;
      ctx.VS2[k][l] = 0.0;
    }
  }

  const hupwt = 1.0;
  let hdcon = 5.0 * hupwt / ctx.HK2 ** 2;
  let hdHk1 = 0.0;
  let hdHk2 = -hdcon * 2.0 / ctx.HK2;

  if (ityp === 3) {
    hdcon = hupwt / ctx.HK2 ** 2;
    hdHk1 = 0.0;
    hdHk2 = -hdcon * 2.0 / ctx.HK2;
  }

  const arg = Math.abs((ctx.HK2 - 1.0) / (ctx.HK1 - 1.0));
  const hl = Math.log(arg);
  const hlHk1 = -1.0 / (ctx.HK1 - 1.0);
  const hlHk2 = 1.0 / (ctx.HK2 - 1.0);

  const hlsq = Math.min(hl ** 2, 15.0);
  const ehh = Math.exp(-hlsq * hdcon);
  const upw = 1.0 - 0.5 * ehh;
  const upwHl = ehh * hl * hdcon;
  const upwHd = 0.5 * ehh * hlsq;

  const upwHk1 = upwHl * hlHk1 + upwHd * hdHk1;
  const upwHk2 = upwHl * hlHk2 + upwHd * hdHk2;

  const upwU1 = upwHk1 * ctx.HK1_U1;
  const upwT1 = upwHk1 * ctx.HK1_T1;
  const upwD1 = upwHk1 * ctx.HK1_D1;
  const upwU2 = upwHk2 * ctx.HK2_U2;
  const upwT2 = upwHk2 * ctx.HK2_T2;
  const upwD2 = upwHk2 * ctx.HK2_D2;
  const upwMs = upwHk1 * ctx.HK1_MS + upwHk2 * ctx.HK2_MS;

  if (ityp === 0) {
    ctx.VS2[1][1] = 1.0;
    ctx.VSR[1] = 0.0;
    ctx.VSREZ[1] = -ctx.AMPL2;
  } else if (ityp === 1) {
    const axOut = axset(ctx.HK1, ctx.T1, ctx.RT1, ctx.AMPL1,
      ctx.HK2, ctx.T2, ctx.RT2, ctx.AMPL2, ctx.AMCRIT, ctx.IDAMPV);
    const ax = axOut.ax;
    const rezc = ctx.AMPL2 - ctx.AMPL1 - ax * (ctx.X2 - ctx.X1);
    const zAx = -(ctx.X2 - ctx.X1);

    ctx.VS1[1][1] = zAx * axOut.axA1 - 1.0;
    ctx.VS1[1][2] = zAx * (axOut.axHk1 * ctx.HK1_T1 + axOut.axT1 + axOut.axRt1 * ctx.RT1_T1);
    ctx.VS1[1][3] = zAx * (axOut.axHk1 * ctx.HK1_D1);
    ctx.VS1[1][4] = zAx * (axOut.axHk1 * ctx.HK1_U1 + axOut.axRt1 * ctx.RT1_U1);
    ctx.VS1[1][5] = ax;
    ctx.VS2[1][1] = zAx * axOut.axA2 + 1.0;
    ctx.VS2[1][2] = zAx * (axOut.axHk2 * ctx.HK2_T2 + axOut.axT2 + axOut.axRt2 * ctx.RT2_T2);
    ctx.VS2[1][3] = zAx * (axOut.axHk2 * ctx.HK2_D2);
    ctx.VS2[1][4] = zAx * (axOut.axHk2 * ctx.HK2_U2 + axOut.axRt2 * ctx.RT2_U2);
    ctx.VS2[1][5] = -ax;
    ctx.VSM[1] = zAx * (axOut.axHk1 * ctx.HK1_MS + axOut.axRt1 * ctx.RT1_MS
      + axOut.axHk2 * ctx.HK2_MS + axOut.axRt2 * ctx.RT2_MS);
    ctx.VSR[1] = zAx * (axOut.axRt1 * ctx.RT1_RE + axOut.axRt2 * ctx.RT2_RE);
    ctx.VSX[1] = 0.0;
    ctx.VSREZ[1] = -rezc;
  } else {
    const sa = (1.0 - upw) * ctx.S1 + upw * ctx.S2;
    const cqa = (1.0 - upw) * ctx.CQ1 + upw * ctx.CQ2;
    const cfa = (1.0 - upw) * ctx.CF1 + upw * ctx.CF2;
    const hka = (1.0 - upw) * ctx.HK1 + upw * ctx.HK2;

    const usa = 0.5 * (ctx.US1 + ctx.US2);
    const rta = 0.5 * (ctx.RT1 + ctx.RT2);
    const dea = 0.5 * (ctx.DE1 + ctx.DE2);
    const da = 0.5 * (ctx.D1 + ctx.D2);

    const ald = ityp === 3 ? ctx.DLCON : 1.0;

    let gcc = 0.0;
    let hkc = hka - 1.0;
    let hkcHka = 1.0;
    let hkcRta = 0.0;
    if (ityp === 2) {
      gcc = ctx.GCCON;
      hkc = hka - 1.0 - gcc / rta;
      hkcHka = 1.0;
      hkcRta = gcc / rta ** 2;
      if (hkc < 0.01) {
        hkc = 0.01;
        hkcHka = 0.0;
        hkcRta = 0.0;
      }
    }

    const hr = hkc / (ctx.GACON * ald * hka);
    const hrHka = hkcHka / (ctx.GACON * ald * hka) - hr / hka;
    const hrRta = hkcRta / (ctx.GACON * ald * hka);

    const uq = (0.5 * cfa - hr ** 2) / (ctx.GBCON * da);
    const uqHka = -2.0 * hr * hrHka / (ctx.GBCON * da);
    const uqRta = -2.0 * hr * hrRta / (ctx.GBCON * da);
    const uqCfa = 0.5 / (ctx.GBCON * da);
    const uqDa = -uq / da;
    const uqUpw = uqCfa * (ctx.CF2 - ctx.CF1) + uqHka * (ctx.HK2 - ctx.HK1);

    let uqT1 = (1.0 - upw) * (uqCfa * ctx.CF1_T1 + uqHka * ctx.HK1_T1) + uqUpw * upwT1;
    let uqD1 = (1.0 - upw) * (uqCfa * ctx.CF1_D1 + uqHka * ctx.HK1_D1) + uqUpw * upwD1;
    let uqU1 = (1.0 - upw) * (uqCfa * ctx.CF1_U1 + uqHka * ctx.HK1_U1) + uqUpw * upwU1;
    let uqT2 = upw * (uqCfa * ctx.CF2_T2 + uqHka * ctx.HK2_T2) + uqUpw * upwT2;
    let uqD2 = upw * (uqCfa * ctx.CF2_D2 + uqHka * ctx.HK2_D2) + uqUpw * upwD2;
    let uqU2 = upw * (uqCfa * ctx.CF2_U2 + uqHka * ctx.HK2_U2) + uqUpw * upwU2;
    let uqMs = (1.0 - upw) * (uqCfa * ctx.CF1_MS + uqHka * ctx.HK1_MS) + uqUpw * upwMs
      + upw * (uqCfa * ctx.CF2_MS + uqHka * ctx.HK2_MS);
    let uqRe = (1.0 - upw) * uqCfa * ctx.CF1_RE + upw * uqCfa * ctx.CF2_RE;

    uqT1 += 0.5 * uqRta * ctx.RT1_T1;
    uqD1 += 0.5 * uqDa;
    uqU1 += 0.5 * uqRta * ctx.RT1_U1;
    uqT2 += 0.5 * uqRta * ctx.RT2_T2;
    uqD2 += 0.5 * uqDa;
    uqU2 += 0.5 * uqRta * ctx.RT2_U2;
    uqMs += 0.5 * uqRta * ctx.RT1_MS + 0.5 * uqRta * ctx.RT2_MS;
    uqRe += 0.5 * uqRta * ctx.RT1_RE + 0.5 * uqRta * ctx.RT2_RE;

    const scc = ctx.SCCON * 1.333 / (1.0 + usa);
    const sccUsa = -scc / (1.0 + usa);

    const sccUs1 = sccUsa * 0.5;
    const sccUs2 = sccUsa * 0.5;

    const slog = Math.log(ctx.S2 / ctx.S1);
    const dxi = ctx.X2 - ctx.X1;

    const rezc = scc * (cqa - sa * ald) * dxi
      - dea * 2.0 * slog
      + dea * 2.0 * (uq * dxi - ulog) * ctx.DUXCON;

    const zCfa = dea * 2.0 * uqCfa * dxi * ctx.DUXCON;
    const zHka = dea * 2.0 * uqHka * dxi * ctx.DUXCON;
    const zDa = dea * 2.0 * uqDa * dxi * ctx.DUXCON;
    const zSl = -dea * 2.0;
    const zUl = -dea * 2.0 * ctx.DUXCON;
    const zDxi = scc * (cqa - sa * ald) + dea * 2.0 * uq * ctx.DUXCON;
    const zUsa = sccUsa * (cqa - sa * ald) * dxi;
    const zCqa = scc * dxi;
    const zSa = -scc * dxi * ald;
    const zDea = 2.0 * ((uq * dxi - ulog) * ctx.DUXCON - slog);

    const zUpw = zCqa * (ctx.CQ2 - ctx.CQ1) + zSa * (ctx.S2 - ctx.S1)
      + zCfa * (ctx.CF2 - ctx.CF1) + zHka * (ctx.HK2 - ctx.HK1);
    const zDe1 = 0.5 * zDea;
    const zDe2 = 0.5 * zDea;
    const zUs1 = 0.5 * zUsa;
    const zUs2 = 0.5 * zUsa;
    const zD1 = 0.5 * zDa;
    const zD2 = 0.5 * zDa;
    const zU1 = -zUl / ctx.U1;
    const zU2 = zUl / ctx.U2;
    const zX1 = -zDxi;
    const zX2 = zDxi;
    const zS1 = (1.0 - upw) * zSa - zSl / ctx.S1;
    const zS2 = upw * zSa + zSl / ctx.S2;
    const zCq1 = (1.0 - upw) * zCqa;
    const zCq2 = upw * zCqa;
    const zCf1 = (1.0 - upw) * zCfa;
    const zCf2 = upw * zCfa;
    const zHk1 = (1.0 - upw) * zHka;
    const zHk2 = upw * zHka;

    ctx.VS1[1][1] = zS1;
    ctx.VS1[1][2] = zUpw * upwT1 + zDe1 * ctx.DE1_T1 + zUs1 * ctx.US1_T1;
    ctx.VS1[1][3] = zD1 + zUpw * upwD1 + zDe1 * ctx.DE1_D1 + zUs1 * ctx.US1_D1;
    ctx.VS1[1][4] = zU1 + zUpw * upwU1 + zDe1 * ctx.DE1_U1 + zUs1 * ctx.US1_U1;
    ctx.VS1[1][5] = zX1;
    ctx.VS2[1][1] = zS2;
    ctx.VS2[1][2] = zUpw * upwT2 + zDe2 * ctx.DE2_T2 + zUs2 * ctx.US2_T2;
    ctx.VS2[1][3] = zD2 + zUpw * upwD2 + zDe2 * ctx.DE2_D2 + zUs2 * ctx.US2_D2;
    ctx.VS2[1][4] = zU2 + zUpw * upwU2 + zDe2 * ctx.DE2_U2 + zUs2 * ctx.US2_U2;
    ctx.VS2[1][5] = zX2;
    ctx.VSM[1] = zUpw * upwMs + zDe1 * ctx.DE1_MS + zUs1 * ctx.US1_MS
      + zDe2 * ctx.DE2_MS + zUs2 * ctx.US2_MS;

    ctx.VS1[1][2] += zCq1 * ctx.CQ1_T1 + zCf1 * ctx.CF1_T1 + zHk1 * ctx.HK1_T1;
    ctx.VS1[1][3] += zCq1 * ctx.CQ1_D1 + zCf1 * ctx.CF1_D1 + zHk1 * ctx.HK1_D1;
    ctx.VS1[1][4] += zCq1 * ctx.CQ1_U1 + zCf1 * ctx.CF1_U1 + zHk1 * ctx.HK1_U1;

    ctx.VS2[1][2] += zCq2 * ctx.CQ2_T2 + zCf2 * ctx.CF2_T2 + zHk2 * ctx.HK2_T2;
    ctx.VS2[1][3] += zCq2 * ctx.CQ2_D2 + zCf2 * ctx.CF2_D2 + zHk2 * ctx.HK2_D2;
    ctx.VS2[1][4] += zCq2 * ctx.CQ2_U2 + zCf2 * ctx.CF2_U2 + zHk2 * ctx.HK2_U2;

    ctx.VSM[1] += zCq1 * ctx.CQ1_MS + zCf1 * ctx.CF1_MS + zHk1 * ctx.HK1_MS
      + zCq2 * ctx.CQ2_MS + zCf2 * ctx.CF2_MS + zHk2 * ctx.HK2_MS;
    ctx.VSR[1] = zCq1 * ctx.CQ1_RE + zCf1 * ctx.CF1_RE + zCq2 * ctx.CQ2_RE + zCf2 * ctx.CF2_RE;
    ctx.VSX[1] = 0.0;
    ctx.VSREZ[1] = -rezc;
  }

  const ha = 0.5 * (ctx.H1 + ctx.H2);
  const ma = 0.5 * (ctx.M1 + ctx.M2);
  const xa = 0.5 * (ctx.X1 + ctx.X2);
  const ta = 0.5 * (ctx.T1 + ctx.T2);
  const hwa = 0.5 * (ctx.DW1 / ctx.T1 + ctx.DW2 / ctx.T2);

  let cfx = 0.50 * ctx.CFM * xa / ta + 0.25 * (ctx.CF1 * ctx.X1 / ctx.T1 + ctx.CF2 * ctx.X2 / ctx.T2);
  const cfxXa = 0.50 * ctx.CFM / ta;
  const cfxTa = -0.50 * ctx.CFM * xa / ta ** 2;

  const cfxX1 = 0.25 * ctx.CF1 / ctx.T1 + cfxXa * 0.5;
  const cfxX2 = 0.25 * ctx.CF2 / ctx.T2 + cfxXa * 0.5;
  const cfxT1 = -0.25 * ctx.CF1 * ctx.X1 / ctx.T1 ** 2 + cfxTa * 0.5;
  const cfxT2 = -0.25 * ctx.CF2 * ctx.X2 / ctx.T2 ** 2 + cfxTa * 0.5;
  const cfxCf1 = 0.25 * ctx.X1 / ctx.T1;
  const cfxCf2 = 0.25 * ctx.X2 / ctx.T2;
  const cfxCfm = 0.50 * xa / ta;

  const btmp = ha + 2.0 - ma + hwa;
  const rezt = tlog + btmp * ulog - xlog * 0.5 * cfx;
  const zCfx = -xlog * 0.5;
  const zHa = ulog;
  const zHwa = ulog;
  const zMa = -ulog;
  const zXl = -ddlog * 0.5 * cfx;
  const zUl = ddlog * btmp;
  const zTl = ddlog;

  const zCfm = zCfx * cfxCfm;
  const zCf1 = zCfx * cfxCf1;
  const zCf2 = zCfx * cfxCf2;

  let zT1 = -zTl / ctx.T1 + zCfx * cfxT1 + zHwa * 0.5 * (-ctx.DW1 / ctx.T1 ** 2);
  let zT2 = zTl / ctx.T2 + zCfx * cfxT2 + zHwa * 0.5 * (-ctx.DW2 / ctx.T2 ** 2);
  let zX1 = -zXl / ctx.X1 + zCfx * cfxX1;
  let zX2 = zXl / ctx.X2 + zCfx * cfxX2;
  const zU1 = -zUl / ctx.U1;
  const zU2 = zUl / ctx.U2;

  ctx.VS1[2][2] = 0.5 * zHa * ctx.H1_T1 + zCfm * ctx.CFM_T1 + zCf1 * ctx.CF1_T1 + zT1;
  ctx.VS1[2][3] = 0.5 * zHa * ctx.H1_D1 + zCfm * ctx.CFM_D1 + zCf1 * ctx.CF1_D1;
  ctx.VS1[2][4] = 0.5 * zMa * ctx.M1_U1 + zCfm * ctx.CFM_U1 + zCf1 * ctx.CF1_U1 + zU1;
  ctx.VS1[2][5] = zX1;
  ctx.VS2[2][2] = 0.5 * zHa * ctx.H2_T2 + zCfm * ctx.CFM_T2 + zCf2 * ctx.CF2_T2 + zT2;
  ctx.VS2[2][3] = 0.5 * zHa * ctx.H2_D2 + zCfm * ctx.CFM_D2 + zCf2 * ctx.CF2_D2;
  ctx.VS2[2][4] = 0.5 * zMa * ctx.M2_U2 + zCfm * ctx.CFM_U2 + zCf2 * ctx.CF2_U2 + zU2;
  ctx.VS2[2][5] = zX2;

  ctx.VSM[2] = 0.5 * zMa * ctx.M1_MS + zCfm * ctx.CFM_MS + zCf1 * ctx.CF1_MS
    + 0.5 * zMa * ctx.M2_MS + zCf2 * ctx.CF2_MS;
  ctx.VSR[2] = zCfm * ctx.CFM_RE + zCf1 * ctx.CF1_RE + zCf2 * ctx.CF2_RE;
  ctx.VSX[2] = 0.0;
  ctx.VSREZ[2] = -rezt;

  const xot1 = ctx.X1 / ctx.T1;
  const xot2 = ctx.X2 / ctx.T2;

  const ha2 = 0.5 * (ctx.H1 + ctx.H2);
  const hsa = 0.5 * (ctx.HS1 + ctx.HS2);
  const hca = 0.5 * (ctx.HC1 + ctx.HC2);
  const hwa2 = 0.5 * (ctx.DW1 / ctx.T1 + ctx.DW2 / ctx.T2);

  const dix = (1.0 - upw) * ctx.DI1 * xot1 + upw * ctx.DI2 * xot2;
  const cfx2 = (1.0 - upw) * ctx.CF1 * xot1 + upw * ctx.CF2 * xot2;
  const dixUpw = ctx.DI2 * xot2 - ctx.DI1 * xot1;
  const cfxUpw = ctx.CF2 * xot2 - ctx.CF1 * xot1;

  const btmp2 = 2.0 * hca / hsa + 1.0 - ha2 - hwa2;
  const rezh = hlog + btmp2 * ulog + xlog * (0.5 * cfx2 - dix);
  const zCfx2 = xlog * 0.5;
  const zDix = -xlog;
  const zHca = 2.0 * ulog / hsa;
  const zHa2 = -ulog;
  const zHwa2 = -ulog;
  const zXl2 = ddlog * (0.5 * cfx2 - dix);
  const zUl2 = ddlog * btmp2;
  const zHl = ddlog;

  const zUpw2 = zCfx2 * cfxUpw + zDix * dixUpw;

  const zHs1 = -hca * ulog / hsa ** 2 - zHl / ctx.HS1;
  const zHs2 = -hca * ulog / hsa ** 2 + zHl / ctx.HS2;

  const zCf1b = (1.0 - upw) * zCfx2 * xot1;
  const zCf2b = upw * zCfx2 * xot2;
  const zDi1 = (1.0 - upw) * zDix * xot1;
  const zDi2 = upw * zDix * xot2;

  let zT1b = (1.0 - upw) * (zCfx2 * ctx.CF1 + zDix * ctx.DI1) * (-xot1 / ctx.T1);
  let zT2b = upw * (zCfx2 * ctx.CF2 + zDix * ctx.DI2) * (-xot2 / ctx.T2);
  let zX1b = (1.0 - upw) * (zCfx2 * ctx.CF1 + zDix * ctx.DI1) / ctx.T1 - zXl2 / ctx.X1;
  let zX2b = upw * (zCfx2 * ctx.CF2 + zDix * ctx.DI2) / ctx.T2 + zXl2 / ctx.X2;
  const zU1b = -zUl2 / ctx.U1;
  const zU2b = zUl2 / ctx.U2;

  zT1b += zHwa2 * 0.5 * (-ctx.DW1 / ctx.T1 ** 2);
  zT2b += zHwa2 * 0.5 * (-ctx.DW2 / ctx.T2 ** 2);

  ctx.VS1[3][1] = zDi1 * ctx.DI1_S1;
  ctx.VS1[3][2] = zHs1 * ctx.HS1_T1 + zCf1b * ctx.CF1_T1 + zDi1 * ctx.DI1_T1 + zT1b;
  ctx.VS1[3][3] = zHs1 * ctx.HS1_D1 + zCf1b * ctx.CF1_D1 + zDi1 * ctx.DI1_D1;
  ctx.VS1[3][4] = zHs1 * ctx.HS1_U1 + zCf1b * ctx.CF1_U1 + zDi1 * ctx.DI1_U1 + zU1b;
  ctx.VS1[3][5] = zX1b;
  ctx.VS2[3][1] = zDi2 * ctx.DI2_S2;
  ctx.VS2[3][2] = zHs2 * ctx.HS2_T2 + zCf2b * ctx.CF2_T2 + zDi2 * ctx.DI2_T2 + zT2b;
  ctx.VS2[3][3] = zHs2 * ctx.HS2_D2 + zCf2b * ctx.CF2_D2 + zDi2 * ctx.DI2_D2;
  ctx.VS2[3][4] = zHs2 * ctx.HS2_U2 + zCf2b * ctx.CF2_U2 + zDi2 * ctx.DI2_U2 + zU2b;
  ctx.VS2[3][5] = zX2b;
  ctx.VSM[3] = zHs1 * ctx.HS1_MS + zCf1b * ctx.CF1_MS + zDi1 * ctx.DI1_MS
    + zHs2 * ctx.HS2_MS + zCf2b * ctx.CF2_MS + zDi2 * ctx.DI2_MS;
  ctx.VSR[3] = zHs1 * ctx.HS1_RE + zCf1b * ctx.CF1_RE + zDi1 * ctx.DI1_RE
    + zHs2 * ctx.HS2_RE + zCf2b * ctx.CF2_RE + zDi2 * ctx.DI2_RE;

  ctx.VS1[3][2] += 0.5 * (zHca * ctx.HC1_T1 + zHa2 * ctx.H1_T1) + zUpw2 * upwT1;
  ctx.VS1[3][3] += 0.5 * (zHca * ctx.HC1_D1 + zHa2 * ctx.H1_D1) + zUpw2 * upwD1;
  ctx.VS1[3][4] += 0.5 * (zHca * ctx.HC1_U1) + zUpw2 * upwU1;
  ctx.VS2[3][2] += 0.5 * (zHca * ctx.HC2_T2 + zHa2 * ctx.H2_T2) + zUpw2 * upwT2;
  ctx.VS2[3][3] += 0.5 * (zHca * ctx.HC2_D2 + zHa2 * ctx.H2_D2) + zUpw2 * upwD2;
  ctx.VS2[3][4] += 0.5 * (zHca * ctx.HC2_U2) + zUpw2 * upwU2;

  ctx.VSM[3] += 0.5 * (zHca * ctx.HC1_MS) + zUpw2 * upwMs
    + 0.5 * (zHca * ctx.HC2_MS);

  ctx.VSX[3] = 0.0;
  ctx.VSREZ[3] = -rezh;
}

// Laminar damping function for transition criteria.
function dampl(hk, th, rt) {
  const dgr = 0.08;
  const hmi = 1.0 / (hk - 1.0);
  const hmiHk = -(hmi ** 2);
  const aa = 2.492 * (hmi ** 0.43);
  const aaHk = (aa / hmi) * 0.43 * hmiHk;
  const bb = Math.tanh(14.0 * hmi - 9.24);
  const bbHk = (1.0 - bb * bb) * 14.0 * hmiHk;
  const grcrit = aa + 0.7 * (bb + 1.0);
  const grcHk = aaHk + 0.7 * bbHk;
  const gr = Math.log10(rt);
  const grRt = 1.0 / (2.3025851 * rt);

  if (gr < grcrit - dgr) {
    return { ax: 0.0, axHk: 0.0, axTh: 0.0, axRt: 0.0 };
  }

  const rnorm = (gr - (grcrit - dgr)) / (2.0 * dgr);
  const rnHk = -grcHk / (2.0 * dgr);
  const rnRt = grRt / (2.0 * dgr);

  let rfac;
  let rfacHk;
  let rfacRt;
  if (rnorm >= 1.0) {
    rfac = 1.0;
    rfacHk = 0.0;
    rfacRt = 0.0;
  } else {
    rfac = 3.0 * rnorm ** 2 - 2.0 * rnorm ** 3;
    const rfacRn = 6.0 * rnorm - 6.0 * rnorm ** 2;
    rfacHk = rfacRn * rnHk;
    rfacRt = rfacRn * rnRt;
  }

  const arg = 3.87 * hmi - 2.52;
  const argHk = 3.87 * hmiHk;
  const ex = Math.exp(-(arg ** 2));
  const exHk = ex * (-2.0 * arg * argHk);
  const dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
  const dadrHk = 0.028 - 0.0345 * exHk;
  const af = -0.05 + 2.7 * hmi - 5.5 * hmi ** 2 + 3.0 * hmi ** 3;
  const afHmi = 2.7 - 11.0 * hmi + 9.0 * hmi ** 2;
  const afHk = afHmi * hmiHk;
  const ax = (af * dadr / th) * rfac;
  const axHk = (afHk * dadr / th + af * dadrHk / th) * rfac + (af * dadr / th) * rfacHk;
  const axTh = -ax / th;
  const axRt = (af * dadr / th) * rfacRt;
  return { ax, axHk, axTh, axRt };
}

// Modified damping for separated flow regions.
function dampl2(hk, th, rt) {
  const dgr = 0.08;
  const hk1 = 3.5;
  const hk2 = 4.0;
  const hmi = 1.0 / (hk - 1.0);
  const hmiHk = -(hmi ** 2);
  const aa = 2.492 * (hmi ** 0.43);
  const aaHk = (aa / hmi) * 0.43 * hmiHk;
  const bb = Math.tanh(14.0 * hmi - 9.24);
  const bbHk = (1.0 - bb * bb) * 14.0 * hmiHk;
  const grc = aa + 0.7 * (bb + 1.0);
  const grcHk = aaHk + 0.7 * bbHk;
  const gr = Math.log10(rt);
  const grRt = 1.0 / (2.3025851 * rt);

  if (gr < grc - dgr) {
    return { ax: 0.0, axHk: 0.0, axTh: 0.0, axRt: 0.0 };
  }

  const rnorm = (gr - (grc - dgr)) / (2.0 * dgr);
  const rnHk = -grcHk / (2.0 * dgr);
  const rnRt = grRt / (2.0 * dgr);

  let rfac;
  let rfacHk;
  let rfacRt;
  if (rnorm >= 1.0) {
    rfac = 1.0;
    rfacHk = 0.0;
    rfacRt = 0.0;
  } else {
    rfac = 3.0 * rnorm ** 2 - 2.0 * rnorm ** 3;
    const rfacRn = 6.0 * rnorm - 6.0 * rnorm ** 2;
    rfacHk = rfacRn * rnHk;
    rfacRt = rfacRn * rnRt;
  }

  const arg = 3.87 * hmi - 2.52;
  const argHk = 3.87 * hmiHk;
  const ex = Math.exp(-(arg ** 2));
  const exHk = ex * (-2.0 * arg * argHk);
  const dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
  const dadrHk = 0.028 - 0.0345 * exHk;
  const brg = -20.0 * hmi;
  const af = -0.05 + 2.7 * hmi - 5.5 * hmi ** 2 + 3.0 * hmi ** 3 + 0.1 * Math.exp(brg);
  const afHmi = 2.7 - 11.0 * hmi + 9.0 * hmi ** 2 - 2.0 * Math.exp(brg);
  const afHk = afHmi * hmiHk;
  let ax = (af * dadr / th) * rfac;
  let axHk = (afHk * dadr / th + af * dadrHk / th) * rfac + (af * dadr / th) * rfacHk;
  let axTh = -ax / th;
  let axRt = (af * dadr / th) * rfacRt;

  if (hk < hk1) {
    return { ax, axHk, axTh, axRt };
  }

  const hnorm = (hk - hk1) / (hk2 - hk1);
  const hnHk = 1.0 / (hk2 - hk1);

  let hfac;
  let hfHk;
  if (hnorm >= 1.0) {
    hfac = 1.0;
    hfHk = 0.0;
  } else {
    hfac = 3.0 * hnorm ** 2 - 2.0 * hnorm ** 3;
    hfHk = (6.0 * hnorm - 6.0 * hnorm ** 2) * hnHk;
  }

  const ax1 = ax;
  const ax1Hk = axHk;
  const ax1Th = axTh;
  const ax1Rt = axRt;

  const gr0 = 0.30 + 0.35 * Math.exp(-0.15 * (hk - 5.0));
  const gr0Hk = -0.35 * Math.exp(-0.15 * (hk - 5.0)) * 0.15;

  const tnr = Math.tanh(1.2 * (gr - gr0));
  const tnrRt = (1.0 - tnr ** 2) * 1.2 * grRt;
  const tnrHk = -(1.0 - tnr ** 2) * 1.2 * gr0Hk;

  let ax2 = (0.086 * tnr - 0.25 / (hk - 1.0) ** 1.5) / th;
  let ax2Hk = (0.086 * tnrHk + 1.5 * 0.25 / (hk - 1.0) ** 2.5) / th;
  let ax2Rt = (0.086 * tnrRt) / th;
  let ax2Th = -ax2 / th;

  if (ax2 < 0.0) {
    ax2 = 0.0;
    ax2Hk = 0.0;
    ax2Rt = 0.0;
    ax2Th = 0.0;
  }

  ax = hfac * ax2 + (1.0 - hfac) * ax1;
  axHk = hfac * ax2Hk + (1.0 - hfac) * ax1Hk + hfHk * (ax2 - ax1);
  axRt = hfac * ax2Rt + (1.0 - hfac) * ax1Rt;
  axTh = hfac * ax2Th + (1.0 - hfac) * ax1Th;

  return { ax, axHk, axTh, axRt };
}

// Kinematic shape factor relation for compressibility correction.
function hkin(h, msq) {
  const hk = (h - 0.29 * msq) / (1.0 + 0.113 * msq);
  const hkH = 1.0 / (1.0 + 0.113 * msq);
  const hkMsq = (-0.29 - 0.113 * hk) / (1.0 + 0.113 * msq);
  return { hk, hkH, hkMsq };
}

function dil(hk, rt) {
  let di;
  let diHk;
  if (hk < 4.0) {
    di = (0.00205 * (4.0 - hk) ** 5.5 + 0.207) / rt;
    diHk = (-0.00205 * 5.5 * (4.0 - hk) ** 4.5) / rt;
  } else {
    const hkb = hk - 4.0;
    const den = 1.0 + 0.02 * hkb ** 2;
    di = (-0.0016 * hkb ** 2 / den + 0.207) / rt;
    diHk = (-0.0016 * 2.0 * hkb * (1.0 / den - 0.02 * hkb ** 2 / den ** 2)) / rt;
  }
  const diRt = -di / rt;
  return { di, diHk, diRt };
}

function dilw(hk, rt) {
  const msq = 0.0;
  const hslRes = hsl(hk, rt, msq);
  const rcd = 1.10 * (1.0 - 1.0 / hk) ** 2 / hk;
  const rcdHk = -1.10 * (1.0 - 1.0 / hk) * 2.0 / hk ** 3 - rcd / hk;
  const di = 2.0 * rcd / (hslRes.hs * rt);
  const diHk = 2.0 * rcdHk / (hslRes.hs * rt) - (di / hslRes.hs) * hslRes.hsHk;
  const diRt = -di / rt - (di / hslRes.hs) * hslRes.hsRt;
  return { di, diHk, diRt };
}

function hsl(hk, rt, msq) {
  let hs;
  let hsHk;
  if (hk < 4.35) {
    const tmp = hk - 4.35;
    hs = 0.0111 * tmp ** 2 / (hk + 1.0)
      - 0.0278 * tmp ** 3 / (hk + 1.0) + 1.528
      - 0.0002 * (tmp * hk) ** 2;
    hsHk = 0.0111 * (2.0 * tmp - tmp ** 2 / (hk + 1.0)) / (hk + 1.0)
      - 0.0278 * (3.0 * tmp ** 2 - tmp ** 3 / (hk + 1.0)) / (hk + 1.0)
      - 0.0002 * 2.0 * tmp * hk * (tmp + hk);
  } else {
    const hs2 = 0.015;
    hs = hs2 * (hk - 4.35) ** 2 / hk + 1.528;
    hsHk = hs2 * 2.0 * (hk - 4.35) / hk - hs2 * (hk - 4.35) ** 2 / hk ** 2;
  }
  return { hs, hsHk, hsRt: 0.0, hsMsq: 0.0 };
}

function cfl(hk, rt, msq) {
  let cf;
  let cfHk;
  if (hk < 5.5) {
    const tmp = (5.5 - hk) ** 3 / (hk + 1.0);
    cf = (0.0727 * tmp - 0.07) / rt;
    cfHk = (-0.0727 * tmp * 3.0 / (5.5 - hk) - 0.0727 * tmp / (hk + 1.0)) / rt;
  } else {
    const tmp = 1.0 - 1.0 / (hk - 4.5);
    cf = (0.015 * tmp ** 2 - 0.07) / rt;
    cfHk = (0.015 * tmp * 2.0 / (hk - 4.5) ** 2) / rt;
  }
  const cfRt = -cf / rt;
  const cfMsq = 0.0;
  return { cf, cfHk, cfRt, cfMsq };
}

function dit(hs, us, cf, st) {
  const di = (0.5 * cf * us + st * st * (1.0 - us)) * 2.0 / hs;
  const diHs = -(0.5 * cf * us + st * st * (1.0 - us)) * 2.0 / hs ** 2;
  const diUs = (0.5 * cf - st * st) * 2.0 / hs;
  const diCf = (0.5 * us) * 2.0 / hs;
  const diSt = (2.0 * st * (1.0 - us)) * 2.0 / hs;
  return { di, diHs, diUs, diCf, diSt };
}

function hst(hk, rt, msq) {
  const hsmin = 1.5;
  const dhsinf = 0.015;
  let ho;
  let hoRt;
  if (rt > 400.0) {
    ho = 3.0 + 400.0 / rt;
    hoRt = -400.0 / rt ** 2;
  } else {
    ho = 4.0;
    hoRt = 0.0;
  }

  let rtz;
  let rtzRt;
  if (rt > 200.0) {
    rtz = rt;
    rtzRt = 1.0;
  } else {
    rtz = 200.0;
    rtzRt = 0.0;
  }

  let hs;
  let hsHk;
  let hsRt;
  if (hk < ho) {
    const hr = (ho - hk) / (ho - 1.0);
    const hrHk = -1.0 / (ho - 1.0);
    const hrRt = (1.0 - hr) / (ho - 1.0) * hoRt;
    hs = (2.0 - hsmin - 4.0 / rtz) * hr ** 2 * 1.5 / (hk + 0.5) + hsmin + 4.0 / rtz;
    hsHk = -(2.0 - hsmin - 4.0 / rtz) * hr ** 2 * 1.5 / (hk + 0.5) ** 2
      + (2.0 - hsmin - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hrHk;
    hsRt = (2.0 - hsmin - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hrRt
      + (hr ** 2 * 1.5 / (hk + 0.5) - 1.0) * 4.0 / rtz ** 2 * rtzRt;
  } else {
    const grt = Math.log(rtz);
    const hdif = hk - ho;
    const rtmp = hk - ho + 4.0 / grt;
    const htmp = 0.007 * grt / rtmp ** 2 + dhsinf / hk;
    const htmpHk = -0.014 * grt / rtmp ** 3 - dhsinf / hk ** 2;
    const htmpRt = -0.014 * grt / rtmp ** 3 * (-hoRt - 4.0 / grt ** 2 / rtz * rtzRt)
      + 0.007 / rtmp ** 2 / rtz * rtzRt;
    hs = hdif ** 2 * htmp + hsmin + 4.0 / rtz;
    hsHk = hdif * 2.0 * htmp + hdif ** 2 * htmpHk;
    hsRt = hdif ** 2 * htmpRt - 4.0 / rtz ** 2 * rtzRt + hdif * 2.0 * htmp * (-hoRt);
  }

  const fm = 1.0 + 0.014 * msq;
  const hsMsq = 0.028 / fm - 0.014 * hs / fm;
  hs = (hs + 0.028 * msq) / fm;
  hsHk /= fm;
  hsRt /= fm;
  return { hs, hsHk, hsRt, hsMsq };
}

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
  const cfMsq = gex * cfo / (fc * grt) * (-0.25 * gm1 / fc ** 2) - 0.25 * gm1 * cf / fc ** 2;
  return { cf, cfHk, cfRt, cfMsq };
}

function hct(hk, msq) {
  const hc = msq * (0.064 / (hk - 0.8) + 0.251);
  const hcHk = msq * (-0.064 / (hk - 0.8) ** 2);
  const hcMsq = 0.064 / (hk - 0.8) + 0.251;
  return { hc, hcHk, hcMsq };
}

export {
  NCOM,
  ensureCtx,
  syncComToVars,
  syncVarsToCom,
  copyCom,
  trchek,
  trchek2,
  axset,
  blsys,
  tesys,
  blprv,
  blkin,
  blvar,
  blmid,
  trdif,
  bldif,
  dampl,
  dampl2,
  hkin,
  dil,
  dilw,
  hsl,
  cfl,
  dit,
  hst,
  cft,
  hct,
};
