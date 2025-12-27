import fs from 'fs';

import { gauss } from '../../js/xsolve.js';
import { blpini, mrchue, xifset } from '../../js/xbl.js';
import {
  blprv,
  blkin,
  blsys,
  tesys,
  trchek,
  hkin,
  blvar,
  blmid,
  syncComToVars,
  ensureCtx,
} from '../../js/xblsys.js';

const noop = () => {};
console.log = noop;
console.warn = noop;

function gauss1(n, z, r) {
  const z0 = new Array(n);
  for (let i = 0; i < n; i += 1) {
    z0[i] = new Float64Array(n);
    for (let j = 0; j < n; j += 1) {
      z0[i][j] = z[i + 1][j + 1];
    }
  }

  const r0 = new Float64Array(n);
  for (let i = 0; i < n; i += 1) {
    r0[i] = r[i + 1];
  }

  gauss(n, z0, r0);

  for (let i = 0; i < n; i += 1) {
    r[i + 1] = r0[i];
    for (let j = 0; j < n; j += 1) {
      z[i + 1][j + 1] = z0[i][j];
    }
  }
}

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const { ctx: ctxIn, bl: blIn } = payload;

function toFloatArray(arr) {
  if (!arr) return null;
  const out = new Float64Array(arr.length);
  for (let i = 0; i < arr.length; i += 1) out[i] = arr[i];
  return out;
}

function toMatrix(arr) {
  if (!arr) return null;
  return arr.map((row) => toFloatArray(row));
}

const ctx = {
  NBL: Int32Array.from(ctxIn.NBL),
  IBLTE: Int32Array.from(ctxIn.IBLTE),
  ITRAN: Int32Array.from(ctxIn.ITRAN),
  XSSITR: Float64Array.from(ctxIn.XSSITR),
  XSTRIP: Float64Array.from(ctxIn.XSTRIP),
  ACRIT: Float64Array.from(ctxIn.ACRIT),
  WGAP: Float64Array.from(ctxIn.WGAP),
  TFORCE: ctxIn.TFORCE.slice(),
  XSSI: toMatrix(ctxIn.XSSI),
  UEDG: toMatrix(ctxIn.UEDG),
  THET: toMatrix(ctxIn.THET),
  DSTR: toMatrix(ctxIn.DSTR),
  CTAU: toMatrix(ctxIn.CTAU),
  MASS: toMatrix(ctxIn.MASS),
  TAU: toMatrix(ctxIn.TAU),
  DIS: toMatrix(ctxIn.DIS),
  CTQ: toMatrix(ctxIn.CTQ),
  DELT: toMatrix(ctxIn.DELT),
  TSTR: toMatrix(ctxIn.TSTR),
  ANTE: ctxIn.ANTE,
  DWTE: blIn.DWTE,
  GAMBL: blIn.GAMBL,
  GM1BL: blIn.GM1BL,
  QINFBL: blIn.QINFBL,
  TKBL: blIn.TKBL,
  TKBL_MS: blIn.TKBL_MS,
  RSTBL: blIn.RSTBL,
  RSTBL_MS: blIn.RSTBL_MS,
  HSTINV: blIn.HSTINV,
  HSTINV_MS: blIn.HSTINV_MS,
  HVRAT: blIn.HVRAT,
  REYBL: blIn.REYBL,
  REYBL_RE: blIn.REYBL_RE,
  REYBL_MS: blIn.REYBL_MS,
  IDAMPV: blIn.IDAMPV,
  blprv,
  blkin,
  blsys,
  tesys,
  trchek,
  hkin,
  blvar,
  blmid,
  syncComToVars,
};

ensureCtx(ctx);
blpini(ctx);
mrchue(ctx);

function firstIterationRlx(is) {
  ctx.AMCRIT = ctx.ACRIT[is];
  xifset(ctx, is);

  const ibl = 2;
  const simi = ibl === 2;
  const wake = ibl > ctx.IBLTE[is];
  ctx.SIMI = simi;
  ctx.WAKE = wake;

  const itrold = ctx.ITRAN[is];
  ctx.TRAN = false;
  ctx.TURB = false;
  ctx.ITRAN[is] = ctx.IBLTE[is];

  let xsi = ctx.XSSI[ibl][is];
  let uei = ctx.UEDG[ibl][is];
  let thi = ctx.THET[ibl][is];
  let dsi = ctx.DSTR[ibl][is];

  let ami = 0.0;
  let cti = 0.03;
  if (ibl < itrold) {
    ami = ctx.CTAU[ibl][is];
  } else {
    cti = ctx.CTAU[ibl][is];
    if (cti <= 0.0) cti = 0.03;
  }

  let dswaki = 0.0;
  if (wake) {
    const iw = ibl - ctx.IBLTE[is];
    dswaki = ctx.WGAP[iw];
  }

  if (ibl <= ctx.IBLTE[is]) {
    dsi = Math.max(dsi - dswaki, 1.02000 * thi) + dswaki;
  }
  if (ibl > ctx.IBLTE[is]) {
    dsi = Math.max(dsi - dswaki, 1.00005 * thi) + dswaki;
  }

  blprv(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
  blkin(ctx);

  if (!simi && !ctx.TURB) {
    trchek(ctx);
    ami = ctx.AMPL2;
    if (ctx.TRAN) ctx.ITRAN[is] = ibl;
    if (!ctx.TRAN) ctx.ITRAN[is] = ibl + 2;
  }

  if (ibl === ctx.IBLTE[is] + 1) {
    const tte = ctx.THET[ctx.IBLTE[1]][1] + ctx.THET[ctx.IBLTE[2]][2];
    const dte = ctx.DSTR[ctx.IBLTE[1]][1] + ctx.DSTR[ctx.IBLTE[2]][2] + ctx.ANTE;
    const cte = (ctx.CTAU[ctx.IBLTE[1]][1] * ctx.THET[ctx.IBLTE[1]][1]
      + ctx.CTAU[ctx.IBLTE[2]][2] * ctx.THET[ctx.IBLTE[2]][2]) / tte;
    tesys(cte, tte, dte, ctx);
  } else {
    blsys(ctx);
  }

  let ueref = ctx.U2;
  let hkref = ctx.HK2;

  if (ibl < ctx.ITRAN[is] && ibl >= itrold) {
    const uem = ctx.UEDG[ibl - 1][is];
    const dsm = ctx.DSTR[ibl - 1][is];
    const thm = ctx.THET[ibl - 1][is];
    const msq = uem * uem * ctx.HSTINV / (ctx.GM1BL * (1.0 - 0.5 * uem * uem * ctx.HSTINV));
    hkref = hkin(dsm / thm, msq).hk;
  }

  if (ibl < itrold) {
    if (ctx.TRAN) ctx.CTAU[ibl][is] = 0.03;
    if (ctx.TURB) ctx.CTAU[ibl][is] = ctx.CTAU[ibl - 1][is];
    if (ctx.TRAN || ctx.TURB) {
      cti = ctx.CTAU[ibl][is];
      ctx.S2 = cti;
    }
  }

  if (simi || ibl === ctx.IBLTE[is] + 1) {
    ctx.VS2[4][1] = 0.0;
    ctx.VS2[4][2] = 0.0;
    ctx.VS2[4][3] = 0.0;
    ctx.VS2[4][4] = ctx.U2_UEI;
    ctx.VSREZ[4] = ueref - ctx.U2;
  } else {
    const vtmp = Array.from({ length: 5 }, () => new Float64Array(5));
    const vztmp = new Float64Array(5);
    for (let k = 1; k <= 4; k += 1) {
      vztmp[k] = ctx.VSREZ[k];
      for (let l = 1; l <= 4; l += 1) {
        vtmp[k][l] = ctx.VS2[k][l];
      }
    }
    vtmp[4][1] = 0.0;
    vtmp[4][2] = ctx.HK2_T2;
    vtmp[4][3] = ctx.HK2_D2;
    vtmp[4][4] = ctx.HK2_U2 * ctx.U2_UEI;
    vztmp[4] = 1.0;
    gauss1(4, vtmp, vztmp);
    const sennew = 1000.0 * vztmp[4] * hkref / ueref;
    const sens = sennew;
    ctx.VS2[4][1] = 0.0;
    ctx.VS2[4][2] = ctx.HK2_T2 * hkref;
    ctx.VS2[4][3] = ctx.HK2_D2 * hkref;
    ctx.VS2[4][4] = (ctx.HK2_U2 * hkref + sens / ueref) * ctx.U2_UEI;
    ctx.VSREZ[4] = -(hkref ** 2) * (ctx.HK2 / hkref - 1.0)
      - sens * (ctx.U2 / ueref - 1.0);
  }

  gauss1(4, ctx.VS2, ctx.VSREZ);

  let dmax = Math.max(Math.abs(ctx.VSREZ[2] / thi),
    Math.abs(ctx.VSREZ[3] / dsi),
    Math.abs(ctx.VSREZ[4] / uei));
  if (ibl >= ctx.ITRAN[is]) {
    dmax = Math.max(dmax, Math.abs(ctx.VSREZ[1] / (10.0 * cti)));
  }

  let rlx = 1.0;
  if (dmax > 0.3) rlx = 0.3 / dmax;

  return {
    ibl,
    itrold,
    simi,
    wake,
    xsi,
    uei,
    thi,
    dsi,
    cti,
    ami,
    ueref,
    hkref,
    dmax,
    rlx,
    VSREZ: [ctx.VSREZ[1], ctx.VSREZ[2], ctx.VSREZ[3], ctx.VSREZ[4]],
  };
}

const results = {
  1: firstIterationRlx(1),
  2: firstIterationRlx(2),
};

process.stdout.write(JSON.stringify({ results }));
