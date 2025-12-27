import fs from 'fs';

import { blpini, mrchue, mrchdu } from '../../js/xbl.js';
import {
  blprv,
  blkin,
  blvar,
  blmid,
  blsys,
  tesys,
  trchek,
  hkin,
  syncComToVars,
  ensureCtx,
} from '../../js/xblsys.js';

const noop = () => {};
console.log = noop;
console.warn = noop;

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
  blvar,
  blmid,
  blsys,
  tesys,
  trchek,
  hkin,
  syncComToVars,
};

ensureCtx(ctx);
blpini(ctx);

mrchue(ctx);
mrchdu(ctx);

const maxNbl = Math.max(ctx.NBL[1], ctx.NBL[2]);
function extractSide(mat, is) {
  const out = new Array(maxNbl + 1).fill(0.0);
  for (let ibl = 1; ibl <= ctx.NBL[is]; ibl += 1) {
    out[ibl] = mat[ibl][is];
  }
  return out;
}

function extractSide1d(arr, is) {
  const out = new Array(maxNbl + 1).fill(0.0);
  const side = arr?.[is];
  if (!side) return out;
  for (let ibl = 1; ibl <= ctx.NBL[is]; ibl += 1) {
    out[ibl] = side[ibl];
  }
  return out;
}

const results = {
  NBL: Array.from(ctx.NBL),
  IBLTE: Array.from(ctx.IBLTE),
  ITRAN: Array.from(ctx.ITRAN),
  XSSITR: Array.from(ctx.XSSITR),
  TFORCE: ctx.TFORCE.slice(),
  THET: { 1: extractSide(ctx.THET, 1), 2: extractSide(ctx.THET, 2) },
  DSTR: { 1: extractSide(ctx.DSTR, 1), 2: extractSide(ctx.DSTR, 2) },
  CTAU: { 1: extractSide(ctx.CTAU, 1), 2: extractSide(ctx.CTAU, 2) },
  UEDG: { 1: extractSide(ctx.UEDG, 1), 2: extractSide(ctx.UEDG, 2) },
  MASS: { 1: extractSide(ctx.MASS, 1), 2: extractSide(ctx.MASS, 2) },
  TAU: { 1: extractSide(ctx.TAU, 1), 2: extractSide(ctx.TAU, 2) },
  DIS: { 1: extractSide(ctx.DIS, 1), 2: extractSide(ctx.DIS, 2) },
  CTQ: { 1: extractSide(ctx.CTQ, 1), 2: extractSide(ctx.CTQ, 2) },
  DELT: { 1: extractSide(ctx.DELT, 1), 2: extractSide(ctx.DELT, 2) },
  TSTR: { 1: extractSide(ctx.TSTR, 1), 2: extractSide(ctx.TSTR, 2) },
  HTARG: { 1: extractSide1d(ctx.HTARG, 1), 2: extractSide1d(ctx.HTARG, 2) },
};

process.stdout.write(JSON.stringify({ results }));
