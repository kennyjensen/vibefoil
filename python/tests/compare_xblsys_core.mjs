import fs from 'fs';

import {
  ensureCtx,
  syncVarsToCom,
  trchek,
  blvar,
  blmid,
  blsys,
  tesys,
} from '../../js/xblsys.js';

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const { base, config } = payload;

function cloneValue(val) {
  if (!Array.isArray(val)) return val;
  if (val.length === 0) return [];
  if (Array.isArray(val[0])) {
    return val.map((row) => row.slice());
  }
  return val.slice();
}

function makeCtx(state) {
  const ctx = {};
  ensureCtx(ctx);
  for (const [key, val] of Object.entries(state)) {
    ctx[key] = cloneValue(val);
  }
  ensureCtx(ctx);
  syncVarsToCom(ctx, 1);
  syncVarsToCom(ctx, 2);
  return ctx;
}

function pickFields(ctx, fields) {
  const out = {};
  for (const name of fields) {
    out[name] = ctx[name];
  }
  return out;
}

function pickArray(ctx, name, rows, cols) {
  const out = [];
  for (let i = 1; i <= rows; i += 1) {
    const row = [];
    for (let j = 1; j <= cols; j += 1) {
      row.push(ctx[name][i][j]);
    }
    out.push(row);
  }
  return out;
}

function pickVector(ctx, name, count) {
  const out = [];
  for (let i = 1; i <= count; i += 1) {
    out.push(ctx[name][i]);
  }
  return out;
}

const results = {};

if (config.trchek) {
  const ctx = makeCtx({ ...base, ...(config.trchek.override || {}) });
  trchek(ctx);
  results.trchek = pickFields(ctx, config.trchek.fields);
}

if (config.blvar) {
  results.blvar = {};
  for (const ityp of config.blvar.ityps) {
    const ctx = makeCtx(base);
    blvar(ityp, ctx);
    results.blvar[String(ityp)] = pickFields(ctx, config.blvar.fields);
  }
}

if (config.blmid) {
  results.blmid = {};
  for (const ityp of config.blmid.ityps) {
    const ctx = makeCtx(base);
    blmid(ityp, ctx);
    results.blmid[String(ityp)] = pickFields(ctx, config.blmid.fields);
  }
}

if (config.blsys) {
  const ctx = makeCtx(base);
  blsys(ctx);
  results.blsys = {
    VS1: pickArray(ctx, 'VS1', 4, 5),
    VS2: pickArray(ctx, 'VS2', 4, 5),
    VSREZ: pickVector(ctx, 'VSREZ', 4),
    VSM: pickVector(ctx, 'VSM', 4),
    VSR: pickVector(ctx, 'VSR', 4),
    VSX: pickVector(ctx, 'VSX', 4),
  };
}

if (config.tesys) {
  const ctx = makeCtx(base);
  tesys(config.tesys.cte, config.tesys.tte, config.tesys.dte, ctx);
  results.tesys = {
    VS1: pickArray(ctx, 'VS1', 4, 5),
    VS2: pickArray(ctx, 'VS2', 4, 5),
    VSREZ: pickVector(ctx, 'VSREZ', 4),
    VSM: pickVector(ctx, 'VSM', 4),
    VSR: pickVector(ctx, 'VSR', 4),
    VSX: pickVector(ctx, 'VSX', 4),
  };
}

process.stdout.write(JSON.stringify({ results }));
