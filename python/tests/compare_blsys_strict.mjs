import fs from 'fs';

import { ensureCtx, syncVarsToCom, blsys } from '../../js/xblsys.js';

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const { base } = payload;

function cloneValue(val) {
  if (!Array.isArray(val)) return val;
  if (val.length === 0) return [];
  if (Array.isArray(val[0])) {
    return val.map((row) => row.slice());
  }
  return val.slice();
}

const ctx = {};
ensureCtx(ctx);
for (const [key, val] of Object.entries(base)) {
  ctx[key] = cloneValue(val);
}
ensureCtx(ctx);
syncVarsToCom(ctx, 1);
syncVarsToCom(ctx, 2);

function pickArray(name, rows, cols) {
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

function pickVector(name, count) {
  const out = [];
  for (let i = 1; i <= count; i += 1) {
    out.push(ctx[name][i]);
  }
  return out;
}

blsys(ctx);

const results = {
  VS1: pickArray('VS1', 4, 5),
  VS2: pickArray('VS2', 4, 5),
  VSREZ: pickVector('VSREZ', 4),
  VSM: pickVector('VSM', 4),
  VSR: pickVector('VSR', 4),
  VSX: pickVector('VSX', 4),
};

process.stdout.write(JSON.stringify({ results }));
