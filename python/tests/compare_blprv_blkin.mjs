import fs from 'fs';

import { blprv, blkin, ensureCtx } from '../../js/xblsys.js';

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const { cases, fields } = payload;

const results = cases.map(({ constants, inputs }) => {
  const ctx = { ...constants };
  ensureCtx(ctx);
  blprv(inputs.xsi, inputs.ami, inputs.cti, inputs.thi, inputs.dsi, inputs.dswaki, inputs.uei, ctx);
  blkin(ctx);

  const out = {};
  for (const field of fields) {
    out[field] = ctx[field];
  }
  return out;
});

process.stdout.write(JSON.stringify({ results }));
