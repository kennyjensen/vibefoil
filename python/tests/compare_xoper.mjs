import fs from 'fs';

import { specal } from '../../js/xoper.js';

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const { qinvu, alfa } = payload;
const n = qinvu.length;

const ctxPanel = {
  N: n,
  NW: 0,
  ALFA: alfa,
  QINVU: qinvu.map((pair) => new Float64Array(pair)),
  QINV: new Float64Array(n + 1),
  QINV_A: new Float64Array(n + 1),
};

const { qinv, qinvA } = specal(ctxPanel, alfa);

const results = {
  qinv: Array.from(qinv),
  qinvA: Array.from(qinvA),
};

process.stdout.write(JSON.stringify({ results }));
