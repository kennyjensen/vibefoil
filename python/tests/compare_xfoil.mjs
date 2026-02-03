import fs from 'fs';

import { tecalc, cpcalc, clcalc } from '../../js/xfoil.js';
import { scalc, segspl } from '../../js/spline.js';

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const { x, y, gam, gamA, q, alfa, minf, qinf, xref, yref } = payload;
const n = x.length;

const s = new Float64Array(n);
scalc(x, y, s, n);
const xp = new Float64Array(n);
const yp = new Float64Array(n);
segspl(x, xp, s, n);
segspl(y, yp, s, n);

const ctx = {
  N: n,
  X: x,
  Y: y,
  XP: xp,
  YP: yp,
  ANTE: 0.0,
  ASTE: 0.0,
  DSTE: 0.0,
  SHARP: false,
};

tecalc(ctx);
const cp = cpcalc(q, qinf, minf);
const clres = clcalc(n, x, y, gam, gamA, alfa, minf, qinf, xref, yref);

const results = {
  te: {
    ante: ctx.ANTE,
    aste: ctx.ASTE,
    dste: ctx.DSTE,
    sharp: ctx.SHARP ? 1.0 : 0.0,
  },
  cp: Array.from(cp),
  cl: clres,
};

process.stdout.write(JSON.stringify({ results }));
