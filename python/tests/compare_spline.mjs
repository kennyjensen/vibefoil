import fs from 'fs';

import { scalc, segspl, seval, deval, d2val, sinvrt } from '../../js/spline.js';

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const { x, y, queries, sinvrtCases } = payload;
const n = x.length;

const s = new Float64Array(n);
scalc(x, y, s, n);

const xs = new Float64Array(n);
segspl(x, xs, s, n);

const evals = queries.map((sq) => ({
  seval: seval(sq, x, xs, s, n),
  deval: deval(sq, x, xs, s, n),
  d2val: d2val(sq, x, xs, s, n),
}));

const sinv = sinvrtCases.map((item) => sinvrt(item.si, item.xi, x, xs, s, n));

const results = {
  scalc: Array.from(s),
  segspl: Array.from(xs),
  evals,
  sinvrt: sinv,
};

process.stdout.write(JSON.stringify({ results }));
