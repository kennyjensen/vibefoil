import fs from 'fs';

import { inside, getxyf, sss } from '../../js/xgdes.js';
import { segspl, scalc } from '../../js/spline.js';

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const { x, y, s, xf, yrel, insidePts, sssCases } = payload;
const n = x.length;

const sArr = new Float64Array(s);
const xp = new Float64Array(n);
const yp = new Float64Array(n);
segspl(x, xp, sArr, n);
segspl(y, yp, sArr, n);

const hinge = getxyf(x, xp, y, yp, sArr, n, xf, yrel);

const insideResults = insidePts.map((pt) => (inside(x, y, n, pt.x, pt.y) ? 1 : 0));

const sssResults = sssCases.map((item) => {
  const out = sss(item.ss, item.del, item.xbf, item.ybf, x, xp, y, yp, sArr, n, item.iside);
  return { s1: out.s1, s2: out.s2 };
});

const results = {
  getxyf: hinge,
  inside: insideResults,
  sss: sssResults,
};

process.stdout.write(JSON.stringify({ results }));
