import fs from 'fs';

import { lefind } from '../../js/xgeom.js';
import { segspl, scalc } from '../../js/spline.js';

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const { x, y, s } = payload;
const n = x.length;

const sArr = new Float64Array(s);
const xp = new Float64Array(n);
const yp = new Float64Array(n);
segspl(x, xp, sArr, n);
segspl(y, yp, sArr, n);

const sle = lefind(x, xp, y, yp, sArr, n);

process.stdout.write(JSON.stringify({ results: { sle } }));
