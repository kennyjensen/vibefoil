import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

import { pangen } from '../../js/xfoil.js';
import '../../js/naca.js';

const input = fs.readFileSync(0, 'utf8');
const payload = JSON.parse(input);

const results = payload.cases.map((caseDef) => {
  const { ides, nside, params } = caseDef;
  const xx = new Float64Array(nside);
  const yt = new Float64Array(nside);
  const yc = new Float64Array(nside);
  const xb = new Float64Array(2 * nside);
  const yb = new Float64Array(2 * nside);

  let res;
  if (ides <= 9999) {
    res = globalThis.Naca.naca4(ides, xx, yt, yc, nside, xb, yb);
  } else {
    res = globalThis.Naca.naca5(ides, xx, yt, yc, nside, xb, yb);
  }

  const nb = res.nb;
  const xb1 = Array(nb + 1).fill(0.0);
  const yb1 = Array(nb + 1).fill(0.0);
  for (let i = 0; i < nb; i += 1) {
    xb1[i + 1] = xb[i];
    yb1[i + 1] = yb[i];
  }

  const panel = pangen(xb.subarray(0, nb), yb.subarray(0, nb), nb, params);
  const x1 = Array(panel.n + 1).fill(0.0);
  const y1 = Array(panel.n + 1).fill(0.0);
  const s1 = Array(panel.n + 1).fill(0.0);
  for (let i = 0; i < panel.n; i += 1) {
    x1[i + 1] = panel.x[i];
    y1[i + 1] = panel.y[i];
    s1[i + 1] = panel.s[i];
  }

  return {
    ides,
    naca: { nb, xb: xb1, yb: yb1 },
    pangen: { n: panel.n, x: x1, y: y1, s: s1 },
  };
});

process.stdout.write(JSON.stringify({ results }));
