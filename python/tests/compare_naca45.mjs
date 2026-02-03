import fs from 'fs';

import '../../js/naca.js';

const input = fs.readFileSync(0, 'utf8');
const payload = JSON.parse(input);

const results = payload.cases.map((caseDef) => {
  const { ides, nside } = caseDef;
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
  return {
    ides,
    nb,
    xb: Array.from(xb.slice(0, nb)),
    yb: Array.from(yb.slice(0, nb)),
  };
});

process.stdout.write(JSON.stringify({ results }));
