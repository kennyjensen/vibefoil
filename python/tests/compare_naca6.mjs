import fs from 'fs';

import '../../js/naca.js';

const input = fs.readFileSync(0, 'utf8');
const payload = JSON.parse(input);

const results = payload.cases.map((caseDef) => {
  const { profile, toc, camber, cl, a, nside } = caseDef;
  const xx = new Float64Array(nside);
  const yt = new Float64Array(nside);
  const yc = new Float64Array(nside);
  const xb = new Float64Array(2 * nside);
  const yb = new Float64Array(2 * nside);

  const res = globalThis.Naca.naca6(
    { profile, toc, camber, cl, a },
    xx,
    yt,
    yc,
    nside,
    xb,
    yb
  );

  const xUpper = new Float64Array(nside);
  const yUpper = new Float64Array(nside);
  for (let i = 0; i < nside; i += 1) {
    xUpper[i] = xb[nside - 1 - i];
    yUpper[i] = yb[nside - 1 - i];
  }
  const xLower = new Float64Array(nside);
  const yLower = new Float64Array(nside);
  xLower[0] = xUpper[0];
  yLower[0] = yUpper[0];
  for (let i = 1; i < nside; i += 1) {
    xLower[i] = xb[nside - 1 + i];
    yLower[i] = yb[nside - 1 + i];
  }

  return {
    profile,
    toc,
    camber,
    cl,
    a,
    nside,
    ok: res.ok,
    x: Array.from(xUpper),
    yu: Array.from(yUpper),
    yl: Array.from(yLower),
  };
});

process.stdout.write(JSON.stringify({ results }));
