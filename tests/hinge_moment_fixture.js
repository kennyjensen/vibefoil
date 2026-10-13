import '../js/naca.js';
import { flap } from '../js/xgdes.js';
import { scalc, segspl } from '../js/spline.js';

export function makeHingeFixture({ foil, sharp, x, yrel, deflection, viscous, symmetric }) {
  const nside = 40;
  const xb = new Float64Array(2 * nside), yb = new Float64Array(2 * nside);
  const generate = foil < 10000 ? globalThis.Naca.naca4 : globalThis.Naca.naca5;
  const { nb } = generate(foil, new Float64Array(nside), new Float64Array(nside),
    new Float64Array(nside), nside, xb, yb, sharp);
  const geometry = flap(xb.slice(0, nb), yb.slice(0, nb), { x, yrel, deflection });
  const n = geometry.nb;
  const panel = { N: n, X: geometry.xb, Y: geometry.yb,
    XP: new Float64Array(n), YP: new Float64Array(n), S: new Float64Array(n) };
  scalc(panel.X, panel.Y, panel.S, n);
  segspl(panel.X, panel.XP, panel.S, n);
  segspl(panel.Y, panel.YP, panel.S, n);
  const cp = new Float64Array(n + 1);
  for (let i = 0; i < n; i++) {
    cp[i + 1] = symmetric ? 1.0 - panel.X[i]
      : 0.2 - 0.7 * panel.X[i] + 0.4 * panel.Y[i]
        + (viscous ? 0.23 : 0.3) * Math.sin(2 * Math.PI * panel.S[i] / panel.S[n - 1]);
  }
  return { panel, cp, hinge: geometry.hinge };
}
