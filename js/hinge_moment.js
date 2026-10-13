// Ported from XFOIL's MHINGE (Mark Drela), GPL-2.0-or-later.
// Geometry is zero-based; Cp retains the worker's one-based indexing.
// Read-only postprocessing, with the same chord reference as CL/CM/CD.
import { seval, sinvrt } from './spline.js';

export function computeHingeMoment(panel, cp, hinge) {
  if (!panel || !cp || !Number.isFinite(hinge?.x) || !Number.isFinite(hinge?.y)) return NaN;
  const { N: n, X: x, Y: y, XP: xp, YP: yp, S: s } = panel;
  if (n < 2 || cp.length <= n) return NaN;
  const xf = hinge.x, yf = hinge.y;
  const tops = sinvrt(xf, xf, x, xp, s, n);
  const bots = sinvrt(s[n - 1] - xf, xf, x, xp, s, n);
  if (!Number.isFinite(tops) || !Number.isFinite(bots)
    || tops < s[0] || bots > s[n - 1] || tops >= bots) return NaN;
  const topx = seval(tops, x, xp, s, n);
  const topy = seval(tops, y, yp, s, n);
  const botx = seval(bots, x, xp, s, n);
  const boty = seval(bots, y, yp, s, n);

  let moment = 0.0;
  const add = (x1, y1, x2, y2, pressure) => {
    const dx = x2 - x1, dy = y2 - y1;
    const xmid = 0.5 * (x2 + x1) - xf;
    const ymid = 0.5 * (y2 + y1) - yf;
    moment += pressure * (xmid * dx + ymid * dy);
  };
  for (let i = 1; i < n; i++) {
    if (s[i - 1] >= tops && s[i] <= bots) continue;
    add(x[i - 1], y[i - 1], x[i], y[i], 0.5 * (cp[i + 1] + cp[i]));
  }

  let i = 1;
  while (i < n && s[i] <= tops) i++;
  if (i === n) return NaN;
  let fraction = s[i] !== s[i - 1] ? (tops - s[i - 1]) / (s[i] - s[i - 1]) : 0.0;
  const topPressure = cp[i + 1] * fraction + cp[i] * (1.0 - fraction);
  let pressure = 0.5 * (topPressure + cp[i]);
  add(x[i - 1], y[i - 1], topx, topy, pressure);
  add(topx, topy, xf, yf, pressure);

  i = n - 1;
  while (i > 0 && s[i - 1] >= bots) i--;
  if (i === 0) return NaN;
  fraction = s[i] !== s[i - 1] ? (bots - s[i - 1]) / (s[i] - s[i - 1]) : 0.0;
  const bottomPressure = cp[i + 1] * fraction + cp[i] * (1.0 - fraction);
  pressure = 0.5 * (bottomPressure + cp[i + 1]);
  add(botx, boty, x[i], y[i], pressure);
  add(xf, yf, botx, boty, pressure);

  // Close the trailing-edge base, including foils with finite TE thickness.
  add(x[n - 1], y[n - 1], x[0], y[0], 0.5 * (cp[1] + cp[n]));
  return Number.isFinite(moment) ? moment : NaN;
}
