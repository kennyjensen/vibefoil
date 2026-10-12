// Derived from XFOIL Fortran source (Mark Drela) and the PSILIN JavaScript port.
// This file remains under the GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original source.

import { psilin } from './xpanel.js';

// The streamline display needs only PSILIN's scalar result, without source or
// geometry linearization. Prepare the panel constants once for the whole grid.
// Keep PSILIN's expression and summation order, including the trailing edge.
// Rebuild this evaluator after changing geometry, circulation, or freestream.
export function createStreamfunctionEvaluator(ctx) {
  if (ctx.LIMAGE) {
    return (xi, yi) => psilin(ctx.N, xi, yi, 0.0, 0.0, false, false, ctx).psi;
  }
  const { N: n, X: x, Y: y, S: s, GAM: gam, PI: pi, QOPI: qopi, HOPI: hopi, QINF: qinf } = ctx;
  const cosa = Math.cos(ctx.ALFA);
  const sina = Math.sin(ctx.ALFA);
  const seps = (s[n - 1] - s[0]) * 1.0e-5;
  const scs = ctx.SHARP ? 1.0 : ctx.ANTE / ctx.DSTE;
  const sds = ctx.SHARP ? 0.0 : ctx.ASTE / ctx.DSTE;
  const panels = [];
  let previousEndpoint = -1;
  for (let jo = 0; jo < n; jo += 1) {
    const jp = jo === n - 1 ? 0 : jo + 1;
    const dx = x[jo] - x[jp];
    const dy = y[jo] - y[jp];
    if (jo === n - 1 && dx * dx + dy * dy < seps * seps) break;
    const dso = Math.sqrt(dx * dx + dy * dy);
    if (dso === 0.0) continue;
    const dsio = 1.0 / dso;
    panels.push({
      x1: x[jo], y1: y[jo], x2: x[jp], y2: y[jp],
      sx: (x[jp] - x[jo]) * dsio,
      sy: (y[jp] - y[jo]) * dsio,
      gsum: gam[jp] + gam[jo], gdif: gam[jp] - gam[jo],
      trailingEdge: jo === n - 1,
      apan: ctx.APANEL[jo],
      sharedStart: previousEndpoint === jo,
    });
    previousEndpoint = jp;
  }
  const sigte = 0.5 * scs * (gam[0] - gam[n - 1]);
  const gamte = -0.5 * sds * (gam[0] - gam[n - 1]);

  return function streamfunction(xi, yi) {
    let psi = 0.0;
    let previousG2 = 0.0;
    for (let j = 0; j < panels.length; j += 1) {
      const panel = panels[j];
      const rx1 = xi - panel.x1;
      const ry1 = yi - panel.y1;
      const rx2 = xi - panel.x2;
      const ry2 = yi - panel.y2;
      const x1 = panel.sx * rx1 + panel.sy * ry1;
      const x2 = panel.sx * rx2 + panel.sy * ry2;
      const yy = panel.sx * ry1 - panel.sy * rx1;
      const rs1 = rx1 * rx1 + ry1 * ry1;
      const rs2 = rx2 * rx2 + ry2 * ry2;
      const sgn = yy >= 0.0 ? 1.0 : -1.0;
      let g1 = 0.0;
      let g2 = 0.0;
      let t1 = 0.0;
      let t2 = 0.0;
      if (rs1 > 0.0) {
        // Adjacent panels use exactly the same endpoint logarithm. Their
        // local angles still need separate work; skipped panels break reuse.
        g1 = panel.sharedStart ? previousG2 : Math.log(rs1);
        t1 = Math.atan2(sgn * x1, sgn * yy) + (0.5 - 0.5 * sgn) * pi;
      }
      if (rs2 > 0.0) {
        g2 = Math.log(rs2);
        t2 = Math.atan2(sgn * x2, sgn * yy) + (0.5 - 0.5 * sgn) * pi;
      }
      previousG2 = g2;
      if (panel.trailingEdge) {
        const psig = 0.5 * yy * (g1 - g2) + x2 * (t2 - panel.apan) - x1 * (t1 - panel.apan);
        const pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
        psi += hopi * (psig * sigte + pgam * gamte);
      } else {
        const dxinv = 1.0 / (x1 - x2);
        const psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
        const psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;
        psi += qopi * (psis * panel.gsum + psid * panel.gdif);
      }
    }
    return psi + qinf * (cosa * yi - sina * xi);
  };
}
