import assert from 'node:assert/strict';
import '../js/naca.js';
import { scalc, segspl } from '../js/spline.js';
import { ncalc, apcalc, psilin } from '../js/xpanel.js';
import { tecalc } from '../js/xfoil.js';
import { createStreamfunctionEvaluator } from '../js/streamfunction.js';

let comparisons = 0;
for (const ides of [12, 2412, 23012]) {
  for (const sharp of [false, true]) {
    for (const rotation of [0, 0.7, Math.PI / 2]) {
      const nside = 40;
      const x = new Float64Array(2 * nside);
      const y = new Float64Array(2 * nside);
      const generate = ides < 10000 ? globalThis.Naca.naca4 : globalThis.Naca.naca5;
      const { nb: n } = generate(ides, new Float64Array(nside), new Float64Array(nside),
        new Float64Array(nside), nside, x, y, sharp);
      for (let i = 0; i < n; i += 1) {
        const xi = x[i], yi = y[i];
        x[i] = xi * Math.cos(rotation) - yi * Math.sin(rotation) + 0.1;
        y[i] = xi * Math.sin(rotation) + yi * Math.cos(rotation) - 0.2;
      }
      // Include an internal zero-length panel as well as closed/sharp TEs.
      if (rotation === 0.7) {
        x[10] = x[9];
        y[10] = y[9];
      }
      const ctx = { N: n, X: x, Y: y, CHORD: 1.0, PI: Math.PI, QOPI: 1 / (4 * Math.PI),
        HOPI: 1 / (2 * Math.PI), QINF: 1.3, ALFA: -0.12, LIMAGE: false, YIMAGE: 2 };
      for (const name of ['S', 'XP', 'YP', 'NX', 'NY', 'APANEL', 'GAM',
        'DZDG', 'DZDN', 'DQDG', 'DZDM', 'DQDM']) ctx[name] = new Float64Array(n);
      ctx.GAMU = Array.from({ length: n }, () => new Float64Array(2));
      for (let i = 0; i < n; i += 1) ctx.GAM[i] = Math.sin(i * 0.21) - 0.03 * i;
      scalc(x, y, ctx.S, n);
      segspl(x, ctx.XP, ctx.S, n);
      segspl(y, ctx.YP, ctx.S, n);
      ncalc(x, y, ctx.S, n, ctx.NX, ctx.NY);
      tecalc(ctx);
      apcalc(ctx);

      for (const image of [false, true]) {
        ctx.LIMAGE = image;
        const evaluate = createStreamfunctionEvaluator(ctx);
        const points = [];
        for (let j = 0; j < 13; j += 1) {
          for (let i = 0; i < 17; i += 1) points.push([-0.4 + 0.1 * i, -0.7 + 0.13 * j]);
        }
        if (!image) {
          for (let i = 0; i < n; i += 1) points.push([x[i], y[i]]);
        }
        for (const [xi, yi] of points) {
          const expected = psilin(n, xi, yi, 0, 0, false, false, ctx).psi;
          assert.ok(Number.isFinite(expected));
          assert.equal(evaluate(xi, yi), expected,
            `NACA ${ides}, sharp=${sharp}, rotation=${rotation}, image=${image}, (${xi},${yi})`);
          comparisons += 1;
        }
      }

      // Carrying a shared endpoint must restart after skipped zero-length
      // panels, including a run of duplicates or a missing first/last panel.
      for (const duplicates of [[1], [n - 1], [8, 9, 10], [1, 2, n - 2, n - 1]]) {
        const changed = structuredClone(ctx);
        changed.LIMAGE = false;
        for (const i of duplicates) {
          changed.X[i] = changed.X[i - 1];
          changed.Y[i] = changed.Y[i - 1];
        }
        // Exercise the panel kernel directly: the spline preprocessor rejects
        // duplicated first/last points. Keep its existing endpoint tangents.
        scalc(changed.X, changed.Y, changed.S, n);
        tecalc(changed);
        apcalc(changed);
        const evaluate = createStreamfunctionEvaluator(changed);
        const points = [[-0.2, 0.3], [1.2, -0.4], [changed.X[0], changed.Y[0]],
          [changed.X[n - 1], changed.Y[n - 1]], [changed.X[10], changed.Y[10]]];
        // Revisit earlier points so no endpoint state can leak across calls.
        for (const [xi, yi] of [...points, ...points.slice().reverse()]) {
          const expected = psilin(n, xi, yi, 0, 0, false, false, changed).psi;
          assert.ok(Number.isFinite(expected));
          assert.equal(evaluate(xi, yi), expected,
            `NACA ${ides}, duplicate endpoints ${duplicates}, (${xi},${yi})`);
          comparisons += 1;
        }
      }

      ctx.LIMAGE = false;
      const evaluate = createStreamfunctionEvaluator(ctx);
      for (const [xi, yi] of [[0, -0], [-0, 0], [Number.MIN_VALUE, -Number.MIN_VALUE],
        [1e-155, -1e-155], [1e155, 1e155], [Infinity, 1], [1, -Infinity], [NaN, 0]]) {
        const expected = psilin(n, xi, yi, 0, 0, false, false, ctx).psi;
        assert.equal(evaluate(xi, yi), expected, 'Extreme evaluation coordinates');
        comparisons += 1;
      }

      if (rotation === 0) {
        // Check the avoided expensive work as well as the numerical result.
        // This is per evaluation point, with no memoized grid approximation.
        const expected = psilin(n, -0.4, 0.7, 0, 0, false, false, ctx).psi;
        const originalLog = Math.log;
        let logarithms = 0;
        try {
          Math.log = value => { logarithms += 1; return originalLog(value); };
          assert.equal(evaluate(-0.4, 0.7), expected);
        } finally {
          Math.log = originalLog;
        }
        assert.ok(logarithms <= n + 1, `At most one logarithm per endpoint: ${logarithms}`);
        comparisons += 1;
      }
    }
  }
}
console.log(`streamfunction.test.js: OK (${comparisons} exact PSILIN comparisons)`);
