import assert from 'node:assert/strict';
import { computeAlphaAxisRanges, yToPy } from '../js/alpha_axes.js';

function nearlyEqual(a, b, eps = 1.0e-6) {
  return Math.abs(a - b) <= eps;
}

// Regression test: when only C_M and C_L/C_D are visible, the y=0 screen
// position must match across the left and right axes.
{
  const points = [
    { cl: -50.0, cm: -0.045, cd: 1.0 },
    { cl: -10.0, cm: -0.035, cd: 1.0 },
    { cl: 20.0, cm: -0.025, cd: 1.0 },
    { cl: 50.0, cm: -0.015, cd: 1.0 },
  ];
  const visibility = {
    cl: false,
    cm: true,
    cd: false,
    ld: true,
  };
  for (const zoom of [1.0, 3.0]) {
    const { yminL, ymaxL, yminR, ymaxR, hasLeftSeries, hasRightSeries } = computeAlphaAxisRanges({
      points,
      visibility,
      zoom,
      zoomCenter: null,
    });

    assert.equal(hasLeftSeries, true, 'expected left axis to be active');
    assert.equal(hasRightSeries, true, 'expected right axis to be active');
    assert.ok(yminL <= 0.0 && ymaxL >= 0.0, 'left axis must include zero');
    assert.ok(yminR <= 0.0 && ymaxR >= 0.0, 'right axis must include zero');

    const top = 30;
    const plotH = 340 - top - 36;
    const y0Left = yToPy(0.0, yminL, ymaxL, top, plotH);
    const y0Right = yToPy(0.0, yminR, ymaxR, top, plotH);

    assert.ok(
      nearlyEqual(y0Left, y0Right),
      `zoom=${zoom} expected aligned zeros, got left=${y0Left} right=${y0Right}`,
    );
  }
}

console.log('alpha_axes.test.js: OK');
