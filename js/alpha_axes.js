function applyZoomToRange(min, max, zoom, clampMin, clampMax, center) {
  if (!Number.isFinite(zoom)) return { min, max };
  if (zoom === 1.0 && !Number.isFinite(center)) return { min, max };
  const mid = Number.isFinite(center) ? center : 0.5 * (min + max);
  const span = Math.max((max - min) / zoom, 1.0e-6);
  let nextMin = mid - span * 0.5;
  let nextMax = mid + span * 0.5;
  if (Number.isFinite(clampMin) && Number.isFinite(clampMax)) {
    if (nextMin < clampMin) {
      nextMax = Math.min(clampMax, nextMax + (clampMin - nextMin));
      nextMin = clampMin;
    }
    if (nextMax > clampMax) {
      nextMin = Math.max(clampMin, nextMin - (nextMax - clampMax));
      nextMax = clampMax;
    }
  }
  return { min: nextMin, max: nextMax };
}

export function computeAlphaAxisRanges({
  points,
  visibility,
  zoom = 1.0,
  zoomCenter = null,
}) {
  const leftVals = points.flatMap((p) => {
    const vals = [];
    if (visibility.cl) vals.push(p.cl);
    if (visibility.cm) vals.push(p.cm);
    return vals;
  }).filter(Number.isFinite);
  const minLeftVal = Math.min(...leftVals);
  const maxLeftVal = Math.max(...leftVals);
  const rightVals = points.flatMap((p) => {
    const vals = [];
    if (visibility.cd) vals.push(p.cd);
    if (
      visibility.ld
      && Number.isFinite(p.cl)
      && Number.isFinite(p.cd)
      && p.cd !== 0.0
    ) {
      vals.push(p.cl / p.cd);
    }
    return vals;
  }).filter(Number.isFinite);

  let yminL = Math.min(...leftVals);
  let ymaxL = Math.max(...leftVals);
  if (!Number.isFinite(yminL) || !Number.isFinite(ymaxL)) {
    yminL = -1.0;
    ymaxL = 1.0;
  } else if (yminL === ymaxL) {
    yminL -= 0.1;
    ymaxL += 0.1;
  }
  const padL = 0.12 * (ymaxL - yminL || 1.0);
  yminL -= padL;
  ymaxL += padL;
  yminL = Math.max(yminL, -2.0);
  ymaxL = Math.min(ymaxL, 2.0);

  const dataMinR = Math.min(...rightVals);
  const dataMaxR = Math.max(...rightVals);
  let yminR = dataMinR;
  let ymaxR = dataMaxR;
  if (!Number.isFinite(yminR) || !Number.isFinite(ymaxR)) {
    yminR = -0.1;
    ymaxR = 0.1;
  } else if (yminR === ymaxR) {
    yminR -= 0.01;
    ymaxR += 0.01;
  }
  const padR = 0.2 * (ymaxR - yminR || 1.0);
  yminR -= padR;
  ymaxR += padR;
  const rightClampAbs = visibility.ld ? 500.0 : 2.0;
  yminR = Math.max(yminR, -rightClampAbs);
  ymaxR = Math.min(ymaxR, rightClampAbs);

  const hasLeftSeries = visibility.cl || visibility.cm;
  const hasRightSeries = rightVals.length > 0;
  const rightIsLdOnly = visibility.ld && !visibility.cd;

  if (hasLeftSeries) {
    if (yminL > 0.0) yminL = 0.0;
    if (ymaxL < 0.0) ymaxL = 0.0;
  }

  if (hasRightSeries) {
    if (!rightIsLdOnly) {
      yminR = 0.0;
    } else if (yminR > 0.0) {
      yminR = 0.0;
    }
    if (ymaxR < 0.0) ymaxR = 0.0;
    yminR = Math.max(yminR, -rightClampAbs);
    ymaxR = Math.min(ymaxR, rightClampAbs);
  }

  // Avoid pinning the left zero to an edge when aligning with CL/CD.
  if (hasLeftSeries && rightIsLdOnly && Number.isFinite(minLeftVal) && Number.isFinite(maxLeftVal)) {
    const leftSpan = Math.max(ymaxL - yminL, 1.0e-6);
    if (maxLeftVal <= 0.0) {
      ymaxL = Math.max(ymaxL, 0.15 * leftSpan);
    } else if (minLeftVal >= 0.0) {
      yminL = Math.min(yminL, -0.15 * leftSpan);
    }
    const zeroFracLeft = (0.0 - yminL) / Math.max(ymaxL - yminL, 1.0e-6);
    const maxZeroFrac = 0.8;
    const targetZeroFrac = 0.7;
    if (zeroFracLeft > maxZeroFrac && yminL < 0.0) {
      const spanForTarget = (-yminL) / targetZeroFrac;
      const targetYmax = yminL + spanForTarget;
      ymaxL = Math.max(ymaxL, targetYmax);
    } else if (zeroFracLeft < (1.0 - maxZeroFrac) && ymaxL > 0.0) {
      const spanForTarget = ymaxL / targetZeroFrac;
      const targetYmin = ymaxL - spanForTarget;
      yminL = Math.min(yminL, targetYmin);
    }
    yminL = Math.max(yminL, -2.0);
    ymaxL = Math.min(ymaxL, 2.0);
  }

  const baseYL = { min: yminL, max: ymaxL };
  const baseYR = { min: yminR, max: ymaxR };

  const zoomYL = applyZoomToRange(
    yminL,
    ymaxL,
    zoom,
    baseYL.min,
    baseYL.max,
    zoomCenter?.yL,
  );
  yminL = zoomYL.min;
  ymaxL = zoomYL.max;

  const zoomYR = applyZoomToRange(
    yminR,
    ymaxR,
    zoom,
    baseYR.min,
    baseYR.max,
    zoomCenter?.yR,
  );
  yminR = zoomYR.min;
  ymaxR = zoomYR.max;

  // Keep zero available for the alignment cases even after zooming.
  if (hasLeftSeries && rightIsLdOnly) {
    if (yminL > 0.0) yminL = 0.0;
    if (ymaxL < 0.0) ymaxL = 0.0;
  }
  if (hasRightSeries) {
    if (rightIsLdOnly) {
      if (yminR > 0.0) yminR = 0.0;
      if (ymaxR < 0.0) ymaxR = 0.0;
    } else {
      // For CD, keep zero pinned to the bottom even when zoomed.
      yminR = 0.0;
      if (ymaxR < 0.0) ymaxR = 0.0;
    }
  }

  // Re-apply left-axis headroom after zoom to keep zero away from the edge.
  if (hasLeftSeries && rightIsLdOnly && Number.isFinite(minLeftVal) && Number.isFinite(maxLeftVal)) {
    const leftSpan = Math.max(ymaxL - yminL, 1.0e-6);
    if (maxLeftVal <= 0.0) {
      ymaxL = Math.max(ymaxL, 0.15 * leftSpan);
    } else if (minLeftVal >= 0.0) {
      yminL = Math.min(yminL, -0.15 * leftSpan);
    }
    const zeroFracLeft = (0.0 - yminL) / Math.max(ymaxL - yminL, 1.0e-6);
    const maxZeroFrac = 0.8;
    const targetZeroFrac = 0.7;
    if (zeroFracLeft > maxZeroFrac && yminL < 0.0) {
      const spanForTarget = (-yminL) / targetZeroFrac;
      const targetYmax = yminL + spanForTarget;
      ymaxL = Math.max(ymaxL, targetYmax);
    } else if (zeroFracLeft < (1.0 - maxZeroFrac) && ymaxL > 0.0) {
      const spanForTarget = ymaxL / targetZeroFrac;
      const targetYmin = ymaxL - spanForTarget;
      yminL = Math.min(yminL, targetYmin);
    }
    yminL = Math.max(yminL, -2.0);
    ymaxL = Math.min(ymaxL, 2.0);
  }

  if (hasLeftSeries && hasRightSeries && rightIsLdOnly) {
    const zeroFracLeft = (0.0 - yminL) / (ymaxL - yminL);
    const f = Math.min(Math.max(zeroFracLeft, 1.0e-3), 1.0 - 1.0e-3);
    const dataSpan = dataMaxR - dataMinR || 1.0;
    const safetyPad = 0.08 * dataSpan;
    const minNeeded = dataMinR - safetyPad;
    const maxNeeded = dataMaxR + safetyPad;
    const spanFromMin = minNeeded < 0.0 ? (-minNeeded) / f : 0.0;
    const spanFromMax = maxNeeded > 0.0 ? maxNeeded / (1.0 - f) : 0.0;
    const spanNeeded = Math.max(spanFromMin, spanFromMax, 1.0e-6);
    yminR = -f * spanNeeded;
    ymaxR = (1.0 - f) * spanNeeded;
    yminR = Math.max(yminR, -rightClampAbs);
    ymaxR = Math.min(ymaxR, rightClampAbs);
  } else if (hasRightSeries && Number.isFinite(dataMinR) && Number.isFinite(dataMaxR)) {
    const dataSpan = dataMaxR - dataMinR || 1.0;
    const safetyPad = 0.08 * dataSpan;
    if (dataMinR < 0.0) {
      yminR = Math.min(yminR, dataMinR - safetyPad);
    } else {
      yminR = 0.0;
    }
    ymaxR = Math.max(ymaxR, dataMaxR + safetyPad);
    yminR = Math.max(yminR, -rightClampAbs);
    ymaxR = Math.min(ymaxR, rightClampAbs);
  }

  return {
    yminL,
    ymaxL,
    yminR,
    ymaxR,
    hasLeftSeries,
    hasRightSeries,
    rightIsLdOnly,
  };
}

export function yToPy(y, ymin, ymax, top, plotH) {
  return top + (1.0 - (y - ymin) / (ymax - ymin)) * plotH;
}
