// Marching-square contour coordinates, with the original edge ordering and
// interpolation arithmetic. Each group of four values describes one segment.
const edgeTable = [
  [],
  [3, 0],
  [0, 1],
  [3, 1],
  [1, 2],
  [3, 0, 1, 2],
  [0, 2],
  [3, 2],
  [2, 3],
  [0, 2],
  [0, 1, 2, 3],
  [1, 2],
  [1, 3],
  [0, 1],
  [3, 0],
  [],
];

function buildContourSegments(bounds, grid, nx, ny, isoValues) {
  const xSpan = bounds.xmax - bounds.xmin;
  const ySpan = bounds.ymax - bounds.ymin;
  const dx = xSpan / (nx - 1);
  const dy = ySpan / (ny - 1);

  const contours = [];
  isoValues.forEach((iso) => {
    const points = [];
    for (let j = 0; j < ny - 1; j += 1) {
      for (let i = 0; i < nx - 1; i += 1) {
        const idx = j * nx + i;
        const v0 = grid[idx];
        const v1 = grid[idx + 1];
        const v2 = grid[idx + 1 + nx];
        const v3 = grid[idx + nx];

        const caseIndex = (v0 > iso ? 1 : 0)
          | (v1 > iso ? 2 : 0)
          | (v2 > iso ? 4 : 0)
          | (v3 > iso ? 8 : 0);

        const edges = edgeTable[caseIndex];
        if (edges.length === 0) {
          continue;
        }

        const x0 = bounds.xmin + i * dx;
        const y0 = bounds.ymin + j * dy;

        for (let e = 0; e < edges.length; e += 1) {
          const edge = edges[e];
          let ax;
          let ay;
          let bx;
          let by;
          let va;
          let vb;

          if (edge === 0) {
            ax = x0;
            ay = y0;
            bx = x0 + dx;
            by = y0;
            va = v0;
            vb = v1;
          } else if (edge === 1) {
            ax = x0 + dx;
            ay = y0;
            bx = x0 + dx;
            by = y0 + dy;
            va = v1;
            vb = v2;
          } else if (edge === 2) {
            ax = x0 + dx;
            ay = y0 + dy;
            bx = x0;
            by = y0 + dy;
            va = v2;
            vb = v3;
          } else {
            ax = x0;
            ay = y0 + dy;
            bx = x0;
            by = y0;
            va = v3;
            vb = v0;
          }

          const denom = vb - va;
          const t = denom !== 0.0 ? (iso - va) / denom : 0.5;
          const px = ax + (bx - ax) * t;
          const py = ay + (by - ay) * t;
          points.push(px, py);
        }
      }
    }
    contours.push(Float64Array.from(points));
  });
  return contours;
}

// Retain one grid's world-space contour segments. Check source values as well
// as identity so in-place edits and nonfinite/signed-zero changes invalidate it.
// Camera transforms are deliberately applied later by the renderer. Returned
// segment arrays are read-only; the renderer never modifies them.
export function createContourCache() {
  let cachedGrid = null;
  let gridValues = null;
  let parameters = null;
  let contours = null;
  return function getContours(bounds, grid, nx, ny, isoValues) {
    const next = [bounds.xmin, bounds.xmax, bounds.ymin, bounds.ymax, nx, ny, isoValues.length];
    isoValues.forEach((value, index) => next.push(index, value));
    let matches = grid === cachedGrid && parameters?.length === next.length
      && gridValues?.length === grid.length;
    if (matches) {
      for (let i = 0; i < next.length; i += 1) {
        if (!Object.is(parameters[i], next[i])) { matches = false; break; }
      }
    }
    if (matches) {
      for (let i = 0; i < grid.length; i += 1) {
        if (!Object.is(gridValues[i], grid[i])) { matches = false; break; }
      }
    }
    if (!matches) {
      const segments = buildContourSegments(bounds, grid, nx, ny, isoValues);
      cachedGrid = grid;
      gridValues = Array.from(grid);
      parameters = next;
      contours = segments;
    }
    return contours;
  };
}
