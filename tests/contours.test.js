import assert from 'node:assert/strict';
import { createContourCache } from '../js/contours.js';

// Capture the old renderer's exact segment order and world coordinates without
// a browser. This retains its per-cell points and two/four-edge branches.
function referenceContours(bounds, grid, nx, ny, isoValues) {
  const contours = [];
  let points;
  const ctx = {
    beginPath() { points = []; },
    moveTo(x, y) { points.push(x, y); },
    lineTo(x, y) { points.push(x, y); },
    stroke() { contours.push(Float64Array.from(points)); },
  };
  const worldToCanvas = (x, y) => ({ x, y });

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

  const xSpan = bounds.xmax - bounds.xmin;
  const ySpan = bounds.ymax - bounds.ymin;
  const dx = xSpan / (nx - 1);
  const dy = ySpan / (ny - 1);

  isoValues.forEach((iso) => {
    ctx.beginPath();
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

        const points = [];

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
          points.push({ x: px, y: py });
        }

        if (points.length === 2) {
          const p0 = worldToCanvas(points[0].x, points[0].y, bounds);
          const p1 = worldToCanvas(points[1].x, points[1].y, bounds);
          ctx.moveTo(p0.x, p0.y);
          ctx.lineTo(p1.x, p1.y);
        } else if (points.length === 4) {
          const p0 = worldToCanvas(points[0].x, points[0].y, bounds);
          const p1 = worldToCanvas(points[1].x, points[1].y, bounds);
          const p2 = worldToCanvas(points[2].x, points[2].y, bounds);
          const p3 = worldToCanvas(points[3].x, points[3].y, bounds);
          ctx.moveTo(p0.x, p0.y);
          ctx.lineTo(p1.x, p1.y);
          ctx.moveTo(p2.x, p2.y);
          ctx.lineTo(p3.x, p3.y);
        }
      }
    }
    ctx.stroke();
  });
  return contours;
}

const getContours = createContourCache();
const bounds = { xmin: -0.37, xmax: 1.42, ymin: -0.28, ymax: 0.83 };
let comparisons = 0;
function check(grid, nx, ny, levels, domain = bounds) {
  const actual = getContours(domain, grid, nx, ny, levels);
  assert.deepStrictEqual(actual, referenceContours(domain, grid, nx, ny, levels), `case=${comparisons}`);
  assert.equal(getContours({ ...domain }, grid, nx, ny, levels.slice()), actual, 'Identical inputs reuse segments');
  comparisons += 1;
  return actual;
}
for (let mask = 0; mask < 16; mask += 1) {
  const corners = [0, 1, 3, 2].map(bit => mask & (1 << bit) ? 1 : -1);
  check(Float64Array.from(corners), 2, 2, [-1, 0, 1]);
}
for (const value of [-0, Number.MIN_VALUE, Number.MAX_VALUE, NaN, Infinity, -Infinity]) {
  for (let corner = 0; corner < 4; corner += 1) {
    const grid = Float64Array.of(-1, 1, 1, -1);
    grid[corner] = value;
    check(grid, 2, 2, [-1, 0, 1]);
  }
}
let seed = 0x24122301;
function random() {
  seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
  return seed / 0x100000000 - 0.5;
}
for (const [nx, ny] of [[1, 1], [2, 7], [7, 2], [13, 9], [60, 32]]) {
  const grid = Float64Array.from({ length: nx * ny }, random);
  const levels = Array.from({ length: 14 }, (_, i) => -0.4 + i * 0.8 / 13);
  check(grid, nx, ny, levels);
  check(grid, nx, ny, []);
}
const grid = Float64Array.from({ length: 48 }, random);
const levels = [-0.25, 0, 0.25];
let previous = check(grid, 8, 6, levels);
// Screen changes must reuse world coordinates; domain changes must rebuild.
assert.equal(getContours({ ...bounds, zoom: 2, panX: 17, angle: 0.2, width: 800, height: 600 }, grid, 8, 6, levels), previous);
for (const key of ['xmin', 'xmax', 'ymin', 'ymax']) {
  const domain = { ...bounds, [key]: bounds[key] + 0.125 };
  const next = check(grid, 8, 6, levels, domain);
  assert.notEqual(next, previous);
  previous = next;
}
for (const [nx, ny] of [[6, 8], [4, 12], [8, 6]]) {
  const next = check(grid, nx, ny, levels);
  assert.notEqual(next, previous);
  previous = next;
}
levels[1] = 0.125;
let next = check(grid, 8, 6, levels);
assert.notEqual(next, previous, 'In-place level changes invalidate the cache');
previous = next;
for (const value of [0, -0, NaN, Infinity, -Infinity, 0.3]) {
  grid[17] = value;
  next = check(grid, 8, 6, levels);
  assert.notEqual(next, previous, 'In-place grid changes invalidate the cache');
  previous = next;
}
check(grid.slice(), 8, 6, Float64Array.from(levels));
const sparseLevels = [-0.25, undefined, 0.25];
previous = check(grid, 8, 6, sparseLevels);
delete sparseLevels[1];
next = check(grid, 8, 6, sparseLevels);
assert.notEqual(next, previous, 'A skipped level differs from an explicit undefined level');
const mutableBounds = { ...bounds };
previous = check(grid, 8, 6, levels, mutableBounds);
mutableBounds.xmin += 0.125;
next = check(grid, 8, 6, levels, mutableBounds);
assert.notEqual(next, previous, 'In-place domain changes invalidate the cache');
console.log(`contours: ${comparisons} exact segment comparisons and cache invalidation checks passed`);
