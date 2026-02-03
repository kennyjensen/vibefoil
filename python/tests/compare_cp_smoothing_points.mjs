import fs from 'fs';

function buildCpEditPoints(points, lePt, tePt) {
  if (!points || !points.length) return null;
  const dxChord = tePt.x - lePt.x;
  const dyChord = tePt.y - lePt.y;
  const chord2 = dxChord * dxChord + dyChord * dyChord || 1.0;
  let leIdx = 0;
  let leDist = Infinity;
  for (let i = 0; i < points.length; i += 1) {
    const p = points[i];
    const dx = p.x - lePt.x;
    const dy = p.y - lePt.y;
    const dist = dx * dx + dy * dy;
    if (dist < leDist) {
      leDist = dist;
      leIdx = i;
    }
  }
  const out = new Array(points.length);
  for (let i = 0; i < points.length; i += 1) {
    const p = points[i];
    const s = ((p.x - lePt.x) * dxChord + (p.y - lePt.y) * dyChord) / chord2;
    const side = i <= leIdx ? 'upper' : 'lower';
    out[i] = { index: i, s, cp: p.cp, side };
  }
  return out;
}

function smoothCpPoints(points) {
  if (!points || points.length < 3) return points;
  const smoothSide = (side) => {
    const idxs = [];
    for (let i = 0; i < points.length; i += 1) {
      if (points[i].side === side) idxs.push(i);
    }
    if (idxs.length < 3) return;
    const nextVals = new Array(idxs.length);
    for (let i = 0; i < idxs.length; i += 1) {
      const idx = idxs[i];
      if (i === 0 || i === idxs.length - 1) {
        nextVals[i] = points[idx].cp;
        continue;
      }
      const prev = points[idxs[i - 1]].cp;
      const cur = points[idx].cp;
      const next = points[idxs[i + 1]].cp;
      nextVals[i] = (prev + 2.0 * cur + next) * 0.25;
    }
    for (let i = 0; i < idxs.length; i += 1) {
      points[idxs[i]].cp = nextVals[i];
    }
  };
  smoothSide('upper');
  smoothSide('lower');
  return points;
}

function calcLeSpike(points, side) {
  const sidePts = points.filter((p) => p.side === side);
  if (sidePts.length < 3) return { spike: null, count: sidePts.length, leIdx: -1 };
  let minIdx = 0;
  for (let i = 1; i < sidePts.length; i += 1) {
    if (sidePts[i].s < sidePts[minIdx].s) minIdx = i;
  }
  const candidates = [minIdx - 1, minIdx, minIdx + 1].filter((i) => i >= 1 && i <= sidePts.length - 2);
  let maxSpike = 0.0;
  for (const idx of candidates) {
    const prev = sidePts[idx - 1].cp;
    const cur = sidePts[idx].cp;
    const next = sidePts[idx + 1].cp;
    const spike = Math.abs(2.0 * cur - prev - next);
    if (spike > maxSpike) maxSpike = spike;
  }
  return { spike: maxSpike, count: sidePts.length, leIdx: minIdx };
}

function main() {
  const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
  const points = payload.points;
  const le = payload.le;
  const te = payload.te;
  const editPoints = buildCpEditPoints(points, le, te);
  const beforeUpper = calcLeSpike(editPoints, 'upper');
  const beforeLower = calcLeSpike(editPoints, 'lower');
  const smoothed = editPoints.map((p) => ({ ...p }));
  smoothCpPoints(smoothed);
  const afterUpper = calcLeSpike(smoothed, 'upper');
  const afterLower = calcLeSpike(smoothed, 'lower');
  const results = {
    before: { upper: beforeUpper, lower: beforeLower },
    after: { upper: afterUpper, lower: afterLower },
  };
  process.stdout.write(JSON.stringify({ results }));
}

main();
