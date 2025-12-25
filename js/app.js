import { scalc, segspl, seval } from './spline.js';
import {
  apcalc,
  tecalc,
  xywake,
  qwcalc,
  ncalc,
  psilin,
  ggcalc,
  qdcalc,
} from './xpanel.js';
import {
  blpini,
  iblsys,
  mrcl,
  comset,
  mrchue,
  setbl,
  blsolv,
  update as blUpdate,
} from './xbl.js';
import {
  NCOM,
  ensureCtx,
  syncComToVars,
  syncVarsToCom,
  blprv,
  blkin,
  blvar,
  blmid,
  trchek,
  tesys,
  blsys,
  hkin,
} from './xblsys.js';

// High-level orchestrator for the XFOIL port: UI, geometry generation,
// inviscid panel solve, viscous BL coupling, and plotting.
const canvas = document.getElementById('plot');
const ctx = canvas.getContext('2d');
const cpCanvas = document.getElementById('cpPlot');
const cpCtx = cpCanvas.getContext('2d');
const alphaCanvas = document.getElementById('alphaPlot');
const alphaCtx = alphaCanvas ? alphaCanvas.getContext('2d') : null;
const polarCanvas = document.getElementById('polarPlot');
const polarCtx = polarCanvas ? polarCanvas.getContext('2d') : null;

// UI controls mirrored from the DOM; values drive geometry and solver setup.
const seriesRadios = Array.from(document.querySelectorAll('input[name="series"]'));
const sourceRadios = Array.from(document.querySelectorAll('input[name="source"]'));
const nacaOptions = document.getElementById('nacaOptions');
const customOptions = document.getElementById('customOptions');
const controls4 = document.getElementById('controls4');
const controls5 = document.getElementById('controls5');
const controls6 = document.getElementById('controls6');
const dataRe = document.getElementById('dataRe');
const dataAlpha = document.getElementById('dataAlpha');
const dataCl = document.getElementById('dataCl');
const dataCm = document.getElementById('dataCm');
const dataCd = document.getElementById('dataCd');
const dataLd = document.getElementById('dataLd');
const dataNcr = document.getElementById('dataNcr');
const dataConverged = document.getElementById('dataConverged');

const mSlider = document.getElementById('m');
const pSlider = document.getElementById('p');
const tSlider = document.getElementById('t');
const t5Slider = document.getElementById('t5');
const series5Select = document.getElementById('series5');
const t6Slider = document.getElementById('t6');
const series6Profile = document.getElementById('series6Profile');
const cl6Input = document.getElementById('cl6');

const mValue = document.getElementById('mValue');
const pValue = document.getElementById('pValue');
const tValue = document.getElementById('tValue');
const t5Value = document.getElementById('t5Value');
const t6Value = document.getElementById('t6Value');
const dataName = document.getElementById('dataName');
const viscousToggle = document.getElementById('viscous');
const machInput = document.getElementById('mach');
const reynoldsInput = document.getElementById('reynolds');
const ncrInput = document.getElementById('ncr');
const loadDatButton = document.getElementById('loadDat');
const datFileInput = document.getElementById('datFile');
const alphaSlider = document.getElementById('alpha');
const alphaValue = document.getElementById('alphaValue');
const alphaMinus = document.getElementById('alphaMinus');
const alphaPlus = document.getElementById('alphaPlus');
const alphaSweep = document.getElementById('alphaSweep');
const sweepStartInput = document.getElementById('sweepStart');
const sweepEndInput = document.getElementById('sweepEnd');
const sweepIncInput = document.getElementById('sweepInc');
const runCaseList = document.getElementById('runCaseList');
const addRunCaseButton = document.getElementById('addRunCase');

// Geometry buffers (airfoil surface and custom inputs); sizes follow XFOIL nside.
const nside = 80;
const xx = new Float64Array(nside);
const yt = new Float64Array(nside);
const yc = new Float64Array(nside);
const xb = new Float64Array(2 * nside);
const yb = new Float64Array(2 * nside);
const xbCustom = new Float64Array(2 * nside);
const ybCustom = new Float64Array(2 * nside);
let customAirfoil = null;
let currentAirfoilName = 'NACA 2412';

// Panel-method working state and run history.
let panelCtx = null;
let panelX = null;
let panelY = null;
let panelXP = null;
let panelYP = null;
const runCases = [];
let activeCaseId = null;
let nextCaseId = 1;
let sweeping = false;

// Matrix helpers mirror Fortran 1-based indexing expectations in BL modules.
function createMatrix(rows, cols) {
  const matrix = new Array(rows);
  for (let i = 0; i < rows; i += 1) {
    matrix[i] = new Float64Array(cols);
  }
  return matrix;
}

function createMatrix1(rows, cols) {
  const matrix = new Array(rows + 1);
  for (let i = 0; i <= rows; i += 1) {
    matrix[i] = new Float64Array(cols + 1);
  }
  return matrix;
}

function createTensor3(d1, d2, d3) {
  const tensor = new Array(d1 + 1);
  for (let i = 0; i <= d1; i += 1) {
    tensor[i] = new Array(d2 + 1);
    for (let j = 0; j <= d2; j += 1) {
      tensor[i][j] = new Float64Array(d3 + 1);
    }
  }
  return tensor;
}

// NACA 4-digit designation formatting for UI/overlay.
function updateLabels(m, p, t) {
  mValue.textContent = `${m}%`;
  pValue.textContent = `${p * 10}%`;
  tValue.textContent = `${t}%`;
  currentAirfoilName = `NACA ${m}${p}${String(t).padStart(2, '0')}`;
}

// Numeric formatting with NaN guard for overlay display.
function formatNum(value, digits = 4) {
  if (!Number.isFinite(value)) return '—';
  return value.toFixed(digits);
}

// Overlay summary: converts solver state to UI-friendly engineering values.
function updateDataBox(alphaRad, coeffs, converged = true) {
  const alphaDeg = alphaRad / (Math.PI / 180.0);
  dataName.textContent = currentAirfoilName || '—';
  dataAlpha.textContent = `${alphaDeg.toFixed(2)}°`;
  dataRe.textContent = coeffs?.REINF1 ? `${Math.round(coeffs.REINF1)}` : '—';
  dataCl.textContent = formatNum(coeffs?.CL, 4);
  dataCm.textContent = formatNum(coeffs?.CM, 4);
  dataCd.textContent = formatNum(coeffs?.CD, 5);
  if (Number.isFinite(coeffs?.CL) && Number.isFinite(coeffs?.CD) && coeffs.CD !== 0.0) {
    dataLd.textContent = formatNum(coeffs.CL / coeffs.CD, 2);
  } else {
    dataLd.textContent = '—';
  }
  const ncr = coeffs?.ACRIT?.[1];
  dataNcr.textContent = Number.isFinite(ncr) ? `${ncr.toFixed(2)}` : '—';
  if (dataConverged) {
    dataConverged.hidden = converged;
  }
}

// 5-digit designation formatting (design camber line + thickness).
function updateLabels5(series, t) {
  t5Value.textContent = `${t}%`;
  currentAirfoilName = `NACA ${series}${String(t).padStart(2, '0')}`;
}

// 6-series designation formatting: profile + design CL digit + thickness.
function updateLabels6(profile, t, camber, cl) {
  t6Value.textContent = `${t}%`;
  const tocCode = String(t).padStart(2, '0');
  let clCode = '0';
  if (camber !== '0' && Number.isFinite(cl)) {
    clCode = Math.min(9, Math.max(0, Math.round(cl * 10))).toString();
  }
  currentAirfoilName = `NACA ${profile}-${clCode}${tocCode}`;
}

// Standard location of minimum pressure for 6-series (second digit).
function defaultSixSeriesA(profile) {
  const digit = parseInt(String(profile).charAt(1), 10);
  if (!Number.isFinite(digit)) return 0.8;
  return Math.min(0.9, Math.max(0.1, digit / 10));
}

// Camber-line selection inferred from CL and 6A suffix.
function inferSixSeriesCamber(profile, cl) {
  if (!Number.isFinite(cl) || Math.abs(cl) < 1.0e-6) return '0';
  const prof = String(profile).toUpperCase();
  return prof.endsWith('A') ? '6A' : '6';
}


function updateAlphaLabel(alpha) {
  alphaValue.textContent = `${alpha.toFixed(1)}°`;
}

// Subscript labeling to mimic textbook notation (C_L, C_D, etc.).
function drawSubLabel(ctx2d, base, sub, x, y) {
  const baseFont = ctx2d.font;
  ctx2d.fillText(base, x, y);
  const baseWidth = ctx2d.measureText(base).width;
  ctx2d.save();
  ctx2d.font = '10px Consolas, "Courier New", monospace';
  ctx2d.fillText(sub, x + baseWidth + 1, y + 4);
  ctx2d.restore();
  ctx2d.font = baseFont;
}

// Clamp and trigger update for alpha changes (slider, buttons, sweep).
function setAlphaValue(alphaDeg) {
  const min = parseFloat(alphaSlider.min);
  const max = parseFloat(alphaSlider.max);
  const clamped = Math.max(min, Math.min(max, alphaDeg));
  alphaSlider.value = `${clamped}`;
  update();
}

// Run-case structure for multi-sweep overlays and color separation.
function createRunCase() {
  const colors = ['#ff4a3d', '#ff9f1a', '#ffd166', '#6dd36f', '#2f7bff', '#9b5de5'];
  const color = colors[(nextCaseId - 1) % colors.length];
  const id = nextCaseId;
  nextCaseId += 1;
  return {
    id,
    name: `Case ${id}`,
    color,
    history: [],
    sweep: null,
  };
}

// Active run case determines which history accumulates.
function getActiveCase() {
  if (runCases.length === 0) return null;
  if (activeCaseId == null) return null;
  return runCases.find((c) => c.id === activeCaseId) || null;
}

// Render run-case list with editable labels and active highlighting.
function renderRunCases() {
  if (!runCaseList) return;
  runCaseList.innerHTML = '';
  runCases.forEach((caseItem) => {
    const item = document.createElement('div');
    item.className = `run-case-item${caseItem.id === activeCaseId ? ' active' : ''}`;
    item.dataset.caseId = `${caseItem.id}`;

    const meta = document.createElement('div');
    meta.className = 'run-case-meta';

    const titleInput = document.createElement('input');
    titleInput.type = 'text';
    titleInput.value = caseItem.name;
    titleInput.className = 'run-case-title';
    titleInput.addEventListener('input', (event) => {
      caseItem.name = event.target.value;
    });
    titleInput.addEventListener('click', (event) => {
      event.stopPropagation();
    });
    titleInput.addEventListener('keydown', (event) => {
      event.stopPropagation();
    });

    meta.appendChild(titleInput);

    const colorInput = document.createElement('input');
    colorInput.type = 'color';
    colorInput.value = caseItem.color;
    colorInput.addEventListener('input', (event) => {
      caseItem.color = event.target.value;
      drawAlphaSweepPlot();
      drawPolarPlot();
    });

    const deleteBtn = document.createElement('button');
    deleteBtn.type = 'button';
    deleteBtn.textContent = '🗑';
    deleteBtn.addEventListener('click', (event) => {
      event.stopPropagation();
      const idx = runCases.findIndex((c) => c.id === caseItem.id);
      if (idx >= 0) runCases.splice(idx, 1);
      if (activeCaseId === caseItem.id) {
        activeCaseId = runCases.length ? runCases[0].id : null;
      }
      if (runCases.length === 0) {
        const newCase = createRunCase();
        runCases.push(newCase);
        activeCaseId = newCase.id;
      }
      renderRunCases();
      drawAlphaSweepPlot();
      drawPolarPlot();
    });

    item.addEventListener('click', () => {
      activeCaseId = activeCaseId === caseItem.id ? null : caseItem.id;
      renderRunCases();
    });

    item.appendChild(meta);
    item.appendChild(colorInput);
    item.appendChild(deleteBtn);
    runCaseList.appendChild(item);
  });
}

// Sequential alpha sweep; preserves deterministic ordering for polar/alpha plots.
async function sweepAlpha() {
  if (sweeping) return;
  sweeping = true;
  const start = parseFloat(sweepStartInput?.value ?? '-10');
  const end = parseFloat(sweepEndInput?.value ?? '10');
  const incRaw = parseFloat(sweepIncInput?.value ?? '1');
  const inc = Math.max(Math.abs(incRaw), 0.1);
  const direction = start <= end ? 1 : -1;
  const active = getActiveCase();
  if (active) {
    active.sweep = {
      start: Number.isFinite(start) ? start : -10,
      end: Number.isFinite(end) ? end : 10,
      inc: Number.isFinite(incRaw) ? incRaw : 1,
    };
    renderRunCases();
  }

  for (let a = start; direction > 0 ? a <= end : a >= end; a += direction * inc) {
    setAlphaValue(a);
    await new Promise((resolve) => setTimeout(resolve, 50));
  }
  sweeping = false;
  activeCaseId = null;
  renderRunCases();
}

// Pixel-density aware canvas sizing to keep text/lines crisp.
function resizeCanvas() {
  const ratio = window.devicePixelRatio || 1;
  const rect = canvas.getBoundingClientRect();
  canvas.width = Math.max(1, Math.floor(rect.width * ratio));
  canvas.height = Math.max(1, Math.floor(rect.height * ratio));
  ctx.setTransform(ratio, 0, 0, ratio, 0, 0);

  if (cpCanvas) {
    const cpRect = cpCanvas.getBoundingClientRect();
    cpCanvas.width = Math.max(1, Math.floor(cpRect.width * ratio));
    cpCanvas.height = Math.max(1, Math.floor(cpRect.height * ratio));
    cpCtx.setTransform(ratio, 0, 0, ratio, 0, 0);
  }

  if (alphaCanvas && alphaCtx) {
    const alphaRect = alphaCanvas.getBoundingClientRect();
    alphaCanvas.width = Math.max(1, Math.floor(alphaRect.width * ratio));
    alphaCanvas.height = Math.max(1, Math.floor(alphaRect.height * ratio));
    alphaCtx.setTransform(ratio, 0, 0, ratio, 0, 0);
  }

  if (polarCanvas && polarCtx) {
    const polarRect = polarCanvas.getBoundingClientRect();
    polarCanvas.width = Math.max(1, Math.floor(polarRect.width * ratio));
    polarCanvas.height = Math.max(1, Math.floor(polarRect.height * ratio));
    polarCtx.setTransform(ratio, 0, 0, ratio, 0, 0);
  }
}

// Insert or replace sweep point (alpha-indexed) to keep plots ordered.
function upsertSweepPoint(alphaDeg, coeffs) {
  const active = getActiveCase();
  if (!active) return;
  if (!Number.isFinite(alphaDeg)
    || !Number.isFinite(coeffs?.CL)
    || !Number.isFinite(coeffs?.CD)
    || !Number.isFinite(coeffs?.CM)) {
    return;
  }

  const history = active.history;
  const eps = 1.0e-6;
  const idx = history.findIndex((p) => Math.abs(p.alpha - alphaDeg) < eps);
  const point = { alpha: alphaDeg, cl: coeffs.CL, cd: coeffs.CD, cm: coeffs.CM };
  if (idx >= 0) {
    history[idx] = point;
  } else {
    history.push(point);
    history.sort((a, b) => a.alpha - b.alpha);
  }
}

// Alpha plot: CL/CM on left axis and CD on right, matching XFOIL conventions.
function drawAlphaSweepPlot() {
  if (!alphaCanvas || !alphaCtx) return;
  alphaCtx.clearRect(0, 0, alphaCanvas.width, alphaCanvas.height);
  const allPoints = runCases.flatMap((c) => c.history);
  if (allPoints.length === 0) return;

  const w = alphaCanvas.width;
  const h = alphaCanvas.height;
  const left = 58;
  const right = 72;
  const top = 18;
  const bottom = 36;
  const plotW = w - left - right;
  const plotH = h - top - bottom;

  const alphas = allPoints.map((p) => p.alpha);
  let xmin = Math.min(...alphas);
  let xmax = Math.max(...alphas);
  if (xmin === xmax) {
    xmin -= 1;
    xmax += 1;
  }

  const clcmVals = allPoints.flatMap((p) => [p.cl, p.cm]).filter(Number.isFinite);
  const cdVals = allPoints.map((p) => p.cd).filter(Number.isFinite);
  let yminL = Math.min(...clcmVals);
  let ymaxL = Math.max(...clcmVals);
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

  let yminR = Math.min(...cdVals);
  let ymaxR = Math.max(...cdVals);
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
  if (yminR <= 0.0 && ymaxR >= 0.0) {
    const zeroLeft = (0.0 - yminL) / (ymaxL - yminL);
    const zeroRight = (0.0 - yminR) / (ymaxR - yminR);
    const span = ymaxR - yminR || 1.0;
    const shift = (zeroRight - zeroLeft) * span;
    yminR += shift;
    ymaxR += shift;
  }

  const xToPx = (x) => left + ((x - xmin) / (xmax - xmin)) * plotW;
  const yToPyLeft = (y) => top + (1.0 - (y - yminL) / (ymaxL - yminL)) * plotH;
  const yToPyRight = (y) => top + (1.0 - (y - yminR) / (ymaxR - yminR)) * plotH;

  alphaCtx.save();
  alphaCtx.strokeStyle = 'rgba(255, 255, 255, 0.35)';
  alphaCtx.lineWidth = 1.3;
  alphaCtx.beginPath();
  alphaCtx.moveTo(left, top);
  alphaCtx.lineTo(left, top + plotH);
  alphaCtx.lineTo(left + plotW, top + plotH);
  alphaCtx.stroke();

  alphaCtx.fillStyle = 'rgba(230, 236, 244, 0.85)';
  alphaCtx.font = '13px Consolas, "Courier New", monospace';
  drawSubLabel(alphaCtx, 'C', 'L', 12, 8);
  drawSubLabel(alphaCtx, 'C', 'M', 12, 22);
  alphaCtx.fillText('α', w - 22, h - 8);
  drawSubLabel(alphaCtx, 'C', 'D', w - 40, 8);

  const ticks = 5;
  alphaCtx.strokeStyle = 'rgba(255, 255, 255, 0.4)';
  for (let i = 0; i < ticks; i += 1) {
    const t = i / (ticks - 1);
    const xv = xmin + t * (xmax - xmin);
    const px = xToPx(xv);
    alphaCtx.beginPath();
    alphaCtx.moveTo(px, top + plotH);
    alphaCtx.lineTo(px, top + plotH + 8);
    alphaCtx.stroke();
    alphaCtx.fillText(xv.toFixed(1), px - 10, top + plotH + 22);
  }

  const niceStep = (range, target) => {
    const rough = range / target;
    const pow = 10 ** Math.floor(Math.log10(Math.max(rough, 1.0e-9)));
    const frac = rough / pow;
    let step = 1.0;
    if (frac <= 1.0) step = 1.0;
    else if (frac <= 2.0) step = 2.0;
    else if (frac <= 2.5) step = 2.5;
    else if (frac <= 5.0) step = 5.0;
    else step = 10.0;
    return step * pow;
  };
  const yStepL = niceStep(ymaxL - yminL, 6);
  const yStartL = Math.ceil(yminL / yStepL) * yStepL;
  for (let yv = yStartL; yv <= ymaxL + 1.0e-6; yv += yStepL) {
    const py = yToPyLeft(yv);
    alphaCtx.beginPath();
    alphaCtx.moveTo(left - 10, py);
    alphaCtx.lineTo(left, py);
    alphaCtx.stroke();
    alphaCtx.fillText(yv.toFixed(2), 6, py + 4);
  }

  const yStepR = niceStep(ymaxR - yminR, 6);
  const yStartR = Math.ceil(yminR / yStepR) * yStepR;
  for (let yv = yStartR; yv <= ymaxR + 1.0e-6; yv += yStepR) {
    const py = yToPyRight(yv);
    alphaCtx.beginPath();
    alphaCtx.moveTo(left + plotW, py);
    alphaCtx.lineTo(left + plotW + 10, py);
    alphaCtx.stroke();
    alphaCtx.fillText(yv.toFixed(3), left + plotW + 12, py + 4);
  }

  if (xmin <= 0.0 && xmax >= 0.0) {
    const px = xToPx(0.0);
    alphaCtx.strokeStyle = 'rgba(255, 255, 255, 0.25)';
    alphaCtx.beginPath();
    alphaCtx.moveTo(px, top);
    alphaCtx.lineTo(px, top + plotH);
    alphaCtx.stroke();
  }
  if (yminL <= 0.0 && ymaxL >= 0.0) {
    const py = yToPyLeft(0.0);
    alphaCtx.strokeStyle = 'rgba(255, 255, 255, 0.25)';
    alphaCtx.beginPath();
    alphaCtx.moveTo(left, py);
    alphaCtx.lineTo(left + plotW, py);
    alphaCtx.stroke();
  }

  const drawMarker = (shape, px, py, color, outlineOnly = false) => {
    alphaCtx.fillStyle = color;
    alphaCtx.strokeStyle = color;
    if (shape === 'square') {
      if (outlineOnly) {
        alphaCtx.lineWidth = 1.8;
        alphaCtx.strokeRect(px - 4, py - 4, 8, 8);
      } else {
        alphaCtx.fillRect(px - 4, py - 4, 8, 8);
      }
      return;
    }
    if (shape === 'triangle') {
      alphaCtx.beginPath();
      alphaCtx.moveTo(px, py - 5);
      alphaCtx.lineTo(px + 5, py + 4);
      alphaCtx.lineTo(px - 5, py + 4);
      alphaCtx.closePath();
      if (outlineOnly) {
        alphaCtx.lineWidth = 1.8;
        alphaCtx.stroke();
      } else {
        alphaCtx.fill();
      }
      return;
    }
    alphaCtx.beginPath();
    alphaCtx.arc(px, py, 4, 0, Math.PI * 2);
    if (outlineOnly) {
      alphaCtx.lineWidth = 1.8;
      alphaCtx.stroke();
    } else {
      alphaCtx.fill();
    }
  };

  runCases.forEach((caseItem) => {
    caseItem.history.forEach((p) => {
      const px = xToPx(p.alpha);
      drawMarker('circle', px, yToPyLeft(p.cl), caseItem.color);
      drawMarker('square', px, yToPyRight(p.cd), caseItem.color, true);
      drawMarker('triangle', px, yToPyLeft(p.cm), caseItem.color);
    });
  });

  const legendX = left + plotW - 70;
  const legendY = top + 14;
  const labelX = legendX + 12;
  alphaCtx.fillStyle = 'rgba(230, 236, 244, 0.85)';
  drawMarker('circle', legendX, legendY, alphaCtx.fillStyle);
  drawMarker('square', legendX, legendY + 18, alphaCtx.fillStyle, true);
  drawMarker('triangle', legendX, legendY + 36, alphaCtx.fillStyle);
  drawSubLabel(alphaCtx, 'C', 'L', labelX, legendY + 4);
  drawSubLabel(alphaCtx, 'C', 'D', labelX, legendY + 22);
  drawSubLabel(alphaCtx, 'C', 'M', labelX, legendY + 40);
  alphaCtx.restore();
}

// Polar plot: CL vs CD with nonnegative CD axis, point-only styling.
function drawPolarPlot() {
  if (!polarCanvas || !polarCtx) return;
  polarCtx.clearRect(0, 0, polarCanvas.width, polarCanvas.height);
  const allPoints = runCases.flatMap((c) => c.history);
  if (allPoints.length === 0) return;

  const w = polarCanvas.width;
  const h = polarCanvas.height;
  const left = 58;
  const right = 16;
  const top = 18;
  const bottom = 36;
  const plotW = w - left - right;
  const plotH = h - top - bottom;

  const points = allPoints.filter((p) => Number.isFinite(p.cd) && Number.isFinite(p.cl));
  if (points.length === 0) return;
  const cds = points.map((p) => Math.max(p.cd, 0.0));
  const cls = points.map((p) => p.cl);
  let xmin = Math.min(...cds, -0.005, 0.0);
  let xmax = Math.max(...cds, 0.0);
  let ymin = Math.min(...cls);
  let ymax = Math.max(...cls);
  if (!Number.isFinite(xmin) || !Number.isFinite(xmax)) {
    xmin = 0.0;
    xmax = 0.05;
  }
  if (!Number.isFinite(ymin) || !Number.isFinite(ymax)) {
    ymin = -1.0;
    ymax = 1.0;
  }
  const xRange = Math.max(xmax - xmin, 1.0e-4);
  const yRange = Math.max(ymax - ymin, 1.0e-4);
  const padX = 0.15 * xRange;
  const padY = 0.15 * yRange;
  xmin = 0.0;
  if (xmax < 0.001) {
    xmax = 0.01;
  } else {
    xmax += padX;
  }
  if (ymax - ymin < 1.0e-4) {
    ymin -= 0.1;
    ymax += 0.1;
  } else {
    ymin -= padY;
    ymax += padY;
  }

  const xToPx = (x) => left + ((x - xmin) / (xmax - xmin)) * plotW;
  const yToPy = (y) => top + (1.0 - (y - ymin) / (ymax - ymin)) * plotH;

  polarCtx.save();
  polarCtx.strokeStyle = 'rgba(255, 255, 255, 0.35)';
  polarCtx.lineWidth = 1.3;
  polarCtx.beginPath();
  polarCtx.moveTo(left, top);
  polarCtx.lineTo(left, top + plotH);
  polarCtx.lineTo(left + plotW, top + plotH);
  polarCtx.stroke();

  polarCtx.fillStyle = 'rgba(230, 236, 244, 0.85)';
  polarCtx.font = '13px Consolas, "Courier New", monospace';
  drawSubLabel(polarCtx, 'C', 'L', 12, 8);
  drawSubLabel(polarCtx, 'C', 'D', w - 28, h - 8);

  const ticks = 5;
  polarCtx.strokeStyle = 'rgba(255, 255, 255, 0.4)';
  for (let i = 0; i < ticks; i += 1) {
    const t = i / (ticks - 1);
    const xv = xmin + t * (xmax - xmin);
    const px = xToPx(xv);
    polarCtx.beginPath();
    polarCtx.moveTo(px, top + plotH);
    polarCtx.lineTo(px, top + plotH + 8);
    polarCtx.stroke();
    polarCtx.fillText(xv.toFixed(4), px - 18, top + plotH + 22);
  }

  for (let i = 0; i < ticks; i += 1) {
    const t = i / (ticks - 1);
    const yv = ymin + t * (ymax - ymin);
    const py = yToPy(yv);
    polarCtx.beginPath();
    polarCtx.moveTo(left - 10, py);
    polarCtx.lineTo(left, py);
    polarCtx.stroke();
    polarCtx.fillText(yv.toFixed(2), 6, py + 4);
  }

  runCases.forEach((caseItem) => {
    polarCtx.fillStyle = caseItem.color;
    caseItem.history.forEach((p) => {
      if (!Number.isFinite(p.cd) || !Number.isFinite(p.cl)) return;
      const px = xToPx(Math.max(p.cd, 0.0));
      const py = yToPy(p.cl);
      polarCtx.beginPath();
      polarCtx.arc(px, py, 4, 0, Math.PI * 2);
      polarCtx.fill();
    });
  });
  polarCtx.restore();
}

function rotatePoint(x, y, angle, ox, oy) {
  const dx = x - ox;
  const dy = y - oy;
  const ca = Math.cos(angle);
  const sa = Math.sin(angle);
  return {
    x: ox + dx * ca - dy * sa,
    y: oy + dx * sa + dy * ca,
  };
}

// Compute rotated bounds and scaling for plot framing (airfoil + margin).
function computeBounds(nb, angle) {
  let xmin = xb[0];
  let xmax = xb[0];
  let ymin = yb[0];
  let ymax = yb[0];

  for (let i = 1; i < nb; i += 1) {
    xmin = Math.min(xmin, xb[i]);
    xmax = Math.max(xmax, xb[i]);
    ymin = Math.min(ymin, yb[i]);
    ymax = Math.max(ymax, yb[i]);
  }

  const chord = xmax - xmin || 1.0;
  const ox = 0.5 * (xmin + xmax);
  const oy = 0.5 * (ymin + ymax);

  let rxmin = xmin;
  let rxmax = xmax;
  let rymin = ymin;
  let rymax = ymax;

  if (angle !== 0.0) {
    rxmin = Infinity;
    rxmax = -Infinity;
    rymin = Infinity;
    rymax = -Infinity;
    for (let i = 0; i < nb; i += 1) {
      const rp = rotatePoint(xb[i], yb[i], angle, ox, oy);
      rxmin = Math.min(rxmin, rp.x);
      rxmax = Math.max(rxmax, rp.x);
      rymin = Math.min(rymin, rp.y);
      rymax = Math.max(rymax, rp.y);
    }
  }

  const marginX = 0.05 * chord;
  const marginY = 0.2 * chord;
  rxmin -= marginX;
  rxmax += marginX;
  rymin -= marginY;
  rymax += marginY;

  const rect = canvas.getBoundingClientRect();
  const width = rect.width;
  const height = rect.height;
  const padding = 32;
  const spanX = rxmax - rxmin || 1.0;
  const spanY = rymax - rymin || 1.0;
  const scale = Math.min((width - padding * 2) / spanX, (height - padding * 2) / spanY);

  return {
    xmin: rxmin,
    xmax: rxmax,
    ymin: rymin,
    ymax: rymax,
    width,
    height,
    padding,
    scale,
    angle,
    ox,
    oy,
  };
}

function worldToCanvas(x, y, bounds) {
  let wx = x;
  let wy = y;
  if (bounds.angle !== 0.0) {
    const rotated = rotatePoint(x, y, bounds.angle, bounds.ox, bounds.oy);
    wx = rotated.x;
    wy = rotated.y;
  }
  return {
    x: bounds.padding + (wx - bounds.xmin) * bounds.scale,
    y: bounds.height - bounds.padding - (wy - bounds.ymin) * bounds.scale,
  };
}

function pointInPolygon(x, y, px, py, n) {
  let inside = false;
  for (let i = 0, j = n - 1; i < n; j = i, i += 1) {
    const xi = px[i];
    const yi = py[i];
    const xj = px[j];
    const yj = py[j];

    const intersect = ((yi > y) !== (yj > y))
      && (x < (xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi);
    if (intersect) {
      inside = !inside;
    }
  }
  return inside;
}

// Panel-method geometry setup: arc length, normals, influence matrices.
function buildPanelContext(nb, alphaRad) {
  const waklen = 1.0;
  const nw = Math.floor(nb / 12) + 10 * Math.floor(waklen);
  const total = nb + nw;
  if (!panelCtx || panelCtx.N !== nb || panelCtx.NW !== nw) {
    panelX = new Float64Array(total);
    panelY = new Float64Array(total);
    panelXP = new Float64Array(total);
    panelYP = new Float64Array(total);
    const n1 = nb + 1;
    panelCtx = {
      N: nb,
      NW: nw,
      WAKLEN: waklen,
      X: panelX,
      Y: panelY,
      XP: panelXP,
      YP: panelYP,
      S: new Float64Array(total),
      NX: new Float64Array(total),
      NY: new Float64Array(total),
      APANEL: new Float64Array(total),
      SHARP: true,
      PI: Math.PI,
      ANTE: 0.0,
      ASTE: 0.0,
      DSTE: 0.0,
      GAM: new Float64Array(n1),
      GAM_A: new Float64Array(n1),
      QINV: new Float64Array(total + 1),
      QINV_A: new Float64Array(total + 1),
      QVIS: new Float64Array(total + 1),
      GAMU: Array.from({ length: n1 }, () => new Float64Array(2)),
      SIG: new Float64Array(total),
      QF0: new Float64Array(nb),
      QF1: new Float64Array(nb),
      QF2: new Float64Array(nb),
      QF3: new Float64Array(nb),
      QINVU: Array.from({ length: total }, () => new Float64Array(2)),
      AIJ: createMatrix(n1, n1),
      BIJ: createMatrix(n1, total),
      AIJPIV: new Int32Array(n1),
      QOPI: 1.0 / (4.0 * Math.PI),
      HOPI: 1.0 / (2.0 * Math.PI),
      ALFA: 0.0,
      QINF: 1.0,
      LIMAGE: false,
      YIMAGE: 0.0,
      XTE: 0.0,
      YTE: 0.0,
      DZDG: new Float64Array(nb),
      DZDN: new Float64Array(nb),
      DQDG: new Float64Array(nb),
      DZDM: new Float64Array(total),
      DQDM: new Float64Array(total),
      LWAKE: false,
      LWDIJ: false,
      LADIJ: false,
      SNEW: new Float64Array(total),
    };
  }

  for (let i = 0; i < nb; i += 1) {
    panelX[i] = xb[i];
    panelY[i] = yb[i];
  }

  scalc(panelX, panelY, panelCtx.S, nb);
  segspl(panelX, panelXP, panelCtx.S, nb);
  segspl(panelY, panelYP, panelCtx.S, nb);
  ncalc(panelX, panelY, panelCtx.S, nb, panelCtx.NX, panelCtx.NY);
  panelCtx.XTE = 0.5 * (panelX[0] + panelX[nb - 1]);
  panelCtx.YTE = 0.5 * (panelY[0] + panelY[nb - 1]);
  tecalc(panelCtx);
  apcalc(panelCtx);
  panelCtx.ALFA = alphaRad;
  ggcalc(panelCtx);

  return panelCtx;
}

// Stagnation point from circulation sign change; used to split upper/lower.
function findStagnation(ctxPanel, nb) {
  const gam = ctxPanel.GAM;
  let ist = Math.floor(nb / 2) - 1;
  for (let i = 0; i < nb - 1; i += 1) {
    if (gam[i] >= 0.0 && gam[i + 1] < 0.0) {
      ist = i;
      break;
    }
  }

  const s = ctxPanel.S;
  const dgam = gam[ist + 1] - gam[ist];
  const ds = s[ist + 1] - s[ist];
  let sst;
  if (dgam !== 0.0) {
    if (gam[ist] < -gam[ist + 1]) {
      sst = s[ist] - ds * (gam[ist] / dgam);
    } else {
      sst = s[ist + 1] - ds * (gam[ist + 1] / dgam);
    }
  } else {
    sst = 0.5 * (s[ist] + s[ist + 1]);
  }

  if (sst <= s[ist]) sst = s[ist] + 1.0e-7;
  if (sst >= s[ist + 1]) sst = s[ist + 1] - 1.0e-7;

  const sstGo = dgam !== 0.0 ? (sst - s[ist + 1]) / dgam : 0.0;
  const sstGp = dgam !== 0.0 ? (s[ist] - sst) / dgam : 0.0;

  return { ist, sst, sstGo, sstGp };
}

function getSurfaceIndices(nb, ctxPanel) {
  const { ist } = findStagnation(ctxPanel, nb);
  const upperIdx = [];
  for (let i = ist; i >= 0; i -= 1) {
    upperIdx.push(i);
  }
  const lowerIdx = [];
  for (let i = ist + 1; i < nb; i += 1) {
    lowerIdx.push(i);
  }
  return { upperIdx, lowerIdx };
}

// Combine QINV basis with alpha rotation to get inviscid edge velocity.
function computeQinv(ctxPanel, alphaRad, nb) {
  const total = nb + (ctxPanel.NW ?? 0);
  const qinv = new Float64Array(total + 1);
  const qinvA = new Float64Array(total + 1);
  const cosA = Math.cos(alphaRad);
  const sinA = Math.sin(alphaRad);
  for (let i = 1; i <= total; i += 1) {
    const q0 = ctxPanel.QINVU[i - 1][0];
    const q90 = ctxPanel.QINVU[i - 1][1];
    const qi = cosA * q0 + sinA * q90;
    qinv[i] = Number.isFinite(qi) ? qi : 0.0;
    const qAi = -sinA * q0 + cosA * q90;
    qinvA[i] = Number.isFinite(qAi) ? qAi : 0.0;
  }
  if (ctxPanel.QINV) {
    ctxPanel.QINV.set(qinv);
  }
  if (ctxPanel.QINV_A) {
    ctxPanel.QINV_A.set(qinvA);
  }
  return { qinv, qinvA };
}

function getChordPoints(nb) {
  let leIdx = 0;
  let teIdx = 0;
  for (let i = 1; i < nb; i += 1) {
    if (xb[i] < xb[leIdx]) leIdx = i;
    if (xb[i] > xb[teIdx]) teIdx = i;
  }
  return {
    le: { x: xb[leIdx], y: yb[leIdx] },
    te: { x: xb[teIdx], y: yb[teIdx] },
  };
}

// Parse XFOIL-style .dat coordinate files (upper->lower or arbitrary order).
function parseDatAirfoil(text) {
  const lines = text.split(/\r?\n/).map((line) => line.trim()).filter((line) => line.length > 0);
  if (lines.length === 0) return null;
  let name = lines[0];
  let start = 1;
  if (/^[+-]?\d/.test(lines[0])) {
    name = 'Loaded Airfoil';
    start = 0;
  }
  const coords = [];
  for (let i = start; i < lines.length; i += 1) {
    const parts = lines[i].split(/[\s,]+/).filter(Boolean);
    if (parts.length < 2) continue;
    const x = parseFloat(parts[0]);
    const y = parseFloat(parts[1]);
    if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
    coords.push({ x, y });
  }
  if (coords.length < 5) return null;
  return { name, coords };
}

// Load custom airfoil coordinates into fixed buffers for panel setup.
function loadCustomAirfoil(data) {
  const count = Math.min(data.coords.length, xbCustom.length);
  for (let i = 0; i < count; i += 1) {
    xbCustom[i] = data.coords[i].x;
    ybCustom[i] = data.coords[i].y;
  }
  customAirfoil = { name: data.name, nb: count };
  currentAirfoilName = data.name;
}

// Map boundary-layer station indices to panel indices (upper/lower + wake).
function iblpan(blCtx, nb, nw = 0) {
  const ist = blCtx.IST;

  let ibl = 1;
  for (let i = ist; i >= 1; i -= 1) {
    ibl += 1;
    blCtx.IPAN[ibl][1] = i;
    blCtx.VTI[ibl][1] = 1.0;
  }
  blCtx.IBLTE[1] = ibl;
  blCtx.NBL[1] = ibl;

  ibl = 1;
  for (let i = ist + 1; i <= nb; i += 1) {
    ibl += 1;
    blCtx.IPAN[ibl][2] = i;
    blCtx.VTI[ibl][2] = -1.0;
  }
  blCtx.IBLTE[2] = ibl;
  for (let iw = 1; iw <= nw; iw += 1) {
    const i = nb + iw;
    const iblw = blCtx.IBLTE[2] + iw;
    blCtx.IPAN[iblw][2] = i;
    blCtx.VTI[iblw][2] = -1.0;
  }
  blCtx.NBL[2] = blCtx.IBLTE[2] + nw;

  for (let iw = 1; iw <= nw; iw += 1) {
    blCtx.IPAN[blCtx.IBLTE[1] + iw][1] = blCtx.IPAN[blCtx.IBLTE[2] + iw][2];
    blCtx.VTI[blCtx.IBLTE[1] + iw][1] = 1.0;
  }

  blCtx.LIPAN = true;
}

function xicalc(blCtx, ctxPanel) {
  const s = ctxPanel.S;
  const nb = ctxPanel.N;
  const nw = ctxPanel.NW ?? 0;
  const x = ctxPanel.X;
  const y = ctxPanel.Y;
  const xp = ctxPanel.XP;
  const yp = ctxPanel.YP;
  const xeps = 1.0e-7 * (s[nb - 1] - s[0]);
  const sAt = (i) => s[i - 1];

  blCtx.XSSI[1][1] = 0.0;
  for (let ibl = 2; ibl <= blCtx.IBLTE[1]; ibl += 1) {
    const i = blCtx.IPAN[ibl][1];
    blCtx.XSSI[ibl][1] = Math.max(blCtx.SST - sAt(i), xeps);
  }

  blCtx.XSSI[1][2] = 0.0;
  for (let ibl = 2; ibl <= blCtx.IBLTE[2]; ibl += 1) {
    const i = blCtx.IPAN[ibl][2];
    blCtx.XSSI[ibl][2] = Math.max(sAt(i) - blCtx.SST, xeps);
  }

  if (nw <= 0) {
    return;
  }

  const ibl1 = blCtx.IBLTE[1] + 1;
  blCtx.XSSI[ibl1][1] = blCtx.XSSI[ibl1 - 1][1];

  const ibl2 = blCtx.IBLTE[2] + 1;
  blCtx.XSSI[ibl2][2] = blCtx.XSSI[ibl2 - 1][2];

  for (let ibl = blCtx.IBLTE[2] + 2; ibl <= blCtx.NBL[2]; ibl += 1) {
    const i = blCtx.IPAN[ibl][2];
    const dxssi = Math.sqrt((x[i - 1] - x[i - 2]) ** 2 + (y[i - 1] - y[i - 2]) ** 2);

    const ibl1w = blCtx.IBLTE[1] + ibl - blCtx.IBLTE[2];
    const ibl2w = blCtx.IBLTE[2] + ibl - blCtx.IBLTE[2];
    blCtx.XSSI[ibl1w][1] = blCtx.XSSI[ibl1w - 1][1] + dxssi;
    blCtx.XSSI[ibl2w][2] = blCtx.XSSI[ibl2w - 1][2] + dxssi;
  }

  const telrat = 2.5;
  const crosp = (xp[0] * yp[nb - 1] - yp[0] * xp[nb - 1])
    / Math.sqrt((xp[0] ** 2 + yp[0] ** 2) * (xp[nb - 1] ** 2 + yp[nb - 1] ** 2));
  let dwdxte = crosp / Math.sqrt(Math.max(1.0 - crosp ** 2, 1.0e-12));
  dwdxte = Math.max(dwdxte, -3.0 / telrat);
  dwdxte = Math.min(dwdxte, 3.0 / telrat);

  const aa = 3.0 + telrat * dwdxte;
  const bb = -2.0 - telrat * dwdxte;

  if (ctxPanel.SHARP) {
    for (let iw = 1; iw <= nw; iw += 1) {
      blCtx.WGAP[iw] = 0.0;
    }
  } else {
    for (let iw = 1; iw <= nw; iw += 1) {
      const ibl = blCtx.IBLTE[2] + iw;
      const zn = 1.0 - (blCtx.XSSI[ibl][2] - blCtx.XSSI[blCtx.IBLTE[2]][2]) / (telrat * ctxPanel.ANTE);
      blCtx.WGAP[iw] = 0.0;
      if (zn >= 0.0) {
        blCtx.WGAP[iw] = ctxPanel.ANTE * (aa + bb * zn) * zn ** 2;
      }
    }
  }
}

// Build BL input arrays from inviscid edge velocities (XFOIL UICALC).
function uicalc(blCtx, qinv, qinvA) {
  for (let is = 1; is <= 2; is += 1) {
    blCtx.UINV[1][is] = 0.0;
    blCtx.UINV_A[1][is] = 0.0;
    for (let ibl = 2; ibl <= blCtx.NBL[is]; ibl += 1) {
      const i = blCtx.IPAN[ibl][is];
      blCtx.UINV[ibl][is] = blCtx.VTI[ibl][is] * qinv[i];
      blCtx.UINV_A[ibl][is] = blCtx.VTI[ibl][is] * qinvA[i];
    }
  }
}

// Compute viscous edge speed array from BL variables (XFOIL QVFUE).
function qvfue(blCtx, qvis) {
  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= blCtx.NBL[is]; ibl += 1) {
      const i = blCtx.IPAN[ibl][is];
      qvis[i] = blCtx.VTI[ibl][is] * blCtx.UEDG[ibl][is];
    }
  }
}

// Update panel circulation from viscous edge velocities (GAMQV).
function gamqv(ctxPanel, qvis, qinvA) {
  const nb = ctxPanel.N;
  for (let i = 1; i <= nb; i += 1) {
    ctxPanel.GAM[i - 1] = qvis[i];
    if (ctxPanel.GAM_A) {
      ctxPanel.GAM_A[i - 1] = qinvA[i];
    }
  }
  if (ctxPanel.QVIS) {
    ctxPanel.QVIS.set(qvis);
  }
}

function stmove(ctxPanel, blCtx, qinv, qinvA) {
  const nb = ctxPanel.N;
  const istOld = blCtx.IST;
  const { ist, sst, sstGo, sstGp } = findStagnation(ctxPanel, nb);
  const istNew = ist + 1;

  blCtx.IST = istNew;
  blCtx.SST = sst;
  blCtx.SST_GO = sstGo;
  blCtx.SST_GP = sstGp;

  if (istNew === istOld) {
    xicalc(blCtx, ctxPanel);
    return;
  }

  iblpan(blCtx, nb, ctxPanel.NW ?? 0);
  uicalc(blCtx, qinv, qinvA);
  xicalc(blCtx, ctxPanel);
  iblsys(blCtx);

  if (istNew > istOld) {
    const idif = istNew - istOld;
    blCtx.ITRAN[1] += idif;
    blCtx.ITRAN[2] -= idif;

    for (let ibl = blCtx.NBL[1]; ibl >= idif + 2; ibl -= 1) {
      blCtx.CTAU[ibl][1] = blCtx.CTAU[ibl - idif][1];
      blCtx.THET[ibl][1] = blCtx.THET[ibl - idif][1];
      blCtx.DSTR[ibl][1] = blCtx.DSTR[ibl - idif][1];
      blCtx.UEDG[ibl][1] = blCtx.UEDG[ibl - idif][1];
    }

    const dudx = blCtx.UEDG[idif + 2][1] / blCtx.XSSI[idif + 2][1];
    for (let ibl = idif + 1; ibl >= 2; ibl -= 1) {
      blCtx.CTAU[ibl][1] = blCtx.CTAU[idif + 2][1];
      blCtx.THET[ibl][1] = blCtx.THET[idif + 2][1];
      blCtx.DSTR[ibl][1] = blCtx.DSTR[idif + 2][1];
      blCtx.UEDG[ibl][1] = dudx * blCtx.XSSI[ibl][1];
    }

    for (let ibl = 2; ibl <= blCtx.NBL[2] - idif; ibl += 1) {
      blCtx.CTAU[ibl][2] = blCtx.CTAU[ibl + idif][2];
      blCtx.THET[ibl][2] = blCtx.THET[ibl + idif][2];
      blCtx.DSTR[ibl][2] = blCtx.DSTR[ibl + idif][2];
      blCtx.UEDG[ibl][2] = blCtx.UEDG[ibl + idif][2];
    }
  } else {
    const idif = istOld - istNew;
    blCtx.ITRAN[1] -= idif;
    blCtx.ITRAN[2] += idif;

    for (let ibl = blCtx.NBL[2]; ibl >= idif + 2; ibl -= 1) {
      blCtx.CTAU[ibl][2] = blCtx.CTAU[ibl - idif][2];
      blCtx.THET[ibl][2] = blCtx.THET[ibl - idif][2];
      blCtx.DSTR[ibl][2] = blCtx.DSTR[ibl - idif][2];
      blCtx.UEDG[ibl][2] = blCtx.UEDG[ibl - idif][2];
    }

    const dudx = blCtx.UEDG[idif + 2][2] / blCtx.XSSI[idif + 2][2];
    for (let ibl = idif + 1; ibl >= 2; ibl -= 1) {
      blCtx.CTAU[ibl][2] = blCtx.CTAU[idif + 2][2];
      blCtx.THET[ibl][2] = blCtx.THET[idif + 2][2];
      blCtx.DSTR[ibl][2] = blCtx.DSTR[idif + 2][2];
      blCtx.UEDG[ibl][2] = dudx * blCtx.XSSI[ibl][2];
    }

    for (let ibl = 2; ibl <= blCtx.NBL[1] - idif; ibl += 1) {
      blCtx.CTAU[ibl][1] = blCtx.CTAU[ibl + idif][1];
      blCtx.THET[ibl][1] = blCtx.THET[ibl + idif][1];
      blCtx.DSTR[ibl][1] = blCtx.DSTR[ibl + idif][1];
      blCtx.UEDG[ibl][1] = blCtx.UEDG[ibl + idif][1];
    }
  }

  const ueps = 1.0e-7;
  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= blCtx.NBL[is]; ibl += 1) {
      if (blCtx.UEDG[ibl][is] <= ueps) {
        blCtx.UEDG[ibl][is] = ueps;
      }
      blCtx.MASS[ibl][is] = blCtx.DSTR[ibl][is] * blCtx.UEDG[ibl][is];
    }
  }

  const ist0 = blCtx.IST - 1;
  const upperIdx = [];
  for (let i = ist0; i >= 0; i -= 1) {
    upperIdx.push(i);
  }
  const lowerIdx = [];
  for (let i = ist0 + 1; i < nb; i += 1) {
    lowerIdx.push(i);
  }
  blCtx.upperIdx = upperIdx;
  blCtx.lowerIdx = lowerIdx;
}

// Cp plot in XFOIL style: viscous Cp with optional inviscid overlay.
function drawCpPlot(nb, cpUpper, cpLower, lePt, tePt, bounds, cpInvAll = [], cpWake = []) {
  if (!cpCanvas) return;
  cpCtx.clearRect(0, 0, cpCanvas.width, cpCanvas.height);

  const dxChord = tePt.x - lePt.x;
  const dyChord = tePt.y - lePt.y;
  const chord2 = dxChord * dxChord + dyChord * dyChord || 1.0;

  const allCp = cpUpper.concat(cpLower).map((p) => p.cp);
  if (allCp.length === 0) return;
  const minCp = Math.min(...allCp);
  let cpMin = -1.0;
  if (minCp < cpMin) {
    cpMin = 0.5 * Math.floor(minCp / 0.5);
  }
  let cpMax = 1.0;

  const w = cpCanvas.width;
  const h = cpCanvas.height;
  const left = 64;
  const right = 18;
  const top = 18;
  const bottom = 36;
  const plotW = w - left - right;
  const plotH = h - top - bottom;

  let leScreen = null;
  let teScreen = null;
  if (bounds) {
    leScreen = worldToCanvas(lePt.x, lePt.y, bounds);
    teScreen = worldToCanvas(tePt.x, tePt.y, bounds);
  }
  const leX = leScreen ? leScreen.x * (w / canvas.width) : left;
  const teX = teScreen ? teScreen.x * (w / canvas.width) : left + plotW;
  const span = teX - leX || plotW;
  const xToPx = (s) => leX + s * span;
  const cpToPy = (cp) => top + ((cp - cpMin) / (cpMax - cpMin)) * plotH;
  const chordFrac = (x, y) => ((x - lePt.x) * dxChord + (y - lePt.y) * dyChord) / chord2;

  cpCtx.save();
  cpCtx.strokeStyle = 'rgba(255, 255, 255, 0.3)';
  cpCtx.lineWidth = 1.2;
  cpCtx.beginPath();
  cpCtx.rect(left, top, plotW, plotH);
  cpCtx.stroke();

  cpCtx.fillStyle = 'rgba(230, 236, 244, 0.85)';
  cpCtx.font = '13px Consolas, "Courier New", monospace';
  cpCtx.fillText('Cp', 12, 8);
  cpCtx.fillText('x/c', w - 40, h - 10);

  const zeroY = cpToPy(0.0);
  if (zeroY >= top && zeroY <= top + plotH) {
    cpCtx.strokeStyle = 'rgba(255, 255, 255, 0.25)';
    cpCtx.beginPath();
    cpCtx.moveTo(left, zeroY);
    cpCtx.lineTo(left + plotW, zeroY);
    cpCtx.stroke();
  }

  const tickCount = Math.round((cpMax - cpMin) / 0.5) + 1;
  cpCtx.strokeStyle = 'rgba(255, 255, 255, 0.35)';
  cpCtx.lineWidth = 1.4;
  cpCtx.beginPath();
  cpCtx.moveTo(left, top);
  cpCtx.lineTo(left, top + plotH);
  cpCtx.lineTo(left + plotW, top + plotH);
  cpCtx.stroke();

  cpCtx.strokeStyle = 'rgba(255, 255, 255, 0.4)';
  cpCtx.fillStyle = 'rgba(230, 236, 244, 0.85)';
  for (let i = 0; i < tickCount; i += 1) {
    const cpVal = cpMax - i * 0.5;
    const y = cpToPy(cpVal);
    cpCtx.beginPath();
    cpCtx.moveTo(left - 10, y);
    cpCtx.lineTo(left, y);
    cpCtx.stroke();
    cpCtx.fillText(cpVal.toFixed(1), 6, y + 4);
  }

  for (let i = 0; i < tickCount; i += 1) {
    const t = i / (tickCount - 1);
    const x = xToPx(t);
    cpCtx.beginPath();
    cpCtx.moveTo(x, top + plotH);
    cpCtx.lineTo(x, top + plotH + 8);
    cpCtx.stroke();
    cpCtx.fillText(t.toFixed(2), x - 10, top + plotH + 22);
  }

  cpCtx.strokeStyle = '#2f7bff';
  cpCtx.lineWidth = 1.6;
  cpCtx.beginPath();
  cpUpper.forEach((p, idx) => {
    const px = xToPx(chordFrac(p.x, p.y));
    const py = cpToPy(p.cp);
    if (idx === 0) {
      cpCtx.moveTo(px, py);
    } else {
      cpCtx.lineTo(px, py);
    }
  });
  cpCtx.stroke();

  cpCtx.strokeStyle = '#ff4a3d';
  cpCtx.beginPath();
  cpLower.forEach((p, idx) => {
    const px = xToPx(chordFrac(p.x, p.y));
    const py = cpToPy(p.cp);
    if (idx === 0) {
      cpCtx.moveTo(px, py);
    } else {
      cpCtx.lineTo(px, py);
    }
  });
  cpCtx.stroke();

  if (cpInvAll.length) {
    cpCtx.save();
    cpCtx.strokeStyle = 'rgba(255, 255, 255, 0.7)';
    cpCtx.lineWidth = 1.1;
    cpCtx.setLineDash([5, 6]);
    cpCtx.beginPath();
    cpInvAll.forEach((p, idx) => {
      const px = xToPx(chordFrac(p.x, p.y));
      const py = cpToPy(p.cp);
      if (idx === 0) {
        cpCtx.moveTo(px, py);
      } else {
        cpCtx.lineTo(px, py);
      }
    });
    cpCtx.stroke();
    cpCtx.restore();
  }

  if (cpWake.length) {
    cpCtx.save();
    cpCtx.strokeStyle = 'rgba(230, 236, 244, 0.55)';
    cpCtx.lineWidth = 1.2;
    cpCtx.beginPath();
    cpWake.forEach((p, idx) => {
      const px = xToPx(chordFrac(p.x, p.y));
      const py = cpToPy(p.cp);
      if (idx === 0) {
        cpCtx.moveTo(px, py);
      } else {
        cpCtx.lineTo(px, py);
      }
    });
    cpCtx.stroke();
    cpCtx.restore();
  }

  cpCtx.restore();
}

// Integrated force/moment coefficients, matching XFOIL's panel conventions.
function computeCoefficients(nb, ctxPanel, blCtx, alphaRad, qinvA, viscous) {
  const x = ctxPanel.X;
  const y = ctxPanel.Y;
  const gam = ctxPanel.GAM;
  const qinf = ctxPanel.QINF ?? 1.0;
  const minf = viscous ? blCtx?.MINF ?? 0.0 : 0.0;
  const xref = 0.25;
  const yref = 0.0;

  const sa = Math.sin(alphaRad);
  const ca = Math.cos(alphaRad);
  const beta = Math.sqrt(Math.max(1.0 - minf ** 2, 0.0));
  const betaMsq = beta === 0.0 ? 0.0 : -0.5 / beta;
  const bfac = 0.5 * minf ** 2 / (1.0 + beta);
  const bfacMsq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * betaMsq;

  let cl = 0.0;
  let cm = 0.0;
  let cdp = 0.0;
  let clAlf = 0.0;
  let clMsq = 0.0;

  let cginc = 1.0 - (gam[0] / qinf) ** 2;
  let cpg1 = cginc / (beta + bfac * cginc);
  let cpg1Msq = -cpg1 / (beta + bfac * cginc) * (betaMsq + bfacMsq * cginc);
  let cpiGam = -2.0 * gam[0] / qinf ** 2;
  let cpcCpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
  let cpg1Alf = cpcCpi * cpiGam * (qinvA?.[1] ?? 0.0);

  for (let i = 0; i < nb; i += 1) {
    const ip = i === nb - 1 ? 0 : i + 1;
    cginc = 1.0 - (gam[ip] / qinf) ** 2;
    const cpg2 = cginc / (beta + bfac * cginc);
    const cpg2Msq = -cpg2 / (beta + bfac * cginc) * (betaMsq + bfacMsq * cginc);
    cpiGam = -2.0 * gam[ip] / qinf ** 2;
    cpcCpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
    const cpg2Alf = cpcCpi * cpiGam * (qinvA?.[ip + 1] ?? 0.0);

    const dx = (x[ip] - x[i]) * ca + (y[ip] - y[i]) * sa;
    const dy = (y[ip] - y[i]) * ca - (x[ip] - x[i]) * sa;
    const dg = cpg2 - cpg1;

    const ax = (0.5 * (x[ip] + x[i]) - xref) * ca + (0.5 * (y[ip] + y[i]) - yref) * sa;
    const ay = (0.5 * (y[ip] + y[i]) - yref) * ca - (0.5 * (x[ip] + x[i]) - xref) * sa;
    const ag = 0.5 * (cpg2 + cpg1);

    const dxAlf = -(x[ip] - x[i]) * sa + (y[ip] - y[i]) * ca;
    const agAlf = 0.5 * (cpg2Alf + cpg1Alf);
    const agMsq = 0.5 * (cpg2Msq + cpg1Msq);

    cl += dx * ag;
    cdp -= dy * ag;
    cm -= dx * (ag * ax + dg * dx / 12.0) + dy * (ag * ay + dg * dy / 12.0);
    clAlf += dx * agAlf + ag * dxAlf;
    clMsq += dx * agMsq;

    cpg1 = cpg2;
    cpg1Alf = cpg2Alf;
    cpg1Msq = cpg2Msq;
  }

  let cd = 0.0;
  let cdf = 0.0;
  if (viscous && blCtx?.LBLINI) {
    const thwake = blCtx.THET[blCtx.NBL[2]][2];
    const urat = blCtx.UEDG[blCtx.NBL[2]][2] / qinf;
    const uewake = blCtx.UEDG[blCtx.NBL[2]][2]
      * (1.0 - blCtx.TKLAM) / (1.0 - blCtx.TKLAM * urat ** 2);
    const shwake = blCtx.DSTR[blCtx.NBL[2]][2] / blCtx.THET[blCtx.NBL[2]][2];
    cd = 2.0 * thwake * (uewake / qinf) ** (0.5 * (5.0 + shwake));

    for (let is = 1; is <= 2; is += 1) {
      for (let ibl = 3; ibl <= blCtx.IBLTE[is]; ibl += 1) {
        const i = blCtx.IPAN[ibl][is] - 1;
        const im = blCtx.IPAN[ibl - 1][is] - 1;
        const dx = (x[i] - x[im]) * ca + (y[i] - y[im]) * sa;
        cdf += 0.5 * (blCtx.TAU[ibl][is] + blCtx.TAU[ibl - 1][is]) * dx * 2.0 / qinf ** 2;
      }
    }
  }

  return { cl, cm, cd, cdf, cdp, clAlf, clMsq };
}

// Compressibility-corrected Cp (CPCALC), matching XFOIL's formulation.
function cpcalc(qvals, qinf, minf) {
  const n = qvals.length - 1;
  const cp = new Float64Array(n + 1);
  const beta = Math.sqrt(Math.max(1.0 - minf ** 2, 0.0));
  const bfac = 0.5 * minf ** 2 / (1.0 + beta);
  for (let i = 1; i <= n; i += 1) {
    const q = qvals[i] ?? 0.0;
    const cpinc = 1.0 - (q / qinf) ** 2;
    const den = beta + bfac * cpinc;
    cp[i] = den !== 0.0 ? cpinc / den : cpinc;
  }
  return cp;
}

// Allocate/initialize BL arrays and geometry mapping used by viscous solver.
function buildBlContext(nb, ctxPanel, alphaRad, ncr) {
  const nw = ctxPanel.NW ?? 0;
  const total = nb + nw;
  const nbl = total + 1;
  const ctxBl = {
    N: nb,
    NW: nw,
    IQX: nb,
    X: new Float64Array(total + 1),
    Y: new Float64Array(total + 1),
    S: new Float64Array(total + 1),
    NBL: new Int32Array(3),
    IBLTE: new Int32Array(3),
    ISYS: createMatrix1(nbl, 2),
    IPAN: createMatrix1(nbl, 2),
    VTI: createMatrix1(nbl, 2),
    XSSI: createMatrix1(nbl, 2),
    UEDG: createMatrix1(nbl, 2),
    UINV: createMatrix1(nbl, 2),
    UINV_A: createMatrix1(nbl, 2),
    THET: createMatrix1(nbl, 2),
    DSTR: createMatrix1(nbl, 2),
    CTAU: createMatrix1(nbl, 2),
    MASS: createMatrix1(nbl, 2),
    TAU: createMatrix1(nbl, 2),
    DIS: createMatrix1(nbl, 2),
    CTQ: createMatrix1(nbl, 2),
    DELT: createMatrix1(nbl, 2),
    TSTR: createMatrix1(nbl, 2),
    WGAP: new Float64Array(nw + 1),
    ACRIT: new Float64Array(3),
    XSTRIP: new Float64Array(3),
    ITRAN: new Int32Array(3),
    TFORCE: new Array(3).fill(false),
    XSSITR: new Float64Array(3),
    W1: new Float64Array(nb + 1),
    W2: new Float64Array(nb + 1),
    W3: new Float64Array(nb + 1),
    W4: new Float64Array(nb + 1),
    DTOR: Math.PI / 180.0,
    ANTE: 0.0,
    GAMBL: 1.4,
    GM1BL: 0.4,
    GAMM1: 0.4,
    MINF: 0.0,
    QINF: 1.0,
    QINFBL: 1.0,
    TKBL: 0.0,
    TKBL_MS: 0.0,
    RSTBL: 1.0,
    RSTBL_MS: 0.0,
    REYBL: 1.0e6,
    REYBL_RE: 1.0,
    REYBL_MS: 0.0,
    HSTINV: 0.0,
    HSTINV_MS: 0.0,
    HVRAT: 0.0,
    IDAMPV: 0,
    NCOM,
    blprv,
    blkin,
    blvar,
    blmid,
    trchek,
    tesys,
    blsys,
    hkin,
    syncComToVars,
    syncVarsToCom,
  };

  for (let i = 1; i <= total; i += 1) {
    ctxBl.X[i] = ctxPanel.X[i - 1];
    ctxBl.Y[i] = ctxPanel.Y[i - 1];
    ctxBl.S[i] = ctxPanel.S[i - 1];
  }

  let leIndex = 0;
  let minX = xb[0];
  for (let i = 1; i < nb; i += 1) {
    if (xb[i] < minX) {
      minX = xb[i];
      leIndex = i;
    }
  }
  const sLE = ctxBl.S[leIndex + 1];
  ctxBl.SLE = sLE;
  const { ist, sst, sstGo, sstGp } = findStagnation(ctxPanel, nb);
  ctxBl.SST = sst;
  ctxBl.SST_GO = sstGo;
  ctxBl.SST_GP = sstGp;
  ctxBl.IST = ist + 1;
  ctxBl.XLE = ctxBl.X[leIndex + 1];
  ctxBl.YLE = ctxBl.Y[leIndex + 1];
  ctxBl.XTE = 0.5 * (ctxBl.X[1] + ctxBl.X[nb]);
  ctxBl.YTE = 0.5 * (ctxBl.Y[1] + ctxBl.Y[nb]);

  const cosA = Math.cos(alphaRad);
  const sinA = Math.sin(alphaRad);
  const qinv = new Float64Array(total + 1);
  const qinvA = new Float64Array(total + 1);
  for (let i = 1; i <= total; i += 1) {
    const q0 = ctxPanel.QINVU[i - 1][0];
    const q90 = ctxPanel.QINVU[i - 1][1];
    const qi = cosA * q0 + sinA * q90;
    qinv[i] = Number.isFinite(qi) ? qi : 0.0;
    const qAi = -sinA * q0 + cosA * q90;
    qinvA[i] = Number.isFinite(qAi) ? qAi : 0.0;
  }

  ctxBl.IVX = nbl;

  const acrit = Number.isFinite(ncr) ? ncr : 9.0;
  ctxBl.ACRIT[1] = acrit;
  ctxBl.ACRIT[2] = acrit;
  ctxBl.XSTRIP[1] = 1.0;
  ctxBl.XSTRIP[2] = 1.0;

  const upperIdx = [];
  for (let i = ist; i >= 0; i -= 1) {
    upperIdx.push(i);
  }
  const lowerIdx = [];
  for (let i = ist + 1; i < nb; i += 1) {
    lowerIdx.push(i);
  }
  ctxBl.upperIdx = upperIdx;
  ctxBl.lowerIdx = lowerIdx;

  iblpan(ctxBl, nb, nw);
  xicalc(ctxBl, ctxPanel);

  for (let is = 1; is <= 2; is += 1) {
    let lastUe = 1.0e-6;
    for (let iblIdx = 2; iblIdx <= ctxBl.NBL[is]; iblIdx += 1) {
      const i = ctxBl.IPAN[iblIdx][is];
      const uinv = ctxBl.VTI[iblIdx][is] * qinv[i];
      const uinvA = ctxBl.VTI[iblIdx][is] * qinvA[i];
      ctxBl.UINV[iblIdx][is] = Number.isFinite(uinv) ? uinv : lastUe;
      ctxBl.UINV_A[iblIdx][is] = Number.isFinite(uinvA) ? uinvA : 0.0;
      const ue = Number.isFinite(uinv) && uinv > 0.0 ? uinv : lastUe;
      ctxBl.UEDG[iblIdx][is] = ue > 0.0 ? ue : 1.0e-6;
      lastUe = ctxBl.UEDG[iblIdx][is];
    }
  }

  ensureCtx(ctxBl);
  blpini(ctxBl);
  iblsys(ctxBl);

  ctxBl.XP = new Float64Array(nb + 1);
  ctxBl.YP = new Float64Array(nb + 1);
  for (let i = 1; i <= nb; i += 1) {
    ctxBl.XP[i] = ctxPanel.XP[i - 1];
    ctxBl.YP[i] = ctxPanel.YP[i - 1];
  }

  const nsys = ctxBl.NSYS;
  ctxBl.VA = createTensor3(3, 2, nsys);
  ctxBl.VB = createTensor3(3, 2, nsys);
  ctxBl.VDEL = createTensor3(3, 2, nsys);
  ctxBl.VM = createTensor3(3, nsys, nsys);
  ctxBl.VZ = createMatrix1(3, 2);
  ctxBl.VACCEL = 0.01;
  ctxBl.XOCTR = new Float64Array(3);
  ctxBl.YOCTR = new Float64Array(3);
  ctxBl.TINDEX = new Float64Array(3);
  ctxBl.CL = 0.0;
  ctxBl.CLSPEC = 0.0;
  ctxBl.LALFA = true;
  ctxBl.IDAMP = 0;
  ctxBl.LBLINI = false;
  ctxBl.GAMMA = 1.4;
  ctxBl.GAMM1 = 0.4;
  ctxBl.MINF1 = ctxBl.MINF;
  ctxBl.REINF1 = ctxBl.REYBL;
  ctxBl.MATYP = 1;
  ctxBl.RETYP = 1;
  ctxBl.MINF_CL = 0.0;
  ctxBl.REINF_CL = 0.0;
  ctxBl.DEBUG_BL = true;
  ctxBl.DEBUG_BL_FAILFAST = false;

  return ctxBl;
}

// Reconstruct viscous edge velocities from BL state when QVIS is absent.
function computeQvisFromUedg(blCtx, nb, qinv) {
  const qvis = new Float64Array(nb);
  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= blCtx.NBL[is]; ibl += 1) {
      const i = blCtx.IPAN[ibl][is];
      const qv = blCtx.VTI[ibl][is] * blCtx.UEDG[ibl][is];
      if (i >= 1 && i <= nb) qvis[i - 1] = qv;
    }
  }
  for (let i = 0; i < nb; i += 1) {
    if (!Number.isFinite(qvis[i]) || qvis[i] === 0.0) {
      qvis[i] = qinv[i + 1] ?? 0.0;
    }
  }
  return qvis;
}

// Initialize BL operating point (alpha, Mach, Re) to match XFOIL's OPER init.
function applyXfoilOperInit(blCtx, alphaRad, reinf1) {
  blCtx.LVISC = true;
  blCtx.LVCONV = false;
  blCtx.LWAKE = blCtx.NW > 0;
  blCtx.LALFA = true;
  blCtx.ALFA = alphaRad;
  blCtx.ADEG = alphaRad / blCtx.DTOR;
  blCtx.AWAKE = alphaRad;
  blCtx.AVISC = alphaRad;
  blCtx.MINF1 = blCtx.MINF;
  blCtx.REINF1 = reinf1;
  blCtx.REINF = reinf1;
  blCtx.MINF_CL = 0.0;
  blCtx.REINF_CL = 0.0;
}

// Reinitialize BL state for retry without changing the governing equations.
function resetBlState(blCtx) {
  blCtx.LBLINI = false;
  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 1; ibl <= blCtx.NBL[is]; ibl += 1) {
      blCtx.CTAU[ibl][is] = 0.0;
      blCtx.THET[ibl][is] = 0.0;
      blCtx.DSTR[ibl][is] = 0.0;
      blCtx.MASS[ibl][is] = 0.0;
      blCtx.TAU[ibl][is] = 0.0;
      blCtx.DIS[ibl][is] = 0.0;
      blCtx.CTQ[ibl][is] = 0.0;
      blCtx.DELT[ibl][is] = 0.0;
      blCtx.TSTR[ibl][is] = 0.0;
      const uinv = blCtx.UINV[ibl][is];
      const ue = Number.isFinite(uinv) && uinv > 0.0 ? uinv : 1.0e-6;
      blCtx.UEDG[ibl][is] = ue;
    }
  }
}

// Iterative viscous solve with optional re-init; mirrors XFOIL convergence logic.
function solveViscous(blCtx, ctxPanel, qinv, qinvA, maxIter = 20, maxRetry = 2) {
  const total = ctxPanel.N + (ctxPanel.NW ?? 0);
  const qvis = new Float64Array(total + 1);
  const eps1 = 1.0e-4;

  uicalc(blCtx, qinv, qinvA);
  let converged = false;
  for (let attempt = 0; attempt <= maxRetry; attempt += 1) {
    if (!blCtx.LBLINI) {
      for (let ibl = 1; ibl <= blCtx.NBL[1]; ibl += 1) {
        blCtx.UEDG[ibl][1] = blCtx.UINV[ibl][1];
      }
      for (let ibl = 1; ibl <= blCtx.NBL[2]; ibl += 1) {
        blCtx.UEDG[ibl][2] = blCtx.UINV[ibl][2];
      }
    }

    for (let iter = 0; iter < maxIter; iter += 1) {
      setbl(blCtx);
      const ok = blsolv(blCtx);
      if (!ok) break;
      blUpdate(blCtx);

      if (blCtx.LALFA) {
        mrcl(blCtx, blCtx.CL);
        comset(blCtx);
      } else {
        uicalc(blCtx, qinv, qinvA);
      }

      qvfue(blCtx, qvis);
      gamqv(ctxPanel, qvis, qinvA);
      stmove(ctxPanel, blCtx, qinv, qinvA);

      if (Number.isFinite(blCtx.RMSBL) && blCtx.RMSBL < eps1) {
        blCtx.LVCONV = true;
        blCtx.AVISC = blCtx.ALFA;
        blCtx.MVISC = blCtx.MINF;
        converged = true;
        break;
      }
    }

    if (converged) break;
    if (attempt < maxRetry) {
      resetBlState(blCtx);
    }
  }

  if (blCtx.DEBUG_BL) {
    for (let is = 1; is <= 2; is += 1) {
      const xt = blCtx.XSSITR[is];
    if (Number.isFinite(xt) && xt > 0.0) {
      const str = is === 1 ? blCtx.SST - xt : blCtx.SST + xt;
      const xtr = seval(str, ctxPanel.X, ctxPanel.XP, ctxPanel.S, ctxPanel.N);
      const ytr = seval(str, ctxPanel.Y, ctxPanel.YP, ctxPanel.S, ctxPanel.N);
      const chx = blCtx.XTE - blCtx.XLE;
      const chy = blCtx.YTE - blCtx.YLE;
      const chsq = chx * chx + chy * chy || 1.0;
      const xoctr = ((xtr - blCtx.XLE) * chx + (ytr - blCtx.YLE) * chy) / chsq;
      if (!Number.isFinite(xtr) || !Number.isFinite(ytr) || !Number.isFinite(xoctr)) {
        console.warn('BL: transition mapping non-finite', {
          is,
          ibl: blCtx.ITRAN[is],
          xt,
          str,
          sst: blCtx.SST,
          sle: blCtx.SLE,
          s1: ctxPanel.S[0],
          sn: ctxPanel.S[ctxPanel.N - 1],
          xle: blCtx.XLE,
          yle: blCtx.YLE,
          xte: blCtx.XTE,
          yte: blCtx.YTE,
          xtr,
          ytr,
          xoctr,
        });
      }
      console.warn(`BL: side ${is} transition x/c = ${Number.isFinite(xoctr) ? xoctr.toFixed(4) : 'NaN'} (ibl=${blCtx.ITRAN[is]})`);
    }
    }
  }

  return { qvis, converged };
}

function drawContours(bounds, grid, nx, ny, isoValues) {
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
}

// Streamline sketch from inviscid flow field, clipped to avoid airfoil interior.
function drawStreamlines(bounds, nb, ctxPanel) {
  const gridX = 60;
  const gridY = 32;
  const xSpan = bounds.xmax - bounds.xmin;
  const ySpan = bounds.ymax - bounds.ymin;
  const dx = xSpan / (gridX - 1);
  const dy = ySpan / (gridY - 1);
  const grid = new Float64Array(gridX * gridY);

  let psiMin = Infinity;
  let psiMax = -Infinity;

  for (let j = 0; j < gridY; j += 1) {
    for (let i = 0; i < gridX; i += 1) {
      const x = bounds.xmin + i * dx;
      const y = bounds.ymin + j * dy;
      const idx = j * gridX + i;

      if (pointInPolygon(x, y, panelX, panelY, nb)) {
        grid[idx] = NaN;
        continue;
      }

      const { psi } = psilin(nb, x, y, 0.0, 0.0, false, false, ctxPanel);
      grid[idx] = psi;
      if (psi < psiMin) psiMin = psi;
      if (psi > psiMax) psiMax = psi;
    }
  }

  if (!Number.isFinite(psiMin) || !Number.isFinite(psiMax)) {
    return;
  }

  const levels = 14;
  const isoValues = [];
  for (let k = 1; k <= levels; k += 1) {
    const frac = k / (levels + 1);
    isoValues.push(psiMin + frac * (psiMax - psiMin));
  }

  ctx.save();
  ctx.strokeStyle = 'rgba(220, 226, 232, 0.4)';
  ctx.lineWidth = 1.0;
  ctx.beginPath();
  ctx.rect(0, 0, canvas.width, canvas.height);
  const start = worldToCanvas(panelX[0], panelY[0], bounds);
  ctx.moveTo(start.x, start.y);
  for (let i = 1; i < nb; i += 1) {
    const p = worldToCanvas(panelX[i], panelY[i], bounds);
    ctx.lineTo(p.x, p.y);
  }
  ctx.closePath();
  ctx.clip('evenodd');
  drawContours(bounds, grid, gridX, gridY, isoValues);
  ctx.restore();
}

// Boundary-layer displacement thickness visualization (upper/lower + wake).
function drawBoundaryLayer(bounds, blCtx, ctxPanel, userScale) {
  ctx.save();
  const scale = Math.max(userScale, 0.05);

  ctx.lineWidth = 3.0;
  ctx.strokeStyle = 'rgba(47, 123, 255, 0.85)';
  ctx.beginPath();
  for (let ibl = 2; ibl <= blCtx.IBLTE[1]; ibl += 1) {
    const i = blCtx.IPAN[ibl][1];
    const dstr = blCtx.DSTR[ibl][1];
    if (!Number.isFinite(dstr)) continue;
    const dVis = dstr * scale;
    const x = ctxPanel.X[i - 1] + ctxPanel.NX[i - 1] * dVis;
    const y = ctxPanel.Y[i - 1] + ctxPanel.NY[i - 1] * dVis;
    const p = worldToCanvas(x, y, bounds);
    if (ibl === 2) {
      ctx.moveTo(p.x, p.y);
    } else {
      ctx.lineTo(p.x, p.y);
    }
  }
  ctx.stroke();

  ctx.strokeStyle = 'rgba(255, 74, 61, 0.85)';
  ctx.beginPath();
  for (let ibl = 2; ibl <= blCtx.IBLTE[2]; ibl += 1) {
    const i = blCtx.IPAN[ibl][2];
    const dstr = blCtx.DSTR[ibl][2];
    if (!Number.isFinite(dstr)) continue;
    const dVis = dstr * scale;
    const x = ctxPanel.X[i - 1] + ctxPanel.NX[i - 1] * dVis;
    const y = ctxPanel.Y[i - 1] + ctxPanel.NY[i - 1] * dVis;
    const p = worldToCanvas(x, y, bounds);
    if (ibl === 2) {
      ctx.moveTo(p.x, p.y);
    } else {
      ctx.lineTo(p.x, p.y);
    }
  }
  ctx.stroke();

  const nw = blCtx.NW ?? 0;
  if (nw > 0) {
    const is = 2;
    const dstrTe = blCtx.DSTR[blCtx.IBLTE[is] + 1][is];
    let dsf1 = 0.5;
    let dsf2 = 0.5;
    if (dstrTe !== 0.0) {
      dsf1 = (blCtx.DSTR[blCtx.IBLTE[1]][1] + 0.5 * ctxPanel.ANTE) / dstrTe;
      dsf2 = (blCtx.DSTR[blCtx.IBLTE[2]][2] + 0.5 * ctxPanel.ANTE) / dstrTe;
    }

    ctx.lineWidth = 2.4;
    ctx.strokeStyle = 'rgba(47, 123, 255, 0.6)';
    ctx.beginPath();
    let first = true;
    for (let ibl = blCtx.IBLTE[is] + 1; ibl <= blCtx.NBL[is]; ibl += 1) {
      const i = blCtx.IPAN[ibl][is];
      const dstr = blCtx.DSTR[ibl][is];
      if (!Number.isFinite(dstr)) continue;
      const dVis = dstr * dsf1 * scale;
      const x = ctxPanel.X[i - 1] - ctxPanel.NX[i - 1] * dVis;
      const y = ctxPanel.Y[i - 1] - ctxPanel.NY[i - 1] * dVis;
      const p = worldToCanvas(x, y, bounds);
      if (first) {
        ctx.moveTo(p.x, p.y);
        first = false;
      } else {
        ctx.lineTo(p.x, p.y);
      }
    }
    ctx.stroke();

    ctx.lineWidth = 2.4;
    ctx.strokeStyle = 'rgba(255, 74, 61, 0.6)';
    ctx.beginPath();
    first = true;
    for (let ibl = blCtx.IBLTE[is] + 1; ibl <= blCtx.NBL[is]; ibl += 1) {
      const i = blCtx.IPAN[ibl][is];
      const dstr = blCtx.DSTR[ibl][is];
      if (!Number.isFinite(dstr)) continue;
      const dVis = dstr * dsf2 * scale;
      const x = ctxPanel.X[i - 1] + ctxPanel.NX[i - 1] * dVis;
      const y = ctxPanel.Y[i - 1] + ctxPanel.NY[i - 1] * dVis;
      const p = worldToCanvas(x, y, bounds);
      if (first) {
        ctx.moveTo(p.x, p.y);
        first = false;
      } else {
        ctx.lineTo(p.x, p.y);
      }
    }
    ctx.stroke();
  }

  ctx.restore();
}

// Airfoil outline and chord line rendering, with custom airfoil handling.
function drawAirfoil(nb, bounds) {
  ctx.save();
  ctx.strokeStyle = '#ffffff';
  ctx.lineWidth = 1;
  ctx.setLineDash([4, 6]);
  const chordStart = worldToCanvas(0.0, 0.0, bounds);
  const chordEnd = worldToCanvas(1.0, 0.0, bounds);
  ctx.beginPath();
  ctx.moveTo(chordStart.x, chordStart.y);
  ctx.lineTo(chordEnd.x, chordEnd.y);
  ctx.stroke();
  ctx.restore();

  const source = sourceRadios.find((radio) => radio.checked)?.value || 'naca';
  const isCustom = source === 'custom' && customAirfoil;

  ctx.save();
  ctx.lineWidth = 2.5;
  ctx.setLineDash([]);
  if (isCustom || nb <= nside) {
    ctx.strokeStyle = '#ffffff';
    ctx.beginPath();
    for (let i = 0; i < nb; i += 1) {
      const p = worldToCanvas(xb[i], yb[i], bounds);
      if (i === 0) {
        ctx.moveTo(p.x, p.y);
      } else {
        ctx.lineTo(p.x, p.y);
      }
    }
    ctx.stroke();
    ctx.restore();
    return;
  }

  ctx.strokeStyle = '#ffffff';
  ctx.beginPath();
  for (let i = 0; i < nside; i += 1) {
    const p = worldToCanvas(xb[i], yb[i], bounds);
    if (i === 0) {
      ctx.moveTo(p.x, p.y);
    } else {
      ctx.lineTo(p.x, p.y);
    }
  }
  ctx.stroke();

  ctx.strokeStyle = '#ffffff';
  ctx.beginPath();
  for (let i = nside - 1; i < nb; i += 1) {
    const p = worldToCanvas(xb[i], yb[i], bounds);
    if (i === nside - 1) {
      ctx.moveTo(p.x, p.y);
    } else {
      ctx.lineTo(p.x, p.y);
    }
  }
  ctx.stroke();
  ctx.restore();
}

// Main compute/render pipeline: geometry -> panel -> BL -> coefficients -> plots.
function update() {
  const mode = seriesRadios.find((radio) => radio.checked)?.value || '4';
  const source = sourceRadios.find((radio) => radio.checked)?.value || 'naca';
  const alphaDeg = parseFloat(alphaSlider.value);
  const alphaRad = alphaDeg * (Math.PI / 180.0);
  const displayAngle = -alphaRad;
  updateAlphaLabel(alphaDeg);

  let nb = 0;
  if (source === 'custom') {
    if (!customAirfoil) {
      currentAirfoilName = 'Load .dat airfoil';
      return;
    }
    nb = customAirfoil.nb;
    xb.set(xbCustom);
    yb.set(ybCustom);
  } else if (mode === '4') {
    const m = parseInt(mSlider.value, 10);
    const p = parseInt(pSlider.value, 10);
    const t = parseInt(tSlider.value, 10);

    updateLabels(m, p, t);

    const ides = m * 1000 + p * 100 + t;
    ({ nb } = globalThis.Naca.naca4(ides, xx, yt, yc, nside, xb, yb));
  } else {
    if (mode === '5') {
      const series = series5Select.value;
      const t = parseInt(t5Slider.value, 10);
      const n5 = parseInt(series.charAt(0), 10);
      const n4 = parseInt(series.charAt(1), 10);
      const n3 = parseInt(series.charAt(2), 10);
      const ides = n5 * 10000 + n4 * 1000 + n3 * 100 + t;

      updateLabels5(series, t);

      const result = globalThis.Naca.naca5(ides, xx, yt, yc, nside, xb, yb);
      if (result.ok) {
        nb = result.nb;
      } else {
        return;
      }
    } else {
      const profile = series6Profile.value;
      const t = parseInt(t6Slider.value, 10);
      const cl = parseFloat(cl6Input.value);
      const camber = inferSixSeriesCamber(profile, cl);
      const fallbackA = defaultSixSeriesA(profile);

      updateLabels6(profile, t, camber, cl);

      const result = globalThis.Naca.naca6(
        {
          profile,
          toc: t / 100,
          camber,
          cl: Number.isFinite(cl) ? cl : 0.0,
          a: fallbackA,
        },
        xx,
        yt,
        yc,
        nside,
        xb,
        yb,
      );
      if (result.ok) {
        nb = result.nb;
        currentAirfoilName = result.name || currentAirfoilName;
      } else {
        return;
      }
    }
  }

  const bounds = computeBounds(nb, displayAngle);

  ctx.clearRect(0, 0, bounds.width, bounds.height);

  const ctxPanel = buildPanelContext(nb, alphaRad);
  if (viscousToggle.checked && ctxPanel) {
    xywake(ctxPanel);
    qwcalc(ctxPanel);
  }
  const { qinv, qinvA } = computeQinv(ctxPanel, alphaRad, nb);

  let blCtx = null;
  if (viscousToggle.checked && ctxPanel) {
    const dij = qdcalc(ctxPanel);
    const ncr = parseFloat(ncrInput.value);
    blCtx = buildBlContext(nb, ctxPanel, alphaRad, ncr);
    blCtx.DIJ = dij;
    const mach = parseFloat(machInput.value);
    const reinf = parseFloat(reynoldsInput.value);
    blCtx.MINF = Number.isFinite(mach) ? mach : 0.0;
    blCtx.MINF1 = blCtx.MINF;
    applyXfoilOperInit(blCtx, alphaRad, Number.isFinite(reinf) ? reinf : 1.0e6);

    const { qvis, converged } = solveViscous(blCtx, ctxPanel, qinv, qinvA, 20);
    for (let i = 0; i < nb; i += 1) {
      if (!Number.isFinite(ctxPanel.GAM[i])) {
        ctxPanel.GAM[i] = qvis[i + 1] ?? qinv[i + 1] ?? 0.0;
      }
    }
    blCtx.LVCONV = converged;

    let uinvMin = Infinity;
    let uinvMax = -Infinity;
    let uedgMin = Infinity;
    let uedgMax = -Infinity;
    let badUedg = 0;
    for (let is = 1; is <= 2; is += 1) {
      for (let ibl = 2; ibl <= blCtx.NBL[is]; ibl += 1) {
        const uinv = blCtx.UINV[ibl][is];
        const uedg = blCtx.UEDG[ibl][is];
        if (Number.isFinite(uinv)) {
          uinvMin = Math.min(uinvMin, uinv);
          uinvMax = Math.max(uinvMax, uinv);
        }
        if (Number.isFinite(uedg)) {
          uedgMin = Math.min(uedgMin, uedg);
          uedgMax = Math.max(uedgMax, uedg);
        } else {
          badUedg += 1;
        }
      }
    }
    console.warn('BL: Ue stats', {
      uinvMin,
      uinvMax,
      uedgMin,
      uedgMax,
      badUedg,
      rmsbl: blCtx.RMSBL,
      rmxbl: blCtx.RMXBL,
    });
  }

  if (ctxPanel) {
    drawStreamlines(bounds, nb, ctxPanel);
  }

  if (ctxPanel && qinv) {
    const qinf = ctxPanel.QINF ?? 1.0;
    const minf = viscousToggle.checked && blCtx ? blCtx.MINF ?? 0.0 : (parseFloat(machInput.value) || 0.0);
    const total = ctxPanel.N + (ctxPanel.NW ?? 0);
    const qvis = (ctxPanel.QVIS && ctxPanel.QVIS.length === total + 1)
      ? ctxPanel.QVIS
      : computeQvisFromUedg(blCtx, nb, qinv);
    const cpInv = cpcalc(qinv, qinf, minf);
    const cpVis = cpcalc(qvis, qinf, minf);

    let cpUpper = [];
    let cpLower = [];
    let cpWake = [];
    if (viscousToggle.checked && blCtx) {
      const ile1 = blCtx.IPAN[2][1] || 0;
      const ile2 = blCtx.IPAN[2][2] || 0;
      for (let i = 1; i <= ile1; i += 1) {
        cpUpper.push({ x: ctxPanel.X[i - 1], y: ctxPanel.Y[i - 1], cp: cpVis[i] });
      }
      for (let i = ile2; i <= nb; i += 1) {
        cpLower.push({ x: ctxPanel.X[i - 1], y: ctxPanel.Y[i - 1], cp: cpVis[i] });
      }
      for (let i = nb + 1; i <= total; i += 1) {
        cpWake.push({ x: ctxPanel.X[i - 1], y: ctxPanel.Y[i - 1], cp: cpVis[i] });
      }
    } else {
      const { upperIdx, lowerIdx } = getSurfaceIndices(nb, ctxPanel);
      for (let k = 0; k < upperIdx.length; k += 1) {
        const i = upperIdx[k] + 1;
        cpUpper.push({ x: ctxPanel.X[i - 1], y: ctxPanel.Y[i - 1], cp: cpInv[i] });
      }
      for (let k = 0; k < lowerIdx.length; k += 1) {
        const i = lowerIdx[k] + 1;
        cpLower.push({ x: ctxPanel.X[i - 1], y: ctxPanel.Y[i - 1], cp: cpInv[i] });
      }
    }

    const cpInvAll = [];
    if (viscousToggle.checked && blCtx) {
      for (let i = 1; i <= total; i += 1) {
        cpInvAll.push({ x: ctxPanel.X[i - 1], y: ctxPanel.Y[i - 1], cp: cpInv[i] });
      }
    }

    let lePt = { x: ctxPanel.XLE, y: ctxPanel.YLE };
    let tePt = { x: ctxPanel.XTE, y: ctxPanel.YTE };
    if (!Number.isFinite(lePt.x) || !Number.isFinite(lePt.y)
      || !Number.isFinite(tePt.x) || !Number.isFinite(tePt.y)) {
      const chordPts = getChordPoints(nb);
      lePt = chordPts.le;
      tePt = chordPts.te;
    }
    drawCpPlot(nb, cpUpper, cpLower, lePt, tePt, bounds, cpInvAll, cpWake);
  } else if (cpCtx) {
    cpCtx.clearRect(0, 0, cpCanvas.width, cpCanvas.height);
  }

  const coeffs = computeCoefficients(nb, ctxPanel, blCtx, alphaRad, qinvA, viscousToggle.checked);
  let displayCoeffs = {
    CL: coeffs.cl,
    CM: coeffs.cm,
    CD: coeffs.cd,
    CDF: coeffs.cdf,
    CDP: coeffs.cdp,
  };
  if (!viscousToggle.checked) {
    displayCoeffs.CD = coeffs.cdp;
  }
  if (blCtx) {
    blCtx.CL = coeffs.cl;
    blCtx.CM = coeffs.cm;
    blCtx.CD = coeffs.cd;
    blCtx.CDF = coeffs.cdf;
    blCtx.CDP = coeffs.cdp;
    displayCoeffs = blCtx;
  }
  const converged = !viscousToggle.checked || (blCtx ? blCtx.LVCONV === true : true);
  updateDataBox(alphaRad, displayCoeffs, converged);
  if (converged) {
    upsertSweepPoint(alphaDeg, displayCoeffs);
  }
  drawAlphaSweepPlot();
  drawPolarPlot();

  if (viscousToggle.checked && ctxPanel && blCtx) {
    let maxDelta = 0.0;
    for (let ibl = 2; ibl <= blCtx.NBL[1]; ibl += 1) {
      const d = blCtx.DELT[ibl][1];
      if (Number.isFinite(d)) maxDelta = Math.max(maxDelta, d);
    }
    for (let ibl = 2; ibl <= blCtx.NBL[2]; ibl += 1) {
      const d = blCtx.DELT[ibl][2];
      if (Number.isFinite(d)) maxDelta = Math.max(maxDelta, d);
    }
    if (maxDelta <= 0.0) {
      console.warn('BL: DELT sample', {
        del1: blCtx.DELT[2][1],
        del2: blCtx.DELT[2][2],
        th1: blCtx.THET[2][1],
        th2: blCtx.THET[2][2],
        ds1: blCtx.DSTR[2][1],
        ds2: blCtx.DSTR[2][2],
        ue1: blCtx.UEDG[2][1],
        ue2: blCtx.UEDG[2][2],
        hk2: blCtx.HK2,
        hs2: blCtx.HS2,
        cf2: blCtx.CF2,
      });
    }
  }

  drawAirfoil(nb, bounds);
  if (viscous && blCtx) {
    drawBoundaryLayer(bounds, blCtx, ctxPanel, 1.0);
  }
}

[mSlider, pSlider, tSlider, t5Slider].forEach((slider) => {
  slider.addEventListener('input', update);
});
alphaSlider.addEventListener('input', update);
machInput.addEventListener('input', update);
reynoldsInput.addEventListener('input', update);
ncrInput.addEventListener('input', update);
if (alphaMinus) {
  alphaMinus.addEventListener('click', () => {
    const current = parseFloat(alphaSlider.value);
    setAlphaValue(current - 1.0);
  });
}
if (alphaPlus) {
  alphaPlus.addEventListener('click', () => {
    const current = parseFloat(alphaSlider.value);
    setAlphaValue(current + 1.0);
  });
}
if (alphaSweep) {
  alphaSweep.addEventListener('click', sweepAlpha);
}
if (addRunCaseButton) {
  addRunCaseButton.addEventListener('click', () => {
    const newCase = createRunCase();
    runCases.push(newCase);
    activeCaseId = newCase.id;
    renderRunCases();
    drawAlphaSweepPlot();
    drawPolarPlot();
  });
}

seriesRadios.forEach((radio) => {
  radio.addEventListener('change', () => {
    const mode = seriesRadios.find((item) => item.checked)?.value || '4';
    const source = sourceRadios.find((item) => item.checked)?.value || 'naca';
    if (customAirfoil && source === 'naca') {
      customAirfoil = null;
    }
    controls4.hidden = mode !== '4';
    controls5.hidden = mode !== '5';
    controls6.hidden = mode !== '6';
    update();
  });
});

sourceRadios.forEach((radio) => {
  radio.addEventListener('change', () => {
    const source = sourceRadios.find((item) => item.checked)?.value || 'naca';
    const isCustom = source === 'custom';
    nacaOptions.hidden = isCustom;
    customOptions.hidden = !isCustom;
    controls4.hidden = isCustom || (seriesRadios.find((item) => item.checked)?.value || '4') !== '4';
    controls5.hidden = isCustom || (seriesRadios.find((item) => item.checked)?.value || '4') !== '5';
    controls6.hidden = isCustom || (seriesRadios.find((item) => item.checked)?.value || '4') !== '6';
    update();
  });
});

viscousToggle.addEventListener('change', update);
series5Select.addEventListener('change', update);
series6Profile.addEventListener('change', () => {
  update();
});
t6Slider.addEventListener('input', update);
cl6Input.addEventListener('input', update);

loadDatButton.addEventListener('click', () => {
  datFileInput.click();
});

datFileInput.addEventListener('change', () => {
  const file = datFileInput.files?.[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = () => {
    const data = parseDatAirfoil(String(reader.result || ''));
    if (!data) {
      console.warn('DAT load failed: no coordinates parsed.');
      return;
    }
    loadCustomAirfoil(data);
    const customRadio = sourceRadios.find((radio) => radio.value === 'custom');
    if (customRadio) {
      customRadio.checked = true;
      nacaOptions.hidden = true;
      customOptions.hidden = false;
      controls4.hidden = true;
      controls5.hidden = true;
    }
    update();
  };
  reader.readAsText(file);
});

window.addEventListener('resize', () => {
  resizeCanvas();
  update();
});

const initSource = sourceRadios.find((item) => item.checked)?.value || 'naca';
nacaOptions.hidden = initSource === 'custom';
customOptions.hidden = initSource !== 'custom';
controls4.hidden = initSource === 'custom' || (seriesRadios.find((item) => item.checked)?.value || '4') !== '4';
controls5.hidden = initSource === 'custom' || (seriesRadios.find((item) => item.checked)?.value || '4') !== '5';
controls6.hidden = initSource === 'custom' || (seriesRadios.find((item) => item.checked)?.value || '4') !== '6';
updateAlphaLabel(parseFloat(alphaSlider.value));
if (runCases.length === 0) {
  const initialCase = createRunCase();
  runCases.push(initialCase);
  activeCaseId = initialCase.id;
}
renderRunCases();
resizeCanvas();
update();
