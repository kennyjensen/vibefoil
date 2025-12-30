import { scalc, segspl } from './spline.js';
import {
  apcalc,
  ncalc,
  psilin,
  ggcalc,
  stfind,
} from './xpanel.js';
import { computeCoefficients, cpcalc, tecalc, pangen } from './xfoil.js';
import { buildBlContext, computeQvisFromUedg, specal, viscal } from './xoper.js';
import { createMatrix } from './arrays.js';

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
const databaseOptions = document.getElementById('databaseOptions');
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
const dataBox = document.getElementById('dataBox');
const DEFAULT_BL_ITER = 20;
const UIUC_BASE_URL = 'https://raw.githubusercontent.com/kennyjensen/uiuc_airfoil_database/master/';
const UIUC_TREE_URL = 'https://api.github.com/repos/kennyjensen/uiuc_airfoil_database/git/trees/master?recursive=1';
const UIUC_FALLBACK_LIST = [
  'e387.dat',
  'e423.dat',
  'fx63-137.dat',
  'fx63-137sm.dat',
  'naca2412.dat',
  'naca4412.dat',
  'naca0012.dat',
  'naca0015.dat',
  'naca6409.dat',
  'sd7003.dat',
  'sd7037.dat',
  's1223.dat',
];
let uiucListLoaded = false;
let uiucListLoading = false;
let lastSolverPayload = null;
let lastViewportSize = {
  width: window.innerWidth,
  height: window.innerHeight,
};

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
const nIterInput = document.getElementById('nIter');
const advancedModeToggle = document.getElementById('advancedMode');
const advancedControls = document.getElementById('advancedControls');
const reuseSolutionToggle = document.getElementById('reuseSolution');
const loadDatButton = document.getElementById('loadDat');
const datFileInput = document.getElementById('datFile');
const uiucNameInput = document.getElementById('uiucName');
const uiucList = document.getElementById('uiucList');
const fetchUiucButton = document.getElementById('fetchUiuc');
const uiucStatus = document.getElementById('uiucStatus');
const flapXInput = document.getElementById('flapX');
const flapYInput = document.getElementById('flapY');
const flapDefInput = document.getElementById('flapDef');
const flapDefValue = document.getElementById('flapDefValue');
const flapXValue = document.getElementById('flapXValue');
const flapYValue = document.getElementById('flapYValue');
const flapDefMinus = document.getElementById('flapDefMinus');
const flapDefPlus = document.getElementById('flapDefPlus');
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
const nside = 123;
const xx = new Float64Array(nside);
const yt = new Float64Array(nside);
const yc = new Float64Array(nside);
const xb = new Float64Array(2 * nside);
const yb = new Float64Array(2 * nside);
const xbBuffer = new Float64Array(2 * nside);
const ybBuffer = new Float64Array(2 * nside);
const xbCustom = new Float64Array(2 * nside);
const ybCustom = new Float64Array(2 * nside);
let customAirfoil = null;
let currentAirfoilName = 'NACA 2412';
let customAirfoilVersion = 0;

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
const panelCache = { ctx: null, key: null };
const blCache = { ctx: null, key: null };
let reuseState = null;
let solverWorker = null;
let solverRequestId = 0;
let latestSolverId = 0;
let solverInFlight = false;
const pendingSolves = new Map();

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

function updateFlapDefLabel(value) {
  if (!flapDefValue) return;
  const num = Number.isFinite(value) ? value : 0.0;
  flapDefValue.textContent = `${num.toFixed(1)}°`;
}

function updateFlapHingeLabels(xValue, yValue) {
  if (flapXValue) {
    const num = Number.isFinite(xValue) ? xValue : 0.0;
    flapXValue.textContent = num.toFixed(2);
  }
  if (flapYValue) {
    const num = Number.isFinite(yValue) ? yValue : 0.0;
    flapYValue.textContent = num.toFixed(2);
  }
}

function adjustFlapDeflection(delta) {
  if (!flapDefInput) return;
  const min = parseFloat(flapDefInput.min);
  const max = parseFloat(flapDefInput.max);
  const current = parseFloat(flapDefInput.value);
  const next = Number.isFinite(current) ? current + delta : delta;
  const clamped = Math.max(min, Math.min(max, next));
  flapDefInput.value = `${clamped}`;
  update();
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
  return update();
}

// Run-case structure for multi-sweep overlays and color separation.
function createRunCase() {
  const colors = ['#ff4a3d', '#ff9f1a', '#ffd166', '#6dd36f', '#2f7bff', '#9b5de5'];
  const color = colors[(nextCaseId - 1) % colors.length];
  const id = nextCaseId;
  nextCaseId += 1;
  const foilLabel = currentAirfoilName || 'Airfoil';
  return {
    id,
    name: `${foilLabel} Run ${id}`,
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
    await setAlphaValue(a);
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

  updateDataBoxScale();
}

function updateDataBoxScale() {
  if (!dataBox || !canvas) return;
  const rect = canvas.getBoundingClientRect();
  if (!rect.width) return;
  const targetWidth = rect.width * 0.25;
  const width = targetWidth;
  const fontSize = width / 9.2;
  const padY = width * 0.06;
  const padX = width * 0.07;
  dataBox.style.setProperty('--data-box-width', `${width.toFixed(1)}px`);
  dataBox.style.setProperty('--data-box-font', `${fontSize.toFixed(1)}px`);
  dataBox.style.setProperty('--data-box-pad-y', `${padY.toFixed(1)}px`);
  dataBox.style.setProperty('--data-box-pad-x', `${padX.toFixed(1)}px`);
}

function reframeBoundsForCanvas(bounds) {
  if (!bounds) return bounds;
  const rect = canvas.getBoundingClientRect();
  if (!rect.width || !rect.height) return bounds;
  const width = rect.width;
  const height = rect.height;
  const minDim = Math.min(width, height);
  const padding = Math.max(8, Math.floor(minDim * 0.025));
  const spanX = bounds.xmax - bounds.xmin || 1.0;
  const spanY = bounds.ymax - bounds.ymin || 1.0;
  const plotWidth = Math.max(1, width - padding * 2);
  const plotHeight = Math.max(1, height - padding * 2);
  const scale = Math.min(plotWidth / spanX, plotHeight / spanY);
  const offsetX = padding + Math.max(0, (plotWidth - spanX * scale) * 0.5);
  const offsetY = padding + Math.max(0, (plotHeight - spanY * scale) * 0.5);
  return {
    ...bounds,
    width,
    height,
    padding,
    offsetX,
    offsetY,
    scale,
  };
}

// Insert or replace sweep point (alpha-indexed) to keep plots ordered.
function upsertSweepPoint(alphaDeg, coeffs, converged = true) {
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
  const point = {
    alpha: alphaDeg,
    cl: coeffs.CL,
    cd: coeffs.CD,
    cm: coeffs.CM,
    converged,
  };
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
  const alphaRect = alphaCanvas.getBoundingClientRect();
  alphaCtx.clearRect(0, 0, alphaRect.width, alphaRect.height);
  const allPoints = runCases.flatMap((c) => c.history);
  if (allPoints.length === 0) return;

  const w = alphaRect.width;
  const h = alphaRect.height;
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

  const staleColor = 'rgba(154, 160, 168, 0.85)';
  runCases.forEach((caseItem) => {
    caseItem.history.forEach((p) => {
      const px = xToPx(p.alpha);
      const color = p.converged === false ? staleColor : caseItem.color;
      drawMarker('circle', px, yToPyLeft(p.cl), color);
      drawMarker('square', px, yToPyRight(p.cd), color, true);
      drawMarker('triangle', px, yToPyLeft(p.cm), color);
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
  const polarRect = polarCanvas.getBoundingClientRect();
  polarCtx.clearRect(0, 0, polarRect.width, polarRect.height);
  const allPoints = runCases.flatMap((c) => c.history);
  if (allPoints.length === 0) return;

  const w = polarRect.width;
  const h = polarRect.height;
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

  const staleColor = 'rgba(154, 160, 168, 0.85)';
  runCases.forEach((caseItem) => {
    caseItem.history.forEach((p) => {
      if (!Number.isFinite(p.cd) || !Number.isFinite(p.cl)) return;
      const px = xToPx(Math.max(p.cd, 0.0));
      const py = yToPy(p.cl);
      polarCtx.fillStyle = p.converged === false ? staleColor : caseItem.color;
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
  const rect = canvas.getBoundingClientRect();
  const width = rect.width;
  const height = rect.height;
  const minDim = Math.min(width, height);
  const padding = Math.max(8, Math.floor(minDim * 0.025));
  const marginY = (minDim < 520 ? 0.04 : 0.12) * chord;
  rxmin -= marginX;
  rxmax += marginX;
  rymin -= marginY;
  rymax += marginY;

  const spanX = rxmax - rxmin || 1.0;
  let spanY = rymax - rymin || 1.0;
  const plotWidth = Math.max(1, width - padding * 2);
  const plotHeight = Math.max(1, height - padding * 2);
  const targetSpanY = spanX * (plotHeight / plotWidth);
  if (spanY < targetSpanY) {
    const expand = targetSpanY - spanY;
    rymin -= 0.5 * expand;
    rymax += 0.5 * expand;
    spanY = targetSpanY;
  }
  const scale = Math.min(plotWidth / spanX, plotHeight / spanY);
  const offsetX = padding + Math.max(0, (plotWidth - spanX * scale) * 0.5);
  const offsetY = padding + Math.max(0, (plotHeight - spanY * scale) * 0.5);

  return {
    xmin: rxmin,
    xmax: rxmax,
    ymin: rymin,
    ymax: rymax,
    width,
    height,
    padding,
    offsetX,
    offsetY,
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
    x: bounds.offsetX + (wx - bounds.xmin) * bounds.scale,
    y: bounds.height - bounds.offsetY - (wy - bounds.ymin) * bounds.scale,
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
function buildPanelContext(nb, alphaRad, opts = {}) {
  const { reusePanel = false, geometryKey = '' } = opts;
  const waklen = 1.0;
  const nw = Math.floor(nb / 12) + 10 * Math.floor(waklen);
  const total = nb + nw;
  if (reusePanel && panelCache.ctx && panelCache.key === geometryKey
    && panelCache.ctx.N === nb && panelCache.ctx.NW === nw) {
    panelCtx = panelCache.ctx;
    panelCtx.ALFA = alphaRad;
    return panelCtx;
  }

  const resetViscous = !reusePanel || panelCache.key !== geometryKey;
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

  if (resetViscous && panelCtx.QVIS) {
    panelCtx.QVIS.fill(0.0);
    panelCtx.LWAKE = false;
    panelCtx.LWDIJ = false;
    panelCtx.LADIJ = false;
  }

  panelCache.ctx = panelCtx;
  panelCache.key = geometryKey;
  return panelCtx;
}

// Stagnation point from circulation sign change; used to split upper/lower.
function getSurfaceIndices(nb, ctxPanel) {
  const { ist } = stfind(ctxPanel, nb);
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
  customAirfoilVersion += 1;
  currentAirfoilName = data.name;
}

function setUiucStatus(message, isError = false) {
  if (!uiucStatus) return;
  uiucStatus.textContent = message;
  uiucStatus.hidden = !message;
  uiucStatus.classList.toggle('error', isError);
}

function normalizeUiucFilename(rawName) {
  const trimmed = String(rawName || '').trim();
  if (!trimmed) return null;
  const base = trimmed.split('/').pop() || '';
  if (!base) return null;
  return base.toLowerCase().endsWith('.dat') ? base : `${base}.dat`;
}

function renderUiucList(items) {
  if (!uiucList) return;
  uiucList.innerHTML = '';
  items.forEach((name) => {
    const option = document.createElement('option');
    option.value = name;
    uiucList.appendChild(option);
  });
}

function parseUiucTree(payload) {
  const tree = Array.isArray(payload?.tree) ? payload.tree : [];
  const names = new Set();
  tree.forEach((item) => {
    if (!item || item.type !== 'blob' || typeof item.path !== 'string') return;
    if (!item.path.toLowerCase().endsWith('.dat')) return;
    const filename = item.path.split('/').pop();
    if (filename) {
      names.add(filename);
    }
  });
  return Array.from(names).sort((a, b) => a.localeCompare(b));
}

async function ensureUiucListLoaded() {
  if (uiucListLoaded || uiucListLoading || !uiucList) return;
  uiucListLoading = true;
  if (!uiucListLoaded) {
    renderUiucList(UIUC_FALLBACK_LIST);
  }
  try {
    const response = await fetch(UIUC_TREE_URL);
    if (!response.ok) {
      throw new Error(`HTTP ${response.status}`);
    }
    const payload = await response.json();
    const names = parseUiucTree(payload);
    if (names.length > 0) {
      renderUiucList(names);
    }
    uiucListLoaded = true;
  } catch (error) {
    uiucListLoaded = true;
    setUiucStatus('Autocomplete using starter list (GitHub API fetch failed).');
  } finally {
    uiucListLoading = false;
  }
}

async function fetchUiucAirfoil() {
  if (!fetchUiucButton || !uiucNameInput) return;
  const filename = normalizeUiucFilename(uiucNameInput.value);
  if (!filename) {
    setUiucStatus('Enter a .dat filename to fetch.', true);
    return;
  }

  const url = `${UIUC_BASE_URL}${encodeURIComponent(filename)}`;
  fetchUiucButton.disabled = true;
  setUiucStatus(`Fetching ${filename}...`);

  try {
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error(`HTTP ${response.status}`);
    }
    const text = await response.text();
    const data = parseDatAirfoil(text);
    if (!data) {
      throw new Error('No coordinates parsed.');
    }
    if (!data.name || data.name === 'Loaded Airfoil') {
      data.name = filename.replace(/\.dat$/i, '');
    }
    loadCustomAirfoil(data);
    const dbRadio = sourceRadios.find((radio) => radio.value === 'database');
    if (dbRadio) {
      dbRadio.checked = true;
    }
    nacaOptions.hidden = true;
    customOptions.hidden = true;
    if (databaseOptions) {
      databaseOptions.hidden = false;
    }
    controls4.hidden = true;
    controls5.hidden = true;
    controls6.hidden = true;
    setUiucStatus(`Loaded ${data.name}.`);
    update();
  } catch (error) {
    setUiucStatus(`Fetch failed: ${error.message}`, true);
  } finally {
    fetchUiucButton.disabled = false;
  }
}


// Cp plot in XFOIL style: viscous Cp with optional inviscid overlay.
function drawCpPlot(nb, cpUpper, cpLower, lePt, tePt, bounds, cpInvAll = [], cpWake = []) {
  if (!cpCanvas) return;
  const cpRect = cpCanvas.getBoundingClientRect();
  cpCtx.clearRect(0, 0, cpRect.width, cpRect.height);

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

  const w = cpRect.width;
  const h = cpRect.height;
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
  const mainRect = canvas.getBoundingClientRect();
  const scaleX = mainRect.width ? w / mainRect.width : 1.0;
  const leX = leScreen ? leScreen.x * scaleX : left;
  const teX = teScreen ? teScreen.x * scaleX : left + plotW;
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
function drawStreamlines(bounds, nb, streamlines) {
  if (!streamlines || !streamlines.grid) return;
  const { grid, gridX, gridY, psiMin, psiMax } = streamlines;
  if (!Number.isFinite(psiMin) || !Number.isFinite(psiMax)) return;

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
  ctx.rect(0, 0, bounds.width, bounds.height);
  const start = worldToCanvas(xb[0], yb[0], bounds);
  ctx.moveTo(start.x, start.y);
  for (let i = 1; i < nb; i += 1) {
    const p = worldToCanvas(xb[i], yb[i], bounds);
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

function drawBoundaryLayerLines(bounds, lines) {
  if (!lines) return;
  const { upper = [], lower = [], wakeUpper = [], wakeLower = [] } = lines;
  const drawLine = (pts, strokeStyle, lineWidth) => {
    if (!pts.length) return;
    ctx.save();
    ctx.lineWidth = lineWidth;
    ctx.strokeStyle = strokeStyle;
    ctx.beginPath();
    pts.forEach((pt, idx) => {
      const p = worldToCanvas(pt.x, pt.y, bounds);
      if (idx === 0) {
        ctx.moveTo(p.x, p.y);
      } else {
        ctx.lineTo(p.x, p.y);
      }
    });
    ctx.stroke();
    ctx.restore();
  };

  drawLine(upper, 'rgba(47, 123, 255, 0.85)', 3.0);
  drawLine(lower, 'rgba(255, 74, 61, 0.85)', 3.0);
  drawLine(wakeUpper, 'rgba(47, 123, 255, 0.6)', 2.4);
  drawLine(wakeLower, 'rgba(255, 74, 61, 0.6)', 2.4);
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
  const isCustom = (source === 'custom' || source === 'database') && customAirfoil;

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

function drawFlapHinge(bounds, hinge) {
  if (!hinge || !Number.isFinite(hinge.x) || !Number.isFinite(hinge.y)) return;
  const p = worldToCanvas(hinge.x, hinge.y, bounds);
  const size = 6;

  ctx.save();
  ctx.strokeStyle = '#9aa0a8';
  ctx.lineWidth = 2.0;
  ctx.beginPath();
  ctx.moveTo(p.x - size, p.y);
  ctx.lineTo(p.x + size, p.y);
  ctx.moveTo(p.x, p.y - size);
  ctx.lineTo(p.x, p.y + size);
  ctx.stroke();

  ctx.fillStyle = 'rgba(154, 160, 168, 0.55)';
  ctx.beginPath();
  ctx.arc(p.x, p.y, 3.5, 0, Math.PI * 2);
  ctx.fill();
  ctx.restore();
}

function applySolverResult(payload) {
  if (!payload?.ok) return;
  lastSolverPayload = payload;
  const {
    nb,
    xb: xbResult,
    yb: ybResult,
    bounds,
    streamlines,
    cpData,
    coeffs,
    converged,
    blLines,
    hinge,
    airfoilName,
    alphaDeg,
    alphaRad,
    viscous,
  } = payload;

  if (!nb || !bounds) return;
  const framedBounds = reframeBoundsForCanvas(bounds);
  if (xbResult && ybResult) {
    xb.set(xbResult);
    yb.set(ybResult);
  }
  if (airfoilName) {
    currentAirfoilName = airfoilName;
  }

  ctx.clearRect(0, 0, framedBounds.width, framedBounds.height);

  if (streamlines) {
    drawStreamlines(framedBounds, nb, streamlines);
  }

  if (cpData?.upper?.length || cpData?.lower?.length) {
    drawCpPlot(
      nb,
      cpData.upper,
      cpData.lower,
      cpData.le,
      cpData.te,
      framedBounds,
      cpData.invAll || [],
      cpData.wake || [],
    );
  } else if (cpCtx) {
    const cpRect = cpCanvas.getBoundingClientRect();
    cpCtx.clearRect(0, 0, cpRect.width, cpRect.height);
  }

  updateDataBox(alphaRad, coeffs, converged);
  upsertSweepPoint(alphaDeg, coeffs, converged);
  drawAlphaSweepPlot();
  drawPolarPlot();

  drawAirfoil(nb, framedBounds);
  drawFlapHinge(framedBounds, hinge);
  if (viscous && blLines) {
    drawBoundaryLayerLines(framedBounds, blLines);
  }
}

function resolvePendingSolve(id, payload) {
  const resolve = pendingSolves.get(id);
  if (resolve) {
    resolve(payload);
    pendingSolves.delete(id);
  }
}

function cancelPendingSolves() {
  pendingSolves.forEach((resolve) => resolve({ canceled: true }));
  pendingSolves.clear();
}

function spawnSolverWorker() {
  solverWorker = new Worker(new URL('./solver_worker.js', import.meta.url), { type: 'module' });
  solverWorker.onmessage = (event) => {
    const payload = event.data;
    resolvePendingSolve(payload.requestId, payload);
    if (payload.requestId !== latestSolverId) return;
    solverInFlight = false;
    if (!payload.ok) {
      if (payload.error) {
        console.warn(`Solver worker error: ${payload.error}`);
      }
      if (payload.errorStack) {
        console.warn(payload.errorStack);
      }
      return;
    }
    if (payload.reuseState) {
      reuseState = payload.reuseState;
    }
    applySolverResult(payload);
  };
  solverWorker.onerror = (event) => {
    solverInFlight = false;
    cancelPendingSolves();
    console.warn('Solver worker error:', event);
  };
  solverWorker.onmessageerror = (event) => {
    solverInFlight = false;
    cancelPendingSolves();
    console.warn('Solver worker message error:', event);
  };
}

function syncAdvancedControls() {
  if (advancedControls && advancedModeToggle) {
    advancedControls.hidden = !advancedModeToggle.checked;
  }
}

// Main compute/render pipeline: geometry -> panel -> BL -> coefficients -> plots.
function update() {
  const mode = seriesRadios.find((radio) => radio.checked)?.value || '4';
  const source = sourceRadios.find((radio) => radio.checked)?.value || 'naca';
  const alphaDeg = parseFloat(alphaSlider.value);
  const alphaRad = alphaDeg * (Math.PI / 180.0);
  updateAlphaLabel(alphaDeg);
  updateFlapDefLabel(parseFloat(flapDefInput?.value));
  updateFlapHingeLabels(
    parseFloat(flapXInput?.value),
    parseFloat(flapYInput?.value),
  );

  if (source === 'custom' || source === 'database') {
    if (!customAirfoil) {
      currentAirfoilName = source === 'database' ? 'Select database airfoil' : 'Load .dat airfoil';
      return Promise.resolve({ skipped: true });
    }
  } else if (mode === '4') {
    const m = parseInt(mSlider.value, 10);
    const p = parseInt(pSlider.value, 10);
    const t = parseInt(tSlider.value, 10);

    updateLabels(m, p, t);
  } else {
    if (mode === '5') {
      const series = series5Select.value;
      const t = parseInt(t5Slider.value, 10);

      updateLabels5(series, t);
    } else {
      const profile = series6Profile.value;
      const t = parseInt(t6Slider.value, 10);
      const cl = parseFloat(cl6Input.value);
      const camber = inferSixSeriesCamber(profile, cl);
      const fallbackA = defaultSixSeriesA(profile);

      updateLabels6(profile, t, camber, cl);
    }
  }

  let geometryKey = '';
  if (source === 'custom' || source === 'database') {
    geometryKey = `${source}:${customAirfoil?.name || 'none'}:${customAirfoilVersion}`;
  } else if (mode === '4') {
    geometryKey = `naca4:${mSlider.value}:${pSlider.value}:${tSlider.value}`;
  } else if (mode === '5') {
    geometryKey = `naca5:${series5Select.value}:${t5Slider.value}`;
  } else {
    geometryKey = `naca6:${series6Profile.value}:${t6Slider.value}:${cl6Input.value}`;
  }
  const flapKey = `${flapXInput?.value ?? '0.75'}:${flapYInput?.value ?? '0.5'}:${flapDefInput?.value ?? '0'}`;
  geometryKey = `${geometryKey}:flap:${flapKey}`;

  const profile = series6Profile.value;
  const cl = parseFloat(cl6Input.value);
  const camber = inferSixSeriesCamber(profile, cl);
  const fallbackA = defaultSixSeriesA(profile);

  const rect = canvas.getBoundingClientRect();
  const advancedMode = !!advancedModeToggle?.checked;
  const reuseEnabled = viscousToggle.checked && !!reuseSolutionToggle?.checked;
  const reuseKeyedState = reuseState?.geometryKey === geometryKey
    ? reuseState
    : null;
  const settings = {
    mode,
    source,
    m: parseInt(mSlider.value, 10),
    p: parseInt(pSlider.value, 10),
    t: parseInt(tSlider.value, 10),
    series5: series5Select.value,
    t5: parseInt(t5Slider.value, 10),
    profile6: profile,
    t6: parseInt(t6Slider.value, 10),
    cl6: Number.isFinite(cl) ? cl : 0.0,
    camber6: camber,
    fallbackA6: fallbackA,
    flap: {
      x: parseFloat(flapXInput?.value),
      yrel: parseFloat(flapYInput?.value),
      deflection: parseFloat(flapDefInput?.value),
    },
    custom: customAirfoil
      ? {
        name: customAirfoil.name,
        nb: customAirfoil.nb,
        x: xbCustom,
        y: ybCustom,
      }
      : null,
    alphaDeg,
    alphaRad,
    geometryKey,
    advancedMode,
    reusePanel: viscousToggle.checked && (advancedMode ? reuseEnabled : true),
    reuseSolution: advancedMode ? reuseEnabled : false,
    reuseState: reuseKeyedState,
    viscous: viscousToggle.checked,
    mach: parseFloat(machInput.value),
    reynolds: parseFloat(reynoldsInput.value),
    ncr: parseFloat(ncrInput.value),
    nIter: parseInt(nIterInput.value, 10),
    canvasWidth: rect.width,
    canvasHeight: rect.height,
  };

  if (!solverWorker) {
    spawnSolverWorker();
  }
  if (solverInFlight && solverWorker) {
    solverWorker.terminate();
    cancelPendingSolves();
    spawnSolverWorker();
    solverInFlight = false;
  }

  solverRequestId += 1;
  latestSolverId = solverRequestId;
  solverInFlight = true;
  solverWorker.postMessage({ requestId: solverRequestId, settings });
  return new Promise((resolve) => {
    pendingSolves.set(solverRequestId, resolve);
  });
}

syncAdvancedControls();

[mSlider, pSlider, tSlider, t5Slider].forEach((slider) => {
  slider.addEventListener('input', update);
});
alphaSlider.addEventListener('input', update);
machInput.addEventListener('input', update);
reynoldsInput.addEventListener('input', update);
ncrInput.addEventListener('input', update);
nIterInput.addEventListener('input', update);
[flapXInput, flapYInput, flapDefInput].forEach((input) => {
  if (input) {
    input.addEventListener('input', update);
  }
});
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
if (flapDefMinus) {
  flapDefMinus.addEventListener('click', () => {
    adjustFlapDeflection(-1.0);
  });
}
if (flapDefPlus) {
  flapDefPlus.addEventListener('click', () => {
    adjustFlapDeflection(1.0);
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
    const isDatabase = source === 'database';
    nacaOptions.hidden = source !== 'naca';
    customOptions.hidden = !isCustom;
    if (databaseOptions) {
      databaseOptions.hidden = !isDatabase;
    }
    controls4.hidden = source !== 'naca' || (seriesRadios.find((item) => item.checked)?.value || '4') !== '4';
    controls5.hidden = source !== 'naca' || (seriesRadios.find((item) => item.checked)?.value || '4') !== '5';
    controls6.hidden = source !== 'naca' || (seriesRadios.find((item) => item.checked)?.value || '4') !== '6';
    update();
  });
});

viscousToggle.addEventListener('change', update);
if (reuseSolutionToggle) {
  reuseSolutionToggle.addEventListener('change', update);
}
if (advancedModeToggle) {
  advancedModeToggle.addEventListener('change', () => {
    syncAdvancedControls();
    update();
  });
}
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
      if (databaseOptions) {
        databaseOptions.hidden = true;
      }
      controls4.hidden = true;
      controls5.hidden = true;
    }
    update();
  };
  reader.readAsText(file);
});

if (fetchUiucButton) {
  fetchUiucButton.addEventListener('click', () => {
    fetchUiucAirfoil();
  });
}

if (uiucNameInput) {
  uiucNameInput.addEventListener('keydown', (event) => {
    if (event.key === 'Enter') {
      event.preventDefault();
      fetchUiucAirfoil();
    }
  });
  uiucNameInput.addEventListener('input', () => {
    ensureUiucListLoaded();
  });
  uiucNameInput.addEventListener('focus', () => {
    ensureUiucListLoaded();
  });
}

const layoutPanel = document.querySelector('.layout');
const pageButtons = Array.from(document.querySelectorAll('.page-indicator button'));

function setActivePageButton(index) {
  if (pageButtons.length === 0) return;
  const clamped = Math.max(0, Math.min(index, pageButtons.length - 1));
  pageButtons.forEach((button, idx) => {
    button.classList.toggle('active', idx === clamped);
  });
}

function updatePageIndicator() {
  if (!layoutPanel) return;
  const pageWidth = layoutPanel.clientWidth;
  if (!pageWidth) return;
  const index = Math.round(layoutPanel.scrollLeft / pageWidth);
  setActivePageButton(index);
}

if (layoutPanel && pageButtons.length > 0) {
  let scrollFrame = null;
  layoutPanel.addEventListener('scroll', () => {
    if (scrollFrame) return;
    scrollFrame = window.requestAnimationFrame(() => {
      scrollFrame = null;
      updatePageIndicator();
    });
  }, { passive: true });

  pageButtons.forEach((button) => {
    button.addEventListener('click', () => {
      const index = Number.parseInt(button.dataset.index || '0', 10);
      const pageWidth = layoutPanel.clientWidth;
      layoutPanel.scrollTo({
        left: pageWidth * index,
        behavior: 'smooth',
      });
      setActivePageButton(index);
    });
  });
}

window.addEventListener('resize', () => {
  const nextWidth = window.innerWidth;
  const nextHeight = window.innerHeight;
  const widthChanged = Math.abs(nextWidth - lastViewportSize.width) > 1;
  const heightChanged = Math.abs(nextHeight - lastViewportSize.height) > 1;
  if (!widthChanged && !heightChanged) return;
  lastViewportSize = { width: nextWidth, height: nextHeight };
  if (!widthChanged) {
    return;
  }
  resizeCanvas();
  if (lastSolverPayload) {
    applySolverResult(lastSolverPayload);
  }
  updatePageIndicator();
});

const initSource = sourceRadios.find((item) => item.checked)?.value || 'naca';
nacaOptions.hidden = initSource !== 'naca';
customOptions.hidden = initSource !== 'custom';
if (databaseOptions) {
  databaseOptions.hidden = initSource !== 'database';
}
controls4.hidden = initSource !== 'naca' || (seriesRadios.find((item) => item.checked)?.value || '4') !== '4';
controls5.hidden = initSource !== 'naca' || (seriesRadios.find((item) => item.checked)?.value || '4') !== '5';
controls6.hidden = initSource !== 'naca' || (seriesRadios.find((item) => item.checked)?.value || '4') !== '6';
updateAlphaLabel(parseFloat(alphaSlider.value));
updateFlapDefLabel(parseFloat(flapDefInput?.value));
updateFlapHingeLabels(
  parseFloat(flapXInput?.value),
  parseFloat(flapYInput?.value),
);
if (runCases.length === 0) {
  const initialCase = createRunCase();
  runCases.push(initialCase);
  activeCaseId = initialCase.id;
}
renderRunCases();
resizeCanvas();
update();
