import { buildBlContext, computeQvisFromUedg, specal, viscal } from './xoper.js';

// High-level orchestrator for the XFOIL port: UI, geometry generation,
// inviscid panel solve, viscous BL coupling, and plotting.
const canvas = document.getElementById('plot');
const ctx = canvas.getContext('2d');
const cpCanvas = document.getElementById('cpPlot');
const cpCtx = cpCanvas.getContext('2d');
const downloadCpButton = document.getElementById('downloadCp');
const downloadBlButton = document.getElementById('downloadBl');
const editAirfoilButton = document.getElementById('editAirfoil');
const editCpButton = document.getElementById('editCp');
const airfoilFrame = document.getElementById('airfoilFrame');
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
const hoverBox = document.getElementById('hoverBox');
const hoverName = document.getElementById('hoverName');
const hoverMach = document.getElementById('hoverMach');
const hoverRe = document.getElementById('hoverRe');
const hoverNcr = document.getElementById('hoverNcr');
const hoverAlpha = document.getElementById('hoverAlpha');
const hoverFlap = document.getElementById('hoverFlap');
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
let lastGeometrySettings = null;
let pendingGeometrySettings = null;
let lastCpPlot = null;
let cpHoverS = null;
let alphaPlotState = null;
let polarPlotState = null;
let alphaHover = null;
let polarHover = null;
let activeHoverSource = null;
let editMode = false;
let editDragIndex = null;
let editHoverIndex = null;
let cpEditMode = false;
let cpEditDragIndex = null;
let cpEditHoverIndex = null;
let cpEditPoints = null;
let cpPanActive = false;
let cpPanStart = null;
let airfoilTouchTimer = null;
let cpTouchTimer = null;
let touchAirfoilDrag = false;
let touchCpDrag = false;
let lastAirfoilBounds = null;
let lastAirfoilNb = 0;
let airfoilZoom = 1.0;
let cpZoom = 1.0;
let alphaZoom = 1.0;
let polarZoom = 1.0;
let cpZoomCenter = null;
let alphaZoomCenter = null;
let polarZoomCenter = null;
let airfoilZoomCenterScreen = null;
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
const runCases = [];
let activeCaseId = null;
let nextCaseId = 1;
let sweeping = false;
const blCache = { ctx: null, key: null };
let reuseState = null;
let solverWorker = null;
let solverRequestId = 0;
let latestSolverId = 0;
let solverInFlight = false;
const pendingSolves = new Map();
const pendingDumps = new Map();
const pendingQdes = new Map();
let nextDumpId = 1;
let nextQdesId = 1;

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

function drawHoverLabel(ctx2d, px, py, lines, bounds) {
  if (!lines.length) return;
  const pad = 6;
  const lineHeight = 14;
  ctx2d.save();
  const baseFont = '12px Consolas, "Courier New", monospace';
  const subFont = '10px Consolas, "Courier New", monospace';
  ctx2d.font = baseFont;
  const measureLine = (line) => {
    if (typeof line === 'string') {
      return ctx2d.measureText(line).width;
    }
    ctx2d.font = baseFont;
    const baseWidth = ctx2d.measureText(line.base).width;
    ctx2d.font = subFont;
    const subWidth = ctx2d.measureText(line.sub).width;
    ctx2d.font = baseFont;
    const valueWidth = ctx2d.measureText(` ${line.value}`).width;
    return baseWidth + 1 + subWidth + valueWidth;
  };
  const widths = lines.map((line) => measureLine(line));
  const boxW = Math.max(...widths, 40) + pad * 2;
  const boxH = lineHeight * lines.length + pad * 2;
  let bx = px + 12;
  let by = py - boxH - 12;
  if (bx + boxW > bounds.w - 6) bx = px - boxW - 12;
  if (bx < 6) bx = 6;
  if (by < 6) by = py + 12;
  if (by + boxH > bounds.h - 6) by = bounds.h - boxH - 6;
  ctx2d.fillStyle = 'rgba(8, 10, 12, 0.82)';
  ctx2d.strokeStyle = 'rgba(230, 236, 244, 0.5)';
  ctx2d.lineWidth = 1.0;
  ctx2d.beginPath();
  ctx2d.rect(bx, by, boxW, boxH);
  ctx2d.fill();
  ctx2d.stroke();
  ctx2d.fillStyle = 'rgba(230, 236, 244, 0.92)';
  lines.forEach((line, idx) => {
    const textY = by + pad + lineHeight * (idx + 1) - 3;
    if (typeof line === 'string') {
      ctx2d.font = baseFont;
      ctx2d.fillText(line, bx + pad, textY);
      return;
    }
    ctx2d.font = baseFont;
    const baseX = bx + pad;
    ctx2d.fillText(line.base, baseX, textY);
    const baseWidth = ctx2d.measureText(line.base).width;
    ctx2d.font = subFont;
    ctx2d.fillText(line.sub, baseX + baseWidth + 1, textY + 3);
    const subWidth = ctx2d.measureText(line.sub).width;
    ctx2d.font = baseFont;
    ctx2d.fillText(` ${line.value}`, baseX + baseWidth + 1 + subWidth, textY);
  });
  ctx2d.restore();
}

function formatHoverValue(value, digits = 2, suffix = '') {
  if (!Number.isFinite(value)) return '—';
  return `${value.toFixed(digits)}${suffix}`;
}

function canvasToWorld(x, y, bounds) {
  let sx = x;
  let sy = y;
  if (bounds.zoom && bounds.zoomCenterScreen) {
    sx = bounds.zoomCenterScreen.x + (sx - bounds.zoomCenterScreen.x) / bounds.zoom;
    sy = bounds.zoomCenterScreen.y + (sy - bounds.zoomCenterScreen.y) / bounds.zoom;
  }
  const wx = (sx - bounds.offsetX) / bounds.scale + bounds.xmin;
  const wy = (bounds.height - bounds.offsetY - sy) / bounds.scale + bounds.ymin;
  if (bounds.angle === 0.0) {
    return { x: wx, y: wy };
  }
  return rotatePoint(wx, wy, -bounds.angle, bounds.ox, bounds.oy);
}

function canvasToRotatedWorld(x, y, bounds) {
  let sx = x;
  let sy = y;
  if (bounds.zoom && bounds.zoomCenterScreen) {
    sx = bounds.zoomCenterScreen.x + (sx - bounds.zoomCenterScreen.x) / bounds.zoom;
    sy = bounds.zoomCenterScreen.y + (sy - bounds.zoomCenterScreen.y) / bounds.zoom;
  }
  const wx = (sx - bounds.offsetX) / bounds.scale + bounds.xmin;
  const wy = (bounds.height - bounds.offsetY - sy) / bounds.scale + bounds.ymin;
  return { x: wx, y: wy };
}

function applyZoomToRange(min, max, zoom, clampMin, clampMax, center) {
  if (!Number.isFinite(zoom) || zoom === 1.0) return { min, max };
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

function updateZoomValue(current, delta) {
  const factor = delta > 0 ? 1 / 1.1 : 1.1;
  const next = Math.min(6.0, Math.max(1.0, current * factor));
  return Math.abs(next - current) > 1.0e-4 ? next : current;
}

function setSourceToCustom() {
  const customRadio = sourceRadios.find((radio) => radio.value === 'custom');
  if (customRadio) {
    customRadio.checked = true;
  }
  nacaOptions.hidden = true;
  customOptions.hidden = false;
  if (databaseOptions) {
    databaseOptions.hidden = true;
  }
  controls4.hidden = true;
  controls5.hidden = true;
  controls6.hidden = true;
}

function drawAirfoilNodes(nb, bounds, activeIndex = null, hoverIndex = null) {
  if (!nb || !bounds) return;
  ctx.save();
  for (let i = 0; i < nb; i += 1) {
    const p = worldToCanvas(xb[i], yb[i], bounds);
    const isActive = i === activeIndex;
    const isHover = i === hoverIndex;
    ctx.beginPath();
    ctx.fillStyle = isActive ? '#2f7bff' : 'rgba(230, 236, 244, 0.85)';
    ctx.strokeStyle = isActive ? '#ffffff' : 'rgba(230, 236, 244, 0.65)';
    ctx.lineWidth = isActive ? 2.0 : 1.0;
    ctx.arc(p.x, p.y, isActive ? 4.8 : 3.4, 0, Math.PI * 2);
    ctx.fill();
    ctx.stroke();
    if (isHover && !isActive) {
      ctx.strokeStyle = '#2f7bff';
      ctx.lineWidth = 2.0;
      ctx.beginPath();
      ctx.arc(p.x, p.y, 6.2, 0, Math.PI * 2);
      ctx.stroke();
    }
  }
  ctx.restore();
}

function isTouchInput() {
  return ('ontouchstart' in window) || (navigator.maxTouchPoints ?? 0) > 0;
}

function getTouchPoint(event, rect) {
  if (!event.touches || event.touches.length === 0) return null;
  const touch = event.touches[0];
  return {
    x: touch.clientX - rect.left,
    y: touch.clientY - rect.top,
  };
}

function renderEditAirfoil() {
  if (!lastAirfoilBounds || !lastAirfoilNb) return;
  ctx.clearRect(0, 0, lastAirfoilBounds.width, lastAirfoilBounds.height);
  if (lastSolverPayload?.streamlines) {
    drawStreamlines(lastAirfoilBounds, lastAirfoilNb, lastSolverPayload.streamlines);
  }
  drawAirfoil(lastAirfoilNb, lastAirfoilBounds);
  if (lastSolverPayload?.hinge) {
    drawFlapHinge(lastAirfoilBounds, lastSolverPayload.hinge);
  }
  drawAirfoilNodes(lastAirfoilNb, lastAirfoilBounds, editDragIndex, editHoverIndex);
}

function findNearestNodeIndex(nb, bounds, mx, my, radius = 12) {
  let best = null;
  let bestDist = radius * radius;
  for (let i = 0; i < nb; i += 1) {
    const p = worldToCanvas(xb[i], yb[i], bounds);
    const dx = mx - p.x;
    const dy = my - p.y;
    const dist = dx * dx + dy * dy;
    if (dist <= bestDist) {
      best = i;
      bestDist = dist;
    }
  }
  return best;
}

function applyEditedAirfoil() {
  if (!lastAirfoilNb) return;
  const coords = [];
  for (let i = 0; i < lastAirfoilNb; i += 1) {
    coords.push({ x: xb[i], y: yb[i] });
  }
  const name = currentAirfoilName || 'Edited Airfoil';
  loadCustomAirfoil({ name, coords });
  setSourceToCustom();
  update();
}

async function applyCpEdit() {
  if (!cpEditPoints || !cpEditPoints.length) return;
  if (solverInFlight) return;
  const cpSpec = cpEditPoints.map((point) => point.cp);
  const result = await requestQdes(cpSpec, 10);
  if (!result?.ok) {
    console.warn(result?.error || 'QDES update failed.');
    return;
  }
  if (result.debug) {
    console.log('[QDES] result', result.debug);
  }
  if (!result.xb || !result.yb || !result.nb) {
    console.warn('QDES update missing airfoil geometry.');
    return;
  }
  const coords = [];
  for (let i = 0; i < result.nb; i += 1) {
    coords.push({ x: result.xb[i], y: result.yb[i] });
  }
  const name = currentAirfoilName || 'QDES Airfoil';
  loadCustomAirfoil({ name, coords });
  setSourceToCustom();
  update();
}

function updateHoverBox(point) {
  if (!hoverBox) return;
  if (!point) {
    if (hoverName) hoverName.textContent = '—';
    if (hoverMach) hoverMach.textContent = '—';
    if (hoverRe) hoverRe.textContent = '—';
    if (hoverNcr) hoverNcr.textContent = '—';
    if (hoverAlpha) hoverAlpha.textContent = '—';
    if (hoverFlap) hoverFlap.textContent = '—';
    return;
  }
  if (hoverName) hoverName.textContent = point.airfoil || '—';
  if (hoverMach) hoverMach.textContent = formatHoverValue(point.mach, 3);
  if (hoverRe) hoverRe.textContent = Number.isFinite(point.re) ? `${Math.round(point.re)}` : '—';
  if (hoverNcr) hoverNcr.textContent = formatHoverValue(point.ncr, 2);
  if (hoverAlpha) hoverAlpha.textContent = formatHoverValue(point.alpha, 2, '°');
  if (hoverFlap) hoverFlap.textContent = formatHoverValue(point.flapDef, 1, '°');
}

function setHoverState(source, point) {
  if (source === 'alpha') {
    alphaHover = point;
  } else {
    polarHover = point;
  }
  if (point) {
    activeHoverSource = source;
    updateHoverBox(point.data);
    return;
  }
  if (activeHoverSource === source) {
    const fallback = source === 'alpha' ? polarHover : alphaHover;
    if (fallback) {
      activeHoverSource = source === 'alpha' ? 'polar' : 'alpha';
      updateHoverBox(fallback.data);
    } else {
      activeHoverSource = null;
      updateHoverBox(null);
    }
  }
}

function findNearestPlotPoint(points, mx, my, radius = 10) {
  if (!points || points.length === 0) return null;
  const r2 = radius * radius;
  let best = null;
  let bestDist = r2;
  for (let i = 0; i < points.length; i += 1) {
    const p = points[i];
    const dx = mx - p.px;
    const dy = my - p.py;
    const d2 = dx * dx + dy * dy;
    if (d2 <= bestDist) {
      best = p;
      bestDist = d2;
    }
  }
  return best;
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
  const machValue = machInput ? parseFloat(machInput.value) : NaN;
  const reValue = Number.isFinite(coeffs?.REINF1)
    ? coeffs.REINF1
    : (reynoldsInput ? parseFloat(reynoldsInput.value) : NaN);
  const ncrValue = Number.isFinite(coeffs?.ACRIT?.[1])
    ? coeffs.ACRIT[1]
    : (ncrInput ? parseFloat(ncrInput.value) : NaN);
  const flapValue = flapDefInput ? parseFloat(flapDefInput.value) : NaN;
  const point = {
    alpha: alphaDeg,
    cl: coeffs.CL,
    cd: coeffs.CD,
    cm: coeffs.CM,
    converged,
    airfoil: currentAirfoilName || '—',
    mach: machValue,
    re: reValue,
    ncr: ncrValue,
    flapDef: flapValue,
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
  if (allPoints.length === 0) {
    alphaPlotState = null;
    return;
  }

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
  const baseX = { min: xmin, max: xmax };

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
  yminL = Math.max(yminL, -2.0);
  ymaxL = Math.min(ymaxL, 2.0);
  const baseYL = { min: yminL, max: ymaxL };

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
  yminR = Math.max(yminR, -2.0);
  ymaxR = Math.min(ymaxR, 2.0);
  const baseYR = { min: yminR, max: ymaxR };
  if (yminR <= 0.0 && ymaxR >= 0.0) {
    const zeroLeft = (0.0 - yminL) / (ymaxL - yminL);
    const zeroRight = (0.0 - yminR) / (ymaxR - yminR);
    const span = ymaxR - yminR || 1.0;
    const shift = (zeroRight - zeroLeft) * span;
    yminR += shift;
    ymaxR += shift;
  }

  const zoomX = applyZoomToRange(
    xmin,
    xmax,
    alphaZoom,
    baseX.min,
    baseX.max,
    alphaZoomCenter?.x,
  );
  xmin = zoomX.min;
  xmax = zoomX.max;
  const zoomYL = applyZoomToRange(
    yminL,
    ymaxL,
    alphaZoom,
    baseYL.min,
    baseYL.max,
    alphaZoomCenter?.yL,
  );
  yminL = zoomYL.min;
  ymaxL = zoomYL.max;
  const zoomYR = applyZoomToRange(
    yminR,
    ymaxR,
    alphaZoom,
    baseYR.min,
    baseYR.max,
    alphaZoomCenter?.yR,
  );
  yminR = zoomYR.min;
  ymaxR = zoomYR.max;

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

  const drawMarker = (shape, px, py, color, outlineOnly = false, size = 4) => {
    alphaCtx.fillStyle = color;
    alphaCtx.strokeStyle = color;
    if (shape === 'square') {
      if (outlineOnly) {
        alphaCtx.lineWidth = 1.8;
        alphaCtx.strokeRect(px - size, py - size, size * 2, size * 2);
      } else {
        alphaCtx.fillRect(px - size, py - size, size * 2, size * 2);
      }
      return;
    }
    if (shape === 'triangle') {
      alphaCtx.beginPath();
      alphaCtx.moveTo(px, py - size - 1);
      alphaCtx.lineTo(px + size + 1, py + size);
      alphaCtx.lineTo(px - size - 1, py + size);
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
    alphaCtx.arc(px, py, size, 0, Math.PI * 2);
    if (outlineOnly) {
      alphaCtx.lineWidth = 1.8;
      alphaCtx.stroke();
    } else {
      alphaCtx.fill();
    }
  };

  const staleColor = 'rgba(154, 160, 168, 0.85)';
  const hoverPoints = [];
  runCases.forEach((caseItem) => {
    caseItem.history.forEach((p) => {
      const px = xToPx(p.alpha);
      const color = p.converged === false ? staleColor : caseItem.color;
      const clPy = yToPyLeft(p.cl);
      const cdPy = yToPyRight(p.cd);
      const cmPy = yToPyLeft(p.cm);
      drawMarker('circle', px, clPy, color);
      drawMarker('square', px, cdPy, color, true);
      drawMarker('triangle', px, cmPy, color);
      const baseKey = `alpha-${caseItem.id}-${p.alpha.toFixed(3)}`;
      hoverPoints.push({
        key: `${baseKey}-cl`,
        px,
        py: clPy,
        type: 'cl',
        shape: 'circle',
        color,
        data: p,
      });
      hoverPoints.push({
        key: `${baseKey}-cd`,
        px,
        py: cdPy,
        type: 'cd',
        shape: 'square',
        color,
        data: p,
      });
      hoverPoints.push({
        key: `${baseKey}-cm`,
        px,
        py: cmPy,
        type: 'cm',
        shape: 'triangle',
        color,
        data: p,
      });
    });
  });

  const legendX = left + plotW - 70;
  const legendY = top + plotH - 50;
  const labelX = legendX + 12;
  alphaCtx.fillStyle = 'rgba(230, 236, 244, 0.85)';
  drawMarker('circle', legendX, legendY, alphaCtx.fillStyle);
  drawMarker('square', legendX, legendY + 18, alphaCtx.fillStyle, true);
  drawMarker('triangle', legendX, legendY + 36, alphaCtx.fillStyle);
  drawSubLabel(alphaCtx, 'C', 'L', labelX, legendY + 4);
  drawSubLabel(alphaCtx, 'C', 'D', labelX, legendY + 22);
  drawSubLabel(alphaCtx, 'C', 'M', labelX, legendY + 40);

  const hoverKey = alphaHover?.key;
  const hoverPoint = hoverKey ? hoverPoints.find((p) => p.key === hoverKey) : null;
  if (hoverPoint) {
    alphaCtx.strokeStyle = '#ffffff';
    alphaCtx.lineWidth = 2.4;
    drawMarker(hoverPoint.shape, hoverPoint.px, hoverPoint.py, '#ffffff', true, 6);
    let valueLine = null;
    if (hoverPoint.type === 'cl') {
      valueLine = { base: 'C', sub: 'L', value: formatNum(hoverPoint.data.cl, 3) };
    } else if (hoverPoint.type === 'cd') {
      valueLine = { base: 'C', sub: 'D', value: formatNum(hoverPoint.data.cd, 4) };
    } else {
      valueLine = { base: 'C', sub: 'M', value: formatNum(hoverPoint.data.cm, 3) };
    }
    const lines = [
      `α ${formatHoverValue(hoverPoint.data.alpha, 2, '°')}`,
      valueLine,
    ];
    drawHoverLabel(alphaCtx, hoverPoint.px, hoverPoint.py, lines, { w, h });
  }
  alphaCtx.restore();
  alphaPlotState = {
    points: hoverPoints,
    w,
    h,
    mapping: {
      left,
      right,
      top,
      bottom,
      plotW,
      plotH,
      xmin,
      xmax,
      yminL,
      ymaxL,
      yminR,
      ymaxR,
    },
  };
}

// Polar plot: CL vs CD with nonnegative CD axis, point-only styling.
function drawPolarPlot() {
  if (!polarCanvas || !polarCtx) return;
  const polarRect = polarCanvas.getBoundingClientRect();
  polarCtx.clearRect(0, 0, polarRect.width, polarRect.height);
  const allPoints = runCases.flatMap((c) => c.history);
  if (allPoints.length === 0) {
    polarPlotState = null;
    return;
  }

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
  xmax = Math.min(xmax, 0.04);
  ymin = Math.max(ymin, -2.0);
  ymax = Math.min(ymax, 2.0);
  const baseX = { min: xmin, max: xmax };
  const baseY = { min: ymin, max: ymax };

  const zoomX = applyZoomToRange(
    xmin,
    xmax,
    polarZoom,
    baseX.min,
    baseX.max,
    polarZoomCenter?.x,
  );
  xmin = zoomX.min;
  xmax = zoomX.max;
  const zoomY = applyZoomToRange(
    ymin,
    ymax,
    polarZoom,
    baseY.min,
    baseY.max,
    polarZoomCenter?.y,
  );
  ymin = zoomY.min;
  ymax = zoomY.max;

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
  const hoverPoints = [];
  runCases.forEach((caseItem) => {
    caseItem.history.forEach((p) => {
      if (!Number.isFinite(p.cd) || !Number.isFinite(p.cl)) return;
      const px = xToPx(Math.max(p.cd, 0.0));
      const py = yToPy(p.cl);
      polarCtx.fillStyle = p.converged === false ? staleColor : caseItem.color;
      polarCtx.beginPath();
      polarCtx.arc(px, py, 4, 0, Math.PI * 2);
      polarCtx.fill();
      hoverPoints.push({
        key: `polar-${caseItem.id}-${p.alpha.toFixed(3)}-${p.cd.toFixed(4)}`,
        px,
        py,
        type: 'polar',
        shape: 'circle',
        color: polarCtx.fillStyle,
        data: p,
      });
    });
  });

  const hoverKey = polarHover?.key;
  const hoverPoint = hoverKey ? hoverPoints.find((p) => p.key === hoverKey) : null;
  if (hoverPoint) {
    polarCtx.strokeStyle = '#ffffff';
    polarCtx.lineWidth = 2.4;
    polarCtx.beginPath();
    polarCtx.arc(hoverPoint.px, hoverPoint.py, 6, 0, Math.PI * 2);
    polarCtx.stroke();
    const lines = [
      `α ${formatHoverValue(hoverPoint.data.alpha, 2, '°')}`,
      { base: 'C', sub: 'L', value: formatNum(hoverPoint.data.cl, 3) },
      { base: 'C', sub: 'D', value: formatNum(hoverPoint.data.cd, 4) },
    ];
    drawHoverLabel(polarCtx, hoverPoint.px, hoverPoint.py, lines, { w, h });
  }
  polarCtx.restore();
  polarPlotState = {
    points: hoverPoints,
    w,
    h,
    mapping: {
      left,
      right,
      top,
      bottom,
      plotW,
      plotH,
      xmin,
      xmax,
      ymin,
      ymax,
    },
  };
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
  let sx = bounds.offsetX + (wx - bounds.xmin) * bounds.scale;
  let sy = bounds.height - bounds.offsetY - (wy - bounds.ymin) * bounds.scale;
  if (bounds.zoom && bounds.zoomCenterScreen) {
    sx = bounds.zoomCenterScreen.x + (sx - bounds.zoomCenterScreen.x) * bounds.zoom;
    sy = bounds.zoomCenterScreen.y + (sy - bounds.zoomCenterScreen.y) * bounds.zoom;
  }
  return { x: sx, y: sy };
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
function buildCpEditPoints(nb, cpData) {
  if (!cpData) return null;
  const invAll = Array.isArray(cpData.invAll) ? cpData.invAll : [];
  const source = invAll.length ? invAll : cpData.upper.concat(cpData.lower);
  if (!source.length) return null;
  const count = Math.min(nb, source.length);
  const lePt = cpData.le;
  const tePt = cpData.te;
  if (!lePt || !tePt) return null;
  const dxChord = tePt.x - lePt.x;
  const dyChord = tePt.y - lePt.y;
  const chord2 = dxChord * dxChord + dyChord * dyChord || 1.0;
  let leIdx = 0;
  let leDist = Infinity;
  for (let i = 0; i < count; i += 1) {
    const p = source[i];
    const dx = p.x - lePt.x;
    const dy = p.y - lePt.y;
    const dist = dx * dx + dy * dy;
    if (dist < leDist) {
      leDist = dist;
      leIdx = i;
    }
  }
  const points = new Array(count);
  for (let i = 0; i < count; i += 1) {
    const p = source[i];
    const s = ((p.x - lePt.x) * dxChord + (p.y - lePt.y) * dyChord) / chord2;
    const side = i <= leIdx ? 'upper' : 'lower';
    points[i] = { index: i, s, cp: p.cp, side };
  }
  return points;
}

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
  const baseCp = { min: cpMin, max: cpMax };

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
  let xMin = 0.0;
  let xMax = 1.0;
  const baseX = { min: xMin, max: xMax };
  const zoomX = applyZoomToRange(
    xMin,
    xMax,
    cpZoom,
    baseX.min,
    baseX.max,
    cpZoomCenter?.x,
  );
  xMin = zoomX.min;
  xMax = zoomX.max;
  cpMin = baseCp.min;
  cpMax = baseCp.max;
  const xToPx = (s) => leX + ((s - xMin) / (xMax - xMin)) * span;
  const cpToPy = (cp) => top + ((cp - cpMin) / (cpMax - cpMin)) * plotH;
  const chordFrac = (x, y) => ((x - lePt.x) * dxChord + (y - lePt.y) * dyChord) / chord2;

  lastCpPlot = {
    nb,
    cpUpper,
    cpLower,
    lePt,
    tePt,
    bounds,
    cpInvAll,
    cpWake,
    mapping: {
      left,
      right,
      top,
      bottom,
      plotW,
      plotH,
      leX,
      teX,
      span,
      xMin,
      xMax,
      cpMin,
      cpMax,
    },
  };

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
    const xv = xMin + t * (xMax - xMin);
    const x = xToPx(xv);
    cpCtx.beginPath();
    cpCtx.moveTo(x, top + plotH);
    cpCtx.lineTo(x, top + plotH + 8);
    cpCtx.stroke();
    cpCtx.fillText(xv.toFixed(2), x - 10, top + plotH + 22);
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

  if (Number.isFinite(cpHoverS)) {
    const s = Math.max(xMin, Math.min(xMax, cpHoverS));
    const hoverX = xToPx(s);
    const findClosest = (arr) => {
      let best = null;
      let bestDist = Infinity;
      arr.forEach((p) => {
        const frac = chordFrac(p.x, p.y);
        const dist = Math.abs(frac - s);
        if (dist < bestDist) {
          bestDist = dist;
          best = { p, frac };
        }
      });
      return best;
    };

    const upperHit = findClosest(cpUpper);
    const lowerHit = findClosest(cpLower);
    cpCtx.save();
    cpCtx.strokeStyle = 'rgba(255, 255, 255, 0.35)';
    cpCtx.lineWidth = 1.0;
    cpCtx.beginPath();
    cpCtx.moveTo(hoverX, top);
    cpCtx.lineTo(hoverX, top + plotH);
    cpCtx.stroke();

    const drawHit = (hit, color, dyText) => {
      if (!hit) return;
      const py = cpToPy(hit.p.cp);
      cpCtx.fillStyle = color;
      cpCtx.beginPath();
      cpCtx.arc(hoverX, py, 4.2, 0, Math.PI * 2);
      cpCtx.fill();
      const label = hit.p.cp.toFixed(3);
      const textX = Math.min(w - 40, Math.max(left + 6, hoverX + 8));
      const textY = Math.min(top + plotH - 6, Math.max(top + 12, py + dyText));
      cpCtx.fillStyle = 'rgba(230, 236, 244, 0.9)';
      cpCtx.fillText(label, textX, textY);
    };

    drawHit(upperHit, '#2f7bff', -8);
    drawHit(lowerHit, '#ff4a3d', 14);
    cpCtx.restore();
  }

  if (cpEditMode && Array.isArray(cpEditPoints)) {
    cpCtx.save();
    cpCtx.lineWidth = 1.2;
    cpEditPoints.forEach((point, idx) => {
      const px = xToPx(point.s);
      const py = cpToPy(point.cp);
      const isActive = idx === cpEditDragIndex;
      const isHover = idx === cpEditHoverIndex;
      const color = point.side === 'lower' ? '#ff4a3d' : '#2f7bff';
      cpCtx.strokeStyle = 'rgba(18, 22, 27, 0.9)';
      cpCtx.fillStyle = color;
      const r = isActive ? 5.5 : 4.0;
      cpCtx.beginPath();
      cpCtx.arc(px, py, r, 0, Math.PI * 2);
      cpCtx.fill();
      cpCtx.stroke();
      if (isHover && !isActive) {
        cpCtx.strokeStyle = color;
        cpCtx.lineWidth = 2.0;
        cpCtx.beginPath();
        cpCtx.arc(px, py, r + 2.4, 0, Math.PI * 2);
        cpCtx.stroke();
        cpCtx.strokeStyle = 'rgba(18, 22, 27, 0.9)';
        cpCtx.lineWidth = 1.2;
      }
    });
    cpCtx.restore();
  }

  cpCtx.restore();
}

function rebuildCpEditPoints(nb, cpData, force = false) {
  if (!cpEditMode) {
    cpEditPoints = null;
    cpEditDragIndex = null;
    cpEditHoverIndex = null;
    return;
  }
  if (!force && cpEditPoints) {
    return;
  }
  cpEditPoints = buildCpEditPoints(nb, cpData);
  cpEditDragIndex = null;
  cpEditHoverIndex = null;
}

function getCpEditPointScreen(point, mapping) {
  if (!mapping) return null;
  const { leX, span, xMin, xMax, top, plotH, cpMin, cpMax } = mapping;
  const px = leX + ((point.s - xMin) / (xMax - xMin)) * span;
  const py = top + ((point.cp - cpMin) / (cpMax - cpMin)) * plotH;
  return { x: px, y: py };
}

function findClosestCpEditPoint(mx, my, mapping, radius = 10) {
  if (!cpEditPoints || !mapping) return null;
  let best = null;
  let bestDist = radius * radius;
  cpEditPoints.forEach((point, idx) => {
    const screen = getCpEditPointScreen(point, mapping);
    if (!screen) return;
    const dx = mx - screen.x;
    const dy = my - screen.y;
    const dist = dx * dx + dy * dy;
    if (dist <= bestDist) {
      best = idx;
      bestDist = dist;
    }
  });
  return best;
}

function screenToCpValue(my, mapping) {
  if (!mapping) return null;
  const { top, plotH, cpMin, cpMax } = mapping;
  const clamped = Math.max(top, Math.min(top + plotH, my));
  return cpMin + ((clamped - top) / plotH) * (cpMax - cpMin);
}

function renderCpPlotFromCache() {
  if (!lastCpPlot) return;
  drawCpPlot(
    lastCpPlot.nb,
    lastCpPlot.cpUpper,
    lastCpPlot.cpLower,
    lastCpPlot.lePt,
    lastCpPlot.tePt,
    lastCpPlot.bounds,
    lastCpPlot.cpInvAll,
    lastCpPlot.cpWake,
  );
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
  const prevPayload = lastSolverPayload;
  const samePayload = prevPayload === payload;
  lastSolverPayload = payload;
  if (pendingGeometrySettings) {
    lastGeometrySettings = pendingGeometrySettings;
  }
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
  const zoomedBounds = {
    ...framedBounds,
    zoom: airfoilZoom,
    zoomCenterScreen: airfoilZoomCenterScreen,
  };
  lastAirfoilBounds = zoomedBounds;
  lastAirfoilNb = nb;
  if (xbResult && ybResult) {
    xb.set(xbResult);
    yb.set(ybResult);
  }
  if (airfoilName) {
    currentAirfoilName = airfoilName;
  }

  ctx.clearRect(0, 0, framedBounds.width, framedBounds.height);

  if (streamlines) {
    drawStreamlines(zoomedBounds, nb, streamlines);
  }

  if (cpData?.upper?.length || cpData?.lower?.length) {
    rebuildCpEditPoints(nb, cpData, !samePayload);
    drawCpPlot(
      nb,
      cpData.upper,
      cpData.lower,
      cpData.le,
      cpData.te,
      zoomedBounds,
      cpData.invAll || [],
      cpData.wake || [],
    );
  } else if (cpCtx) {
    const cpRect = cpCanvas.getBoundingClientRect();
    cpCtx.clearRect(0, 0, cpRect.width, cpRect.height);
    lastCpPlot = null;
    rebuildCpEditPoints(0, null);
  }

  updateDataBox(alphaRad, coeffs, converged);
  upsertSweepPoint(alphaDeg, coeffs, converged);
  drawAlphaSweepPlot();
  drawPolarPlot();

  drawAirfoil(nb, zoomedBounds);
  drawFlapHinge(zoomedBounds, hinge);
  if (viscous && blLines) {
    drawBoundaryLayerLines(zoomedBounds, blLines);
  }
  if (editMode) {
    drawAirfoilNodes(nb, zoomedBounds, editDragIndex, editHoverIndex);
  }

  updateDownloadButtons();
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

function resolvePendingQdes(id, payload) {
  const resolve = pendingQdes.get(id);
  if (resolve) {
    resolve(payload);
    pendingQdes.delete(id);
  }
}

function cancelPendingQdes() {
  pendingQdes.forEach((resolve) => resolve({ canceled: true }));
  pendingQdes.clear();
}

function resolvePendingDump(id, payload) {
  const resolve = pendingDumps.get(id);
  if (resolve) {
    resolve(payload);
    pendingDumps.delete(id);
  }
}

function cancelPendingDumps() {
  pendingDumps.forEach((resolve) => resolve({ canceled: true }));
  pendingDumps.clear();
}

function spawnSolverWorker() {
  solverWorker = new Worker(new URL('./solver_worker.js', import.meta.url), { type: 'module' });
  solverWorker.onmessage = (event) => {
    const payload = event.data;
    if (payload?.type === 'dump') {
      resolvePendingDump(payload.requestId, payload);
      return;
    }
    if (payload?.type === 'qdes') {
      resolvePendingQdes(payload.requestId, payload);
      return;
    }
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
    cancelPendingQdes();
    cancelPendingDumps();
    console.warn('Solver worker error:', event);
  };
  solverWorker.onmessageerror = (event) => {
    solverInFlight = false;
    cancelPendingSolves();
    cancelPendingQdes();
    cancelPendingDumps();
    console.warn('Solver worker message error:', event);
  };
}

function downloadTextFile(filename, content) {
  const blob = new Blob([content], { type: 'text/plain;charset=ascii' });
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement('a');
  anchor.href = url;
  anchor.download = filename;
  document.body.appendChild(anchor);
  anchor.click();
  document.body.removeChild(anchor);
  URL.revokeObjectURL(url);
}

function requestDump(kind) {
  if (!solverWorker) {
    return Promise.resolve({ ok: false, error: 'Solver worker not ready.' });
  }
  const requestId = nextDumpId;
  nextDumpId += 1;
  solverWorker.postMessage({ requestId, action: 'dump', kind, kdelim: 1 });
  return new Promise((resolve) => {
    pendingDumps.set(requestId, resolve);
  });
}

function requestQdes(cpSpec, niter = 10) {
  if (!solverWorker) {
    return Promise.resolve({ ok: false, error: 'Solver worker not ready.' });
  }
  if (!Array.isArray(cpSpec) || cpSpec.length === 0) {
    return Promise.resolve({ ok: false, error: 'Cp edit points are missing.' });
  }
  const requestId = nextQdesId;
  nextQdesId += 1;
  solverWorker.postMessage({ requestId, action: 'qdes', cpSpec, niter });
  return new Promise((resolve) => {
    pendingQdes.set(requestId, resolve);
  });
}

function updateDownloadButtons() {
  const enabled = !!lastSolverPayload?.ok;
  if (downloadCpButton) {
    downloadCpButton.disabled = !enabled;
  }
  if (downloadBlButton) {
    downloadBlButton.disabled = !enabled;
  }
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

  if (mode === '4') {
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

  const profile = series6Profile.value;
  const cl = parseFloat(cl6Input.value);
  const camber = inferSixSeriesCamber(profile, cl);
  const fallbackA = defaultSixSeriesA(profile);

  const uiGeometry = {
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
    custom: customAirfoil
      ? {
        name: customAirfoil.name,
        nb: customAirfoil.nb,
        x: xbCustom,
        y: ybCustom,
      }
      : null,
  };

  let geometry = uiGeometry;
  if ((source === 'custom' || source === 'database') && !customAirfoil) {
    if (lastGeometrySettings) {
      geometry = lastGeometrySettings;
    } else {
      currentAirfoilName = source === 'database' ? 'Select database airfoil' : 'Load .dat airfoil';
      return Promise.resolve({ skipped: true });
    }
  }

  if ((geometry.source === 'custom' || geometry.source === 'database') && !geometry.custom) {
    return Promise.resolve({ skipped: true });
  }

  let geometryKey = '';
  if (geometry.source === 'custom' || geometry.source === 'database') {
    geometryKey = `${geometry.source}:${geometry.custom?.name || 'none'}:${customAirfoilVersion}`;
  } else if (geometry.mode === '4') {
    geometryKey = `naca4:${geometry.m}:${geometry.p}:${geometry.t}`;
  } else if (geometry.mode === '5') {
    geometryKey = `naca5:${geometry.series5}:${geometry.t5}`;
  } else {
    geometryKey = `naca6:${geometry.profile6}:${geometry.t6}:${geometry.cl6}`;
  }
  const flapKey = `${flapXInput?.value ?? '0.75'}:${flapYInput?.value ?? '0.5'}:${flapDefInput?.value ?? '0'}`;
  geometryKey = `${geometryKey}:flap:${flapKey}`;

  const rect = canvas.getBoundingClientRect();
  const advancedMode = !!advancedModeToggle?.checked;
  const reuseEnabled = viscousToggle.checked && !!reuseSolutionToggle?.checked;
  const reuseKeyedState = reuseState?.geometryKey === geometryKey
    ? reuseState
    : null;
  const settings = {
    mode: geometry.mode,
    source: geometry.source,
    m: geometry.m,
    p: geometry.p,
    t: geometry.t,
    series5: geometry.series5,
    t5: geometry.t5,
    profile6: geometry.profile6,
    t6: geometry.t6,
    cl6: geometry.cl6,
    camber6: geometry.camber6,
    fallbackA6: geometry.fallbackA6,
    flap: {
      x: parseFloat(flapXInput?.value),
      yrel: parseFloat(flapYInput?.value),
      deflection: parseFloat(flapDefInput?.value),
    },
    custom: geometry.custom,
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
    cancelPendingDumps();
    spawnSolverWorker();
    solverInFlight = false;
  }

  solverRequestId += 1;
  latestSolverId = solverRequestId;
  solverInFlight = true;
  pendingGeometrySettings = geometry;
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

if (cpCanvas) {
  cpCanvas.addEventListener('wheel', (event) => {
    event.preventDefault();
    let cpAnchor = null;
    if (lastCpPlot?.mapping) {
      const rect = cpCanvas.getBoundingClientRect();
      const mx = event.clientX - rect.left;
      const my = event.clientY - rect.top;
      const {
        leX,
        top,
        span,
        plotH,
        xMin,
        xMax,
      } = lastCpPlot.mapping;
      if (mx >= Math.min(leX, leX + span)
        && mx <= Math.max(leX, leX + span)
        && my >= top
        && my <= top + plotH) {
        const x = xMin + ((mx - leX) / span) * (xMax - xMin);
        cpZoomCenter = { x };
        cpAnchor = { x, mx };
        if (lastCpPlot.lePt && lastCpPlot.tePt && lastAirfoilBounds) {
          const dx = lastCpPlot.tePt.x - lastCpPlot.lePt.x;
          const dy = lastCpPlot.tePt.y - lastCpPlot.lePt.y;
          const chord2 = dx * dx + dy * dy || 1.0;
          const world = {
            x: lastCpPlot.lePt.x + dx * x,
            y: lastCpPlot.lePt.y + dy * x,
          };
          const screen = worldToCanvas(world.x, world.y, lastAirfoilBounds);
          airfoilZoomCenterScreen = { x: screen.x, y: screen.y };
        }
      }
    }
    cpZoom = updateZoomValue(cpZoom, event.deltaY);
    airfoilZoom = cpZoom;
    if (lastSolverPayload) {
      applySolverResult(lastSolverPayload);
      if (cpAnchor && lastCpPlot?.mapping) {
        const { leX, span } = lastCpPlot.mapping;
        const f = span !== 0 ? (cpAnchor.mx - leX) / span : 0.5;
        const center = cpAnchor.x + (0.5 - f) / cpZoom;
        cpZoomCenter = { x: Math.max(0.0, Math.min(1.0, center)) };
        renderCpPlotFromCache();
      }
    } else {
      renderCpPlotFromCache();
    }
  }, { passive: false });

  cpCanvas.addEventListener('mousedown', (event) => {
    if (event.ctrlKey && lastCpPlot?.mapping) {
      const rect = cpCanvas.getBoundingClientRect();
      cpPanActive = true;
      cpPanStart = {
        x: event.clientX - rect.left,
        center: cpZoomCenter?.x ?? null,
      };
      event.preventDefault();
      return;
    }
    if (!cpEditMode || !lastCpPlot?.mapping) return;
    const rect = cpCanvas.getBoundingClientRect();
    const x = event.clientX - rect.left;
    const y = event.clientY - rect.top;
    const idx = findClosestCpEditPoint(x, y, lastCpPlot.mapping, 12);
    if (idx == null) return;
    event.preventDefault();
    cpEditDragIndex = idx;
    cpEditHoverIndex = idx;
    renderCpPlotFromCache();
  });

  cpCanvas.addEventListener('mousemove', (event) => {
    if (!lastCpPlot?.mapping) return;
    const rect = cpCanvas.getBoundingClientRect();
    const x = event.clientX - rect.left;
    const y = event.clientY - rect.top;
    const { top, plotH, leX, teX, xMin, xMax } = lastCpPlot.mapping;
    const minX = Math.min(leX, teX);
    const maxX = Math.max(leX, teX);
    if (cpPanActive && cpPanStart) {
      const span = teX - leX || 1.0;
      const dx = x - cpPanStart.x;
      const baseCenter = cpPanStart.center ?? 0.5 * (xMin + xMax);
      const delta = (dx / span) * (xMax - xMin);
      cpZoomCenter = { x: baseCenter - delta };
      renderCpPlotFromCache();
      return;
    }
    if (cpEditMode && cpEditPoints) {
      if (cpEditDragIndex != null) {
        const cpVal = screenToCpValue(y, lastCpPlot.mapping);
        if (Number.isFinite(cpVal)) {
          cpEditPoints[cpEditDragIndex].cp = cpVal;
        }
        renderCpPlotFromCache();
        return;
      }
      const hoverIdx = findClosestCpEditPoint(x, y, lastCpPlot.mapping, 12);
      if (hoverIdx !== cpEditHoverIndex) {
        cpEditHoverIndex = hoverIdx;
        renderCpPlotFromCache();
      }
    }
    if (x < minX || x > maxX || y < top || y > top + plotH) {
      if (cpHoverS !== null) {
        cpHoverS = null;
        renderCpPlotFromCache();
      }
      return;
    }
    const span = teX - leX || 1.0;
    cpHoverS = xMin + ((x - leX) / span) * (xMax - xMin);
    renderCpPlotFromCache();
  });

  cpCanvas.addEventListener('mouseleave', () => {
    if (cpPanActive) {
      cpPanActive = false;
      cpPanStart = null;
    }
    if (cpEditMode && cpEditDragIndex != null) {
      cpEditDragIndex = null;
      cpEditHoverIndex = null;
      applyCpEdit();
      return;
    }
    if (cpEditHoverIndex != null) {
      cpEditHoverIndex = null;
      renderCpPlotFromCache();
    }
    if (cpHoverS !== null) {
      cpHoverS = null;
      renderCpPlotFromCache();
    }
  });

  cpCanvas.addEventListener('mouseup', () => {
    if (cpPanActive) {
      cpPanActive = false;
      cpPanStart = null;
      return;
    }
    if (!cpEditMode || cpEditDragIndex == null) return;
    cpEditDragIndex = null;
    cpEditHoverIndex = null;
    applyCpEdit();
  });

  if (isTouchInput()) {
    cpCanvas.addEventListener('touchstart', (event) => {
      if (!cpEditMode || !lastCpPlot?.mapping) return;
      const rect = cpCanvas.getBoundingClientRect();
      const point = getTouchPoint(event, rect);
      if (!point) return;
      const idx = findClosestCpEditPoint(point.x, point.y, lastCpPlot.mapping, 18);
      if (idx == null) return;
      event.preventDefault();
      cpEditHoverIndex = idx;
      renderCpPlotFromCache();
      if (cpTouchTimer) clearTimeout(cpTouchTimer);
      cpTouchTimer = setTimeout(() => {
        cpEditDragIndex = idx;
        touchCpDrag = true;
        renderCpPlotFromCache();
      }, 250);
    }, { passive: false });

    cpCanvas.addEventListener('touchmove', (event) => {
      if (!cpEditMode || !lastCpPlot?.mapping) return;
      const rect = cpCanvas.getBoundingClientRect();
      const point = getTouchPoint(event, rect);
      if (!point) return;
      if (touchCpDrag && cpEditDragIndex != null) {
        event.preventDefault();
        const cpVal = screenToCpValue(point.y, lastCpPlot.mapping);
        if (Number.isFinite(cpVal)) {
          cpEditPoints[cpEditDragIndex].cp = cpVal;
        }
        renderCpPlotFromCache();
      }
    }, { passive: false });

    cpCanvas.addEventListener('touchend', () => {
      if (cpTouchTimer) {
        clearTimeout(cpTouchTimer);
        cpTouchTimer = null;
      }
      if (!cpEditMode) return;
      if (touchCpDrag && cpEditDragIndex != null) {
        cpEditDragIndex = null;
        cpEditHoverIndex = null;
        touchCpDrag = false;
        applyCpEdit();
      }
    });
  }
}

if (alphaCanvas) {
  alphaCanvas.addEventListener('wheel', (event) => {
    event.preventDefault();
    if (alphaPlotState?.mapping) {
      const rect = alphaCanvas.getBoundingClientRect();
      const mx = event.clientX - rect.left;
      const my = event.clientY - rect.top;
      const {
        left,
        top,
        plotW,
        plotH,
        xmin,
        xmax,
        yminL,
        ymaxL,
        yminR,
        ymaxR,
      } = alphaPlotState.mapping;
      if (mx >= left && mx <= left + plotW && my >= top && my <= top + plotH) {
        const x = xmin + ((mx - left) / plotW) * (xmax - xmin);
        const yL = ymaxL - ((my - top) / plotH) * (ymaxL - yminL);
        const yR = ymaxR - ((my - top) / plotH) * (ymaxR - yminR);
        alphaZoomCenter = { x, yL, yR };
      }
    }
    alphaZoom = updateZoomValue(alphaZoom, event.deltaY);
    drawAlphaSweepPlot();
  }, { passive: false });

  alphaCanvas.addEventListener('mousemove', (event) => {
    if (!alphaPlotState?.points) return;
    const rect = alphaCanvas.getBoundingClientRect();
    const mx = event.clientX - rect.left;
    const my = event.clientY - rect.top;
    const hit = findNearestPlotPoint(alphaPlotState.points, mx, my, 10);
    const nextKey = hit?.key || null;
    if (nextKey !== alphaHover?.key) {
      setHoverState('alpha', hit);
      drawAlphaSweepPlot();
    }
  });

  alphaCanvas.addEventListener('mouseleave', () => {
    if (alphaHover) {
      setHoverState('alpha', null);
      drawAlphaSweepPlot();
    }
  });
}

if (polarCanvas) {
  polarCanvas.addEventListener('wheel', (event) => {
    event.preventDefault();
    if (polarPlotState?.mapping) {
      const rect = polarCanvas.getBoundingClientRect();
      const mx = event.clientX - rect.left;
      const my = event.clientY - rect.top;
      const {
        left,
        top,
        plotW,
        plotH,
        xmin,
        xmax,
        ymin,
        ymax,
      } = polarPlotState.mapping;
      if (mx >= left && mx <= left + plotW && my >= top && my <= top + plotH) {
        const x = xmin + ((mx - left) / plotW) * (xmax - xmin);
        const y = ymax - ((my - top) / plotH) * (ymax - ymin);
        polarZoomCenter = { x, y };
      }
    }
    polarZoom = updateZoomValue(polarZoom, event.deltaY);
    drawPolarPlot();
  }, { passive: false });

  polarCanvas.addEventListener('mousemove', (event) => {
    if (!polarPlotState?.points) return;
    const rect = polarCanvas.getBoundingClientRect();
    const mx = event.clientX - rect.left;
    const my = event.clientY - rect.top;
    const hit = findNearestPlotPoint(polarPlotState.points, mx, my, 10);
    const nextKey = hit?.key || null;
    if (nextKey !== polarHover?.key) {
      setHoverState('polar', hit);
      drawPolarPlot();
    }
  });

  polarCanvas.addEventListener('mouseleave', () => {
    if (polarHover) {
      setHoverState('polar', null);
      drawPolarPlot();
    }
  });
}

if (downloadCpButton) {
  downloadCpButton.addEventListener('click', async () => {
    const payload = await requestDump('cp');
    if (!payload?.ok || !payload.content) {
      console.warn(payload?.error || 'CPDUMP unavailable.');
      return;
    }
    downloadTextFile(payload.filename || 'xfoil.cp', payload.content);
  });
}

if (downloadBlButton) {
  downloadBlButton.addEventListener('click', async () => {
    const payload = await requestDump('bl');
    if (!payload?.ok || !payload.content) {
      console.warn(payload?.error || 'BLDUMP unavailable.');
      return;
    }
    downloadTextFile(payload.filename || 'xfoil.bl', payload.content);
  });
}

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

if (editAirfoilButton) {
  editAirfoilButton.addEventListener('click', () => {
    editMode = !editMode;
    editAirfoilButton.classList.toggle('active', editMode);
    editAirfoilButton.setAttribute('aria-pressed', editMode ? 'true' : 'false');
    if (airfoilFrame) {
      airfoilFrame.classList.toggle('editing', editMode);
    }
    editDragIndex = null;
    if (editMode) {
      renderEditAirfoil();
    } else if (lastSolverPayload) {
      applySolverResult(lastSolverPayload);
    }
  });
}

if (editCpButton) {
  editCpButton.addEventListener('click', () => {
    cpEditMode = !cpEditMode;
    editCpButton.classList.toggle('active', cpEditMode);
    editCpButton.setAttribute('aria-pressed', cpEditMode ? 'true' : 'false');
    if (cpCanvas) {
      cpCanvas.parentElement?.classList.toggle('editing', cpEditMode);
    }
    rebuildCpEditPoints(lastCpPlot?.nb ?? 0, lastSolverPayload?.cpData ?? null);
    renderCpPlotFromCache();
  });
}

if (canvas) {
  canvas.addEventListener('wheel', (event) => {
    event.preventDefault();
    if (lastAirfoilBounds) {
      const rect = canvas.getBoundingClientRect();
      const mx = event.clientX - rect.left;
      const my = event.clientY - rect.top;
      airfoilZoomCenterScreen = { x: mx, y: my };
      if (lastCpPlot?.lePt && lastCpPlot?.tePt) {
        const world = canvasToWorld(mx, my, lastAirfoilBounds);
        const dx = lastCpPlot.tePt.x - lastCpPlot.lePt.x;
        const dy = lastCpPlot.tePt.y - lastCpPlot.lePt.y;
        const chord2 = dx * dx + dy * dy || 1.0;
        const s = ((world.x - lastCpPlot.lePt.x) * dx + (world.y - lastCpPlot.lePt.y) * dy) / chord2;
        cpZoomCenter = { x: s };
      }
    }
    airfoilZoom = updateZoomValue(airfoilZoom, event.deltaY);
    cpZoom = airfoilZoom;
    if (lastSolverPayload) {
      applySolverResult(lastSolverPayload);
    }
  }, { passive: false });

  canvas.addEventListener('mousedown', (event) => {
    if (!editMode || !lastAirfoilBounds) return;
    const rect = canvas.getBoundingClientRect();
    const mx = event.clientX - rect.left;
    const my = event.clientY - rect.top;
    const idx = findNearestNodeIndex(lastAirfoilNb, lastAirfoilBounds, mx, my, 12);
    if (idx == null) return;
    event.preventDefault();
    editDragIndex = idx;
    editHoverIndex = idx;
    renderEditAirfoil();
  });

  canvas.addEventListener('mousemove', (event) => {
    if (!editMode || !lastAirfoilBounds) return;
    const rect = canvas.getBoundingClientRect();
    const mx = event.clientX - rect.left;
    const my = event.clientY - rect.top;
    if (editDragIndex != null) {
      const { x, y } = canvasToWorld(mx, my, lastAirfoilBounds);
      xb[editDragIndex] = x;
      yb[editDragIndex] = y;
      renderEditAirfoil();
      return;
    }
    const hoverIdx = findNearestNodeIndex(lastAirfoilNb, lastAirfoilBounds, mx, my, 12);
    if (hoverIdx !== editHoverIndex) {
      editHoverIndex = hoverIdx;
      renderEditAirfoil();
    }
  });

  const finishDrag = () => {
    if (!editMode || editDragIndex == null) return;
    editDragIndex = null;
    editHoverIndex = null;
    applyEditedAirfoil();
  };

  canvas.addEventListener('mouseup', finishDrag);
  canvas.addEventListener('mouseleave', () => {
    if (!editMode) return;
    if (editDragIndex != null) {
      finishDrag();
      return;
    }
    if (editHoverIndex != null) {
      editHoverIndex = null;
      renderEditAirfoil();
    }
  });

  if (isTouchInput()) {
    canvas.addEventListener('touchstart', (event) => {
      if (!editMode || !lastAirfoilBounds) return;
      const rect = canvas.getBoundingClientRect();
      const point = getTouchPoint(event, rect);
      if (!point) return;
      const idx = findNearestNodeIndex(lastAirfoilNb, lastAirfoilBounds, point.x, point.y, 18);
      if (idx == null) return;
      event.preventDefault();
      editHoverIndex = idx;
      renderEditAirfoil();
      if (airfoilTouchTimer) clearTimeout(airfoilTouchTimer);
      airfoilTouchTimer = setTimeout(() => {
        editDragIndex = idx;
        touchAirfoilDrag = true;
        renderEditAirfoil();
      }, 250);
    }, { passive: false });

    canvas.addEventListener('touchmove', (event) => {
      if (!editMode || !lastAirfoilBounds) return;
      const rect = canvas.getBoundingClientRect();
      const point = getTouchPoint(event, rect);
      if (!point) return;
      if (touchAirfoilDrag && editDragIndex != null) {
        event.preventDefault();
        const { x, y } = canvasToWorld(point.x, point.y, lastAirfoilBounds);
        xb[editDragIndex] = x;
        yb[editDragIndex] = y;
        renderEditAirfoil();
      }
    }, { passive: false });

    canvas.addEventListener('touchend', () => {
      if (airfoilTouchTimer) {
        clearTimeout(airfoilTouchTimer);
        airfoilTouchTimer = null;
      }
      if (!editMode) return;
      if (touchAirfoilDrag && editDragIndex != null) {
        editDragIndex = null;
        editHoverIndex = null;
        touchAirfoilDrag = false;
        applyEditedAirfoil();
      }
    });
  }
}

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
updateDownloadButtons();
update();
