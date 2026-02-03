import fs from 'fs';

import {
  apcalc,
  xywake,
  qwcalc,
  ncalc,
  ggcalc,
  psilin,
  pswlin,
  qdcalc,
  stfind,
  iblpan,
  xicalc,
  uicalc,
  qvfue,
  qiset,
  gamqv,
  stmove,
} from '../../js/xpanel.js';

const input = fs.readFileSync(0, 'utf8');
const payload = JSON.parse(input);

function buildPanelCtx(state) {
  const n = state.n ?? 0;
  const nw = state.nw ?? 0;
  const total = n + nw;
  const makeArray = (src, len, fill = 0.0) =>
    Float64Array.from(src ?? new Array(len).fill(fill));
  return {
    N: n,
    NW: nw,
    WAKLEN: state.waklen ?? 1.0,
    CHORD: state.chord ?? 1.0,
    X: makeArray(state.x, total),
    Y: makeArray(state.y, total),
    XP: makeArray(state.xp, total),
    YP: makeArray(state.yp, total),
    S: makeArray(state.s, total),
    NX: makeArray(state.nx, total),
    NY: makeArray(state.ny, total),
    APANEL: makeArray(state.apanel, total),
    SHARP: !!state.sharp,
    PI: state.pi ?? Math.PI,
    ANTE: state.ante ?? 0.0,
    ASTE: state.aste ?? 0.0,
    DSTE: state.dste ?? 1.0,
    GAM: makeArray(state.gam, n),
    GAM_A: makeArray(state.gamA, n),
    QINV: makeArray(state.qinv, total + 1),
    QINV_A: makeArray(state.qinvA, total + 1),
    QVIS: makeArray(state.qvis, total + 1),
    GAMU: Array.from({ length: n + 1 }, () => new Float64Array(2)),
    SIG: makeArray(state.sig, total),
    QF0: makeArray(state.qf0, n),
    QF1: makeArray(state.qf1, n),
    QF2: makeArray(state.qf2, n),
    QF3: makeArray(state.qf3, n),
    QINVU: Array.from({ length: total }, (_, idx) => Float64Array.from(state.qinvu?.[idx] ?? [0.0, 0.0])),
    AIJ: Array.from({ length: n + 1 }, () => new Float64Array(n + 1)),
    BIJ: Array.from({ length: n + 1 }, () => new Float64Array(total)),
    AIJPIV: new Int32Array(n + 1),
    QOPI: state.qopi ?? 1.0 / (4.0 * Math.PI),
    HOPI: state.hopi ?? 1.0 / (2.0 * Math.PI),
    ALFA: state.alfa ?? 0.0,
    QINF: state.qinf ?? 1.0,
    LIMAGE: !!state.limage,
    YIMAGE: state.yimage ?? 0.0,
    XTE: state.xte ?? 0.0,
    YTE: state.yte ?? 0.0,
    DZDG: new Float64Array(n),
    DZDN: new Float64Array(n),
    DQDG: new Float64Array(n),
    DZDM: new Float64Array(total),
    DQDM: new Float64Array(total),
    LWAKE: false,
    LWDIJ: false,
    LADIJ: false,
    SNEW: new Float64Array(total),
  };
}

function buildBlCtx(state, panelCtx) {
  const nbl1 = state.nbl1 ?? 0;
  const nbl2 = state.nbl2 ?? 0;
  const nb = panelCtx?.N ?? state.n ?? 0;
  const nw = state.nw ?? panelCtx?.NW ?? 0;
  const maxNbl = Math.max(nbl1, nbl2, nb) + nw + 2;
  const make2D = (rows, cols, fill = 0.0) =>
    Array.from({ length: rows }, () => {
      const arr = new Float64Array(cols + 1);
      if (fill !== 0.0) arr.fill(fill);
      return arr;
    });
  const ipan = Array.from({ length: maxNbl + 1 }, () => new Int32Array(3));
  const vti = make2D(maxNbl + 1, 2);
  const uinv = make2D(maxNbl + 1, 2);
  const uinvA = make2D(maxNbl + 1, 2);
  const uedg = make2D(maxNbl + 1, 2);
  const ctau = make2D(maxNbl + 1, 2);
  const thet = make2D(maxNbl + 1, 2);
  const dstr = make2D(maxNbl + 1, 2);
  const xssi = make2D(maxNbl + 1, 2);
  const mass = make2D(maxNbl + 1, 2);
  const wgap = new Float64Array((state.nw ?? 0) + 2);
  const isys = Array.from({ length: maxNbl + 1 }, () => new Int32Array(3));
  const itran = new Int32Array(3);
  itran[1] = state.itran1 ?? 0;
  itran[2] = state.itran2 ?? 0;

  if (state.ipan1) {
    state.ipan1.forEach((val, idx) => {
      if (idx + 1 < ipan.length) ipan[idx + 1][1] = val;
    });
  }
  if (state.ipan2) {
    state.ipan2.forEach((val, idx) => {
      if (idx + 1 < ipan.length) ipan[idx + 1][2] = val;
    });
  }
  if (state.vti1) {
    state.vti1.forEach((val, idx) => {
      if (idx + 1 < vti.length) vti[idx + 1][1] = val;
    });
  }
  if (state.vti2) {
    state.vti2.forEach((val, idx) => {
      if (idx + 1 < vti.length) vti[idx + 1][2] = val;
    });
  }
  if (state.uedg1) {
    state.uedg1.forEach((val, idx) => {
      if (idx + 1 < uedg.length) uedg[idx + 1][1] = val;
    });
  }
  if (state.uedg2) {
    state.uedg2.forEach((val, idx) => {
      if (idx + 1 < uedg.length) uedg[idx + 1][2] = val;
    });
  }
  if (state.ctau1) {
    state.ctau1.forEach((val, idx) => {
      if (idx + 1 < ctau.length) ctau[idx + 1][1] = val;
    });
  }
  if (state.ctau2) {
    state.ctau2.forEach((val, idx) => {
      if (idx + 1 < ctau.length) ctau[idx + 1][2] = val;
    });
  }
  if (state.thet1) {
    state.thet1.forEach((val, idx) => {
      if (idx + 1 < thet.length) thet[idx + 1][1] = val;
    });
  }
  if (state.thet2) {
    state.thet2.forEach((val, idx) => {
      if (idx + 1 < thet.length) thet[idx + 1][2] = val;
    });
  }
  if (state.dstr1) {
    state.dstr1.forEach((val, idx) => {
      if (idx + 1 < dstr.length) dstr[idx + 1][1] = val;
    });
  }
  if (state.dstr2) {
    state.dstr2.forEach((val, idx) => {
      if (idx + 1 < dstr.length) dstr[idx + 1][2] = val;
    });
  }
  if (state.xssi1) {
    state.xssi1.forEach((val, idx) => {
      if (idx + 1 < xssi.length) xssi[idx + 1][1] = val;
    });
  }
  if (state.xssi2) {
    state.xssi2.forEach((val, idx) => {
      if (idx + 1 < xssi.length) xssi[idx + 1][2] = val;
    });
  }

  const blCtx = {
    IST: state.ist ?? 1,
    SST: state.sst ?? 0.0,
    SST_GO: state.sstGo ?? 0.0,
    SST_GP: state.sstGp ?? 0.0,
    IBLTE: new Int32Array([0, state.iblte1 ?? 0, state.iblte2 ?? 0]),
    NBL: new Int32Array([0, nbl1, nbl2]),
    IPAN: ipan,
    VTI: vti,
    UINV: uinv,
    UINV_A: uinvA,
    UEDG: uedg,
    CTAU: ctau,
    THET: thet,
    DSTR: dstr,
    XSSI: xssi,
    MASS: mass,
    WGAP: wgap,
    ISYS: isys,
    ITRAN: itran,
    upperIdx: [],
    lowerIdx: [],
  };
  return blCtx;
}

const results = payload.cases.map((caseDef) => {
  const state = caseDef.state ?? {};
  if (caseDef.kind === 'ncalc') {
    const { x, y, s } = caseDef;
    const n = x.length;
    const xArr = Float64Array.from(x);
    const yArr = Float64Array.from(y);
    const sArr = Float64Array.from(s);
    const xn = new Float64Array(n);
    const yn = new Float64Array(n);
    ncalc(xArr, yArr, sArr, n, xn, yn);
    return {
      kind: 'ncalc',
      xn: Array.from(xn),
      yn: Array.from(yn),
    };
  }

  const panelCtx = buildPanelCtx(state);
  const blCtx = buildBlCtx(state, panelCtx);

  switch (caseDef.kind) {
    case 'apcalc':
      apcalc(panelCtx);
      return { kind: 'apcalc', apanel: Array.from(panelCtx.APANEL) };
    case 'psilin': {
      const res = psilin(
        caseDef.i,
        caseDef.xi,
        caseDef.yi,
        caseDef.nxi,
        caseDef.nyi,
        caseDef.geolin,
        caseDef.siglin,
        panelCtx
      );
      return {
        kind: 'psilin',
        result: res,
        dzdg: Array.from(panelCtx.DZDG),
        dzdn: Array.from(panelCtx.DZDN),
        dzdm: Array.from(panelCtx.DZDM),
        dqdg: Array.from(panelCtx.DQDG),
        dqdm: Array.from(panelCtx.DQDM),
      };
    }
    case 'pswlin': {
      const res = pswlin(
        caseDef.i,
        caseDef.xi,
        caseDef.yi,
        caseDef.nxi,
        caseDef.nyi,
        panelCtx
      );
      return {
        kind: 'pswlin',
        result: res,
        dzdm: Array.from(panelCtx.DZDM),
        dqdm: Array.from(panelCtx.DQDM),
      };
    }
    case 'ggcalc':
      ggcalc(panelCtx);
      return {
        kind: 'ggcalc',
        qinvu: panelCtx.QINVU.map((row) => [row[0], row[1]]),
      };
    case 'qwcalc':
      qwcalc(panelCtx);
      return {
        kind: 'qwcalc',
        qinvu: panelCtx.QINVU.map((row) => [row[0], row[1]]),
      };
    case 'qdcalc':
      ggcalc(panelCtx);
      qdcalc(panelCtx);
      return {
        kind: 'qdcalc',
        dij: panelCtx.DIJ.map((row) => Array.from(row)),
      };
    case 'xywake':
      xywake(panelCtx);
      return {
        kind: 'xywake',
        x: Array.from(panelCtx.X),
        y: Array.from(panelCtx.Y),
        s: Array.from(panelCtx.S),
        nx: Array.from(panelCtx.NX),
        ny: Array.from(panelCtx.NY),
        apanel: Array.from(panelCtx.APANEL),
      };
    case 'stfind': {
      const res = stfind(panelCtx, panelCtx.N);
      return { kind: 'stfind', result: res };
    }
    case 'iblpan':
      iblpan(blCtx, panelCtx.N, panelCtx.NW ?? 0);
      return {
        kind: 'iblpan',
        iblte: [blCtx.IBLTE[1], blCtx.IBLTE[2]],
        nbl: [blCtx.NBL[1], blCtx.NBL[2]],
        ipan1: Array.from(blCtx.IPAN).slice(1, blCtx.NBL[1] + 1).map((row) => row[1]),
        ipan2: Array.from(blCtx.IPAN).slice(1, blCtx.NBL[2] + 1).map((row) => row[2]),
        vti1: Array.from(blCtx.VTI).slice(1, blCtx.NBL[1] + 1).map((row) => row[1]),
        vti2: Array.from(blCtx.VTI).slice(1, blCtx.NBL[2] + 1).map((row) => row[2]),
      };
    case 'xicalc':
      xicalc(blCtx, panelCtx);
      return {
        kind: 'xicalc',
        xssi1: Array.from(blCtx.XSSI).slice(1, blCtx.NBL[1] + 1).map((row) => row[1]),
        xssi2: Array.from(blCtx.XSSI).slice(1, blCtx.NBL[2] + 1).map((row) => row[2]),
        wgap: Array.from(blCtx.WGAP ?? new Float64Array(panelCtx.NW ?? 0)),
      };
    case 'uicalc':
      uicalc(blCtx, panelCtx.QINV, panelCtx.QINV_A);
      return {
        kind: 'uicalc',
        uinv1: Array.from(blCtx.UINV).slice(1, blCtx.NBL[1] + 1).map((row) => row[1]),
        uinv2: Array.from(blCtx.UINV).slice(1, blCtx.NBL[2] + 1).map((row) => row[2]),
        uinvA1: Array.from(blCtx.UINV_A).slice(1, blCtx.NBL[1] + 1).map((row) => row[1]),
        uinvA2: Array.from(blCtx.UINV_A).slice(1, blCtx.NBL[2] + 1).map((row) => row[2]),
      };
    case 'qvfue':
      qvfue(blCtx, panelCtx.QVIS);
      return {
        kind: 'qvfue',
        qvis: Array.from(panelCtx.QVIS),
      };
    case 'qiset':
      qiset(panelCtx, panelCtx.ALFA);
      return {
        kind: 'qiset',
        qinv: Array.from(panelCtx.QINV),
        qinvA: Array.from(panelCtx.QINV_A),
      };
    case 'gamqv':
      gamqv(panelCtx, panelCtx.QVIS, panelCtx.QINV_A);
      return {
        kind: 'gamqv',
        gam: Array.from(panelCtx.GAM),
        gamA: Array.from(panelCtx.GAM_A),
      };
    case 'stmove':
      stmove(panelCtx, blCtx, panelCtx.QINV, panelCtx.QINV_A);
      return {
        kind: 'stmove',
        ist: blCtx.IST,
        sst: blCtx.SST,
        sstGo: blCtx.SST_GO,
        sstGp: blCtx.SST_GP,
        nbl: [blCtx.NBL[1], blCtx.NBL[2]],
        uegd1: Array.from(blCtx.UEDG).slice(1, blCtx.NBL[1] + 1).map((row) => row[1]),
        uegd2: Array.from(blCtx.UEDG).slice(1, blCtx.NBL[2] + 1).map((row) => row[2]),
        ctau1: Array.from(blCtx.CTAU).slice(1, blCtx.NBL[1] + 1).map((row) => row[1]),
        ctau2: Array.from(blCtx.CTAU).slice(1, blCtx.NBL[2] + 1).map((row) => row[2]),
        thet1: Array.from(blCtx.THET).slice(1, blCtx.NBL[1] + 1).map((row) => row[1]),
        thet2: Array.from(blCtx.THET).slice(1, blCtx.NBL[2] + 1).map((row) => row[2]),
        dstr1: Array.from(blCtx.DSTR).slice(1, blCtx.NBL[1] + 1).map((row) => row[1]),
        dstr2: Array.from(blCtx.DSTR).slice(1, blCtx.NBL[2] + 1).map((row) => row[2]),
        xssi1: Array.from(blCtx.XSSI).slice(1, blCtx.NBL[1] + 1).map((row) => row[1]),
        xssi2: Array.from(blCtx.XSSI).slice(1, blCtx.NBL[2] + 1).map((row) => row[2]),
        mass1: Array.from(blCtx.MASS).slice(1, blCtx.NBL[1] + 1).map((row) => row[1]),
        mass2: Array.from(blCtx.MASS).slice(1, blCtx.NBL[2] + 1).map((row) => row[2]),
      };
    default:
      return { kind: caseDef.kind, error: 'unknown' };
  }
});

process.stdout.write(JSON.stringify({ results }));
