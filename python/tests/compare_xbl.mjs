import fs from 'fs';

import {
  dslim,
  blpini,
  iblsys,
  mrcl,
  comset,
  mrchue,
  mrchdu,
  ueset,
  xifset,
  update,
} from '../../js/xbl.js';
import {
  hkin,
  ensureCtx,
  syncComToVars,
  blprv,
  blkin,
  blsys,
  tesys,
  blvar,
  blmid,
  trchek,
} from '../../js/xblsys.js';

const input = fs.readFileSync(0, 'utf8');
const payload = JSON.parse(input);

function make2D(rows, cols, fill = 0.0) {
  return Array.from({ length: rows }, () => {
    const row = new Float64Array(cols + 1);
    if (fill !== 0.0) row.fill(fill);
    return row;
  });
}

function buildCtx(state) {
  const n = state.n ?? 0;
  const nbl1 = state.nbl1 ?? 0;
  const nbl2 = state.nbl2 ?? 0;
  const maxNbl = Math.max(nbl1, nbl2) + 2;
  const ctx = {
    N: n,
    NBL: new Int32Array([0, nbl1, nbl2]),
    IBLTE: new Int32Array([0, state.iblte1 ?? 0, state.iblte2 ?? 0]),
    ISYS: Array.from({ length: maxNbl + 1 }, () => new Int32Array(3)),
    IVX: state.ivx ?? 200,
    IQX: state.iqx ?? state.n ?? 0,
    MATYP: state.matyp ?? 1,
    RETYP: state.retyp ?? 1,
    MINF1: state.minf1 ?? 0.0,
    REINF1: state.reinf1 ?? 1.0,
    MINF: state.minf ?? 0.0,
    TKLAM: 0.0,
    TK_MSQ: 0.0,
    hkin,
    XLE: state.xle ?? 0.0,
    YLE: state.yle ?? 0.0,
    XTE: state.xte ?? 1.0,
    YTE: state.yte ?? 0.0,
    SLE: state.sle ?? 0.0,
    SST: state.sst ?? 0.0,
    XSTRIP: new Float64Array([0.0, state.xstrip1 ?? 1.0, state.xstrip2 ?? 1.0]),
    ITRAN: new Int32Array([0, state.itran1 ?? 0, state.itran2 ?? 0]),
    ACRIT: new Float64Array([0.0, state.acrit1 ?? 9.0, state.acrit2 ?? 9.0]),
    TFORCE: new Array(3).fill(false),
    XSSITR: new Float64Array(3),
    REYBL: state.reybl ?? 1.0e6,
    GM1BL: state.gm1bl ?? 0.4,
    HSTINV: state.hstinv ?? 0.0,
    HVRAT: state.hvrat ?? 0.0,
    ANTE: state.ante ?? 0.0,
  };

  if (n > 0) {
    const arr = (src, len, fill = 0.0) => {
      const out = new Float64Array(len);
      if (src) {
        const max = Math.min(src.length, len);
        for (let i = 0; i < max; i += 1) out[i] = src[i];
        if (max < len && fill !== 0.0) out.fill(fill, max);
      } else if (fill !== 0.0) {
        out.fill(fill);
      }
      return out;
    };
    ctx.X = arr(state.x, n + 1);
    ctx.Y = arr(state.y, n + 1);
    ctx.S = arr(state.s, n + 1);
    ctx.W1 = arr(state.w1, n + 1);
    ctx.W2 = arr(state.w2, n + 1);
    ctx.W3 = arr(state.w3, n + 1);
    ctx.W4 = arr(state.w4, n + 1);
  }

  if (nbl1 > 0 || nbl2 > 0) {
    ctx.IPAN = Array.from({ length: maxNbl + 1 }, () => new Int32Array(3));
    ctx.VTI = make2D(maxNbl + 1, 2);
    ctx.UINV = make2D(maxNbl + 1, 2);
    ctx.UINV_A = make2D(maxNbl + 1, 2);
    ctx.MASS = make2D(maxNbl + 1, 2);
    ctx.UEDG = make2D(maxNbl + 1, 2);
    ctx.XSSI = make2D(maxNbl + 1, 2);
    ctx.THET = make2D(maxNbl + 1, 2);
    ctx.DSTR = make2D(maxNbl + 1, 2);
    ctx.CTAU = make2D(maxNbl + 1, 2);
    ctx.TAU = make2D(maxNbl + 1, 2);
    ctx.DIS = make2D(maxNbl + 1, 2);
    ctx.CTQ = make2D(maxNbl + 1, 2);
    ctx.DELT = make2D(maxNbl + 1, 2);
    ctx.TSTR = make2D(maxNbl + 1, 2);
    ctx.WGAP = new Float64Array((state.nw ?? 0) + 2);
    ctx.HTARG = [null, new Float64Array(maxNbl + 1), new Float64Array(maxNbl + 1)];
  }

  if (state.ipan1) {
    state.ipan1.forEach((val, idx) => {
      ctx.IPAN[idx + 2][1] = val;
    });
  }
  if (state.ipan2) {
    state.ipan2.forEach((val, idx) => {
      ctx.IPAN[idx + 2][2] = val;
    });
  }
  if (state.vti1) {
    state.vti1.forEach((val, idx) => {
      ctx.VTI[idx + 2][1] = val;
    });
  }
  if (state.vti2) {
    state.vti2.forEach((val, idx) => {
      ctx.VTI[idx + 2][2] = val;
    });
  }
  if (state.uinv1) {
    state.uinv1.forEach((val, idx) => {
      ctx.UINV[idx + 2][1] = val;
    });
  }
  if (state.uinv2) {
    state.uinv2.forEach((val, idx) => {
      ctx.UINV[idx + 2][2] = val;
    });
  }
  if (state.uinvA1) {
    state.uinvA1.forEach((val, idx) => {
      ctx.UINV_A[idx + 2][1] = val;
    });
  }
  if (state.uinvA2) {
    state.uinvA2.forEach((val, idx) => {
      ctx.UINV_A[idx + 2][2] = val;
    });
  }
  if (state.mass1) {
    state.mass1.forEach((val, idx) => {
      ctx.MASS[idx + 2][1] = val;
    });
  }
  if (state.mass2) {
    state.mass2.forEach((val, idx) => {
      ctx.MASS[idx + 2][2] = val;
    });
  }
  if (state.uedg1) {
    state.uedg1.forEach((val, idx) => {
      ctx.UEDG[idx + 2][1] = val;
    });
  }
  if (state.uedg2) {
    state.uedg2.forEach((val, idx) => {
      ctx.UEDG[idx + 2][2] = val;
    });
  }
  if (state.thet1) {
    state.thet1.forEach((val, idx) => {
      ctx.THET[idx + 2][1] = val;
    });
  }
  if (state.thet2) {
    state.thet2.forEach((val, idx) => {
      ctx.THET[idx + 2][2] = val;
    });
  }
  if (state.dstr1) {
    state.dstr1.forEach((val, idx) => {
      ctx.DSTR[idx + 2][1] = val;
    });
  }
  if (state.dstr2) {
    state.dstr2.forEach((val, idx) => {
      ctx.DSTR[idx + 2][2] = val;
    });
  }
  if (state.ctau1) {
    state.ctau1.forEach((val, idx) => {
      ctx.CTAU[idx + 2][1] = val;
    });
  }
  if (state.ctau2) {
    state.ctau2.forEach((val, idx) => {
      ctx.CTAU[idx + 2][2] = val;
    });
  }
  if (state.xssi1) {
    state.xssi1.forEach((val, idx) => {
      ctx.XSSI[idx + 1][1] = val;
    });
  }
  if (state.xssi2) {
    state.xssi2.forEach((val, idx) => {
      ctx.XSSI[idx + 1][2] = val;
    });
  }

  if (state.itran1 != null || state.itran2 != null) {
    ctx.ITRAN = new Int32Array([0, state.itran1 ?? 0, state.itran2 ?? 0]);
  }

  if (state.dij) {
    ctx.DIJ = Array.from({ length: n + 1 }, (_, i) => {
      const row = new Float64Array(n + 1);
      if (state.dij[i]) {
        state.dij[i].forEach((val, j) => {
          if (j <= n) row[j] = val;
        });
      }
      return row;
    });
  }

  if (state.vdel) {
    const vdel = Array.from({ length: 4 }, () =>
      Array.from({ length: 3 }, () => new Float64Array((state.nsys ?? 0) + 1)));
    for (let k = 1; k <= 3; k += 1) {
      for (let j = 1; j <= 2; j += 1) {
        const row = state.vdel?.[k]?.[j] ?? [];
        for (let iv = 1; iv <= (state.nsys ?? 0); iv += 1) {
          vdel[k][j][iv] = row[iv] ?? 0.0;
        }
      }
    }
    ctx.VDEL = vdel;
  }

  if (state.isys1 || state.isys2) {
    if (!ctx.ISYS) ctx.ISYS = Array.from({ length: maxNbl + 1 }, () => new Int32Array(3));
    if (state.isys1) {
      state.isys1.forEach((val, idx) => {
        ctx.ISYS[idx + 2][1] = val;
      });
    }
    if (state.isys2) {
      state.isys2.forEach((val, idx) => {
        ctx.ISYS[idx + 2][2] = val;
      });
    }
  }

  if (state.nsys != null) {
    ctx.NSYS = state.nsys;
  }

  ctx.MATYP = state.matyp ?? ctx.MATYP ?? 1;
  ctx.LALFA = state.lalfa ?? ctx.LALFA ?? false;
  ctx.CL = state.cl ?? ctx.CL ?? 0.0;
  ctx.CLSPEC = state.clspec ?? ctx.CLSPEC ?? 0.0;
  ctx.ALFA = state.alfa ?? ctx.ALFA ?? 0.0;
  ctx.DTOR = state.dtor ?? ctx.DTOR ?? Math.PI / 180.0;
  ctx.MINF = state.minf ?? ctx.MINF ?? 0.0;
  ctx.QINF = state.qinf ?? ctx.QINF ?? 1.0;
  ctx.GAMMA = state.gamma ?? ctx.GAMMA ?? 1.4;
  ctx.GAMM1 = state.gamm1 ?? ctx.GAMM1 ?? 0.4;
  ctx.MINF_CL = state.minfCl ?? ctx.MINF_CL ?? 0.0;
  ctx.HVRAT = state.hvrat ?? ctx.HVRAT ?? 0.0;
  ctx.QINV_A = ctx.UINV_A;

  ensureCtx(ctx);
  ctx.syncComToVars = syncComToVars;
  ctx.blprv = blprv;
  ctx.blkin = blkin;
  ctx.blsys = blsys;
  ctx.tesys = tesys;
  ctx.blvar = blvar;
  ctx.blmid = blmid;
  ctx.trchek = trchek;
  if (Number.isFinite(ctx.MINF)) {
    comset(ctx);
    ctx.GAMBL = ctx.GAMMA;
    ctx.GM1BL = ctx.GAMM1;
    ctx.QINFBL = ctx.QINF;
    ctx.TKBL = ctx.TKLAM;
    ctx.TKBL_MS = ctx.TK_MSQ;
    const minf = ctx.MINF;
    const gm1 = ctx.GM1BL;
    const qinf = ctx.QINFBL;
    if (!Number.isFinite(ctx.RSTBL)) {
      ctx.RSTBL = (1.0 + 0.5 * gm1 * minf ** 2) ** (1.0 / gm1);
      ctx.RSTBL_MS = 0.5 * ctx.RSTBL / (1.0 + 0.5 * gm1 * minf ** 2);
    }
    if (!Number.isFinite(ctx.HSTINV) || ctx.HSTINV === 0.0) {
      ctx.HSTINV = gm1 * (minf / qinf) ** 2 / (1.0 + 0.5 * gm1 * minf ** 2);
    }
    ctx.HSTINV_MS = gm1 * (1.0 / qinf) ** 2 / (1.0 + 0.5 * gm1 * minf ** 2)
      - 0.5 * gm1 * ctx.HSTINV / (1.0 + 0.5 * gm1 * minf ** 2);
    const herat = 1.0 - 0.5 * qinf ** 2 * ctx.HSTINV;
    const heratMs = -0.5 * qinf ** 2 * ctx.HSTINV_MS;
    const reinf = state.reinf ?? ctx.REINF ?? ctx.REINF1 ?? 1.0e6;
    if (!Number.isFinite(ctx.REYBL) || ctx.REYBL === 0.0) {
      ctx.REYBL = reinf * Math.sqrt(herat ** 3) * (1.0 + ctx.HVRAT) / (herat + ctx.HVRAT);
    }
    ctx.REYBL_RE = Math.sqrt(herat ** 3) * (1.0 + ctx.HVRAT) / (herat + ctx.HVRAT);
    ctx.REYBL_MS = ctx.REYBL * (1.5 / herat - 1.0 / (herat + ctx.HVRAT)) * heratMs;
  }
  return ctx;
}

const results = payload.cases.map((caseDef) => {
  const state = caseDef.state ?? {};
  const ctx = buildCtx(state);

  switch (caseDef.kind) {
    case 'dslim': {
      const dstr = dslim(ctx, caseDef.dstr, caseDef.thet, caseDef.uedg, caseDef.msq, caseDef.hklim);
      return { kind: 'dslim', dstr };
    }
    case 'blpini':
      blpini(ctx);
      return {
        kind: 'blpini',
        constants: {
          SCCON: ctx.SCCON,
          GACON: ctx.GACON,
          GBCON: ctx.GBCON,
          GCCON: ctx.GCCON,
          DLCON: ctx.DLCON,
          CTRCON: ctx.CTRCON,
          CTRCEX: ctx.CTRCEX,
          DUXCON: ctx.DUXCON,
          CTCON: ctx.CTCON,
          CFFAC: ctx.CFFAC,
        },
      };
    case 'mrcl': {
      const res = mrcl(ctx, caseDef.cls);
      return { kind: 'mrcl', minf: ctx.MINF, reinf: ctx.REINF, mCls: res.mCls, rCls: res.rCls };
    }
    case 'comset':
      comset(ctx);
      return { kind: 'comset', tklam: ctx.TKLAM, tkMsq: ctx.TK_MSQ };
    case 'iblsys':
      iblsys(ctx);
      return {
        kind: 'iblsys',
        nsys: ctx.NSYS,
        isys1: Array.from(ctx.ISYS).slice(2, ctx.NBL[1] + 1).map((row) => row[1]),
        isys2: Array.from(ctx.ISYS).slice(2, ctx.NBL[2] + 1).map((row) => row[2]),
      };
    case 'xifset':
      xifset(ctx, caseDef.is);
      return { kind: 'xifset', xiforc: ctx.XIFORC };
    case 'ueset':
      ueset(ctx);
      return {
        kind: 'ueset',
        uedg1: Array.from(ctx.UEDG).slice(2, ctx.NBL[1] + 1).map((row) => row[1]),
        uedg2: Array.from(ctx.UEDG).slice(2, ctx.NBL[2] + 1).map((row) => row[2]),
      };
    case 'mrchue':
      blpini(ctx);
      {
        const prevLog = console.log;
        console.log = () => {};
        try {
          mrchue(ctx);
        } finally {
          console.log = prevLog;
        }
      }
      return {
        kind: 'mrchue',
        itran1: ctx.ITRAN?.[1] ?? 0,
        itran2: ctx.ITRAN?.[2] ?? 0,
        thet1: Array.from(ctx.THET).slice(2, ctx.NBL[1] + 1).map((row) => row[1]),
        thet2: Array.from(ctx.THET).slice(2, ctx.NBL[2] + 1).map((row) => row[2]),
        dstr1: Array.from(ctx.DSTR).slice(2, ctx.NBL[1] + 1).map((row) => row[1]),
        dstr2: Array.from(ctx.DSTR).slice(2, ctx.NBL[2] + 1).map((row) => row[2]),
        ctau1: Array.from(ctx.CTAU).slice(2, ctx.NBL[1] + 1).map((row) => row[1]),
        ctau2: Array.from(ctx.CTAU).slice(2, ctx.NBL[2] + 1).map((row) => row[2]),
        mass1: Array.from(ctx.MASS).slice(2, ctx.NBL[1] + 1).map((row) => row[1]),
        mass2: Array.from(ctx.MASS).slice(2, ctx.NBL[2] + 1).map((row) => row[2]),
      };
    case 'mrchdu':
      blpini(ctx);
      {
        const prevLog = console.log;
        console.log = () => {};
        try {
          mrchdu(ctx);
        } finally {
          console.log = prevLog;
        }
      }
      return {
        kind: 'mrchdu',
        itran1: ctx.ITRAN?.[1] ?? 0,
        itran2: ctx.ITRAN?.[2] ?? 0,
        thet1: Array.from(ctx.THET).slice(2, ctx.NBL[1] + 1).map((row) => row[1]),
        thet2: Array.from(ctx.THET).slice(2, ctx.NBL[2] + 1).map((row) => row[2]),
        dstr1: Array.from(ctx.DSTR).slice(2, ctx.NBL[1] + 1).map((row) => row[1]),
        dstr2: Array.from(ctx.DSTR).slice(2, ctx.NBL[2] + 1).map((row) => row[2]),
        ctau1: Array.from(ctx.CTAU).slice(2, ctx.NBL[1] + 1).map((row) => row[1]),
        ctau2: Array.from(ctx.CTAU).slice(2, ctx.NBL[2] + 1).map((row) => row[2]),
        mass1: Array.from(ctx.MASS).slice(2, ctx.NBL[1] + 1).map((row) => row[1]),
        mass2: Array.from(ctx.MASS).slice(2, ctx.NBL[2] + 1).map((row) => row[2]),
      };
    case 'update':
      update(ctx);
      return {
        kind: 'update',
        cl: ctx.CL,
        alfa: ctx.ALFA,
        rmsbl: ctx.RMSBL,
        rlx: ctx.RLX,
        uedg1: Array.from(ctx.UEDG).slice(2, ctx.NBL[1] + 1).map((row) => row[1]),
        uedg2: Array.from(ctx.UEDG).slice(2, ctx.NBL[2] + 1).map((row) => row[2]),
        thet1: Array.from(ctx.THET).slice(2, ctx.NBL[1] + 1).map((row) => row[1]),
        thet2: Array.from(ctx.THET).slice(2, ctx.NBL[2] + 1).map((row) => row[2]),
        dstr1: Array.from(ctx.DSTR).slice(2, ctx.NBL[1] + 1).map((row) => row[1]),
        dstr2: Array.from(ctx.DSTR).slice(2, ctx.NBL[2] + 1).map((row) => row[2]),
        ctau1: Array.from(ctx.CTAU).slice(2, ctx.NBL[1] + 1).map((row) => row[1]),
        ctau2: Array.from(ctx.CTAU).slice(2, ctx.NBL[2] + 1).map((row) => row[2]),
        mass1: Array.from(ctx.MASS).slice(2, ctx.NBL[1] + 1).map((row) => row[1]),
        mass2: Array.from(ctx.MASS).slice(2, ctx.NBL[2] + 1).map((row) => row[2]),
      };
    default:
      return { kind: caseDef.kind, error: 'unknown' };
  }
});

process.stdout.write(JSON.stringify({ results }));
