import fs from 'fs';

import { SPLQSP, SMOOQ, SYMQSP, GAMQSP, QINCOM } from '../../js/xqdes.js';
import { createMatrix } from '../../js/arrays.js';

const noop = () => {};
console.log = noop;
console.warn = noop;

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));
const {
  nsp,
  sspec,
  qspec,
  qgamm,
  kq1,
  kq2,
  lqslop,
  algam,
  clgam,
  cmgam,
  liqset,
  qincomCases,
} = payload;

function baseCtx() {
  const ctx = {
    NSP: nsp,
    SSPEC: new Float64Array(nsp),
    QSPEC: createMatrix(nsp, 1),
    QSPECP: createMatrix(nsp, 1),
    W1: new Float64Array(nsp),
    W2: new Float64Array(nsp),
    W3: new Float64Array(nsp),
    W8: new Float64Array(nsp),
    LQSLOP: lqslop,
    QDOF0: 0.0,
    QDOF1: 0.0,
    QDOF2: 0.0,
    QDOF3: 0.0,
    ALQSP: new Float64Array(1),
    CLQSP: new Float64Array(1),
    CMQSP: new Float64Array(1),
    QGAMM: new Float64Array(nsp),
    ALGAM: algam,
    CLGAM: clgam,
    CMGAM: cmgam,
    LIQSET: liqset,
    IQ1: 0,
    IQ2: nsp - 1,
  };

  for (let i = 0; i < nsp; i += 1) {
    ctx.SSPEC[i] = sspec[i];
    ctx.QSPEC[i][0] = qspec[i];
    ctx.QGAMM[i] = qgamm[i];
  }

  return ctx;
}

const ctxSpl = baseCtx();
SPLQSP(ctxSpl, 0);
const splqsp = {
  qspecp: Array.from({ length: nsp }, (_, i) => ctxSpl.QSPECP[i][0]),
};

const ctxSmo = baseCtx();
SMOOQ(ctxSmo, kq1, kq2, 0);
const smooq = {
  qspec: Array.from({ length: nsp }, (_, i) => ctxSmo.QSPEC[i][0]),
};

const ctxSym = baseCtx();
SYMQSP(ctxSym, 0);
const symqsp = {
  sspec: Array.from({ length: nsp }, (_, i) => ctxSym.SSPEC[i]),
  qspec: Array.from({ length: nsp }, (_, i) => ctxSym.QSPEC[i][0]),
};

const ctxGam = baseCtx();
GAMQSP(ctxGam, 0);
const gamqsp = {
  alqsp: ctxGam.ALQSP[0],
  clqsp: ctxGam.CLQSP[0],
  cmqsp: ctxGam.CMQSP[0],
  qspec: Array.from({ length: nsp }, (_, i) => ctxGam.QSPEC[i][0]),
  qdof0: ctxGam.QDOF0,
  qdof1: ctxGam.QDOF1,
  qdof2: ctxGam.QDOF2,
  qdof3: ctxGam.QDOF3,
  iq1: ctxGam.IQ1,
  iq2: ctxGam.IQ2,
};

const qincom = qincomCases.map((item) => QINCOM(item.qc, item.qinf, item.tklam));

process.stdout.write(JSON.stringify({
  results: { splqsp, smooq, symqsp, gamqsp, qincom },
}));
