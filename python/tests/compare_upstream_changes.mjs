import fs from 'node:fs';
import '../../js/naca.js';
import { cft } from '../../js/xblsys.js';
import { cft as profileCft } from '../../js/blu.js';
import { segspld, segspl, scalc } from '../../js/spline.js';
import { lefind } from '../../js/xgeom.js';
import { apcalc, ncalc, psilin } from '../../js/xpanel.js';
import { tecalc } from '../../js/xfoil.js';
import { mrcl, comset } from '../../js/xbl.js';

const p = JSON.parse(fs.readFileSync(0, 'utf8'));
let result;
if (p.mode === 1) {
  const out = cft(...p.args);
  const profile = profileCft(...p.args);
  if (JSON.stringify(out) !== JSON.stringify(profile)) throw new Error('CFT copies disagree');
  result = [out.cf, out.cfHk, out.cfRt, out.cfMsq];
} else if (p.mode === 2) {
  const [ides, n, sharp] = p.args;
  const xx = new Float64Array(n), yt = new Float64Array(n), yc = new Float64Array(n);
  const x = new Float64Array(2 * n), y = new Float64Array(2 * n);
  const generate = ides <= 9999 ? globalThis.Naca.naca4 : globalThis.Naca.naca5;
  const { nb } = generate(ides, xx, yt, yc, n, x, y, sharp);
  result = Array.from({ length: nb }, (_, i) => [x[i], y[i]]).flat();
} else if (p.mode === 3) {
  const xs = new Float64Array(p.x.length);
  segspld(p.x, xs, p.s, p.x.length, ...p.slopes);
  result = Array.from(xs);
} else if (p.mode === 4) {
  const n = p.x.length, s = new Float64Array(n), xp = new Float64Array(n), yp = new Float64Array(n);
  scalc(p.x, p.y, s, n);
  segspl(p.x, xp, s, n);
  segspl(p.y, yp, s, n);
  result = [lefind(p.x, xp, p.y, yp, s, n)];
} else if (p.mode === 5) {
  const n = p.x.length;
  const [ip, sharp, limage, alfa, yimage, xi, yi, nxi, nyi] = p.args;
  const ctx = { N: n, CHORD: 1.0, PI: Math.PI, HOPI: 0.5 / Math.PI, QOPI: 0.25 / Math.PI,
    QINF: 1.0, ALFA: alfa, LIMAGE: Boolean(limage), YIMAGE: yimage,
    X: p.x, Y: p.y, GAMU: Array.from({ length: n }, (_, i) => [0.2 * (i + 1), -0.1 * (i + 1)]) };
  for (const key of ['S', 'XP', 'YP', 'NX', 'NY', 'APANEL', 'DZDG', 'DZDN', 'DQDG', 'DZDM', 'DQDM']) ctx[key] = new Float64Array(n);
  for (const [key, factor] of [['GAM', 0.1], ['SIG', 0.005], ['QF0', 0.03], ['QF1', -0.02], ['QF2', 0.01], ['QF3', -0.04]])
    ctx[key] = Float64Array.from({ length: n }, (_, i) => factor * (i + 1 - (key === 'GAM' ? 3 : 0)));
  scalc(ctx.X, ctx.Y, ctx.S, n);
  segspl(ctx.X, ctx.XP, ctx.S, n);
  segspl(ctx.Y, ctx.YP, ctx.S, n);
  ncalc(ctx.X, ctx.Y, ctx.S, n, ctx.NX, ctx.NY);
  tecalc(ctx);
  ctx.SHARP = Boolean(sharp);
  apcalc(ctx);
  const value = psilin(ip - 1, xi, yi, nxi, nyi, true, true, ctx);
  result = ['psi', 'psiNi', 'qtan1', 'qtan2', 'zQinf', 'zAlfa', 'zQdof0', 'zQdof1', 'zQdof2', 'zQdof3'].map(key => value[key]);
  for (let j = 0; j < n; j += 1) for (const key of ['DZDG', 'DZDN', 'DQDG', 'DZDM', 'DQDM']) result.push(ctx[key][j]);
} else if (p.mode === 6) {
  const [matyp, retyp, mach, reynolds, cl] = p.args;
  const ctx = { MATYP: matyp, RETYP: retyp, MINF1: mach, REINF1: reynolds,
    QINF: 1.0, GAMMA: 1.4, GAMM1: 1.4 - 1.0 };
  mrcl(ctx, cl);
  comset(ctx);
  result = ['MATYP', 'RETYP', 'MINF', 'REINF', 'MINF_CL', 'REINF_CL',
    'TKLAM', 'TKL_MSQ', 'CPSTAR', 'QSTAR'].map(key => ctx[key]);
} else if (p.mode === 7) {
  const n = p.x.length;
  const ctx = { N: n, X: p.x, Y: p.y, S: new Float64Array(n),
    XP: new Float64Array(n), YP: new Float64Array(n) };
  scalc(ctx.X, ctx.Y, ctx.S, n);
  segspl(ctx.X, ctx.XP, ctx.S, n);
  segspl(ctx.Y, ctx.YP, ctx.S, n);
  tecalc(ctx);
  result = [ctx.ANTE, ctx.ASTE, ctx.DSTE, Number(ctx.SHARP)];
}
process.stdout.write(JSON.stringify(result));
