import fs from 'fs';

import {
  hkin,
  dil,
  dilw,
  hsl,
  cfl,
  dit,
  hst,
  cft,
  hct,
} from '../../js/xblsys.js';
import { dslim } from '../../js/xbl.js';

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));

const ctx = { hkin };

const results = {
  hkin: payload.hkin.map(({ h, msq }) => hkin(h, msq)),
  dil: payload.dil.map(({ hk, rt }) => dil(hk, rt)),
  dilw: payload.dilw.map(({ hk, rt }) => dilw(hk, rt)),
  hsl: payload.hsl.map(({ hk, rt, msq }) => hsl(hk, rt, msq)),
  cfl: payload.cfl.map(({ hk, rt, msq }) => cfl(hk, rt, msq)),
  dit: payload.dit.map(({ hs, us, cf, st }) => dit(hs, us, cf, st)),
  hst: payload.hst.map(({ hk, rt, msq }) => hst(hk, rt, msq)),
  cft: payload.cft.map(({ hk, rt, msq, cffac }) => cft(hk, rt, msq, cffac)),
  hct: payload.hct.map(({ hk, msq }) => hct(hk, msq)),
  dslim: payload.dslim.map(({ dstr, thet, uedg, msq, hklim }) => dslim(ctx, dstr, thet, uedg, msq, hklim)),
};

process.stdout.write(JSON.stringify({ results }));
