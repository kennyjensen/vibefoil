import fs from 'fs';

import { atanc, setexp } from '../../js/xutils.js';

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));

const results = {
  atanc: payload.atanc.map((item) => atanc(item.y, item.x, item.thold)),
  setexp: payload.setexp.map((item) => {
    const s = new Float64Array(item.nn);
    setexp(s, item.ds1, item.smax, item.nn);
    return Array.from(s);
  }),
};

process.stdout.write(JSON.stringify({ results }));
