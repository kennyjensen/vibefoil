import fs from 'fs';

import { cft } from '../../js/blu.js';

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));

const results = {
  cft: payload.cft.map((item) => cft(item.hk, item.rt, item.msq, item.cffac ?? 1.0)),
};

process.stdout.write(JSON.stringify({ results }));
