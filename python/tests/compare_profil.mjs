import fs from 'fs';

import { prwall, uwall, fs as fsProfile } from '../../js/profil.js';

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));

const prwallResults = payload.prwall.map((item) => prwall(item));

const uwallResults = payload.uwall.map((item) => {
  const yy = new Array(item.n + 1).fill(0.0);
  const xx = new Array(item.n + 1).fill(0.0);
  const out = uwall({ ...item, yy, xx });
  return {
    y: out.y,
    u: out.u,
  };
});

const fsResults = payload.fs.map((item) => {
  const eta = new Array(item.n + 1).fill(0.0);
  const f = new Array(item.n + 1).fill(0.0);
  const u = new Array(item.n + 1).fill(0.0);
  const s = new Array(item.n + 1).fill(0.0);
  const out = fsProfile({ ...item, eta, f, u, s });
  return {
    eta,
    f,
    u,
    s,
    delta: out.delta,
  };
});

const results = {
  prwall: prwallResults,
  uwall: uwallResults,
  fs: fsResults,
};

process.stdout.write(JSON.stringify({ results }));
