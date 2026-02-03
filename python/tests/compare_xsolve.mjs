import fs from 'fs';

import { gauss, ludcmp, baksub, blsolv } from '../../js/xsolve.js';
import { createMatrix1, createTensor3 } from '../../js/arrays.js';

const payload = JSON.parse(fs.readFileSync(0, 'utf8'));

const results = {};

if (payload.gauss) {
  const { n, nrhs, z, r } = payload.gauss;
  const zWork = z.map((row) => Float64Array.from(row));
  const rWork = r.map((row) => row.slice());
  gauss(n, zWork, rWork, nrhs);
  results.gauss = {
    n,
    nrhs,
    z: zWork.map((row) => Array.from(row)),
    r: rWork.map((row) => Array.from(row)),
  };
}

if (payload.ludcmp) {
  const { n, a, b } = payload.ludcmp;
  const aWork = a.map((row) => Float64Array.from(row));
  const bWork = Float64Array.from(b);
  const indx = new Int32Array(n);
  ludcmp(n, aWork, indx);
  baksub(n, aWork, indx, bWork);
  results.ludcmp = {
    n,
    a: aWork.map((row) => Array.from(row)),
    indx: Array.from(indx),
    b: Array.from(bWork),
  };
}

if (payload.blsolv) {
  const {
    nsys,
    iblte1,
    iblte2,
    ivte1,
    ivz,
    vaccel,
    s1,
    s2,
    va,
    vb,
    vm,
    vdel,
    vz,
  } = payload.blsolv;

  const ctx = {
    NSYS: nsys,
    N: 2,
    VACCEL: vaccel,
    S: new Float64Array(3),
    IBLTE: new Int32Array(3),
    ISYS: new Array(nsys + 2),
    VA: createTensor3(3, 2, nsys),
    VB: createTensor3(3, 2, nsys),
    VM: createTensor3(3, nsys, nsys),
    VDEL: createTensor3(3, 2, nsys),
    VZ: createMatrix1(3, 2),
  };

  ctx.S[1] = s1;
  ctx.S[2] = s2;
  ctx.IBLTE[1] = iblte1;
  ctx.IBLTE[2] = iblte2;
  for (let i = 0; i < ctx.ISYS.length; i += 1) {
    ctx.ISYS[i] = new Int32Array(3);
  }
  if (iblte1 >= 1 && iblte1 <= nsys) {
    ctx.ISYS[iblte1][1] = ivte1;
  }
  if (iblte2 + 1 >= 1 && iblte2 + 1 <= nsys) {
    ctx.ISYS[iblte2 + 1][2] = ivz;
  }

  let idx = 0;
  for (let k = 1; k <= 3; k += 1) {
    for (let j = 1; j <= 2; j += 1) {
      for (let i = 1; i <= nsys; i += 1) {
        ctx.VA[k][j][i] = va[idx];
        idx += 1;
      }
    }
  }

  idx = 0;
  for (let k = 1; k <= 3; k += 1) {
    for (let j = 1; j <= 2; j += 1) {
      for (let i = 1; i <= nsys; i += 1) {
        ctx.VB[k][j][i] = vb[idx];
        idx += 1;
      }
    }
  }

  idx = 0;
  for (let k = 1; k <= 3; k += 1) {
    for (let l = 1; l <= nsys; l += 1) {
      for (let i = 1; i <= nsys; i += 1) {
        ctx.VM[k][l][i] = vm[idx];
        idx += 1;
      }
    }
  }

  idx = 0;
  for (let k = 1; k <= 3; k += 1) {
    for (let j = 1; j <= 2; j += 1) {
      for (let i = 1; i <= nsys; i += 1) {
        ctx.VDEL[k][j][i] = vdel[idx];
        idx += 1;
      }
    }
  }

  idx = 0;
  for (let k = 1; k <= 3; k += 1) {
    for (let j = 1; j <= 2; j += 1) {
      ctx.VZ[k][j] = vz[idx];
      idx += 1;
    }
  }

  if (!blsolv(ctx)) {
    throw new Error('blsolv failed');
  }

  const flattenTensor3 = (arr, d1, d2, d3) => {
    const data = new Array(d1 * d2 * d3);
    let out = 0;
    for (let k = 1; k <= d1; k += 1) {
      for (let j = 1; j <= d2; j += 1) {
        for (let i = 1; i <= d3; i += 1) {
          data[out] = arr[k][j][i] ?? 0.0;
          out += 1;
        }
      }
    }
    return data;
  };

  results.blsolv = {
    nsys,
    va: flattenTensor3(ctx.VA, 3, 2, nsys),
    vm: flattenTensor3(ctx.VM, 3, nsys, nsys),
    vdel: flattenTensor3(ctx.VDEL, 3, 2, nsys),
  };
}

process.stdout.write(JSON.stringify({ results }));
