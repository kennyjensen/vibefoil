import { createMatrix1, createTensor3Flat } from '../js/arrays.js';

let seed = 0x24122301;
function random() {
  seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
  return seed / 0x100000000;
}

export function makeContext(nsys, vaccel) {
  const ctx = {
    NSYS: nsys, N: 2, VACCEL: vaccel,
    S: Float64Array.of(0, 0.2, 2.7),
    IBLTE: Int32Array.of(0, 1, 2),
    ISYS: createMatrix1(3, 2), VZ: createMatrix1(3, 2),
  };
  const trailingEdge = Math.max(1, Math.floor(nsys / 3));
  ctx.ISYS[1][1] = trailingEdge;
  ctx.ISYS[3][2] = Math.min(nsys, trailingEdge + 2);
  for (const name of ['VA', 'VB', 'VM', 'VDEL']) {
    const { flat, view } = createTensor3Flat(3, name === 'VM' ? nsys : 2, nsys);
    // Nonzero padding also detects unintended writes outside active entries.
    for (let i = 0; i < flat.length; i += 1) flat[i] = 0.1 * (random() - 0.5);
    ctx[`${name}F`] = flat;
    ctx[name] = view;
  }
  for (let iv = 1; iv <= nsys; iv += 1) {
    ctx.VA[1][1][iv] += 3;
    ctx.VA[2][2][iv] += 3;
    ctx.VM[3][iv][iv] += 4;
  }
  for (let k = 1; k <= 3; k += 1) {
    ctx.VZ[k][1] = 0.1 * (random() - 0.5);
    ctx.VZ[k][2] = 0.1 * (random() - 0.5);
  }
  return ctx;
}

