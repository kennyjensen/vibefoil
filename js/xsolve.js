// Port of xsolve.f LU decomposition and back-substitution.
// Numerical linear algebra used by the panel/BL solvers (classic pivoting).

const NVX = 500;

// In-place Gaussian elimination with partial pivoting for multiple RHS (GAUSS).
// Solves general NxN system in NN unknowns with NRHS right-hand sides.
// Z is the coefficient matrix (destroyed during solution); R is overwritten
// by the solution vector(s).
function gauss(n, z, r, nrhs = 1) {
  const rhsIsVector = !Array.isArray(r[0]);

  const getR = (i, l) => (rhsIsVector ? r[i] : r[i][l]);
  const setR = (i, l, value) => {
    if (rhsIsVector) {
      r[i] = value;
    } else {
      r[i][l] = value;
    }
  };

  for (let np = 0; np < n - 1; np += 1) {
    const np1 = np + 1;

    // Find max pivot index NX.
    let nx = np;
    for (let nn = np1; nn < n; nn += 1) {
      if (Math.abs(z[nn][np]) > Math.abs(z[nx][np])) {
        nx = nn;
      }
    }

    const pivot = 1.0 / z[nx][np];

    // Switch pivots (copy NP diagonal into NX row without swapping the NP diagonal).
    z[nx][np] = z[np][np];

    // Switch rows & normalize pivot row.
    for (let l = np1; l < n; l += 1) {
      const temp = z[nx][l] * pivot;
      z[nx][l] = z[np][l];
      z[np][l] = temp;
    }

    // Switch RHS rows & normalize.
    for (let l = 0; l < nrhs; l += 1) {
      const temp = getR(nx, l) * pivot;
      setR(nx, l, getR(np, l));
      setR(np, l, temp);
    }

    // Forward eliminate everything.
    for (let k = np1; k < n; k += 1) {
      const ztmp = z[k][np];
      for (let l = np1; l < n; l += 1) {
        z[k][l] = z[k][l] - ztmp * z[np][l];
      }
      for (let l = 0; l < nrhs; l += 1) {
        setR(k, l, getR(k, l) - ztmp * getR(np, l));
      }
    }
  }

  // Solve for last row.
  for (let l = 0; l < nrhs; l += 1) {
    setR(n - 1, l, getR(n - 1, l) / z[n - 1][n - 1]);
  }

  // Back substitute everything.
  for (let np = n - 2; np >= 0; np -= 1) {
    const np1 = np + 1;
    for (let l = 0; l < nrhs; l += 1) {
      let value = getR(np, l);
      for (let k = np1; k < n; k += 1) {
        value -= z[np][k] * getR(k, l);
      }
      setR(np, l, value);
    }
  }
}

// LU decomposition with implicit scaling (LUDCMP).
function ludcmp(n, a, indx) {
  if (n > NVX) {
    throw new Error('LUDCMP: Array overflow. Increase NVX.');
  }

  const vv = new Float64Array(n);

  for (let i = 0; i < n; i += 1) {
    let aamax = 0.0;
    for (let j = 0; j < n; j += 1) {
      aamax = Math.max(Math.abs(a[i][j]), aamax);
    }
    vv[i] = 1.0 / aamax;
  }

  for (let j = 0; j < n; j += 1) {
    for (let i = 0; i < j; i += 1) {
      let sum = a[i][j];
      for (let k = 0; k < i; k += 1) {
        sum -= a[i][k] * a[k][j];
      }
      a[i][j] = sum;
    }

    let aamax = 0.0;
    let imax = j;
    for (let i = j; i < n; i += 1) {
      let sum = a[i][j];
      for (let k = 0; k < j; k += 1) {
        sum -= a[i][k] * a[k][j];
      }
      a[i][j] = sum;

      const dum = vv[i] * Math.abs(sum);
      if (dum >= aamax) {
        imax = i;
        aamax = dum;
      }
    }

    if (j !== imax) {
      const rowTmp = a[imax];
      a[imax] = a[j];
      a[j] = rowTmp;
      vv[imax] = vv[j];
    }

    indx[j] = imax;

    if (j !== n - 1) {
      const dum = 1.0 / a[j][j];
      for (let i = j + 1; i < n; i += 1) {
        a[i][j] *= dum;
      }
    }
  }
}

// Back substitution for LU factorization (BAKSUB).
function baksub(n, a, indx, b) {
  let ii = -1;

  for (let i = 0; i < n; i += 1) {
    const ll = indx[i];
    let sum = b[ll];
    b[ll] = b[i];
    if (ii !== -1) {
      for (let j = ii; j <= i - 1; j += 1) {
        sum -= a[i][j] * b[j];
      }
    } else if (sum !== 0.0) {
      ii = i;
    }
    b[i] = sum;
  }

  for (let i = n - 1; i >= 0; i -= 1) {
    let sum = b[i];
    if (i < n - 1) {
      for (let j = i + 1; j < n; j += 1) {
        sum -= a[i][j] * b[j];
      }
    }
    b[i] = sum / a[i][i];
  }
}

// Custom solver for coupled viscous-inviscid Newton system (XFOIL BLSOLV).
// Solves the block-tridiagonal BL system:
//   A d = R - dRe * S
// where A/B/Z are 3x3 blocks and d is the Newton delta for (Ctau, Theta, m).
function blsolv(ctx) {
  const nsys = ctx.NSYS;
  const ivte1 = ctx.ISYS[ctx.IBLTE[1]][1];
  const VA = ctx.VA;
  const VB = ctx.VB;
  const VM = ctx.VM;
  const VDEL = ctx.VDEL;
  const vacc1 = ctx.VACCEL;
  const vacc2 = ctx.VACCEL * 2.0 / (ctx.S[ctx.N] - ctx.S[1]);
  const vacc3 = ctx.VACCEL * 2.0 / (ctx.S[ctx.N] - ctx.S[1]);

  for (let iv = 1; iv <= nsys; iv += 1) {
    const ivp = iv + 1;

    const va11 = VA[1][1][iv];
    const va22 = VA[2][2][iv];
    const vm33 = VM[3][iv][iv];
    if (!Number.isFinite(va11) || !Number.isFinite(va22) || !Number.isFinite(vm33)
      || va11 === 0.0 || va22 === 0.0 || vm33 === 0.0) {
      console.warn('BLSOLV: singular block', { iv, va11, va22, vm33 });
      return false;
    }

    let pivot = 1.0 / va11;
    VA[1][2][iv] *= pivot;
    for (let l = iv; l <= nsys; l += 1) {
      VM[1][l][iv] *= pivot;
    }
    VDEL[1][1][iv] *= pivot;
    VDEL[1][2][iv] *= pivot;

    for (let k = 2; k <= 3; k += 1) {
      const vtmp = VA[k][1][iv];
      VA[k][2][iv] -= vtmp * VA[1][2][iv];
      for (let l = iv; l <= nsys; l += 1) {
        VM[k][l][iv] -= vtmp * VM[1][l][iv];
      }
      VDEL[k][1][iv] -= vtmp * VDEL[1][1][iv];
      VDEL[k][2][iv] -= vtmp * VDEL[1][2][iv];
    }

    pivot = 1.0 / VA[2][2][iv];
    for (let l = iv; l <= nsys; l += 1) {
      VM[2][l][iv] *= pivot;
    }
    VDEL[2][1][iv] *= pivot;
    VDEL[2][2][iv] *= pivot;

    {
      const k = 3;
      const vtmp = VA[k][2][iv];
      for (let l = iv; l <= nsys; l += 1) {
        VM[k][l][iv] -= vtmp * VM[2][l][iv];
      }
      VDEL[k][1][iv] -= vtmp * VDEL[2][1][iv];
      VDEL[k][2][iv] -= vtmp * VDEL[2][2][iv];
    }

    pivot = 1.0 / VM[3][iv][iv];
    for (let l = ivp; l <= nsys; l += 1) {
      VM[3][l][iv] *= pivot;
    }
    VDEL[3][1][iv] *= pivot;
    VDEL[3][2][iv] *= pivot;

    {
      const vtmp1 = VM[1][iv][iv];
      const vtmp2 = VM[2][iv][iv];
      for (let l = ivp; l <= nsys; l += 1) {
        VM[1][l][iv] -= vtmp1 * VM[3][l][iv];
        VM[2][l][iv] -= vtmp2 * VM[3][l][iv];
      }
      VDEL[1][1][iv] -= vtmp1 * VDEL[3][1][iv];
      VDEL[2][1][iv] -= vtmp2 * VDEL[3][1][iv];
      VDEL[1][2][iv] -= vtmp1 * VDEL[3][2][iv];
      VDEL[2][2][iv] -= vtmp2 * VDEL[3][2][iv];
    }

    {
      const vtmp = VA[1][2][iv];
      for (let l = ivp; l <= nsys; l += 1) {
        VM[1][l][iv] -= vtmp * VM[2][l][iv];
      }
      VDEL[1][1][iv] -= vtmp * VDEL[2][1][iv];
      VDEL[1][2][iv] -= vtmp * VDEL[2][2][iv];
    }

    if (iv === nsys) continue;

    for (let k = 1; k <= 3; k += 1) {
      const vtmp1 = VB[k][1][ivp];
      const vtmp2 = VB[k][2][ivp];
      const vtmp3 = VM[k][iv][ivp];
      for (let l = ivp; l <= nsys; l += 1) {
        VM[k][l][ivp] -= vtmp1 * VM[1][l][iv]
          + vtmp2 * VM[2][l][iv]
          + vtmp3 * VM[3][l][iv];
      }
      VDEL[k][1][ivp] -= vtmp1 * VDEL[1][1][iv]
        + vtmp2 * VDEL[2][1][iv]
        + vtmp3 * VDEL[3][1][iv];
      VDEL[k][2][ivp] -= vtmp1 * VDEL[1][2][iv]
        + vtmp2 * VDEL[2][2][iv]
        + vtmp3 * VDEL[3][2][iv];
    }

    if (iv === ivte1) {
      const ivz = ctx.ISYS[ctx.IBLTE[2] + 1][2];
      for (let k = 1; k <= 3; k += 1) {
        const vtmp1 = ctx.VZ[k][1];
        const vtmp2 = ctx.VZ[k][2];
        for (let l = ivp; l <= nsys; l += 1) {
          VM[k][l][ivz] -= vtmp1 * VM[1][l][iv]
            + vtmp2 * VM[2][l][iv];
        }
        VDEL[k][1][ivz] -= vtmp1 * VDEL[1][1][iv]
          + vtmp2 * VDEL[2][1][iv];
        VDEL[k][2][ivz] -= vtmp1 * VDEL[1][2][iv]
          + vtmp2 * VDEL[2][2][iv];
      }
    }

    if (ivp === nsys) continue;

    for (let kv = iv + 2; kv <= nsys; kv += 1) {
      const vtmp1 = VM[1][iv][kv];
      const vtmp2 = VM[2][iv][kv];
      const vtmp3 = VM[3][iv][kv];

      if (Math.abs(vtmp1) > vacc1) {
        for (let l = ivp; l <= nsys; l += 1) {
          VM[1][l][kv] -= vtmp1 * VM[3][l][iv];
        }
        VDEL[1][1][kv] -= vtmp1 * VDEL[3][1][iv];
        VDEL[1][2][kv] -= vtmp1 * VDEL[3][2][iv];
      }
      if (Math.abs(vtmp2) > vacc2) {
        for (let l = ivp; l <= nsys; l += 1) {
          VM[2][l][kv] -= vtmp2 * VM[3][l][iv];
        }
        VDEL[2][1][kv] -= vtmp2 * VDEL[3][1][iv];
        VDEL[2][2][kv] -= vtmp2 * VDEL[3][2][iv];
      }
      if (Math.abs(vtmp3) > vacc3) {
        for (let l = ivp; l <= nsys; l += 1) {
          VM[3][l][kv] -= vtmp3 * VM[3][l][iv];
        }
        VDEL[3][1][kv] -= vtmp3 * VDEL[3][1][iv];
        VDEL[3][2][kv] -= vtmp3 * VDEL[3][2][iv];
      }
    }
  }

  for (let iv = nsys; iv >= 2; iv -= 1) {
    let vtmp = VDEL[3][1][iv];
    for (let kv = iv - 1; kv >= 1; kv -= 1) {
      VDEL[1][1][kv] -= VM[1][iv][kv] * vtmp;
      VDEL[2][1][kv] -= VM[2][iv][kv] * vtmp;
      VDEL[3][1][kv] -= VM[3][iv][kv] * vtmp;
    }
    vtmp = VDEL[3][2][iv];
    for (let kv = iv - 1; kv >= 1; kv -= 1) {
      VDEL[1][2][kv] -= VM[1][iv][kv] * vtmp;
      VDEL[2][2][kv] -= VM[2][iv][kv] * vtmp;
      VDEL[3][2][kv] -= VM[3][iv][kv] * vtmp;
    }
  }

  return true;
}

export { gauss, ludcmp, baksub, blsolv };
