// Port of xsolve.f LU decomposition and back-substitution.
// Numerical linear algebra used by the panel/BL solvers (classic pivoting).

const NVX = 500;

// In-place Gaussian elimination with partial pivoting for multiple RHS.
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

    let nx = np;
    for (let nn = np1; nn < n; nn += 1) {
      if (Math.abs(z[nn][np]) > Math.abs(z[nx][np])) {
        nx = nn;
      }
    }

    const pivot = 1.0 / z[nx][np];

    const tempPivot = z[nx][np];
    z[nx][np] = z[np][np];
    z[np][np] = tempPivot;

    for (let l = np1; l < n; l += 1) {
      const temp = z[nx][l] * pivot;
      z[nx][l] = z[np][l];
      z[np][l] = temp;
    }

    for (let l = 0; l < nrhs; l += 1) {
      const temp = getR(nx, l) * pivot;
      setR(nx, l, getR(np, l));
      setR(np, l, temp);
    }

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

  for (let l = 0; l < nrhs; l += 1) {
    setR(n - 1, l, getR(n - 1, l) / z[n - 1][n - 1]);
  }

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

// LU decomposition with implicit scaling (Numerical Recipes style).
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

// Back substitution for LU factorization (permute then triangular solve).
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

export { gauss, ludcmp, baksub };
