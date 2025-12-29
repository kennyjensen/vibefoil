// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

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
  const VAF = ctx.VAF;
  const VBF = ctx.VBF;
  const VMF = ctx.VMF;
  const VDELF = ctx.VDELF;
  if (VAF && VBF && VMF && VDELF) {
    const ISYS = ctx.ISYS;
    const IBLTE = ctx.IBLTE;
    const ivte1 = ISYS[IBLTE[1]][1];
    const VZ = ctx.VZ;
    const s1 = ctx.S[1];
    const sN = ctx.S[ctx.N];
    const vacc1 = ctx.VACCEL;
    const vaccScale = ctx.VACCEL * 2.0 / (sN - s1);
    const vacc2 = vaccScale;
    const vacc3 = vaccScale;
    const vStride3 = nsys + 1;
    const vStride2 = 3 * vStride3;
    const vmStride3 = vStride3;
    const vmStride2 = (nsys + 1) * vmStride3;
    const vaBase1 = vStride2;
    const vaBase2 = 2 * vStride2;
    const vaBase3 = 3 * vStride2;
    const vbBase1 = vStride2;
    const vbBase2 = 2 * vStride2;
    const vbBase3 = 3 * vStride2;
    const vdelBase1 = vStride2;
    const vdelBase2 = 2 * vStride2;
    const vdelBase3 = 3 * vStride2;
    const vmBase1 = vmStride2;
    const vmBase2 = 2 * vmStride2;
    const vmBase3 = 3 * vmStride2;
    const va11Base = vaBase1 + vStride3;
    const va12Base = vaBase1 + 2 * vStride3;
    const va21Base = vaBase2 + vStride3;
    const va22Base = vaBase2 + 2 * vStride3;
    const va31Base = vaBase3 + vStride3;
    const va32Base = vaBase3 + 2 * vStride3;
    const vb11Base = vbBase1 + vStride3;
    const vb12Base = vbBase1 + 2 * vStride3;
    const vb21Base = vbBase2 + vStride3;
    const vb22Base = vbBase2 + 2 * vStride3;
    const vb31Base = vbBase3 + vStride3;
    const vb32Base = vbBase3 + 2 * vStride3;
    const vdel11Base = vdelBase1 + vStride3;
    const vdel12Base = vdelBase1 + 2 * vStride3;
    const vdel21Base = vdelBase2 + vStride3;
    const vdel22Base = vdelBase2 + 2 * vStride3;
    const vdel31Base = vdelBase3 + vStride3;
    const vdel32Base = vdelBase3 + 2 * vStride3;

    for (let iv = 1; iv <= nsys; iv += 1) {
      const ivp = iv + 1;

      const va11 = VAF[va11Base + iv];
      const va22 = VAF[va22Base + iv];
      const vm33 = VMF[vmBase3 + iv * vmStride3 + iv];
      if (!Number.isFinite(va11) || !Number.isFinite(va22) || !Number.isFinite(vm33)
        || va11 === 0.0 || va22 === 0.0 || vm33 === 0.0) {
        console.warn('BLSOLV: singular block', { iv, va11, va22, vm33 });
        return false;
      }

      let pivot = 1.0 / va11;
      VAF[va12Base + iv] *= pivot;
      for (let l = iv; l <= nsys; l += 1) {
        VMF[vmBase1 + l * vmStride3 + iv] *= pivot;
      }
      VDELF[vdel11Base + iv] *= pivot;
      VDELF[vdel12Base + iv] *= pivot;

      for (let k = 2; k <= 3; k += 1) {
        const vaBase = k === 2 ? vaBase2 : vaBase3;
        const vdelBase = k === 2 ? vdelBase2 : vdelBase3;
        const vmBase = k === 2 ? vmBase2 : vmBase3;
        const vtmp = VAF[vaBase + vStride3 + iv];
        VAF[vaBase + 2 * vStride3 + iv] -= vtmp * VAF[va12Base + iv];
        for (let l = iv; l <= nsys; l += 1) {
          const vmkl = vmBase + l * vmStride3 + iv;
          const vm1l = vmBase1 + l * vmStride3 + iv;
          VMF[vmkl] -= vtmp * VMF[vm1l];
        }
        VDELF[vdelBase + vStride3 + iv] -= vtmp * VDELF[vdel11Base + iv];
        VDELF[vdelBase + 2 * vStride3 + iv] -= vtmp * VDELF[vdel12Base + iv];
      }

      pivot = 1.0 / VAF[va22Base + iv];
      for (let l = iv; l <= nsys; l += 1) {
        VMF[vmBase2 + l * vmStride3 + iv] *= pivot;
      }
      VDELF[vdel21Base + iv] *= pivot;
      VDELF[vdel22Base + iv] *= pivot;

      {
        const vtmp = VAF[va32Base + iv];
        for (let l = iv; l <= nsys; l += 1) {
          const vm3l = vmBase3 + l * vmStride3 + iv;
          const vm2l = vmBase2 + l * vmStride3 + iv;
          VMF[vm3l] -= vtmp * VMF[vm2l];
        }
        VDELF[vdel31Base + iv] -= vtmp * VDELF[vdel21Base + iv];
        VDELF[vdel32Base + iv] -= vtmp * VDELF[vdel22Base + iv];
      }

      pivot = 1.0 / VMF[vmBase3 + iv * vmStride3 + iv];
      for (let l = ivp; l <= nsys; l += 1) {
        VMF[vmBase3 + l * vmStride3 + iv] *= pivot;
      }
      VDELF[vdel31Base + iv] *= pivot;
      VDELF[vdel32Base + iv] *= pivot;

      {
        const vtmp1 = VMF[vmBase1 + iv * vmStride3 + iv];
        const vtmp2 = VMF[vmBase2 + iv * vmStride3 + iv];
        for (let l = ivp; l <= nsys; l += 1) {
          const vm1l = vmBase1 + l * vmStride3 + iv;
          const vm2l = vmBase2 + l * vmStride3 + iv;
          const vm3l = vmBase3 + l * vmStride3 + iv;
          VMF[vm1l] -= vtmp1 * VMF[vm3l];
          VMF[vm2l] -= vtmp2 * VMF[vm3l];
        }
        VDELF[vdel11Base + iv] -= vtmp1 * VDELF[vdel31Base + iv];
        VDELF[vdel21Base + iv] -= vtmp2 * VDELF[vdel31Base + iv];
        VDELF[vdel12Base + iv] -= vtmp1 * VDELF[vdel32Base + iv];
        VDELF[vdel22Base + iv] -= vtmp2 * VDELF[vdel32Base + iv];
      }

      {
        const vtmp = VAF[va12Base + iv];
        for (let l = ivp; l <= nsys; l += 1) {
          const vm1l = vmBase1 + l * vmStride3 + iv;
          const vm2l = vmBase2 + l * vmStride3 + iv;
          VMF[vm1l] -= vtmp * VMF[vm2l];
        }
        VDELF[vdel11Base + iv] -= vtmp * VDELF[vdel21Base + iv];
        VDELF[vdel12Base + iv] -= vtmp * VDELF[vdel22Base + iv];
      }

      if (iv === nsys) continue;

      for (let k = 1; k <= 3; k += 1) {
        const vbBase = k === 1 ? vbBase1 : (k === 2 ? vbBase2 : vbBase3);
        const vdelBase = k === 1 ? vdelBase1 : (k === 2 ? vdelBase2 : vdelBase3);
        const vmBase = k === 1 ? vmBase1 : (k === 2 ? vmBase2 : vmBase3);
        const vtmp1 = VBF[vbBase + vStride3 + ivp];
        const vtmp2 = VBF[vbBase + 2 * vStride3 + ivp];
        const vtmp3 = VMF[vmBase + iv * vmStride3 + ivp];
        for (let l = ivp; l <= nsys; l += 1) {
          const vmkl = vmBase + l * vmStride3 + ivp;
          const vm1l = vmBase1 + l * vmStride3 + iv;
          const vm2l = vmBase2 + l * vmStride3 + iv;
          const vm3l = vmBase3 + l * vmStride3 + iv;
          VMF[vmkl] -= vtmp1 * VMF[vm1l]
            + vtmp2 * VMF[vm2l]
            + vtmp3 * VMF[vm3l];
        }
        VDELF[vdelBase + vStride3 + ivp] -= vtmp1 * VDELF[vdel11Base + iv]
          + vtmp2 * VDELF[vdel21Base + iv]
          + vtmp3 * VDELF[vdel31Base + iv];
        VDELF[vdelBase + 2 * vStride3 + ivp] -= vtmp1 * VDELF[vdel12Base + iv]
          + vtmp2 * VDELF[vdel22Base + iv]
          + vtmp3 * VDELF[vdel32Base + iv];
      }

      if (iv === ivte1) {
        const ivz = ISYS[IBLTE[2] + 1][2];
        for (let k = 1; k <= 3; k += 1) {
          const vtmp1 = VZ[k][1];
          const vtmp2 = VZ[k][2];
          const vdelBase = k === 1 ? vdelBase1 : (k === 2 ? vdelBase2 : vdelBase3);
          const vmBase = k === 1 ? vmBase1 : (k === 2 ? vmBase2 : vmBase3);
          for (let l = ivp; l <= nsys; l += 1) {
            const vmkl = vmBase + l * vmStride3 + ivz;
            const vm1l = vmBase1 + l * vmStride3 + iv;
            const vm2l = vmBase2 + l * vmStride3 + iv;
            VMF[vmkl] -= vtmp1 * VMF[vm1l]
              + vtmp2 * VMF[vm2l];
          }
          VDELF[vdelBase + vStride3 + ivz] -= vtmp1 * VDELF[vdel11Base + iv]
            + vtmp2 * VDELF[vdel21Base + iv];
          VDELF[vdelBase + 2 * vStride3 + ivz] -= vtmp1 * VDELF[vdel12Base + iv]
            + vtmp2 * VDELF[vdel22Base + iv];
        }
      }

      if (ivp === nsys) continue;

      for (let kv = iv + 2; kv <= nsys; kv += 1) {
        const vtmp1 = VMF[vmBase1 + iv * vmStride3 + kv];
        const vtmp2 = VMF[vmBase2 + iv * vmStride3 + kv];
        const vtmp3 = VMF[vmBase3 + iv * vmStride3 + kv];

        if (Math.abs(vtmp1) > vacc1) {
          for (let l = ivp; l <= nsys; l += 1) {
            const vm1l = vmBase1 + l * vmStride3 + kv;
            const vm3l = vmBase3 + l * vmStride3 + iv;
            VMF[vm1l] -= vtmp1 * VMF[vm3l];
          }
          VDELF[vdel11Base + kv] -= vtmp1 * VDELF[vdel31Base + iv];
          VDELF[vdel12Base + kv] -= vtmp1 * VDELF[vdel32Base + iv];
        }
        if (Math.abs(vtmp2) > vacc2) {
          for (let l = ivp; l <= nsys; l += 1) {
            const vm2l = vmBase2 + l * vmStride3 + kv;
            const vm3l = vmBase3 + l * vmStride3 + iv;
            VMF[vm2l] -= vtmp2 * VMF[vm3l];
          }
          VDELF[vdel21Base + kv] -= vtmp2 * VDELF[vdel31Base + iv];
          VDELF[vdel22Base + kv] -= vtmp2 * VDELF[vdel32Base + iv];
        }
        if (Math.abs(vtmp3) > vacc3) {
          for (let l = ivp; l <= nsys; l += 1) {
            const vm3l = vmBase3 + l * vmStride3 + kv;
            const vm3lIv = vmBase3 + l * vmStride3 + iv;
            VMF[vm3l] -= vtmp3 * VMF[vm3lIv];
          }
          VDELF[vdel31Base + kv] -= vtmp3 * VDELF[vdel31Base + iv];
          VDELF[vdel32Base + kv] -= vtmp3 * VDELF[vdel32Base + iv];
        }
      }
    }

    for (let iv = nsys; iv >= 2; iv -= 1) {
      let vtmp = VDELF[vdel31Base + iv];
      for (let kv = iv - 1; kv >= 1; kv -= 1) {
        VDELF[vdel11Base + kv] -= VMF[vmBase1 + iv * vmStride3 + kv] * vtmp;
        VDELF[vdel21Base + kv] -= VMF[vmBase2 + iv * vmStride3 + kv] * vtmp;
        VDELF[vdel31Base + kv] -= VMF[vmBase3 + iv * vmStride3 + kv] * vtmp;
      }
      vtmp = VDELF[vdel32Base + iv];
      for (let kv = iv - 1; kv >= 1; kv -= 1) {
        VDELF[vdel12Base + kv] -= VMF[vmBase1 + iv * vmStride3 + kv] * vtmp;
        VDELF[vdel22Base + kv] -= VMF[vmBase2 + iv * vmStride3 + kv] * vtmp;
        VDELF[vdel32Base + kv] -= VMF[vmBase3 + iv * vmStride3 + kv] * vtmp;
      }
    }

    return true;
  }
  const ISYS = ctx.ISYS;
  const IBLTE = ctx.IBLTE;
  const ivte1 = ISYS[IBLTE[1]][1];
  const VA = ctx.VA;
  const VB = ctx.VB;
  const VM = ctx.VM;
  const VDEL = ctx.VDEL;
  const VZ = ctx.VZ;
  const VA1 = VA[1];
  const VA2 = VA[2];
  const VA3 = VA[3];
  const VM1 = VM[1];
  const VM2 = VM[2];
  const VM3 = VM[3];
  const VDEL1 = VDEL[1];
  const VDEL2 = VDEL[2];
  const VDEL3 = VDEL[3];
  const s1 = ctx.S[1];
  const sN = ctx.S[ctx.N];
  const vacc1 = ctx.VACCEL;
  const vaccScale = ctx.VACCEL * 2.0 / (sN - s1);
  const vacc2 = vaccScale;
  const vacc3 = vaccScale;

  for (let iv = 1; iv <= nsys; iv += 1) {
    const ivp = iv + 1;

    const va11 = VA1[1][iv];
    const va22 = VA2[2][iv];
    const vm33 = VM3[iv][iv];
    if (!Number.isFinite(va11) || !Number.isFinite(va22) || !Number.isFinite(vm33)
      || va11 === 0.0 || va22 === 0.0 || vm33 === 0.0) {
      console.warn('BLSOLV: singular block', { iv, va11, va22, vm33 });
      return false;
    }

    let pivot = 1.0 / va11;
    VA1[2][iv] *= pivot;
    for (let l = iv; l <= nsys; l += 1) {
      const vm1l = VM1[l];
      vm1l[iv] *= pivot;
    }
    VDEL1[1][iv] *= pivot;
    VDEL1[2][iv] *= pivot;

    for (let k = 2; k <= 3; k += 1) {
      const vtmp = VA[k][1][iv];
      VA[k][2][iv] -= vtmp * VA1[2][iv];
      const VMk = VM[k];
      const VDELk = VDEL[k];
      for (let l = iv; l <= nsys; l += 1) {
        const vmkl = VMk[l];
        const vm1l = VM1[l];
        vmkl[iv] -= vtmp * vm1l[iv];
      }
      VDELk[1][iv] -= vtmp * VDEL1[1][iv];
      VDELk[2][iv] -= vtmp * VDEL1[2][iv];
    }

    pivot = 1.0 / VA2[2][iv];
    for (let l = iv; l <= nsys; l += 1) {
      const vm2l = VM2[l];
      vm2l[iv] *= pivot;
    }
    VDEL2[1][iv] *= pivot;
    VDEL2[2][iv] *= pivot;

    {
      const k = 3;
      const vtmp = VA3[2][iv];
      for (let l = iv; l <= nsys; l += 1) {
        const vm3l = VM3[l];
        const vm2l = VM2[l];
        vm3l[iv] -= vtmp * vm2l[iv];
      }
      VDEL3[1][iv] -= vtmp * VDEL2[1][iv];
      VDEL3[2][iv] -= vtmp * VDEL2[2][iv];
    }

    pivot = 1.0 / VM3[iv][iv];
    for (let l = ivp; l <= nsys; l += 1) {
      const vm3l = VM3[l];
      vm3l[iv] *= pivot;
    }
    VDEL3[1][iv] *= pivot;
    VDEL3[2][iv] *= pivot;

    {
      const vtmp1 = VM1[iv][iv];
      const vtmp2 = VM2[iv][iv];
      for (let l = ivp; l <= nsys; l += 1) {
        const vm1l = VM1[l];
        const vm2l = VM2[l];
        const vm3l = VM3[l];
        vm1l[iv] -= vtmp1 * vm3l[iv];
        vm2l[iv] -= vtmp2 * vm3l[iv];
      }
      VDEL1[1][iv] -= vtmp1 * VDEL3[1][iv];
      VDEL2[1][iv] -= vtmp2 * VDEL3[1][iv];
      VDEL1[2][iv] -= vtmp1 * VDEL3[2][iv];
      VDEL2[2][iv] -= vtmp2 * VDEL3[2][iv];
    }

    {
      const vtmp = VA1[2][iv];
      for (let l = ivp; l <= nsys; l += 1) {
        const vm1l = VM1[l];
        const vm2l = VM2[l];
        vm1l[iv] -= vtmp * vm2l[iv];
      }
      VDEL1[1][iv] -= vtmp * VDEL2[1][iv];
      VDEL1[2][iv] -= vtmp * VDEL2[2][iv];
    }

    if (iv === nsys) continue;

    for (let k = 1; k <= 3; k += 1) {
      const vtmp1 = VB[k][1][ivp];
      const vtmp2 = VB[k][2][ivp];
      const vtmp3 = VM[k][iv][ivp];
      const VMk = VM[k];
      const VDELk = VDEL[k];
      for (let l = ivp; l <= nsys; l += 1) {
        const vmkl = VMk[l];
        const vm1l = VM1[l];
        const vm2l = VM2[l];
        const vm3l = VM3[l];
        vmkl[ivp] -= vtmp1 * vm1l[iv]
          + vtmp2 * vm2l[iv]
          + vtmp3 * vm3l[iv];
      }
      VDELk[1][ivp] -= vtmp1 * VDEL1[1][iv]
        + vtmp2 * VDEL2[1][iv]
        + vtmp3 * VDEL3[1][iv];
      VDELk[2][ivp] -= vtmp1 * VDEL1[2][iv]
        + vtmp2 * VDEL2[2][iv]
        + vtmp3 * VDEL3[2][iv];
    }

    if (iv === ivte1) {
      const ivz = ISYS[IBLTE[2] + 1][2];
      for (let k = 1; k <= 3; k += 1) {
        const vtmp1 = VZ[k][1];
        const vtmp2 = VZ[k][2];
        const VMk = VM[k];
        const VDELk = VDEL[k];
        for (let l = ivp; l <= nsys; l += 1) {
          const vmkl = VMk[l];
          const vm1l = VM1[l];
          const vm2l = VM2[l];
          vmkl[ivz] -= vtmp1 * vm1l[iv]
            + vtmp2 * vm2l[iv];
        }
        VDELk[1][ivz] -= vtmp1 * VDEL1[1][iv]
          + vtmp2 * VDEL2[1][iv];
        VDELk[2][ivz] -= vtmp1 * VDEL1[2][iv]
          + vtmp2 * VDEL2[2][iv];
      }
    }

    if (ivp === nsys) continue;

    for (let kv = iv + 2; kv <= nsys; kv += 1) {
      const vtmp1 = VM1[iv][kv];
      const vtmp2 = VM2[iv][kv];
      const vtmp3 = VM3[iv][kv];

      if (Math.abs(vtmp1) > vacc1) {
        for (let l = ivp; l <= nsys; l += 1) {
          const vm1l = VM1[l];
          const vm3l = VM3[l];
          vm1l[kv] -= vtmp1 * vm3l[iv];
        }
        VDEL1[1][kv] -= vtmp1 * VDEL3[1][iv];
        VDEL1[2][kv] -= vtmp1 * VDEL3[2][iv];
      }
      if (Math.abs(vtmp2) > vacc2) {
        for (let l = ivp; l <= nsys; l += 1) {
          const vm2l = VM2[l];
          const vm3l = VM3[l];
          vm2l[kv] -= vtmp2 * vm3l[iv];
        }
        VDEL2[1][kv] -= vtmp2 * VDEL3[1][iv];
        VDEL2[2][kv] -= vtmp2 * VDEL3[2][iv];
      }
      if (Math.abs(vtmp3) > vacc3) {
        for (let l = ivp; l <= nsys; l += 1) {
          const vm3l = VM3[l];
          vm3l[kv] -= vtmp3 * vm3l[iv];
        }
        VDEL3[1][kv] -= vtmp3 * VDEL3[1][iv];
        VDEL3[2][kv] -= vtmp3 * VDEL3[2][iv];
      }
    }
  }

  for (let iv = nsys; iv >= 2; iv -= 1) {
    let vtmp = VDEL3[1][iv];
    const vm1iv = VM1[iv];
    const vm2iv = VM2[iv];
    const vm3iv = VM3[iv];
    for (let kv = iv - 1; kv >= 1; kv -= 1) {
      VDEL1[1][kv] -= vm1iv[kv] * vtmp;
      VDEL2[1][kv] -= vm2iv[kv] * vtmp;
      VDEL3[1][kv] -= vm3iv[kv] * vtmp;
    }
    vtmp = VDEL3[2][iv];
    for (let kv = iv - 1; kv >= 1; kv -= 1) {
      VDEL1[2][kv] -= vm1iv[kv] * vtmp;
      VDEL2[2][kv] -= vm2iv[kv] * vtmp;
      VDEL3[2][kv] -= vm3iv[kv] * vtmp;
    }
  }

  return true;
}

export { gauss, ludcmp, baksub, blsolv };
