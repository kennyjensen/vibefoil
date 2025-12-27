// Port of xbl.f (partial). Marching and utility routines for the integral BL.
// Integral boundary-layer method: marching, closures, and coupling terms.

import { splind, sinvrt, seval } from './spline.js';
import { gauss } from './xsolve.js';

// 1-based to 0-based wrapper for Gaussian elimination (Fortran compatibility).
function gauss1(n, z, r) {
  const z0 = new Array(n);
  for (let i = 0; i < n; i += 1) {
    z0[i] = new Float64Array(n);
    for (let j = 0; j < n; j += 1) {
      z0[i][j] = z[i + 1][j + 1];
    }
  }

  const r0 = new Float64Array(n);
  for (let i = 0; i < n; i += 1) {
    r0[i] = r[i + 1];
  }

  gauss(n, z0, r0);

  for (let i = 0; i < n; i += 1) {
    r[i + 1] = r0[i];
    for (let j = 0; j < n; j += 1) {
      z[i + 1][j + 1] = z0[i][j];
    }
  }
}

// Displacement thickness limiter based on shape factor constraint (DSLIM).
function dslim(ctx, dstr, thet, uedg, msq, hklim) {
  const h = dstr / thet;
  const { hk, hkH } = ctx.hkin(h, msq);
  const dh = Math.max(0.0, hklim - hk) / hkH;
  return dstr + dh * thet;
}

// Initialize BL empirical constants (BLPINI).
function blpini(ctx) {
  ctx.SCCON = 5.6;
  ctx.GACON = 6.70;
  ctx.GBCON = 0.75;
  ctx.GCCON = 18.0;
  ctx.DLCON = 0.9;

  ctx.CTRCON = 1.8;
  ctx.CTRCEX = 3.3;

  ctx.DUXCON = 1.0;

  ctx.CTCON = 0.5 / (ctx.GACON ** 2 * ctx.GBCON);

  ctx.CFFAC = 1.0;
}

// Mach/Re updates for fixed CL or alpha families (MRCL).
function mrcl(ctx, cls) {
  const cla = Math.max(cls, 1.0e-6);
  let minf;
  let mCls = 0.0;
  let reinf;
  let rCls = 0.0;

  const matyp = ctx.MATYP ?? 1;
  const retyp = ctx.RETYP ?? 1;

  if (matyp === 1) {
    minf = ctx.MINF1;
    mCls = 0.0;
  } else if (matyp === 2) {
    minf = ctx.MINF1 / Math.sqrt(cla);
    mCls = -0.5 * minf / cla;
  } else {
    minf = ctx.MINF1;
    mCls = 0.0;
  }

  if (retyp === 1) {
    reinf = ctx.REINF1;
    rCls = 0.0;
  } else if (retyp === 2) {
    reinf = ctx.REINF1 / Math.sqrt(cla);
    rCls = -0.5 * reinf / cla;
  } else {
    reinf = ctx.REINF1 / cla;
    rCls = -reinf / cla;
  }

  ctx.MINF = minf;
  ctx.REINF = reinf;
  return { mCls, rCls };
}

// Precompute coefficients and intermediate BL parameters for marching (COMSET).
function comset(ctx) {
  const beta = Math.sqrt(Math.max(1.0 - ctx.MINF ** 2, 0.0));
  const betaMsq = beta === 0.0 ? 0.0 : -0.5 / beta;
  ctx.TKLAM = ctx.MINF ** 2 / (1.0 + beta) ** 2;
  ctx.TK_MSQ = 1.0 / (1.0 + beta) ** 2 - 2.0 * ctx.TKLAM / (1.0 + beta) * betaMsq;
}

// Initialize BL system matrices and coupling terms for marching (IBLSYS).
function iblsys(ctx) {
  let iv = 0;
  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
      iv += 1;
      ctx.ISYS[ibl][is] = iv;
    }
  }

  ctx.NSYS = iv;
  if (ctx.NSYS > 2 * ctx.IVX) {
    throw new Error('*** IBLSYS: BL system array overflow. ***');
  }
}

// March Ue and related quantities along the surface and wake (MRCHUE).
// Laminar/turbulent switching and wake handling follow XFOIL's marching scheme.
function mrchue(ctx) {
  const hlmax = 3.8;
  const htmax = 2.5;
  if (!ctx.HTARG || ctx.HTARG.length < 3
    || ctx.HTARG[1].length <= ctx.NBL[1]
    || ctx.HTARG[2].length <= ctx.NBL[2]) {
    ctx.HTARG = [
      null,
      new Float64Array(ctx.NBL[1] + 1),
      new Float64Array(ctx.NBL[2] + 1),
    ];
  } else {
    ctx.HTARG[1].fill(0.0);
    ctx.HTARG[2].fill(0.0);
  }

  function resetThwaites(xsi, uei) {
    const ucon = uei / (xsi ** ctx.BULE);
    const tsq = 0.45 / (ucon * (5.0 * ctx.BULE + 1.0) * ctx.REYBL)
      * (xsi ** (1.0 - ctx.BULE));
    const thi = Math.sqrt(tsq);
    const dsi = 2.2 * thi;
    return { thi, dsi };
  }

  for (let is = 1; is <= 2; is += 1) {
    ctx.AMCRIT = ctx.ACRIT[is];

    xifset(ctx, is);

    let ibl = 2;
    let xsi = ctx.XSSI[ibl][is];
    let uei = ctx.UEDG[ibl][is];
    ctx.BULE = 1.0;
    const init = resetThwaites(xsi, uei);
    let thi = init.thi;
    let dsi = init.dsi;
    let ami = 0.0;

    let cti = 0.03;

    // Initialize laminar start on each side.
    ctx.TRAN = false;
    ctx.TURB = false;
    ctx.ITRAN[is] = ctx.IBLTE[is];
    for (ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
      const ibm = ibl - 1;

      const simi = ibl === 2;
      const wake = ibl > ctx.IBLTE[is];
      ctx.TRAN = ibl === ctx.ITRAN[is];
      ctx.TURB = ibl > ctx.ITRAN[is];
      ctx.SIMI = simi;
      ctx.WAKE = wake;

      xsi = ctx.XSSI[ibl][is];
      uei = ctx.UEDG[ibl][is];

      let dswaki = 0.0;
      if (wake) {
        const iw = ibl - ctx.IBLTE[is];
        dswaki = ctx.WGAP[iw];
      }

      let direct = true;
      let htarg = 0.0;
      let dmax = 0.0;
      let msq = 0.0;
      let hklim = 0.0;

      for (let itbl = 1; itbl <= 25; itbl += 1) {
        ctx.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
        ctx.blkin(ctx);

        if (!simi && !ctx.TURB) {
          ctx.trchek(ctx);
          ami = ctx.AMPL2;

          if (ctx.TRAN) {
            ctx.ITRAN[is] = ibl;
            if (cti <= 0.0) {
              cti = 0.03;
              ctx.S2 = cti;
            }
          } else {
            ctx.ITRAN[is] = ibl + 2;
          }
        }

      if (ibl === ctx.IBLTE[is] + 1) {
        const tte = ctx.THET[ctx.IBLTE[1]][1] + ctx.THET[ctx.IBLTE[2]][2];
        const dte = ctx.DSTR[ctx.IBLTE[1]][1] + ctx.DSTR[ctx.IBLTE[2]][2] + ctx.ANTE;
        const cte = (ctx.CTAU[ctx.IBLTE[1]][1] * ctx.THET[ctx.IBLTE[1]][1]
          + ctx.CTAU[ctx.IBLTE[2]][2] * ctx.THET[ctx.IBLTE[2]][2]) / tte;
        ctx.CTE = cte;
        ctx.TTE = tte;
        ctx.DTE = dte;
        ctx.tesys(cte, tte, dte, ctx);
      } else {
        ctx.blsys(ctx);
      }


        if (direct) {
          ctx.VS2[4][1] = 0.0;
          ctx.VS2[4][2] = 0.0;
          ctx.VS2[4][3] = 0.0;
          ctx.VS2[4][4] = 1.0;
          ctx.VSREZ[4] = 0.0;

          gauss1(4, ctx.VS2, ctx.VSREZ);

          dmax = Math.max(Math.abs(ctx.VSREZ[2] / thi), Math.abs(ctx.VSREZ[3] / dsi));
          if (ibl < ctx.ITRAN[is]) dmax = Math.max(dmax, Math.abs(ctx.VSREZ[1] / 10.0));
          if (ibl >= ctx.ITRAN[is]) dmax = Math.max(dmax, Math.abs(ctx.VSREZ[1] / cti));

          let rlx = 1.0;
          if (dmax > 0.3) rlx = 0.3 / dmax;

          let hktest = null;
          let hmax = null;
          if (ibl !== ctx.IBLTE[is] + 1) {
            msq = uei * uei * ctx.HSTINV / (ctx.GM1BL * (1.0 - 0.5 * uei * uei * ctx.HSTINV));
            const htest = (dsi + rlx * ctx.VSREZ[3]) / (thi + rlx * ctx.VSREZ[2]);
            const hkin = ctx.hkin(htest, msq);
            hktest = hkin.hk;

            hmax = htmax;
            if (ibl < ctx.ITRAN[is]) hmax = hlmax;
            direct = hktest < hmax;
          }

        if (direct) {
          if (ibl >= ctx.ITRAN[is]) cti = cti + rlx * ctx.VSREZ[1];
          thi = thi + rlx * ctx.VSREZ[2];
          dsi = dsi + rlx * ctx.VSREZ[3];
        } else {
            if (ibl < ctx.ITRAN[is]) {
              htarg = ctx.HK1 + 0.03 * (ctx.X2 - ctx.X1) / ctx.T1;
            } else if (ibl === ctx.ITRAN[is]) {
              htarg = ctx.HK1 + (0.03 * (ctx.XT - ctx.X1) - 0.15 * (ctx.X2 - ctx.XT)) / ctx.T1;
            } else if (wake) {
              const constant = 0.03 * (ctx.X2 - ctx.X1) / ctx.T1;
              let hk2 = ctx.HK1;
              hk2 = hk2 - (hk2 + constant * (hk2 - 1.0) ** 3 - ctx.HK1)
                / (1.0 + 3.0 * constant * (hk2 - 1.0) ** 2);
              hk2 = hk2 - (hk2 + constant * (hk2 - 1.0) ** 3 - ctx.HK1)
                / (1.0 + 3.0 * constant * (hk2 - 1.0) ** 2);
              hk2 = hk2 - (hk2 + constant * (hk2 - 1.0) ** 3 - ctx.HK1)
                / (1.0 + 3.0 * constant * (hk2 - 1.0) ** 2);
              htarg = hk2;
          } else {
            htarg = ctx.HK1 - 0.15 * (ctx.X2 - ctx.X1) / ctx.T1;
          }

          if (wake) {
            htarg = Math.max(htarg, 1.01);
          } else {
            const hmax = ibl < ctx.ITRAN[is] ? hlmax : htmax;
            htarg = Math.max(htarg, hmax);
          }

          ctx.HTARG[is][ibl] = htarg;
          console.log(` MRCHUE: Inverse mode at ${String(ibl).padStart(4)}     Hk =${htarg.toFixed(3).padStart(8)}`);

          // Inverse mode: retry with prescribed Hk.
          direct = false;
          continue;
          }
        } else {
          ctx.VS2[4][1] = 0.0;
          ctx.VS2[4][2] = ctx.HK2_T2;
          ctx.VS2[4][3] = ctx.HK2_D2;
          ctx.VS2[4][4] = ctx.HK2_U2;
          ctx.VSREZ[4] = htarg - ctx.HK2;

          gauss1(4, ctx.VS2, ctx.VSREZ);

          dmax = Math.max(Math.abs(ctx.VSREZ[2] / thi),
            Math.abs(ctx.VSREZ[3] / dsi),
            Math.abs(ctx.VSREZ[4] / uei));
          if (ibl >= ctx.ITRAN[is]) dmax = Math.max(dmax, Math.abs(ctx.VSREZ[1] / cti));

          let rlx = 1.0;
          if (dmax > 0.3) rlx = 0.3 / dmax;

          if (ibl >= ctx.ITRAN[is]) cti = cti + rlx * ctx.VSREZ[1];
          thi = thi + rlx * ctx.VSREZ[2];
          dsi = dsi + rlx * ctx.VSREZ[3];
          uei = uei + rlx * ctx.VSREZ[4];
        }

        if (ibl >= ctx.ITRAN[is]) {
          cti = Math.min(cti, 0.30);
          cti = Math.max(cti, 0.0000001);
        }

        if (ibl <= ctx.IBLTE[is]) {
          hklim = 1.02;
        } else {
          hklim = 1.00005;
        }
        msq = uei * uei * ctx.HSTINV / (ctx.GM1BL * (1.0 - 0.5 * uei * uei * ctx.HSTINV));
        let dsw = dsi - dswaki;
        dsw = dslim(ctx, dsw, thi, uei, msq, hklim);
        dsi = dsw + dswaki;

        if (dmax <= 1.0e-5) {
          break;
        }
      }

      if (dmax > 1.0e-5) {
        console.log(` MRCHUE: Convergence failed at${String(ibl).padStart(4)}  side${String(is).padStart(3)}    Res =${dmax.toExponential(4).padStart(12)}`);
        if (dmax > 0.1 && ibl > 3) {
          if (ibl <= ctx.IBLTE[is]) {
            thi = ctx.THET[ibm][is] * (ctx.XSSI[ibl][is] / ctx.XSSI[ibm][is]) ** 0.5;
            dsi = ctx.DSTR[ibm][is] * (ctx.XSSI[ibl][is] / ctx.XSSI[ibm][is]) ** 0.5;
          } else if (ibl === ctx.IBLTE[is] + 1) {
            cti = ctx.CTE;
            thi = ctx.TTE;
            dsi = ctx.DTE;
          } else {
            thi = ctx.THET[ibm][is];
            const ratlen = (ctx.XSSI[ibl][is] - ctx.XSSI[ibm][is]) / (10.0 * ctx.DSTR[ibm][is]);
            dsi = (ctx.DSTR[ibm][is] + thi * ratlen) / (1.0 + ratlen);
          }
          if (ibl === ctx.ITRAN[is]) cti = 0.05;
          if (ibl > ctx.ITRAN[is]) cti = ctx.CTAU[ibm][is];

          uei = ctx.UEDG[ibl][is];
          if (ibl > 2 && ibl < ctx.NBL[is]) {
            uei = 0.5 * (ctx.UEDG[ibl - 1][is] + ctx.UEDG[ibl + 1][is]);
          }
        }

        ctx.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
        ctx.blkin(ctx);
        if (!simi && !ctx.TURB) {
          ctx.trchek(ctx);
          ami = ctx.AMPL2;
          if (ctx.TRAN) ctx.ITRAN[is] = ibl;
          if (!ctx.TRAN) ctx.ITRAN[is] = ibl + 2;
        }
        if (ibl < ctx.ITRAN[is]) ctx.blvar(1, ctx);
        if (ibl >= ctx.ITRAN[is]) ctx.blvar(2, ctx);
        if (wake) ctx.blvar(3, ctx);

        if (ibl < ctx.ITRAN[is]) ctx.blmid(1, ctx);
        if (ibl >= ctx.ITRAN[is]) ctx.blmid(2, ctx);
        if (wake) ctx.blmid(3, ctx);
      }

      if (ibl < ctx.ITRAN[is]) ctx.CTAU[ibl][is] = ami;
      if (ibl >= ctx.ITRAN[is]) ctx.CTAU[ibl][is] = cti;
      ctx.THET[ibl][is] = thi;
      ctx.DSTR[ibl][is] = dsi;
      ctx.UEDG[ibl][is] = uei;
      ctx.MASS[ibl][is] = dsi * uei;
      ctx.TAU[ibl][is] = 0.5 * ctx.R2 * ctx.U2 * ctx.U2 * ctx.CF2;
      ctx.DIS[ibl][is] = ctx.R2 * ctx.U2 * ctx.U2 * ctx.U2 * ctx.DI2 * ctx.HS2 * 0.5;
      ctx.CTQ[ibl][is] = ctx.CQ2;
      ctx.DELT[ibl][is] = ctx.DE2;
      ctx.TSTR[ibl][is] = ctx.HS2 * ctx.T2;


      ctx.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
      ctx.blkin(ctx);
      for (let icom = 1; icom <= ctx.NCOM; icom += 1) {
        ctx.COM1[icom] = ctx.COM2[icom];
      }
      if (typeof ctx.syncComToVars === 'function') {
        ctx.syncComToVars(ctx, 1);
      }

    // Set transition state and force turbulence at TE.
    if (ctx.TRAN || ibl === ctx.IBLTE[is]) {
      ctx.TURB = true;
      ctx.TFORCE[is] = ctx.TRFORC;
      ctx.XSSITR[is] = ctx.XT;
    }

      ctx.TRAN = false;

      if (ibl === ctx.IBLTE[is]) {
        thi = ctx.THET[ctx.IBLTE[1]][1] + ctx.THET[ctx.IBLTE[2]][2];
        dsi = ctx.DSTR[ctx.IBLTE[1]][1] + ctx.DSTR[ctx.IBLTE[2]][2] + ctx.ANTE;
      }
    }
  }
}

// March Ue and BL variables with transition logic (MRCHDU).
function mrchdu(ctx) {
  // Fortran comments (MRCHDU) highlight:
  // - March with transition logic and sensitivity weighting.
  // - Update CTAU/THET/DSTR/UEDG based on branch (lam/turb/wake).
  const deps = 5.0e-6;
  const senswt = 1000.0;

  for (let is = 1; is <= 2; is += 1) {
    ctx.AMCRIT = ctx.ACRIT[is];
    xifset(ctx, is);

    let ibl = 2;
    let xsi = ctx.XSSI[ibl][is];
    let uei = ctx.UEDG[ibl][is];
    ctx.BULE = 1.0;

    // Save previous transition index.
    const itrold = ctx.ITRAN[is];

    ctx.TRAN = false;
    ctx.TURB = false;
    ctx.ITRAN[is] = ctx.IBLTE[is];

    let sens = 0.0;
    let ami = 0.0;

    for (ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
      const ibm = ibl - 1;

      const simi = ibl === 2;
      const wake = ibl > ctx.IBLTE[is];
      ctx.SIMI = simi;
      ctx.WAKE = wake;
      ctx.TRAN = ibl === ctx.ITRAN[is];
      ctx.TURB = ibl > ctx.ITRAN[is];

      xsi = ctx.XSSI[ibl][is];
      uei = ctx.UEDG[ibl][is];
      let thi = ctx.THET[ibl][is];
      let dsi = ctx.DSTR[ibl][is];

      let cti = 0.03;
      if (ibl < itrold) {
        ami = ctx.CTAU[ibl][is];
      } else {
        cti = ctx.CTAU[ibl][is];
        if (cti <= 0.0) cti = 0.03;
      }

      let dswaki = 0.0;
      if (wake) {
        const iw = ibl - ctx.IBLTE[is];
        dswaki = ctx.WGAP[iw];
      }

      if (ibl <= ctx.IBLTE[is]) {
        dsi = Math.max(dsi - dswaki, 1.02000 * thi) + dswaki;
      }
      if (ibl > ctx.IBLTE[is]) {
        dsi = Math.max(dsi - dswaki, 1.00005 * thi) + dswaki;
      }

      let ueref = 0.0;
      let hkref = 0.0;
      let sennew = sens;
      let dmax = 0.0;

      for (let itbl = 1; itbl <= 25; itbl += 1) {
        ctx.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
        ctx.blkin(ctx);

        if (!simi && !ctx.TURB) {
          ctx.trchek(ctx);
          ami = ctx.AMPL2;
          if (ctx.TRAN) ctx.ITRAN[is] = ibl;
          if (!ctx.TRAN) ctx.ITRAN[is] = ibl + 2;
        }

        if (ibl === ctx.IBLTE[is] + 1) {
          const tte = ctx.THET[ctx.IBLTE[1]][1] + ctx.THET[ctx.IBLTE[2]][2];
          const dte = ctx.DSTR[ctx.IBLTE[1]][1] + ctx.DSTR[ctx.IBLTE[2]][2] + ctx.ANTE;
          const cte = (ctx.CTAU[ctx.IBLTE[1]][1] * ctx.THET[ctx.IBLTE[1]][1]
            + ctx.CTAU[ctx.IBLTE[2]][2] * ctx.THET[ctx.IBLTE[2]][2]) / tte;
          ctx.CTE = cte;
          ctx.TTE = tte;
          ctx.DTE = dte;
          ctx.tesys(cte, tte, dte, ctx);
        } else {
          ctx.blsys(ctx);
        }

        if (itbl === 1) {
          ueref = ctx.U2;
          hkref = ctx.HK2;

          if (ibl < ctx.ITRAN[is] && ibl >= itrold) {
            const uem = ctx.UEDG[ibl - 1][is];
            const dsm = ctx.DSTR[ibl - 1][is];
            const thm = ctx.THET[ibl - 1][is];
            const msq = uem * uem * ctx.HSTINV / (ctx.GM1BL * (1.0 - 0.5 * uem * uem * ctx.HSTINV));
            hkref = ctx.hkin(dsm / thm, msq).hk;
          }

          if (ibl < itrold) {
            if (ctx.TRAN) ctx.CTAU[ibl][is] = 0.03;
            if (ctx.TURB) ctx.CTAU[ibl][is] = ctx.CTAU[ibl - 1][is];
        // Turbulent or transition branch.
        if (ctx.TRAN || ctx.TURB) {
              cti = ctx.CTAU[ibl][is];
              ctx.S2 = cti;
            }
          }
        }

        if (simi || ibl === ctx.IBLTE[is] + 1) {
          ctx.VS2[4][1] = 0.0;
          ctx.VS2[4][2] = 0.0;
          ctx.VS2[4][3] = 0.0;
          ctx.VS2[4][4] = ctx.U2_UEI;
          ctx.VSREZ[4] = ueref - ctx.U2;
        } else {
          const vtmp = Array.from({ length: 5 }, () => new Float64Array(5));
          const vztmp = new Float64Array(5);

          for (let k = 1; k <= 4; k += 1) {
            vztmp[k] = ctx.VSREZ[k];
            for (let l = 1; l <= 4; l += 1) {
              vtmp[k][l] = ctx.VS2[k][l];
            }
          }

          vtmp[4][1] = 0.0;
          vtmp[4][2] = ctx.HK2_T2;
          vtmp[4][3] = ctx.HK2_D2;
          vtmp[4][4] = ctx.HK2_U2 * ctx.U2_UEI;
          vztmp[4] = 1.0;

          gauss1(4, vtmp, vztmp);

          sennew = senswt * vztmp[4] * hkref / ueref;
          if (itbl <= 5) {
            sens = sennew;
          } else if (itbl <= 15) {
            sens = 0.5 * (sens + sennew);
          }

          ctx.VS2[4][1] = 0.0;
          ctx.VS2[4][2] = ctx.HK2_T2 * hkref;
          ctx.VS2[4][3] = ctx.HK2_D2 * hkref;
          ctx.VS2[4][4] = (ctx.HK2_U2 * hkref + sens / ueref) * ctx.U2_UEI;
          ctx.VSREZ[4] = -(hkref ** 2) * (ctx.HK2 / hkref - 1.0)
            - sens * (ctx.U2 / ueref - 1.0);
        }

        gauss1(4, ctx.VS2, ctx.VSREZ);

        dmax = Math.max(Math.abs(ctx.VSREZ[2] / thi),
          Math.abs(ctx.VSREZ[3] / dsi),
          Math.abs(ctx.VSREZ[4] / uei));
        if (ibl >= ctx.ITRAN[is]) dmax = Math.max(dmax, Math.abs(ctx.VSREZ[1] / (10.0 * cti)));

        let rlx = 1.0;
        if (dmax > 0.3) rlx = 0.3 / dmax;

        if (ibl < ctx.ITRAN[is]) ami = ami + rlx * ctx.VSREZ[1];
        if (ibl >= ctx.ITRAN[is]) cti = cti + rlx * ctx.VSREZ[1];
        thi = thi + rlx * ctx.VSREZ[2];
        dsi = dsi + rlx * ctx.VSREZ[3];
        uei = uei + rlx * ctx.VSREZ[4];

        if (ibl >= ctx.ITRAN[is]) {
          cti = Math.min(cti, 0.30);
          cti = Math.max(cti, 0.0000001);
        }

        let hklim = 0.0;
        if (ibl <= ctx.IBLTE[is]) {
          hklim = 1.02;
        } else {
          hklim = 1.00005;
        }
        const msq = uei * uei * ctx.HSTINV / (ctx.GM1BL * (1.0 - 0.5 * uei * uei * ctx.HSTINV));
        let dsw = dsi - dswaki;
        dsw = dslim(ctx, dsw, thi, uei, msq, hklim);
        dsi = dsw + dswaki;

        if (dmax <= deps) {
          break;
        }
      }

      if (dmax > deps && dmax > 0.1 && ibl > 3) {
        if (ibl <= ctx.IBLTE[is]) {
          thi = ctx.THET[ibm][is] * (ctx.XSSI[ibl][is] / ctx.XSSI[ibm][is]) ** 0.5;
          dsi = ctx.DSTR[ibm][is] * (ctx.XSSI[ibl][is] / ctx.XSSI[ibm][is]) ** 0.5;
          uei = ctx.UEDG[ibm][is];
        } else if (ibl === ctx.IBLTE[is] + 1) {
          cti = ctx.CTE;
          thi = ctx.TTE;
          dsi = ctx.DTE;
          uei = ctx.UEDG[ibm][is];
        } else {
          thi = ctx.THET[ibm][is];
          const ratlen = (ctx.XSSI[ibl][is] - ctx.XSSI[ibm][is]) / (10.0 * ctx.DSTR[ibm][is]);
          dsi = (ctx.DSTR[ibm][is] + thi * ratlen) / (1.0 + ratlen);
          uei = ctx.UEDG[ibm][is];
        }
        if (ibl === ctx.ITRAN[is]) cti = 0.05;
        if (ibl > ctx.ITRAN[is]) cti = ctx.CTAU[ibm][is];
      }

      if (dmax > deps) {
        ctx.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
        ctx.blkin(ctx);

        if (!simi && !ctx.TURB) {
          ctx.trchek(ctx);
          ami = ctx.AMPL2;
          if (ctx.TRAN) ctx.ITRAN[is] = ibl;
          if (!ctx.TRAN) ctx.ITRAN[is] = ibl + 2;
        }

        if (ibl < ctx.ITRAN[is]) ctx.blvar(1, ctx);
        if (ibl >= ctx.ITRAN[is]) ctx.blvar(2, ctx);
        if (wake) ctx.blvar(3, ctx);

        if (ibl < ctx.ITRAN[is]) ctx.blmid(1, ctx);
        if (ibl >= ctx.ITRAN[is]) ctx.blmid(2, ctx);
        if (wake) ctx.blmid(3, ctx);
      }

      if (ibl < ctx.ITRAN[is]) ctx.CTAU[ibl][is] = ami;
      if (ibl >= ctx.ITRAN[is]) ctx.CTAU[ibl][is] = cti;
      ctx.THET[ibl][is] = thi;
      ctx.DSTR[ibl][is] = dsi;
      ctx.UEDG[ibl][is] = uei;
      ctx.MASS[ibl][is] = dsi * uei;
      ctx.TAU[ibl][is] = 0.5 * ctx.R2 * ctx.U2 * ctx.U2 * ctx.CF2;
      ctx.DIS[ibl][is] = ctx.R2 * ctx.U2 * ctx.U2 * ctx.U2 * ctx.DI2 * ctx.HS2 * 0.5;
      ctx.CTQ[ibl][is] = ctx.CQ2;
      ctx.DELT[ibl][is] = ctx.DE2;
      ctx.TSTR[ibl][is] = ctx.HS2 * ctx.T2;

      ctx.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
      ctx.blkin(ctx);
      for (let icom = 1; icom <= ctx.NCOM; icom += 1) {
        ctx.COM1[icom] = ctx.COM2[icom];
      }
      if (typeof ctx.syncComToVars === 'function') {
        ctx.syncComToVars(ctx, 1);
      }

      // Force transition at TE or if transition detected.
      if (ctx.TRAN || ibl === ctx.IBLTE[is]) {
        ctx.TURB = true;
        ctx.TFORCE[is] = ctx.TRFORC;
        ctx.XSSITR[is] = ctx.XT;
      }

      ctx.TRAN = false;
    }
  }
}

// Assemble current BL state into system form prior to Newton iteration (SETBL).
function setbl(ctx) {
  // Fortran comments (SETBL) highlight:
  // - Assemble system matrices and residuals for Newton iteration.
  // - Set transition indices, compute TE coupling, and fill V arrays.
  const nsys = ctx.NSYS;
  const COM1 = ctx.COM1;
  const COM2 = ctx.COM2;
  const VA = ctx.VA;
  const VB = ctx.VB;
  const VDEL = ctx.VDEL;
  const VM = ctx.VM;
  // No tracing; setbl assembles VA/VB/VDEL/VM directly.
  // Clear COM arrays before assembling system.
  COM1.fill(0.0, 1);
  COM2.fill(0.0, 1);
  // Index maps for system line -> side/BL index.
  if (!ctx.IV_TO_IS || ctx.IV_TO_IS.length < nsys + 1) {
    ctx.IV_TO_IS = new Int32Array(nsys + 1);
    ctx.IV_TO_IBL = new Int32Array(nsys + 1);
  } else {
    ctx.IV_TO_IS.fill(0);
    ctx.IV_TO_IBL.fill(0);
  }
  const ivToIs = ctx.IV_TO_IS;
  const ivToIbl = ctx.IV_TO_IBL;
  const NBL = ctx.NBL;
  const ISYS = ctx.ISYS;
  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= NBL[is]; ibl += 1) {
      const iv = ISYS[ibl][is];
      ivToIs[iv] = is;
      ivToIbl[iv] = ibl;
    }
  }
  for (let k = 1; k <= 3; k += 1) {
    for (let j = 1; j <= 2; j += 1) {
      VA[k][j].fill(0.0, 1, nsys + 1);
      VB[k][j].fill(0.0, 1, nsys + 1);
      VDEL[k][j].fill(0.0, 1, nsys + 1);
    }
    for (let j = 1; j <= nsys; j += 1) {
      VM[k][j].fill(0.0, 1, nsys + 1);
    }
  }

  const clmr = ctx.LALFA ? ctx.CL : ctx.CLSPEC;
  const { mCls: maClmr, rCls: reClmr } = mrcl(ctx, clmr);
  const msqClmr = 2.0 * ctx.MINF * maClmr;

  comset(ctx);

  ctx.GAMBL = ctx.GAMMA;
  ctx.GM1BL = ctx.GAMM1;
  ctx.QINFBL = ctx.QINF;
  ctx.TKBL = ctx.TKLAM;
  ctx.TKBL_MS = ctx.TK_MSQ;

  ctx.RSTBL = (1.0 + 0.5 * ctx.GM1BL * ctx.MINF ** 2) ** (1.0 / ctx.GM1BL);
  ctx.RSTBL_MS = 0.5 * ctx.RSTBL / (1.0 + 0.5 * ctx.GM1BL * ctx.MINF ** 2);

  ctx.HSTINV = ctx.GM1BL * (ctx.MINF / ctx.QINFBL) ** 2 / (1.0 + 0.5 * ctx.GM1BL * ctx.MINF ** 2);
  ctx.HSTINV_MS = ctx.GM1BL * (1.0 / ctx.QINFBL) ** 2 / (1.0 + 0.5 * ctx.GM1BL * ctx.MINF ** 2)
    - 0.5 * ctx.GM1BL * ctx.HSTINV / (1.0 + 0.5 * ctx.GM1BL * ctx.MINF ** 2);

  const herat = 1.0 - 0.5 * ctx.QINFBL ** 2 * ctx.HSTINV;
  const heratMs = -0.5 * ctx.QINFBL ** 2 * ctx.HSTINV_MS;

  ctx.REYBL = ctx.REINF * Math.sqrt(herat ** 3) * (1.0 + ctx.HVRAT) / (herat + ctx.HVRAT);
  ctx.REYBL_RE = Math.sqrt(herat ** 3) * (1.0 + ctx.HVRAT) / (herat + ctx.HVRAT);
  ctx.REYBL_MS = ctx.REYBL * (1.5 / herat - 1.0 / (herat + ctx.HVRAT)) * heratMs;

  ctx.IDAMPV = ctx.IDAMP ?? ctx.IDAMPV ?? 0;
  ctx.DWTE = ctx.WGAP[1] ?? 0.0;

  if (!ctx.LBLINI) {
    mrchue(ctx);
    ctx.LBLINI = true;
  }

  mrchdu(ctx);

  const maxNbl = Math.max(ctx.NBL[1], ctx.NBL[2]);
  const usav = new Array(maxNbl + 1);
  for (let i = 0; i <= maxNbl; i += 1) {
    usav[i] = new Float64Array(3);
  }

  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
      usav[ibl][is] = ctx.UEDG[ibl][is];
    }
  }

  ueset(ctx);

  let badUeSet = 0;
  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
      if (!Number.isFinite(ctx.UEDG[ibl][is])) {
        badUeSet += 1;
      }
    }
  }

  if (badUeSet === 0) {
    for (let is = 1; is <= 2; is += 1) {
      for (let ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
        const temp = usav[ibl][is];
        usav[ibl][is] = ctx.UEDG[ibl][is];
        ctx.UEDG[ibl][is] = temp;
      }
    }
  } else {
    for (let is = 1; is <= 2; is += 1) {
      for (let ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
        ctx.UEDG[ibl][is] = usav[ibl][is];
      }
    }
  }

  const ile1 = ctx.IPAN[2][1];
  const ile2 = ctx.IPAN[2][2];
  const ite1 = ctx.IPAN[ctx.IBLTE[1]][1];
  const ite2 = ctx.IPAN[ctx.IBLTE[2]][2];
  const jvte1 = ctx.ISYS[ctx.IBLTE[1]][1];
  const jvte2 = ctx.ISYS[ctx.IBLTE[2]][2];
  let dule1 = ctx.UEDG[2][1] - usav[2][1];
  let dule2 = ctx.UEDG[2][2] - usav[2][2];

  const u1M = new Float64Array(nsys + 1);
  const u2M = new Float64Array(nsys + 1);
  const d1M = new Float64Array(nsys + 1);
  const d2M = new Float64Array(nsys + 1);
  const ule1M = new Float64Array(nsys + 1);
  const ule2M = new Float64Array(nsys + 1);
  const ute1M = new Float64Array(nsys + 1);
  const ute2M = new Float64Array(nsys + 1);

  for (let js = 1; js <= 2; js += 1) {
    for (let jbl = 2; jbl <= ctx.NBL[js]; jbl += 1) {
      const j = ctx.IPAN[jbl][js];
      const jv = ctx.ISYS[jbl][js];
      ule1M[jv] = -ctx.VTI[2][1] * ctx.VTI[jbl][js] * ctx.DIJ[ile1][j];
      ule2M[jv] = -ctx.VTI[2][2] * ctx.VTI[jbl][js] * ctx.DIJ[ile2][j];
      ute1M[jv] = -ctx.VTI[ctx.IBLTE[1]][1] * ctx.VTI[jbl][js] * ctx.DIJ[ite1][j];
      ute2M[jv] = -ctx.VTI[ctx.IBLTE[2]][2] * ctx.VTI[jbl][js] * ctx.DIJ[ite2][j];
    }
  }

  const ule1A = ctx.UINV_A[2][1];
  const ule2A = ctx.UINV_A[2][2];

  ctx.TINDEX[1] = 0.0;
  ctx.TINDEX[2] = 0.0;

  for (let is = 1; is <= 2; is += 1) {
    for (let js = 1; js <= 2; js += 1) {
      for (let jbl = 2; jbl <= ctx.NBL[js]; jbl += 1) {
        const jv = ctx.ISYS[jbl][js];
        u1M[jv] = 0.0;
        d1M[jv] = 0.0;
      }
    }
    let u1A = 0.0;
    let d1A = 0.0;
    let due1 = 0.0;
    let dds1 = 0.0;

    ctx.BULE = 1.0;
    ctx.AMCRIT = ctx.ACRIT[is];
    xifset(ctx, is);

    ctx.TRAN = false;
    ctx.TURB = false;

    let ami = 0.0;
    let cti = 0.0;
    for (let ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
      const iv = ctx.ISYS[ibl][is];

      ctx.SIMI = ibl === 2;
      ctx.WAKE = ibl > ctx.IBLTE[is];
      ctx.TRAN = ibl === ctx.ITRAN[is];
      ctx.TURB = ibl > ctx.ITRAN[is];

      const i = ctx.IPAN[ibl][is];

      const xsi = ctx.XSSI[ibl][is];
      if (ibl < ctx.ITRAN[is]) ami = ctx.CTAU[ibl][is];
      if (ibl >= ctx.ITRAN[is]) cti = ctx.CTAU[ibl][is];
      const uei = ctx.UEDG[ibl][is];
      const thi = ctx.THET[ibl][is];
      const mdi = ctx.MASS[ibl][is];
      const dsi = mdi / uei;

      let dswaki = 0.0;
      if (ctx.WAKE) {
        const iw = ibl - ctx.IBLTE[is];
        dswaki = ctx.WGAP[iw];
      }

      const d2M2 = 1.0 / uei;
      const d2U2 = -dsi / uei;

      for (let js = 1; js <= 2; js += 1) {
        for (let jbl = 2; jbl <= ctx.NBL[js]; jbl += 1) {
          const j = ctx.IPAN[jbl][js];
          const jv = ctx.ISYS[jbl][js];
          u2M[jv] = -ctx.VTI[ibl][is] * ctx.VTI[jbl][js] * ctx.DIJ[i][j];
          d2M[jv] = d2U2 * u2M[jv];
        }
      }
      d2M[iv] += d2M2;

      const u2A = ctx.UINV_A[ibl][is];
      const d2A = d2U2 * u2A;

      const due2 = ctx.UEDG[ibl][is] - usav[ibl][is];
      const dds2 = d2U2 * due2;

      ctx.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
      ctx.blkin(ctx);

      if (ctx.TRAN) {
        ctx.trchek(ctx);
        ami = ctx.AMPL2;
      }
      if (ibl === ctx.ITRAN[is] && !ctx.TRAN) {
        // no-op: log suppressed
      }

      const xiUle1 = is === 1 ? ctx.SST_GO : -ctx.SST_GO;
      const xiUle2 = is === 1 ? -ctx.SST_GP : ctx.SST_GP;

      if (ibl === ctx.IBLTE[is] + 1) {
        const tte = ctx.THET[ctx.IBLTE[1]][1] + ctx.THET[ctx.IBLTE[2]][2];
        const dte = ctx.DSTR[ctx.IBLTE[1]][1] + ctx.DSTR[ctx.IBLTE[2]][2] + ctx.ANTE;
        const cte = (ctx.CTAU[ctx.IBLTE[1]][1] * ctx.THET[ctx.IBLTE[1]][1]
          + ctx.CTAU[ctx.IBLTE[2]][2] * ctx.THET[ctx.IBLTE[2]][2]) / tte;
        ctx.tesys(cte, tte, dte, ctx);

        const tteTte1 = 1.0;
        const tteTte2 = 1.0;
        const dteMte1 = 1.0 / ctx.UEDG[ctx.IBLTE[1]][1];
        const dteUte1 = -ctx.DSTR[ctx.IBLTE[1]][1] / ctx.UEDG[ctx.IBLTE[1]][1];
        const dteMte2 = 1.0 / ctx.UEDG[ctx.IBLTE[2]][2];
        const dteUte2 = -ctx.DSTR[ctx.IBLTE[2]][2] / ctx.UEDG[ctx.IBLTE[2]][2];
        const cteCte1 = ctx.THET[ctx.IBLTE[1]][1] / tte;
        const cteCte2 = ctx.THET[ctx.IBLTE[2]][2] / tte;
        const cteTte1 = (ctx.CTAU[ctx.IBLTE[1]][1] - cte) / tte;
        const cteTte2 = (ctx.CTAU[ctx.IBLTE[2]][2] - cte) / tte;

        for (let js = 1; js <= 2; js += 1) {
          for (let jbl = 2; jbl <= ctx.NBL[js]; jbl += 1) {
            const j = ctx.IPAN[jbl][js];
            const jv = ctx.ISYS[jbl][js];
            d1M[jv] = dteUte1 * ute1M[jv] + dteUte2 * ute2M[jv];
          }
        }
        d1M[jvte1] += dteMte1;
        d1M[jvte2] += dteMte2;

        due1 = 0.0;
        dds1 = dteUte1 * (ctx.UEDG[ctx.IBLTE[1]][1] - usav[ctx.IBLTE[1]][1])
          + dteUte2 * (ctx.UEDG[ctx.IBLTE[2]][2] - usav[ctx.IBLTE[2]][2]);

        for (let jv = 1; jv <= nsys; jv += 1) {
          ctx.VM[1][jv][iv] = ctx.VS1[1][3] * d1M[jv] + ctx.VS1[1][4] * u1M[jv]
            + ctx.VS2[1][3] * d2M[jv] + ctx.VS2[1][4] * u2M[jv]
            + (ctx.VS1[1][5] + ctx.VS2[1][5] + ctx.VSX[1])
            * (xiUle1 * ule1M[jv] + xiUle2 * ule2M[jv]);
        }

        ctx.VB[1][1][iv] = ctx.VS1[1][1];
        ctx.VB[1][2][iv] = ctx.VS1[1][2];
        ctx.VA[1][1][iv] = ctx.VS2[1][1];
        ctx.VA[1][2][iv] = ctx.VS2[1][2];

        ctx.VDEL[1][2][iv] = ctx.LALFA
          ? ctx.VSR[1] * reClmr + ctx.VSM[1] * msqClmr
          : (ctx.VS1[1][4] * u1A + ctx.VS1[1][3] * d1A)
            + (ctx.VS2[1][4] * u2A + ctx.VS2[1][3] * d2A);
        ctx.VDEL[1][2][iv] += (ctx.VS1[1][5] + ctx.VS2[1][5] + ctx.VSX[1])
          * (xiUle1 * ule1A + xiUle2 * ule2A);
        ctx.VDEL[1][1][iv] = ctx.VSREZ[1]
          + (ctx.VS1[1][4] * due1 + ctx.VS1[1][3] * dds1)
          + (ctx.VS2[1][4] * due2 + ctx.VS2[1][3] * dds2)
          + (ctx.VS1[1][5] + ctx.VS2[1][5] + ctx.VSX[1])
          * (xiUle1 * dule1 + xiUle2 * dule2);

        for (let jv = 1; jv <= nsys; jv += 1) {
          ctx.VM[2][jv][iv] = ctx.VS1[2][3] * d1M[jv] + ctx.VS1[2][4] * u1M[jv]
            + ctx.VS2[2][3] * d2M[jv] + ctx.VS2[2][4] * u2M[jv]
            + (ctx.VS1[2][5] + ctx.VS2[2][5] + ctx.VSX[2])
            * (xiUle1 * ule1M[jv] + xiUle2 * ule2M[jv]);
        }
        ctx.VB[2][1][iv] = ctx.VS1[2][1];
        ctx.VB[2][2][iv] = ctx.VS1[2][2];
        ctx.VA[2][1][iv] = ctx.VS2[2][1];
        ctx.VA[2][2][iv] = ctx.VS2[2][2];
        ctx.VDEL[2][2][iv] = ctx.LALFA
          ? ctx.VSR[2] * reClmr + ctx.VSM[2] * msqClmr
          : (ctx.VS1[2][4] * u1A + ctx.VS1[2][3] * d1A)
            + (ctx.VS2[2][4] * u2A + ctx.VS2[2][3] * d2A);
        ctx.VDEL[2][2][iv] += (ctx.VS1[2][5] + ctx.VS2[2][5] + ctx.VSX[2])
          * (xiUle1 * ule1A + xiUle2 * ule2A);
        ctx.VDEL[2][1][iv] = ctx.VSREZ[2]
          + (ctx.VS1[2][4] * due1 + ctx.VS1[2][3] * dds1)
          + (ctx.VS2[2][4] * due2 + ctx.VS2[2][3] * dds2)
          + (ctx.VS1[2][5] + ctx.VS2[2][5] + ctx.VSX[2])
          * (xiUle1 * dule1 + xiUle2 * dule2);

        for (let jv = 1; jv <= nsys; jv += 1) {
          ctx.VM[3][jv][iv] = ctx.VS1[3][3] * d1M[jv] + ctx.VS1[3][4] * u1M[jv]
            + ctx.VS2[3][3] * d2M[jv] + ctx.VS2[3][4] * u2M[jv]
            + (ctx.VS1[3][5] + ctx.VS2[3][5] + ctx.VSX[3])
            * (xiUle1 * ule1M[jv] + xiUle2 * ule2M[jv]);
        }
        ctx.VB[3][1][iv] = ctx.VS1[3][1];
        ctx.VB[3][2][iv] = ctx.VS1[3][2];
        ctx.VA[3][1][iv] = ctx.VS2[3][1];
        ctx.VA[3][2][iv] = ctx.VS2[3][2];
        ctx.VDEL[3][2][iv] = ctx.LALFA
          ? ctx.VSR[3] * reClmr + ctx.VSM[3] * msqClmr
          : (ctx.VS1[3][4] * u1A + ctx.VS1[3][3] * d1A)
            + (ctx.VS2[3][4] * u2A + ctx.VS2[3][3] * d2A);
        ctx.VDEL[3][2][iv] += (ctx.VS1[3][5] + ctx.VS2[3][5] + ctx.VSX[3])
          * (xiUle1 * ule1A + xiUle2 * ule2A);
        ctx.VDEL[3][1][iv] = ctx.VSREZ[3]
          + (ctx.VS1[3][4] * due1 + ctx.VS1[3][3] * dds1)
          + (ctx.VS2[3][4] * due2 + ctx.VS2[3][3] * dds2)
          + (ctx.VS1[3][5] + ctx.VS2[3][5] + ctx.VSX[3])
          * (xiUle1 * dule1 + xiUle2 * dule2);

        ctx.VZ[1][1] = ctx.VS1[1][1] * cteCte1;
        ctx.VZ[1][2] = ctx.VS1[1][1] * cteTte1 + ctx.VS1[1][2] * tteTte1;
        ctx.VB[1][1][iv] = ctx.VS1[1][1] * cteCte2;
        ctx.VB[1][2][iv] = ctx.VS1[1][1] * cteTte2 + ctx.VS1[1][2] * tteTte2;

        ctx.VZ[2][1] = ctx.VS1[2][1] * cteCte1;
        ctx.VZ[2][2] = ctx.VS1[2][1] * cteTte1 + ctx.VS1[2][2] * tteTte1;
        ctx.VB[2][1][iv] = ctx.VS1[2][1] * cteCte2;
        ctx.VB[2][2][iv] = ctx.VS1[2][1] * cteTte2 + ctx.VS1[2][2] * tteTte2;

        ctx.VZ[3][1] = ctx.VS1[3][1] * cteCte1;
        ctx.VZ[3][2] = ctx.VS1[3][1] * cteTte1 + ctx.VS1[3][2] * tteTte1;
        ctx.VB[3][1][iv] = ctx.VS1[3][1] * cteCte2;
        ctx.VB[3][2][iv] = ctx.VS1[3][1] * cteTte2 + ctx.VS1[3][2] * tteTte2;
      } else {
        ctx.blsys(ctx);

        for (let jv = 1; jv <= nsys; jv += 1) {
          ctx.VM[1][jv][iv] = ctx.VS1[1][3] * d1M[jv] + ctx.VS1[1][4] * u1M[jv]
            + ctx.VS2[1][3] * d2M[jv] + ctx.VS2[1][4] * u2M[jv]
            + (ctx.VS1[1][5] + ctx.VS2[1][5] + ctx.VSX[1])
            * (xiUle1 * ule1M[jv] + xiUle2 * ule2M[jv]);
        }
        ctx.VB[1][1][iv] = ctx.VS1[1][1];
        ctx.VB[1][2][iv] = ctx.VS1[1][2];
        ctx.VA[1][1][iv] = ctx.VS2[1][1];
        ctx.VA[1][2][iv] = ctx.VS2[1][2];
        ctx.VDEL[1][2][iv] = ctx.LALFA
          ? ctx.VSR[1] * reClmr + ctx.VSM[1] * msqClmr
          : (ctx.VS1[1][4] * u1A + ctx.VS1[1][3] * d1A)
            + (ctx.VS2[1][4] * u2A + ctx.VS2[1][3] * d2A)
            + (ctx.VS1[1][5] + ctx.VS2[1][5] + ctx.VSX[1])
            * (xiUle1 * ule1A + xiUle2 * ule2A);
        ctx.VDEL[1][1][iv] = ctx.VSREZ[1]
          + (ctx.VS1[1][4] * due1 + ctx.VS1[1][3] * dds1)
          + (ctx.VS2[1][4] * due2 + ctx.VS2[1][3] * dds2)
          + (ctx.VS1[1][5] + ctx.VS2[1][5] + ctx.VSX[1])
          * (xiUle1 * dule1 + xiUle2 * dule2);

        for (let jv = 1; jv <= nsys; jv += 1) {
          ctx.VM[2][jv][iv] = ctx.VS1[2][3] * d1M[jv] + ctx.VS1[2][4] * u1M[jv]
            + ctx.VS2[2][3] * d2M[jv] + ctx.VS2[2][4] * u2M[jv]
            + (ctx.VS1[2][5] + ctx.VS2[2][5] + ctx.VSX[2])
            * (xiUle1 * ule1M[jv] + xiUle2 * ule2M[jv]);
        }
        ctx.VB[2][1][iv] = ctx.VS1[2][1];
        ctx.VB[2][2][iv] = ctx.VS1[2][2];
        ctx.VA[2][1][iv] = ctx.VS2[2][1];
        ctx.VA[2][2][iv] = ctx.VS2[2][2];
        ctx.VDEL[2][2][iv] = ctx.LALFA
          ? ctx.VSR[2] * reClmr + ctx.VSM[2] * msqClmr
          : (ctx.VS1[2][4] * u1A + ctx.VS1[2][3] * d1A)
            + (ctx.VS2[2][4] * u2A + ctx.VS2[2][3] * d2A)
            + (ctx.VS1[2][5] + ctx.VS2[2][5] + ctx.VSX[2])
            * (xiUle1 * ule1A + xiUle2 * ule2A);
        ctx.VDEL[2][1][iv] = ctx.VSREZ[2]
          + (ctx.VS1[2][4] * due1 + ctx.VS1[2][3] * dds1)
          + (ctx.VS2[2][4] * due2 + ctx.VS2[2][3] * dds2)
          + (ctx.VS1[2][5] + ctx.VS2[2][5] + ctx.VSX[2])
          * (xiUle1 * dule1 + xiUle2 * dule2);

        for (let jv = 1; jv <= nsys; jv += 1) {
          ctx.VM[3][jv][iv] = ctx.VS1[3][3] * d1M[jv] + ctx.VS1[3][4] * u1M[jv]
            + ctx.VS2[3][3] * d2M[jv] + ctx.VS2[3][4] * u2M[jv]
            + (ctx.VS1[3][5] + ctx.VS2[3][5] + ctx.VSX[3])
            * (xiUle1 * ule1M[jv] + xiUle2 * ule2M[jv]);
        }
        ctx.VB[3][1][iv] = ctx.VS1[3][1];
        ctx.VB[3][2][iv] = ctx.VS1[3][2];
        ctx.VA[3][1][iv] = ctx.VS2[3][1];
        ctx.VA[3][2][iv] = ctx.VS2[3][2];
        ctx.VDEL[3][2][iv] = ctx.LALFA
          ? ctx.VSR[3] * reClmr + ctx.VSM[3] * msqClmr
          : (ctx.VS1[3][4] * u1A + ctx.VS1[3][3] * d1A)
            + (ctx.VS2[3][4] * u2A + ctx.VS2[3][3] * d2A)
            + (ctx.VS1[3][5] + ctx.VS2[3][5] + ctx.VSX[3])
            * (xiUle1 * ule1A + xiUle2 * ule2A);
        ctx.VDEL[3][1][iv] = ctx.VSREZ[3]
          + (ctx.VS1[3][4] * due1 + ctx.VS1[3][3] * dds1)
          + (ctx.VS2[3][4] * due2 + ctx.VS2[3][3] * dds2)
          + (ctx.VS1[3][5] + ctx.VS2[3][5] + ctx.VSX[3])
          * (xiUle1 * dule1 + xiUle2 * dule2);

      }

      if (ctx.TRAN) {
        ctx.TURB = true;
        ctx.ITRAN[is] = ibl;
        ctx.TFORCE[is] = ctx.TRFORC;
        ctx.XSSITR[is] = ctx.XT;

        const str = is === 1 ? ctx.SST - ctx.XT : ctx.SST + ctx.XT;
        const chx = ctx.XTE - ctx.XLE;
        const chy = ctx.YTE - ctx.YLE;
        const chsq = chx ** 2 + chy ** 2;
        const s = ctx.S.subarray(1, ctx.N + 1);
        const x = ctx.X.subarray(1, ctx.N + 1);
        const y = ctx.Y.subarray(1, ctx.N + 1);
        const xp = ctx.XP.subarray(1, ctx.N + 1);
        const yp = ctx.YP.subarray(1, ctx.N + 1);
        const xtr = seval(str, x, xp, s, ctx.N);
        const ytr = seval(str, y, yp, s, ctx.N);
        ctx.XOCTR[is] = ((xtr - ctx.XLE) * chx + (ytr - ctx.YLE) * chy) / chsq;
        ctx.YOCTR[is] = ((ytr - ctx.YLE) * chx - (xtr - ctx.XLE) * chy) / chsq;
      }

      ctx.TRAN = false;

      if (ibl === ctx.IBLTE[is]) {
        ctx.TURB = true;
        ctx.WAKE = true;
        ctx.blvar(3, ctx);
        ctx.blmid(3, ctx);
      }

      for (let js = 1; js <= 2; js += 1) {
        for (let jbl = 2; jbl <= ctx.NBL[js]; jbl += 1) {
          const jv = ctx.ISYS[jbl][js];
          u1M[jv] = u2M[jv];
          d1M[jv] = d2M[jv];
        }
      }
      u1A = u2A;
      d1A = d2A;
      due1 = due2;
      dds1 = dds2;

      if (ibl === ctx.ITRAN[is] && ctx.X2 > ctx.X1) {
        if (is === 1) {
          ctx.TINDEX[is] = (ctx.IST - ctx.ITRAN[is] + 3) - (ctx.XT - ctx.X1) / (ctx.X2 - ctx.X1);
        } else {
          ctx.TINDEX[is] = (ctx.IST + ctx.ITRAN[is] - 2) + (ctx.XT - ctx.X1) / (ctx.X2 - ctx.X1);
        }
      }

      for (let icom = 1; icom <= ctx.NCOM; icom += 1) {
        ctx.COM1[icom] = ctx.COM2[icom];
      }
      if (typeof ctx.syncComToVars === 'function') {
        ctx.syncComToVars(ctx, 1);
      }
    }
  }

  // No additional diagnostics.
}

// Forced transition position on each side (XIFSET).
// Map x/c to arc-length position along each surface.
function xifset(ctx, is) {
  if (ctx.XSTRIP[is] >= 1.0) {
    ctx.XIFORC = ctx.XSSI[ctx.IBLTE[is]][is];
    return;
  }

  const chx = ctx.XTE - ctx.XLE;
  const chy = ctx.YTE - ctx.YLE;
  const chsq = chx ** 2 + chy ** 2;

  let dxNan = 0;
  let agNan = 0;
  for (let i = 1; i <= ctx.N; i += 1) {
    ctx.W1[i] = ((ctx.X[i] - ctx.XLE) * chx + (ctx.Y[i] - ctx.YLE) * chy) / chsq;
    ctx.W2[i] = ((ctx.Y[i] - ctx.YLE) * chx - (ctx.X[i] - ctx.XLE) * chy) / chsq;
  }

  const s1 = ctx.S.subarray(1, ctx.N + 1);
  const w1 = ctx.W1.subarray(1, ctx.N + 1);
  const w2 = ctx.W2.subarray(1, ctx.N + 1);
  const w3 = ctx.W3.subarray(1, ctx.N + 1);
  const w4 = ctx.W4.subarray(1, ctx.N + 1);
  splind(w1, w3, s1, ctx.N, -999.0, -999.0);
  splind(w2, w4, s1, ctx.N, -999.0, -999.0);

  if (is === 1) {
    let str = ctx.SLE + (ctx.S[1] - ctx.SLE) * ctx.XSTRIP[is];
    str = sinvrt(str, ctx.XSTRIP[is], w1, w3, s1, ctx.N);
    ctx.XIFORC = Math.min((ctx.SST - str), ctx.XSSI[ctx.IBLTE[is]][is]);
  } else {
    let str = ctx.SLE + (ctx.S[ctx.N] - ctx.SLE) * ctx.XSTRIP[is];
    str = sinvrt(str, ctx.XSTRIP[is], w1, w3, s1, ctx.N);
    ctx.XIFORC = Math.min((str - ctx.SST), ctx.XSSI[ctx.IBLTE[is]][is]);
  }

  if (ctx.XIFORC < 0.0) {
    ctx.XIFORC = ctx.XSSI[ctx.IBLTE[is]][is];
  }
}

// Finalize BL step: update state variables and derived quantities (UPDATE).
// Computes Newton deltas for (Ctau, Theta, m) and applies them to the BL state.
function update(ctx) {
  // Fortran comments (UPDATE) highlight:
  // - Solve for alpha/CL correction (DAC) with relaxation.
  // - Update BL variables and enforce limits.
  const dalmax = 0.5 * ctx.DTOR;
  const dalmin = -0.5 * ctx.DTOR;

  // Allowable alpha/CL step limits.
  let dclmin = -0.5;
  const dclmax = 0.5;
  if (ctx.MATYP !== 1) dclmin = Math.max(-0.5, -0.9 * ctx.CL);

  const hstinv = ctx.GAMM1 * (ctx.MINF / ctx.QINF) ** 2 / (1.0 + 0.5 * ctx.GAMM1 * ctx.MINF ** 2);

  const unew = Array.from({ length: ctx.IVX + 1 }, () => new Float64Array(3));
  const uAc = Array.from({ length: ctx.IVX + 1 }, () => new Float64Array(3));
  const qnew = new Float64Array(ctx.IQX + 1);
  const qAc = new Float64Array(ctx.IQX + 1);

  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
      const i = ctx.IPAN[ibl][is];

      let dui = 0.0;
      let duiAc = 0.0;
      for (let js = 1; js <= 2; js += 1) {
        for (let jbl = 2; jbl <= ctx.NBL[js]; jbl += 1) {
          const j = ctx.IPAN[jbl][js];
          const jv = ctx.ISYS[jbl][js];
          const ueM = -ctx.VTI[ibl][is] * ctx.VTI[jbl][js] * ctx.DIJ[i][j];
          dui += ueM * (ctx.MASS[jbl][js] + ctx.VDEL[3][1][jv]);
          duiAc += ueM * (-ctx.VDEL[3][2][jv]);
        }
      }

      let uinvAc = 0.0;
      if (!ctx.LALFA) {
        uinvAc = ctx.UINV_A[ibl][is];
      }

      unew[ibl][is] = ctx.UINV[ibl][is] + dui;
      uAc[ibl][is] = uinvAc + duiAc;
    }
  }

  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= ctx.IBLTE[is]; ibl += 1) {
      const i = ctx.IPAN[ibl][is];
      qnew[i] = ctx.VTI[ibl][is] * unew[ibl][is];
      qAc[i] = ctx.VTI[ibl][is] * uAc[ibl][is];
    }
  }

  const sa = Math.sin(ctx.ALFA);
  const ca = Math.cos(ctx.ALFA);

  const beta = Math.sqrt(1.0 - ctx.MINF ** 2);
  const betaMsq = -0.5 / beta;

  const bfac = 0.5 * ctx.MINF ** 2 / (1.0 + beta);
  const bfacMsq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * betaMsq;

  let clnew = 0.0;
  let clA = 0.0;
  let clMs = 0.0;
  let clAc = 0.0;

  let cpg1;
  let cpg1Ms;
  let cpg1Ac;

  {
    const i = 1;
    const cginc = 1.0 - (qnew[i] / ctx.QINF) ** 2;
    const den1 = beta + bfac * cginc;
    cpg1 = cginc / den1;
    cpg1Ms = -cpg1 / (beta + bfac * cginc) * (betaMsq + bfacMsq * cginc);
    const cpiQ = -2.0 * qnew[i] / ctx.QINF ** 2;
    const cpcCpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
    cpg1Ac = cpcCpi * cpiQ * qAc[i];
  }
  for (let i = 1; i <= ctx.N; i += 1) {
    let ip = i + 1;
    if (i === ctx.N) ip = 1;

    const cginc = 1.0 - (qnew[ip] / ctx.QINF) ** 2;
    const den2 = beta + bfac * cginc;
    const cpg2 = cginc / den2;
    const cpg2Ms = -cpg2 / (beta + bfac * cginc) * (betaMsq + bfacMsq * cginc);

    const cpiQ = -2.0 * qnew[ip] / ctx.QINF ** 2;
    const cpcCpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
    const cpg2Ac = cpcCpi * cpiQ * qAc[ip];

    const dx = (ctx.X[ip] - ctx.X[i]) * ca + (ctx.Y[ip] - ctx.Y[i]) * sa;
    const dxA = -(ctx.X[ip] - ctx.X[i]) * sa + (ctx.Y[ip] - ctx.Y[i]) * ca;

    const ag = 0.5 * (cpg2 + cpg1);
    const agMs = 0.5 * (cpg2Ms + cpg1Ms);
    const agAc = 0.5 * (cpg2Ac + cpg1Ac);

    clnew += dx * ag;
    clA += dxA * ag;
    clMs += dx * agMs;
    clAc += dx * agAc;

    cpg1 = cpg2;
    cpg1Ms = cpg2Ms;
    cpg1Ac = cpg2Ac;
  }

  let rlx = 1.0;
  let dac = 0.0;
  let denom = 0.0;

  if (ctx.LALFA) {
    denom = 1.0 - clAc - clMs * 2.0 * ctx.MINF * ctx.MINF_CL;
    dac = (clnew - ctx.CL) / denom;
    if (rlx * dac > dclmax) rlx = dclmax / dac;
    if (rlx * dac < dclmin) rlx = dclmin / dac;
  } else {
    denom = 0.0 - clAc - clA;
    dac = (clnew - ctx.CLSPEC) / denom;
    if (rlx * dac > dalmax) rlx = dalmax / dac;
    if (rlx * dac < dalmin) rlx = dalmin / dac;
  }
  ctx.DAC = dac;

  ctx.RMSBL = 0.0;
  ctx.RMXBL = 0.0;

  const dhi = 1.5;
  const dlo = -0.5;

  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
      const iv = ctx.ISYS[ibl][is];

      const dctau = ctx.VDEL[1][1][iv] - dac * ctx.VDEL[1][2][iv];
      const dthet = ctx.VDEL[2][1][iv] - dac * ctx.VDEL[2][2][iv];
      const dmass = ctx.VDEL[3][1][iv] - dac * ctx.VDEL[3][2][iv];
      const duedg = unew[ibl][is] + dac * uAc[ibl][is] - ctx.UEDG[ibl][is];
      const ddstr = (dmass - ctx.DSTR[ibl][is] * duedg) / ctx.UEDG[ibl][is];

      let dn1 = dctau / 10.0;
      if (ibl >= ctx.ITRAN[is]) dn1 = dctau / ctx.CTAU[ibl][is];
      const dn2 = dthet / ctx.THET[ibl][is];
      const dn3 = ddstr / ctx.DSTR[ibl][is];
      const dn4 = Math.abs(duedg) / 0.25;

      ctx.RMSBL += dn1 ** 2 + dn2 ** 2 + dn3 ** 2 + dn4 ** 2;

      const rdn1 = rlx * dn1;
      if (Math.abs(dn1) > Math.abs(ctx.RMXBL)) {
        ctx.RMXBL = dn1;
        ctx.VMXBL = ibl < ctx.ITRAN[is] ? 'n' : 'C';
        ctx.IMXBL = ibl;
        ctx.ISMXBL = is;
      }
      if (rdn1 > dhi) rlx = dhi / dn1;
      if (rdn1 < dlo) rlx = dlo / dn1;

      const rdn2 = rlx * dn2;
      if (Math.abs(dn2) > Math.abs(ctx.RMXBL)) {
        ctx.RMXBL = dn2;
        ctx.VMXBL = 'T';
        ctx.IMXBL = ibl;
        ctx.ISMXBL = is;
      }
      if (rdn2 > dhi) rlx = dhi / dn2;
      if (rdn2 < dlo) rlx = dlo / dn2;

      const rdn3 = rlx * dn3;
      if (Math.abs(dn3) > Math.abs(ctx.RMXBL)) {
        ctx.RMXBL = dn3;
        ctx.VMXBL = 'D';
        ctx.IMXBL = ibl;
        ctx.ISMXBL = is;
      }
      if (rdn3 > dhi) rlx = dhi / dn3;
      if (rdn3 < dlo) rlx = dlo / dn3;

      const rdn4 = rlx * dn4;
      if (Math.abs(dn4) > Math.abs(ctx.RMXBL)) {
        ctx.RMXBL = duedg;
        ctx.VMXBL = 'U';
        ctx.IMXBL = ibl;
        ctx.ISMXBL = is;
      }
      if (rdn4 > dhi) rlx = dhi / dn4;
      if (rdn4 < dlo) rlx = dlo / dn4;
    }
  }

  ctx.RMSBL = Math.sqrt(ctx.RMSBL / (4.0 * (ctx.NBL[1] + ctx.NBL[2])));
  ctx.RLX = rlx;

  if (ctx.LALFA) {
    ctx.CL = ctx.CL + rlx * dac;
  } else {
    ctx.ALFA = ctx.ALFA + rlx * dac;
    ctx.ADEG = ctx.ALFA / ctx.DTOR;
  }

  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
      const iv = ctx.ISYS[ibl][is];

      const dctau = ctx.VDEL[1][1][iv] - dac * ctx.VDEL[1][2][iv];
      const dthet = ctx.VDEL[2][1][iv] - dac * ctx.VDEL[2][2][iv];
      const dmass = ctx.VDEL[3][1][iv] - dac * ctx.VDEL[3][2][iv];
      const duedg = unew[ibl][is] + dac * uAc[ibl][is] - ctx.UEDG[ibl][is];
      const ddstr = (dmass - ctx.DSTR[ibl][is] * duedg) / ctx.UEDG[ibl][is];

      ctx.CTAU[ibl][is] = ctx.CTAU[ibl][is] + rlx * dctau;
      ctx.THET[ibl][is] = ctx.THET[ibl][is] + rlx * dthet;
      ctx.DSTR[ibl][is] = ctx.DSTR[ibl][is] + rlx * ddstr;
      ctx.UEDG[ibl][is] = ctx.UEDG[ibl][is] + rlx * duedg;

      let dswaki = 0.0;
      if (ibl > ctx.IBLTE[is]) {
        const iw = ibl - ctx.IBLTE[is];
        dswaki = ctx.WGAP[iw];
      }

      if (ibl >= ctx.ITRAN[is]) {
        ctx.CTAU[ibl][is] = Math.min(ctx.CTAU[ibl][is], 0.25);
      }

      let hklim = 0.0;
      if (ibl <= ctx.IBLTE[is]) {
        hklim = 1.02;
      } else {
        hklim = 1.00005;
      }
      const msq = ctx.UEDG[ibl][is] ** 2 * hstinv / (ctx.GAMM1 * (1.0 - 0.5 * ctx.UEDG[ibl][is] ** 2 * hstinv));
      let dsw = ctx.DSTR[ibl][is] - dswaki;
      dsw = dslim(ctx, dsw, ctx.THET[ibl][is], ctx.UEDG[ibl][is], msq, hklim);
      ctx.DSTR[ibl][is] = dsw + dswaki;

      ctx.MASS[ibl][is] = ctx.DSTR[ibl][is] * ctx.UEDG[ibl][is];
    }

    for (let ibl = 3; ibl <= ctx.IBLTE[is]; ibl += 1) {
      if (ctx.UEDG[ibl - 1][is] > 0.0 && ctx.UEDG[ibl][is] <= 0.0) {
        ctx.UEDG[ibl][is] = ctx.UEDG[ibl - 1][is];
        ctx.MASS[ibl][is] = ctx.DSTR[ibl][is] * ctx.UEDG[ibl][is];
      }
    }
  }

  for (let kbl = 1; kbl <= ctx.NBL[2] - ctx.IBLTE[2]; kbl += 1) {
    ctx.CTAU[ctx.IBLTE[1] + kbl][1] = ctx.CTAU[ctx.IBLTE[2] + kbl][2];
    ctx.THET[ctx.IBLTE[1] + kbl][1] = ctx.THET[ctx.IBLTE[2] + kbl][2];
    ctx.DSTR[ctx.IBLTE[1] + kbl][1] = ctx.DSTR[ctx.IBLTE[2] + kbl][2];
    ctx.UEDG[ctx.IBLTE[1] + kbl][1] = ctx.UEDG[ctx.IBLTE[2] + kbl][2];
    ctx.TAU[ctx.IBLTE[1] + kbl][1] = ctx.TAU[ctx.IBLTE[2] + kbl][2];
    ctx.DIS[ctx.IBLTE[1] + kbl][1] = ctx.DIS[ctx.IBLTE[2] + kbl][2];
    ctx.CTQ[ctx.IBLTE[1] + kbl][1] = ctx.CTQ[ctx.IBLTE[2] + kbl][2];
    ctx.DELT[ctx.IBLTE[1] + kbl][1] = ctx.DELT[ctx.IBLTE[2] + kbl][2];
    ctx.TSTR[ctx.IBLTE[1] + kbl][1] = ctx.TSTR[ctx.IBLTE[2] + kbl][2];
  }
}

function ueset(ctx) {
  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
      const i = ctx.IPAN[ibl][is];
      let dui = 0.0;
      for (let js = 1; js <= 2; js += 1) {
        for (let jbl = 2; jbl <= ctx.NBL[js]; jbl += 1) {
          const j = ctx.IPAN[jbl][js];
          const ueM = -ctx.VTI[ibl][is] * ctx.VTI[jbl][js] * ctx.DIJ[i][j];
          dui += ueM * ctx.MASS[jbl][js];
        }
      }
      ctx.UEDG[ibl][is] = ctx.UINV[ibl][is] + dui;
    }
  }
}

export {
  dslim,
  blpini,
  iblsys,
  mrcl,
  comset,
  mrchue,
  mrchdu,
  ueset,
  setbl,
  xifset,
  update,
};
