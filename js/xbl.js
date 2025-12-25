// Port of xbl.f (partial). Marching and utility routines for the integral BL.
// Integral boundary-layer method: marching, closures, and coupling terms.

import { splind, sinvrt, seval } from './spline.js';
import { gauss } from './xsolve.js';

// Debug guard for non-finite BL state; mirrors XFOIL's diagnostic checks.
function blCheck(ctx, label, extra) {
  if (!ctx.DEBUG_BL) return;
  const fields = {
    x1: ctx.X1,
    x2: ctx.X2,
    u1: ctx.U1,
    u2: ctx.U2,
    t1: ctx.T1,
    t2: ctx.T2,
    d1: ctx.D1,
    d2: ctx.D2,
    hk1: ctx.HK1,
    hk2: ctx.HK2,
    rt1: ctx.RT1,
    rt2: ctx.RT2,
    ampl1: ctx.AMPL1,
    ampl2: ctx.AMPL2,
  };
  const bad = Object.entries(fields).filter(([, v]) => !Number.isFinite(v));
  if (bad.length === 0) return;
  const record = { label, bad: bad.map(([k]) => k), ...fields, ...extra };
  if (!ctx.BL_TRACE) ctx.BL_TRACE = [];
  if (ctx.BL_TRACE.length < 50) ctx.BL_TRACE.push(record);
  console.warn('BLCHK:', record);
  if (ctx.DEBUG_BL_FAILFAST) {
    throw new Error(`BLCHK ${label}: non-finite ${bad.map(([k]) => k).join(', ')}`);
  }
}

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

// Displacement thickness limiter based on shape factor constraint.
function dslim(ctx, dstr, thet, uedg, msq, hklim) {
  const h = dstr / thet;
  const { hk, hkH } = ctx.hkin(h, msq);
  const dh = Math.max(0.0, hklim - hk) / hkH;
  return dstr + dh * thet;
}

// Initialize BL empirical constants (XFOIL defaults).
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

// Mach/Re updates for fixed CL or alpha families (operating condition map).
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

// Precompute coefficients and intermediate BL parameters for marching.
function comset(ctx) {
  const beta = Math.sqrt(Math.max(1.0 - ctx.MINF ** 2, 0.0));
  const betaMsq = beta === 0.0 ? 0.0 : -0.5 / beta;
  ctx.TKLAM = ctx.MINF ** 2 / (1.0 + beta) ** 2;
  ctx.TK_MSQ = 1.0 / (1.0 + beta) ** 2 - 2.0 * ctx.TKLAM / (1.0 + beta) * betaMsq;
}

// Initialize BL system matrices and coupling terms for marching.
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

// March Ue and related quantities along the surface and wake.
function mrchue(ctx) {
  const hlmax = 3.8;
  const htmax = 2.5;

  function resetThwaites(xsi, uei) {
    const ucon = uei / (xsi ** ctx.BULE);
    const safeUcon = ucon > 0.0 ? ucon : 1.0e-6;
    const tsq = 0.45 / (safeUcon * (5.0 * ctx.BULE + 1.0) * ctx.REYBL)
      * (xsi ** (1.0 - ctx.BULE));
    const thi = Math.sqrt(Math.max(tsq, 1.0e-12));
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

    ctx.TRAN = false;
    ctx.TURB = false;
    ctx.ITRAN[is] = ctx.IBLTE[is];
    ctx.THET[2][is] = thi;
    ctx.DSTR[2][is] = dsi;
    ctx.CTAU[2][is] = ami;
    ctx.MASS[2][is] = dsi * uei;

    for (ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
      const ibm = ibl - 1;

      const simi = ibl === 2;
      const wake = ibl > ctx.IBLTE[is];
      ctx.SIMI = simi;
      ctx.WAKE = wake;

      if (ibl > 2) {
        const xsiPrev = ctx.XSSI[ibm][is];
        const ueiPrev = ctx.UEDG[ibm][is];
        const thiPrev = ctx.THET[ibm][is];
        const dsiPrev = ctx.DSTR[ibm][is];
        if (ctx.DEBUG_BL && (!Number.isFinite(thiPrev) || !Number.isFinite(dsiPrev))) {
          console.warn('MRCHUE:prev-input non-finite', {
            is,
            ibl,
            ibm,
            xsiPrev,
            ueiPrev,
            thiPrev,
            dsiPrev,
          });
          if (ctx.DEBUG_BL_FAILFAST) {
            throw new Error('MRCHUE:prev-input non-finite');
          }
        }
        let amiPrev = ibm < ctx.ITRAN[is] ? ctx.CTAU[ibm][is] : 0.0;
        let ctiPrev = ibm >= ctx.ITRAN[is] ? ctx.CTAU[ibm][is] : 0.0;
        if (!Number.isFinite(amiPrev)) amiPrev = 0.0;
        if (!Number.isFinite(ctiPrev) || ctiPrev <= 0.0) ctiPrev = 0.03;
        let dswPrev = 0.0;
        if (ibm > ctx.IBLTE[is]) {
          const iwPrev = ibm - ctx.IBLTE[is];
          dswPrev = ctx.WGAP[iwPrev];
        }
        ctx.blprv(xsiPrev, amiPrev, ctiPrev, thiPrev, dsiPrev, dswPrev, ueiPrev, ctx);
        ctx.blkin(ctx);
        blCheck(ctx, 'MRCHUE:prev', { is, ibl, ibm });
        for (let icom = 1; icom <= ctx.NCOM; icom += 1) {
          ctx.COM1[icom] = ctx.COM2[icom];
        }
        if (typeof ctx.syncComToVars === 'function') {
          ctx.syncComToVars(ctx, 1);
        }
      }

      xsi = ctx.XSSI[ibl][is];
      uei = ctx.UEDG[ibl]?.[is];
      if (!Number.isFinite(uei) || uei <= 0.0) {
        const fallback = ctx.UEDG[ibm]?.[is];
        uei = Number.isFinite(fallback) && fallback > 0.0 ? fallback : 1.0e-6;
        if (ctx.UEDG[ibl]) ctx.UEDG[ibl][is] = uei;
      }

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
        if (!Number.isFinite(uei) || uei <= 0.0) {
          const uFromArray = ctx.UEDG[ibl]?.[is];
          const uFromPrev = ctx.UEDG[ibm]?.[is];
          if (Number.isFinite(uFromArray) && uFromArray > 0.0) {
            uei = uFromArray;
          } else if (Number.isFinite(uFromPrev) && uFromPrev > 0.0) {
            uei = uFromPrev;
          } else {
            uei = 1.0e-6;
          }
        }
        if (!Number.isFinite(thi) || !Number.isFinite(dsi) || thi <= 0.0 || dsi <= 0.0) {
          const reset = resetThwaites(xsi, uei);
          thi = reset.thi;
          dsi = reset.dsi;
        }
        ctx.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
        ctx.blkin(ctx);
        if (!Number.isFinite(ctx.HK2) || !Number.isFinite(ctx.RT2)) {
          console.warn('BL: non-finite after BLKIN', {
            is,
            ibl,
            xsi,
            uei,
            thi,
            dsi,
            hk2: ctx.HK2,
            rt2: ctx.RT2,
            u2: ctx.U2,
            t2: ctx.T2,
            d2: ctx.D2,
          });
          break;
        }

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

        if (ctx.DEBUG_BL && ibl === 2 && !ctx._dbgIbl2Checked) {
          let badVs = false;
          for (let k = 1; k <= 4; k += 1) {
            if (!Number.isFinite(ctx.VSREZ[k])) badVs = true;
            for (let l = 1; l <= 4; l += 1) {
              if (!Number.isFinite(ctx.VS2[k][l])) badVs = true;
            }
          }
          if (badVs) {
            console.warn('MRCHUE:ibl2 bad VS', {
              is,
              ibl,
              itbl,
              thi,
              dsi,
              uei,
              hk2: ctx.HK2,
              rt2: ctx.RT2,
              vsrez: [ctx.VSREZ[1], ctx.VSREZ[2], ctx.VSREZ[3], ctx.VSREZ[4]],
              vs2: [
                [ctx.VS2[1][1], ctx.VS2[1][2], ctx.VS2[1][3], ctx.VS2[1][4]],
                [ctx.VS2[2][1], ctx.VS2[2][2], ctx.VS2[2][3], ctx.VS2[2][4]],
                [ctx.VS2[3][1], ctx.VS2[3][2], ctx.VS2[3][3], ctx.VS2[3][4]],
                [ctx.VS2[4][1], ctx.VS2[4][2], ctx.VS2[4][3], ctx.VS2[4][4]],
              ],
            });
            if (ctx.DEBUG_BL_FAILFAST) {
              throw new Error('MRCHUE:ibl2 bad VS');
            }
          }
          ctx._dbgIbl2Checked = true;
        }

        if (direct) {
          ctx.VS2[4][1] = 0.0;
          ctx.VS2[4][2] = 0.0;
          ctx.VS2[4][3] = 0.0;
          ctx.VS2[4][4] = 1.0;
          ctx.VSREZ[4] = 0.0;

          gauss1(4, ctx.VS2, ctx.VSREZ);

          if (ctx.DEBUG_BL && ibl === 2 && !Number.isFinite(ctx.VSREZ[3])) {
            console.warn('MRCHUE:ibl2 bad VSREZ', {
              is,
              ibl,
              itbl,
              thi,
              dsi,
              uei,
              vsrez: [ctx.VSREZ[1], ctx.VSREZ[2], ctx.VSREZ[3], ctx.VSREZ[4]],
            });
            if (ctx.DEBUG_BL_FAILFAST) {
              throw new Error('MRCHUE:ibl2 bad VSREZ');
            }
          }

          dmax = Math.max(Math.abs(ctx.VSREZ[2] / thi), Math.abs(ctx.VSREZ[3] / dsi));
          if (ibl < ctx.ITRAN[is]) dmax = Math.max(dmax, Math.abs(ctx.VSREZ[1] / 10.0));
          if (ibl >= ctx.ITRAN[is]) dmax = Math.max(dmax, Math.abs(ctx.VSREZ[1] / cti));

          let rlx = 1.0;
          if (dmax > 0.3) rlx = 0.3 / dmax;

          if (ibl !== ctx.IBLTE[is] + 1) {
            msq = uei * uei * ctx.HSTINV / (ctx.GM1BL * (1.0 - 0.5 * uei * uei * ctx.HSTINV));
            const htest = (dsi + rlx * ctx.VSREZ[3]) / (thi + rlx * ctx.VSREZ[2]);
            const hkin = ctx.hkin(htest, msq);
            const hktest = hkin.hk;

            let hmax = htmax;
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

      if (ctx.DEBUG_BL && ibl === 2) {
        const th2 = ctx.THET[ibl][is];
        const ds2 = ctx.DSTR[ibl][is];
        if (!Number.isFinite(th2) || !Number.isFinite(ds2)) {
          console.warn('MRCHUE:ibl2 non-finite', { is, ibl, th2, ds2 });
          if (ctx.DEBUG_BL_FAILFAST) {
            throw new Error('MRCHUE:ibl2 non-finite');
          }
        }
      }

      ctx.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, ctx);
      ctx.blkin(ctx);
      for (let icom = 1; icom <= ctx.NCOM; icom += 1) {
        ctx.COM1[icom] = ctx.COM2[icom];
      }
      if (typeof ctx.syncComToVars === 'function') {
        ctx.syncComToVars(ctx, 1);
      }

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

function mrchdu(ctx) {
  const deps = 5.0e-6;
  const senswt = 1000.0;

  for (let is = 1; is <= 2; is += 1) {
    ctx.AMCRIT = ctx.ACRIT[is];
    xifset(ctx, is);

    let ibl = 2;
    let xsi = ctx.XSSI[ibl][is];
    let uei = ctx.UEDG[ibl][is];
    ctx.BULE = 1.0;

    const itrold = ctx.ITRAN[is];

    ctx.TRAN = false;
    ctx.TURB = false;
    ctx.ITRAN[is] = ctx.IBLTE[is];

    let sens = 0.0;

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
      if (ctx.DEBUG_BL && ibl === 2 && !Number.isFinite(uei)) {
        console.warn('MRCHDU:ibl2 non-finite UEDG', { is, ibl, uei });
        if (ctx.DEBUG_BL_FAILFAST) {
          throw new Error('MRCHDU:ibl2 non-finite UEDG');
        }
      }
      let thi = ctx.THET[ibl][is];
      let dsi = ctx.DSTR[ibl][is];
      if (ctx.DEBUG_BL && ibl === 2 && (!Number.isFinite(thi) || !Number.isFinite(dsi))) {
        console.warn('MRCHDU:ibl2 non-finite THET/DSTR', { is, ibl, thi, dsi });
        if (ctx.DEBUG_BL_FAILFAST) {
          throw new Error('MRCHDU:ibl2 non-finite THET/DSTR');
        }
      }

      let ami = 0.0;
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

      if (ctx.TRAN || ibl === ctx.IBLTE[is]) {
        ctx.TURB = true;
        ctx.TFORCE[is] = ctx.TRFORC;
        ctx.XSSITR[is] = ctx.XT;
      }

      ctx.TRAN = false;
    }
  }
}

// Assemble current BL state into system form prior to Newton iteration.
function setbl(ctx) {
  const nsys = ctx.NSYS;
  const COM1 = ctx.COM1;
  const COM2 = ctx.COM2;
  const VA = ctx.VA;
  const VB = ctx.VB;
  const VDEL = ctx.VDEL;
  const VM = ctx.VM;
  COM1.fill(0.0, 1);
  COM2.fill(0.0, 1);
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

  const UINV = ctx.UINV;
  const UEDG = ctx.UEDG;
  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= NBL[is]; ibl += 1) {
      const uinv = UINV[ibl][is];
      let ue = UEDG[ibl][is];
      if (!Number.isFinite(ue) || ue <= 0.0) {
        ue = Number.isFinite(uinv) ? Math.abs(uinv) : 1.0e-6;
        UEDG[ibl][is] = ue > 0.0 ? ue : 1.0e-6;
      }
    }
  }

  if (!ctx.LBLINI) {
    mrchue(ctx);
    ctx.LBLINI = true;
  }

  mrchdu(ctx);

  const CTAU = ctx.CTAU;
  const ITRAN = ctx.ITRAN;
  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= NBL[is]; ibl += 1) {
      if (!Number.isFinite(CTAU[ibl][is])) {
        CTAU[ibl][is] = ibl < ITRAN[is] ? 0.0 : 0.03;
      }
    }
  }

  for (let is = 1; is <= 2; is += 1) {
    for (let ibl = 2; ibl <= NBL[is]; ibl += 1) {
      let uei = UEDG[ibl][is];
      if (!Number.isFinite(uei) || uei <= 0.0) {
        const uinv = UINV[ibl][is];
        uei = Number.isFinite(uinv) ? Math.abs(uinv) : 1.0e-6;
        UEDG[ibl][is] = uei > 0.0 ? uei : 1.0e-6;
      }
      let thi = ctx.THET[ibl][is];
      let dsi = ctx.DSTR[ibl][is];
      if (!Number.isFinite(thi) || !Number.isFinite(dsi) || thi <= 0.0 || dsi <= 0.0) {
        const xsi = ctx.XSSI[ibl][is];
        const ucon = uei / (xsi ** 1.0);
        const tsq = 0.45 / (Math.max(ucon, 1.0e-6) * 6.0 * ctx.REYBL) * xsi ** 0.0;
        thi = Math.sqrt(Math.max(tsq, 1.0e-12));
        dsi = 2.2 * thi;
        ctx.THET[ibl][is] = thi;
        ctx.DSTR[ibl][is] = dsi;
      }
      ctx.MASS[ibl][is] = ctx.DSTR[ibl][is] * ctx.UEDG[ibl][is];
    }
  }

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
    console.warn('SETBL: UESET produced non-finite UEDG', { badUeSet });
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

    for (let ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
      const iv = ctx.ISYS[ibl][is];

      ctx.SIMI = ibl === 2;
      ctx.WAKE = ibl > ctx.IBLTE[is];
      ctx.TRAN = ibl === ctx.ITRAN[is];
      ctx.TURB = ibl > ctx.ITRAN[is];

      const i = ctx.IPAN[ibl][is];

      const xsi = ctx.XSSI[ibl][is];
      let ami = ibl < ctx.ITRAN[is] ? ctx.CTAU[ibl][is] : 0.0;
      let cti = ibl >= ctx.ITRAN[is] ? ctx.CTAU[ibl][is] : 0.0;
      if (!Number.isFinite(ami)) ami = 0.0;
      if (!Number.isFinite(cti)) cti = 0.03;
      let uei = ctx.UEDG[ibl][is];
      let thi = ctx.THET[ibl][is];
      let mdi = ctx.MASS[ibl][is];
      let ueiSafe = uei;
      if (!Number.isFinite(ueiSafe) || ueiSafe <= 0.0) {
        const uinv = ctx.UINV[ibl][is];
        ueiSafe = Number.isFinite(uinv) ? Math.abs(uinv) : 1.0e-6;
        ctx.UEDG[ibl][is] = ueiSafe;
      }
      const ueiUse = Number.isFinite(ueiSafe) && ueiSafe > 0.0 ? ueiSafe : 1.0e-6;
      uei = ueiUse;
      let dsi = mdi / uei;
      if (!Number.isFinite(thi) || !Number.isFinite(dsi) || thi <= 0.0 || dsi <= 0.0) {
        const xsi = ctx.XSSI[ibl][is];
        const ucon = ueiSafe / (xsi ** 1.0);
        const tsq = 0.45 / (Math.max(ucon, 1.0e-6) * 6.0 * ctx.REYBL) * xsi ** 0.0;
        thi = Math.sqrt(Math.max(tsq, 1.0e-12));
        dsi = 2.2 * thi;
        ctx.THET[ibl][is] = thi;
        ctx.DSTR[ibl][is] = dsi;
      }
      if (!Number.isFinite(mdi)) {
        mdi = dsi * ueiSafe;
        ctx.MASS[ibl][is] = mdi;
      }
      dsi = mdi / uei;

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

      ctx.blprv(xsi, ami, cti, thi, dsi, dswaki, ueiSafe, ctx);
      ctx.blkin(ctx);
      blCheck(ctx, 'SETBL:cur', { is, ibl });

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

        if (!Number.isFinite(ctx.VS2[1][1])) {
          console.warn('SETBL: non-finite VS2(1,1)', {
            is,
            ibl,
            simi: ctx.SIMI,
            tran: ctx.TRAN,
            turb: ctx.TURB,
            x1: ctx.X1,
            x2: ctx.X2,
            t1: ctx.T1,
            t2: ctx.T2,
            d1: ctx.D1,
            d2: ctx.D2,
            u1: ctx.U1,
            u2: ctx.U2,
            hk1: ctx.HK1,
            hk2: ctx.HK2,
            rt1: ctx.RT1,
            rt2: ctx.RT2,
            ampl1: ctx.AMPL1,
            ampl2: ctx.AMPL2,
            s1: ctx.S1,
            s2: ctx.S2,
          });
        }

        if (!Number.isFinite(ctx.VS2[1][1])) {
          console.warn('SETBL: non-finite VS2(1,1)', {
            is,
            ibl,
            simi: ctx.SIMI,
            tran: ctx.TRAN,
            turb: ctx.TURB,
            x1: ctx.X1,
            x2: ctx.X2,
            t1: ctx.T1,
            t2: ctx.T2,
            d1: ctx.D1,
            d2: ctx.D2,
            u1: ctx.U1,
            u2: ctx.U2,
            hk1: ctx.HK1,
            hk2: ctx.HK2,
            rt1: ctx.RT1,
            rt2: ctx.RT2,
            ampl1: ctx.AMPL1,
            ampl2: ctx.AMPL2,
            s1: ctx.S1,
            s2: ctx.S2,
          });
        }

        if (!Number.isFinite(ctx.VS2[1][1])) {
          console.warn('SETBL: non-finite VS2(1,1)', {
            is,
            ibl,
            simi: ctx.SIMI,
            tran: ctx.TRAN,
            turb: ctx.TURB,
            x1: ctx.X1,
            x2: ctx.X2,
            t1: ctx.T1,
            t2: ctx.T2,
            d1: ctx.D1,
            d2: ctx.D2,
            u1: ctx.U1,
            u2: ctx.U2,
            hk1: ctx.HK1,
            hk2: ctx.HK2,
            rt1: ctx.RT1,
            rt2: ctx.RT2,
            ampl1: ctx.AMPL1,
            ampl2: ctx.AMPL2,
            s1: ctx.S1,
            s2: ctx.S2,
          });
        }

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

        if (!Number.isFinite(ctx.VS2[2][2]) || !Number.isFinite(ctx.HK2) || !Number.isFinite(ctx.RT2)) {
          console.warn('SETBL: non-finite system terms', {
            is,
            ibl,
            xsi: ctx.X2,
            uei: ctx.U2,
            t2: ctx.T2,
            d2: ctx.D2,
            hk2: ctx.HK2,
            rt2: ctx.RT2,
            vs22: ctx.VS2[2][2],
          });
        }
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
        const xtr = seval(str, ctx.X, ctx.XP, ctx.S, ctx.N);
        const ytr = seval(str, ctx.Y, ctx.YP, ctx.S, ctx.N);
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

  let badVa = 0;
  for (let iv = 1; iv <= nsys; iv += 1) {
    if (!Number.isFinite(VA[1][1][iv])) {
      badVa += 1;
      if (badVa <= 3) {
        console.warn('SETBL: non-finite VA11', {
          iv,
          is: ivToIs[iv],
          ibl: ivToIbl[iv],
          va11: VA[1][1][iv],
          vb11: VB[1][1][iv],
          va22: VA[2][2][iv],
          vm33: VM[3][iv][iv],
        });
      }
    }
  }
  if (badVa > 0) {
    console.warn('SETBL: non-finite VA11 count', { badVa });
  }
}

// Solve BL system with Newton iterations (laminar/turbulent switching).
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

function xifset(ctx, is) {
  if (ctx.XSTRIP[is] >= 1.0) {
    ctx.XIFORC = ctx.XSSI[ctx.IBLTE[is]][is];
    return;
  }

  const chx = ctx.XTE - ctx.XLE;
  const chy = ctx.YTE - ctx.YLE;
  const chsq = chx ** 2 + chy ** 2;

  for (let i = 1; i <= ctx.N; i += 1) {
    ctx.W1[i] = ((ctx.X[i] - ctx.XLE) * chx + (ctx.Y[i] - ctx.YLE) * chy) / chsq;
    ctx.W2[i] = ((ctx.Y[i] - ctx.YLE) * chx - (ctx.X[i] - ctx.XLE) * chy) / chsq;
  }

  splind(ctx.W1, ctx.W3, ctx.S, ctx.N, -999.0, -999.0);
  splind(ctx.W2, ctx.W4, ctx.S, ctx.N, -999.0, -999.0);

  if (is === 1) {
    let str = ctx.SLE + (ctx.S[1] - ctx.SLE) * ctx.XSTRIP[is];
    str = sinvrt(str, ctx.XSTRIP[is], ctx.W1, ctx.W3, ctx.S, ctx.N);
    ctx.XIFORC = Math.min((ctx.SST - str), ctx.XSSI[ctx.IBLTE[is]][is]);
  } else {
    let str = ctx.SLE + (ctx.S[ctx.N] - ctx.SLE) * ctx.XSTRIP[is];
    str = sinvrt(str, ctx.XSTRIP[is], ctx.W1, ctx.W3, ctx.S, ctx.N);
    ctx.XIFORC = Math.min((str - ctx.SST), ctx.XSSI[ctx.IBLTE[is]][is]);
  }

  if (ctx.XIFORC < 0.0) {
    ctx.XIFORC = ctx.XSSI[ctx.IBLTE[is]][is];
  }
}

// Finalize BL step: update state variables and derived quantities.
function update(ctx) {
  const dalmax = 0.5 * ctx.DTOR;
  const dalmin = -0.5 * ctx.DTOR;

  let dclmin = -0.5;
  const dclmax = 0.5;
  if (ctx.MATYP !== 1) dclmin = Math.max(-0.5, -0.9 * ctx.CL);

  const hstinv = ctx.GAMM1 * (ctx.MINF / ctx.QINF) ** 2 / (1.0 + 0.5 * ctx.GAMM1 * ctx.MINF ** 2);

  const unew = Array.from({ length: ctx.IVX + 1 }, () => new Float64Array(3));
  const uAc = Array.from({ length: ctx.IVX + 1 }, () => new Float64Array(3));
  const qnew = new Float64Array(ctx.IQX + 1);
  const qAc = new Float64Array(ctx.IQX + 1);

  if (ctx.DEBUG_BL) {
    let badMass = 0;
    let badTh = 0;
    let badDs = 0;
    let badUe = 0;
    for (let is = 1; is <= 2; is += 1) {
      for (let ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
        if (!Number.isFinite(ctx.MASS[ibl][is])) badMass += 1;
        if (!Number.isFinite(ctx.THET[ibl][is])) badTh += 1;
        if (!Number.isFinite(ctx.DSTR[ibl][is])) badDs += 1;
        if (!Number.isFinite(ctx.UEDG[ibl][is])) badUe += 1;
      }
    }
    if (badMass || badTh || badDs || badUe) {
      console.warn('BL: non-finite state pre-update', {
        badMass,
        badTh,
        badDs,
        badUe,
      });
      if (ctx.DEBUG_BL_FAILFAST) {
        throw new Error('BL: non-finite state pre-update');
      }
    }
  }

  if (ctx.DEBUG_BL) {
    let badVdel = 0;
    for (let is = 1; is <= 2; is += 1) {
      for (let ibl = 2; ibl <= ctx.NBL[is]; ibl += 1) {
        const iv = ctx.ISYS[ibl][is];
        if (!Number.isFinite(ctx.VDEL[1][1][iv]) || !Number.isFinite(ctx.VDEL[2][1][iv])
          || !Number.isFinite(ctx.VDEL[3][1][iv])) {
          badVdel += 1;
        }
      }
    }
    if (badVdel > 0) {
      console.warn('BL: non-finite VDEL', { badVdel });
      if (ctx.DEBUG_BL_FAILFAST) {
        throw new Error('BL: non-finite VDEL');
      }
    }
  }

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
    cpg1 = cginc / (beta + bfac * cginc);
    cpg1Ms = -cpg1 / (beta + bfac * cginc) * (betaMsq + bfacMsq * cginc);
    const cpiQ = -2.0 * qnew[i] / ctx.QINF ** 2;
    const cpcCpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
    cpg1Ac = cpcCpi * cpiQ * qAc[i];
  }

  for (let i = 1; i <= ctx.N; i += 1) {
    let ip = i + 1;
    if (i === ctx.N) ip = 1;

    const cginc = 1.0 - (qnew[ip] / ctx.QINF) ** 2;
    const cpg2 = cginc / (beta + bfac * cginc);
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

  if (ctx.LALFA) {
    dac = (clnew - ctx.CL) / (1.0 - clAc - clMs * 2.0 * ctx.MINF * ctx.MINF_CL);
    if (rlx * dac > dclmax) rlx = dclmax / dac;
    if (rlx * dac < dclmin) rlx = dclmin / dac;
  } else {
    dac = (clnew - ctx.CLSPEC) / (0.0 - clAc - clA);
    if (rlx * dac > dalmax) rlx = dalmax / dac;
    if (rlx * dac < dalmin) rlx = dalmin / dac;
  }

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
      let ue = ctx.UINV[ibl][is] + dui;
      if (!Number.isFinite(ue) || ue <= 0.0) ue = 1.0e-6;
      ctx.UEDG[ibl][is] = ue;
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
  blsolv,
  xifset,
  update,
};
