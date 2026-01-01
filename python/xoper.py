# Ported from XFOIL Fortran source (Mark Drela).
# This file is a derived work and remains under the terms of the
# GNU General Public License v2 or later.
# See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

import math
from .spline import sinvrt, seval
from .xfoil import cpcalc, clcalc, cdcalc, comset, mrcl
from .xgdes import getxyf
from .xpanel import gamqv, iblpan, qdcalc, qiset, qvfue, qwcalc, stfind, stmove, uicalc, xicalc, xywake
from .xsolve import blsolv
from .xbl import iblsys, setbl, update
from .xblsys import hkin


def _format_fixed(value, width, decimals):
    text = f"{value:.{decimals}f}"
    if len(text) < width:
        text = (" " * (width - len(text))) + text
    return text


def _bstrip(line):
    stripped = line.strip()
    return stripped.replace(" ", "")


def cpdump(ctx, fname1, kdelim=1):
    if not fname1 or not str(fname1).strip():
        raise ValueError("CPDUMP requires a filename.")
    delim = " "
    if kdelim == 1:
        delim = ","
    elif kdelim == 2:
        delim = "\t"
    elif kdelim != 0:
        print("? Illegal delimiter.  Using blank.")
        delim = " "

    with open(fname1, "w", encoding="ascii") as lu:
        if kdelim == 0:
            lu.write("#      x          Cp  \n")
        else:
            lu.write(f"#x{delim}Cp\n")

        comset(ctx)
        beta = math.sqrt(1.0 - ctx.MINF**2)
        bfac = 0.5 * ctx.MINF**2 / (1.0 + beta)

        for i in range(1, ctx.N + 1):
            cpinc = 1.0 - (ctx.GAM[i] / ctx.QINF) ** 2
            den = beta + bfac * cpinc
            cpcom = cpinc / den
            if kdelim == 0:
                line = f" {_format_fixed(ctx.X[i], 11, 5)}{_format_fixed(cpcom, 11, 5)}"
                lu.write(f"{line}\n")
            else:
                line = f" {_format_fixed(ctx.X[i], 11, 5)}{delim}{_format_fixed(cpcom, 11, 5)}{delim}"
                lu.write(f"{_bstrip(line)}\n")


def bldump(ctx, fname1, kdelim=1):
    if not fname1 or not str(fname1).strip():
        raise ValueError("BLDUMP requires a filename.")
    delim = " "
    if kdelim == 1:
        delim = ","
    elif kdelim == 2:
        delim = "\t"
    elif kdelim != 0:
        print("? Illegal delimiter.  Using blank.")
        delim = " "

    with open(fname1, "w", encoding="ascii") as lu:
        if kdelim == 0:
            lu.write(
                "#    s        x        y     Ue/Vinf    Dstar     Theta      Cf       H       H*        P         m          K          tau         Di\n"
            )
        else:
            lu.write(f"#s{delim}x{delim}y{delim}Ue/Vinf{delim}Dstar{delim}Theta{delim}Cf{delim}H\n")

        comset(ctx)
        hstinv = ctx.GAMM1 * (ctx.MINF / ctx.QINF) ** 2 / (1.0 + 0.5 * ctx.GAMM1 * ctx.MINF**2)

        for i in range(1, ctx.N + 1):
            is_ = 1
            if ctx.GAM[i] < 0.0:
                is_ = 2

            if ctx.LIPAN and ctx.LVISC:
                if is_ == 1:
                    ibl = ctx.IBLTE[is_] - i + 1
                else:
                    ibl = ctx.IBLTE[is_] + i - ctx.N
                ds = ctx.DSTR[ibl][is_]
                th = ctx.THET[ibl][is_]
                ts = ctx.TSTR[ibl][is_]
                cf = ctx.TAU[ibl][is_] / (0.5 * ctx.QINF**2)
                if th == 0.0:
                    h = 1.0
                    hs = 1.0
                else:
                    h = ds / th
                    hs = ts / th
            else:
                ds = 0.0
                th = 0.0
                ts = 0.0
                cf = 0.0
                h = 1.0
                hs = 2.0

            ue = (ctx.GAM[i] / ctx.QINF) * (1.0 - ctx.TKLAM) / (1.0 - ctx.TKLAM * (ctx.GAM[i] / ctx.QINF) ** 2)
            amsq = ue * ue * hstinv / (ctx.GAMM1 * (1.0 - 0.5 * ue * ue * hstinv))
            hk, _, _ = hkin(h, amsq)

            if kdelim == 0:
                line = (
                    f" {_format_fixed(ctx.S[i], 9, 5)}"
                    f"{_format_fixed(ctx.X[i], 9, 5)}"
                    f"{_format_fixed(ctx.Y[i], 9, 5)}"
                    f"{_format_fixed(ue, 9, 5)}"
                    f"{_format_fixed(ds, 10, 6)}"
                    f"{_format_fixed(th, 10, 6)}"
                    f"{_format_fixed(cf, 10, 6)}"
                    f"{_format_fixed(hk, 10, 4)}"
                    f"{_format_fixed(hs, 10, 4)}"
                    f"{_format_fixed(th * ue**2, 9, 5)}"
                    f"{_format_fixed(ds * ue, 9, 5)}"
                    f"{_format_fixed(ts * ue**3, 9, 5)}"
                )
                lu.write(f"{line}\n")
            else:
                line = (
                    f" {_format_fixed(ctx.S[i], 9, 5)}{delim}"
                    f"{_format_fixed(ctx.X[i], 9, 5)}{delim}"
                    f"{_format_fixed(ctx.Y[i], 9, 5)}{delim}"
                    f"{_format_fixed(ue, 9, 5)}{delim}"
                    f"{_format_fixed(ds, 10, 6)}{delim}"
                    f"{_format_fixed(th, 10, 6)}{delim}"
                    f"{_format_fixed(cf, 10, 6)}{delim}"
                    f"{_format_fixed(hk, 10, 4)}"
                )
                lu.write(f"{_bstrip(line)}\n")

        if ctx.LWAKE:
            is_ = 2
            for i in range(ctx.N + 1, ctx.N + ctx.NW + 1):
                ibl = ctx.IBLTE[is_] + i - ctx.N
                ds = ctx.DSTR[ibl][is_]
                th = ctx.THET[ibl][is_]
                h = ds / th
                cf = 0.0
                ui = ctx.UEDG[ibl][is_]
                ue = (ui / ctx.QINF) * (1.0 - ctx.TKLAM) / (1.0 - ctx.TKLAM * (ui / ctx.QINF) ** 2)
                amsq = ue * ue * hstinv / (ctx.GAMM1 * (1.0 - 0.5 * ue * ue * hstinv))
                hk, _, _ = hkin(h, amsq)

                if kdelim == 0:
                    line = (
                        f" {_format_fixed(ctx.S[i], 9, 5)}"
                        f"{_format_fixed(ctx.X[i], 9, 5)}"
                        f"{_format_fixed(ctx.Y[i], 9, 5)}"
                        f"{_format_fixed(ue, 9, 5)}"
                        f"{_format_fixed(ds, 10, 6)}"
                        f"{_format_fixed(th, 10, 6)}"
                        f"{_format_fixed(cf, 10, 6)}"
                        f"{_format_fixed(hk, 10, 4)}"
                    )
                    lu.write(f"{line}\n")
                else:
                    line = (
                        f" {_format_fixed(ctx.S[i], 9, 5)}{delim}"
                        f"{_format_fixed(ctx.X[i], 9, 5)}{delim}"
                        f"{_format_fixed(ctx.Y[i], 9, 5)}{delim}"
                        f"{_format_fixed(ue, 9, 5)}{delim}"
                        f"{_format_fixed(ds, 10, 6)}{delim}"
                        f"{_format_fixed(th, 10, 6)}{delim}"
                        f"{_format_fixed(cf, 10, 6)}{delim}"
                        f"{_format_fixed(hk, 10, 4)}"
                    )
                    lu.write(f"{_bstrip(line)}\n")


def mhinge(ctx):
    if not ctx.LFLAP:
        tops, bots, ctx.XOF, ctx.YOF = getxyf(ctx.X, ctx.XP, ctx.Y, ctx.YP, ctx.S, ctx.N, 0.0, 0.0, ctx.XOF, ctx.YOF)
        ctx.LFLAP = True
    else:
        tops = ctx.XOF
        bots = ctx.S[ctx.N] - ctx.XOF
        tops = sinvrt(tops, ctx.XOF, ctx.X, ctx.XP, ctx.S, ctx.N)
        bots = sinvrt(bots, ctx.XOF, ctx.X, ctx.XP, ctx.S, ctx.N)

    topx = seval(tops, ctx.X, ctx.XP, ctx.S, ctx.N)
    topy = seval(tops, ctx.Y, ctx.YP, ctx.S, ctx.N)
    botx = seval(bots, ctx.X, ctx.XP, ctx.S, ctx.N)
    boty = seval(bots, ctx.Y, ctx.YP, ctx.S, ctx.N)

    ctx.HMOM = 0.0
    ctx.HFX = 0.0
    ctx.HFY = 0.0

    for i in range(2, ctx.N + 1):
        if ctx.S[i - 1] >= tops and ctx.S[i] <= bots:
            continue

        dx = ctx.X[i] - ctx.X[i - 1]
        dy = ctx.Y[i] - ctx.Y[i - 1]
        xmid = 0.5 * (ctx.X[i] + ctx.X[i - 1]) - ctx.XOF
        ymid = 0.5 * (ctx.Y[i] + ctx.Y[i - 1]) - ctx.YOF
        if ctx.LVISC:
            pmid = 0.5 * (ctx.CPV[i] + ctx.CPV[i - 1])
        else:
            pmid = 0.5 * (ctx.CPI[i] + ctx.CPI[i - 1])
        ctx.HMOM = ctx.HMOM + pmid * (xmid * dx + ymid * dy)
        ctx.HFX = ctx.HFX - pmid * dy
        ctx.HFY = ctx.HFY + pmid * dx

    for i in range(2, ctx.N + 1):
        if ctx.S[i] > tops:
            break

    dx = topx - ctx.X[i - 1]
    dy = topy - ctx.Y[i - 1]
    xmid = 0.5 * (topx + ctx.X[i - 1]) - ctx.XOF
    ymid = 0.5 * (topy + ctx.Y[i - 1]) - ctx.YOF
    if ctx.S[i] != ctx.S[i - 1]:
        frac = (tops - ctx.S[i - 1]) / (ctx.S[i] - ctx.S[i - 1])
    else:
        frac = 0.0
    if ctx.LVISC:
        topp = ctx.CPV[i] * frac + ctx.CPV[i - 1] * (1.0 - frac)
        pmid = 0.5 * (topp + ctx.CPV[i - 1])
    else:
        topp = ctx.CPI[i] * frac + ctx.CPI[i - 1] * (1.0 - frac)
        pmid = 0.5 * (topp + ctx.CPI[i - 1])
    ctx.HMOM = ctx.HMOM + pmid * (xmid * dx + ymid * dy)
    ctx.HFX = ctx.HFX - pmid * dy
    ctx.HFY = ctx.HFY + pmid * dx

    dx = ctx.XOF - topx
    dy = ctx.YOF - topy
    xmid = 0.5 * (topx + ctx.XOF) - ctx.XOF
    ymid = 0.5 * (topy + ctx.YOF) - ctx.YOF
    ctx.HMOM = ctx.HMOM + pmid * (xmid * dx + ymid * dy)
    ctx.HFX = ctx.HFX - pmid * dy
    ctx.HFY = ctx.HFY + pmid * dx

    for i in range(ctx.N, 1, -1):
        if ctx.S[i - 1] < bots:
            break

    dx = ctx.X[i] - botx
    dy = ctx.Y[i] - boty
    xmid = 0.5 * (ctx.X[i] + botx) - ctx.XOF
    ymid = 0.5 * (ctx.Y[i] + boty) - ctx.YOF
    if ctx.S[i] != ctx.S[i - 1]:
        frac = (ctx.S[i] - bots) / (ctx.S[i] - ctx.S[i - 1])
    else:
        frac = 0.0
    if ctx.LVISC:
        botp = ctx.CPV[i - 1] * frac + ctx.CPV[i] * (1.0 - frac)
        pmid = 0.5 * (botp + ctx.CPV[i])
    else:
        botp = ctx.CPI[i - 1] * frac + ctx.CPI[i] * (1.0 - frac)
        pmid = 0.5 * (botp + ctx.CPI[i])
    ctx.HMOM = ctx.HMOM + pmid * (xmid * dx + ymid * dy)
    ctx.HFX = ctx.HFX - pmid * dy
    ctx.HFY = ctx.HFY + pmid * dx

    dx = botx - ctx.XOF
    dy = boty - ctx.YOF
    xmid = 0.5 * (botx + ctx.XOF) - ctx.XOF
    ymid = 0.5 * (boty + ctx.YOF) - ctx.YOF
    ctx.HMOM = ctx.HMOM + pmid * (xmid * dx + ymid * dy)
    ctx.HFX = ctx.HFX - pmid * dy
    ctx.HFY = ctx.HFY + pmid * dx


def viscal(ctx, bl, niter1):
    eps1 = 1.0e-4

    niter = niter1

    if not ctx.LWAKE:
        xywake(ctx)

    qwcalc(ctx)
    qiset(ctx)

    if ctx.LALFA:
        ctx.CL, ctx.CM, ctx.CDP, ctx.CL_ALF, ctx.CL_MSQ = clcalc(
            ctx.N, ctx.X, ctx.Y, ctx.GAM, ctx.GAM_A, ctx.ALFA, ctx.MINF, ctx.QINF, ctx.XCMREF, ctx.YCMREF
        )

    if not ctx.LIPAN:
        if ctx.LBLINI:
            gamqv(ctx)

        stfind(ctx)
        iblpan(ctx)
        xicalc(ctx)
        iblsys(ctx)

    uicalc(ctx)

    if not ctx.LBLINI:
        for ibl in range(1, ctx.NBL[1] + 1):
            ctx.UEDG[ibl][1] = ctx.UINV[ibl][1]
        for ibl in range(1, ctx.NBL[2] + 1):
            ctx.UEDG[ibl][2] = ctx.UINV[ibl][2]

    if ctx.LVCONV:
        qvfue(ctx)
        if ctx.LVISC:
            cpcalc(ctx.N + ctx.NW, ctx.QVIS, ctx.QINF, ctx.MINF, ctx.CPV)
            cpcalc(ctx.N + ctx.NW, ctx.QINV, ctx.QINF, ctx.MINF, ctx.CPI)
        else:
            cpcalc(ctx.N, ctx.QINV, ctx.QINF, ctx.MINF, ctx.CPI)
        gamqv(ctx)
        ctx.CL, ctx.CM, ctx.CDP, ctx.CL_ALF, ctx.CL_MSQ = clcalc(
            ctx.N, ctx.X, ctx.Y, ctx.GAM, ctx.GAM_A, ctx.ALFA, ctx.MINF, ctx.QINF, ctx.XCMREF, ctx.YCMREF
        )
        cdcalc(ctx)

    if not ctx.LWDIJ or not ctx.LADIJ:
        qdcalc(ctx)

    if niter == 0:
        raise RuntimeError("VISCAL: NITER=0 not supported.")

    print()
    print("Solving BL system ...")
    for iter_ in range(1, niter + 1):
        setbl(ctx, bl)
        blsolv(ctx)
        update(ctx, bl)

        if ctx.LALFA:
            mrcl(ctx, ctx.CL)
            comset(ctx)
        else:
            qiset(ctx)
            uicalc(ctx)

        qvfue(ctx)
        gamqv(ctx)
        stmove(ctx)

        ctx.CL, ctx.CM, ctx.CDP, ctx.CL_ALF, ctx.CL_MSQ = clcalc(
            ctx.N, ctx.X, ctx.Y, ctx.GAM, ctx.GAM_A, ctx.ALFA, ctx.MINF, ctx.QINF, ctx.XCMREF, ctx.YCMREF
        )
        cdcalc(ctx)

        flags = ""

        if ctx.RLX < 1.0:
            print(
                f"\n{iter_:3d}   rms: {ctx.RMSBL:10.4E}   max: {ctx.RMXBL:10.4E}   {ctx.VMXBL} at {ctx.IMXBL:4d}{ctx.ISMXBL:3d}   RLX:{ctx.RLX:6.3f}{flags}"
            )
        if ctx.RLX == 1.0:
            print(
                f"\n{iter_:3d}   rms: {ctx.RMSBL:10.4E}   max: {ctx.RMXBL:10.4E}   {ctx.VMXBL} at {ctx.IMXBL:4d}{ctx.ISMXBL:3d}{flags}"
            )
        cdpdif = ctx.CD - ctx.CDF
        print(
            f"    a ={ctx.ALFA/ctx.DTOR:7.3f}      CL ={ctx.CL:8.4f}\n"
            f"   Cm ={ctx.CM:8.4f}     CD ={ctx.CD:9.5f}   =>   CDf ={ctx.CDF:9.5f}    CDp ={cdpdif:9.5f}"
        )

        if ctx.RMSBL < eps1:
            ctx.LVCONV = True
            ctx.AVISC = ctx.ALFA
            ctx.MVISC = ctx.MINF
            break
    else:
        print("VISCAL:  Convergence failed")

    cpcalc(ctx.N + ctx.NW, ctx.QINV, ctx.QINF, ctx.MINF, ctx.CPI)
    cpcalc(ctx.N + ctx.NW, ctx.QVIS, ctx.QINF, ctx.MINF, ctx.CPV)
    if ctx.LFLAP:
        mhinge(ctx)

    is_ = 1
    hkmax = 0.0
    hkm = 0.0
    psep = 0.0
    patt = 0.0
    for ibl in range(2, ctx.IBLTE[is_] + 1):
        hki = ctx.DSTR[ibl][is_] / ctx.THET[ibl][is_]
        hkmax = max(hki, hkmax)
        if hkm < 4.0 and hki >= 4.0:
            hfrac = (4.0 - hkm) / (hki - hkm)
            pdefm = ctx.UEDG[ibl - 1][is_] ** 2 * ctx.THET[ibl - 1][is_]
            pdefi = ctx.UEDG[ibl][is_] ** 2 * ctx.THET[ibl][is_]
            psep = pdefm * (1.0 - hfrac) + pdefi * hfrac
        if hkm > 4.0 and hki < 4.0:
            hfrac = (4.0 - hkm) / (hki - hkm)
            pdefm = ctx.UEDG[ibl - 1][is_] ** 2 * ctx.THET[ibl - 1][is_]
            pdefi = ctx.UEDG[ibl][is_] ** 2 * ctx.THET[ibl][is_]
            patt = pdefm * (1.0 - hfrac) + pdefi * hfrac
        hkm = hki
    delp = patt - psep

    print(f" {ctx.ACRIT[is_]:10.3f}{hkmax:10.4f}{ctx.CD:11.6f}{2.0*psep:11.6f}{2.0*patt:11.6f}{2.0*delp:11.6f}{ctx.XOCTR[is_]:10.4f}     #")

    fnum = ctx.XSTRIP[is_] * 100.0
    iten = int(fnum / 9.99999)
    ione = int((fnum - float(10 * iten)) / 0.99999)
    idec = int((fnum - float(10 * iten) - float(ione)) / 0.09999)
    fname = f"{iten:d}{ione:d}{idec:d}.bl"

    with open(fname, "w", encoding="ascii") as lu:
        lu.write("#       s         ue          H          P         K         x    -m du/dx\n")
        for ibl in range(2, ctx.IBLTE[is_] + 1):
            iblm = max(ibl - 1, 2)
            iblp = min(ibl + 1, ctx.IBLTE[is_])
            i = ctx.IPAN[ibl][is_]
            hk = ctx.DSTR[ibl][is_] / ctx.THET[ibl][is_]
            ddef = ctx.DSTR[ibl][is_] * ctx.UEDG[ibl][is_]
            pdef = ctx.THET[ibl][is_] * ctx.UEDG[ibl][is_] ** 2
            edef = ctx.TSTR[ibl][is_] * ctx.UEDG[ibl][is_] ** 3 * 0.5
            duds = (ctx.UEDG[iblp][is_] - ctx.UEDG[iblm][is_]) / (ctx.XSSI[iblp][is_] - ctx.XSSI[iblm][is_])
            dpds = -ddef * duds
            lu.write(
                f" {ctx.XSSI[ibl][is_]:11.4f}{ctx.UEDG[ibl][is_]:11.4f}{hk:11.4f}{pdef:11.6f}{edef:11.6f}{ctx.X[i]:11.3f}{dpds:14.6e}\n"
            )
