import math
from .spline import sinvrt, seval
from .xfoil import cpcalc, clcalc, cdcalc, comset, mrcl
from .xgdes import getxyf
from .xpanel import gamqv, iblpan, qdcalc, qiset, qvfue, qwcalc, stfind, stmove, uicalc, xicalc, xywake
from .xsolve import blsolv
from .xbl import iblsys, setbl, update


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
