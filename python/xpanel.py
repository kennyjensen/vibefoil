# Ported from XFOIL Fortran source (Mark Drela).
# This file is a derived work and remains under the terms of the
# GNU General Public License v2 or later.
# See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

import math

from .xsolve import ludcmp, baksub
from .xutils import atanc, setexp


def apcalc(ctx):
    for i in range(1, ctx.N):
        sx = ctx.X[i + 1] - ctx.X[i]
        sy = ctx.Y[i + 1] - ctx.Y[i]
        if sx == 0.0 and sy == 0.0:
            ctx.APANEL[i] = math.atan2(-ctx.NY[i], -ctx.NX[i])
        else:
            ctx.APANEL[i] = math.atan2(sx, -sy)

    i = ctx.N
    ip = 1
    if ctx.SHARP:
        ctx.APANEL[i] = ctx.PI
    else:
        sx = ctx.X[ip] - ctx.X[i]
        sy = ctx.Y[ip] - ctx.Y[i]
        ctx.APANEL[i] = math.atan2(-sx, sy) + ctx.PI


def ncalc(x, y, s, n, xn, yn):
    from .spline import segspl

    if n <= 1:
        return

    segspl(x, xn, s, n)
    segspl(y, yn, s, n)
    for i in range(1, n + 1):
        sx = yn[i]
        sy = -xn[i]
        smod = math.sqrt(sx * sx + sy * sy)
        if smod == 0.0:
            xn[i] = -1.0
            yn[i] = 0.0
        else:
            xn[i] = sx / smod
            yn[i] = sy / smod

    for i in range(1, n):
        if s[i] == s[i + 1]:
            sx = 0.5 * (xn[i] + xn[i + 1])
            sy = 0.5 * (yn[i] + yn[i + 1])
            smod = math.sqrt(sx * sx + sy * sy)
            if smod == 0.0:
                xn[i] = -1.0
                yn[i] = 0.0
                xn[i + 1] = -1.0
                yn[i + 1] = 0.0
            else:
                xn[i] = sx / smod
                yn[i] = sy / smod
                xn[i + 1] = sx / smod
                yn[i + 1] = sy / smod


def psilin(ctx, i, xi, yi, nxi, nyi, geolin, siglin):
    seps = (ctx.S[ctx.N] - ctx.S[1]) * 1.0e-5
    io = i

    cosa = math.cos(ctx.ALFA)
    sina = math.sin(ctx.ALFA)

    for jo in range(1, ctx.N + 1):
        ctx.DZDG[jo] = 0.0
        ctx.DZDN[jo] = 0.0
        ctx.DQDG[jo] = 0.0

    for jo in range(1, ctx.N + 1):
        ctx.DZDM[jo] = 0.0
        ctx.DQDM[jo] = 0.0

    ctx.Z_QINF = 0.0
    ctx.Z_ALFA = 0.0
    ctx.Z_QDOF0 = 0.0
    ctx.Z_QDOF1 = 0.0
    ctx.Z_QDOF2 = 0.0
    ctx.Z_QDOF3 = 0.0

    psi = 0.0
    psi_ni = 0.0

    qt1 = 0.0
    qt2 = 0.0
    qtanm = 0.0

    if ctx.SHARP:
        scs = 1.0
        sds = 0.0
    else:
        scs = ctx.ANTE / ctx.DSTE
        sds = ctx.ASTE / ctx.DSTE

    for jo in range(1, ctx.N + 1):
        jp = jo + 1
        jm = jo - 1
        jq = jp + 1

        if jo == 1:
            jm = jo
        elif jo == ctx.N - 1:
            jq = jp
        elif jo == ctx.N:
            jp = 1
            if (ctx.X[jo] - ctx.X[jp]) ** 2 + (ctx.Y[jo] - ctx.Y[jp]) ** 2 < seps**2:
                continue

        dso = math.sqrt((ctx.X[jo] - ctx.X[jp]) ** 2 + (ctx.Y[jo] - ctx.Y[jp]) ** 2)
        if dso == 0.0:
            continue

        dsio = 1.0 / dso
        apan = ctx.APANEL[jo]

        rx1 = xi - ctx.X[jo]
        ry1 = yi - ctx.Y[jo]
        rx2 = xi - ctx.X[jp]
        ry2 = yi - ctx.Y[jp]

        sx = (ctx.X[jp] - ctx.X[jo]) * dsio
        sy = (ctx.Y[jp] - ctx.Y[jo]) * dsio

        x1 = sx * rx1 + sy * ry1
        x2 = sx * rx2 + sy * ry2
        yy = sx * ry1 - sy * rx1

        rs1 = rx1 * rx1 + ry1 * ry1
        rs2 = rx2 * rx2 + ry2 * ry2

        if io >= 1 and io <= ctx.N:
            sgn = 1.0
        else:
            sgn = math.copysign(1.0, yy)

        if io != jo and rs1 > 0.0:
            g1 = math.log(rs1)
            t1 = math.atan2(sgn * x1, sgn * yy) + (0.5 - 0.5 * sgn) * ctx.PI
        else:
            g1 = 0.0
            t1 = 0.0

        if io != jp and rs2 > 0.0:
            g2 = math.log(rs2)
            t2 = math.atan2(sgn * x2, sgn * yy) + (0.5 - 0.5 * sgn) * ctx.PI
        else:
            g2 = 0.0
            t2 = 0.0

        x1i = sx * nxi + sy * nyi
        x2i = sx * nxi + sy * nyi
        yyi = sx * nyi - sy * nxi

        if geolin:
            nxo = ctx.NX[jo]
            nyo = ctx.NY[jo]
            nxp = ctx.NX[jp]
            nyp = ctx.NY[jp]

            x1o = -((rx1 - x1 * sx) * nxo + (ry1 - x1 * sy) * nyo) * dsio - (sx * nxo + sy * nyo)
            x1p = ((rx1 - x1 * sx) * nxp + (ry1 - x1 * sy) * nyp) * dsio
            x2o = -((rx2 - x2 * sx) * nxo + (ry2 - x2 * sy) * nyo) * dsio
            x2p = ((rx2 - x2 * sx) * nxp + (ry2 - x2 * sy) * nyp) * dsio - (sx * nxp + sy * nyp)
            yyo = ((rx1 + x1 * sy) * nyo - (ry1 - x1 * sx) * nxo) * dsio - (sx * nyo - sy * nxo)
            yyp = -((rx1 - x1 * sy) * nyp - (ry1 + x1 * sx) * nxp) * dsio

        if jo == ctx.N:
            continue

        if siglin:
            x0 = 0.5 * (x1 + x2)
            rs0 = x0 * x0 + yy * yy
            g0 = math.log(rs0)
            t0 = math.atan2(sgn * x0, sgn * yy) + (0.5 - 0.5 * sgn) * ctx.PI

            dxinv = 1.0 / (x1 - x0)
            psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0)
            pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy) * dxinv

            psx1 = -(t1 - apan)
            psx0 = t0 - apan
            psyy = 0.5 * (g1 - g0)

            pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv
            pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv
            pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv

            dsm = math.sqrt((ctx.X[jp] - ctx.X[jm]) ** 2 + (ctx.Y[jp] - ctx.Y[jm]) ** 2)
            dsim = 1.0 / dsm

            ssum = (ctx.SIG[jp] - ctx.SIG[jo]) * dsio + (ctx.SIG[jp] - ctx.SIG[jm]) * dsim
            sdif = (ctx.SIG[jp] - ctx.SIG[jo]) * dsio - (ctx.SIG[jp] - ctx.SIG[jm]) * dsim

            psi = psi + ctx.QOPI * (psum * ssum + pdif * sdif)

            ctx.DZDM[jm] = ctx.DZDM[jm] + ctx.QOPI * (-psum * dsim + pdif * dsim)
            ctx.DZDM[jo] = ctx.DZDM[jo] + ctx.QOPI * (-psum * dsio - pdif * dsio)
            ctx.DZDM[jp] = ctx.DZDM[jp] + ctx.QOPI * (psum * (dsio + dsim) + pdif * (dsio - dsim))

            psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi
            pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi
            psi_ni = psi_ni + ctx.QOPI * (psni * ssum + pdni * sdif)

            qtanm = qtanm + ctx.QOPI * (psni * ssum + pdni * sdif)

            ctx.DQDM[jm] = ctx.DQDM[jm] + ctx.QOPI * (-psni * dsim + pdni * dsim)
            ctx.DQDM[jo] = ctx.DQDM[jo] + ctx.QOPI * (-psni * dsio - pdni * dsio)
            ctx.DQDM[jp] = ctx.DQDM[jp] + ctx.QOPI * (psni * (dsio + dsim) + pdni * (dsio - dsim))

            dxinv = 1.0 / (x0 - x2)
            psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2)
            pdif = ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy) * dxinv

            psx0 = -(t0 - apan)
            psx2 = t2 - apan
            psyy = 0.5 * (g0 - g2)

            pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv
            pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv
            pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv

            dsp = math.sqrt((ctx.X[jq] - ctx.X[jo]) ** 2 + (ctx.Y[jq] - ctx.Y[jo]) ** 2)
            dsip = 1.0 / dsp

            ssum = (ctx.SIG[jq] - ctx.SIG[jo]) * dsip + (ctx.SIG[jp] - ctx.SIG[jo]) * dsio
            sdif = (ctx.SIG[jq] - ctx.SIG[jo]) * dsip - (ctx.SIG[jp] - ctx.SIG[jo]) * dsio

            psi = psi + ctx.QOPI * (psum * ssum + pdif * sdif)

            ctx.DZDM[jo] = ctx.DZDM[jo] + ctx.QOPI * (-psum * (dsip + dsio) - pdif * (dsip - dsio))
            ctx.DZDM[jp] = ctx.DZDM[jp] + ctx.QOPI * (psum * dsio - pdif * dsio)
            ctx.DZDM[jq] = ctx.DZDM[jq] + ctx.QOPI * (psum * dsip + pdif * dsip)

            psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi
            pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi
            psi_ni = psi_ni + ctx.QOPI * (psni * ssum + pdni * sdif)

            qtanm = qtanm + ctx.QOPI * (psni * ssum + pdni * sdif)

            ctx.DQDM[jo] = ctx.DQDM[jo] + ctx.QOPI * (-psni * (dsip + dsio) - pdni * (dsip - dsio))
            ctx.DQDM[jp] = ctx.DQDM[jp] + ctx.QOPI * (psni * dsio - pdni * dsio)
            ctx.DQDM[jq] = ctx.DQDM[jq] + ctx.QOPI * (psni * dsip + pdni * dsip)

        dxinv = 1.0 / (x1 - x2)
        psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2)
        psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv

        psx1 = 0.5 * g1
        psx2 = -0.5 * g2
        psyy = t1 - t2

        pdx1 = ((x1 + x2) * psx1 + psis - x1 * g1 - psid) * dxinv
        pdx2 = ((x1 + x2) * psx2 + psis + x2 * g2 + psid) * dxinv
        pdyy = ((x1 + x2) * psyy - yy * (g1 - g2)) * dxinv

        gsum1 = ctx.GAMU[jp][1] + ctx.GAMU[jo][1]
        gsum2 = ctx.GAMU[jp][2] + ctx.GAMU[jo][2]
        gdif1 = ctx.GAMU[jp][1] - ctx.GAMU[jo][1]
        gdif2 = ctx.GAMU[jp][2] - ctx.GAMU[jo][2]

        gsum = ctx.GAM[jp] + ctx.GAM[jo]
        gdif = ctx.GAM[jp] - ctx.GAM[jo]

        psi = psi + ctx.QOPI * (psis * gsum + psid * gdif)

        ctx.DZDG[jo] = ctx.DZDG[jo] + ctx.QOPI * (psis - psid)
        ctx.DZDG[jp] = ctx.DZDG[jp] + ctx.QOPI * (psis + psid)

        psni = psx1 * x1i + psx2 * x2i + psyy * yyi
        pdni = pdx1 * x1i + pdx2 * x2i + pdyy * yyi
        psi_ni = psi_ni + ctx.QOPI * (gsum * psni + gdif * pdni)

        qt1 = qt1 + ctx.QOPI * (gsum1 * psni + gdif1 * pdni)
        qt2 = qt2 + ctx.QOPI * (gsum2 * psni + gdif2 * pdni)

        ctx.DQDG[jo] = ctx.DQDG[jo] + ctx.QOPI * (psni - pdni)
        ctx.DQDG[jp] = ctx.DQDG[jp] + ctx.QOPI * (psni + pdni)

        if geolin:
            ctx.DZDN[jo] = ctx.DZDN[jo] + ctx.QOPI * gsum * (psx1 * x1o + psx2 * x2o + psyy * yyo) + ctx.QOPI * gdif * (
                pdx1 * x1o + pdx2 * x2o + pdyy * yyo
            )
            ctx.DZDN[jp] = ctx.DZDN[jp] + ctx.QOPI * gsum * (psx1 * x1p + psx2 * x2p + psyy * yyp) + ctx.QOPI * gdif * (
                pdx1 * x1p + pdx2 * x2p + pdyy * yyp
            )
            ctx.Z_QDOF0 = ctx.Z_QDOF0 + ctx.QOPI * ((psis - psid) * ctx.QF0[jo] + (psis + psid) * ctx.QF0[jp])
            ctx.Z_QDOF1 = ctx.Z_QDOF1 + ctx.QOPI * ((psis - psid) * ctx.QF1[jo] + (psis + psid) * ctx.QF1[jp])
            ctx.Z_QDOF2 = ctx.Z_QDOF2 + ctx.QOPI * ((psis - psid) * ctx.QF2[jo] + (psis + psid) * ctx.QF2[jp])
            ctx.Z_QDOF3 = ctx.Z_QDOF3 + ctx.QOPI * ((psis - psid) * ctx.QF3[jo] + (psis + psid) * ctx.QF3[jp])

    psig = 0.5 * yy * (g1 - g2) + x2 * (t2 - apan) - x1 * (t1 - apan)
    pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2)

    psigx1 = -(t1 - apan)
    psigx2 = t2 - apan
    psigyy = 0.5 * (g1 - g2)
    pgamx1 = 0.5 * g1
    pgamx2 = -0.5 * g2
    pgamyy = t1 - t2

    psigni = psigx1 * x1i + psigx2 * x2i + psigyy * yyi
    pgamni = pgamx1 * x1i + pgamx2 * x2i + pgamyy * yyi

    sigte1 = 0.5 * scs * (ctx.GAMU[jp][1] - ctx.GAMU[jo][1])
    sigte2 = 0.5 * scs * (ctx.GAMU[jp][2] - ctx.GAMU[jo][2])
    gamte1 = -0.5 * sds * (ctx.GAMU[jp][1] - ctx.GAMU[jo][1])
    gamte2 = -0.5 * sds * (ctx.GAMU[jp][2] - ctx.GAMU[jo][2])

    ctx.SIGTE = 0.5 * scs * (ctx.GAM[jp] - ctx.GAM[jo])
    ctx.GAMTE = -0.5 * sds * (ctx.GAM[jp] - ctx.GAM[jo])

    psi = psi + ctx.HOPI * (psig * ctx.SIGTE + pgam * ctx.GAMTE)

    ctx.DZDG[jo] = ctx.DZDG[jo] - ctx.HOPI * psig * scs * 0.5
    ctx.DZDG[jp] = ctx.DZDG[jp] + ctx.HOPI * psig * scs * 0.5

    ctx.DZDG[jo] = ctx.DZDG[jo] + ctx.HOPI * pgam * sds * 0.5
    ctx.DZDG[jp] = ctx.DZDG[jp] - ctx.HOPI * pgam * sds * 0.5

    psi_ni = psi_ni + ctx.HOPI * (psigni * ctx.SIGTE + pgamni * ctx.GAMTE)

    qt1 = qt1 + ctx.HOPI * (psigni * sigte1 + pgamni * gamte1)
    qt2 = qt2 + ctx.HOPI * (psigni * sigte2 + pgamni * gamte2)

    ctx.DQDG[jo] = ctx.DQDG[jo] - ctx.HOPI * (psigni * 0.5 * scs - pgamni * 0.5 * sds)
    ctx.DQDG[jp] = ctx.DQDG[jp] + ctx.HOPI * (psigni * 0.5 * scs - pgamni * 0.5 * sds)

    if geolin:
        ctx.DZDN[jo] = ctx.DZDN[jo] + ctx.HOPI * (psigx1 * x1o + psigx2 * x2o + psigyy * yyo) * ctx.SIGTE + ctx.HOPI * (
            pgamx1 * x1o + pgamx2 * x2o + pgamyy * yyo
        ) * ctx.GAMTE
        ctx.DZDN[jp] = ctx.DZDN[jp] + ctx.HOPI * (psigx1 * x1p + psigx2 * x2p + psigyy * yyp) * ctx.SIGTE + ctx.HOPI * (
            pgamx1 * x1p + pgamx2 * x2p + pgamyy * yyp
        ) * ctx.GAMTE

        ctx.Z_QDOF0 = ctx.Z_QDOF0 + ctx.HOPI * psig * 0.5 * (ctx.QF0[jp] - ctx.QF0[jo]) * scs - ctx.HOPI * pgam * 0.5 * (
            ctx.QF0[jp] - ctx.QF0[jo]
        ) * sds
        ctx.Z_QDOF1 = ctx.Z_QDOF1 + ctx.HOPI * psig * 0.5 * (ctx.QF1[jp] - ctx.QF1[jo]) * scs - ctx.HOPI * pgam * 0.5 * (
            ctx.QF1[jp] - ctx.QF1[jo]
        ) * sds
        ctx.Z_QDOF2 = ctx.Z_QDOF2 + ctx.HOPI * psig * 0.5 * (ctx.QF2[jp] - ctx.QF2[jo]) * scs - ctx.HOPI * pgam * 0.5 * (
            ctx.QF2[jp] - ctx.QF2[jo]
        ) * sds
        ctx.Z_QDOF3 = ctx.Z_QDOF3 + ctx.HOPI * psig * 0.5 * (ctx.QF3[jp] - ctx.QF3[jo]) * scs - ctx.HOPI * pgam * 0.5 * (
            ctx.QF3[jp] - ctx.QF3[jo]
        ) * sds

    psi = psi + ctx.QINF * (cosa * yi - sina * xi)
    psi_ni = psi_ni + ctx.QINF * (cosa * nyi - sina * nxi)

    qt1 = qt1 + ctx.QINF * nyi
    qt2 = qt2 - ctx.QINF * nxi

    ctx.Z_QINF = ctx.Z_QINF + (cosa * yi - sina * xi)
    ctx.Z_ALFA = ctx.Z_ALFA - ctx.QINF * (sina * yi + cosa * xi)

    if ctx.LIMAGE:
        for jo in range(1, ctx.N + 1):
            jp = jo + 1
            jm = jo - 1
            jq = jp + 1

            if jo == 1:
                jm = jo
            elif jo == ctx.N - 1:
                jq = jp
            elif jo == ctx.N:
                jp = 1
                if (ctx.X[jo] - ctx.X[jp]) ** 2 + (ctx.Y[jo] - ctx.Y[jp]) ** 2 < seps**2:
                    continue

            dso = math.sqrt((ctx.X[jo] - ctx.X[jp]) ** 2 + (ctx.Y[jo] - ctx.Y[jp]) ** 2)
            if dso == 0.0:
                continue

            dsio = 1.0 / dso
            apan = ctx.PI - ctx.APANEL[jo] + 2.0 * ctx.ALFA

            xjo = ctx.X[jo] + 2.0 * (ctx.YIMAGE + ctx.Y[jo]) * sina
            yjo = ctx.Y[jo] - 2.0 * (ctx.YIMAGE + ctx.Y[jo]) * cosa
            xjp = ctx.X[jp] + 2.0 * (ctx.YIMAGE + ctx.Y[jp]) * sina
            yjp = ctx.Y[jp] - 2.0 * (ctx.YIMAGE + ctx.Y[jp]) * cosa

            rx1 = xi - xjo
            ry1 = yi - yjo
            rx2 = xi - xjp
            ry2 = yi - yjp

            sx = (xjp - xjo) * dsio
            sy = (yjp - yjo) * dsio

            x1 = sx * rx1 + sy * ry1
            x2 = sx * rx2 + sy * ry2
            yy = sx * ry1 - sy * rx1

            rs1 = rx1 * rx1 + ry1 * ry1
            rs2 = rx2 * rx2 + ry2 * ry2

            if io >= 1 and io <= ctx.N:
                sgn = 1.0
            else:
                sgn = math.copysign(1.0, yy)

            if io != jo and rs1 > 0.0:
                g1 = math.log(rs1)
                t1 = math.atan2(sgn * x1, sgn * yy) + (0.5 - 0.5 * sgn) * ctx.PI
            else:
                g1 = 0.0
                t1 = 0.0

            if io != jp and rs2 > 0.0:
                g2 = math.log(rs2)
                t2 = math.atan2(sgn * x2, sgn * yy) + (0.5 - 0.5 * sgn) * ctx.PI
            else:
                g2 = 0.0
                t2 = 0.0

            x1i = sx * nxi + sy * nyi
            x2i = sx * nxi + sy * nyi
            yyi = sx * nyi - sy * nxi

            if geolin:
                nxo = ctx.NX[jo]
                nyo = ctx.NY[jo]
                nxp = ctx.NX[jp]
                nyp = ctx.NY[jp]

                x1o = -((rx1 - x1 * sx) * nxo + (ry1 - x1 * sy) * nyo) * dsio - (sx * nxo + sy * nyo)
                x1p = ((rx1 - x1 * sx) * nxp + (ry1 - x1 * sy) * nyp) * dsio
                x2o = -((rx2 - x2 * sx) * nxo + (ry2 - x2 * sy) * nyo) * dsio
                x2p = ((rx2 - x2 * sx) * nxp + (ry2 - x2 * sy) * nyp) * dsio - (sx * nxp + sy * nyp)
                yyo = ((rx1 + x1 * sy) * nyo - (ry1 - x1 * sx) * nxo) * dsio - (sx * nyo - sy * nxo)
                yyp = -((rx1 - x1 * sy) * nyp - (ry1 + x1 * sx) * nxp) * dsio

            if jo == ctx.N:
                continue

            if siglin:
                x0 = 0.5 * (x1 + x2)
                rs0 = x0 * x0 + yy * yy
                g0 = math.log(rs0)
                t0 = math.atan2(sgn * x0, sgn * yy) + (0.5 - 0.5 * sgn) * ctx.PI

                dxinv = 1.0 / (x1 - x0)
                psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0)
                pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy) * dxinv

                psx1 = -(t1 - apan)
                psx0 = t0 - apan
                psyy = 0.5 * (g1 - g0)

                pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv
                pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv
                pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv

                dsm = math.sqrt((ctx.X[jp] - ctx.X[jm]) ** 2 + (ctx.Y[jp] - ctx.Y[jm]) ** 2)
                dsim = 1.0 / dsm

                ssum = (ctx.SIG[jp] - ctx.SIG[jo]) * dsio + (ctx.SIG[jp] - ctx.SIG[jm]) * dsim
                sdif = (ctx.SIG[jp] - ctx.SIG[jo]) * dsio - (ctx.SIG[jp] - ctx.SIG[jm]) * dsim

                psi = psi + ctx.QOPI * (psum * ssum + pdif * sdif)

                ctx.DZDM[jm] = ctx.DZDM[jm] + ctx.QOPI * (-psum * dsim + pdif * dsim)
                ctx.DZDM[jo] = ctx.DZDM[jo] + ctx.QOPI * (-psum * dsio - pdif * dsio)
                ctx.DZDM[jp] = ctx.DZDM[jp] + ctx.QOPI * (psum * (dsio + dsim) + pdif * (dsio - dsim))

                psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi
                pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi
                psi_ni = psi_ni + ctx.QOPI * (psni * ssum + pdni * sdif)

                ctx.DQDM[jm] = ctx.DQDM[jm] + ctx.QOPI * (-psni * dsim + pdni * dsim)
                ctx.DQDM[jo] = ctx.DQDM[jo] + ctx.QOPI * (-psni * dsio - pdni * dsio)
                ctx.DQDM[jp] = ctx.DQDM[jp] + ctx.QOPI * (psni * (dsio + dsim) + pdni * (dsio - dsim))

                dxinv = 1.0 / (x0 - x2)
                psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2)
                pdif = ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy) * dxinv

                psx0 = -(t0 - apan)
                psx2 = t2 - apan
                psyy = 0.5 * (g0 - g2)

                pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv
                pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv
                pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv

                dsp = math.sqrt((ctx.X[jq] - ctx.X[jo]) ** 2 + (ctx.Y[jq] - ctx.Y[jo]) ** 2)
                dsip = 1.0 / dsp

                ssum = (ctx.SIG[jq] - ctx.SIG[jo]) * dsip + (ctx.SIG[jp] - ctx.SIG[jo]) * dsio
                sdif = (ctx.SIG[jq] - ctx.SIG[jo]) * dsip - (ctx.SIG[jp] - ctx.SIG[jo]) * dsio

                psi = psi + ctx.QOPI * (psum * ssum + pdif * sdif)

                ctx.DZDM[jo] = ctx.DZDM[jo] + ctx.QOPI * (-psum * (dsip + dsio) - pdif * (dsip - dsio))
                ctx.DZDM[jp] = ctx.DZDM[jp] + ctx.QOPI * (psum * dsio - pdif * dsio)
                ctx.DZDM[jq] = ctx.DZDM[jq] + ctx.QOPI * (psum * dsip + pdif * dsip)

                psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi
                pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi
                psi_ni = psi_ni + ctx.QOPI * (psni * ssum + pdni * sdif)

                ctx.DQDM[jo] = ctx.DQDM[jo] + ctx.QOPI * (-psni * (dsip + dsio) - pdni * (dsip - dsio))
                ctx.DQDM[jp] = ctx.DQDM[jp] + ctx.QOPI * (psni * dsio - pdni * dsio)
                ctx.DQDM[jq] = ctx.DQDM[jq] + ctx.QOPI * (psni * dsip + pdni * dsip)

            dxinv = 1.0 / (x1 - x2)
            psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2)
            psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv

            psx1 = 0.5 * g1
            psx2 = -0.5 * g2
            psyy = t1 - t2

            pdx1 = ((x1 + x2) * psx1 + psis - x1 * g1 - psid) * dxinv
            pdx2 = ((x1 + x2) * psx2 + psis + x2 * g2 + psid) * dxinv
            pdyy = ((x1 + x2) * psyy - yy * (g1 - g2)) * dxinv

            gsum1 = ctx.GAMU[jp][1] + ctx.GAMU[jo][1]
            gsum2 = ctx.GAMU[jp][2] + ctx.GAMU[jo][2]
            gdif1 = ctx.GAMU[jp][1] - ctx.GAMU[jo][1]
            gdif2 = ctx.GAMU[jp][2] - ctx.GAMU[jo][2]

            gsum = ctx.GAM[jp] + ctx.GAM[jo]
            gdif = ctx.GAM[jp] - ctx.GAM[jo]

            psi = psi - ctx.QOPI * (psis * gsum + psid * gdif)

            ctx.DZDG[jo] = ctx.DZDG[jo] - ctx.QOPI * (psis - psid)
            ctx.DZDG[jp] = ctx.DZDG[jp] - ctx.QOPI * (psis + psid)

            psni = psx1 * x1i + psx2 * x2i + psyy * yyi
            pdni = pdx1 * x1i + pdx2 * x2i + pdyy * yyi
            psi_ni = psi_ni - ctx.QOPI * (gsum * psni + gdif * pdni)

            qt1 = qt1 - ctx.QOPI * (gsum1 * psni + gdif1 * pdni)
            qt2 = qt2 - ctx.QOPI * (gsum2 * psni + gdif2 * pdni)

            ctx.DQDG[jo] = ctx.DQDG[jo] - ctx.QOPI * (psni - pdni)
            ctx.DQDG[jp] = ctx.DQDG[jp] - ctx.QOPI * (psni + pdni)

            if geolin:
                ctx.DZDN[jo] = ctx.DZDN[jo] - ctx.QOPI * gsum * (psx1 * x1o + psx2 * x2o + psyy * yyo) - ctx.QOPI * gdif * (
                    pdx1 * x1o + pdx2 * x2o + pdyy * yyo
                )
                ctx.DZDN[jp] = ctx.DZDN[jp] - ctx.QOPI * gsum * (psx1 * x1p + psx2 * x2p + psyy * yyp) - ctx.QOPI * gdif * (
                    pdx1 * x1p + pdx2 * x2p + pdyy * yyp
                )
                ctx.Z_QDOF0 = ctx.Z_QDOF0 - ctx.QOPI * ((psis - psid) * ctx.QF0[jo] + (psis + psid) * ctx.QF0[jp])
                ctx.Z_QDOF1 = ctx.Z_QDOF1 - ctx.QOPI * ((psis - psid) * ctx.QF1[jo] + (psis + psid) * ctx.QF1[jp])
                ctx.Z_QDOF2 = ctx.Z_QDOF2 - ctx.QOPI * ((psis - psid) * ctx.QF2[jo] + (psis + psid) * ctx.QF2[jp])
                ctx.Z_QDOF3 = ctx.Z_QDOF3 - ctx.QOPI * ((psis - psid) * ctx.QF3[jo] + (psis + psid) * ctx.QF3[jp])

        psig = 0.5 * yy * (g1 - g2) + x2 * (t2 - apan) - x1 * (t1 - apan)
        pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2)

        psigx1 = -(t1 - apan)
        psigx2 = t2 - apan
        psigyy = 0.5 * (g1 - g2)
        pgamx1 = 0.5 * g1
        pgamx2 = -0.5 * g2
        pgamyy = t1 - t2

        psigni = psigx1 * x1i + psigx2 * x2i + psigyy * yyi
        pgamni = pgamx1 * x1i + pgamx2 * x2i + pgamyy * yyi

        sigte1 = 0.5 * scs * (ctx.GAMU[jp][1] - ctx.GAMU[jo][1])
        sigte2 = 0.5 * scs * (ctx.GAMU[jp][2] - ctx.GAMU[jo][2])
        gamte1 = -0.5 * sds * (ctx.GAMU[jp][1] - ctx.GAMU[jo][1])
        gamte2 = -0.5 * sds * (ctx.GAMU[jp][2] - ctx.GAMU[jo][2])

        ctx.SIGTE = 0.5 * scs * (ctx.GAM[jp] - ctx.GAM[jo])
        ctx.GAMTE = -0.5 * sds * (ctx.GAM[jp] - ctx.GAM[jo])

        psi = psi + ctx.HOPI * (psig * ctx.SIGTE - pgam * ctx.GAMTE)

        ctx.DZDG[jo] = ctx.DZDG[jo] - ctx.HOPI * psig * scs * 0.5
        ctx.DZDG[jp] = ctx.DZDG[jp] + ctx.HOPI * psig * scs * 0.5

        ctx.DZDG[jo] = ctx.DZDG[jo] - ctx.HOPI * pgam * sds * 0.5
        ctx.DZDG[jp] = ctx.DZDG[jp] + ctx.HOPI * pgam * sds * 0.5

        psi_ni = psi_ni + ctx.HOPI * (psigni * ctx.SIGTE - pgamni * ctx.GAMTE)

        qt1 = qt1 + ctx.HOPI * (psigni * sigte1 - pgamni * gamte1)
        qt2 = qt2 + ctx.HOPI * (psigni * sigte2 - pgamni * gamte2)

        ctx.DQDG[jo] = ctx.DQDG[jo] - ctx.HOPI * (psigni * 0.5 * scs + pgamni * 0.5 * sds)
        ctx.DQDG[jp] = ctx.DQDG[jp] + ctx.HOPI * (psigni * 0.5 * scs + pgamni * 0.5 * sds)

        if geolin:
            ctx.DZDN[jo] = ctx.DZDN[jo] + ctx.HOPI * (psigx1 * x1o + psigx2 * x2o + psigyy * yyo) * ctx.SIGTE - ctx.HOPI * (
                pgamx1 * x1o + pgamx2 * x2o + pgamyy * yyo
            ) * ctx.GAMTE
            ctx.DZDN[jp] = ctx.DZDN[jp] + ctx.HOPI * (psigx1 * x1p + psigx2 * x2p + psigyy * yyp) * ctx.SIGTE - ctx.HOPI * (
                pgamx1 * x1p + pgamx2 * x2p + pgamyy * yyp
            ) * ctx.GAMTE
            ctx.Z_QDOF0 = ctx.Z_QDOF0 + ctx.HOPI * psig * 0.5 * (ctx.QF0[jp] - ctx.QF0[jo]) * scs + ctx.HOPI * pgam * 0.5 * (
                ctx.QF0[jp] - ctx.QF0[jo]
            ) * sds
            ctx.Z_QDOF1 = ctx.Z_QDOF1 + ctx.HOPI * psig * 0.5 * (ctx.QF1[jp] - ctx.QF1[jo]) * scs + ctx.HOPI * pgam * 0.5 * (
                ctx.QF1[jp] - ctx.QF1[jo]
            ) * sds
            ctx.Z_QDOF2 = ctx.Z_QDOF2 + ctx.HOPI * psig * 0.5 * (ctx.QF2[jp] - ctx.QF2[jo]) * scs + ctx.HOPI * pgam * 0.5 * (
                ctx.QF2[jp] - ctx.QF2[jo]
            ) * sds
            ctx.Z_QDOF3 = ctx.Z_QDOF3 + ctx.HOPI * psig * 0.5 * (ctx.QF3[jp] - ctx.QF3[jo]) * scs + ctx.HOPI * pgam * 0.5 * (
                ctx.QF3[jp] - ctx.QF3[jo]
            ) * sds

    return psi, psi_ni, qt1, qt2, qtanm


def pswlin(ctx, i, xi, yi, nxi, nyi):
    io = i
    cosa = math.cos(ctx.ALFA)
    sina = math.sin(ctx.ALFA)

    for jo in range(ctx.N + 1, ctx.N + ctx.NW + 1):
        ctx.DZDM[jo] = 0.0
        ctx.DQDM[jo] = 0.0

    psi = 0.0
    psi_ni = 0.0

    for jo in range(ctx.N + 1, ctx.N + ctx.NW):
        jp = jo + 1

        jm = jo - 1
        jq = jp + 1
        if jo == ctx.N + 1:
            jm = jo
        elif jo == ctx.N + ctx.NW - 1:
            jq = jp

        dso = math.sqrt((ctx.X[jo] - ctx.X[jp]) ** 2 + (ctx.Y[jo] - ctx.Y[jp]) ** 2)
        dsio = 1.0 / dso

        apan = ctx.APANEL[jo]

        rx1 = xi - ctx.X[jo]
        ry1 = yi - ctx.Y[jo]
        rx2 = xi - ctx.X[jp]
        ry2 = yi - ctx.Y[jp]

        sx = (ctx.X[jp] - ctx.X[jo]) * dsio
        sy = (ctx.Y[jp] - ctx.Y[jo]) * dsio

        x1 = sx * rx1 + sy * ry1
        x2 = sx * rx2 + sy * ry2
        yy = sx * ry1 - sy * rx1

        rs1 = rx1 * rx1 + ry1 * ry1
        rs2 = rx2 * rx2 + ry2 * ry2

        if io >= ctx.N + 1 and io <= ctx.N + ctx.NW:
            sgn = 1.0
        else:
            sgn = math.copysign(1.0, yy)

        if io != jo and rs1 > 0.0:
            g1 = math.log(rs1)
            t1 = math.atan2(sgn * x1, sgn * yy) - (0.5 - 0.5 * sgn) * ctx.PI
        else:
            g1 = 0.0
            t1 = 0.0

        if io != jp and rs2 > 0.0:
            g2 = math.log(rs2)
            t2 = math.atan2(sgn * x2, sgn * yy) - (0.5 - 0.5 * sgn) * ctx.PI
        else:
            g2 = 0.0
            t2 = 0.0

        x1i = sx * nxi + sy * nyi
        x2i = sx * nxi + sy * nyi
        yyi = sx * nyi - sy * nxi

        x0 = 0.5 * (x1 + x2)
        rs0 = x0 * x0 + yy * yy
        g0 = math.log(rs0)
        t0 = math.atan2(sgn * x0, sgn * yy) - (0.5 - 0.5 * sgn) * ctx.PI

        dxinv = 1.0 / (x1 - x0)
        psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0)
        pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy) * dxinv

        psx1 = -(t1 - apan)
        psx0 = t0 - apan
        psyy = 0.5 * (g1 - g0)

        pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv
        pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv
        pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv

        dsm = math.sqrt((ctx.X[jp] - ctx.X[jm]) ** 2 + (ctx.Y[jp] - ctx.Y[jm]) ** 2)
        dsim = 1.0 / dsm

        ssum = (ctx.SIG[jp] - ctx.SIG[jo]) * dsio + (ctx.SIG[jp] - ctx.SIG[jm]) * dsim
        sdif = (ctx.SIG[jp] - ctx.SIG[jo]) * dsio - (ctx.SIG[jp] - ctx.SIG[jm]) * dsim

        psi = psi + ctx.QOPI * (psum * ssum + pdif * sdif)

        ctx.DZDM[jm] = ctx.DZDM[jm] + ctx.QOPI * (-psum * dsim + pdif * dsim)
        ctx.DZDM[jo] = ctx.DZDM[jo] + ctx.QOPI * (-psum * dsio - pdif * dsio)
        ctx.DZDM[jp] = ctx.DZDM[jp] + ctx.QOPI * (psum * (dsio + dsim) + pdif * (dsio - dsim))

        psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi
        pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi
        psi_ni = psi_ni + ctx.QOPI * (psni * ssum + pdni * sdif)

        ctx.DQDM[jm] = ctx.DQDM[jm] + ctx.QOPI * (-psni * dsim + pdni * dsim)
        ctx.DQDM[jo] = ctx.DQDM[jo] + ctx.QOPI * (-psni * dsio - pdni * dsio)
        ctx.DQDM[jp] = ctx.DQDM[jp] + ctx.QOPI * (psni * (dsio + dsim) + pdni * (dsio - dsim))

        dxinv = 1.0 / (x0 - x2)
        psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2)
        pdif = ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy) * dxinv

        psx0 = -(t0 - apan)
        psx2 = t2 - apan
        psyy = 0.5 * (g0 - g2)

        pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv
        pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv
        pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv

        dsp = math.sqrt((ctx.X[jq] - ctx.X[jo]) ** 2 + (ctx.Y[jq] - ctx.Y[jo]) ** 2)
        dsip = 1.0 / dsp

        ssum = (ctx.SIG[jq] - ctx.SIG[jo]) * dsip + (ctx.SIG[jp] - ctx.SIG[jo]) * dsio
        sdif = (ctx.SIG[jq] - ctx.SIG[jo]) * dsip - (ctx.SIG[jp] - ctx.SIG[jo]) * dsio

        psi = psi + ctx.QOPI * (psum * ssum + pdif * sdif)

        ctx.DZDM[jo] = ctx.DZDM[jo] + ctx.QOPI * (-psum * (dsip + dsio) - pdif * (dsip - dsio))
        ctx.DZDM[jp] = ctx.DZDM[jp] + ctx.QOPI * (psum * dsio - pdif * dsio)
        ctx.DZDM[jq] = ctx.DZDM[jq] + ctx.QOPI * (psum * dsip + pdif * dsip)

        psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi
        pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi
        psi_ni = psi_ni + ctx.QOPI * (psni * ssum + pdni * sdif)

        ctx.DQDM[jo] = ctx.DQDM[jo] + ctx.QOPI * (-psni * (dsip + dsio) - pdni * (dsip - dsio))
        ctx.DQDM[jp] = ctx.DQDM[jp] + ctx.QOPI * (psni * dsio - pdni * dsio)
        ctx.DQDM[jq] = ctx.DQDM[jq] + ctx.QOPI * (psni * dsip + pdni * dsip)

    psi = psi + ctx.QINF * (cosa * yi - sina * xi)
    psi_ni = psi_ni + ctx.QINF * (cosa * nyi - sina * nxi)

    return psi, psi_ni


def ggcalc(ctx):
    bwt = 0.1

    print("Calculating unit vorticity distributions ...")

    for i in range(1, ctx.N + 1):
        ctx.GAM[i] = 0.0
        ctx.GAMU[i][1] = 0.0
        ctx.GAMU[i][2] = 0.0
    psio = 0.0

    for i in range(1, ctx.N + 1):
        psi, psi_n, _, _, _ = psilin(ctx, i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i], False, True)

        psiinf = ctx.QINF * (math.cos(ctx.ALFA) * ctx.Y[i] - math.sin(ctx.ALFA) * ctx.X[i])

        res1 = ctx.QINF * ctx.Y[i]
        res2 = -ctx.QINF * ctx.X[i]

        for j in range(1, ctx.N + 1):
            ctx.AIJ[i][j] = ctx.DZDG[j]

        for j in range(1, ctx.N + 1):
            ctx.BIJ[i][j] = -ctx.DZDM[j]

        ctx.AIJ[i][ctx.N + 1] = -1.0

        ctx.GAMU[i][1] = -res1
        ctx.GAMU[i][2] = -res2

    res = 0.0

    for j in range(1, ctx.N + 2):
        ctx.AIJ[ctx.N + 1][j] = 0.0

    ctx.AIJ[ctx.N + 1][1] = 1.0
    ctx.AIJ[ctx.N + 1][ctx.N] = 1.0

    ctx.GAMU[ctx.N + 1][1] = -res
    ctx.GAMU[ctx.N + 1][2] = -res

    for j in range(1, ctx.N + 1):
        ctx.BIJ[ctx.N + 1][j] = 0.0

    if ctx.SHARP:
        ag1 = math.atan2(-ctx.YP[1], -ctx.XP[1])
        ag2 = atanc(ctx.YP[ctx.N], ctx.XP[ctx.N], ag1)
        abis = 0.5 * (ag1 + ag2)
        cbis = math.cos(abis)
        sbis = math.sin(abis)

        ds1 = math.sqrt((ctx.X[1] - ctx.X[2]) ** 2 + (ctx.Y[1] - ctx.Y[2]) ** 2)
        ds2 = math.sqrt((ctx.X[ctx.N] - ctx.X[ctx.N - 1]) ** 2 + (ctx.Y[ctx.N] - ctx.Y[ctx.N - 1]) ** 2)
        dsmin = min(ds1, ds2)

        xbis = ctx.XTE - bwt * dsmin * cbis
        ybis = ctx.YTE - bwt * dsmin * sbis

        psi, qbis, _, _, _ = psilin(ctx, 0, xbis, ybis, -sbis, cbis, False, True)

        res = qbis

        for j in range(1, ctx.N + 1):
            ctx.AIJ[ctx.N][j] = ctx.DQDG[j]

        for j in range(1, ctx.N + 1):
            ctx.BIJ[ctx.N][j] = -ctx.DQDM[j]

        ctx.AIJ[ctx.N][ctx.N + 1] = 0.0

        ctx.GAMU[ctx.N][1] = -cbis
        ctx.GAMU[ctx.N][2] = -sbis

    ludcmp(len(ctx.AIJ) - 1, ctx.N + 1, ctx.AIJ, ctx.AIJPIV)
    ctx.LQAIJ = True

    col1 = [0.0] * (ctx.N + 2)
    col2 = [0.0] * (ctx.N + 2)
    for i in range(1, ctx.N + 2):
        col1[i] = ctx.GAMU[i][1]
        col2[i] = ctx.GAMU[i][2]
    baksub(len(ctx.AIJ) - 1, ctx.N + 1, ctx.AIJ, ctx.AIJPIV, col1)
    baksub(len(ctx.AIJ) - 1, ctx.N + 1, ctx.AIJ, ctx.AIJPIV, col2)
    for i in range(1, ctx.N + 2):
        ctx.GAMU[i][1] = col1[i]
        ctx.GAMU[i][2] = col2[i]

    for i in range(1, ctx.N + 1):
        ctx.QINVU[i][1] = ctx.GAMU[i][1]
        ctx.QINVU[i][2] = ctx.GAMU[i][2]

    ctx.LGAMU = True


def qwcalc(ctx):
    ctx.QINVU[ctx.N + 1][1] = ctx.QINVU[ctx.N][1]
    ctx.QINVU[ctx.N + 1][2] = ctx.QINVU[ctx.N][2]

    for i in range(ctx.N + 2, ctx.N + ctx.NW + 1):
        psi, psi_ni, qt1, qt2, _ = psilin(ctx, i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i], False, False)
        ctx.QINVU[i][1] = qt1
        ctx.QINVU[i][2] = qt2


def qdcalc(ctx):
    print("Calculating source influence matrix ...")

    if not ctx.LADIJ:
        for j in range(1, ctx.N + 1):
            col = [0.0] * (ctx.N + 2)
            for i in range(1, ctx.N + 2):
                col[i] = ctx.BIJ[i][j]
            baksub(len(ctx.AIJ) - 1, ctx.N + 1, ctx.AIJ, ctx.AIJPIV, col)
            for i in range(1, ctx.N + 2):
                ctx.BIJ[i][j] = col[i]
            for i in range(1, ctx.N + 1):
                ctx.DIJ[i][j] = ctx.BIJ[i][j]
        ctx.LADIJ = True

    for i in range(1, ctx.N + 1):
        psi, psi_n = pswlin(ctx, i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i])
        for j in range(ctx.N + 1, ctx.N + ctx.NW + 1):
            ctx.BIJ[i][j] = -ctx.DZDM[j]

    for j in range(ctx.N + 1, ctx.N + ctx.NW + 1):
        ctx.BIJ[ctx.N + 1][j] = 0.0

    if ctx.SHARP:
        bwt = 0.1
        ag1 = math.atan2(-ctx.YP[1], -ctx.XP[1])
        ag2 = atanc(ctx.YP[ctx.N], ctx.XP[ctx.N], ag1)
        abis = 0.5 * (ag1 + ag2)
        cbis = math.cos(abis)
        sbis = math.sin(abis)
        ds1 = math.sqrt((ctx.X[1] - ctx.X[2]) ** 2 + (ctx.Y[1] - ctx.Y[2]) ** 2)
        ds2 = math.sqrt((ctx.X[ctx.N] - ctx.X[ctx.N - 1]) ** 2 + (ctx.Y[ctx.N] - ctx.Y[ctx.N - 1]) ** 2)
        dsmin = min(ds1, ds2)
        xbis = ctx.XTE - bwt * dsmin * cbis
        ybis = ctx.YTE - bwt * dsmin * sbis
        pswlin(ctx, 0, xbis, ybis, -sbis, cbis)
        for j in range(ctx.N + 1, ctx.N + ctx.NW + 1):
            ctx.BIJ[ctx.N][j] = -ctx.DQDM[j]

    for j in range(ctx.N + 1, ctx.N + ctx.NW + 1):
        col = [0.0] * (ctx.N + 2)
        for i in range(1, ctx.N + 2):
            col[i] = ctx.BIJ[i][j]
        baksub(len(ctx.AIJ) - 1, ctx.N + 1, ctx.AIJ, ctx.AIJPIV, col)
        for i in range(1, ctx.N + 2):
            ctx.BIJ[i][j] = col[i]

    for i in range(1, ctx.N + 1):
        for j in range(ctx.N + 1, ctx.N + ctx.NW + 1):
            ctx.DIJ[i][j] = ctx.BIJ[i][j]

    for i in range(ctx.N + 1, ctx.N + ctx.NW + 1):
        iw = i - ctx.N
        psi, psi_n, _, _, _ = psilin(ctx, i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i], False, True)
        for j in range(1, ctx.N + 1):
            ctx.CIJ[iw][j] = ctx.DQDG[j]
        for j in range(1, ctx.N + 1):
            ctx.DIJ[i][j] = ctx.DQDM[j]

        psi, psi_n = pswlin(ctx, i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i])
        for j in range(ctx.N + 1, ctx.N + ctx.NW + 1):
            ctx.DIJ[i][j] = ctx.DQDM[j]

    for i in range(ctx.N + 1, ctx.N + ctx.NW + 1):
        iw = i - ctx.N
        for j in range(1, ctx.N + 1):
            summ = 0.0
            for k in range(1, ctx.N + 1):
                summ = summ + ctx.CIJ[iw][k] * ctx.DIJ[k][j]
            ctx.DIJ[i][j] = ctx.DIJ[i][j] + summ

        for j in range(ctx.N + 1, ctx.N + ctx.NW + 1):
            summ = 0.0
            for k in range(1, ctx.N + 1):
                summ = summ + ctx.CIJ[iw][k] * ctx.BIJ[k][j]
            ctx.DIJ[i][j] = ctx.DIJ[i][j] + summ

    for j in range(1, ctx.N + ctx.NW + 1):
        ctx.DIJ[ctx.N + 1][j] = ctx.DIJ[ctx.N][j]

    ctx.LWDIJ = True


def xywake(ctx):
    print("Calculating wake trajectory ...")

    ctx.NW = ctx.N // 12 + 10 * int(ctx.WAKLEN)
    if ctx.NW > len(ctx.WGAP) - 1:
        print()
        print("Array size (IWX) too small.  Last wake point index reduced.")
        ctx.NW = len(ctx.WGAP) - 1

    ds1 = 0.5 * (ctx.S[2] - ctx.S[1] + ctx.S[ctx.N] - ctx.S[ctx.N - 1])
    s_tmp = [0.0] * (ctx.NW + 1)
    setexp(s_tmp, ds1, ctx.WAKLEN * ctx.CHORD, ctx.NW)
    for i in range(1, ctx.NW + 1):
        ctx.SNEW[ctx.N + i] = s_tmp[i]

    ctx.XTE = 0.5 * (ctx.X[1] + ctx.X[ctx.N])
    ctx.YTE = 0.5 * (ctx.Y[1] + ctx.Y[ctx.N])

    i = ctx.N + 1
    sx = 0.5 * (ctx.YP[ctx.N] - ctx.YP[1])
    sy = 0.5 * (ctx.XP[1] - ctx.XP[ctx.N])
    smod = math.sqrt(sx**2 + sy**2)
    ctx.NX[i] = sx / smod
    ctx.NY[i] = sy / smod
    ctx.X[i] = ctx.XTE - 0.0001 * ctx.NY[i]
    ctx.Y[i] = ctx.YTE + 0.0001 * ctx.NX[i]
    ctx.S[i] = ctx.S[ctx.N]

    psi, psi_x, _, _, _ = psilin(ctx, i, ctx.X[i], ctx.Y[i], 1.0, 0.0, False, False)
    psi, psi_y, _, _, _ = psilin(ctx, i, ctx.X[i], ctx.Y[i], 0.0, 1.0, False, False)

    ctx.NX[i + 1] = -psi_x / math.sqrt(psi_x**2 + psi_y**2)
    ctx.NY[i + 1] = -psi_y / math.sqrt(psi_x**2 + psi_y**2)

    ctx.APANEL[i] = math.atan2(psi_y, psi_x)

    for i in range(ctx.N + 2, ctx.N + ctx.NW + 1):
        ds = ctx.SNEW[i] - ctx.SNEW[i - 1]
        ctx.X[i] = ctx.X[i - 1] - ds * ctx.NY[i]
        ctx.Y[i] = ctx.Y[i - 1] + ds * ctx.NX[i]
        ctx.S[i] = ctx.S[i - 1] + ds

        if i == ctx.N + ctx.NW:
            continue

        psi, psi_x, _, _, _ = psilin(ctx, i, ctx.X[i], ctx.Y[i], 1.0, 0.0, False, False)
        psi, psi_y, _, _, _ = psilin(ctx, i, ctx.X[i], ctx.Y[i], 0.0, 1.0, False, False)

        ctx.NX[i + 1] = -psi_x / math.sqrt(psi_x**2 + psi_y**2)
        ctx.NY[i + 1] = -psi_y / math.sqrt(psi_x**2 + psi_y**2)

        ctx.APANEL[i] = math.atan2(psi_y, psi_x)

    ctx.LWAKE = True
    ctx.AWAKE = ctx.ALFA
    ctx.LWDIJ = False


def stfind(ctx):
    for i in range(1, ctx.N):
        if ctx.GAM[i] >= 0.0 and ctx.GAM[i + 1] < 0.0:
            break
    else:
        print("STFIND: Stagnation point not found. Continuing ...")
        i = ctx.N // 2

    ctx.IST = i
    dgam = ctx.GAM[i + 1] - ctx.GAM[i]
    ds = ctx.S[i + 1] - ctx.S[i]

    if ctx.GAM[i] < -ctx.GAM[i + 1]:
        ctx.SST = ctx.S[i] - ds * (ctx.GAM[i] / dgam)
    else:
        ctx.SST = ctx.S[i + 1] - ds * (ctx.GAM[i + 1] / dgam)

    if ctx.SST <= ctx.S[i]:
        ctx.SST = ctx.S[i] + 1.0e-7
    if ctx.SST >= ctx.S[i + 1]:
        ctx.SST = ctx.S[i + 1] - 1.0e-7

    ctx.SST_GO = (ctx.SST - ctx.S[i + 1]) / dgam
    ctx.SST_GP = (ctx.S[i] - ctx.SST) / dgam


def iblpan(ctx):
    is_ = 1
    ibl = 1
    for i in range(ctx.IST, 0, -1):
        ibl += 1
        ctx.IPAN[ibl][is_] = i
        ctx.VTI[ibl][is_] = 1.0

    ctx.IBLTE[is_] = ibl
    ctx.NBL[is_] = ibl

    is_ = 2
    ibl = 1
    for i in range(ctx.IST + 1, ctx.N + 1):
        ibl += 1
        ctx.IPAN[ibl][is_] = i
        ctx.VTI[ibl][is_] = -1.0

    ctx.IBLTE[is_] = ibl

    for iw in range(1, ctx.NW + 1):
        i = ctx.N + iw
        ibl = ctx.IBLTE[is_] + iw
        ctx.IPAN[ibl][is_] = i
        ctx.VTI[ibl][is_] = -1.0

    ctx.NBL[is_] = ctx.IBLTE[is_] + ctx.NW

    for iw in range(1, ctx.NW + 1):
        ctx.IPAN[ctx.IBLTE[1] + iw][1] = ctx.IPAN[ctx.IBLTE[2] + iw][2]
        ctx.VTI[ctx.IBLTE[1] + iw][1] = 1.0

    iblmax = max(ctx.IBLTE[1], ctx.IBLTE[2]) + ctx.NW
    if iblmax > len(ctx.XSSI) - 1:
        raise RuntimeError(f" ***  BL array overflow.  Increase IVX to at least {iblmax}")

    ctx.LIPAN = True


def xicalc(ctx):
    xfeps = 1.0e-7
    xeps = xfeps * (ctx.S[ctx.N] - ctx.S[1])

    is_ = 1
    ctx.XSSI[1][is_] = 0.0
    for ibl in range(2, ctx.IBLTE[is_] + 1):
        i = ctx.IPAN[ibl][is_]
        ctx.XSSI[ibl][is_] = max(ctx.SST - ctx.S[i], xeps)

    is_ = 2
    ctx.XSSI[1][is_] = 0.0
    for ibl in range(2, ctx.IBLTE[is_] + 1):
        i = ctx.IPAN[ibl][is_]
        ctx.XSSI[ibl][is_] = max(ctx.S[i] - ctx.SST, xeps)

    is1 = 1
    is2 = 2
    ibl1 = ctx.IBLTE[is1] + 1
    ctx.XSSI[ibl1][is1] = ctx.XSSI[ibl1 - 1][is1]

    ibl2 = ctx.IBLTE[is2] + 1
    ctx.XSSI[ibl2][is2] = ctx.XSSI[ibl2 - 1][is2]

    for ibl in range(ctx.IBLTE[is_] + 2, ctx.NBL[is_] + 1):
        i = ctx.IPAN[ibl][is_]
        dxssi = math.sqrt((ctx.X[i] - ctx.X[i - 1]) ** 2 + (ctx.Y[i] - ctx.Y[i - 1]) ** 2)

        ibl1 = ctx.IBLTE[is1] + ibl - ctx.IBLTE[is_]
        ibl2 = ctx.IBLTE[is2] + ibl - ctx.IBLTE[is_]
        ctx.XSSI[ibl1][is1] = ctx.XSSI[ibl1 - 1][is1] + dxssi
        ctx.XSSI[ibl2][is2] = ctx.XSSI[ibl2 - 1][is2] + dxssi

    telrat = 2.50

    crosp = (ctx.XP[1] * ctx.YP[ctx.N] - ctx.YP[1] * ctx.XP[ctx.N]) / math.sqrt(
        (ctx.XP[1] ** 2 + ctx.YP[1] ** 2) * (ctx.XP[ctx.N] ** 2 + ctx.YP[ctx.N] ** 2)
    )
    dwdxte = crosp / math.sqrt(1.0 - crosp**2)

    dwdxte = max(dwdxte, -3.0 / telrat)
    dwdxte = min(dwdxte, 3.0 / telrat)

    aa = 3.0 + telrat * dwdxte
    bb = -2.0 - telrat * dwdxte

    if ctx.SHARP:
        for iw in range(1, ctx.NW + 1):
            ctx.WGAP[iw] = 0.0
    else:
        is_ = 2
        for iw in range(1, ctx.NW + 1):
            ibl = ctx.IBLTE[is_] + iw
            zn = 1.0 - (ctx.XSSI[ibl][is_] - ctx.XSSI[ctx.IBLTE[is_]][is_]) / (telrat * ctx.ANTE)
            ctx.WGAP[iw] = 0.0
            if zn >= 0.0:
                ctx.WGAP[iw] = ctx.ANTE * (aa + bb * zn) * zn**2


def uicalc(ctx):
    for is_ in range(1, 3):
        ctx.UINV[1][is_] = 0.0
        ctx.UINV_A[1][is_] = 0.0
        for ibl in range(2, ctx.NBL[is_] + 1):
            i = ctx.IPAN[ibl][is_]
            ctx.UINV[ibl][is_] = ctx.VTI[ibl][is_] * ctx.QINV[i]
            ctx.UINV_A[ibl][is_] = ctx.VTI[ibl][is_] * ctx.QINV_A[i]


def qvfue(ctx):
    for is_ in range(1, 3):
        for ibl in range(2, ctx.NBL[is_] + 1):
            i = ctx.IPAN[ibl][is_]
            ctx.QVIS[i] = ctx.VTI[ibl][is_] * ctx.UEDG[ibl][is_]


def qiset(ctx):
    ctx.COSA = math.cos(ctx.ALFA)
    ctx.SINA = math.sin(ctx.ALFA)

    for i in range(1, ctx.N + ctx.NW + 1):
        ctx.QINV[i] = ctx.COSA * ctx.QINVU[i][1] + ctx.SINA * ctx.QINVU[i][2]
        ctx.QINV_A[i] = -ctx.SINA * ctx.QINVU[i][1] + ctx.COSA * ctx.QINVU[i][2]


def gamqv(ctx):
    for i in range(1, ctx.N + 1):
        ctx.GAM[i] = ctx.QVIS[i]
        ctx.GAM_A[i] = ctx.QINV_A[i]


def stmove(ctx):
    istold = ctx.IST
    stfind(ctx)

    if istold == ctx.IST:
        xicalc(ctx)
    else:
        iblpan(ctx)
        uicalc(ctx)
        xicalc(ctx)
        from .xbl import iblsys

        iblsys(ctx)

        if ctx.IST > istold:
            idif = ctx.IST - istold

            ctx.ITRAN[1] = ctx.ITRAN[1] + idif
            ctx.ITRAN[2] = ctx.ITRAN[2] - idif

            for ibl in range(ctx.NBL[1], idif + 1, -1):
                ctx.CTAU[ibl][1] = ctx.CTAU[ibl - idif][1]
                ctx.THET[ibl][1] = ctx.THET[ibl - idif][1]
                ctx.DSTR[ibl][1] = ctx.DSTR[ibl - idif][1]
                ctx.UEDG[ibl][1] = ctx.UEDG[ibl - idif][1]

            dudx = ctx.UEDG[idif + 2][1] / ctx.XSSI[idif + 2][1]
            for ibl in range(idif + 1, 1, -1):
                ctx.CTAU[ibl][1] = ctx.CTAU[idif + 2][1]
                ctx.THET[ibl][1] = ctx.THET[idif + 2][1]
                ctx.DSTR[ibl][1] = ctx.DSTR[idif + 2][1]
                ctx.UEDG[ibl][1] = dudx * ctx.XSSI[ibl][1]

            for ibl in range(2, ctx.NBL[2] + 1):
                ctx.CTAU[ibl][2] = ctx.CTAU[ibl + idif][2]
                ctx.THET[ibl][2] = ctx.THET[ibl + idif][2]
                ctx.DSTR[ibl][2] = ctx.DSTR[ibl + idif][2]
                ctx.UEDG[ibl][2] = ctx.UEDG[ibl + idif][2]

        else:
            idif = istold - ctx.IST

            ctx.ITRAN[1] = ctx.ITRAN[1] - idif
            ctx.ITRAN[2] = ctx.ITRAN[2] + idif

            for ibl in range(ctx.NBL[2], idif + 1, -1):
                ctx.CTAU[ibl][2] = ctx.CTAU[ibl - idif][2]
                ctx.THET[ibl][2] = ctx.THET[ibl - idif][2]
                ctx.DSTR[ibl][2] = ctx.DSTR[ibl - idif][2]
                ctx.UEDG[ibl][2] = ctx.UEDG[ibl - idif][2]

            dudx = ctx.UEDG[idif + 2][2] / ctx.XSSI[idif + 2][2]
            for ibl in range(idif + 1, 1, -1):
                ctx.CTAU[ibl][2] = ctx.CTAU[idif + 2][2]
                ctx.THET[ibl][2] = ctx.THET[idif + 2][2]
                ctx.DSTR[ibl][2] = ctx.DSTR[idif + 2][2]
                ctx.UEDG[ibl][2] = dudx * ctx.XSSI[ibl][2]

            for ibl in range(2, ctx.NBL[1] + 1):
                ctx.CTAU[ibl][1] = ctx.CTAU[ibl + idif][1]
                ctx.THET[ibl][1] = ctx.THET[ibl + idif][1]
                ctx.DSTR[ibl][1] = ctx.DSTR[ibl + idif][1]
                ctx.UEDG[ibl][1] = ctx.UEDG[ibl + idif][1]

        ueps = 1.0e-7
        for is_ in range(1, 3):
            for ibl in range(2, ctx.NBL[is_] + 1):
                i = ctx.IPAN[ibl][is_]
                if ctx.UEDG[ibl][is_] <= ueps:
                    ctx.UEDG[ibl][is_] = ueps
                    if 1 <= i <= ctx.N + ctx.NW:
                        ctx.QVIS[i] = ctx.VTI[ibl][is_] * ueps
                    if 1 <= i <= ctx.N:
                        ctx.GAM[i] = ctx.VTI[ibl][is_] * ueps

    for is_ in range(1, 3):
        for ibl in range(2, ctx.NBL[is_] + 1):
            ctx.MASS[ibl][is_] = ctx.DSTR[ibl][is_] * ctx.UEDG[ibl][is_]


def ueset(ctx):
    for is_ in range(1, 3):
        for ibl in range(2, ctx.NBL[is_] + 1):
            i = ctx.IPAN[ibl][is_]

            dui = 0.0
            for js in range(1, 3):
                for jbl in range(2, ctx.NBL[js] + 1):
                    j = ctx.IPAN[jbl][js]
                    ue_m = -ctx.VTI[ibl][is_] * ctx.VTI[jbl][js] * ctx.DIJ[i][j]
                    dui = dui + ue_m * ctx.MASS[jbl][js]

            ctx.UEDG[ibl][is_] = ctx.UINV[ibl][is_] + dui


def dsset(ctx):
    for is_ in range(1, 3):
        for ibl in range(2, ctx.NBL[is_] + 1):
            ctx.DSTR[ibl][is_] = ctx.MASS[ibl][is_] / ctx.UEDG[ibl][is_]
