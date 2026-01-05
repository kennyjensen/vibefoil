# Ported from XFOIL Fortran source (Mark Drela).
# This file is a derived work and remains under the terms of the
# GNU General Public License v2 or later.
# See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

import math

from .spline import scalc, splind, trisol
from .xpanel import ncalc, psilin
from .xfoil import clcalc, tecalc
from .xsolve import gauss
from .xutils import atanc
from .xbl import IQX


def QDES(_ctx):
    raise NotImplementedError("QDES interactive command loop is not ported.")


def NEWPLOTQ(_ctx):
    return None


def QPLINI(_ctx, _ldef):
    return None


def QSPLOT(_ctx):
    return None


def QSPPLT(_ctx, _iqspl1, _iqspl2, _kqsp, _nt):
    return None


def IQSGET(_ctx):
    raise NotImplementedError("IQSGET cursor selection is not ported.")


def SPLQSP(ctx, kqsp):
    if ctx.NSP < 2:
        return

    nmid = ctx.NSP - 2
    if nmid > 0:
        q_local = [0.0] * (nmid + 1)
        s_local = [0.0] * (nmid + 1)
        qp_local = [0.0] * (nmid + 1)
        for i in range(1, nmid + 1):
            idx = i + 1
            q_local[i] = ctx.QSPEC[idx][kqsp]
            s_local[i] = ctx.SSPEC[idx]
        splind(q_local, qp_local, s_local, nmid, -999.0, -999.0)
        for i in range(1, nmid + 1):
            idx = i + 1
            ctx.QSPECP[idx][kqsp] = qp_local[i]

    i = 1
    q_local = [0.0] * 3
    s_local = [0.0] * 3
    qp_local = [0.0] * 3
    q_local[1] = ctx.QSPEC[i][kqsp]
    q_local[2] = ctx.QSPEC[i + 1][kqsp]
    s_local[1] = ctx.SSPEC[i]
    s_local[2] = ctx.SSPEC[i + 1]
    splind(q_local, qp_local, s_local, 2, -999.0, ctx.QSPECP[i + 1][kqsp])
    ctx.QSPECP[i][kqsp] = qp_local[1]

    i = ctx.NSP - 1
    q_local[1] = ctx.QSPEC[i][kqsp]
    q_local[2] = ctx.QSPEC[i + 1][kqsp]
    s_local[1] = ctx.SSPEC[i]
    s_local[2] = ctx.SSPEC[i + 1]
    splind(q_local, qp_local, s_local, 2, ctx.QSPECP[i][kqsp], -999.0)
    ctx.QSPECP[i][kqsp] = qp_local[1]


def SMOOQ(ctx, kq1, kq2, kqsp):
    for i in range(1, ctx.NSP + 1):
        ctx.W8[i] = ctx.SSPEC[i]

    if kq2 - kq1 < 2:
        return

    smool = 0.002 * (ctx.W8[ctx.NSP] - ctx.W8[1])
    smoosq = smool * smool

    for i in range(kq1 + 1, kq2):
        dsm = ctx.W8[i] - ctx.W8[i - 1]
        dsp = ctx.W8[i + 1] - ctx.W8[i]
        dso = 0.5 * (ctx.W8[i + 1] - ctx.W8[i - 1])
        ctx.W1[i] = smoosq * (-1.0 / dsm) / dso
        ctx.W2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0
        ctx.W3[i] = smoosq * (-1.0 / dsp) / dso

    ctx.W2[kq1] = 1.0
    ctx.W3[kq1] = 0.0
    ctx.W1[kq2] = 0.0
    ctx.W2[kq2] = 1.0

    if ctx.LQSLOP:
        i = kq1 + 1
        dsm = ctx.W8[i] - ctx.W8[i - 1]
        dsp = ctx.W8[i + 1] - ctx.W8[i]
        ds = ctx.W8[i + 1] - ctx.W8[i - 1]
        ctx.W1[i] = -1.0 / dsm - (dsm / ds) / dsm
        ctx.W2[i] = 1.0 / dsm + (dsm / ds) / dsm + (dsm / ds) / dsp
        ctx.W3[i] = -(dsm / ds) / dsp
        qsp1 = (
            ctx.W1[i] * ctx.QSPEC[i - 1][kqsp]
            + ctx.W2[i] * ctx.QSPEC[i][kqsp]
            + ctx.W3[i] * ctx.QSPEC[i + 1][kqsp]
        )

        i = kq2 - 1
        dsm = ctx.W8[i] - ctx.W8[i - 1]
        dsp = ctx.W8[i + 1] - ctx.W8[i]
        ds = ctx.W8[i + 1] - ctx.W8[i - 1]
        ctx.W1[i] = (dsp / ds) / dsm
        ctx.W2[i] = -1.0 / dsp - (dsp / ds) / dsp - (dsp / ds) / dsm
        ctx.W3[i] = 1.0 / dsp + (dsp / ds) / dsp
        qsp2 = (
            ctx.W1[i] * ctx.QSPEC[i - 1][kqsp]
            + ctx.W2[i] * ctx.QSPEC[i][kqsp]
            + ctx.W3[i] * ctx.QSPEC[i + 1][kqsp]
        )

        ctx.QSPEC[kq1 + 1][kqsp] = qsp1
        ctx.QSPEC[kq2 - 1][kqsp] = qsp2

    kk = kq2 - kq1 + 1
    a = [0.0] * (kk + 1)
    b = [0.0] * (kk + 1)
    c = [0.0] * (kk + 1)
    d = [0.0] * (kk + 1)
    for i in range(1, kk + 1):
        idx = kq1 + i - 1
        a[i] = ctx.W2[idx]
        b[i] = ctx.W1[idx]
        c[i] = ctx.W3[idx]
        d[i] = ctx.QSPEC[idx][kqsp]

    trisol(a, b, c, d, kk)

    for i in range(1, kk + 1):
        idx = kq1 + i - 1
        ctx.QSPEC[idx][kqsp] = d[i]


def QINCOM(qc, qinf, tklam):
    if tklam < 1.0e-4 or abs(qc) < 1.0e-4:
        return qc / (1.0 - tklam)
    tmp = 0.5 * (1.0 - tklam) * qinf / (qc * tklam)
    return qinf * tmp * (math.sqrt(1.0 + 1.0 / (tklam * tmp * tmp)) - 1.0)


def GAMQSP(ctx, kqsp):
    ctx.ALQSP[kqsp] = ctx.ALGAM
    ctx.CLQSP[kqsp] = ctx.CLGAM
    ctx.CMQSP[kqsp] = ctx.CMGAM

    for i in range(1, ctx.NSP + 1):
        ctx.QSPEC[i][kqsp] = ctx.QGAMM[i]

    ctx.QDOF0 = 0.0
    ctx.QDOF1 = 0.0
    ctx.QDOF2 = 0.0
    ctx.QDOF3 = 0.0

    SPLQSP(ctx, kqsp)

    if not ctx.LIQSET:
        ctx.IQ1 = 1
        ctx.IQ2 = ctx.NSP


def SYMQSP(ctx, kqsp):
    ctx.ALQSP[kqsp] = 0.0
    ctx.CLQSP[kqsp] = 0.0
    ctx.CMQSP[kqsp] = 0.0

    sspmid = 0.5 * (ctx.SSPEC[ctx.NSP] - ctx.SSPEC[1])
    for i in range(1, (ctx.NSP + 1) // 2 + 1):
        ctx.SSPEC[i] = sspmid + 0.5 * (ctx.SSPEC[i] - ctx.SSPEC[ctx.NSP - i + 1])
        ctx.QSPEC[i][kqsp] = 0.5 * (ctx.QSPEC[i][kqsp] - ctx.QSPEC[ctx.NSP - i + 1][kqsp])

    for i in range((ctx.NSP + 1) // 2 + 1, ctx.NSP + 1):
        ctx.SSPEC[i] = -ctx.SSPEC[ctx.NSP - i + 1] + 2.0 * sspmid
        ctx.QSPEC[i][kqsp] = -ctx.QSPEC[ctx.NSP - i + 1][kqsp]

    ctx.QDOF0 = 0.0
    ctx.QDOF1 = 0.0
    ctx.QDOF2 = 0.0
    ctx.QDOF3 = 0.0

    SPLQSP(ctx, kqsp)


def MIXED(ctx, kqsp, niterq):
    bwt = 0.1

    ctx.COSA = math.cos(ctx.ALFA)
    ctx.SINA = math.sin(ctx.ALFA)
    scalc(ctx.X, ctx.Y, ctx.S, ctx.N)

    for i in range(1, ctx.N + 1):
        ctx.QF0[i] = 0.0
        ctx.QF1[i] = 0.0
        ctx.QF2[i] = 0.0
        ctx.QF3[i] = 0.0

    for i in range(ctx.IQ1, ctx.IQ2 + 1):
        fs = (ctx.S[i] - ctx.S[ctx.IQ1]) / (ctx.S[ctx.IQ2] - ctx.S[ctx.IQ1])
        ctx.QF0[i] = 1.0 - fs
        ctx.QF1[i] = fs
        if ctx.LCPXX:
            ctx.QF2[i] = math.exp(-5.0 * fs)
            ctx.QF3[i] = math.exp(-5.0 * (1.0 - fs))
        else:
            ctx.QF2[i] = 0.0
            ctx.QF3[i] = 0.0
        ctx.GAM[i] = (
            ctx.QSPEC[i][kqsp]
            + ctx.QDOF0 * ctx.QF0[i]
            + ctx.QDOF1 * ctx.QF1[i]
            + ctx.QDOF2 * ctx.QF2[i]
            + ctx.QDOF3 * ctx.QF3[i]
        )

    for _ in range(1, niterq + 1):
        for i in range(1, ctx.N + 6):
            for j in range(1, ctx.N + 6):
                ctx.Q[i][j] = 0.0

        ncalc(ctx.X, ctx.Y, ctx.S, ctx.N, ctx.NX, ctx.NY)

        for i in range(1, ctx.N + 1):
            psi, psi_n, _, _, _ = psilin(ctx, i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i], True, False)
            ctx.DZDN[i] = ctx.DZDN[i] + psi_n

            for j in range(1, ctx.IQ1):
                ctx.Q[i][j] = ctx.Q[i][j] + ctx.DZDG[j]
            for j in range(ctx.IQ1, ctx.IQ2 + 1):
                ctx.Q[i][j] = ctx.Q[i][j] + ctx.DZDN[j]
            for j in range(ctx.IQ2 + 1, ctx.N + 1):
                ctx.Q[i][j] = ctx.Q[i][j] + ctx.DZDG[j]

            ctx.DQ[i] = ctx.PSIO - psi

            ctx.Q[i][ctx.N + 1] = ctx.Q[i][ctx.N + 1] - 1.0
            ctx.Q[i][ctx.N + 2] = ctx.Q[i][ctx.N + 2] + ctx.Z_QDOF0
            ctx.Q[i][ctx.N + 3] = ctx.Q[i][ctx.N + 3] + ctx.Z_QDOF1
            ctx.Q[i][ctx.N + 4] = ctx.Q[i][ctx.N + 4] + ctx.Z_QDOF2
            ctx.Q[i][ctx.N + 5] = ctx.Q[i][ctx.N + 5] + ctx.Z_QDOF3

        ctx.DQ[ctx.N + 1] = -(ctx.GAM[1] + ctx.GAM[ctx.N])
        GAMLIN(ctx, ctx.N + 1, 1, 1.0)
        GAMLIN(ctx, ctx.N + 1, ctx.N, 1.0)

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
            _, qbis, _, _, _ = psilin(ctx, 0, xbis, ybis, -sbis, cbis, False, True)
            res = qbis

            for j in range(1, ctx.N + 6):
                ctx.Q[ctx.N][j] = 0.0

            for j in range(1, ctx.N + 1):
                GAMLIN(ctx, ctx.N, j, ctx.DQDG[j])
                ctx.Q[ctx.N][j] = ctx.DQDG[j]

            ctx.Q[ctx.N][ctx.N + 1] = 0.0
            ctx.DQ[ctx.N] = -res

        ctx.Q[ctx.N + 2][ctx.IQ1] = 1.0
        ctx.DQ[ctx.N + 2] = 0.0
        ctx.Q[ctx.N + 3][ctx.IQ2] = 1.0
        ctx.DQ[ctx.N + 3] = 0.0

        if ctx.IQ1 > 1 and ctx.LCPXX:
            res = (
                ctx.GAM[ctx.IQ1 - 1]
                - 2.0 * ctx.GAM[ctx.IQ1]
                + ctx.GAM[ctx.IQ1 + 1]
                - (
                    ctx.QSPEC[ctx.IQ1 - 1][kqsp]
                    - 2.0 * ctx.QSPEC[ctx.IQ1][kqsp]
                    + ctx.QSPEC[ctx.IQ1 + 1][kqsp]
                )
            )
            GAMLIN(ctx, ctx.N + 4, ctx.IQ1 - 1, 1.0)
            GAMLIN(ctx, ctx.N + 4, ctx.IQ1, -2.0)
            GAMLIN(ctx, ctx.N + 4, ctx.IQ1 + 1, 1.0)
            ctx.DQ[ctx.N + 4] = -res
        else:
            ctx.Q[ctx.N + 4][ctx.N + 4] = 1.0
            ctx.DQ[ctx.N + 4] = -ctx.QDOF2

        if ctx.IQ2 < ctx.N and ctx.LCPXX:
            res = (
                ctx.GAM[ctx.IQ2 - 1]
                - 2.0 * ctx.GAM[ctx.IQ2]
                + ctx.GAM[ctx.IQ2 + 1]
                - (
                    ctx.QSPEC[ctx.IQ2 - 1][kqsp]
                    - 2.0 * ctx.QSPEC[ctx.IQ2][kqsp]
                    + ctx.QSPEC[ctx.IQ2 + 1][kqsp]
                )
            )
            GAMLIN(ctx, ctx.N + 5, ctx.IQ2 - 1, 1.0)
            GAMLIN(ctx, ctx.N + 5, ctx.IQ2, -2.0)
            GAMLIN(ctx, ctx.N + 5, ctx.IQ2 + 1, 1.0)
            ctx.DQ[ctx.N + 5] = -res
        else:
            ctx.Q[ctx.N + 5][ctx.N + 5] = 1.0
            ctx.DQ[ctx.N + 5] = -ctx.QDOF3

        dq_mat = [[0.0, 0.0] for _ in range(ctx.N + 6)]
        for i in range(1, ctx.N + 6):
            dq_mat[i][1] = ctx.DQ[i]

        gauss(IQX, ctx.N + 5, ctx.Q, dq_mat, 1)

        for i in range(1, ctx.N + 6):
            ctx.DQ[i] = dq_mat[i][1]

        inmax = 0
        igmax = 0
        dnmax = 0.0
        dgmax = 0.0

        for i in range(1, ctx.IQ1):
            ctx.GAM[i] = ctx.GAM[i] + ctx.DQ[i]
            if abs(ctx.DQ[i]) > abs(dgmax):
                dgmax = ctx.DQ[i]
                igmax = i

        for i in range(ctx.IQ1, ctx.IQ2 + 1):
            ctx.X[i] = ctx.X[i] + ctx.NX[i] * ctx.DQ[i]
            ctx.Y[i] = ctx.Y[i] + ctx.NY[i] * ctx.DQ[i]
            if abs(ctx.DQ[i]) > abs(dnmax):
                dnmax = ctx.DQ[i]
                inmax = i

        for i in range(ctx.IQ2 + 1, ctx.N + 1):
            ctx.GAM[i] = ctx.GAM[i] + ctx.DQ[i]
            if abs(ctx.DQ[i]) > abs(dgmax):
                dgmax = ctx.DQ[i]
                igmax = i

        ctx.PSIO = ctx.PSIO + ctx.DQ[ctx.N + 1]
        ctx.QDOF0 = ctx.QDOF0 + ctx.DQ[ctx.N + 2]
        ctx.QDOF1 = ctx.QDOF1 + ctx.DQ[ctx.N + 3]
        ctx.QDOF2 = ctx.QDOF2 + ctx.DQ[ctx.N + 4]
        ctx.QDOF3 = ctx.QDOF3 + ctx.DQ[ctx.N + 5]

        ctx.COSA = math.cos(ctx.ALFA)
        ctx.SINA = math.sin(ctx.ALFA)
        scalc(ctx.X, ctx.Y, ctx.S, ctx.N)

        for i in range(ctx.IQ1, ctx.IQ2 + 1):
            ctx.GAM[i] = (
                ctx.QSPEC[i][kqsp]
                + ctx.QDOF0 * ctx.QF0[i]
                + ctx.QDOF1 * ctx.QF1[i]
                + ctx.QDOF2 * ctx.QF2[i]
                + ctx.QDOF3 * ctx.QF3[i]
            )

        tecalc(ctx)
        cl, cm, cdp, cl_alf, cl_msq = clcalc(
            ctx.N,
            ctx.X,
            ctx.Y,
            ctx.GAM,
            ctx.GAM_A,
            ctx.ALFA,
            ctx.MINF,
            ctx.QINF,
            ctx.XCMREF,
            ctx.YCMREF,
        )
        ctx.CL = cl
        ctx.CM = cm
        ctx.CDP = cdp
        ctx.CL_ALF = cl_alf
        ctx.CL_MSQ = cl_msq

        if abs(dnmax) < 5.0e-5 and abs(dgmax) < 5.0e-4:
            return


def GAMLIN(ctx, i, j, coef):
    if ctx.IQ1 <= j <= ctx.IQ2:
        ctx.Q[i][ctx.N + 2] = ctx.Q[i][ctx.N + 2] + coef * ctx.QF0[j]
        ctx.Q[i][ctx.N + 3] = ctx.Q[i][ctx.N + 3] + coef * ctx.QF1[j]
        ctx.Q[i][ctx.N + 4] = ctx.Q[i][ctx.N + 4] + coef * ctx.QF2[j]
        ctx.Q[i][ctx.N + 5] = ctx.Q[i][ctx.N + 5] + coef * ctx.QF3[j]
    else:
        ctx.Q[i][j] = ctx.Q[i][j] + coef
