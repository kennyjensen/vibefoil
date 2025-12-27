import math

from .naca import naca4, naca5
from .spline import curv, deval, scalc, seval, segspl, trisol
from .userio import strip_string
from .xgeom import geopar, lefind
from .xpanel import apcalc, ncalc


def mrcl(ctx, cls):
    cla = max(cls, 0.000001)

    if ctx.RETYP < 1 or ctx.RETYP > 3:
        print("MRCL:  Illegal Re(CL) dependence trigger.")
        print("       Setting fixed Re.")
        ctx.RETYP = 1
    if ctx.MATYP < 1 or ctx.MATYP > 3:
        print("MRCL:  Illegal Mach(CL) dependence trigger.")
        print("       Setting fixed Mach.")
        ctx.MATYP = 1

    if ctx.MATYP == 1:
        ctx.MINF = ctx.MINF1
        m_cls = 0.0
    elif ctx.MATYP == 2:
        ctx.MINF = ctx.MINF1 / math.sqrt(cla)
        m_cls = -0.5 * ctx.MINF / cla
    else:
        ctx.MINF = ctx.MINF1
        m_cls = 0.0

    if ctx.RETYP == 1:
        ctx.REINF = ctx.REINF1
        r_cls = 0.0
    elif ctx.RETYP == 2:
        ctx.REINF = ctx.REINF1 / math.sqrt(cla)
        r_cls = -0.5 * ctx.REINF / cla
    else:
        ctx.REINF = ctx.REINF1 / cla
        r_cls = -ctx.REINF / cla

    if ctx.MINF >= 0.99:
        print()
        print("MRCL: CL too low for chosen Mach(CL) dependence")
        print("      Aritificially limiting Mach to  0.99")
        ctx.MINF = 0.99
        m_cls = 0.0

    rrat = 1.0
    if ctx.REINF1 > 0.0:
        rrat = ctx.REINF / ctx.REINF1

    if rrat > 100.0:
        print()
        print("MRCL: CL too low for chosen Re(CL) dependence")
        print("      Aritificially limiting Re to ", ctx.REINF1 * 100.0)
        ctx.REINF = ctx.REINF1 * 100.0
        r_cls = 0.0

    ctx.MINF_CL = m_cls
    ctx.REINF_CL = r_cls
    return m_cls, r_cls


def comset(ctx):
    beta = math.sqrt(1.0 - ctx.MINF**2)
    beta_msq = -0.5 / beta

    ctx.TKLAM = ctx.MINF**2 / (1.0 + beta) ** 2
    ctx.TKL_MSQ = 1.0 / (1.0 + beta) ** 2 - 2.0 * ctx.TKLAM / (1.0 + beta) * beta_msq

    if ctx.MINF == 0.0:
        ctx.CPSTAR = -999.0
        ctx.QSTAR = 999.0
    else:
        ctx.CPSTAR = (
            2.0
            / (ctx.GAMMA * ctx.MINF**2)
            * (
                ((1.0 + 0.5 * ctx.GAMM1 * ctx.MINF**2) / (1.0 + 0.5 * ctx.GAMM1)) ** (ctx.GAMMA / ctx.GAMM1)
                - 1.0
            )
        )
        ctx.QSTAR = ctx.QINF / ctx.MINF * math.sqrt(
            (1.0 + 0.5 * ctx.GAMM1 * ctx.MINF**2) / (1.0 + 0.5 * ctx.GAMM1)
        )


def cpcalc(n, q, qinf, minf, cp):
    beta = math.sqrt(1.0 - minf**2)
    bfac = 0.5 * minf**2 / (1.0 + beta)

    denneg = False

    for i in range(1, n + 1):
        cpinc = 1.0 - (q[i] / qinf) ** 2
        den = beta + bfac * cpinc
        cp[i] = cpinc / den
        if den <= 0.0:
            denneg = True

    if denneg:
        print()
        print("CPCALC: Local speed too large. Compressibility corrections invalid.")


def clcalc(n, x, y, gam, gam_a, alfa, minf, qinf, xref, yref):
    sa = math.sin(alfa)
    ca = math.cos(alfa)

    beta = math.sqrt(1.0 - minf**2)
    beta_msq = -0.5 / beta

    bfac = 0.5 * minf**2 / (1.0 + beta)
    bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq

    cl = 0.0
    cm = 0.0
    cdp = 0.0
    cl_alf = 0.0
    cl_msq = 0.0

    i = 1
    cginc = 1.0 - (gam[i] / qinf) ** 2
    cpg1 = cginc / (beta + bfac * cginc)
    cpg1_msq = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc)

    cpi_gam = -2.0 * gam[i] / qinf**2
    cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc)
    cpg1_alf = cpc_cpi * cpi_gam * gam_a[i]

    for i in range(1, n + 1):
        ip = i + 1
        if i == n:
            ip = 1

        cginc = 1.0 - (gam[ip] / qinf) ** 2
        cpg2 = cginc / (beta + bfac * cginc)
        cpg2_msq = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc)

        cpi_gam = -2.0 * gam[ip] / qinf**2
        cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc)
        cpg2_alf = cpc_cpi * cpi_gam * gam_a[ip]

        dx = (x[ip] - x[i]) * ca + (y[ip] - y[i]) * sa
        dy = (y[ip] - y[i]) * ca - (x[ip] - x[i]) * sa
        dg = cpg2 - cpg1

        ax = (0.5 * (x[ip] + x[i]) - xref) * ca + (0.5 * (y[ip] + y[i]) - yref) * sa
        ay = (0.5 * (y[ip] + y[i]) - yref) * ca - (0.5 * (x[ip] + x[i]) - xref) * sa
        ag = 0.5 * (cpg2 + cpg1)

        dx_alf = -(x[ip] - x[i]) * sa + (y[ip] - y[i]) * ca
        ag_alf = 0.5 * (cpg2_alf + cpg1_alf)
        ag_msq = 0.5 * (cpg2_msq + cpg1_msq)

        cl = cl + dx * ag
        cdp = cdp - dy * ag
        cm = cm - dx * (ag * ax + dg * dx / 12.0) - dy * (ag * ay + dg * dy / 12.0)

        cl_alf = cl_alf + dx * ag_alf + ag * dx_alf
        cl_msq = cl_msq + dx * ag_msq

        cpg1 = cpg2
        cpg1_alf = cpg2_alf
        cpg1_msq = cpg2_msq

    return cl, cm, cdp, cl_alf, cl_msq


def cdcalc(ctx):
    sa = math.sin(ctx.ALFA)
    ca = math.cos(ctx.ALFA)

    if ctx.LVISC and ctx.LBLINI:
        thwake = ctx.THET[ctx.NBL[2]][2]
        urat = ctx.UEDG[ctx.NBL[2]][2] / ctx.QINF
        uewake = ctx.UEDG[ctx.NBL[2]][2] * (1.0 - ctx.TKLAM) / (1.0 - ctx.TKLAM * urat**2)
        shwake = ctx.DSTR[ctx.NBL[2]][2] / ctx.THET[ctx.NBL[2]][2]
        ctx.CD = 2.0 * thwake * (uewake / ctx.QINF) ** (0.5 * (5.0 + shwake))
    else:
        ctx.CD = 0.0

    ctx.CDF = 0.0
    for is_ in range(1, 3):
        for ibl in range(3, ctx.IBLTE[is_] + 1):
            i = ctx.IPAN[ibl][is_]
            im = ctx.IPAN[ibl - 1][is_]
            dx = (ctx.X[i] - ctx.X[im]) * ca + (ctx.Y[i] - ctx.Y[im]) * sa
            ctx.CDF = ctx.CDF + 0.5 * (ctx.TAU[ibl][is_] + ctx.TAU[ibl - 1][is_]) * dx * 2.0 / ctx.QINF**2


def tecalc(ctx):
    dxte = ctx.X[1] - ctx.X[ctx.N]
    dyte = ctx.Y[1] - ctx.Y[ctx.N]
    dxs = 0.5 * (-ctx.XP[1] + ctx.XP[ctx.N])
    dys = 0.5 * (-ctx.YP[1] + ctx.YP[ctx.N])

    ctx.ANTE = dxs * dyte - dys * dxte
    ctx.ASTE = dxs * dxte + dys * dyte

    ctx.DSTE = math.sqrt(dxte**2 + dyte**2)

    ctx.SHARP = ctx.DSTE < 0.0001 * ctx.CHORD

    if ctx.SHARP:
        scs = 1.0
        sds = 0.0
    else:
        scs = ctx.ANTE / ctx.DSTE
        sds = ctx.ASTE / ctx.DSTE

    ctx.SIGTE = 0.5 * (ctx.GAM[1] - ctx.GAM[ctx.N]) * scs
    ctx.GAMTE = -0.5 * (ctx.GAM[1] - ctx.GAM[ctx.N]) * sds

    ctx.SIGTE_A = 0.5 * (ctx.GAM_A[1] - ctx.GAM_A[ctx.N]) * scs
    ctx.GAMTE_A = -0.5 * (ctx.GAM_A[1] - ctx.GAM_A[ctx.N]) * sds


def naca(ctx, ides1):
    iqx = (len(ctx.W1) - 1) // 6
    nside = iqx // 3

    if ides1 <= 0:
        raise RuntimeError("NACA: IDES must be specified.")
    else:
        ides = ides1

    itype = 0
    if ides <= 25099:
        itype = 5
    if ides <= 9999:
        itype = 4

    if itype == 0:
        print("This designation not implemented.")
        return

    xx = ctx.W1
    yt = ctx.W2
    yc = ctx.W3

    if itype == 4:
        nb, name = naca4(ides, xx, yt, yc, nside, ctx.XB, ctx.YB)
    else:
        nb, name = naca5(ides, xx, yt, yc, nside, ctx.XB, ctx.YB)

    ctx.NB = nb
    ctx.NAME = name
    ctx.NAME, ctx.NNAME = strip_string(ctx.NAME)

    if ctx.NB == 0:
        return

    ctx.LCLOCK = False

    ctx.XBF = 0.0
    ctx.YBF = 0.0
    ctx.LBFLAP = False

    scalc(ctx.XB, ctx.YB, ctx.SB, ctx.NB)
    segspl(ctx.XB, ctx.XBP, ctx.SB, ctx.NB)
    segspl(ctx.YB, ctx.YBP, ctx.SB, ctx.NB)

    gp = geopar(ctx.XB, ctx.XBP, ctx.YB, ctx.YBP, ctx.SB, ctx.NB, ctx.W1)
    ctx.SBLE = gp["sle"]
    ctx.CHORDB = gp["chord"]
    ctx.AREAB = gp["area"]
    ctx.RADBLE = gp["radle"]
    ctx.ANGBTE = gp["angte"]
    ctx.EI11BA = gp["ei11a"]
    ctx.EI22BA = gp["ei22a"]
    ctx.APX1BA = gp["apx1a"]
    ctx.APX2BA = gp["apx2a"]
    ctx.EI11BT = gp["ei11t"]
    ctx.EI22BT = gp["ei22t"]
    ctx.APX1BT = gp["apx1t"]
    ctx.APX2BT = gp["apx2t"]
    ctx.THICKB = gp["thick"]
    ctx.CAMBRB = gp["cambr"]

    print(f"\n Buffer airfoil set using{ctx.NB:4d} points")

    pangen(ctx, True)


def pangen(ctx, shopar):
    if ctx.NB < 2:
        print("PANGEN: Buffer airfoil not available.")
        ctx.N = 0
        return

    ipfac = 5

    ctx.N = ctx.NPAN

    scalc(ctx.XB, ctx.YB, ctx.SB, ctx.NB)
    segspl(ctx.XB, ctx.XBP, ctx.SB, ctx.NB)
    segspl(ctx.YB, ctx.YBP, ctx.SB, ctx.NB)

    sbref = 0.5 * (ctx.SB[ctx.NB] - ctx.SB[1])

    for i in range(1, ctx.NB + 1):
        ctx.W5[i] = abs(curv(ctx.SB[i], ctx.XB, ctx.XBP, ctx.YB, ctx.YBP, ctx.SB, ctx.NB)) * sbref

    ctx.SBLE = lefind(ctx.XB, ctx.XBP, ctx.YB, ctx.YBP, ctx.SB, ctx.NB)
    cvle = abs(curv(ctx.SBLE, ctx.XB, ctx.XBP, ctx.YB, ctx.YBP, ctx.SB, ctx.NB)) * sbref

    ible = 0
    for i in range(1, ctx.NB):
        if ctx.SBLE == ctx.SB[i] and ctx.SBLE == ctx.SB[i + 1]:
            ible = i
            print()
            print("Sharp leading edge")
            break

    xble = seval(ctx.SBLE, ctx.XB, ctx.XBP, ctx.SB, ctx.NB)
    yble = seval(ctx.SBLE, ctx.YB, ctx.YBP, ctx.SB, ctx.NB)
    xbte = 0.5 * (ctx.XB[1] + ctx.XB[ctx.NB])
    ybte = 0.5 * (ctx.YB[1] + ctx.YB[ctx.NB])
    chbsq = (xbte - xble) ** 2 + (ybte - yble) ** 2

    nk = 3
    cvsum = 0.0
    for k in range(-nk, nk + 1):
        frac = float(k) / float(nk)
        sbk = ctx.SBLE + frac * sbref / max(cvle, 20.0)
        cvk = abs(curv(sbk, ctx.XB, ctx.XBP, ctx.YB, ctx.YBP, ctx.SB, ctx.NB)) * sbref
        cvsum = cvsum + cvk
    cvavg = cvsum / float(2 * nk + 1)

    if ible != 0:
        cvavg = 10.0

    cc = 6.0 * ctx.CVPAR

    cvte = cvavg * ctx.CTERAT
    ctx.W5[1] = cvte
    ctx.W5[ctx.NB] = cvte

    smool = max(1.0 / max(cvavg, 20.0), 0.25 / float(ctx.NPAN // 2))
    smoosq = (smool * sbref) ** 2

    ctx.W2[1] = 1.0
    ctx.W3[1] = 0.0
    for i in range(2, ctx.NB):
        dsm = ctx.SB[i] - ctx.SB[i - 1]
        dsp = ctx.SB[i + 1] - ctx.SB[i]
        dso = 0.5 * (ctx.SB[i + 1] - ctx.SB[i - 1])

        if dsm == 0.0 or dsp == 0.0:
            ctx.W1[i] = 0.0
            ctx.W2[i] = 1.0
            ctx.W3[i] = 0.0
        else:
            ctx.W1[i] = smoosq * (-1.0 / dsm) / dso
            ctx.W2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0
            ctx.W3[i] = smoosq * (-1.0 / dsp) / dso
    ctx.W1[ctx.NB] = 0.0
    ctx.W2[ctx.NB] = 1.0

    for i in range(2, ctx.NB):
        if ctx.SB[i] == ctx.SBLE or i == ible or i == ible + 1:
            ctx.W1[i] = 0.0
            ctx.W2[i] = 1.0
            ctx.W3[i] = 0.0
            ctx.W5[i] = cvle
        elif ctx.SB[i - 1] < ctx.SBLE and ctx.SB[i] > ctx.SBLE:
            dsm = ctx.SB[i - 1] - ctx.SB[i - 2]
            dsp = ctx.SBLE - ctx.SB[i - 1]
            dso = 0.5 * (ctx.SBLE - ctx.SB[i - 2])

            ctx.W1[i - 1] = smoosq * (-1.0 / dsm) / dso
            ctx.W2[i - 1] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0
            ctx.W3[i - 1] = 0.0
            ctx.W5[i - 1] = ctx.W5[i - 1] + smoosq * cvle / (dsp * dso)

            dsm = ctx.SB[i] - ctx.SBLE
            dsp = ctx.SB[i + 1] - ctx.SB[i]
            dso = 0.5 * (ctx.SB[i + 1] - ctx.SBLE)
            ctx.W1[i] = 0.0
            ctx.W2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0
            ctx.W3[i] = smoosq * (-1.0 / dsp) / dso
            ctx.W5[i] = ctx.W5[i] + smoosq * cvle / (dsm * dso)
            break

    for i in range(2, ctx.NB):
        xoc = ((ctx.XB[i] - xble) * (xbte - xble) + (ctx.YB[i] - yble) * (ybte - yble)) / chbsq

        if ctx.SB[i] < ctx.SBLE:
            if xoc > ctx.XSREF1 and xoc < ctx.XSREF2:
                ctx.W1[i] = 0.0
                ctx.W2[i] = 1.0
                ctx.W3[i] = 0.0
                ctx.W5[i] = cvle * ctx.CTRRAT
        else:
            if xoc > ctx.XPREF1 and xoc < ctx.XPREF2:
                ctx.W1[i] = 0.0
                ctx.W2[i] = 1.0
                ctx.W3[i] = 0.0
                ctx.W5[i] = cvle * ctx.CTRRAT

    if ible == 0:
        trisol(ctx.W2, ctx.W1, ctx.W3, ctx.W5, ctx.NB)
    else:
        nn1 = ible
        a = [0.0] * (nn1 + 1)
        b = [0.0] * (nn1 + 1)
        c = [0.0] * (nn1 + 1)
        d = [0.0] * (nn1 + 1)
        for i in range(1, nn1 + 1):
            a[i] = ctx.W2[i]
            b[i] = ctx.W1[i]
            c[i] = ctx.W3[i]
            d[i] = ctx.W5[i]
        trisol(a, b, c, d, nn1)
        for i in range(1, nn1 + 1):
            ctx.W5[i] = d[i]

        nn2 = ctx.NB - ible
        a = [0.0] * (nn2 + 1)
        b = [0.0] * (nn2 + 1)
        c = [0.0] * (nn2 + 1)
        d = [0.0] * (nn2 + 1)
        for i in range(1, nn2 + 1):
            idx = ible + i
            a[i] = ctx.W2[idx]
            b[i] = ctx.W1[idx]
            c[i] = ctx.W3[idx]
            d[i] = ctx.W5[idx]
        trisol(a, b, c, d, nn2)
        for i in range(1, nn2 + 1):
            ctx.W5[ible + i] = d[i]

    cvmax = 0.0
    for i in range(1, ctx.NB + 1):
        cvmax = max(cvmax, abs(ctx.W5[i]))

    for i in range(1, ctx.NB + 1):
        ctx.W5[i] = ctx.W5[i] / cvmax

    segspl(ctx.W5, ctx.W6, ctx.SB, ctx.NB)

    nn = ipfac * (ctx.N - 1) + 1

    rdste = 0.667
    rtf = (rdste - 1.0) * 2.0 + 1.0

    if ible == 0:
        dsavg = (ctx.SB[ctx.NB] - ctx.SB[1]) / (float(nn - 3) + 2.0 * rtf)
        ctx.SNEW[1] = ctx.SB[1]
        for i in range(2, nn):
            ctx.SNEW[i] = ctx.SB[1] + dsavg * (float(i - 2) + rtf)
        ctx.SNEW[nn] = ctx.SB[ctx.NB]
    else:
        nfrac1 = (ctx.N * ible) // ctx.NB
        nn1 = ipfac * (nfrac1 - 1) + 1
        dsavg1 = (ctx.SBLE - ctx.SB[1]) / (float(nn1 - 2) + rtf)
        ctx.SNEW[1] = ctx.SB[1]
        for i in range(2, nn1 + 1):
            ctx.SNEW[i] = ctx.SB[1] + dsavg1 * (float(i - 2) + rtf)

        nn2 = nn - nn1 + 1
        dsavg2 = (ctx.SB[ctx.NB] - ctx.SBLE) / (float(nn2 - 2) + rtf)
        for i in range(2, nn2):
            ctx.SNEW[i - 1 + nn1] = ctx.SBLE + dsavg2 * (float(i - 2) + rtf)
        ctx.SNEW[nn] = ctx.SB[ctx.NB]

    for _ in range(1, 21):
        cv1 = seval(ctx.SNEW[1], ctx.W5, ctx.W6, ctx.SB, ctx.NB)
        cv2 = seval(ctx.SNEW[2], ctx.W5, ctx.W6, ctx.SB, ctx.NB)
        cvs1 = deval(ctx.SNEW[1], ctx.W5, ctx.W6, ctx.SB, ctx.NB)
        cvs2 = deval(ctx.SNEW[2], ctx.W5, ctx.W6, ctx.SB, ctx.NB)

        cavm = math.sqrt(cv1**2 + cv2**2)
        if cavm == 0.0:
            cavm_s1 = 0.0
            cavm_s2 = 0.0
        else:
            cavm_s1 = cvs1 * cv1 / cavm
            cavm_s2 = cvs2 * cv2 / cavm

        for i in range(2, nn):
            dsm = ctx.SNEW[i] - ctx.SNEW[i - 1]
            dsp = ctx.SNEW[i] - ctx.SNEW[i + 1]
            cv3 = seval(ctx.SNEW[i + 1], ctx.W5, ctx.W6, ctx.SB, ctx.NB)
            cvs3 = deval(ctx.SNEW[i + 1], ctx.W5, ctx.W6, ctx.SB, ctx.NB)

            cavp = math.sqrt(cv3**2 + cv2**2)
            if cavp == 0.0:
                cavp_s2 = 0.0
                cavp_s3 = 0.0
            else:
                cavp_s2 = cvs2 * cv2 / cavp
                cavp_s3 = cvs3 * cv3 / cavp

            fm = cc * cavm + 1.0
            fp = cc * cavp + 1.0

            rez = dsp * fp + dsm * fm

            ctx.W1[i] = -fm + cc * dsm * cavm_s1
            ctx.W2[i] = fp + fm + cc * (dsp * cavp_s2 + dsm * cavm_s2)
            ctx.W3[i] = -fp + cc * dsp * cavp_s3

            ctx.W4[i] = -rez

            cv1 = cv2
            cv2 = cv3
            cvs1 = cvs2
            cvs2 = cvs3
            cavm = cavp
            cavm_s1 = cavp_s2
            cavm_s2 = cavp_s3

        ctx.W2[1] = 1.0
        ctx.W3[1] = 0.0
        ctx.W4[1] = 0.0
        ctx.W1[nn] = 0.0
        ctx.W2[nn] = 1.0
        ctx.W4[nn] = 0.0

        if rtf != 1.0:
            i = 2
            ctx.W4[i] = -((ctx.SNEW[i] - ctx.SNEW[i - 1]) + rtf * (ctx.SNEW[i] - ctx.SNEW[i + 1]))
            ctx.W1[i] = -1.0
            ctx.W2[i] = 1.0 + rtf
            ctx.W3[i] = -rtf

            i = nn - 1
            ctx.W4[i] = -((ctx.SNEW[i] - ctx.SNEW[i + 1]) + rtf * (ctx.SNEW[i] - ctx.SNEW[i - 1]))
            ctx.W3[i] = -1.0
            ctx.W2[i] = 1.0 + rtf
            ctx.W1[i] = -rtf

        if ible != 0:
            i = nn1
            ctx.W1[i] = 0.0
            ctx.W2[i] = 1.0
            ctx.W3[i] = 0.0
            ctx.W4[i] = ctx.SBLE - ctx.SNEW[i]

        trisol(ctx.W2, ctx.W1, ctx.W3, ctx.W4, nn)

        rlx = 1.0
        dmax = 0.0
        for i in range(1, nn):
            ds = ctx.SNEW[i + 1] - ctx.SNEW[i]
            dds = ctx.W4[i + 1] - ctx.W4[i]
            dsrat = 1.0 + rlx * dds / ds
            if dsrat > 4.0:
                rlx = (4.0 - 1.0) * ds / dds
            if dsrat < 0.2:
                rlx = (0.2 - 1.0) * ds / dds
            dmax = max(abs(ctx.W4[i]), dmax)

        for i in range(2, nn):
            ctx.SNEW[i] = ctx.SNEW[i] + rlx * ctx.W4[i]

        if abs(dmax) < 1.0e-3:
            break
    else:
        print("Paneling convergence failed.  Continuing anyway...")

    for i in range(1, ctx.N + 1):
        ind = ipfac * (i - 1) + 1
        ctx.S[i] = ctx.SNEW[ind]
        ctx.X[i] = seval(ctx.SNEW[ind], ctx.XB, ctx.XBP, ctx.SB, ctx.NB)
        ctx.Y[i] = seval(ctx.SNEW[ind], ctx.YB, ctx.YBP, ctx.SB, ctx.NB)

    for ib in range(1, ctx.NB):
        if ctx.SB[ib] == ctx.SB[ib + 1]:
            xbcorn = ctx.XB[ib]
            ybcorn = ctx.YB[ib]
            sbcorn = ctx.SB[ib]

            for i in range(1, ctx.N + 1):
                if ctx.S[i] <= sbcorn:
                    continue

                for j in range(ctx.N, i - 1, -1):
                    ctx.X[j + 1] = ctx.X[j]
                    ctx.Y[j + 1] = ctx.Y[j]
                    ctx.S[j + 1] = ctx.S[j]
                ctx.N = ctx.N + 1

                if ctx.N > (len(ctx.X) - 1) - 1:
                    raise RuntimeError("PANEL: Too many panels. Increase IQX in XFOIL.INC")

                ctx.X[i] = xbcorn
                ctx.Y[i] = ybcorn
                ctx.S[i] = sbcorn

                if i - 2 >= 1:
                    ctx.S[i - 1] = 0.5 * (ctx.S[i] + ctx.S[i - 2])
                    ctx.X[i - 1] = seval(ctx.S[i - 1], ctx.XB, ctx.XBP, ctx.SB, ctx.NB)
                    ctx.Y[i - 1] = seval(ctx.S[i - 1], ctx.YB, ctx.YBP, ctx.SB, ctx.NB)

                if i + 2 <= ctx.N:
                    ctx.S[i + 1] = 0.5 * (ctx.S[i] + ctx.S[i + 2])
                    ctx.X[i + 1] = seval(ctx.S[i + 1], ctx.XB, ctx.XBP, ctx.SB, ctx.NB)
                    ctx.Y[i + 1] = seval(ctx.S[i + 1], ctx.YB, ctx.YBP, ctx.SB, ctx.NB)

                break

    scalc(ctx.X, ctx.Y, ctx.S, ctx.N)
    segspl(ctx.X, ctx.XP, ctx.S, ctx.N)
    segspl(ctx.Y, ctx.YP, ctx.S, ctx.N)
    ctx.SLE = lefind(ctx.X, ctx.XP, ctx.Y, ctx.YP, ctx.S, ctx.N)

    ctx.XLE = seval(ctx.SLE, ctx.X, ctx.XP, ctx.S, ctx.N)
    ctx.YLE = seval(ctx.SLE, ctx.Y, ctx.YP, ctx.S, ctx.N)
    ctx.XTE = 0.5 * (ctx.X[1] + ctx.X[ctx.N])
    ctx.YTE = 0.5 * (ctx.Y[1] + ctx.Y[ctx.N])
    ctx.CHORD = math.sqrt((ctx.XTE - ctx.XLE) ** 2 + (ctx.YTE - ctx.YLE) ** 2)

    dsmin = 1000.0
    dsmax = -1000.0
    for i in range(1, ctx.N):
        ds = ctx.S[i + 1] - ctx.S[i]
        if ds == 0.0:
            continue
        dsmin = min(dsmin, ds)
        dsmax = max(dsmax, ds)

    dsmin = dsmin * float(ctx.N - 1) / ctx.S[ctx.N]
    dsmax = dsmax * float(ctx.N - 1) / ctx.S[ctx.N]

    ctx.LGAMU = False
    ctx.LQINU = False
    ctx.LWAKE = False
    ctx.LQAIJ = False
    ctx.LADIJ = False
    ctx.LWDIJ = False
    ctx.LIPAN = False
    ctx.LBLINI = False
    ctx.LVCONV = False
    ctx.LSCINI = False
    ctx.LQSPEC = False
    ctx.LGSAME = False

    if ctx.LBFLAP:
        ctx.XOF = ctx.XBF
        ctx.YOF = ctx.YBF
        ctx.LFLAP = True

    tecalc(ctx)

    ncalc(ctx.X, ctx.Y, ctx.S, ctx.N, ctx.NX, ctx.NY)
    apcalc(ctx)

    if ctx.SHARP:
        print("\nSharp trailing edge")
    else:
        gap = math.sqrt((ctx.X[1] - ctx.X[ctx.N]) ** 2 + (ctx.Y[1] - ctx.Y[ctx.N]) ** 2)
        print(f"\nBlunt trailing edge.  Gap ={gap:9.5f}")

    if shopar:
        print(
            "\n Paneling parameters used..."
            f"\n   Number of panel nodes      {ctx.NPAN:4d}"
            f"\n   Panel bunching parameter   {ctx.CVPAR:6.3f}"
            f"\n   TE/LE panel density ratio  {ctx.CTERAT:6.3f}"
            f"\n   Refined-area/LE panel density ratio   {ctx.CTRRAT:6.3f}"
            f"\n   Top    side refined area x/c limits  {ctx.XSREF1:6.3f}{ctx.XSREF2:6.3f}"
            f"\n   Bottom side refined area x/c limits  {ctx.XPREF1:6.3f}{ctx.XPREF2:6.3f}"
        )
