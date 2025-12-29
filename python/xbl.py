# Ported from XFOIL Fortran source (Mark Drela).
# This file is a derived work and remains under the terms of the
# GNU General Public License v2 or later.
# See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

import math

from .xblcom import NCOM, COM1_NAMES, COM2_NAMES, sync_com_to_vars, sync_vars_to_com
from .xblsys import blprv, blkin, trchek, tesys, blsys, blvar, blmid, hkin
from .xfoil import mrcl, comset
from .xpanel import ueset
from .xsolve import gauss
from .spline import splind, sinvrt, seval

IQX = 370
IPX = 5
ISX = 2
IWX = IQX // 8 + 2
IBX = 4 * IQX
IZX = IQX + IWX
IVX = IQX // 2 + IWX + 50


def make_1d(n, fill=0.0):
    return [fill] * (n + 1)


def make_2d(n1, n2, fill=0.0):
    return [[fill] * (n2 + 1) for _ in range(n1 + 1)]


def make_3d(n1, n2, n3, fill=0.0):
    return [[[fill] * (n3 + 1) for _ in range(n2 + 1)] for _ in range(n1 + 1)]


class XFoilState:
    def __init__(self):
        self.LALFA = False
        self.LBLINI = False

        self.CL = 0.0
        self.CM = 0.0
        self.CD = 0.0
        self.CDP = 0.0
        self.CDF = 0.0
        self.CL_ALF = 0.0
        self.CL_MSQ = 0.0
        self.CLSPEC = 0.0
        self.MINF = 0.0
        self.MINF1 = 0.0
        self.MINF_CL = 0.0
        self.QINF = 1.0
        self.TKLAM = 0.0
        self.TKL_MSQ = 0.0
        self.GAMMA = 1.4
        self.GAMM1 = self.GAMMA - 1.0
        self.HVRAT = 0.0
        self.REINF = 0.0
        self.REINF1 = 0.0
        self.REINF_CL = 0.0
        self.RETYP = 1
        self.MATYP = 1
        self.CPSTAR = 0.0
        self.QSTAR = 0.0
        self.IDAMP = 0

        self.SLE = 0.0
        self.XLE = 0.0
        self.YLE = 0.0
        self.XTE = 0.0
        self.YTE = 0.0
        self.SST = 0.0
        self.SST_GO = 0.0
        self.SST_GP = 0.0

        self.ALFA = 0.0
        self.ADEG = 0.0
        self.DTOR = math.pi / 180.0
        self.AVISC = 0.0
        self.MVISC = 0.0
        self.XCMREF = 0.0
        self.YCMREF = 0.0

        self.N = 0
        self.IST = 0
        self.NW = 0
        self.NPAN = 0
        self.NB = 0
        self.LCLOCK = False

        self.X = make_1d(IZX)
        self.Y = make_1d(IZX)
        self.XP = make_1d(IZX)
        self.YP = make_1d(IZX)
        self.S = make_1d(IZX)
        self.SNEW = make_1d(5 * IBX)

        self.W1 = make_1d(6 * IQX)
        self.W2 = make_1d(6 * IQX)
        self.W3 = make_1d(6 * IQX)
        self.W4 = make_1d(6 * IQX)
        self.W5 = make_1d(6 * IQX)
        self.W6 = make_1d(6 * IQX)
        self.W7 = make_1d(6 * IQX)
        self.W8 = make_1d(6 * IQX)

        self.XB = make_1d(IBX)
        self.YB = make_1d(IBX)
        self.XBP = make_1d(IBX)
        self.YBP = make_1d(IBX)
        self.SB = make_1d(IBX)
        self.SBLE = 0.0
        self.CHORDB = 0.0
        self.AREAB = 0.0
        self.RADBLE = 0.0
        self.ANGBTE = 0.0
        self.EI11BA = 0.0
        self.EI22BA = 0.0
        self.APX1BA = 0.0
        self.APX2BA = 0.0
        self.EI11BT = 0.0
        self.EI22BT = 0.0
        self.APX1BT = 0.0
        self.APX2BT = 0.0
        self.THICKB = 0.0
        self.CAMBRB = 0.0

        self.XSSI = make_2d(IVX, ISX)
        self.UEDG = make_2d(IVX, ISX)
        self.UINV = make_2d(IVX, ISX)
        self.UINV_A = make_2d(IVX, ISX)
        self.MASS = make_2d(IVX, ISX)
        self.THET = make_2d(IVX, ISX)
        self.DSTR = make_2d(IVX, ISX)
        self.CTAU = make_2d(IVX, ISX)
        self.DELT = make_2d(IVX, ISX)
        self.TSTR = make_2d(IVX, ISX)
        self.USLP = make_2d(IVX, ISX)
        self.GUXQ = make_2d(IVX, ISX)
        self.GUXD = make_2d(IVX, ISX)
        self.TAU = make_2d(IVX, ISX)
        self.DIS = make_2d(IVX, ISX)
        self.CTQ = make_2d(IVX, ISX)
        self.VTI = make_2d(IVX, ISX)

        self.ACRIT = make_1d(ISX)
        self.XSTRIP = make_1d(ISX)
        self.XOCTR = make_1d(ISX)
        self.YOCTR = make_1d(ISX)
        self.XSSITR = make_1d(ISX)
        self.TINDEX = make_1d(ISX)

        self.IBLTE = [0] * (ISX + 1)
        self.NBL = [0] * (ISX + 1)
        self.IPAN = [[0] * (ISX + 1) for _ in range(IVX + 1)]
        self.ISYS = [[0] * (ISX + 1) for _ in range(IVX + 1)]
        self.NSYS = 0
        self.ITRAN = [0] * (ISX + 1)
        self.TFORCE = [False] * (ISX + 1)

        self.WGAP = make_1d(IWX)
        self.DWTE = 0.0
        self.ANTE = 0.0
        self.DSTE = 0.0
        self.ASTE = 0.0
        self.WAKLEN = 0.0
        self.CHORD = 0.0
        self.YIMAGE = 0.0
        self.SHARP = False

        self.XSTRIP[1] = 1.0
        self.XSTRIP[2] = 1.0
        self.XOCTR[1] = 1.0
        self.XOCTR[2] = 1.0
        self.YOCTR[1] = 0.0
        self.YOCTR[2] = 0.0

        self.DIJ = make_2d(IZX, IZX)
        self.AIJ = make_2d(IQX, IQX)
        self.BIJ = make_2d(IQX, IZX)
        self.CIJ = make_2d(IWX, IQX)
        self.VM = make_3d(3, IZX, IZX)
        self.VA = make_3d(3, 2, IZX)
        self.VB = make_3d(3, 2, IZX)
        self.VDEL = make_3d(3, 2, IZX)
        self.VZ = make_2d(3, 2)
        self.AIJPIV = [0] * (IQX + 1)

        self.QINV = make_1d(IZX)
        self.QVIS = make_1d(IZX)
        self.CPI = make_1d(IZX)
        self.CPV = make_1d(IZX)
        self.QINVU = make_2d(IZX, 2)
        self.QINV_A = make_1d(IZX)

        self.GAM = make_1d(IQX)
        self.GAMU = make_2d(IQX, 2)
        self.GAM_A = make_1d(IQX)
        self.SIG = make_1d(IZX)
        self.NX = make_1d(IZX)
        self.NY = make_1d(IZX)
        self.APANEL = make_1d(IZX)
        self.SST = 0.0
        self.SST_GO = 0.0
        self.SST_GP = 0.0
        self.GAMTE = 0.0
        self.GAMTE_A = 0.0
        self.SIGTE = 0.0
        self.SIGTE_A = 0.0

        self.DZDG = make_1d(IQX)
        self.DZDN = make_1d(IQX)
        self.DZDM = make_1d(IZX)
        self.DQDG = make_1d(IQX)
        self.DQDM = make_1d(IZX)
        self.QTAN1 = 0.0
        self.QTAN2 = 0.0
        self.Z_QINF = 0.0
        self.Z_ALFA = 0.0
        self.Z_QDOF0 = 0.0
        self.Z_QDOF1 = 0.0
        self.Z_QDOF2 = 0.0
        self.Z_QDOF3 = 0.0

        self.QF0 = make_1d(IQX)
        self.QF1 = make_1d(IQX)
        self.QF2 = make_1d(IQX)
        self.QF3 = make_1d(IQX)

        self.PI = math.pi
        self.HOPI = 0.5 / self.PI
        self.QOPI = 0.25 / self.PI

        self.RMSBL = 0.0
        self.RMXBL = 0.0
        self.RLX = 0.0
        self.VACCEL = 0.0
        self.IMXBL = 0
        self.ISMXBL = 0
        self.VMXBL = " "

        self.LQAIJ = False
        self.LADIJ = False
        self.LWDIJ = False
        self.LWAKE = False
        self.LGAMU = False
        self.LVISC = False
        self.LVCONV = False
        self.LFLAP = False
        self.LIMAGE = False
        self.LQINU = False
        self.LQSPEC = False
        self.LGSAME = False
        self.LSCINI = False
        self.AWAKE = 0.0

        self.XOF = 0.0
        self.YOF = 0.0
        self.HMOM = 0.0
        self.HFX = 0.0
        self.HFY = 0.0

        self.CVPAR = 0.0
        self.CTERAT = 0.0
        self.CTRRAT = 0.0
        self.XSREF1 = 0.0
        self.XSREF2 = 0.0
        self.XPREF1 = 0.0
        self.XPREF2 = 0.0

        self.XBF = 0.0
        self.YBF = 0.0
        self.LBFLAP = False

        self.NAME = ""
        self.NNAME = 0


class XBlState:
    def __init__(self):
        self.COM1 = make_1d(NCOM)
        self.COM2 = make_1d(NCOM)
        self.C1SAV = make_1d(NCOM)
        self.C2SAV = make_1d(NCOM)

        self.SIMI = False
        self.TRAN = False
        self.TURB = False
        self.WAKE = False
        self.TRFORC = False
        self.TRFREE = False
        self.IDAMPV = 0

        self.VS1 = make_2d(4, 5)
        self.VS2 = make_2d(4, 5)
        self.VSREZ = make_1d(4)
        self.VSR = make_1d(4)
        self.VSM = make_1d(4)
        self.VSX = make_1d(4)

        self.SCCON = 0.0
        self.GACON = 0.0
        self.GBCON = 0.0
        self.GCCON = 0.0
        self.DLCON = 0.0
        self.CTRCON = 0.0
        self.CTRCEX = 0.0
        self.DUXCON = 0.0
        self.CTCON = 0.0
        self.CFFAC = 1.0

        for name in COM1_NAMES[1:] + COM2_NAMES[1:]:
            setattr(self, name, 0.0)

        for name in [
            "CFM",
            "CFM_MS",
            "CFM_RE",
            "CFM_U1",
            "CFM_T1",
            "CFM_D1",
            "CFM_U2",
            "CFM_T2",
            "CFM_D2",
            "XT",
            "XT_A1",
            "XT_MS",
            "XT_RE",
            "XT_XF",
            "XT_X1",
            "XT_T1",
            "XT_D1",
            "XT_U1",
            "XT_X2",
            "XT_T2",
            "XT_D2",
            "XT_U2",
            "DWTE",
            "QINFBL",
            "TKBL",
            "TKBL_MS",
            "RSTBL",
            "RSTBL_MS",
            "HSTINV",
            "HSTINV_MS",
            "REYBL",
            "REYBL_MS",
            "REYBL_RE",
            "GAMBL",
            "GM1BL",
            "HVRAT",
            "BULE",
            "XIFORC",
            "AMCRIT",
        ]:
            setattr(self, name, 0.0)


def setbl(ctx, bl):
    # No tracing; setbl assembles VA/VB/VDEL/VM directly.

    if ctx.LALFA:
        clmr = ctx.CL
    else:
        clmr = ctx.CLSPEC

    ma_clmr, re_clmr = mrcl(ctx, clmr)
    msq_clmr = 2.0 * ctx.MINF * ma_clmr

    comset(ctx)

    bl.GAMBL = ctx.GAMMA
    bl.GM1BL = ctx.GAMM1

    bl.QINFBL = ctx.QINF
    bl.TKBL = ctx.TKLAM
    bl.TKBL_MS = ctx.TKL_MSQ

    bl.RSTBL = (1.0 + 0.5 * bl.GM1BL * ctx.MINF**2) ** (1.0 / bl.GM1BL)
    bl.RSTBL_MS = 0.5 * bl.RSTBL / (1.0 + 0.5 * bl.GM1BL * ctx.MINF**2)

    bl.HSTINV = (
        bl.GM1BL * (ctx.MINF / bl.QINFBL) ** 2 / (1.0 + 0.5 * bl.GM1BL * ctx.MINF**2)
    )
    bl.HSTINV_MS = (
        bl.GM1BL * (1.0 / bl.QINFBL) ** 2 / (1.0 + 0.5 * bl.GM1BL * ctx.MINF**2)
        - 0.5 * bl.GM1BL * bl.HSTINV / (1.0 + 0.5 * bl.GM1BL * ctx.MINF**2)
    )

    herat = 1.0 - 0.5 * bl.QINFBL**2 * bl.HSTINV
    herat_ms = -0.5 * bl.QINFBL**2 * bl.HSTINV_MS

    bl.HVRAT = ctx.HVRAT
    bl.REYBL = ctx.REINF * math.sqrt(herat**3) * (1.0 + bl.HVRAT) / (herat + bl.HVRAT)
    bl.REYBL_RE = math.sqrt(herat**3) * (1.0 + bl.HVRAT) / (herat + bl.HVRAT)
    bl.REYBL_MS = bl.REYBL * (1.5 / herat - 1.0 / (herat + bl.HVRAT)) * herat_ms

    bl.IDAMPV = ctx.IDAMP
    bl.DWTE = ctx.WGAP[1]


    if not ctx.LBLINI:
        print()
        print("Initializing BL ...")
        mrchue(ctx, bl)
        ctx.LBLINI = True

    print()

    mrchdu(ctx, bl)

    usav = make_2d(IVX, 2)
    for is_ in range(1, 3):
        for ibl in range(2, ctx.NBL[is_] + 1):
            usav[ibl][is_] = ctx.UEDG[ibl][is_]

    ueset(ctx)

    for is_ in range(1, 3):
        for ibl in range(2, ctx.NBL[is_] + 1):
            temp = usav[ibl][is_]
            usav[ibl][is_] = ctx.UEDG[ibl][is_]
            ctx.UEDG[ibl][is_] = temp

    ile1 = ctx.IPAN[2][1]
    ile2 = ctx.IPAN[2][2]
    ite1 = ctx.IPAN[ctx.IBLTE[1]][1]
    ite2 = ctx.IPAN[ctx.IBLTE[2]][2]

    jvte1 = ctx.ISYS[ctx.IBLTE[1]][1]
    jvte2 = ctx.ISYS[ctx.IBLTE[2]][2]

    dule1 = ctx.UEDG[2][1] - usav[2][1]
    dule2 = ctx.UEDG[2][2] - usav[2][2]

    ule1_m = make_1d(2 * IVX)
    ule2_m = make_1d(2 * IVX)
    ute1_m = make_1d(2 * IVX)
    ute2_m = make_1d(2 * IVX)

    for js in range(1, 3):
        for jbl in range(2, ctx.NBL[js] + 1):
            j = ctx.IPAN[jbl][js]
            jv = ctx.ISYS[jbl][js]
            ule1_m[jv] = -ctx.VTI[2][1] * ctx.VTI[jbl][js] * ctx.DIJ[ile1][j]
            ule2_m[jv] = -ctx.VTI[2][2] * ctx.VTI[jbl][js] * ctx.DIJ[ile2][j]
            ute1_m[jv] = -ctx.VTI[ctx.IBLTE[1]][1] * ctx.VTI[jbl][js] * ctx.DIJ[ite1][j]
            ute2_m[jv] = -ctx.VTI[ctx.IBLTE[2]][2] * ctx.VTI[jbl][js] * ctx.DIJ[ite2][j]

    ule1_a = ctx.UINV_A[2][1]
    ule2_a = ctx.UINV_A[2][2]

    ctx.TINDEX[1] = 0.0
    ctx.TINDEX[2] = 0.0

    u1_m = make_1d(2 * IVX)
    u2_m = make_1d(2 * IVX)
    d1_m = make_1d(2 * IVX)
    d2_m = make_1d(2 * IVX)

    for is_ in range(1, 3):
        for js in range(1, 3):
            for jbl in range(2, ctx.NBL[js] + 1):
                jv = ctx.ISYS[jbl][js]
                u1_m[jv] = 0.0
                d1_m[jv] = 0.0
        u1_a = 0.0
        d1_a = 0.0

        due1 = 0.0
        dds1 = 0.0

        ibl = 2
        bl.BULE = 1.0
        bl.AMCRIT = ctx.ACRIT[is_]

        xifset(ctx, bl, is_)

        bl.TRAN = False
        bl.TURB = False

        cti = 0.0
        for ibl in range(2, ctx.NBL[is_] + 1):
            iv = ctx.ISYS[ibl][is_]

            bl.SIMI = ibl == 2
            bl.WAKE = ibl > ctx.IBLTE[is_]
            bl.TRAN = ibl == ctx.ITRAN[is_]
            bl.TURB = ibl > ctx.ITRAN[is_]

            i = ctx.IPAN[ibl][is_]

            xsi = ctx.XSSI[ibl][is_]
            if ibl < ctx.ITRAN[is_]:
                ami = ctx.CTAU[ibl][is_]
            if ibl >= ctx.ITRAN[is_]:
                cti = ctx.CTAU[ibl][is_]
            uei = ctx.UEDG[ibl][is_]
            thi = ctx.THET[ibl][is_]
            mdi = ctx.MASS[ibl][is_]

            dsi = mdi / uei

            if bl.WAKE:
                iw = ibl - ctx.IBLTE[is_]
                dswaki = ctx.WGAP[iw]
            else:
                dswaki = 0.0

            d2_m2 = 1.0 / uei
            d2_u2 = -dsi / uei

            for js in range(1, 3):
                for jbl in range(2, ctx.NBL[js] + 1):
                    j = ctx.IPAN[jbl][js]
                    jv = ctx.ISYS[jbl][js]
                    u2_m[jv] = -ctx.VTI[ibl][is_] * ctx.VTI[jbl][js] * ctx.DIJ[i][j]
                    d2_m[jv] = d2_u2 * u2_m[jv]
            d2_m[iv] = d2_m[iv] + d2_m2

            u2_a = ctx.UINV_A[ibl][is_]
            d2_a = d2_u2 * u2_a

            due2 = ctx.UEDG[ibl][is_] - usav[ibl][is_]
            dds2 = d2_u2 * due2

            blprv(bl, xsi, ami, cti, thi, dsi, dswaki, uei)
            blkin(bl)

            if bl.TRAN:
                trchek(bl)
                ami = bl.AMPL2
            if ibl == ctx.ITRAN[is_] and not bl.TRAN:
                print("SETBL: Xtr???  n1 n2: ", bl.AMPL1, bl.AMPL2)

            if ibl == ctx.IBLTE[is_] + 1:
                tte = ctx.THET[ctx.IBLTE[1]][1] + ctx.THET[ctx.IBLTE[2]][2]
                dte = ctx.DSTR[ctx.IBLTE[1]][1] + ctx.DSTR[ctx.IBLTE[2]][2] + ctx.ANTE
                cte = (
                    ctx.CTAU[ctx.IBLTE[1]][1] * ctx.THET[ctx.IBLTE[1]][1]
                    + ctx.CTAU[ctx.IBLTE[2]][2] * ctx.THET[ctx.IBLTE[2]][2]
                ) / tte
                tesys(bl, cte, tte, dte)

                tte_tte1 = 1.0
                tte_tte2 = 1.0
                dte_mte1 = 1.0 / ctx.UEDG[ctx.IBLTE[1]][1]
                dte_ute1 = -ctx.DSTR[ctx.IBLTE[1]][1] / ctx.UEDG[ctx.IBLTE[1]][1]
                dte_mte2 = 1.0 / ctx.UEDG[ctx.IBLTE[2]][2]
                dte_ute2 = -ctx.DSTR[ctx.IBLTE[2]][2] / ctx.UEDG[ctx.IBLTE[2]][2]
                cte_cte1 = ctx.THET[ctx.IBLTE[1]][1] / tte
                cte_cte2 = ctx.THET[ctx.IBLTE[2]][2] / tte
                cte_tte1 = (ctx.CTAU[ctx.IBLTE[1]][1] - cte) / tte
                cte_tte2 = (ctx.CTAU[ctx.IBLTE[2]][2] - cte) / tte

                for js in range(1, 3):
                    for jbl in range(2, ctx.NBL[js] + 1):
                        j = ctx.IPAN[jbl][js]
                        jv = ctx.ISYS[jbl][js]
                        d1_m[jv] = dte_ute1 * ute1_m[jv] + dte_ute2 * ute2_m[jv]
                d1_m[jvte1] = d1_m[jvte1] + dte_mte1
                d1_m[jvte2] = d1_m[jvte2] + dte_mte2

                due1 = 0.0
                dds1 = dte_ute1 * (ctx.UEDG[ctx.IBLTE[1]][1] - usav[ctx.IBLTE[1]][1]) + dte_ute2 * (
                    ctx.UEDG[ctx.IBLTE[2]][2] - usav[ctx.IBLTE[2]][2]
                )
            else:
                blsys(bl)

            ctx.TAU[ibl][is_] = 0.5 * bl.R2 * bl.U2 * bl.U2 * bl.CF2
            ctx.DIS[ibl][is_] = bl.R2 * bl.U2 * bl.U2 * bl.U2 * bl.DI2 * bl.HS2 * 0.5
            ctx.CTQ[ibl][is_] = bl.CQ2
            ctx.DELT[ibl][is_] = bl.DE2
            ctx.USLP[ibl][is_] = 1.60 / (1.0 + bl.US2)

            if is_ == 1:
                xi_ule1 = ctx.SST_GO
                xi_ule2 = -ctx.SST_GP
            else:
                xi_ule1 = -ctx.SST_GO
                xi_ule2 = ctx.SST_GP

            for jv in range(1, ctx.NSYS + 1):
                ctx.VM[1][jv][iv] = (
                    bl.VS1[1][3] * d1_m[jv]
                    + bl.VS1[1][4] * u1_m[jv]
                    + bl.VS2[1][3] * d2_m[jv]
                    + bl.VS2[1][4] * u2_m[jv]
                    + (bl.VS1[1][5] + bl.VS2[1][5] + bl.VSX[1]) * (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv])
                )
            ctx.VB[1][1][iv] = bl.VS1[1][1]
            ctx.VB[1][2][iv] = bl.VS1[1][2]

            ctx.VA[1][1][iv] = bl.VS2[1][1]
            ctx.VA[1][2][iv] = bl.VS2[1][2]

            if ctx.LALFA:
                ctx.VDEL[1][2][iv] = bl.VSR[1] * re_clmr + bl.VSM[1] * msq_clmr
            else:
                ctx.VDEL[1][2][iv] = (
                    bl.VS1[1][4] * u1_a
                    + bl.VS1[1][3] * d1_a
                    + bl.VS2[1][4] * u2_a
                    + bl.VS2[1][3] * d2_a
                    + (bl.VS1[1][5] + bl.VS2[1][5] + bl.VSX[1]) * (xi_ule1 * ule1_a + xi_ule2 * ule2_a)
                )

            ctx.VDEL[1][1][iv] = (
                bl.VSREZ[1]
                + (bl.VS1[1][4] * due1 + bl.VS1[1][3] * dds1)
                + (bl.VS2[1][4] * due2 + bl.VS2[1][3] * dds2)
                + (bl.VS1[1][5] + bl.VS2[1][5] + bl.VSX[1]) * (xi_ule1 * dule1 + xi_ule2 * dule2)
            )

            for jv in range(1, ctx.NSYS + 1):
                ctx.VM[2][jv][iv] = (
                    bl.VS1[2][3] * d1_m[jv]
                    + bl.VS1[2][4] * u1_m[jv]
                    + bl.VS2[2][3] * d2_m[jv]
                    + bl.VS2[2][4] * u2_m[jv]
                    + (bl.VS1[2][5] + bl.VS2[2][5] + bl.VSX[2]) * (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv])
                )
            ctx.VB[2][1][iv] = bl.VS1[2][1]
            ctx.VB[2][2][iv] = bl.VS1[2][2]
            ctx.VA[2][1][iv] = bl.VS2[2][1]
            ctx.VA[2][2][iv] = bl.VS2[2][2]

            if ctx.LALFA:
                ctx.VDEL[2][2][iv] = bl.VSR[2] * re_clmr + bl.VSM[2] * msq_clmr
            else:
                ctx.VDEL[2][2][iv] = (
                    bl.VS1[2][4] * u1_a
                    + bl.VS1[2][3] * d1_a
                    + bl.VS2[2][4] * u2_a
                    + bl.VS2[2][3] * d2_a
                    + (bl.VS1[2][5] + bl.VS2[2][5] + bl.VSX[2]) * (xi_ule1 * ule1_a + xi_ule2 * ule2_a)
                )

            ctx.VDEL[2][1][iv] = (
                bl.VSREZ[2]
                + (bl.VS1[2][4] * due1 + bl.VS1[2][3] * dds1)
                + (bl.VS2[2][4] * due2 + bl.VS2[2][3] * dds2)
                + (bl.VS1[2][5] + bl.VS2[2][5] + bl.VSX[2]) * (xi_ule1 * dule1 + xi_ule2 * dule2)
            )

            for jv in range(1, ctx.NSYS + 1):
                ctx.VM[3][jv][iv] = (
                    bl.VS1[3][3] * d1_m[jv]
                    + bl.VS1[3][4] * u1_m[jv]
                    + bl.VS2[3][3] * d2_m[jv]
                    + bl.VS2[3][4] * u2_m[jv]
                    + (bl.VS1[3][5] + bl.VS2[3][5] + bl.VSX[3]) * (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv])
                )
            ctx.VB[3][1][iv] = bl.VS1[3][1]
            ctx.VB[3][2][iv] = bl.VS1[3][2]
            ctx.VA[3][1][iv] = bl.VS2[3][1]
            ctx.VA[3][2][iv] = bl.VS2[3][2]

            if ctx.LALFA:
                ctx.VDEL[3][2][iv] = bl.VSR[3] * re_clmr + bl.VSM[3] * msq_clmr
            else:
                ctx.VDEL[3][2][iv] = (
                    bl.VS1[3][4] * u1_a
                    + bl.VS1[3][3] * d1_a
                    + bl.VS2[3][4] * u2_a
                    + bl.VS2[3][3] * d2_a
                    + (bl.VS1[3][5] + bl.VS2[3][5] + bl.VSX[3]) * (xi_ule1 * ule1_a + xi_ule2 * ule2_a)
                )

            ctx.VDEL[3][1][iv] = (
                bl.VSREZ[3]
                + (bl.VS1[3][4] * due1 + bl.VS1[3][3] * dds1)
                + (bl.VS2[3][4] * due2 + bl.VS2[3][3] * dds2)
                + (bl.VS1[3][5] + bl.VS2[3][5] + bl.VSX[3]) * (xi_ule1 * dule1 + xi_ule2 * dule2)
            )

            if ibl == ctx.IBLTE[is_] + 1:
                ctx.VZ[1][1] = bl.VS1[1][1] * cte_cte1
                ctx.VZ[1][2] = bl.VS1[1][1] * cte_tte1 + bl.VS1[1][2] * tte_tte1
                ctx.VB[1][1][iv] = bl.VS1[1][1] * cte_cte2
                ctx.VB[1][2][iv] = bl.VS1[1][1] * cte_tte2 + bl.VS1[1][2] * tte_tte2

                ctx.VZ[2][1] = bl.VS1[2][1] * cte_cte1
                ctx.VZ[2][2] = bl.VS1[2][1] * cte_tte1 + bl.VS1[2][2] * tte_tte1
                ctx.VB[2][1][iv] = bl.VS1[2][1] * cte_cte2
                ctx.VB[2][2][iv] = bl.VS1[2][1] * cte_tte2 + bl.VS1[2][2] * tte_tte2

                ctx.VZ[3][1] = bl.VS1[3][1] * cte_cte1
                ctx.VZ[3][2] = bl.VS1[3][1] * cte_tte1 + bl.VS1[3][2] * tte_tte1
                ctx.VB[3][1][iv] = bl.VS1[3][1] * cte_cte2
                ctx.VB[3][2][iv] = bl.VS1[3][1] * cte_tte2 + bl.VS1[3][2] * tte_tte2

            if bl.TRAN:
                bl.TURB = True
                ctx.ITRAN[is_] = ibl
                ctx.TFORCE[is_] = bl.TRFORC
                ctx.XSSITR[is_] = bl.XT

                if is_ == 1:
                    str_ = ctx.SST - bl.XT
                else:
                    str_ = ctx.SST + bl.XT
                chx = ctx.XTE - ctx.XLE
                chy = ctx.YTE - ctx.YLE
                chsq = chx**2 + chy**2
                xtr = seval(str_, ctx.X, ctx.XP, ctx.S, ctx.N)
                ytr = seval(str_, ctx.Y, ctx.YP, ctx.S, ctx.N)
                ctx.XOCTR[is_] = ((xtr - ctx.XLE) * chx + (ytr - ctx.YLE) * chy) / chsq
                ctx.YOCTR[is_] = ((ytr - ctx.YLE) * chx - (xtr - ctx.XLE) * chy) / chsq

            bl.TRAN = False

            if ibl == ctx.IBLTE[is_]:
                bl.TURB = True
                bl.WAKE = True
                blvar(bl, 3)
                blmid(bl, 3)

            for js in range(1, 3):
                for jbl in range(2, ctx.NBL[js] + 1):
                    jv = ctx.ISYS[jbl][js]
                    u1_m[jv] = u2_m[jv]
                    d1_m[jv] = d2_m[jv]

            u1_a = u2_a
            d1_a = d2_a
            due1 = due2
            dds1 = dds2

            if ibl == ctx.ITRAN[is_] and bl.X2 > bl.X1:
                if is_ == 1:
                    ctx.TINDEX[is_] = float(ctx.IST - ctx.ITRAN[is_] + 3) - (bl.XT - bl.X1) / (bl.X2 - bl.X1)
                else:
                    ctx.TINDEX[is_] = float(ctx.IST + ctx.ITRAN[is_] - 2) + (bl.XT - bl.X1) / (bl.X2 - bl.X1)


            for icom in range(1, NCOM + 1):
                bl.COM1[icom] = bl.COM2[icom]
            sync_com_to_vars(bl, 1)

        if ctx.TFORCE[is_]:
            print(f" Side{is_:2d} forced transition at x/c = {ctx.XOCTR[is_]:7.4f}{ctx.ITRAN[is_]:5d}")
        else:
            print(f" Side{is_:2d}  free  transition at x/c = {ctx.XOCTR[is_]:7.4f}{ctx.ITRAN[is_]:5d}")


def iblsys(ctx):
    iv = 0
    for is_ in range(1, 3):
        for ibl in range(2, ctx.NBL[is_] + 1):
            iv += 1
            ctx.ISYS[ibl][is_] = iv

    ctx.NSYS = iv
    if ctx.NSYS > 2 * IVX:
        raise RuntimeError("*** IBLSYS: BL system array overflow. ***")


def mrchue(ctx, bl):
    hlmax = 3.8
    htmax = 2.5
    max_nbl = max(ctx.NBL[1], ctx.NBL[2])
    if not hasattr(ctx, "HTARG") or len(ctx.HTARG) < 3 or len(ctx.HTARG[1]) <= max_nbl:
        ctx.HTARG = [[0.0] * (max_nbl + 1) for _ in range(3)]
    else:
        for is_ in range(1, 3):
            for ibl in range(1, max_nbl + 1):
                ctx.HTARG[is_][ibl] = 0.0

    for is_ in range(1, 3):
        print("   side ", is_, " ...")

        bl.AMCRIT = ctx.ACRIT[is_]

        xifset(ctx, bl, is_)

        ibl = 2
        xsi = ctx.XSSI[ibl][is_]
        uei = ctx.UEDG[ibl][is_]
        bl.BULE = 1.0
        ucon = uei / xsi**bl.BULE
        tsq = 0.45 / (ucon * (5.0 * bl.BULE + 1.0) * bl.REYBL) * xsi ** (1.0 - bl.BULE)
        thi = math.sqrt(tsq)
        dsi = 2.2 * thi
        ami = 0.0

        cti = 0.03

        bl.TRAN = False
        bl.TURB = False
        ctx.ITRAN[is_] = ctx.IBLTE[is_]

        sens = 0.0
        for ibl in range(2, ctx.NBL[is_] + 1):
            ibm = ibl - 1

            bl.SIMI = ibl == 2
            bl.WAKE = ibl > ctx.IBLTE[is_]

            xsi = ctx.XSSI[ibl][is_]
            uei = ctx.UEDG[ibl][is_]

            if bl.WAKE:
                iw = ibl - ctx.IBLTE[is_]
                dswaki = ctx.WGAP[iw]
            else:
                dswaki = 0.0

            direct = True

            for itbl in range(1, 26):
                sennew = sens
                blprv(bl, xsi, ami, cti, thi, dsi, dswaki, uei)
                blkin(bl)

                if (not bl.SIMI) and (not bl.TURB):
                    trchek(bl)
                    ami = bl.AMPL2

                    if bl.TRAN:
                        ctx.ITRAN[is_] = ibl
                        if cti <= 0.0:
                            cti = 0.03
                            bl.S2 = cti
                    else:
                        ctx.ITRAN[is_] = ibl + 2

                if ibl == ctx.IBLTE[is_] + 1:
                    tte = ctx.THET[ctx.IBLTE[1]][1] + ctx.THET[ctx.IBLTE[2]][2]
                    dte = ctx.DSTR[ctx.IBLTE[1]][1] + ctx.DSTR[ctx.IBLTE[2]][2] + ctx.ANTE
                    cte = (
                        ctx.CTAU[ctx.IBLTE[1]][1] * ctx.THET[ctx.IBLTE[1]][1]
                        + ctx.CTAU[ctx.IBLTE[2]][2] * ctx.THET[ctx.IBLTE[2]][2]
                    ) / tte
                    tesys(bl, cte, tte, dte)
                else:
                    blsys(bl)

                if direct:
                    bl.VS2[4][1] = 0.0
                    bl.VS2[4][2] = 0.0
                    bl.VS2[4][3] = 0.0
                    bl.VS2[4][4] = 1.0
                    bl.VSREZ[4] = 0.0

                    gauss(4, 4, bl.VS2, bl.VSREZ, 1)

                    dmax = max(abs(bl.VSREZ[2] / thi), abs(bl.VSREZ[3] / dsi))
                    if ibl < ctx.ITRAN[is_]:
                        dmax = max(dmax, abs(bl.VSREZ[1] / 10.0))
                    if ibl >= ctx.ITRAN[is_]:
                        dmax = max(dmax, abs(bl.VSREZ[1] / cti))

                    rlx = 1.0
                    if dmax > 0.3:
                        rlx = 0.3 / dmax

                    if ibl != ctx.IBLTE[is_] + 1:
                        msq = uei * uei * bl.HSTINV / (
                            bl.GM1BL * (1.0 - 0.5 * uei * uei * bl.HSTINV)
                        )
                        htest = (dsi + rlx * bl.VSREZ[3]) / (thi + rlx * bl.VSREZ[2])
                        hktest, _, _ = hkin(htest, msq)

                        if ibl < ctx.ITRAN[is_]:
                            hmax = hlmax
                        if ibl >= ctx.ITRAN[is_]:
                            hmax = htmax
                        direct = hktest < hmax

                    if direct:
                        if ibl >= ctx.ITRAN[is_]:
                            cti = cti + rlx * bl.VSREZ[1]
                        thi = thi + rlx * bl.VSREZ[2]
                        dsi = dsi + rlx * bl.VSREZ[3]
                    else:
                        if ibl < ctx.ITRAN[is_]:
                            htarg = bl.HK1 + 0.03 * (bl.X2 - bl.X1) / bl.T1
                        elif ibl == ctx.ITRAN[is_]:
                            htarg = bl.HK1 + (0.03 * (bl.XT - bl.X1) - 0.15 * (bl.X2 - bl.XT)) / bl.T1
                        elif bl.WAKE:
                            const = 0.03 * (bl.X2 - bl.X1) / bl.T1
                            hk2 = bl.HK1
                            hk2 = hk2 - (hk2 + const * (hk2 - 1.0) ** 3 - bl.HK1) / (
                                1.0 + 3.0 * const * (hk2 - 1.0) ** 2
                            )
                            hk2 = hk2 - (hk2 + const * (hk2 - 1.0) ** 3 - bl.HK1) / (
                                1.0 + 3.0 * const * (hk2 - 1.0) ** 2
                            )
                            hk2 = hk2 - (hk2 + const * (hk2 - 1.0) ** 3 - bl.HK1) / (
                                1.0 + 3.0 * const * (hk2 - 1.0) ** 2
                            )
                            htarg = hk2
                        else:
                            htarg = bl.HK1 - 0.15 * (bl.X2 - bl.X1) / bl.T1

                        if bl.WAKE:
                            htarg = max(htarg, 1.01)
                        else:
                            htarg = max(htarg, hmax)

                        print(f" MRCHUE: Inverse mode at{ibl:4d}     Hk ={htarg:8.3f}")
                        ctx.HTARG[is_][ibl] = htarg

                        continue
                else:
                    bl.VS2[4][1] = 0.0
                    bl.VS2[4][2] = bl.HK2_T2
                    bl.VS2[4][3] = bl.HK2_D2
                    bl.VS2[4][4] = bl.HK2_U2
                    bl.VSREZ[4] = htarg - bl.HK2

                    gauss(4, 4, bl.VS2, bl.VSREZ, 1)

                    dmax = max(abs(bl.VSREZ[2] / thi), abs(bl.VSREZ[3] / dsi), abs(bl.VSREZ[4] / uei))
                    if ibl >= ctx.ITRAN[is_]:
                        dmax = max(dmax, abs(bl.VSREZ[1] / cti))

                    rlx = 1.0
                    if dmax > 0.3:
                        rlx = 0.3 / dmax

                    if ibl >= ctx.ITRAN[is_]:
                        cti = cti + rlx * bl.VSREZ[1]
                    thi = thi + rlx * bl.VSREZ[2]
                    dsi = dsi + rlx * bl.VSREZ[3]
                    uei = uei + rlx * bl.VSREZ[4]

                if ibl >= ctx.ITRAN[is_]:
                    cti = min(cti, 0.30)
                    cti = max(cti, 0.0000001)

                if ibl <= ctx.IBLTE[is_]:
                    hklim = 1.02
                else:
                    hklim = 1.00005
                msq = uei * uei * bl.HSTINV / (bl.GM1BL * (1.0 - 0.5 * uei * uei * bl.HSTINV))
                dsw = dsi - dswaki
                dsw = dslim(dsw, thi, uei, msq, hklim)
                dsi = dsw + dswaki

                if dmax <= 1.0e-5:
                    break
            else:
                print(f" MRCHUE: Convergence failed at{ibl:4d}  side{is_:2d}    Res ={dmax:12.4E}")
                if dmax <= 0.1:
                    pass
                else:
                    if ibl > 3:
                        if ibl <= ctx.IBLTE[is_]:
                            thi = ctx.THET[ibm][is_] * (ctx.XSSI[ibl][is_] / ctx.XSSI[ibm][is_]) ** 0.5
                            dsi = ctx.DSTR[ibm][is_] * (ctx.XSSI[ibl][is_] / ctx.XSSI[ibm][is_]) ** 0.5
                        elif ibl == ctx.IBLTE[is_] + 1:
                            cti = cte
                            thi = tte
                            dsi = dte
                        else:
                            thi = ctx.THET[ibm][is_]
                            ratlen = (ctx.XSSI[ibl][is_] - ctx.XSSI[ibm][is_]) / (10.0 * ctx.DSTR[ibm][is_])
                            dsi = (ctx.DSTR[ibm][is_] + thi * ratlen) / (1.0 + ratlen)
                        if ibl == ctx.ITRAN[is_]:
                            cti = 0.05
                        if ibl > ctx.ITRAN[is_]:
                            cti = ctx.CTAU[ibm][is_]

                        uei = ctx.UEDG[ibl][is_]
                        if ibl > 2 and ibl < ctx.NBL[is_]:
                            uei = 0.5 * (ctx.UEDG[ibl - 1][is_] + ctx.UEDG[ibl + 1][is_])

                blprv(bl, xsi, ami, cti, thi, dsi, dswaki, uei)
                blkin(bl)

                if (not bl.SIMI) and (not bl.TURB):
                    trchek(bl)
                    ami = bl.AMPL2
                    if bl.TRAN:
                        ctx.ITRAN[is_] = ibl
                    if not bl.TRAN:
                        ctx.ITRAN[is_] = ibl + 2

                if ibl < ctx.ITRAN[is_]:
                    blvar(bl, 1)
                if ibl >= ctx.ITRAN[is_]:
                    blvar(bl, 2)
                if bl.WAKE:
                    blvar(bl, 3)

                if ibl < ctx.ITRAN[is_]:
                    blmid(bl, 1)
                if ibl >= ctx.ITRAN[is_]:
                    blmid(bl, 2)
                if bl.WAKE:
                    blmid(bl, 3)

            if ibl < ctx.ITRAN[is_]:
                ctx.CTAU[ibl][is_] = ami
            if ibl >= ctx.ITRAN[is_]:
                ctx.CTAU[ibl][is_] = cti
            ctx.THET[ibl][is_] = thi
            ctx.DSTR[ibl][is_] = dsi
            ctx.UEDG[ibl][is_] = uei
            ctx.MASS[ibl][is_] = dsi * uei
            ctx.TAU[ibl][is_] = 0.5 * bl.R2 * bl.U2 * bl.U2 * bl.CF2
            ctx.DIS[ibl][is_] = bl.R2 * bl.U2 * bl.U2 * bl.U2 * bl.DI2 * bl.HS2 * 0.5
            ctx.CTQ[ibl][is_] = bl.CQ2
            ctx.DELT[ibl][is_] = bl.DE2
            ctx.TSTR[ibl][is_] = bl.HS2 * bl.T2

            blprv(bl, xsi, ami, cti, thi, dsi, dswaki, uei)
            blkin(bl)
            for icom in range(1, NCOM + 1):
                bl.COM1[icom] = bl.COM2[icom]
            sync_com_to_vars(bl, 1)

            if bl.TRAN or ibl == ctx.IBLTE[is_]:
                bl.TURB = True
                ctx.TFORCE[is_] = bl.TRFORC
                ctx.XSSITR[is_] = bl.XT

            bl.TRAN = False

            if ibl == ctx.IBLTE[is_]:
                thi = ctx.THET[ctx.IBLTE[1]][1] + ctx.THET[ctx.IBLTE[2]][2]
                dsi = ctx.DSTR[ctx.IBLTE[1]][1] + ctx.DSTR[ctx.IBLTE[2]][2] + ctx.ANTE


def mrchdu(ctx, bl):
    deps = 5.0e-6
    senswt = 1000.0

    for is_ in range(1, 3):
        bl.AMCRIT = ctx.ACRIT[is_]

        xifset(ctx, bl, is_)

        ibl = 2
        xsi = ctx.XSSI[ibl][is_]
        uei = ctx.UEDG[ibl][is_]
        bl.BULE = 1.0

        itrold = ctx.ITRAN[is_]

        bl.TRAN = False
        bl.TURB = False
        ctx.ITRAN[is_] = ctx.IBLTE[is_]

        sens = 0.0
        for ibl in range(2, ctx.NBL[is_] + 1):
            ibm = ibl - 1

            bl.SIMI = ibl == 2
            bl.WAKE = ibl > ctx.IBLTE[is_]

            xsi = ctx.XSSI[ibl][is_]
            uei = ctx.UEDG[ibl][is_]
            thi = ctx.THET[ibl][is_]
            dsi = ctx.DSTR[ibl][is_]

            if ibl < itrold:
                ami = ctx.CTAU[ibl][is_]
                cti = 0.03
            else:
                cti = ctx.CTAU[ibl][is_]
                if cti <= 0.0:
                    cti = 0.03

            if bl.WAKE:
                iw = ibl - ctx.IBLTE[is_]
                dswaki = ctx.WGAP[iw]
            else:
                dswaki = 0.0

            if ibl <= ctx.IBLTE[is_]:
                dsi = max(dsi - dswaki, 1.02000 * thi) + dswaki
            if ibl > ctx.IBLTE[is_]:
                dsi = max(dsi - dswaki, 1.00005 * thi) + dswaki

            sennew = sens
            for itbl in range(1, 26):
                blprv(bl, xsi, ami, cti, thi, dsi, dswaki, uei)
                blkin(bl)

                if (not bl.SIMI) and (not bl.TURB):
                    trchek(bl)
                    ami = bl.AMPL2
                    if bl.TRAN:
                        ctx.ITRAN[is_] = ibl
                    if not bl.TRAN:
                        ctx.ITRAN[is_] = ibl + 2

                if ibl == ctx.IBLTE[is_] + 1:
                    tte = ctx.THET[ctx.IBLTE[1]][1] + ctx.THET[ctx.IBLTE[2]][2]
                    dte = ctx.DSTR[ctx.IBLTE[1]][1] + ctx.DSTR[ctx.IBLTE[2]][2] + ctx.ANTE
                    cte = (
                        ctx.CTAU[ctx.IBLTE[1]][1] * ctx.THET[ctx.IBLTE[1]][1]
                        + ctx.CTAU[ctx.IBLTE[2]][2] * ctx.THET[ctx.IBLTE[2]][2]
                    ) / tte
                    tesys(bl, cte, tte, dte)
                else:
                    blsys(bl)

                if itbl == 1:
                    ueref = bl.U2
                    hkref = bl.HK2

                    if ibl < ctx.ITRAN[is_] and ibl >= itrold:
                        uem = ctx.UEDG[ibl - 1][is_]
                        dsm = ctx.DSTR[ibl - 1][is_]
                        thm = ctx.THET[ibl - 1][is_]
                        msq = uem * uem * bl.HSTINV / (
                            bl.GM1BL * (1.0 - 0.5 * uem * uem * bl.HSTINV)
                        )
                        hkref, _, _ = hkin(dsm / thm, msq)

                    if ibl < itrold:
                        if bl.TRAN:
                            ctx.CTAU[ibl][is_] = 0.03
                        if bl.TURB:
                            ctx.CTAU[ibl][is_] = ctx.CTAU[ibl - 1][is_]
                        if bl.TRAN or bl.TURB:
                            cti = ctx.CTAU[ibl][is_]
                            bl.S2 = cti

                if bl.SIMI or ibl == ctx.IBLTE[is_] + 1:
                    bl.VS2[4][1] = 0.0
                    bl.VS2[4][2] = 0.0
                    bl.VS2[4][3] = 0.0
                    bl.VS2[4][4] = bl.U2_UEI
                    bl.VSREZ[4] = ueref - bl.U2
                else:
                    vtmp = make_2d(4, 5)
                    vztmp = make_1d(4)
                    for k in range(1, 5):
                        vztmp[k] = bl.VSREZ[k]
                        for l in range(1, 6):
                            vtmp[k][l] = bl.VS2[k][l]

                    vtmp[4][1] = 0.0
                    vtmp[4][2] = bl.HK2_T2
                    vtmp[4][3] = bl.HK2_D2
                    vtmp[4][4] = bl.HK2_U2 * bl.U2_UEI
                    vztmp[4] = 1.0

                    gauss(4, 4, vtmp, vztmp, 1)

                    sennew = senswt * vztmp[4] * hkref / ueref
                    if itbl <= 5:
                        sens = sennew
                    elif itbl <= 15:
                        sens = 0.5 * (sens + sennew)

                    bl.VS2[4][1] = 0.0
                    bl.VS2[4][2] = bl.HK2_T2 * hkref
                    bl.VS2[4][3] = bl.HK2_D2 * hkref
                    bl.VS2[4][4] = (bl.HK2_U2 * hkref + sens / ueref) * bl.U2_UEI
                    bl.VSREZ[4] = -(hkref**2) * (bl.HK2 / hkref - 1.0) - sens * (bl.U2 / ueref - 1.0)

                gauss(4, 4, bl.VS2, bl.VSREZ, 1)

                dmax = max(abs(bl.VSREZ[2] / thi), abs(bl.VSREZ[3] / dsi), abs(bl.VSREZ[4] / uei))
                if ibl >= ctx.ITRAN[is_]:
                    dmax = max(dmax, abs(bl.VSREZ[1] / (10.0 * cti)))

                rlx = 1.0
                if dmax > 0.3:
                    rlx = 0.3 / dmax

                if ibl < ctx.ITRAN[is_]:
                    ami = ami + rlx * bl.VSREZ[1]
                if ibl >= ctx.ITRAN[is_]:
                    cti = cti + rlx * bl.VSREZ[1]
                thi = thi + rlx * bl.VSREZ[2]
                dsi = dsi + rlx * bl.VSREZ[3]
                uei = uei + rlx * bl.VSREZ[4]

                if ibl >= ctx.ITRAN[is_]:
                    cti = min(cti, 0.30)
                    cti = max(cti, 0.0000001)

                if ibl <= ctx.IBLTE[is_]:
                    hklim = 1.02
                else:
                    hklim = 1.00005
                msq = uei * uei * bl.HSTINV / (bl.GM1BL * (1.0 - 0.5 * uei * uei * bl.HSTINV))
                dsw = dsi - dswaki
                dsw = dslim(dsw, thi, uei, msq, hklim)
                dsi = dsw + dswaki

                if dmax <= deps:
                    break
            else:
                print(f" MRCHDU: Convergence failed at{ibl:4d}  side{is_:2d}    Res ={dmax:12.4E}")
                if dmax <= 0.1:
                    pass
                else:
                    if ibl > 3:
                        if ibl <= ctx.IBLTE[is_]:
                            thi = ctx.THET[ibm][is_] * (ctx.XSSI[ibl][is_] / ctx.XSSI[ibm][is_]) ** 0.5
                            dsi = ctx.DSTR[ibm][is_] * (ctx.XSSI[ibl][is_] / ctx.XSSI[ibm][is_]) ** 0.5
                            uei = ctx.UEDG[ibm][is_]
                        elif ibl == ctx.IBLTE[is_] + 1:
                            cti = cte
                            thi = tte
                            dsi = dte
                            uei = ctx.UEDG[ibm][is_]
                        else:
                            thi = ctx.THET[ibm][is_]
                            ratlen = (ctx.XSSI[ibl][is_] - ctx.XSSI[ibm][is_]) / (10.0 * ctx.DSTR[ibm][is_])
                            dsi = (ctx.DSTR[ibm][is_] + thi * ratlen) / (1.0 + ratlen)
                            uei = ctx.UEDG[ibm][is_]
                        if ibl == ctx.ITRAN[is_]:
                            cti = 0.05
                        if ibl > ctx.ITRAN[is_]:
                            cti = ctx.CTAU[ibm][is_]

                blprv(bl, xsi, ami, cti, thi, dsi, dswaki, uei)
                blkin(bl)

                if (not bl.SIMI) and (not bl.TURB):
                    trchek(bl)
                    ami = bl.AMPL2
                    if bl.TRAN:
                        ctx.ITRAN[is_] = ibl
                    if not bl.TRAN:
                        ctx.ITRAN[is_] = ibl + 2

                if ibl < ctx.ITRAN[is_]:
                    blvar(bl, 1)
                if ibl >= ctx.ITRAN[is_]:
                    blvar(bl, 2)
                if bl.WAKE:
                    blvar(bl, 3)

                if ibl < ctx.ITRAN[is_]:
                    blmid(bl, 1)
                if ibl >= ctx.ITRAN[is_]:
                    blmid(bl, 2)
                if bl.WAKE:
                    blmid(bl, 3)

            sens = sennew

            if ibl < ctx.ITRAN[is_]:
                ctx.CTAU[ibl][is_] = ami
            if ibl >= ctx.ITRAN[is_]:
                ctx.CTAU[ibl][is_] = cti
            ctx.THET[ibl][is_] = thi
            ctx.DSTR[ibl][is_] = dsi
            ctx.UEDG[ibl][is_] = uei
            ctx.MASS[ibl][is_] = dsi * uei
            ctx.TAU[ibl][is_] = 0.5 * bl.R2 * bl.U2 * bl.U2 * bl.CF2
            ctx.DIS[ibl][is_] = bl.R2 * bl.U2 * bl.U2 * bl.U2 * bl.DI2 * bl.HS2 * 0.5
            ctx.CTQ[ibl][is_] = bl.CQ2
            ctx.DELT[ibl][is_] = bl.DE2
            ctx.TSTR[ibl][is_] = bl.HS2 * bl.T2

            blprv(bl, xsi, ami, cti, thi, dsi, dswaki, uei)
            blkin(bl)
            for icom in range(1, NCOM + 1):
                bl.COM1[icom] = bl.COM2[icom]
            sync_com_to_vars(bl, 1)

            if bl.TRAN or ibl == ctx.IBLTE[is_]:
                bl.TURB = True
                ctx.TFORCE[is_] = bl.TRFORC
                ctx.XSSITR[is_] = bl.XT

            bl.TRAN = False


def xifset(ctx, bl, is_):
    if ctx.XSTRIP[is_] >= 1.0:
        bl.XIFORC = ctx.XSSI[ctx.IBLTE[is_]][is_]
        return

    chx = ctx.XTE - ctx.XLE
    chy = ctx.YTE - ctx.YLE
    chsq = chx**2 + chy**2

    for i in range(1, ctx.N + 1):
        ctx.W1[i] = ((ctx.X[i] - ctx.XLE) * chx + (ctx.Y[i] - ctx.YLE) * chy) / chsq
        ctx.W2[i] = ((ctx.Y[i] - ctx.YLE) * chx - (ctx.X[i] - ctx.XLE) * chy) / chsq

    splind(ctx.W1, ctx.W3, ctx.S, ctx.N, -999.0, -999.0)
    splind(ctx.W2, ctx.W4, ctx.S, ctx.N, -999.0, -999.0)

    if is_ == 1:
        str_ = ctx.SLE + (ctx.S[1] - ctx.SLE) * ctx.XSTRIP[is_]
        str_ = sinvrt(str_, ctx.XSTRIP[is_], ctx.W1, ctx.W3, ctx.S, ctx.N)
        bl.XIFORC = min((ctx.SST - str_), ctx.XSSI[ctx.IBLTE[is_]][is_])
    else:
        str_ = ctx.SLE + (ctx.S[ctx.N] - ctx.SLE) * ctx.XSTRIP[is_]
        str_ = sinvrt(str_, ctx.XSTRIP[is_], ctx.W1, ctx.W3, ctx.S, ctx.N)
        bl.XIFORC = min((str_ - ctx.SST), ctx.XSSI[ctx.IBLTE[is_]][is_])

    if bl.XIFORC < 0.0:
        print()
        print(f" ***  Stagnation point is past trip on side{is_:2d}  ***")
        bl.XIFORC = ctx.XSSI[ctx.IBLTE[is_]][is_]


def update(ctx, bl):
    unew = make_2d(IVX, 2)
    u_ac = make_2d(IVX, 2)
    qnew = make_1d(IQX)
    q_ac = make_1d(IQX)

    dalmax = 0.5 * ctx.DTOR
    dalmin = -0.5 * ctx.DTOR

    dclmax = 0.5
    dclmin = -0.5
    if ctx.MATYP != 1:
        dclmin = max(-0.5, -0.9 * ctx.CL)

    hstinv = ctx.GAMM1 * (ctx.MINF / ctx.QINF) ** 2 / (1.0 + 0.5 * ctx.GAMM1 * ctx.MINF**2)

    for is_ in range(1, 3):
        for ibl in range(2, ctx.NBL[is_] + 1):
            i = ctx.IPAN[ibl][is_]

            dui = 0.0
            dui_ac = 0.0
            for js in range(1, 3):
                for jbl in range(2, ctx.NBL[js] + 1):
                    j = ctx.IPAN[jbl][js]
                    jv = ctx.ISYS[jbl][js]
                    ue_m = -ctx.VTI[ibl][is_] * ctx.VTI[jbl][js] * ctx.DIJ[i][j]
                    dui = dui + ue_m * (ctx.MASS[jbl][js] + ctx.VDEL[3][1][jv])
                    dui_ac = dui_ac + ue_m * (0.0 - ctx.VDEL[3][2][jv])

            if ctx.LALFA:
                uinv_ac = 0.0
            else:
                uinv_ac = ctx.UINV_A[ibl][is_]

            unew[ibl][is_] = ctx.UINV[ibl][is_] + dui
            u_ac[ibl][is_] = uinv_ac + dui_ac

    for is_ in range(1, 3):
        for ibl in range(2, ctx.IBLTE[is_] + 1):
            i = ctx.IPAN[ibl][is_]
            qnew[i] = ctx.VTI[ibl][is_] * unew[ibl][is_]
            q_ac[i] = ctx.VTI[ibl][is_] * u_ac[ibl][is_]

    sa = math.sin(ctx.ALFA)
    ca = math.cos(ctx.ALFA)

    beta = math.sqrt(1.0 - ctx.MINF**2)
    beta_msq = -0.5 / beta

    bfac = 0.5 * ctx.MINF**2 / (1.0 + beta)
    bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq

    clnew = 0.0
    cl_a = 0.0
    cl_ms = 0.0
    cl_ac = 0.0

    i = 1
    cginc = 1.0 - (qnew[i] / ctx.QINF) ** 2
    cpg1 = cginc / (beta + bfac * cginc)
    cpg1_ms = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc)

    cpi_q = -2.0 * qnew[i] / ctx.QINF**2
    cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc)
    cpg1_ac = cpc_cpi * cpi_q * q_ac[i]

    for i in range(1, ctx.N + 1):
        ip = i + 1
        if i == ctx.N:
            ip = 1

        cginc = 1.0 - (qnew[ip] / ctx.QINF) ** 2
        cpg2 = cginc / (beta + bfac * cginc)
        cpg2_ms = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc)

        cpi_q = -2.0 * qnew[ip] / ctx.QINF**2
        cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc)
        cpg2_ac = cpc_cpi * cpi_q * q_ac[ip]

        dx = (ctx.X[ip] - ctx.X[i]) * ca + (ctx.Y[ip] - ctx.Y[i]) * sa
        dx_a = -(ctx.X[ip] - ctx.X[i]) * sa + (ctx.Y[ip] - ctx.Y[i]) * ca

        ag = 0.5 * (cpg2 + cpg1)
        ag_ms = 0.5 * (cpg2_ms + cpg1_ms)
        ag_ac = 0.5 * (cpg2_ac + cpg1_ac)

        clnew = clnew + dx * ag
        cl_a = cl_a + dx_a * ag
        cl_ms = cl_ms + dx * ag_ms
        cl_ac = cl_ac + dx * ag_ac

        cpg1 = cpg2
        cpg1_ms = cpg2_ms
        cpg1_ac = cpg2_ac

    rlx = 1.0

    if ctx.LALFA:
        dac = (clnew - ctx.CL) / (1.0 - cl_ac - cl_ms * 2.0 * ctx.MINF * ctx.MINF_CL)

        if rlx * dac > dclmax:
            rlx = dclmax / dac
        if rlx * dac < dclmin:
            rlx = dclmin / dac
    else:
        dac = (clnew - ctx.CLSPEC) / (0.0 - cl_ac - cl_a)

        if rlx * dac > dalmax:
            rlx = dalmax / dac
        if rlx * dac < dalmin:
            rlx = dalmin / dac

    ctx.RMSBL = 0.0
    ctx.RMXBL = 0.0

    dhi = 1.5
    dlo = -0.5

    for is_ in range(1, 3):
        for ibl in range(2, ctx.NBL[is_] + 1):
            iv = ctx.ISYS[ibl][is_]

            dctau = ctx.VDEL[1][1][iv] - dac * ctx.VDEL[1][2][iv]
            dthet = ctx.VDEL[2][1][iv] - dac * ctx.VDEL[2][2][iv]
            dmass = ctx.VDEL[3][1][iv] - dac * ctx.VDEL[3][2][iv]
            duedg = unew[ibl][is_] + dac * u_ac[ibl][is_] - ctx.UEDG[ibl][is_]
            ddstr = (dmass - ctx.DSTR[ibl][is_] * duedg) / ctx.UEDG[ibl][is_]

            if ibl < ctx.ITRAN[is_]:
                dn1 = dctau / 10.0
            if ibl >= ctx.ITRAN[is_]:
                dn1 = dctau / ctx.CTAU[ibl][is_]
            dn2 = dthet / ctx.THET[ibl][is_]
            dn3 = ddstr / ctx.DSTR[ibl][is_]
            dn4 = abs(duedg) / 0.25

            ctx.RMSBL = ctx.RMSBL + dn1**2 + dn2**2 + dn3**2 + dn4**2

            rdn1 = rlx * dn1
            if abs(dn1) > abs(ctx.RMXBL):
                ctx.RMXBL = dn1
                if ibl < ctx.ITRAN[is_]:
                    ctx.VMXBL = "n"
                if ibl >= ctx.ITRAN[is_]:
                    ctx.VMXBL = "C"
                ctx.IMXBL = ibl
                ctx.ISMXBL = is_
            if rdn1 > dhi:
                rlx = dhi / dn1
            if rdn1 < dlo:
                rlx = dlo / dn1

            rdn2 = rlx * dn2
            if abs(dn2) > abs(ctx.RMXBL):
                ctx.RMXBL = dn2
                ctx.VMXBL = "T"
                ctx.IMXBL = ibl
                ctx.ISMXBL = is_
            if rdn2 > dhi:
                rlx = dhi / dn2
            if rdn2 < dlo:
                rlx = dlo / dn2

            rdn3 = rlx * dn3
            if abs(dn3) > abs(ctx.RMXBL):
                ctx.RMXBL = dn3
                ctx.VMXBL = "D"
                ctx.IMXBL = ibl
                ctx.ISMXBL = is_
            if rdn3 > dhi:
                rlx = dhi / dn3
            if rdn3 < dlo:
                rlx = dlo / dn3

            rdn4 = rlx * dn4
            if abs(dn4) > abs(ctx.RMXBL):
                ctx.RMXBL = duedg
                ctx.VMXBL = "U"
                ctx.IMXBL = ibl
                ctx.ISMXBL = is_
            if rdn4 > dhi:
                rlx = dhi / dn4
            if rdn4 < dlo:
                rlx = dlo / dn4

    ctx.RMSBL = math.sqrt(ctx.RMSBL / (4.0 * float(ctx.NBL[1] + ctx.NBL[2])))
    ctx.RLX = rlx

    if ctx.LALFA:
        ctx.CL = ctx.CL + rlx * dac
    else:
        ctx.ALFA = ctx.ALFA + rlx * dac
        ctx.ADEG = ctx.ALFA / ctx.DTOR

    for is_ in range(1, 3):
        for ibl in range(2, ctx.NBL[is_] + 1):
            iv = ctx.ISYS[ibl][is_]

            dctau = ctx.VDEL[1][1][iv] - dac * ctx.VDEL[1][2][iv]
            dthet = ctx.VDEL[2][1][iv] - dac * ctx.VDEL[2][2][iv]
            dmass = ctx.VDEL[3][1][iv] - dac * ctx.VDEL[3][2][iv]
            duedg = unew[ibl][is_] + dac * u_ac[ibl][is_] - ctx.UEDG[ibl][is_]
            ddstr = (dmass - ctx.DSTR[ibl][is_] * duedg) / ctx.UEDG[ibl][is_]

            ctx.CTAU[ibl][is_] = ctx.CTAU[ibl][is_] + rlx * dctau
            ctx.THET[ibl][is_] = ctx.THET[ibl][is_] + rlx * dthet
            ctx.DSTR[ibl][is_] = ctx.DSTR[ibl][is_] + rlx * ddstr
            ctx.UEDG[ibl][is_] = ctx.UEDG[ibl][is_] + rlx * duedg

            if ibl > ctx.IBLTE[is_]:
                iw = ibl - ctx.IBLTE[is_]
                dswaki = ctx.WGAP[iw]
            else:
                dswaki = 0.0

            if ibl >= ctx.ITRAN[is_]:
                ctx.CTAU[ibl][is_] = min(ctx.CTAU[ibl][is_], 0.25)

            if ibl <= ctx.IBLTE[is_]:
                hklim = 1.02
            else:
                hklim = 1.00005
            msq = ctx.UEDG[ibl][is_] ** 2 * hstinv / (
                ctx.GAMM1 * (1.0 - 0.5 * ctx.UEDG[ibl][is_] ** 2 * hstinv)
            )
            dsw = ctx.DSTR[ibl][is_] - dswaki
            dsw = dslim(dsw, ctx.THET[ibl][is_], ctx.UEDG[ibl][is_], msq, hklim)
            ctx.DSTR[ibl][is_] = dsw + dswaki

            ctx.MASS[ibl][is_] = ctx.DSTR[ibl][is_] * ctx.UEDG[ibl][is_]

        for ibl in range(3, ctx.IBLTE[is_] + 1):
            if ctx.UEDG[ibl - 1][is_] > 0.0 and ctx.UEDG[ibl][is_] <= 0.0:
                ctx.UEDG[ibl][is_] = ctx.UEDG[ibl - 1][is_]
                ctx.MASS[ibl][is_] = ctx.DSTR[ibl][is_] * ctx.UEDG[ibl][is_]

    for kbl in range(1, ctx.NBL[2] - ctx.IBLTE[2] + 1):
        ctx.CTAU[ctx.IBLTE[1] + kbl][1] = ctx.CTAU[ctx.IBLTE[2] + kbl][2]
        ctx.THET[ctx.IBLTE[1] + kbl][1] = ctx.THET[ctx.IBLTE[2] + kbl][2]
        ctx.DSTR[ctx.IBLTE[1] + kbl][1] = ctx.DSTR[ctx.IBLTE[2] + kbl][2]
        ctx.UEDG[ctx.IBLTE[1] + kbl][1] = ctx.UEDG[ctx.IBLTE[2] + kbl][2]
        ctx.TAU[ctx.IBLTE[1] + kbl][1] = ctx.TAU[ctx.IBLTE[2] + kbl][2]
        ctx.DIS[ctx.IBLTE[1] + kbl][1] = ctx.DIS[ctx.IBLTE[2] + kbl][2]
        ctx.CTQ[ctx.IBLTE[1] + kbl][1] = ctx.CTQ[ctx.IBLTE[2] + kbl][2]
        ctx.DELT[ctx.IBLTE[1] + kbl][1] = ctx.DELT[ctx.IBLTE[2] + kbl][2]
        ctx.TSTR[ctx.IBLTE[1] + kbl][1] = ctx.TSTR[ctx.IBLTE[2] + kbl][2]


def dslim(dstr, thet, uedg, msq, hklim):
    h = dstr / thet
    hk, hk_h, _ = hkin(h, msq)
    dh = max(0.0, hklim - hk) / hk_h
    dstr = dstr + dh * thet
    return dstr


def blpini(bl):
    bl.SCCON = 5.6
    bl.GACON = 6.70
    bl.GBCON = 0.75
    bl.GCCON = 18.0
    bl.DLCON = 0.9

    bl.CTRCON = 1.8
    bl.CTRCEX = 3.3

    bl.DUXCON = 1.0

    bl.CTCON = 0.5 / (bl.GACON**2 * bl.GBCON)

    bl.CFFAC = 1.0
