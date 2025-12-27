#include "xbl.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <tuple>

#include "spline.h"
#include "xblcom.h"
#include "xblsys.h"
#include "xfoil.h"
#include "xpanel.h"
#include "xsolve.h"

XFoilState::XFoilState() {
    LALFA = false;
    LBLINI = false;

    CL = 0.0;
    CM = 0.0;
    CD = 0.0;
    CDP = 0.0;
    CDF = 0.0;
    CL_ALF = 0.0;
    CL_MSQ = 0.0;
    CLSPEC = 0.0;
    MINF = 0.0;
    MINF1 = 0.0;
    MINF_CL = 0.0;
    QINF = 1.0;
    TKLAM = 0.0;
    TKL_MSQ = 0.0;
    GAMMA = 1.4;
    GAMM1 = GAMMA - 1.0;
    HVRAT = 0.0;
    REINF = 0.0;
    REINF1 = 0.0;
    REINF_CL = 0.0;
    RETYP = 1;
    MATYP = 1;
    CPSTAR = 0.0;
    QSTAR = 0.0;
    IDAMP = 0;

    SLE = 0.0;
    XLE = 0.0;
    YLE = 0.0;
    XTE = 0.0;
    YTE = 0.0;
    SST = 0.0;
    SST_GO = 0.0;
    SST_GP = 0.0;

    ALFA = 0.0;
    ADEG = 0.0;
    COSA = 1.0;
    SINA = 0.0;
    DTOR = std::acos(-1.0) / 180.0;
    AVISC = 0.0;
    MVISC = 0.0;
    XCMREF = 0.0;
    YCMREF = 0.0;

    N = 0;
    IST = 0;
    NW = 0;
    NPAN = 0;
    NB = 0;
    LCLOCK = false;

    X = make_1d(IZX);
    Y = make_1d(IZX);
    XP = make_1d(IZX);
    YP = make_1d(IZX);
    S = make_1d(IZX);
    SNEW = make_1d(5 * IBX);

    W1 = make_1d(6 * IQX);
    W2 = make_1d(6 * IQX);
    W3 = make_1d(6 * IQX);
    W4 = make_1d(6 * IQX);
    W5 = make_1d(6 * IQX);
    W6 = make_1d(6 * IQX);
    W7 = make_1d(6 * IQX);
    W8 = make_1d(6 * IQX);

    XB = make_1d(IBX);
    YB = make_1d(IBX);
    XBP = make_1d(IBX);
    YBP = make_1d(IBX);
    SB = make_1d(IBX);
    SBLE = 0.0;
    CHORDB = 0.0;
    AREAB = 0.0;
    RADBLE = 0.0;
    ANGBTE = 0.0;
    EI11BA = 0.0;
    EI22BA = 0.0;
    APX1BA = 0.0;
    APX2BA = 0.0;
    EI11BT = 0.0;
    EI22BT = 0.0;
    APX1BT = 0.0;
    APX2BT = 0.0;
    THICKB = 0.0;
    CAMBRB = 0.0;

    XSSI = make_2d(IVX, ISX);
    UEDG = make_2d(IVX, ISX);
    UINV = make_2d(IVX, ISX);
    UINV_A = make_2d(IVX, ISX);
    MASS = make_2d(IVX, ISX);
    THET = make_2d(IVX, ISX);
    DSTR = make_2d(IVX, ISX);
    CTAU = make_2d(IVX, ISX);
    DELT = make_2d(IVX, ISX);
    TSTR = make_2d(IVX, ISX);
    USLP = make_2d(IVX, ISX);
    GUXQ = make_2d(IVX, ISX);
    GUXD = make_2d(IVX, ISX);
    TAU = make_2d(IVX, ISX);
    DIS = make_2d(IVX, ISX);
    CTQ = make_2d(IVX, ISX);
    VTI = make_2d(IVX, ISX);

    ACRIT = make_1d(ISX);
    XSTRIP = make_1d(ISX);
    XOCTR = make_1d(ISX);
    YOCTR = make_1d(ISX);
    XSSITR = make_1d(ISX);
    TINDEX = make_1d(ISX);

    IBLTE = std::vector<int>(ISX + 1, 0);
    NBL = std::vector<int>(ISX + 1, 0);
    IPAN = std::vector<std::vector<int>>(IVX + 1, std::vector<int>(ISX + 1, 0));
    ISYS = std::vector<std::vector<int>>(IVX + 1, std::vector<int>(ISX + 1, 0));
    NSYS = 0;
    ITRAN = std::vector<int>(ISX + 1, 0);
    TFORCE = std::vector<bool>(ISX + 1, false);

    WGAP = make_1d(IWX);
    DWTE = 0.0;
    ANTE = 0.0;
    DSTE = 0.0;
    ASTE = 0.0;
    WAKLEN = 0.0;
    CHORD = 0.0;
    YIMAGE = 0.0;
    SHARP = false;

    DIJ = make_2d(IZX, IZX);
    AIJ = make_2d(IQX, IQX);
    BIJ = make_2d(IQX, IZX);
    CIJ = make_2d(IWX, IQX);
    VM = make_3d(3, IZX, IZX);
    VA = make_3d(3, 2, IZX);
    VB = make_3d(3, 2, IZX);
    VDEL = make_3d(3, 2, IZX);
    VZ = make_2d(3, 2);
    AIJPIV = std::vector<int>(IQX + 1, 0);

    QINV = make_1d(IZX);
    QVIS = make_1d(IZX);
    CPI = make_1d(IZX);
    CPV = make_1d(IZX);
    QINVU = make_2d(IZX, 2);
    QINV_A = make_1d(IZX);

    GAM = make_1d(IQX);
    GAMU = make_2d(IQX, 2);
    GAM_A = make_1d(IQX);
    SIG = make_1d(IZX);
    NX = make_1d(IZX);
    NY = make_1d(IZX);
    APANEL = make_1d(IZX);
    GAMTE = 0.0;
    GAMTE_A = 0.0;
    SIGTE = 0.0;
    SIGTE_A = 0.0;

    DZDG = make_1d(IQX);
    DZDN = make_1d(IQX);
    DZDM = make_1d(IZX);
    DQDG = make_1d(IQX);
    DQDM = make_1d(IZX);
    QTAN1 = 0.0;
    QTAN2 = 0.0;
    Z_QINF = 0.0;
    Z_ALFA = 0.0;
    Z_QDOF0 = 0.0;
    Z_QDOF1 = 0.0;
    Z_QDOF2 = 0.0;
    Z_QDOF3 = 0.0;

    QF0 = make_1d(IQX);
    QF1 = make_1d(IQX);
    QF2 = make_1d(IQX);
    QF3 = make_1d(IQX);

    PI = std::acos(-1.0);
    HOPI = 0.5 / PI;
    QOPI = 0.25 / PI;

    RMSBL = 0.0;
    RMXBL = 0.0;
    RLX = 0.0;
    VACCEL = 0.0;
    IMXBL = 0;
    ISMXBL = 0;
    VMXBL = " ";

    LQAIJ = false;
    LADIJ = false;
    LWDIJ = false;
    LWAKE = false;
    LGAMU = false;
    LVISC = false;
    LVCONV = false;
    LFLAP = false;
    LIMAGE = false;
    LQINU = false;
    LQSPEC = false;
    LGSAME = false;
    LSCINI = false;
    LIPAN = false;
    AWAKE = 0.0;

    XOF = 0.0;
    YOF = 0.0;
    HMOM = 0.0;
    HFX = 0.0;
    HFY = 0.0;

    CVPAR = 0.0;
    CTERAT = 0.0;
    CTRRAT = 0.0;
    XSREF1 = 0.0;
    XSREF2 = 0.0;
    XPREF1 = 0.0;
    XPREF2 = 0.0;

    XBF = 0.0;
    YBF = 0.0;
    LBFLAP = false;

    NAME = "";
    NNAME = 0;

    HTARG = std::vector<std::vector<double>>();
}

XBlState::XBlState() {
    COM1 = make_1d(NCOM);
    COM2 = make_1d(NCOM);
    C1SAV = make_1d(NCOM);
    C2SAV = make_1d(NCOM);

    SIMI = false;
    TRAN = false;
    TURB = false;
    WAKE = false;
    TRFORC = false;
    TRFREE = false;
    IDAMPV = 0;

    VS1 = make_2d(4, 5);
    VS2 = make_2d(4, 5);
    VSREZ = make_1d(4);
    VSR = make_1d(4);
    VSM = make_1d(4);
    VSX = make_1d(4);

    SCCON = 0.0;
    GACON = 0.0;
    GBCON = 0.0;
    GCCON = 0.0;
    DLCON = 0.0;
    CTRCON = 0.0;
    CTRCEX = 0.0;
    DUXCON = 0.0;
    CTCON = 0.0;
    CFFAC = 1.0;

    X1 = U1 = T1 = D1 = S1 = AMPL1 = U1_UEI = U1_MS = DW1 = H1 = H1_T1 = H1_D1 = 0.0;
    M1 = M1_U1 = M1_MS = R1 = R1_U1 = R1_MS = V1 = V1_U1 = V1_MS = V1_RE = 0.0;
    HK1 = HK1_U1 = HK1_T1 = HK1_D1 = HK1_MS = 0.0;
    HS1 = HS1_U1 = HS1_T1 = HS1_D1 = HS1_MS = HS1_RE = 0.0;
    HC1 = HC1_U1 = HC1_T1 = HC1_D1 = HC1_MS = 0.0;
    RT1 = RT1_U1 = RT1_T1 = RT1_MS = RT1_RE = 0.0;
    CF1 = CF1_U1 = CF1_T1 = CF1_D1 = CF1_MS = CF1_RE = 0.0;
    DI1 = DI1_U1 = DI1_T1 = DI1_D1 = DI1_S1 = DI1_MS = DI1_RE = 0.0;
    US1 = US1_U1 = US1_T1 = US1_D1 = US1_MS = US1_RE = 0.0;
    CQ1 = CQ1_U1 = CQ1_T1 = CQ1_D1 = CQ1_MS = CQ1_RE = 0.0;
    DE1 = DE1_U1 = DE1_T1 = DE1_D1 = DE1_MS = 0.0;

    X2 = U2 = T2 = D2 = S2 = AMPL2 = U2_UEI = U2_MS = DW2 = H2 = H2_T2 = H2_D2 = 0.0;
    M2 = M2_U2 = M2_MS = R2 = R2_U2 = R2_MS = V2 = V2_U2 = V2_MS = V2_RE = 0.0;
    HK2 = HK2_U2 = HK2_T2 = HK2_D2 = HK2_MS = 0.0;
    HS2 = HS2_U2 = HS2_T2 = HS2_D2 = HS2_MS = HS2_RE = 0.0;
    HC2 = HC2_U2 = HC2_T2 = HC2_D2 = HC2_MS = 0.0;
    RT2 = RT2_U2 = RT2_T2 = RT2_MS = RT2_RE = 0.0;
    CF2 = CF2_HK2 = CF2_M2 = CF2_RT2 = CF2_U2 = CF2_T2 = CF2_D2 = CF2_MS = CF2_RE = 0.0;
    DI2 = DI2_U2 = DI2_T2 = DI2_D2 = DI2_S2 = DI2_MS = DI2_RE = 0.0;
    US2 = US2_U2 = US2_T2 = US2_D2 = US2_MS = US2_RE = 0.0;
    CQ2 = CQ2_U2 = CQ2_T2 = CQ2_D2 = CQ2_MS = CQ2_RE = 0.0;
    DE2 = DE2_U2 = DE2_T2 = DE2_D2 = DE2_MS = 0.0;

    CFM = CFM_HKA = CFM_MA = CFM_MS = CFM_RE = CFM_RTA = CFM_U1 = CFM_T1 = CFM_D1 = CFM_U2 = CFM_T2 = CFM_D2 = 0.0;
    XT = XT_A1 = XT_A2 = XT_MS = XT_RE = XT_XF = XT_X1 = XT_T1 = XT_D1 = XT_U1 = XT_X2 = XT_T2 = XT_D2 = XT_U2 = 0.0;
    DWTE = QINFBL = TKBL = TKBL_MS = RSTBL = RSTBL_MS = HSTINV = HSTINV_MS = 0.0;
    REYBL = REYBL_MS = REYBL_RE = GAMBL = GM1BL = HVRAT = BULE = XIFORC = AMCRIT = 0.0;

    com1_fields = std::vector<double *>(NCOM + 1, nullptr);
    com2_fields = std::vector<double *>(NCOM + 1, nullptr);

    com1_fields[1] = &X1;
    com1_fields[2] = &U1;
    com1_fields[3] = &T1;
    com1_fields[4] = &D1;
    com1_fields[5] = &S1;
    com1_fields[6] = &AMPL1;
    com1_fields[7] = &U1_UEI;
    com1_fields[8] = &U1_MS;
    com1_fields[9] = &DW1;
    com1_fields[10] = &H1;
    com1_fields[11] = &H1_T1;
    com1_fields[12] = &H1_D1;
    com1_fields[13] = &M1;
    com1_fields[14] = &M1_U1;
    com1_fields[15] = &M1_MS;
    com1_fields[16] = &R1;
    com1_fields[17] = &R1_U1;
    com1_fields[18] = &R1_MS;
    com1_fields[19] = &V1;
    com1_fields[20] = &V1_U1;
    com1_fields[21] = &V1_MS;
    com1_fields[22] = &V1_RE;
    com1_fields[23] = &HK1;
    com1_fields[24] = &HK1_U1;
    com1_fields[25] = &HK1_T1;
    com1_fields[26] = &HK1_D1;
    com1_fields[27] = &HK1_MS;
    com1_fields[28] = &HS1;
    com1_fields[29] = &HS1_U1;
    com1_fields[30] = &HS1_T1;
    com1_fields[31] = &HS1_D1;
    com1_fields[32] = &HS1_MS;
    com1_fields[33] = &HS1_RE;
    com1_fields[34] = &HC1;
    com1_fields[35] = &HC1_U1;
    com1_fields[36] = &HC1_T1;
    com1_fields[37] = &HC1_D1;
    com1_fields[38] = &HC1_MS;
    com1_fields[39] = &RT1;
    com1_fields[40] = &RT1_U1;
    com1_fields[41] = &RT1_T1;
    com1_fields[42] = &RT1_MS;
    com1_fields[43] = &RT1_RE;
    com1_fields[44] = &CF1;
    com1_fields[45] = &CF1_U1;
    com1_fields[46] = &CF1_T1;
    com1_fields[47] = &CF1_D1;
    com1_fields[48] = &CF1_MS;
    com1_fields[49] = &CF1_RE;
    com1_fields[50] = &DI1;
    com1_fields[51] = &DI1_U1;
    com1_fields[52] = &DI1_T1;
    com1_fields[53] = &DI1_D1;
    com1_fields[54] = &DI1_S1;
    com1_fields[55] = &DI1_MS;
    com1_fields[56] = &DI1_RE;
    com1_fields[57] = &US1;
    com1_fields[58] = &US1_U1;
    com1_fields[59] = &US1_T1;
    com1_fields[60] = &US1_D1;
    com1_fields[61] = &US1_MS;
    com1_fields[62] = &US1_RE;
    com1_fields[63] = &CQ1;
    com1_fields[64] = &CQ1_U1;
    com1_fields[65] = &CQ1_T1;
    com1_fields[66] = &CQ1_D1;
    com1_fields[67] = &CQ1_MS;
    com1_fields[68] = &CQ1_RE;
    com1_fields[69] = &DE1;
    com1_fields[70] = &DE1_U1;
    com1_fields[71] = &DE1_T1;
    com1_fields[72] = &DE1_D1;
    com1_fields[73] = &DE1_MS;

    com2_fields[1] = &X2;
    com2_fields[2] = &U2;
    com2_fields[3] = &T2;
    com2_fields[4] = &D2;
    com2_fields[5] = &S2;
    com2_fields[6] = &AMPL2;
    com2_fields[7] = &U2_UEI;
    com2_fields[8] = &U2_MS;
    com2_fields[9] = &DW2;
    com2_fields[10] = &H2;
    com2_fields[11] = &H2_T2;
    com2_fields[12] = &H2_D2;
    com2_fields[13] = &M2;
    com2_fields[14] = &M2_U2;
    com2_fields[15] = &M2_MS;
    com2_fields[16] = &R2;
    com2_fields[17] = &R2_U2;
    com2_fields[18] = &R2_MS;
    com2_fields[19] = &V2;
    com2_fields[20] = &V2_U2;
    com2_fields[21] = &V2_MS;
    com2_fields[22] = &V2_RE;
    com2_fields[23] = &HK2;
    com2_fields[24] = &HK2_U2;
    com2_fields[25] = &HK2_T2;
    com2_fields[26] = &HK2_D2;
    com2_fields[27] = &HK2_MS;
    com2_fields[28] = &HS2;
    com2_fields[29] = &HS2_U2;
    com2_fields[30] = &HS2_T2;
    com2_fields[31] = &HS2_D2;
    com2_fields[32] = &HS2_MS;
    com2_fields[33] = &HS2_RE;
    com2_fields[34] = &HC2;
    com2_fields[35] = &HC2_U2;
    com2_fields[36] = &HC2_T2;
    com2_fields[37] = &HC2_D2;
    com2_fields[38] = &HC2_MS;
    com2_fields[39] = &RT2;
    com2_fields[40] = &RT2_U2;
    com2_fields[41] = &RT2_T2;
    com2_fields[42] = &RT2_MS;
    com2_fields[43] = &RT2_RE;
    com2_fields[44] = &CF2;
    com2_fields[45] = &CF2_U2;
    com2_fields[46] = &CF2_T2;
    com2_fields[47] = &CF2_D2;
    com2_fields[48] = &CF2_MS;
    com2_fields[49] = &CF2_RE;
    com2_fields[50] = &DI2;
    com2_fields[51] = &DI2_U2;
    com2_fields[52] = &DI2_T2;
    com2_fields[53] = &DI2_D2;
    com2_fields[54] = &DI2_S2;
    com2_fields[55] = &DI2_MS;
    com2_fields[56] = &DI2_RE;
    com2_fields[57] = &US2;
    com2_fields[58] = &US2_U2;
    com2_fields[59] = &US2_T2;
    com2_fields[60] = &US2_D2;
    com2_fields[61] = &US2_MS;
    com2_fields[62] = &US2_RE;
    com2_fields[63] = &CQ2;
    com2_fields[64] = &CQ2_U2;
    com2_fields[65] = &CQ2_T2;
    com2_fields[66] = &CQ2_D2;
    com2_fields[67] = &CQ2_MS;
    com2_fields[68] = &CQ2_RE;
    com2_fields[69] = &DE2;
    com2_fields[70] = &DE2_U2;
    com2_fields[71] = &DE2_T2;
    com2_fields[72] = &DE2_D2;
    com2_fields[73] = &DE2_MS;
}

void setbl(XFoilState &ctx, XBlState &bl) {
    const double clmr = ctx.LALFA ? ctx.CL : ctx.CLSPEC;

    const auto mr = mrcl(ctx, clmr);
    const double ma_clmr = mr.first;
    const double re_clmr = mr.second;
    const double msq_clmr = 2.0 * ctx.MINF * ma_clmr;

    comset(ctx);

    bl.GAMBL = ctx.GAMMA;
    bl.GM1BL = ctx.GAMM1;

    bl.QINFBL = ctx.QINF;
    bl.TKBL = ctx.TKLAM;
    bl.TKBL_MS = ctx.TKL_MSQ;

    bl.RSTBL = std::pow(1.0 + 0.5 * bl.GM1BL * ctx.MINF * ctx.MINF, 1.0 / bl.GM1BL);
    bl.RSTBL_MS = 0.5 * bl.RSTBL / (1.0 + 0.5 * bl.GM1BL * ctx.MINF * ctx.MINF);

    bl.HSTINV = bl.GM1BL * (ctx.MINF / bl.QINFBL) * (ctx.MINF / bl.QINFBL)
                / (1.0 + 0.5 * bl.GM1BL * ctx.MINF * ctx.MINF);
    bl.HSTINV_MS = bl.GM1BL * (1.0 / bl.QINFBL) * (1.0 / bl.QINFBL)
                   / (1.0 + 0.5 * bl.GM1BL * ctx.MINF * ctx.MINF)
                   - 0.5 * bl.GM1BL * bl.HSTINV / (1.0 + 0.5 * bl.GM1BL * ctx.MINF * ctx.MINF);

    const double herat = 1.0 - 0.5 * bl.QINFBL * bl.QINFBL * bl.HSTINV;
    const double herat_ms = -0.5 * bl.QINFBL * bl.QINFBL * bl.HSTINV_MS;

    bl.HVRAT = ctx.HVRAT;
    bl.REYBL = ctx.REINF * std::sqrt(herat * herat * herat) * (1.0 + bl.HVRAT) / (herat + bl.HVRAT);
    bl.REYBL_RE = std::sqrt(herat * herat * herat) * (1.0 + bl.HVRAT) / (herat + bl.HVRAT);
    bl.REYBL_MS = bl.REYBL * (1.5 / herat - 1.0 / (herat + bl.HVRAT)) * herat_ms;

    bl.IDAMPV = ctx.IDAMP;
    bl.DWTE = ctx.WGAP[1];

    if (!ctx.LBLINI) {
        std::cout << std::endl;
        std::cout << "Initializing BL ..." << std::endl;
        mrchue(ctx, bl);
        ctx.LBLINI = true;
    }

    std::cout << std::endl;

    mrchdu(ctx, bl);

    std::vector<std::vector<double>> usav = make_2d(IVX, 2);
    for (int is_ = 1; is_ <= 2; ++is_) {
        for (int ibl = 2; ibl <= ctx.NBL[is_]; ++ibl) {
            usav[ibl][is_] = ctx.UEDG[ibl][is_];
        }
    }

    ueset(ctx);

    for (int is_ = 1; is_ <= 2; ++is_) {
        for (int ibl = 2; ibl <= ctx.NBL[is_]; ++ibl) {
            const double temp = usav[ibl][is_];
            usav[ibl][is_] = ctx.UEDG[ibl][is_];
            ctx.UEDG[ibl][is_] = temp;
        }
    }

    const int ile1 = ctx.IPAN[2][1];
    const int ile2 = ctx.IPAN[2][2];
    const int ite1 = ctx.IPAN[ctx.IBLTE[1]][1];
    const int ite2 = ctx.IPAN[ctx.IBLTE[2]][2];

    const int jvte1 = ctx.ISYS[ctx.IBLTE[1]][1];
    const int jvte2 = ctx.ISYS[ctx.IBLTE[2]][2];

    const double dule1 = ctx.UEDG[2][1] - usav[2][1];
    const double dule2 = ctx.UEDG[2][2] - usav[2][2];

    std::vector<double> ule1_m = make_1d(2 * IVX);
    std::vector<double> ule2_m = make_1d(2 * IVX);
    std::vector<double> ute1_m = make_1d(2 * IVX);
    std::vector<double> ute2_m = make_1d(2 * IVX);

    for (int js = 1; js <= 2; ++js) {
        for (int jbl = 2; jbl <= ctx.NBL[js]; ++jbl) {
            const int j = ctx.IPAN[jbl][js];
            const int jv = ctx.ISYS[jbl][js];
            ule1_m[jv] = -ctx.VTI[2][1] * ctx.VTI[jbl][js] * ctx.DIJ[ile1][j];
            ule2_m[jv] = -ctx.VTI[2][2] * ctx.VTI[jbl][js] * ctx.DIJ[ile2][j];
            ute1_m[jv] = -ctx.VTI[ctx.IBLTE[1]][1] * ctx.VTI[jbl][js] * ctx.DIJ[ite1][j];
            ute2_m[jv] = -ctx.VTI[ctx.IBLTE[2]][2] * ctx.VTI[jbl][js] * ctx.DIJ[ite2][j];
        }
    }

    const double ule1_a = ctx.UINV_A[2][1];
    const double ule2_a = ctx.UINV_A[2][2];

    ctx.TINDEX[1] = 0.0;
    ctx.TINDEX[2] = 0.0;

    std::vector<double> u1_m = make_1d(2 * IVX);
    std::vector<double> u2_m = make_1d(2 * IVX);
    std::vector<double> d1_m = make_1d(2 * IVX);
    std::vector<double> d2_m = make_1d(2 * IVX);

    for (int is_ = 1; is_ <= 2; ++is_) {
        for (int js = 1; js <= 2; ++js) {
            for (int jbl = 2; jbl <= ctx.NBL[js]; ++jbl) {
                const int jv = ctx.ISYS[jbl][js];
                u1_m[jv] = 0.0;
                d1_m[jv] = 0.0;
            }
        }
        double u1_a = 0.0;
        double d1_a = 0.0;

        double due1 = 0.0;
        double dds1 = 0.0;

        int ibl = 2;
        bl.BULE = 1.0;
        bl.AMCRIT = ctx.ACRIT[is_];

        xifset(ctx, bl, is_);

        bl.TRAN = false;
        bl.TURB = false;

        double cti = 0.0;
        double ami = 0.0;
        for (ibl = 2; ibl <= ctx.NBL[is_]; ++ibl) {
            const int iv = ctx.ISYS[ibl][is_];

            bl.SIMI = (ibl == 2);
            bl.WAKE = (ibl > ctx.IBLTE[is_]);
            bl.TRAN = (ibl == ctx.ITRAN[is_]);
            bl.TURB = (ibl > ctx.ITRAN[is_]);

            const int i = ctx.IPAN[ibl][is_];

            const double xsi = ctx.XSSI[ibl][is_];
            if (ibl < ctx.ITRAN[is_]) {
                ami = ctx.CTAU[ibl][is_];
            }
            if (ibl >= ctx.ITRAN[is_]) {
                cti = ctx.CTAU[ibl][is_];
            }
            const double uei = ctx.UEDG[ibl][is_];
            const double thi = ctx.THET[ibl][is_];
            const double mdi = ctx.MASS[ibl][is_];

            double dsi = mdi / uei;

            double dswaki = 0.0;
            if (bl.WAKE) {
                const int iw = ibl - ctx.IBLTE[is_];
                dswaki = ctx.WGAP[iw];
            }

            const double d2_m2 = 1.0 / uei;
            const double d2_u2 = -dsi / uei;

            for (int js = 1; js <= 2; ++js) {
                for (int jbl = 2; jbl <= ctx.NBL[js]; ++jbl) {
                    const int j = ctx.IPAN[jbl][js];
                    const int jv = ctx.ISYS[jbl][js];
                    u2_m[jv] = -ctx.VTI[ibl][is_] * ctx.VTI[jbl][js] * ctx.DIJ[i][j];
                    d2_m[jv] = d2_u2 * u2_m[jv];
                }
            }
            d2_m[iv] = d2_m[iv] + d2_m2;

            const double u2_a = ctx.UINV_A[ibl][is_];
            const double d2_a = d2_u2 * u2_a;

            const double due2 = ctx.UEDG[ibl][is_] - usav[ibl][is_];
            const double dds2 = d2_u2 * due2;

            blprv(bl, xsi, ami, cti, thi, dsi, dswaki, uei);
            blkin(bl);

            if (bl.TRAN) {
                trchek(bl);
                ami = bl.AMPL2;
            }
            if (ibl == ctx.ITRAN[is_] && !bl.TRAN) {
                std::cout << "SETBL: Xtr???  n1 n2: " << bl.AMPL1 << " " << bl.AMPL2 << std::endl;
            }

            double tte_tte1 = 0.0;
            double tte_tte2 = 0.0;
            double cte_cte1 = 0.0;
            double cte_cte2 = 0.0;
            double cte_tte1 = 0.0;
            double cte_tte2 = 0.0;

            if (ibl == ctx.IBLTE[is_] + 1) {
                const double tte = ctx.THET[ctx.IBLTE[1]][1] + ctx.THET[ctx.IBLTE[2]][2];
                const double dte = ctx.DSTR[ctx.IBLTE[1]][1] + ctx.DSTR[ctx.IBLTE[2]][2] + ctx.ANTE;
                const double cte = (ctx.CTAU[ctx.IBLTE[1]][1] * ctx.THET[ctx.IBLTE[1]][1]
                                    + ctx.CTAU[ctx.IBLTE[2]][2] * ctx.THET[ctx.IBLTE[2]][2])
                                   / tte;
                tesys(bl, cte, tte, dte);

                tte_tte1 = 1.0;
                tte_tte2 = 1.0;
                const double dte_mte1 = 1.0 / ctx.UEDG[ctx.IBLTE[1]][1];
                const double dte_ute1 = -ctx.DSTR[ctx.IBLTE[1]][1] / ctx.UEDG[ctx.IBLTE[1]][1];
                const double dte_mte2 = 1.0 / ctx.UEDG[ctx.IBLTE[2]][2];
                const double dte_ute2 = -ctx.DSTR[ctx.IBLTE[2]][2] / ctx.UEDG[ctx.IBLTE[2]][2];
                cte_cte1 = ctx.THET[ctx.IBLTE[1]][1] / tte;
                cte_cte2 = ctx.THET[ctx.IBLTE[2]][2] / tte;
                cte_tte1 = (ctx.CTAU[ctx.IBLTE[1]][1] - cte) / tte;
                cte_tte2 = (ctx.CTAU[ctx.IBLTE[2]][2] - cte) / tte;

                for (int js = 1; js <= 2; ++js) {
                    for (int jbl = 2; jbl <= ctx.NBL[js]; ++jbl) {
                        const int j = ctx.IPAN[jbl][js];
                        const int jv = ctx.ISYS[jbl][js];
                        d1_m[jv] = dte_ute1 * ute1_m[jv] + dte_ute2 * ute2_m[jv];
                    }
                }
                d1_m[jvte1] = d1_m[jvte1] + dte_mte1;
                d1_m[jvte2] = d1_m[jvte2] + dte_mte2;

                due1 = 0.0;
                dds1 = dte_ute1 * (ctx.UEDG[ctx.IBLTE[1]][1] - usav[ctx.IBLTE[1]][1])
                       + dte_ute2 * (ctx.UEDG[ctx.IBLTE[2]][2] - usav[ctx.IBLTE[2]][2]);

            } else {
                blsys(bl);
            }

            ctx.TAU[ibl][is_] = 0.5 * bl.R2 * bl.U2 * bl.U2 * bl.CF2;
            ctx.DIS[ibl][is_] = bl.R2 * bl.U2 * bl.U2 * bl.U2 * bl.DI2 * bl.HS2 * 0.5;
            ctx.CTQ[ibl][is_] = bl.CQ2;
            ctx.DELT[ibl][is_] = bl.DE2;
            ctx.USLP[ibl][is_] = 1.60 / (1.0 + bl.US2);

            double xi_ule1 = 0.0;
            double xi_ule2 = 0.0;
            if (is_ == 1) {
                xi_ule1 = ctx.SST_GO;
                xi_ule2 = -ctx.SST_GP;
            } else {
                xi_ule1 = -ctx.SST_GO;
                xi_ule2 = ctx.SST_GP;
            }

            for (int jv = 1; jv <= ctx.NSYS; ++jv) {
                ctx.VM[1][jv][iv] = bl.VS1[1][3] * d1_m[jv] + bl.VS1[1][4] * u1_m[jv] + bl.VS2[1][3] * d2_m[jv]
                                   + bl.VS2[1][4] * u2_m[jv]
                                   + (bl.VS1[1][5] + bl.VS2[1][5] + bl.VSX[1])
                                         * (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
            }
            ctx.VB[1][1][iv] = bl.VS1[1][1];
            ctx.VB[1][2][iv] = bl.VS1[1][2];

            ctx.VA[1][1][iv] = bl.VS2[1][1];
            ctx.VA[1][2][iv] = bl.VS2[1][2];

            if (ctx.LALFA) {
                ctx.VDEL[1][2][iv] = bl.VSR[1] * re_clmr + bl.VSM[1] * msq_clmr;
            } else {
                ctx.VDEL[1][2][iv] = bl.VS1[1][4] * u1_a + bl.VS1[1][3] * d1_a + bl.VS2[1][4] * u2_a
                                     + bl.VS2[1][3] * d2_a
                                     + (bl.VS1[1][5] + bl.VS2[1][5] + bl.VSX[1]) * (xi_ule1 * ule1_a + xi_ule2 * ule2_a);
            }

            ctx.VDEL[1][1][iv] = bl.VSREZ[1] + (bl.VS1[1][4] * due1 + bl.VS1[1][3] * dds1)
                                 + (bl.VS2[1][4] * due2 + bl.VS2[1][3] * dds2)
                                 + (bl.VS1[1][5] + bl.VS2[1][5] + bl.VSX[1]) * (xi_ule1 * dule1 + xi_ule2 * dule2);

            for (int jv = 1; jv <= ctx.NSYS; ++jv) {
                ctx.VM[2][jv][iv] = bl.VS1[2][3] * d1_m[jv] + bl.VS1[2][4] * u1_m[jv] + bl.VS2[2][3] * d2_m[jv]
                                   + bl.VS2[2][4] * u2_m[jv]
                                   + (bl.VS1[2][5] + bl.VS2[2][5] + bl.VSX[2])
                                         * (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
            }
            ctx.VB[2][1][iv] = bl.VS1[2][1];
            ctx.VB[2][2][iv] = bl.VS1[2][2];
            ctx.VA[2][1][iv] = bl.VS2[2][1];
            ctx.VA[2][2][iv] = bl.VS2[2][2];

            if (ctx.LALFA) {
                ctx.VDEL[2][2][iv] = bl.VSR[2] * re_clmr + bl.VSM[2] * msq_clmr;
            } else {
                ctx.VDEL[2][2][iv] = bl.VS1[2][4] * u1_a + bl.VS1[2][3] * d1_a + bl.VS2[2][4] * u2_a
                                     + bl.VS2[2][3] * d2_a
                                     + (bl.VS1[2][5] + bl.VS2[2][5] + bl.VSX[2]) * (xi_ule1 * ule1_a + xi_ule2 * ule2_a);
            }

            ctx.VDEL[2][1][iv] = bl.VSREZ[2] + (bl.VS1[2][4] * due1 + bl.VS1[2][3] * dds1)
                                 + (bl.VS2[2][4] * due2 + bl.VS2[2][3] * dds2)
                                 + (bl.VS1[2][5] + bl.VS2[2][5] + bl.VSX[2]) * (xi_ule1 * dule1 + xi_ule2 * dule2);

            for (int jv = 1; jv <= ctx.NSYS; ++jv) {
                ctx.VM[3][jv][iv] = bl.VS1[3][3] * d1_m[jv] + bl.VS1[3][4] * u1_m[jv] + bl.VS2[3][3] * d2_m[jv]
                                   + bl.VS2[3][4] * u2_m[jv]
                                   + (bl.VS1[3][5] + bl.VS2[3][5] + bl.VSX[3])
                                         * (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
            }
            ctx.VB[3][1][iv] = bl.VS1[3][1];
            ctx.VB[3][2][iv] = bl.VS1[3][2];
            ctx.VA[3][1][iv] = bl.VS2[3][1];
            ctx.VA[3][2][iv] = bl.VS2[3][2];

            if (ctx.LALFA) {
                ctx.VDEL[3][2][iv] = bl.VSR[3] * re_clmr + bl.VSM[3] * msq_clmr;
            } else {
                ctx.VDEL[3][2][iv] = bl.VS1[3][4] * u1_a + bl.VS1[3][3] * d1_a + bl.VS2[3][4] * u2_a
                                     + bl.VS2[3][3] * d2_a
                                     + (bl.VS1[3][5] + bl.VS2[3][5] + bl.VSX[3]) * (xi_ule1 * ule1_a + xi_ule2 * ule2_a);
            }

            ctx.VDEL[3][1][iv] = bl.VSREZ[3] + (bl.VS1[3][4] * due1 + bl.VS1[3][3] * dds1)
                                 + (bl.VS2[3][4] * due2 + bl.VS2[3][3] * dds2)
                                 + (bl.VS1[3][5] + bl.VS2[3][5] + bl.VSX[3]) * (xi_ule1 * dule1 + xi_ule2 * dule2);

            if (ibl == ctx.IBLTE[is_] + 1) {
                ctx.VZ[1][1] = bl.VS1[1][1] * cte_cte1;
                ctx.VZ[1][2] = bl.VS1[1][1] * cte_tte1 + bl.VS1[1][2] * tte_tte1;
                ctx.VB[1][1][iv] = bl.VS1[1][1] * cte_cte2;
                ctx.VB[1][2][iv] = bl.VS1[1][1] * cte_tte2 + bl.VS1[1][2] * tte_tte2;

                ctx.VZ[2][1] = bl.VS1[2][1] * cte_cte1;
                ctx.VZ[2][2] = bl.VS1[2][1] * cte_tte1 + bl.VS1[2][2] * tte_tte1;
                ctx.VB[2][1][iv] = bl.VS1[2][1] * cte_cte2;
                ctx.VB[2][2][iv] = bl.VS1[2][1] * cte_tte2 + bl.VS1[2][2] * tte_tte2;

                ctx.VZ[3][1] = bl.VS1[3][1] * cte_cte1;
                ctx.VZ[3][2] = bl.VS1[3][1] * cte_tte1 + bl.VS1[3][2] * tte_tte1;
                ctx.VB[3][1][iv] = bl.VS1[3][1] * cte_cte2;
                ctx.VB[3][2][iv] = bl.VS1[3][1] * cte_tte2 + bl.VS1[3][2] * tte_tte2;
            }

            if (bl.TRAN) {
                bl.TURB = true;
                ctx.ITRAN[is_] = ibl;
                ctx.TFORCE[is_] = bl.TRFORC;
                ctx.XSSITR[is_] = bl.XT;

                double str_ = 0.0;
                if (is_ == 1) {
                    str_ = ctx.SST - bl.XT;
                } else {
                    str_ = ctx.SST + bl.XT;
                }
                const double chx = ctx.XTE - ctx.XLE;
                const double chy = ctx.YTE - ctx.YLE;
                const double chsq = chx * chx + chy * chy;
                const double xtr = seval(str_, ctx.X, ctx.XP, ctx.S, ctx.N);
                const double ytr = seval(str_, ctx.Y, ctx.YP, ctx.S, ctx.N);
                ctx.XOCTR[is_] = ((xtr - ctx.XLE) * chx + (ytr - ctx.YLE) * chy) / chsq;
                ctx.YOCTR[is_] = ((ytr - ctx.YLE) * chx - (xtr - ctx.XLE) * chy) / chsq;
            }

            bl.TRAN = false;

            if (ibl == ctx.IBLTE[is_]) {
                bl.TURB = true;
                bl.WAKE = true;
                blvar(bl, 3);
                blmid(bl, 3);
            }

            for (int js = 1; js <= 2; ++js) {
                for (int jbl = 2; jbl <= ctx.NBL[js]; ++jbl) {
                    const int jv = ctx.ISYS[jbl][js];
                    u1_m[jv] = u2_m[jv];
                    d1_m[jv] = d2_m[jv];
                }
            }

            u1_a = u2_a;
            d1_a = d2_a;
            due1 = due2;
            dds1 = dds2;

            if (ibl == ctx.ITRAN[is_] && bl.X2 > bl.X1) {
                if (is_ == 1) {
                    ctx.TINDEX[is_] = static_cast<double>(ctx.IST - ctx.ITRAN[is_] + 3) - (bl.XT - bl.X1) / (bl.X2 - bl.X1);
                } else {
                    ctx.TINDEX[is_] = static_cast<double>(ctx.IST + ctx.ITRAN[is_] - 2) + (bl.XT - bl.X1) / (bl.X2 - bl.X1);
                }
            }

            for (int icom = 1; icom <= NCOM; ++icom) {
                bl.COM1[icom] = bl.COM2[icom];
            }
            sync_com_to_vars(bl, 1);
        }

        if (ctx.TFORCE[is_]) {
            std::cout << " Side" << std::setw(2) << is_ << " forced transition at x/c = " << std::setw(7) << std::fixed
                      << std::setprecision(4) << ctx.XOCTR[is_] << std::setw(5) << ctx.ITRAN[is_] << std::defaultfloat << std::endl;
        } else {
            std::cout << " Side" << std::setw(2) << is_ << "  free  transition at x/c = " << std::setw(7) << std::fixed
                      << std::setprecision(4) << ctx.XOCTR[is_] << std::setw(5) << ctx.ITRAN[is_] << std::defaultfloat << std::endl;
        }
    }
}

void iblsys(XFoilState &ctx) {
    int iv = 0;
    for (int is_ = 1; is_ <= 2; ++is_) {
        for (int ibl = 2; ibl <= ctx.NBL[is_]; ++ibl) {
            iv += 1;
            ctx.ISYS[ibl][is_] = iv;
        }
    }

    ctx.NSYS = iv;
    if (ctx.NSYS > 2 * IVX) {
        throw std::runtime_error("*** IBLSYS: BL system array overflow. ***");
    }
}

void mrchue(XFoilState &ctx, XBlState &bl) {
    const double hlmax = 3.8;
    const double htmax = 2.5;

    const int max_nbl = std::max(ctx.NBL[1], ctx.NBL[2]);
    if (ctx.HTARG.size() < 3 || ctx.HTARG[1].size() <= static_cast<size_t>(max_nbl)) {
        ctx.HTARG = std::vector<std::vector<double>>(3, std::vector<double>(static_cast<size_t>(max_nbl + 1), 0.0));
    } else {
        for (int is_ = 1; is_ <= 2; ++is_) {
            for (int ibl = 1; ibl <= max_nbl; ++ibl) {
                ctx.HTARG[is_][ibl] = 0.0;
            }
        }
    }

    for (int is_ = 1; is_ <= 2; ++is_) {
        std::cout << "   side " << is_ << " ..." << std::endl;

        bl.AMCRIT = ctx.ACRIT[is_];

        xifset(ctx, bl, is_);

        int ibl = 2;
        double xsi = ctx.XSSI[ibl][is_];
        double uei = ctx.UEDG[ibl][is_];
        bl.BULE = 1.0;
        const double ucon = uei / std::pow(xsi, bl.BULE);
        const double tsq = 0.45 / (ucon * (5.0 * bl.BULE + 1.0) * bl.REYBL) * std::pow(xsi, 1.0 - bl.BULE);
        double thi = std::sqrt(tsq);
        double dsi = 2.2 * thi;
        double ami = 0.0;

        double cti = 0.03;

        bl.TRAN = false;
        bl.TURB = false;
        ctx.ITRAN[is_] = ctx.IBLTE[is_];

        double sens = 0.0;
        for (ibl = 2; ibl <= ctx.NBL[is_]; ++ibl) {
            const int ibm = ibl - 1;

            bl.SIMI = (ibl == 2);
            bl.WAKE = (ibl > ctx.IBLTE[is_]);

            xsi = ctx.XSSI[ibl][is_];
            uei = ctx.UEDG[ibl][is_];

            double dswaki = 0.0;
            if (bl.WAKE) {
                const int iw = ibl - ctx.IBLTE[is_];
                dswaki = ctx.WGAP[iw];
            }

            bool direct = true;
            double htarg = 0.0;
            double dmax = 0.0;
            double sennew = sens;
            double tte = 0.0;
            double dte = 0.0;
            double cte = 0.0;

            int itbl = 0;
            for (itbl = 1; itbl <= 25; ++itbl) {
                sennew = sens;
                blprv(bl, xsi, ami, cti, thi, dsi, dswaki, uei);
                blkin(bl);

                if (!bl.SIMI && !bl.TURB) {
                    trchek(bl);
                    ami = bl.AMPL2;

                    if (bl.TRAN) {
                        ctx.ITRAN[is_] = ibl;
                        if (cti <= 0.0) {
                            cti = 0.03;
                            bl.S2 = cti;
                        }
                    } else {
                        ctx.ITRAN[is_] = ibl + 2;
                    }
                }

                if (ibl == ctx.IBLTE[is_] + 1) {
                    tte = ctx.THET[ctx.IBLTE[1]][1] + ctx.THET[ctx.IBLTE[2]][2];
                    dte = ctx.DSTR[ctx.IBLTE[1]][1] + ctx.DSTR[ctx.IBLTE[2]][2] + ctx.ANTE;
                    cte = (ctx.CTAU[ctx.IBLTE[1]][1] * ctx.THET[ctx.IBLTE[1]][1]
                           + ctx.CTAU[ctx.IBLTE[2]][2] * ctx.THET[ctx.IBLTE[2]][2])
                          / tte;
                    tesys(bl, cte, tte, dte);
                } else {
                    blsys(bl);
                }

                if (direct) {
                    bl.VS2[4][1] = 0.0;
                    bl.VS2[4][2] = 0.0;
                    bl.VS2[4][3] = 0.0;
                    bl.VS2[4][4] = 1.0;
                    bl.VSREZ[4] = 0.0;

                    gauss(4, 4, bl.VS2, bl.VSREZ, 1);

                    dmax = std::max(std::abs(bl.VSREZ[2] / thi), std::abs(bl.VSREZ[3] / dsi));
                    if (ibl < ctx.ITRAN[is_]) {
                        dmax = std::max(dmax, std::abs(bl.VSREZ[1] / 10.0));
                    }
                    if (ibl >= ctx.ITRAN[is_]) {
                        dmax = std::max(dmax, std::abs(bl.VSREZ[1] / cti));
                    }

                    double rlx = 1.0;
                    if (dmax > 0.3) {
                        rlx = 0.3 / dmax;
                    }

                    if (ibl != ctx.IBLTE[is_] + 1) {
                        const double msq = uei * uei * bl.HSTINV / (bl.GM1BL * (1.0 - 0.5 * uei * uei * bl.HSTINV));
                        const double htest = (dsi + rlx * bl.VSREZ[3]) / (thi + rlx * bl.VSREZ[2]);
                        double hktest = 0.0;
                        double dummy1 = 0.0;
                        double dummy2 = 0.0;
                        std::tie(hktest, dummy1, dummy2) = hkin(htest, msq);

                        double hmax = htmax;
                        if (ibl < ctx.ITRAN[is_]) {
                            hmax = hlmax;
                        }
                        if (ibl >= ctx.ITRAN[is_]) {
                            hmax = htmax;
                        }
                        direct = (hktest < hmax);
                    }

                    if (direct) {
                        if (ibl >= ctx.ITRAN[is_]) {
                            cti = cti + rlx * bl.VSREZ[1];
                        }
                        thi = thi + rlx * bl.VSREZ[2];
                        dsi = dsi + rlx * bl.VSREZ[3];
                    } else {
                        if (ibl < ctx.ITRAN[is_]) {
                            htarg = bl.HK1 + 0.03 * (bl.X2 - bl.X1) / bl.T1;
                        } else if (ibl == ctx.ITRAN[is_]) {
                            htarg = bl.HK1 + (0.03 * (bl.XT - bl.X1) - 0.15 * (bl.X2 - bl.XT)) / bl.T1;
                        } else if (bl.WAKE) {
                            const double cst = 0.03 * (bl.X2 - bl.X1) / bl.T1;
                            double hk2 = bl.HK1;
                            hk2 = hk2 - (hk2 + cst * std::pow(hk2 - 1.0, 3.0) - bl.HK1)
                                          / (1.0 + 3.0 * cst * std::pow(hk2 - 1.0, 2.0));
                            hk2 = hk2 - (hk2 + cst * std::pow(hk2 - 1.0, 3.0) - bl.HK1)
                                          / (1.0 + 3.0 * cst * std::pow(hk2 - 1.0, 2.0));
                            hk2 = hk2 - (hk2 + cst * std::pow(hk2 - 1.0, 3.0) - bl.HK1)
                                          / (1.0 + 3.0 * cst * std::pow(hk2 - 1.0, 2.0));
                            htarg = hk2;
                        } else {
                            htarg = bl.HK1 - 0.15 * (bl.X2 - bl.X1) / bl.T1;
                        }

                        if (bl.WAKE) {
                            htarg = std::max(htarg, 1.01);
                        } else {
                            htarg = std::max(htarg, (ibl < ctx.ITRAN[is_]) ? hlmax : htmax);
                        }

                        std::cout << " MRCHUE: Inverse mode at" << std::setw(4) << ibl << "     Hk =" << std::setw(8)
                                  << std::fixed << std::setprecision(3) << htarg << std::defaultfloat << std::endl;
                        ctx.HTARG[is_][ibl] = htarg;

                        continue;
                    }
                } else {
                    bl.VS2[4][1] = 0.0;
                    bl.VS2[4][2] = bl.HK2_T2;
                    bl.VS2[4][3] = bl.HK2_D2;
                    bl.VS2[4][4] = bl.HK2_U2;
                    bl.VSREZ[4] = htarg - bl.HK2;

                    gauss(4, 4, bl.VS2, bl.VSREZ, 1);

                    dmax = std::max(std::abs(bl.VSREZ[2] / thi), std::abs(bl.VSREZ[3] / dsi));
                    dmax = std::max(dmax, std::abs(bl.VSREZ[4] / uei));
                    if (ibl >= ctx.ITRAN[is_]) {
                        dmax = std::max(dmax, std::abs(bl.VSREZ[1] / cti));
                    }

                    double rlx = 1.0;
                    if (dmax > 0.3) {
                        rlx = 0.3 / dmax;
                    }

                    if (ibl >= ctx.ITRAN[is_]) {
                        cti = cti + rlx * bl.VSREZ[1];
                    }
                    thi = thi + rlx * bl.VSREZ[2];
                    dsi = dsi + rlx * bl.VSREZ[3];
                    uei = uei + rlx * bl.VSREZ[4];
                }

                if (ibl >= ctx.ITRAN[is_]) {
                    cti = std::min(cti, 0.30);
                    cti = std::max(cti, 0.0000001);
                }

                const double hklim = (ibl <= ctx.IBLTE[is_]) ? 1.02 : 1.00005;
                const double msq = uei * uei * bl.HSTINV / (bl.GM1BL * (1.0 - 0.5 * uei * uei * bl.HSTINV));
                double dsw = dsi - dswaki;
                dsw = dslim(dsw, thi, uei, msq, hklim);
                dsi = dsw + dswaki;

                if (dmax <= 1.0e-5) {
                    break;
                }
            }

            if (itbl > 25) {
                std::cout << " MRCHUE: Convergence failed at" << std::setw(4) << ibl << "  side" << std::setw(2) << is_
                          << "    Res =" << std::setw(12) << std::scientific << dmax << std::defaultfloat << std::endl;
                if (dmax > 0.1) {
                    if (ibl > 3) {
                        if (ibl <= ctx.IBLTE[is_]) {
                            thi = ctx.THET[ibm][is_] * std::sqrt(ctx.XSSI[ibl][is_] / ctx.XSSI[ibm][is_]);
                            dsi = ctx.DSTR[ibm][is_] * std::sqrt(ctx.XSSI[ibl][is_] / ctx.XSSI[ibm][is_]);
                        } else if (ibl == ctx.IBLTE[is_] + 1) {
                            cti = cte;
                            thi = tte;
                            dsi = dte;
                        } else {
                            thi = ctx.THET[ibm][is_];
                            const double ratlen = (ctx.XSSI[ibl][is_] - ctx.XSSI[ibm][is_]) / (10.0 * ctx.DSTR[ibm][is_]);
                            dsi = (ctx.DSTR[ibm][is_] + thi * ratlen) / (1.0 + ratlen);
                        }
                        if (ibl == ctx.ITRAN[is_]) {
                            cti = 0.05;
                        }
                        if (ibl > ctx.ITRAN[is_]) {
                            cti = ctx.CTAU[ibm][is_];
                        }

                        uei = ctx.UEDG[ibl][is_];
                        if (ibl > 2 && ibl < ctx.NBL[is_]) {
                            uei = 0.5 * (ctx.UEDG[ibl - 1][is_] + ctx.UEDG[ibl + 1][is_]);
                        }
                    }
                }

                blprv(bl, xsi, ami, cti, thi, dsi, dswaki, uei);
                blkin(bl);

                if (!bl.SIMI && !bl.TURB) {
                    trchek(bl);
                    ami = bl.AMPL2;
                    if (bl.TRAN) {
                        ctx.ITRAN[is_] = ibl;
                    }
                    if (!bl.TRAN) {
                        ctx.ITRAN[is_] = ibl + 2;
                    }
                }

                if (ibl < ctx.ITRAN[is_]) {
                    blvar(bl, 1);
                }
                if (ibl >= ctx.ITRAN[is_]) {
                    blvar(bl, 2);
                }
                if (bl.WAKE) {
                    blvar(bl, 3);
                }

                if (ibl < ctx.ITRAN[is_]) {
                    blmid(bl, 1);
                }
                if (ibl >= ctx.ITRAN[is_]) {
                    blmid(bl, 2);
                }
                if (bl.WAKE) {
                    blmid(bl, 3);
                }
            }

            sens = sennew;

            if (ibl < ctx.ITRAN[is_]) {
                ctx.CTAU[ibl][is_] = ami;
            }
            if (ibl >= ctx.ITRAN[is_]) {
                ctx.CTAU[ibl][is_] = cti;
            }
            ctx.THET[ibl][is_] = thi;
            ctx.DSTR[ibl][is_] = dsi;
            ctx.UEDG[ibl][is_] = uei;
            ctx.MASS[ibl][is_] = dsi * uei;
            ctx.TAU[ibl][is_] = 0.5 * bl.R2 * bl.U2 * bl.U2 * bl.CF2;
            ctx.DIS[ibl][is_] = bl.R2 * bl.U2 * bl.U2 * bl.U2 * bl.DI2 * bl.HS2 * 0.5;
            ctx.CTQ[ibl][is_] = bl.CQ2;
            ctx.DELT[ibl][is_] = bl.DE2;
            ctx.TSTR[ibl][is_] = bl.HS2 * bl.T2;

            blprv(bl, xsi, ami, cti, thi, dsi, dswaki, uei);
            blkin(bl);
            for (int icom = 1; icom <= NCOM; ++icom) {
                bl.COM1[icom] = bl.COM2[icom];
            }
            sync_com_to_vars(bl, 1);

            if (bl.TRAN || ibl == ctx.IBLTE[is_]) {
                bl.TURB = true;
                ctx.TFORCE[is_] = bl.TRFORC;
                ctx.XSSITR[is_] = bl.XT;
            }

            bl.TRAN = false;

            if (ibl == ctx.IBLTE[is_]) {
                thi = ctx.THET[ctx.IBLTE[1]][1] + ctx.THET[ctx.IBLTE[2]][2];
                dsi = ctx.DSTR[ctx.IBLTE[1]][1] + ctx.DSTR[ctx.IBLTE[2]][2] + ctx.ANTE;
            }
        }
    }
}

void mrchdu(XFoilState &ctx, XBlState &bl) {
    const double deps = 5.0e-6;
    const double senswt = 1000.0;

    for (int is_ = 1; is_ <= 2; ++is_) {
        bl.AMCRIT = ctx.ACRIT[is_];

        xifset(ctx, bl, is_);

        int ibl = 2;
        double xsi = ctx.XSSI[ibl][is_];
        double uei = ctx.UEDG[ibl][is_];
        bl.BULE = 1.0;

        const int itrold = ctx.ITRAN[is_];

        bl.TRAN = false;
        bl.TURB = false;
        ctx.ITRAN[is_] = ctx.IBLTE[is_];

        double sens = 0.0;
        double ami = 0.0;
        for (ibl = 2; ibl <= ctx.NBL[is_]; ++ibl) {
            const int ibm = ibl - 1;

            bl.SIMI = (ibl == 2);
            bl.WAKE = (ibl > ctx.IBLTE[is_]);

            xsi = ctx.XSSI[ibl][is_];
            uei = ctx.UEDG[ibl][is_];
            double thi = ctx.THET[ibl][is_];
            double dsi = ctx.DSTR[ibl][is_];

            double cti = 0.0;
            if (ibl < itrold) {
                ami = ctx.CTAU[ibl][is_];
                cti = 0.03;
            } else {
                cti = ctx.CTAU[ibl][is_];
                if (cti <= 0.0) {
                    cti = 0.03;
                }
            }

            double dswaki = 0.0;
            if (bl.WAKE) {
                const int iw = ibl - ctx.IBLTE[is_];
                dswaki = ctx.WGAP[iw];
            }

            if (ibl <= ctx.IBLTE[is_]) {
                dsi = std::max(dsi - dswaki, 1.02000 * thi) + dswaki;
            }
            if (ibl > ctx.IBLTE[is_]) {
                dsi = std::max(dsi - dswaki, 1.00005 * thi) + dswaki;
            }

            double sennew = sens;
            double ueref = 0.0;
            double hkref = 0.0;
            double tte = 0.0;
            double dte = 0.0;
            double cte = 0.0;
            double dmax = 0.0;

            int itbl = 0;
            for (itbl = 1; itbl <= 25; ++itbl) {
                blprv(bl, xsi, ami, cti, thi, dsi, dswaki, uei);
                blkin(bl);

                if (!bl.SIMI && !bl.TURB) {
                    trchek(bl);
                    ami = bl.AMPL2;
                    if (bl.TRAN) {
                        ctx.ITRAN[is_] = ibl;
                    }
                    if (!bl.TRAN) {
                        ctx.ITRAN[is_] = ibl + 2;
                    }
                }

                if (ibl == ctx.IBLTE[is_] + 1) {
                    tte = ctx.THET[ctx.IBLTE[1]][1] + ctx.THET[ctx.IBLTE[2]][2];
                    dte = ctx.DSTR[ctx.IBLTE[1]][1] + ctx.DSTR[ctx.IBLTE[2]][2] + ctx.ANTE;
                    cte = (ctx.CTAU[ctx.IBLTE[1]][1] * ctx.THET[ctx.IBLTE[1]][1]
                           + ctx.CTAU[ctx.IBLTE[2]][2] * ctx.THET[ctx.IBLTE[2]][2])
                          / tte;
                    tesys(bl, cte, tte, dte);
                } else {
                    blsys(bl);
                }

                if (itbl == 1) {
                    ueref = bl.U2;
                    hkref = bl.HK2;

                    if (ibl < ctx.ITRAN[is_] && ibl >= itrold) {
                        const double uem = ctx.UEDG[ibl - 1][is_];
                        const double dsm = ctx.DSTR[ibl - 1][is_];
                        const double thm = ctx.THET[ibl - 1][is_];
                        const double msq = uem * uem * bl.HSTINV / (bl.GM1BL * (1.0 - 0.5 * uem * uem * bl.HSTINV));
                        double dummy1 = 0.0;
                        double dummy2 = 0.0;
                        std::tie(hkref, dummy1, dummy2) = hkin(dsm / thm, msq);
                    }

                    if (ibl < itrold) {
                        if (bl.TRAN) {
                            ctx.CTAU[ibl][is_] = 0.03;
                        }
                        if (bl.TURB) {
                            ctx.CTAU[ibl][is_] = ctx.CTAU[ibl - 1][is_];
                        }
                        if (bl.TRAN || bl.TURB) {
                            cti = ctx.CTAU[ibl][is_];
                            bl.S2 = cti;
                        }
                    }
                }

                if (bl.SIMI || ibl == ctx.IBLTE[is_] + 1) {
                    bl.VS2[4][1] = 0.0;
                    bl.VS2[4][2] = 0.0;
                    bl.VS2[4][3] = 0.0;
                    bl.VS2[4][4] = bl.U2_UEI;
                    bl.VSREZ[4] = ueref - bl.U2;
                } else {
                    std::vector<std::vector<double>> vtmp = make_2d(4, 5);
                    std::vector<double> vztmp = make_1d(4);
                    for (int k = 1; k <= 4; ++k) {
                        vztmp[k] = bl.VSREZ[k];
                        for (int l = 1; l <= 5; ++l) {
                            vtmp[k][l] = bl.VS2[k][l];
                        }
                    }

                    vtmp[4][1] = 0.0;
                    vtmp[4][2] = bl.HK2_T2;
                    vtmp[4][3] = bl.HK2_D2;
                    vtmp[4][4] = bl.HK2_U2 * bl.U2_UEI;
                    vztmp[4] = 1.0;

                    gauss(4, 4, vtmp, vztmp, 1);

                    sennew = senswt * vztmp[4] * hkref / ueref;
                    if (itbl <= 5) {
                        sens = sennew;
                    } else if (itbl <= 15) {
                        sens = 0.5 * (sens + sennew);
                    }

                    bl.VS2[4][1] = 0.0;
                    bl.VS2[4][2] = bl.HK2_T2 * hkref;
                    bl.VS2[4][3] = bl.HK2_D2 * hkref;
                    bl.VS2[4][4] = (bl.HK2_U2 * hkref + sens / ueref) * bl.U2_UEI;
                    bl.VSREZ[4] = -(hkref * hkref) * (bl.HK2 / hkref - 1.0) - sens * (bl.U2 / ueref - 1.0);
                }

                gauss(4, 4, bl.VS2, bl.VSREZ, 1);

                dmax = std::max(std::abs(bl.VSREZ[2] / thi), std::abs(bl.VSREZ[3] / dsi));
                dmax = std::max(dmax, std::abs(bl.VSREZ[4] / uei));
                if (ibl >= ctx.ITRAN[is_]) {
                    dmax = std::max(dmax, std::abs(bl.VSREZ[1] / (10.0 * cti)));
                }

                double rlx = 1.0;
                if (dmax > 0.3) {
                    rlx = 0.3 / dmax;
                }

                if (ibl < ctx.ITRAN[is_]) {
                    ami = ami + rlx * bl.VSREZ[1];
                }
                if (ibl >= ctx.ITRAN[is_]) {
                    cti = cti + rlx * bl.VSREZ[1];
                }
                thi = thi + rlx * bl.VSREZ[2];
                dsi = dsi + rlx * bl.VSREZ[3];
                uei = uei + rlx * bl.VSREZ[4];

                if (ibl >= ctx.ITRAN[is_]) {
                    cti = std::min(cti, 0.30);
                    cti = std::max(cti, 0.0000001);
                }

                const double hklim = (ibl <= ctx.IBLTE[is_]) ? 1.02 : 1.00005;
                const double msq = uei * uei * bl.HSTINV / (bl.GM1BL * (1.0 - 0.5 * uei * uei * bl.HSTINV));
                double dsw = dsi - dswaki;
                dsw = dslim(dsw, thi, uei, msq, hklim);
                dsi = dsw + dswaki;

                if (dmax <= deps) {
                    break;
                }
            }

            if (itbl > 25) {
                std::cout << " MRCHDU: Convergence failed at" << std::setw(4) << ibl << "  side" << std::setw(2) << is_
                          << "    Res =" << std::setw(12) << std::scientific << dmax << std::defaultfloat << std::endl;
                if (dmax > 0.1) {
                    if (ibl > 3) {
                        if (ibl <= ctx.IBLTE[is_]) {
                            thi = ctx.THET[ibm][is_] * std::sqrt(ctx.XSSI[ibl][is_] / ctx.XSSI[ibm][is_]);
                            dsi = ctx.DSTR[ibm][is_] * std::sqrt(ctx.XSSI[ibl][is_] / ctx.XSSI[ibm][is_]);
                            uei = ctx.UEDG[ibm][is_];
                        } else if (ibl == ctx.IBLTE[is_] + 1) {
                            cti = cte;
                            thi = tte;
                            dsi = dte;
                            uei = ctx.UEDG[ibm][is_];
                        } else {
                            thi = ctx.THET[ibm][is_];
                            const double ratlen = (ctx.XSSI[ibl][is_] - ctx.XSSI[ibm][is_]) / (10.0 * ctx.DSTR[ibm][is_]);
                            dsi = (ctx.DSTR[ibm][is_] + thi * ratlen) / (1.0 + ratlen);
                            uei = ctx.UEDG[ibm][is_];
                        }
                        if (ibl == ctx.ITRAN[is_]) {
                            cti = 0.05;
                        }
                        if (ibl > ctx.ITRAN[is_]) {
                            cti = ctx.CTAU[ibm][is_];
                        }
                    }
                }

                blprv(bl, xsi, ami, cti, thi, dsi, dswaki, uei);
                blkin(bl);

                if (!bl.SIMI && !bl.TURB) {
                    trchek(bl);
                    ami = bl.AMPL2;
                    if (bl.TRAN) {
                        ctx.ITRAN[is_] = ibl;
                    }
                    if (!bl.TRAN) {
                        ctx.ITRAN[is_] = ibl + 2;
                    }
                }

                if (ibl < ctx.ITRAN[is_]) {
                    blvar(bl, 1);
                }
                if (ibl >= ctx.ITRAN[is_]) {
                    blvar(bl, 2);
                }
                if (bl.WAKE) {
                    blvar(bl, 3);
                }

                if (ibl < ctx.ITRAN[is_]) {
                    blmid(bl, 1);
                }
                if (ibl >= ctx.ITRAN[is_]) {
                    blmid(bl, 2);
                }
                if (bl.WAKE) {
                    blmid(bl, 3);
                }
            }

            sens = sennew;

            if (ibl < ctx.ITRAN[is_]) {
                ctx.CTAU[ibl][is_] = ami;
            }
            if (ibl >= ctx.ITRAN[is_]) {
                ctx.CTAU[ibl][is_] = cti;
            }
            ctx.THET[ibl][is_] = thi;
            ctx.DSTR[ibl][is_] = dsi;
            ctx.UEDG[ibl][is_] = uei;
            ctx.MASS[ibl][is_] = dsi * uei;
            ctx.TAU[ibl][is_] = 0.5 * bl.R2 * bl.U2 * bl.U2 * bl.CF2;
            ctx.DIS[ibl][is_] = bl.R2 * bl.U2 * bl.U2 * bl.U2 * bl.DI2 * bl.HS2 * 0.5;
            ctx.CTQ[ibl][is_] = bl.CQ2;
            ctx.DELT[ibl][is_] = bl.DE2;
            ctx.TSTR[ibl][is_] = bl.HS2 * bl.T2;

            blprv(bl, xsi, ami, cti, thi, dsi, dswaki, uei);
            blkin(bl);
            for (int icom = 1; icom <= NCOM; ++icom) {
                bl.COM1[icom] = bl.COM2[icom];
            }
            sync_com_to_vars(bl, 1);

            if (bl.TRAN || ibl == ctx.IBLTE[is_]) {
                bl.TURB = true;
                ctx.TFORCE[is_] = bl.TRFORC;
                ctx.XSSITR[is_] = bl.XT;
            }

            bl.TRAN = false;
        }
    }
}

void xifset(XFoilState &ctx, XBlState &bl, int is_) {
    if (ctx.XSTRIP[is_] >= 1.0) {
        bl.XIFORC = ctx.XSSI[ctx.IBLTE[is_]][is_];
        return;
    }

    const double chx = ctx.XTE - ctx.XLE;
    const double chy = ctx.YTE - ctx.YLE;
    const double chsq = chx * chx + chy * chy;

    for (int i = 1; i <= ctx.N; ++i) {
        ctx.W1[i] = ((ctx.X[i] - ctx.XLE) * chx + (ctx.Y[i] - ctx.YLE) * chy) / chsq;
        ctx.W2[i] = ((ctx.Y[i] - ctx.YLE) * chx - (ctx.X[i] - ctx.XLE) * chy) / chsq;
    }

    splind(ctx.W1, ctx.W3, ctx.S, ctx.N, -999.0, -999.0);
    splind(ctx.W2, ctx.W4, ctx.S, ctx.N, -999.0, -999.0);

    double str_ = 0.0;
    if (is_ == 1) {
        str_ = ctx.SLE + (ctx.S[1] - ctx.SLE) * ctx.XSTRIP[is_];
        str_ = sinvrt(str_, ctx.XSTRIP[is_], ctx.W1, ctx.W3, ctx.S, ctx.N);
        bl.XIFORC = std::min((ctx.SST - str_), ctx.XSSI[ctx.IBLTE[is_]][is_]);
    } else {
        str_ = ctx.SLE + (ctx.S[ctx.N] - ctx.SLE) * ctx.XSTRIP[is_];
        str_ = sinvrt(str_, ctx.XSTRIP[is_], ctx.W1, ctx.W3, ctx.S, ctx.N);
        bl.XIFORC = std::min((str_ - ctx.SST), ctx.XSSI[ctx.IBLTE[is_]][is_]);
    }

    if (bl.XIFORC < 0.0) {
        std::cout << std::endl;
        std::cout << " ***  Stagnation point is past trip on side" << std::setw(2) << is_ << "  ***" << std::endl;
        bl.XIFORC = ctx.XSSI[ctx.IBLTE[is_]][is_];
    }
}

void update(XFoilState &ctx, XBlState &bl) {
    std::vector<std::vector<double>> unew = make_2d(IVX, 2);
    std::vector<std::vector<double>> u_ac = make_2d(IVX, 2);
    std::vector<double> qnew = make_1d(IQX);
    std::vector<double> q_ac = make_1d(IQX);

    const double dalmax = 0.5 * ctx.DTOR;
    const double dalmin = -0.5 * ctx.DTOR;

    const double dclmax = 0.5;
    double dclmin = -0.5;
    if (ctx.MATYP != 1) {
        dclmin = std::max(-0.5, -0.9 * ctx.CL);
    }

    const double hstinv = ctx.GAMM1 * (ctx.MINF / ctx.QINF) * (ctx.MINF / ctx.QINF)
                          / (1.0 + 0.5 * ctx.GAMM1 * ctx.MINF * ctx.MINF);

    for (int is_ = 1; is_ <= 2; ++is_) {
        for (int ibl = 2; ibl <= ctx.NBL[is_]; ++ibl) {
            const int i = ctx.IPAN[ibl][is_];

            double dui = 0.0;
            double dui_ac = 0.0;
            for (int js = 1; js <= 2; ++js) {
                for (int jbl = 2; jbl <= ctx.NBL[js]; ++jbl) {
                    const int j = ctx.IPAN[jbl][js];
                    const int jv = ctx.ISYS[jbl][js];
                    const double ue_m = -ctx.VTI[ibl][is_] * ctx.VTI[jbl][js] * ctx.DIJ[i][j];
                    dui = dui + ue_m * (ctx.MASS[jbl][js] + ctx.VDEL[3][1][jv]);
                    dui_ac = dui_ac + ue_m * (0.0 - ctx.VDEL[3][2][jv]);
                }
            }

            const double uinv_ac = ctx.LALFA ? 0.0 : ctx.UINV_A[ibl][is_];

            unew[ibl][is_] = ctx.UINV[ibl][is_] + dui;
            u_ac[ibl][is_] = uinv_ac + dui_ac;
        }
    }

    for (int is_ = 1; is_ <= 2; ++is_) {
        for (int ibl = 2; ibl <= ctx.IBLTE[is_]; ++ibl) {
            const int i = ctx.IPAN[ibl][is_];
            qnew[i] = ctx.VTI[ibl][is_] * unew[ibl][is_];
            q_ac[i] = ctx.VTI[ibl][is_] * u_ac[ibl][is_];
        }
    }

    const double sa = std::sin(ctx.ALFA);
    const double ca = std::cos(ctx.ALFA);

    const double beta = std::sqrt(1.0 - ctx.MINF * ctx.MINF);
    const double beta_msq = -0.5 / beta;

    const double bfac = 0.5 * ctx.MINF * ctx.MINF / (1.0 + beta);
    const double bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

    double clnew = 0.0;
    double cl_a = 0.0;
    double cl_ms = 0.0;
    double cl_ac = 0.0;

    int i = 1;
    double cginc = 1.0 - (qnew[i] / ctx.QINF) * (qnew[i] / ctx.QINF);
    double cpg1 = cginc / (beta + bfac * cginc);
    double cpg1_ms = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

    double cpi_q = -2.0 * qnew[i] / (ctx.QINF * ctx.QINF);
    double cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
    double cpg1_ac = cpc_cpi * cpi_q * q_ac[i];

    for (i = 1; i <= ctx.N; ++i) {
        int ip = i + 1;
        if (i == ctx.N) {
            ip = 1;
        }

        cginc = 1.0 - (qnew[ip] / ctx.QINF) * (qnew[ip] / ctx.QINF);
        const double cpg2 = cginc / (beta + bfac * cginc);
        const double cpg2_ms = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

        cpi_q = -2.0 * qnew[ip] / (ctx.QINF * ctx.QINF);
        cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
        const double cpg2_ac = cpc_cpi * cpi_q * q_ac[ip];

        const double dx = (ctx.X[ip] - ctx.X[i]) * ca + (ctx.Y[ip] - ctx.Y[i]) * sa;
        const double dx_a = -(ctx.X[ip] - ctx.X[i]) * sa + (ctx.Y[ip] - ctx.Y[i]) * ca;

        const double ag = 0.5 * (cpg2 + cpg1);
        const double ag_ms = 0.5 * (cpg2_ms + cpg1_ms);
        const double ag_ac = 0.5 * (cpg2_ac + cpg1_ac);

        clnew = clnew + dx * ag;
        cl_a = cl_a + dx_a * ag;
        cl_ms = cl_ms + dx * ag_ms;
        cl_ac = cl_ac + dx * ag_ac;

        cpg1 = cpg2;
        cpg1_ms = cpg2_ms;
        cpg1_ac = cpg2_ac;
    }

    double rlx = 1.0;
    double dac = 0.0;

    if (ctx.LALFA) {
        dac = (clnew - ctx.CL) / (1.0 - cl_ac - cl_ms * 2.0 * ctx.MINF * ctx.MINF_CL);

        if (rlx * dac > dclmax) {
            rlx = dclmax / dac;
        }
        if (rlx * dac < dclmin) {
            rlx = dclmin / dac;
        }
    } else {
        dac = (clnew - ctx.CLSPEC) / (0.0 - cl_ac - cl_a);

        if (rlx * dac > dalmax) {
            rlx = dalmax / dac;
        }
        if (rlx * dac < dalmin) {
            rlx = dalmin / dac;
        }
    }

    ctx.RMSBL = 0.0;
    ctx.RMXBL = 0.0;

    const double dhi = 1.5;
    const double dlo = -0.5;

    for (int is_ = 1; is_ <= 2; ++is_) {
        for (int ibl = 2; ibl <= ctx.NBL[is_]; ++ibl) {
            const int iv = ctx.ISYS[ibl][is_];

            const double dctau = ctx.VDEL[1][1][iv] - dac * ctx.VDEL[1][2][iv];
            const double dthet = ctx.VDEL[2][1][iv] - dac * ctx.VDEL[2][2][iv];
            const double dmass = ctx.VDEL[3][1][iv] - dac * ctx.VDEL[3][2][iv];
            const double duedg = unew[ibl][is_] + dac * u_ac[ibl][is_] - ctx.UEDG[ibl][is_];
            const double ddstr = (dmass - ctx.DSTR[ibl][is_] * duedg) / ctx.UEDG[ibl][is_];

            double dn1 = 0.0;
            if (ibl < ctx.ITRAN[is_]) {
                dn1 = dctau / 10.0;
            }
            if (ibl >= ctx.ITRAN[is_]) {
                dn1 = dctau / ctx.CTAU[ibl][is_];
            }
            const double dn2 = dthet / ctx.THET[ibl][is_];
            const double dn3 = ddstr / ctx.DSTR[ibl][is_];
            const double dn4 = std::abs(duedg) / 0.25;

            ctx.RMSBL = ctx.RMSBL + dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;

            const double rdn1 = rlx * dn1;
            if (std::abs(dn1) > std::abs(ctx.RMXBL)) {
                ctx.RMXBL = dn1;
                ctx.VMXBL = (ibl < ctx.ITRAN[is_]) ? "n" : "C";
                ctx.IMXBL = ibl;
                ctx.ISMXBL = is_;
            }
            if (rdn1 > dhi) {
                rlx = dhi / dn1;
            }
            if (rdn1 < dlo) {
                rlx = dlo / dn1;
            }

            const double rdn2 = rlx * dn2;
            if (std::abs(dn2) > std::abs(ctx.RMXBL)) {
                ctx.RMXBL = dn2;
                ctx.VMXBL = "T";
                ctx.IMXBL = ibl;
                ctx.ISMXBL = is_;
            }
            if (rdn2 > dhi) {
                rlx = dhi / dn2;
            }
            if (rdn2 < dlo) {
                rlx = dlo / dn2;
            }

            const double rdn3 = rlx * dn3;
            if (std::abs(dn3) > std::abs(ctx.RMXBL)) {
                ctx.RMXBL = dn3;
                ctx.VMXBL = "D";
                ctx.IMXBL = ibl;
                ctx.ISMXBL = is_;
            }
            if (rdn3 > dhi) {
                rlx = dhi / dn3;
            }
            if (rdn3 < dlo) {
                rlx = dlo / dn3;
            }

            const double rdn4 = rlx * dn4;
            if (std::abs(dn4) > std::abs(ctx.RMXBL)) {
                ctx.RMXBL = duedg;
                ctx.VMXBL = "U";
                ctx.IMXBL = ibl;
                ctx.ISMXBL = is_;
            }
            if (rdn4 > dhi) {
                rlx = dhi / dn4;
            }
            if (rdn4 < dlo) {
                rlx = dlo / dn4;
            }
        }
    }

    ctx.RMSBL = std::sqrt(ctx.RMSBL / (4.0 * static_cast<double>(ctx.NBL[1] + ctx.NBL[2])));
    ctx.RLX = rlx;

    if (ctx.LALFA) {
        ctx.CL = ctx.CL + rlx * dac;
    } else {
        ctx.ALFA = ctx.ALFA + rlx * dac;
        ctx.ADEG = ctx.ALFA / ctx.DTOR;
    }

    for (int is_ = 1; is_ <= 2; ++is_) {
        for (int ibl = 2; ibl <= ctx.NBL[is_]; ++ibl) {
            const int iv = ctx.ISYS[ibl][is_];

            const double dctau = ctx.VDEL[1][1][iv] - dac * ctx.VDEL[1][2][iv];
            const double dthet = ctx.VDEL[2][1][iv] - dac * ctx.VDEL[2][2][iv];
            const double dmass = ctx.VDEL[3][1][iv] - dac * ctx.VDEL[3][2][iv];
            const double duedg = unew[ibl][is_] + dac * u_ac[ibl][is_] - ctx.UEDG[ibl][is_];
            const double ddstr = (dmass - ctx.DSTR[ibl][is_] * duedg) / ctx.UEDG[ibl][is_];

            ctx.CTAU[ibl][is_] = ctx.CTAU[ibl][is_] + rlx * dctau;
            ctx.THET[ibl][is_] = ctx.THET[ibl][is_] + rlx * dthet;
            ctx.DSTR[ibl][is_] = ctx.DSTR[ibl][is_] + rlx * ddstr;
            ctx.UEDG[ibl][is_] = ctx.UEDG[ibl][is_] + rlx * duedg;

            double dswaki = 0.0;
            if (ibl > ctx.IBLTE[is_]) {
                const int iw = ibl - ctx.IBLTE[is_];
                dswaki = ctx.WGAP[iw];
            }

            if (ibl >= ctx.ITRAN[is_]) {
                ctx.CTAU[ibl][is_] = std::min(ctx.CTAU[ibl][is_], 0.25);
            }

            double hklim = 0.0;
            if (ibl <= ctx.IBLTE[is_]) {
                hklim = 1.02;
            } else {
                hklim = 1.00005;
            }
            const double msq = ctx.UEDG[ibl][is_] * ctx.UEDG[ibl][is_] * hstinv
                               / (ctx.GAMM1 * (1.0 - 0.5 * ctx.UEDG[ibl][is_] * ctx.UEDG[ibl][is_] * hstinv));
            double dsw = ctx.DSTR[ibl][is_] - dswaki;
            dsw = dslim(dsw, ctx.THET[ibl][is_], ctx.UEDG[ibl][is_], msq, hklim);
            ctx.DSTR[ibl][is_] = dsw + dswaki;

            ctx.MASS[ibl][is_] = ctx.DSTR[ibl][is_] * ctx.UEDG[ibl][is_];
        }

        for (int ibl = 3; ibl <= ctx.IBLTE[is_]; ++ibl) {
            if (ctx.UEDG[ibl - 1][is_] > 0.0 && ctx.UEDG[ibl][is_] <= 0.0) {
                ctx.UEDG[ibl][is_] = ctx.UEDG[ibl - 1][is_];
                ctx.MASS[ibl][is_] = ctx.DSTR[ibl][is_] * ctx.UEDG[ibl][is_];
            }
        }
    }

    for (int kbl = 1; kbl <= ctx.NBL[2] - ctx.IBLTE[2]; ++kbl) {
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

double dslim(double dstr, double thet, double /*uedg*/, double msq, double hklim) {
    double hk = 0.0;
    double hk_h = 0.0;
    double dummy = 0.0;
    std::tie(hk, hk_h, dummy) = hkin(dstr / thet, msq);
    const double dh = std::max(0.0, hklim - hk) / hk_h;
    dstr = dstr + dh * thet;
    return dstr;
}

void blpini(XBlState &bl) {
    bl.SCCON = 5.6;
    bl.GACON = 6.70;
    bl.GBCON = 0.75;
    bl.GCCON = 18.0;
    bl.DLCON = 0.9;

    bl.CTRCON = 1.8;
    bl.CTRCEX = 3.3;

    bl.DUXCON = 1.0;

    bl.CTCON = 0.5 / (bl.GACON * bl.GACON * bl.GBCON);

    bl.CFFAC = 1.0;
}
