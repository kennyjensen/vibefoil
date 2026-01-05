// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

#ifndef WASM_XBL_H
#define WASM_XBL_H

#include <string>
#include <vector>

constexpr int IQX = 370;
constexpr int IPX = 5;
constexpr int ISX = 2;
constexpr int IWX = IQX / 8 + 2;
constexpr int IBX = 4 * IQX;
constexpr int IZX = IQX + IWX;
constexpr int IVX = IQX / 2 + IWX + 50;

inline std::vector<double> make_1d(int n, double fill = 0.0) {
    return std::vector<double>(static_cast<size_t>(n + 1), fill);
}

inline std::vector<std::vector<double>> make_2d(int n1, int n2, double fill = 0.0) {
    return std::vector<std::vector<double>>(static_cast<size_t>(n1 + 1), std::vector<double>(static_cast<size_t>(n2 + 1), fill));
}

inline std::vector<std::vector<std::vector<double>>> make_3d(int n1, int n2, int n3, double fill = 0.0) {
    return std::vector<std::vector<std::vector<double>>>(
        static_cast<size_t>(n1 + 1),
        std::vector<std::vector<double>>(static_cast<size_t>(n2 + 1), std::vector<double>(static_cast<size_t>(n3 + 1), fill)));
}

struct XFoilState {
    XFoilState();

    bool LALFA;
    bool LBLINI;

    double CL;
    double CM;
    double CD;
    double CDP;
    double CDF;
    double CL_ALF;
    double CL_MSQ;
    double CLSPEC;
    double MINF;
    double MINF1;
    double MINF_CL;
    double QINF;
    double TKLAM;
    double TKL_MSQ;
    double GAMMA;
    double GAMM1;
    double HVRAT;
    double REINF;
    double REINF1;
    double REINF_CL;
    int RETYP;
    int MATYP;
    double CPSTAR;
    double QSTAR;
    int IDAMP;

    double SLE;
    double XLE;
    double YLE;
    double XTE;
    double YTE;
    double SST;
    double SST_GO;
    double SST_GP;

    double ALFA;
    double ADEG;
    double COSA;
    double SINA;
    double DTOR;
    double AVISC;
    double MVISC;
    double XCMREF;
    double YCMREF;
    double PSIO;
    double CIRC;

    int N;
    int IST;
    int NW;
    int NPAN;
    int NB;
    bool LCLOCK;

    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> XP;
    std::vector<double> YP;
    std::vector<double> S;
    std::vector<double> SNEW;

    std::vector<double> W1;
    std::vector<double> W2;
    std::vector<double> W3;
    std::vector<double> W4;
    std::vector<double> W5;
    std::vector<double> W6;
    std::vector<double> W7;
    std::vector<double> W8;

    std::vector<double> XB;
    std::vector<double> YB;
    std::vector<double> XBP;
    std::vector<double> YBP;
    std::vector<double> SB;
    double SBLE;
    double CHORDB;
    double AREAB;
    double RADBLE;
    double ANGBTE;
    double EI11BA;
    double EI22BA;
    double APX1BA;
    double APX2BA;
    double EI11BT;
    double EI22BT;
    double APX1BT;
    double APX2BT;
    double THICKB;
    double CAMBRB;

    std::vector<std::vector<double>> XSSI;
    std::vector<std::vector<double>> UEDG;
    std::vector<std::vector<double>> UINV;
    std::vector<std::vector<double>> UINV_A;
    std::vector<std::vector<double>> MASS;
    std::vector<std::vector<double>> THET;
    std::vector<std::vector<double>> DSTR;
    std::vector<std::vector<double>> CTAU;
    std::vector<std::vector<double>> DELT;
    std::vector<std::vector<double>> TSTR;
    std::vector<std::vector<double>> USLP;
    std::vector<std::vector<double>> GUXQ;
    std::vector<std::vector<double>> GUXD;
    std::vector<std::vector<double>> TAU;
    std::vector<std::vector<double>> DIS;
    std::vector<std::vector<double>> CTQ;
    std::vector<std::vector<double>> VTI;

    std::vector<double> ACRIT;
    std::vector<double> XSTRIP;
    std::vector<double> XOCTR;
    std::vector<double> YOCTR;
    std::vector<double> XSSITR;
    std::vector<double> TINDEX;

    std::vector<int> IBLTE;
    std::vector<int> NBL;
    std::vector<std::vector<int>> IPAN;
    std::vector<std::vector<int>> ISYS;
    int NSYS;
    std::vector<int> ITRAN;
    std::vector<bool> TFORCE;

    std::vector<double> WGAP;
    double DWTE;
    double ANTE;
    double DSTE;
    double ASTE;
    double WAKLEN;
    double CHORD;
    double YIMAGE;
    bool SHARP;

    std::vector<std::vector<double>> DIJ;
    std::vector<std::vector<double>> AIJ;
    std::vector<std::vector<double>> BIJ;
    std::vector<std::vector<double>> CIJ;
    std::vector<std::vector<std::vector<double>>> VM;
    std::vector<std::vector<std::vector<double>>> VA;
    std::vector<std::vector<std::vector<double>>> VB;
    std::vector<std::vector<std::vector<double>>> VDEL;
    std::vector<std::vector<double>> VZ;
    std::vector<int> AIJPIV;
    std::vector<std::vector<double>> Q;
    std::vector<double> DQ;

    std::vector<double> QINV;
    std::vector<double> QVIS;
    std::vector<double> CPI;
    std::vector<double> CPV;
    std::vector<std::vector<double>> QINVU;
    std::vector<double> QINV_A;
    double SSPLE;
    std::vector<double> SSPEC;
    std::vector<double> XSPOC;
    std::vector<double> YSPOC;
    std::vector<double> QGAMM;
    std::vector<std::vector<double>> QSPEC;
    std::vector<std::vector<double>> QSPECP;
    double ALGAM;
    double CLGAM;
    double CMGAM;
    std::vector<double> ALQSP;
    std::vector<double> CLQSP;
    std::vector<double> CMQSP;
    double QDOF0;
    double QDOF1;
    double QDOF2;
    double QDOF3;
    double CLSPEC;
    double FFILT;

    std::vector<double> GAM;
    std::vector<std::vector<double>> GAMU;
    std::vector<double> GAM_A;
    std::vector<double> SIG;
    std::vector<double> NX;
    std::vector<double> NY;
    std::vector<double> APANEL;
    double GAMTE;
    double GAMTE_A;
    double SIGTE;
    double SIGTE_A;

    std::vector<double> DZDG;
    std::vector<double> DZDN;
    std::vector<double> DZDM;
    std::vector<double> DQDG;
    std::vector<double> DQDM;
    double QTAN1;
    double QTAN2;
    double Z_QINF;
    double Z_ALFA;
    double Z_QDOF0;
    double Z_QDOF1;
    double Z_QDOF2;
    double Z_QDOF3;

    std::vector<double> QF0;
    std::vector<double> QF1;
    std::vector<double> QF2;
    std::vector<double> QF3;

    double PI;
    double HOPI;
    double QOPI;

    double RMSBL;
    double RMXBL;
    double RLX;
    double VACCEL;
    int IMXBL;
    int ISMXBL;
    std::string VMXBL;

    bool LQAIJ;
    bool LADIJ;
    bool LWDIJ;
    bool LWAKE;
    bool LGAMU;
    bool LVISC;
    bool LVCONV;
    bool LFLAP;
    bool LIMAGE;
    bool LQINU;
    bool LQSPEC;
    bool LGSAME;
    bool LSCINI;
    bool LCPXX;
    bool LQVDES;
    bool LQREFL;
    bool LQSYM;
    bool LQSLOP;
    bool LQSPPL;
    bool LIQSET;
    bool LIPAN;
    double AWAKE;

    double XOF;
    double YOF;
    double HMOM;
    double HFX;
    double HFY;

    double CVPAR;
    double CTERAT;
    double CTRRAT;
    double XSREF1;
    double XSREF2;
    double XPREF1;
    double XPREF2;

    double XBF;
    double YBF;
    bool LBFLAP;

    std::string NAME;
    int NNAME;

    std::vector<std::vector<double>> HTARG;

    int IQ1;
    int IQ2;
    int NSP;
    int NQSP;
    int KQTARG;
    int IACQSP;
    int NC1;
};

struct XBlState {
    XBlState();

    std::vector<double> COM1;
    std::vector<double> COM2;
    std::vector<double> C1SAV;
    std::vector<double> C2SAV;

    bool SIMI;
    bool TRAN;
    bool TURB;
    bool WAKE;
    bool TRFORC;
    bool TRFREE;
    int IDAMPV;

    std::vector<std::vector<double>> VS1;
    std::vector<std::vector<double>> VS2;
    std::vector<double> VSREZ;
    std::vector<double> VSR;
    std::vector<double> VSM;
    std::vector<double> VSX;

    double SCCON;
    double GACON;
    double GBCON;
    double GCCON;
    double DLCON;
    double CTRCON;
    double CTRCEX;
    double DUXCON;
    double CTCON;
    double CFFAC;

    double X1;
    double U1;
    double T1;
    double D1;
    double S1;
    double AMPL1;
    double U1_UEI;
    double U1_MS;
    double DW1;
    double H1;
    double H1_T1;
    double H1_D1;
    double M1;
    double M1_U1;
    double M1_MS;
    double R1;
    double R1_U1;
    double R1_MS;
    double V1;
    double V1_U1;
    double V1_MS;
    double V1_RE;
    double HK1;
    double HK1_U1;
    double HK1_T1;
    double HK1_D1;
    double HK1_MS;
    double HS1;
    double HS1_U1;
    double HS1_T1;
    double HS1_D1;
    double HS1_MS;
    double HS1_RE;
    double HC1;
    double HC1_U1;
    double HC1_T1;
    double HC1_D1;
    double HC1_MS;
    double RT1;
    double RT1_U1;
    double RT1_T1;
    double RT1_MS;
    double RT1_RE;
    double CF1;
    double CF1_U1;
    double CF1_T1;
    double CF1_D1;
    double CF1_MS;
    double CF1_RE;
    double DI1;
    double DI1_U1;
    double DI1_T1;
    double DI1_D1;
    double DI1_S1;
    double DI1_MS;
    double DI1_RE;
    double US1;
    double US1_U1;
    double US1_T1;
    double US1_D1;
    double US1_MS;
    double US1_RE;
    double CQ1;
    double CQ1_U1;
    double CQ1_T1;
    double CQ1_D1;
    double CQ1_MS;
    double CQ1_RE;
    double DE1;
    double DE1_U1;
    double DE1_T1;
    double DE1_D1;
    double DE1_MS;

    double X2;
    double U2;
    double T2;
    double D2;
    double S2;
    double AMPL2;
    double U2_UEI;
    double U2_MS;
    double DW2;
    double H2;
    double H2_T2;
    double H2_D2;
    double M2;
    double M2_U2;
    double M2_MS;
    double R2;
    double R2_U2;
    double R2_MS;
    double V2;
    double V2_U2;
    double V2_MS;
    double V2_RE;
    double HK2;
    double HK2_U2;
    double HK2_T2;
    double HK2_D2;
    double HK2_MS;
    double HS2;
    double HS2_U2;
    double HS2_T2;
    double HS2_D2;
    double HS2_MS;
    double HS2_RE;
    double HC2;
    double HC2_U2;
    double HC2_T2;
    double HC2_D2;
    double HC2_MS;
    double RT2;
    double RT2_U2;
    double RT2_T2;
    double RT2_MS;
    double RT2_RE;
    double CF2;
    double CF2_HK2;
    double CF2_M2;
    double CF2_RT2;
    double CF2_U2;
    double CF2_T2;
    double CF2_D2;
    double CF2_MS;
    double CF2_RE;
    double DI2;
    double DI2_U2;
    double DI2_T2;
    double DI2_D2;
    double DI2_S2;
    double DI2_MS;
    double DI2_RE;
    double US2;
    double US2_U2;
    double US2_T2;
    double US2_D2;
    double US2_MS;
    double US2_RE;
    double CQ2;
    double CQ2_U2;
    double CQ2_T2;
    double CQ2_D2;
    double CQ2_MS;
    double CQ2_RE;
    double DE2;
    double DE2_U2;
    double DE2_T2;
    double DE2_D2;
    double DE2_MS;

    double CFM;
    double CFM_HKA;
    double CFM_MA;
    double CFM_MS;
    double CFM_RE;
    double CFM_RTA;
    double CFM_U1;
    double CFM_T1;
    double CFM_D1;
    double CFM_U2;
    double CFM_T2;
    double CFM_D2;
    double XT;
    double XT_A1;
    double XT_A2;
    double XT_MS;
    double XT_RE;
    double XT_XF;
    double XT_X1;
    double XT_T1;
    double XT_D1;
    double XT_U1;
    double XT_X2;
    double XT_T2;
    double XT_D2;
    double XT_U2;
    double DWTE;
    double QINFBL;
    double TKBL;
    double TKBL_MS;
    double RSTBL;
    double RSTBL_MS;
    double HSTINV;
    double HSTINV_MS;
    double REYBL;
    double REYBL_MS;
    double REYBL_RE;
    double GAMBL;
    double GM1BL;
    double HVRAT;
    double BULE;
    double XIFORC;
    double AMCRIT;

    std::vector<double *> com1_fields;
    std::vector<double *> com2_fields;
};

void setbl(XFoilState &ctx, XBlState &bl);
void iblsys(XFoilState &ctx);
void mrchue(XFoilState &ctx, XBlState &bl);
void mrchdu(XFoilState &ctx, XBlState &bl);
void xifset(XFoilState &ctx, XBlState &bl, int is_);
void update(XFoilState &ctx, XBlState &bl);

double dslim(double dstr, double thet, double uedg, double msq, double hklim);
void blpini(XBlState &bl);

#endif  // WASM_XBL_H
