#include "xfoil.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "naca.h"
#include "spline.h"
#include "userio.h"
#include "xbl.h"
#include "xgeom.h"
#include "xpanel.h"

std::pair<double, double> mrcl(XFoilState &ctx, double cls) {
    const double cla = std::max(cls, 0.000001);

    if (ctx.RETYP < 1 || ctx.RETYP > 3) {
        std::cout << "MRCL:  Illegal Re(CL) dependence trigger." << std::endl;
        std::cout << "       Setting fixed Re." << std::endl;
        ctx.RETYP = 1;
    }
    if (ctx.MATYP < 1 || ctx.MATYP > 3) {
        std::cout << "MRCL:  Illegal Mach(CL) dependence trigger." << std::endl;
        std::cout << "       Setting fixed Mach." << std::endl;
        ctx.MATYP = 1;
    }

    double m_cls = 0.0;
    if (ctx.MATYP == 1) {
        ctx.MINF = ctx.MINF1;
        m_cls = 0.0;
    } else if (ctx.MATYP == 2) {
        ctx.MINF = ctx.MINF1 / std::sqrt(cla);
        m_cls = -0.5 * ctx.MINF / cla;
    } else {
        ctx.MINF = ctx.MINF1;
        m_cls = 0.0;
    }

    double r_cls = 0.0;
    if (ctx.RETYP == 1) {
        ctx.REINF = ctx.REINF1;
        r_cls = 0.0;
    } else if (ctx.RETYP == 2) {
        ctx.REINF = ctx.REINF1 / std::sqrt(cla);
        r_cls = -0.5 * ctx.REINF / cla;
    } else {
        ctx.REINF = ctx.REINF1 / cla;
        r_cls = -ctx.REINF / cla;
    }

    if (ctx.MINF >= 0.99) {
        std::cout << std::endl;
        std::cout << "MRCL: CL too low for chosen Mach(CL) dependence" << std::endl;
        std::cout << "      Aritificially limiting Mach to  0.99" << std::endl;
        ctx.MINF = 0.99;
        m_cls = 0.0;
    }

    double rrat = 1.0;
    if (ctx.REINF1 > 0.0) {
        rrat = ctx.REINF / ctx.REINF1;
    }

    if (rrat > 100.0) {
        std::cout << std::endl;
        std::cout << "MRCL: CL too low for chosen Re(CL) dependence" << std::endl;
        std::cout << "      Aritificially limiting Re to  " << ctx.REINF1 * 100.0 << std::endl;
        ctx.REINF = ctx.REINF1 * 100.0;
        r_cls = 0.0;
    }

    ctx.MINF_CL = m_cls;
    ctx.REINF_CL = r_cls;
    return {m_cls, r_cls};
}

void comset(XFoilState &ctx) {
    const double beta = std::sqrt(1.0 - ctx.MINF * ctx.MINF);
    const double beta_msq = -0.5 / beta;

    ctx.TKLAM = ctx.MINF * ctx.MINF / std::pow(1.0 + beta, 2.0);
    ctx.TKL_MSQ = 1.0 / std::pow(1.0 + beta, 2.0) - 2.0 * ctx.TKLAM / (1.0 + beta) * beta_msq;

    if (ctx.MINF == 0.0) {
        ctx.CPSTAR = -999.0;
        ctx.QSTAR = 999.0;
    } else {
        ctx.CPSTAR = 2.0 / (ctx.GAMMA * ctx.MINF * ctx.MINF)
                     * (std::pow((1.0 + 0.5 * ctx.GAMM1 * ctx.MINF * ctx.MINF) / (1.0 + 0.5 * ctx.GAMM1),
                                 ctx.GAMMA / ctx.GAMM1)
                        - 1.0);
        ctx.QSTAR = ctx.QINF / ctx.MINF
                    * std::sqrt((1.0 + 0.5 * ctx.GAMM1 * ctx.MINF * ctx.MINF) / (1.0 + 0.5 * ctx.GAMM1));
    }
}

void cpcalc(int n, const std::vector<double> &q, double qinf, double minf, std::vector<double> &cp) {
    const double beta = std::sqrt(1.0 - minf * minf);
    const double bfac = 0.5 * minf * minf / (1.0 + beta);

    bool denneg = false;

    for (int i = 1; i <= n; ++i) {
        const double cpinc = 1.0 - std::pow(q[i] / qinf, 2.0);
        const double den = beta + bfac * cpinc;
        cp[i] = cpinc / den;
        if (den <= 0.0) {
            denneg = true;
        }
    }

    if (denneg) {
        std::cout << std::endl;
        std::cout << "CPCALC: Local speed too large. Compressibility corrections invalid." << std::endl;
    }
}

std::tuple<double, double, double, double, double> clcalc(int n, const std::vector<double> &x, const std::vector<double> &y,
                                                          const std::vector<double> &gam, const std::vector<double> &gam_a,
                                                          double alfa, double minf, double qinf, double xref, double yref) {
    const double sa = std::sin(alfa);
    const double ca = std::cos(alfa);

    const double beta = std::sqrt(1.0 - minf * minf);
    const double beta_msq = -0.5 / beta;

    const double bfac = 0.5 * minf * minf / (1.0 + beta);
    const double bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

    double cl = 0.0;
    double cm = 0.0;
    double cdp = 0.0;
    double cl_alf = 0.0;
    double cl_msq = 0.0;

    int i = 1;
    double cginc = 1.0 - std::pow(gam[i] / qinf, 2.0);
    double cpg1 = cginc / (beta + bfac * cginc);
    double cpg1_msq = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

    double cpi_gam = -2.0 * gam[i] / (qinf * qinf);
    double cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
    double cpg1_alf = cpc_cpi * cpi_gam * gam_a[i];

    for (i = 1; i <= n; ++i) {
        int ip = i + 1;
        if (i == n) {
            ip = 1;
        }

        cginc = 1.0 - std::pow(gam[ip] / qinf, 2.0);
        const double cpg2 = cginc / (beta + bfac * cginc);
        const double cpg2_msq = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);

        cpi_gam = -2.0 * gam[ip] / (qinf * qinf);
        cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
        const double cpg2_alf = cpc_cpi * cpi_gam * gam_a[ip];

        const double dx = (x[ip] - x[i]) * ca + (y[ip] - y[i]) * sa;
        const double dy = (y[ip] - y[i]) * ca - (x[ip] - x[i]) * sa;
        const double dg = cpg2 - cpg1;

        const double ax = (0.5 * (x[ip] + x[i]) - xref) * ca + (0.5 * (y[ip] + y[i]) - yref) * sa;
        const double ay = (0.5 * (y[ip] + y[i]) - yref) * ca - (0.5 * (x[ip] + x[i]) - xref) * sa;
        const double ag = 0.5 * (cpg2 + cpg1);

        const double dx_alf = -(x[ip] - x[i]) * sa + (y[ip] - y[i]) * ca;
        const double ag_alf = 0.5 * (cpg2_alf + cpg1_alf);
        const double ag_msq = 0.5 * (cpg2_msq + cpg1_msq);

        cl = cl + dx * ag;
        cdp = cdp - dy * ag;
        cm = cm - dx * (ag * ax + dg * dx / 12.0) - dy * (ag * ay + dg * dy / 12.0);

        cl_alf = cl_alf + dx * ag_alf + ag * dx_alf;
        cl_msq = cl_msq + dx * ag_msq;

        cpg1 = cpg2;
        cpg1_alf = cpg2_alf;
        cpg1_msq = cpg2_msq;
    }

    return {cl, cm, cdp, cl_alf, cl_msq};
}

void cdcalc(XFoilState &ctx) {
    const double sa = std::sin(ctx.ALFA);
    const double ca = std::cos(ctx.ALFA);

    if (ctx.LVISC && ctx.LBLINI) {
        const double thwake = ctx.THET[ctx.NBL[2]][2];
        const double urat = ctx.UEDG[ctx.NBL[2]][2] / ctx.QINF;
        const double uewake = ctx.UEDG[ctx.NBL[2]][2] * (1.0 - ctx.TKLAM) / (1.0 - ctx.TKLAM * urat * urat);
        const double shwake = ctx.DSTR[ctx.NBL[2]][2] / ctx.THET[ctx.NBL[2]][2];
        ctx.CD = 2.0 * thwake * std::pow(uewake / ctx.QINF, 0.5 * (5.0 + shwake));
    } else {
        ctx.CD = 0.0;
    }

    ctx.CDF = 0.0;
    for (int is_ = 1; is_ <= 2; ++is_) {
        for (int ibl = 3; ibl <= ctx.IBLTE[is_]; ++ibl) {
            const int i = ctx.IPAN[ibl][is_];
            const int im = ctx.IPAN[ibl - 1][is_];
            const double dx = (ctx.X[i] - ctx.X[im]) * ca + (ctx.Y[i] - ctx.Y[im]) * sa;
            ctx.CDF = ctx.CDF + 0.5 * (ctx.TAU[ibl][is_] + ctx.TAU[ibl - 1][is_]) * dx * 2.0 / (ctx.QINF * ctx.QINF);
        }
    }
}

void tecalc(XFoilState &ctx) {
    const double dxte = ctx.X[1] - ctx.X[ctx.N];
    const double dyte = ctx.Y[1] - ctx.Y[ctx.N];
    const double dxs = 0.5 * (-ctx.XP[1] + ctx.XP[ctx.N]);
    const double dys = 0.5 * (-ctx.YP[1] + ctx.YP[ctx.N]);

    ctx.ANTE = dxs * dyte - dys * dxte;
    ctx.ASTE = dxs * dxte + dys * dyte;

    ctx.DSTE = std::sqrt(dxte * dxte + dyte * dyte);

    ctx.SHARP = ctx.DSTE < 0.0001 * ctx.CHORD;

    double scs = 1.0;
    double sds = 0.0;
    if (!ctx.SHARP) {
        scs = ctx.ANTE / ctx.DSTE;
        sds = ctx.ASTE / ctx.DSTE;
    }

    ctx.SIGTE = 0.5 * (ctx.GAM[1] - ctx.GAM[ctx.N]) * scs;
    ctx.GAMTE = -0.5 * (ctx.GAM[1] - ctx.GAM[ctx.N]) * sds;

    ctx.SIGTE_A = 0.5 * (ctx.GAM_A[1] - ctx.GAM_A[ctx.N]) * scs;
    ctx.GAMTE_A = -0.5 * (ctx.GAM_A[1] - ctx.GAM_A[ctx.N]) * sds;
}

void naca(XFoilState &ctx, int ides1) {
    const int iqx = (static_cast<int>(ctx.W1.size()) - 1) / 6;
    const int nside = iqx / 3;

    if (ides1 <= 0) {
        throw std::runtime_error("NACA: IDES must be specified.");
    }
    int ides = ides1;

    int itype = 0;
    if (ides <= 25099) {
        itype = 5;
    }
    if (ides <= 9999) {
        itype = 4;
    }

    if (itype == 0) {
        std::cout << "This designation not implemented." << std::endl;
        return;
    }

    std::vector<double> &xx = ctx.W1;
    std::vector<double> &yt = ctx.W2;
    std::vector<double> &yc = ctx.W3;

    int nb = 0;
    std::string name;
    if (itype == 4) {
        auto out = naca4(ides, xx, yt, yc, nside, ctx.XB, ctx.YB);
        nb = out.first;
        name = out.second;
    } else {
        auto out = naca5(ides, xx, yt, yc, nside, ctx.XB, ctx.YB);
        nb = out.first;
        name = out.second;
    }

    ctx.NB = nb;
    ctx.NAME = name;
    auto stripped = strip_string(ctx.NAME);
    ctx.NAME = stripped.first;
    ctx.NNAME = stripped.second;

    if (ctx.NB == 0) {
        return;
    }

    ctx.LCLOCK = false;

    ctx.XBF = 0.0;
    ctx.YBF = 0.0;
    ctx.LBFLAP = false;

    scalc(ctx.XB, ctx.YB, ctx.SB, ctx.NB);
    segspl(ctx.XB, ctx.XBP, ctx.SB, ctx.NB);
    segspl(ctx.YB, ctx.YBP, ctx.SB, ctx.NB);

    const GeoparResults gp = geopar(ctx.XB, ctx.XBP, ctx.YB, ctx.YBP, ctx.SB, ctx.NB, ctx.W1);
    ctx.SBLE = gp.sle;
    ctx.CHORDB = gp.chord;
    ctx.AREAB = gp.area;
    ctx.RADBLE = gp.radle;
    ctx.ANGBTE = gp.angte;
    ctx.EI11BA = gp.ei11a;
    ctx.EI22BA = gp.ei22a;
    ctx.APX1BA = gp.apx1a;
    ctx.APX2BA = gp.apx2a;
    ctx.EI11BT = gp.ei11t;
    ctx.EI22BT = gp.ei22t;
    ctx.APX1BT = gp.apx1t;
    ctx.APX2BT = gp.apx2t;
    ctx.THICKB = gp.thick;
    ctx.CAMBRB = gp.cambr;

    std::cout << "\n Buffer airfoil set using" << std::setw(4) << ctx.NB << " points" << std::endl;

    pangen(ctx, true);
}

void pangen(XFoilState &ctx, bool shopar) {
    if (ctx.NB < 2) {
        std::cout << "PANGEN: Buffer airfoil not available." << std::endl;
        ctx.N = 0;
        return;
    }

    const int ipfac = 5;

    ctx.N = ctx.NPAN;

    scalc(ctx.XB, ctx.YB, ctx.SB, ctx.NB);
    segspl(ctx.XB, ctx.XBP, ctx.SB, ctx.NB);
    segspl(ctx.YB, ctx.YBP, ctx.SB, ctx.NB);

    const double sbref = 0.5 * (ctx.SB[ctx.NB] - ctx.SB[1]);

    for (int i = 1; i <= ctx.NB; ++i) {
        ctx.W5[i] = std::abs(curv(ctx.SB[i], ctx.XB, ctx.XBP, ctx.YB, ctx.YBP, ctx.SB, ctx.NB)) * sbref;
    }

    ctx.SBLE = lefind(ctx.XB, ctx.XBP, ctx.YB, ctx.YBP, ctx.SB, ctx.NB);
    const double cvle = std::abs(curv(ctx.SBLE, ctx.XB, ctx.XBP, ctx.YB, ctx.YBP, ctx.SB, ctx.NB)) * sbref;

    int ible = 0;
    for (int i = 1; i <= ctx.NB - 1; ++i) {
        if (ctx.SBLE == ctx.SB[i] && ctx.SBLE == ctx.SB[i + 1]) {
            ible = i;
            std::cout << std::endl;
            std::cout << "Sharp leading edge" << std::endl;
            break;
        }
    }

    const double xble = seval(ctx.SBLE, ctx.XB, ctx.XBP, ctx.SB, ctx.NB);
    const double yble = seval(ctx.SBLE, ctx.YB, ctx.YBP, ctx.SB, ctx.NB);
    const double xbte = 0.5 * (ctx.XB[1] + ctx.XB[ctx.NB]);
    const double ybte = 0.5 * (ctx.YB[1] + ctx.YB[ctx.NB]);
    const double chbsq = (xbte - xble) * (xbte - xble) + (ybte - yble) * (ybte - yble);

    const int nk = 3;
    double cvsum = 0.0;
    for (int k = -nk; k <= nk; ++k) {
        const double frac = static_cast<double>(k) / static_cast<double>(nk);
        const double sbk = ctx.SBLE + frac * sbref / std::max(cvle, 20.0);
        const double cvk = std::abs(curv(sbk, ctx.XB, ctx.XBP, ctx.YB, ctx.YBP, ctx.SB, ctx.NB)) * sbref;
        cvsum = cvsum + cvk;
    }
    double cvavg = cvsum / static_cast<double>(2 * nk + 1);

    if (ible != 0) {
        cvavg = 10.0;
    }

    const double cc = 6.0 * ctx.CVPAR;

    const double cvte = cvavg * ctx.CTERAT;
    ctx.W5[1] = cvte;
    ctx.W5[ctx.NB] = cvte;

    const double smool = std::max(1.0 / std::max(cvavg, 20.0), 0.25 / static_cast<double>(ctx.NPAN / 2));
    const double smoosq = (smool * sbref) * (smool * sbref);

    ctx.W2[1] = 1.0;
    ctx.W3[1] = 0.0;
    for (int i = 2; i <= ctx.NB - 1; ++i) {
        const double dsm = ctx.SB[i] - ctx.SB[i - 1];
        const double dsp = ctx.SB[i + 1] - ctx.SB[i];
        const double dso = 0.5 * (ctx.SB[i + 1] - ctx.SB[i - 1]);

        if (dsm == 0.0 || dsp == 0.0) {
            ctx.W1[i] = 0.0;
            ctx.W2[i] = 1.0;
            ctx.W3[i] = 0.0;
        } else {
            ctx.W1[i] = smoosq * (-1.0 / dsm) / dso;
            ctx.W2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            ctx.W3[i] = smoosq * (-1.0 / dsp) / dso;
        }
    }
    ctx.W1[ctx.NB] = 0.0;
    ctx.W2[ctx.NB] = 1.0;

    for (int i = 2; i <= ctx.NB - 1; ++i) {
        if (ctx.SB[i] == ctx.SBLE || i == ible || i == ible + 1) {
            ctx.W1[i] = 0.0;
            ctx.W2[i] = 1.0;
            ctx.W3[i] = 0.0;
            ctx.W5[i] = cvle;
        } else if (ctx.SB[i - 1] < ctx.SBLE && ctx.SB[i] > ctx.SBLE) {
            double dsm = ctx.SB[i - 1] - ctx.SB[i - 2];
            double dsp = ctx.SBLE - ctx.SB[i - 1];
            double dso = 0.5 * (ctx.SBLE - ctx.SB[i - 2]);

            ctx.W1[i - 1] = smoosq * (-1.0 / dsm) / dso;
            ctx.W2[i - 1] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            ctx.W3[i - 1] = 0.0;
            ctx.W5[i - 1] = ctx.W5[i - 1] + smoosq * cvle / (dsp * dso);

            dsm = ctx.SB[i] - ctx.SBLE;
            dsp = ctx.SB[i + 1] - ctx.SB[i];
            dso = 0.5 * (ctx.SB[i + 1] - ctx.SBLE);
            ctx.W1[i] = 0.0;
            ctx.W2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            ctx.W3[i] = smoosq * (-1.0 / dsp) / dso;
            ctx.W5[i] = ctx.W5[i] + smoosq * cvle / (dsm * dso);
            break;
        }
    }

    for (int i = 2; i <= ctx.NB - 1; ++i) {
        const double xoc = ((ctx.XB[i] - xble) * (xbte - xble) + (ctx.YB[i] - yble) * (ybte - yble)) / chbsq;

        if (ctx.SB[i] < ctx.SBLE) {
            if (xoc > ctx.XSREF1 && xoc < ctx.XSREF2) {
                ctx.W1[i] = 0.0;
                ctx.W2[i] = 1.0;
                ctx.W3[i] = 0.0;
                ctx.W5[i] = cvle * ctx.CTRRAT;
            }
        } else {
            if (xoc > ctx.XPREF1 && xoc < ctx.XPREF2) {
                ctx.W1[i] = 0.0;
                ctx.W2[i] = 1.0;
                ctx.W3[i] = 0.0;
                ctx.W5[i] = cvle * ctx.CTRRAT;
            }
        }
    }

    if (ible == 0) {
        trisol(ctx.W2, ctx.W1, ctx.W3, ctx.W5, ctx.NB);
    } else {
        const int nn1 = ible;
        std::vector<double> a(static_cast<size_t>(nn1 + 1), 0.0);
        std::vector<double> b(static_cast<size_t>(nn1 + 1), 0.0);
        std::vector<double> c(static_cast<size_t>(nn1 + 1), 0.0);
        std::vector<double> d(static_cast<size_t>(nn1 + 1), 0.0);
        for (int i = 1; i <= nn1; ++i) {
            a[i] = ctx.W2[i];
            b[i] = ctx.W1[i];
            c[i] = ctx.W3[i];
            d[i] = ctx.W5[i];
        }
        trisol(a, b, c, d, nn1);
        for (int i = 1; i <= nn1; ++i) {
            ctx.W5[i] = d[i];
        }

        const int nn2 = ctx.NB - ible;
        a.assign(static_cast<size_t>(nn2 + 1), 0.0);
        b.assign(static_cast<size_t>(nn2 + 1), 0.0);
        c.assign(static_cast<size_t>(nn2 + 1), 0.0);
        d.assign(static_cast<size_t>(nn2 + 1), 0.0);
        for (int i = 1; i <= nn2; ++i) {
            const int idx = ible + i;
            a[i] = ctx.W2[idx];
            b[i] = ctx.W1[idx];
            c[i] = ctx.W3[idx];
            d[i] = ctx.W5[idx];
        }
        trisol(a, b, c, d, nn2);
        for (int i = 1; i <= nn2; ++i) {
            ctx.W5[ible + i] = d[i];
        }
    }

    double cvmax = 0.0;
    for (int i = 1; i <= ctx.NB; ++i) {
        cvmax = std::max(cvmax, std::abs(ctx.W5[i]));
    }

    for (int i = 1; i <= ctx.NB; ++i) {
        ctx.W5[i] = ctx.W5[i] / cvmax;
    }

    segspl(ctx.W5, ctx.W6, ctx.SB, ctx.NB);

    const int nn = ipfac * (ctx.N - 1) + 1;
    int nn1 = 0;

    const double rdste = 0.667;
    const double rtf = (rdste - 1.0) * 2.0 + 1.0;

    if (ible == 0) {
        const double dsavg = (ctx.SB[ctx.NB] - ctx.SB[1]) / (static_cast<double>(nn - 3) + 2.0 * rtf);
        ctx.SNEW[1] = ctx.SB[1];
        for (int i = 2; i <= nn - 1; ++i) {
            ctx.SNEW[i] = ctx.SB[1] + dsavg * (static_cast<double>(i - 2) + rtf);
        }
        ctx.SNEW[nn] = ctx.SB[ctx.NB];
    } else {
        const int nfrac1 = (ctx.N * ible) / ctx.NB;
        nn1 = ipfac * (nfrac1 - 1) + 1;
        const double dsavg1 = (ctx.SBLE - ctx.SB[1]) / (static_cast<double>(nn1 - 2) + rtf);
        ctx.SNEW[1] = ctx.SB[1];
        for (int i = 2; i <= nn1; ++i) {
            ctx.SNEW[i] = ctx.SB[1] + dsavg1 * (static_cast<double>(i - 2) + rtf);
        }

        const int nn2 = nn - nn1 + 1;
        const double dsavg2 = (ctx.SB[ctx.NB] - ctx.SBLE) / (static_cast<double>(nn2 - 2) + rtf);
        for (int i = 2; i <= nn2 - 1; ++i) {
            ctx.SNEW[i - 1 + nn1] = ctx.SBLE + dsavg2 * (static_cast<double>(i - 2) + rtf);
        }
        ctx.SNEW[nn] = ctx.SB[ctx.NB];
    }

    for (int iter = 1; iter <= 20; ++iter) {
        double cv1 = seval(ctx.SNEW[1], ctx.W5, ctx.W6, ctx.SB, ctx.NB);
        double cv2 = seval(ctx.SNEW[2], ctx.W5, ctx.W6, ctx.SB, ctx.NB);
        double cvs1 = deval(ctx.SNEW[1], ctx.W5, ctx.W6, ctx.SB, ctx.NB);
        double cvs2 = deval(ctx.SNEW[2], ctx.W5, ctx.W6, ctx.SB, ctx.NB);

        double cavm = std::sqrt(cv1 * cv1 + cv2 * cv2);
        double cavm_s1 = 0.0;
        double cavm_s2 = 0.0;
        if (cavm != 0.0) {
            cavm_s1 = cvs1 * cv1 / cavm;
            cavm_s2 = cvs2 * cv2 / cavm;
        }

        for (int i = 2; i <= nn - 1; ++i) {
            const double dsm = ctx.SNEW[i] - ctx.SNEW[i - 1];
            const double dsp = ctx.SNEW[i] - ctx.SNEW[i + 1];
            const double cv3 = seval(ctx.SNEW[i + 1], ctx.W5, ctx.W6, ctx.SB, ctx.NB);
            const double cvs3 = deval(ctx.SNEW[i + 1], ctx.W5, ctx.W6, ctx.SB, ctx.NB);

            double cavp = std::sqrt(cv3 * cv3 + cv2 * cv2);
            double cavp_s2 = 0.0;
            double cavp_s3 = 0.0;
            if (cavp != 0.0) {
                cavp_s2 = cvs2 * cv2 / cavp;
                cavp_s3 = cvs3 * cv3 / cavp;
            }

            const double fm = cc * cavm + 1.0;
            const double fp = cc * cavp + 1.0;

            const double rez = dsp * fp + dsm * fm;

            ctx.W1[i] = -fm + cc * dsm * cavm_s1;
            ctx.W2[i] = fp + fm + cc * (dsp * cavp_s2 + dsm * cavm_s2);
            ctx.W3[i] = -fp + cc * dsp * cavp_s3;

            ctx.W4[i] = -rez;

            cv1 = cv2;
            cv2 = cv3;
            cvs1 = cvs2;
            cvs2 = cvs3;
            cavm = cavp;
            cavm_s1 = cavp_s2;
            cavm_s2 = cavp_s3;
        }

        ctx.W2[1] = 1.0;
        ctx.W3[1] = 0.0;
        ctx.W4[1] = 0.0;
        ctx.W1[nn] = 0.0;
        ctx.W2[nn] = 1.0;
        ctx.W4[nn] = 0.0;

        if (rtf != 1.0) {
            int i = 2;
            ctx.W4[i] = -((ctx.SNEW[i] - ctx.SNEW[i - 1]) + rtf * (ctx.SNEW[i] - ctx.SNEW[i + 1]));
            ctx.W1[i] = -1.0;
            ctx.W2[i] = 1.0 + rtf;
            ctx.W3[i] = -rtf;

            i = nn - 1;
            ctx.W4[i] = -((ctx.SNEW[i] - ctx.SNEW[i + 1]) + rtf * (ctx.SNEW[i] - ctx.SNEW[i - 1]));
            ctx.W3[i] = -1.0;
            ctx.W2[i] = 1.0 + rtf;
            ctx.W1[i] = -rtf;
        }

        if (ible != 0) {
            const int i = nn1;
            ctx.W1[i] = 0.0;
            ctx.W2[i] = 1.0;
            ctx.W3[i] = 0.0;
            ctx.W4[i] = ctx.SBLE - ctx.SNEW[i];
        }

        trisol(ctx.W2, ctx.W1, ctx.W3, ctx.W4, nn);

        double rlx = 1.0;
        double dmax = 0.0;
        for (int i = 1; i <= nn - 1; ++i) {
            const double ds = ctx.SNEW[i + 1] - ctx.SNEW[i];
            const double dds = ctx.W4[i + 1] - ctx.W4[i];
            const double dsrat = 1.0 + rlx * dds / ds;
            if (dsrat > 4.0) {
                rlx = (4.0 - 1.0) * ds / dds;
            }
            if (dsrat < 0.2) {
                rlx = (0.2 - 1.0) * ds / dds;
            }
            dmax = std::max(std::abs(ctx.W4[i]), dmax);
        }

        for (int i = 2; i <= nn - 1; ++i) {
            ctx.SNEW[i] = ctx.SNEW[i] + rlx * ctx.W4[i];
        }

        if (std::abs(dmax) < 1.0e-3) {
            break;
        }
    }

    for (int i = 1; i <= ctx.N; ++i) {
        const int ind = ipfac * (i - 1) + 1;
        ctx.S[i] = ctx.SNEW[ind];
        ctx.X[i] = seval(ctx.SNEW[ind], ctx.XB, ctx.XBP, ctx.SB, ctx.NB);
        ctx.Y[i] = seval(ctx.SNEW[ind], ctx.YB, ctx.YBP, ctx.SB, ctx.NB);
    }

    for (int ib = 1; ib <= ctx.NB - 1; ++ib) {
        if (ctx.SB[ib] == ctx.SB[ib + 1]) {
            const double xbcorn = ctx.XB[ib];
            const double ybcorn = ctx.YB[ib];
            const double sbcorn = ctx.SB[ib];

            for (int i = 1; i <= ctx.N; ++i) {
                if (ctx.S[i] <= sbcorn) {
                    continue;
                }

                for (int j = ctx.N; j >= i; --j) {
                    ctx.X[j + 1] = ctx.X[j];
                    ctx.Y[j + 1] = ctx.Y[j];
                    ctx.S[j + 1] = ctx.S[j];
                }
                ctx.N = ctx.N + 1;

                if (ctx.N > static_cast<int>(ctx.X.size()) - 2) {
                    throw std::runtime_error("PANEL: Too many panels. Increase IQX in XFOIL.INC");
                }

                ctx.X[i] = xbcorn;
                ctx.Y[i] = ybcorn;
                ctx.S[i] = sbcorn;

                if (i - 2 >= 1) {
                    ctx.S[i - 1] = 0.5 * (ctx.S[i] + ctx.S[i - 2]);
                    ctx.X[i - 1] = seval(ctx.S[i - 1], ctx.XB, ctx.XBP, ctx.SB, ctx.NB);
                    ctx.Y[i - 1] = seval(ctx.S[i - 1], ctx.YB, ctx.YBP, ctx.SB, ctx.NB);
                }

                if (i + 2 <= ctx.N) {
                    ctx.S[i + 1] = 0.5 * (ctx.S[i] + ctx.S[i + 2]);
                    ctx.X[i + 1] = seval(ctx.S[i + 1], ctx.XB, ctx.XBP, ctx.SB, ctx.NB);
                    ctx.Y[i + 1] = seval(ctx.S[i + 1], ctx.YB, ctx.YBP, ctx.SB, ctx.NB);
                }

                break;
            }
        }
    }

    scalc(ctx.X, ctx.Y, ctx.S, ctx.N);
    segspl(ctx.X, ctx.XP, ctx.S, ctx.N);
    segspl(ctx.Y, ctx.YP, ctx.S, ctx.N);
    ctx.SLE = lefind(ctx.X, ctx.XP, ctx.Y, ctx.YP, ctx.S, ctx.N);

    ctx.XLE = seval(ctx.SLE, ctx.X, ctx.XP, ctx.S, ctx.N);
    ctx.YLE = seval(ctx.SLE, ctx.Y, ctx.YP, ctx.S, ctx.N);
    ctx.XTE = 0.5 * (ctx.X[1] + ctx.X[ctx.N]);
    ctx.YTE = 0.5 * (ctx.Y[1] + ctx.Y[ctx.N]);
    ctx.CHORD = std::sqrt((ctx.XTE - ctx.XLE) * (ctx.XTE - ctx.XLE) + (ctx.YTE - ctx.YLE) * (ctx.YTE - ctx.YLE));

    double dsmin = 1000.0;
    double dsmax = -1000.0;
    for (int i = 1; i <= ctx.N - 1; ++i) {
        const double ds = ctx.S[i + 1] - ctx.S[i];
        if (ds == 0.0) {
            continue;
        }
        dsmin = std::min(dsmin, ds);
        dsmax = std::max(dsmax, ds);
    }

    dsmin = dsmin * static_cast<double>(ctx.N - 1) / ctx.S[ctx.N];
    dsmax = dsmax * static_cast<double>(ctx.N - 1) / ctx.S[ctx.N];

    ctx.LGAMU = false;
    ctx.LQINU = false;
    ctx.LWAKE = false;
    ctx.LQAIJ = false;
    ctx.LADIJ = false;
    ctx.LWDIJ = false;
    ctx.LIPAN = false;
    ctx.LBLINI = false;
    ctx.LVCONV = false;
    ctx.LSCINI = false;
    ctx.LQSPEC = false;
    ctx.LGSAME = false;

    if (ctx.LBFLAP) {
        ctx.XOF = ctx.XBF;
        ctx.YOF = ctx.YBF;
        ctx.LFLAP = true;
    }

    tecalc(ctx);

    ncalc(ctx.X, ctx.Y, ctx.S, ctx.N, ctx.NX, ctx.NY);
    apcalc(ctx);

    if (ctx.SHARP) {
        std::cout << "\nSharp trailing edge" << std::endl;
    } else {
        const double gap = std::sqrt((ctx.X[1] - ctx.X[ctx.N]) * (ctx.X[1] - ctx.X[ctx.N])
                                     + (ctx.Y[1] - ctx.Y[ctx.N]) * (ctx.Y[1] - ctx.Y[ctx.N]));
        std::cout << "\nBlunt trailing edge.  Gap =" << std::setw(9) << std::fixed << std::setprecision(5) << gap
                  << std::endl;
    }

    if (shopar) {
        std::cout << "\n Paneling parameters used..."
                  << "\n   Number of panel nodes      " << std::setw(4) << ctx.NPAN
                  << "\n   Panel bunching parameter   " << std::setw(6) << std::fixed << std::setprecision(3) << ctx.CVPAR
                  << "\n   TE/LE panel density ratio  " << std::setw(6) << std::fixed << std::setprecision(3) << ctx.CTERAT
                  << "\n   Refined-area/LE panel density ratio   " << std::setw(6) << std::fixed << std::setprecision(3) << ctx.CTRRAT
                  << "\n   Top    side refined area x/c limits  " << std::setw(6) << std::fixed << std::setprecision(3) << ctx.XSREF1
                  << std::setw(6) << std::fixed << std::setprecision(3) << ctx.XSREF2
                  << "\n   Bottom side refined area x/c limits  " << std::setw(6) << std::fixed << std::setprecision(3) << ctx.XPREF1
                  << std::setw(6) << std::fixed << std::setprecision(3) << ctx.XPREF2 << std::endl;
    }
}
