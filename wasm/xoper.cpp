#include "xoper.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "spline.h"
#include "xbl.h"
#include "xblsys.h"
#include "xfoil.h"
#include "xgdes.h"
#include "xpanel.h"
#include "xsolve.h"
#include "xutils.h"

void mhinge(XFoilState &ctx) {
    double tops = 0.0;
    double bots = 0.0;
    if (!ctx.LFLAP) {
        std::tie(tops, bots, ctx.XOF, ctx.YOF) = getxyf(ctx.X, ctx.XP, ctx.Y, ctx.YP, ctx.S, ctx.N, 0.0, 0.0, ctx.XOF, ctx.YOF);
        ctx.LFLAP = true;
    } else {
        tops = ctx.XOF;
        bots = ctx.S[ctx.N] - ctx.XOF;
        tops = sinvrt(tops, ctx.XOF, ctx.X, ctx.XP, ctx.S, ctx.N);
        bots = sinvrt(bots, ctx.XOF, ctx.X, ctx.XP, ctx.S, ctx.N);
    }

    const double topx = seval(tops, ctx.X, ctx.XP, ctx.S, ctx.N);
    const double topy = seval(tops, ctx.Y, ctx.YP, ctx.S, ctx.N);
    const double botx = seval(bots, ctx.X, ctx.XP, ctx.S, ctx.N);
    const double boty = seval(bots, ctx.Y, ctx.YP, ctx.S, ctx.N);

    ctx.HMOM = 0.0;
    ctx.HFX = 0.0;
    ctx.HFY = 0.0;

    for (int i = 2; i <= ctx.N; ++i) {
        if (ctx.S[i - 1] >= tops && ctx.S[i] <= bots) {
            continue;
        }

        const double dx = ctx.X[i] - ctx.X[i - 1];
        const double dy = ctx.Y[i] - ctx.Y[i - 1];
        const double xmid = 0.5 * (ctx.X[i] + ctx.X[i - 1]) - ctx.XOF;
        const double ymid = 0.5 * (ctx.Y[i] + ctx.Y[i - 1]) - ctx.YOF;
        const double pmid = ctx.LVISC ? 0.5 * (ctx.CPV[i] + ctx.CPV[i - 1]) : 0.5 * (ctx.CPI[i] + ctx.CPI[i - 1]);
        ctx.HMOM = ctx.HMOM + pmid * (xmid * dx + ymid * dy);
        ctx.HFX = ctx.HFX - pmid * dy;
        ctx.HFY = ctx.HFY + pmid * dx;
    }

    int i = 2;
    for (i = 2; i <= ctx.N; ++i) {
        if (ctx.S[i] > tops) {
            break;
        }
    }

    double dx = topx - ctx.X[i - 1];
    double dy = topy - ctx.Y[i - 1];
    double xmid = 0.5 * (topx + ctx.X[i - 1]) - ctx.XOF;
    double ymid = 0.5 * (topy + ctx.Y[i - 1]) - ctx.YOF;
    double frac = 0.0;
    if (ctx.S[i] != ctx.S[i - 1]) {
        frac = (tops - ctx.S[i - 1]) / (ctx.S[i] - ctx.S[i - 1]);
    }
    double pmid = 0.0;
    if (ctx.LVISC) {
        const double topp = ctx.CPV[i] * frac + ctx.CPV[i - 1] * (1.0 - frac);
        pmid = 0.5 * (topp + ctx.CPV[i - 1]);
    } else {
        const double topp = ctx.CPI[i] * frac + ctx.CPI[i - 1] * (1.0 - frac);
        pmid = 0.5 * (topp + ctx.CPI[i - 1]);
    }
    ctx.HMOM = ctx.HMOM + pmid * (xmid * dx + ymid * dy);
    ctx.HFX = ctx.HFX - pmid * dy;
    ctx.HFY = ctx.HFY + pmid * dx;

    dx = ctx.XOF - topx;
    dy = ctx.YOF - topy;
    xmid = 0.5 * (topx + ctx.XOF) - ctx.XOF;
    ymid = 0.5 * (topy + ctx.YOF) - ctx.YOF;
    ctx.HMOM = ctx.HMOM + pmid * (xmid * dx + ymid * dy);
    ctx.HFX = ctx.HFX - pmid * dy;
    ctx.HFY = ctx.HFY + pmid * dx;

    for (i = ctx.N; i >= 2; --i) {
        if (ctx.S[i - 1] < bots) {
            break;
        }
    }

    dx = ctx.X[i] - botx;
    dy = ctx.Y[i] - boty;
    xmid = 0.5 * (ctx.X[i] + botx) - ctx.XOF;
    ymid = 0.5 * (ctx.Y[i] + boty) - ctx.YOF;
    if (ctx.S[i] != ctx.S[i - 1]) {
        frac = (ctx.S[i] - bots) / (ctx.S[i] - ctx.S[i - 1]);
    } else {
        frac = 0.0;
    }
    if (ctx.LVISC) {
        const double botp = ctx.CPV[i - 1] * frac + ctx.CPV[i] * (1.0 - frac);
        pmid = 0.5 * (botp + ctx.CPV[i]);
    } else {
        const double botp = ctx.CPI[i - 1] * frac + ctx.CPI[i] * (1.0 - frac);
        pmid = 0.5 * (botp + ctx.CPI[i]);
    }
    ctx.HMOM = ctx.HMOM + pmid * (xmid * dx + ymid * dy);
    ctx.HFX = ctx.HFX - pmid * dy;
    ctx.HFY = ctx.HFY + pmid * dx;

    dx = botx - ctx.XOF;
    dy = boty - ctx.YOF;
    xmid = 0.5 * (botx + ctx.XOF) - ctx.XOF;
    ymid = 0.5 * (boty + ctx.YOF) - ctx.YOF;
    ctx.HMOM = ctx.HMOM + pmid * (xmid * dx + ymid * dy);
    ctx.HFX = ctx.HFX - pmid * dy;
    ctx.HFY = ctx.HFY + pmid * dx;
}

void viscal(XFoilState &ctx, XBlState &bl, int niter1) {
    const double eps1 = 1.0e-4;
    const int niter = niter1;

    if (!ctx.LWAKE) {
        xywake(ctx);
    }

    qwcalc(ctx);
    qiset(ctx);

    if (!ctx.LIPAN) {
        if (ctx.LBLINI) {
            gamqv(ctx);
        }

        stfind(ctx);
        iblpan(ctx);
        xicalc(ctx);
        iblsys(ctx);
    }

    uicalc(ctx);

    if (!ctx.LBLINI) {
        for (int ibl = 1; ibl <= ctx.NBL[1]; ++ibl) {
            ctx.UEDG[ibl][1] = ctx.UINV[ibl][1];
        }
        for (int ibl = 1; ibl <= ctx.NBL[2]; ++ibl) {
            ctx.UEDG[ibl][2] = ctx.UINV[ibl][2];
        }
    }

    if (ctx.LVCONV) {
        qvfue(ctx);
        if (ctx.LVISC) {
            cpcalc(ctx.N + ctx.NW, ctx.QVIS, ctx.QINF, ctx.MINF, ctx.CPV);
            cpcalc(ctx.N + ctx.NW, ctx.QINV, ctx.QINF, ctx.MINF, ctx.CPI);
        } else {
            cpcalc(ctx.N, ctx.QINV, ctx.QINF, ctx.MINF, ctx.CPI);
        }
        gamqv(ctx);
        std::tie(ctx.CL, ctx.CM, ctx.CDP, ctx.CL_ALF, ctx.CL_MSQ) = clcalc(ctx.N, ctx.X, ctx.Y, ctx.GAM, ctx.GAM_A, ctx.ALFA,
                                                                           ctx.MINF, ctx.QINF, ctx.XCMREF, ctx.YCMREF);
        cdcalc(ctx);
    }

    if (!ctx.LWDIJ || !ctx.LADIJ) {
        qdcalc(ctx);
    }

    if (niter == 0) {
        throw std::runtime_error("VISCAL: NITER=0 not supported.");
    }

    std::cout << std::endl;
    std::cout << "Solving BL system ..." << std::endl;
    for (int iter_ = 1; iter_ <= niter; ++iter_) {
        setbl(ctx, bl);
        blsolv(ctx);
        update(ctx, bl);

        if (ctx.LALFA) {
            mrcl(ctx, ctx.CL);
            comset(ctx);
        } else {
            qiset(ctx);
            uicalc(ctx);
        }

        qvfue(ctx);
        gamqv(ctx);
        stmove(ctx);

        std::tie(ctx.CL, ctx.CM, ctx.CDP, ctx.CL_ALF, ctx.CL_MSQ) = clcalc(ctx.N, ctx.X, ctx.Y, ctx.GAM, ctx.GAM_A, ctx.ALFA,
                                                                           ctx.MINF, ctx.QINF, ctx.XCMREF, ctx.YCMREF);
        cdcalc(ctx);

        if (ctx.RLX < 1.0) {
            std::cout << "\n" << std::setw(3) << iter_ << "   rms: " << std::setw(10) << std::scientific << ctx.RMSBL
                      << "   max: " << std::setw(10) << ctx.RMXBL << std::defaultfloat << "   " << ctx.VMXBL << " at "
                      << std::setw(4) << ctx.IMXBL << std::setw(3) << ctx.ISMXBL << "   RLX:" << std::setw(6)
                      << std::fixed << std::setprecision(3) << ctx.RLX << std::defaultfloat << std::endl;
        }
        if (ctx.RLX == 1.0) {
            std::cout << "\n" << std::setw(3) << iter_ << "   rms: " << std::setw(10) << std::scientific << ctx.RMSBL
                      << "   max: " << std::setw(10) << ctx.RMXBL << std::defaultfloat << "   " << ctx.VMXBL << " at "
                      << std::setw(4) << ctx.IMXBL << std::setw(3) << ctx.ISMXBL << std::endl;
        }
        const double cdpdif = ctx.CD - ctx.CDF;
        std::cout << "    a =" << std::setw(7) << std::fixed << std::setprecision(3) << ctx.ALFA / ctx.DTOR
                  << "      CL =" << std::setw(8) << std::setprecision(4) << ctx.CL << std::defaultfloat << "\n"
                  << "   Cm =" << std::setw(8) << std::fixed << std::setprecision(4) << ctx.CM << std::defaultfloat
                  << "     CD =" << std::setw(9) << std::fixed << std::setprecision(5) << ctx.CD << std::defaultfloat
                  << "   =>   CDf =" << std::setw(9) << std::fixed << std::setprecision(5) << ctx.CDF << std::defaultfloat
                  << "    CDp =" << std::setw(9) << std::fixed << std::setprecision(5) << cdpdif << std::defaultfloat
                  << std::endl;

        if (ctx.RMSBL < eps1) {
            ctx.LVCONV = true;
            ctx.AVISC = ctx.ALFA;
            ctx.MVISC = ctx.MINF;
            break;
        }
    }

    if (!ctx.LVCONV) {
        std::cout << "VISCAL:  Convergence failed" << std::endl;
    }

    cpcalc(ctx.N + ctx.NW, ctx.QINV, ctx.QINF, ctx.MINF, ctx.CPI);
    cpcalc(ctx.N + ctx.NW, ctx.QVIS, ctx.QINF, ctx.MINF, ctx.CPV);
    if (ctx.LFLAP) {
        mhinge(ctx);
    }

    int is_ = 1;
    double hkmax = 0.0;
    double hkm = 0.0;
    double psep = 0.0;
    double patt = 0.0;
    for (int ibl = 2; ibl <= ctx.IBLTE[is_]; ++ibl) {
        const double hki = ctx.DSTR[ibl][is_] / ctx.THET[ibl][is_];
        hkmax = std::max(hki, hkmax);
        if (hkm < 4.0 && hki >= 4.0) {
            const double hfrac = (4.0 - hkm) / (hki - hkm);
            const double pdefm = ctx.UEDG[ibl - 1][is_] * ctx.UEDG[ibl - 1][is_] * ctx.THET[ibl - 1][is_];
            const double pdefi = ctx.UEDG[ibl][is_] * ctx.UEDG[ibl][is_] * ctx.THET[ibl][is_];
            psep = pdefm * (1.0 - hfrac) + pdefi * hfrac;
        }
        if (hkm > 4.0 && hki < 4.0) {
            const double hfrac = (4.0 - hkm) / (hki - hkm);
            const double pdefm = ctx.UEDG[ibl - 1][is_] * ctx.UEDG[ibl - 1][is_] * ctx.THET[ibl - 1][is_];
            const double pdefi = ctx.UEDG[ibl][is_] * ctx.UEDG[ibl][is_] * ctx.THET[ibl][is_];
            patt = pdefm * (1.0 - hfrac) + pdefi * hfrac;
        }
        hkm = hki;
    }
    const double delp = patt - psep;

    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << ctx.ACRIT[is_] << std::setw(10) << std::setprecision(4)
              << hkmax << std::setw(11) << std::setprecision(6) << ctx.CD << std::setw(11) << std::setprecision(6) << 2.0 * psep
              << std::setw(11) << std::setprecision(6) << 2.0 * patt << std::setw(11) << std::setprecision(6) << 2.0 * delp
              << std::setw(10) << std::setprecision(4) << ctx.XOCTR[is_] << std::defaultfloat << "     #" << std::endl;

    const double fnum = ctx.XSTRIP[is_] * 100.0;
    const int iten = static_cast<int>(fnum / 9.99999);
    const int ione = static_cast<int>((fnum - static_cast<double>(10 * iten)) / 0.99999);
    const int idec = static_cast<int>((fnum - static_cast<double>(10 * iten) - static_cast<double>(ione)) / 0.09999);
    const std::string fname = std::to_string(iten) + std::to_string(ione) + std::to_string(idec) + ".bl";

    std::ofstream lu(fname, std::ios::out);
    lu << "#       s         ue          H          P         K         x    -m du/dx\n";
    for (int ibl = 2; ibl <= ctx.IBLTE[is_]; ++ibl) {
        const int iblm = std::max(ibl - 1, 2);
        const int iblp = std::min(ibl + 1, ctx.IBLTE[is_]);
        const int i = ctx.IPAN[ibl][is_];
        const double hk = ctx.DSTR[ibl][is_] / ctx.THET[ibl][is_];
        const double ddef = ctx.DSTR[ibl][is_] * ctx.UEDG[ibl][is_];
        const double pdef = ctx.THET[ibl][is_] * ctx.UEDG[ibl][is_] * ctx.UEDG[ibl][is_];
        const double edef = ctx.TSTR[ibl][is_] * ctx.UEDG[ibl][is_] * ctx.UEDG[ibl][is_] * ctx.UEDG[ibl][is_] * 0.5;
        const double duds = (ctx.UEDG[iblp][is_] - ctx.UEDG[iblm][is_]) / (ctx.XSSI[iblp][is_] - ctx.XSSI[iblm][is_]);
        const double dpds = -ddef * duds;
        lu << " " << std::setw(11) << std::fixed << std::setprecision(4) << ctx.XSSI[ibl][is_] << std::setw(11)
           << std::setprecision(4) << ctx.UEDG[ibl][is_] << std::setw(11) << std::setprecision(4) << hk << std::setw(11)
           << std::setprecision(6) << pdef << std::setw(11) << std::setprecision(6) << edef << std::setw(11)
           << std::setprecision(3) << ctx.X[i] << std::setw(14) << std::scientific << std::setprecision(6) << dpds
           << std::defaultfloat << "\n";
    }
}
