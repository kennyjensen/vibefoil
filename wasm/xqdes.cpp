// Ported from XFOIL Fortran source (Mark Drela).
// This file is a derived work and remains under the terms of the
// GNU General Public License v2 or later.
// See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

#include "xqdes.h"

#include <cmath>
#include <stdexcept>

#include "spline.h"
#include "xbl.h"
#include "xfoil.h"
#include "xpanel.h"
#include "xsolve.h"
#include "xutils.h"

void QDES(XFoilState & /*ctx*/) {
    throw std::runtime_error("QDES interactive command loop is not ported.");
}

void NEWPLOTQ(XFoilState & /*ctx*/) {
}

void QPLINI(XFoilState & /*ctx*/, bool /*ldef*/) {
}

void QSPLOT(XFoilState & /*ctx*/) {
}

void QSPPLT(XFoilState & /*ctx*/, int /*iqspl1*/, int /*iqspl2*/, int /*kqsp*/, int /*nt*/) {
}

void IQSGET(XFoilState & /*ctx*/) {
    throw std::runtime_error("IQSGET cursor selection is not ported.");
}

static void splind_offset(const std::vector<double> &x, std::vector<double> &xs, const std::vector<double> &s,
                          int start, int n, double xs1, double xs2) {
    if (n <= 0) {
        return;
    }
    std::vector<double> x_local(static_cast<size_t>(n + 1), 0.0);
    std::vector<double> xs_local(static_cast<size_t>(n + 1), 0.0);
    std::vector<double> s_local(static_cast<size_t>(n + 1), 0.0);
    for (int i = 1; i <= n; ++i) {
        const int idx = start + i - 1;
        x_local[i] = x[idx];
        s_local[i] = s[idx];
    }
    splind(x_local, xs_local, s_local, n, xs1, xs2);
    for (int i = 1; i <= n; ++i) {
        const int idx = start + i - 1;
        xs[idx] = xs_local[i];
    }
}

void SPLQSP(XFoilState &ctx, int kqsp) {
    if (ctx.NSP < 2) {
        return;
    }

    const int nmid = ctx.NSP - 2;
    if (nmid > 0) {
        std::vector<double> q_local(static_cast<size_t>(nmid + 1), 0.0);
        std::vector<double> s_local(static_cast<size_t>(nmid + 1), 0.0);
        std::vector<double> qp_local(static_cast<size_t>(nmid + 1), 0.0);
        for (int i = 1; i <= nmid; ++i) {
            const int idx = i + 1;
            q_local[i] = ctx.QSPEC[idx][kqsp];
            s_local[i] = ctx.SSPEC[idx];
        }
        splind(q_local, qp_local, s_local, nmid, -999.0, -999.0);
        for (int i = 1; i <= nmid; ++i) {
            const int idx = i + 1;
            ctx.QSPECP[idx][kqsp] = qp_local[i];
        }
    }

    std::vector<double> q_local(3, 0.0);
    std::vector<double> s_local(3, 0.0);
    std::vector<double> qp_local(3, 0.0);
    q_local[1] = ctx.QSPEC[1][kqsp];
    q_local[2] = ctx.QSPEC[2][kqsp];
    s_local[1] = ctx.SSPEC[1];
    s_local[2] = ctx.SSPEC[2];
    splind(q_local, qp_local, s_local, 2, -999.0, ctx.QSPECP[2][kqsp]);
    ctx.QSPECP[1][kqsp] = qp_local[1];

    const int i = ctx.NSP - 1;
    q_local[1] = ctx.QSPEC[i][kqsp];
    q_local[2] = ctx.QSPEC[i + 1][kqsp];
    s_local[1] = ctx.SSPEC[i];
    s_local[2] = ctx.SSPEC[i + 1];
    splind(q_local, qp_local, s_local, 2, ctx.QSPECP[i][kqsp], -999.0);
    ctx.QSPECP[i][kqsp] = qp_local[1];
}

void SMOOQ(XFoilState &ctx, int kq1, int kq2, int kqsp) {
    for (int i = 1; i <= ctx.NSP; ++i) {
        ctx.W8[i] = ctx.SSPEC[i];
    }

    if (kq2 - kq1 < 2) {
        return;
    }

    const double smool = 0.002 * (ctx.W8[ctx.NSP] - ctx.W8[1]);
    const double smoosq = smool * smool;

    for (int i = kq1 + 1; i <= kq2 - 1; ++i) {
        const double dsm = ctx.W8[i] - ctx.W8[i - 1];
        const double dsp = ctx.W8[i + 1] - ctx.W8[i];
        const double dso = 0.5 * (ctx.W8[i + 1] - ctx.W8[i - 1]);
        ctx.W1[i] = smoosq * (-1.0 / dsm) / dso;
        ctx.W2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
        ctx.W3[i] = smoosq * (-1.0 / dsp) / dso;
    }

    ctx.W2[kq1] = 1.0;
    ctx.W3[kq1] = 0.0;
    ctx.W1[kq2] = 0.0;
    ctx.W2[kq2] = 1.0;

    if (ctx.LQSLOP) {
        int i = kq1 + 1;
        double dsm = ctx.W8[i] - ctx.W8[i - 1];
        double dsp = ctx.W8[i + 1] - ctx.W8[i];
        double ds = ctx.W8[i + 1] - ctx.W8[i - 1];
        ctx.W1[i] = -1.0 / dsm - (dsm / ds) / dsm;
        ctx.W2[i] = 1.0 / dsm + (dsm / ds) / dsm + (dsm / ds) / dsp;
        ctx.W3[i] = -(dsm / ds) / dsp;
        const double qsp1 = ctx.W1[i] * ctx.QSPEC[i - 1][kqsp]
            + ctx.W2[i] * ctx.QSPEC[i][kqsp]
            + ctx.W3[i] * ctx.QSPEC[i + 1][kqsp];

        i = kq2 - 1;
        dsm = ctx.W8[i] - ctx.W8[i - 1];
        dsp = ctx.W8[i + 1] - ctx.W8[i];
        ds = ctx.W8[i + 1] - ctx.W8[i - 1];
        ctx.W1[i] = (dsp / ds) / dsm;
        ctx.W2[i] = -1.0 / dsp - (dsp / ds) / dsp - (dsp / ds) / dsm;
        ctx.W3[i] = 1.0 / dsp + (dsp / ds) / dsp;
        const double qsp2 = ctx.W1[i] * ctx.QSPEC[i - 1][kqsp]
            + ctx.W2[i] * ctx.QSPEC[i][kqsp]
            + ctx.W3[i] * ctx.QSPEC[i + 1][kqsp];

        ctx.QSPEC[kq1 + 1][kqsp] = qsp1;
        ctx.QSPEC[kq2 - 1][kqsp] = qsp2;
    }

    const int kk = kq2 - kq1 + 1;
    std::vector<double> a(static_cast<size_t>(kk + 1), 0.0);
    std::vector<double> b(static_cast<size_t>(kk + 1), 0.0);
    std::vector<double> c(static_cast<size_t>(kk + 1), 0.0);
    std::vector<double> d(static_cast<size_t>(kk + 1), 0.0);
    for (int i = 1; i <= kk; ++i) {
        const int idx = kq1 + i - 1;
        a[i] = ctx.W2[idx];
        b[i] = ctx.W1[idx];
        c[i] = ctx.W3[idx];
        d[i] = ctx.QSPEC[idx][kqsp];
    }

    trisol(a, b, c, d, kk);

    for (int i = 1; i <= kk; ++i) {
        const int idx = kq1 + i - 1;
        ctx.QSPEC[idx][kqsp] = d[i];
    }
}

double QINCOM(double qc, double qinf, double tklam) {
    if (tklam < 1.0e-4 || std::abs(qc) < 1.0e-4) {
        return qc / (1.0 - tklam);
    }
    const double tmp = 0.5 * (1.0 - tklam) * qinf / (qc * tklam);
    return qinf * tmp * (std::sqrt(1.0 + 1.0 / (tklam * tmp * tmp)) - 1.0);
}

void GAMQSP(XFoilState &ctx, int kqsp) {
    ctx.ALQSP[kqsp] = ctx.ALGAM;
    ctx.CLQSP[kqsp] = ctx.CLGAM;
    ctx.CMQSP[kqsp] = ctx.CMGAM;

    for (int i = 1; i <= ctx.NSP; ++i) {
        ctx.QSPEC[i][kqsp] = ctx.QGAMM[i];
    }

    ctx.QDOF0 = 0.0;
    ctx.QDOF1 = 0.0;
    ctx.QDOF2 = 0.0;
    ctx.QDOF3 = 0.0;

    SPLQSP(ctx, kqsp);

    if (!ctx.LIQSET) {
        ctx.IQ1 = 1;
        ctx.IQ2 = ctx.NSP;
    }
}

void SYMQSP(XFoilState &ctx, int kqsp) {
    ctx.ALQSP[kqsp] = 0.0;
    ctx.CLQSP[kqsp] = 0.0;
    ctx.CMQSP[kqsp] = 0.0;

    const double sspmid = 0.5 * (ctx.SSPEC[ctx.NSP] - ctx.SSPEC[1]);
    for (int i = 1; i <= (ctx.NSP + 1) / 2; ++i) {
        ctx.SSPEC[i] = sspmid + 0.5 * (ctx.SSPEC[i] - ctx.SSPEC[ctx.NSP - i + 1]);
        ctx.QSPEC[i][kqsp] = 0.5 * (ctx.QSPEC[i][kqsp] - ctx.QSPEC[ctx.NSP - i + 1][kqsp]);
    }
    for (int i = (ctx.NSP + 1) / 2 + 1; i <= ctx.NSP; ++i) {
        ctx.SSPEC[i] = -ctx.SSPEC[ctx.NSP - i + 1] + 2.0 * sspmid;
        ctx.QSPEC[i][kqsp] = -ctx.QSPEC[ctx.NSP - i + 1][kqsp];
    }

    ctx.QDOF0 = 0.0;
    ctx.QDOF1 = 0.0;
    ctx.QDOF2 = 0.0;
    ctx.QDOF3 = 0.0;

    SPLQSP(ctx, kqsp);
}

void MIXED(XFoilState &ctx, int kqsp, int niterq) {
    const double bwt = 0.1;

    ctx.COSA = std::cos(ctx.ALFA);
    ctx.SINA = std::sin(ctx.ALFA);
    scalc(ctx.X, ctx.Y, ctx.S, ctx.N);

    for (int i = 1; i <= ctx.N; ++i) {
        ctx.QF0[i] = 0.0;
        ctx.QF1[i] = 0.0;
        ctx.QF2[i] = 0.0;
        ctx.QF3[i] = 0.0;
    }

    for (int i = ctx.IQ1; i <= ctx.IQ2; ++i) {
        const double fs = (ctx.S[i] - ctx.S[ctx.IQ1]) / (ctx.S[ctx.IQ2] - ctx.S[ctx.IQ1]);
        ctx.QF0[i] = 1.0 - fs;
        ctx.QF1[i] = fs;
        if (ctx.LCPXX) {
            ctx.QF2[i] = std::exp(-5.0 * fs);
            ctx.QF3[i] = std::exp(-5.0 * (1.0 - fs));
        } else {
            ctx.QF2[i] = 0.0;
            ctx.QF3[i] = 0.0;
        }
        ctx.GAM[i] = ctx.QSPEC[i][kqsp]
            + ctx.QDOF0 * ctx.QF0[i]
            + ctx.QDOF1 * ctx.QF1[i]
            + ctx.QDOF2 * ctx.QF2[i]
            + ctx.QDOF3 * ctx.QF3[i];
    }

    for (int iter = 1; iter <= niterq; ++iter) {
        for (int i = 1; i <= ctx.N + 5; ++i) {
            for (int j = 1; j <= ctx.N + 5; ++j) {
                ctx.Q[i][j] = 0.0;
            }
        }

        ncalc(ctx.X, ctx.Y, ctx.S, ctx.N, ctx.NX, ctx.NY);

        for (int i = 1; i <= ctx.N; ++i) {
            const auto res = psilin(ctx, i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i], true, false);
            ctx.DZDN[i] = ctx.DZDN[i] + res.psiNi;

            for (int j = 1; j <= ctx.IQ1 - 1; ++j) {
                ctx.Q[i][j] = ctx.Q[i][j] + ctx.DZDG[j];
            }
            for (int j = ctx.IQ1; j <= ctx.IQ2; ++j) {
                ctx.Q[i][j] = ctx.Q[i][j] + ctx.DZDN[j];
            }
            for (int j = ctx.IQ2 + 1; j <= ctx.N; ++j) {
                ctx.Q[i][j] = ctx.Q[i][j] + ctx.DZDG[j];
            }

            ctx.DQ[i] = ctx.PSIO - res.psi;

            ctx.Q[i][ctx.N + 1] = ctx.Q[i][ctx.N + 1] - 1.0;
            ctx.Q[i][ctx.N + 2] = ctx.Q[i][ctx.N + 2] + ctx.Z_QDOF0;
            ctx.Q[i][ctx.N + 3] = ctx.Q[i][ctx.N + 3] + ctx.Z_QDOF1;
            ctx.Q[i][ctx.N + 4] = ctx.Q[i][ctx.N + 4] + ctx.Z_QDOF2;
            ctx.Q[i][ctx.N + 5] = ctx.Q[i][ctx.N + 5] + ctx.Z_QDOF3;
        }

        ctx.DQ[ctx.N + 1] = -(ctx.GAM[1] + ctx.GAM[ctx.N]);
        GAMLIN(ctx, ctx.N + 1, 1, 1.0);
        GAMLIN(ctx, ctx.N + 1, ctx.N, 1.0);

        if (ctx.SHARP) {
            const double ag1 = std::atan2(-ctx.YP[1], -ctx.XP[1]);
            const double ag2 = atanc(ctx.YP[ctx.N], ctx.XP[ctx.N], ag1);
            const double abis = 0.5 * (ag1 + ag2);
            const double cbis = std::cos(abis);
            const double sbis = std::sin(abis);

            const double ds1 = std::sqrt((ctx.X[1] - ctx.X[2]) * (ctx.X[1] - ctx.X[2])
                                         + (ctx.Y[1] - ctx.Y[2]) * (ctx.Y[1] - ctx.Y[2]));
            const double ds2 = std::sqrt((ctx.X[ctx.N] - ctx.X[ctx.N - 1]) * (ctx.X[ctx.N] - ctx.X[ctx.N - 1])
                                         + (ctx.Y[ctx.N] - ctx.Y[ctx.N - 1]) * (ctx.Y[ctx.N] - ctx.Y[ctx.N - 1]));
            const double dsmin = std::min(ds1, ds2);

            const double xbis = ctx.XTE - bwt * dsmin * cbis;
            const double ybis = ctx.YTE - bwt * dsmin * sbis;
            const auto bis = psilin(ctx, 0, xbis, ybis, -sbis, cbis, false, true);
            const double res = bis.psiNi;

            for (int j = 1; j <= ctx.N + 5; ++j) {
                ctx.Q[ctx.N][j] = 0.0;
            }

            for (int j = 1; j <= ctx.N; ++j) {
                GAMLIN(ctx, ctx.N, j, ctx.DQDG[j]);
                ctx.Q[ctx.N][j] = ctx.DQDG[j];
            }

            ctx.Q[ctx.N][ctx.N + 1] = 0.0;
            ctx.DQ[ctx.N] = -res;
        }

        ctx.Q[ctx.N + 2][ctx.IQ1] = 1.0;
        ctx.DQ[ctx.N + 2] = 0.0;
        ctx.Q[ctx.N + 3][ctx.IQ2] = 1.0;
        ctx.DQ[ctx.N + 3] = 0.0;

        if (ctx.IQ1 > 1 && ctx.LCPXX) {
            const double res = ctx.GAM[ctx.IQ1 - 1] - 2.0 * ctx.GAM[ctx.IQ1] + ctx.GAM[ctx.IQ1 + 1]
                - (ctx.QSPEC[ctx.IQ1 - 1][kqsp] - 2.0 * ctx.QSPEC[ctx.IQ1][kqsp] + ctx.QSPEC[ctx.IQ1 + 1][kqsp]);
            GAMLIN(ctx, ctx.N + 4, ctx.IQ1 - 1, 1.0);
            GAMLIN(ctx, ctx.N + 4, ctx.IQ1, -2.0);
            GAMLIN(ctx, ctx.N + 4, ctx.IQ1 + 1, 1.0);
            ctx.DQ[ctx.N + 4] = -res;
        } else {
            ctx.Q[ctx.N + 4][ctx.N + 4] = 1.0;
            ctx.DQ[ctx.N + 4] = -ctx.QDOF2;
        }

        if (ctx.IQ2 < ctx.N && ctx.LCPXX) {
            const double res = ctx.GAM[ctx.IQ2 - 1] - 2.0 * ctx.GAM[ctx.IQ2] + ctx.GAM[ctx.IQ2 + 1]
                - (ctx.QSPEC[ctx.IQ2 - 1][kqsp] - 2.0 * ctx.QSPEC[ctx.IQ2][kqsp] + ctx.QSPEC[ctx.IQ2 + 1][kqsp]);
            GAMLIN(ctx, ctx.N + 5, ctx.IQ2 - 1, 1.0);
            GAMLIN(ctx, ctx.N + 5, ctx.IQ2, -2.0);
            GAMLIN(ctx, ctx.N + 5, ctx.IQ2 + 1, 1.0);
            ctx.DQ[ctx.N + 5] = -res;
        } else {
            ctx.Q[ctx.N + 5][ctx.N + 5] = 1.0;
            ctx.DQ[ctx.N + 5] = -ctx.QDOF3;
        }

        gauss(IQX, ctx.N + 5, ctx.Q, ctx.DQ, 1);

        double dnmax = 0.0;
        double dgmax = 0.0;

        for (int i = 1; i <= ctx.IQ1 - 1; ++i) {
            ctx.GAM[i] = ctx.GAM[i] + ctx.DQ[i];
            if (std::abs(ctx.DQ[i]) > std::abs(dgmax)) {
                dgmax = ctx.DQ[i];
            }
        }

        for (int i = ctx.IQ1; i <= ctx.IQ2; ++i) {
            ctx.X[i] = ctx.X[i] + ctx.NX[i] * ctx.DQ[i];
            ctx.Y[i] = ctx.Y[i] + ctx.NY[i] * ctx.DQ[i];
            if (std::abs(ctx.DQ[i]) > std::abs(dnmax)) {
                dnmax = ctx.DQ[i];
            }
        }

        for (int i = ctx.IQ2 + 1; i <= ctx.N; ++i) {
            ctx.GAM[i] = ctx.GAM[i] + ctx.DQ[i];
            if (std::abs(ctx.DQ[i]) > std::abs(dgmax)) {
                dgmax = ctx.DQ[i];
            }
        }

        ctx.PSIO = ctx.PSIO + ctx.DQ[ctx.N + 1];
        ctx.QDOF0 = ctx.QDOF0 + ctx.DQ[ctx.N + 2];
        ctx.QDOF1 = ctx.QDOF1 + ctx.DQ[ctx.N + 3];
        ctx.QDOF2 = ctx.QDOF2 + ctx.DQ[ctx.N + 4];
        ctx.QDOF3 = ctx.QDOF3 + ctx.DQ[ctx.N + 5];

        ctx.COSA = std::cos(ctx.ALFA);
        ctx.SINA = std::sin(ctx.ALFA);
        scalc(ctx.X, ctx.Y, ctx.S, ctx.N);

        for (int i = ctx.IQ1; i <= ctx.IQ2; ++i) {
            ctx.GAM[i] = ctx.QSPEC[i][kqsp]
                + ctx.QDOF0 * ctx.QF0[i]
                + ctx.QDOF1 * ctx.QF1[i]
                + ctx.QDOF2 * ctx.QF2[i]
                + ctx.QDOF3 * ctx.QF3[i];
        }

        tecalc(ctx);
        const auto coeffs = clcalc(ctx.N, ctx.X, ctx.Y, ctx.GAM, ctx.GAM_A, ctx.ALFA, ctx.MINF, ctx.QINF, ctx.XCMREF, ctx.YCMREF);
        ctx.CL = std::get<0>(coeffs);
        ctx.CM = std::get<1>(coeffs);
        ctx.CDP = std::get<2>(coeffs);
        ctx.CL_ALF = std::get<3>(coeffs);
        ctx.CL_MSQ = std::get<4>(coeffs);

        if (std::abs(dnmax) < 5.0e-5 && std::abs(dgmax) < 5.0e-4) {
            return;
        }
    }
}

void GAMLIN(XFoilState &ctx, int i, int j, double coef) {
    if (j >= ctx.IQ1 && j <= ctx.IQ2) {
        ctx.Q[i][ctx.N + 2] = ctx.Q[i][ctx.N + 2] + coef * ctx.QF0[j];
        ctx.Q[i][ctx.N + 3] = ctx.Q[i][ctx.N + 3] + coef * ctx.QF1[j];
        ctx.Q[i][ctx.N + 4] = ctx.Q[i][ctx.N + 4] + coef * ctx.QF2[j];
        ctx.Q[i][ctx.N + 5] = ctx.Q[i][ctx.N + 5] + coef * ctx.QF3[j];
    } else {
        ctx.Q[i][j] = ctx.Q[i][j] + coef;
    }
}
