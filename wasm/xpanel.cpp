#include "xpanel.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "spline.h"
#include "xbl.h"
#include "xsolve.h"
#include "xutils.h"

void apcalc(XFoilState &ctx) {
    for (int i = 1; i <= ctx.N - 1; ++i) {
        const double sx = ctx.X[i + 1] - ctx.X[i];
        const double sy = ctx.Y[i + 1] - ctx.Y[i];
        if (sx == 0.0 && sy == 0.0) {
            ctx.APANEL[i] = std::atan2(-ctx.NY[i], -ctx.NX[i]);
        } else {
            ctx.APANEL[i] = std::atan2(sx, -sy);
        }
    }

    const int i = ctx.N;
    const int ip = 1;
    if (ctx.SHARP) {
        ctx.APANEL[i] = ctx.PI;
    } else {
        const double sx = ctx.X[ip] - ctx.X[i];
        const double sy = ctx.Y[ip] - ctx.Y[i];
        ctx.APANEL[i] = std::atan2(-sx, sy) + ctx.PI;
    }
}

void ncalc(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s, int n,
           std::vector<double> &xn, std::vector<double> &yn) {
    if (n <= 1) {
        return;
    }

    segspl(x, xn, s, n);
    segspl(y, yn, s, n);
    for (int i = 1; i <= n; ++i) {
        const double sx = yn[i];
        const double sy = -xn[i];
        const double smod = std::sqrt(sx * sx + sy * sy);
        if (smod == 0.0) {
            xn[i] = -1.0;
            yn[i] = 0.0;
        } else {
            xn[i] = sx / smod;
            yn[i] = sy / smod;
        }
    }

    for (int i = 1; i <= n - 1; ++i) {
        if (s[i] == s[i + 1]) {
            const double sx = 0.5 * (xn[i] + xn[i + 1]);
            const double sy = 0.5 * (yn[i] + yn[i + 1]);
            const double smod = std::sqrt(sx * sx + sy * sy);
            if (smod == 0.0) {
                xn[i] = -1.0;
                yn[i] = 0.0;
                xn[i + 1] = -1.0;
                yn[i + 1] = 0.0;
            } else {
                xn[i] = sx / smod;
                yn[i] = sy / smod;
                xn[i + 1] = sx / smod;
                yn[i + 1] = sy / smod;
            }
        }
    }
}

PsinResult psilin(XFoilState &ctx, int i, double xi, double yi, double nxi, double nyi, bool geolin, bool siglin) {
    const double seps = (ctx.S[ctx.N] - ctx.S[1]) * 1.0e-5;
    const int io = i;

    const double cosa = std::cos(ctx.ALFA);
    const double sina = std::sin(ctx.ALFA);

    for (int jo = 1; jo <= ctx.N; ++jo) {
        ctx.DZDG[jo] = 0.0;
        ctx.DZDN[jo] = 0.0;
        ctx.DQDG[jo] = 0.0;
    }

    for (int jo = 1; jo <= ctx.N; ++jo) {
        ctx.DZDM[jo] = 0.0;
        ctx.DQDM[jo] = 0.0;
    }

    ctx.Z_QINF = 0.0;
    ctx.Z_ALFA = 0.0;
    ctx.Z_QDOF0 = 0.0;
    ctx.Z_QDOF1 = 0.0;
    ctx.Z_QDOF2 = 0.0;
    ctx.Z_QDOF3 = 0.0;

    double psi = 0.0;
    double psi_ni = 0.0;

    double qt1 = 0.0;
    double qt2 = 0.0;
    double qtanm = 0.0;

    double scs = 1.0;
    double sds = 0.0;
    if (!ctx.SHARP) {
        scs = ctx.ANTE / ctx.DSTE;
        sds = ctx.ASTE / ctx.DSTE;
    }

    int jo_last = 0;
    int jp_last = 0;
    double apan = 0.0;
    double x1 = 0.0;
    double x2 = 0.0;
    double yy = 0.0;
    double g1 = 0.0;
    double g2 = 0.0;
    double t1 = 0.0;
    double t2 = 0.0;
    double x1i = 0.0;
    double x2i = 0.0;
    double yyi = 0.0;
    double x1o = 0.0;
    double x1p = 0.0;
    double x2o = 0.0;
    double x2p = 0.0;
    double yyo = 0.0;
    double yyp = 0.0;

    for (int jo = 1; jo <= ctx.N; ++jo) {
        int jp = jo + 1;
        int jm = jo - 1;
        int jq = jp + 1;

        if (jo == 1) {
            jm = jo;
        } else if (jo == ctx.N - 1) {
            jq = jp;
        } else if (jo == ctx.N) {
            jp = 1;
            if ((ctx.X[jo] - ctx.X[jp]) * (ctx.X[jo] - ctx.X[jp]) + (ctx.Y[jo] - ctx.Y[jp]) * (ctx.Y[jo] - ctx.Y[jp])
                < seps * seps) {
                continue;
            }
        }

        const double dso = std::sqrt((ctx.X[jo] - ctx.X[jp]) * (ctx.X[jo] - ctx.X[jp])
                                     + (ctx.Y[jo] - ctx.Y[jp]) * (ctx.Y[jo] - ctx.Y[jp]));
        if (dso == 0.0) {
            continue;
        }

        const double dsio = 1.0 / dso;
        apan = ctx.APANEL[jo];
        jo_last = jo;
        jp_last = jp;

        const double rx1 = xi - ctx.X[jo];
        const double ry1 = yi - ctx.Y[jo];
        const double rx2 = xi - ctx.X[jp];
        const double ry2 = yi - ctx.Y[jp];

        const double sx = (ctx.X[jp] - ctx.X[jo]) * dsio;
        const double sy = (ctx.Y[jp] - ctx.Y[jo]) * dsio;

        x1 = sx * rx1 + sy * ry1;
        x2 = sx * rx2 + sy * ry2;
        yy = sx * ry1 - sy * rx1;

        const double rs1 = rx1 * rx1 + ry1 * ry1;
        const double rs2 = rx2 * rx2 + ry2 * ry2;

        double sgn = 0.0;
        if (io >= 1 && io <= ctx.N) {
            sgn = 1.0;
        } else {
            sgn = std::copysign(1.0, yy);
        }

        g1 = 0.0;
        t1 = 0.0;
        if (io != jo && rs1 > 0.0) {
            g1 = std::log(rs1);
            t1 = std::atan2(sgn * x1, sgn * yy) + (0.5 - 0.5 * sgn) * ctx.PI;
        }

        g2 = 0.0;
        t2 = 0.0;
        if (io != jp && rs2 > 0.0) {
            g2 = std::log(rs2);
            t2 = std::atan2(sgn * x2, sgn * yy) + (0.5 - 0.5 * sgn) * ctx.PI;
        }

        x1i = sx * nxi + sy * nyi;
        x2i = sx * nxi + sy * nyi;
        yyi = sx * nyi - sy * nxi;

        if (geolin) {
            const double nxo = ctx.NX[jo];
            const double nyo = ctx.NY[jo];
            const double nxp = ctx.NX[jp];
            const double nyp = ctx.NY[jp];

            x1o = -((rx1 - x1 * sx) * nxo + (ry1 - x1 * sy) * nyo) * dsio - (sx * nxo + sy * nyo);
            x1p = ((rx1 - x1 * sx) * nxp + (ry1 - x1 * sy) * nyp) * dsio;
            x2o = -((rx2 - x2 * sx) * nxo + (ry2 - x2 * sy) * nyo) * dsio;
            x2p = ((rx2 - x2 * sx) * nxp + (ry2 - x2 * sy) * nyp) * dsio - (sx * nxp + sy * nyp);
            yyo = ((rx1 + x1 * sy) * nyo - (ry1 - x1 * sx) * nxo) * dsio - (sx * nyo - sy * nxo);
            yyp = -((rx1 - x1 * sy) * nyp - (ry1 + x1 * sx) * nxp) * dsio;
        }

        if (jo == ctx.N) {
            continue;
        }

        if (siglin) {
            const double x0 = 0.5 * (x1 + x2);
            const double rs0 = x0 * x0 + yy * yy;
            const double g0 = std::log(rs0);
            const double t0 = std::atan2(sgn * x0, sgn * yy) + (0.5 - 0.5 * sgn) * ctx.PI;

            double dxinv = 1.0 / (x1 - x0);
            double psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
            double pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy) * dxinv;

            double psx1 = -(t1 - apan);
            double psx0 = t0 - apan;
            double psyy = 0.5 * (g1 - g0);

            double pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
            double pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
            double pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;

            const double dsm = std::sqrt((ctx.X[jp] - ctx.X[jm]) * (ctx.X[jp] - ctx.X[jm])
                                         + (ctx.Y[jp] - ctx.Y[jm]) * (ctx.Y[jp] - ctx.Y[jm]));
            const double dsim = 1.0 / dsm;

            double ssum = (ctx.SIG[jp] - ctx.SIG[jo]) * dsio + (ctx.SIG[jp] - ctx.SIG[jm]) * dsim;
            double sdif = (ctx.SIG[jp] - ctx.SIG[jo]) * dsio - (ctx.SIG[jp] - ctx.SIG[jm]) * dsim;

            psi = psi + ctx.QOPI * (psum * ssum + pdif * sdif);

            ctx.DZDM[jm] = ctx.DZDM[jm] + ctx.QOPI * (-psum * dsim + pdif * dsim);
            ctx.DZDM[jo] = ctx.DZDM[jo] + ctx.QOPI * (-psum * dsio - pdif * dsio);
            ctx.DZDM[jp] = ctx.DZDM[jp] + ctx.QOPI * (psum * (dsio + dsim) + pdif * (dsio - dsim));

            double psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
            double pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
            psi_ni = psi_ni + ctx.QOPI * (psni * ssum + pdni * sdif);

            qtanm = qtanm + ctx.QOPI * (psni * ssum + pdni * sdif);

            ctx.DQDM[jm] = ctx.DQDM[jm] + ctx.QOPI * (-psni * dsim + pdni * dsim);
            ctx.DQDM[jo] = ctx.DQDM[jo] + ctx.QOPI * (-psni * dsio - pdni * dsio);
            ctx.DQDM[jp] = ctx.DQDM[jp] + ctx.QOPI * (psni * (dsio + dsim) + pdni * (dsio - dsim));

            dxinv = 1.0 / (x0 - x2);
            psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
            pdif = ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy) * dxinv;

            psx0 = -(t0 - apan);
            double psx2 = t2 - apan;
            psyy = 0.5 * (g0 - g2);

            pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
            double pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
            pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv;

            const double dsp = std::sqrt((ctx.X[jq] - ctx.X[jo]) * (ctx.X[jq] - ctx.X[jo])
                                         + (ctx.Y[jq] - ctx.Y[jo]) * (ctx.Y[jq] - ctx.Y[jo]));
            const double dsip = 1.0 / dsp;

            ssum = (ctx.SIG[jq] - ctx.SIG[jo]) * dsip + (ctx.SIG[jp] - ctx.SIG[jo]) * dsio;
            sdif = (ctx.SIG[jq] - ctx.SIG[jo]) * dsip - (ctx.SIG[jp] - ctx.SIG[jo]) * dsio;

            psi = psi + ctx.QOPI * (psum * ssum + pdif * sdif);

            ctx.DZDM[jo] = ctx.DZDM[jo] + ctx.QOPI * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
            ctx.DZDM[jp] = ctx.DZDM[jp] + ctx.QOPI * (psum * dsio - pdif * dsio);
            ctx.DZDM[jq] = ctx.DZDM[jq] + ctx.QOPI * (psum * dsip + pdif * dsip);

            psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
            pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
            psi_ni = psi_ni + ctx.QOPI * (psni * ssum + pdni * sdif);

            qtanm = qtanm + ctx.QOPI * (psni * ssum + pdni * sdif);

            ctx.DQDM[jo] = ctx.DQDM[jo] + ctx.QOPI * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
            ctx.DQDM[jp] = ctx.DQDM[jp] + ctx.QOPI * (psni * dsio - pdni * dsio);
            ctx.DQDM[jq] = ctx.DQDM[jq] + ctx.QOPI * (psni * dsip + pdni * dsip);
        }

        const double dxinv = 1.0 / (x1 - x2);
        const double psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
        const double psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;

        const double psx1 = 0.5 * g1;
        const double psx2 = -0.5 * g2;
        const double psyy = t1 - t2;

        const double pdx1 = ((x1 + x2) * psx1 + psis - x1 * g1 - psid) * dxinv;
        const double pdx2 = ((x1 + x2) * psx2 + psis + x2 * g2 + psid) * dxinv;
        const double pdyy = ((x1 + x2) * psyy - yy * (g1 - g2)) * dxinv;

        const double gsum1 = ctx.GAMU[jp][1] + ctx.GAMU[jo][1];
        const double gsum2 = ctx.GAMU[jp][2] + ctx.GAMU[jo][2];
        const double gdif1 = ctx.GAMU[jp][1] - ctx.GAMU[jo][1];
        const double gdif2 = ctx.GAMU[jp][2] - ctx.GAMU[jo][2];

        const double gsum = ctx.GAM[jp] + ctx.GAM[jo];
        const double gdif = ctx.GAM[jp] - ctx.GAM[jo];

        psi = psi + ctx.QOPI * (psis * gsum + psid * gdif);

        ctx.DZDG[jo] = ctx.DZDG[jo] + ctx.QOPI * (psis - psid);
        ctx.DZDG[jp] = ctx.DZDG[jp] + ctx.QOPI * (psis + psid);

        const double psni = psx1 * x1i + psx2 * x2i + psyy * yyi;
        const double pdni = pdx1 * x1i + pdx2 * x2i + pdyy * yyi;
        psi_ni = psi_ni + ctx.QOPI * (gsum * psni + gdif * pdni);

        qt1 = qt1 + ctx.QOPI * (gsum1 * psni + gdif1 * pdni);
        qt2 = qt2 + ctx.QOPI * (gsum2 * psni + gdif2 * pdni);

        ctx.DQDG[jo] = ctx.DQDG[jo] + ctx.QOPI * (psni - pdni);
        ctx.DQDG[jp] = ctx.DQDG[jp] + ctx.QOPI * (psni + pdni);

        if (geolin) {
            ctx.DZDN[jo] = ctx.DZDN[jo] + ctx.QOPI * gsum * (psx1 * x1o + psx2 * x2o + psyy * yyo)
                           + ctx.QOPI * gdif * (pdx1 * x1o + pdx2 * x2o + pdyy * yyo);
            ctx.DZDN[jp] = ctx.DZDN[jp] + ctx.QOPI * gsum * (psx1 * x1p + psx2 * x2p + psyy * yyp)
                           + ctx.QOPI * gdif * (pdx1 * x1p + pdx2 * x2p + pdyy * yyp);
            ctx.Z_QDOF0 = ctx.Z_QDOF0 + ctx.QOPI * ((psis - psid) * ctx.QF0[jo] + (psis + psid) * ctx.QF0[jp]);
            ctx.Z_QDOF1 = ctx.Z_QDOF1 + ctx.QOPI * ((psis - psid) * ctx.QF1[jo] + (psis + psid) * ctx.QF1[jp]);
            ctx.Z_QDOF2 = ctx.Z_QDOF2 + ctx.QOPI * ((psis - psid) * ctx.QF2[jo] + (psis + psid) * ctx.QF2[jp]);
            ctx.Z_QDOF3 = ctx.Z_QDOF3 + ctx.QOPI * ((psis - psid) * ctx.QF3[jo] + (psis + psid) * ctx.QF3[jp]);
        }
    }

    const double psig = 0.5 * yy * (g1 - g2) + x2 * (t2 - apan) - x1 * (t1 - apan);
    const double pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);

    const double psigx1 = -(t1 - apan);
    const double psigx2 = t2 - apan;
    const double psigyy = 0.5 * (g1 - g2);
    const double pgamx1 = 0.5 * g1;
    const double pgamx2 = -0.5 * g2;
    const double pgamyy = t1 - t2;

    const double psigni = psigx1 * x1i + psigx2 * x2i + psigyy * yyi;
    const double pgamni = pgamx1 * x1i + pgamx2 * x2i + pgamyy * yyi;

    const double sigte1 = 0.5 * scs * (ctx.GAMU[jp_last][1] - ctx.GAMU[jo_last][1]);
    const double sigte2 = 0.5 * scs * (ctx.GAMU[jp_last][2] - ctx.GAMU[jo_last][2]);
    const double gamte1 = -0.5 * sds * (ctx.GAMU[jp_last][1] - ctx.GAMU[jo_last][1]);
    const double gamte2 = -0.5 * sds * (ctx.GAMU[jp_last][2] - ctx.GAMU[jo_last][2]);

    ctx.SIGTE = 0.5 * scs * (ctx.GAM[jp_last] - ctx.GAM[jo_last]);
    ctx.GAMTE = -0.5 * sds * (ctx.GAM[jp_last] - ctx.GAM[jo_last]);

    psi = psi + ctx.HOPI * (psig * ctx.SIGTE + pgam * ctx.GAMTE);

    ctx.DZDG[jo_last] = ctx.DZDG[jo_last] - ctx.HOPI * psig * scs * 0.5;
    ctx.DZDG[jp_last] = ctx.DZDG[jp_last] + ctx.HOPI * psig * scs * 0.5;

    ctx.DZDG[jo_last] = ctx.DZDG[jo_last] + ctx.HOPI * pgam * sds * 0.5;
    ctx.DZDG[jp_last] = ctx.DZDG[jp_last] - ctx.HOPI * pgam * sds * 0.5;

    psi_ni = psi_ni + ctx.HOPI * (psigni * ctx.SIGTE + pgamni * ctx.GAMTE);

    qt1 = qt1 + ctx.HOPI * (psigni * sigte1 + pgamni * gamte1);
    qt2 = qt2 + ctx.HOPI * (psigni * sigte2 + pgamni * gamte2);

    ctx.DQDG[jo_last] = ctx.DQDG[jo_last] - ctx.HOPI * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
    ctx.DQDG[jp_last] = ctx.DQDG[jp_last] + ctx.HOPI * (psigni * 0.5 * scs - pgamni * 0.5 * sds);

    if (geolin) {
        ctx.DZDN[jo_last] = ctx.DZDN[jo_last] + ctx.HOPI * (psigx1 * x1o + psigx2 * x2o + psigyy * yyo) * ctx.SIGTE
                            + ctx.HOPI * (pgamx1 * x1o + pgamx2 * x2o + pgamyy * yyo) * ctx.GAMTE;
        ctx.DZDN[jp_last] = ctx.DZDN[jp_last] + ctx.HOPI * (psigx1 * x1p + psigx2 * x2p + psigyy * yyp) * ctx.SIGTE
                            + ctx.HOPI * (pgamx1 * x1p + pgamx2 * x2p + pgamyy * yyp) * ctx.GAMTE;

        ctx.Z_QDOF0 = ctx.Z_QDOF0 + ctx.HOPI * psig * 0.5 * (ctx.QF0[jp_last] - ctx.QF0[jo_last]) * scs
                      - ctx.HOPI * pgam * 0.5 * (ctx.QF0[jp_last] - ctx.QF0[jo_last]) * sds;
        ctx.Z_QDOF1 = ctx.Z_QDOF1 + ctx.HOPI * psig * 0.5 * (ctx.QF1[jp_last] - ctx.QF1[jo_last]) * scs
                      - ctx.HOPI * pgam * 0.5 * (ctx.QF1[jp_last] - ctx.QF1[jo_last]) * sds;
        ctx.Z_QDOF2 = ctx.Z_QDOF2 + ctx.HOPI * psig * 0.5 * (ctx.QF2[jp_last] - ctx.QF2[jo_last]) * scs
                      - ctx.HOPI * pgam * 0.5 * (ctx.QF2[jp_last] - ctx.QF2[jo_last]) * sds;
        ctx.Z_QDOF3 = ctx.Z_QDOF3 + ctx.HOPI * psig * 0.5 * (ctx.QF3[jp_last] - ctx.QF3[jo_last]) * scs
                      - ctx.HOPI * pgam * 0.5 * (ctx.QF3[jp_last] - ctx.QF3[jo_last]) * sds;
    }

    psi = psi + ctx.QINF * (cosa * yi - sina * xi);
    psi_ni = psi_ni + ctx.QINF * (cosa * nyi - sina * nxi);

    qt1 = qt1 + ctx.QINF * nyi;
    qt2 = qt2 - ctx.QINF * nxi;

    ctx.Z_QINF = ctx.Z_QINF + (cosa * yi - sina * xi);
    ctx.Z_ALFA = ctx.Z_ALFA - ctx.QINF * (sina * yi + cosa * xi);

    if (ctx.LIMAGE) {
        for (int jo = 1; jo <= ctx.N; ++jo) {
            int jp = jo + 1;
            int jm = jo - 1;
            int jq = jp + 1;

            if (jo == 1) {
                jm = jo;
            } else if (jo == ctx.N - 1) {
                jq = jp;
            } else if (jo == ctx.N) {
                jp = 1;
                if ((ctx.X[jo] - ctx.X[jp]) * (ctx.X[jo] - ctx.X[jp]) + (ctx.Y[jo] - ctx.Y[jp]) * (ctx.Y[jo] - ctx.Y[jp])
                    < seps * seps) {
                    continue;
                }
            }

            const double dso = std::sqrt((ctx.X[jo] - ctx.X[jp]) * (ctx.X[jo] - ctx.X[jp])
                                         + (ctx.Y[jo] - ctx.Y[jp]) * (ctx.Y[jo] - ctx.Y[jp]));
            if (dso == 0.0) {
                continue;
            }

            const double dsio = 1.0 / dso;
            const double apan = ctx.PI - ctx.APANEL[jo] + 2.0 * ctx.ALFA;

            const double xjo = ctx.X[jo] + 2.0 * (ctx.YIMAGE + ctx.Y[jo]) * sina;
            const double yjo = ctx.Y[jo] - 2.0 * (ctx.YIMAGE + ctx.Y[jo]) * cosa;
            const double xjp = ctx.X[jp] + 2.0 * (ctx.YIMAGE + ctx.Y[jp]) * sina;
            const double yjp = ctx.Y[jp] - 2.0 * (ctx.YIMAGE + ctx.Y[jp]) * cosa;

            const double rx1 = xi - xjo;
            const double ry1 = yi - yjo;
            const double rx2 = xi - xjp;
            const double ry2 = yi - yjp;

            const double sx = (xjp - xjo) * dsio;
            const double sy = (yjp - yjo) * dsio;

            const double x1 = sx * rx1 + sy * ry1;
            const double x2 = sx * rx2 + sy * ry2;
            const double yy = sx * ry1 - sy * rx1;

            const double rs1 = rx1 * rx1 + ry1 * ry1;
            const double rs2 = rx2 * rx2 + ry2 * ry2;

            double sgn = 0.0;
            if (io >= 1 && io <= ctx.N) {
                sgn = 1.0;
            } else {
                sgn = std::copysign(1.0, yy);
            }

            double g1 = 0.0;
            double t1 = 0.0;
            if (io != jo && rs1 > 0.0) {
                g1 = std::log(rs1);
                t1 = std::atan2(sgn * x1, sgn * yy) + (0.5 - 0.5 * sgn) * ctx.PI;
            }

            double g2 = 0.0;
            double t2 = 0.0;
            if (io != jp && rs2 > 0.0) {
                g2 = std::log(rs2);
                t2 = std::atan2(sgn * x2, sgn * yy) + (0.5 - 0.5 * sgn) * ctx.PI;
            }

            const double x1i = sx * nxi + sy * nyi;
            const double x2i = sx * nxi + sy * nyi;
            const double yyi = sx * nyi - sy * nxi;

            double x1o = 0.0;
            double x1p = 0.0;
            double x2o = 0.0;
            double x2p = 0.0;
            double yyo = 0.0;
            double yyp = 0.0;

            if (geolin) {
                const double nxo = ctx.NX[jo];
                const double nyo = ctx.NY[jo];
                const double nxp = ctx.NX[jp];
                const double nyp = ctx.NY[jp];

                x1o = -((rx1 - x1 * sx) * nxo + (ry1 - x1 * sy) * nyo) * dsio - (sx * nxo + sy * nyo);
                x1p = ((rx1 - x1 * sx) * nxp + (ry1 - x1 * sy) * nyp) * dsio;
                x2o = -((rx2 - x2 * sx) * nxo + (ry2 - x2 * sy) * nyo) * dsio;
                x2p = ((rx2 - x2 * sx) * nxp + (ry2 - x2 * sy) * nyp) * dsio - (sx * nxp + sy * nyp);
                yyo = ((rx1 + x1 * sy) * nyo - (ry1 - x1 * sx) * nxo) * dsio - (sx * nyo - sy * nxo);
                yyp = -((rx1 - x1 * sy) * nyp - (ry1 + x1 * sx) * nxp) * dsio;
            }

            if (jo == ctx.N) {
                continue;
            }

            if (siglin) {
                const double x0 = 0.5 * (x1 + x2);
                const double rs0 = x0 * x0 + yy * yy;
                const double g0 = std::log(rs0);
                const double t0 = std::atan2(sgn * x0, sgn * yy) + (0.5 - 0.5 * sgn) * ctx.PI;

                double dxinv = 1.0 / (x1 - x0);
                double psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
                double pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy) * dxinv;

                double psx1 = -(t1 - apan);
                double psx0 = t0 - apan;
                double psyy = 0.5 * (g1 - g0);

                double pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
                double pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
                double pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;

                const double dsm = std::sqrt((ctx.X[jp] - ctx.X[jm]) * (ctx.X[jp] - ctx.X[jm])
                                             + (ctx.Y[jp] - ctx.Y[jm]) * (ctx.Y[jp] - ctx.Y[jm]));
                const double dsim = 1.0 / dsm;

                double ssum = (ctx.SIG[jp] - ctx.SIG[jo]) * dsio + (ctx.SIG[jp] - ctx.SIG[jm]) * dsim;
                double sdif = (ctx.SIG[jp] - ctx.SIG[jo]) * dsio - (ctx.SIG[jp] - ctx.SIG[jm]) * dsim;

                psi = psi + ctx.QOPI * (psum * ssum + pdif * sdif);

                ctx.DZDM[jm] = ctx.DZDM[jm] + ctx.QOPI * (-psum * dsim + pdif * dsim);
                ctx.DZDM[jo] = ctx.DZDM[jo] + ctx.QOPI * (-psum * dsio - pdif * dsio);
                ctx.DZDM[jp] = ctx.DZDM[jp] + ctx.QOPI * (psum * (dsio + dsim) + pdif * (dsio - dsim));

                double psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
                double pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
                psi_ni = psi_ni + ctx.QOPI * (psni * ssum + pdni * sdif);

                ctx.DQDM[jm] = ctx.DQDM[jm] + ctx.QOPI * (-psni * dsim + pdni * dsim);
                ctx.DQDM[jo] = ctx.DQDM[jo] + ctx.QOPI * (-psni * dsio - pdni * dsio);
                ctx.DQDM[jp] = ctx.DQDM[jp] + ctx.QOPI * (psni * (dsio + dsim) + pdni * (dsio - dsim));

                dxinv = 1.0 / (x0 - x2);
                psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
                pdif = ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy) * dxinv;

                psx0 = -(t0 - apan);
                double psx2 = t2 - apan;
                psyy = 0.5 * (g0 - g2);

                pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
                double pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
                pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv;

                const double dsp = std::sqrt((ctx.X[jq] - ctx.X[jo]) * (ctx.X[jq] - ctx.X[jo])
                                             + (ctx.Y[jq] - ctx.Y[jo]) * (ctx.Y[jq] - ctx.Y[jo]));
                const double dsip = 1.0 / dsp;

                ssum = (ctx.SIG[jq] - ctx.SIG[jo]) * dsip + (ctx.SIG[jp] - ctx.SIG[jo]) * dsio;
                sdif = (ctx.SIG[jq] - ctx.SIG[jo]) * dsip - (ctx.SIG[jp] - ctx.SIG[jo]) * dsio;

                psi = psi + ctx.QOPI * (psum * ssum + pdif * sdif);

                ctx.DZDM[jo] = ctx.DZDM[jo] + ctx.QOPI * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
                ctx.DZDM[jp] = ctx.DZDM[jp] + ctx.QOPI * (psum * dsio - pdif * dsio);
                ctx.DZDM[jq] = ctx.DZDM[jq] + ctx.QOPI * (psum * dsip + pdif * dsip);

                psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
                pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
                psi_ni = psi_ni + ctx.QOPI * (psni * ssum + pdni * sdif);

                ctx.DQDM[jo] = ctx.DQDM[jo] + ctx.QOPI * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
                ctx.DQDM[jp] = ctx.DQDM[jp] + ctx.QOPI * (psni * dsio - pdni * dsio);
                ctx.DQDM[jq] = ctx.DQDM[jq] + ctx.QOPI * (psni * dsip + pdni * dsip);
            }

            const double dxinv2 = 1.0 / (x1 - x2);
            const double psis2 = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
            const double psid2 = ((x1 + x2) * psis2 + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv2;

            const double psx1b = 0.5 * g1;
            const double psx2b = -0.5 * g2;
            const double psyyb = t1 - t2;

            const double pdx1b = ((x1 + x2) * psx1b + psis2 - x1 * g1 - psid2) * dxinv2;
            const double pdx2b = ((x1 + x2) * psx2b + psis2 + x2 * g2 + psid2) * dxinv2;
            const double pdyyb = ((x1 + x2) * psyyb - yy * (g1 - g2)) * dxinv2;

            const double gsum1 = ctx.GAMU[jp][1] + ctx.GAMU[jo][1];
            const double gsum2 = ctx.GAMU[jp][2] + ctx.GAMU[jo][2];
            const double gdif1 = ctx.GAMU[jp][1] - ctx.GAMU[jo][1];
            const double gdif2 = ctx.GAMU[jp][2] - ctx.GAMU[jo][2];

            const double gsum = ctx.GAM[jp] + ctx.GAM[jo];
            const double gdif = ctx.GAM[jp] - ctx.GAM[jo];

            psi = psi + ctx.QOPI * (psis2 * gsum + psid2 * gdif);

            ctx.DZDG[jo] = ctx.DZDG[jo] + ctx.QOPI * (psis2 - psid2);
            ctx.DZDG[jp] = ctx.DZDG[jp] + ctx.QOPI * (psis2 + psid2);

            const double psni2 = psx1b * x1i + psx2b * x2i + psyyb * yyi;
            const double pdni2 = pdx1b * x1i + pdx2b * x2i + pdyyb * yyi;
            psi_ni = psi_ni + ctx.QOPI * (gsum * psni2 + gdif * pdni2);

            qt1 = qt1 + ctx.QOPI * (gsum1 * psni2 + gdif1 * pdni2);
            qt2 = qt2 + ctx.QOPI * (gsum2 * psni2 + gdif2 * pdni2);

            ctx.DQDG[jo] = ctx.DQDG[jo] + ctx.QOPI * (psni2 - pdni2);
            ctx.DQDG[jp] = ctx.DQDG[jp] + ctx.QOPI * (psni2 + pdni2);

            if (geolin) {
                ctx.DZDN[jo] = ctx.DZDN[jo] + ctx.QOPI * gsum * (psx1b * x1o + psx2b * x2o + psyyb * yyo)
                               + ctx.QOPI * gdif * (pdx1b * x1o + pdx2b * x2o + pdyyb * yyo);
                ctx.DZDN[jp] = ctx.DZDN[jp] + ctx.QOPI * gsum * (psx1b * x1p + psx2b * x2p + psyyb * yyp)
                               + ctx.QOPI * gdif * (pdx1b * x1p + pdx2b * x2p + pdyyb * yyp);
                ctx.Z_QDOF0 = ctx.Z_QDOF0 + ctx.QOPI * ((psis2 - psid2) * ctx.QF0[jo] + (psis2 + psid2) * ctx.QF0[jp]);
                ctx.Z_QDOF1 = ctx.Z_QDOF1 + ctx.QOPI * ((psis2 - psid2) * ctx.QF1[jo] + (psis2 + psid2) * ctx.QF1[jp]);
                ctx.Z_QDOF2 = ctx.Z_QDOF2 + ctx.QOPI * ((psis2 - psid2) * ctx.QF2[jo] + (psis2 + psid2) * ctx.QF2[jp]);
                ctx.Z_QDOF3 = ctx.Z_QDOF3 + ctx.QOPI * ((psis2 - psid2) * ctx.QF3[jo] + (psis2 + psid2) * ctx.QF3[jp]);
            }
        }
    }

    return {psi, psi_ni, qt1, qt2, qtanm};
}

std::pair<double, double> pswlin(XFoilState &ctx, int i, double xi, double yi, double nxi, double nyi) {
    const bool siglin = true;
    const bool geolin = false;
    const PsinResult res = psilin(ctx, i, xi, yi, nxi, nyi, geolin, siglin);
    return {res.psi, res.psi_ni};
}

void ggcalc(XFoilState &ctx) {
    const double bwt = 0.1;

    std::cout << "Calculating unit vorticity distributions ..." << std::endl;

    for (int i = 1; i <= ctx.N; ++i) {
        ctx.GAM[i] = 0.0;
        ctx.GAMU[i][1] = 0.0;
        ctx.GAMU[i][2] = 0.0;
    }

    double psio = 0.0;

    for (int i = 1; i <= ctx.N; ++i) {
        PsinResult res = psilin(ctx, i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i], false, true);
        const double psiinf = ctx.QINF * (std::cos(ctx.ALFA) * ctx.Y[i] - std::sin(ctx.ALFA) * ctx.X[i]);

        const double res1 = ctx.QINF * ctx.Y[i];
        const double res2 = -ctx.QINF * ctx.X[i];

        for (int j = 1; j <= ctx.N; ++j) {
            ctx.AIJ[i][j] = ctx.DZDG[j];
        }

        for (int j = 1; j <= ctx.N; ++j) {
            ctx.BIJ[i][j] = -ctx.DZDM[j];
        }

        ctx.AIJ[i][ctx.N + 1] = -1.0;

        ctx.GAMU[i][1] = -res1;
        ctx.GAMU[i][2] = -res2;

        psio = res.psi;
        (void)psiinf;
    }

    double res = 0.0;

    for (int j = 1; j <= ctx.N + 1; ++j) {
        ctx.AIJ[ctx.N + 1][j] = 0.0;
    }

    ctx.AIJ[ctx.N + 1][1] = 1.0;
    ctx.AIJ[ctx.N + 1][ctx.N] = 1.0;

    ctx.GAMU[ctx.N + 1][1] = -res;
    ctx.GAMU[ctx.N + 1][2] = -res;

    for (int j = 1; j <= ctx.N; ++j) {
        ctx.BIJ[ctx.N + 1][j] = 0.0;
    }

    if (ctx.SHARP) {
        const double ag1 = std::atan2(-ctx.YP[1], -ctx.XP[1]);
        const double ag2 = atanc(ctx.YP[ctx.N], ctx.XP[ctx.N], ag1);
        const double abis = 0.5 * (ag1 + ag2);
        const double cbis = std::cos(abis);
        const double sbis = std::sin(abis);

        const double ds1 = std::sqrt((ctx.X[1] - ctx.X[2]) * (ctx.X[1] - ctx.X[2]) + (ctx.Y[1] - ctx.Y[2]) * (ctx.Y[1] - ctx.Y[2]));
        const double ds2 = std::sqrt((ctx.X[ctx.N] - ctx.X[ctx.N - 1]) * (ctx.X[ctx.N] - ctx.X[ctx.N - 1])
                                     + (ctx.Y[ctx.N] - ctx.Y[ctx.N - 1]) * (ctx.Y[ctx.N] - ctx.Y[ctx.N - 1]));
        const double dsmin = std::min(ds1, ds2);

        const double xbis = ctx.XTE - bwt * dsmin * cbis;
        const double ybis = ctx.YTE - bwt * dsmin * sbis;

        PsinResult resb = psilin(ctx, 0, xbis, ybis, -sbis, cbis, false, true);
        const double qbis = resb.psi_ni;

        res = qbis;

        for (int j = 1; j <= ctx.N; ++j) {
            ctx.AIJ[ctx.N][j] = ctx.DQDG[j];
        }

        for (int j = 1; j <= ctx.N; ++j) {
            ctx.BIJ[ctx.N][j] = -ctx.DQDM[j];
        }

        ctx.AIJ[ctx.N][ctx.N + 1] = 0.0;

        ctx.GAMU[ctx.N][1] = -cbis;
        ctx.GAMU[ctx.N][2] = -sbis;
    }

    ludcmp(static_cast<int>(ctx.AIJ.size()) - 1, ctx.N + 1, ctx.AIJ, ctx.AIJPIV);
    ctx.LQAIJ = true;

    std::vector<double> col1(ctx.N + 2, 0.0);
    std::vector<double> col2(ctx.N + 2, 0.0);
    for (int i = 1; i <= ctx.N + 1; ++i) {
        col1[i] = ctx.GAMU[i][1];
        col2[i] = ctx.GAMU[i][2];
    }
    baksub(static_cast<int>(ctx.AIJ.size()) - 1, ctx.N + 1, ctx.AIJ, ctx.AIJPIV, col1);
    baksub(static_cast<int>(ctx.AIJ.size()) - 1, ctx.N + 1, ctx.AIJ, ctx.AIJPIV, col2);
    for (int i = 1; i <= ctx.N + 1; ++i) {
        ctx.GAMU[i][1] = col1[i];
        ctx.GAMU[i][2] = col2[i];
    }

    for (int i = 1; i <= ctx.N; ++i) {
        ctx.QINVU[i][1] = ctx.GAMU[i][1];
        ctx.QINVU[i][2] = ctx.GAMU[i][2];
    }

    ctx.LGAMU = true;
    (void)psio;
}

void qwcalc(XFoilState &ctx) {
    ctx.QINVU[ctx.N + 1][1] = ctx.QINVU[ctx.N][1];
    ctx.QINVU[ctx.N + 1][2] = ctx.QINVU[ctx.N][2];

    for (int i = ctx.N + 2; i <= ctx.N + ctx.NW; ++i) {
        PsinResult res = psilin(ctx, i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i], false, false);
        ctx.QINVU[i][1] = res.qt1;
        ctx.QINVU[i][2] = res.qt2;
    }
}

void qdcalc(XFoilState &ctx) {
    std::cout << "Calculating source influence matrix ..." << std::endl;

    if (!ctx.LADIJ) {
        for (int j = 1; j <= ctx.N; ++j) {
            std::vector<double> col(ctx.N + 2, 0.0);
            for (int i = 1; i <= ctx.N + 1; ++i) {
                col[i] = ctx.BIJ[i][j];
            }
            baksub(static_cast<int>(ctx.AIJ.size()) - 1, ctx.N + 1, ctx.AIJ, ctx.AIJPIV, col);
            for (int i = 1; i <= ctx.N + 1; ++i) {
                ctx.BIJ[i][j] = col[i];
            }
            for (int i = 1; i <= ctx.N; ++i) {
                ctx.DIJ[i][j] = ctx.BIJ[i][j];
            }
        }
        ctx.LADIJ = true;
    }

    for (int i = 1; i <= ctx.N; ++i) {
        auto res = pswlin(ctx, i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i]);
        (void)res;
        for (int j = ctx.N + 1; j <= ctx.N + ctx.NW; ++j) {
            ctx.BIJ[i][j] = -ctx.DZDM[j];
        }
    }

    for (int j = ctx.N + 1; j <= ctx.N + ctx.NW; ++j) {
        ctx.BIJ[ctx.N + 1][j] = 0.0;
    }

    if (ctx.SHARP) {
        for (int j = ctx.N + 1; j <= ctx.N + ctx.NW; ++j) {
            ctx.BIJ[ctx.N][j] = 0.0;
        }
    }

    for (int j = ctx.N + 1; j <= ctx.N + ctx.NW; ++j) {
        std::vector<double> col(ctx.N + 2, 0.0);
        for (int i = 1; i <= ctx.N + 1; ++i) {
            col[i] = ctx.BIJ[i][j];
        }
        baksub(static_cast<int>(ctx.AIJ.size()) - 1, ctx.N + 1, ctx.AIJ, ctx.AIJPIV, col);
        for (int i = 1; i <= ctx.N + 1; ++i) {
            ctx.BIJ[i][j] = col[i];
        }
    }

    for (int i = 1; i <= ctx.N; ++i) {
        for (int j = ctx.N + 1; j <= ctx.N + ctx.NW; ++j) {
            ctx.DIJ[i][j] = ctx.BIJ[i][j];
        }
    }

    for (int i = ctx.N + 1; i <= ctx.N + ctx.NW; ++i) {
        const int iw = i - ctx.N;
        PsinResult res = psilin(ctx, i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i], false, true);
        (void)res;
        for (int j = 1; j <= ctx.N; ++j) {
            ctx.CIJ[iw][j] = ctx.DQDG[j];
        }
        for (int j = 1; j <= ctx.N; ++j) {
            ctx.DIJ[i][j] = ctx.DQDM[j];
        }

        auto res2 = pswlin(ctx, i, ctx.X[i], ctx.Y[i], ctx.NX[i], ctx.NY[i]);
        (void)res2;
        for (int j = ctx.N + 1; j <= ctx.N + ctx.NW; ++j) {
            ctx.DIJ[i][j] = ctx.DQDM[j];
        }
    }

    for (int i = ctx.N + 1; i <= ctx.N + ctx.NW; ++i) {
        const int iw = i - ctx.N;
        for (int j = 1; j <= ctx.N; ++j) {
            double summ = 0.0;
            for (int k = 1; k <= ctx.N; ++k) {
                summ = summ + ctx.CIJ[iw][k] * ctx.DIJ[k][j];
            }
            ctx.DIJ[i][j] = ctx.DIJ[i][j] + summ;
        }

        for (int j = ctx.N + 1; j <= ctx.N + ctx.NW; ++j) {
            double summ = 0.0;
            for (int k = 1; k <= ctx.N; ++k) {
                summ = summ + ctx.CIJ[iw][k] * ctx.BIJ[k][j];
            }
            ctx.DIJ[i][j] = ctx.DIJ[i][j] + summ;
        }
    }

    for (int j = 1; j <= ctx.N + ctx.NW; ++j) {
        ctx.DIJ[ctx.N + 1][j] = ctx.DIJ[ctx.N][j];
    }

    ctx.LWDIJ = true;
}

void xywake(XFoilState &ctx) {
    std::cout << "Calculating wake trajectory ..." << std::endl;

    ctx.NW = ctx.N / 12 + 10 * static_cast<int>(ctx.WAKLEN);
    if (ctx.NW > static_cast<int>(ctx.WGAP.size()) - 1) {
        std::cout << "\nArray size (IWX) too small.  Last wake point index reduced." << std::endl;
        ctx.NW = static_cast<int>(ctx.WGAP.size()) - 1;
    }

    const double ds1 = 0.5 * (ctx.S[2] - ctx.S[1] + ctx.S[ctx.N] - ctx.S[ctx.N - 1]);
    std::vector<double> s_tmp(ctx.NW + 1, 0.0);
    setexp(s_tmp, ds1, ctx.WAKLEN * ctx.CHORD, ctx.NW);
    for (int i = 1; i <= ctx.NW; ++i) {
        ctx.SNEW[ctx.N + i] = s_tmp[i];
    }

    ctx.XTE = 0.5 * (ctx.X[1] + ctx.X[ctx.N]);
    ctx.YTE = 0.5 * (ctx.Y[1] + ctx.Y[ctx.N]);

    int i = ctx.N + 1;
    double sx = 0.5 * (ctx.YP[ctx.N] - ctx.YP[1]);
    double sy = 0.5 * (ctx.XP[1] - ctx.XP[ctx.N]);
    double smod = std::sqrt(sx * sx + sy * sy);
    ctx.NX[i] = sx / smod;
    ctx.NY[i] = sy / smod;
    ctx.X[i] = ctx.XTE - 0.0001 * ctx.NY[i];
    ctx.Y[i] = ctx.YTE + 0.0001 * ctx.NX[i];
    ctx.S[i] = ctx.S[ctx.N];

    PsinResult resx = psilin(ctx, i, ctx.X[i], ctx.Y[i], 1.0, 0.0, false, false);
    PsinResult resy = psilin(ctx, i, ctx.X[i], ctx.Y[i], 0.0, 1.0, false, false);

    ctx.NX[i + 1] = -resx.psi_ni / std::sqrt(resx.psi_ni * resx.psi_ni + resy.psi_ni * resy.psi_ni);
    ctx.NY[i + 1] = -resy.psi_ni / std::sqrt(resx.psi_ni * resx.psi_ni + resy.psi_ni * resy.psi_ni);

    ctx.APANEL[i] = std::atan2(resy.psi_ni, resx.psi_ni);

    for (i = ctx.N + 2; i <= ctx.N + ctx.NW; ++i) {
        const double ds = ctx.SNEW[i] - ctx.SNEW[i - 1];
        ctx.X[i] = ctx.X[i - 1] - ds * ctx.NY[i];
        ctx.Y[i] = ctx.Y[i - 1] + ds * ctx.NX[i];
        ctx.S[i] = ctx.S[i - 1] + ds;

        if (i == ctx.N + ctx.NW) {
            continue;
        }

        PsinResult resx2 = psilin(ctx, i, ctx.X[i], ctx.Y[i], 1.0, 0.0, false, false);
        PsinResult resy2 = psilin(ctx, i, ctx.X[i], ctx.Y[i], 0.0, 1.0, false, false);

        ctx.NX[i + 1] = -resx2.psi_ni / std::sqrt(resx2.psi_ni * resx2.psi_ni + resy2.psi_ni * resy2.psi_ni);
        ctx.NY[i + 1] = -resy2.psi_ni / std::sqrt(resx2.psi_ni * resx2.psi_ni + resy2.psi_ni * resy2.psi_ni);

        ctx.APANEL[i] = std::atan2(resy2.psi_ni, resx2.psi_ni);
    }

    ctx.LWAKE = true;
    ctx.AWAKE = ctx.ALFA;
    ctx.LWDIJ = false;
}

void stfind(XFoilState &ctx) {
    int i = 0;
    for (i = 1; i <= ctx.N - 1; ++i) {
        if (ctx.GAM[i] >= 0.0 && ctx.GAM[i + 1] < 0.0) {
            break;
        }
    }
    if (i > ctx.N - 1) {
        std::cout << "STFIND: Stagnation point not found. Continuing ..." << std::endl;
        i = ctx.N / 2;
    }

    ctx.IST = i;
    const double dgam = ctx.GAM[i + 1] - ctx.GAM[i];
    const double ds = ctx.S[i + 1] - ctx.S[i];

    if (ctx.GAM[i] < -ctx.GAM[i + 1]) {
        ctx.SST = ctx.S[i] - ds * (ctx.GAM[i] / dgam);
    } else {
        ctx.SST = ctx.S[i + 1] - ds * (ctx.GAM[i + 1] / dgam);
    }

    if (ctx.SST <= ctx.S[i]) {
        ctx.SST = ctx.S[i] + 1.0e-7;
    }
    if (ctx.SST >= ctx.S[i + 1]) {
        ctx.SST = ctx.S[i + 1] - 1.0e-7;
    }

    ctx.SST_GO = (ctx.SST - ctx.S[i + 1]) / dgam;
    ctx.SST_GP = (ctx.S[i] - ctx.SST) / dgam;
}

void iblpan(XFoilState &ctx) {
    int is_ = 1;
    int ibl = 1;
    for (int i = ctx.IST; i >= 1; --i) {
        ibl += 1;
        ctx.IPAN[ibl][is_] = i;
        ctx.VTI[ibl][is_] = 1.0;
    }

    ctx.IBLTE[is_] = ibl;
    ctx.NBL[is_] = ibl;

    is_ = 2;
    ibl = 1;
    for (int i = ctx.IST + 1; i <= ctx.N; ++i) {
        ibl += 1;
        ctx.IPAN[ibl][is_] = i;
        ctx.VTI[ibl][is_] = -1.0;
    }

    ctx.IBLTE[is_] = ibl;

    for (int iw = 1; iw <= ctx.NW; ++iw) {
        const int i = ctx.N + iw;
        ibl = ctx.IBLTE[is_] + iw;
        ctx.IPAN[ibl][is_] = i;
        ctx.VTI[ibl][is_] = -1.0;
    }

    ctx.NBL[is_] = ctx.IBLTE[is_] + ctx.NW;

    for (int iw = 1; iw <= ctx.NW; ++iw) {
        ctx.IPAN[ctx.IBLTE[1] + iw][1] = ctx.IPAN[ctx.IBLTE[2] + iw][2];
        ctx.VTI[ctx.IBLTE[1] + iw][1] = 1.0;
    }

    const int iblmax = std::max(ctx.IBLTE[1], ctx.IBLTE[2]) + ctx.NW;
    if (iblmax > static_cast<int>(ctx.XSSI.size()) - 1) {
        throw std::runtime_error(" ***  BL array overflow.  Increase IVX to at least " + std::to_string(iblmax));
    }

    ctx.LIPAN = true;
}

void xicalc(XFoilState &ctx) {
    const double xfeps = 1.0e-7;
    const double xeps = xfeps * (ctx.S[ctx.N] - ctx.S[1]);

    int is_ = 1;
    ctx.XSSI[1][is_] = 0.0;
    for (int ibl = 2; ibl <= ctx.IBLTE[is_]; ++ibl) {
        const int i = ctx.IPAN[ibl][is_];
        ctx.XSSI[ibl][is_] = std::max(ctx.SST - ctx.S[i], xeps);
    }

    is_ = 2;
    ctx.XSSI[1][is_] = 0.0;
    for (int ibl = 2; ibl <= ctx.IBLTE[is_]; ++ibl) {
        const int i = ctx.IPAN[ibl][is_];
        ctx.XSSI[ibl][is_] = std::max(ctx.S[i] - ctx.SST, xeps);
    }

    const int is1 = 1;
    const int is2 = 2;
    int ibl1 = ctx.IBLTE[is1] + 1;
    ctx.XSSI[ibl1][is1] = ctx.XSSI[ibl1 - 1][is1];

    int ibl2 = ctx.IBLTE[is2] + 1;
    ctx.XSSI[ibl2][is2] = ctx.XSSI[ibl2 - 1][is2];

    for (int ibl = ctx.IBLTE[is_] + 2; ibl <= ctx.NBL[is_]; ++ibl) {
        const int i = ctx.IPAN[ibl][is_];
        const double dxssi = std::sqrt((ctx.X[i] - ctx.X[i - 1]) * (ctx.X[i] - ctx.X[i - 1])
                                       + (ctx.Y[i] - ctx.Y[i - 1]) * (ctx.Y[i] - ctx.Y[i - 1]));

        ibl1 = ctx.IBLTE[is1] + ibl - ctx.IBLTE[is_];
        ibl2 = ctx.IBLTE[is2] + ibl - ctx.IBLTE[is_];
        ctx.XSSI[ibl1][is1] = ctx.XSSI[ibl1 - 1][is1] + dxssi;
        ctx.XSSI[ibl2][is2] = ctx.XSSI[ibl2 - 1][is2] + dxssi;
    }

    const double telrat = 2.50;

    const double crosp = (ctx.XP[1] * ctx.YP[ctx.N] - ctx.YP[1] * ctx.XP[ctx.N])
                         / std::sqrt((ctx.XP[1] * ctx.XP[1] + ctx.YP[1] * ctx.YP[1])
                                     * (ctx.XP[ctx.N] * ctx.XP[ctx.N] + ctx.YP[ctx.N] * ctx.YP[ctx.N]));
    double dwdxte = crosp / std::sqrt(1.0 - crosp * crosp);

    dwdxte = std::max(dwdxte, -3.0 / telrat);
    dwdxte = std::min(dwdxte, 3.0 / telrat);

    const double aa = 3.0 + telrat * dwdxte;
    const double bb = -2.0 - telrat * dwdxte;

    if (ctx.SHARP) {
        for (int iw = 1; iw <= ctx.NW; ++iw) {
            ctx.WGAP[iw] = 0.0;
        }
    } else {
        is_ = 2;
        for (int iw = 1; iw <= ctx.NW; ++iw) {
            const int ibl = ctx.IBLTE[is_] + iw;
            const double zn = 1.0 - (ctx.XSSI[ibl][is_] - ctx.XSSI[ctx.IBLTE[is_]][is_]) / (telrat * ctx.ANTE);
            ctx.WGAP[iw] = 0.0;
            if (zn >= 0.0) {
                ctx.WGAP[iw] = ctx.ANTE * (aa + bb * zn) * zn * zn;
            }
        }
    }
}

void uicalc(XFoilState &ctx) {
    for (int is_ = 1; is_ <= 2; ++is_) {
        ctx.UINV[1][is_] = 0.0;
        ctx.UINV_A[1][is_] = 0.0;
        for (int ibl = 2; ibl <= ctx.NBL[is_]; ++ibl) {
            const int i = ctx.IPAN[ibl][is_];
            ctx.UINV[ibl][is_] = ctx.VTI[ibl][is_] * ctx.QINV[i];
            ctx.UINV_A[ibl][is_] = ctx.VTI[ibl][is_] * ctx.QINV_A[i];
        }
    }
}

void qvfue(XFoilState &ctx) {
    for (int is_ = 1; is_ <= 2; ++is_) {
        for (int ibl = 2; ibl <= ctx.NBL[is_]; ++ibl) {
            const int i = ctx.IPAN[ibl][is_];
            ctx.QVIS[i] = ctx.VTI[ibl][is_] * ctx.UEDG[ibl][is_];
        }
    }
}

void qiset(XFoilState &ctx) {
    ctx.COSA = std::cos(ctx.ALFA);
    ctx.SINA = std::sin(ctx.ALFA);

    for (int i = 1; i <= ctx.N + ctx.NW; ++i) {
        ctx.QINV[i] = ctx.COSA * ctx.QINVU[i][1] + ctx.SINA * ctx.QINVU[i][2];
        ctx.QINV_A[i] = -ctx.SINA * ctx.QINVU[i][1] + ctx.COSA * ctx.QINVU[i][2];
    }
}

void gamqv(XFoilState &ctx) {
    for (int i = 1; i <= ctx.N; ++i) {
        ctx.GAM[i] = ctx.QVIS[i];
        ctx.GAM_A[i] = ctx.QINV_A[i];
    }
}

void stmove(XFoilState &ctx) {
    const int istold = ctx.IST;
    stfind(ctx);

    if (istold == ctx.IST) {
        xicalc(ctx);
    } else {
        iblpan(ctx);
        uicalc(ctx);
        xicalc(ctx);
        iblsys(ctx);

        if (ctx.IST > istold) {
            const int idif = ctx.IST - istold;

            ctx.ITRAN[1] = ctx.ITRAN[1] + idif;
            ctx.ITRAN[2] = ctx.ITRAN[2] - idif;

            for (int ibl = ctx.NBL[1]; ibl >= idif + 1; --ibl) {
                ctx.CTAU[ibl][1] = ctx.CTAU[ibl - idif][1];
                ctx.THET[ibl][1] = ctx.THET[ibl - idif][1];
                ctx.DSTR[ibl][1] = ctx.DSTR[ibl - idif][1];
                ctx.UEDG[ibl][1] = ctx.UEDG[ibl - idif][1];
            }

            const double dudx = ctx.UEDG[idif + 2][1] / ctx.XSSI[idif + 2][1];
            for (int ibl = idif + 1; ibl >= 1; --ibl) {
                ctx.CTAU[ibl][1] = ctx.CTAU[idif + 2][1];
                ctx.THET[ibl][1] = ctx.THET[idif + 2][1];
                ctx.DSTR[ibl][1] = ctx.DSTR[idif + 2][1];
                ctx.UEDG[ibl][1] = dudx * ctx.XSSI[ibl][1];
            }

            for (int ibl = 2; ibl <= ctx.NBL[2]; ++ibl) {
                ctx.CTAU[ibl][2] = ctx.CTAU[ibl + idif][2];
                ctx.THET[ibl][2] = ctx.THET[ibl + idif][2];
                ctx.DSTR[ibl][2] = ctx.DSTR[ibl + idif][2];
                ctx.UEDG[ibl][2] = ctx.UEDG[ibl + idif][2];
            }

        } else {
            const int idif = istold - ctx.IST;

            ctx.ITRAN[1] = ctx.ITRAN[1] - idif;
            ctx.ITRAN[2] = ctx.ITRAN[2] + idif;

            for (int ibl = ctx.NBL[2]; ibl >= idif + 1; --ibl) {
                ctx.CTAU[ibl][2] = ctx.CTAU[ibl - idif][2];
                ctx.THET[ibl][2] = ctx.THET[ibl - idif][2];
                ctx.DSTR[ibl][2] = ctx.DSTR[ibl - idif][2];
                ctx.UEDG[ibl][2] = ctx.UEDG[ibl - idif][2];
            }

            const double dudx = ctx.UEDG[idif + 2][2] / ctx.XSSI[idif + 2][2];
            for (int ibl = idif + 1; ibl >= 1; --ibl) {
                ctx.CTAU[ibl][2] = ctx.CTAU[idif + 2][2];
                ctx.THET[ibl][2] = ctx.THET[idif + 2][2];
                ctx.DSTR[ibl][2] = ctx.DSTR[idif + 2][2];
                ctx.UEDG[ibl][2] = dudx * ctx.XSSI[ibl][2];
            }

            for (int ibl = 2; ibl <= ctx.NBL[1]; ++ibl) {
                ctx.CTAU[ibl][1] = ctx.CTAU[ibl + idif][1];
                ctx.THET[ibl][1] = ctx.THET[ibl + idif][1];
                ctx.DSTR[ibl][1] = ctx.DSTR[ibl + idif][1];
                ctx.UEDG[ibl][1] = ctx.UEDG[ibl + idif][1];
            }
        }

        const double ueps = 1.0e-7;
        for (int is_ = 1; is_ <= 2; ++is_) {
            for (int ibl = 2; ibl <= ctx.NBL[is_]; ++ibl) {
                const int i = ctx.IPAN[ibl][is_];
                if (ctx.UEDG[ibl][is_] <= ueps) {
                    ctx.UEDG[ibl][is_] = ueps;
                    ctx.QVIS[i] = ctx.VTI[ibl][is_] * ueps;
                    ctx.GAM[i] = ctx.VTI[ibl][is_] * ueps;
                }
            }
        }
    }

    for (int is_ = 1; is_ <= 2; ++is_) {
        for (int ibl = 2; ibl <= ctx.NBL[is_]; ++ibl) {
            ctx.MASS[ibl][is_] = ctx.DSTR[ibl][is_] * ctx.UEDG[ibl][is_];
        }
    }
}

void ueset(XFoilState &ctx) {
    for (int is_ = 1; is_ <= 2; ++is_) {
        for (int ibl = 2; ibl <= ctx.NBL[is_]; ++ibl) {
            const int i = ctx.IPAN[ibl][is_];

            double dui = 0.0;
            for (int js = 1; js <= 2; ++js) {
                for (int jbl = 2; jbl <= ctx.NBL[js]; ++jbl) {
                    const int j = ctx.IPAN[jbl][js];
                    const double ue_m = -ctx.VTI[ibl][is_] * ctx.VTI[jbl][js] * ctx.DIJ[i][j];
                    dui = dui + ue_m * ctx.MASS[jbl][js];
                }
            }

            ctx.UEDG[ibl][is_] = ctx.UINV[ibl][is_] + dui;
        }
    }
}

void dsset(XFoilState &ctx) {
    for (int is_ = 1; is_ <= 2; ++is_) {
        for (int ibl = 2; ibl <= ctx.NBL[is_]; ++ibl) {
            ctx.DSTR[ibl][is_] = ctx.MASS[ibl][is_] / ctx.UEDG[ibl][is_];
        }
    }
}
