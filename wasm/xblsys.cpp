#include "xblsys.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <vector>

#include "xbl.h"
#include "xblcom.h"

static std::vector<std::vector<double>> create2d(int rows, int cols, double fill = 0.0) {
    return std::vector<std::vector<double>>(static_cast<size_t>(rows + 1), std::vector<double>(static_cast<size_t>(cols + 1), fill));
}

static std::tuple<double, double, double, double> dampl(double hk, double th, double rt);
static std::tuple<double, double, double, double> dampl2(double hk, double th, double rt);
static std::tuple<double, double, double> hkin_impl(double h, double msq);
static std::tuple<double, double, double> dil_impl(double hk, double rt);
static std::tuple<double, double, double> dilw_impl(double hk, double rt);
static std::tuple<double, double, double, double> hsl(double hk, double rt, double msq);
static std::tuple<double, double, double, double> cfl(double hk, double rt, double msq);
static std::tuple<double, double, double, double, double> dit(double hs, double us, double cf, double st);
static std::tuple<double, double, double, double> hst(double hk, double rt, double msq);
static std::tuple<double, double, double, double> cft(double hk, double rt, double msq, double cffac);
static std::tuple<double, double, double> hct(double hk, double msq);

static void trchek2(XBlState &bl);
static void trdif(XBlState &bl);
static void bldif(XBlState &bl, int ityp);

static std::tuple<double, double, double, double, double, double, double, double, double>
axset(double hk1, double t1, double rt1, double a1, double hk2, double t2, double rt2, double a2, double acrit, int idampv) {
    double ax1 = 0.0;
    double ax1_hk1 = 0.0;
    double ax1_t1 = 0.0;
    double ax1_rt1 = 0.0;
    double ax2 = 0.0;
    double ax2_hk2 = 0.0;
    double ax2_t2 = 0.0;
    double ax2_rt2 = 0.0;

    if (idampv == 0) {
        auto out1 = dampl(hk1, t1, rt1);
        ax1 = std::get<0>(out1);
        ax1_hk1 = std::get<1>(out1);
        ax1_t1 = std::get<2>(out1);
        ax1_rt1 = std::get<3>(out1);

        auto out2 = dampl(hk2, t2, rt2);
        ax2 = std::get<0>(out2);
        ax2_hk2 = std::get<1>(out2);
        ax2_t2 = std::get<2>(out2);
        ax2_rt2 = std::get<3>(out2);
    } else {
        auto out1 = dampl2(hk1, t1, rt1);
        ax1 = std::get<0>(out1);
        ax1_hk1 = std::get<1>(out1);
        ax1_t1 = std::get<2>(out1);
        ax1_rt1 = std::get<3>(out1);

        auto out2 = dampl2(hk2, t2, rt2);
        ax2 = std::get<0>(out2);
        ax2_hk2 = std::get<1>(out2);
        ax2_t2 = std::get<2>(out2);
        ax2_rt2 = std::get<3>(out2);
    }

    const double axsq = 0.5 * (ax1 * ax1 + ax2 * ax2);
    double axa = 0.0;
    double axa_ax1 = 0.0;
    double axa_ax2 = 0.0;
    if (axsq <= 0.0) {
        axa = 0.0;
        axa_ax1 = 0.0;
        axa_ax2 = 0.0;
    } else {
        axa = std::sqrt(axsq);
        axa_ax1 = 0.5 * ax1 / axa;
        axa_ax2 = 0.5 * ax2 / axa;
    }

    const double arg = std::min(20.0 * (acrit - 0.5 * (a1 + a2)), 20.0);
    double exn = 0.0;
    double exn_a1 = 0.0;
    double exn_a2 = 0.0;
    if (arg <= 0.0) {
        exn = 1.0;
        exn_a1 = 0.0;
        exn_a2 = 0.0;
    } else {
        exn = std::exp(-arg);
        exn_a1 = 20.0 * 0.5 * exn;
        exn_a2 = 20.0 * 0.5 * exn;
    }

    const double dax = exn * 0.002 / (t1 + t2);
    const double dax_a1 = exn_a1 * 0.002 / (t1 + t2);
    const double dax_a2 = exn_a2 * 0.002 / (t1 + t2);
    const double dax_t1 = -dax / (t1 + t2);
    const double dax_t2 = -dax / (t1 + t2);

    const double ax = axa + dax;

    const double ax_hk1 = axa_ax1 * ax1_hk1;
    const double ax_t1 = axa_ax1 * ax1_t1 + dax_t1;
    const double ax_rt1 = axa_ax1 * ax1_rt1;
    const double ax_a1 = dax_a1;

    const double ax_hk2 = axa_ax2 * ax2_hk2;
    const double ax_t2 = axa_ax2 * ax2_t2 + dax_t2;
    const double ax_rt2 = axa_ax2 * ax2_rt2;
    const double ax_a2 = dax_a2;

    return {ax, ax_hk1, ax_t1, ax_rt1, ax_a1, ax_hk2, ax_t2, ax_rt2, ax_a2};
}

void trchek(XBlState &bl) {
    // wrapper for trchek2
    trchek2(bl);
}

void trchek2(XBlState &bl) {
    const double daeps = 5.0e-5;

    for (int icom = 1; icom <= NCOM; ++icom) {
        bl.C2SAV[icom] = bl.COM2[icom];
    }

    auto ax_out = axset(bl.HK1, bl.T1, bl.RT1, bl.AMPL1, bl.HK2, bl.T2, bl.RT2, bl.AMPL2, bl.AMCRIT, bl.IDAMPV);
    double ax = std::get<0>(ax_out);
    bl.AMPL2 = bl.AMPL1 + ax * (bl.X2 - bl.X1);

    bool converged = false;
    double last_da2 = 0.0;
    double last_ax = ax;
    double amplt = 0.0;
    double amplt_a2 = 0.0;
    double sfa = 0.0;
    double sfa_a1 = 0.0;
    double sfa_a2 = 0.0;
    double wf1 = 0.0;
    double wf1_a1 = 0.0;
    double wf1_a2 = 0.0;
    double wf1_x1 = 0.0;
    double wf1_x2 = 0.0;
    double wf1_xf = 0.0;
    double wf2 = 0.0;
    double wf2_a1 = 0.0;
    double wf2_a2 = 0.0;
    double wf2_x1 = 0.0;
    double wf2_x2 = 0.0;
    double wf2_xf = 0.0;
    double xt = 0.0;
    double xt_a2 = 0.0;
    double tt_a2 = 0.0;
    double dt_a2 = 0.0;
    double ut_a2 = 0.0;
    double hkt_tt = 0.0;
    double hkt_dt = 0.0;
    double hkt_ut = 0.0;
    double hkt_ms = 0.0;
    double rtt_tt = 0.0;
    double rtt_ut = 0.0;
    double rtt_ms = 0.0;
    double rtt_re = 0.0;
    for (int iter = 0; iter < 30; ++iter) {
        amplt = 0.0;
        amplt_a2 = 0.0;
        sfa = 0.0;
        sfa_a1 = 0.0;
        sfa_a2 = 0.0;
        if (bl.AMPL2 <= bl.AMCRIT) {
            amplt = bl.AMPL2;
            amplt_a2 = 1.0;
            sfa = 1.0;
            sfa_a1 = 0.0;
            sfa_a2 = 0.0;
        } else {
            amplt = bl.AMCRIT;
            amplt_a2 = 0.0;
            sfa = (amplt - bl.AMPL1) / (bl.AMPL2 - bl.AMPL1);
            sfa_a1 = (sfa - 1.0) / (bl.AMPL2 - bl.AMPL1);
            sfa_a2 = (-sfa) / (bl.AMPL2 - bl.AMPL1);
        }

        double sfx = 0.0;
        double sfx_x1 = 0.0;
        double sfx_x2 = 0.0;
        double sfx_xf = 0.0;
        if (bl.XIFORC < bl.X2) {
            sfx = (bl.XIFORC - bl.X1) / (bl.X2 - bl.X1);
            sfx_x1 = (sfx - 1.0) / (bl.X2 - bl.X1);
            sfx_x2 = (-sfx) / (bl.X2 - bl.X1);
            sfx_xf = 1.0 / (bl.X2 - bl.X1);
        } else {
            sfx = 1.0;
            sfx_x1 = 0.0;
            sfx_x2 = 0.0;
            sfx_xf = 0.0;
        }

        wf2 = 0.0;
        wf2_a1 = 0.0;
        wf2_a2 = 0.0;
        wf2_x1 = 0.0;
        wf2_x2 = 0.0;
        wf2_xf = 0.0;
        if (sfa < sfx) {
            wf2 = sfa;
            wf2_a1 = sfa_a1;
            wf2_a2 = sfa_a2;
            wf2_x1 = 0.0;
            wf2_x2 = 0.0;
            wf2_xf = 0.0;
        } else {
            wf2 = sfx;
            wf2_a1 = 0.0;
            wf2_a2 = 0.0;
            wf2_x1 = sfx_x1;
            wf2_x2 = sfx_x2;
            wf2_xf = sfx_xf;
        }

        wf1 = 1.0 - wf2;
        wf1_a1 = -wf2_a1;
        wf1_a2 = -wf2_a2;
        wf1_x1 = -wf2_x1;
        wf1_x2 = -wf2_x2;
        wf1_xf = -wf2_xf;

        xt = bl.X1 * wf1 + bl.X2 * wf2;
        const double tt = bl.T1 * wf1 + bl.T2 * wf2;
        const double dt = bl.D1 * wf1 + bl.D2 * wf2;
        const double ut = bl.U1 * wf1 + bl.U2 * wf2;

        xt_a2 = bl.X1 * wf1_a2 + bl.X2 * wf2_a2;
        tt_a2 = bl.T1 * wf1_a2 + bl.T2 * wf2_a2;
        dt_a2 = bl.D1 * wf1_a2 + bl.D2 * wf2_a2;
        ut_a2 = bl.U1 * wf1_a2 + bl.U2 * wf2_a2;

        bl.X2 = xt;
        bl.T2 = tt;
        bl.D2 = dt;
        bl.U2 = ut;

        blkin(bl);

        const double hkt = bl.HK2;
        hkt_tt = bl.HK2_T2;
        hkt_dt = bl.HK2_D2;
        hkt_ut = bl.HK2_U2;
        hkt_ms = bl.HK2_MS;

        const double rtt = bl.RT2;
        rtt_tt = bl.RT2_T2;
        rtt_ut = bl.RT2_U2;
        rtt_ms = bl.RT2_MS;
        rtt_re = bl.RT2_RE;

        const double amsave = bl.AMPL2;
        for (int icom = 1; icom <= NCOM; ++icom) {
            bl.COM2[icom] = bl.C2SAV[icom];
        }
        sync_com_to_vars(bl, 2);
        bl.AMPL2 = amsave;

        ax_out = axset(bl.HK1, bl.T1, bl.RT1, bl.AMPL1, hkt, tt, rtt, amplt, bl.AMCRIT, bl.IDAMPV);
        ax = std::get<0>(ax_out);
        last_ax = ax;

        if (ax <= 0.0) {
            converged = true;
            break;
        }

        const double ax_a2 = (std::get<5>(ax_out) * hkt_tt + std::get<6>(ax_out) + std::get<7>(ax_out) * rtt_tt) * tt_a2
                             + (std::get<5>(ax_out) * hkt_dt) * dt_a2
                             + (std::get<5>(ax_out) * hkt_ut + std::get<7>(ax_out) * rtt_ut) * ut_a2
                             + std::get<8>(ax_out) * amplt_a2;

        const double res = bl.AMPL2 - bl.AMPL1 - ax * (bl.X2 - bl.X1);
        const double res_a2 = 1.0 - ax_a2 * (bl.X2 - bl.X1);
        const double da2 = -res / res_a2;
        last_da2 = da2;

        double rlx = 1.0;
        const double dxt = xt_a2 * da2;
        if (rlx * std::abs(dxt / (bl.X2 - bl.X1)) > 0.05) {
            rlx = 0.05 * std::abs((bl.X2 - bl.X1) / dxt);
        }
        if (rlx * std::abs(da2) > 1.0) {
            rlx = 1.0 * std::abs(1.0 / da2);
        }

        if (std::abs(da2) < daeps) {
            converged = true;
            break;
        }

        if ((bl.AMPL2 > bl.AMCRIT && bl.AMPL2 + rlx * da2 < bl.AMCRIT)
            || (bl.AMPL2 < bl.AMCRIT && bl.AMPL2 + rlx * da2 > bl.AMCRIT)) {
            bl.AMPL2 = bl.AMCRIT;
        } else {
            bl.AMPL2 = bl.AMPL2 + rlx * da2;
        }
    }

    if (!converged) {
        std::cout << "TRCHEK2: N2 convergence failed." << std::endl;
        std::cout << " x:" << std::setw(9) << std::fixed << std::setprecision(5) << bl.X1 << std::setw(9) << bl.XT << std::setw(9)
                  << bl.X2 << "  N:" << std::setw(7) << std::setprecision(3) << bl.AMPL1 << std::setw(7) << bl.AMCRIT
                  << std::setw(7) << bl.AMPL2 << "  Nx:" << std::setw(8) << std::setprecision(3) << last_ax
                  << "   dN:" << std::setw(10) << std::scientific << last_da2 << std::defaultfloat << std::endl;
    }

    bl.XT = xt;
    bl.XT_A2 = xt_a2;

    bl.TRFREE = bl.AMPL2 >= bl.AMCRIT;
    bl.TRFORC = bl.XIFORC > bl.X1 && bl.XIFORC <= bl.X2;
    bl.TRAN = bl.TRFORC || bl.TRFREE;
    if (!bl.TRAN) {
        return;
    }

    if (bl.TRFREE && bl.TRFORC) {
        bl.TRFORC = bl.XIFORC < xt;
        bl.TRFREE = bl.XIFORC >= xt;
    }

    if (bl.TRFORC) {
        bl.XT = bl.XIFORC;
        bl.XT_A1 = 0.0;
        bl.XT_X1 = 0.0;
        bl.XT_T1 = 0.0;
        bl.XT_D1 = 0.0;
        bl.XT_U1 = 0.0;
        bl.XT_X2 = 0.0;
        bl.XT_T2 = 0.0;
        bl.XT_D2 = 0.0;
        bl.XT_U2 = 0.0;
        bl.XT_MS = 0.0;
        bl.XT_RE = 0.0;
        bl.XT_XF = 1.0;
        return;
    }

    bl.XT = bl.X1 * wf1 + bl.X2 * wf2;
    bl.XT_A2 = xt_a2;
    const double tt = bl.T1 * wf1 + bl.T2 * wf2;
    const double dt = bl.D1 * wf1 + bl.D2 * wf2;
    const double ut = bl.U1 * wf1 + bl.U2 * wf2;

    bl.XT_X1 = wf1;
    double tt_t1 = wf1;
    double dt_d1 = wf1;
    double ut_u1 = wf1;

    bl.XT_X2 = wf2;
    double tt_t2 = wf2;
    double dt_d2 = wf2;
    double ut_u2 = wf2;

    bl.XT_A1 = bl.X1 * wf1_a1 + bl.X2 * wf2_a1;
    double tt_a1 = bl.T1 * wf1_a1 + bl.T2 * wf2_a1;
    double dt_a1 = bl.D1 * wf1_a1 + bl.D2 * wf2_a1;
    double ut_a1 = bl.U1 * wf1_a1 + bl.U2 * wf2_a1;

    bl.XT_X1 = bl.X1 * wf1_x1 + bl.X2 * wf2_x1 + bl.XT_X1;
    double tt_x1 = bl.T1 * wf1_x1 + bl.T2 * wf2_x1;
    double dt_x1 = bl.D1 * wf1_x1 + bl.D2 * wf2_x1;
    double ut_x1 = bl.U1 * wf1_x1 + bl.U2 * wf2_x1;

    bl.XT_X2 = bl.X1 * wf1_x2 + bl.X2 * wf2_x2 + bl.XT_X2;
    double tt_x2 = bl.T1 * wf1_x2 + bl.T2 * wf2_x2;
    double dt_x2 = bl.D1 * wf1_x2 + bl.D2 * wf2_x2;
    double ut_x2 = bl.U1 * wf1_x2 + bl.U2 * wf2_x2;

    bl.XT_XF = bl.X1 * wf1_xf + bl.X2 * wf2_xf;
    double tt_xf = bl.T1 * wf1_xf + bl.T2 * wf2_xf;
    double dt_xf = bl.D1 * wf1_xf + bl.D2 * wf2_xf;
    double ut_xf = bl.U1 * wf1_xf + bl.U2 * wf2_xf;

    const double ax_t1 = std::get<1>(ax_out) * bl.HK1_T1 + std::get<2>(ax_out) + std::get<3>(ax_out) * bl.RT1_T1
                         + (std::get<5>(ax_out) * hkt_tt + std::get<6>(ax_out) + std::get<7>(ax_out) * rtt_tt) * tt_t1;
    const double ax_d1 = std::get<1>(ax_out) * bl.HK1_D1 + (std::get<5>(ax_out) * hkt_dt) * dt_d1;
    const double ax_u1 = std::get<1>(ax_out) * bl.HK1_U1 + std::get<3>(ax_out) * bl.RT1_U1
                         + (std::get<5>(ax_out) * hkt_ut + std::get<7>(ax_out) * rtt_ut) * ut_u1;
    const double ax_a1 = std::get<4>(ax_out)
                         + (std::get<5>(ax_out) * hkt_tt + std::get<6>(ax_out) + std::get<7>(ax_out) * rtt_tt) * tt_a1
                         + (std::get<5>(ax_out) * hkt_dt) * dt_a1
                         + (std::get<5>(ax_out) * hkt_ut + std::get<7>(ax_out) * rtt_ut) * ut_a1;
    const double ax_x1 = (std::get<5>(ax_out) * hkt_tt + std::get<6>(ax_out) + std::get<7>(ax_out) * rtt_tt) * tt_x1
                         + (std::get<5>(ax_out) * hkt_dt) * dt_x1
                         + (std::get<5>(ax_out) * hkt_ut + std::get<7>(ax_out) * rtt_ut) * ut_x1;

    const double ax_t2 = (std::get<5>(ax_out) * hkt_tt + std::get<6>(ax_out) + std::get<7>(ax_out) * rtt_tt) * tt_t2;
    const double ax_d2 = (std::get<5>(ax_out) * hkt_dt) * dt_d2;
    const double ax_u2 = (std::get<5>(ax_out) * hkt_ut + std::get<7>(ax_out) * rtt_ut) * ut_u2;
    const double ax_a2 = std::get<8>(ax_out) * amplt_a2
                         + (std::get<5>(ax_out) * hkt_tt + std::get<6>(ax_out) + std::get<7>(ax_out) * rtt_tt) * tt_a2
                         + (std::get<5>(ax_out) * hkt_dt) * dt_a2
                         + (std::get<5>(ax_out) * hkt_ut + std::get<7>(ax_out) * rtt_ut) * ut_a2;
    const double ax_x2 = (std::get<5>(ax_out) * hkt_tt + std::get<6>(ax_out) + std::get<7>(ax_out) * rtt_tt) * tt_x2
                         + (std::get<5>(ax_out) * hkt_dt) * dt_x2
                         + (std::get<5>(ax_out) * hkt_ut + std::get<7>(ax_out) * rtt_ut) * ut_x2;

    const double ax_xf = (std::get<5>(ax_out) * hkt_tt + std::get<6>(ax_out) + std::get<7>(ax_out) * rtt_tt) * tt_xf
                         + (std::get<5>(ax_out) * hkt_dt) * dt_xf
                         + (std::get<5>(ax_out) * hkt_ut + std::get<7>(ax_out) * rtt_ut) * ut_xf;

    const double ax_ms = std::get<5>(ax_out) * hkt_ms + std::get<7>(ax_out) * rtt_ms + std::get<1>(ax_out) * bl.HK1_MS
                         + std::get<3>(ax_out) * bl.RT1_MS;
    const double ax_re = std::get<7>(ax_out) * rtt_re + std::get<3>(ax_out) * bl.RT1_RE;

    const double z_ax = -(bl.X2 - bl.X1);
    const double z_a1 = z_ax * ax_a1 - 1.0;
    const double z_t1 = z_ax * ax_t1;
    const double z_d1 = z_ax * ax_d1;
    const double z_u1 = z_ax * ax_u1;
    const double z_x1 = z_ax * ax_x1 + ax;
    const double z_a2 = z_ax * ax_a2 + 1.0;
    const double z_t2 = z_ax * ax_t2;
    const double z_d2 = z_ax * ax_d2;
    const double z_u2 = z_ax * ax_u2;
    const double z_x2 = z_ax * ax_x2 - ax;
    const double z_xf = z_ax * ax_xf;
    const double z_ms = z_ax * ax_ms;
    const double z_re = z_ax * ax_re;

    bl.XT_A1 = bl.XT_A1 - (bl.XT_A2 / z_a2) * z_a1;
    bl.XT_T1 = -(bl.XT_A2 / z_a2) * z_t1;
    bl.XT_D1 = -(bl.XT_A2 / z_a2) * z_d1;
    bl.XT_U1 = -(bl.XT_A2 / z_a2) * z_u1;
    bl.XT_X1 = bl.XT_X1 - (bl.XT_A2 / z_a2) * z_x1;
    bl.XT_T2 = -(bl.XT_A2 / z_a2) * z_t2;
    bl.XT_D2 = -(bl.XT_A2 / z_a2) * z_d2;
    bl.XT_U2 = -(bl.XT_A2 / z_a2) * z_u2;
    bl.XT_X2 = bl.XT_X2 - (bl.XT_A2 / z_a2) * z_x2;
    bl.XT_MS = -(bl.XT_A2 / z_a2) * z_ms;
    bl.XT_RE = -(bl.XT_A2 / z_a2) * z_re;
    bl.XT_XF = 0.0;
}

void blsys(XBlState &bl) {
    if (bl.WAKE) {
        blvar(bl, 3);
        blmid(bl, 3);
    } else if (bl.TURB || bl.TRAN) {
        blvar(bl, 2);
        blmid(bl, 2);
    } else {
        blvar(bl, 1);
        blmid(bl, 1);
    }

    if (bl.SIMI) {
        for (int icom = 1; icom <= NCOM; ++icom) {
            bl.COM1[icom] = bl.COM2[icom];
        }
        sync_com_to_vars(bl, 1);
    }

    if (bl.TRAN) {
        trdif(bl);
    } else if (bl.SIMI) {
        bldif(bl, 0);
    } else if (!bl.TURB) {
        bldif(bl, 1);
    } else if (bl.WAKE) {
        bldif(bl, 3);
    } else {
        bldif(bl, 2);
    }

    if (bl.SIMI) {
        for (int k = 1; k <= 4; ++k) {
            for (int l = 1; l <= 5; ++l) {
                bl.VS2[k][l] = bl.VS1[k][l] + bl.VS2[k][l];
                bl.VS1[k][l] = 0.0;
            }
        }
    }

    for (int k = 1; k <= 4; ++k) {
        const double res_u1 = bl.VS1[k][4];
        const double res_u2 = bl.VS2[k][4];
        const double res_ms = bl.VSM[k];
        bl.VS1[k][4] = res_u1 * bl.U1_UEI;
        bl.VS2[k][4] = res_u2 * bl.U2_UEI;
        bl.VSM[k] = res_u1 * bl.U1_MS + res_u2 * bl.U2_MS + res_ms;
    }
}

void tesys(XBlState &bl, double cte, double tte, double dte) {
    for (int k = 1; k <= 4; ++k) {
        bl.VSREZ[k] = 0.0;
        bl.VSM[k] = 0.0;
        bl.VSR[k] = 0.0;
        bl.VSX[k] = 0.0;
        for (int l = 1; l <= 5; ++l) {
            bl.VS1[k][l] = 0.0;
            bl.VS2[k][l] = 0.0;
        }
    }

    blvar(bl, 3);

    bl.VS1[1][1] = -1.0;
    bl.VS2[1][1] = 1.0;
    bl.VSREZ[1] = cte - bl.S2;

    bl.VS1[2][2] = -1.0;
    bl.VS2[2][2] = 1.0;
    bl.VSREZ[2] = tte - bl.T2;

    bl.VS1[3][3] = -1.0;
    bl.VS2[3][3] = 1.0;
    bl.VSREZ[3] = dte - bl.D2 - bl.DW2;
}

void blprv(XBlState &bl, double xsi, double ami, double cti, double thi, double dsi, double dswaki, double uei) {
    bl.X2 = xsi;
    bl.AMPL2 = ami;
    bl.S2 = cti;
    bl.T2 = thi;
    bl.D2 = dsi - dswaki;
    bl.DW2 = dswaki;

    bl.U2 = uei * (1.0 - bl.TKBL) / (1.0 - bl.TKBL * std::pow(uei / bl.QINFBL, 2));
    bl.U2_UEI = (1.0 + bl.TKBL * (2.0 * bl.U2 * uei / (bl.QINFBL * bl.QINFBL) - 1.0))
               / (1.0 - bl.TKBL * std::pow(uei / bl.QINFBL, 2));
    bl.U2_MS = (bl.U2 * std::pow(uei / bl.QINFBL, 2) - uei) * bl.TKBL_MS / (1.0 - bl.TKBL * std::pow(uei / bl.QINFBL, 2));
    sync_vars_to_com(bl, 2);
}

void blkin(XBlState &bl) {
    bl.M2 = bl.U2 * bl.U2 * bl.HSTINV / (bl.GM1BL * (1.0 - 0.5 * bl.U2 * bl.U2 * bl.HSTINV));
    const double tr2 = 1.0 + 0.5 * bl.GM1BL * bl.M2;
    bl.M2_U2 = 2.0 * bl.M2 * tr2 / bl.U2;
    bl.M2_MS = bl.U2 * bl.U2 * tr2 / (bl.GM1BL * (1.0 - 0.5 * bl.U2 * bl.U2 * bl.HSTINV)) * bl.HSTINV_MS;

    bl.R2 = bl.RSTBL * std::pow(tr2, (-1.0 / bl.GM1BL));
    bl.R2_U2 = -bl.R2 / tr2 * 0.5 * bl.M2_U2;
    bl.R2_MS = -bl.R2 / tr2 * 0.5 * bl.M2_MS + bl.RSTBL_MS * std::pow(tr2, (-1.0 / bl.GM1BL));

    bl.H2 = bl.D2 / bl.T2;
    bl.H2_D2 = 1.0 / bl.T2;
    bl.H2_T2 = -bl.H2 / bl.T2;

    const double herat = 1.0 - 0.5 * bl.U2 * bl.U2 * bl.HSTINV;
    const double he_u2 = -bl.U2 * bl.HSTINV;
    const double he_ms = -0.5 * bl.U2 * bl.U2 * bl.HSTINV_MS;

    bl.V2 = std::sqrt(std::pow(herat, 3)) * (1.0 + bl.HVRAT) / (herat + bl.HVRAT) / bl.REYBL;
    const double v2_he = bl.V2 * (1.5 / herat - 1.0 / (herat + bl.HVRAT));

    bl.V2_U2 = v2_he * he_u2;
    bl.V2_MS = -bl.V2 / bl.REYBL * bl.REYBL_MS + v2_he * he_ms;
    bl.V2_RE = -bl.V2 / bl.REYBL * bl.REYBL_RE;

    auto hk_out = hkin_impl(bl.H2, bl.M2);
    bl.HK2 = std::get<0>(hk_out);
    const double hk2_h2 = std::get<1>(hk_out);
    const double hk2_m2 = std::get<2>(hk_out);

    bl.HK2_U2 = hk2_m2 * bl.M2_U2;
    bl.HK2_T2 = hk2_h2 * bl.H2_T2;
    bl.HK2_D2 = hk2_h2 * bl.H2_D2;
    bl.HK2_MS = hk2_m2 * bl.M2_MS;

    bl.RT2 = bl.R2 * bl.U2 * bl.T2 / bl.V2;
    bl.RT2_U2 = bl.RT2 * (1.0 / bl.U2 + bl.R2_U2 / bl.R2 - bl.V2_U2 / bl.V2);
    bl.RT2_T2 = bl.RT2 / bl.T2;
    bl.RT2_MS = bl.RT2 * (bl.R2_MS / bl.R2 - bl.V2_MS / bl.V2);
    bl.RT2_RE = bl.RT2 * (-bl.V2_RE / bl.V2);
    sync_vars_to_com(bl, 2);
}

void blvar(XBlState &bl, int ityp) {
    if (ityp == 3) {
        bl.HK2 = std::max(bl.HK2, 1.00005);
    }
    if (ityp != 3) {
        bl.HK2 = std::max(bl.HK2, 1.05000);
    }

    auto hc_out = hct(bl.HK2, bl.M2);
    bl.HC2 = std::get<0>(hc_out);
    const double hc2_hk2 = std::get<1>(hc_out);
    const double hc2_m2 = std::get<2>(hc_out);
    bl.HC2_U2 = hc2_hk2 * bl.HK2_U2 + hc2_m2 * bl.M2_U2;
    bl.HC2_T2 = hc2_hk2 * bl.HK2_T2;
    bl.HC2_D2 = hc2_hk2 * bl.HK2_D2;
    bl.HC2_MS = hc2_hk2 * bl.HK2_MS + hc2_m2 * bl.M2_MS;

    double hs2 = 0.0;
    double hs2_hk2 = 0.0;
    double hs2_rt2 = 0.0;
    double hs2_m2 = 0.0;
    if (ityp == 1) {
        auto hs_out = hsl(bl.HK2, bl.RT2, bl.M2);
        hs2 = std::get<0>(hs_out);
        hs2_hk2 = std::get<1>(hs_out);
        hs2_rt2 = std::get<2>(hs_out);
        hs2_m2 = std::get<3>(hs_out);
    } else {
        auto hs_out = hst(bl.HK2, bl.RT2, bl.M2);
        hs2 = std::get<0>(hs_out);
        hs2_hk2 = std::get<1>(hs_out);
        hs2_rt2 = std::get<2>(hs_out);
        hs2_m2 = std::get<3>(hs_out);
    }
    bl.HS2 = hs2;

    bl.HS2_U2 = hs2_hk2 * bl.HK2_U2 + hs2_rt2 * bl.RT2_U2 + hs2_m2 * bl.M2_U2;
    bl.HS2_T2 = hs2_hk2 * bl.HK2_T2 + hs2_rt2 * bl.RT2_T2;
    bl.HS2_D2 = hs2_hk2 * bl.HK2_D2;
    bl.HS2_MS = hs2_hk2 * bl.HK2_MS + hs2_rt2 * bl.RT2_MS + hs2_m2 * bl.M2_MS;
    bl.HS2_RE = hs2_rt2 * bl.RT2_RE;

    bl.US2 = 0.5 * bl.HS2 * (1.0 - (bl.HK2 - 1.0) / (bl.GBCON * bl.H2));
    const double us2_hs2 = 0.5 * (1.0 - (bl.HK2 - 1.0) / (bl.GBCON * bl.H2));
    const double us2_hk2 = 0.5 * bl.HS2 * (-1.0 / (bl.GBCON * bl.H2));
    const double us2_h2 = 0.5 * bl.HS2 * (bl.HK2 - 1.0) / (bl.GBCON * bl.H2 * bl.H2);

    bl.US2_U2 = us2_hs2 * bl.HS2_U2 + us2_hk2 * bl.HK2_U2;
    bl.US2_T2 = us2_hs2 * bl.HS2_T2 + us2_hk2 * bl.HK2_T2 + us2_h2 * bl.H2_T2;
    bl.US2_D2 = us2_hs2 * bl.HS2_D2 + us2_hk2 * bl.HK2_D2 + us2_h2 * bl.H2_D2;
    bl.US2_MS = us2_hs2 * bl.HS2_MS + us2_hk2 * bl.HK2_MS;
    bl.US2_RE = us2_hs2 * bl.HS2_RE;

    if (ityp <= 2 && bl.US2 > 0.95) {
        bl.US2 = 0.98;
        bl.US2_U2 = 0.0;
        bl.US2_T2 = 0.0;
        bl.US2_D2 = 0.0;
        bl.US2_MS = 0.0;
        bl.US2_RE = 0.0;
    }

    if (ityp == 3 && bl.US2 > 0.99995) {
        bl.US2 = 0.99995;
        bl.US2_U2 = 0.0;
        bl.US2_T2 = 0.0;
        bl.US2_D2 = 0.0;
        bl.US2_MS = 0.0;
        bl.US2_RE = 0.0;
    }

    double hkc = 0.0;
    double hkc_hk2 = 0.0;
    double hkc_rt2 = 0.0;
    if (ityp == 2) {
        const double gcc = bl.GCCON;
        hkc = bl.HK2 - 1.0 - gcc / bl.RT2;
        hkc_hk2 = 1.0;
        hkc_rt2 = gcc / (bl.RT2 * bl.RT2);
        if (hkc < 0.01) {
            hkc = 0.01;
            hkc_hk2 = 0.0;
            hkc_rt2 = 0.0;
        }
    } else {
        hkc = bl.HK2 - 1.0;
        hkc_hk2 = 1.0;
        hkc_rt2 = 0.0;
    }

    const double hkb = bl.HK2 - 1.0;
    const double usb = 1.0 - bl.US2;
    bl.CQ2 = std::sqrt(bl.CTCON * bl.HS2 * hkb * hkc * hkc / (usb * bl.H2 * bl.HK2 * bl.HK2));
    const double cq2_hs2 = bl.CTCON * hkb * hkc * hkc / (usb * bl.H2 * bl.HK2 * bl.HK2) * 0.5 / bl.CQ2;
    const double cq2_us2 = bl.CTCON * bl.HS2 * hkb * hkc * hkc / (usb * bl.H2 * bl.HK2 * bl.HK2) / usb * 0.5 / bl.CQ2;
    const double cq2_hk2 = bl.CTCON * bl.HS2 * hkc * hkc / (usb * bl.H2 * bl.HK2 * bl.HK2) * 0.5 / bl.CQ2
                           - bl.CTCON * bl.HS2 * hkb * hkc * hkc / (usb * bl.H2 * std::pow(bl.HK2, 3)) * 2.0 * 0.5 / bl.CQ2
                           + bl.CTCON * bl.HS2 * hkb * hkc / (usb * bl.H2 * bl.HK2 * bl.HK2) * 2.0 * 0.5 / bl.CQ2 * hkc_hk2;
    const double cq2_h2 = -bl.CTCON * bl.HS2 * hkb * hkc * hkc / (usb * bl.H2 * bl.HK2 * bl.HK2) / bl.H2 * 0.5 / bl.CQ2;
    const double cq2_rt2 = bl.CTCON * bl.HS2 * hkb * hkc / (usb * bl.H2 * bl.HK2 * bl.HK2) * 2.0 * 0.5 / bl.CQ2 * hkc_rt2;

    bl.CQ2_U2 = cq2_hs2 * bl.HS2_U2 + cq2_us2 * bl.US2_U2 + cq2_hk2 * bl.HK2_U2;
    bl.CQ2_T2 = cq2_hs2 * bl.HS2_T2 + cq2_us2 * bl.US2_T2 + cq2_hk2 * bl.HK2_T2 + cq2_h2 * bl.H2_T2;
    bl.CQ2_D2 = cq2_hs2 * bl.HS2_D2 + cq2_us2 * bl.US2_D2 + cq2_hk2 * bl.HK2_D2 + cq2_h2 * bl.H2_D2;
    bl.CQ2_MS = cq2_hs2 * bl.HS2_MS + cq2_us2 * bl.US2_MS + cq2_hk2 * bl.HK2_MS;
    bl.CQ2_RE = cq2_hs2 * bl.HS2_RE + cq2_us2 * bl.US2_RE + cq2_rt2 * bl.RT2_RE;

    bl.CQ2_U2 = bl.CQ2_U2 + cq2_rt2 * bl.RT2_U2;
    bl.CQ2_T2 = bl.CQ2_T2 + cq2_rt2 * bl.RT2_T2;
    bl.CQ2_MS = bl.CQ2_MS + cq2_rt2 * bl.RT2_MS;

    if (ityp == 3) {
        bl.CF2 = 0.0;
        bl.CF2_HK2 = 0.0;
        bl.CF2_RT2 = 0.0;
        bl.CF2_M2 = 0.0;
    } else if (ityp == 1) {
        auto cf_out = cfl(bl.HK2, bl.RT2, bl.M2);
        bl.CF2 = std::get<0>(cf_out);
        bl.CF2_HK2 = std::get<1>(cf_out);
        bl.CF2_RT2 = std::get<2>(cf_out);
        bl.CF2_M2 = std::get<3>(cf_out);
    } else {
        auto cf_out = cft(bl.HK2, bl.RT2, bl.M2, bl.CFFAC);
        bl.CF2 = std::get<0>(cf_out);
        bl.CF2_HK2 = std::get<1>(cf_out);
        bl.CF2_RT2 = std::get<2>(cf_out);
        bl.CF2_M2 = std::get<3>(cf_out);
        auto cfl_out = cfl(bl.HK2, bl.RT2, bl.M2);
        if (std::get<0>(cfl_out) > bl.CF2) {
            bl.CF2 = std::get<0>(cfl_out);
            bl.CF2_HK2 = std::get<1>(cfl_out);
            bl.CF2_RT2 = std::get<2>(cfl_out);
            bl.CF2_M2 = std::get<3>(cfl_out);
        }
    }

    bl.CF2_U2 = bl.CF2_HK2 * bl.HK2_U2 + bl.CF2_RT2 * bl.RT2_U2 + bl.CF2_M2 * bl.M2_U2;
    bl.CF2_T2 = bl.CF2_HK2 * bl.HK2_T2 + bl.CF2_RT2 * bl.RT2_T2;
    bl.CF2_D2 = bl.CF2_HK2 * bl.HK2_D2;
    bl.CF2_MS = bl.CF2_HK2 * bl.HK2_MS + bl.CF2_RT2 * bl.RT2_MS + bl.CF2_M2 * bl.M2_MS;
    bl.CF2_RE = bl.CF2_RT2 * bl.RT2_RE;

    if (ityp == 1) {
        auto di_out = dil_impl(bl.HK2, bl.RT2);
        bl.DI2 = std::get<0>(di_out);
        const double di2_hk2 = std::get<1>(di_out);
        const double di2_rt2 = std::get<2>(di_out);
        bl.DI2_U2 = di2_hk2 * bl.HK2_U2 + di2_rt2 * bl.RT2_U2;
        bl.DI2_T2 = di2_hk2 * bl.HK2_T2 + di2_rt2 * bl.RT2_T2;
        bl.DI2_D2 = di2_hk2 * bl.HK2_D2;
        bl.DI2_S2 = 0.0;
        bl.DI2_MS = di2_hk2 * bl.HK2_MS + di2_rt2 * bl.RT2_MS;
        bl.DI2_RE = di2_rt2 * bl.RT2_RE;
    } else if (ityp == 2) {
        auto cf2t_out = cft(bl.HK2, bl.RT2, bl.M2, bl.CFFAC);
        const double cf2t = std::get<0>(cf2t_out);
        const double cf2t_hk = std::get<1>(cf2t_out);
        const double cf2t_rt = std::get<2>(cf2t_out);
        const double cf2t_ms = std::get<3>(cf2t_out);

        const double cf2t_u2 = cf2t_hk * bl.HK2_U2 + cf2t_rt * bl.RT2_U2 + cf2t_ms * bl.M2_U2;
        const double cf2t_t2 = cf2t_hk * bl.HK2_T2 + cf2t_rt * bl.RT2_T2;
        const double cf2t_d2 = cf2t_hk * bl.HK2_D2;
        const double cf2t_ms2 = cf2t_hk * bl.HK2_MS + cf2t_rt * bl.RT2_MS + cf2t_ms * bl.M2_MS;
        const double cf2t_re = cf2t_rt * bl.RT2_RE;

        bl.DI2 = (0.5 * cf2t * bl.US2) * 2.0 / bl.HS2;
        const double di2_hs2 = -(0.5 * cf2t * bl.US2) * 2.0 / (bl.HS2 * bl.HS2);
        const double di2_us2 = (0.5 * cf2t) * 2.0 / bl.HS2;
        const double di2_cf2t = (0.5 * bl.US2) * 2.0 / bl.HS2;

        bl.DI2_S2 = 0.0;
        bl.DI2_U2 = di2_hs2 * bl.HS2_U2 + di2_us2 * bl.US2_U2 + di2_cf2t * cf2t_u2;
        bl.DI2_T2 = di2_hs2 * bl.HS2_T2 + di2_us2 * bl.US2_T2 + di2_cf2t * cf2t_t2;
        bl.DI2_D2 = di2_hs2 * bl.HS2_D2 + di2_us2 * bl.US2_D2 + di2_cf2t * cf2t_d2;
        bl.DI2_MS = di2_hs2 * bl.HS2_MS + di2_us2 * bl.US2_MS + di2_cf2t * cf2t_ms2;
        bl.DI2_RE = di2_hs2 * bl.HS2_RE + di2_us2 * bl.US2_RE + di2_cf2t * cf2t_re;

        const double grt = std::log(bl.RT2);
        const double hmin = 1.0 + 2.1 / grt;
        const double hm_rt2 = -(2.1 / (grt * grt)) / bl.RT2;

        const double fl = (bl.HK2 - 1.0) / (hmin - 1.0);
        const double fl_hk2 = 1.0 / (hmin - 1.0);
        const double fl_rt2 = (-fl / (hmin - 1.0)) * hm_rt2;

        const double tfl = std::tanh(fl);
        const double dfac = 0.5 + 0.5 * tfl;
        const double df_fl = 0.5 * (1.0 - tfl * tfl);

        const double df_hk2 = df_fl * fl_hk2;
        const double df_rt2 = df_fl * fl_rt2;

        bl.DI2_S2 = bl.DI2_S2 * dfac;
        bl.DI2_U2 = bl.DI2_U2 * dfac + bl.DI2 * (df_hk2 * bl.HK2_U2 + df_rt2 * bl.RT2_U2);
        bl.DI2_T2 = bl.DI2_T2 * dfac + bl.DI2 * (df_hk2 * bl.HK2_T2 + df_rt2 * bl.RT2_T2);
        bl.DI2_D2 = bl.DI2_D2 * dfac + bl.DI2 * (df_hk2 * bl.HK2_D2);
        bl.DI2_MS = bl.DI2_MS * dfac + bl.DI2 * (df_hk2 * bl.HK2_MS + df_rt2 * bl.RT2_MS);
        bl.DI2_RE = bl.DI2_RE * dfac + bl.DI2 * (df_rt2 * bl.RT2_RE);
        bl.DI2 = bl.DI2 * dfac;
    } else {
        bl.DI2 = 0.0;
        bl.DI2_S2 = 0.0;
        bl.DI2_U2 = 0.0;
        bl.DI2_T2 = 0.0;
        bl.DI2_D2 = 0.0;
        bl.DI2_MS = 0.0;
        bl.DI2_RE = 0.0;
    }

    if (ityp != 1) {
        double dd = bl.S2 * bl.S2 * (0.995 - bl.US2) * 2.0 / bl.HS2;
        const double dd_hs2 = -(bl.S2 * bl.S2) * (0.995 - bl.US2) * 2.0 / (bl.HS2 * bl.HS2);
        const double dd_us2 = -(bl.S2 * bl.S2) * 2.0 / bl.HS2;
        const double dd_s2 = bl.S2 * 2.0 * (0.995 - bl.US2) * 2.0 / bl.HS2;

        bl.DI2 += dd;
        bl.DI2_S2 = dd_s2;
        bl.DI2_U2 += dd_hs2 * bl.HS2_U2 + dd_us2 * bl.US2_U2;
        bl.DI2_T2 += dd_hs2 * bl.HS2_T2 + dd_us2 * bl.US2_T2;
        bl.DI2_D2 += dd_hs2 * bl.HS2_D2 + dd_us2 * bl.US2_D2;
        bl.DI2_MS += dd_hs2 * bl.HS2_MS + dd_us2 * bl.US2_MS;
        bl.DI2_RE += dd_hs2 * bl.HS2_RE + dd_us2 * bl.US2_RE;

        dd = 0.15 * std::pow(0.995 - bl.US2, 2) / bl.RT2 * 2.0 / bl.HS2;
        const double dd_us2b = -0.15 * (0.995 - bl.US2) * 2.0 / bl.RT2 * 2.0 / bl.HS2;
        const double dd_hs2b = -dd / bl.HS2;
        const double dd_rt2 = -dd / bl.RT2;

        bl.DI2 += dd;
        bl.DI2_U2 += dd_hs2b * bl.HS2_U2 + dd_us2b * bl.US2_U2 + dd_rt2 * bl.RT2_U2;
        bl.DI2_T2 += dd_hs2b * bl.HS2_T2 + dd_us2b * bl.US2_T2 + dd_rt2 * bl.RT2_T2;
        bl.DI2_D2 += dd_hs2b * bl.HS2_D2 + dd_us2b * bl.US2_D2;
        bl.DI2_MS += dd_hs2b * bl.HS2_MS + dd_us2b * bl.US2_MS + dd_rt2 * bl.RT2_MS;
        bl.DI2_RE += dd_hs2b * bl.HS2_RE + dd_us2b * bl.US2_RE + dd_rt2 * bl.RT2_RE;
    }

    if (ityp == 2) {
        auto di2l_out = dil_impl(bl.HK2, bl.RT2);
        if (std::get<0>(di2l_out) > bl.DI2) {
            bl.DI2 = std::get<0>(di2l_out);
            bl.DI2_S2 = 0.0;
            bl.DI2_U2 = std::get<1>(di2l_out) * bl.HK2_U2 + std::get<2>(di2l_out) * bl.RT2_U2;
            bl.DI2_T2 = std::get<1>(di2l_out) * bl.HK2_T2 + std::get<2>(di2l_out) * bl.RT2_T2;
            bl.DI2_D2 = std::get<1>(di2l_out) * bl.HK2_D2;
            bl.DI2_MS = std::get<1>(di2l_out) * bl.HK2_MS + std::get<2>(di2l_out) * bl.RT2_MS;
            bl.DI2_RE = std::get<2>(di2l_out) * bl.RT2_RE;
        }
    }

    if (ityp == 3) {
        auto di2l_out = dilw_impl(bl.HK2, bl.RT2);
        if (std::get<0>(di2l_out) > bl.DI2) {
            bl.DI2 = std::get<0>(di2l_out);
            bl.DI2_S2 = 0.0;
            bl.DI2_U2 = std::get<1>(di2l_out) * bl.HK2_U2 + std::get<2>(di2l_out) * bl.RT2_U2;
            bl.DI2_T2 = std::get<1>(di2l_out) * bl.HK2_T2 + std::get<2>(di2l_out) * bl.RT2_T2;
            bl.DI2_D2 = std::get<1>(di2l_out) * bl.HK2_D2;
            bl.DI2_MS = std::get<1>(di2l_out) * bl.HK2_MS + std::get<2>(di2l_out) * bl.RT2_MS;
            bl.DI2_RE = std::get<2>(di2l_out) * bl.RT2_RE;
        }
    }

    if (ityp == 3) {
        bl.DI2 *= 2.0;
        bl.DI2_S2 *= 2.0;
        bl.DI2_U2 *= 2.0;
        bl.DI2_T2 *= 2.0;
        bl.DI2_D2 *= 2.0;
        bl.DI2_MS *= 2.0;
        bl.DI2_RE *= 2.0;
    }

    bl.DE2 = (3.15 + 1.72 / (bl.HK2 - 1.0)) * bl.T2 + bl.D2;
    const double de2_hk2 = (-1.72 / std::pow(bl.HK2 - 1.0, 2)) * bl.T2;

    bl.DE2_U2 = de2_hk2 * bl.HK2_U2;
    bl.DE2_T2 = de2_hk2 * bl.HK2_T2 + (3.15 + 1.72 / (bl.HK2 - 1.0));
    bl.DE2_D2 = de2_hk2 * bl.HK2_D2 + 1.0;
    bl.DE2_MS = de2_hk2 * bl.HK2_MS;

    const double hdmax = 12.0;
    if (bl.DE2 > hdmax * bl.T2) {
        bl.DE2 = hdmax * bl.T2;
        bl.DE2_U2 = 0.0;
        bl.DE2_T2 = hdmax;
        bl.DE2_D2 = 0.0;
        bl.DE2_MS = 0.0;
    }
    sync_vars_to_com(bl, 2);
}

void blmid(XBlState &bl, int ityp) {
    if (bl.SIMI) {
        bl.HK1 = bl.HK2;
        bl.HK1_T1 = bl.HK2_T2;
        bl.HK1_D1 = bl.HK2_D2;
        bl.HK1_U1 = bl.HK2_U2;
        bl.HK1_MS = bl.HK2_MS;
        bl.RT1 = bl.RT2;
        bl.RT1_T1 = bl.RT2_T2;
        bl.RT1_U1 = bl.RT2_U2;
        bl.RT1_MS = bl.RT2_MS;
        bl.RT1_RE = bl.RT2_RE;
        bl.M1 = bl.M2;
        bl.M1_U1 = bl.M2_U2;
        bl.M1_MS = bl.M2_MS;
        sync_vars_to_com(bl, 1);
    }

    const double hka = 0.5 * (bl.HK1 + bl.HK2);
    const double rta = 0.5 * (bl.RT1 + bl.RT2);
    const double ma = 0.5 * (bl.M1 + bl.M2);

    if (ityp == 3) {
        bl.CFM = 0.0;
        bl.CFM_HKA = 0.0;
        bl.CFM_RTA = 0.0;
        bl.CFM_MA = 0.0;
        bl.CFM_MS = 0.0;
    } else if (ityp == 1) {
        auto cfm_out = cfl(hka, rta, ma);
        bl.CFM = std::get<0>(cfm_out);
        bl.CFM_HKA = std::get<1>(cfm_out);
        bl.CFM_RTA = std::get<2>(cfm_out);
        bl.CFM_MA = std::get<3>(cfm_out);
    } else {
        auto cfm_out = cft(hka, rta, ma, bl.CFFAC);
        bl.CFM = std::get<0>(cfm_out);
        bl.CFM_HKA = std::get<1>(cfm_out);
        bl.CFM_RTA = std::get<2>(cfm_out);
        bl.CFM_MA = std::get<3>(cfm_out);
        auto cfml = cfl(hka, rta, ma);
        if (std::get<0>(cfml) > bl.CFM) {
            bl.CFM = std::get<0>(cfml);
            bl.CFM_HKA = std::get<1>(cfml);
            bl.CFM_RTA = std::get<2>(cfml);
            bl.CFM_MA = std::get<3>(cfml);
        }
    }

    bl.CFM_U1 = 0.5 * (bl.CFM_HKA * bl.HK1_U1 + bl.CFM_MA * bl.M1_U1 + bl.CFM_RTA * bl.RT1_U1);
    bl.CFM_T1 = 0.5 * (bl.CFM_HKA * bl.HK1_T1 + bl.CFM_RTA * bl.RT1_T1);
    bl.CFM_D1 = 0.5 * (bl.CFM_HKA * bl.HK1_D1);

    bl.CFM_U2 = 0.5 * (bl.CFM_HKA * bl.HK2_U2 + bl.CFM_MA * bl.M2_U2 + bl.CFM_RTA * bl.RT2_U2);
    bl.CFM_T2 = 0.5 * (bl.CFM_HKA * bl.HK2_T2 + bl.CFM_RTA * bl.RT2_T2);
    bl.CFM_D2 = 0.5 * (bl.CFM_HKA * bl.HK2_D2);

    bl.CFM_MS = 0.5
                * (bl.CFM_HKA * bl.HK1_MS + bl.CFM_MA * bl.M1_MS + bl.CFM_RTA * bl.RT1_MS + bl.CFM_HKA * bl.HK2_MS
                   + bl.CFM_MA * bl.M2_MS + bl.CFM_RTA * bl.RT2_MS);
    bl.CFM_RE = 0.5 * (bl.CFM_RTA * bl.RT1_RE + bl.CFM_RTA * bl.RT2_RE);
}

void trdif(XBlState &bl) {
    auto bl1 = create2d(4, 5);
    auto bl2 = create2d(4, 5);
    std::vector<double> blrez(static_cast<size_t>(5), 0.0);
    std::vector<double> blm(static_cast<size_t>(5), 0.0);
    std::vector<double> blr(static_cast<size_t>(5), 0.0);
    std::vector<double> blx(static_cast<size_t>(5), 0.0);
    auto bt1 = create2d(4, 5);
    auto bt2 = create2d(4, 5);
    std::vector<double> btrez(static_cast<size_t>(5), 0.0);
    std::vector<double> btm(static_cast<size_t>(5), 0.0);
    std::vector<double> btr(static_cast<size_t>(5), 0.0);
    std::vector<double> btx(static_cast<size_t>(5), 0.0);

    for (int icom = 1; icom <= NCOM; ++icom) {
        bl.C1SAV[icom] = bl.COM1[icom];
        bl.C2SAV[icom] = bl.COM2[icom];
    }

    const double wf2 = (bl.XT - bl.X1) / (bl.X2 - bl.X1);
    const double wf2_xt = 1.0 / (bl.X2 - bl.X1);

    const double wf2_a1 = wf2_xt * bl.XT_A1;
    const double wf2_x1 = wf2_xt * bl.XT_X1 + (wf2 - 1.0) / (bl.X2 - bl.X1);
    const double wf2_x2 = wf2_xt * bl.XT_X2 - wf2 / (bl.X2 - bl.X1);
    const double wf2_t1 = wf2_xt * bl.XT_T1;
    const double wf2_t2 = wf2_xt * bl.XT_T2;
    const double wf2_d1 = wf2_xt * bl.XT_D1;
    const double wf2_d2 = wf2_xt * bl.XT_D2;
    const double wf2_u1 = wf2_xt * bl.XT_U1;
    const double wf2_u2 = wf2_xt * bl.XT_U2;
    const double wf2_ms = wf2_xt * bl.XT_MS;
    const double wf2_re = wf2_xt * bl.XT_RE;
    const double wf2_xf = wf2_xt * bl.XT_XF;

    const double wf1 = 1.0 - wf2;
    const double wf1_a1 = -wf2_a1;
    const double wf1_x1 = -wf2_x1;
    const double wf1_x2 = -wf2_x2;
    const double wf1_t1 = -wf2_t1;
    const double wf1_t2 = -wf2_t2;
    const double wf1_d1 = -wf2_d1;
    const double wf1_d2 = -wf2_d2;
    const double wf1_u1 = -wf2_u1;
    const double wf1_u2 = -wf2_u2;
    const double wf1_ms = -wf2_ms;
    const double wf1_re = -wf2_re;
    const double wf1_xf = -wf2_xf;

    const double tt = bl.T1 * wf1 + bl.T2 * wf2;
    const double tt_a1 = bl.T1 * wf1_a1 + bl.T2 * wf2_a1;
    const double tt_x1 = bl.T1 * wf1_x1 + bl.T2 * wf2_x1;
    const double tt_x2 = bl.T1 * wf1_x2 + bl.T2 * wf2_x2;
    const double tt_t1 = bl.T1 * wf1_t1 + bl.T2 * wf2_t1 + wf1;
    const double tt_t2 = bl.T1 * wf1_t2 + bl.T2 * wf2_t2 + wf2;
    const double tt_d1 = bl.T1 * wf1_d1 + bl.T2 * wf2_d1;
    const double tt_d2 = bl.T1 * wf1_d2 + bl.T2 * wf2_d2;
    const double tt_u1 = bl.T1 * wf1_u1 + bl.T2 * wf2_u1;
    const double tt_u2 = bl.T1 * wf1_u2 + bl.T2 * wf2_u2;
    const double tt_ms = bl.T1 * wf1_ms + bl.T2 * wf2_ms;
    const double tt_re = bl.T1 * wf1_re + bl.T2 * wf2_re;
    const double tt_xf = bl.T1 * wf1_xf + bl.T2 * wf2_xf;

    const double dt = bl.D1 * wf1 + bl.D2 * wf2;
    const double dt_a1 = bl.D1 * wf1_a1 + bl.D2 * wf2_a1;
    const double dt_x1 = bl.D1 * wf1_x1 + bl.D2 * wf2_x1;
    const double dt_x2 = bl.D1 * wf1_x2 + bl.D2 * wf2_x2;
    const double dt_t1 = bl.D1 * wf1_t1 + bl.D2 * wf2_t1;
    const double dt_t2 = bl.D1 * wf1_t2 + bl.D2 * wf2_t2;
    const double dt_d1 = bl.D1 * wf1_d1 + bl.D2 * wf2_d1 + wf1;
    const double dt_d2 = bl.D1 * wf1_d2 + bl.D2 * wf2_d2 + wf2;
    const double dt_u1 = bl.D1 * wf1_u1 + bl.D2 * wf2_u1;
    const double dt_u2 = bl.D1 * wf1_u2 + bl.D2 * wf2_u2;
    const double dt_ms = bl.D1 * wf1_ms + bl.D2 * wf2_ms;
    const double dt_re = bl.D1 * wf1_re + bl.D2 * wf2_re;
    const double dt_xf = bl.D1 * wf1_xf + bl.D2 * wf2_xf;

    const double ut = bl.U1 * wf1 + bl.U2 * wf2;
    const double ut_a1 = bl.U1 * wf1_a1 + bl.U2 * wf2_a1;
    const double ut_x1 = bl.U1 * wf1_x1 + bl.U2 * wf2_x1;
    const double ut_x2 = bl.U1 * wf1_x2 + bl.U2 * wf2_x2;
    const double ut_t1 = bl.U1 * wf1_t1 + bl.U2 * wf2_t1;
    const double ut_t2 = bl.U1 * wf1_t2 + bl.U2 * wf2_t2;
    const double ut_d1 = bl.U1 * wf1_d1 + bl.U2 * wf2_d1;
    const double ut_d2 = bl.U1 * wf1_d2 + bl.U2 * wf2_d2;
    const double ut_u1 = bl.U1 * wf1_u1 + bl.U2 * wf2_u1 + wf1;
    const double ut_u2 = bl.U1 * wf1_u2 + bl.U2 * wf2_u2 + wf2;
    const double ut_ms = bl.U1 * wf1_ms + bl.U2 * wf2_ms;
    const double ut_re = bl.U1 * wf1_re + bl.U2 * wf2_re;
    const double ut_xf = bl.U1 * wf1_xf + bl.U2 * wf2_xf;

    bl.X2 = bl.XT;
    bl.T2 = tt;
    bl.D2 = dt;
    bl.U2 = ut;

    bl.AMPL2 = bl.AMCRIT;
    bl.S2 = 0.0;
    sync_vars_to_com(bl, 2);

    blkin(bl);
    blvar(bl, 1);
    blmid(bl, 1);
    bldif(bl, 1);

    for (int k = 2; k <= 3; ++k) {
        blrez[k] = bl.VSREZ[k];
        blm[k] = bl.VSM[k] + bl.VS2[k][2] * tt_ms + bl.VS2[k][3] * dt_ms + bl.VS2[k][4] * ut_ms + bl.VS2[k][5] * bl.XT_MS;
        blr[k] = bl.VSR[k] + bl.VS2[k][2] * tt_re + bl.VS2[k][3] * dt_re + bl.VS2[k][4] * ut_re + bl.VS2[k][5] * bl.XT_RE;
        blx[k] = bl.VSX[k] + bl.VS2[k][2] * tt_xf + bl.VS2[k][3] * dt_xf + bl.VS2[k][4] * ut_xf + bl.VS2[k][5] * bl.XT_XF;

        bl1[k][1] = bl.VS1[k][1] + bl.VS2[k][2] * tt_a1 + bl.VS2[k][3] * dt_a1 + bl.VS2[k][4] * ut_a1 + bl.VS2[k][5] * bl.XT_A1;
        bl1[k][2] = bl.VS1[k][2] + bl.VS2[k][2] * tt_t1 + bl.VS2[k][3] * dt_t1 + bl.VS2[k][4] * ut_t1 + bl.VS2[k][5] * bl.XT_T1;
        bl1[k][3] = bl.VS1[k][3] + bl.VS2[k][2] * tt_d1 + bl.VS2[k][3] * dt_d1 + bl.VS2[k][4] * ut_d1 + bl.VS2[k][5] * bl.XT_D1;
        bl1[k][4] = bl.VS1[k][4] + bl.VS2[k][2] * tt_u1 + bl.VS2[k][3] * dt_u1 + bl.VS2[k][4] * ut_u1 + bl.VS2[k][5] * bl.XT_U1;
        bl1[k][5] = bl.VS1[k][5] + bl.VS2[k][2] * tt_x1 + bl.VS2[k][3] * dt_x1 + bl.VS2[k][4] * ut_x1 + bl.VS2[k][5] * bl.XT_X1;

        bl2[k][1] = 0.0;
        bl2[k][2] = bl.VS2[k][2] * tt_t2 + bl.VS2[k][3] * dt_t2 + bl.VS2[k][4] * ut_t2 + bl.VS2[k][5] * bl.XT_T2;
        bl2[k][3] = bl.VS2[k][2] * tt_d2 + bl.VS2[k][3] * dt_d2 + bl.VS2[k][4] * ut_d2 + bl.VS2[k][5] * bl.XT_D2;
        bl2[k][4] = bl.VS2[k][2] * tt_u2 + bl.VS2[k][3] * dt_u2 + bl.VS2[k][4] * ut_u2 + bl.VS2[k][5] * bl.XT_U2;
        bl2[k][5] = bl.VS2[k][2] * tt_x2 + bl.VS2[k][3] * dt_x2 + bl.VS2[k][4] * ut_x2 + bl.VS2[k][5] * bl.XT_X2;
    }

    blvar(bl, 2);
    const double ctr = bl.CTRCON * std::exp(-bl.CTRCEX / (bl.HK2 - 1.0));
    const double ctr_hk2 = ctr * bl.CTRCEX / std::pow(bl.HK2 - 1.0, 2);

    const double st = ctr * bl.CQ2;
    const double st_tt = ctr * bl.CQ2_T2 + bl.CQ2 * ctr_hk2 * bl.HK2_T2;
    const double st_dt = ctr * bl.CQ2_D2 + bl.CQ2 * ctr_hk2 * bl.HK2_D2;
    const double st_ut = ctr * bl.CQ2_U2 + bl.CQ2 * ctr_hk2 * bl.HK2_U2;
    const double st_ms = ctr * bl.CQ2_MS + bl.CQ2 * ctr_hk2 * bl.HK2_MS;
    const double st_re = ctr * bl.CQ2_RE;

    const double st_a1 = st_tt * tt_a1 + st_dt * dt_a1 + st_ut * ut_a1;
    const double st_x1 = st_tt * tt_x1 + st_dt * dt_x1 + st_ut * ut_x1;
    const double st_x2 = st_tt * tt_x2 + st_dt * dt_x2 + st_ut * ut_x2;
    const double st_t1 = st_tt * tt_t1 + st_dt * dt_t1 + st_ut * ut_t1;
    const double st_t2 = st_tt * tt_t2 + st_dt * dt_t2 + st_ut * ut_t2;
    const double st_d1 = st_tt * tt_d1 + st_dt * dt_d1 + st_ut * ut_d1;
    const double st_d2 = st_tt * tt_d2 + st_dt * dt_d2 + st_ut * ut_d2;
    const double st_u1 = st_tt * tt_u1 + st_dt * dt_u1 + st_ut * ut_u1;
    const double st_u2 = st_tt * tt_u2 + st_dt * dt_u2 + st_ut * ut_u2;
    const double st_ms2 = st_tt * tt_ms + st_dt * dt_ms + st_ut * ut_ms + st_ms;
    const double st_re2 = st_tt * tt_re + st_dt * dt_re + st_ut * ut_re + st_re;
    const double st_xf = st_tt * tt_xf + st_dt * dt_xf + st_ut * ut_xf;

    bl.AMPL2 = 0.0;
    bl.S2 = st;

    blvar(bl, 2);

    for (int icom = 1; icom <= NCOM; ++icom) {
        bl.COM1[icom] = bl.COM2[icom];
        bl.COM2[icom] = bl.C2SAV[icom];
    }
    sync_com_to_vars(bl, 1);
    sync_com_to_vars(bl, 2);

    blmid(bl, 2);
    bldif(bl, 2);

    for (int k = 1; k <= 3; ++k) {
        btrez[k] = bl.VSREZ[k];
        btm[k] = bl.VSM[k] + bl.VS1[k][1] * st_ms2 + bl.VS1[k][2] * tt_ms + bl.VS1[k][3] * dt_ms + bl.VS1[k][4] * ut_ms
                 + bl.VS1[k][5] * bl.XT_MS;
        btr[k] = bl.VSR[k] + bl.VS1[k][1] * st_re2 + bl.VS1[k][2] * tt_re + bl.VS1[k][3] * dt_re + bl.VS1[k][4] * ut_re
                 + bl.VS1[k][5] * bl.XT_RE;
        btx[k] = bl.VSX[k] + bl.VS1[k][1] * st_xf + bl.VS1[k][2] * tt_xf + bl.VS1[k][3] * dt_xf + bl.VS1[k][4] * ut_xf
                 + bl.VS1[k][5] * bl.XT_XF;

        bt1[k][1] = bl.VS1[k][1] * st_a1 + bl.VS1[k][2] * tt_a1 + bl.VS1[k][3] * dt_a1 + bl.VS1[k][4] * ut_a1
                    + bl.VS1[k][5] * bl.XT_A1;
        bt1[k][2] = bl.VS1[k][1] * st_t1 + bl.VS1[k][2] * tt_t1 + bl.VS1[k][3] * dt_t1 + bl.VS1[k][4] * ut_t1
                    + bl.VS1[k][5] * bl.XT_T1;
        bt1[k][3] = bl.VS1[k][1] * st_d1 + bl.VS1[k][2] * tt_d1 + bl.VS1[k][3] * dt_d1 + bl.VS1[k][4] * ut_d1
                    + bl.VS1[k][5] * bl.XT_D1;
        bt1[k][4] = bl.VS1[k][1] * st_u1 + bl.VS1[k][2] * tt_u1 + bl.VS1[k][3] * dt_u1 + bl.VS1[k][4] * ut_u1
                    + bl.VS1[k][5] * bl.XT_U1;
        bt1[k][5] = bl.VS1[k][1] * st_x1 + bl.VS1[k][2] * tt_x1 + bl.VS1[k][3] * dt_x1 + bl.VS1[k][4] * ut_x1
                    + bl.VS1[k][5] * bl.XT_X1;

        bt2[k][1] = bl.VS2[k][1];
        bt2[k][2] = bl.VS2[k][2] + bl.VS1[k][1] * st_t2 + bl.VS1[k][2] * tt_t2 + bl.VS1[k][3] * dt_t2 + bl.VS1[k][4] * ut_t2
                    + bl.VS1[k][5] * bl.XT_T2;
        bt2[k][3] = bl.VS2[k][3] + bl.VS1[k][1] * st_d2 + bl.VS1[k][2] * tt_d2 + bl.VS1[k][3] * dt_d2 + bl.VS1[k][4] * ut_d2
                    + bl.VS1[k][5] * bl.XT_D2;
        bt2[k][4] = bl.VS2[k][4] + bl.VS1[k][1] * st_u2 + bl.VS1[k][2] * tt_u2 + bl.VS1[k][3] * dt_u2 + bl.VS1[k][4] * ut_u2
                    + bl.VS1[k][5] * bl.XT_U2;
        bt2[k][5] = bl.VS2[k][5] + bl.VS1[k][1] * st_x2 + bl.VS1[k][2] * tt_x2 + bl.VS1[k][3] * dt_x2 + bl.VS1[k][4] * ut_x2
                    + bl.VS1[k][5] * bl.XT_X2;
    }

    bl.VSREZ[1] = btrez[1];
    bl.VSREZ[2] = blrez[2] + btrez[2];
    bl.VSREZ[3] = blrez[3] + btrez[3];
    bl.VSM[1] = btm[1];
    bl.VSM[2] = blm[2] + btm[2];
    bl.VSM[3] = blm[3] + btm[3];
    bl.VSR[1] = btr[1];
    bl.VSR[2] = blr[2] + btr[2];
    bl.VSR[3] = blr[3] + btr[3];
    bl.VSX[1] = btx[1];
    bl.VSX[2] = blx[2] + btx[2];
    bl.VSX[3] = blx[3] + btx[3];
    for (int l = 1; l <= 5; ++l) {
        bl.VS1[1][l] = bt1[1][l];
        bl.VS2[1][l] = bt2[1][l];
        bl.VS1[2][l] = bl1[2][l] + bt1[2][l];
        bl.VS2[2][l] = bl2[2][l] + bt2[2][l];
        bl.VS1[3][l] = bl1[3][l] + bt1[3][l];
        bl.VS2[3][l] = bl2[3][l] + bt2[3][l];
    }

    for (int icom = 1; icom <= NCOM; ++icom) {
        bl.COM1[icom] = bl.C1SAV[icom];
    }
    sync_com_to_vars(bl, 1);
}

void bldif(XBlState &bl, int ityp) {
    double xlog = 0.0;
    double ulog = 0.0;
    double tlog = 0.0;
    double hlog = 0.0;
    double ddlog = 0.0;
    if (ityp == 0) {
        xlog = 1.0;
        ulog = bl.BULE;
        tlog = 0.5 * (1.0 - bl.BULE);
        hlog = 0.0;
        ddlog = 0.0;
    } else {
        xlog = std::log(bl.X2 / bl.X1);
        ulog = std::log(bl.U2 / bl.U1);
        tlog = std::log(bl.T2 / bl.T1);
        hlog = std::log(bl.HS2 / bl.HS1);
        ddlog = 1.0;
    }

    for (int k = 1; k <= 4; ++k) {
        bl.VSREZ[k] = 0.0;
        bl.VSM[k] = 0.0;
        bl.VSR[k] = 0.0;
        bl.VSX[k] = 0.0;
        for (int l = 1; l <= 5; ++l) {
            bl.VS1[k][l] = 0.0;
            bl.VS2[k][l] = 0.0;
        }
    }

    const double hupwt = 1.0;

    double hdcon = 5.0 * hupwt / (bl.HK2 * bl.HK2);
    double hd_hk1 = 0.0;
    double hd_hk2 = -hdcon * 2.0 / bl.HK2;

    if (ityp == 3) {
        hdcon = hupwt / (bl.HK2 * bl.HK2);
        hd_hk1 = 0.0;
        hd_hk2 = -hdcon * 2.0 / bl.HK2;
    }

    const double arg = std::abs((bl.HK2 - 1.0) / (bl.HK1 - 1.0));
    const double hl = std::log(arg);
    const double hl_hk1 = -1.0 / (bl.HK1 - 1.0);
    const double hl_hk2 = 1.0 / (bl.HK2 - 1.0);

    const double hlsq = std::min(hl * hl, 15.0);
    const double ehh = std::exp(-hlsq * hdcon);
    const double upw = 1.0 - 0.5 * ehh;
    const double upw_hl = ehh * hl * hdcon;
    const double upw_hd = 0.5 * ehh * hlsq;

    const double upw_hk1 = upw_hl * hl_hk1 + upw_hd * hd_hk1;
    const double upw_hk2 = upw_hl * hl_hk2 + upw_hd * hd_hk2;

    const double upw_u1 = upw_hk1 * bl.HK1_U1;
    const double upw_t1 = upw_hk1 * bl.HK1_T1;
    const double upw_d1 = upw_hk1 * bl.HK1_D1;
    const double upw_u2 = upw_hk2 * bl.HK2_U2;
    const double upw_t2 = upw_hk2 * bl.HK2_T2;
    const double upw_d2 = upw_hk2 * bl.HK2_D2;
    const double upw_ms = upw_hk1 * bl.HK1_MS + upw_hk2 * bl.HK2_MS;

    if (ityp == 0) {
        bl.VS2[1][1] = 1.0;
        bl.VSR[1] = 0.0;
        bl.VSREZ[1] = -bl.AMPL2;
    } else if (ityp == 1) {
        auto ax_out = axset(bl.HK1, bl.T1, bl.RT1, bl.AMPL1, bl.HK2, bl.T2, bl.RT2, bl.AMPL2, bl.AMCRIT, bl.IDAMPV);
        const double ax = std::get<0>(ax_out);
        const double rezc = bl.AMPL2 - bl.AMPL1 - ax * (bl.X2 - bl.X1);
        const double z_ax = -(bl.X2 - bl.X1);

        bl.VS1[1][1] = z_ax * std::get<4>(ax_out) - 1.0;
        bl.VS1[1][2] = z_ax * (std::get<1>(ax_out) * bl.HK1_T1 + std::get<2>(ax_out) + std::get<3>(ax_out) * bl.RT1_T1);
        bl.VS1[1][3] = z_ax * (std::get<1>(ax_out) * bl.HK1_D1);
        bl.VS1[1][4] = z_ax * (std::get<1>(ax_out) * bl.HK1_U1 + std::get<3>(ax_out) * bl.RT1_U1);
        bl.VS1[1][5] = ax;
        bl.VS2[1][1] = z_ax * std::get<8>(ax_out) + 1.0;
        bl.VS2[1][2] = z_ax * (std::get<5>(ax_out) * bl.HK2_T2 + std::get<6>(ax_out) + std::get<7>(ax_out) * bl.RT2_T2);
        bl.VS2[1][3] = z_ax * (std::get<5>(ax_out) * bl.HK2_D2);
        bl.VS2[1][4] = z_ax * (std::get<5>(ax_out) * bl.HK2_U2 + std::get<7>(ax_out) * bl.RT2_U2);
        bl.VS2[1][5] = -ax;
        bl.VSM[1] = z_ax * (std::get<1>(ax_out) * bl.HK1_MS + std::get<3>(ax_out) * bl.RT1_MS + std::get<5>(ax_out) * bl.HK2_MS
                            + std::get<7>(ax_out) * bl.RT2_MS);
        bl.VSR[1] = z_ax * (std::get<3>(ax_out) * bl.RT1_RE + std::get<7>(ax_out) * bl.RT2_RE);
        bl.VSX[1] = 0.0;
        bl.VSREZ[1] = -rezc;
    } else {
        const double sa = (1.0 - upw) * bl.S1 + upw * bl.S2;
        const double cqa = (1.0 - upw) * bl.CQ1 + upw * bl.CQ2;
        const double cfa = (1.0 - upw) * bl.CF1 + upw * bl.CF2;
        const double hka = (1.0 - upw) * bl.HK1 + upw * bl.HK2;

        const double usa = 0.5 * (bl.US1 + bl.US2);
        const double rta = 0.5 * (bl.RT1 + bl.RT2);
        const double dea = 0.5 * (bl.DE1 + bl.DE2);
        const double da = 0.5 * (bl.D1 + bl.D2);

        const double ald = (ityp == 3) ? bl.DLCON : 1.0;

        double hkc = 0.0;
        double hkc_hka = 0.0;
        double hkc_rta = 0.0;
        if (ityp == 2) {
            const double gcc = bl.GCCON;
            hkc = hka - 1.0 - gcc / rta;
            hkc_hka = 1.0;
            hkc_rta = gcc / (rta * rta);
            if (hkc < 0.01) {
                hkc = 0.01;
                hkc_hka = 0.0;
                hkc_rta = 0.0;
            }
        } else {
            hkc = hka - 1.0;
            hkc_hka = 1.0;
            hkc_rta = 0.0;
        }

        const double hr = hkc / (bl.GACON * ald * hka);
        const double hr_hka = hkc_hka / (bl.GACON * ald * hka) - hr / hka;
        const double hr_rta = hkc_rta / (bl.GACON * ald * hka);

        const double uq = (0.5 * cfa - hr * hr) / (bl.GBCON * da);
        const double uq_hka = -2.0 * hr * hr_hka / (bl.GBCON * da);
        const double uq_rta = -2.0 * hr * hr_rta / (bl.GBCON * da);
        const double uq_cfa = 0.5 / (bl.GBCON * da);
        const double uq_da = -uq / da;
        const double uq_upw = uq_cfa * (bl.CF2 - bl.CF1) + uq_hka * (bl.HK2 - bl.HK1);

        double uq_t1 = (1.0 - upw) * (uq_cfa * bl.CF1_T1 + uq_hka * bl.HK1_T1) + uq_upw * upw_t1;
        double uq_d1 = (1.0 - upw) * (uq_cfa * bl.CF1_D1 + uq_hka * bl.HK1_D1) + uq_upw * upw_d1;
        double uq_u1 = (1.0 - upw) * (uq_cfa * bl.CF1_U1 + uq_hka * bl.HK1_U1) + uq_upw * upw_u1;
        double uq_t2 = upw * (uq_cfa * bl.CF2_T2 + uq_hka * bl.HK2_T2) + uq_upw * upw_t2;
        double uq_d2 = upw * (uq_cfa * bl.CF2_D2 + uq_hka * bl.HK2_D2) + uq_upw * upw_d2;
        double uq_u2 = upw * (uq_cfa * bl.CF2_U2 + uq_hka * bl.HK2_U2) + uq_upw * upw_u2;
        double uq_ms = (1.0 - upw) * (uq_cfa * bl.CF1_MS + uq_hka * bl.HK1_MS) + uq_upw * upw_ms
                       + upw * (uq_cfa * bl.CF2_MS + uq_hka * bl.HK2_MS);
        double uq_re = (1.0 - upw) * uq_cfa * bl.CF1_RE + upw * uq_cfa * bl.CF2_RE;

        uq_t1 = uq_t1 + 0.5 * uq_rta * bl.RT1_T1;
        uq_d1 = uq_d1 + 0.5 * uq_da;
        uq_u1 = uq_u1 + 0.5 * uq_rta * bl.RT1_U1;
        uq_t2 = uq_t2 + 0.5 * uq_rta * bl.RT2_T2;
        uq_d2 = uq_d2 + 0.5 * uq_da;
        uq_u2 = uq_u2 + 0.5 * uq_rta * bl.RT2_U2;
        uq_ms = uq_ms + 0.5 * uq_rta * bl.RT1_MS + 0.5 * uq_rta * bl.RT2_MS;
        uq_re = uq_re + 0.5 * uq_rta * bl.RT1_RE + 0.5 * uq_rta * bl.RT2_RE;

        const double scc = bl.SCCON * 1.333 / (1.0 + usa);
        const double scc_usa = -scc / (1.0 + usa);

        const double scc_us1 = scc_usa * 0.5;
        const double scc_us2 = scc_usa * 0.5;

        const double slog = std::log(bl.S2 / bl.S1);
        const double dxi = bl.X2 - bl.X1;

        const double rezc = scc * (cqa - sa * ald) * dxi - dea * 2.0 * slog + dea * 2.0 * (uq * dxi - ulog) * bl.DUXCON;

        const double z_cfa = dea * 2.0 * uq_cfa * dxi * bl.DUXCON;
        const double z_hka = dea * 2.0 * uq_hka * dxi * bl.DUXCON;
        const double z_da = dea * 2.0 * uq_da * dxi * bl.DUXCON;
        const double z_sl = -dea * 2.0;
        const double z_ul = -dea * 2.0 * bl.DUXCON;
        const double z_dxi = scc * (cqa - sa * ald) + dea * 2.0 * uq * bl.DUXCON;
        const double z_usa = scc_usa * (cqa - sa * ald) * dxi;
        const double z_cqa = scc * dxi;
        const double z_sa = -scc * dxi * ald;
        const double z_dea = 2.0 * ((uq * dxi - ulog) * bl.DUXCON - slog);

        const double z_upw = z_cqa * (bl.CQ2 - bl.CQ1) + z_sa * (bl.S2 - bl.S1) + z_cfa * (bl.CF2 - bl.CF1)
                             + z_hka * (bl.HK2 - bl.HK1);
        const double z_de1 = 0.5 * z_dea;
        const double z_de2 = 0.5 * z_dea;
        const double z_us1 = 0.5 * z_usa;
        const double z_us2 = 0.5 * z_usa;
        const double z_d1 = 0.5 * z_da;
        const double z_d2 = 0.5 * z_da;
        const double z_u1 = -z_ul / bl.U1;
        const double z_u2 = z_ul / bl.U2;
        const double z_x1 = -z_dxi;
        const double z_x2 = z_dxi;
        const double z_s1 = (1.0 - upw) * z_sa - z_sl / bl.S1;
        const double z_s2 = upw * z_sa + z_sl / bl.S2;
        const double z_cq1 = (1.0 - upw) * z_cqa;
        const double z_cq2 = upw * z_cqa;
        const double z_cf1 = (1.0 - upw) * z_cfa;
        const double z_cf2 = upw * z_cfa;
        const double z_hk1 = (1.0 - upw) * z_hka;
        const double z_hk2 = upw * z_hka;

        bl.VS1[1][1] = z_s1;
        bl.VS1[1][2] = z_upw * upw_t1 + z_de1 * bl.DE1_T1 + z_us1 * bl.US1_T1;
        bl.VS1[1][3] = z_d1 + z_upw * upw_d1 + z_de1 * bl.DE1_D1 + z_us1 * bl.US1_D1;
        bl.VS1[1][4] = z_u1 + z_upw * upw_u1 + z_de1 * bl.DE1_U1 + z_us1 * bl.US1_U1;
        bl.VS1[1][5] = z_x1;
        bl.VS2[1][1] = z_s2;
        bl.VS2[1][2] = z_upw * upw_t2 + z_de2 * bl.DE2_T2 + z_us2 * bl.US2_T2;
        bl.VS2[1][3] = z_d2 + z_upw * upw_d2 + z_de2 * bl.DE2_D2 + z_us2 * bl.US2_D2;
        bl.VS2[1][4] = z_u2 + z_upw * upw_u2 + z_de2 * bl.DE2_U2 + z_us2 * bl.US2_U2;
        bl.VS2[1][5] = z_x2;
        bl.VSM[1] = z_upw * upw_ms + z_de1 * bl.DE1_MS + z_us1 * bl.US1_MS + z_de2 * bl.DE2_MS + z_us2 * bl.US2_MS;

        bl.VS1[1][2] = bl.VS1[1][2] + z_cq1 * bl.CQ1_T1 + z_cf1 * bl.CF1_T1 + z_hk1 * bl.HK1_T1;
        bl.VS1[1][3] = bl.VS1[1][3] + z_cq1 * bl.CQ1_D1 + z_cf1 * bl.CF1_D1 + z_hk1 * bl.HK1_D1;
        bl.VS1[1][4] = bl.VS1[1][4] + z_cq1 * bl.CQ1_U1 + z_cf1 * bl.CF1_U1 + z_hk1 * bl.HK1_U1;

        bl.VS2[1][2] = bl.VS2[1][2] + z_cq2 * bl.CQ2_T2 + z_cf2 * bl.CF2_T2 + z_hk2 * bl.HK2_T2;
        bl.VS2[1][3] = bl.VS2[1][3] + z_cq2 * bl.CQ2_D2 + z_cf2 * bl.CF2_D2 + z_hk2 * bl.HK2_D2;
        bl.VS2[1][4] = bl.VS2[1][4] + z_cq2 * bl.CQ2_U2 + z_cf2 * bl.CF2_U2 + z_hk2 * bl.HK2_U2;

        bl.VSM[1] = bl.VSM[1] + z_cq1 * bl.CQ1_MS + z_cf1 * bl.CF1_MS + z_hk1 * bl.HK1_MS + z_cq2 * bl.CQ2_MS
                    + z_cf2 * bl.CF2_MS + z_hk2 * bl.HK2_MS;
        bl.VSR[1] = z_cq1 * bl.CQ1_RE + z_cf1 * bl.CF1_RE + z_cq2 * bl.CQ2_RE + z_cf2 * bl.CF2_RE;
        bl.VSX[1] = 0.0;
        bl.VSREZ[1] = -rezc;
    }

    const double ha = 0.5 * (bl.H1 + bl.H2);
    const double ma = 0.5 * (bl.M1 + bl.M2);
    const double xa = 0.5 * (bl.X1 + bl.X2);
    const double ta = 0.5 * (bl.T1 + bl.T2);
    const double hwa = 0.5 * (bl.DW1 / bl.T1 + bl.DW2 / bl.T2);

    const double cfx = 0.50 * bl.CFM * xa / ta + 0.25 * (bl.CF1 * bl.X1 / bl.T1 + bl.CF2 * bl.X2 / bl.T2);
    const double cfx_xa = 0.50 * bl.CFM / ta;
    const double cfx_ta = -0.50 * bl.CFM * xa / (ta * ta);

    const double cfx_x1 = 0.25 * bl.CF1 / bl.T1 + cfx_xa * 0.5;
    const double cfx_x2 = 0.25 * bl.CF2 / bl.T2 + cfx_xa * 0.5;
    const double cfx_t1 = -0.25 * bl.CF1 * bl.X1 / (bl.T1 * bl.T1) + cfx_ta * 0.5;
    const double cfx_t2 = -0.25 * bl.CF2 * bl.X2 / (bl.T2 * bl.T2) + cfx_ta * 0.5;
    const double cfx_cf1 = 0.25 * bl.X1 / bl.T1;
    const double cfx_cf2 = 0.25 * bl.X2 / bl.T2;
    const double cfx_cfm = 0.50 * xa / ta;

    const double btmp = ha + 2.0 - ma + hwa;

    const double rezt = tlog + btmp * ulog - xlog * 0.5 * cfx;
    const double z_cfx = -xlog * 0.5;
    const double z_ha = ulog;
    const double z_hwa = ulog;
    const double z_ma = -ulog;
    const double z_xl = -ddlog * 0.5 * cfx;
    const double z_ul = ddlog * btmp;
    const double z_tl = ddlog;

    const double z_cfm = z_cfx * cfx_cfm;
    const double z_cf1 = z_cfx * cfx_cf1;
    const double z_cf2 = z_cfx * cfx_cf2;

    double z_t1 = -z_tl / bl.T1 + z_cfx * cfx_t1 + z_hwa * 0.5 * (-bl.DW1 / (bl.T1 * bl.T1));
    double z_t2 = z_tl / bl.T2 + z_cfx * cfx_t2 + z_hwa * 0.5 * (-bl.DW2 / (bl.T2 * bl.T2));
    const double z_x1 = -z_xl / bl.X1 + z_cfx * cfx_x1;
    const double z_x2 = z_xl / bl.X2 + z_cfx * cfx_x2;
    const double z_u1 = -z_ul / bl.U1;
    const double z_u2 = z_ul / bl.U2;

    bl.VS1[2][2] = 0.5 * z_ha * bl.H1_T1 + z_cfm * bl.CFM_T1 + z_cf1 * bl.CF1_T1 + z_t1;
    bl.VS1[2][3] = 0.5 * z_ha * bl.H1_D1 + z_cfm * bl.CFM_D1 + z_cf1 * bl.CF1_D1;
    bl.VS1[2][4] = 0.5 * z_ma * bl.M1_U1 + z_cfm * bl.CFM_U1 + z_cf1 * bl.CF1_U1 + z_u1;
    bl.VS1[2][5] = z_x1;
    bl.VS2[2][2] = 0.5 * z_ha * bl.H2_T2 + z_cfm * bl.CFM_T2 + z_cf2 * bl.CF2_T2 + z_t2;
    bl.VS2[2][3] = 0.5 * z_ha * bl.H2_D2 + z_cfm * bl.CFM_D2 + z_cf2 * bl.CF2_D2;
    bl.VS2[2][4] = 0.5 * z_ma * bl.M2_U2 + z_cfm * bl.CFM_U2 + z_cf2 * bl.CF2_U2 + z_u2;
    bl.VS2[2][5] = z_x2;

    bl.VSM[2] = 0.5 * z_ma * bl.M1_MS + z_cfm * bl.CFM_MS + z_cf1 * bl.CF1_MS + 0.5 * z_ma * bl.M2_MS
                + z_cf2 * bl.CF2_MS;
    bl.VSR[2] = z_cfm * bl.CFM_RE + z_cf1 * bl.CF1_RE + z_cf2 * bl.CF2_RE;
    bl.VSX[2] = 0.0;
    bl.VSREZ[2] = -rezt;

    const double xot1 = bl.X1 / bl.T1;
    const double xot2 = bl.X2 / bl.T2;

    const double hsa = 0.5 * (bl.HS1 + bl.HS2);
    const double hca = 0.5 * (bl.HC1 + bl.HC2);
    const double hwa2 = 0.5 * (bl.DW1 / bl.T1 + bl.DW2 / bl.T2);

    const double dix = (1.0 - upw) * bl.DI1 * xot1 + upw * bl.DI2 * xot2;
    const double cfx2 = (1.0 - upw) * bl.CF1 * xot1 + upw * bl.CF2 * xot2;
    const double dix_upw = bl.DI2 * xot2 - bl.DI1 * xot1;
    const double cfx_upw = bl.CF2 * xot2 - bl.CF1 * xot1;

    const double btmp2 = 2.0 * hca / hsa + 1.0 - ha - hwa2;

    const double rezh = hlog + btmp2 * ulog + xlog * (0.5 * cfx2 - dix);
    const double z_cfx2 = xlog * 0.5;
    const double z_dix = -xlog;
    const double z_hca = 2.0 * ulog / hsa;
    const double z_ha2 = -ulog;
    const double z_hwa2 = -ulog;
    const double z_xl2 = ddlog * (0.5 * cfx2 - dix);
    const double z_ul2 = ddlog * btmp2;
    const double z_hl = ddlog;

    const double z_upw2 = z_cfx2 * cfx_upw + z_dix * dix_upw;

    const double z_hs1 = -hca * ulog / (hsa * hsa) - z_hl / bl.HS1;
    const double z_hs2 = -hca * ulog / (hsa * hsa) + z_hl / bl.HS2;

    const double z_cf1b = (1.0 - upw) * z_cfx2 * xot1;
    const double z_cf2b = upw * z_cfx2 * xot2;
    const double z_di1 = (1.0 - upw) * z_dix * xot1;
    const double z_di2 = upw * z_dix * xot2;

    double z_t1b = (1.0 - upw) * (z_cfx2 * bl.CF1 + z_dix * bl.DI1) * (-xot1 / bl.T1);
    double z_t2b = upw * (z_cfx2 * bl.CF2 + z_dix * bl.DI2) * (-xot2 / bl.T2);
    const double z_x1b = (1.0 - upw) * (z_cfx2 * bl.CF1 + z_dix * bl.DI1) / bl.T1 - z_xl2 / bl.X1;
    const double z_x2b = upw * (z_cfx2 * bl.CF2 + z_dix * bl.DI2) / bl.T2 + z_xl2 / bl.X2;
    const double z_u1b = -z_ul2 / bl.U1;
    const double z_u2b = z_ul2 / bl.U2;

    z_t1b = z_t1b + z_hwa2 * 0.5 * (-bl.DW1 / (bl.T1 * bl.T1));
    z_t2b = z_t2b + z_hwa2 * 0.5 * (-bl.DW2 / (bl.T2 * bl.T2));

    bl.VS1[3][1] = z_di1 * bl.DI1_S1;
    bl.VS1[3][2] = z_hs1 * bl.HS1_T1 + z_cf1b * bl.CF1_T1 + z_di1 * bl.DI1_T1 + z_t1b;
    bl.VS1[3][3] = z_hs1 * bl.HS1_D1 + z_cf1b * bl.CF1_D1 + z_di1 * bl.DI1_D1;
    bl.VS1[3][4] = z_hs1 * bl.HS1_U1 + z_cf1b * bl.CF1_U1 + z_di1 * bl.DI1_U1 + z_u1b;
    bl.VS1[3][5] = z_x1b;
    bl.VS2[3][1] = z_di2 * bl.DI2_S2;
    bl.VS2[3][2] = z_hs2 * bl.HS2_T2 + z_cf2b * bl.CF2_T2 + z_di2 * bl.DI2_T2 + z_t2b;
    bl.VS2[3][3] = z_hs2 * bl.HS2_D2 + z_cf2b * bl.CF2_D2 + z_di2 * bl.DI2_D2;
    bl.VS2[3][4] = z_hs2 * bl.HS2_U2 + z_cf2b * bl.CF2_U2 + z_di2 * bl.DI2_U2 + z_u2b;
    bl.VS2[3][5] = z_x2b;
    bl.VSM[3] = z_hs1 * bl.HS1_MS + z_cf1b * bl.CF1_MS + z_di1 * bl.DI1_MS + z_hs2 * bl.HS2_MS + z_cf2b * bl.CF2_MS
                + z_di2 * bl.DI2_MS;
    bl.VSR[3] = z_hs1 * bl.HS1_RE + z_cf1b * bl.CF1_RE + z_di1 * bl.DI1_RE + z_hs2 * bl.HS2_RE + z_cf2b * bl.CF2_RE
                + z_di2 * bl.DI2_RE;

    bl.VS1[3][2] = bl.VS1[3][2] + 0.5 * (z_hca * bl.HC1_T1 + z_ha2 * bl.H1_T1) + z_upw2 * upw_t1;
    bl.VS1[3][3] = bl.VS1[3][3] + 0.5 * (z_hca * bl.HC1_D1 + z_ha2 * bl.H1_D1) + z_upw2 * upw_d1;
    bl.VS1[3][4] = bl.VS1[3][4] + 0.5 * (z_hca * bl.HC1_U1) + z_upw2 * upw_u1;
    bl.VS2[3][2] = bl.VS2[3][2] + 0.5 * (z_hca * bl.HC2_T2 + z_ha2 * bl.H2_T2) + z_upw2 * upw_t2;
    bl.VS2[3][3] = bl.VS2[3][3] + 0.5 * (z_hca * bl.HC2_D2 + z_ha2 * bl.H2_D2) + z_upw2 * upw_d2;
    bl.VS2[3][4] = bl.VS2[3][4] + 0.5 * (z_hca * bl.HC2_U2) + z_upw2 * upw_u2;

    bl.VSM[3] = bl.VSM[3] + 0.5 * (z_hca * bl.HC1_MS) + z_upw2 * upw_ms + 0.5 * (z_hca * bl.HC2_MS);

    bl.VSX[3] = 0.0;
    bl.VSREZ[3] = -rezh;
}

std::tuple<double, double, double> hkin(double h, double msq) {
    return hkin_impl(h, msq);
}

static std::tuple<double, double, double, double> dampl(double hk, double th, double rt) {
    const double dgr = 0.08;
    const double hmi = 1.0 / (hk - 1.0);
    const double hmi_hk = -(hmi * hmi);

    const double aa = 2.492 * std::pow(hmi, 0.43);
    const double aa_hk = (aa / hmi) * 0.43 * hmi_hk;

    const double bb = std::tanh(14.0 * hmi - 9.24);
    const double bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;

    const double grcrit = aa + 0.7 * (bb + 1.0);
    const double grc_hk = aa_hk + 0.7 * bb_hk;

    const double gr = std::log10(rt);
    const double gr_rt = 1.0 / (2.3025851 * rt);

    if (gr < grcrit - dgr) {
        return {0.0, 0.0, 0.0, 0.0};
    }

    const double rnorm = (gr - (grcrit - dgr)) / (2.0 * dgr);
    const double rn_hk = -grc_hk / (2.0 * dgr);
    const double rn_rt = gr_rt / (2.0 * dgr);

    double rfac = 0.0;
    double rfac_hk = 0.0;
    double rfac_rt = 0.0;
    if (rnorm >= 1.0) {
        rfac = 1.0;
        rfac_hk = 0.0;
        rfac_rt = 0.0;
    } else {
        rfac = 3.0 * rnorm * rnorm - 2.0 * rnorm * rnorm * rnorm;
        const double rfac_rn = 6.0 * rnorm - 6.0 * rnorm * rnorm;
        rfac_hk = rfac_rn * rn_hk;
        rfac_rt = rfac_rn * rn_rt;
    }

    const double arg = 3.87 * hmi - 2.52;
    const double arg_hk = 3.87 * hmi_hk;
    const double ex = std::exp(-(arg * arg));
    const double ex_hk = ex * (-2.0 * arg * arg_hk);

    const double dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
    const double dadr_hk = 0.028 - 0.0345 * ex_hk;

    const double af = -0.05 + 2.7 * hmi - 5.5 * hmi * hmi + 3.0 * hmi * hmi * hmi;
    const double af_hmi = 2.7 - 11.0 * hmi + 9.0 * hmi * hmi;
    const double af_hk = af_hmi * hmi_hk;

    const double ax = (af * dadr / th) * rfac;
    const double ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac + (af * dadr / th) * rfac_hk;
    const double ax_th = -ax / th;
    const double ax_rt = (af * dadr / th) * rfac_rt;

    return {ax, ax_hk, ax_th, ax_rt};
}

static std::tuple<double, double, double, double> dampl2(double hk, double th, double rt) {
    const double dgr = 0.08;
    const double hk1 = 3.5;
    const double hk2 = 4.0;
    const double hmi = 1.0 / (hk - 1.0);
    const double hmi_hk = -(hmi * hmi);

    const double aa = 2.492 * std::pow(hmi, 0.43);
    const double aa_hk = (aa / hmi) * 0.43 * hmi_hk;

    const double bb = std::tanh(14.0 * hmi - 9.24);
    const double bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;

    const double grc = aa + 0.7 * (bb + 1.0);
    const double grc_hk = aa_hk + 0.7 * bb_hk;

    const double gr = std::log10(rt);
    const double gr_rt = 1.0 / (2.3025851 * rt);

    if (gr < grc - dgr) {
        return {0.0, 0.0, 0.0, 0.0};
    }

    const double rnorm = (gr - (grc - dgr)) / (2.0 * dgr);
    const double rn_hk = -grc_hk / (2.0 * dgr);
    const double rn_rt = gr_rt / (2.0 * dgr);

    double rfac = 0.0;
    double rfac_hk = 0.0;
    double rfac_rt = 0.0;
    if (rnorm >= 1.0) {
        rfac = 1.0;
        rfac_hk = 0.0;
        rfac_rt = 0.0;
    } else {
        rfac = 3.0 * rnorm * rnorm - 2.0 * rnorm * rnorm * rnorm;
        const double rfac_rn = 6.0 * rnorm - 6.0 * rnorm * rnorm;
        rfac_hk = rfac_rn * rn_hk;
        rfac_rt = rfac_rn * rn_rt;
    }

    const double arg = 3.87 * hmi - 2.52;
    const double arg_hk = 3.87 * hmi_hk;
    const double ex = std::exp(-(arg * arg));
    const double ex_hk = ex * (-2.0 * arg * arg_hk);

    const double dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
    const double dadr_hk = 0.028 - 0.0345 * ex_hk;

    const double brg = -20.0 * hmi;
    const double af = -0.05 + 2.7 * hmi - 5.5 * hmi * hmi + 3.0 * hmi * hmi * hmi + 0.1 * std::exp(brg);
    const double af_hmi = 2.7 - 11.0 * hmi + 9.0 * hmi * hmi - 2.0 * std::exp(brg);
    const double af_hk = af_hmi * hmi_hk;

    double ax = (af * dadr / th) * rfac;
    double ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac + (af * dadr / th) * rfac_hk;
    double ax_th = -ax / th;
    double ax_rt = (af * dadr / th) * rfac_rt;

    if (hk < hk1) {
        return {ax, ax_hk, ax_th, ax_rt};
    }

    const double hnorm = (hk - hk1) / (hk2 - hk1);
    const double hn_hk = 1.0 / (hk2 - hk1);

    double hfac = 0.0;
    double hf_hk = 0.0;
    if (hnorm >= 1.0) {
        hfac = 1.0;
        hf_hk = 0.0;
    } else {
        hfac = 3.0 * hnorm * hnorm - 2.0 * hnorm * hnorm * hnorm;
        hf_hk = (6.0 * hnorm - 6.0 * hnorm * hnorm) * hn_hk;
    }

    const double ax1 = ax;
    const double ax1_hk = ax_hk;
    const double ax1_th = ax_th;
    const double ax1_rt = ax_rt;

    const double gr0 = 0.30 + 0.35 * std::exp(-0.15 * (hk - 5.0));
    const double gr0_hk = -0.35 * std::exp(-0.15 * (hk - 5.0)) * 0.15;

    const double tnr = std::tanh(1.2 * (gr - gr0));
    const double tnr_rt = (1.0 - tnr * tnr) * 1.2 * gr_rt;
    const double tnr_hk = -(1.0 - tnr * tnr) * 1.2 * gr0_hk;

    double ax2 = (0.086 * tnr - 0.25 / std::pow(hk - 1.0, 1.5)) / th;
    double ax2_hk = (0.086 * tnr_hk + 1.5 * 0.25 / std::pow(hk - 1.0, 2.5)) / th;
    double ax2_rt = (0.086 * tnr_rt) / th;
    double ax2_th = -ax2 / th;

    if (ax2 < 0.0) {
        ax2 = 0.0;
        ax2_hk = 0.0;
        ax2_rt = 0.0;
        ax2_th = 0.0;
    }

    ax = hfac * ax2 + (1.0 - hfac) * ax1;
    ax_hk = hfac * ax2_hk + (1.0 - hfac) * ax1_hk + hf_hk * (ax2 - ax1);
    ax_rt = hfac * ax2_rt + (1.0 - hfac) * ax1_rt;
    ax_th = hfac * ax2_th + (1.0 - hfac) * ax1_th;

    return {ax, ax_hk, ax_th, ax_rt};
}

static std::tuple<double, double, double> hkin_impl(double h, double msq) {
    const double hk = (h - 0.29 * msq) / (1.0 + 0.113 * msq);
    const double hk_h = 1.0 / (1.0 + 0.113 * msq);
    const double hk_msq = (-0.29 - 0.113 * hk) / (1.0 + 0.113 * msq);
    return {hk, hk_h, hk_msq};
}

static std::tuple<double, double, double> dil_impl(double hk, double rt) {
    double di = 0.0;
    double di_hk = 0.0;
    if (hk < 4.0) {
        di = (0.00205 * std::pow(4.0 - hk, 5.5) + 0.207) / rt;
        di_hk = (-0.00205 * 5.5 * std::pow(4.0 - hk, 4.5)) / rt;
    } else {
        const double hkb = hk - 4.0;
        const double den = 1.0 + 0.02 * hkb * hkb;
        di = (-0.0016 * hkb * hkb / den + 0.207) / rt;
        di_hk = (-0.0016 * 2.0 * hkb * (1.0 / den - 0.02 * hkb * hkb / (den * den))) / rt;
    }
    const double di_rt = -di / rt;
    return {di, di_hk, di_rt};
}

static std::tuple<double, double, double> dilw_impl(double hk, double rt) {
    const double msq = 0.0;
    auto hs_out = hsl(hk, rt, msq);
    const double hs = std::get<0>(hs_out);
    const double hs_hk = std::get<1>(hs_out);
    const double hs_rt = std::get<2>(hs_out);

    const double rcd = 1.10 * std::pow(1.0 - 1.0 / hk, 2) / hk;
    const double rcd_hk = -1.10 * (1.0 - 1.0 / hk) * 2.0 / std::pow(hk, 3) - rcd / hk;

    const double di = 2.0 * rcd / (hs * rt);
    const double di_hk = 2.0 * rcd_hk / (hs * rt) - (di / hs) * hs_hk;
    const double di_rt = -di / rt - (di / hs) * hs_rt;
    return {di, di_hk, di_rt};
}

static std::tuple<double, double, double, double> hsl(double hk, double rt, double msq) {
    double hs = 0.0;
    double hs_hk = 0.0;
    if (hk < 4.35) {
        const double tmp = hk - 4.35;
        hs = 0.0111 * tmp * tmp / (hk + 1.0) - 0.0278 * tmp * tmp * tmp / (hk + 1.0) + 1.528
             - 0.0002 * std::pow(tmp * hk, 2);
        hs_hk = 0.0111 * (2.0 * tmp - tmp * tmp / (hk + 1.0)) / (hk + 1.0)
                - 0.0278 * (3.0 * tmp * tmp - tmp * tmp * tmp / (hk + 1.0)) / (hk + 1.0)
                - 0.0002 * 2.0 * tmp * hk * (tmp + hk);
    } else {
        const double hs2 = 0.015;
        hs = hs2 * std::pow(hk - 4.35, 2) / hk + 1.528;
        hs_hk = hs2 * 2.0 * (hk - 4.35) / hk - hs2 * std::pow(hk - 4.35, 2) / (hk * hk);
    }

    const double hs_rt = 0.0;
    const double hs_msq = 0.0;
    return {hs, hs_hk, hs_rt, hs_msq};
}

static std::tuple<double, double, double, double> cfl(double hk, double rt, double msq) {
    double cf = 0.0;
    double cf_hk = 0.0;
    if (hk < 5.5) {
        const double tmp = std::pow(5.5 - hk, 3) / (hk + 1.0);
        cf = (0.0727 * tmp - 0.07) / rt;
        cf_hk = (-0.0727 * tmp * 3.0 / (5.5 - hk) - 0.0727 * tmp / (hk + 1.0)) / rt;
    } else {
        const double tmp = 1.0 - 1.0 / (hk - 4.5);
        cf = (0.015 * tmp * tmp - 0.07) / rt;
        cf_hk = (0.015 * tmp * 2.0 / std::pow(hk - 4.5, 2)) / rt;
    }
    const double cf_rt = -cf / rt;
    const double cf_msq = 0.0;
    return {cf, cf_hk, cf_rt, cf_msq};
}

static std::tuple<double, double, double, double, double> dit(double hs, double us, double cf, double st) {
    const double di = (0.5 * cf * us + st * st * (1.0 - us)) * 2.0 / hs;
    const double di_hs = -(0.5 * cf * us + st * st * (1.0 - us)) * 2.0 / (hs * hs);
    const double di_us = (0.5 * cf - st * st) * 2.0 / hs;
    const double di_cf = (0.5 * us) * 2.0 / hs;
    const double di_st = (2.0 * st * (1.0 - us)) * 2.0 / hs;
    return {di, di_hs, di_us, di_cf, di_st};
}

static std::tuple<double, double, double, double> hst(double hk, double rt, double msq) {
    const double hsmin = 1.500;
    const double dhsinf = 0.015;

    double ho = 0.0;
    double ho_rt = 0.0;
    if (rt > 400.0) {
        ho = 3.0 + 400.0 / rt;
        ho_rt = -400.0 / (rt * rt);
    } else {
        ho = 4.0;
        ho_rt = 0.0;
    }

    double rtz = 0.0;
    double rtz_rt = 0.0;
    if (rt > 200.0) {
        rtz = rt;
        rtz_rt = 1.0;
    } else {
        rtz = 200.0;
        rtz_rt = 0.0;
    }

    double hs = 0.0;
    double hs_hk = 0.0;
    double hs_rt = 0.0;
    if (hk < ho) {
        const double hr = (ho - hk) / (ho - 1.0);
        const double hr_hk = -1.0 / (ho - 1.0);
        const double hr_rt = (1.0 - hr) / (ho - 1.0) * ho_rt;
        hs = (2.0 - hsmin - 4.0 / rtz) * hr * hr * 1.5 / (hk + 0.5) + hsmin + 4.0 / rtz;
        hs_hk = -(2.0 - hsmin - 4.0 / rtz) * hr * hr * 1.5 / std::pow(hk + 0.5, 2)
                + (2.0 - hsmin - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_hk;
        hs_rt = (2.0 - hsmin - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_rt
                + (hr * hr * 1.5 / (hk + 0.5) - 1.0) * 4.0 / (rtz * rtz) * rtz_rt;
    } else {
        const double grt = std::log(rtz);
        const double hdif = hk - ho;
        const double rtmp = hk - ho + 4.0 / grt;
        const double htmp = 0.007 * grt / (rtmp * rtmp) + dhsinf / hk;
        const double htmp_hk = -0.014 * grt / (rtmp * rtmp * rtmp) - dhsinf / (hk * hk);
        const double htmp_rt = -0.014 * grt / (rtmp * rtmp * rtmp) * (-ho_rt - 4.0 / (grt * grt) / rtz * rtz_rt)
                                + 0.007 / (rtmp * rtmp) / rtz * rtz_rt;
        hs = hdif * hdif * htmp + hsmin + 4.0 / rtz;
        hs_hk = hdif * 2.0 * htmp + hdif * hdif * htmp_hk;
        hs_rt = hdif * hdif * htmp_rt - 4.0 / (rtz * rtz) * rtz_rt + hdif * 2.0 * htmp * (-ho_rt);
    }

    const double fm = 1.0 + 0.014 * msq;
    hs = (hs + 0.028 * msq) / fm;
    hs_hk = hs_hk / fm;
    hs_rt = hs_rt / fm;
    const double hs_msq = 0.028 / fm - 0.014 * hs / fm;
    return {hs, hs_hk, hs_rt, hs_msq};
}

static std::tuple<double, double, double, double> cft(double hk, double rt, double msq, double cffac) {
    const double gam = 1.4;
    const double gm1 = gam - 1.0;
    const double fc = std::sqrt(1.0 + 0.5 * gm1 * msq);
    double grt = std::log(rt / fc);
    grt = std::max(grt, 3.0);

    const double gex = -1.74 - 0.31 * hk;

    double arg = -1.33 * hk;
    arg = std::max(-20.0, arg);

    const double thk = std::tanh(4.0 - hk / 0.875);

    const double cfo = cffac * 0.3 * std::exp(arg) * std::pow(grt / 2.3026, gex);
    const double cf = (cfo + 1.1e-4 * (thk - 1.0)) / fc;
    const double cf_hk = (-1.33 * cfo - 0.31 * std::log(grt / 2.3026) * cfo - 1.1e-4 * (1.0 - thk * thk) / 0.875) / fc;
    const double cf_rt = gex * cfo / (fc * grt) / rt;
    const double cf_msq = gex * cfo / (fc * grt) * (-0.25 * gm1 / (fc * fc)) - 0.25 * gm1 * cf / (fc * fc);

    return {cf, cf_hk, cf_rt, cf_msq};
}

static std::tuple<double, double, double> hct(double hk, double msq) {
    const double hc = msq * (0.064 / (hk - 0.8) + 0.251);
    const double hc_hk = msq * (-0.064 / std::pow(hk - 0.8, 2));
    const double hc_msq = 0.064 / (hk - 0.8) + 0.251;
    return {hc, hc_hk, hc_msq};
}
