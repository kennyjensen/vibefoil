# Ported from XFOIL Fortran source (Mark Drela).
# This file is a derived work and remains under the terms of the
# GNU General Public License v2 or later.
# See https://web.mit.edu/drela/Public/web/xfoil/ for the original code and license text.

import math

from .xblcom import NCOM, sync_com_to_vars, sync_vars_to_com


def create2d(rows, cols, fill=0.0):
    return [[fill] * (cols + 1) for _ in range(rows + 1)]


def trchek(bl):
    trchek2(bl)


def axset(hk1, t1, rt1, a1, hk2, t2, rt2, a2, acrit, idampv):
    if idampv == 0:
        ax1, ax1_hk1, ax1_t1, ax1_rt1 = dampl(hk1, t1, rt1)
        ax2, ax2_hk2, ax2_t2, ax2_rt2 = dampl(hk2, t2, rt2)
    else:
        ax1, ax1_hk1, ax1_t1, ax1_rt1 = dampl2(hk1, t1, rt1)
        ax2, ax2_hk2, ax2_t2, ax2_rt2 = dampl2(hk2, t2, rt2)

    axsq = 0.5 * (ax1**2 + ax2**2)
    if axsq <= 0.0:
        axa = 0.0
        axa_ax1 = 0.0
        axa_ax2 = 0.0
    else:
        axa = math.sqrt(axsq)
        axa_ax1 = 0.5 * ax1 / axa
        axa_ax2 = 0.5 * ax2 / axa

    arg = min(20.0 * (acrit - 0.5 * (a1 + a2)), 20.0)
    if arg <= 0.0:
        exn = 1.0
        exn_a1 = 0.0
        exn_a2 = 0.0
    else:
        exn = math.exp(-arg)
        exn_a1 = 20.0 * 0.5 * exn
        exn_a2 = 20.0 * 0.5 * exn

    dax = exn * 0.002 / (t1 + t2)
    dax_a1 = exn_a1 * 0.002 / (t1 + t2)
    dax_a2 = exn_a2 * 0.002 / (t1 + t2)
    dax_t1 = -dax / (t1 + t2)
    dax_t2 = -dax / (t1 + t2)

    ax = axa + dax

    ax_hk1 = axa_ax1 * ax1_hk1
    ax_t1 = axa_ax1 * ax1_t1 + dax_t1
    ax_rt1 = axa_ax1 * ax1_rt1
    ax_a1 = dax_a1

    ax_hk2 = axa_ax2 * ax2_hk2
    ax_t2 = axa_ax2 * ax2_t2 + dax_t2
    ax_rt2 = axa_ax2 * ax2_rt2
    ax_a2 = dax_a2

    return ax, ax_hk1, ax_t1, ax_rt1, ax_a1, ax_hk2, ax_t2, ax_rt2, ax_a2


def trchek2(bl):
    daeps = 5.0e-5

    for icom in range(1, NCOM + 1):
        bl.C2SAV[icom] = bl.COM2[icom]

    ax_out = axset(
        bl.HK1,
        bl.T1,
        bl.RT1,
        bl.AMPL1,
        bl.HK2,
        bl.T2,
        bl.RT2,
        bl.AMPL2,
        bl.AMCRIT,
        bl.IDAMPV,
    )
    ax = ax_out[0]
    bl.AMPL2 = bl.AMPL1 + ax * (bl.X2 - bl.X1)

    converged = False
    last_da2 = 0.0
    last_ax = ax
    for _ in range(30):
        if bl.AMPL2 <= bl.AMCRIT:
            amplt = bl.AMPL2
            amplt_a2 = 1.0
            sfa = 1.0
            sfa_a1 = 0.0
            sfa_a2 = 0.0
        else:
            amplt = bl.AMCRIT
            amplt_a2 = 0.0
            sfa = (amplt - bl.AMPL1) / (bl.AMPL2 - bl.AMPL1)
            sfa_a1 = (sfa - 1.0) / (bl.AMPL2 - bl.AMPL1)
            sfa_a2 = (-sfa) / (bl.AMPL2 - bl.AMPL1)

        xif_eps = 1.0e-12 * max(1.0, abs(bl.X2))
        if bl.XIFORC >= bl.X2 - xif_eps and bl.XIFORC <= bl.X2 + xif_eps:
            xif = bl.X2
        else:
            xif = bl.XIFORC
        if xif < bl.X2:
            sfx = (xif - bl.X1) / (bl.X2 - bl.X1)
            sfx_x1 = (sfx - 1.0) / (bl.X2 - bl.X1)
            sfx_x2 = (-sfx) / (bl.X2 - bl.X1)
            sfx_xf = 1.0 / (bl.X2 - bl.X1)
        else:
            sfx = 1.0
            sfx_x1 = 0.0
            sfx_x2 = 0.0
            sfx_xf = 0.0

        if sfa < sfx:
            wf2 = sfa
            wf2_a1 = sfa_a1
            wf2_a2 = sfa_a2
            wf2_x1 = 0.0
            wf2_x2 = 0.0
            wf2_xf = 0.0
        else:
            wf2 = sfx
            wf2_a1 = 0.0
            wf2_a2 = 0.0
            wf2_x1 = sfx_x1
            wf2_x2 = sfx_x2
            wf2_xf = sfx_xf

        wf1 = 1.0 - wf2
        wf1_a1 = -wf2_a1
        wf1_a2 = -wf2_a2
        wf1_x1 = -wf2_x1
        wf1_x2 = -wf2_x2
        wf1_xf = -wf2_xf

        xt = bl.X1 * wf1 + bl.X2 * wf2
        tt = bl.T1 * wf1 + bl.T2 * wf2
        dt = bl.D1 * wf1 + bl.D2 * wf2
        ut = bl.U1 * wf1 + bl.U2 * wf2

        xt_a2 = bl.X1 * wf1_a2 + bl.X2 * wf2_a2
        tt_a2 = bl.T1 * wf1_a2 + bl.T2 * wf2_a2
        dt_a2 = bl.D1 * wf1_a2 + bl.D2 * wf2_a2
        ut_a2 = bl.U1 * wf1_a2 + bl.U2 * wf2_a2

        bl.X2 = xt
        bl.T2 = tt
        bl.D2 = dt
        bl.U2 = ut

        blkin(bl)

        hkt = bl.HK2
        hkt_tt = bl.HK2_T2
        hkt_dt = bl.HK2_D2
        hkt_ut = bl.HK2_U2
        hkt_ms = bl.HK2_MS

        rtt = bl.RT2
        rtt_tt = bl.RT2_T2
        rtt_ut = bl.RT2_U2
        rtt_ms = bl.RT2_MS
        rtt_re = bl.RT2_RE

        amsave = bl.AMPL2
        for icom in range(1, NCOM + 1):
            bl.COM2[icom] = bl.C2SAV[icom]
        sync_com_to_vars(bl, 2)
        bl.AMPL2 = amsave

        ax_out = axset(
            bl.HK1,
            bl.T1,
            bl.RT1,
            bl.AMPL1,
            hkt,
            tt,
            rtt,
            amplt,
            bl.AMCRIT,
            bl.IDAMPV,
        )
        ax = ax_out[0]
        last_ax = ax

        if ax <= 0.0:
            converged = True
            break

        ax_a2 = (
            (ax_out[5] * hkt_tt + ax_out[6] + ax_out[7] * rtt_tt) * tt_a2
            + (ax_out[5] * hkt_dt) * dt_a2
            + (ax_out[5] * hkt_ut + ax_out[7] * rtt_ut) * ut_a2
            + ax_out[8] * amplt_a2
        )

        res = bl.AMPL2 - bl.AMPL1 - ax * (bl.X2 - bl.X1)
        res_a2 = 1.0 - ax_a2 * (bl.X2 - bl.X1)
        da2 = -res / res_a2
        last_da2 = da2

        rlx = 1.0
        dxt = xt_a2 * da2
        if rlx * abs(dxt / (bl.X2 - bl.X1)) > 0.05:
            rlx = 0.05 * abs((bl.X2 - bl.X1) / dxt)
        if rlx * abs(da2) > 1.0:
            rlx = 1.0 * abs(1.0 / da2)

        if abs(da2) < daeps:
            converged = True
            break

        if (bl.AMPL2 > bl.AMCRIT and bl.AMPL2 + rlx * da2 < bl.AMCRIT) or (
            bl.AMPL2 < bl.AMCRIT and bl.AMPL2 + rlx * da2 > bl.AMCRIT
        ):
            bl.AMPL2 = bl.AMCRIT
        else:
            bl.AMPL2 = bl.AMPL2 + rlx * da2

    if not converged:
        print("TRCHEK2: N2 convergence failed.")
        print(
            f" x:{bl.X1:9.5f}{bl.XT:9.5f}{bl.X2:9.5f}  N:{bl.AMPL1:7.3f}{amplt:7.3f}{bl.AMPL2:7.3f}  Nx:{last_ax:8.3f}   dN:{last_da2:10.3E}"
        )

    bl.XT = xt
    bl.XT_A2 = xt_a2

    bl.TRFREE = bl.AMPL2 >= bl.AMCRIT
    xif_eps_final = 1.0e-12 * max(1.0, abs(bl.X2))
    if bl.XIFORC >= bl.X2 - xif_eps_final and bl.XIFORC <= bl.X2 + xif_eps_final:
        xif_final = bl.X2
    else:
        xif_final = bl.XIFORC
    bl.TRFORC = xif_final > bl.X1 and xif_final <= bl.X2
    bl.TRAN = bl.TRFORC or bl.TRFREE
    if not bl.TRAN:
        return

    if bl.TRFREE and bl.TRFORC:
        bl.TRFORC = xif_final < xt
        bl.TRFREE = xif_final >= xt

    if bl.TRFORC:
        bl.XT = xif_final
        bl.XT_A1 = 0.0
        bl.XT_X1 = 0.0
        bl.XT_T1 = 0.0
        bl.XT_D1 = 0.0
        bl.XT_U1 = 0.0
        bl.XT_X2 = 0.0
        bl.XT_T2 = 0.0
        bl.XT_D2 = 0.0
        bl.XT_U2 = 0.0
        bl.XT_MS = 0.0
        bl.XT_RE = 0.0
        bl.XT_XF = 1.0
        return

    bl.XT = bl.X1 * wf1 + bl.X2 * wf2
    bl.XT_A2 = xt_a2
    tt = bl.T1 * wf1 + bl.T2 * wf2
    dt = bl.D1 * wf1 + bl.D2 * wf2
    ut = bl.U1 * wf1 + bl.U2 * wf2

    bl.XT_X1 = wf1
    tt_t1 = wf1
    dt_d1 = wf1
    ut_u1 = wf1

    bl.XT_X2 = wf2
    tt_t2 = wf2
    dt_d2 = wf2
    ut_u2 = wf2

    bl.XT_A1 = bl.X1 * wf1_a1 + bl.X2 * wf2_a1
    tt_a1 = bl.T1 * wf1_a1 + bl.T2 * wf2_a1
    dt_a1 = bl.D1 * wf1_a1 + bl.D2 * wf2_a1
    ut_a1 = bl.U1 * wf1_a1 + bl.U2 * wf2_a1

    bl.XT_X1 = bl.X1 * wf1_x1 + bl.X2 * wf2_x1 + bl.XT_X1
    tt_x1 = bl.T1 * wf1_x1 + bl.T2 * wf2_x1
    dt_x1 = bl.D1 * wf1_x1 + bl.D2 * wf2_x1
    ut_x1 = bl.U1 * wf1_x1 + bl.U2 * wf2_x1

    bl.XT_X2 = bl.X1 * wf1_x2 + bl.X2 * wf2_x2 + bl.XT_X2
    tt_x2 = bl.T1 * wf1_x2 + bl.T2 * wf2_x2
    dt_x2 = bl.D1 * wf1_x2 + bl.D2 * wf2_x2
    ut_x2 = bl.U1 * wf1_x2 + bl.U2 * wf2_x2

    bl.XT_XF = bl.X1 * wf1_xf + bl.X2 * wf2_xf
    tt_xf = bl.T1 * wf1_xf + bl.T2 * wf2_xf
    dt_xf = bl.D1 * wf1_xf + bl.D2 * wf2_xf
    ut_xf = bl.U1 * wf1_xf + bl.U2 * wf2_xf

    ax_t1 = ax_out[1] * bl.HK1_T1 + ax_out[2] + ax_out[3] * bl.RT1_T1 + (
        ax_out[5] * hkt_tt + ax_out[6] + ax_out[7] * rtt_tt
    ) * tt_t1
    ax_d1 = ax_out[1] * bl.HK1_D1 + (ax_out[5] * hkt_dt) * dt_d1
    ax_u1 = ax_out[1] * bl.HK1_U1 + ax_out[3] * bl.RT1_U1 + (
        ax_out[5] * hkt_ut + ax_out[7] * rtt_ut
    ) * ut_u1
    ax_a1 = ax_out[4] + (ax_out[5] * hkt_tt + ax_out[6] + ax_out[7] * rtt_tt) * tt_a1 + (
        ax_out[5] * hkt_dt
    ) * dt_a1 + (ax_out[5] * hkt_ut + ax_out[7] * rtt_ut) * ut_a1
    ax_x1 = (ax_out[5] * hkt_tt + ax_out[6] + ax_out[7] * rtt_tt) * tt_x1 + (
        ax_out[5] * hkt_dt
    ) * dt_x1 + (ax_out[5] * hkt_ut + ax_out[7] * rtt_ut) * ut_x1

    ax_t2 = (ax_out[5] * hkt_tt + ax_out[6] + ax_out[7] * rtt_tt) * tt_t2
    ax_d2 = (ax_out[5] * hkt_dt) * dt_d2
    ax_u2 = (ax_out[5] * hkt_ut + ax_out[7] * rtt_ut) * ut_u2
    ax_a2 = ax_out[8] * amplt_a2 + (ax_out[5] * hkt_tt + ax_out[6] + ax_out[7] * rtt_tt) * tt_a2 + (
        ax_out[5] * hkt_dt
    ) * dt_a2 + (ax_out[5] * hkt_ut + ax_out[7] * rtt_ut) * ut_a2
    ax_x2 = (ax_out[5] * hkt_tt + ax_out[6] + ax_out[7] * rtt_tt) * tt_x2 + (
        ax_out[5] * hkt_dt
    ) * dt_x2 + (ax_out[5] * hkt_ut + ax_out[7] * rtt_ut) * ut_x2

    ax_xf = (ax_out[5] * hkt_tt + ax_out[6] + ax_out[7] * rtt_tt) * tt_xf + (
        ax_out[5] * hkt_dt
    ) * dt_xf + (ax_out[5] * hkt_ut + ax_out[7] * rtt_ut) * ut_xf

    ax_ms = ax_out[5] * hkt_ms + ax_out[7] * rtt_ms + ax_out[1] * bl.HK1_MS + ax_out[3] * bl.RT1_MS
    ax_re = ax_out[7] * rtt_re + ax_out[3] * bl.RT1_RE

    z_ax = -(bl.X2 - bl.X1)
    z_a1 = z_ax * ax_a1 - 1.0
    z_t1 = z_ax * ax_t1
    z_d1 = z_ax * ax_d1
    z_u1 = z_ax * ax_u1
    z_x1 = z_ax * ax_x1 + ax
    z_a2 = z_ax * ax_a2 + 1.0
    z_t2 = z_ax * ax_t2
    z_d2 = z_ax * ax_d2
    z_u2 = z_ax * ax_u2
    z_x2 = z_ax * ax_x2 - ax
    z_xf = z_ax * ax_xf
    z_ms = z_ax * ax_ms
    z_re = z_ax * ax_re

    bl.XT_A1 = bl.XT_A1 - (bl.XT_A2 / z_a2) * z_a1
    bl.XT_T1 = -(bl.XT_A2 / z_a2) * z_t1
    bl.XT_D1 = -(bl.XT_A2 / z_a2) * z_d1
    bl.XT_U1 = -(bl.XT_A2 / z_a2) * z_u1
    bl.XT_X1 = bl.XT_X1 - (bl.XT_A2 / z_a2) * z_x1
    bl.XT_T2 = -(bl.XT_A2 / z_a2) * z_t2
    bl.XT_D2 = -(bl.XT_A2 / z_a2) * z_d2
    bl.XT_U2 = -(bl.XT_A2 / z_a2) * z_u2
    bl.XT_X2 = bl.XT_X2 - (bl.XT_A2 / z_a2) * z_x2
    bl.XT_MS = -(bl.XT_A2 / z_a2) * z_ms
    bl.XT_RE = -(bl.XT_A2 / z_a2) * z_re
    bl.XT_XF = 0.0


def blsys(bl):
    if bl.WAKE:
        blvar(bl, 3)
        blmid(bl, 3)
    elif bl.TURB or bl.TRAN:
        blvar(bl, 2)
        blmid(bl, 2)
    else:
        blvar(bl, 1)
        blmid(bl, 1)

    if bl.SIMI:
        for icom in range(1, NCOM + 1):
            bl.COM1[icom] = bl.COM2[icom]
        sync_com_to_vars(bl, 1)

    if bl.TRAN:
        trdif(bl)
    elif bl.SIMI:
        bldif(bl, 0)
    elif not bl.TURB:
        bldif(bl, 1)
    elif bl.WAKE:
        bldif(bl, 3)
    else:
        bldif(bl, 2)

    if bl.SIMI:
        for k in range(1, 5):
            for l in range(1, 6):
                bl.VS2[k][l] = bl.VS1[k][l] + bl.VS2[k][l]
                bl.VS1[k][l] = 0.0

    for k in range(1, 5):
        res_u1 = bl.VS1[k][4]
        res_u2 = bl.VS2[k][4]
        res_ms = bl.VSM[k]
        bl.VS1[k][4] = res_u1 * bl.U1_UEI
        bl.VS2[k][4] = res_u2 * bl.U2_UEI
        bl.VSM[k] = res_u1 * bl.U1_MS + res_u2 * bl.U2_MS + res_ms


def tesys(bl, cte, tte, dte):
    for k in range(1, 5):
        bl.VSREZ[k] = 0.0
        bl.VSM[k] = 0.0
        bl.VSR[k] = 0.0
        bl.VSX[k] = 0.0
        for l in range(1, 6):
            bl.VS1[k][l] = 0.0
            bl.VS2[k][l] = 0.0

    blvar(bl, 3)

    bl.VS1[1][1] = -1.0
    bl.VS2[1][1] = 1.0
    bl.VSREZ[1] = cte - bl.S2

    bl.VS1[2][2] = -1.0
    bl.VS2[2][2] = 1.0
    bl.VSREZ[2] = tte - bl.T2

    bl.VS1[3][3] = -1.0
    bl.VS2[3][3] = 1.0
    bl.VSREZ[3] = dte - bl.D2 - bl.DW2


def blprv(bl, xsi, ami, cti, thi, dsi, dswaki, uei):
    bl.X2 = xsi
    bl.AMPL2 = ami
    bl.S2 = cti
    bl.T2 = thi
    bl.D2 = dsi - dswaki
    bl.DW2 = dswaki

    bl.U2 = uei * (1.0 - bl.TKBL) / (1.0 - bl.TKBL * (uei / bl.QINFBL) ** 2)
    bl.U2_UEI = (1.0 + bl.TKBL * (2.0 * bl.U2 * uei / bl.QINFBL**2 - 1.0)) / (
        1.0 - bl.TKBL * (uei / bl.QINFBL) ** 2
    )
    bl.U2_MS = (bl.U2 * (uei / bl.QINFBL) ** 2 - uei) * bl.TKBL_MS / (
        1.0 - bl.TKBL * (uei / bl.QINFBL) ** 2
    )
    sync_vars_to_com(bl, 2)


def blkin(bl):
    bl.M2 = bl.U2 * bl.U2 * bl.HSTINV / (bl.GM1BL * (1.0 - 0.5 * bl.U2 * bl.U2 * bl.HSTINV))
    tr2 = 1.0 + 0.5 * bl.GM1BL * bl.M2
    bl.M2_U2 = 2.0 * bl.M2 * tr2 / bl.U2
    bl.M2_MS = (
        bl.U2 * bl.U2 * tr2 / (bl.GM1BL * (1.0 - 0.5 * bl.U2 * bl.U2 * bl.HSTINV)) * bl.HSTINV_MS
    )

    bl.R2 = bl.RSTBL * tr2 ** (-1.0 / bl.GM1BL)
    bl.R2_U2 = -bl.R2 / tr2 * 0.5 * bl.M2_U2
    bl.R2_MS = -bl.R2 / tr2 * 0.5 * bl.M2_MS + bl.RSTBL_MS * tr2 ** (-1.0 / bl.GM1BL)

    bl.H2 = bl.D2 / bl.T2
    bl.H2_D2 = 1.0 / bl.T2
    bl.H2_T2 = -bl.H2 / bl.T2

    herat = 1.0 - 0.5 * bl.U2 * bl.U2 * bl.HSTINV
    he_u2 = -bl.U2 * bl.HSTINV
    he_ms = -0.5 * bl.U2 * bl.U2 * bl.HSTINV_MS

    bl.V2 = math.sqrt((herat) ** 3) * (1.0 + bl.HVRAT) / (herat + bl.HVRAT) / bl.REYBL
    v2_he = bl.V2 * (1.5 / herat - 1.0 / (herat + bl.HVRAT))

    bl.V2_U2 = v2_he * he_u2
    bl.V2_MS = -bl.V2 / bl.REYBL * bl.REYBL_MS + v2_he * he_ms
    bl.V2_RE = -bl.V2 / bl.REYBL * bl.REYBL_RE

    bl.HK2, hk2_h2, hk2_m2 = hkin(bl.H2, bl.M2)

    bl.HK2_U2 = hk2_m2 * bl.M2_U2
    bl.HK2_T2 = hk2_h2 * bl.H2_T2
    bl.HK2_D2 = hk2_h2 * bl.H2_D2
    bl.HK2_MS = hk2_m2 * bl.M2_MS

    bl.RT2 = bl.R2 * bl.U2 * bl.T2 / bl.V2
    bl.RT2_U2 = bl.RT2 * (1.0 / bl.U2 + bl.R2_U2 / bl.R2 - bl.V2_U2 / bl.V2)
    bl.RT2_T2 = bl.RT2 / bl.T2
    bl.RT2_MS = bl.RT2 * (bl.R2_MS / bl.R2 - bl.V2_MS / bl.V2)
    bl.RT2_RE = bl.RT2 * (-bl.V2_RE / bl.V2)
    sync_vars_to_com(bl, 2)


def blvar(bl, ityp):
    if ityp == 3:
        bl.HK2 = max(bl.HK2, 1.00005)
    if ityp != 3:
        bl.HK2 = max(bl.HK2, 1.05000)

    bl.HC2, hc2_hk2, hc2_m2 = hct(bl.HK2, bl.M2)
    bl.HC2_U2 = hc2_hk2 * bl.HK2_U2 + hc2_m2 * bl.M2_U2
    bl.HC2_T2 = hc2_hk2 * bl.HK2_T2
    bl.HC2_D2 = hc2_hk2 * bl.HK2_D2
    bl.HC2_MS = hc2_hk2 * bl.HK2_MS + hc2_m2 * bl.M2_MS

    if ityp == 1:
        bl.HS2, hs2_hk2, hs2_rt2, hs2_m2 = hsl(bl.HK2, bl.RT2, bl.M2)
    else:
        bl.HS2, hs2_hk2, hs2_rt2, hs2_m2 = hst(bl.HK2, bl.RT2, bl.M2)

    bl.HS2_U2 = hs2_hk2 * bl.HK2_U2 + hs2_rt2 * bl.RT2_U2 + hs2_m2 * bl.M2_U2
    bl.HS2_T2 = hs2_hk2 * bl.HK2_T2 + hs2_rt2 * bl.RT2_T2
    bl.HS2_D2 = hs2_hk2 * bl.HK2_D2
    bl.HS2_MS = hs2_hk2 * bl.HK2_MS + hs2_rt2 * bl.RT2_MS + hs2_m2 * bl.M2_MS
    bl.HS2_RE = hs2_rt2 * bl.RT2_RE

    bl.US2 = 0.5 * bl.HS2 * (1.0 - (bl.HK2 - 1.0) / (bl.GBCON * bl.H2))
    us2_hs2 = 0.5 * (1.0 - (bl.HK2 - 1.0) / (bl.GBCON * bl.H2))
    us2_hk2 = 0.5 * bl.HS2 * (-1.0 / (bl.GBCON * bl.H2))
    us2_h2 = 0.5 * bl.HS2 * (bl.HK2 - 1.0) / (bl.GBCON * bl.H2**2)

    bl.US2_U2 = us2_hs2 * bl.HS2_U2 + us2_hk2 * bl.HK2_U2
    bl.US2_T2 = us2_hs2 * bl.HS2_T2 + us2_hk2 * bl.HK2_T2 + us2_h2 * bl.H2_T2
    bl.US2_D2 = us2_hs2 * bl.HS2_D2 + us2_hk2 * bl.HK2_D2 + us2_h2 * bl.H2_D2
    bl.US2_MS = us2_hs2 * bl.HS2_MS + us2_hk2 * bl.HK2_MS
    bl.US2_RE = us2_hs2 * bl.HS2_RE

    if ityp <= 2 and bl.US2 > 0.95:
        bl.US2 = 0.98
        bl.US2_U2 = 0.0
        bl.US2_T2 = 0.0
        bl.US2_D2 = 0.0
        bl.US2_MS = 0.0
        bl.US2_RE = 0.0

    if ityp == 3 and bl.US2 > 0.99995:
        bl.US2 = 0.99995
        bl.US2_U2 = 0.0
        bl.US2_T2 = 0.0
        bl.US2_D2 = 0.0
        bl.US2_MS = 0.0
        bl.US2_RE = 0.0

    if ityp == 2:
        gcc = bl.GCCON
        hkc = bl.HK2 - 1.0 - gcc / bl.RT2
        hkc_hk2 = 1.0
        hkc_rt2 = gcc / bl.RT2**2
        if hkc < 0.01:
            hkc = 0.01
            hkc_hk2 = 0.0
            hkc_rt2 = 0.0
    else:
        hkc = bl.HK2 - 1.0
        hkc_hk2 = 1.0
        hkc_rt2 = 0.0

    hkb = bl.HK2 - 1.0
    usb = 1.0 - bl.US2
    bl.CQ2 = math.sqrt(bl.CTCON * bl.HS2 * hkb * hkc**2 / (usb * bl.H2 * bl.HK2**2))
    cq2_hs2 = bl.CTCON * hkb * hkc**2 / (usb * bl.H2 * bl.HK2**2) * 0.5 / bl.CQ2
    cq2_us2 = bl.CTCON * bl.HS2 * hkb * hkc**2 / (usb * bl.H2 * bl.HK2**2) / usb * 0.5 / bl.CQ2
    cq2_hk2 = (
        bl.CTCON * bl.HS2 * hkc**2 / (usb * bl.H2 * bl.HK2**2) * 0.5 / bl.CQ2
        - bl.CTCON * bl.HS2 * hkb * hkc**2 / (usb * bl.H2 * bl.HK2**3) * 2.0 * 0.5 / bl.CQ2
        + bl.CTCON * bl.HS2 * hkb * hkc / (usb * bl.H2 * bl.HK2**2) * 2.0 * 0.5 / bl.CQ2 * hkc_hk2
    )
    cq2_h2 = -bl.CTCON * bl.HS2 * hkb * hkc**2 / (usb * bl.H2 * bl.HK2**2) / bl.H2 * 0.5 / bl.CQ2
    cq2_rt2 = bl.CTCON * bl.HS2 * hkb * hkc / (usb * bl.H2 * bl.HK2**2) * 2.0 * 0.5 / bl.CQ2 * hkc_rt2

    bl.CQ2_U2 = cq2_hs2 * bl.HS2_U2 + cq2_us2 * bl.US2_U2 + cq2_hk2 * bl.HK2_U2
    bl.CQ2_T2 = cq2_hs2 * bl.HS2_T2 + cq2_us2 * bl.US2_T2 + cq2_hk2 * bl.HK2_T2 + cq2_h2 * bl.H2_T2
    bl.CQ2_D2 = cq2_hs2 * bl.HS2_D2 + cq2_us2 * bl.US2_D2 + cq2_hk2 * bl.HK2_D2 + cq2_h2 * bl.H2_D2
    bl.CQ2_MS = cq2_hs2 * bl.HS2_MS + cq2_us2 * bl.US2_MS + cq2_hk2 * bl.HK2_MS
    bl.CQ2_RE = cq2_hs2 * bl.HS2_RE + cq2_us2 * bl.US2_RE + cq2_rt2 * bl.RT2_RE

    bl.CQ2_U2 = bl.CQ2_U2 + cq2_rt2 * bl.RT2_U2
    bl.CQ2_T2 = bl.CQ2_T2 + cq2_rt2 * bl.RT2_T2
    bl.CQ2_MS = bl.CQ2_MS + cq2_rt2 * bl.RT2_MS

    if ityp == 3:
        bl.CF2 = 0.0
        bl.CF2_HK2 = 0.0
        bl.CF2_RT2 = 0.0
        bl.CF2_M2 = 0.0
    elif ityp == 1:
        bl.CF2, bl.CF2_HK2, bl.CF2_RT2, bl.CF2_M2 = cfl(bl.HK2, bl.RT2, bl.M2)
    else:
        bl.CF2, bl.CF2_HK2, bl.CF2_RT2, bl.CF2_M2 = cft(bl.HK2, bl.RT2, bl.M2, bl.CFFAC)
        cfl_cf, cfl_hk, cfl_rt, cfl_ms = cfl(bl.HK2, bl.RT2, bl.M2)
        if cfl_cf > bl.CF2:
            bl.CF2 = cfl_cf
            bl.CF2_HK2 = cfl_hk
            bl.CF2_RT2 = cfl_rt
            bl.CF2_M2 = cfl_ms

    bl.CF2_U2 = bl.CF2_HK2 * bl.HK2_U2 + bl.CF2_RT2 * bl.RT2_U2 + bl.CF2_M2 * bl.M2_U2
    bl.CF2_T2 = bl.CF2_HK2 * bl.HK2_T2 + bl.CF2_RT2 * bl.RT2_T2
    bl.CF2_D2 = bl.CF2_HK2 * bl.HK2_D2
    bl.CF2_MS = bl.CF2_HK2 * bl.HK2_MS + bl.CF2_RT2 * bl.RT2_MS + bl.CF2_M2 * bl.M2_MS
    bl.CF2_RE = bl.CF2_RT2 * bl.RT2_RE

    if ityp == 1:
        bl.DI2, di2_hk2, di2_rt2 = dil(bl.HK2, bl.RT2)
        bl.DI2_U2 = di2_hk2 * bl.HK2_U2 + di2_rt2 * bl.RT2_U2
        bl.DI2_T2 = di2_hk2 * bl.HK2_T2 + di2_rt2 * bl.RT2_T2
        bl.DI2_D2 = di2_hk2 * bl.HK2_D2
        bl.DI2_S2 = 0.0
        bl.DI2_MS = di2_hk2 * bl.HK2_MS + di2_rt2 * bl.RT2_MS
        bl.DI2_RE = di2_rt2 * bl.RT2_RE
    elif ityp == 2:
        cf2t, cf2t_hk, cf2t_rt, cf2t_ms = cft(bl.HK2, bl.RT2, bl.M2, bl.CFFAC)
        cf2t_u2 = cf2t_hk * bl.HK2_U2 + cf2t_rt * bl.RT2_U2 + cf2t_ms * bl.M2_U2
        cf2t_t2 = cf2t_hk * bl.HK2_T2 + cf2t_rt * bl.RT2_T2
        cf2t_d2 = cf2t_hk * bl.HK2_D2
        cf2t_ms2 = cf2t_hk * bl.HK2_MS + cf2t_rt * bl.RT2_MS + cf2t_ms * bl.M2_MS
        cf2t_re = cf2t_rt * bl.RT2_RE

        bl.DI2 = (0.5 * cf2t * bl.US2) * 2.0 / bl.HS2
        di2_hs2 = -(0.5 * cf2t * bl.US2) * 2.0 / bl.HS2**2
        di2_us2 = (0.5 * cf2t) * 2.0 / bl.HS2
        di2_cf2t = (0.5 * bl.US2) * 2.0 / bl.HS2

        bl.DI2_S2 = 0.0
        bl.DI2_U2 = di2_hs2 * bl.HS2_U2 + di2_us2 * bl.US2_U2 + di2_cf2t * cf2t_u2
        bl.DI2_T2 = di2_hs2 * bl.HS2_T2 + di2_us2 * bl.US2_T2 + di2_cf2t * cf2t_t2
        bl.DI2_D2 = di2_hs2 * bl.HS2_D2 + di2_us2 * bl.US2_D2 + di2_cf2t * cf2t_d2
        bl.DI2_MS = di2_hs2 * bl.HS2_MS + di2_us2 * bl.US2_MS + di2_cf2t * cf2t_ms2
        bl.DI2_RE = di2_hs2 * bl.HS2_RE + di2_us2 * bl.US2_RE + di2_cf2t * cf2t_re

        grt = math.log(bl.RT2)
        hmin = 1.0 + 2.1 / grt
        hm_rt2 = -(2.1 / grt**2) / bl.RT2

        fl = (bl.HK2 - 1.0) / (hmin - 1.0)
        fl_hk2 = 1.0 / (hmin - 1.0)
        fl_rt2 = (-fl / (hmin - 1.0)) * hm_rt2

        tfl = math.tanh(fl)
        dfac = 0.5 + 0.5 * tfl
        df_fl = 0.5 * (1.0 - tfl**2)

        df_hk2 = df_fl * fl_hk2
        df_rt2 = df_fl * fl_rt2

        bl.DI2_S2 = bl.DI2_S2 * dfac
        bl.DI2_U2 = bl.DI2_U2 * dfac + bl.DI2 * (df_hk2 * bl.HK2_U2 + df_rt2 * bl.RT2_U2)
        bl.DI2_T2 = bl.DI2_T2 * dfac + bl.DI2 * (df_hk2 * bl.HK2_T2 + df_rt2 * bl.RT2_T2)
        bl.DI2_D2 = bl.DI2_D2 * dfac + bl.DI2 * (df_hk2 * bl.HK2_D2)
        bl.DI2_MS = bl.DI2_MS * dfac + bl.DI2 * (df_hk2 * bl.HK2_MS + df_rt2 * bl.RT2_MS)
        bl.DI2_RE = bl.DI2_RE * dfac + bl.DI2 * (df_rt2 * bl.RT2_RE)
        bl.DI2 = bl.DI2 * dfac
    else:
        bl.DI2 = 0.0
        bl.DI2_S2 = 0.0
        bl.DI2_U2 = 0.0
        bl.DI2_T2 = 0.0
        bl.DI2_D2 = 0.0
        bl.DI2_MS = 0.0
        bl.DI2_RE = 0.0

    if ityp != 1:
        dd = bl.S2**2 * (0.995 - bl.US2) * 2.0 / bl.HS2
        dd_hs2 = -(bl.S2**2) * (0.995 - bl.US2) * 2.0 / bl.HS2**2
        dd_us2 = -(bl.S2**2) * 2.0 / bl.HS2
        dd_s2 = bl.S2 * 2.0 * (0.995 - bl.US2) * 2.0 / bl.HS2

        bl.DI2 += dd
        bl.DI2_S2 = dd_s2
        bl.DI2_U2 += dd_hs2 * bl.HS2_U2 + dd_us2 * bl.US2_U2
        bl.DI2_T2 += dd_hs2 * bl.HS2_T2 + dd_us2 * bl.US2_T2
        bl.DI2_D2 += dd_hs2 * bl.HS2_D2 + dd_us2 * bl.US2_D2
        bl.DI2_MS += dd_hs2 * bl.HS2_MS + dd_us2 * bl.US2_MS
        bl.DI2_RE += dd_hs2 * bl.HS2_RE + dd_us2 * bl.US2_RE

        dd = 0.15 * (0.995 - bl.US2) ** 2 / bl.RT2 * 2.0 / bl.HS2
        dd_us2b = -0.15 * (0.995 - bl.US2) * 2.0 / bl.RT2 * 2.0 / bl.HS2
        dd_hs2b = -dd / bl.HS2
        dd_rt2 = -dd / bl.RT2

        bl.DI2 += dd
        bl.DI2_U2 += dd_hs2b * bl.HS2_U2 + dd_us2b * bl.US2_U2 + dd_rt2 * bl.RT2_U2
        bl.DI2_T2 += dd_hs2b * bl.HS2_T2 + dd_us2b * bl.US2_T2 + dd_rt2 * bl.RT2_T2
        bl.DI2_D2 += dd_hs2b * bl.HS2_D2 + dd_us2b * bl.US2_D2
        bl.DI2_MS += dd_hs2b * bl.HS2_MS + dd_us2b * bl.US2_MS + dd_rt2 * bl.RT2_MS
        bl.DI2_RE += dd_hs2b * bl.HS2_RE + dd_us2b * bl.US2_RE + dd_rt2 * bl.RT2_RE

    if ityp == 2:
        di2l, di2l_hk, di2l_rt = dil(bl.HK2, bl.RT2)
        if di2l > bl.DI2:
            bl.DI2 = di2l
            bl.DI2_S2 = 0.0
            bl.DI2_U2 = di2l_hk * bl.HK2_U2 + di2l_rt * bl.RT2_U2
            bl.DI2_T2 = di2l_hk * bl.HK2_T2 + di2l_rt * bl.RT2_T2
            bl.DI2_D2 = di2l_hk * bl.HK2_D2
            bl.DI2_MS = di2l_hk * bl.HK2_MS + di2l_rt * bl.RT2_MS
            bl.DI2_RE = di2l_rt * bl.RT2_RE

    if ityp == 3:
        di2l, di2l_hk, di2l_rt = dilw(bl.HK2, bl.RT2)
        if di2l > bl.DI2:
            bl.DI2 = di2l
            bl.DI2_S2 = 0.0
            bl.DI2_U2 = di2l_hk * bl.HK2_U2 + di2l_rt * bl.RT2_U2
            bl.DI2_T2 = di2l_hk * bl.HK2_T2 + di2l_rt * bl.RT2_T2
            bl.DI2_D2 = di2l_hk * bl.HK2_D2
            bl.DI2_MS = di2l_hk * bl.HK2_MS + di2l_rt * bl.RT2_MS
            bl.DI2_RE = di2l_rt * bl.RT2_RE

    if ityp == 3:
        bl.DI2 *= 2.0
        bl.DI2_S2 *= 2.0
        bl.DI2_U2 *= 2.0
        bl.DI2_T2 *= 2.0
        bl.DI2_D2 *= 2.0
        bl.DI2_MS *= 2.0
        bl.DI2_RE *= 2.0

    bl.DE2 = (3.15 + 1.72 / (bl.HK2 - 1.0)) * bl.T2 + bl.D2
    de2_hk2 = (-1.72 / (bl.HK2 - 1.0) ** 2) * bl.T2

    bl.DE2_U2 = de2_hk2 * bl.HK2_U2
    bl.DE2_T2 = de2_hk2 * bl.HK2_T2 + (3.15 + 1.72 / (bl.HK2 - 1.0))
    bl.DE2_D2 = de2_hk2 * bl.HK2_D2 + 1.0
    bl.DE2_MS = de2_hk2 * bl.HK2_MS

    hdmax = 12.0
    if bl.DE2 > hdmax * bl.T2:
        bl.DE2 = hdmax * bl.T2
        bl.DE2_U2 = 0.0
        bl.DE2_T2 = hdmax
        bl.DE2_D2 = 0.0
        bl.DE2_MS = 0.0
    sync_vars_to_com(bl, 2)


def blmid(bl, ityp):
    if bl.SIMI:
        bl.HK1 = bl.HK2
        bl.HK1_T1 = bl.HK2_T2
        bl.HK1_D1 = bl.HK2_D2
        bl.HK1_U1 = bl.HK2_U2
        bl.HK1_MS = bl.HK2_MS
        bl.RT1 = bl.RT2
        bl.RT1_T1 = bl.RT2_T2
        bl.RT1_U1 = bl.RT2_U2
        bl.RT1_MS = bl.RT2_MS
        bl.RT1_RE = bl.RT2_RE
        bl.M1 = bl.M2
        bl.M1_U1 = bl.M2_U2
        bl.M1_MS = bl.M2_MS
        sync_vars_to_com(bl, 1)

    hka = 0.5 * (bl.HK1 + bl.HK2)
    rta = 0.5 * (bl.RT1 + bl.RT2)
    ma = 0.5 * (bl.M1 + bl.M2)

    if ityp == 3:
        bl.CFM = 0.0
        bl.CFM_HKA = 0.0
        bl.CFM_RTA = 0.0
        bl.CFM_MA = 0.0
        bl.CFM_MS = 0.0
    elif ityp == 1:
        bl.CFM, bl.CFM_HKA, bl.CFM_RTA, bl.CFM_MA = cfl(hka, rta, ma)
    else:
        bl.CFM, bl.CFM_HKA, bl.CFM_RTA, bl.CFM_MA = cft(hka, rta, ma, bl.CFFAC)
        cfml, cfml_hka, cfml_rta, cfml_ma = cfl(hka, rta, ma)
        if cfml > bl.CFM:
            bl.CFM = cfml
            bl.CFM_HKA = cfml_hka
            bl.CFM_RTA = cfml_rta
            bl.CFM_MA = cfml_ma

    bl.CFM_U1 = 0.5 * (bl.CFM_HKA * bl.HK1_U1 + bl.CFM_MA * bl.M1_U1 + bl.CFM_RTA * bl.RT1_U1)
    bl.CFM_T1 = 0.5 * (bl.CFM_HKA * bl.HK1_T1 + bl.CFM_RTA * bl.RT1_T1)
    bl.CFM_D1 = 0.5 * (bl.CFM_HKA * bl.HK1_D1)

    bl.CFM_U2 = 0.5 * (bl.CFM_HKA * bl.HK2_U2 + bl.CFM_MA * bl.M2_U2 + bl.CFM_RTA * bl.RT2_U2)
    bl.CFM_T2 = 0.5 * (bl.CFM_HKA * bl.HK2_T2 + bl.CFM_RTA * bl.RT2_T2)
    bl.CFM_D2 = 0.5 * (bl.CFM_HKA * bl.HK2_D2)

    bl.CFM_MS = 0.5 * (
        bl.CFM_HKA * bl.HK1_MS
        + bl.CFM_MA * bl.M1_MS
        + bl.CFM_RTA * bl.RT1_MS
        + bl.CFM_HKA * bl.HK2_MS
        + bl.CFM_MA * bl.M2_MS
        + bl.CFM_RTA * bl.RT2_MS
    )
    bl.CFM_RE = 0.5 * (bl.CFM_RTA * bl.RT1_RE + bl.CFM_RTA * bl.RT2_RE)


def trdif(bl):
    bl1 = create2d(4, 5)
    bl2 = create2d(4, 5)
    blrez = [0.0] * 5
    blm = [0.0] * 5
    blr = [0.0] * 5
    blx = [0.0] * 5
    bt1 = create2d(4, 5)
    bt2 = create2d(4, 5)
    btrez = [0.0] * 5
    btm = [0.0] * 5
    btr = [0.0] * 5
    btx = [0.0] * 5

    for icom in range(1, NCOM + 1):
        bl.C1SAV[icom] = bl.COM1[icom]
        bl.C2SAV[icom] = bl.COM2[icom]

    wf2 = (bl.XT - bl.X1) / (bl.X2 - bl.X1)
    wf2_xt = 1.0 / (bl.X2 - bl.X1)

    wf2_a1 = wf2_xt * bl.XT_A1
    wf2_x1 = wf2_xt * bl.XT_X1 + (wf2 - 1.0) / (bl.X2 - bl.X1)
    wf2_x2 = wf2_xt * bl.XT_X2 - wf2 / (bl.X2 - bl.X1)
    wf2_t1 = wf2_xt * bl.XT_T1
    wf2_t2 = wf2_xt * bl.XT_T2
    wf2_d1 = wf2_xt * bl.XT_D1
    wf2_d2 = wf2_xt * bl.XT_D2
    wf2_u1 = wf2_xt * bl.XT_U1
    wf2_u2 = wf2_xt * bl.XT_U2
    wf2_ms = wf2_xt * bl.XT_MS
    wf2_re = wf2_xt * bl.XT_RE
    wf2_xf = wf2_xt * bl.XT_XF

    wf1 = 1.0 - wf2
    wf1_a1 = -wf2_a1
    wf1_x1 = -wf2_x1
    wf1_x2 = -wf2_x2
    wf1_t1 = -wf2_t1
    wf1_t2 = -wf2_t2
    wf1_d1 = -wf2_d1
    wf1_d2 = -wf2_d2
    wf1_u1 = -wf2_u1
    wf1_u2 = -wf2_u2
    wf1_ms = -wf2_ms
    wf1_re = -wf2_re
    wf1_xf = -wf2_xf

    tt = bl.T1 * wf1 + bl.T2 * wf2
    tt_a1 = bl.T1 * wf1_a1 + bl.T2 * wf2_a1
    tt_x1 = bl.T1 * wf1_x1 + bl.T2 * wf2_x1
    tt_x2 = bl.T1 * wf1_x2 + bl.T2 * wf2_x2
    tt_t1 = bl.T1 * wf1_t1 + bl.T2 * wf2_t1 + wf1
    tt_t2 = bl.T1 * wf1_t2 + bl.T2 * wf2_t2 + wf2
    tt_d1 = bl.T1 * wf1_d1 + bl.T2 * wf2_d1
    tt_d2 = bl.T1 * wf1_d2 + bl.T2 * wf2_d2
    tt_u1 = bl.T1 * wf1_u1 + bl.T2 * wf2_u1
    tt_u2 = bl.T1 * wf1_u2 + bl.T2 * wf2_u2
    tt_ms = bl.T1 * wf1_ms + bl.T2 * wf2_ms
    tt_re = bl.T1 * wf1_re + bl.T2 * wf2_re
    tt_xf = bl.T1 * wf1_xf + bl.T2 * wf2_xf

    dt = bl.D1 * wf1 + bl.D2 * wf2
    dt_a1 = bl.D1 * wf1_a1 + bl.D2 * wf2_a1
    dt_x1 = bl.D1 * wf1_x1 + bl.D2 * wf2_x1
    dt_x2 = bl.D1 * wf1_x2 + bl.D2 * wf2_x2
    dt_t1 = bl.D1 * wf1_t1 + bl.D2 * wf2_t1
    dt_t2 = bl.D1 * wf1_t2 + bl.D2 * wf2_t2
    dt_d1 = bl.D1 * wf1_d1 + bl.D2 * wf2_d1 + wf1
    dt_d2 = bl.D1 * wf1_d2 + bl.D2 * wf2_d2 + wf2
    dt_u1 = bl.D1 * wf1_u1 + bl.D2 * wf2_u1
    dt_u2 = bl.D1 * wf1_u2 + bl.D2 * wf2_u2
    dt_ms = bl.D1 * wf1_ms + bl.D2 * wf2_ms
    dt_re = bl.D1 * wf1_re + bl.D2 * wf2_re
    dt_xf = bl.D1 * wf1_xf + bl.D2 * wf2_xf

    ut = bl.U1 * wf1 + bl.U2 * wf2
    ut_a1 = bl.U1 * wf1_a1 + bl.U2 * wf2_a1
    ut_x1 = bl.U1 * wf1_x1 + bl.U2 * wf2_x1
    ut_x2 = bl.U1 * wf1_x2 + bl.U2 * wf2_x2
    ut_t1 = bl.U1 * wf1_t1 + bl.U2 * wf2_t1
    ut_t2 = bl.U1 * wf1_t2 + bl.U2 * wf2_t2
    ut_d1 = bl.U1 * wf1_d1 + bl.U2 * wf2_d1
    ut_d2 = bl.U1 * wf1_d2 + bl.U2 * wf2_d2
    ut_u1 = bl.U1 * wf1_u1 + bl.U2 * wf2_u1 + wf1
    ut_u2 = bl.U1 * wf1_u2 + bl.U2 * wf2_u2 + wf2
    ut_ms = bl.U1 * wf1_ms + bl.U2 * wf2_ms
    ut_re = bl.U1 * wf1_re + bl.U2 * wf2_re
    ut_xf = bl.U1 * wf1_xf + bl.U2 * wf2_xf

    bl.X2 = bl.XT
    bl.T2 = tt
    bl.D2 = dt
    bl.U2 = ut

    bl.AMPL2 = bl.AMCRIT
    bl.S2 = 0.0
    sync_vars_to_com(bl, 2)

    blkin(bl)
    blvar(bl, 1)
    blmid(bl, 1)
    bldif(bl, 1)

    for k in range(2, 4):
        blrez[k] = bl.VSREZ[k]
        blm[k] = bl.VSM[k] + bl.VS2[k][2] * tt_ms + bl.VS2[k][3] * dt_ms + bl.VS2[k][4] * ut_ms + bl.VS2[k][
            5
        ] * bl.XT_MS
        blr[k] = bl.VSR[k] + bl.VS2[k][2] * tt_re + bl.VS2[k][3] * dt_re + bl.VS2[k][4] * ut_re + bl.VS2[k][
            5
        ] * bl.XT_RE
        blx[k] = bl.VSX[k] + bl.VS2[k][2] * tt_xf + bl.VS2[k][3] * dt_xf + bl.VS2[k][4] * ut_xf + bl.VS2[k][
            5
        ] * bl.XT_XF

        bl1[k][1] = bl.VS1[k][1] + bl.VS2[k][2] * tt_a1 + bl.VS2[k][3] * dt_a1 + bl.VS2[k][4] * ut_a1 + bl.VS2[
            k
        ][5] * bl.XT_A1
        bl1[k][2] = bl.VS1[k][2] + bl.VS2[k][2] * tt_t1 + bl.VS2[k][3] * dt_t1 + bl.VS2[k][4] * ut_t1 + bl.VS2[
            k
        ][5] * bl.XT_T1
        bl1[k][3] = bl.VS1[k][3] + bl.VS2[k][2] * tt_d1 + bl.VS2[k][3] * dt_d1 + bl.VS2[k][4] * ut_d1 + bl.VS2[
            k
        ][5] * bl.XT_D1
        bl1[k][4] = bl.VS1[k][4] + bl.VS2[k][2] * tt_u1 + bl.VS2[k][3] * dt_u1 + bl.VS2[k][4] * ut_u1 + bl.VS2[
            k
        ][5] * bl.XT_U1
        bl1[k][5] = bl.VS1[k][5] + bl.VS2[k][2] * tt_x1 + bl.VS2[k][3] * dt_x1 + bl.VS2[k][4] * ut_x1 + bl.VS2[
            k
        ][5] * bl.XT_X1

        bl2[k][1] = 0.0
        bl2[k][2] = bl.VS2[k][2] * tt_t2 + bl.VS2[k][3] * dt_t2 + bl.VS2[k][4] * ut_t2 + bl.VS2[k][5] * bl.XT_T2
        bl2[k][3] = bl.VS2[k][2] * tt_d2 + bl.VS2[k][3] * dt_d2 + bl.VS2[k][4] * ut_d2 + bl.VS2[k][5] * bl.XT_D2
        bl2[k][4] = bl.VS2[k][2] * tt_u2 + bl.VS2[k][3] * dt_u2 + bl.VS2[k][4] * ut_u2 + bl.VS2[k][5] * bl.XT_U2
        bl2[k][5] = bl.VS2[k][2] * tt_x2 + bl.VS2[k][3] * dt_x2 + bl.VS2[k][4] * ut_x2 + bl.VS2[k][5] * bl.XT_X2

    blvar(bl, 2)
    ctr = bl.CTRCON * math.exp(-bl.CTRCEX / (bl.HK2 - 1.0))
    ctr_hk2 = ctr * bl.CTRCEX / (bl.HK2 - 1.0) ** 2

    st = ctr * bl.CQ2
    st_tt = ctr * bl.CQ2_T2 + bl.CQ2 * ctr_hk2 * bl.HK2_T2
    st_dt = ctr * bl.CQ2_D2 + bl.CQ2 * ctr_hk2 * bl.HK2_D2
    st_ut = ctr * bl.CQ2_U2 + bl.CQ2 * ctr_hk2 * bl.HK2_U2
    st_ms = ctr * bl.CQ2_MS + bl.CQ2 * ctr_hk2 * bl.HK2_MS
    st_re = ctr * bl.CQ2_RE

    st_a1 = st_tt * tt_a1 + st_dt * dt_a1 + st_ut * ut_a1
    st_x1 = st_tt * tt_x1 + st_dt * dt_x1 + st_ut * ut_x1
    st_x2 = st_tt * tt_x2 + st_dt * dt_x2 + st_ut * ut_x2
    st_t1 = st_tt * tt_t1 + st_dt * dt_t1 + st_ut * ut_t1
    st_t2 = st_tt * tt_t2 + st_dt * dt_t2 + st_ut * ut_t2
    st_d1 = st_tt * tt_d1 + st_dt * dt_d1 + st_ut * ut_d1
    st_d2 = st_tt * tt_d2 + st_dt * dt_d2 + st_ut * ut_d2
    st_u1 = st_tt * tt_u1 + st_dt * dt_u1 + st_ut * ut_u1
    st_u2 = st_tt * tt_u2 + st_dt * dt_u2 + st_ut * ut_u2
    st_ms2 = st_tt * tt_ms + st_dt * dt_ms + st_ut * ut_ms + st_ms
    st_re2 = st_tt * tt_re + st_dt * dt_re + st_ut * ut_re + st_re
    st_xf = st_tt * tt_xf + st_dt * dt_xf + st_ut * ut_xf

    bl.AMPL2 = 0.0
    bl.S2 = st

    blvar(bl, 2)

    for icom in range(1, NCOM + 1):
        bl.COM1[icom] = bl.COM2[icom]
        bl.COM2[icom] = bl.C2SAV[icom]
    sync_com_to_vars(bl, 1)
    sync_com_to_vars(bl, 2)

    blmid(bl, 2)
    bldif(bl, 2)

    for k in range(1, 4):
        btrez[k] = bl.VSREZ[k]
        btm[k] = (
            bl.VSM[k]
            + bl.VS1[k][1] * st_ms2
            + bl.VS1[k][2] * tt_ms
            + bl.VS1[k][3] * dt_ms
            + bl.VS1[k][4] * ut_ms
            + bl.VS1[k][5] * bl.XT_MS
        )
        btr[k] = (
            bl.VSR[k]
            + bl.VS1[k][1] * st_re2
            + bl.VS1[k][2] * tt_re
            + bl.VS1[k][3] * dt_re
            + bl.VS1[k][4] * ut_re
            + bl.VS1[k][5] * bl.XT_RE
        )
        btx[k] = (
            bl.VSX[k]
            + bl.VS1[k][1] * st_xf
            + bl.VS1[k][2] * tt_xf
            + bl.VS1[k][3] * dt_xf
            + bl.VS1[k][4] * ut_xf
            + bl.VS1[k][5] * bl.XT_XF
        )

        bt1[k][1] = (
            bl.VS1[k][1] * st_a1
            + bl.VS1[k][2] * tt_a1
            + bl.VS1[k][3] * dt_a1
            + bl.VS1[k][4] * ut_a1
            + bl.VS1[k][5] * bl.XT_A1
        )
        bt1[k][2] = (
            bl.VS1[k][1] * st_t1
            + bl.VS1[k][2] * tt_t1
            + bl.VS1[k][3] * dt_t1
            + bl.VS1[k][4] * ut_t1
            + bl.VS1[k][5] * bl.XT_T1
        )
        bt1[k][3] = (
            bl.VS1[k][1] * st_d1
            + bl.VS1[k][2] * tt_d1
            + bl.VS1[k][3] * dt_d1
            + bl.VS1[k][4] * ut_d1
            + bl.VS1[k][5] * bl.XT_D1
        )
        bt1[k][4] = (
            bl.VS1[k][1] * st_u1
            + bl.VS1[k][2] * tt_u1
            + bl.VS1[k][3] * dt_u1
            + bl.VS1[k][4] * ut_u1
            + bl.VS1[k][5] * bl.XT_U1
        )
        bt1[k][5] = (
            bl.VS1[k][1] * st_x1
            + bl.VS1[k][2] * tt_x1
            + bl.VS1[k][3] * dt_x1
            + bl.VS1[k][4] * ut_x1
            + bl.VS1[k][5] * bl.XT_X1
        )

        bt2[k][1] = bl.VS2[k][1]
        bt2[k][2] = (
            bl.VS2[k][2]
            + bl.VS1[k][1] * st_t2
            + bl.VS1[k][2] * tt_t2
            + bl.VS1[k][3] * dt_t2
            + bl.VS1[k][4] * ut_t2
            + bl.VS1[k][5] * bl.XT_T2
        )
        bt2[k][3] = (
            bl.VS2[k][3]
            + bl.VS1[k][1] * st_d2
            + bl.VS1[k][2] * tt_d2
            + bl.VS1[k][3] * dt_d2
            + bl.VS1[k][4] * ut_d2
            + bl.VS1[k][5] * bl.XT_D2
        )
        bt2[k][4] = (
            bl.VS2[k][4]
            + bl.VS1[k][1] * st_u2
            + bl.VS1[k][2] * tt_u2
            + bl.VS1[k][3] * dt_u2
            + bl.VS1[k][4] * ut_u2
            + bl.VS1[k][5] * bl.XT_U2
        )
        bt2[k][5] = (
            bl.VS2[k][5]
            + bl.VS1[k][1] * st_x2
            + bl.VS1[k][2] * tt_x2
            + bl.VS1[k][3] * dt_x2
            + bl.VS1[k][4] * ut_x2
            + bl.VS1[k][5] * bl.XT_X2
        )

    bl.VSREZ[1] = btrez[1]
    bl.VSREZ[2] = blrez[2] + btrez[2]
    bl.VSREZ[3] = blrez[3] + btrez[3]
    bl.VSM[1] = btm[1]
    bl.VSM[2] = blm[2] + btm[2]
    bl.VSM[3] = blm[3] + btm[3]
    bl.VSR[1] = btr[1]
    bl.VSR[2] = blr[2] + btr[2]
    bl.VSR[3] = blr[3] + btr[3]
    bl.VSX[1] = btx[1]
    bl.VSX[2] = blx[2] + btx[2]
    bl.VSX[3] = blx[3] + btx[3]
    for l in range(1, 6):
        bl.VS1[1][l] = bt1[1][l]
        bl.VS2[1][l] = bt2[1][l]
        bl.VS1[2][l] = bl1[2][l] + bt1[2][l]
        bl.VS2[2][l] = bl2[2][l] + bt2[2][l]
        bl.VS1[3][l] = bl1[3][l] + bt1[3][l]
        bl.VS2[3][l] = bl2[3][l] + bt2[3][l]

    for icom in range(1, NCOM + 1):
        bl.COM1[icom] = bl.C1SAV[icom]
    sync_com_to_vars(bl, 1)


def bldif(bl, ityp):
    if ityp == 0:
        xlog = 1.0
        ulog = bl.BULE
        tlog = 0.5 * (1.0 - bl.BULE)
        hlog = 0.0
        ddlog = 0.0
    else:
        xlog = math.log(bl.X2 / bl.X1)
        ulog = math.log(bl.U2 / bl.U1)
        tlog = math.log(bl.T2 / bl.T1)
        hlog = math.log(bl.HS2 / bl.HS1)
        ddlog = 1.0

    for k in range(1, 5):
        bl.VSREZ[k] = 0.0
        bl.VSM[k] = 0.0
        bl.VSR[k] = 0.0
        bl.VSX[k] = 0.0
        for l in range(1, 6):
            bl.VS1[k][l] = 0.0
            bl.VS2[k][l] = 0.0

    hupwt = 1.0

    hdcon = 5.0 * hupwt / bl.HK2**2
    hd_hk1 = 0.0
    hd_hk2 = -hdcon * 2.0 / bl.HK2

    if ityp == 3:
        hdcon = hupwt / bl.HK2**2
        hd_hk1 = 0.0
        hd_hk2 = -hdcon * 2.0 / bl.HK2

    arg = abs((bl.HK2 - 1.0) / (bl.HK1 - 1.0))
    hl = math.log(arg)
    hl_hk1 = -1.0 / (bl.HK1 - 1.0)
    hl_hk2 = 1.0 / (bl.HK2 - 1.0)

    hlsq = min(hl**2, 15.0)
    ehh = math.exp(-hlsq * hdcon)
    upw = 1.0 - 0.5 * ehh
    upw_hl = ehh * hl * hdcon
    upw_hd = 0.5 * ehh * hlsq

    upw_hk1 = upw_hl * hl_hk1 + upw_hd * hd_hk1
    upw_hk2 = upw_hl * hl_hk2 + upw_hd * hd_hk2

    upw_u1 = upw_hk1 * bl.HK1_U1
    upw_t1 = upw_hk1 * bl.HK1_T1
    upw_d1 = upw_hk1 * bl.HK1_D1
    upw_u2 = upw_hk2 * bl.HK2_U2
    upw_t2 = upw_hk2 * bl.HK2_T2
    upw_d2 = upw_hk2 * bl.HK2_D2
    upw_ms = upw_hk1 * bl.HK1_MS + upw_hk2 * bl.HK2_MS

    if ityp == 0:
        bl.VS2[1][1] = 1.0
        bl.VSR[1] = 0.0
        bl.VSREZ[1] = -bl.AMPL2
    elif ityp == 1:
        ax, ax_hk1, ax_t1, ax_rt1, ax_a1, ax_hk2, ax_t2, ax_rt2, ax_a2 = axset(
            bl.HK1, bl.T1, bl.RT1, bl.AMPL1, bl.HK2, bl.T2, bl.RT2, bl.AMPL2, bl.AMCRIT, bl.IDAMPV
        )
        rezc = bl.AMPL2 - bl.AMPL1 - ax * (bl.X2 - bl.X1)
        z_ax = -(bl.X2 - bl.X1)

        bl.VS1[1][1] = z_ax * ax_a1 - 1.0
        bl.VS1[1][2] = z_ax * (ax_hk1 * bl.HK1_T1 + ax_t1 + ax_rt1 * bl.RT1_T1)
        bl.VS1[1][3] = z_ax * (ax_hk1 * bl.HK1_D1)
        bl.VS1[1][4] = z_ax * (ax_hk1 * bl.HK1_U1 + ax_rt1 * bl.RT1_U1)
        bl.VS1[1][5] = ax
        bl.VS2[1][1] = z_ax * ax_a2 + 1.0
        bl.VS2[1][2] = z_ax * (ax_hk2 * bl.HK2_T2 + ax_t2 + ax_rt2 * bl.RT2_T2)
        bl.VS2[1][3] = z_ax * (ax_hk2 * bl.HK2_D2)
        bl.VS2[1][4] = z_ax * (ax_hk2 * bl.HK2_U2 + ax_rt2 * bl.RT2_U2)
        bl.VS2[1][5] = -ax
        bl.VSM[1] = z_ax * (ax_hk1 * bl.HK1_MS + ax_rt1 * bl.RT1_MS + ax_hk2 * bl.HK2_MS + ax_rt2 * bl.RT2_MS)
        bl.VSR[1] = z_ax * (ax_rt1 * bl.RT1_RE + ax_rt2 * bl.RT2_RE)
        bl.VSX[1] = 0.0
        bl.VSREZ[1] = -rezc
    else:
        sa = (1.0 - upw) * bl.S1 + upw * bl.S2
        cqa = (1.0 - upw) * bl.CQ1 + upw * bl.CQ2
        cfa = (1.0 - upw) * bl.CF1 + upw * bl.CF2
        hka = (1.0 - upw) * bl.HK1 + upw * bl.HK2

        usa = 0.5 * (bl.US1 + bl.US2)
        rta = 0.5 * (bl.RT1 + bl.RT2)
        dea = 0.5 * (bl.DE1 + bl.DE2)
        da = 0.5 * (bl.D1 + bl.D2)

        if ityp == 3:
            ald = bl.DLCON
        else:
            ald = 1.0

        if ityp == 2:
            gcc = bl.GCCON
            hkc = hka - 1.0 - gcc / rta
            hkc_hka = 1.0
            hkc_rta = gcc / rta**2
            if hkc < 0.01:
                hkc = 0.01
                hkc_hka = 0.0
                hkc_rta = 0.0
        else:
            hkc = hka - 1.0
            hkc_hka = 1.0
            hkc_rta = 0.0

        hr = hkc / (bl.GACON * ald * hka)
        hr_hka = hkc_hka / (bl.GACON * ald * hka) - hr / hka
        hr_rta = hkc_rta / (bl.GACON * ald * hka)

        uq = (0.5 * cfa - hr**2) / (bl.GBCON * da)
        uq_hka = -2.0 * hr * hr_hka / (bl.GBCON * da)
        uq_rta = -2.0 * hr * hr_rta / (bl.GBCON * da)
        uq_cfa = 0.5 / (bl.GBCON * da)
        uq_da = -uq / da
        uq_upw = uq_cfa * (bl.CF2 - bl.CF1) + uq_hka * (bl.HK2 - bl.HK1)

        uq_t1 = (1.0 - upw) * (uq_cfa * bl.CF1_T1 + uq_hka * bl.HK1_T1) + uq_upw * upw_t1
        uq_d1 = (1.0 - upw) * (uq_cfa * bl.CF1_D1 + uq_hka * bl.HK1_D1) + uq_upw * upw_d1
        uq_u1 = (1.0 - upw) * (uq_cfa * bl.CF1_U1 + uq_hka * bl.HK1_U1) + uq_upw * upw_u1
        uq_t2 = upw * (uq_cfa * bl.CF2_T2 + uq_hka * bl.HK2_T2) + uq_upw * upw_t2
        uq_d2 = upw * (uq_cfa * bl.CF2_D2 + uq_hka * bl.HK2_D2) + uq_upw * upw_d2
        uq_u2 = upw * (uq_cfa * bl.CF2_U2 + uq_hka * bl.HK2_U2) + uq_upw * upw_u2
        uq_ms = (1.0 - upw) * (uq_cfa * bl.CF1_MS + uq_hka * bl.HK1_MS) + uq_upw * upw_ms + upw * (
            uq_cfa * bl.CF2_MS + uq_hka * bl.HK2_MS
        )
        uq_re = (1.0 - upw) * uq_cfa * bl.CF1_RE + upw * uq_cfa * bl.CF2_RE

        uq_t1 = uq_t1 + 0.5 * uq_rta * bl.RT1_T1
        uq_d1 = uq_d1 + 0.5 * uq_da
        uq_u1 = uq_u1 + 0.5 * uq_rta * bl.RT1_U1
        uq_t2 = uq_t2 + 0.5 * uq_rta * bl.RT2_T2
        uq_d2 = uq_d2 + 0.5 * uq_da
        uq_u2 = uq_u2 + 0.5 * uq_rta * bl.RT2_U2
        uq_ms = uq_ms + 0.5 * uq_rta * bl.RT1_MS + 0.5 * uq_rta * bl.RT2_MS
        uq_re = uq_re + 0.5 * uq_rta * bl.RT1_RE + 0.5 * uq_rta * bl.RT2_RE

        scc = bl.SCCON * 1.333 / (1.0 + usa)
        scc_usa = -scc / (1.0 + usa)

        scc_us1 = scc_usa * 0.5
        scc_us2 = scc_usa * 0.5

        slog = math.log(bl.S2 / bl.S1)
        dxi = bl.X2 - bl.X1

        rezc = scc * (cqa - sa * ald) * dxi - dea * 2.0 * slog + dea * 2.0 * (uq * dxi - ulog) * bl.DUXCON

        z_cfa = dea * 2.0 * uq_cfa * dxi * bl.DUXCON
        z_hka = dea * 2.0 * uq_hka * dxi * bl.DUXCON
        z_da = dea * 2.0 * uq_da * dxi * bl.DUXCON
        z_sl = -dea * 2.0
        z_ul = -dea * 2.0 * bl.DUXCON
        z_dxi = scc * (cqa - sa * ald) + dea * 2.0 * uq * bl.DUXCON
        z_usa = scc_usa * (cqa - sa * ald) * dxi
        z_cqa = scc * dxi
        z_sa = -scc * dxi * ald
        z_dea = 2.0 * ((uq * dxi - ulog) * bl.DUXCON - slog)

        z_upw = z_cqa * (bl.CQ2 - bl.CQ1) + z_sa * (bl.S2 - bl.S1) + z_cfa * (bl.CF2 - bl.CF1) + z_hka * (
            bl.HK2 - bl.HK1
        )
        z_de1 = 0.5 * z_dea
        z_de2 = 0.5 * z_dea
        z_us1 = 0.5 * z_usa
        z_us2 = 0.5 * z_usa
        z_d1 = 0.5 * z_da
        z_d2 = 0.5 * z_da
        z_u1 = -z_ul / bl.U1
        z_u2 = z_ul / bl.U2
        z_x1 = -z_dxi
        z_x2 = z_dxi
        z_s1 = (1.0 - upw) * z_sa - z_sl / bl.S1
        z_s2 = upw * z_sa + z_sl / bl.S2
        z_cq1 = (1.0 - upw) * z_cqa
        z_cq2 = upw * z_cqa
        z_cf1 = (1.0 - upw) * z_cfa
        z_cf2 = upw * z_cfa
        z_hk1 = (1.0 - upw) * z_hka
        z_hk2 = upw * z_hka

        bl.VS1[1][1] = z_s1
        bl.VS1[1][2] = z_upw * upw_t1 + z_de1 * bl.DE1_T1 + z_us1 * bl.US1_T1
        bl.VS1[1][3] = z_d1 + z_upw * upw_d1 + z_de1 * bl.DE1_D1 + z_us1 * bl.US1_D1
        bl.VS1[1][4] = z_u1 + z_upw * upw_u1 + z_de1 * bl.DE1_U1 + z_us1 * bl.US1_U1
        bl.VS1[1][5] = z_x1
        bl.VS2[1][1] = z_s2
        bl.VS2[1][2] = z_upw * upw_t2 + z_de2 * bl.DE2_T2 + z_us2 * bl.US2_T2
        bl.VS2[1][3] = z_d2 + z_upw * upw_d2 + z_de2 * bl.DE2_D2 + z_us2 * bl.US2_D2
        bl.VS2[1][4] = z_u2 + z_upw * upw_u2 + z_de2 * bl.DE2_U2 + z_us2 * bl.US2_U2
        bl.VS2[1][5] = z_x2
        bl.VSM[1] = z_upw * upw_ms + z_de1 * bl.DE1_MS + z_us1 * bl.US1_MS + z_de2 * bl.DE2_MS + z_us2 * bl.US2_MS

        bl.VS1[1][2] = bl.VS1[1][2] + z_cq1 * bl.CQ1_T1 + z_cf1 * bl.CF1_T1 + z_hk1 * bl.HK1_T1
        bl.VS1[1][3] = bl.VS1[1][3] + z_cq1 * bl.CQ1_D1 + z_cf1 * bl.CF1_D1 + z_hk1 * bl.HK1_D1
        bl.VS1[1][4] = bl.VS1[1][4] + z_cq1 * bl.CQ1_U1 + z_cf1 * bl.CF1_U1 + z_hk1 * bl.HK1_U1

        bl.VS2[1][2] = bl.VS2[1][2] + z_cq2 * bl.CQ2_T2 + z_cf2 * bl.CF2_T2 + z_hk2 * bl.HK2_T2
        bl.VS2[1][3] = bl.VS2[1][3] + z_cq2 * bl.CQ2_D2 + z_cf2 * bl.CF2_D2 + z_hk2 * bl.HK2_D2
        bl.VS2[1][4] = bl.VS2[1][4] + z_cq2 * bl.CQ2_U2 + z_cf2 * bl.CF2_U2 + z_hk2 * bl.HK2_U2

        bl.VSM[1] = bl.VSM[1] + z_cq1 * bl.CQ1_MS + z_cf1 * bl.CF1_MS + z_hk1 * bl.HK1_MS + z_cq2 * bl.CQ2_MS + z_cf2 * bl.CF2_MS + z_hk2 * bl.HK2_MS
        bl.VSR[1] = z_cq1 * bl.CQ1_RE + z_cf1 * bl.CF1_RE + z_cq2 * bl.CQ2_RE + z_cf2 * bl.CF2_RE
        bl.VSX[1] = 0.0
        bl.VSREZ[1] = -rezc

    ha = 0.5 * (bl.H1 + bl.H2)
    ma = 0.5 * (bl.M1 + bl.M2)
    xa = 0.5 * (bl.X1 + bl.X2)
    ta = 0.5 * (bl.T1 + bl.T2)
    hwa = 0.5 * (bl.DW1 / bl.T1 + bl.DW2 / bl.T2)

    cfx = 0.50 * bl.CFM * xa / ta + 0.25 * (bl.CF1 * bl.X1 / bl.T1 + bl.CF2 * bl.X2 / bl.T2)
    cfx_xa = 0.50 * bl.CFM / ta
    cfx_ta = -0.50 * bl.CFM * xa / ta**2

    cfx_x1 = 0.25 * bl.CF1 / bl.T1 + cfx_xa * 0.5
    cfx_x2 = 0.25 * bl.CF2 / bl.T2 + cfx_xa * 0.5
    cfx_t1 = -0.25 * bl.CF1 * bl.X1 / bl.T1**2 + cfx_ta * 0.5
    cfx_t2 = -0.25 * bl.CF2 * bl.X2 / bl.T2**2 + cfx_ta * 0.5
    cfx_cf1 = 0.25 * bl.X1 / bl.T1
    cfx_cf2 = 0.25 * bl.X2 / bl.T2
    cfx_cfm = 0.50 * xa / ta

    btmp = ha + 2.0 - ma + hwa

    rezt = tlog + btmp * ulog - xlog * 0.5 * cfx
    z_cfx = -xlog * 0.5
    z_ha = ulog
    z_hwa = ulog
    z_ma = -ulog
    z_xl = -ddlog * 0.5 * cfx
    z_ul = ddlog * btmp
    z_tl = ddlog

    z_cfm = z_cfx * cfx_cfm
    z_cf1 = z_cfx * cfx_cf1
    z_cf2 = z_cfx * cfx_cf2

    z_t1 = -z_tl / bl.T1 + z_cfx * cfx_t1 + z_hwa * 0.5 * (-bl.DW1 / bl.T1**2)
    z_t2 = z_tl / bl.T2 + z_cfx * cfx_t2 + z_hwa * 0.5 * (-bl.DW2 / bl.T2**2)
    z_x1 = -z_xl / bl.X1 + z_cfx * cfx_x1
    z_x2 = z_xl / bl.X2 + z_cfx * cfx_x2
    z_u1 = -z_ul / bl.U1
    z_u2 = z_ul / bl.U2

    bl.VS1[2][2] = 0.5 * z_ha * bl.H1_T1 + z_cfm * bl.CFM_T1 + z_cf1 * bl.CF1_T1 + z_t1
    bl.VS1[2][3] = 0.5 * z_ha * bl.H1_D1 + z_cfm * bl.CFM_D1 + z_cf1 * bl.CF1_D1
    bl.VS1[2][4] = 0.5 * z_ma * bl.M1_U1 + z_cfm * bl.CFM_U1 + z_cf1 * bl.CF1_U1 + z_u1
    bl.VS1[2][5] = z_x1
    bl.VS2[2][2] = 0.5 * z_ha * bl.H2_T2 + z_cfm * bl.CFM_T2 + z_cf2 * bl.CF2_T2 + z_t2
    bl.VS2[2][3] = 0.5 * z_ha * bl.H2_D2 + z_cfm * bl.CFM_D2 + z_cf2 * bl.CF2_D2
    bl.VS2[2][4] = 0.5 * z_ma * bl.M2_U2 + z_cfm * bl.CFM_U2 + z_cf2 * bl.CF2_U2 + z_u2
    bl.VS2[2][5] = z_x2

    bl.VSM[2] = 0.5 * z_ma * bl.M1_MS + z_cfm * bl.CFM_MS + z_cf1 * bl.CF1_MS + 0.5 * z_ma * bl.M2_MS + z_cf2 * bl.CF2_MS
    bl.VSR[2] = z_cfm * bl.CFM_RE + z_cf1 * bl.CF1_RE + z_cf2 * bl.CF2_RE
    bl.VSX[2] = 0.0
    bl.VSREZ[2] = -rezt

    xot1 = bl.X1 / bl.T1
    xot2 = bl.X2 / bl.T2

    ha = 0.5 * (bl.H1 + bl.H2)
    hsa = 0.5 * (bl.HS1 + bl.HS2)
    hca = 0.5 * (bl.HC1 + bl.HC2)
    hwa = 0.5 * (bl.DW1 / bl.T1 + bl.DW2 / bl.T2)

    dix = (1.0 - upw) * bl.DI1 * xot1 + upw * bl.DI2 * xot2
    cfx = (1.0 - upw) * bl.CF1 * xot1 + upw * bl.CF2 * xot2
    dix_upw = bl.DI2 * xot2 - bl.DI1 * xot1
    cfx_upw = bl.CF2 * xot2 - bl.CF1 * xot1

    btmp = 2.0 * hca / hsa + 1.0 - ha - hwa

    rezh = hlog + btmp * ulog + xlog * (0.5 * cfx - dix)
    z_cfx = xlog * 0.5
    z_dix = -xlog
    z_hca = 2.0 * ulog / hsa
    z_ha = -ulog
    z_hwa = -ulog
    z_xl = ddlog * (0.5 * cfx - dix)
    z_ul = ddlog * btmp
    z_hl = ddlog

    z_upw = z_cfx * cfx_upw + z_dix * dix_upw

    z_hs1 = -hca * ulog / hsa**2 - z_hl / bl.HS1
    z_hs2 = -hca * ulog / hsa**2 + z_hl / bl.HS2

    z_cf1 = (1.0 - upw) * z_cfx * xot1
    z_cf2 = upw * z_cfx * xot2
    z_di1 = (1.0 - upw) * z_dix * xot1
    z_di2 = upw * z_dix * xot2

    z_t1 = (1.0 - upw) * (z_cfx * bl.CF1 + z_dix * bl.DI1) * (-xot1 / bl.T1)
    z_t2 = upw * (z_cfx * bl.CF2 + z_dix * bl.DI2) * (-xot2 / bl.T2)
    z_x1 = (1.0 - upw) * (z_cfx * bl.CF1 + z_dix * bl.DI1) / bl.T1 - z_xl / bl.X1
    z_x2 = upw * (z_cfx * bl.CF2 + z_dix * bl.DI2) / bl.T2 + z_xl / bl.X2
    z_u1 = -z_ul / bl.U1
    z_u2 = z_ul / bl.U2

    z_t1 = z_t1 + z_hwa * 0.5 * (-bl.DW1 / bl.T1**2)
    z_t2 = z_t2 + z_hwa * 0.5 * (-bl.DW2 / bl.T2**2)

    bl.VS1[3][1] = z_di1 * bl.DI1_S1
    bl.VS1[3][2] = z_hs1 * bl.HS1_T1 + z_cf1 * bl.CF1_T1 + z_di1 * bl.DI1_T1 + z_t1
    bl.VS1[3][3] = z_hs1 * bl.HS1_D1 + z_cf1 * bl.CF1_D1 + z_di1 * bl.DI1_D1
    bl.VS1[3][4] = z_hs1 * bl.HS1_U1 + z_cf1 * bl.CF1_U1 + z_di1 * bl.DI1_U1 + z_u1
    bl.VS1[3][5] = z_x1
    bl.VS2[3][1] = z_di2 * bl.DI2_S2
    bl.VS2[3][2] = z_hs2 * bl.HS2_T2 + z_cf2 * bl.CF2_T2 + z_di2 * bl.DI2_T2 + z_t2
    bl.VS2[3][3] = z_hs2 * bl.HS2_D2 + z_cf2 * bl.CF2_D2 + z_di2 * bl.DI2_D2
    bl.VS2[3][4] = z_hs2 * bl.HS2_U2 + z_cf2 * bl.CF2_U2 + z_di2 * bl.DI2_U2 + z_u2
    bl.VS2[3][5] = z_x2
    bl.VSM[3] = z_hs1 * bl.HS1_MS + z_cf1 * bl.CF1_MS + z_di1 * bl.DI1_MS + z_hs2 * bl.HS2_MS + z_cf2 * bl.CF2_MS + z_di2 * bl.DI2_MS
    bl.VSR[3] = z_hs1 * bl.HS1_RE + z_cf1 * bl.CF1_RE + z_di1 * bl.DI1_RE + z_hs2 * bl.HS2_RE + z_cf2 * bl.CF2_RE + z_di2 * bl.DI2_RE

    bl.VS1[3][2] = bl.VS1[3][2] + 0.5 * (z_hca * bl.HC1_T1 + z_ha * bl.H1_T1) + z_upw * upw_t1
    bl.VS1[3][3] = bl.VS1[3][3] + 0.5 * (z_hca * bl.HC1_D1 + z_ha * bl.H1_D1) + z_upw * upw_d1
    bl.VS1[3][4] = bl.VS1[3][4] + 0.5 * (z_hca * bl.HC1_U1) + z_upw * upw_u1
    bl.VS2[3][2] = bl.VS2[3][2] + 0.5 * (z_hca * bl.HC2_T2 + z_ha * bl.H2_T2) + z_upw * upw_t2
    bl.VS2[3][3] = bl.VS2[3][3] + 0.5 * (z_hca * bl.HC2_D2 + z_ha * bl.H2_D2) + z_upw * upw_d2
    bl.VS2[3][4] = bl.VS2[3][4] + 0.5 * (z_hca * bl.HC2_U2) + z_upw * upw_u2

    bl.VSM[3] = bl.VSM[3] + 0.5 * (z_hca * bl.HC1_MS) + z_upw * upw_ms + 0.5 * (z_hca * bl.HC2_MS)

    bl.VSX[3] = 0.0
    bl.VSREZ[3] = -rezh


def dampl(hk, th, rt):
    dgr = 0.08
    hmi = 1.0 / (hk - 1.0)
    hmi_hk = -(hmi**2)

    aa = 2.492 * hmi**0.43
    aa_hk = (aa / hmi) * 0.43 * hmi_hk

    bb = math.tanh(14.0 * hmi - 9.24)
    bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk

    grcrit = aa + 0.7 * (bb + 1.0)
    grc_hk = aa_hk + 0.7 * bb_hk

    gr = math.log10(rt)
    gr_rt = 1.0 / (2.3025851 * rt)

    if gr < grcrit - dgr:
        return 0.0, 0.0, 0.0, 0.0

    rnorm = (gr - (grcrit - dgr)) / (2.0 * dgr)
    rn_hk = -grc_hk / (2.0 * dgr)
    rn_rt = gr_rt / (2.0 * dgr)

    if rnorm >= 1.0:
        rfac = 1.0
        rfac_hk = 0.0
        rfac_rt = 0.0
    else:
        rfac = 3.0 * rnorm**2 - 2.0 * rnorm**3
        rfac_rn = 6.0 * rnorm - 6.0 * rnorm**2
        rfac_hk = rfac_rn * rn_hk
        rfac_rt = rfac_rn * rn_rt

    arg = 3.87 * hmi - 2.52
    arg_hk = 3.87 * hmi_hk
    ex = math.exp(-(arg**2))
    ex_hk = ex * (-2.0 * arg * arg_hk)

    dadr = 0.028 * (hk - 1.0) - 0.0345 * ex
    dadr_hk = 0.028 - 0.0345 * ex_hk

    af = -0.05 + 2.7 * hmi - 5.5 * hmi**2 + 3.0 * hmi**3
    af_hmi = 2.7 - 11.0 * hmi + 9.0 * hmi**2
    af_hk = af_hmi * hmi_hk

    ax = (af * dadr / th) * rfac
    ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac + (af * dadr / th) * rfac_hk
    ax_th = -ax / th
    ax_rt = (af * dadr / th) * rfac_rt

    return ax, ax_hk, ax_th, ax_rt


def dampl2(hk, th, rt):
    dgr = 0.08
    hk1 = 3.5
    hk2 = 4.0
    hmi = 1.0 / (hk - 1.0)
    hmi_hk = -(hmi**2)

    aa = 2.492 * hmi**0.43
    aa_hk = (aa / hmi) * 0.43 * hmi_hk

    bb = math.tanh(14.0 * hmi - 9.24)
    bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk

    grc = aa + 0.7 * (bb + 1.0)
    grc_hk = aa_hk + 0.7 * bb_hk

    gr = math.log10(rt)
    gr_rt = 1.0 / (2.3025851 * rt)

    if gr < grc - dgr:
        return 0.0, 0.0, 0.0, 0.0

    rnorm = (gr - (grc - dgr)) / (2.0 * dgr)
    rn_hk = -grc_hk / (2.0 * dgr)
    rn_rt = gr_rt / (2.0 * dgr)

    if rnorm >= 1.0:
        rfac = 1.0
        rfac_hk = 0.0
        rfac_rt = 0.0
    else:
        rfac = 3.0 * rnorm**2 - 2.0 * rnorm**3
        rfac_rn = 6.0 * rnorm - 6.0 * rnorm**2
        rfac_hk = rfac_rn * rn_hk
        rfac_rt = rfac_rn * rn_rt

    arg = 3.87 * hmi - 2.52
    arg_hk = 3.87 * hmi_hk
    ex = math.exp(-(arg**2))
    ex_hk = ex * (-2.0 * arg * arg_hk)

    dadr = 0.028 * (hk - 1.0) - 0.0345 * ex
    dadr_hk = 0.028 - 0.0345 * ex_hk

    brg = -20.0 * hmi
    af = -0.05 + 2.7 * hmi - 5.5 * hmi**2 + 3.0 * hmi**3 + 0.1 * math.exp(brg)
    af_hmi = 2.7 - 11.0 * hmi + 9.0 * hmi**2 - 2.0 * math.exp(brg)
    af_hk = af_hmi * hmi_hk

    ax = (af * dadr / th) * rfac
    ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac + (af * dadr / th) * rfac_hk
    ax_th = -ax / th
    ax_rt = (af * dadr / th) * rfac_rt

    if hk < hk1:
        return ax, ax_hk, ax_th, ax_rt

    hnorm = (hk - hk1) / (hk2 - hk1)
    hn_hk = 1.0 / (hk2 - hk1)

    if hnorm >= 1.0:
        hfac = 1.0
        hf_hk = 0.0
    else:
        hfac = 3.0 * hnorm**2 - 2.0 * hnorm**3
        hf_hk = (6.0 * hnorm - 6.0 * hnorm**2) * hn_hk

    ax1 = ax
    ax1_hk = ax_hk
    ax1_th = ax_th
    ax1_rt = ax_rt

    gr0 = 0.30 + 0.35 * math.exp(-0.15 * (hk - 5.0))
    gr0_hk = -0.35 * math.exp(-0.15 * (hk - 5.0)) * 0.15

    tnr = math.tanh(1.2 * (gr - gr0))
    tnr_rt = (1.0 - tnr**2) * 1.2 * gr_rt
    tnr_hk = -(1.0 - tnr**2) * 1.2 * gr0_hk

    ax2 = (0.086 * tnr - 0.25 / (hk - 1.0) ** 1.5) / th
    ax2_hk = (0.086 * tnr_hk + 1.5 * 0.25 / (hk - 1.0) ** 2.5) / th
    ax2_rt = (0.086 * tnr_rt) / th
    ax2_th = -ax2 / th

    if ax2 < 0.0:
        ax2 = 0.0
        ax2_hk = 0.0
        ax2_rt = 0.0
        ax2_th = 0.0

    ax = hfac * ax2 + (1.0 - hfac) * ax1
    ax_hk = hfac * ax2_hk + (1.0 - hfac) * ax1_hk + hf_hk * (ax2 - ax1)
    ax_rt = hfac * ax2_rt + (1.0 - hfac) * ax1_rt
    ax_th = hfac * ax2_th + (1.0 - hfac) * ax1_th

    return ax, ax_hk, ax_th, ax_rt


def hkin(h, msq):
    hk = (h - 0.29 * msq) / (1.0 + 0.113 * msq)
    hk_h = 1.0 / (1.0 + 0.113 * msq)
    hk_msq = (-0.29 - 0.113 * hk) / (1.0 + 0.113 * msq)
    return hk, hk_h, hk_msq


def dil(hk, rt):
    if hk < 4.0:
        di = (0.00205 * (4.0 - hk) ** 5.5 + 0.207) / rt
        di_hk = (-0.00205 * 5.5 * (4.0 - hk) ** 4.5) / rt
    else:
        hkb = hk - 4.0
        den = 1.0 + 0.02 * hkb**2
        di = (-0.0016 * hkb**2 / den + 0.207) / rt
        di_hk = (-0.0016 * 2.0 * hkb * (1.0 / den - 0.02 * hkb**2 / den**2)) / rt
    di_rt = -di / rt
    return di, di_hk, di_rt


def dilw(hk, rt):
    msq = 0.0
    hs, hs_hk, hs_rt, _ = hsl(hk, rt, msq)

    rcd = 1.10 * (1.0 - 1.0 / hk) ** 2 / hk
    rcd_hk = -1.10 * (1.0 - 1.0 / hk) * 2.0 / hk**3 - rcd / hk

    di = 2.0 * rcd / (hs * rt)
    di_hk = 2.0 * rcd_hk / (hs * rt) - (di / hs) * hs_hk
    di_rt = -di / rt - (di / hs) * hs_rt
    return di, di_hk, di_rt


def hsl(hk, rt, msq):
    if hk < 4.35:
        tmp = hk - 4.35
        hs = 0.0111 * tmp**2 / (hk + 1.0) - 0.0278 * tmp**3 / (hk + 1.0) + 1.528 - 0.0002 * (tmp * hk) ** 2
        hs_hk = (
            0.0111 * (2.0 * tmp - tmp**2 / (hk + 1.0)) / (hk + 1.0)
            - 0.0278 * (3.0 * tmp**2 - tmp**3 / (hk + 1.0)) / (hk + 1.0)
            - 0.0002 * 2.0 * tmp * hk * (tmp + hk)
        )
    else:
        hs2 = 0.015
        hs = hs2 * (hk - 4.35) ** 2 / hk + 1.528
        hs_hk = hs2 * 2.0 * (hk - 4.35) / hk - hs2 * (hk - 4.35) ** 2 / hk**2

    hs_rt = 0.0
    hs_msq = 0.0
    return hs, hs_hk, hs_rt, hs_msq


def cfl(hk, rt, msq):
    if hk < 5.5:
        tmp = (5.5 - hk) ** 3 / (hk + 1.0)
        cf = (0.0727 * tmp - 0.07) / rt
        cf_hk = (-0.0727 * tmp * 3.0 / (5.5 - hk) - 0.0727 * tmp / (hk + 1.0)) / rt
    else:
        tmp = 1.0 - 1.0 / (hk - 4.5)
        cf = (0.015 * tmp**2 - 0.07) / rt
        cf_hk = (0.015 * tmp * 2.0 / (hk - 4.5) ** 2) / rt
    cf_rt = -cf / rt
    cf_msq = 0.0
    return cf, cf_hk, cf_rt, cf_msq


def dit(hs, us, cf, st):
    di = (0.5 * cf * us + st * st * (1.0 - us)) * 2.0 / hs
    di_hs = -(0.5 * cf * us + st * st * (1.0 - us)) * 2.0 / hs**2
    di_us = (0.5 * cf - st * st) * 2.0 / hs
    di_cf = (0.5 * us) * 2.0 / hs
    di_st = (2.0 * st * (1.0 - us)) * 2.0 / hs
    return di, di_hs, di_us, di_cf, di_st


def hst(hk, rt, msq):
    hsmin = 1.500
    dhsinf = 0.015

    if rt > 400.0:
        ho = 3.0 + 400.0 / rt
        ho_rt = -400.0 / rt**2
    else:
        ho = 4.0
        ho_rt = 0.0

    if rt > 200.0:
        rtz = rt
        rtz_rt = 1.0
    else:
        rtz = 200.0
        rtz_rt = 0.0

    if hk < ho:
        hr = (ho - hk) / (ho - 1.0)
        hr_hk = -1.0 / (ho - 1.0)
        hr_rt = (1.0 - hr) / (ho - 1.0) * ho_rt
        hs = (2.0 - hsmin - 4.0 / rtz) * hr**2 * 1.5 / (hk + 0.5) + hsmin + 4.0 / rtz
        hs_hk = (
            -(2.0 - hsmin - 4.0 / rtz) * hr**2 * 1.5 / (hk + 0.5) ** 2
            + (2.0 - hsmin - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_hk
        )
        hs_rt = (
            (2.0 - hsmin - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_rt
            + (hr**2 * 1.5 / (hk + 0.5) - 1.0) * 4.0 / rtz**2 * rtz_rt
        )
    else:
        grt = math.log(rtz)
        hdif = hk - ho
        rtmp = hk - ho + 4.0 / grt
        htmp = 0.007 * grt / rtmp**2 + dhsinf / hk
        htmp_hk = -0.014 * grt / rtmp**3 - dhsinf / hk**2
        htmp_rt = -0.014 * grt / rtmp**3 * (-ho_rt - 4.0 / grt**2 / rtz * rtz_rt) + 0.007 / rtmp**2 / rtz * rtz_rt
        hs = hdif**2 * htmp + hsmin + 4.0 / rtz
        hs_hk = hdif * 2.0 * htmp + hdif**2 * htmp_hk
        hs_rt = hdif**2 * htmp_rt - 4.0 / rtz**2 * rtz_rt + hdif * 2.0 * htmp * (-ho_rt)

    fm = 1.0 + 0.014 * msq
    hs = (hs + 0.028 * msq) / fm
    hs_hk = hs_hk / fm
    hs_rt = hs_rt / fm
    hs_msq = 0.028 / fm - 0.014 * hs / fm
    return hs, hs_hk, hs_rt, hs_msq


def cft(hk, rt, msq, cffac):
    gam = 1.4
    gm1 = gam - 1.0
    fc = math.sqrt(1.0 + 0.5 * gm1 * msq)
    grt = math.log(rt / fc)
    grt = max(grt, 3.0)

    gex = -1.74 - 0.31 * hk

    arg = -1.33 * hk
    arg = max(-20.0, arg)

    thk = math.tanh(4.0 - hk / 0.875)

    cfo = cffac * 0.3 * math.exp(arg) * (grt / 2.3026) ** gex
    cf = (cfo + 1.1e-4 * (thk - 1.0)) / fc
    cf_hk = (-1.33 * cfo - 0.31 * math.log(grt / 2.3026) * cfo - 1.1e-4 * (1.0 - thk**2) / 0.875) / fc
    cf_rt = gex * cfo / (fc * grt) / rt
    cf_msq = gex * cfo / (fc * grt) * (-0.25 * gm1 / fc**2) - 0.25 * gm1 * cf / fc**2
    return cf, cf_hk, cf_rt, cf_msq


def hct(hk, msq):
    hc = msq * (0.064 / (hk - 0.8) + 0.251)
    hc_hk = msq * (-0.064 / (hk - 0.8) ** 2)
    hc_msq = 0.064 / (hk - 0.8) + 0.251
    return hc, hc_hk, hc_msq
