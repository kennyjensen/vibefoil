import json
import math
import pathlib
import shutil
import subprocess
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import python.xbl as xbl_mod
from python.xbl import XBlState, blpini, setbl, iblsys
from python.xpanel import (
    xywake,
    qwcalc,
    qiset,
    stfind,
    iblpan,
    xicalc,
    uicalc,
    qdcalc,
    qvfue,
    gamqv,
    stmove,
)
from python.xsolve import gauss as gauss_base, blsolv
from python.tests.test_compare_viscal import build_viscal_context
from python.tests.test_compare_viscal_subfuncs import build_payload_ctx


def compute_update_metrics(ctx):
    max_nbl = max(ctx.NBL[1], ctx.NBL[2])
    unew = [[0.0] * 3 for _ in range(max_nbl + 2)]
    u_ac = [[0.0] * 3 for _ in range(max_nbl + 2)]
    qnew = [0.0] * (ctx.N + ctx.NW + 1)
    q_ac = [0.0] * (ctx.N + ctx.NW + 1)

    dclmin = -0.5
    dclmax = 0.5
    if ctx.MATYP != 1:
        dclmin = max(-0.5, -0.9 * ctx.CL)

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

            uinv_ac = 0.0 if ctx.LALFA else ctx.UINV_A[ibl][is_]
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
    dac = 0.0
    if ctx.LALFA:
        dac = (clnew - ctx.CL) / (1.0 - cl_ac - cl_ms * 2.0 * ctx.MINF * ctx.MINF_CL)
        if rlx * dac > dclmax:
            rlx = dclmax / dac
        if rlx * dac < dclmin:
            rlx = dclmin / dac
    else:
        dalmax = 0.5 * ctx.DTOR
        dalmin = -0.5 * ctx.DTOR
        dac = (clnew - ctx.CLSPEC) / (0.0 - cl_ac - cl_a)
        if rlx * dac > dalmax:
            rlx = dalmax / dac
        if rlx * dac < dalmin:
            rlx = dalmin / dac

    rlx_dac = rlx
    rmsbl = 0.0
    rmxbl = 0.0
    vmxbl = " "
    imxbl = 0
    ismxbl = 0
    dn_max = [0.0, 0.0, 0.0, 0.0]

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

            rmsbl = rmsbl + dn1**2 + dn2**2 + dn3**2 + dn4**2

            rdn1 = rlx * dn1
            if abs(dn1) > abs(rmxbl):
                rmxbl = dn1
                vmxbl = "n" if ibl < ctx.ITRAN[is_] else "C"
                imxbl = ibl
                ismxbl = is_
                dn_max = [dn1, dn2, dn3, dn4]
            if rdn1 > dhi:
                rlx = dhi / dn1
            if rdn1 < dlo:
                rlx = dlo / dn1

            rdn2 = rlx * dn2
            if abs(dn2) > abs(rmxbl):
                rmxbl = dn2
                vmxbl = "T"
                imxbl = ibl
                ismxbl = is_
                dn_max = [dn1, dn2, dn3, dn4]
            if rdn2 > dhi:
                rlx = dhi / dn2
            if rdn2 < dlo:
                rlx = dlo / dn2

            rdn3 = rlx * dn3
            if abs(dn3) > abs(rmxbl):
                rmxbl = dn3
                vmxbl = "D"
                imxbl = ibl
                ismxbl = is_
                dn_max = [dn1, dn2, dn3, dn4]
            if rdn3 > dhi:
                rlx = dhi / dn3
            if rdn3 < dlo:
                rlx = dlo / dn3

            rdn4 = rlx * dn4
            if abs(dn4) > abs(rmxbl):
                rmxbl = duedg
                vmxbl = "U"
                imxbl = ibl
                ismxbl = is_
                dn_max = [dn1, dn2, dn3, dn4]
            if rdn4 > dhi:
                rlx = dhi / dn4
            if rdn4 < dlo:
                rlx = dlo / dn4

    rmsbl = math.sqrt(rmsbl / (4.0 * float(ctx.NBL[1] + ctx.NBL[2])))

    return {
        "clnew": clnew,
        "clAc": cl_ac,
        "clMs": cl_ms,
        "dac": dac,
        "rlxDac": rlx_dac,
        "rlxFinal": rlx,
        "rmsbl": rmsbl,
        "rmxbl": rmxbl,
        "vmxbl": vmxbl,
        "imxbl": imxbl,
        "ismxbl": ismxbl,
        "dnMax": dn_max,
    }


class TestUpdateRlxCase(unittest.TestCase):
    def test_update_rlx_case(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        def gauss1(nsiz, nn, z, r, nrhs):
            rmat = [[0.0] * (nrhs + 1) for _ in range(nn + 1)]
            for i in range(1, nn + 1):
                rmat[i][1] = r[i]
            gauss_base(nsiz, nn, z, rmat, nrhs)
            for i in range(1, nn + 1):
                r[i] = rmat[i][1]

        xbl_mod.gauss = gauss1

        case = {"ides": 2412, "minf": 0.0, "reinf": 1.0e6, "alfa": -10.0 * math.pi / 180.0, "waklen": 1.0}
        ctx = build_viscal_context(case["ides"], case["minf"], case["reinf"], case["alfa"], case["waklen"])
        cosa = math.cos(ctx.ALFA)
        sina = math.sin(ctx.ALFA)
        for i in range(1, ctx.N + 1):
            ctx.GAM[i] = cosa * ctx.GAMU[i][1] + sina * ctx.GAMU[i][2]
        bl = XBlState()
        blpini(bl)

        xywake(ctx)
        qwcalc(ctx)
        qiset(ctx)
        stfind(ctx)
        iblpan(ctx)
        xicalc(ctx)
        iblsys(ctx)
        uicalc(ctx)
        for is_ in (1, 2):
            for ibl in range(1, ctx.NBL[is_] + 1):
                ctx.UEDG[ibl][is_] = ctx.UINV[ibl][is_]
        qdcalc(ctx)
        qvfue(ctx)
        gamqv(ctx)
        stmove(ctx)
        setbl(ctx, bl)
        blsolv(ctx)

        py_results = compute_update_metrics(ctx)

        payload = build_payload_ctx(ctx, case["waklen"], 9.0)
        script = pathlib.Path(__file__).with_name("compare_update_rlx_case.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        tol = 1.0e-9
        for key in ("clnew", "clAc", "clMs", "dac", "rlxDac", "rlxFinal", "rmsbl", "rmxbl"):
            self.assertLessEqual(abs(py_results[key] - js_results[key]), tol, msg=key)
        self.assertEqual(py_results["vmxbl"], js_results["vmxbl"])
        self.assertEqual(py_results["imxbl"], js_results["imxbl"])
        self.assertEqual(py_results["ismxbl"], js_results["ismxbl"])
        for i in range(4):
            self.assertLessEqual(abs(py_results["dnMax"][i] - js_results["dnMax"][i]), tol, msg=f"dnMax[{i}]")
