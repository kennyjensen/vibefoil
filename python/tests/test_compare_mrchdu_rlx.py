import json
import math
import pathlib
import shutil
import subprocess
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from python.tests.test_compare_mrchue import build_context, build_bl, build_payload
import python.xbl as xbl_mod
from python.xbl import mrchue, xifset
from python.xblsys import blprv, blkin, blsys, tesys, trchek, hkin
from python.xsolve import gauss as gauss_base


def compute_first_rlx(ctx, bl, is_):
    bl.AMCRIT = ctx.ACRIT[is_]
    xifset(ctx, bl, is_)

    ibl = 2
    bl.SIMI = ibl == 2
    bl.WAKE = ibl > ctx.IBLTE[is_]

    itrold = ctx.ITRAN[is_]
    bl.TRAN = False
    bl.TURB = False
    ctx.ITRAN[is_] = ctx.IBLTE[is_]

    xsi = ctx.XSSI[ibl][is_]
    uei = ctx.UEDG[ibl][is_]
    thi = ctx.THET[ibl][is_]
    dsi = ctx.DSTR[ibl][is_]

    ami = 0.0
    cti = 0.03
    if ibl < itrold:
        ami = ctx.CTAU[ibl][is_]
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

    ueref = bl.U2
    hkref = bl.HK2

    if ibl < ctx.ITRAN[is_] and ibl >= itrold:
        uem = ctx.UEDG[ibl - 1][is_]
        dsm = ctx.DSTR[ibl - 1][is_]
        thm = ctx.THET[ibl - 1][is_]
        msq = uem * uem * bl.HSTINV / (bl.GM1BL * (1.0 - 0.5 * uem * uem * bl.HSTINV))
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
        vtmp = [[0.0] * 5 for _ in range(5)]
        vztmp = [0.0] * 5
        for k in range(1, 5):
            vztmp[k] = bl.VSREZ[k]
            for l in range(1, 5):
                vtmp[k][l] = bl.VS2[k][l]

        vtmp[4][1] = 0.0
        vtmp[4][2] = bl.HK2_T2
        vtmp[4][3] = bl.HK2_D2
        vtmp[4][4] = bl.HK2_U2 * bl.U2_UEI
        vztmp[4] = 1.0

        xbl_mod.gauss(4, 4, vtmp, vztmp, 1)

        sennew = 1000.0 * vztmp[4] * hkref / ueref
        sens = sennew

        bl.VS2[4][1] = 0.0
        bl.VS2[4][2] = bl.HK2_T2 * hkref
        bl.VS2[4][3] = bl.HK2_D2 * hkref
        bl.VS2[4][4] = (bl.HK2_U2 * hkref + sens / ueref) * bl.U2_UEI
        bl.VSREZ[4] = -(hkref**2) * (bl.HK2 / hkref - 1.0) - sens * (bl.U2 / ueref - 1.0)

    xbl_mod.gauss(4, 4, bl.VS2, bl.VSREZ, 1)

    dmax = max(abs(bl.VSREZ[2] / thi), abs(bl.VSREZ[3] / dsi), abs(bl.VSREZ[4] / uei))
    if ibl >= ctx.ITRAN[is_]:
        dmax = max(dmax, abs(bl.VSREZ[1] / (10.0 * cti)))

    rlx = 1.0
    if dmax > 0.3:
        rlx = 0.3 / dmax

    return {
        "ibl": ibl,
        "itrold": itrold,
        "simi": bl.SIMI,
        "wake": bl.WAKE,
        "xsi": xsi,
        "uei": uei,
        "thi": thi,
        "dsi": dsi,
        "cti": cti,
        "ami": ami,
        "ueref": ueref,
        "hkref": hkref,
        "dmax": dmax,
        "rlx": rlx,
        "VSREZ": [bl.VSREZ[1], bl.VSREZ[2], bl.VSREZ[3], bl.VSREZ[4]],
    }


class TestMrchduRlxParity(unittest.TestCase):
    def test_mrchdu_rlx_first_iteration(self):
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

        case = {
            "ides": 2412,
            "minf": 0.0,
            "reinf": 1.0e6,
            "alfa": -10.0 * math.pi / 180.0,
            "waklen": 1.0,
        }
        ctx = build_context(
            case["ides"],
            case["minf"],
            case["reinf"],
            case["alfa"],
            case["waklen"],
        )
        bl = build_bl(ctx)

        payload = build_payload(ctx, bl)
        script = pathlib.Path(__file__).with_name("compare_mrchdu_rlx.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        mrchue(ctx, bl)

        tol = 1.0e-9
        for is_ in (1, 2):
            py_res = compute_first_rlx(ctx, bl, is_)
            js_res = js_results[str(is_)]

            for key in ("xsi", "uei", "thi", "dsi", "cti", "ami", "ueref", "hkref", "dmax", "rlx"):
                self.assertLessEqual(abs(py_res[key] - js_res[key]), tol, msg=f"{is_} {key}")

            for i, val in enumerate(py_res["VSREZ"]):
                self.assertLessEqual(abs(val - js_res["VSREZ"][i]), tol, msg=f"{is_} VSREZ[{i+1}]")


if __name__ == "__main__":
    unittest.main()
