import copy
import math
import json
import pathlib
import shutil
import subprocess
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from python.xbl import XFoilState, XBlState, iblsys, mrchue, blpini
import python.xbl as xbl_mod
from python.xsolve import gauss as gauss_base
from python.xfoil import comset, naca
from python.xpanel import ggcalc, xywake, qwcalc, qiset, stfind, iblpan, xicalc, uicalc


def build_context(ides, minf, reinf, alfa, waklen):
    ctx = XFoilState()

    ctx.NPAN = 160
    ctx.CVPAR = 1.0
    ctx.CTERAT = 0.15
    ctx.CTRRAT = 0.2
    ctx.XSREF1 = 1.0
    ctx.XSREF2 = 1.0
    ctx.XPREF1 = 1.0
    ctx.XPREF2 = 1.0

    ctx.WAKLEN = waklen
    ctx.ALFA = alfa
    ctx.ADEG = alfa / ctx.DTOR
    ctx.QINF = 1.0
    ctx.MINF = minf
    ctx.MINF1 = minf
    ctx.REINF = reinf
    ctx.REINF1 = reinf
    ctx.LALFA = True

    ctx.ACRIT[1] = 9.0
    ctx.ACRIT[2] = 9.0
    ctx.XSTRIP[1] = 1.0
    ctx.XSTRIP[2] = 1.0

    naca(ctx, ides)
    comset(ctx)
    ggcalc(ctx)
    xywake(ctx)
    qwcalc(ctx)
    qiset(ctx)
    for i in range(1, ctx.N + 1):
        ctx.GAM[i] = 1.0 if i <= ctx.N // 2 else -1.0
    stfind(ctx)
    iblpan(ctx)
    xicalc(ctx)
    iblsys(ctx)
    uicalc(ctx)

    for is_ in range(1, 3):
        for ibl in range(1, ctx.NBL[is_] + 1):
            ctx.UEDG[ibl][is_] = abs(ctx.UINV[ibl][is_])

    return ctx


def build_bl(ctx):
    bl = XBlState()
    blpini(bl)

    bl.GAMBL = ctx.GAMMA
    bl.GM1BL = ctx.GAMM1
    bl.QINFBL = ctx.QINF
    bl.TKBL = ctx.TKLAM
    bl.TKBL_MS = ctx.TKL_MSQ

    bl.RSTBL = (1.0 + 0.5 * bl.GM1BL * ctx.MINF**2) ** (1.0 / bl.GM1BL)
    bl.RSTBL_MS = 0.5 * bl.RSTBL / (1.0 + 0.5 * bl.GM1BL * ctx.MINF**2)

    bl.HSTINV = bl.GM1BL * (ctx.MINF / bl.QINFBL) ** 2 / (1.0 + 0.5 * bl.GM1BL * ctx.MINF**2)
    bl.HSTINV_MS = (
        bl.GM1BL * (1.0 / bl.QINFBL) ** 2 / (1.0 + 0.5 * bl.GM1BL * ctx.MINF**2)
        - 0.5 * bl.GM1BL * bl.HSTINV / (1.0 + 0.5 * bl.GM1BL * ctx.MINF**2)
    )

    herat = 1.0 - 0.5 * bl.QINFBL**2 * bl.HSTINV
    herat_ms = -0.5 * bl.QINFBL**2 * bl.HSTINV_MS

    bl.HVRAT = ctx.HVRAT
    bl.REYBL = ctx.REINF * (herat ** 1.5) * (1.0 + bl.HVRAT) / (herat + bl.HVRAT)
    bl.REYBL_RE = (herat ** 1.5) * (1.0 + bl.HVRAT) / (herat + bl.HVRAT)
    bl.REYBL_MS = bl.REYBL * (1.5 / herat - 1.0 / (herat + bl.HVRAT)) * herat_ms

    bl.IDAMPV = ctx.IDAMP
    bl.DWTE = ctx.WGAP[1]

    return bl


def serialize_side_matrix(mat, max_nbl):
    out = []
    for ibl in range(0, max_nbl + 1):
        row = [0.0, 0.0, 0.0]
        row[1] = mat[ibl][1]
        row[2] = mat[ibl][2]
        out.append(row)
    return out


def serialize_vector(vec):
    return [vec[0], vec[1], vec[2]]


def build_payload(ctx, bl):
    max_nbl = max(ctx.NBL[1], ctx.NBL[2])
    nw = ctx.NBL[2] - ctx.IBLTE[2]
    wgap = [0.0] * (nw + 1)
    for i in range(1, nw + 1):
        wgap[i] = ctx.WGAP[i]

    payload = {
        "ctx": {
            "NBL": serialize_vector(ctx.NBL),
            "IBLTE": serialize_vector(ctx.IBLTE),
            "ITRAN": serialize_vector(ctx.ITRAN),
            "XSSITR": serialize_vector(ctx.XSSITR),
            "XSTRIP": serialize_vector(ctx.XSTRIP),
            "ACRIT": serialize_vector(ctx.ACRIT),
            "WGAP": wgap,
            "TFORCE": [False, False, False],
            "XSSI": serialize_side_matrix(ctx.XSSI, max_nbl),
            "UEDG": serialize_side_matrix(ctx.UEDG, max_nbl),
            "THET": serialize_side_matrix(ctx.THET, max_nbl),
            "DSTR": serialize_side_matrix(ctx.DSTR, max_nbl),
            "CTAU": serialize_side_matrix(ctx.CTAU, max_nbl),
            "MASS": serialize_side_matrix(ctx.MASS, max_nbl),
            "TAU": serialize_side_matrix(ctx.TAU, max_nbl),
            "DIS": serialize_side_matrix(ctx.DIS, max_nbl),
            "CTQ": serialize_side_matrix(ctx.CTQ, max_nbl),
            "DELT": serialize_side_matrix(ctx.DELT, max_nbl),
            "TSTR": serialize_side_matrix(ctx.TSTR, max_nbl),
            "ANTE": ctx.ANTE,
        },
        "bl": {
            "GAMBL": bl.GAMBL,
            "GM1BL": bl.GM1BL,
            "QINFBL": bl.QINFBL,
            "TKBL": bl.TKBL,
            "TKBL_MS": bl.TKBL_MS,
            "RSTBL": bl.RSTBL,
            "RSTBL_MS": bl.RSTBL_MS,
            "HSTINV": bl.HSTINV,
            "HSTINV_MS": bl.HSTINV_MS,
            "HVRAT": bl.HVRAT,
            "REYBL": bl.REYBL,
            "REYBL_RE": bl.REYBL_RE,
            "REYBL_MS": bl.REYBL_MS,
            "IDAMPV": bl.IDAMPV,
            "DWTE": bl.DWTE,
        },
    }
    return payload


def max_abs_list_diff(a, b):
    if len(a) != len(b):
        return float("inf")
    return max(abs(ai - bi) for ai, bi in zip(a, b))


class TestMrchueParity(unittest.TestCase):
    def test_mrchue(self):
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

        cases = [
            {"ides": 12, "minf": 0.0, "reinf": 1.0e6, "alfa": 0.0, "waklen": 1.0},
            {"ides": 2412, "minf": 0.1, "reinf": 3.0e6, "alfa": 0.0, "waklen": 1.0},
            {"ides": 23012, "minf": 0.2, "reinf": 5.0e6, "alfa": 0.0, "waklen": 2.0},
        ]

        script = pathlib.Path(__file__).with_name("compare_mrchue.mjs")
        tol = 1.0e-9

        for case in cases:
            with self.subTest(case=case):
                ctx = build_context(
                    case["ides"],
                    case["minf"],
                    case["reinf"],
                    case["alfa"],
                    case["waklen"],
                )
                bl = build_bl(ctx)

                payload = build_payload(ctx, bl)
                proc = subprocess.run(
                    ["node", str(script)],
                    input=json.dumps(payload),
                    text=True,
                    capture_output=True,
                    check=True,
                )
                js_results = json.loads(proc.stdout)["results"]

                mrchue(ctx, bl)

                self.assertEqual(ctx.NBL[1], js_results["NBL"][1])
                self.assertEqual(ctx.NBL[2], js_results["NBL"][2])
                self.assertEqual(ctx.IBLTE[1], js_results["IBLTE"][1])
                self.assertEqual(ctx.IBLTE[2], js_results["IBLTE"][2])
                self.assertEqual(ctx.ITRAN[1], js_results["ITRAN"][1])
                self.assertEqual(ctx.ITRAN[2], js_results["ITRAN"][2])

                self.assertLessEqual(abs(ctx.XSSITR[1] - js_results["XSSITR"][1]), tol)
                self.assertLessEqual(abs(ctx.XSSITR[2] - js_results["XSSITR"][2]), tol)

                max_nbl = max(ctx.NBL[1], ctx.NBL[2])

                for is_ in (1, 2):
                    th_py = [0.0] * (max_nbl + 1)
                    ds_py = [0.0] * (max_nbl + 1)
                    ct_py = [0.0] * (max_nbl + 1)
                    ue_py = [0.0] * (max_nbl + 1)
                    ht_py = [0.0] * (max_nbl + 1)
                    for ibl in range(1, ctx.NBL[is_] + 1):
                        th_py[ibl] = ctx.THET[ibl][is_]
                        ds_py[ibl] = ctx.DSTR[ibl][is_]
                        ct_py[ibl] = ctx.CTAU[ibl][is_]
                        ue_py[ibl] = ctx.UEDG[ibl][is_]
                        ht_py[ibl] = ctx.HTARG[is_][ibl]

                    th_js = js_results["THET"][str(is_)]
                    ds_js = js_results["DSTR"][str(is_)]
                    ct_js = js_results["CTAU"][str(is_)]
                    ue_js = js_results["UEDG"][str(is_)]
                    ht_js = js_results["HTARG"][str(is_)]

                    self.assertLessEqual(max_abs_list_diff(th_py, th_js), tol)
                    self.assertLessEqual(max_abs_list_diff(ds_py, ds_js), tol)
                    self.assertLessEqual(max_abs_list_diff(ct_py, ct_js), tol)
                    self.assertLessEqual(max_abs_list_diff(ue_py, ue_js), tol)
                    self.assertLessEqual(max_abs_list_diff(ht_py, ht_js), tol)


if __name__ == "__main__":
    unittest.main()
