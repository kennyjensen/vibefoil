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
from python.xbl import XFoilState, XBlState, blpini
from python.xfoil import comset, naca
from python.xpanel import ggcalc
from python.xoper import viscal
from python.xsolve import gauss as gauss_base


def build_viscal_context(ides, minf, reinf, alfa, waklen):
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
    ctx.VACCEL = 0.01

    ctx.ACRIT[1] = 9.0
    ctx.ACRIT[2] = 9.0
    ctx.XSTRIP[1] = 1.0
    ctx.XSTRIP[2] = 1.0

    naca(ctx, ides)
    comset(ctx)
    ggcalc(ctx)
    for i in range(1, ctx.N + 1):
        ctx.GAM[i] = 1.0 if i <= ctx.N // 2 else -1.0

    return ctx


def build_payload(ctx, minf, reinf, alfa, waklen):
    xb = [ctx.XB[i] for i in range(1, ctx.NB + 1)]
    yb = [ctx.YB[i] for i in range(1, ctx.NB + 1)]
    params = {
        "npan": ctx.NPAN,
        "cvpar": ctx.CVPAR,
        "cterat": ctx.CTERAT,
        "ctrrat": ctx.CTRRAT,
        "xsref1": ctx.XSREF1,
        "xsref2": ctx.XSREF2,
        "xpref1": ctx.XPREF1,
        "xpref2": ctx.XPREF2,
    }
    return {
        "xb": xb,
        "yb": yb,
        "nb": ctx.NB,
        "params": params,
        "alphaRad": alfa,
        "reinf": reinf,
        "minf": minf,
        "waklen": waklen,
        "ncrit": 9.0,
    }


def max_abs_list_diff(a, b):
    if len(a) != len(b):
        return float("inf")
    return max(abs(ai - bi) for ai, bi in zip(a, b))


class TestViscalParity(unittest.TestCase):
    def test_viscal_one_iter(self):
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

        script = pathlib.Path(__file__).with_name("compare_viscal.mjs")
        tol_coeff = 5.0e-1
        tol_arr = 5.0e-2

        for case in cases:
            with self.subTest(case=case):
                ctx = build_viscal_context(
                    case["ides"],
                    case["minf"],
                    case["reinf"],
                    case["alfa"],
                    case["waklen"],
                )
                bl = XBlState()
                blpini(bl)

                payload = build_payload(ctx, case["minf"], case["reinf"], case["alfa"], case["waklen"])
                proc = subprocess.run(
                    ["node", str(script)],
                    input=json.dumps(payload),
                    text=True,
                    capture_output=True,
                    check=True,
                )
                js_results = json.loads(proc.stdout)["results"]

                viscal(ctx, bl, 1)

                self.assertLessEqual(abs(ctx.CL - js_results["CL"]), tol_coeff)
                self.assertLessEqual(abs(ctx.CM - js_results["CM"]), tol_coeff)
                self.assertLessEqual(abs(ctx.CD - js_results["CD"]), tol_coeff)
                self.assertLessEqual(abs(ctx.CDF - js_results["CDF"]), tol_coeff)
                self.assertLessEqual(abs(ctx.CDP - js_results["CDP"]), tol_coeff)

                self.assertEqual(ctx.NBL[1], js_results["NBL"][1])
                self.assertEqual(ctx.NBL[2], js_results["NBL"][2])
                self.assertEqual(ctx.IBLTE[1], js_results["IBLTE"][1])
                self.assertEqual(ctx.IBLTE[2], js_results["IBLTE"][2])
                self.assertEqual(ctx.ITRAN[1], js_results["ITRAN"][1])
                self.assertEqual(ctx.ITRAN[2], js_results["ITRAN"][2])
                self.assertLessEqual(abs(ctx.XSSITR[1] - js_results["XSSITR"][1]), tol_arr)
                self.assertLessEqual(abs(ctx.XSSITR[2] - js_results["XSSITR"][2]), tol_arr)


if __name__ == "__main__":
    unittest.main()
