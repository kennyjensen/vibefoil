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
from python.xbl import XBlState, blpini
from python.xoper import viscal
from python.xsolve import gauss as gauss_base
from python.tests.test_compare_viscal import build_viscal_context, build_payload


class TestViscalSummaryCase(unittest.TestCase):
    def test_viscal_summary_case(self):
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

        payload = build_payload(ctx, case["minf"], case["reinf"], case["alfa"], case["waklen"])
        script = pathlib.Path(__file__).with_name("compare_viscal_summary.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        viscal(ctx, bl, 1)

        expected = {
            "XOCTR1": 1.0000,
            "XOCTR2": 0.0196,
            "ITRAN1": 68,
            "ITRAN2": 28,
            "TFORCE1": True,
            "TFORCE2": False,
            "RMSBL": 0.2309,
            "RMXBL": 1.233,
            "IMXBL": 98,
            "ISMXBL": 2,
            "VMXBL": "D",
            "RLX": 0.810,
            "ALFA_DEG": -10.000,
            "CL": -0.9933,
            "CM": -0.0296,
            "CD": 0.01369,
            "CDF": 0.00531,
            "CDP": 0.00838,
        }

        tol_transition = 5.0e-4
        tol_summary = 5.0e-4
        tol_coeff = 5.0e-4
        tol_rlx = 5.0e-4

        self.assertLessEqual(abs(ctx.XOCTR[1] - expected["XOCTR1"]), tol_transition, "XOCTR1")
        self.assertLessEqual(abs(ctx.XOCTR[2] - expected["XOCTR2"]), tol_transition, "XOCTR2")
        self.assertEqual(ctx.ITRAN[1], expected["ITRAN1"], "ITRAN1")
        self.assertEqual(ctx.ITRAN[2], expected["ITRAN2"], "ITRAN2")
        self.assertEqual(bool(ctx.TFORCE[1]), expected["TFORCE1"], "TFORCE1")
        self.assertEqual(bool(ctx.TFORCE[2]), expected["TFORCE2"], "TFORCE2")

        self.assertLessEqual(abs(ctx.RMSBL - expected["RMSBL"]), tol_summary, "RMSBL")
        self.assertLessEqual(abs(ctx.RMXBL - expected["RMXBL"]), tol_summary, "RMXBL")
        self.assertEqual(ctx.IMXBL, expected["IMXBL"], "IMXBL")
        self.assertEqual(ctx.ISMXBL, expected["ISMXBL"], "ISMXBL")
        self.assertEqual(ctx.VMXBL.strip(), expected["VMXBL"], "VMXBL")
        self.assertLessEqual(abs(ctx.RLX - expected["RLX"]), tol_rlx, "RLX")

        self.assertLessEqual(abs(ctx.ALFA / ctx.DTOR - expected["ALFA_DEG"]), tol_coeff, "ALFA")
        self.assertLessEqual(abs(ctx.CL - expected["CL"]), tol_coeff, "CL")
        self.assertLessEqual(abs(ctx.CM - expected["CM"]), tol_coeff, "CM")
        self.assertLessEqual(abs(ctx.CD - expected["CD"]), tol_coeff, "CD")
        self.assertLessEqual(abs(ctx.CDF - expected["CDF"]), tol_coeff, "CDF")
        self.assertLessEqual(abs(ctx.CD - ctx.CDF - expected["CDP"]), tol_coeff, "CDP")

        tol_js = 1.0e-6
        self.assertLessEqual(abs(js_results["ALFA_DEG"] - ctx.ALFA / ctx.DTOR), tol_js, "JS ALFA")
        self.assertLessEqual(abs(js_results["CL"] - ctx.CL), tol_js, "JS CL")
        self.assertLessEqual(abs(js_results["CM"] - ctx.CM), tol_js, "JS CM")
        self.assertLessEqual(abs(js_results["CD"] - ctx.CD), tol_js, "JS CD")
        self.assertLessEqual(abs(js_results["CDF"] - ctx.CDF), tol_js, "JS CDF")
        self.assertLessEqual(abs(js_results["CDP"] - (ctx.CD - ctx.CDF)), tol_js, "JS CDP")
        self.assertLessEqual(abs(js_results["RMSBL"] - ctx.RMSBL), tol_js, "JS RMSBL")
        self.assertLessEqual(abs(js_results["RMXBL"] - ctx.RMXBL), tol_js, "JS RMXBL")
        self.assertEqual(js_results["IMXBL"], ctx.IMXBL, "JS IMXBL")
        self.assertEqual(js_results["ISMXBL"], ctx.ISMXBL, "JS ISMXBL")
        self.assertEqual(js_results["VMXBL"], ctx.VMXBL, "JS VMXBL")
        self.assertLessEqual(abs(js_results["RLX"] - ctx.RLX), tol_js, "JS RLX")
        self.assertLessEqual(abs(js_results["XOCTR"][1] - ctx.XOCTR[1]), tol_js, "JS XOCTR1")
        self.assertLessEqual(abs(js_results["XOCTR"][2] - ctx.XOCTR[2]), tol_js, "JS XOCTR2")
        self.assertEqual(js_results["ITRAN"][1], ctx.ITRAN[1], "JS ITRAN1")
        self.assertEqual(js_results["ITRAN"][2], ctx.ITRAN[2], "JS ITRAN2")
        self.assertEqual(bool(js_results["TFORCE"][1]), bool(ctx.TFORCE[1]), "JS TFORCE1")
        self.assertEqual(bool(js_results["TFORCE"][2]), bool(ctx.TFORCE[2]), "JS TFORCE2")


if __name__ == "__main__":
    unittest.main()
