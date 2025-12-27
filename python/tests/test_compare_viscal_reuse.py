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
from python.tests.test_compare_viscal import build_payload, build_viscal_context
from python.xoper import viscal
from python.xsolve import gauss as gauss_base


class TestViscalReuseParity(unittest.TestCase):
    def test_viscal_reuse_two_calls(self):
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

        script = pathlib.Path(__file__).with_name("compare_viscal_reuse.mjs")
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
                cosa = math.cos(ctx.ALFA)
                sina = math.sin(ctx.ALFA)
                for i in range(1, ctx.N + 1):
                    ctx.GAM[i] = cosa * ctx.GAMU[i][1] + sina * ctx.GAMU[i][2]
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

                viscal(ctx, bl, 5)
                viscal(ctx, bl, 5)

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
