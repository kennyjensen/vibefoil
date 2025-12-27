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
from python.xbl import mrchue, mrchdu
from python.xsolve import gauss as gauss_base


def max_abs_list_diff(a, b):
    if len(a) != len(b):
        return float("inf")
    return max(abs(ai - bi) for ai, bi in zip(a, b))


class TestMrchduParity(unittest.TestCase):
    def test_mrchdu(self):
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

        script = pathlib.Path(__file__).with_name("compare_mrchdu.mjs")
        tol = 1.0e-8

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
                mrchdu(ctx, bl)

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
                    mass_py = [0.0] * (max_nbl + 1)
                    tau_py = [0.0] * (max_nbl + 1)
                    dis_py = [0.0] * (max_nbl + 1)
                    ctq_py = [0.0] * (max_nbl + 1)
                    delt_py = [0.0] * (max_nbl + 1)
                    tstr_py = [0.0] * (max_nbl + 1)
                    ht_py = [0.0] * (max_nbl + 1)
                    for ibl in range(1, ctx.NBL[is_] + 1):
                        th_py[ibl] = ctx.THET[ibl][is_]
                        ds_py[ibl] = ctx.DSTR[ibl][is_]
                        ct_py[ibl] = ctx.CTAU[ibl][is_]
                        ue_py[ibl] = ctx.UEDG[ibl][is_]
                        mass_py[ibl] = ctx.MASS[ibl][is_]
                        tau_py[ibl] = ctx.TAU[ibl][is_]
                        dis_py[ibl] = ctx.DIS[ibl][is_]
                        ctq_py[ibl] = ctx.CTQ[ibl][is_]
                        delt_py[ibl] = ctx.DELT[ibl][is_]
                        tstr_py[ibl] = ctx.TSTR[ibl][is_]
                        ht_py[ibl] = ctx.HTARG[is_][ibl]

                    th_js = js_results["THET"][str(is_)]
                    ds_js = js_results["DSTR"][str(is_)]
                    ct_js = js_results["CTAU"][str(is_)]
                    ue_js = js_results["UEDG"][str(is_)]
                    mass_js = js_results["MASS"][str(is_)]
                    tau_js = js_results["TAU"][str(is_)]
                    dis_js = js_results["DIS"][str(is_)]
                    ctq_js = js_results["CTQ"][str(is_)]
                    delt_js = js_results["DELT"][str(is_)]
                    tstr_js = js_results["TSTR"][str(is_)]
                    ht_js = js_results["HTARG"][str(is_)]

                    self.assertLessEqual(max_abs_list_diff(th_py, th_js), tol)
                    self.assertLessEqual(max_abs_list_diff(ds_py, ds_js), tol)
                    self.assertLessEqual(max_abs_list_diff(ct_py, ct_js), tol)
                    self.assertLessEqual(max_abs_list_diff(ue_py, ue_js), tol)
                    self.assertLessEqual(max_abs_list_diff(mass_py, mass_js), tol)
                    self.assertLessEqual(max_abs_list_diff(tau_py, tau_js), tol)
                    self.assertLessEqual(max_abs_list_diff(dis_py, dis_js), tol)
                    self.assertLessEqual(max_abs_list_diff(ctq_py, ctq_js), tol)
                    self.assertLessEqual(max_abs_list_diff(delt_py, delt_js), tol)
                    self.assertLessEqual(max_abs_list_diff(tstr_py, tstr_js), tol)
                    self.assertLessEqual(max_abs_list_diff(ht_py, ht_js), tol)


if __name__ == "__main__":
    unittest.main()
