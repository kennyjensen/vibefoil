import json
import pathlib
import shutil
import subprocess
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from python.naca import naca4, naca5
from python.xbl import XFoilState
from python.xfoil import pangen


def max_abs_diff(a, b):
    if len(a) != len(b):
        return float("inf")
    return max(abs(ai - bi) for ai, bi in zip(a, b))


class TestNacaPangenParity(unittest.TestCase):
    def test_naca_and_pangen(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        cases = [
            {"ides": 12, "nside": 123},
            {"ides": 2412, "nside": 123},
            {"ides": 23012, "nside": 123},
            {"ides": 25012, "nside": 123},
        ]

        params = {
            "npan": 160,
            "cvpar": 1.0,
            "cterat": 0.15,
            "ctrrat": 0.2,
            "xsref1": 1.0,
            "xsref2": 1.0,
            "xpref1": 1.0,
            "xpref2": 1.0,
        }

        payload = {"cases": [{**case, "params": params} for case in cases]}
        script = pathlib.Path(__file__).with_name("compare_naca_pangen.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = {item["ides"]: item for item in json.loads(proc.stdout)["results"]}

        tol = 1.0e-10

        for case in cases:
            ides = case["ides"]
            js_case = js_results[ides]

            ctx = XFoilState()
            ctx.NPAN = params["npan"]
            ctx.CVPAR = params["cvpar"]
            ctx.CTERAT = params["cterat"]
            ctx.CTRRAT = params["ctrrat"]
            ctx.XSREF1 = params["xsref1"]
            ctx.XSREF2 = params["xsref2"]
            ctx.XPREF1 = params["xpref1"]
            ctx.XPREF2 = params["xpref2"]

            nside = case["nside"]
            xx = [0.0] * (nside + 1)
            yt = [0.0] * (nside + 1)
            yc = [0.0] * (nside + 1)

            if ides <= 9999:
                nb, _ = naca4(ides, xx, yt, yc, nside, ctx.XB, ctx.YB)
            else:
                nb, _ = naca5(ides, xx, yt, yc, nside, ctx.XB, ctx.YB)

            self.assertEqual(nb, js_case["naca"]["nb"])

            xb_py = ctx.XB[: nb + 1]
            yb_py = ctx.YB[: nb + 1]
            xb_js = js_case["naca"]["xb"]
            yb_js = js_case["naca"]["yb"]
            self.assertLessEqual(max_abs_diff(xb_py, xb_js), tol)
            self.assertLessEqual(max_abs_diff(yb_py, yb_js), tol)

            ctx.NB = nb
            pangen(ctx, False)

            self.assertEqual(ctx.N, js_case["pangen"]["n"])

            x_py = ctx.X[: ctx.N + 1]
            y_py = ctx.Y[: ctx.N + 1]
            s_py = ctx.S[: ctx.N + 1]
            x_js = js_case["pangen"]["x"]
            y_js = js_case["pangen"]["y"]
            s_js = js_case["pangen"]["s"]
            self.assertLessEqual(max_abs_diff(x_py, x_js), tol)
            self.assertLessEqual(max_abs_diff(y_py, y_js), tol)
            self.assertLessEqual(max_abs_diff(s_py, s_js), tol)


if __name__ == "__main__":
    unittest.main()
