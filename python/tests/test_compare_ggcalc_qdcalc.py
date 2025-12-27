import json
import math
import pathlib
import shutil
import subprocess
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from python.xbl import XFoilState
from python.xfoil import comset, naca, pangen
from python.xpanel import ggcalc, qdcalc, xywake


def build_context(ides, params, alpha_rad, minf, waklen, nside, run_ggcalc_qdcalc=True):
    ctx = XFoilState()
    ctx.NPAN = params["npan"]
    ctx.CVPAR = params["cvpar"]
    ctx.CTERAT = params["cterat"]
    ctx.CTRRAT = params["ctrrat"]
    ctx.XSREF1 = params["xsref1"]
    ctx.XSREF2 = params["xsref2"]
    ctx.XPREF1 = params["xpref1"]
    ctx.XPREF2 = params["xpref2"]

    ctx.ALFA = alpha_rad
    ctx.ADEG = alpha_rad / ctx.DTOR
    ctx.MINF = minf
    ctx.MINF1 = minf
    ctx.QINF = 1.0
    ctx.WAKLEN = waklen

    naca(ctx, ides)
    pangen(ctx, False)
    comset(ctx)

    for i in range(1, ctx.N + 1):
        ctx.GAM[i] = 1.0 if i <= ctx.N // 2 else -1.0

    xywake(ctx)
    if run_ggcalc_qdcalc:
        ggcalc(ctx)
        qdcalc(ctx)

    return ctx


def build_payload(ctx):
    total = ctx.N + ctx.NW
    def slice_array(arr):
        return [arr[i] for i in range(1, total + 1)]

    return {
        "ctx": {
            "N": ctx.N,
            "NW": ctx.NW,
            "WAKLEN": ctx.WAKLEN,
            "ALFA": ctx.ALFA,
            "QINF": ctx.QINF,
            "SHARP": ctx.SHARP,
            "ANTE": ctx.ANTE,
            "ASTE": ctx.ASTE,
            "DSTE": ctx.DSTE,
            "XTE": ctx.XTE,
            "YTE": ctx.YTE,
            "X": slice_array(ctx.X),
            "Y": slice_array(ctx.Y),
            "XP": slice_array(ctx.XP),
            "YP": slice_array(ctx.YP),
            "S": slice_array(ctx.S),
            "NX": slice_array(ctx.NX),
            "NY": slice_array(ctx.NY),
            "APANEL": slice_array(ctx.APANEL),
        }
    }


def metrics_matrix(mat, n1, n2, samples):
    total = 0.0
    sumsq = 0.0
    maxabs = 0.0
    for i in range(1, n1 + 1):
        for j in range(1, n2 + 1):
            val = mat[i][j]
            total += val
            sumsq += val * val
            maxabs = max(maxabs, abs(val))
    sample_vals = []
    for ri, ci in samples:
        sample_vals.append(mat[ri][ci])
    return {"sum": total, "sumsq": sumsq, "maxabs": maxabs, "samples": sample_vals}


class TestGgcalcQdcalcParity(unittest.TestCase):
    def test_ggcalc_qdcalc(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

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
        cases = [
            {"ides": 12, "minf": 0.0, "waklen": 1.0},
            {"ides": 2412, "minf": 0.1, "waklen": 1.0},
            {"ides": 23012, "minf": 0.2, "waklen": 2.0},
        ]
        nside = 123
        script = pathlib.Path(__file__).with_name("compare_ggcalc_qdcalc.mjs")

        tol = 1.0e-7

        for case in cases:
            with self.subTest(case=case):
                ctx = build_context(case["ides"], params, 0.0, case["minf"], case["waklen"], nside, False)
                payload = build_payload(ctx)
                ggcalc(ctx)
                qdcalc(ctx)

                proc = subprocess.run(
                    ["node", str(script)],
                    input=json.dumps(payload),
                    text=True,
                    capture_output=True,
                    check=True,
                )
                js_results = json.loads(proc.stdout)["results"]

                n = ctx.N
                total = ctx.N + ctx.NW
                ai_samples = [
                    (1, 1),
                    (1, n),
                    (n, 1),
                    (n, n),
                    (n + 1, 1),
                    (1, n + 1),
                    (n + 1, n + 1),
                    (max(1, (n + 1) // 2), max(1, (n + 1) // 2)),
                ]
                di_samples = [
                    (1, 1),
                    (1, total),
                    (total, 1),
                    (total, total),
                    (max(1, (total + 1) // 2), max(1, (total + 1) // 2)),
                ]

                ai_metrics = metrics_matrix(ctx.AIJ, n + 1, n + 1, ai_samples)
                di_metrics = metrics_matrix(ctx.DIJ, total, total, di_samples)

                self.assertLessEqual(abs(ai_metrics["sum"] - js_results["ai"]["sum"]), tol)
                self.assertLessEqual(abs(ai_metrics["sumsq"] - js_results["ai"]["sumsq"]), tol)
                self.assertLessEqual(abs(ai_metrics["maxabs"] - js_results["ai"]["maxabs"]), tol)
                for py_val, js_val in zip(ai_metrics["samples"], js_results["ai"]["samples"]):
                    self.assertLessEqual(abs(py_val - js_val), tol)

                self.assertLessEqual(abs(di_metrics["sum"] - js_results["di"]["sum"]), tol)
                self.assertLessEqual(abs(di_metrics["sumsq"] - js_results["di"]["sumsq"]), tol)
                self.assertLessEqual(abs(di_metrics["maxabs"] - js_results["di"]["maxabs"]), tol)
                for py_val, js_val in zip(di_metrics["samples"], js_results["di"]["samples"]):
                    self.assertLessEqual(abs(py_val - js_val), tol)


if __name__ == "__main__":
    unittest.main()
