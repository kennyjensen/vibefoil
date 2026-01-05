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
from python.xfoil import comset, naca, pangen, tecalc
from python.spline import scalc, segspl
from python.xpanel import apcalc, ncalc
from python.xqdes import SPLQSP, SMOOQ, SYMQSP, GAMQSP, QINCOM, MIXED


def build_spline_data(nsp):
    sspec = []
    qspec = []
    qgamm = []
    for i in range(1, nsp + 1):
        s = float(i - 1) / float(nsp - 1)
        sspec.append(s)
        qspec.append(math.sin(2.0 * math.pi * s) + 0.3 * s)
        qgamm.append(0.5 * math.cos(2.0 * math.pi * s) - 0.1 * s)
    return sspec, qspec, qgamm


def assert_array_close(testcase, arr_py, arr_js, tol, label):
    testcase.assertEqual(len(arr_py), len(arr_js), f"length mismatch for {label}")
    for idx, (py_val, js_val) in enumerate(zip(arr_py, arr_js)):
        testcase.assertLessEqual(abs(py_val - js_val), tol, f"{label}[{idx}] mismatch")


def metrics_array(values, samples):
    total = 0.0
    sumsq = 0.0
    maxabs = 0.0
    for val in values:
        total += val
        sumsq += val * val
        maxabs = max(maxabs, abs(val))
    sample_vals = [values[i - 1] for i in samples]
    return {"sum": total, "sumsq": sumsq, "maxabs": maxabs, "samples": sample_vals}


class TestXqdesParity(unittest.TestCase):
    def test_qdes_spline_parity(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        nsp = 12
        kqsp = 1
        kq1 = 3
        kq2 = nsp - 2
        lqslop = True
        algam = 0.12
        clgam = 0.65
        cmgam = -0.02
        liqset = False

        sspec, qspec, qgamm = build_spline_data(nsp)

        ctx_spl = XFoilState()
        ctx_spl.NSP = nsp
        ctx_spl.LQSLOP = lqslop
        for i in range(1, nsp + 1):
            ctx_spl.SSPEC[i] = sspec[i - 1]
            ctx_spl.QSPEC[i][kqsp] = qspec[i - 1]
        SPLQSP(ctx_spl, kqsp)
        py_spl = [ctx_spl.QSPECP[i][kqsp] for i in range(1, nsp + 1)]

        ctx_smo = XFoilState()
        ctx_smo.NSP = nsp
        ctx_smo.LQSLOP = lqslop
        for i in range(1, nsp + 1):
            ctx_smo.SSPEC[i] = sspec[i - 1]
            ctx_smo.QSPEC[i][kqsp] = qspec[i - 1]
        SMOOQ(ctx_smo, kq1, kq2, kqsp)
        py_smo = [ctx_smo.QSPEC[i][kqsp] for i in range(1, nsp + 1)]

        ctx_sym = XFoilState()
        ctx_sym.NSP = nsp
        for i in range(1, nsp + 1):
            ctx_sym.SSPEC[i] = sspec[i - 1]
            ctx_sym.QSPEC[i][kqsp] = qspec[i - 1]
        SYMQSP(ctx_sym, kqsp)
        py_sym_s = [ctx_sym.SSPEC[i] for i in range(1, nsp + 1)]
        py_sym_q = [ctx_sym.QSPEC[i][kqsp] for i in range(1, nsp + 1)]

        ctx_gam = XFoilState()
        ctx_gam.NSP = nsp
        ctx_gam.LIQSET = liqset
        ctx_gam.ALGAM = algam
        ctx_gam.CLGAM = clgam
        ctx_gam.CMGAM = cmgam
        for i in range(1, nsp + 1):
            ctx_gam.SSPEC[i] = sspec[i - 1]
            ctx_gam.QGAMM[i] = qgamm[i - 1]
        GAMQSP(ctx_gam, kqsp)
        py_gam_qspec = [ctx_gam.QSPEC[i][kqsp] for i in range(1, nsp + 1)]

        qincom_cases = [
            {"qc": 0.05, "qinf": 1.0, "tklam": 0.2},
            {"qc": -0.2, "qinf": 0.9, "tklam": 0.05},
            {"qc": 0.15, "qinf": 1.1, "tklam": 0.4},
        ]
        py_qincom = [QINCOM(case["qc"], case["qinf"], case["tklam"]) for case in qincom_cases]

        payload = {
            "nsp": nsp,
            "sspec": sspec,
            "qspec": qspec,
            "qgamm": qgamm,
            "kq1": kq1 - 1,
            "kq2": kq2 - 1,
            "lqslop": lqslop,
            "algam": algam,
            "clgam": clgam,
            "cmgam": cmgam,
            "liqset": liqset,
            "qincomCases": qincom_cases,
        }

        script = pathlib.Path(__file__).with_name("compare_xqdes_spline.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        results = json.loads(proc.stdout)["results"]

        tol = 1.0e-7
        assert_array_close(self, py_spl, results["splqsp"]["qspecp"], tol, "splqsp.qspecp")
        assert_array_close(self, py_smo, results["smooq"]["qspec"], tol, "smooq.qspec")
        assert_array_close(self, py_sym_s, results["symqsp"]["sspec"], tol, "symqsp.sspec")
        assert_array_close(self, py_sym_q, results["symqsp"]["qspec"], tol, "symqsp.qspec")

        gam_results = results["gamqsp"]
        self.assertLessEqual(abs(ctx_gam.ALQSP[kqsp] - gam_results["alqsp"]), tol)
        self.assertLessEqual(abs(ctx_gam.CLQSP[kqsp] - gam_results["clqsp"]), tol)
        self.assertLessEqual(abs(ctx_gam.CMQSP[kqsp] - gam_results["cmqsp"]), tol)
        assert_array_close(self, py_gam_qspec, gam_results["qspec"], tol, "gamqsp.qspec")
        self.assertLessEqual(abs(ctx_gam.QDOF0 - gam_results["qdof0"]), tol)
        self.assertLessEqual(abs(ctx_gam.QDOF1 - gam_results["qdof1"]), tol)
        self.assertLessEqual(abs(ctx_gam.QDOF2 - gam_results["qdof2"]), tol)
        self.assertLessEqual(abs(ctx_gam.QDOF3 - gam_results["qdof3"]), tol)
        self.assertEqual(ctx_gam.IQ1, gam_results["iq1"] + 1)
        self.assertEqual(ctx_gam.IQ2, gam_results["iq2"] + 1)

        assert_array_close(self, py_qincom, results["qincom"], tol, "qincom")

    def test_qdes_mixed_parity(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        ctx = XFoilState()
        ctx.NPAN = 64
        ctx.CVPAR = 1.0
        ctx.CTERAT = 0.15
        ctx.CTRRAT = 0.2
        ctx.XSREF1 = 1.0
        ctx.XSREF2 = 1.0
        ctx.XPREF1 = 1.0
        ctx.XPREF2 = 1.0

        alpha_deg = 2.0
        ctx.ALFA = alpha_deg * ctx.DTOR
        ctx.ADEG = alpha_deg
        ctx.MINF = 0.1
        ctx.MINF1 = ctx.MINF
        ctx.QINF = 1.0
        ctx.XCMREF = 0.0
        ctx.YCMREF = 0.0

        naca(ctx, 2412)
        pangen(ctx, False)
        comset(ctx)

        scalc(ctx.X, ctx.Y, ctx.S, ctx.N)
        segspl(ctx.X, ctx.XP, ctx.S, ctx.N)
        segspl(ctx.Y, ctx.YP, ctx.S, ctx.N)
        ncalc(ctx.X, ctx.Y, ctx.S, ctx.N, ctx.NX, ctx.NY)
        apcalc(ctx)
        tecalc(ctx)

        ctx.NSP = ctx.N
        ctx.IQ1 = 2
        ctx.IQ2 = ctx.N - 1
        ctx.LCPXX = True
        ctx.PSIO = 0.0
        ctx.QDOF0 = 0.01
        ctx.QDOF1 = -0.015
        ctx.QDOF2 = 0.004
        ctx.QDOF3 = -0.003
        ctx.LIMAGE = False
        ctx.YIMAGE = 0.0

        for i in range(1, ctx.N + 1):
            ctx.SSPEC[i] = ctx.S[i]
            sfrac = (ctx.S[i] - ctx.S[1]) / (ctx.S[ctx.N] - ctx.S[1])
            ctx.QSPEC[i][1] = 1.0 + 0.1 * math.sin(2.0 * math.pi * sfrac)
            ctx.GAM[i] = ctx.QSPEC[i][1]
            ctx.SIG[i] = 0.0
            ctx.GAM_A[i] = 0.0

        samples = [1, (ctx.N + 1) // 2, ctx.N]
        payload = {
            "n": ctx.N,
            "nsp": ctx.NSP,
            "x": [ctx.X[i] for i in range(1, ctx.N + 1)],
            "y": [ctx.Y[i] for i in range(1, ctx.N + 1)],
            "s": [ctx.S[i] for i in range(1, ctx.N + 1)],
            "xp": [ctx.XP[i] for i in range(1, ctx.N + 1)],
            "yp": [ctx.YP[i] for i in range(1, ctx.N + 1)],
            "nx": [ctx.NX[i] for i in range(1, ctx.N + 1)],
            "ny": [ctx.NY[i] for i in range(1, ctx.N + 1)],
            "apanel": [ctx.APANEL[i] for i in range(1, ctx.N + 1)],
            "qspec": [ctx.QSPEC[i][1] for i in range(1, ctx.N + 1)],
            "sspec": [ctx.SSPEC[i] for i in range(1, ctx.N + 1)],
            "gam": [ctx.GAM[i] for i in range(1, ctx.N + 1)],
            "sig": [ctx.SIG[i] for i in range(1, ctx.N + 1)],
            "alfa": ctx.ALFA,
            "minf": ctx.MINF,
            "qinf": ctx.QINF,
            "xcmref": ctx.XCMREF,
            "ycmref": ctx.YCMREF,
            "psio": ctx.PSIO,
            "qdof0": ctx.QDOF0,
            "qdof1": ctx.QDOF1,
            "qdof2": ctx.QDOF2,
            "qdof3": ctx.QDOF3,
            "iq1": ctx.IQ1 - 1,
            "iq2": ctx.IQ2 - 1,
            "lcpXX": ctx.LCPXX,
            "limage": ctx.LIMAGE,
            "yimage": ctx.YIMAGE,
            "sharp": ctx.SHARP,
            "ante": ctx.ANTE,
            "aste": ctx.ASTE,
            "dste": ctx.DSTE,
            "xte": ctx.XTE,
            "yte": ctx.YTE,
            "niterq": 2,
            "samples": samples,
        }

        script = pathlib.Path(__file__).with_name("compare_xqdes_mixed.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )

        MIXED(ctx, 1, 2)

        py_metrics = {
            "x": metrics_array([ctx.X[i] for i in range(1, ctx.N + 1)], samples),
            "y": metrics_array([ctx.Y[i] for i in range(1, ctx.N + 1)], samples),
            "gam": metrics_array([ctx.GAM[i] for i in range(1, ctx.N + 1)], samples),
        }

        results = json.loads(proc.stdout)["results"]
        tol = 1.0e-6

        for key in ["x", "y", "gam"]:
            self.assertLessEqual(abs(py_metrics[key]["sum"] - results["metrics"][key]["sum"]), tol)
            self.assertLessEqual(abs(py_metrics[key]["sumsq"] - results["metrics"][key]["sumsq"]), tol)
            self.assertLessEqual(abs(py_metrics[key]["maxabs"] - results["metrics"][key]["maxabs"]), tol)
            assert_array_close(self, py_metrics[key]["samples"], results["metrics"][key]["samples"], tol, f"mixed.{key}.samples")

        self.assertLessEqual(abs(ctx.PSIO - results["psio"]), tol)
        self.assertLessEqual(abs(ctx.QDOF0 - results["qdof0"]), tol)
        self.assertLessEqual(abs(ctx.QDOF1 - results["qdof1"]), tol)
        self.assertLessEqual(abs(ctx.QDOF2 - results["qdof2"]), tol)
        self.assertLessEqual(abs(ctx.QDOF3 - results["qdof3"]), tol)
        self.assertLessEqual(abs(ctx.CL - results["cl"]), tol)
        self.assertLessEqual(abs(ctx.CM - results["cm"]), tol)
        self.assertLessEqual(abs(ctx.CDP - results["cdp"]), tol)
        self.assertLessEqual(abs(ctx.CL_ALF - results["clAlf"]), tol)
        self.assertLessEqual(abs(ctx.CL_MSQ - results["clMsq"]), tol)


if __name__ == "__main__":
    unittest.main()
