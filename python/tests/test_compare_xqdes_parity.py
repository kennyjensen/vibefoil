import json
import math
import pathlib
import shutil
import subprocess
import sys
import unittest
import tempfile

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
    @classmethod
    def setUpClass(cls):
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for Fortran parity checks")

        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "xqdes_driver.f"
        stub_src = ROOT / "python" / "tests" / "xqdes_stubs.f"
        srcs = [
            driver_src,
            stub_src,
            ROOT / "python" / "tests" / "xfoil_subs.f",
            ROOT / "python" / "tests" / "xgeom_subs.f",
            ROOT / "third_party" / "Xfoil" / "src" / "xqdes.f",
            ROOT / "third_party" / "Xfoil" / "src" / "xpanel.f",
            ROOT / "third_party" / "Xfoil" / "src" / "xsolve.f",
            ROOT / "third_party" / "Xfoil" / "src" / "spline.f",
            ROOT / "third_party" / "Xfoil" / "src" / "xutils.f",
            ROOT / "third_party" / "Xfoil" / "src" / "userio.f",
        ]
        driver_path = pathlib.Path(cls._tmpdir.name) / "xqdes_driver"
        subprocess.run(
            [
                "gfortran",
                "-O2",
                "-ffixed-form",
                "-I",
                str(ROOT / "third_party" / "Xfoil" / "src"),
                *map(str, srcs),
                "-o",
                str(driver_path),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        cls.driver_path = driver_path

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, "_tmpdir"):
            cls._tmpdir.cleanup()
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

        # Run Fortran reference for MIXED
        input_lines = []
        input_lines.append(f"{ctx.N} {ctx.NSP} {ctx.IQ1} {ctx.IQ2} {2}")
        input_lines.append(
            f"{ctx.ALFA} {ctx.MINF} {ctx.QINF} {ctx.XCMREF} {ctx.YCMREF} {ctx.PSIO} {ctx.QDOF0} {ctx.QDOF1} {ctx.QDOF2} {ctx.QDOF3}"
        )
        input_lines.append(f"{1 if ctx.LCPXX else 0} {1 if ctx.LIMAGE else 0} {1 if ctx.SHARP else 0}")
        input_lines.append(f"{ctx.ANTE} {ctx.ASTE} {ctx.DSTE} {ctx.XTE} {ctx.YTE} {ctx.YIMAGE}")
        input_lines.extend(str(ctx.X[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.Y[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.S[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.XP[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.YP[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.NX[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.NY[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.APANEL[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.SSPEC[i]) for i in range(1, ctx.NSP + 1))
        input_lines.extend(str(ctx.QSPEC[i][1]) for i in range(1, ctx.NSP + 1))
        input_lines.extend(str(ctx.GAM[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.SIG[i]) for i in range(1, ctx.N + 1))

        proc_f = subprocess.run(
            [str(self.driver_path)],
            input="\n".join(input_lines) + "\n",
            text=True,
            capture_output=True,
        )
        if proc_f.returncode != 0:
            self.fail(
                f"Fortran driver failed (code {proc_f.returncode}).\\n"
                f"stdout:\\n{proc_f.stdout}\\n"
                f"stderr:\\n{proc_f.stderr}"
            )
        lines = proc_f.stdout.strip().splitlines()
        if not lines:
            self.fail("Empty output from Fortran driver")
        start_idx = None
        for i, line in enumerate(lines):
            parts = line.strip().split()
            if len(parts) == 1:
                try:
                    n_out = int(parts[0])
                except ValueError:
                    continue
                start_idx = i
                break
        if start_idx is None:
            self.fail("Fortran output missing panel count line")
        self.assertEqual(n_out, ctx.N)
        x_f = []
        y_f = []
        gam_f = []
        idx = start_idx + 1
        for _ in range(ctx.N):
            parts = lines[idx].strip().split()
            if len(parts) < 3:
                self.fail("Fortran output row missing fields")
            x_f.append(float(parts[0]))
            y_f.append(float(parts[1]))
            gam_f.append(float(parts[2]))
            idx += 1
        coeff_parts = lines[idx].strip().split()
        if len(coeff_parts) < 10:
            self.fail("Fortran output missing coefficients")
        psio_f = float(coeff_parts[0])
        qdof0_f = float(coeff_parts[1])
        qdof1_f = float(coeff_parts[2])
        qdof2_f = float(coeff_parts[3])
        qdof3_f = float(coeff_parts[4])
        cl_f = float(coeff_parts[5])
        cm_f = float(coeff_parts[6])
        cdp_f = float(coeff_parts[7])
        clalf_f = float(coeff_parts[8])
        clmsq_f = float(coeff_parts[9])

        results = json.loads(proc.stdout)["results"]
        # Fortran is single-precision; allow small FP drift vs JS double.
        tol = 2.0e-2

        f_metrics = {
            "x": metrics_array(x_f, samples),
            "y": metrics_array(y_f, samples),
            "gam": metrics_array(gam_f, samples),
        }

        for key in ["x", "y", "gam"]:
            self.assertLessEqual(abs(f_metrics[key]["maxabs"] - results["metrics"][key]["maxabs"]), tol)
            assert_array_close(self, f_metrics[key]["samples"], results["metrics"][key]["samples"], tol, f"mixed.{key}.samples")

        self.assertLessEqual(abs(psio_f - results["psio"]), tol)
        self.assertLessEqual(abs(qdof0_f - results["qdof0"]), tol)
        self.assertLessEqual(abs(qdof1_f - results["qdof1"]), tol)
        self.assertLessEqual(abs(qdof2_f - results["qdof2"]), tol)
        self.assertLessEqual(abs(qdof3_f - results["qdof3"]), tol)
        self.assertLessEqual(abs(cl_f - results["cl"]), tol)
        self.assertLessEqual(abs(cm_f - results["cm"]), tol)
        self.assertLessEqual(abs(cdp_f - results["cdp"]), tol)
        self.assertLessEqual(abs(clalf_f - results["clAlf"]), tol)
        self.assertLessEqual(abs(clmsq_f - results["clMsq"]), tol)


if __name__ == "__main__":
    unittest.main()
