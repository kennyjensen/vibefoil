import json
import math
import pathlib
import shutil
import subprocess
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]

import sys

sys.path.insert(0, str(ROOT))

from python.xbl import XFoilState
from python.xfoil import naca, pangen


class TestXfoilParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for Fortran parity checks")

        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "xfoil_driver.f"
        srcs = [
            driver_src,
            ROOT / "python" / "tests" / "xfoil_subs.f",
            ROOT / "third_party" / "Xfoil" / "src" / "spline.f",
        ]
        driver_path = pathlib.Path(cls._tmpdir.name) / "xfoil_driver"
        subprocess.run(
            ["gfortran", "-O2", "-ffixed-form", "-I", str(ROOT / "third_party" / "Xfoil" / "src"), *map(str, srcs), "-o", str(driver_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        cls.driver_path = driver_path

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, "_tmpdir"):
            cls._tmpdir.cleanup()

    def test_xfoil_functions(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        ctx = XFoilState()
        ctx.NPAN = 60
        ctx.CVPAR = 1.0
        ctx.CTERAT = 0.15
        ctx.CTRRAT = 0.2
        ctx.XSREF1 = 1.0
        ctx.XSREF2 = 1.0
        ctx.XPREF1 = 1.0
        ctx.XPREF2 = 1.0

        naca(ctx, 12)
        pangen(ctx, False)

        n = ctx.N
        x = [ctx.X[i] for i in range(1, n + 1)]
        y = [ctx.Y[i] for i in range(1, n + 1)]

        gam = []
        gam_a = []
        q = [0.0]
        for i in range(n):
            sfrac = i / (n - 1)
            g = 1.0 + 0.1 * math.sin(2.0 * math.pi * sfrac)
            gam.append(g)
            gam_a.append(0.0)
            q.append(1.0 + 0.05 * math.cos(2.0 * math.pi * sfrac))

        alfa = 2.0 * math.pi / 180.0
        minf = 0.1
        qinf = 1.0
        xref = 0.25
        yref = 0.0

        payload = {
            "x": x,
            "y": y,
            "gam": gam,
            "gamA": gam_a,
            "q": q,
            "alfa": alfa,
            "minf": minf,
            "qinf": qinf,
            "xref": xref,
            "yref": yref,
        }

        script = pathlib.Path(__file__).with_name("compare_xfoil.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        input_lines = [str(n)]
        input_lines.append(f"{alfa} {minf} {qinf} {xref} {yref}")
        input_lines.extend(str(v) for v in x)
        input_lines.extend(str(v) for v in y)
        input_lines.extend(str(v) for v in gam)
        input_lines.extend(str(v) for v in gam_a)
        input_lines.extend(str(v) for v in q[1:])

        proc_f = subprocess.run(
            [str(self.driver_path)],
            input="\n".join(input_lines) + "\n",
            text=True,
            capture_output=True,
            check=True,
        )
        lines = [line.strip() for line in proc_f.stdout.splitlines() if line.strip()]

        out_te = None
        out_cp = None
        out_cl = None
        for line in lines:
            if line.startswith("TE"):
                out_te = [float(v) for v in line.split()[1:]]
            elif line.startswith("CP"):
                out_cp = [float(v) for v in line.split()[1:]]
            elif line.startswith("CL"):
                out_cl = [float(v) for v in line.split()[1:]]

        tol = 1.0e-6
        self.assertIsNotNone(out_te)
        self.assertIsNotNone(out_cp)
        self.assertIsNotNone(out_cl)

        self.assertLessEqual(abs(out_te[0] - js_results["te"]["ante"]), tol)
        self.assertLessEqual(abs(out_te[1] - js_results["te"]["aste"]), tol)
        self.assertLessEqual(abs(out_te[2] - js_results["te"]["dste"]), tol)
        self.assertLessEqual(abs(out_te[3] - js_results["te"]["sharp"]), tol)

        for a, b in zip(out_cp, js_results["cp"][1:]):
            self.assertLessEqual(abs(a - b), tol)

        self.assertLessEqual(abs(out_cl[0] - js_results["cl"]["cl"]), tol)
        self.assertLessEqual(abs(out_cl[1] - js_results["cl"]["cm"]), tol)
        self.assertLessEqual(abs(out_cl[2] - js_results["cl"]["cdp"]), tol)
        self.assertLessEqual(abs(out_cl[3] - js_results["cl"]["clAlf"]), tol)
        self.assertLessEqual(abs(out_cl[4] - js_results["cl"]["clMsq"]), tol)


if __name__ == "__main__":
    unittest.main()
