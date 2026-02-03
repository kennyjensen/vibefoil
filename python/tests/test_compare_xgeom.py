import json
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
from python.spline import scalc, segspl


class TestXgeomParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for Fortran parity checks")

        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "xgeom_driver.f"
        srcs = [
            driver_src,
            ROOT / "python" / "tests" / "xgeom_subs.f",
            ROOT / "third_party" / "Xfoil" / "src" / "spline.f",
        ]
        driver_path = pathlib.Path(cls._tmpdir.name) / "xgeom_driver"
        subprocess.run(
            ["gfortran", "-O2", "-ffixed-form", *map(str, srcs), "-o", str(driver_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        cls.driver_path = driver_path

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, "_tmpdir"):
            cls._tmpdir.cleanup()

    def test_lefind_parity(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        ctx = XFoilState()
        ctx.NPAN = 80
        ctx.CVPAR = 1.0
        ctx.CTERAT = 0.15
        ctx.CTRRAT = 0.2
        ctx.XSREF1 = 1.0
        ctx.XSREF2 = 1.0
        ctx.XPREF1 = 1.0
        ctx.XPREF2 = 1.0

        naca(ctx, 2412)
        pangen(ctx, False)

        scalc(ctx.X, ctx.Y, ctx.S, ctx.N)
        segspl(ctx.X, ctx.XP, ctx.S, ctx.N)
        segspl(ctx.Y, ctx.YP, ctx.S, ctx.N)

        payload = {
            "x": [ctx.X[i] for i in range(1, ctx.N + 1)],
            "y": [ctx.Y[i] for i in range(1, ctx.N + 1)],
            "s": [ctx.S[i] for i in range(1, ctx.N + 1)],
        }

        script = pathlib.Path(__file__).with_name("compare_xgeom.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        input_lines = [str(ctx.N)]
        input_lines.extend(str(ctx.X[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.Y[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.S[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.XP[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.YP[i]) for i in range(1, ctx.N + 1))

        proc_f = subprocess.run(
            [str(self.driver_path)],
            input="\n".join(input_lines) + "\n",
            text=True,
            capture_output=True,
            check=True,
        )
        sle_f = float(proc_f.stdout.strip().splitlines()[-1].strip())

        tol = 1.0e-6
        self.assertLessEqual(abs(sle_f - js_results["sle"]), tol)


if __name__ == "__main__":
    unittest.main()
