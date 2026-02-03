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
from python.spline import scalc, segspl, sinvrt, seval


class TestXgdesParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for Fortran parity checks")

        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "xgdes_driver.f"
        srcs = [
            driver_src,
            ROOT / "python" / "tests" / "xgdes_subs.f",
            ROOT / "python" / "tests" / "xgeom_subs.f",
            ROOT / "third_party" / "Xfoil" / "src" / "spline.f",
            ROOT / "third_party" / "Xfoil" / "src" / "userio.f",
        ]
        driver_path = pathlib.Path(cls._tmpdir.name) / "xgdes_driver"
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

    def test_xgdes_helpers(self):
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

        xf = 0.7
        yrel = 0.4
        tops = ctx.S[1] + (ctx.X[1] - xf)
        bots = ctx.S[ctx.N] - (ctx.X[ctx.N] - xf)
        tops = sinvrt(tops, xf, ctx.X, ctx.XP, ctx.S, ctx.N)
        bots = sinvrt(bots, xf, ctx.X, ctx.XP, ctx.S, ctx.N)
        topy = seval(tops, ctx.Y, ctx.YP, ctx.S, ctx.N)
        boty = seval(bots, ctx.Y, ctx.YP, ctx.S, ctx.N)
        ybf = topy * yrel + boty * (1.0 - yrel)

        inside_pts = [
            {"x": xf, "y": ybf},
            {"x": 1.5, "y": 0.0},
        ]

        sss_cases = [
            {"ss": tops, "del": 0.1, "xbf": xf, "ybf": ybf, "iside": 1},
            {"ss": bots, "del": 0.1, "xbf": xf, "ybf": ybf, "iside": 2},
        ]

        payload = {
            "x": [ctx.X[i] for i in range(1, ctx.N + 1)],
            "y": [ctx.Y[i] for i in range(1, ctx.N + 1)],
            "s": [ctx.S[i] for i in range(1, ctx.N + 1)],
            "xf": xf,
            "yrel": yrel,
            "insidePts": inside_pts,
            "sssCases": sss_cases,
        }

        script = pathlib.Path(__file__).with_name("compare_xgdes.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        input_lines = [
            f"{ctx.N} {len(inside_pts)} {len(sss_cases)}",
        ]
        input_lines.extend(str(ctx.X[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.Y[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.S[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.XP[i]) for i in range(1, ctx.N + 1))
        input_lines.extend(str(ctx.YP[i]) for i in range(1, ctx.N + 1))
        input_lines.append(f"{xf} 999.0")
        input_lines.append(str(yrel))
        input_lines.extend(f"{pt['x']} {pt['y']}" for pt in inside_pts)
        input_lines.extend(f"{case['ss']} {case['del']} {case['xbf']} {case['ybf']} {case['iside']}" for case in sss_cases)

        proc_f = subprocess.run(
            [str(self.driver_path)],
            input="\n".join(input_lines) + "\n",
            text=True,
            capture_output=True,
            check=True,
        )
        lines = [line.strip() for line in proc_f.stdout.splitlines() if line.strip()]

        out_getxyf = None
        out_inside = []
        out_sss = []
        for line in lines:
            if "GETXYF" in line:
                tail = line.split("GETXYF", 1)[1].strip()
                out_getxyf = [float(v) for v in tail.split()]
                continue
            parts = line.split()
            if not parts:
                continue
            tag = parts[0].lower()
            if tag == "inside":
                out_inside.append(int(parts[1]))
            elif tag == "sss":
                out_sss.append([float(parts[1]), float(parts[2])])

        tol = 1.0e-6
        self.assertIsNotNone(out_getxyf)
        self.assertLessEqual(abs(out_getxyf[0] - js_results["getxyf"]["tops"]), tol)
        self.assertLessEqual(abs(out_getxyf[1] - js_results["getxyf"]["bots"]), tol)
        self.assertLessEqual(abs(out_getxyf[2] - js_results["getxyf"]["xf"]), tol)
        self.assertLessEqual(abs(out_getxyf[3] - js_results["getxyf"]["yf"]), tol)

        for idx, val in enumerate(js_results["inside"]):
            self.assertEqual(out_inside[idx], val)

        for idx, val in enumerate(js_results["sss"]):
            self.assertLessEqual(abs(out_sss[idx][0] - val["s1"]), tol)
            self.assertLessEqual(abs(out_sss[idx][1] - val["s2"]), tol)


if __name__ == "__main__":
    unittest.main()
