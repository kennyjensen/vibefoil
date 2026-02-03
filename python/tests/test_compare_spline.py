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



class TestSplineParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for Fortran parity checks")

        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "spline_driver.f"
        srcs = [
            driver_src,
            ROOT / "third_party" / "Xfoil" / "src" / "spline.f",
        ]
        driver_path = pathlib.Path(cls._tmpdir.name) / "spline_driver"
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

    def test_spline_functions(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        n = 30
        x = [i / (n - 1) for i in range(n)]
        y = [0.1 * math.sin(2.0 * math.pi * xi) for xi in x]

        s = [0.0 for _ in range(n)]
        for i in range(1, n):
            dx = x[i] - x[i - 1]
            dy = y[i] - y[i - 1]
            s[i] = s[i - 1] + math.sqrt(dx * dx + dy * dy)
        s_total = s[-1]

        queries = [s[0], s[n // 2], 0.73 * s_total]
        sinvrt_cases = [
            {"si": s[n // 5], "xi": x[n // 5]},
            {"si": s[n // 2], "xi": x[n // 2]},
            {"si": s[(4 * n) // 5], "xi": x[(4 * n) // 5]},
        ]

        payload = {
            "x": x,
            "y": y,
            "queries": queries,
            "sinvrtCases": sinvrt_cases,
        }

        script = pathlib.Path(__file__).with_name("compare_spline.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        input_lines = [f"{n} {len(queries)} {len(sinvrt_cases)}"]
        input_lines.extend(str(val) for val in x)
        input_lines.extend(str(val) for val in y)
        input_lines.extend(str(val) for val in queries)
        input_lines.extend(f"{case['si']} {case['xi']}" for case in sinvrt_cases)

        proc_f = subprocess.run(
            [str(self.driver_path)],
            input="\n".join(input_lines) + "\n",
            text=True,
            capture_output=True,
            check=True,
        )
        lines = [line.strip() for line in proc_f.stdout.splitlines() if line.strip()]

        out_scalc = []
        out_segspl = []
        out_evals = []
        out_sinvrt = []
        for line in lines:
            parts = line.split()
            tag = parts[0].lower()
            vals = [float(v) for v in parts[1:]]
            if tag == "scalc":
                out_scalc = vals
            elif tag == "segspl":
                out_segspl = vals
            elif tag == "eval":
                out_evals.append(vals)
            elif tag == "sinvrt":
                out_sinvrt.append(vals[0])

        tol = 1.0e-5
        for a, b in zip(out_scalc, js_results["scalc"]):
            self.assertLessEqual(abs(a - b), tol)

        for a, b in zip(out_segspl, js_results["segspl"]):
            self.assertLessEqual(abs(a - b), tol)

        for idx, evals in enumerate(js_results["evals"]):
            self.assertLessEqual(abs(out_evals[idx][0] - evals["seval"]), tol)
            self.assertLessEqual(abs(out_evals[idx][1] - evals["deval"]), tol)
            self.assertLessEqual(abs(out_evals[idx][2] - evals["d2val"]), tol)

        for idx, val in enumerate(js_results["sinvrt"]):
            self.assertLessEqual(abs(out_sinvrt[idx] - val), tol)


if __name__ == "__main__":
    unittest.main()
