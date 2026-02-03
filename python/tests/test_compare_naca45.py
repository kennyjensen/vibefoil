import json
import pathlib
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def max_abs_diff(a, b):
    if len(a) != len(b):
        return float("inf")
    return max(abs(ai - bi) for ai, bi in zip(a, b))


class TestNaca45Parity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("node"):
            raise unittest.SkipTest("node is required for JS/Fortran parity checks")
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for NACA4/5 reference")

        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "naca45_driver.f90"
        srcs = [
            driver_src,
            ROOT / "third_party" / "Xfoil" / "src" / "naca.f",
        ]
        driver_path = pathlib.Path(cls._tmpdir.name) / "naca45_driver"
        subprocess.run(
            ["gfortran", "-O2", *map(str, srcs), "-o", str(driver_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        cls.driver_path = driver_path

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, "_tmpdir"):
            cls._tmpdir.cleanup()

    def run_fortran(self, ides, nside):
        proc = subprocess.run(
            [str(self.driver_path)],
            input=f"{ides} {nside}\n",
            text=True,
            capture_output=True,
            check=True,
        )
        lines = proc.stdout.strip().splitlines()
        if not lines:
            self.fail("Empty output from Fortran driver")
        count = int(lines[0].strip())
        x = []
        y = []
        for line in lines[1:]:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            x.append(float(parts[0]))
            y.append(float(parts[1]))
        self.assertEqual(len(x), count)
        return count, x, y

    def test_naca4_naca5_against_xfoil(self):
        cases = [
            {"ides": 12, "nside": 123},
            {"ides": 2412, "nside": 123},
            {"ides": 23012, "nside": 123},
            {"ides": 25012, "nside": 123},
        ]

        payload = {"cases": cases}
        script = pathlib.Path(__file__).with_name("compare_naca45.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = {item["ides"]: item for item in json.loads(proc.stdout)["results"]}

        tol = 5.0e-6

        for case in cases:
            ides = case["ides"]
            js_case = js_results[ides]
            count, x_ref, y_ref = self.run_fortran(ides, case["nside"])

            self.assertEqual(count, js_case["nb"])
            self.assertLessEqual(max_abs_diff(x_ref, js_case["xb"]), tol)
            self.assertLessEqual(max_abs_diff(y_ref, js_case["yb"]), tol)


if __name__ == "__main__":
    unittest.main()
