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


class TestNaca6Parity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("node"):
            raise unittest.SkipTest("node is required for JS/Fortran parity checks")
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for NACA456 reference")

        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "naca6_driver.f90"
        srcs = [
            driver_src,
            ROOT / "third_party" / "naca456" / "splprocs.f90",
            ROOT / "third_party" / "naca456" / "epspsi.f90",
            ROOT / "third_party" / "naca456" / "nacax.f90",
        ]
        driver_path = pathlib.Path(cls._tmpdir.name) / "naca6_driver"
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

    def run_fortran(self, family, a, cl, toc, nside):
        proc = subprocess.run(
            [str(self.driver_path)],
            input=f"{family} {a} {cl} {toc} {nside}\n",
            text=True,
            capture_output=True,
            check=True,
        )
        lines = proc.stdout.strip().splitlines()
        if not lines:
            self.fail("Empty output from Fortran driver")
        count = int(lines[0].strip())
        self.assertEqual(count, nside)
        x = []
        yu = []
        yl = []
        for line in lines[1:]:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            x.append(float(parts[0]))
            yu.append(float(parts[1]))
            yl.append(float(parts[2]))
        self.assertEqual(len(x), nside)
        return x, yu, yl

    def test_naca6_against_naca456(self):
        cases = [
            {"profile": "63", "family": 1, "toc": 0.12, "camber": "6", "cl": 0.2, "a": 1.0, "nside": 121},
            {"profile": "65", "family": 3, "toc": 0.10, "camber": "6", "cl": 0.3, "a": 0.8, "nside": 121},
            {"profile": "63A", "family": 6, "toc": 0.12, "camber": "6", "cl": 0.25, "a": 0.5, "nside": 121},
            {"profile": "64A", "family": 7, "toc": 0.15, "camber": "0", "cl": 0.0, "a": 0.8, "nside": 121},
        ]

        payload = {"cases": cases}
        script = pathlib.Path(__file__).with_name("compare_naca6.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = {item["profile"]: item for item in json.loads(proc.stdout)["results"]}

        tol = 5.0e-6

        for case in cases:
            js_case = js_results[case["profile"]]
            x_ref, yu_ref, yl_ref = self.run_fortran(
                case["family"],
                case["a"],
                case["cl"],
                case["toc"],
                case["nside"],
            )

            self.assertLessEqual(max_abs_diff(x_ref, js_case["x"]), tol)
            self.assertLessEqual(max_abs_diff(yu_ref, js_case["yu"]), tol)
            self.assertLessEqual(max_abs_diff(yl_ref, js_case["yl"]), tol)


if __name__ == "__main__":
    unittest.main()
