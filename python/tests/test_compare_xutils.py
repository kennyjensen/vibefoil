import json
import pathlib
import shutil
import subprocess
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]


class TestXutilsParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for Fortran parity checks")

        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "xutils_driver.f"
        srcs = [
            driver_src,
            ROOT / "third_party" / "Xfoil" / "src" / "xutils.f",
        ]
        driver_path = pathlib.Path(cls._tmpdir.name) / "xutils_driver"
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

    def test_xutils_functions(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        payload = {
            "atanc": [
                {"y": 0.2, "x": 0.8, "thold": 0.1},
                {"y": -0.5, "x": 0.1, "thold": 1.2},
                {"y": 0.7, "x": -0.3, "thold": -2.4},
            ],
            "setexp": [
                {"ds1": 0.02, "smax": 1.0, "nn": 12},
                {"ds1": 0.01, "smax": 0.5, "nn": 20},
            ],
        }

        script = pathlib.Path(__file__).with_name("compare_xutils.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        input_lines = [
            f"{len(payload['atanc'])} {len(payload['setexp'])}",
        ]
        input_lines.extend(f"{case['y']} {case['x']} {case['thold']}" for case in payload["atanc"])
        input_lines.extend(f"{case['ds1']} {case['smax']} {case['nn']}" for case in payload["setexp"])

        proc_f = subprocess.run(
            [str(self.driver_path)],
            input="\n".join(input_lines) + "\n",
            text=True,
            capture_output=True,
            check=True,
        )
        lines = [line.strip() for line in proc_f.stdout.splitlines() if line.strip()]
        out_atanc = []
        out_setexp = []
        for line in lines:
            parts = line.split()
            tag = parts[0].lower()
            vals = [float(v) for v in parts[1:]]
            if tag == "atanc":
                out_atanc.append(vals[0])
            elif tag == "setexp":
                out_setexp.append(vals)

        tol = 1.0e-6
        for idx, val in enumerate(js_results["atanc"]):
            self.assertLessEqual(abs(out_atanc[idx] - val), tol)

        for idx, arr in enumerate(js_results["setexp"]):
            self.assertEqual(len(out_setexp[idx]), len(arr))
            for a, b in zip(out_setexp[idx], arr):
                self.assertLessEqual(abs(a - b), tol)


if __name__ == "__main__":
    unittest.main()
