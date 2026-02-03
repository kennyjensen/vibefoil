import json
import pathlib
import shutil
import subprocess
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]


class TestBluParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for Fortran parity checks")

        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "blu_driver.f"
        srcs = [
            driver_src,
            ROOT / "python" / "tests" / "blu_subs.f",
            ROOT / "python" / "tests" / "blu_cft.f",
        ]
        driver_path = pathlib.Path(cls._tmpdir.name) / "blu_driver"
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

    def test_cft_parity(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        payload = {
            "cft": [
                {"hk": 3.8, "rt": 200.0, "msq": 0.0, "cffac": 1.0},
                {"hk": 5.0, "rt": 800.0, "msq": 0.2, "cffac": 1.0},
                {"hk": 6.5, "rt": 2000.0, "msq": 0.3, "cffac": 1.0},
            ],
        }

        script = pathlib.Path(__file__).with_name("compare_blu.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        input_lines = [str(len(payload["cft"]))]
        input_lines.extend(f"{c['hk']} {c['rt']} {c['msq']}" for c in payload["cft"])

        proc_f = subprocess.run(
            [str(self.driver_path)],
            input="\n".join(input_lines) + "\n",
            text=True,
            capture_output=True,
            check=True,
        )
        lines = [line.strip() for line in proc_f.stdout.splitlines() if line.strip()]
        out = []
        for line in lines:
            parts = line.split()
            if parts and parts[0].lower() == "cft":
                out.append([float(v) for v in parts[1:]])

        tol = 1.0e-6
        for idx, js in enumerate(js_results["cft"]):
            vals = out[idx]
            self.assertLessEqual(abs(vals[0] - js["cf"]), tol)
            self.assertLessEqual(abs(vals[1] - js["cfHk"]), tol)
            self.assertLessEqual(abs(vals[2] - js["cfRt"]), tol)
            self.assertLessEqual(abs(vals[3] - js["cfMsq"]), tol)


if __name__ == "__main__":
    unittest.main()
