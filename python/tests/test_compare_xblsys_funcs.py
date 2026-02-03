import json
import pathlib
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def max_abs_list_diff(a, b):
    if len(a) != len(b):
        return float("inf")
    return max(abs(ai - bi) for ai, bi in zip(a, b))


class TestXblsysParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for Fortran parity checks")

        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "xblsys_driver.f"
        srcs = [
            driver_src,
            ROOT / "python" / "tests" / "xbl_subs.f",
            ROOT / "third_party" / "Xfoil" / "src" / "xblsys.f",
        ]
        driver_path = pathlib.Path(cls._tmpdir.name) / "xblsys_driver"
        subprocess.run(
            [
                "gfortran",
                "-O2",
                "-ffixed-form",
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

    def test_xblsys_functions(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        payload = {
            "hkin": [
                {"h": 1.5, "msq": 0.0},
                {"h": 2.2, "msq": 0.2},
                {"h": 4.5, "msq": 0.5},
            ],
            "dil": [
                {"hk": 3.0, "rt": 150.0},
                {"hk": 4.0, "rt": 500.0},
                {"hk": 6.5, "rt": 1500.0},
            ],
            "dilw": [
                {"hk": 2.5, "rt": 200.0},
                {"hk": 4.2, "rt": 800.0},
                {"hk": 6.0, "rt": 2000.0},
            ],
            "hsl": [
                {"hk": 3.9, "rt": 120.0, "msq": 0.0},
                {"hk": 4.35, "rt": 500.0, "msq": 0.2},
                {"hk": 5.2, "rt": 900.0, "msq": 0.3},
            ],
            "cfl": [
                {"hk": 3.9, "rt": 120.0, "msq": 0.0},
                {"hk": 5.5, "rt": 500.0, "msq": 0.2},
                {"hk": 6.2, "rt": 900.0, "msq": 0.3},
            ],
            "dit": [
                {"hs": 1.5, "us": 0.2, "cf": 0.003, "st": 0.02},
                {"hs": 2.0, "us": 0.4, "cf": 0.001, "st": 0.03},
                {"hs": 1.8, "us": 0.7, "cf": 0.004, "st": 0.01},
            ],
            "hst": [
                {"hk": 3.0, "rt": 150.0, "msq": 0.0},
                {"hk": 4.0, "rt": 500.0, "msq": 0.2},
                {"hk": 6.5, "rt": 1500.0, "msq": 0.5},
            ],
            "cft": [
                {"hk": 3.8, "rt": 200.0, "msq": 0.0, "cffac": 1.0},
                {"hk": 5.0, "rt": 800.0, "msq": 0.2, "cffac": 1.0},
                {"hk": 6.5, "rt": 2000.0, "msq": 0.3, "cffac": 1.0},
                {"hk": 2.2, "rt": 50.0, "msq": 0.1, "cffac": 1.0},
                {"hk": 20.0, "rt": 1000.0, "msq": 0.2, "cffac": 1.0},
            ],
            "hct": [
                {"hk": 2.0, "msq": 0.1},
                {"hk": 3.5, "msq": 0.3},
                {"hk": 6.0, "msq": 0.5},
            ],
            "dslim": [
                {"dstr": 0.004, "thet": 0.002, "uedg": 0.9, "msq": 0.2, "hklim": 1.02},
                {"dstr": 0.002, "thet": 0.001, "uedg": 1.1, "msq": 0.4, "hklim": 1.00005},
                {"dstr": 0.006, "thet": 0.003, "uedg": 0.7, "msq": 0.1, "hklim": 1.02},
                {"dstr": 0.01, "thet": 0.002, "uedg": 1.0, "msq": 0.2, "hklim": 2.0},
            ],
        }

        script = pathlib.Path(__file__).with_name("compare_xblsys_funcs.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        # Run Fortran driver
        counts = [
            len(payload["hkin"]),
            len(payload["dil"]),
            len(payload["dilw"]),
            len(payload["hsl"]),
            len(payload["cfl"]),
            len(payload["dit"]),
            len(payload["hst"]),
            len(payload["cft"]),
            len(payload["hct"]),
            len(payload["dslim"]),
        ]
        input_lines = [" ".join(str(c) for c in counts)]
        input_lines.extend(f"{case['h']} {case['msq']}" for case in payload["hkin"])
        input_lines.extend(f"{case['hk']} {case['rt']}" for case in payload["dil"])
        input_lines.extend(f"{case['hk']} {case['rt']}" for case in payload["dilw"])
        input_lines.extend(f"{case['hk']} {case['rt']} {case['msq']}" for case in payload["hsl"])
        input_lines.extend(f"{case['hk']} {case['rt']} {case['msq']}" for case in payload["cfl"])
        input_lines.extend(f"{case['hs']} {case['us']} {case['cf']} {case['st']}" for case in payload["dit"])
        input_lines.extend(f"{case['hk']} {case['rt']} {case['msq']}" for case in payload["hst"])
        input_lines.extend(f"{case['hk']} {case['rt']} {case['msq']} {case['cffac']}" for case in payload["cft"])
        input_lines.extend(f"{case['hk']} {case['msq']}" for case in payload["hct"])
        input_lines.extend(
            f"{case['dstr']} {case['thet']} {case['uedg']} {case['msq']} {case['hklim']}"
            for case in payload["dslim"]
        )

        proc_f = subprocess.run(
            [str(self.driver_path)],
            input="\n".join(input_lines) + "\n",
            text=True,
            capture_output=True,
            check=True,
        )
        lines = [line.strip() for line in proc_f.stdout.splitlines() if line.strip()]
        out = {
            "hkin": [],
            "dil": [],
            "dilw": [],
            "hsl": [],
            "cfl": [],
            "dit": [],
            "hst": [],
            "cft": [],
            "hct": [],
            "dslim": [],
        }
        for line in lines:
            parts = line.split()
            tag = parts[0].lower()
            vals = [float(v) for v in parts[1:]]
            out[tag].append(vals)

        tol = 1.0e-3
        tol_cft = 1.0e-2

        for i, js in enumerate(js_results["hkin"]):
            self.assertLessEqual(max_abs_list_diff(out["hkin"][i], [js["hk"], js["hkH"], js["hkMsq"]]), tol)

        for i, js in enumerate(js_results["dil"]):
            self.assertLessEqual(max_abs_list_diff(out["dil"][i], [js["di"], js["diHk"], js["diRt"]]), tol)

        for i, js in enumerate(js_results["dilw"]):
            self.assertLessEqual(max_abs_list_diff(out["dilw"][i], [js["di"], js["diHk"], js["diRt"]]), tol)

        for i, js in enumerate(js_results["hsl"]):
            self.assertLessEqual(max_abs_list_diff(out["hsl"][i], [js["hs"], js["hsHk"], js["hsRt"], js["hsMsq"]]), tol)

        for i, js in enumerate(js_results["cfl"]):
            self.assertLessEqual(max_abs_list_diff(out["cfl"][i], [js["cf"], js["cfHk"], js["cfRt"], js["cfMsq"]]), tol)

        for i, js in enumerate(js_results["dit"]):
            self.assertLessEqual(
                max_abs_list_diff(out["dit"][i], [js["di"], js["diHs"], js["diUs"], js["diCf"], js["diSt"]]),
                tol,
            )

        for i, js in enumerate(js_results["hst"]):
            self.assertLessEqual(max_abs_list_diff(out["hst"][i], [js["hs"], js["hsHk"], js["hsRt"], js["hsMsq"]]), tol)

        for i, js in enumerate(js_results["cft"]):
            self.assertLessEqual(
                max_abs_list_diff(out["cft"][i], [js["cf"], js["cfHk"], js["cfRt"], js["cfMsq"]]),
                tol_cft,
            )

        for i, js in enumerate(js_results["hct"]):
            self.assertLessEqual(max_abs_list_diff(out["hct"][i], [js["hc"], js["hcHk"], js["hcMsq"]]), tol)

        for i, js_val in enumerate(js_results["dslim"]):
            self.assertLessEqual(abs(out["dslim"][i][0] - js_val), tol)


if __name__ == "__main__":
    unittest.main()
