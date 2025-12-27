import json
import pathlib
import shutil
import subprocess
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from python.xbl import dslim as py_dslim
from python.xblsys import cfl, cft, dil, dilw, dit, hct, hkin, hsl, hst


def max_abs_list_diff(a, b):
    if len(a) != len(b):
        return float("inf")
    return max(abs(ai - bi) for ai, bi in zip(a, b))


class TestXblsysParity(unittest.TestCase):
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

        tol = 1.0e-12

        for i, case in enumerate(payload["hkin"]):
            hk, hk_h, hk_msq = hkin(case["h"], case["msq"])
            js = js_results["hkin"][i]
            self.assertLessEqual(max_abs_list_diff([hk, hk_h, hk_msq], [js["hk"], js["hkH"], js["hkMsq"]]), tol)

        for i, case in enumerate(payload["dil"]):
            di, di_hk, di_rt = dil(case["hk"], case["rt"])
            js = js_results["dil"][i]
            self.assertLessEqual(max_abs_list_diff([di, di_hk, di_rt], [js["di"], js["diHk"], js["diRt"]]), tol)

        for i, case in enumerate(payload["dilw"]):
            di, di_hk, di_rt = dilw(case["hk"], case["rt"])
            js = js_results["dilw"][i]
            self.assertLessEqual(max_abs_list_diff([di, di_hk, di_rt], [js["di"], js["diHk"], js["diRt"]]), tol)

        for i, case in enumerate(payload["hsl"]):
            hs, hs_hk, hs_rt, hs_msq = hsl(case["hk"], case["rt"], case["msq"])
            js = js_results["hsl"][i]
            self.assertLessEqual(
                max_abs_list_diff([hs, hs_hk, hs_rt, hs_msq], [js["hs"], js["hsHk"], js["hsRt"], js["hsMsq"]]),
                tol,
            )

        for i, case in enumerate(payload["cfl"]):
            cf, cf_hk, cf_rt, cf_msq = cfl(case["hk"], case["rt"], case["msq"])
            js = js_results["cfl"][i]
            self.assertLessEqual(
                max_abs_list_diff([cf, cf_hk, cf_rt, cf_msq], [js["cf"], js["cfHk"], js["cfRt"], js["cfMsq"]]),
                tol,
            )

        for i, case in enumerate(payload["dit"]):
            di, di_hs, di_us, di_cf, di_st = dit(case["hs"], case["us"], case["cf"], case["st"])
            js = js_results["dit"][i]
            self.assertLessEqual(
                max_abs_list_diff([di, di_hs, di_us, di_cf, di_st], [js["di"], js["diHs"], js["diUs"], js["diCf"], js["diSt"]]),
                tol,
            )

        for i, case in enumerate(payload["hst"]):
            hs, hs_hk, hs_rt, hs_msq = hst(case["hk"], case["rt"], case["msq"])
            js = js_results["hst"][i]
            self.assertLessEqual(
                max_abs_list_diff([hs, hs_hk, hs_rt, hs_msq], [js["hs"], js["hsHk"], js["hsRt"], js["hsMsq"]]),
                tol,
            )

        for i, case in enumerate(payload["cft"]):
            cf, cf_hk, cf_rt, cf_msq = cft(case["hk"], case["rt"], case["msq"], case["cffac"])
            js = js_results["cft"][i]
            self.assertLessEqual(
                max_abs_list_diff([cf, cf_hk, cf_rt, cf_msq], [js["cf"], js["cfHk"], js["cfRt"], js["cfMsq"]]),
                tol,
            )

        for i, case in enumerate(payload["hct"]):
            hc, hc_hk, hc_msq = hct(case["hk"], case["msq"])
            js = js_results["hct"][i]
            self.assertLessEqual(max_abs_list_diff([hc, hc_hk, hc_msq], [js["hc"], js["hcHk"], js["hcMsq"]]), tol)

        for i, case in enumerate(payload["dslim"]):
            dstr = py_dslim(case["dstr"], case["thet"], case["uedg"], case["msq"], case["hklim"])
            js = js_results["dslim"][i]
            self.assertLessEqual(abs(dstr - js), tol)


if __name__ == "__main__":
    unittest.main()
