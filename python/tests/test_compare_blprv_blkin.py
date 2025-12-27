import json
import pathlib
import shutil
import subprocess
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from python.xbl import XBlState
from python.xblsys import blkin, blprv


def build_cases():
    return [
        {
            "constants": {
                "TKBL": 0.1,
                "TKBL_MS": 0.005,
                "QINFBL": 1.0,
                "HSTINV": 0.2,
                "HSTINV_MS": 0.01,
                "GM1BL": 0.4,
                "RSTBL": 1.1,
                "RSTBL_MS": 0.02,
                "HVRAT": 0.2,
                "REYBL": 1.0e6,
                "REYBL_MS": 1.0e4,
                "REYBL_RE": 1.0,
            },
            "inputs": {
                "xsi": 0.1,
                "ami": 0.02,
                "cti": 0.003,
                "thi": 0.0015,
                "dsi": 0.002,
                "dswaki": 0.0001,
                "uei": 0.8,
            },
        },
        {
            "constants": {
                "TKBL": 0.05,
                "TKBL_MS": 0.002,
                "QINFBL": 1.2,
                "HSTINV": 0.15,
                "HSTINV_MS": 0.006,
                "GM1BL": 0.4,
                "RSTBL": 1.05,
                "RSTBL_MS": 0.015,
                "HVRAT": 0.1,
                "REYBL": 2.0e6,
                "REYBL_MS": 5.0e4,
                "REYBL_RE": 1.5,
            },
            "inputs": {
                "xsi": 0.25,
                "ami": 0.015,
                "cti": 0.0025,
                "thi": 0.0012,
                "dsi": 0.0018,
                "dswaki": 0.0002,
                "uei": 1.0,
            },
        },
        {
            "constants": {
                "TKBL": 0.2,
                "TKBL_MS": 0.008,
                "QINFBL": 0.9,
                "HSTINV": 0.25,
                "HSTINV_MS": 0.012,
                "GM1BL": 0.4,
                "RSTBL": 1.08,
                "RSTBL_MS": 0.018,
                "HVRAT": 0.25,
                "REYBL": 1.5e6,
                "REYBL_MS": 7.5e4,
                "REYBL_RE": 0.85,
            },
            "inputs": {
                "xsi": 0.4,
                "ami": 0.01,
                "cti": 0.0035,
                "thi": 0.0018,
                "dsi": 0.0022,
                "dswaki": 0.00015,
                "uei": 0.7,
            },
        },
        {
            "constants": {
                "TKBL": 0.3,
                "TKBL_MS": 0.01,
                "QINFBL": 1.1,
                "HSTINV": 0.18,
                "HSTINV_MS": 0.007,
                "GM1BL": 0.4,
                "RSTBL": 1.2,
                "RSTBL_MS": 0.03,
                "HVRAT": 0.15,
                "REYBL": 8.0e5,
                "REYBL_MS": 2.5e4,
                "REYBL_RE": 0.95,
            },
            "inputs": {
                "xsi": 0.6,
                "ami": 0.03,
                "cti": 0.004,
                "thi": 0.0021,
                "dsi": 0.0026,
                "dswaki": 0.00012,
                "uei": 0.85,
            },
        },
        {
            "constants": {
                "TKBL": 0.12,
                "TKBL_MS": 0.004,
                "QINFBL": 0.95,
                "HSTINV": 0.22,
                "HSTINV_MS": 0.009,
                "GM1BL": 0.4,
                "RSTBL": 1.03,
                "RSTBL_MS": 0.012,
                "HVRAT": 0.3,
                "REYBL": 2.5e6,
                "REYBL_MS": 6.5e4,
                "REYBL_RE": 1.2,
            },
            "inputs": {
                "xsi": 0.75,
                "ami": 0.012,
                "cti": 0.0022,
                "thi": 0.0014,
                "dsi": 0.0016,
                "dswaki": 0.00018,
                "uei": 0.95,
            },
        },
    ]


class TestBlprvBlkinParity(unittest.TestCase):
    def test_blprv_blkin(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        fields = [
            "X2",
            "AMPL2",
            "S2",
            "T2",
            "D2",
            "DW2",
            "U2",
            "U2_UEI",
            "U2_MS",
            "M2",
            "M2_U2",
            "M2_MS",
            "R2",
            "R2_U2",
            "R2_MS",
            "H2",
            "H2_D2",
            "H2_T2",
            "V2",
            "V2_U2",
            "V2_MS",
            "V2_RE",
            "HK2",
            "HK2_U2",
            "HK2_T2",
            "HK2_D2",
            "HK2_MS",
            "RT2",
            "RT2_U2",
            "RT2_T2",
            "RT2_MS",
            "RT2_RE",
        ]

        cases = build_cases()
        payload = {"cases": cases, "fields": fields}

        script = pathlib.Path(__file__).with_name("compare_blprv_blkin.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        tol = 1.0e-12

        for idx, case in enumerate(cases):
            bl = XBlState()
            for key, value in case["constants"].items():
                setattr(bl, key, value)

            inp = case["inputs"]
            blprv(bl, inp["xsi"], inp["ami"], inp["cti"], inp["thi"], inp["dsi"], inp["dswaki"], inp["uei"])
            blkin(bl)

            js_case = js_results[idx]
            for field in fields:
                py_val = getattr(bl, field)
                js_val = js_case[field]
                self.assertLessEqual(abs(py_val - js_val), tol, msg=f"Case {idx} field {field}")


if __name__ == "__main__":
    unittest.main()
