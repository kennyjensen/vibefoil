import copy
import json
import pathlib
import shutil
import subprocess
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from python.xblsys import blsys
from python.tests.test_compare_xblsys_core import prepare_bl_state, serialize_bl


def pick_array(mat, rows, cols):
    return [[mat[i][j] for j in range(1, cols + 1)] for i in range(1, rows + 1)]


def pick_vector(vec, count):
    return [vec[i] for i in range(1, count + 1)]


class TestBlsysStrictParity(unittest.TestCase):
    def test_blsys_strict(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        _, bl = prepare_bl_state()

        payload = {"base": serialize_bl(bl)}
        script = pathlib.Path(__file__).with_name("compare_blsys_strict.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        bl_py = copy.deepcopy(bl)
        blsys(bl_py)

        py_results = {
            "VS1": pick_array(bl_py.VS1, 4, 5),
            "VS2": pick_array(bl_py.VS2, 4, 5),
            "VSREZ": pick_vector(bl_py.VSREZ, 4),
            "VSM": pick_vector(bl_py.VSM, 4),
            "VSR": pick_vector(bl_py.VSR, 4),
            "VSX": pick_vector(bl_py.VSX, 4),
        }

        tol = 1.0e-9
        for name in ("VS1", "VS2"):
            for i in range(4):
                for j in range(5):
                    self.assertLessEqual(
                        abs(py_results[name][i][j] - js_results[name][i][j]),
                        tol,
                        msg=f"blsys {name}[{i+1}][{j+1}]",
                    )
        for name in ("VSREZ", "VSM", "VSR", "VSX"):
            for i in range(4):
                self.assertLessEqual(
                    abs(py_results[name][i] - js_results[name][i]),
                    tol,
                    msg=f"blsys {name}[{i+1}]",
                )


if __name__ == "__main__":
    unittest.main()
