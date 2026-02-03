import json
import math
import pathlib
import shutil
import subprocess
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]


class TestXoperParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for Fortran parity checks")

        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "xoper_driver.f"
        srcs = [
            driver_src,
            ROOT / "python" / "tests" / "xpanel_qiset.f",
        ]
        driver_path = pathlib.Path(cls._tmpdir.name) / "xoper_driver"
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

    def test_specal_qiset_parity(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        n = 12
        alfa = 5.0 * math.pi / 180.0
        qinvu = []
        for i in range(n):
            s = i / (n - 1)
            q0 = 1.0 + 0.2 * math.sin(2.0 * math.pi * s)
            q90 = 0.5 + 0.1 * math.cos(2.0 * math.pi * s)
            qinvu.append([q0, q90])

        payload = {
            "qinvu": qinvu,
            "alfa": alfa,
        }

        script = pathlib.Path(__file__).with_name("compare_xoper.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        input_lines = [f"{n} 0", str(alfa)]
        input_lines.extend(f"{pair[0]} {pair[1]}" for pair in qinvu)

        proc_f = subprocess.run(
            [str(self.driver_path)],
            input="\n".join(input_lines) + "\n",
            text=True,
            capture_output=True,
            check=True,
        )
        lines = [line.strip() for line in proc_f.stdout.splitlines() if line.strip()]
        out_qinv = []
        out_qinva = []
        for line in lines:
            if line.startswith("QINV "):
                out_qinv = [float(v) for v in line.split()[1:]]
            elif line.startswith("QINVA"):
                out_qinva = [float(v) for v in line.split()[1:]]

        tol = 1.0e-6
        # JS arrays are length n+1 with index 0 unused.
        for i in range(n):
            self.assertLessEqual(abs(out_qinv[i] - js_results["qinv"][i + 1]), tol)
            self.assertLessEqual(abs(out_qinva[i] - js_results["qinvA"][i + 1]), tol)


if __name__ == "__main__":
    unittest.main()
