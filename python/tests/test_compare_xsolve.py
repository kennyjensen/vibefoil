import json
import math
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


def flatten_matrix(mat):
    return [item for row in mat for item in row]


class TestXsolveParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("node"):
            raise unittest.SkipTest("node is required for JS/Fortran parity checks")
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for Fortran parity checks")

        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "xsolve_driver.f90"
        srcs = [
            driver_src,
            ROOT / "third_party" / "Xfoil" / "src" / "xsolve.f",
        ]
        driver_path = pathlib.Path(cls._tmpdir.name) / "xsolve_driver"
        subprocess.run(
            [
                "gfortran",
                "-O2",
                "-ffixed-form",
                "-I",
                str(ROOT / "third_party" / "Xfoil" / "src"),
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

    def run_fortran(self, input_str):
        proc = subprocess.run(
            [str(self.driver_path)],
            input=input_str,
            text=True,
            capture_output=True,
            check=True,
        )
        return proc.stdout.strip().splitlines()

    def run_js(self, payload):
        script = pathlib.Path(__file__).with_name("compare_xsolve.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        return json.loads(proc.stdout)["results"]

    def test_gauss_parity(self):
        n = 4
        nrhs = 2
        z = []
        r = []
        for i in range(n):
            row = []
            for j in range(n):
                val = (4.0 + 0.2 * (i + 1)) if i == j else 0.1 * (i + j + 1)
                row.append(val)
            z.append(row)
            r.append([1.0 + 0.1 * (i + 1), 0.5 + 0.05 * (i + 1)])

        payload = {"gauss": {"n": n, "nrhs": nrhs, "z": z, "r": r}}
        js_results = self.run_js(payload)["gauss"]

        lines = self.run_fortran(
            "\n".join(
                [
                    "1",
                    f"{n} {n} {nrhs}",
                    *[str(val) for val in flatten_matrix(z)],
                    *[str(val) for val in flatten_matrix(r)],
                ]
            )
            + "\n"
        )
        header = lines[0].strip().split()
        self.assertEqual(int(header[0]), n)
        self.assertEqual(int(header[1]), nrhs)
        z_out = []
        r_out = []
        idx = 1
        for _ in range(n):
            parts = lines[idx].strip().split()
            z_out.append([float(p) for p in parts])
            idx += 1
        r_flat = []
        for _ in range(n):
            parts = lines[idx].strip().split()
            r_flat.extend(float(p) for p in parts)
            idx += 1
        expected_len = n * nrhs
        r_flat = r_flat[:expected_len]
        for i in range(n):
            r_out.append(r_flat[i * nrhs : (i + 1) * nrhs])

        tol = 1.0e-6
        self.assertLessEqual(max_abs_diff(flatten_matrix(z_out), flatten_matrix(js_results["z"])), tol)
        self.assertLessEqual(max_abs_diff(flatten_matrix(r_out), flatten_matrix(js_results["r"])), tol)

    def test_ludcmp_baksub_parity(self):
        n = 5
        a = []
        b = []
        for i in range(n):
            row = []
            for j in range(n):
                val = (3.0 + 0.3 * (i + 1)) if i == j else 0.05 * (i + j + 1)
                row.append(val)
            a.append(row)
            b.append(0.2 + 0.1 * (i + 1))

        payload = {"ludcmp": {"n": n, "a": a, "b": b}}
        js_results = self.run_js(payload)["ludcmp"]

        lines = self.run_fortran(
            "\n".join(
                [
                    "2",
                    f"{n} {n}",
                    *[str(val) for val in flatten_matrix(a)],
                    *[str(val) for val in b],
                ]
            )
            + "\n"
        )
        n_out = int(lines[0].strip())
        self.assertEqual(n_out, n)
        idx = 1
        a_out = []
        for _ in range(n):
            parts = lines[idx].strip().split()
            a_out.append([float(p) for p in parts])
            idx += 1
        indx_out = [int(val) - 1 for val in lines[idx].strip().split()]
        idx += 1
        b_out = [float(val) for val in lines[idx].strip().split()]

        tol = 1.0e-6
        self.assertLessEqual(max_abs_diff(flatten_matrix(a_out), flatten_matrix(js_results["a"])), tol)
        self.assertEqual(indx_out, js_results["indx"])
        self.assertLessEqual(max_abs_diff(b_out, js_results["b"]), tol)

    def test_blsolv_parity(self):
        nsys = 4
        iblte1 = 1
        iblte2 = 1
        ivte1 = 2
        ivz = 3
        vaccel = 0.01
        s1 = 0.0
        s2 = 1.0

        def va_val(k, j, i):
            if k == 1 and j == 1:
                return 2.0 + 0.1 * i
            if k == 2 and j == 2:
                return 1.5 + 0.07 * i
            return 0.05 * k + 0.02 * j + 0.01 * i

        def vb_val(k, j, i):
            return 0.02 * k + 0.01 * j + 0.005 * i

        def vm_val(k, l, i):
            if k == 3 and l == i:
                return 2.2 + 0.08 * i
            return 0.01 * k + 0.002 * l + 0.003 * i

        def vdel_val(k, j, i):
            return 0.03 * k + 0.01 * j + 0.004 * i

        def vz_val(k, j):
            return 0.02 * k + 0.01 * j

        va = [
            va_val(k, j, i)
            for k in range(1, 4)
            for j in range(1, 3)
            for i in range(1, nsys + 1)
        ]
        vb = [
            vb_val(k, j, i)
            for k in range(1, 4)
            for j in range(1, 3)
            for i in range(1, nsys + 1)
        ]
        vm = [
            vm_val(k, l, i)
            for k in range(1, 4)
            for l in range(1, nsys + 1)
            for i in range(1, nsys + 1)
        ]
        vdel = [
            vdel_val(k, j, i)
            for k in range(1, 4)
            for j in range(1, 3)
            for i in range(1, nsys + 1)
        ]
        vz = [vz_val(k, j) for k in range(1, 4) for j in range(1, 3)]

        payload = {
            "blsolv": {
                "nsys": nsys,
                "iblte1": iblte1,
                "iblte2": iblte2,
                "ivte1": ivte1,
                "ivz": ivz,
                "vaccel": vaccel,
                "s1": s1,
                "s2": s2,
                "va": va,
                "vb": vb,
                "vm": vm,
                "vdel": vdel,
                "vz": vz,
            }
        }
        js_results = self.run_js(payload)["blsolv"]

        lines = self.run_fortran(
            "\n".join(
                [
                    "3",
                    f"{nsys}",
                    f"{iblte1} {iblte2} {ivte1} {ivz}",
                    f"{vaccel}",
                    f"{s1} {s2}",
                    *[str(val) for val in va],
                    *[str(val) for val in vb],
                    *[str(val) for val in vm],
                    *[str(val) for val in vdel],
                    *[str(val) for val in vz],
                ]
            )
            + "\n"
        )
        nsys_out = int(lines[0].strip())
        self.assertEqual(nsys_out, nsys)
        idx = 1
        nva = int(lines[idx].strip())
        idx += 1
        va_out = [float(lines[idx + i].strip()) for i in range(nva)]
        idx += nva
        nvm = int(lines[idx].strip())
        idx += 1
        vm_out = [float(lines[idx + i].strip()) for i in range(nvm)]
        idx += nvm
        nvdel = int(lines[idx].strip())
        idx += 1
        vdel_out = [float(lines[idx + i].strip()) for i in range(nvdel)]

        tol = 1.0e-6
        self.assertLessEqual(max_abs_diff(va_out, js_results["va"]), tol)
        self.assertLessEqual(max_abs_diff(vm_out, js_results["vm"]), tol)
        self.assertLessEqual(max_abs_diff(vdel_out, js_results["vdel"]), tol)


if __name__ == "__main__":
    unittest.main()
