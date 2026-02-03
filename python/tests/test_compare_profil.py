import json
import pathlib
import shutil
import subprocess
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]


class TestProfilParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for Fortran parity checks")

        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "profil_driver.f"
        srcs = [
            driver_src,
            ROOT / "python" / "tests" / "blu_subs.f",
            ROOT / "python" / "tests" / "blu_cft.f",
            ROOT / "third_party" / "Xfoil" / "src" / "profil.f",
        ]
        driver_path = pathlib.Path(cls._tmpdir.name) / "profil_driver"
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

    def test_profil_functions(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        prwall_cases = [
            {"dstar": 0.02, "theta": 0.01, "uo": 1.0, "rt": 500.0, "ms": 0.0, "ct": 0.0, "cffac": 1.0},
            {"dstar": 0.03, "theta": 0.012, "uo": 1.0, "rt": 1200.0, "ms": 0.2, "ct": 0.0, "cffac": 1.0},
        ]

        fs_cases = [
            {"inorm": 3, "ispec": 2, "bspec": 0.0, "hspec": 2.8, "n": 30, "etae": 6.0, "geo": 1.0},
            {"inorm": 3, "ispec": 2, "bspec": 0.0, "hspec": 3.6, "n": 40, "etae": 7.0, "geo": 1.0},
        ]

        # Use JS prwall outputs to seed UWALL comparison inputs.
        payload = {
            "prwall": prwall_cases,
            "uwall": [],
            "fs": fs_cases,
        }
        script = pathlib.Path(__file__).with_name("compare_profil.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        uwall_cases = []
        for idx, case in enumerate(prwall_cases):
            pr = js_results["prwall"][idx]
            uwall_cases.append(
                {
                    "th": case["theta"],
                    "uo": case["uo"],
                    "de": pr["de"],
                    "us": pr["us"],
                    "rt": case["rt"],
                    "cf": pr["cf"],
                    "bb": pr["bb"],
                    "n": 32,
                }
            )

        payload["uwall"] = uwall_cases
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        input_lines = [f"{len(prwall_cases)} {len(uwall_cases)} {len(fs_cases)}"]
        input_lines.extend(
            f"{c['dstar']} {c['theta']} {c['uo']} {c['rt']} {c['ms']} {c['ct']} {c['cffac']}"
            for c in prwall_cases
        )
        input_lines.extend(
            f"{c['th']} {c['uo']} {c['de']} {c['us']} {c['rt']} {c['cf']} {c['bb']} {c['n']}"
            for c in uwall_cases
        )
        input_lines.extend(
            f"{c['inorm']} {c['ispec']} {c['bspec']} {c['hspec']} {c['n']} {c['etae']} {c['geo']}"
            for c in fs_cases
        )

        proc_f = subprocess.run(
            [str(self.driver_path)],
            input="\n".join(input_lines) + "\n",
            text=True,
            capture_output=True,
            check=True,
        )
        lines = [line.strip() for line in proc_f.stdout.splitlines() if line.strip()]

        out_prwall = []
        out_uwall = []
        out_fs = []
        out_fs_delta = []

        for line in lines:
            if line.startswith("PRWALL"):
                parts = line.split()
                out_prwall.append([float(v) for v in parts[1:]])
            elif line.startswith("UWALL"):
                parts = line.split()[1:]
                vals = [float(v) for v in parts]
                out_uwall.append(vals)
            elif line.startswith("FSDELTA"):
                out_fs_delta.append(float(line.split()[1]))
            elif line.startswith("FS "):
                parts = line.split()[1:]
                vals = [float(v) for v in parts]
                out_fs.append(vals)
            else:
                # Skip informational lines like "FS: Convergence failed"
                continue

        tol = 1.0e-6
        for idx, js in enumerate(js_results["prwall"]):
            vals = out_prwall[idx]
            self.assertLessEqual(abs(vals[0] - js["de"]), tol)
            self.assertLessEqual(abs(vals[1] - js["us"]), tol)
            self.assertLessEqual(abs(vals[2] - js["cf"]), tol)
            self.assertLessEqual(abs(vals[3] - js["bb"]), tol)

        for idx, js in enumerate(js_results["uwall"]):
            vals = out_uwall[idx]
            # Fortran outputs N pairs; JS arrays are 1-based with index 0 unused.
            y_js = js["y"][1:]
            u_js = js["u"][1:]
            for j in range(len(y_js)):
                y_f = vals[2 * j]
                u_f = vals[2 * j + 1]
                self.assertLessEqual(abs(y_f - y_js[j]), tol)
                self.assertLessEqual(abs(u_f - u_js[j]), tol)

        for idx, js in enumerate(js_results["fs"]):
            vals = out_fs[idx]
            self.assertLessEqual(abs(out_fs_delta[idx] - js["delta"]), tol)
            eta_js = js["eta"][1:]
            f_js = js["f"][1:]
            u_js = js["u"][1:]
            s_js = js["s"][1:]
            for j in range(len(eta_js)):
                base = 4 * j
                self.assertLessEqual(abs(vals[base] - eta_js[j]), tol)
                self.assertLessEqual(abs(vals[base + 1] - f_js[j]), tol)
                self.assertLessEqual(abs(vals[base + 2] - u_js[j]), tol)
                self.assertLessEqual(abs(vals[base + 3] - s_js[j]), tol)


if __name__ == "__main__":
    unittest.main()
