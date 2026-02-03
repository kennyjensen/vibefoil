import json
import pathlib
import shutil
import subprocess
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]


def spike_delta(metrics, side):
    before = metrics["before"][side]["spike"]
    after = metrics["after"][side]["spike"]
    if before is None or after is None:
        return None
    return after - before


class TestCpSmoothing2412Fortran(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for Fortran parity checks")
        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "cp_smoothing_fortran_driver.f"
        srcs = [
            driver_src,
            ROOT / "python" / "tests" / "xfoil_subs.f",
            ROOT / "third_party" / "Xfoil" / "src" / "xpanel.f",
            ROOT / "third_party" / "Xfoil" / "src" / "xsolve.f",
            ROOT / "third_party" / "Xfoil" / "src" / "xutils.f",
            ROOT / "third_party" / "Xfoil" / "src" / "spline.f",
            ROOT / "python" / "tests" / "xbl_ref.f",
        ]
        driver_path = pathlib.Path(cls._tmpdir.name) / "cp_smoothing_fortran_driver"
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

    def test_fortran_cp_smoothing_spike(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS smoothing checks")

        build_script = pathlib.Path(__file__).with_name("build_cp_smoothing_input_2412.mjs")
        proc = subprocess.run(
            ["node", str(build_script)],
            text=True,
            capture_output=True,
            check=True,
        )
        payload = json.loads(proc.stdout)

        arrays = payload["arrays"]
        n = payload["nb"]
        lines = [
            str(n),
            f"{payload['alpha']} {payload['minf']} {payload['qinf']}",
            str(payload["sharp"]),
            f"{payload['ante']} {payload['aste']} {payload['dste']}",
            f"{payload['qopi']} {payload['hopi']} {payload['pi']}",
        ]
        for x, y, s, xp, yp, nx, ny in zip(
            arrays["x"], arrays["y"], arrays["s"], arrays["xp"], arrays["yp"], arrays["nx"], arrays["ny"]
        ):
            lines.append(f"{x} {y} {s} {xp} {yp} {nx} {ny}")
        for val in arrays["sig"]:
            lines.append(str(val))
        for q0, q1, q2, q3 in zip(arrays["qf0"], arrays["qf1"], arrays["qf2"], arrays["qf3"]):
            lines.append(f"{q0} {q1} {q2} {q3}")

        proc_f = subprocess.run(
            [str(self.driver_path)],
            input="\n".join(lines) + "\n",
            text=True,
            capture_output=True,
            check=True,
        )
        out_lines = [line.strip() for line in proc_f.stdout.splitlines() if line.strip()]
        numeric_lines = []
        for line in out_lines:
            parts = line.split()
            ok = True
            for tok in parts:
                try:
                    float(tok)
                except ValueError:
                    ok = False
                    break
            if ok:
                numeric_lines.append(line)
        count = int(float(numeric_lines[0]))
        self.assertEqual(count, n)
        points = []
        for line in numeric_lines[1:]:
            parts = line.split()
            if len(parts) < 3:
                continue
            points.append({"x": float(parts[0]), "y": float(parts[1]), "cp": float(parts[2])})
        self.assertEqual(len(points), n)

        compare_script = pathlib.Path(__file__).with_name("compare_cp_smoothing_points.mjs")
        js_metrics = subprocess.run(
            ["node", str(compare_script)],
            input=json.dumps({"points": payload["jsPoints"], "le": payload["le"], "te": payload["te"]}),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(js_metrics.stdout)["results"]

        f_metrics = subprocess.run(
            ["node", str(compare_script)],
            input=json.dumps({"points": points, "le": payload["le"], "te": payload["te"]}),
            text=True,
            capture_output=True,
            check=True,
        )
        f_results = json.loads(f_metrics.stdout)["results"]

        for side in ("upper", "lower"):
            js_delta = spike_delta(js_results, side)
            f_delta = spike_delta(f_results, side)
            if js_delta is None or f_delta is None:
                self.fail(f"Not enough points to compute spike for {side} side")
            js_sign = 0 if abs(js_delta) < 1.0e-6 else (1 if js_delta > 0 else -1)
            f_sign = 0 if abs(f_delta) < 1.0e-6 else (1 if f_delta > 0 else -1)
            self.assertEqual(
                js_sign,
                f_sign,
                msg=(
                    f"Spike change sign mismatch on {side}: "
                    f"js_delta={js_delta}, fortran_delta={f_delta}"
                ),
            )


if __name__ == "__main__":
    unittest.main()
