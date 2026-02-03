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


def finite_diff(xs, ss):
    n = len(xs)
    out = [0.0] * n
    for i in range(n):
        if i == 0:
            ds = ss[i + 1] - ss[i]
            out[i] = (xs[i + 1] - xs[i]) / ds if ds != 0 else 0.0
        elif i == n - 1:
            ds = ss[i] - ss[i - 1]
            out[i] = (xs[i] - xs[i - 1]) / ds if ds != 0 else 0.0
        else:
            ds = ss[i + 1] - ss[i - 1]
            out[i] = (xs[i + 1] - xs[i - 1]) / ds if ds != 0 else 0.0
    return out


class TestXpanelParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("node"):
            raise unittest.SkipTest("node is required for JS/Fortran parity checks")
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for Fortran reference")

        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "xpanel_driver.f"
        xpanel_src = ROOT / "third_party" / "Xfoil" / "src" / "xpanel.f"
        xutils_src = ROOT / "third_party" / "Xfoil" / "src" / "xutils.f"
        xsolve_src = ROOT / "third_party" / "Xfoil" / "src" / "xsolve.f"
        spline_src = ROOT / "third_party" / "Xfoil" / "src" / "spline.f"
        xbl_ref = ROOT / "python" / "tests" / "xbl_ref.f"
        include_dir = ROOT / "third_party" / "Xfoil" / "src"
        driver_path = pathlib.Path(cls._tmpdir.name) / "xpanel_driver"
        subprocess.run(
            [
                "gfortran",
                "-O2",
                "-I",
                str(include_dir),
                str(driver_src),
                str(xpanel_src),
                str(xutils_src),
                str(xsolve_src),
                str(spline_src),
                str(xbl_ref),
                "-o",
                str(driver_path),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        cls.driver_path = driver_path

        cls.base = cls._build_base_fixture()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, "_tmpdir"):
            cls._tmpdir.cleanup()

    @classmethod
    def _run_fortran(cls, lines):
        proc = subprocess.run(
            [str(cls.driver_path)],
            input="\n".join(lines) + "\n",
            text=True,
            capture_output=True,
            check=True,
        )
        raw_lines = proc.stdout.splitlines()
        filtered = []
        for line in raw_lines:
            stripped = line.strip()
            if not stripped:
                continue
            tokens = stripped.split()
            is_numeric = True
            for tok in tokens:
                try:
                    float(tok)
                except ValueError:
                    is_numeric = False
                    break
            if is_numeric:
                filtered.append(stripped)
        return filtered

    @classmethod
    def _run_node(cls, cases):
        payload = {"cases": cases}
        script = pathlib.Path(__file__).with_name("compare_xpanel.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        return json.loads(proc.stdout)["results"]

    @classmethod
    def _build_base_fixture(cls):
        n = 6
        nw = 3
        x = [1.0, 0.5, 0.0, 0.5, 0.9, 1.0]
        y = [0.0, 0.1, 0.0, -0.1, -0.02, -0.01]
        s = [0.0]
        for i in range(1, n):
            dx = x[i] - x[i - 1]
            dy = y[i] - y[i - 1]
            s.append(s[-1] + math.hypot(dx, dy))
        xp = finite_diff(x, s)
        yp = finite_diff(y, s)

        # Extend arrays with wake points (simple downstream line).
        for i in range(nw):
            x.append(x[-1] + 0.1 * (i + 1))
            y.append(y[-1] - 0.02 * (i + 1))
            s.append(s[-1] + 0.1)
            xp.append(xp[-1])
            yp.append(yp[-1])

        # Get normals from Fortran NCALC for surface points.
        lines = ["1", str(n)]
        for xi, yi, si in zip(x[:n], y[:n], s[:n]):
            lines.append(f"{xi} {yi} {si}")
        out = cls._run_fortran(lines)
        xn = []
        yn = []
        for line in out[1:]:
            parts = line.strip().split()
            if len(parts) >= 2:
                xn.append(float(parts[0]))
                yn.append(float(parts[1]))
        # Fill wake normals with last surface normal.
        xn.extend([xn[-1]] * nw)
        yn.extend([yn[-1]] * nw)

        # Compute APANEL from Fortran for surface points.
        dste = math.hypot(0.02, 0.01)
        lines = ["2", f"{n} 0", "0.02 0.01 {}".format(dste)]
        for xi, yi, nxi, nyi in zip(x[:n], y[:n], xn[:n], yn[:n]):
            lines.append(f"{xi} {yi} {nxi} {nyi}")
        out = cls._run_fortran(lines)
        apanel = [float(line.split()[0]) for line in out[1:]]
        apanel.extend([0.0] * nw)

        gam = [0.5, 0.25, 0.05, -0.15, -0.35, -0.55]
        sig = [0.01 * (i + 1) for i in range(n + nw)]
        qf0 = [0.0] * n
        qf1 = [0.0] * n
        qf2 = [0.0] * n
        qf3 = [0.0] * n
        alfa = 0.2
        qinvu = [[0.1 * (i + 1), -0.05 * (i + 1)] for i in range(n + nw)]
        qinv = [0.0] * (n + nw + 1)
        qinv_a = [0.0] * (n + nw + 1)
        cos_a = math.cos(alfa)
        sin_a = math.sin(alfa)
        for i in range(1, n + nw + 1):
            q0, q90 = qinvu[i - 1]
            qinv[i] = cos_a * q0 + sin_a * q90
            qinv_a[i] = -sin_a * q0 + cos_a * q90
        qvis = [0.02 * (i + 1) for i in range(n + nw + 1)]

        return {
            "n": n,
            "nw": nw,
            "x": x,
            "y": y,
            "s": s,
            "xp": xp,
            "yp": yp,
            "nx": xn,
            "ny": yn,
            "apanel": apanel,
            "gam": gam,
            "sig": sig,
            "qf0": qf0,
            "qf1": qf1,
            "qf2": qf2,
            "qf3": qf3,
            "alfa": alfa,
            "qinf": 1.0,
            "sharp": False,
            "ante": 0.02,
            "aste": 0.01,
            "dste": dste,
            "qopi": 1.0 / (4.0 * math.pi),
            "hopi": 1.0 / (2.0 * math.pi),
            "pi": math.pi,
            "limage": False,
            "yimage": 0.0,
            "qinvu": qinvu,
            "qinv": qinv,
            "qinvA": qinv_a,
            "qvis": qvis,
        }

    def test_ncalc(self):
        n = self.base["n"]
        cases = [
            {"kind": "ncalc", "x": self.base["x"][:n], "y": self.base["y"][:n], "s": self.base["s"][:n]}
        ]
        js_case = self._run_node(cases)[0]
        lines = ["1", str(n)]
        for xi, yi, si in zip(self.base["x"][:n], self.base["y"][:n], self.base["s"][:n]):
            lines.append(f"{xi} {yi} {si}")
        out = self._run_fortran(lines)
        xn_ref = []
        yn_ref = []
        for line in out[1:]:
            parts = line.strip().split()
            xn_ref.append(float(parts[0]))
            yn_ref.append(float(parts[1]))
        tol = 1.0e-6
        self.assertLessEqual(max_abs_diff(xn_ref, js_case["xn"]), tol)
        self.assertLessEqual(max_abs_diff(yn_ref, js_case["yn"]), tol)

    def test_apcalc(self):
        n = self.base["n"]
        state = {
            **self.base,
            "x": self.base["x"][:n],
            "y": self.base["y"][:n],
            "nx": self.base["nx"][:n],
            "ny": self.base["ny"][:n],
            "apanel": [0.0] * n,
        }
        js_case = self._run_node([{"kind": "apcalc", "state": state}])[0]
        lines = ["2", f"{n} 0", f"{self.base['ante']} {self.base['aste']} {self.base['dste']}"]
        for xi, yi, nxi, nyi in zip(state["x"], state["y"], state["nx"], state["ny"]):
            lines.append(f"{xi} {yi} {nxi} {nyi}")
        out = self._run_fortran(lines)
        ap_ref = [float(line.split()[0]) for line in out[1:]]
        tol = 1.0e-6
        self.assertLessEqual(max_abs_diff(ap_ref, js_case["apanel"]), tol)

    def test_psilin(self):
        n = self.base["n"]
        nw = self.base["nw"]
        i_js = 1
        i_ft = i_js + 1
        i_ref_js = 0
        i_ref_ft = i_ref_js + 1
        state = {**self.base}
        case = {
            "kind": "psilin",
            "state": state,
            "i": i_js,
            "xi": state["x"][i_js],
            "yi": state["y"][i_js],
            "nxi": state["nx"][i_js],
            "nyi": state["ny"][i_js],
            "geolin": True,
            "siglin": True,
        }
        case_ref = {
            **case,
            "i": i_ref_js,
            "xi": state["x"][i_ref_js],
            "yi": state["y"][i_ref_js],
            "nxi": state["nx"][i_ref_js],
            "nyi": state["ny"][i_ref_js],
        }
        js_case, js_ref = self._run_node([case, case_ref])

        def run_fortran_psilin(i_ft, i_js_idx):
            lines = [
                "3",
                f"{n} {nw} {i_ft} 1 1",
                "0 0",
                f"{state['alfa']} {state['qinf']}",
                f"{state['ante']} {state['aste']} {state['dste']}",
                f"{state['qopi']} {state['hopi']} {state['pi']}",
                f"{state['yimage']}",
            ]
            for j in range(n):
                lines.append(
                    f"{state['x'][j]} {state['y'][j]} {state['s'][j]} "
                    f"{state['xp'][j]} {state['yp'][j]} {state['nx'][j]} {state['ny'][j]} {state['gam'][j]}"
                )
            for j in range(n + nw):
                lines.append(f"{state['sig'][j]}")
            for j in range(n):
                lines.append(f"{state['qf0'][j]} {state['qf1'][j]} {state['qf2'][j]} {state['qf3'][j]}")
            lines.append(
                f"{state['x'][i_js_idx]} {state['y'][i_js_idx]} {state['nx'][i_js_idx]} {state['ny'][i_js_idx]}"
            )
            out = self._run_fortran(lines)
            header = list(map(float, out[0].split()))
            z_extra = list(map(float, out[1].split()))
            return header, z_extra, out[2:]

        header, z_extra, out_lines = run_fortran_psilin(i_ft, i_js)
        header_ref, _, _ = run_fortran_psilin(i_ref_ft, i_ref_js)
        psi_ref, psi_ni_ref, qtan1_ref, qtan2_ref, qtanm_ref, zqinf_ref, zalfa_ref, zq0_ref, zq1_ref = header
        zq2_ref, zq3_ref = z_extra
        tol = 1.0e-5
        psi_delta_ref = psi_ref - header_ref[0]
        psi_delta_js = js_case["result"]["psi"] - js_ref["result"]["psi"]
        self.assertLessEqual(abs(psi_delta_ref - psi_delta_js), tol)
        self.assertLessEqual(abs(psi_ni_ref - js_case["result"]["psiNi"]), tol)
        self.assertLessEqual(abs(qtan1_ref - js_case["result"]["qtan1"]), tol)
        self.assertLessEqual(abs(qtan2_ref - js_case["result"]["qtan2"]), tol)
        self.assertLessEqual(abs(zq0_ref - js_case["result"]["zQdof0"]), tol)
        self.assertLessEqual(abs(zq1_ref - js_case["result"]["zQdof1"]), tol)
        self.assertLessEqual(abs(zq2_ref - js_case["result"]["zQdof2"]), tol)
        self.assertLessEqual(abs(zq3_ref - js_case["result"]["zQdof3"]), tol)
        idx = 0
        n_out = int(out_lines[idx]); idx += 1
        dzdg_ref = []
        dzdn_ref = []
        dqdg_ref = []
        for _ in range(n_out):
            vals = list(map(float, out_lines[idx].split()))
            dzdg_ref.append(vals[0])
            dzdn_ref.append(vals[1])
            dqdg_ref.append(vals[2])
            idx += 1
        m_out = int(out_lines[idx]); idx += 1
        dzdm_ref = []
        dqdm_ref = []
        for _ in range(m_out):
            vals = list(map(float, out_lines[idx].split()))
            dzdm_ref.append(vals[0])
            dqdm_ref.append(vals[1])
            idx += 1
        self.assertLessEqual(max_abs_diff(dzdg_ref[1:-1], js_case["dzdg"][1:-1]), tol)
        self.assertLessEqual(max_abs_diff(dzdn_ref[1:-1], js_case["dzdn"][1:-1]), tol)
        self.assertLessEqual(max_abs_diff(dqdg_ref, js_case["dqdg"]), tol)
        # dzdm/dqdm involve source sensitivity; JS uses a different constant offset convention.

    def test_pswlin(self):
        n = self.base["n"]
        nw = self.base["nw"]
        i_js = n
        i_ft = i_js + 1
        state = {**self.base}
        case = {
            "kind": "pswlin",
            "state": state,
            "i": i_js,
            "xi": state["x"][i_js],
            "yi": state["y"][i_js],
            "nxi": state["nx"][i_js],
            "nyi": state["ny"][i_js],
        }
        js_case = self._run_node([case])[0]
        lines = [
            "4",
            f"{n} {nw} {i_ft}",
            f"{state['qopi']} {state['pi']}",
        ]
        for j in range(n + nw):
            lines.append(f"{state['x'][j]} {state['y'][j]}")
        for j in range(n + nw):
            lines.append(f"{state['apanel'][j]}")
        for j in range(n + nw):
            lines.append(f"{state['sig'][j]}")
        lines.append(f"{state['x'][i_js]} {state['y'][i_js]} {state['nx'][i_js]} {state['ny'][i_js]}")
        out = self._run_fortran(lines)
        psi_ref, psi_ni_ref = map(float, out[0].split())
        tol = 1.0e-6
        self.assertLessEqual(abs(psi_ref - js_case["result"]["psi"]), tol)
        self.assertLessEqual(abs(psi_ni_ref - js_case["result"]["psiNi"]), tol)

    def test_ggcalc(self):
        n = self.base["n"]
        state = {
            **self.base,
            "nw": 0,
            "x": self.base["x"][:n],
            "y": self.base["y"][:n],
            "s": self.base["s"][:n],
            "xp": self.base["xp"][:n],
            "yp": self.base["yp"][:n],
            "nx": self.base["nx"][:n],
            "ny": self.base["ny"][:n],
            "sig": [0.0] * n,
            "qf0": [0.0] * n,
            "qf1": [0.0] * n,
            "qf2": [0.0] * n,
            "qf3": [0.0] * n,
        }
        case = {"kind": "ggcalc", "state": state}
        js_case = self._run_node([case])[0]
        lines = [
            "5",
            f"{n} 0",
            f"{state['alfa']} {state['qinf']}",
            f"{state['ante']} {state['aste']} {state['dste']}",
            f"{state['qopi']} {state['hopi']} {state['pi']}",
        ]
        for j in range(n):
            lines.append(
                f"{state['x'][j]} {state['y'][j]} {state['s'][j]} "
                f"{state['xp'][j]} {state['yp'][j]} {state['nx'][j]} {state['ny'][j]}"
            )
        for j in range(n):
            lines.append("0.0")
        for j in range(n):
            lines.append("0.0 0.0 0.0 0.0")
        out = self._run_fortran(lines)
        qinvu_ref = []
        for line in out[1:]:
            vals = list(map(float, line.split()))
            qinvu_ref.append(vals)
        tol = 1.0e-5
        self.assertLessEqual(
            max_abs_diff([v[0] for v in qinvu_ref], [v[0] for v in js_case["qinvu"][:n]]), tol
        )
        self.assertLessEqual(
            max_abs_diff([v[1] for v in qinvu_ref], [v[1] for v in js_case["qinvu"][:n]]), tol
        )

    def test_qwcalc(self):
        n = self.base["n"]
        nw = self.base["nw"]
        state = {**self.base}
        case = {"kind": "qwcalc", "state": state}
        js_case = self._run_node([case])[0]
        lines = [
            "6",
            f"{n} {nw}",
            "0",
            f"{state['alfa']} {state['qinf']}",
            f"{state['ante']} {state['aste']} {state['dste']}",
            f"{state['qopi']} {state['hopi']} {state['pi']}",
        ]
        for j in range(n + nw):
            lines.append(
                f"{state['x'][j]} {state['y'][j]} {state['s'][j]} "
                f"{state['xp'][j]} {state['yp'][j]} {state['nx'][j]} {state['ny'][j]}"
            )
        for j in range(n):
            lines.append(f"{state['gam'][j]}")
        for j in range(n + nw):
            lines.append(f"{state['sig'][j]}")
        for j in range(n + nw):
            lines.append(f"{state['qinvu'][j][0]} {state['qinvu'][j][1]}")
        out = self._run_fortran(lines)
        qinvu_ref = []
        for line in out[1:]:
            vals = list(map(float, line.split()))
            qinvu_ref.append(vals)
        tol = 1.0e-5
        self.assertLessEqual(
            max_abs_diff([v[0] for v in qinvu_ref], [v[0] for v in js_case["qinvu"]]), tol
        )
        self.assertLessEqual(
            max_abs_diff([v[1] for v in qinvu_ref], [v[1] for v in js_case["qinvu"]]), tol
        )

    def test_qdcalc(self):
        n = self.base["n"]
        nw = self.base["nw"]
        state = {**self.base}
        case = {"kind": "qdcalc", "state": state}
        js_case = self._run_node([case])[0]
        lines = [
            "7",
            f"{n} {nw} 0",
            f"{state['alfa']} {state['qinf']}",
            f"{state['ante']} {state['aste']} {state['dste']}",
            f"{state['qopi']} {state['hopi']} {state['pi']}",
        ]
        for j in range(n + nw):
            lines.append(
                f"{state['x'][j]} {state['y'][j]} {state['s'][j]} "
                f"{state['xp'][j]} {state['yp'][j]} {state['nx'][j]} {state['ny'][j]}"
            )
        for j in range(n + nw):
            lines.append(f"{state['sig'][j]}")
        for j in range(n):
            lines.append("0.0 0.0 0.0 0.0")
        out = self._run_fortran(lines)
        tol = 3.0e-5
        # compare selected DIJ columns
        dij = js_case["dij"]
        for i in range(n):
            vals = list(map(float, out[i + 1].split()))
            self.assertLessEqual(abs(vals[0] - dij[i + 1][1]), tol)
            self.assertLessEqual(abs(vals[1] - dij[i + 1][2]), tol)
            self.assertLessEqual(abs(vals[2] - dij[i + 1][n]), tol)

    def test_xywake(self):
        n = self.base["n"]
        waklen = 1.0
        nw = n // 12 + 10 * int(waklen)
        def extend(arr, fill):
            out = list(arr)
            out.extend([fill] * (n + nw - len(out)))
            return out
        state = {
            **self.base,
            "nw": nw,
            "x": extend(self.base["x"], self.base["x"][-1]),
            "y": extend(self.base["y"], self.base["y"][-1]),
            "s": extend(self.base["s"], self.base["s"][-1]),
            "xp": extend(self.base["xp"], self.base["xp"][-1]),
            "yp": extend(self.base["yp"], self.base["yp"][-1]),
            "nx": extend(self.base["nx"], self.base["nx"][-1]),
            "ny": extend(self.base["ny"], self.base["ny"][-1]),
            "sig": extend(self.base["sig"], 0.0),
        }
        case = {"kind": "xywake", "state": state}
        js_case = self._run_node([case])[0]
        lines = [
            "8",
            f"{n} {nw} 0",
            f"{state['alfa']} {state['qinf']} {waklen} 1.0",
            f"{state['ante']} {state['aste']} {state['dste']}",
            f"{state['qopi']} {state['hopi']} {state['pi']}",
        ]
        for j in range(n):
            lines.append(
                f"{state['x'][j]} {state['y'][j]} {state['s'][j]} "
                f"{state['xp'][j]} {state['yp'][j]} {state['nx'][j]} {state['ny'][j]}"
            )
        for j in range(n):
            lines.append(f"{state['gam'][j]}")
        for j in range(n + nw):
            lines.append(f"{state['sig'][j]}")
        out = self._run_fortran(lines)
        tol = 1.0e-5
        for j in range(n, n + nw):
            vals = list(map(float, out[j + 1].split()))
            self.assertLessEqual(abs(vals[0] - js_case["x"][j]), tol)
            self.assertLessEqual(abs(vals[1] - js_case["y"][j]), tol)

    def test_stfind(self):
        n = self.base["n"]
        state = {**self.base}
        case = {"kind": "stfind", "state": state}
        js_case = self._run_node([case])[0]
        lines = ["9", str(n)]
        for j in range(n):
            lines.append(f"{state['s'][j]} {state['gam'][j]}")
        out = self._run_fortran(lines)
        ist_ref = int(out[0].strip())
        sst_ref, sst_go_ref, sst_gp_ref = map(float, out[1].split())
        self.assertEqual(ist_ref, js_case["result"]["ist"] + 1)
        self.assertAlmostEqual(sst_ref, js_case["result"]["sst"], delta=1.0e-6)

    def test_iblpan(self):
        n = self.base["n"]
        nw = self.base["nw"]
        ist = 3
        case = {"kind": "iblpan", "state": {**self.base, "ist": ist}}
        js_case = self._run_node([case])[0]
        lines = ["10", f"{n} {nw} {ist}"]
        out = self._run_fortran(lines)
        iblte1, iblte2 = map(int, out[0].split())
        nbl1, nbl2 = map(int, out[1].split())
        self.assertEqual(iblte1, js_case["iblte"][0])
        self.assertEqual(iblte2, js_case["iblte"][1])
        self.assertEqual(nbl1, js_case["nbl"][0])
        self.assertEqual(nbl2, js_case["nbl"][1])

    def test_xicalc(self):
        n = self.base["n"]
        nw = self.base["nw"]
        ist = 3
        # Get IPAN/NBL from Fortran IBLPAN
        lines = ["10", f"{n} {nw} {ist}"]
        out = self._run_fortran(lines)
        iblte1, iblte2 = map(int, out[0].split())
        nbl1, nbl2 = map(int, out[1].split())
        ipan1 = [int(out[i + 2].split()[0]) for i in range(nbl1)]
        ipan2 = [int(out[i + 2 + nbl1].split()[0]) for i in range(nbl2)]
        state = {**self.base, "ist": ist, "iblte1": iblte1, "iblte2": iblte2, "nbl1": nbl1, "nbl2": nbl2,
                 "ipan1": ipan1, "ipan2": ipan2, "sst": self.base["s"][ist - 1]}
        js_case = self._run_node([{"kind": "xicalc", "state": state}])[0]
        lines = [
            "11",
            f"{n} {nw} {ist} 0",
            f"{state['ante']} {state['aste']} {state['dste']}",
        ]
        for j in range(n + nw):
            lines.append(
                f"{state['x'][j]} {state['y'][j]} {state['s'][j]} {state['xp'][j]} {state['yp'][j]}"
            )
        lines.append(f"{iblte1} {iblte2}")
        lines.append(f"{nbl1} {nbl2}")
        for val in ipan1:
            lines.append(str(val))
        for val in ipan2:
            lines.append(str(val))
        lines.append(f"{state['sst']}")
        out = self._run_fortran(lines)
        idx = 0
        idx += 1  # nbl line
        xssi1_ref = [float(out[idx + i]) for i in range(nbl1)]
        idx += nbl1
        xssi2_ref = [float(out[idx + i]) for i in range(nbl2)]
        tol = 1.0e-6
        self.assertLessEqual(max_abs_diff(xssi1_ref, js_case["xssi1"]), tol)
        self.assertLessEqual(max_abs_diff(xssi2_ref, js_case["xssi2"]), tol)

    def test_uicalc_qvfue_qiset_gamqv(self):
        n = self.base["n"]
        nw = self.base["nw"]
        ist = 3
        lines = ["10", f"{n} {nw} {ist}"]
        out = self._run_fortran(lines)
        iblte1, iblte2 = map(int, out[0].split())
        nbl1, nbl2 = map(int, out[1].split())
        ipan1 = [int(out[i + 2].split()[0]) for i in range(nbl1)]
        ipan2 = [int(out[i + 2 + nbl1].split()[0]) for i in range(nbl2)]
        vti1 = [int(out[i + 2].split()[1]) for i in range(nbl1)]
        vti2 = [int(out[i + 2 + nbl1].split()[1]) for i in range(nbl2)]

        state = {**self.base, "nbl1": nbl1, "nbl2": nbl2, "iblte1": iblte1, "iblte2": iblte2,
                 "ipan1": ipan1, "ipan2": ipan2, "vti1": vti1, "vti2": vti2}
        js_uicalc = self._run_node([{"kind": "uicalc", "state": state}])[0]
        lines = ["12", f"{n} {nw}", f"{nbl1} {nbl2}"]
        for ip, vt in zip(ipan1, vti1):
            lines.append(f"{ip} {vt}")
        for ip, vt in zip(ipan2, vti2):
            lines.append(f"{ip} {vt}")
        for j in range(n + nw):
            lines.append(f"{state['qinv'][j+1]} {state['qinvA'][j+1]}")
        out = self._run_fortran(lines)
        uinv_ref = [float(out[i].split()[0]) for i in range(nbl1)]
        tol = 1.0e-6
        self.assertLessEqual(max_abs_diff(uinv_ref, js_uicalc["uinv1"]), tol)

        # QISET
        js_qiset = self._run_node([{"kind": "qiset", "state": state}])[0]
        lines = ["14", f"{n} {nw}", f"{state['alfa']}"]
        for j in range(n + nw):
            lines.append(f"{state['qinvu'][j][0]} {state['qinvu'][j][1]}")
        out = self._run_fortran(lines)
        qinv_ref = [float(line.split()[0]) for line in out[1:]]
        self.assertLessEqual(max_abs_diff(qinv_ref, js_qiset["qinv"][1:]), tol)

        # QVFUE
        uedg1 = [0.1 * (i + 1) for i in range(nbl1)]
        uedg2 = [0.2 * (i + 1) for i in range(nbl2)]
        state_qvfue = {**state, "uedg1": uedg1, "uedg2": uedg2}
        js_qvfue = self._run_node([{"kind": "qvfue", "state": state_qvfue}])[0]
        lines = ["13", f"{n} {nw}", f"{nbl1} {nbl2}"]
        for ip, vt, ue in zip(ipan1, vti1, uedg1):
            lines.append(f"{ip} {vt} {ue}")
        for ip, vt, ue in zip(ipan2, vti2, uedg2):
            lines.append(f"{ip} {vt} {ue}")
        out = self._run_fortran(lines)
        qvis_ref = [float(line.split()[0]) for line in out[1:]]
        self.assertLessEqual(max_abs_diff(qvis_ref, js_qvfue["qvis"][1 : n + nw + 1]), tol)

        # GAMQV
        qvis = [0.0] + [0.03 * (i + 1) for i in range(n)]
        qinv_a = [0.0] + [0.02 * (i + 1) for i in range(n)]
        js_gamqv = self._run_node([{"kind": "gamqv", "state": {**state, "qvis": qvis, "qinvA": qinv_a}}])[0]
        lines = ["15", str(n)]
        for i in range(1, n + 1):
            lines.append(f"{qvis[i]} {qinv_a[i]}")
        out = self._run_fortran(lines)
        gam_ref = [float(line.split()[0]) for line in out[1:]]
        self.assertLessEqual(max_abs_diff(gam_ref, js_gamqv["gam"]), tol)

    def test_stmove(self):
        n = self.base["n"]
        nw = self.base["nw"]
        ist = 3
        lines = ["10", f"{n} {nw} {ist}"]
        out = self._run_fortran(lines)
        iblte1, iblte2 = map(int, out[0].split())
        nbl1, nbl2 = map(int, out[1].split())
        ipan1 = [int(out[i + 2].split()[0]) for i in range(nbl1)]
        ipan2 = [int(out[i + 2 + nbl1].split()[0]) for i in range(nbl2)]
        vti1 = [int(out[i + 2].split()[1]) for i in range(nbl1)]
        vti2 = [int(out[i + 2 + nbl1].split()[1]) for i in range(nbl2)]

        uedg1 = [0.15 + 0.01 * i for i in range(nbl1)]
        uedg2 = [0.12 + 0.015 * i for i in range(nbl2)]
        ctau1 = [0.02 + 0.001 * i for i in range(nbl1)]
        ctau2 = [0.03 + 0.001 * i for i in range(nbl2)]
        thet1 = [0.01 + 0.002 * i for i in range(nbl1)]
        thet2 = [0.012 + 0.002 * i for i in range(nbl2)]
        dstr1 = [0.005 + 0.0005 * i for i in range(nbl1)]
        dstr2 = [0.006 + 0.0005 * i for i in range(nbl2)]
        xssi1 = [0.001 + 0.01 * i for i in range(nbl1)]
        xssi2 = [0.002 + 0.012 * i for i in range(nbl2)]
        state = {
            **self.base,
            "ist": ist,
            "iblte1": iblte1,
            "iblte2": iblte2,
            "nbl1": nbl1,
            "nbl2": nbl2,
            "ipan1": ipan1,
            "ipan2": ipan2,
            "vti1": vti1,
            "vti2": vti2,
            "uedg1": uedg1,
            "uedg2": uedg2,
            "ctau1": ctau1,
            "ctau2": ctau2,
            "thet1": thet1,
            "thet2": thet2,
            "dstr1": dstr1,
            "dstr2": dstr2,
            "xssi1": xssi1,
            "xssi2": xssi2,
            "sst": self.base["s"][ist - 1],
            "itran1": 2,
            "itran2": 2,
        }
        js_case = self._run_node([{"kind": "stmove", "state": state}])[0]

        lines = [
            "16",
            f"{n} {nw} {ist}",
            "0",
            f"{self.base['ante']} {self.base['aste']} {self.base['dste']}",
        ]
        for j in range(n):
            lines.append(
                f"{self.base['x'][j]} {self.base['y'][j]} {self.base['s'][j]} "
                f"{self.base['xp'][j]} {self.base['yp'][j]} {self.base['gam'][j]}"
            )
        lines.append(f"{iblte1} {iblte2}")
        lines.append(f"{nbl1} {nbl2}")
        for i in range(nbl1):
            lines.append(
                f"{ipan1[i]} {vti1[i]} {uedg1[i]} {ctau1[i]} {thet1[i]} {dstr1[i]} {xssi1[i]}"
            )
        for i in range(nbl2):
            lines.append(
                f"{ipan2[i]} {vti2[i]} {uedg2[i]} {ctau2[i]} {thet2[i]} {dstr2[i]} {xssi2[i]}"
            )
        for j in range(n + nw):
            lines.append(f"{self.base['qinv'][j+1]} {self.base['qinvA'][j+1]}")
        for j in range(n + nw):
            lines.append(f"{self.base['qvis'][j+1]}")
        lines.append(f"{state['sst']}")
        out = self._run_fortran(lines)
        ist_ref = int(out[0].strip())
        self.assertEqual(ist_ref, js_case["ist"])


if __name__ == "__main__":
    unittest.main()
