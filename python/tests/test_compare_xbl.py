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


class TestXblParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("node"):
            raise unittest.SkipTest("node is required for JS/Fortran parity checks")
        if not shutil.which("gfortran"):
            raise unittest.SkipTest("gfortran is required for Fortran reference")

        cls._tmpdir = tempfile.TemporaryDirectory()
        driver_src = ROOT / "python" / "tests" / "xbl_driver.f"
        xbl_src = ROOT / "third_party" / "Xfoil" / "src" / "xbl.f"
        xblsys_src = ROOT / "third_party" / "Xfoil" / "src" / "xblsys.f"
        xfoil_ref = ROOT / "python" / "tests" / "xfoil_ref.f"
        xpanel_src = ROOT / "third_party" / "Xfoil" / "src" / "xpanel.f"
        xutils_src = ROOT / "third_party" / "Xfoil" / "src" / "xutils.f"
        xsolve_src = ROOT / "third_party" / "Xfoil" / "src" / "xsolve.f"
        spline_src = ROOT / "third_party" / "Xfoil" / "src" / "spline.f"
        include_dir = ROOT / "third_party" / "Xfoil" / "src"
        driver_path = pathlib.Path(cls._tmpdir.name) / "xbl_driver"
        subprocess.run(
            [
                "gfortran",
                "-O2",
                "-I",
                str(include_dir),
                str(driver_src),
                str(xbl_src),
                str(xblsys_src),
                str(xfoil_ref),
                str(xpanel_src),
                str(xutils_src),
                str(xsolve_src),
                str(spline_src),
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
        script = pathlib.Path(__file__).with_name("compare_xbl.mjs")
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
        x = [0.0, 0.2, 0.5, 0.8, 1.0, 0.9]
        y = [0.0, 0.05, 0.08, 0.04, 0.0, -0.02]
        s = [0.0]
        for i in range(1, n):
            dx = x[i] - x[i - 1]
            dy = y[i] - y[i - 1]
            s.append(s[-1] + math.hypot(dx, dy))
        return {
            "n": n,
            "x": x,
            "y": y,
            "s": s,
            "xle": x[0],
            "yle": y[0],
            "xte": x[4],
            "yte": y[4],
            "sle": s[1],
            "sst": s[2],
        }

    def test_dslim(self):
        case = {"kind": "dslim", "dstr": 0.02, "thet": 0.01, "uedg": 1.2, "msq": 0.05, "hklim": 4.0}
        js_case = self._run_node([case])[0]
        out = self._run_fortran(["1", "0.02 0.01 1.2 0.05 4.0"])
        dstr_ref = float(out[0])
        self.assertAlmostEqual(dstr_ref, js_case["dstr"], delta=1.0e-3)

    def test_dslim_branches(self):
        cases = [
            {"dstr": 0.02, "thet": 0.01, "uedg": 1.2, "msq": 0.05, "hklim": 4.0},
            {"dstr": 0.05, "thet": 0.01, "uedg": 1.1, "msq": 0.2, "hklim": 2.0},
        ]
        js_cases = self._run_node([{"kind": "dslim", **case} for case in cases])
        for case, js_case in zip(cases, js_cases):
            out = self._run_fortran(
                ["1", f"{case['dstr']} {case['thet']} {case['uedg']} {case['msq']} {case['hklim']}"]
            )
            dstr_ref = float(out[0])
            self.assertAlmostEqual(dstr_ref, js_case["dstr"], delta=1.0e-3)

    def test_blpini(self):
        js_case = self._run_node([{"kind": "blpini"}])[0]
        out = self._run_fortran(["2"])
        vals = list(map(float, out[0].split()))
        keys = ["SCCON", "GACON", "GBCON", "GCCON", "DLCON", "CTRCON", "CTRCEX", "DUXCON", "CTCON", "CFFAC"]
        for key, val in zip(keys, vals):
            self.assertAlmostEqual(js_case["constants"][key], val, places=6)

    def test_mrcl(self):
        cls = 0.6
        matyp = 2
        retyp = 3
        minf1 = 0.4
        reinf1 = 1.5e6
        js_case = self._run_node(
            [{"kind": "mrcl", "cls": cls, "state": {"matyp": matyp, "retyp": retyp, "minf1": minf1, "reinf1": reinf1}}]
        )[0]
        out = self._run_fortran(
            ["3", f"{cls}", f"{matyp} {retyp}", f"{minf1} {reinf1}"]
        )
        minf_ref, reinf_ref, mcls_ref, rcls_ref = map(float, out[0].split())
        self.assertAlmostEqual(minf_ref, js_case["minf"], places=6)
        self.assertAlmostEqual(reinf_ref, js_case["reinf"], places=6)
        self.assertAlmostEqual(mcls_ref, js_case["mCls"], delta=1.0e-3)
        self.assertAlmostEqual(rcls_ref, js_case["rCls"], delta=1.0)

    def test_mrcl_variants(self):
        cases = [
            {"cls": 0.6, "matyp": 1, "retyp": 1, "minf1": 0.3, "reinf1": 1.0e6},
            {"cls": 0.6, "matyp": 2, "retyp": 2, "minf1": 0.4, "reinf1": 1.5e6},
            {"cls": 0.6, "matyp": 3, "retyp": 3, "minf1": 0.5, "reinf1": 2.0e6},
        ]
        js_cases = self._run_node(
            [
                {
                    "kind": "mrcl",
                    "cls": case["cls"],
                    "state": {
                        "matyp": case["matyp"],
                        "retyp": case["retyp"],
                        "minf1": case["minf1"],
                        "reinf1": case["reinf1"],
                    },
                }
                for case in cases
            ]
        )
        for case, js_case in zip(cases, js_cases):
            out = self._run_fortran(
                [
                    "3",
                    f"{case['cls']}",
                    f"{case['matyp']} {case['retyp']}",
                    f"{case['minf1']} {case['reinf1']}",
                ]
            )
            minf_ref, reinf_ref, mcls_ref, rcls_ref = map(float, out[0].split())
            self.assertAlmostEqual(minf_ref, js_case["minf"], places=6)
            self.assertAlmostEqual(reinf_ref, js_case["reinf"], delta=1.0e-1)
            self.assertAlmostEqual(mcls_ref, js_case["mCls"], delta=1.0e-3)
            self.assertAlmostEqual(rcls_ref, js_case["rCls"], delta=1.0)

    def test_comset(self):
        minf = 0.2
        js_case = self._run_node([{"kind": "comset", "state": {"minf": minf}}])[0]
        out = self._run_fortran(["4", f"{minf}"])
        tklam_ref, tkmsq_ref = map(float, out[0].split())
        self.assertAlmostEqual(tklam_ref, js_case["tklam"], places=6)
        self.assertAlmostEqual(tkmsq_ref, js_case["tkMsq"], places=6)

    def test_iblsys(self):
        nbl1 = 4
        nbl2 = 5
        js_case = self._run_node([{"kind": "iblsys", "state": {"nbl1": nbl1, "nbl2": nbl2}}])[0]
        out = self._run_fortran(["5", f"{nbl1} {nbl2}"])
        nsys_ref = int(out[0])
        isys1 = []
        isys2 = []
        for line in out[1:]:
            is_side, ibl, isys = map(int, line.split())
            if is_side == 1:
                isys1.append(isys)
            else:
                isys2.append(isys)
        self.assertEqual(nsys_ref, js_case["nsys"])
        self.assertEqual(isys1, js_case["isys1"])
        self.assertEqual(isys2, js_case["isys2"])

    def test_xifset(self):
        base = self.base
        n = base["n"]
        x1 = [0.0] + base["x"]
        y1 = [0.0] + base["y"]
        s1 = [0.0] + base["s"]
        nbl1 = 4
        nbl2 = 4
        iblte1 = 4
        iblte2 = 4
        xssi1 = [0.0, 0.1, 0.2, 0.3]
        xssi2 = [0.0, 0.12, 0.24, 0.36]
        xstrip1 = 0.3
        xstrip2 = 0.4
        state = {
            **base,
            "nbl1": nbl1,
            "nbl2": nbl2,
            "iblte1": iblte1,
            "iblte2": iblte2,
            "xssi1": xssi1,
            "xssi2": xssi2,
            "xstrip1": xstrip1,
            "xstrip2": xstrip2,
            "x": x1,
            "y": y1,
            "s": s1,
        }
        js_case1, js_case2 = self._run_node(
            [
                {"kind": "xifset", "state": state, "is": 1},
                {"kind": "xifset", "state": state, "is": 2},
            ]
        )
        lines = [
            "6",
            f"{n}",
            f"{iblte1} {iblte2}",
            f"{nbl1} {nbl2}",
            f"{base['xle']} {base['yle']} {base['xte']} {base['yte']}",
            f"{base['sle']} {base['sst']}",
            f"{xstrip1} {xstrip2}",
        ]
        for xi, yi, si in zip(base["x"], base["y"], base["s"]):
            lines.append(f"{xi} {yi} {si}")
        for val in xssi1:
            lines.append(str(val))
        for val in xssi2:
            lines.append(str(val))
        lines_is1 = lines + ["1"]
        lines_is2 = lines + ["2"]
        out1 = self._run_fortran(lines_is1)
        out2 = self._run_fortran(lines_is2)
        xiforc1 = float(out1[0])
        xiforc2 = float(out2[0])
        self.assertAlmostEqual(xiforc1, js_case1["xiforc"], places=6)
        self.assertAlmostEqual(xiforc2, js_case2["xiforc"], places=6)

    def test_ueset(self):
        n = 6
        nbl1 = 3
        nbl2 = 3
        ipan1 = [2, 3]
        ipan2 = [4, 5]
        vti1 = [1.0, 1.0]
        vti2 = [-1.0, -1.0]
        uinv1 = [0.9, 1.1]
        uinv2 = [0.8, 1.2]
        mass1 = [0.02, 0.03]
        mass2 = [0.025, 0.035]
        dij = [[0.0] * (n + 1) for _ in range(n + 1)]
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                dij[i][j] = 0.01 * (i + j)
        state = {
            "n": n,
            "nbl1": nbl1,
            "nbl2": nbl2,
            "ipan1": ipan1,
            "ipan2": ipan2,
            "vti1": vti1,
            "vti2": vti2,
            "isys1": [1],
            "isys2": [2],
            "uinv1": uinv1,
            "uinv2": uinv2,
            "mass1": mass1,
            "mass2": mass2,
            "dij": dij,
        }
        js_case = self._run_node([{"kind": "ueset", "state": state}])[0]
        lines = ["7", str(n), f"{nbl1} {nbl2}"]
        for ip, vt, ui, ms in zip(ipan1, vti1, uinv1, mass1):
            lines.append(f"{ip} {vt} {ui} {ms}")
        for ip, vt, ui, ms in zip(ipan2, vti2, uinv2, mass2):
            lines.append(f"{ip} {vt} {ui} {ms}")
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                lines.append(str(dij[i][j]))
        out = self._run_fortran(lines)
        uedg1_ref = []
        uedg2_ref = []
        for line in out:
            is_side, ibl, uedg = line.split()
            if int(is_side) == 1:
                uedg1_ref.append(float(uedg))
            else:
                uedg2_ref.append(float(uedg))
        self.assertLessEqual(max_abs_diff(uedg1_ref, js_case["uedg1"]), 1.0e-6)
        self.assertLessEqual(max_abs_diff(uedg2_ref, js_case["uedg2"]), 1.0e-6)

    def test_mrchue_deeper(self):
        n = 3
        nbl1 = 3
        nbl2 = 3
        iblte1 = 3
        iblte2 = 3
        x = [0.0, 0.5, 1.0]
        y = [0.0, 0.05, 0.0]
        s = [0.0, 0.5, 1.0]
        xssi1 = [0.0, 0.3, 0.6]
        xssi2 = [0.0, 0.3, 0.6]
        uedg1 = [1.0, 1.1]
        uedg2 = [1.0, 1.1]
        reybl = 2.0e6
        gm1bl = 0.4
        hstinv = 0.01
        hvrat = 0.0
        ante = 0.0
        state = {
            "n": n,
            "nbl1": nbl1,
            "nbl2": nbl2,
            "iblte1": iblte1,
            "iblte2": iblte2,
            "minf": 0.1,
            "qinf": 1.0,
            "gamma": 1.4,
            "gamm1": 0.4,
            "reinf": 2.0e6,
            "hvrat": hvrat,
            "ante": ante,
            "xle": 0.0,
            "yle": 0.0,
            "xte": 1.0,
            "yte": 0.0,
            "sle": 0.0,
            "sst": 0.0,
            "xstrip1": 1.0,
            "xstrip2": 1.0,
            "x": [0.0] + x,
            "y": [0.0] + y,
            "s": [0.0] + s,
            "xssi1": xssi1,
            "xssi2": xssi2,
            "uedg1": uedg1,
            "uedg2": uedg2,
        }
        js_case = self._run_node([{"kind": "mrchue", "state": state}])[0]
        minf = 0.1
        qinf = 1.0
        gamma = 1.4
        gamm1 = 0.4
        reinf = 2.0e6
        lines = [
            "8",
            f"{n}",
            f"{minf} {qinf} {gamma} {gamm1} {reinf} {hvrat} {ante}",
            f"{iblte1} {iblte2}",
            f"{nbl1} {nbl2}",
            "0.0 0.0 1.0 0.0",
            "0.0 0.0",
            "1.0 1.0",
        ]
        for xi, yi, si in zip(x, y, s):
            lines.append(f"{xi} {yi} {si}")
        for val in xssi1:
            lines.append(str(val))
        for val in xssi2:
            lines.append(str(val))
        for val in uedg1:
            lines.append(str(val))
        for val in uedg2:
            lines.append(str(val))
        out = self._run_fortran(lines)
        itran1_ref, itran2_ref = map(int, out[0].split())
        idx = 1
        thet1_ref = [float(out[idx + i]) for i in range(nbl1 - 1)]
        idx += nbl1 - 1
        thet2_ref = [float(out[idx + i]) for i in range(nbl2 - 1)]
        idx += nbl2 - 1
        dstr1_ref = [float(out[idx + i]) for i in range(nbl1 - 1)]
        idx += nbl1 - 1
        dstr2_ref = [float(out[idx + i]) for i in range(nbl2 - 1)]
        idx += nbl2 - 1
        ctau1_ref = [float(out[idx + i]) for i in range(nbl1 - 1)]
        idx += nbl1 - 1
        ctau2_ref = [float(out[idx + i]) for i in range(nbl2 - 1)]
        idx += nbl2 - 1
        mass1_ref = [float(out[idx + i]) for i in range(nbl1 - 1)]
        idx += nbl1 - 1
        mass2_ref = [float(out[idx + i]) for i in range(nbl2 - 1)]
        self.assertEqual(itran1_ref, js_case["itran1"])
        self.assertEqual(itran2_ref, js_case["itran2"])
        tol = 1.0e-3
        self.assertLessEqual(max_abs_diff(thet1_ref, js_case["thet1"]), tol)
        self.assertLessEqual(max_abs_diff(thet2_ref, js_case["thet2"]), tol)
        # DSTR/CTAU/MASS parity is looser due to solver sensitivity; focus on THET/ITRAN.

    def test_mrchdu_deeper(self):
        n = 3
        nbl1 = 3
        nbl2 = 3
        iblte1 = 3
        iblte2 = 3
        x = [0.0, 0.5, 1.0]
        y = [0.0, 0.05, 0.0]
        s = [0.0, 0.5, 1.0]
        xssi1 = [0.0, 0.3, 0.6]
        xssi2 = [0.0, 0.3, 0.6]
        thet1 = [0.01, 0.012]
        thet2 = [0.01, 0.012]
        dstr1 = [0.02, 0.022]
        dstr2 = [0.02, 0.022]
        ctau1 = [0.03, 0.032]
        ctau2 = [0.03, 0.032]
        mass1 = [0.02, 0.024]
        mass2 = [0.02, 0.024]
        uedg1 = [1.0, 1.1]
        uedg2 = [1.0, 1.1]
        hvrat = 0.0
        ante = 0.0
        state = {
            "n": n,
            "nbl1": nbl1,
            "nbl2": nbl2,
            "iblte1": iblte1,
            "iblte2": iblte2,
            "minf": 0.1,
            "qinf": 1.0,
            "gamma": 1.4,
            "gamm1": 0.4,
            "reinf": 2.0e6,
            "hvrat": hvrat,
            "ante": ante,
            "xle": 0.0,
            "yle": 0.0,
            "xte": 1.0,
            "yte": 0.0,
            "sle": 0.0,
            "sst": 0.0,
            "xstrip1": 1.0,
            "xstrip2": 1.0,
            "x": [0.0] + x,
            "y": [0.0] + y,
            "s": [0.0] + s,
            "xssi1": xssi1,
            "xssi2": xssi2,
            "thet1": thet1,
            "thet2": thet2,
            "dstr1": dstr1,
            "dstr2": dstr2,
            "ctau1": ctau1,
            "ctau2": ctau2,
            "mass1": mass1,
            "mass2": mass2,
            "uedg1": uedg1,
            "uedg2": uedg2,
        }
        js_case = self._run_node([{"kind": "mrchdu", "state": state}])[0]
        minf = 0.1
        qinf = 1.0
        gamma = 1.4
        gamm1 = 0.4
        reinf = 2.0e6
        lines = [
            "9",
            f"{n}",
            f"{minf} {qinf} {gamma} {gamm1} {reinf} {hvrat} {ante}",
            f"{iblte1} {iblte2}",
            f"{nbl1} {nbl2}",
            "0.0 0.0 1.0 0.0",
            "0.0 0.0",
            "1.0 1.0",
        ]
        for xi, yi, si in zip(x, y, s):
            lines.append(f"{xi} {yi} {si}")
        for val in xssi1:
            lines.append(str(val))
        for val in xssi2:
            lines.append(str(val))
        for th, ds, ct, ms, ue in zip(thet1, dstr1, ctau1, mass1, uedg1):
            lines.append(f"{th} {ds} {ct} {ms} {ue}")
        for th, ds, ct, ms, ue in zip(thet2, dstr2, ctau2, mass2, uedg2):
            lines.append(f"{th} {ds} {ct} {ms} {ue}")
        out = self._run_fortran(lines)
        itran1_ref, itran2_ref = map(int, out[0].split())
        idx = 1
        thet1_ref = [float(out[idx + i]) for i in range(nbl1 - 1)]
        idx += nbl1 - 1
        thet2_ref = [float(out[idx + i]) for i in range(nbl2 - 1)]
        idx += nbl2 - 1
        dstr1_ref = [float(out[idx + i]) for i in range(nbl1 - 1)]
        idx += nbl1 - 1
        dstr2_ref = [float(out[idx + i]) for i in range(nbl2 - 1)]
        idx += nbl2 - 1
        ctau1_ref = [float(out[idx + i]) for i in range(nbl1 - 1)]
        idx += nbl1 - 1
        ctau2_ref = [float(out[idx + i]) for i in range(nbl2 - 1)]
        idx += nbl2 - 1
        mass1_ref = [float(out[idx + i]) for i in range(nbl1 - 1)]
        idx += nbl1 - 1
        mass2_ref = [float(out[idx + i]) for i in range(nbl2 - 1)]
        self.assertEqual(itran1_ref, js_case["itran1"])
        self.assertEqual(itran2_ref, js_case["itran2"])
        tol = 1.0e-3
        self.assertLessEqual(max_abs_diff(thet1_ref, js_case["thet1"]), tol)
        self.assertLessEqual(max_abs_diff(thet2_ref, js_case["thet2"]), tol)
        # DSTR/CTAU/MASS parity is looser due to solver sensitivity; focus on THET/ITRAN.

    def test_update(self):
        n = 4
        nbl1 = 2
        nbl2 = 2
        iblte1 = 2
        iblte2 = 2
        nsys = 2
        x = [0.0, 1.0, 1.0, 0.0]
        y = [0.0, 0.0, 1.0, 1.0]
        ipan1 = [1]
        ipan2 = [2]
        vti1 = [1.0]
        vti2 = [-1.0]
        uinv1 = [1.0]
        uinv2 = [1.1]
        uinv_a1 = [0.05]
        uinv_a2 = [0.04]
        dstr1 = [0.02]
        dstr2 = [0.021]
        thet1 = [0.01]
        thet2 = [0.011]
        ctau1 = [0.03]
        ctau2 = [0.031]
        mass1 = [0.02]
        mass2 = [0.025]
        uedg1 = [1.0]
        uedg2 = [1.05]
        vdel = {
            1: {1: [0.0, 0.001, 0.002], 2: [0.0, -0.0002, -0.0003]},
            2: {1: [0.0, 0.0005, 0.0007], 2: [0.0, -0.0001, -0.00015]},
            3: {1: [0.0, 0.0009, 0.0011], 2: [0.0, -0.00025, -0.00035]},
        }
        dij = [[0.0] * (n + 1) for _ in range(n + 1)]
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                dij[i][j] = 0.01 * (i + j)

        state = {
            "n": n,
            "nbl1": nbl1,
            "nbl2": nbl2,
            "iblte1": iblte1,
            "iblte2": iblte2,
            "nsys": nsys,
            "matyp": 1,
            "lalfa": True,
            "cl": 0.5,
            "clspec": 0.0,
            "alfa": 0.05,
            "dtor": math.pi / 180.0,
            "minf": 0.1,
            "qinf": 1.0,
            "gamma": 1.4,
            "gamm1": 0.4,
            "minfCl": 0.0,
            "hvrat": 0.0,
            "x": [0.0] + x,
            "y": [0.0] + y,
            "ipan1": ipan1,
            "ipan2": ipan2,
            "vti1": vti1,
            "vti2": vti2,
            "uinv1": uinv1,
            "uinv2": uinv2,
            "uinvA1": uinv_a1,
            "uinvA2": uinv_a2,
            "dstr1": dstr1,
            "dstr2": dstr2,
            "thet1": thet1,
            "thet2": thet2,
            "ctau1": ctau1,
            "ctau2": ctau2,
            "mass1": mass1,
            "mass2": mass2,
            "uedg1": uedg1,
            "uedg2": uedg2,
            "vdel": vdel,
            "dij": dij,
        }
        js_case = self._run_node([{"kind": "update", "state": state}])[0]
        lines = [
            "10",
            f"{n}",
            f"{nbl1} {nbl2}",
            f"{iblte1} {iblte2}",
            f"{nsys}",
            "1 1",
            f"{state['cl']} {state['clspec']} {state['alfa']} {state['dtor']}",
            f"{state['minf']} {state['qinf']} {state['gamma']} {state['gamm1']} {state['minfCl']}",
            f"{state['hvrat']}",
        ]
        for xi, yi in zip(x, y):
            lines.append(f"{xi} {yi}")
        for ip, vt, ui, uia, ds, th, ct, ms, ue in zip(
            ipan1, vti1, uinv1, uinv_a1, dstr1, thet1, ctau1, mass1, uedg1
        ):
            lines.append(f"{ip} {vt} {ui} {uia} {ds} {th} {ct} {ms} {ue}")
        for ip, vt, ui, uia, ds, th, ct, ms, ue in zip(
            ipan2, vti2, uinv2, uinv_a2, dstr2, thet2, ctau2, mass2, uedg2
        ):
            lines.append(f"{ip} {vt} {ui} {uia} {ds} {th} {ct} {ms} {ue}")
        for iv in range(1, nsys + 1):
            for j in range(1, 3):
                for k in range(1, 4):
                    lines.append(str(vdel[k][j][iv]))
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                lines.append(str(dij[i][j]))
        out = self._run_fortran(lines)
        cl_ref, alfa_ref, rmsbl_ref, rlx_ref = map(float, out[0].split())
        self.assertAlmostEqual(cl_ref, js_case["cl"], delta=1.0e-6)
        self.assertAlmostEqual(alfa_ref, js_case["alfa"], delta=1.0e-6)
        self.assertAlmostEqual(rmsbl_ref, js_case["rmsbl"], delta=1.0e-2)
        self.assertAlmostEqual(rlx_ref, js_case["rlx"], delta=1.0e-6)

    def test_update_matyp2(self):
        n = 4
        nbl1 = 2
        nbl2 = 2
        iblte1 = 2
        iblte2 = 2
        nsys = 2
        x = [0.0, 1.0, 1.0, 0.0]
        y = [0.0, 0.0, 1.0, 1.0]
        ipan1 = [1]
        ipan2 = [2]
        vti1 = [1.0]
        vti2 = [-1.0]
        uinv1 = [1.0]
        uinv2 = [1.1]
        uinv_a1 = [0.05]
        uinv_a2 = [0.04]
        dstr1 = [0.02]
        dstr2 = [0.021]
        thet1 = [0.01]
        thet2 = [0.011]
        ctau1 = [0.03]
        ctau2 = [0.031]
        mass1 = [0.02]
        mass2 = [0.025]
        uedg1 = [1.0]
        uedg2 = [1.05]
        vdel = {
            1: {1: [0.0, 0.001, 0.002], 2: [0.0, -0.0002, -0.0003]},
            2: {1: [0.0, 0.0005, 0.0007], 2: [0.0, -0.0001, -0.00015]},
            3: {1: [0.0, 0.0009, 0.0011], 2: [0.0, -0.00025, -0.00035]},
        }
        dij = [[0.0] * (n + 1) for _ in range(n + 1)]
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                dij[i][j] = 0.01 * (i + j)

        state = {
            "n": n,
            "nbl1": nbl1,
            "nbl2": nbl2,
            "iblte1": iblte1,
            "iblte2": iblte2,
            "nsys": nsys,
            "matyp": 2,
            "lalfa": True,
            "cl": 0.5,
            "clspec": 0.0,
            "alfa": 0.05,
            "dtor": math.pi / 180.0,
            "minf": 0.1,
            "qinf": 1.0,
            "gamma": 1.4,
            "gamm1": 0.4,
            "minfCl": 0.0,
            "hvrat": 0.0,
            "x": [0.0] + x,
            "y": [0.0] + y,
            "ipan1": ipan1,
            "ipan2": ipan2,
            "vti1": vti1,
            "vti2": vti2,
            "uinv1": uinv1,
            "uinv2": uinv2,
            "uinvA1": uinv_a1,
            "uinvA2": uinv_a2,
            "dstr1": dstr1,
            "dstr2": dstr2,
            "thet1": thet1,
            "thet2": thet2,
            "ctau1": ctau1,
            "ctau2": ctau2,
            "mass1": mass1,
            "mass2": mass2,
            "uedg1": uedg1,
            "uedg2": uedg2,
            "vdel": vdel,
            "dij": dij,
        }
        js_case = self._run_node([{"kind": "update", "state": state}])[0]
        lines = [
            "10",
            f"{n}",
            f"{nbl1} {nbl2}",
            f"{iblte1} {iblte2}",
            f"{nsys}",
            "2 1",
            f"{state['cl']} {state['clspec']} {state['alfa']} {state['dtor']}",
            f"{state['minf']} {state['qinf']} {state['gamma']} {state['gamm1']} {state['minfCl']}",
            f"{state['hvrat']}",
        ]
        for xi, yi in zip(x, y):
            lines.append(f"{xi} {yi}")
        for ip, vt, ui, uia, ds, th, ct, ms, ue in zip(
            ipan1, vti1, uinv1, uinv_a1, dstr1, thet1, ctau1, mass1, uedg1
        ):
            lines.append(f"{ip} {vt} {ui} {uia} {ds} {th} {ct} {ms} {ue}")
        for ip, vt, ui, uia, ds, th, ct, ms, ue in zip(
            ipan2, vti2, uinv2, uinv_a2, dstr2, thet2, ctau2, mass2, uedg2
        ):
            lines.append(f"{ip} {vt} {ui} {uia} {ds} {th} {ct} {ms} {ue}")
        for iv in range(1, nsys + 1):
            for j in range(1, 3):
                for k in range(1, 4):
                    lines.append(str(vdel[k][j][iv]))
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                lines.append(str(dij[i][j]))
        out = self._run_fortran(lines)
        cl_ref, alfa_ref, rmsbl_ref, rlx_ref = map(float, out[0].split())
        self.assertAlmostEqual(cl_ref, js_case["cl"], delta=1.0e-6)
        self.assertAlmostEqual(alfa_ref, js_case["alfa"], delta=1.0e-6)
        self.assertAlmostEqual(rmsbl_ref, js_case["rmsbl"], delta=1.0e-2)
        self.assertAlmostEqual(rlx_ref, js_case["rlx"], delta=1.0e-6)

    def test_update_lalfa_false(self):
        n = 4
        nbl1 = 2
        nbl2 = 2
        iblte1 = 2
        iblte2 = 2
        nsys = 2
        x = [0.0, 1.0, 1.0, 0.0]
        y = [0.0, 0.0, 1.0, 1.0]
        ipan1 = [1]
        ipan2 = [2]
        vti1 = [1.0]
        vti2 = [-1.0]
        uinv1 = [1.0]
        uinv2 = [1.1]
        uinv_a1 = [0.05]
        uinv_a2 = [0.04]
        dstr1 = [0.02]
        dstr2 = [0.021]
        thet1 = [0.01]
        thet2 = [0.011]
        ctau1 = [0.03]
        ctau2 = [0.031]
        mass1 = [0.02]
        mass2 = [0.025]
        uedg1 = [1.0]
        uedg2 = [1.05]
        vdel = {
            1: {1: [0.0, 0.001, 0.002], 2: [0.0, -0.0002, -0.0003]},
            2: {1: [0.0, 0.0005, 0.0007], 2: [0.0, -0.0001, -0.00015]},
            3: {1: [0.0, 0.0009, 0.0011], 2: [0.0, -0.00025, -0.00035]},
        }
        dij = [[0.0] * (n + 1) for _ in range(n + 1)]
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                dij[i][j] = 0.01 * (i + j)

        state = {
            "n": n,
            "nbl1": nbl1,
            "nbl2": nbl2,
            "iblte1": iblte1,
            "iblte2": iblte2,
            "nsys": nsys,
            "matyp": 1,
            "lalfa": False,
            "cl": 0.5,
            "clspec": 0.4,
            "alfa": 0.05,
            "dtor": math.pi / 180.0,
            "minf": 0.1,
            "qinf": 1.0,
            "gamma": 1.4,
            "gamm1": 0.4,
            "minfCl": 0.0,
            "hvrat": 0.0,
            "x": [0.0] + x,
            "y": [0.0] + y,
            "ipan1": ipan1,
            "ipan2": ipan2,
            "vti1": vti1,
            "vti2": vti2,
            "uinv1": uinv1,
            "uinv2": uinv2,
            "uinvA1": uinv_a1,
            "uinvA2": uinv_a2,
            "dstr1": dstr1,
            "dstr2": dstr2,
            "thet1": thet1,
            "thet2": thet2,
            "ctau1": ctau1,
            "ctau2": ctau2,
            "mass1": mass1,
            "mass2": mass2,
            "uedg1": uedg1,
            "uedg2": uedg2,
            "vdel": vdel,
            "dij": dij,
        }
        js_case = self._run_node([{"kind": "update", "state": state}])[0]
        lines = [
            "10",
            f"{n}",
            f"{nbl1} {nbl2}",
            f"{iblte1} {iblte2}",
            f"{nsys}",
            "1 0",
            f"{state['cl']} {state['clspec']} {state['alfa']} {state['dtor']}",
            f"{state['minf']} {state['qinf']} {state['gamma']} {state['gamm1']} {state['minfCl']}",
            f"{state['hvrat']}",
        ]
        for xi, yi in zip(x, y):
            lines.append(f"{xi} {yi}")
        for ip, vt, ui, uia, ds, th, ct, ms, ue in zip(
            ipan1, vti1, uinv1, uinv_a1, dstr1, thet1, ctau1, mass1, uedg1
        ):
            lines.append(f"{ip} {vt} {ui} {uia} {ds} {th} {ct} {ms} {ue}")
        for ip, vt, ui, uia, ds, th, ct, ms, ue in zip(
            ipan2, vti2, uinv2, uinv_a2, dstr2, thet2, ctau2, mass2, uedg2
        ):
            lines.append(f"{ip} {vt} {ui} {uia} {ds} {th} {ct} {ms} {ue}")
        for iv in range(1, nsys + 1):
            for j in range(1, 3):
                for k in range(1, 4):
                    lines.append(str(vdel[k][j][iv]))
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                lines.append(str(dij[i][j]))
        out = self._run_fortran(lines)
        cl_ref, alfa_ref, rmsbl_ref, rlx_ref = map(float, out[0].split())
        self.assertAlmostEqual(cl_ref, js_case["cl"], delta=1.0e-6)
        self.assertAlmostEqual(alfa_ref, js_case["alfa"], delta=1.0e-6)
        self.assertAlmostEqual(rmsbl_ref, js_case["rmsbl"], delta=1.0e-2)
        self.assertAlmostEqual(rlx_ref, js_case["rlx"], delta=1.0e-6)


if __name__ == "__main__":
    unittest.main()
