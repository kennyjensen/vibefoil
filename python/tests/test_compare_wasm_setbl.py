import copy
import json
import pathlib
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import python.xbl as xbl_mod
from python.xbl import XBlState, blpini, setbl, iblsys
from python.xpanel import (
    xywake,
    qwcalc,
    qiset,
    stfind,
    iblpan,
    xicalc,
    uicalc,
    qdcalc,
    qvfue,
    gamqv,
    stmove,
)
from python.xsolve import gauss as gauss_base
from python.tests.test_compare_viscal_subfuncs import build_context


def gauss1(nsiz, nn, z, r, nrhs):
    rmat = [[0.0] * (nrhs + 1) for _ in range(nn + 1)]
    for i in range(1, nn + 1):
        rmat[i][1] = r[i]
    gauss_base(nsiz, nn, z, rmat, nrhs)
    for i in range(1, nn + 1):
        r[i] = rmat[i][1]


def setup_state(ides, minf, reinf, alfa, waklen):
    ctx = build_context(ides, minf, reinf, alfa, waklen)
    bl = XBlState()
    blpini(bl)

    xywake(ctx)
    qwcalc(ctx)
    qiset(ctx)
    stfind(ctx)
    iblpan(ctx)
    xicalc(ctx)
    iblsys(ctx)
    uicalc(ctx)

    for is_ in (1, 2):
        for ibl in range(1, ctx.NBL[is_] + 1):
            ctx.UEDG[ibl][is_] = ctx.UINV[ibl][is_]

    qdcalc(ctx)
    qvfue(ctx)
    gamqv(ctx)
    stmove(ctx)

    return ctx, bl


def metrics_nested(value):
    total = 0.0
    sumsq = 0.0
    maxabs = 0.0
    count = 0
    stack = [value]
    while stack:
        item = stack.pop()
        if isinstance(item, (int, float, bool)):
            val = float(item)
            total += val
            sumsq += val * val
            maxabs = max(maxabs, abs(val))
            count += 1
        elif isinstance(item, list):
            stack.extend(item)
    return {"sum": total, "sumsq": sumsq, "maxabs": maxabs, "count": count}


class TestWasmSetblParity(unittest.TestCase):
    def run_wasm_setbl(self, ctx, bl):
        if not shutil.which("node"):
            self.skipTest("node is required for wasm parity checks")

        compare_js = ROOT / "wasm" / "build-wasm" / "compare_xbl.js"
        compare_wasm = ROOT / "wasm" / "build-wasm" / "compare_xbl.wasm"
        if not compare_js.exists() or not compare_wasm.exists():
            self.skipTest("compare_xbl.js not built; run emcmake cmake and build wasm first")

        payload = {"ctx": ctx.__dict__, "bl": bl.__dict__}

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)
            input_path = tmpdir / "state.json"
            output_path = tmpdir / "out.json"
            runner_js = tmpdir / "compare_xbl.cjs"
            runner_wasm = tmpdir / "compare_xbl.wasm"
            with input_path.open("w") as handle:
                json.dump(payload, handle)

            runner_js.write_text(compare_js.read_text())
            runner_wasm.write_bytes(compare_wasm.read_bytes())

            proc = subprocess.run(
                ["node", str(runner_js), str(input_path), "--out", str(output_path)],
                text=True,
                capture_output=True,
            )
            if proc.returncode != 0:
                raise RuntimeError(
                    "compare_xbl.js failed\n"
                    f"stdout:\n{proc.stdout}\n"
                    f"stderr:\n{proc.stderr}\n"
                )

            with output_path.open() as handle:
                return json.load(handle)

    def assert_close(self, name, got, expected, tol=1.0e-5, rel_tol=1.0e-8):
        scale = max(1.0, abs(expected))
        allowed = max(tol, rel_tol * scale)
        self.assertLessEqual(abs(got - expected), allowed, f"{name} {got} vs {expected}")

    def assert_metrics_close(self, label, got, expected):
        self.assertEqual(got["count"], expected["count"], f"{label} count")
        self.assert_close(f"{label} sum", got["sum"], expected["sum"], tol=1.0e-4, rel_tol=2.0e-8)
        self.assert_close(f"{label} sumsq", got["sumsq"], expected["sumsq"], tol=1.0e-4, rel_tol=2.0e-8)
        self.assert_close(f"{label} maxabs", got["maxabs"], expected["maxabs"], tol=1.0e-4, rel_tol=2.0e-8)

    def run_case(self, ides, minf, reinf, alfa, waklen):
        orig_gauss = xbl_mod.gauss
        xbl_mod.gauss = gauss1
        try:
            ctx, bl = setup_state(ides, minf, reinf, alfa, waklen)
            ctx_in = copy.deepcopy(ctx)
            bl_in = copy.deepcopy(bl)
            setbl(ctx, bl)
        finally:
            xbl_mod.gauss = orig_gauss

        wasm_out = self.run_wasm_setbl(ctx_in, bl_in)
        ctx_out = wasm_out["ctx"]
        bl_out = wasm_out["bl"]

        scalar_fields = [
            "CL",
            "CM",
            "CD",
            "CDP",
            "CDF",
            "CL_ALF",
            "CL_MSQ",
            "MINF",
            "REINF",
            "QINF",
            "TKLAM",
            "CPSTAR",
            "QSTAR",
            "SLE",
            "XLE",
            "YLE",
            "XTE",
            "YTE",
            "SST",
            "ALFA",
            "ADEG",
            "XCMREF",
            "YCMREF",
            "CHORD",
            "DSTE",
            "ANTE",
            "ASTE",
        ]
        for name in scalar_fields:
            self.assert_close(f"ctx.{name}", ctx_out[name], getattr(ctx, name))

        bl_scalar_fields = [
            "XT",
            "XT_A1",
            "XT_A2",
            "XT_X1",
            "XT_X2",
            "XT_T1",
            "XT_T2",
            "XT_D1",
            "XT_D2",
            "XT_U1",
            "XT_U2",
            "XT_MS",
            "XT_RE",
            "AMPL1",
            "AMPL2",
            "AMCRIT",
            "HK1",
            "HK2",
            "RT1",
            "RT2",
            "T1",
            "T2",
            "D1",
            "D2",
            "U1",
            "U2",
        ]
        for name in bl_scalar_fields:
            self.assert_close(f"bl.{name}", bl_out[name], getattr(bl, name))

        array_fields = [
            "VA",
            "VB",
            "VDEL",
            "VM",
            "UEDG",
            "THET",
            "DSTR",
            "CTAU",
            "MASS",
            "TAU",
            "DIS",
            "CTQ",
            "DELT",
            "TSTR",
        ]
        for name in array_fields:
            metrics_py = metrics_nested(getattr(ctx, name))
            metrics_wasm = metrics_nested(ctx_out[name])
            self.assert_metrics_close(f"ctx.{name}", metrics_wasm, metrics_py)

        bl_array_fields = ["VS1", "VS2", "VSREZ"]
        for name in bl_array_fields:
            metrics_py = metrics_nested(getattr(bl, name))
            metrics_wasm = metrics_nested(bl_out[name])
            self.assert_metrics_close(f"bl.{name}", metrics_wasm, metrics_py)

    def test_setbl_parity_naca2412(self):
        self.run_case(2412, 0.1, 3.0e6, 0.0, 1.0)

    def test_setbl_parity_naca0012(self):
        self.run_case(12, 0.08, 6.0e6, 2.0 * (3.141592653589793 / 180.0), 1.0)
