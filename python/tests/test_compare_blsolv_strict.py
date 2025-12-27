import copy
import json
import math
import pathlib
import shutil
import subprocess
import sys
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
from python.xsolve import gauss as gauss_base, blsolv
from python.tests.test_compare_viscal_subfuncs import build_context, build_payload_ctx


def setup_state():
    def gauss1(nsiz, nn, z, r, nrhs):
        rmat = [[0.0] * (nrhs + 1) for _ in range(nn + 1)]
        for i in range(1, nn + 1):
            rmat[i][1] = r[i]
        gauss_base(nsiz, nn, z, rmat, nrhs)
        for i in range(1, nn + 1):
            r[i] = rmat[i][1]

    xbl_mod.gauss = gauss1

    ctx = build_context(2412, 0.1, 3.0e6, 0.0, 1.0)
    bl = XBlState()
    blpini(bl)

    xywake(ctx)
    qwcalc(ctx)
    ctx_payload = copy.deepcopy(ctx)
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

    return ctx, bl, ctx_payload


def flatten_tensor3(mat, d1, d2, d3):
    data = []
    for k in range(1, d1 + 1):
        for j in range(1, d2 + 1):
            for i in range(1, d3 + 1):
                data.append(mat[k][j][i])
    return data


class TestBlsolvStrictParity(unittest.TestCase):
    def assert_array_close(self, py_data, js_data, abs_tol, rel_tol, label):
        self.assertEqual(len(py_data), len(js_data), f"{label} length")
        for idx, (py_val, js_val) in enumerate(zip(py_data, js_data)):
            if not (math.isfinite(py_val) and math.isfinite(js_val)):
                self.fail(f"{label} idx {idx}: non-finite {py_val} vs {js_val}")
            tol = max(abs_tol, rel_tol * max(1.0, abs(js_val)))
            if abs(py_val - js_val) > tol:
                self.fail(f"{label} idx {idx}: {py_val} vs {js_val} (tol {tol})")

    def test_blsolv_strict_parity(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        ctx, bl, ctx_payload = setup_state()
        setbl(ctx, bl)
        blsolv(ctx)

        script = pathlib.Path(__file__).with_name("compare_blsolv_strict.mjs")
        payload = build_payload_ctx(ctx_payload, 1.0, 9.0)
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        self.assertEqual(ctx.NSYS, js_results["nsys"], "NSYS")

        abs_tol = 1.0e-6
        rel_tol = 1.0e-12

        va_py = flatten_tensor3(ctx.VA, 3, 2, ctx.NSYS)
        vdel_py = flatten_tensor3(ctx.VDEL, 3, 2, ctx.NSYS)
        vm_py = flatten_tensor3(ctx.VM, 3, ctx.NSYS, ctx.NSYS)

        self.assert_array_close(va_py, js_results["va"]["data"], abs_tol, rel_tol, "VA")
        self.assert_array_close(vdel_py, js_results["vdel"]["data"], abs_tol, rel_tol, "VDEL")
        self.assert_array_close(vm_py, js_results["vm"]["data"], abs_tol, rel_tol, "VM")
