import copy
import json
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
from python.xsolve import gauss as gauss_base
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


def metrics_3d(mat, nsys):
    total = 0.0
    sumsq = 0.0
    maxabs = 0.0
    for k in range(1, 4):
        for j in range(1, 3):
            for iv in range(1, nsys + 1):
                val = mat[k][j][iv]
                total += val
                sumsq += val * val
                maxabs = max(maxabs, abs(val))
    return {"sum": total, "sumsq": sumsq, "maxabs": maxabs}


class TestSetblParity(unittest.TestCase):
    def test_setbl_metrics(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        ctx, bl, ctx_payload = setup_state()
        setbl(ctx, bl)

        py_metrics = {
            "VA": metrics_3d(ctx.VA, ctx.NSYS),
            "VB": metrics_3d(ctx.VB, ctx.NSYS),
            "VDEL": metrics_3d(ctx.VDEL, ctx.NSYS),
            "VM": metrics_3d(ctx.VM, ctx.NSYS),
        }

        script = pathlib.Path(__file__).with_name("compare_viscal_subfuncs.mjs")
        payload = build_payload_ctx(ctx_payload, 1.0, 9.0)
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]["setbl"]

        tol = 1.0e-3
        rel_tol = 1.0e-12
        for key in ("VA", "VB", "VDEL", "VM"):
            py_m = py_metrics[key]
            js_m = js_results[key]
            sum_tol = max(tol, rel_tol * max(1.0, abs(js_m["sum"])))
            sumsq_tol = max(tol, rel_tol * max(1.0, abs(js_m["sumsq"])))
            maxabs_tol = max(tol, rel_tol * max(1.0, abs(js_m["maxabs"])))
            self.assertLessEqual(abs(py_m["sum"] - js_m["sum"]), sum_tol, f"{key} sum")
            self.assertLessEqual(abs(py_m["sumsq"] - js_m["sumsq"]), sumsq_tol, f"{key} sumsq")
            self.assertLessEqual(abs(py_m["maxabs"] - js_m["maxabs"]), maxabs_tol, f"{key} maxabs")
