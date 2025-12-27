import math
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
from python.xbl import XBlState, blpini
from python.xoper import viscal
from python.xsolve import gauss as gauss_base
from python.tests.test_compare_viscal import build_viscal_context


def gauss1(nsiz, nn, z, r, nrhs):
    rmat = [[0.0] * (nrhs + 1) for _ in range(nn + 1)]
    for i in range(1, nn + 1):
        rmat[i][1] = r[i]
    gauss_base(nsiz, nn, z, rmat, nrhs)
    for i in range(1, nn + 1):
        r[i] = rmat[i][1]


class TestWasmViscalParity(unittest.TestCase):
    def run_wasm_viscal(self, ctx, bl, niter):
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
                [
                    "node",
                    str(runner_js),
                    str(input_path),
                    "--out",
                    str(output_path),
                    "--mode",
                    "viscal",
                    "--niter",
                    str(niter),
                ],
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

    def assert_close(self, name, got, expected, tol=1.0e-2, rel_tol=1.0e-6):
        scale = max(1.0, abs(expected))
        allowed = max(tol, rel_tol * scale)
        self.assertLessEqual(abs(got - expected), allowed, f"{name} {got} vs {expected}")

    def test_viscal_one_iter(self):
        orig_gauss = xbl_mod.gauss
        xbl_mod.gauss = gauss1
        try:
            cases = [
                {"ides": 12, "minf": 0.0, "reinf": 1.0e6, "alfa": 0.0, "waklen": 1.0},
                {"ides": 2412, "minf": 0.1, "reinf": 3.0e6, "alfa": 0.0, "waklen": 1.0},
            ]

            for case in cases:
                with self.subTest(case=case):
                    ctx = build_viscal_context(
                        case["ides"],
                        case["minf"],
                        case["reinf"],
                        case["alfa"],
                        case["waklen"],
                    )
                    cosa = math.cos(ctx.ALFA)
                    sina = math.sin(ctx.ALFA)
                    for i in range(1, ctx.N + 1):
                        ctx.GAM[i] = cosa * ctx.GAMU[i][1] + sina * ctx.GAMU[i][2]

                    bl = XBlState()
                    blpini(bl)

                    wasm_out = self.run_wasm_viscal(ctx, bl, 1)
                    viscal(ctx, bl, 1)

                    ctx_out = wasm_out["ctx"]

                    self.assert_close("ctx.CL", ctx_out["CL"], ctx.CL)
                    self.assert_close("ctx.CM", ctx_out["CM"], ctx.CM)
                    self.assert_close("ctx.CD", ctx_out["CD"], ctx.CD)
                    self.assert_close("ctx.CDF", ctx_out["CDF"], ctx.CDF)
                    self.assert_close("ctx.CDP", ctx_out["CDP"], ctx.CDP, tol=1.0e-3, rel_tol=1.0e-6)

                    self.assertEqual(ctx.ITRAN[1], ctx_out["ITRAN"][1])
                    self.assertEqual(ctx.ITRAN[2], ctx_out["ITRAN"][2])
                    self.assert_close("ctx.XSSITR[1]", ctx_out["XSSITR"][1], ctx.XSSITR[1], tol=1.0e-3, rel_tol=1.0e-6)
                    self.assert_close("ctx.XSSITR[2]", ctx_out["XSSITR"][2], ctx.XSSITR[2], tol=1.0e-3, rel_tol=1.0e-6)
        finally:
            xbl_mod.gauss = orig_gauss


if __name__ == "__main__":
    unittest.main()
