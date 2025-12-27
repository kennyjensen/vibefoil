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
from python.xbl import mrchdu, mrchue
from python.xsolve import gauss as gauss_base
from python.tests.test_compare_mrchue import build_bl, build_context


def gauss1(nsiz, nn, z, r, nrhs):
    rmat = [[0.0] * (nrhs + 1) for _ in range(nn + 1)]
    for i in range(1, nn + 1):
        rmat[i][1] = r[i]
    gauss_base(nsiz, nn, z, rmat, nrhs)
    for i in range(1, nn + 1):
        r[i] = rmat[i][1]


def max_abs_list_diff(a, b):
    if len(a) != len(b):
        return float("inf")
    return max(abs(ai - bi) for ai, bi in zip(a, b))


class TestWasmMrchueParity(unittest.TestCase):
    def run_wasm_mode(self, ctx, bl, mode):
        if not shutil.which("node"):
            self.skipTest("node is required for wasm parity checks")

        compare_js = ROOT / "wasm" / "build-wasm" / "compare_xbl.js"
        compare_wasm = ROOT / "wasm" / "build-wasm" / "compare_xbl.wasm"
        if not compare_js.exists() or not compare_wasm.exists():
            self.skipTest("compare_xbl.js not built; run emcmake cmake and build wasm first")

        ctx_payload = ctx.__dict__ if hasattr(ctx, "__dict__") else ctx
        bl_payload = bl.__dict__ if hasattr(bl, "__dict__") else bl
        if "HTARG" in ctx_payload:
            ctx_payload = dict(ctx_payload)
            ctx_payload.pop("HTARG", None)
        payload = {"ctx": ctx_payload, "bl": bl_payload}

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
                    mode,
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

    def assert_close(self, name, got, expected, tol=1.0e-6, rel_tol=1.0e-9):
        scale = max(1.0, abs(expected))
        allowed = max(tol, rel_tol * scale)
        self.assertLessEqual(abs(got - expected), allowed, f"{name} {got} vs {expected}")

    def compare_side_arrays(self, ctx_py, ctx_out, key, max_nbl, is_, tol):
        src = getattr(ctx_py, key)
        py_values = [0.0] * (max_nbl + 1)
        js_values = [0.0] * (max_nbl + 1)
        for ibl in range(1, max_nbl + 1):
            if key == "HTARG":
                py_values[ibl] = src[is_][ibl]
                js_values[ibl] = ctx_out[key][is_][ibl]
            else:
                py_values[ibl] = src[ibl][is_]
                js_values[ibl] = ctx_out[key][ibl][is_]
        self.assertLessEqual(max_abs_list_diff(py_values, js_values), tol, f"{key} side {is_}")

    def test_mrchue_and_mrchdu(self):
        orig_gauss = xbl_mod.gauss
        xbl_mod.gauss = gauss1
        try:
            cases = [
                {"ides": 12, "minf": 0.0, "reinf": 1.0e6, "alfa": 0.0, "waklen": 1.0},
                {"ides": 2412, "minf": 0.1, "reinf": 3.0e6, "alfa": 0.0, "waklen": 1.0},
            ]

            tol = 1.0e-6

            for case in cases:
                with self.subTest(case=case):
                    ctx = build_context(
                        case["ides"],
                        case["minf"],
                        case["reinf"],
                        case["alfa"],
                        case["waklen"],
                    )
                    bl = build_bl(ctx)

                    wasm_mrchue = self.run_wasm_mode(ctx, bl, "mrchue")
                    mrchue(ctx, bl)

                    ctx_out = wasm_mrchue["ctx"]
                    bl_out = wasm_mrchue["bl"]

                    self.assertEqual(ctx.ITRAN[1], ctx_out["ITRAN"][1])
                    self.assertEqual(ctx.ITRAN[2], ctx_out["ITRAN"][2])
                    self.assert_close("ctx.XSSITR[1]", ctx_out["XSSITR"][1], ctx.XSSITR[1], tol=tol)
                    self.assert_close("ctx.XSSITR[2]", ctx_out["XSSITR"][2], ctx.XSSITR[2], tol=tol)

                    max_nbl = max(ctx.NBL[1], ctx.NBL[2])
                    for is_ in (1, 2):
                        self.compare_side_arrays(ctx, ctx_out, "THET", max_nbl, is_, tol)
                        self.compare_side_arrays(ctx, ctx_out, "DSTR", max_nbl, is_, tol)
                        self.compare_side_arrays(ctx, ctx_out, "CTAU", max_nbl, is_, tol)
                        self.compare_side_arrays(ctx, ctx_out, "UEDG", max_nbl, is_, tol)
                        if "HTARG" in ctx_out:
                            self.compare_side_arrays(ctx, ctx_out, "HTARG", max_nbl, is_, tol)

                    self.assert_close("bl.AMPL1", bl_out["AMPL1"], bl.AMPL1, tol=tol)
                    self.assert_close("bl.AMPL2", bl_out["AMPL2"], bl.AMPL2, tol=tol)

                    wasm_mrchdu = self.run_wasm_mode(ctx_out, bl_out, "mrchdu")
                    mrchdu(ctx, bl)

                    ctx_out2 = wasm_mrchdu["ctx"]
                    bl_out2 = wasm_mrchdu["bl"]

                    self.assertEqual(ctx.ITRAN[1], ctx_out2["ITRAN"][1])
                    self.assertEqual(ctx.ITRAN[2], ctx_out2["ITRAN"][2])
                    self.assert_close("ctx.XSSITR[1]", ctx_out2["XSSITR"][1], ctx.XSSITR[1], tol=tol)
                    self.assert_close("ctx.XSSITR[2]", ctx_out2["XSSITR"][2], ctx.XSSITR[2], tol=tol)

                    for is_ in (1, 2):
                        self.compare_side_arrays(ctx, ctx_out2, "THET", max_nbl, is_, tol)
                        self.compare_side_arrays(ctx, ctx_out2, "DSTR", max_nbl, is_, tol)
                        self.compare_side_arrays(ctx, ctx_out2, "CTAU", max_nbl, is_, tol)
                        self.compare_side_arrays(ctx, ctx_out2, "UEDG", max_nbl, is_, tol)
                        self.compare_side_arrays(ctx, ctx_out2, "MASS", max_nbl, is_, tol)
                        self.compare_side_arrays(ctx, ctx_out2, "TAU", max_nbl, is_, tol)
                        self.compare_side_arrays(ctx, ctx_out2, "DIS", max_nbl, is_, tol)
                        self.compare_side_arrays(ctx, ctx_out2, "CTQ", max_nbl, is_, tol)
                        self.compare_side_arrays(ctx, ctx_out2, "DELT", max_nbl, is_, tol)
                        self.compare_side_arrays(ctx, ctx_out2, "TSTR", max_nbl, is_, tol)

                    self.assert_close("bl.AMPL1", bl_out2["AMPL1"], bl.AMPL1, tol=tol)
                    self.assert_close("bl.AMPL2", bl_out2["AMPL2"], bl.AMPL2, tol=tol)
        finally:
            xbl_mod.gauss = orig_gauss


if __name__ == "__main__":
    unittest.main()
