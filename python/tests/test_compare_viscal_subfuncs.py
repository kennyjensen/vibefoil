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
from python.xbl import XFoilState, XBlState, blpini, iblsys
from python.xfoil import comset, naca, cpcalc, clcalc, cdcalc
from python.xpanel import (
    ggcalc,
    qdcalc,
    qwcalc,
    qiset,
    stfind,
    iblpan,
    xicalc,
    uicalc,
    qvfue,
    gamqv,
    stmove,
    xywake,
)
from python.xsolve import gauss as gauss_base, blsolv
from python.xbl import setbl, update


def build_context(ides, minf, reinf, alfa, waklen):
    ctx = XFoilState()
    ctx.NPAN = 160
    ctx.CVPAR = 1.0
    ctx.CTERAT = 0.15
    ctx.CTRRAT = 0.2
    ctx.XSREF1 = 1.0
    ctx.XSREF2 = 1.0
    ctx.XPREF1 = 1.0
    ctx.XPREF2 = 1.0

    ctx.WAKLEN = waklen
    ctx.ALFA = alfa
    ctx.ADEG = alfa / ctx.DTOR
    ctx.QINF = 1.0
    ctx.MINF = minf
    ctx.MINF1 = minf
    ctx.REINF = reinf
    ctx.REINF1 = reinf
    ctx.LALFA = True
    ctx.VACCEL = 0.01
    ctx.LVISC = True

    ctx.ACRIT[1] = 9.0
    ctx.ACRIT[2] = 9.0
    ctx.XSTRIP[1] = 1.0
    ctx.XSTRIP[2] = 1.0
    ctx.XCMREF = 0.25

    naca(ctx, ides)
    comset(ctx)
    ggcalc(ctx)
    for i in range(1, ctx.N + 1):
        ctx.GAM[i] = 1.0 if i <= ctx.N // 2 else -1.0

    return ctx


def build_payload_ctx(ctx, waklen, ncrit):
    nw = ctx.N // 12 + 10 * int(waklen)
    total = ctx.N + nw
    ctx.NW = nw

    def slice_array(arr):
        return [arr[i] for i in range(1, total + 1)]

    gam = [0.0] * (ctx.N + 1)
    for i in range(1, ctx.N + 1):
        gam[i] = ctx.GAM[i]

    qinvu1 = [0.0] * (total + 1)
    qinvu2 = [0.0] * (total + 1)
    for i in range(1, total + 1):
        qinvu1[i] = ctx.QINVU[i][1]
        qinvu2[i] = ctx.QINVU[i][2]

    return {
        "ctx": {
            "N": ctx.N,
            "NW": nw,
            "WAKLEN": waklen,
            "ALFA": ctx.ALFA,
            "MINF": ctx.MINF,
            "REINF": ctx.REINF,
            "QINF": ctx.QINF,
            "NCRIT": ncrit,
            "SHARP": ctx.SHARP,
            "ANTE": ctx.ANTE,
            "ASTE": ctx.ASTE,
            "DSTE": ctx.DSTE,
            "XTE": ctx.XTE,
            "YTE": ctx.YTE,
            "CHORD": ctx.CHORD,
            "XCMREF": ctx.XCMREF,
            "YCMREF": ctx.YCMREF,
            "X": slice_array(ctx.X),
            "Y": slice_array(ctx.Y),
            "XP": slice_array(ctx.XP),
            "YP": slice_array(ctx.YP),
            "S": slice_array(ctx.S),
            "NX": slice_array(ctx.NX),
            "NY": slice_array(ctx.NY),
            "APANEL": slice_array(ctx.APANEL),
            "GAM": gam,
            "QINVU1": qinvu1,
            "QINVU2": qinvu2,
        }
    }


def sample_indices(length):
    if length <= 0:
        return [1]
    mid = max(1, (length + 1) // 2)
    return [1, mid, length]


def metrics_1d(arr, length):
    total = 0.0
    sumsq = 0.0
    maxabs = 0.0
    for i in range(1, length + 1):
        val = arr[i]
        total += val
        sumsq += val * val
        maxabs = max(maxabs, abs(val))
    samples = [arr[i] for i in sample_indices(length)]
    return {"sum": total, "sumsq": sumsq, "maxabs": maxabs, "samples": samples}


def metrics_1d_panel(arr, length):
    total = 0.0
    sumsq = 0.0
    maxabs = 0.0
    for i in range(1, length + 1):
        val = arr[i]
        total += val
        sumsq += val * val
        maxabs = max(maxabs, abs(val))
    samples = [arr[i] for i in sample_indices(length)]
    return {"sum": total, "sumsq": sumsq, "maxabs": maxabs, "samples": samples}


def metrics_2d(mat, nbl1, nbl2):
    total = 0.0
    sumsq = 0.0
    maxabs = 0.0
    for is_ in (1, 2):
        nbl = nbl1 if is_ == 1 else nbl2
        for ibl in range(1, nbl + 1):
            val = mat[ibl][is_]
            total += val
            sumsq += val * val
            maxabs = max(maxabs, abs(val))
    s1 = min(nbl1, max(1, (nbl1 + 1) // 2))
    s2 = min(nbl2, max(1, (nbl2 + 1) // 2))
    samples = [
        mat[min(2, nbl1)][1],
        mat[nbl1][1],
        mat[min(2, nbl2)][2],
        mat[nbl2][2],
        mat[s1][1],
        mat[s2][2],
    ]
    return {"sum": total, "sumsq": sumsq, "maxabs": maxabs, "samples": samples}


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
    mid = max(1, (nsys + 1) // 2)
    samples = [
        mat[1][1][1],
        mat[2][1][mid],
        mat[3][2][nsys],
    ]
    return {"sum": total, "sumsq": sumsq, "maxabs": maxabs, "samples": samples}


class TestViscalSubfuncParity(unittest.TestCase):
    def test_viscal_subfunctions(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        def gauss1(nsiz, nn, z, r, nrhs):
            rmat = [[0.0] * (nrhs + 1) for _ in range(nn + 1)]
            for i in range(1, nn + 1):
                rmat[i][1] = r[i]
            gauss_base(nsiz, nn, z, rmat, nrhs)
            for i in range(1, nn + 1):
                r[i] = rmat[i][1]

        xbl_mod.gauss = gauss1

        case = {"ides": 2412, "minf": 0.1, "reinf": 3.0e6, "alfa": 0.0, "waklen": 1.0}
        ctx = build_context(case["ides"], case["minf"], case["reinf"], case["alfa"], case["waklen"])
        bl = XBlState()
        blpini(bl)

        script = pathlib.Path(__file__).with_name("compare_viscal_subfuncs.mjs")
        payload = build_payload_ctx(ctx, case["waklen"], 9.0)

        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        xywake(ctx)
        total = ctx.N + ctx.NW
        py_results = {
            "xywake": {
                "X": metrics_1d_panel(ctx.X, total),
                "Y": metrics_1d_panel(ctx.Y, total),
                "S": metrics_1d_panel(ctx.S, total),
                "NX": metrics_1d_panel(ctx.NX, total),
                "NY": metrics_1d_panel(ctx.NY, total),
                "APANEL": metrics_1d_panel(ctx.APANEL, total),
            }
        }

        qwcalc(ctx)
        qinvu1 = [0.0] * (total + 1)
        qinvu2 = [0.0] * (total + 1)
        for i in range(1, total + 1):
            qinvu1[i] = ctx.QINVU[i][1]
            qinvu2[i] = ctx.QINVU[i][2]
        py_results["qwcalc"] = {
            "QINVU1": metrics_1d(qinvu1, total),
            "QINVU2": metrics_1d(qinvu2, total),
        }

        qiset(ctx)
        py_results["qiset"] = {
            "QINV": metrics_1d(ctx.QINV, total),
            "QINV_A": metrics_1d(ctx.QINV_A, total),
        }

        stfind(ctx)
        py_results["stfind"] = {"IST": ctx.IST, "SST": ctx.SST}

        iblpan(ctx)
        py_results["iblpan"] = {
            "NBL": [ctx.NBL[0], ctx.NBL[1], ctx.NBL[2]],
            "IBLTE": [ctx.IBLTE[0], ctx.IBLTE[1], ctx.IBLTE[2]],
            "IPAN": metrics_2d(ctx.IPAN, ctx.NBL[1], ctx.NBL[2]),
            "VTI": metrics_2d(ctx.VTI, ctx.NBL[1], ctx.NBL[2]),
        }

        xicalc(ctx)
        py_results["xicalc"] = {"XSSI": metrics_2d(ctx.XSSI, ctx.NBL[1], ctx.NBL[2])}

        iblsys(ctx)
        py_results["iblsys"] = {
            "NSYS": ctx.NSYS,
            "ISYS": metrics_2d(ctx.ISYS, ctx.NBL[1], ctx.NBL[2]),
        }

        uicalc(ctx)
        py_results["uicalc"] = {
            "UINV": metrics_2d(ctx.UINV, ctx.NBL[1], ctx.NBL[2]),
            "UINV_A": metrics_2d(ctx.UINV_A, ctx.NBL[1], ctx.NBL[2]),
        }

        for is_ in (1, 2):
            for ibl in range(1, ctx.NBL[is_] + 1):
                ctx.UEDG[ibl][is_] = ctx.UINV[ibl][is_]

        qdcalc(ctx)

        qvfue(ctx)
        py_results["qvfue"] = {"QVIS": metrics_1d(ctx.QVIS, total)}

        gamqv(ctx)
        py_results["gamqv"] = {
            "GAM": metrics_1d_panel(ctx.GAM, ctx.N),
            "GAM_A": metrics_1d_panel(ctx.GAM_A, ctx.N),
        }

        stmove(ctx)
        py_results["stmove"] = {
            "IST": ctx.IST,
            "XSSI": metrics_2d(ctx.XSSI, ctx.NBL[1], ctx.NBL[2]),
        }

        setbl(ctx, bl)
        py_results["setbl"] = {
            "REYBL": bl.REYBL,
            "RSTBL": bl.RSTBL,
            "HSTINV": bl.HSTINV,
            "VA": metrics_3d(ctx.VA, ctx.NSYS),
            "VB": metrics_3d(ctx.VB, ctx.NSYS),
            "VDEL": metrics_3d(ctx.VDEL, ctx.NSYS),
            "VM": metrics_3d(ctx.VM, ctx.NSYS),
        }

        blsolv(ctx)
        py_results["blsolv"] = {"VDEL": metrics_3d(ctx.VDEL, ctx.NSYS)}

        update(ctx, bl)
        py_results["update"] = {
            "CTAU": metrics_2d(ctx.CTAU, ctx.NBL[1], ctx.NBL[2]),
            "THET": metrics_2d(ctx.THET, ctx.NBL[1], ctx.NBL[2]),
            "DSTR": metrics_2d(ctx.DSTR, ctx.NBL[1], ctx.NBL[2]),
            "UEDG": metrics_2d(ctx.UEDG, ctx.NBL[1], ctx.NBL[2]),
            "CL": ctx.CL,
            "RMSBL": ctx.RMSBL,
            "RMXBL": ctx.RMXBL,
            "RLX": ctx.RLX,
        }

        cpcalc(total, ctx.QINV, ctx.QINF, ctx.MINF, ctx.CPI)
        cl_res = clcalc(
            ctx.N,
            ctx.X,
            ctx.Y,
            ctx.GAM,
            ctx.GAM_A,
            ctx.ALFA,
            ctx.MINF,
            ctx.QINF,
            ctx.XCMREF,
            ctx.YCMREF,
        )
        ctx.CL, ctx.CM, ctx.CDP, ctx.CL_ALF, ctx.CL_MSQ = cl_res
        cdcalc(ctx)
        py_results["cpcalc"] = {"CPI": metrics_1d(ctx.CPI, total)}
        py_results["clcalc"] = {"CL": ctx.CL, "CM": ctx.CM, "CDP": ctx.CDP}
        py_results["cdcalc"] = {"CD": ctx.CD, "CDF": ctx.CDF}

        tol = 1.0e-3
        tol_scalar = 1.0e-3

        def assert_metrics(py_m, js_m, tol_val=tol, rel_tol=None):
            if js_m["sum"] is None or js_m["sumsq"] is None or js_m["maxabs"] is None:
                self.fail("JS metrics contain non-finite values")
            if any(val is None for val in js_m["samples"]):
                self.fail("JS metrics contain non-finite sample values")
            sum_tol = tol_val
            sumsq_tol = tol_val
            maxabs_tol = tol_val
            if rel_tol is not None:
                sum_tol = max(sum_tol, rel_tol * max(1.0, abs(js_m["sum"])))
                sumsq_tol = max(sumsq_tol, rel_tol * max(1.0, abs(js_m["sumsq"])))
                maxabs_tol = max(maxabs_tol, rel_tol * max(1.0, abs(js_m["maxabs"])))
            self.assertLessEqual(abs(py_m["sum"] - js_m["sum"]), sum_tol)
            self.assertLessEqual(abs(py_m["sumsq"] - js_m["sumsq"]), sumsq_tol)
            self.assertLessEqual(abs(py_m["maxabs"] - js_m["maxabs"]), maxabs_tol)
            for py_v, js_v in zip(py_m["samples"], js_m["samples"]):
                self.assertLessEqual(abs(py_v - js_v), tol_val)

        for key in ("X", "Y", "S", "NX", "NY", "APANEL"):
            with self.subTest(step="xywake", field=key):
                assert_metrics(py_results["xywake"][key], js_results["xywake"][key])

        with self.subTest(step="qwcalc", field="QINVU1"):
            assert_metrics(py_results["qwcalc"]["QINVU1"], js_results["qwcalc"]["QINVU1"])
        with self.subTest(step="qwcalc", field="QINVU2"):
            assert_metrics(py_results["qwcalc"]["QINVU2"], js_results["qwcalc"]["QINVU2"])
        with self.subTest(step="qiset", field="QINV"):
            assert_metrics(py_results["qiset"]["QINV"], js_results["qiset"]["QINV"])
        with self.subTest(step="qiset", field="QINV_A"):
            assert_metrics(py_results["qiset"]["QINV_A"], js_results["qiset"]["QINV_A"])

        with self.subTest(step="stfind", field="IST"):
            self.assertEqual(py_results["stfind"]["IST"], js_results["stfind"]["IST"])
        with self.subTest(step="stfind", field="SST"):
            self.assertLessEqual(abs(py_results["stfind"]["SST"] - js_results["stfind"]["SST"]), tol_scalar)

        with self.subTest(step="iblpan", field="NBL"):
            self.assertEqual(py_results["iblpan"]["NBL"], js_results["iblpan"]["NBL"])
        with self.subTest(step="iblpan", field="IBLTE"):
            self.assertEqual(py_results["iblpan"]["IBLTE"], js_results["iblpan"]["IBLTE"])
        with self.subTest(step="iblpan", field="IPAN"):
            assert_metrics(py_results["iblpan"]["IPAN"], js_results["iblpan"]["IPAN"], tol_scalar)
        with self.subTest(step="iblpan", field="VTI"):
            assert_metrics(py_results["iblpan"]["VTI"], js_results["iblpan"]["VTI"], tol_scalar)

        with self.subTest(step="xicalc", field="XSSI"):
            assert_metrics(py_results["xicalc"]["XSSI"], js_results["xicalc"]["XSSI"])

        with self.subTest(step="iblsys", field="NSYS"):
            self.assertEqual(py_results["iblsys"]["NSYS"], js_results["iblsys"]["NSYS"])
        with self.subTest(step="iblsys", field="ISYS"):
            assert_metrics(py_results["iblsys"]["ISYS"], js_results["iblsys"]["ISYS"], tol_scalar)

        with self.subTest(step="uicalc", field="UINV"):
            assert_metrics(py_results["uicalc"]["UINV"], js_results["uicalc"]["UINV"])
        with self.subTest(step="uicalc", field="UINV_A"):
            assert_metrics(py_results["uicalc"]["UINV_A"], js_results["uicalc"]["UINV_A"])

        with self.subTest(step="qvfue", field="QVIS"):
            assert_metrics(py_results["qvfue"]["QVIS"], js_results["qvfue"]["QVIS"])
        with self.subTest(step="gamqv", field="GAM"):
            assert_metrics(py_results["gamqv"]["GAM"], js_results["gamqv"]["GAM"])
        with self.subTest(step="gamqv", field="GAM_A"):
            assert_metrics(py_results["gamqv"]["GAM_A"], js_results["gamqv"]["GAM_A"])

        with self.subTest(step="stmove", field="IST"):
            self.assertEqual(py_results["stmove"]["IST"], js_results["stmove"]["IST"])
        with self.subTest(step="stmove", field="XSSI"):
            assert_metrics(py_results["stmove"]["XSSI"], js_results["stmove"]["XSSI"])

        with self.subTest(step="setbl", field="REYBL"):
            if js_results["setbl"]["REYBL"] is None:
                self.fail("JS setbl scalars contain non-finite values")
            self.assertLessEqual(abs(py_results["setbl"]["REYBL"] - js_results["setbl"]["REYBL"]), tol_scalar)
        with self.subTest(step="setbl", field="RSTBL"):
            if js_results["setbl"]["RSTBL"] is None:
                self.fail("JS setbl scalars contain non-finite values")
            self.assertLessEqual(abs(py_results["setbl"]["RSTBL"] - js_results["setbl"]["RSTBL"]), tol_scalar)
        with self.subTest(step="setbl", field="HSTINV"):
            if js_results["setbl"]["HSTINV"] is None:
                self.fail("JS setbl scalars contain non-finite values")
            self.assertLessEqual(abs(py_results["setbl"]["HSTINV"] - js_results["setbl"]["HSTINV"]), tol_scalar)
        with self.subTest(step="setbl", field="VA"):
            assert_metrics(py_results["setbl"]["VA"], js_results["setbl"]["VA"], tol_scalar)
        with self.subTest(step="setbl", field="VB"):
            assert_metrics(py_results["setbl"]["VB"], js_results["setbl"]["VB"], tol_scalar)
        with self.subTest(step="setbl", field="VDEL"):
            assert_metrics(py_results["setbl"]["VDEL"], js_results["setbl"]["VDEL"], tol_scalar)
        with self.subTest(step="setbl", field="VM"):
            assert_metrics(py_results["setbl"]["VM"], js_results["setbl"]["VM"], tol_scalar, rel_tol=1.0e-12)

        with self.subTest(step="blsolv", field="VDEL"):
            assert_metrics(py_results["blsolv"]["VDEL"], js_results["blsolv"]["VDEL"], tol_scalar)

        with self.subTest(step="update", field="CTAU"):
            assert_metrics(py_results["update"]["CTAU"], js_results["update"]["CTAU"], tol_scalar)
        with self.subTest(step="update", field="THET"):
            assert_metrics(py_results["update"]["THET"], js_results["update"]["THET"], tol_scalar)
        with self.subTest(step="update", field="DSTR"):
            assert_metrics(py_results["update"]["DSTR"], js_results["update"]["DSTR"], tol_scalar)
        with self.subTest(step="update", field="UEDG"):
            assert_metrics(py_results["update"]["UEDG"], js_results["update"]["UEDG"], tol_scalar)
        with self.subTest(step="update", field="CL"):
            if js_results["update"]["CL"] is None:
                self.fail("JS update scalars contain non-finite values")
            self.assertLessEqual(abs(py_results["update"]["CL"] - js_results["update"]["CL"]), tol_scalar)
        with self.subTest(step="update", field="RMSBL"):
            if js_results["update"]["RMSBL"] is None:
                self.fail("JS update scalars contain non-finite values")
            self.assertLessEqual(abs(py_results["update"]["RMSBL"] - js_results["update"]["RMSBL"]), tol_scalar)
        with self.subTest(step="update", field="RMXBL"):
            if js_results["update"]["RMXBL"] is None:
                self.fail("JS update scalars contain non-finite values")
            self.assertLessEqual(abs(py_results["update"]["RMXBL"] - js_results["update"]["RMXBL"]), tol_scalar)
        with self.subTest(step="update", field="RLX"):
            if js_results["update"]["RLX"] is None:
                self.fail("JS update scalars contain non-finite values")
            self.assertLessEqual(abs(py_results["update"]["RLX"] - js_results["update"]["RLX"]), tol_scalar)

        with self.subTest(step="cpcalc", field="CPI"):
            assert_metrics(py_results["cpcalc"]["CPI"], js_results["cpcalc"]["CPI"], tol_scalar)
        with self.subTest(step="clcalc", field="CL"):
            if js_results["clcalc"]["CL"] is None:
                self.fail("JS clcalc scalars contain non-finite values")
            self.assertLessEqual(abs(py_results["clcalc"]["CL"] - js_results["clcalc"]["CL"]), tol_scalar)
        with self.subTest(step="clcalc", field="CM"):
            if js_results["clcalc"]["CM"] is None:
                self.fail("JS clcalc scalars contain non-finite values")
            self.assertLessEqual(abs(py_results["clcalc"]["CM"] - js_results["clcalc"]["CM"]), tol_scalar)
        with self.subTest(step="clcalc", field="CDP"):
            if js_results["clcalc"]["CDP"] is None:
                self.fail("JS clcalc scalars contain non-finite values")
            self.assertLessEqual(abs(py_results["clcalc"]["CDP"] - js_results["clcalc"]["CDP"]), tol_scalar)
        with self.subTest(step="cdcalc", field="CD"):
            if js_results["cdcalc"]["CD"] is None:
                self.fail("JS cdcalc scalars contain non-finite values")
            self.assertLessEqual(abs(py_results["cdcalc"]["CD"] - js_results["cdcalc"]["CD"]), tol_scalar)
        with self.subTest(step="cdcalc", field="CDF"):
            if js_results["cdcalc"]["CDF"] is None:
                self.fail("JS cdcalc scalars contain non-finite values")
            self.assertLessEqual(abs(py_results["cdcalc"]["CDF"] - js_results["cdcalc"]["CDF"]), tol_scalar)


if __name__ == "__main__":
    unittest.main()
