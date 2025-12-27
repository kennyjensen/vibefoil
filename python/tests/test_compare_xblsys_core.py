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
from python.xbl import XFoilState, XBlState, blpini, iblsys
from python.xfoil import comset, naca
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
from python.xsolve import gauss as gauss_base
from python.xbl import setbl
from python.xblsys import blsys, blvar, blmid, trchek, tesys
from python.xblcom import sync_vars_to_com


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


def prepare_bl_state():
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
    setbl(ctx, bl)

    return ctx, bl


def serialize_bl(bl):
    data = {}
    for key, val in bl.__dict__.items():
        if isinstance(val, list):
            data[key] = val
        elif isinstance(val, (int, float, bool)):
            data[key] = val
    return data


def pick_fields(bl, fields):
    return {name: getattr(bl, name) for name in fields}


def pick_array(mat, rows, cols):
    return [[mat[i][j] for j in range(1, cols + 1)] for i in range(1, rows + 1)]


def pick_vector(vec, count):
    return [vec[i] for i in range(1, count + 1)]


class TestXblsysCoreParity(unittest.TestCase):
    def test_xblsys_core(self):
        if not shutil.which("node"):
            self.skipTest("node is required for JS/Python parity checks")

        _, bl = prepare_bl_state()
        trchek_override = {}
        if abs(bl.X2 - bl.X1) < 1.0e-12:
            trchek_override["X1"] = 0.0
            trchek_override["X2"] = 1.0

        trchek_fields = [
            "AMPL2",
            "XT",
            "XT_A1",
            "XT_A2",
            "XT_T1",
            "XT_D1",
            "XT_U1",
            "XT_X1",
            "XT_T2",
            "XT_D2",
            "XT_U2",
            "XT_X2",
            "XT_MS",
            "XT_RE",
            "XT_XF",
        ]

        blvar_fields = [
            "HK2",
            "HK2_U2",
            "HK2_T2",
            "HK2_D2",
            "HK2_MS",
            "HC2",
            "HC2_U2",
            "HC2_T2",
            "HC2_D2",
            "HC2_MS",
            "HS2",
            "HS2_U2",
            "HS2_T2",
            "HS2_D2",
            "HS2_MS",
            "HS2_RE",
            "US2",
            "US2_U2",
            "US2_T2",
            "US2_D2",
            "US2_MS",
            "US2_RE",
            "CF2",
            "CF2_U2",
            "CF2_T2",
            "CF2_D2",
            "CF2_MS",
            "CF2_RE",
            "DI2",
            "DI2_U2",
            "DI2_T2",
            "DI2_D2",
            "DI2_S2",
            "DI2_MS",
            "DI2_RE",
            "CQ2",
            "CQ2_U2",
            "CQ2_T2",
            "CQ2_D2",
            "CQ2_MS",
            "CQ2_RE",
            "DE2",
            "DE2_U2",
            "DE2_T2",
            "DE2_D2",
            "DE2_MS",
        ]

        blmid_fields = [
            "CFM",
            "CFM_HKA",
            "CFM_RTA",
            "CFM_MA",
            "CFM_MS",
            "CFM_RE",
            "CFM_U1",
            "CFM_T1",
            "CFM_D1",
            "CFM_U2",
            "CFM_T2",
            "CFM_D2",
        ]

        cte = bl.S2 + 0.01
        tte = bl.T2 + 0.002
        dte = bl.D2 + bl.DW2 + 0.001

        payload = {
            "base": serialize_bl(bl),
            "config": {
                "trchek": {"fields": trchek_fields, "override": trchek_override},
                "blvar": {"ityps": [1, 2], "fields": blvar_fields},
                "blmid": {"ityps": [1, 2], "fields": blmid_fields},
                "blsys": True,
                "tesys": {"cte": cte, "tte": tte, "dte": dte},
            },
        }

        script = pathlib.Path(__file__).with_name("compare_xblsys_core.mjs")
        proc = subprocess.run(
            ["node", str(script)],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            check=True,
        )
        js_results = json.loads(proc.stdout)["results"]

        tol = 1.0e-9

        bl_tr = copy.deepcopy(bl)
        if "X1" in trchek_override:
            bl_tr.X1 = trchek_override["X1"]
        if "X2" in trchek_override:
            bl_tr.X2 = trchek_override["X2"]
        if trchek_override:
            sync_vars_to_com(bl_tr, 1)
            sync_vars_to_com(bl_tr, 2)
        trchek(bl_tr)
        py_tr = pick_fields(bl_tr, trchek_fields)
        for name in trchek_fields:
            self.assertLessEqual(abs(py_tr[name] - js_results["trchek"][name]), tol, msg=f"trchek {name}")

        for ityp in (1, 2):
            bl_v = copy.deepcopy(bl)
            blvar(bl_v, ityp)
            py_v = pick_fields(bl_v, blvar_fields)
            js_v = js_results["blvar"][str(ityp)]
            for name in blvar_fields:
                self.assertLessEqual(abs(py_v[name] - js_v[name]), tol, msg=f"blvar {ityp} {name}")

        for ityp in (1, 2):
            bl_m = copy.deepcopy(bl)
            blmid(bl_m, ityp)
            py_m = pick_fields(bl_m, blmid_fields)
            js_m = js_results["blmid"][str(ityp)]
            for name in blmid_fields:
                self.assertLessEqual(abs(py_m[name] - js_m[name]), tol, msg=f"blmid {ityp} {name}")

        bl_s = copy.deepcopy(bl)
        blsys(bl_s)
        py_blsys = {
            "VS1": pick_array(bl_s.VS1, 4, 5),
            "VS2": pick_array(bl_s.VS2, 4, 5),
            "VSREZ": pick_vector(bl_s.VSREZ, 4),
            "VSM": pick_vector(bl_s.VSM, 4),
            "VSR": pick_vector(bl_s.VSR, 4),
            "VSX": pick_vector(bl_s.VSX, 4),
        }
        js_blsys = js_results["blsys"]
        for name in ("VS1", "VS2"):
            for i in range(4):
                for j in range(5):
                    self.assertLessEqual(
                        abs(py_blsys[name][i][j] - js_blsys[name][i][j]),
                        tol,
                        msg=f"blsys {name}[{i+1}][{j+1}]",
                    )
        for name in ("VSREZ", "VSM", "VSR", "VSX"):
            for i in range(4):
                self.assertLessEqual(
                    abs(py_blsys[name][i] - js_blsys[name][i]),
                    tol,
                    msg=f"blsys {name}[{i+1}]",
                )

        bl_t = copy.deepcopy(bl)
        tesys(bl_t, cte, tte, dte)
        py_tesys = {
            "VS1": pick_array(bl_t.VS1, 4, 5),
            "VS2": pick_array(bl_t.VS2, 4, 5),
            "VSREZ": pick_vector(bl_t.VSREZ, 4),
            "VSM": pick_vector(bl_t.VSM, 4),
            "VSR": pick_vector(bl_t.VSR, 4),
            "VSX": pick_vector(bl_t.VSX, 4),
        }
        js_tesys = js_results["tesys"]
        for name in ("VS1", "VS2"):
            for i in range(4):
                for j in range(5):
                    self.assertLessEqual(
                        abs(py_tesys[name][i][j] - js_tesys[name][i][j]),
                        tol,
                        msg=f"tesys {name}[{i+1}][{j+1}]",
                    )
        for name in ("VSREZ", "VSM", "VSR", "VSX"):
            for i in range(4):
                self.assertLessEqual(
                    abs(py_tesys[name][i] - js_tesys[name][i]),
                    tol,
                    msg=f"tesys {name}[{i+1}]",
                )


if __name__ == "__main__":
    unittest.main()
