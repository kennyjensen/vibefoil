"""Complete operating points against the checked-in double-precision XFOIL."""
import contextlib
import io
import json
import math
import pathlib
import shutil
import subprocess
import tempfile
import unittest

from python.tests.fortran_reference import build_viscal_reference, run_viscal_reference
from python.tests.test_compare_viscal import ROOT, build_viscal_context, build_payload
from python.tests import test_compare_wasm_viscal as wasm_test
from python.xbl import XBlState, blpini
from python.xfoil import clcalc
from python.xoper import viscal

BL_FIELDS = ('UEDG', 'THET', 'DSTR', 'CTAU', 'MASS', 'TSTR', 'DSSR', 'MDEF', 'PDEF', 'KDEF')
SCALARS = ('CL', 'CM', 'CD', 'CDF', 'CLF', 'RMSBL', 'RMXBL', 'RLX', 'IMXBL', 'ISMXBL', 'VMXBL', 'LVCONV', 'N', 'NW')


def result_from_state(ctx):
    out = {key: ctx[key] for key in SCALARS}
    out.update({key: ctx[key] for key in ('XOCTR', 'ITRAN', 'TFORCE')})
    out['ALFA_DEG'] = ctx['ALFA'] / ctx['DTOR']
    out['CDP_SURFACE'] = ctx['CDP']
    out['CDP'] = ctx['CD'] - ctx['CDF']
    out['panel'] = [[ctx[key][i] for key in ('X', 'Y', 'QINV', 'QVIS', 'CPV')]
                    for i in range(1, ctx['N'] + ctx['NW'] + 1)]
    out['bl'] = [[[ctx[key][i][side] for key in BL_FIELDS]
                  for i in range(1, ctx['NBL'][side] + 1)] for side in (1, 2)]
    return out


class TestUpstreamViscal(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran') or not shutil.which('node'):
            raise unittest.SkipTest('gfortran and node are required for upstream parity')
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.reference = build_viscal_reference(cls.tmp.name)

    def assert_result_close(self, actual, expected, language):
        def check(a, b, label):
            if isinstance(b, list):
                self.assertEqual(len(a), len(b), label)
                for i, (ai, bi) in enumerate(zip(a, b)):
                    check(ai, bi, f'{label}[{i}]')
            elif isinstance(b, (bool, str, int)):
                self.assertEqual(a, b, label)
            else:
                self.assertIsNotNone(a, label)
                self.assertTrue(math.isfinite(a) and math.isfinite(b), label)
                tol = 1e-8 if label.startswith(('bl[', 'panel[')) else 1e-9
                # V8 and libm round log/atan2 differently. Panel cancellation
                # and the coupled Newton solve amplify those last-bit changes.
                # Replaying Python's elementary results in JS makes the input
                # influence matrices bit-identical; retain tighter native bounds.
                if language == 'JavaScript':
                    tol = 1e-7 if label.startswith(('bl[', 'panel[')) else 1e-8
                self.assertLessEqual(abs(a - b), tol * max(1.0, abs(b)), f'{label}: {a} vs {b}')
        for key, value in expected.items():
            check(actual[key], value, key)

    def test_operating_points(self):
        cases = [
            # Current default panel settings; symmetric and cambered profiles.
            dict(ides=12, minf=0.0, reinf=1e6, alpha_deg=0.0, npan=240, iterations=1),
            dict(ides=2412, minf=0.1, reinf=3e6, alpha_deg=4.0, npan=120, iterations=20),
            dict(ides=23012, minf=0.2, reinf=5e6, alpha_deg=2.0, npan=120, iterations=20),
            # Exercise the new sharp trailing-edge generation through the solver.
            dict(ides=12, minf=0.0, reinf=1e6, alpha_deg=3.0, npan=120, iterations=1, sharp=True),
        ]
        wasm_available = (ROOT / 'wasm/build-wasm/compare_xbl.js').exists()
        for case in cases:
            with self.subTest(case=case):
                reference = run_viscal_reference(self.reference, **case)
                with contextlib.redirect_stdout(io.StringIO()):
                    ctx = build_viscal_context(case['ides'], case['minf'], case['reinf'],
                        math.radians(case['alpha_deg']), 1.0, npan=case['npan'],
                        cterat=0.3, sharp=case.get('sharp', False))
                    for i in range(1, ctx.N + 1):
                        ctx.GAM[i] = math.cos(ctx.ALFA) * ctx.GAMU[i][1] + math.sin(ctx.ALFA) * ctx.GAMU[i][2]
                        ctx.GAM_A[i] = -math.sin(ctx.ALFA) * ctx.GAMU[i][1] + math.cos(ctx.ALFA) * ctx.GAMU[i][2]
                    ctx.CL, ctx.CM, ctx.CDP, ctx.CL_ALF, ctx.CL_MSQ = clcalc(
                        ctx.N, ctx.X, ctx.Y, ctx.GAM, ctx.GAM_A, ctx.ALFA, ctx.MINF, ctx.QINF, ctx.XCMREF, ctx.YCMREF)
                    bl = XBlState()
                    blpini(bl)
                    payload = build_payload(ctx, ctx.MINF, ctx.REINF, ctx.ALFA, 1.0)
                    payload['iterations'] = case['iterations']
                    if wasm_available:
                        wasm = wasm_test.TestWasmViscalParity.run_wasm_viscal(self, ctx, bl, case['iterations'])
                    viscal(ctx, bl, case['iterations'])
                js = subprocess.run(['node', str(pathlib.Path(__file__).with_name('compare_viscal_summary.mjs'))],
                    input=json.dumps(payload), text=True, capture_output=True, check=True)
                outputs = {'Python': result_from_state(ctx.__dict__), 'JavaScript': json.loads(js.stdout)['results']}
                if wasm_available:
                    outputs['WASM'] = result_from_state(wasm['ctx'])
                for language, result in outputs.items():
                    with self.subTest(language=language):
                        self.assert_result_close(result, reference, language)

    def test_wasm_coverage(self):
        if not (ROOT / 'wasm/build-wasm/compare_xbl.js').exists():
            self.skipTest('Build wasm/compare_xbl before running the full three-port audit')
