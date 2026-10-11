"""Strict binary64 oracles for numerical changes in XFOIL 6.996."""
import json
import contextlib
import io
import math
import pathlib
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from python import naca, spline, xgeom, xblsys, xpanel, xfoil
from python.xbl import XFoilState
from python.tests.fortran_reference import SOURCE_DIR, reference_source


def run(command, **kwargs):
    result = subprocess.run(command, capture_output=True, text=True, **kwargs)
    if result.returncode:
        raise RuntimeError(f"{command}\n{result.stdout}\n{result.stderr}")
    return result.stdout


class TestUpstreamNumericalChanges(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('gfortran') or not shutil.which('node'):
            raise unittest.SkipTest('gfortran and node are required for upstream parity')
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        temp = pathlib.Path(cls.tmp.name)
        cls.fortran = temp / 'reference'
        run(['gfortran', '-O2', '-ffp-contract=off', '-fdefault-real-8', '-fdefault-double-8',
             '-ffunction-sections', '-fdata-sections', '-Wl,--gc-sections',
             '-I', str(SOURCE_DIR), str(ROOT / 'python/tests/upstream_changes_driver.f90'),
             str(ROOT / 'python/tests/upstream_changes_common.f'),
             str(SOURCE_DIR / 'naca.f'), str(SOURCE_DIR / 'spline.f'),
             str(SOURCE_DIR / 'xpanel.f'), str(SOURCE_DIR / 'xutils.f'),
             str(reference_source(temp, 'xfoil.f', ('TECALC', 'MRCL', 'COMSET'))),
             str(reference_source(temp, 'xgeom.f', ('LEFIND',))),
             str(reference_source(temp, 'xblsys.f', ('CFT',))), '-o', str(cls.fortran)])
        cls.ports = {}
        sources = [ROOT / 'python/tests/upstream_changes_driver.cpp'] + [
            ROOT / 'wasm' / name for name in ['naca.cpp', 'spline.cpp', 'xgeom.cpp', 'xutils.cpp',
                'xfoil.cpp', 'xbl.cpp', 'xpanel.cpp']]
        for compiler, name, suffix in [('c++', 'C++', ''), ('em++', 'WASM', '.cjs')]:
            if not shutil.which(compiler):
                continue
            target = temp / ('port' + suffix)
            flags = ['-sENVIRONMENT=node', '-sALLOW_MEMORY_GROWTH=1'] if suffix else []
            run([compiler, '-std=c++17', '-O2', '-ffp-contract=off', '-ffunction-sections',
                 '-fdata-sections', '-Wl,--gc-sections', *flags, *map(str, sources), '-o', str(target)])
            cls.ports[name] = ['node', str(target)] if suffix else [str(target)]

    def compare(self, payload, text, python_values):
        expected = [float(value) for value in run([str(self.fortran)], input=text).split('REFERENCE_RESULT')[-1].split()]
        values = {'Python': python_values, 'JavaScript': json.loads(run(
            ['node', str(ROOT / 'python/tests/compare_upstream_changes.mjs')], input=json.dumps(payload)))}
        for name, command in self.ports.items():
            values[name] = [float(value) for value in run(command, input=text).split('REFERENCE_RESULT')[-1].split()]
        for name, actual in values.items():
            with self.subTest(port=name, payload=payload):
                self.assertEqual(len(actual), len(expected))
                for index, (a, b) in enumerate(zip(actual, expected)):
                    self.assertTrue(math.isfinite(a) and math.isfinite(b), (index, a, b))
                    self.assertLessEqual(abs(a - b), 2e-12 * max(1.0, abs(b)), (index, a, b))
        return expected

    def test_cft_limiter_derivatives(self):
        for rt in [0.5, 5.0, 20.0, math.exp(3) * (1 - 1e-8), math.exp(3) * (1 + 1e-8), 25.0, 2000.0]:
            for msq in [0.0, 0.49]:
                args = [2.3, rt, msq, 1.0]
                expected = self.compare({'mode': 1, 'args': args}, '1\n' + ' '.join(map(str, args)), xblsys.cft(*args))
                if math.log(rt / math.sqrt(1 + 0.2 * msq)) < 3:
                    self.assertEqual(expected[2], 0.0)

    def test_naca_standard_and_sharp(self):
        for ides in [12, 2412, 4415, 21010, 23012, 25018]:
            for sharp in [False, True]:
                n = 200
                xx, yt, yc = ([0.0] * (n + 1) for _ in range(3))
                x, y = ([0.0] * (2 * n + 1) for _ in range(2))
                generate = naca.naca4 if ides <= 9999 else naca.naca5
                nb, _ = generate(ides, xx, yt, yc, n, x, y, sharp)
                values = [v for pair in zip(x[1:nb + 1], y[1:nb + 1]) for v in pair]
                self.compare({'mode': 2, 'args': [ides, n, sharp]}, f'2\n{ides} {n} {"T" if sharp else "F"}\n', values)
                if sharp:
                    self.assertLess(abs(y[1] - y[nb]), 1e-15)

    def test_segmented_endpoint_slopes(self):
        for s, x in [
            ([0, 1, 2, 2, 3, 4], [0, 0.5, 1, 1, 0.5, 0]),
            ([0, 1, 1, 2, 3, 3, 4], [0, 1, 1, 1.3, 2, 2, 3]),
            ([0, 1, 2, 3], [1, 0, 1, 0]),
            (list(range(1201)), [math.sin(i / 200) for i in range(1201)]),
        ]:
            n = len(s)
            first, last = 0.7, -0.3
            xs = [0.0] * (n + 1)
            spline.segspld([0.0] + x, xs, [0.0] + s, n, first, last)
            text = f'3\n{n} {first} {last}\n' + '\n'.join(f'{a} {b}' for a, b in zip(s, x))
            self.compare({'mode': 3, 's': s, 'x': x, 'slopes': [first, last]}, text, xs[1:])

    def test_leading_edge_search(self):
        for x, y in [
            ([1, .8, .5, 0, .5, .8, 1], [0, .06, .1, 0, -.05, -.04, 0]),
            ([1, .8, .5, .51, .2, 0, .3, .8, 1], [0, .03, .04, .2, .1, 0, -.07, -.03, 0]),
        ]:
            n = len(x)
            x1, y1 = [0.0] + x, [0.0] + y
            s, xp, yp = ([0.0] * (n + 1) for _ in range(3))
            spline.scalc(x1, y1, s, n)
            spline.segspl(x1, xp, s, n)
            spline.segspl(y1, yp, s, n)
            value = xgeom.lefind(x1, xp, y1, yp, s, n)
            self.compare({'mode': 4, 'x': x, 'y': y}, f'4\n{n}\n' + '\n'.join(f'{a} {b}' for a, b in zip(x, y)), [value])

    def test_native_and_wasm_coverage(self):
        for name in ['C++', 'WASM']:
            with self.subTest(port=name):
                if name not in self.ports:
                    self.skipTest(f'{name} compiler is required to cover this port')

    def test_mach_reynolds_scaling_and_limits(self):
        for matyp, retyp, mach, reynolds, cl in [
            (1, 1, 0.0, 1e6, 0.0), (1, 1, .3, 3e6, .5),
            (2, 2, .2, 1e6, .25), (3, 3, .2, 1e6, .25),
            (2, 2, .2, 1e6, -.1), (2, 3, .1, 1e6, .0001),
            (0, 4, .1, 1e6, .5), (1, 1, 1.1, 0.0, 0.0),
        ]:
            ctx = XFoilState()
            ctx.MATYP, ctx.RETYP, ctx.MINF1, ctx.REINF1 = matyp, retyp, mach, reynolds
            with contextlib.redirect_stdout(io.StringIO()):
                xfoil.mrcl(ctx, cl)
                xfoil.comset(ctx)
            values = [getattr(ctx, key) for key in ['MATYP', 'RETYP', 'MINF', 'REINF',
                'MINF_CL', 'REINF_CL', 'TKLAM', 'TKL_MSQ', 'CPSTAR', 'QSTAR']]
            args = [matyp, retyp, mach, reynolds, cl]
            self.compare({'mode': 6, 'args': args}, '6\n' + ' '.join(map(str, args)), values)

    def test_trailing_edge_classification_rotated(self):
        for angle in [0.0, 0.7, math.pi / 2]:
            for gap in [0.00008, 0.00012]:
                xy = [(1, gap / 2), (.5, .1), (0, 0), (.5, -.1), (1, -gap / 2)]
                x = [a * math.cos(angle) - b * math.sin(angle) for a, b in xy]
                y = [a * math.sin(angle) + b * math.cos(angle) for a, b in xy]
                ctx = XFoilState()
                ctx.N = len(x)
                ctx.X[1:ctx.N + 1], ctx.Y[1:ctx.N + 1] = x, y
                spline.scalc(ctx.X, ctx.Y, ctx.S, ctx.N)
                spline.segspl(ctx.X, ctx.XP, ctx.S, ctx.N)
                spline.segspl(ctx.Y, ctx.YP, ctx.S, ctx.N)
                sle = xgeom.lefind(ctx.X, ctx.XP, ctx.Y, ctx.YP, ctx.S, ctx.N)
                dx = .5 * (x[0] + x[-1]) - spline.seval(sle, ctx.X, ctx.XP, ctx.S, ctx.N)
                dy = .5 * (y[0] + y[-1]) - spline.seval(sle, ctx.Y, ctx.YP, ctx.S, ctx.N)
                ctx.CHORD = math.sqrt(dx * dx + dy * dy)
                xfoil.tecalc(ctx)
                values = [ctx.ANTE, ctx.ASTE, ctx.DSTE, int(ctx.SHARP)]
                self.compare({'mode': 7, 'x': x, 'y': y}, f'7\n{ctx.N}\n' +
                    '\n'.join(f'{a} {b}' for a, b in zip(x, y)), values)

    def test_panel_influence_sharp_blunt_and_reflected(self):
        for sharp in [False, True]:
            for reflected in [False, True]:
                for ip, xi, yi in [(0, 0.3, 0.2), (1, 1.0, 0.0)]:
                    x = [1.0, 0.5, 0.0, 0.5, 0.9, 1.0]
                    y = [0.0, 0.1, 0.0, -0.1, -0.02, 0.0 if sharp else -0.01]
                    ctx = XFoilState()
                    ctx.N = len(x)
                    ctx.CHORD = 1.0
                    ctx.ALFA = 0.2
                    ctx.LIMAGE = reflected
                    ctx.YIMAGE = -0.75
                    for j, (xx, yy) in enumerate(zip(x, y), 1):
                        ctx.X[j], ctx.Y[j] = xx, yy
                        ctx.GAM[j] = 0.1 * (j - 3)
                        ctx.GAMU[j][1], ctx.GAMU[j][2] = 0.2 * j, -0.1 * j
                        for key, factor in [('SIG', 0.005), ('QF0', .03), ('QF1', -.02), ('QF2', .01), ('QF3', -.04)]:
                            getattr(ctx, key)[j] = factor * j
                    spline.scalc(ctx.X, ctx.Y, ctx.S, ctx.N)
                    spline.segspl(ctx.X, ctx.XP, ctx.S, ctx.N)
                    spline.segspl(ctx.Y, ctx.YP, ctx.S, ctx.N)
                    xpanel.ncalc(ctx.X, ctx.Y, ctx.S, ctx.N, ctx.NX, ctx.NY)
                    xfoil.tecalc(ctx)
                    ctx.SHARP = sharp
                    xpanel.apcalc(ctx)
                    values = list(xpanel.psilin(ctx, ip, xi, yi, 0.6, 0.8, True, True)[:4])
                    values += [getattr(ctx, key) for key in ['Z_QINF', 'Z_ALFA', 'Z_QDOF0', 'Z_QDOF1', 'Z_QDOF2', 'Z_QDOF3']]
                    values += [getattr(ctx, key)[j] for j in range(1, ctx.N + 1)
                               for key in ['DZDG', 'DZDN', 'DQDG', 'DZDM', 'DQDM']]
                    args = [ip, int(sharp), int(reflected), .2, -.75, xi, yi, .6, .8]
                    text = '5\n' + ' '.join(map(str, [ctx.N] + args)) + '\n'
                    text += '\n'.join(f'{a} {b}' for a, b in zip(x, y))
                    self.compare({'mode': 5, 'x': x, 'y': y, 'args': args}, text, values)
