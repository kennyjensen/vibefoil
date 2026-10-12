"""Exercise the accuracy gate itself, including deliberate numerical regressions."""
import contextlib
import copy
import gzip
import io
import json
import math
import pathlib
import re
import shutil
import tempfile
import types
import subprocess
import unittest
from unittest.mock import patch

import run as canonical
from python.tests.fortran_reference import extract_routines, set_residual_tolerance


class ComparatorTests(unittest.TestCase):
    def setUp(self):
        self.manifest = canonical.load_manifest()
        self.expected = {
            'N': 1, 'NW': 0, 'converged': True,
            'coefficients': {'CL': 0.8, 'CM': -0.05, 'CD': 0.007},
            'panel': [[0.4, 0.06, 1.5, 1.4, -1.1, -0.9]],
            'transition': {'x': [0.4, 0.8], 'index': [50, 80], 'forced': [False, True]},
            'bl': [[[1, 0.1, 1.2, 0.0004, 0.0007, 0.03, 0.0008, 0.0009, 0.001, 0.002, 0.03, 0.04]]],
        }

    def compare(self, actual, expected=None):
        return canonical.compare(actual, self.expected if expected is None else expected,
                                 {}, self.manifest['tolerances'])

    def test_exact_and_last_bit_roundoff_pass(self):
        actual = copy.deepcopy(self.expected)
        self.assertTrue(self.compare(actual)['passed'])
        actual['coefficients']['CL'] = math.nextafter(0.8, math.inf)
        actual['bl'][0][0][3] = math.nextafter(0.0004, -math.inf)
        self.assertTrue(self.compare(actual)['passed'])

    def test_field_specific_bounds_catch_real_regressions(self):
        for path, group in [
            (('coefficients', 'CL'), 'coefficient'),
            (('coefficients', 'CD'), 'drag'),
            (('panel', 0, 0), 'geometry'),
            (('panel', 0, 5), 'distribution'),
            (('bl', 0, 0, 3), 'thickness'),
            (('transition', 'x', 0), 'transition'),
        ]:
            for multiple, passed in [(0.5, True), (2.0, False)]:
                with self.subTest(path=path, multiple=multiple):
                    actual = copy.deepcopy(self.expected)
                    parent = actual
                    for key in path[:-1]:
                        parent = parent[key]
                    tol = self.manifest['tolerances'][group]
                    parent[path[-1]] += multiple * (tol['atol'] + tol['rtol'] * abs(parent[path[-1]]))
                    result = self.compare(actual)
                    self.assertEqual(result['passed'], passed)
                    if not passed:
                        self.assertEqual(result['failures'][0]['path'], '.'.join(map(str, path)))

    def test_discrete_values_and_schema_are_exact(self):
        mutations = [
            lambda a: a.update(N=1.0),
            lambda a: a.update(converged=1),
            lambda a: a.update(converged=False),
            lambda a: a['transition']['index'].__setitem__(0, 51),
            lambda a: a['transition']['forced'].__setitem__(0, True),
            lambda a: a['panel'].append(a['panel'][0]),
            lambda a: a['panel'][0].pop(),
            lambda a: a['coefficients'].pop('CD'),
            lambda a: a.update(unexpected=0),
        ]
        for mutate in mutations:
            with self.subTest(mutation=mutate):
                actual = copy.deepcopy(self.expected)
                mutate(actual)
                self.assertFalse(self.compare(actual)['passed'])

    def test_invalid_actual_or_reference_is_never_a_pass(self):
        for value in (math.nan, math.inf, -math.inf, None, True, '0.8'):
            with self.subTest(value=value):
                actual = copy.deepcopy(self.expected)
                actual['coefficients']['CL'] = value
                self.assertFalse(self.compare(actual)['passed'])
        for value in (math.nan, math.inf, -math.inf, None):
            expected = copy.deepcopy(self.expected)
            expected['coefficients']['CL'] = value
            self.assertFalse(self.compare(expected, expected)['passed'])

    def test_json_artifacts_preserve_nonfinite_failure_evidence(self):
        safe = canonical.json_safe({'values': [math.nan, math.inf, -math.inf]})
        self.assertEqual(json.loads(json.dumps(safe, allow_nan=False)),
                         {'values': ['nan', 'inf', '-inf']})


class ManifestTests(unittest.TestCase):
    def test_grid_coverage_and_unique_ids(self):
        manifest = canonical.load_manifest()
        full = canonical.cases(manifest)
        quick = canonical.cases(manifest, True)
        self.assertEqual(len(full), 615)
        self.assertEqual(len({case['id'] for case in full}), 615)
        self.assertEqual(len(quick), 25)
        self.assertEqual(manifest['residual_tolerance'], 1e-8)
        self.assertTrue(all(case['residual_tolerance'] == 1e-8 for case in full))
        self.assertLessEqual({case['id'] for case in quick}, {case['id'] for case in full})
        for foil in manifest['airfoils']:
            subset = [case for case in full if case['foil'] == foil]
            self.assertEqual(sum(case['viscous'] for case in subset), 108)
            for field, values in manifest['viscous_grid'].items():
                self.assertEqual({case[field] for case in subset if case['viscous']}, set(values))

    def test_reference_source_drift_is_rejected(self):
        with patch.object(canonical, 'source_digest', return_value='changed'):
            with self.assertRaisesRegex(ValueError, 'Original XFOIL source changed'):
                canonical.load_manifest()

    def test_coordinate_drift_is_rejected(self):
        original = pathlib.Path.read_bytes

        def read(path):
            return b'changed' if path.name == 'e387.dat' else original(path)

        with patch.object(pathlib.Path, 'read_bytes', read):
            with self.assertRaisesRegex(ValueError, 'Canonical coordinates changed'):
                canonical.load_manifest()


class RunnerTests(unittest.TestCase):
    def test_reference_threshold_override_only_changes_stopping_constant(self):
        source = extract_routines((canonical.SOURCE_DIR / 'xoper.f').read_text(), ('VISCAL',))
        changed = set_residual_tolerance(source, 1e-8)
        differences = [(a, b) for a, b in zip(source.splitlines(), changed.splitlines(), strict=True) if a != b]
        self.assertEqual(differences, [('      DATA EPS1 / 1.0E-4 /',
                                      '      DATA EPS1 / 1.00000000000000002E-08 /')])
        for value in (0, -1, math.nan, math.inf):
            with self.assertRaisesRegex(ValueError, 'finite and positive'):
                set_residual_tolerance(source, value)
        with self.assertRaisesRegex(ValueError, 'exactly one'):
            set_residual_tolerance(changed, 1e-8)

    @unittest.skipUnless(shutil.which('node'), 'Node required for worker convergence test')
    def test_worker_honors_tighter_residual_without_changing_default(self):
        case = next(case for case in canonical.cases(canonical.load_manifest(), True)
                    if case['id'] == 'naca0012-visc-a8-m0.15-re500000-n9')
        def solve(selected):
            output = subprocess.run(['node', str(canonical.HERE / 'worker.mjs'), 'javascript'],
                                    input=json.dumps(selected) + '\n', text=True, capture_output=True,
                                    timeout=30, check=True, cwd=canonical.ROOT)
            value = json.loads(output.stdout)
            self.assertNotIn('error', value)
            iterations = [line.split() for line in value['logs']
                          if re.match(r'^\s*\d+\s', line) and ' at ' in line]
            return value, iterations
        default_case = dict(case)
        del default_case['residual_tolerance']
        default, default_iterations = solve(default_case)
        explicit_default, _ = solve(dict(case, residual_tolerance=1e-4))
        self.assertEqual(default['result'], explicit_default['result'])
        tight, tight_iterations = solve(case)
        self.assertTrue(tight['result']['converged'])
        self.assertLess(float(tight_iterations[-1][1]), 1e-8)
        self.assertGreater(len(tight_iterations), len(default_iterations))
        limited, _ = solve(dict(case, iterations=len(default_iterations)))
        self.assertFalse(limited['result']['converged'])

    def test_sensitivity_perturbs_only_one_input_by_one_ulp(self):
        case = canonical.cases(canonical.load_manifest(), True)[1]
        case = dict(case, alpha=4, cterat=0.3)
        original = copy.deepcopy(case)
        reference = {'diagnostics': {'arithmetic_flags': {'invalid': False, 'overflow': False, 'divide_by_zero': False}},
                     'result': {'N': 0, 'NW': 0, 'converged': True,
                               'coefficients': {'CL': 0.8}, 'panel': []}}
        with patch.object(canonical, 'run_reference', return_value=reference) as run:
            diagnostic = canonical.diagnose_reference_sensitivity(
                pathlib.Path('/stub'), case, reference, canonical.load_manifest()['tolerances'])
        self.assertEqual(case, original)
        self.assertEqual(run.call_count, 10)
        self.assertFalse(diagnostic['proven'])
        self.assertFalse(diagnostic['convergence_changed'])
        self.assertEqual(diagnostic['finite_budget_failures'], 0)
        for call, variant in zip(run.call_args_list, diagnostic['variants']):
            altered = call.args[1]
            changed = [key for key in case if case[key] != altered[key]]
            self.assertEqual(changed, [variant['field']])
            self.assertIn(altered[changed[0]], [math.nextafter(case[changed[0]], -math.inf),
                                              math.nextafter(case[changed[0]], math.inf)])

    def test_sensitivity_records_convergence_changes_without_replacing_reference(self):
        case = canonical.cases(canonical.load_manifest(), True)[0]
        reference = {'diagnostics': {'arithmetic_flags': {'invalid': False, 'overflow': False, 'divide_by_zero': False}},
                     'result': {'N': 0, 'NW': 0, 'converged': True,
                               'coefficients': {'CL': 0.8}, 'panel': []}}
        altered = copy.deepcopy(reference)
        altered['result']['converged'] = False
        altered['result']['coefficients']['CL'] = 1.8
        with patch.object(canonical, 'run_reference', return_value=altered):
            diagnostic = canonical.diagnose_reference_sensitivity(
                pathlib.Path('/stub'), case, reference, canonical.load_manifest()['tolerances'])
        self.assertTrue(diagnostic['convergence_changed'])
        self.assertEqual(len(diagnostic['variants']), 6)  # Mach zero is skipped.
        self.assertEqual(diagnostic['finite_budget_failures'], 6)
        self.assertTrue(reference['result']['converged'])
        self.assertEqual(reference['result']['coefficients']['CL'], 0.8)

    def test_instability_needs_evidence_and_both_guards(self):
        passed, failed = {'passed': True}, {'passed': False}
        self.assertEqual(canonical.classify(passed, None, None, None), 'match')
        for proof, initial, fixed in [(False, passed, passed), (True, failed, passed),
                                      (True, passed, failed), (True, None, passed)]:
            self.assertEqual(canonical.classify(failed, {'proven': proof}, initial, fixed), 'mismatch')
        self.assertEqual(canonical.classify(failed, {'proven': True}, passed, passed), 'reference_unstable')
        self.assertFalse(canonical.proves_instability({'failures': [
            {'reason': 'nonfinite/non-numeric value', 'actual': math.nan, 'reference': math.nan}]}))

    def test_instability_does_not_excuse_geometry_or_operating_errors(self):
        expected = {'N': 1, 'NW': 1, 'panel': [[0.2, 0.04, 1.2], [1.2, 0.1, 1.1]],
                    'operating': {'mach': 0.3, 'reynolds': 1e6}}
        tol = canonical.load_manifest()['tolerances']
        actual = copy.deepcopy(expected)
        actual['panel'][1][1] = 5.0  # Wake trajectory is solution-dependent.
        self.assertTrue(canonical.fixed_inputs_check(actual, expected, {}, tol)['passed'])
        actual['panel'][0][1] += 1e-5
        self.assertFalse(canonical.fixed_inputs_check(actual, expected, {}, tol)['passed'])
        actual = copy.deepcopy(expected)
        actual['operating']['mach'] = 0.0
        self.assertFalse(canonical.fixed_inputs_check(actual, expected, {}, tol)['passed'])

    def test_invalid_reference_needs_clean_early_comparison_and_fixed_inputs(self):
        passed, failed = {'passed': True}, {'passed': False}
        # Even identical final values cannot validate a failed arithmetic solve.
        for final in (passed, failed):
            self.assertEqual(canonical.classify(final, None, passed, passed,
                invalid_reference=True, valid_initial_reference=True), 'reference_invalid')
            for initial, fixed, valid in [(failed, passed, True), (passed, failed, True),
                                           (passed, passed, False), (None, passed, True)]:
                self.assertEqual(canonical.classify(final, None, initial, fixed,
                    invalid_reference=True, valid_initial_reference=valid), 'mismatch')
        self.assertEqual(canonical.classify(failed, {'proven': True}, passed, passed,
            valid_initial_reference=False), 'mismatch')

    def test_reference_arithmetic_flags_required_and_parsed(self):
        case = dict(canonical.cases(canonical.load_manifest(), True)[0], viscous=False)
        body = 'CANONICAL_RESULT\n0 0 F\n' + '0 ' * 11 + '\n0 0 0\nCANONICAL_END\n'
        with tempfile.TemporaryDirectory() as directory:
            for flags in ('F F F', 'T F F', 'F T F', 'F F T'):
                stdout = 'CANONICAL_FLAGS ' + flags + '\n' + body
                with patch.object(canonical.subprocess, 'run', return_value=types.SimpleNamespace(stdout=stdout)):
                    reference = canonical.run_reference(pathlib.Path(directory) / 'oracle', case)
                self.assertEqual(list(reference['diagnostics']['arithmetic_flags'].values()),
                                 [flag == 'T' for flag in flags.split()])
                self.assertEqual(canonical.arithmetic_failed(reference), 'T' in flags)
            for prefix in ('', 'CANONICAL_FLAGS F F\n', 'CANONICAL_FLAGS F F F\n' * 2):
                with patch.object(canonical.subprocess, 'run', return_value=types.SimpleNamespace(stdout=prefix + body)):
                    with self.assertRaisesRegex(ValueError, 'arithmetic flags'):
                        canonical.run_reference(pathlib.Path(directory) / 'oracle', case)

    def test_coordinate_probe_moves_one_interior_coordinate_and_preserves_closure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            path = root / 'foil.dat'
            path.write_text('foil\n1 0\n0.5 0.1\n0 0\n0.5 -0.1\n1 0\n')
            case = dict(canonical.cases(canonical.load_manifest(), True)[0],
                        foil={'coordinates': str(path), 'id': 'test'})
            ref = {'diagnostics': {'arithmetic_flags': {'invalid': False, 'overflow': False, 'divide_by_zero': False}},
                     'result': {'N': 0, 'NW': 0, 'converged': True,
                              'coefficients': {'CL': 0.8}, 'panel': []}}

            def solve(executable, altered, timeout):
                result = copy.deepcopy(ref)
                if altered['foil']['coordinates'] != str(path):
                    result['result']['coefficients']['CL'] = 0.9
                return result

            with patch.object(canonical, 'run_reference', side_effect=solve):
                evidence = canonical.diagnose_reference_sensitivity(
                    pathlib.Path('/stub'), case, ref, canonical.load_manifest()['tolerances'],
                    probe_directory=root / 'probes')
            self.assertTrue(evidence['proven'])
            last = evidence['variants'][-1]
            self.assertEqual(last['field'], 'coordinates')
            before = path.read_text().splitlines()
            after = pathlib.Path(last['input']['coordinates']).read_text().splitlines()
            self.assertEqual(before[0:2], after[0:2])
            self.assertEqual(before[3:], after[3:])
            old_values = list(map(float, before[2].split()))
            new_values = list(map(float, after[2].split()))
            self.assertEqual(new_values[1], old_values[1])
            self.assertEqual(new_values[0], math.nextafter(old_values[0], -math.inf))

    def test_load_uses_short_name_in_long_workspace_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            source = root / ('long-coordinate-name-' * 6 + '.dat')
            source.write_bytes(b'foil\n1 0\n0 0\n1 0\n')
            case = dict(canonical.cases(canonical.load_manifest(), True)[0],
                        foil={'coordinates': str(source)}, viscous=False)
            stdout = 'CANONICAL_FLAGS F F F\nCANONICAL_RESULT\n0 0 F\n' + '0 ' * 11 + '\n0 0 0\nCANONICAL_END\n'
            with patch.object(canonical.subprocess, 'run', return_value=types.SimpleNamespace(stdout=stdout)) as run:
                canonical.run_reference(root / 'oracle', case)
            name = run.call_args.kwargs['input'].splitlines()[1]
            self.assertLess(len(name), 80)
            self.assertEqual((root / name).read_bytes(), source.read_bytes())
            self.assertEqual(run.call_args.kwargs['cwd'], root)

    def test_proven_instability_reported_separately_and_strict_mode_fails(self):
        case = canonical.cases(canonical.load_manifest(), True)[0]
        ref = {'diagnostics': {'arithmetic_flags': {'invalid': False, 'overflow': False, 'divide_by_zero': False}},
                     'result': {'N': 0, 'NW': 0, 'converged': True,
                          'coefficients': {'CL': 0.8}, 'panel': []}}
        actual = copy.deepcopy(ref)
        actual['result']['coefficients']['CL'] = 0.9
        for strict in (False, True):
            with self.subTest(strict=strict), tempfile.TemporaryDirectory() as directory:
                with patch.object(canonical.shutil, 'which', return_value='/stub'), \
                     patch.object(canonical.subprocess, 'check_output', return_value='test'), \
                     patch.object(canonical, 'cases', return_value=[case]), \
                     patch.object(canonical, 'build_viscal_reference', return_value=pathlib.Path('/stub')), \
                     patch.object(canonical, 'run_reference', return_value=ref), \
                     patch.object(canonical, 'diagnose_reference_sensitivity', return_value={'proven': True}), \
                     patch.object(canonical, 'Worker') as worker, \
                     contextlib.redirect_stdout(io.StringIO()):
                    worker.return_value.run.side_effect = [actual, ref]
                    args = ['--quick', '--backend', 'auto', '--diagnose-roundoff', '--output', directory]
                    if strict:
                        args.append('--strict-reference')
                    code = canonical.main(args)
                summary = json.loads((pathlib.Path(directory) / 'summary.json').read_text())
                self.assertEqual(code, int(strict))
                self.assertEqual(summary['passed'], not strict)
                self.assertEqual(summary['reference_unstable_cases'], [case['id']])
                self.assertEqual(summary['reference_unstable_comparisons'], 1)
                self.assertEqual(summary['matched_comparisons'], 0)
                self.assertEqual(summary['initial_comparisons'], 1)

    def test_invalid_reference_reported_separately_only_when_unconverged(self):
        case = canonical.cases(canonical.load_manifest(), True)[0]
        clean = {'diagnostics': {'arithmetic_flags': {'invalid': False, 'overflow': False, 'divide_by_zero': False}},
                 'result': {'N': 0, 'NW': 0, 'converged': False,
                            'coefficients': {'CL': 0.8}, 'panel': []}}
        for converged, strict in [(False, False), (False, True), (True, False)]:
            with self.subTest(converged=converged, strict=strict), tempfile.TemporaryDirectory() as directory:
                reference = copy.deepcopy(clean)
                reference['result']['converged'] = converged
                reference['diagnostics']['arithmetic_flags']['invalid'] = True
                actual = copy.deepcopy(reference)
                actual['result']['coefficients']['CL'] = None
                with patch.object(canonical.shutil, 'which', return_value='/stub'), \
                     patch.object(canonical.subprocess, 'check_output', return_value='test'), \
                     patch.object(canonical, 'cases', return_value=[case]), \
                     patch.object(canonical, 'build_viscal_reference', return_value=pathlib.Path('/stub')), \
                     patch.object(canonical, 'run_reference', side_effect=[reference, clean]), \
                     patch.object(canonical, 'diagnose_reference_sensitivity', return_value={'proven': False}), \
                     patch.object(canonical, 'Worker') as worker, \
                     contextlib.redirect_stdout(io.StringIO()):
                    worker.return_value.run.side_effect = [actual, clean]
                    args = ['--quick', '--backend', 'auto', '--output', directory]
                    if strict:
                        args.append('--strict-reference')
                    code = canonical.main(args)
                summary = json.loads((pathlib.Path(directory) / 'summary.json').read_text())
                self.assertEqual(code, int(strict or converged))
                self.assertEqual(summary['reference_invalid_comparisons'], int(not converged))
                self.assertEqual(summary['matched_comparisons'], 0)
                with gzip.open(pathlib.Path(directory) / 'results.jsonl.gz', 'rt') as records:
                    record = json.loads(next(records))
                self.assertFalse(record['checks']['auto']['passed'])
                self.assertEqual(record['classifications']['auto'],
                                 'mismatch' if converged else 'reference_invalid')

    def test_completed_run_reports_success_and_numerical_failure(self):
        case = canonical.cases(canonical.load_manifest(), True)[0]
        reference = {'diagnostics': {'arithmetic_flags': {'invalid': False, 'overflow': False, 'divide_by_zero': False}},
                     'result': {'N': 0, 'NW': 0, 'converged': True,
                               'coefficients': {'CL': 0.8}, 'panel': []}}
        for coefficient, passed in [(0.8, True), (0.9, False)]:
            with self.subTest(coefficient=coefficient), tempfile.TemporaryDirectory() as directory:
                actual = copy.deepcopy(reference)
                actual['result']['coefficients']['CL'] = coefficient
                with patch.object(canonical.shutil, 'which', return_value='/stub'), \
                     patch.object(canonical.subprocess, 'check_output', return_value='test'), \
                     patch.object(canonical, 'cases', return_value=[case]), \
                     patch.object(canonical, 'build_viscal_reference', return_value=pathlib.Path('/stub')), \
                     patch.object(canonical, 'run_reference', return_value=reference), \
                     patch.object(canonical, 'Worker') as worker, \
                     contextlib.redirect_stdout(io.StringIO()):
                    worker.return_value.run.return_value = actual
                    code = canonical.main(['--quick', '--backend', 'auto', '--output', directory])
                summary = json.loads((pathlib.Path(directory) / 'summary.json').read_text())
                self.assertEqual(code, 0 if passed else 1)
                self.assertEqual(summary['passed'], passed)
                self.assertTrue(summary['completed'])
                self.assertFalse(summary['full_gate'])
                self.assertEqual(summary['comparisons'], summary['expected_comparisons'])

    def test_missing_compiler_is_error_not_skip(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as error:
                canonical.main(['--compiler', 'canonical-nonexistent-compiler'])
        self.assertEqual(error.exception.code, 2)

    def test_build_failure_cannot_leave_a_passing_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(canonical.shutil, 'which', return_value='/stub'), \
                 patch.object(canonical.subprocess, 'check_output', return_value='test'), \
                 patch.object(canonical, 'build_viscal_reference', side_effect=RuntimeError('build failed')), \
                 contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                code = canonical.main(['--quick', '--output', directory])
            summary = json.loads((pathlib.Path(directory) / 'summary.json').read_text())
            self.assertEqual(code, 1)
            self.assertFalse(summary['passed'])
            self.assertFalse(summary['completed'])
            self.assertEqual(summary['status'], 'error')
            self.assertIn('build failed', summary['error'])

    def test_runtime_failure_after_success_cannot_pass(self):
        case = canonical.cases(canonical.load_manifest(), True)[0]
        reference = {'diagnostics': {'arithmetic_flags': {'invalid': False, 'overflow': False, 'divide_by_zero': False}},
                     'result': {'N': 0, 'NW': 0, 'converged': True,
                               'coefficients': {'CL': 0.8}, 'panel': []}}
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(canonical.shutil, 'which', return_value='/stub'), \
                 patch.object(canonical.subprocess, 'check_output', return_value='test'), \
                 patch.object(canonical, 'cases', return_value=[case, dict(case, id='second')]), \
                 patch.object(canonical, 'build_viscal_reference', return_value=pathlib.Path('/stub')), \
                 patch.object(canonical, 'run_reference', side_effect=[reference, RuntimeError('crashed')]), \
                 patch.object(canonical, 'Worker') as worker, \
                 contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                worker.return_value.run.return_value = reference
                code = canonical.main(['--quick', '--backend', 'auto', '--output', directory])
            summary = json.loads((pathlib.Path(directory) / 'summary.json').read_text())
            self.assertEqual(code, 1)
            self.assertFalse(summary['passed'])
            self.assertFalse(summary['completed'])
            self.assertEqual(summary['comparisons'], 1)
            self.assertEqual(summary['completed_cases'], 1)
            self.assertEqual(summary['current_case'], 'second')
            worker.return_value.close.assert_called_once()
            with gzip.open(pathlib.Path(directory) / 'results.jsonl.gz', 'rt') as records:
                self.assertEqual(len(list(records)), 1)

    @unittest.skipUnless(shutil.which('node'), 'Node required for worker protocol tests')
    def test_worker_timeout_and_early_exit(self):
        for script, error in [
            ('process.stdin.resume();', TimeoutError),
            ('process.exit(1);', RuntimeError),
        ]:
            with self.subTest(script=script), tempfile.TemporaryDirectory() as directory:
                output = pathlib.Path(directory)
                (output / 'worker.mjs').write_text(script)
                with patch.object(canonical, 'HERE', output):
                    worker = canonical.Worker('node', 'auto', output, timeout=0.5)
                try:
                    with self.assertRaises(error):
                        worker.run({'id': 'timeout', 'foil': {}})
                finally:
                    worker.close()
                self.assertIsNotNone(worker.proc.returncode)


if __name__ == '__main__':
    unittest.main()
