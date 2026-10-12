#!/usr/bin/env python3
"""Compile original XFOIL and compare the website worker over the canonical grid."""
import argparse
import csv
import datetime
import fnmatch
import gzip
import hashlib
import itertools
import json
import math
import pathlib
import queue
import re
import shutil
import subprocess
import sys
import threading
import time
import traceback

ROOT = pathlib.Path(__file__).resolve().parents[2]
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from python.tests.fortran_reference import SOURCE_DIR, REFERENCE_FLAGS, build_viscal_reference

BL_FIELDS = ['IPAN', 'XSSI', 'UEDG', 'THET', 'DSTR', 'CTAU', 'MASS', 'TSTR', 'TAU', 'DIS', 'CTQ', 'DELT']
PANEL_FIELDS = ['X', 'Y', 'QINV', 'QVIS', 'CPI', 'CPV']


def source_digest():
    digest = hashlib.sha256()
    for path in sorted(SOURCE_DIR.iterdir()):
        if path.suffix.lower() in ('.f', '.inc'):
            digest.update(path.name.encode() + b'\0' + path.read_bytes() + b'\0')
    return digest.hexdigest()


def load_manifest():
    manifest = json.loads((HERE / 'cases.json').read_text())
    if not math.isfinite(manifest['residual_tolerance']) or manifest['residual_tolerance'] <= 0:
        raise ValueError('Residual tolerance must be finite and positive')
    if manifest['reference']['source_sha256'] != source_digest():
        raise ValueError('Original XFOIL source changed; review the reference revision and manifest fingerprint')
    for foil in manifest['airfoils']:
        if 'coordinates' in foil:
            path = ROOT / foil['coordinates']
            if hashlib.sha256(path.read_bytes()).hexdigest() != foil['sha256']:
                raise ValueError(f'Canonical coordinates changed: {path}')
    return manifest


def cases(manifest, quick=False):
    result = []
    grid = manifest['viscous_grid']
    for foil in manifest['airfoils']:
        for alpha, mach, re, ncrit in itertools.product(*(grid[key] for key in ('alpha', 'mach', 'reynolds', 'ncrit'))):
            # A fixed diagnostic subset; the default gate always runs the full grid.
            if quick and (alpha, mach, re, ncrit) not in [
                (-4, 0, 500000, 5), (0, 0.15, 1000000, 9),
                (4, 0.3, 5000000, 12), (8, 0.15, 500000, 9),
            ]:
                continue
            result.append(dict(foil=foil, alpha=alpha, mach=mach, reynolds=re,
                               ncrit=ncrit, viscous=True, iterations=manifest['iterations']))
        for alpha, mach in itertools.product(*(manifest['inviscid_grid'][key] for key in ('alpha', 'mach'))):
            if quick and (alpha, mach) != (4, 0.3):
                continue
            result.append(dict(foil=foil, alpha=alpha, mach=mach, reynolds=1000000,
                               ncrit=9, viscous=False, iterations=0))
    for case in result:
        case['residual_tolerance'] = manifest['residual_tolerance']
        case['id'] = (f"{case['foil']['id']}-{'visc' if case['viscous'] else 'inv'}"
                      f"-a{case['alpha']:g}-m{case['mach']:g}"
                      f"-re{case['reynolds']:g}-n{case['ncrit']:g}")
    return result


def run_reference(executable, case, timeout=120):
    foil = case['foil']
    args = [foil.get('naca', 0), case['mach'], case['reynolds'], case['alpha'],
            case['ncrit'], case['iterations'], int(case['viscous']), case.get('cterat', 0.3)]
    coordinate = '-'
    if 'coordinates' in foil:
        data = (ROOT / foil['coordinates']).read_bytes()
        # Upstream LOAD truncates long file names. Use a short, content-addressed
        # copy in the oracle's working directory, preserving the input bytes.
        coordinate = 'foil-' + hashlib.sha256(data).hexdigest()[:20] + '.dat'
        local = executable.parent / coordinate
        if not local.exists():
            local.write_bytes(data)
        elif local.read_bytes() != data:
            raise ValueError('Reference coordinate-copy fingerprint collision')
    proc = subprocess.run([str(executable)], input=' '.join(map(str, args)) + '\n' + coordinate + '\n',
                          text=True, capture_output=True, timeout=timeout, check=True, cwd=executable.parent)
    prefix, body = proc.stdout.split('CANONICAL_RESULT\n', 1)
    rows = iter(body.splitlines())
    n, nw, converged = next(rows).split()
    values = list(map(float, next(rows).split()))
    if len(values) != 11:
        raise ValueError('Invalid reference scalar record')
    cl, cm, cd, cdf, cdp, top, bottom, mach, reynolds, nct, ncb = values
    diagnostics = dict(zip(('RMSBL', 'RMXBL', 'RLX'), map(float, next(rows).split())))
    diagnostics['iterations'] = len(re.findall(r'^\s*\d+\s+rms:', prefix, re.M))
    flags = re.findall(r'^\s*CANONICAL_FLAGS\s+([TF])\s+([TF])\s+([TF])\s*$', prefix, re.M)
    if len(flags) != 1:
        raise ValueError('Missing or invalid reference arithmetic flags')
    diagnostics['arithmetic_flags'] = dict(zip(
        ('invalid', 'overflow', 'divide_by_zero'), (flag == 'T' for flag in flags[0])))
    result = dict(N=int(n), NW=int(nw), converged=converged == 'T' if case['viscous'] else True,
                  coefficients=dict(CL=cl, CM=cm, CDP_SURFACE=cdp))
    if case['viscous']:
        n1, n2, e1, e2, t1, t2 = map(int, next(rows).split())
        forced = [value == 'T' for value in next(rows).split()]
        result.update(NBL=[n1, n2], IBLTE=[e1, e2],
                      transition=dict(x=[top, bottom], index=[t1, t2], forced=forced),
                      operating=dict(mach=mach, reynolds=reynolds, ncrit=[nct, ncb]))
        result['coefficients'].update(CD=cd, CDF=cdf, CDP=cd - cdf)
    panel = [list(map(float, next(rows).split())) for _ in range(result['N'] + result['NW'])]
    if any(len(row) != 6 for row in panel):
        raise ValueError('Invalid reference panel record')
    result['panel'] = panel if case['viscous'] else [[row[0], row[1], row[4]] for row in panel]
    if case['viscous']:
        result['bl'] = []
        for count in result['NBL']:
            side = []
            for _ in range(count):
                fields = next(rows).split()
                if len(fields) != len(BL_FIELDS):
                    raise ValueError('Invalid reference boundary-layer record')
                side.append([int(fields[0]), *map(float, fields[1:])])
            result['bl'].append(side)
    if next(rows) != 'CANONICAL_END' or list(rows):
        raise ValueError('Incomplete or unexpected reference output')
    return dict(result=result, diagnostics=diagnostics, logs=prefix.splitlines())


def tolerance_group(path, case):
    if path[0] == 'coefficients':
        return 'drag' if path[1] in ('CD', 'CDF', 'CDP') else 'coefficient'
    if path[0] == 'panel':
        # Wake coordinates depend on the solution; surface coordinates do not.
        if path[2] < 2 and path[1] < case['surface_count']:
            return 'geometry'
        return 'distribution'
    if path[0] == 'bl':
        return 'thickness' if BL_FIELDS[path[3]] in ('THET', 'DSTR', 'MASS', 'TSTR', 'TAU', 'DIS') else 'distribution'
    return 'transition'


def compare(actual, expected, case, tolerances):
    """Strict schema, finite floats, exact flags/indices; abs + relative bounds."""
    failures = []
    worst = {}
    count = 0
    case = dict(case, surface_count=expected.get('N', 0))

    def visit(a, b, path):
        nonlocal count
        label = '.'.join(map(str, path))
        if isinstance(b, dict):
            if not isinstance(a, dict) or a.keys() != b.keys():
                failures.append(dict(path=label, reason='object keys differ'))
                return
            for key in b:
                visit(a[key], b[key], (*path, key))
        elif isinstance(b, list):
            if not isinstance(a, list) or len(a) != len(b):
                failures.append(dict(path=label, reason='array length/type differs'))
                return
            for i, (ai, bi) in enumerate(zip(a, b)):
                visit(ai, bi, (*path, i))
        elif isinstance(b, (bool, int, str)):
            if type(a) is not type(b) or a != b:
                failures.append(dict(path=label, actual=a, reference=b, reason='discrete value differs'))
        else:
            count += 1
            if (not isinstance(b, float) or not isinstance(a, (float, int))
                    or isinstance(a, bool) or not math.isfinite(a) or not math.isfinite(b)):
                failures.append(dict(path=label, actual=a, reference=b, reason='nonfinite/non-numeric value'))
                return
            group = tolerance_group(path, case)
            tol = tolerances[group]
            bound = tol['atol'] + tol['rtol'] * abs(b)
            error = abs(a - b)
            ratio = error / bound
            item = dict(path=label, actual=a, reference=b, absolute_error=error, bound=bound, ratio=ratio)
            if ratio > worst.get(group, {}).get('ratio', -1):
                worst[group] = item
            if ratio > 1:
                failures.append(item)

    visit(actual, expected, ())
    return dict(passed=not failures, numeric_values=count, failures=failures, worst=worst)


def proves_instability(check):
    # Nonfinite values shared by two failed solves are not proof of sensitivity.
    return any('absolute_error' in item or item.get('reason') in (
        'discrete value differs', 'array length/type differs', 'object keys differ')
        for item in check['failures'])


def diagnose_reference_sensitivity(executable, case, reference, tolerances, timeout=120,
                                   probe_directory=None):
    """Probe the original solver with adjacent binary64 inputs, retaining evidence.

    Change one nonzero operating input to its adjacent binary64 value. Zero is
    skipped because its adjacent subnormal generally underflows during setup.
    These runs use unchanged source and a fresh process, just like the baseline.
    """
    variants = []

    def probe(altered, **details):
        result = run_reference(executable, altered, timeout)
        check = compare(result['result'], reference['result'], case, tolerances)
        variants.append(dict(**details, output=result, comparison=check))

    for field in ('alpha', 'mach', 'reynolds', 'ncrit'):
        if case[field] == 0:
            continue
        for direction in (-math.inf, math.inf):
            value = math.nextafter(case[field], direction)
            altered = dict(case, **{field: value})
            probe(altered, field=field, original=case[field], perturbed=value)
    if not any(proves_instability(item['comparison']) for item in variants) and 'naca' in case['foil']:
        original = case.get('cterat', 0.3)
        for direction in (-math.inf, math.inf):
            value = math.nextafter(original, direction)
            probe(dict(case, cterat=value), field='cterat', original=original, perturbed=value)
    if (not any(proves_instability(item['comparison']) for item in variants)
            and 'coordinates' in case['foil'] and probe_directory is not None):
        source = (ROOT / case['foil']['coordinates']).read_text().splitlines(keepends=True)
        probe_directory.mkdir(parents=True, exist_ok=True)
        found = False
        # Preserve the two trailing-edge endpoints and any closure condition.
        # Move one nonzero coordinate of an interior node, never the whole foil.
        for row in range(2, len(source) - 1):
            fields = source[row].split()
            if len(fields) != 2:
                raise ValueError('Expected two coordinates per canonical input node')
            for column in (0, 1):
                original = float(fields[column])
                if original == 0:
                    continue
                for direction in (-math.inf, math.inf):
                    value = math.nextafter(original, direction)
                    changed, values = source.copy(), fields.copy()
                    values[column] = format(value, '.17g')
                    changed[row] = ' '.join(values) + '\n'
                    path = probe_directory / f'node-{row}-{column}-{direction > 0}.dat'
                    path.write_text(''.join(changed))
                    foil = dict(case['foil'], coordinates=str(path.resolve()),
                                sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                    probe(dict(case, foil=foil), field='coordinates', line=row + 1, column=column,
                          original=original, perturbed=value, input=foil)
                    if proves_instability(variants[-1]['comparison']):
                        found = True
                        break
                if found:
                    break
            if found:
                break
    return dict(proven=any(proves_instability(item['comparison']) for item in variants),
                perturbation='one input by one binary64 ULP',
                convergence_changed=any(item['output']['result']['converged'] !=
                                        reference['result']['converged'] for item in variants),
                finite_budget_failures=sum(
                    'absolute_error' in failure
                    for item in variants for failure in item['comparison']['failures']),
                variants=variants)


def fixed_inputs_check(actual, expected, case, tolerances):
    """Reference instability never excuses the wrong geometry or conditions."""
    def fixed(result):
        return dict(N=result['N'], NW=result['NW'], operating=result.get('operating', {}),
                    panel=[row[:2] for row in result['panel'][:result['N']]])
    return compare(fixed(actual), fixed(expected), case, tolerances)


def arithmetic_failed(reference):
    return any(reference['diagnostics']['arithmetic_flags'].values())


def classify(check, sensitivity, initial_check, fixed_check, *,
             invalid_reference=False, valid_initial_reference=True):
    if check['passed'] and not invalid_reference:
        return 'match'
    if (valid_initial_reference and initial_check and initial_check['passed']
            and fixed_check and fixed_check['passed']):
        if invalid_reference:
            return 'reference_invalid'
        if sensitivity and sensitivity['proven']:
            return 'reference_unstable'
    return 'mismatch'


class Worker:
    def __init__(self, node, backend, output, timeout=120):
        self.timeout = timeout
        self.failure_path = output / f'{backend}-failure.json'
        self.stderr = (output / f'{backend}-stderr.log').open('w')
        try:
            self.proc = subprocess.Popen([node, str(HERE / 'worker.mjs'), backend],
                                     stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                     stderr=self.stderr, text=True, bufsize=1, cwd=ROOT)
        except BaseException:
            self.stderr.close()
            raise
        self.responses = queue.Queue()
        self.reader = threading.Thread(target=self._read, daemon=True)
        self.reader.start()

    def _read(self):
        try:
            for line in self.proc.stdout:
                self.responses.put(line)
        finally:
            self.responses.put(None)

    def run(self, case):
        case = dict(case, foil=dict(case['foil']))
        if 'coordinates' in case['foil']:
            case['foil']['coordinates'] = str(ROOT / case['foil']['coordinates'])
        self.proc.stdin.write(json.dumps(case) + '\n')
        self.proc.stdin.flush()
        try:
            line = self.responses.get(timeout=self.timeout)
        except queue.Empty:
            raise TimeoutError(f'Website worker exceeded {self.timeout}s for {case["id"]}') from None
        if not line:
            raise RuntimeError('Website worker exited; see its stderr log')
        result = json.loads(line)
        if result.get('error'):
            save_json(self.failure_path, result)
            raise RuntimeError(result['error'])
        if result['id'] != case['id']:
            raise RuntimeError('Worker response case ID differs')
        return result

    def close(self):
        try:
            self.proc.stdin.close()
        except BrokenPipeError:
            pass
        try:
            self.proc.wait(timeout=2)
        except subprocess.TimeoutExpired:
            self.proc.terminate()
            try:
                self.proc.wait(timeout=2)
            except subprocess.TimeoutExpired:
                self.proc.kill()
                self.proc.wait()
        self.reader.join(timeout=2)
        self.proc.stdout.close()
        self.stderr.close()


def json_safe(value):
    """Keep failed nonfinite results readable in strict JSON artifacts."""
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def save_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(json_safe(value), indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def input_fingerprints():
    paths = [*sorted((ROOT / 'js').glob('*.js')), *sorted(HERE.glob('*.py')),
             HERE / 'worker.mjs', HERE / 'reference_driver.f', HERE / 'cases.json',
             ROOT / 'python/tests/fortran_reference.py']
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='25-case diagnostic subset; not the canonical full gate')
    parser.add_argument('--case', action='append', default=[], help='Case ID glob for diagnosis')
    parser.add_argument('--backend', choices=['both', 'auto', 'javascript'], default='both')
    parser.add_argument('--compiler', default='gfortran', help='Fortran compiler executable')
    parser.add_argument('--node', default='node')
    parser.add_argument('--timeout', type=float, default=120, help='Seconds per solver operating point')
    parser.add_argument('--diagnose-roundoff', action='store_true',
                        help='Retain sensitivity diagnostics even with --strict-reference')
    parser.add_argument('--strict-reference', action='store_true',
                        help='Fail on unstable or invalid references too (default: report them separately)')
    parser.add_argument('--output', type=pathlib.Path)
    args = parser.parse_args(argv)
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error('--timeout must be finite and positive')
    for command in (args.compiler, args.node):
        if not shutil.which(command):
            parser.error(f'{command} is required; the canonical gate never skips missing toolchains')
    manifest = load_manifest()
    selected = cases(manifest, args.quick)
    if args.case:
        selected = [case for case in selected if any(fnmatch.fnmatchcase(case['id'], pattern) for pattern in args.case)]
    if not selected:
        parser.error('No cases selected')
    timestamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    output = (args.output or ROOT / 'artifacts' / 'canonical' / timestamp).resolve()
    if output.exists() and any(output.iterdir()):
        parser.error('Output directory must be empty to preserve earlier results')
    output.mkdir(parents=True, exist_ok=True)
    build = output / 'build'
    build.mkdir()
    metadata = dict(schema_version=1, manifest=manifest, cases=selected,
                    compiler=subprocess.check_output([args.compiler, '--version'], text=True, timeout=15),
                    node=subprocess.check_output([args.node, '--version'], text=True, timeout=15).strip(),
                    driver_sha256=hashlib.sha256((HERE / 'reference_driver.f').read_bytes()).hexdigest(),
                    flags=list(REFERENCE_FLAGS), input_sha256=input_fingerprints(),
                    reference_adjustments=dict(viscal_residual_tolerance=manifest['residual_tolerance'],
                                               upstream_residual_tolerance=1e-4),
                    command=argv if argv is not None else sys.argv[1:],
                    panel_fields=PANEL_FIELDS, inviscid_panel_fields=['X', 'Y', 'CPI'], bl_fields=BL_FIELDS)
    save_json(output / 'metadata.json', metadata)
    print(f'Compiling original XFOIL; {len(selected)} operating points; output: {output}', flush=True)
    backends = ['auto', 'javascript'] if args.backend == 'both' else [args.backend]
    workers = {}
    summary = dict(passed=False, completed=False, status='running', cases=len(selected),
                   completed_cases=0, comparisons=0, expected_comparisons=len(selected) * len(backends), failed_cases=[],
                   unconverged_reference=[], worst={}, numeric_values=0, backends=backends,
                   matched_comparisons=0, reference_unstable_comparisons=0, reference_unstable_cases=[],
                   reference_invalid_comparisons=0, reference_invalid_cases=[],
                   initial_comparisons=0, strict_reference=args.strict_reference,
                   full_gate=not args.quick and not args.case and args.backend == 'both')
    start = time.monotonic()
    save_json(output / 'summary.json', summary)
    try:
        executable = build_viscal_reference(build, driver=HERE / 'reference_driver.f',
                                            load_coordinates=True, compiler=args.compiler,
                                            residual_tolerance=manifest['residual_tolerance'])
        for backend in backends:
            workers[backend] = Worker(args.node, backend, output, args.timeout)
        with gzip.open(output / 'results.jsonl.gz', 'wt', encoding='utf8') as records, (output / 'polar.csv').open('w', newline='') as polar:
            columns = ['case', 'implementation', 'airfoil', 'alpha', 'mach', 'reynolds', 'ncrit', 'viscous',
                       'converged', 'CL', 'CM', 'CD', 'CDF', 'CDP', 'CDP_SURFACE', 'XTRT', 'XTRB']
            writer = csv.DictWriter(polar, fieldnames=columns)
            writer.writeheader()
            for index, case in enumerate(selected):
                summary['current_case'] = case['id']
                reference = run_reference(executable, case, args.timeout)
                if not reference['result']['converged']:
                    summary['unconverged_reference'].append(case['id'])
                outputs = {'fortran': reference}
                checks = {}
                for backend, worker in workers.items():
                    actual = worker.run(case)
                    check = compare(actual['result'], reference['result'], case, manifest['tolerances'])
                    checks[backend] = check
                    outputs[backend] = actual
                    summary['comparisons'] += 1
                    summary['numeric_values'] += check['numeric_values']
                    for group, item in check['worst'].items():
                        if item['ratio'] > summary['worst'].get(group, {}).get('ratio', -1):
                            summary['worst'][group] = dict(item, case=case['id'], backend=backend)
                record = dict(case=case, outputs=outputs, checks=checks)
                invalid_reference = not reference['result']['converged'] and arithmetic_failed(reference)
                if (not invalid_reference and (not args.strict_reference or args.diagnose_roundoff)
                        and any(not check['passed'] for check in checks.values())):
                    record['reference_sensitivity'] = diagnose_reference_sensitivity(
                        executable, case, reference, manifest['tolerances'], args.timeout,
                        output / 'probes' / case['id'])
                if (invalid_reference or record.get('reference_sensitivity', {}).get('proven')) and case['viscous']:
                    initial_case = dict(case, id=case['id'] + '-first-iteration', iterations=1)
                    initial_reference = run_reference(executable, initial_case, args.timeout)
                    initial = dict(case=initial_case, outputs={'fortran': initial_reference}, checks={},
                                   valid_reference=not arithmetic_failed(initial_reference))
                    for backend, worker in workers.items():
                        initial_actual = worker.run(initial_case)
                        initial['outputs'][backend] = initial_actual
                        initial['checks'][backend] = compare(initial_actual['result'],
                            initial_reference['result'], initial_case, manifest['tolerances'])
                        summary['initial_comparisons'] += 1
                    record['first_iteration'] = initial
                record['classifications'] = {}
                record['fixed_input_checks'] = {}
                for backend, check in checks.items():
                    invariant = fixed_inputs_check(outputs[backend]['result'], reference['result'],
                                                   case, manifest['tolerances'])
                    record['fixed_input_checks'][backend] = invariant
                    initial_check = record.get('first_iteration', {}).get('checks', {}).get(backend)
                    category = classify(check, record.get('reference_sensitivity'), initial_check, invariant,
                                        invalid_reference=invalid_reference,
                                        valid_initial_reference=record.get('first_iteration', {}).get('valid_reference', False))
                    record['classifications'][backend] = category
                    if category == 'match':
                        summary['matched_comparisons'] += 1
                    elif category == 'reference_unstable':
                        summary['reference_unstable_comparisons'] += 1
                    elif category == 'reference_invalid':
                        summary['reference_invalid_comparisons'] += 1
                    if category == 'mismatch' or (args.strict_reference and category != 'match'):
                        summary['failed_cases'].append(dict(case=case['id'], backend=backend,
                                                           failures=len(check['failures'])))
                if 'reference_unstable' in record['classifications'].values():
                    summary['reference_unstable_cases'].append(case['id'])
                if 'reference_invalid' in record['classifications'].values():
                    summary['reference_invalid_cases'].append(case['id'])
                records.write(json.dumps(json_safe(record), allow_nan=False) + '\n')
                records.flush()
                for implementation, value in outputs.items():
                    r = value['result']
                    row = dict(case=case['id'], implementation=implementation, airfoil=case['foil']['id'],
                               **{key: case[key] for key in ('alpha', 'mach', 'reynolds', 'ncrit', 'viscous')},
                               converged=r['converged'], **r['coefficients'])
                    if 'transition' in r:
                        row.update(XTRT=r['transition']['x'][0], XTRB=r['transition']['x'][1])
                    writer.writerow(row)
                polar.flush()
                categories = record['classifications'].values()
                label = ('FAIL' if 'mismatch' in categories else
                         'REFERENCE_INVALID' if 'reference_invalid' in categories else
                         'REFERENCE_UNSTABLE' if 'reference_unstable' in categories else 'PASS')
                if label in ('REFERENCE_UNSTABLE', 'REFERENCE_INVALID') and args.strict_reference:
                    label += ' (strict failure)'
                print(f'[{index + 1}/{len(selected)}] {case["id"]}: {label}', flush=True)
                summary['elapsed_seconds'] = time.monotonic() - start
                summary['completed_cases'] += 1
                save_json(output / 'summary.json', summary)
        summary['completed'] = True
        summary['status'] = 'completed'
        summary['passed'] = (not summary['failed_cases']
                             and summary['comparisons'] == summary['expected_comparisons'])
    except (Exception, KeyboardInterrupt) as error:
        summary['status'] = 'error'
        summary['error'] = f'{type(error).__name__}: {error}'
        (output / 'error.log').write_text(traceback.format_exc())
        if isinstance(error, (subprocess.CalledProcessError, subprocess.TimeoutExpired)):
            for name in ('stdout', 'stderr'):
                value = getattr(error, name, None)
                if value:
                    (output / f'failed-process-{name}.log').write_text(
                        value.decode(errors='replace') if isinstance(value, bytes) else value)
        print(summary['error'], file=sys.stderr)
    finally:
        for worker in workers.values():
            worker.close()
        summary['elapsed_seconds'] = time.monotonic() - start
        save_json(output / 'summary.json', summary)
    print(f'{"PASS" if summary["passed"] else "FAIL"}: {summary["comparisons"]} comparisons, '
          f'{summary["matched_comparisons"]} match, '
          f'{summary["reference_unstable_comparisons"]} reference unstable, '
          f'{summary["reference_invalid_comparisons"]} reference invalid, '
          f'{len(summary["failed_cases"])} failures; see {output / "summary.json"}')
    return 0 if summary['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
