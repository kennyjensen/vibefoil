"""Compile test oracles from the checked-in upstream source, never stale copies."""

import pathlib
import math
import re
import subprocess

SOURCE_DIR = pathlib.Path(__file__).resolve().parents[2] / "third_party" / "Xfoil" / "src"
REFERENCE_FLAGS = ('-O2', '-fdefault-real-8', '-fdefault-double-8',
                   '-ffp-contract=off', '-ffunction-sections', '-fdata-sections',
                   '-Wl,--gc-sections')
_START = re.compile(r"^      \s*(?:(?:REAL|COMPLEX|INTEGER|LOGICAL)\s+)?(?:SUBROUTINE|FUNCTION)\s+(\w+)", re.I)
_END = re.compile(r"^      \s*END\s*(?:!.*)?$", re.I)


def extract_routines(source, names):
    """Return complete fixed-form program units without rewriting their math."""
    units = {}
    name = None
    lines = []
    for line in source.splitlines(keepends=True):
        match = _START.match(line)
        if match:
            if name is not None:
                raise ValueError(f"Unterminated Fortran routine {name}")
            name = match[1].upper()
            lines = []
        if name is not None:
            lines.append(line)
            if _END.match(line.rstrip('\n')):
                units[name] = ''.join(lines)
                name = None
    return '\n'.join(units[name.upper()] for name in names)


def reference_source(output_dir, source_name, names):
    source = (SOURCE_DIR / source_name).read_text()
    output = pathlib.Path(output_dir) / (pathlib.Path(source_name).stem + '_' + '_'.join(names) + '.f')
    output.write_text(extract_routines(source, names))
    return output


def set_residual_tolerance(source, tolerance):
    """Change only VISCAL's stopping threshold in a generated build copy."""
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError('Residual tolerance must be finite and positive')
    original = '      DATA EPS1 / 1.0E-4 /'
    if source.count(original) != 1:
        raise ValueError('Expected exactly one original VISCAL stopping threshold')
    return source.replace(original, f'      DATA EPS1 / {tolerance:.17E} /')


def build_viscal_reference(output_dir, *, driver=None, load_coordinates=False,
                           compiler='gfortran', residual_tolerance=None):
    """Build the current upstream solver without its interactive/plot UI."""
    sources = [pathlib.Path(driver) if driver else
               pathlib.Path(__file__).with_name('viscal_reference_driver.f')]
    sources += [SOURCE_DIR / name for name in (
        'naca.f', 'spline.f', 'xgeom.f', 'xpanel.f', 'xbl.f', 'xblsys.f',
        'xsolve.f', 'xutils.f', 'userio.f')]
    routines_by_file = {
        'xfoil.f': ('INIT', 'MRCL', 'COMSET', 'CPCALC', 'CLCALC', 'CDCALC',
                    'NACA', 'PANGEN', 'TECALC'),
        'xoper.f': ('SPECAL', 'VISCAL', 'MHINGE'),
        'xgdes.f': ('GETXYF', 'INSIDE'),
    }
    if load_coordinates:
        sources.append(SOURCE_DIR / 'aread.f')
        routines_by_file['xfoil.f'] += ('LOAD',)
        routines_by_file['xgdes.f'] += ('ABCOPY',)
    for filename, routines in routines_by_file.items():
        source = reference_source(output_dir, filename, routines)
        if filename == 'xoper.f' and residual_tolerance is not None:
            source.write_text(set_residual_tolerance(source.read_text(), residual_tolerance))
        sources.append(source)
    target = pathlib.Path(output_dir) / 'viscal_reference'
    result = subprocess.run([
        compiler, *REFERENCE_FLAGS, '-I', str(SOURCE_DIR), *map(str, sources),
        '-o', str(target)], capture_output=True, text=True, timeout=300)
    (pathlib.Path(output_dir) / 'compiler.log').write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    return target


def run_viscal_reference(executable, ides, minf, reinf, alpha_deg, *,
                         npan=240, cterat=0.3, waklen=1.0, iterations=1, sharp=False):
    args = [ides, npan, cterat, minf, reinf, alpha_deg, waklen, iterations, int(sharp)]
    result = subprocess.run([str(executable)], input=' '.join(map(str, args)) + '\n',
                            capture_output=True, text=True, check=True)
    rows = iter(result.stdout.split('REFERENCE_RESULT\n', 1)[1].splitlines())
    out = dict(zip(('ALFA_DEG', 'CL', 'CM', 'CD', 'CDF', 'CDP_SURFACE', 'CLF', 'RMSBL', 'RMXBL'),
                   map(float, next(rows).split())))
    out['RLX'], top, bottom = map(float, next(rows).split())
    out['XOCTR'] = [0.0, top, bottom]
    itop, ibottom, out['IMXBL'], out['ISMXBL'], nbl1, nbl2 = map(int, next(rows).split())
    out['ITRAN'] = [0, itop, ibottom]
    top, bottom, converged, out['VMXBL'] = next(rows).split()
    out['TFORCE'] = [False, top == 'T', bottom == 'T']
    out['LVCONV'] = converged == 'T'
    out['CDP'] = out['CD'] - out['CDF']
    out['N'], out['NW'] = map(int, next(rows).split())
    out['panel'] = [list(map(float, next(rows).split())) for _ in range(out['N'] + out['NW'])]
    out['bl'] = [[list(map(float, next(rows).split())) for _ in range(nbl)] for nbl in [nbl1, nbl2]]
    return out
