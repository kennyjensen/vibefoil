"""Refresh MHINGE fixtures using unmodified, checked-in original XFOIL."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from python.tests.fortran_reference import SOURCE_DIR, REFERENCE_FLAGS, reference_source

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--compiler', default='gfortran')
args = parser.parse_args()
reference = ROOT / 'tests/hinge_moment_reference.json'
data = json.loads(reference.read_text())
script = '''
import { readFileSync } from 'node:fs';
import { makeHingeFixture } from './tests/hinge_moment_fixture.js';
const cases = JSON.parse(readFileSync('tests/hinge_moment_reference.json')).cases;
console.log(JSON.stringify(cases.map(c => {
  const {panel:p, cp, hinge:h} = makeHingeFixture(c);
  return [p.N, h.x, h.y, +c.viscous].join(' ') + '\\n'
    + Array.from({length:p.N}, (_,i) => [p.X[i],p.Y[i],p.XP[i],p.YP[i],p.S[i],cp[i+1]].join(' ')).join('\\n') + '\\n';
})));
'''
inputs = json.loads(subprocess.check_output(['node', '--input-type=module', '-e', script], cwd=ROOT, text=True))
with tempfile.TemporaryDirectory(prefix='hinge-reference-') as directory:
    work = Path(directory)
    driver = work / 'driver.f'
    driver.write_text('''      PROGRAM HINGE_REFERENCE
      INCLUDE 'XFOIL.INC'
      INTEGER IVISC
      REAL PRESSURE
      READ(*,*) N,XOF,YOF,IVISC
      LFLAP = .TRUE.
      LVISC = IVISC.EQ.1
      DO I=1,N
        READ(*,*) X(I),Y(I),XP(I),YP(I),S(I),PRESSURE
        CPI(I) = -99.
        CPV(I) = -99.
        IF(LVISC) THEN
          CPV(I) = PRESSURE
        ELSE
          CPI(I) = PRESSURE
        ENDIF
      ENDDO
      CALL MHINGE
      WRITE(*,'(A,ES26.17E3)') 'HINGE_RESULT ',HMOM
      END
''')
    sources = [driver, reference_source(work, 'xoper.f', ('MHINGE',)),
               reference_source(work, 'xgdes.f', ('GETXYF',)),
               SOURCE_DIR / 'spline.f', SOURCE_DIR / 'userio.f']
    executable = work / 'hinge_reference'
    subprocess.run([args.compiler, *REFERENCE_FLAGS, '-I', str(SOURCE_DIR),
                    *map(str, sources), '-o', str(executable)], check=True)
    for case, source in zip(data['cases'], inputs):
        output = subprocess.check_output([str(executable)], input=source, text=True)
        case['ch'] = float(output.split('HINGE_RESULT ')[1].split()[0])
data['sourceSha256'] = hashlib.sha256((SOURCE_DIR / 'xoper.f').read_bytes()).hexdigest()
reference.write_text(json.dumps(data, indent=2) + '\n')
print(f'Updated {len(data["cases"])} original-XFOIL hinge-moment fixtures')
