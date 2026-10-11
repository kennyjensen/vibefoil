# Python Tests

All Python tests in this repository use the standard library `unittest` framework.

Run the full suite from the repo root:

```sh
python3 -m unittest discover -s python/tests -p "test_*.py"
```

Run a single test module:

```sh
python3 -m unittest python.tests.test_compare_blsolv_strict
```

## XFOIL numerical audit

The numerical reference is the XFOIL 6.996 source checked into
`third_party/Xfoil` by commit `32ff3c2`, including the September 2020 and January
2023 addenda. The tests compile that source directly, or extract complete
subroutines with `fortran_reference.py`. They no longer use copied solver
snapshots. Plotting and interactive input are stubbed where necessary; numerical
solver routines under comparison are not replaced. The isolated BLU tests retain
a minimal BLPINI setup that initializes CFFAC; the complete operating-point
reference uses upstream initialization.

The Fortran reference uses `-fdefault-real-8 -fdefault-double-8
-ffp-contract=off` to match the ports' binary64 arithmetic. Comparing a default
single-precision XFOIL executable would measure both precision and algorithm
differences. The C++ build also disables contraction of multiply/add operations.

### Corrections covered by the audit

| Area | Corrections |
| --- | --- |
| Updated upstream algorithms | Skin-friction limiter derivatives; optional sharp NACA 4/5-digit sections; leading-edge search; segmented spline endpoint constraints. |
| Updated defaults and capacity | 240 panels, trailing-edge density ratio 0.3, 200 NACA points per side, IQX 601, spline workspace 2500, LU workspace 515. Explicit caller settings still take precedence. |
| Boundary-layer diagnostics | Density-defect thickness, mass/momentum/kinetic-energy defects, wake copying, and friction lift coefficient CLF. |
| Panel influence | Python/C++ skip the zero-length trailing-edge base panel; reflected-panel self-influence is retained; C++ reflected vortex signs and the missing reflected base contribution are corrected. JavaScript uses spline chord length for sharp-edge classification and wake length. |
| Coupled solver | Python accepts vector right-hand sides in production Gaussian elimination; JavaScript accepts typed-array matrix right-hand sides, saves updated shear/diagnostics during SETBL, applies the wake angle derivative only in the appropriate solve mode, and refreshes mass defects when the stagnation index stays unchanged. Mach/Reynolds limiting and cached operating values follow MRCL. |
| C++ build and operating setup | Include inverse-design code in CMake, repair its field names, initialize the operating-point lift before viscous iteration, and match numerical state defaults. |
| Test fixtures | Repair the CFT call signature, initialize panel angles and boundary-layer settings consistently, remove an uninitialized diagnostic from driver output, and exercise the production linear solver without test monkeypatches. |

The obsolete unconditional `1000.bl` diagnostic write was also removed from
Python/C++, matching the upstream disabled block. Artificial zero-pivot
replacement was removed from LU factorization to match the reference algorithm;
singular inputs remain outside the valid solver domain.

### Coverage and tolerances

`test_upstream_numerical_changes.py` checks Python, JavaScript, native C++, and
actual WebAssembly against Fortran. It covers limiter branches, standard and
sharp NACA 4/5-digit shapes, segmented splines including 1201 points, leading-edge
search, Mach/Reynolds scaling and clipping, rotated trailing edges, and all panel
influence derivatives for sharp/blunt and reflected/unreflected geometries.
These checks use `2e-12 * max(1, abs(reference))`.

`test_upstream_viscal.py` runs complete operating points through the upstream
NACA, PANGEN, SPECAL, and VISCAL routines:

| Airfoil | Mach | Reynolds | Alpha | Panels | Iterations |
| --- | --- | --- | --- | --- | --- |
| NACA 0012 | 0 | 1,000,000 | 0 degrees | 240 | 1 |
| NACA 2412 | 0.1 | 3,000,000 | 4 degrees | 120 | 20 |
| NACA 23012 | 0.2 | 5,000,000 | 2 degrees | 120 | 20 |
| Sharp NACA 0012 | 0 | 1,000,000 | 3 degrees | 120 | 1 |

It compares coefficients, residuals, transition indices, convergence flags,
surface and wake coordinates/velocities/Cp, and boundary-layer distributions. Python/WASM
use scaled bounds of `1e-9` for scalars and `1e-8` for distributions; JavaScript
uses `1e-8` and `1e-7`. Discrete flags and indices must match exactly. These are
regression tolerances for the listed cases, not global error bounds.

V8 and the native math library round `log` and `atan2` differently. Replaying
Python's elementary-function results in JavaScript made the NACA 2412 test's
influence matrices bit-identical; native V8 results differed by up to about
`1.4e-11` after factorization. Cancellation and Newton iteration amplify that
roundoff. The strict SETBL test therefore shares the source-influence matrix to
isolate boundary-layer assembly, while the complete operating-point tests retain
each port's independent panel calculation. Bit-for-bit identity across runtimes
is not claimed.

The existing suite additionally checks panel/wake construction, boundary-layer
closures and marching, matrix assembly/solution, stagnation relocation, viscous
reuse, pressure smoothing, geometry helpers, coefficient integration, and dump
formats. The inverse-design test compares full MIXED output, including all node
coordinates and vorticity, against Fortran in all four ports, using a scaled
`1e-7` bound. NACA 6-series generation is an extension tested against the separate
PDAS `third_party/naca456` source.

### Reproducing the complete checks

Validated on 2026-10-06: all 76 tests passed with no skips (447 seconds), the
Release WebAssembly build succeeded, and both `npm test` suites passed. The run
used Python 3.14.4, Node 24.16.0, GNU Fortran 15.2.0, Clang 19.1.7, and Emscripten
3.1.69. Changed JavaScript files also passed `node --check`.

Install Python 3, Node, gfortran, a C++17 compiler exposed as `c++`, CMake, and
Emscripten with `em++` on PATH. Build the WASM operating-point driver before
running the suite:

```sh
emcmake cmake -S wasm -B wasm/build-wasm -DCMAKE_BUILD_TYPE=Release
cmake --build wasm/build-wasm -j2
python3 -m unittest discover -s python/tests -p "test_*.py"
npm test
```

The focused kernel and inverse-design tests compile their native and WASM
drivers themselves. Tests requiring unavailable toolchains are skipped; a run
with skips does not establish coverage of those ports. The generated WASM build
is ignored by Git and must be rebuilt after changing C++ sources.

### Scope of equivalence

This is an audit of the implemented numerical routines and the exercised
operating paths. It is not a proof for every geometry, transition event,
separated-flow case, or compiler/runtime, nor a complete port of the XFOIL CLI.
The ports do not expose all upstream geometry commands (including NORM, XLFIND,
CANG, SQTE, TECIRC and the new TGAP/SYTB controls), full inverse MDES, plotting,
SVG/DXF output, polar-file management, or BLDUMP2/BLDUMP3. The JavaScript SPECAL
entry point is a fixed-Mach velocity-superposition wrapper, not the upstream
general CL-dependent Mach iteration. The application uses the fixed-Mach path;
MRCL's scaling kernels are tested separately.

Passing these checks supports numerical agreement for the covered paths.
It does not establish complete feature equivalence with original XFOIL.
